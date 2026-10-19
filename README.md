# Figure 3: pooled myofiber RNA-seq revision analysis

This package contains the inputs, analysis code, frozen reference results, and figure-generation code for the revised pooled-myofiber RNA-seq analysis in DMD rats.

Open **`Fig3_pooled_myofiber_RNAseq_revision.ipynb`** for the annotated workflow. The corresponding **`Fig3_pooled_myofiber_RNAseq_revision.py`** provides a script entry point.

## Experimental and statistical units

RNA-seq libraries were prepared from pooled isolated myofibers: 10 GFP-positive or 30 GFP-negative fibers per flexor digitorum brevis (FDB) muscle. The analysis includes all 17 available libraries, comprising eight GFP-negative and nine GFP-positive samples from five rats.

Available left- and right-side samples are combined within each rat and GFP group after the input conversion described below. The resulting paired comparison uses **five independent rats**, with one GFP-positive and one GFP-negative count vector per rat. Multiple FDB samples from the same rat are not treated as independent animals.

## Analysis

### Input conversion and differential expression

The source matrices are gene-level Salmon estimated counts and TPM values exported from RaNA-seq. The original Salmon commands are retained in `data/salmon_commands.txt`. This workflow starts from those gene-level matrices; it does not rerun Salmon on sequencing reads.

For each library, the Python code applies the scaledTPM formula:

```text
scaled_count[g, s] = TPM[g, s] / sum_g(TPM[g, s]) * sum_g(estimated_count[g, s])
```

The conversion is performed separately for each library, before summing available left- and right-side values within each rat and GFP group. Values are rounded to integers once, after summation. TPM values are therefore not supplied directly to PyDESeq2. The tximport package is not invoked, and no transcript-length offset is added after this conversion.

- Differential expression: **PyDESeq2 0.5.4**, design `~ animal + condition`, parametric dispersion fit, and a Wald comparison of GFP-positive versus GFP-negative samples.
- Gene filter: total aggregated count of at least 10 across the ten rat/group count vectors.
- DEG threshold: Benjamini–Hochberg (BH)-adjusted P < 0.05. No fold-change threshold defines significance.
- PCA: variance-stabilizing transformation without using the design, followed by PCA of the 500 most variable genes.

### Preranked GSEA

Genes with finite Wald statistics are ranked by the signed statistic. Positive values indicate higher expression in GFP-positive samples. The ranked list is not restricted to significant DEGs.

Frozen rat Reactome and KEGG gene sets are supplied in `references/`. Sets must contain 15–500 genes represented in the ranked list. GSEA uses **GSEApy 1.3.1**, the multilevel preranked method, `weight=1`, `sample_size=1001`, `eps=1e-10`, `seed=20260915`, and `permutation_num=1000`. The configured thread count is four.

Nominal enrichment P values are BH-adjusted jointly across **all 1,310 tested gene sets from both databases**. Gene-set significance is defined by `BH_padj_all_pathways < 0.05`.

### Frozen reference results and figure display

| Item | Reference result |
|---|---:|
| Libraries included | 17 |
| Paired rats | 5 |
| Retained genes | 23,015 |
| DEGs at BH-adjusted P < 0.05 | 2 |
| Tested gene sets | 1,310 |
| Significant gene sets at BH-adjusted P < 0.05 | 49 |
| Selected significant gene sets displayed in Fig. 3i | 10 |

These values describe the frozen reference analysis; they are not targets used to select analysis parameters.

Fig. 3i displays an illustrative subset of ten significant gene sets. Its display filter is applied **after** correction across all tested sets. The complete tested and significant results, together with a plot of all 49 significant sets, are also retained.

PCA ellipses are descriptive data contours, not confidence intervals or tests of group separation. The volcano plot labels both DEGs; its vertical lines at log2 fold changes of ±1 are visual guides. A full-range volcano plot is supplied alongside the main display.

## Package contents

Keep the following files and directories together:

```text
Fig3_pooled_myofiber_RNAseq_revision.ipynb
Fig3_pooled_myofiber_RNAseq_revision.py
README.md
requirements.txt
analysis_config.json
manifest.json
data/
  counts_all_samples.csv.gz
  tpm_all_samples.csv.gz
  sample_metadata.csv
  salmon_commands.txt
references/
  gene_annotation.csv
  rat_gene_sets.json
  rat_pathway_names.json
  sources.json
snapshot/
  DE_annotated.csv
  GSEA_all_pathways.csv
  VST_counts.csv.gz
  normalized_counts.csv.gz
  gene_diagnostics.csv
audit/
comparison/
results_reanalyze/
validation.json
release_manifest.json
```

`manifest.json` records SHA256 hashes for the analysis inputs, references, and frozen numerical results. The code checks these hashes before proceeding. Both execution modes require the complete package, including `data/`, `references/`, and `snapshot/`; the count matrix and metadata alone are insufficient.

Original sample identifiers are retained. `PoorQC` is a legacy sample-name label and does not determine inclusion; all 17 libraries are used.

`results_reanalyze/` contains the saved primary result tables and figures. `analysis_config.json` records the fixed analysis settings. `audit/` contains input-provenance and numerical-validation records. `comparison/` retains detailed numerical comparisons with the earlier counts-input version of the paired analysis.

Once dependencies are installed, the analysis uses the supplied files and does not query external pathway databases.

## Running the workflow

Use Python 3.12 and install the packages listed in `requirements.txt`. Start from the directory containing the notebook, `manifest.json`, and the package subdirectories.

For a virtual environment on macOS or Linux:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
jupyter lab
```

Open `Fig3_pooled_myofiber_RNAseq_revision.ipynb`, select the virtual environment's Python kernel, set `MODE` in the setup cell, and run the cells in order.

| Mode | Behavior | Output directory |
|---|---|---|
| `MODE = "reanalyze"` (default) | Repeats input conversion, paired differential expression, GSEA, PCA, and figure generation from the supplied counts and TPM values. | `results_reanalyze/` |
| `MODE = "replot"` | Loads the frozen DE, GSEA, and VST results and regenerates figures and output tables without rerunning the primary statistical inference. | `results_replot/` |

To use the script, set `MODE` near its top and run it from the same package directory:

```bash
python Fig3_pooled_myofiber_RNAseq_revision.py
```

Outputs include DE tables, significant genes, all tested GSEA results, significant gene sets, PCA coordinates, figure data, software versions, and `run_summary.json`. Figures are saved as PNG, SVG, and PDF. The summary reports numerical differences from the frozen reference results. PCA component signs may reverse across numerical environments without changing distances among samples.

## Optional numerical-convergence audit

The final code cell contains `RUN_CONVERGENCE_SENSITIVITY = False`. Setting it to `True` runs a separate GSEA calculation after removing genes flagged for non-convergence in any recorded PyDESeq2 fitting step. This optional calculation runs when enabled in either execution mode; it does not replace the primary analysis.

In the stored audit, removing 169 flagged genes yielded 56 significant gene sets, including **46 of the 49** significant sets from the primary analysis. Neither primary DEG was among the flagged genes. Audit results are retained in `audit/`.

This check examines sensitivity to numerical fitting flags. It is not a test of biological robustness or of sensitivity to excluding individual libraries.

## Reproducibility and interpretation

Reproducing stored numbers demonstrates numerical reproducibility of this fixed calculation. It does not establish independent biological validation, rule out sample-dependent or technical effects, or assign the RNA signals exclusively to myonuclei. GSEA results are exploratory associations and do not establish pathway activation, cellular senescence, or a causal effect of p16.

Packaging checks are recorded in `validation.json`. For this English-language release, all 12 code cells were executed sequentially in Python 3.12 in `replot` mode. A Jupyter kernel was not used, and differential-expression and GSEA inference were not rerun. The statistical code, inputs, and saved primary results are unchanged; the original numerical-validation records are retained.
