# MRBIGR2 Tutorial

**Project:** `MRBIGR2`  
**Repository:** `https://github.com/CrazyHsu/MRBIGR2`  
**Original MRBIGR Repository:** `https://gitee.com/crazyhsu/MRBIGR`  
**Original MRBIGR Tutorial:** `https://mrbigr.github.io/`  
**Version:** `1.0.0`  
**Author:** `CrazyHsu <crazyhsu9527@gmail.com>`

---

## Contents

1. [Module Overview](#1-module-overview)
2. [Genotype Analysis](#2-genotype-analysis)
3. [Phenotype Analysis](#3-phenotype-analysis)
4. [GWAS Analysis](#4-gwas-analysis)
5. [Annotation and QTL Analysis](#5-annotation-and-qtl-analysis)
6. [Peak and Haplotype Analysis](#6-peak-and-haplotype-analysis)
7. [Mendelian Randomization](#7-mendelian-randomization)
8. [GO, KEGG, and GSEA](#8-go-kegg-and-gsea)
9. [Network Analysis](#9-network-analysis)
10. [Visualization](#10-visualization)
11. [Utility Helpers](#11-utility-helpers)
12. [Quick Start Checklist](#12-quick-start-checklist)

---

## 1. Module Overview

MRBIGR2 is a pure Python and agent-oriented reorganization of the MRBIGR workflow. It keeps the core analysis logic in Python and provides both CLI and MCP access.

Main modules:

| Module | Role |
| --- | --- |
| `geno.py` | genotype conversion, QC, subset, PCA, kinship, IBD, pruning, clumping |
| `pheno.py` | filtering, scaling, imputation, outlier handling, BLUP, BLUE |
| `gwas.py` | LM, LMM, PLINK GWAS, SNP ranking, lambda, clumping preparation |
| `anno.py` | GTF parsing, SNP annotation, QTL annotation, variant effect prediction |
| `qtl.py` | QTL region detection, lead SNP selection, QTL summaries, candidate gene mapping |
| `peak.py` | haplotype-oriented testing and local region plotting |
| `mr.py` | MR estimation, heterogeneity, pleiotropy, QTL target analysis |
| `go.py` | GO enrichment, GSEA, KEGG enrichment, plotting, simplification, report export |
| `net.py` | edge weighting, module detection, hub identification, network metrics |
| `vis.py` | Manhattan, QQ, PCA, phenotype and LD plots |
| `parallel.py` | parallel execution and utility helpers (internal, not MCP-exposed) |

## 2. Genotype Analysis

Typical genotype workflows start from PLINK-format data.

### Common tasks

- convert VCF or HapMap to PLINK
- run SNP QC
- build subsets for fast testing
- calculate PCA
- calculate kinship or IBD
- perform pruning or clumping

### Python examples

```python
from mrbigr.core import geno

geno.snp_qc("data/geno", "output/geno_qc", maf=0.05, missing_rate=0.2, mind=0.2)
pc_df, variance_ratio = geno.calculate_pca("output/geno_qc", n_components=10)
stats = geno.get_snp_stats("output/geno_qc")
```

### Notes

- For fast testing, create a random SNP subset before running expensive analyses.
- PCA and kinship are often required downstream for GWAS and trait correction.

## 3. Phenotype Analysis

Phenotype workflows typically include filtering, scaling, missing-value handling, and mixed-model summaries.

### Common tasks

- filter low-abundance or high-missing entries
- apply log or Z-score transformations
- impute missing values
- detect outliers
- run BLUP or BLUE

### Python examples

```python
import pandas as pd
from mrbigr.core import pheno

pheno_df = pd.read_csv("data/pheno.csv", index_col=0)
scaled = pheno.scale_wrapper(pheno_df, method="zscore")
blup_result = pheno.blup(pheno_df, method="matrix")
```

### Notes

- Use BLUP or BLUE when the phenotype design includes replicated environments.
- Keep sample IDs synchronized with the genotype dataset before launching GWAS.

## 4. GWAS Analysis

The project provides three primary GWAS paths:

- linear model (`gwas_lm`)
- linear mixed model (`gwas_lmm`)
- PLINK linear regression (`gwas_plink`)

### Recommended practice

- Prefer `gwas_lmm` when relatedness or structure may matter.
- Use `gwas_lm` mainly for simple baselines or method comparison.
- Use `gwas_plink` when PLINK-based validation is required.

### Python examples

```python
from mrbigr.core import gwas

gwas.gwas_lmm("data/pheno.csv", "data/chr_HAMP", num_threads=4, out_dir="output/gwas_lmm")
top_hits = gwas.get_top_snps("output/gwas_lmm/trait.assoc.txt", n=10)
lambda_gc = gwas.calculate_lambda("output/gwas_lmm/trait.assoc.txt")
```

### Interpretation notes

- Lambda values close to 1 usually indicate better calibration.
- Strong inflation in LM but not LMM often suggests population structure effects rather than an implementation bug.

## 5. Annotation and QTL Analysis

Annotation and QTL workflows are closely linked in practice.

### Annotation tasks

- parse a GTF file
- annotate SNPs
- query genes in a region
- annotate QTL intervals
- create a local annotation database
- predict variant effects

### QTL tasks

- detect QTL regions from GWAS output
- identify lead SNPs
- summarize QTLs
- map QTLs to genes
- export BED intervals

### Python examples

```python
from mrbigr.core import anno, qtl

gtf_df = anno.parse_gtf("data/genes.gtf.gz")
annotated = anno.annotate_snps_simple(["1:46746", "1:51053"], "data/genes.gtf.gz")
qtl_df = qtl.detect_qtl("output/trait.assoc.txt", p1=1e-7, p2=1e-5, p2n=5)
lead = qtl.get_lead_snp("output/trait.assoc.txt", "1", 100000, 300000)
```

### Notes

- Use real annotation files whenever possible; QTL interpretation quality depends heavily on annotation quality.
- Region overlap logic and point-QTL logic should both be validated during testing.

## 6. Peak and Haplotype Analysis

The peak module focuses on local post-QTL interpretation.

### Common tasks

- region-level boxplots
- haplotype tests
- gene haplotype tests
- SNP extraction for local windows

### Notes

- These analyses are most useful after lead SNPs or QTL intervals are already known.
- Keep phenotype and genotype sample ordering aligned before haplotype testing.

## 7. Mendelian Randomization

The MR module supports summary-statistic workflows.

### Common tasks

- run MR estimation
- calculate causal estimates
- test pleiotropy
- test heterogeneity
- format QTL tables for MR
- run QTL target analysis

### Notes

- Treat instrument validity checks as part of the main workflow, not as optional extras.
- Heterogeneity and pleiotropy outputs are part of the interpretation layer, not just diagnostics.

## 8. GO, KEGG, and GSEA

The GO module supports both local and online workflows.

### GO enrichment modes

- `local`: preferred when species-specific local resources are available
- `online`: uses remote Enrichr-style services
- `auto`: prefers local input when gene-set resources are provided

### Common tasks

- run GO enrichment
- run KEGG enrichment
- run GSEA
- draw GO barplots or dotplots
- draw GSEA curves, trace plots, or NES summary plots
- simplify GO terms
- export reports

### Notes

- For maize workflows, local mode is often the most reliable option.
- Local mode benefits from a consistent `gene -> GO` file and a matching `go-basic.obo`.
- Standard GSEA curves, trace plots, and NES barplots serve different interpretation needs.

## 9. Network Analysis

The network module supports edge weighting, module detection, and hub identification.

### Important distinction

- If Java is available, the bundled ClusterONE jar can be used.
- If Java is unavailable, the workflow falls back to a NetworkX-based implementation.

### Notes

- The fallback path is operationally useful, but should not be treated as numerically identical to ClusterONE.
- Hub results depend on the module definition, so module-path differences can propagate into biological interpretation.

## 10. Visualization

Visualization is implemented in `mrbigr.core.vis` and is used both directly and downstream of analysis modules.

### Available plot classes

- Manhattan plots
- QQ plots
- PCA plots
- LD heatmaps
- phenotype histograms
- phenotype boxplots
- phenotype correlation plots

### Notes

- For GWAS workflows, Manhattan and QQ plots should normally be treated as part of the default reporting bundle.
- For GO workflows, choose either split-ontology or combined plotting based on reporting needs.

## 11. Parallel Execution (Internal)

`parallel.py` (formerly `multi.py`) provides common infrastructure for parallel execution and thread selection. It is **internal infrastructure** used by other core modules (e.g. `gwas.py`) and is not exposed as MCP tools.

### Key functions

- `run_cmd` — execute a shell command
- `parallel_run` — run multiple commands in parallel
- `parallel_map` / `parallel_starmap` — parallel mapping
- `get_optimal_threads` — recommended thread count

## 12. Quick Start Checklist

Use the following checklist for a new environment:

1. Clone the repository from `https://github.com/CrazyHsu/MRBIGR2`.
2. Run `./install.sh`.
3. If network module parity with ClusterONE matters, rerun `./install.sh --with-java`.
4. Verify imports with a minimal smoke test.
5. Register MCP servers with `mrbigr install-all`, or start the aggregated server with `python src/server.py`.
6. Keep all derived outputs in a dedicated output directory and never overwrite original source data.

## Author

- **CrazyHsu**
- **Email:** `crazyhsu9527@gmail.com`
