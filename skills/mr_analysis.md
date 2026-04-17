---
name: MR Analysis
description: Mendelian Randomization workflow — format instruments, estimate causal effects, test pleiotropy and heterogeneity.
tags: [mr, causal, mendelian-randomization]
tools:
  - mr_mcp
  - qtl_mcp
  - go_mcp
---

# Mendelian Randomization Analysis

Two-sample MR workflow to estimate causal effects between an exposure and an outcome trait.

## Inputs

- `EXPOSURE_GWAS`: GWAS summary statistics for the exposure trait
- `OUTCOME_GWAS`: GWAS summary statistics for the outcome trait
- `QTL_FILE` (optional): QTL file to format as MR instruments

## Steps

### 1. Format instruments (if starting from QTL)

```
format_qtl_for_mr(qtl_file=QTL_FILE, output_file="output/mr_instruments.csv")
```

### 2. Run MR analysis (IVW)

```
run_mr_analysis(exposure_gwas=EXPOSURE_GWAS, outcome_gwas=OUTCOME_GWAS, method="ivw")
```

### 3. Causal estimate

```
calculate_mr_causal_estimate(
    exposure_gwas=EXPOSURE_GWAS, outcome_gwas=OUTCOME_GWAS,
    snp_col="rs", beta_col="beta", se_col="se", method="ivw"
)
```

### 4. Sensitivity analyses

#### Pleiotropy test (MR-Egger intercept)

```
test_mr_pleiotropy(exposure_gwas=EXPOSURE_GWAS, outcome_gwas=OUTCOME_GWAS)
```

#### Heterogeneity test (Cochran's Q)

```
test_mr_heterogeneity(exposure_gwas=EXPOSURE_GWAS, outcome_gwas=OUTCOME_GWAS)
```

### 5. QTL-target analysis (optional)

```
run_qtl_target_analysis(qtl_df=qtl_data, tf_genes=tf_list, target_genes=target_list, window=500000)
```

### 6. Functional enrichment of causal genes (optional)

```
run_go_enrichment(gene_list=causal_genes, organism="Mouse", mode="local")
run_kegg_enrichment(gene_list=causal_genes, organism="mmu")
```

## Interpretation

| Result | Meaning |
|--------|---------|
| IVW p < 0.05 | Evidence for causal effect |
| MR-Egger intercept p > 0.05 | No directional pleiotropy |
| Cochran's Q p > 0.05 | No significant heterogeneity |

## Output

- Causal effect estimates (beta, SE, p-value)
- Pleiotropy and heterogeneity test results
- Enrichment results for causal candidate genes
