---
name: GWAS Pipeline
description: End-to-end GWAS workflow — genotype QC, phenotype prep, association analysis, visualization, and QTL detection.
tags: [gwas, genomics, pipeline]
tools:
  - geno_mcp
  - pheno_mcp
  - gwas_mcp
  - vis_mcp
  - qtl_mcp
  - anno_mcp
---

# GWAS Pipeline

Run a complete Genome-Wide Association Study from raw genotype + phenotype files to annotated QTL regions.

## Inputs

- `GENO_PREFIX`: PLINK bed/bim/fam prefix (e.g., `data/chr_HAMP`)
- `PHENO_FILE`: phenotype CSV with sample IDs in the first column
- `GTF_FILE` (optional): GTF annotation for gene mapping

## Steps

### 1. Genotype QC

```
run_snp_qc(input_prefix=GENO_PREFIX, output_prefix="output/qc", maf=0.05, missing_rate=0.2)
```

### 2. PCA for population structure

```
run_genotype_pca(input_prefix="output/qc", n_components=5)
```

### 3. Phenotype preparation

```
filter_missing(d=pheno_data, missing_ratio=0.3)
remove_outliers(d=filtered, method="zscore")
scale_phenotype(d=cleaned, method="zscore")
```

### 4. Kinship matrix

```
calculate_kinship(input_prefix="output/qc", output_prefix="output/kinship")
```

### 5. GWAS (Linear Mixed Model)

```
run_gwas_lmm(phe=scaled_pheno, geno_prefix="output/qc", output_dir="output/gwas")
```

### 6. Visualization

```
plot_manhattan(gwas_file="output/gwas/result.assoc.txt", output_file="output/manhattan.png")
plot_qq(gwas_file="output/gwas/result.assoc.txt", output_file="output/qq.png")
```

### 7. QTL detection

```
detect_qtl_regions(gwas_file="output/gwas/result.assoc.txt", p1=1e-7, p2=1e-5)
```

### 8. Gene annotation (if GTF provided)

```
map_qtl_to_genes(qtl_df=qtl_result, annotation_file=GTF_FILE)
```

## Output

- QC'd genotype files in `output/`
- GWAS results (`.assoc.txt`)
- Manhattan and QQ plots (`.png`)
- QTL region table with annotated genes
