---
name: QTL Mapping
description: Detect QTL regions from GWAS results, extract haplotypes, test allelic effects, and annotate candidate genes.
tags: [qtl, gwas, haplotype]
tools:
  - qtl_mcp
  - peak_mcp
  - anno_mcp
  - geno_mcp
---

# QTL Mapping

Post-GWAS workflow for QTL detection, haplotype analysis, and candidate gene identification.

## Inputs

- `GWAS_FILE`: GWAS result file (`.assoc.txt`)
- `GENO_PREFIX`: PLINK bed/bim/fam prefix
- `PHENO_FILE`: phenotype CSV
- `GTF_FILE`: GTF annotation file

## Steps

### 1. Detect QTL regions

```
detect_qtl_regions(gwas_file=GWAS_FILE, p1=1e-7, p2=1e-5, window=500000)
```

### 2. Identify peak SNPs

```
identify_peak_snps(gwas_dir="output/gwas/", p_threshold=1e-5, max_peaks=50)
```

### 3. Extract QTL genotypes

```
extract_qtl_genotypes(geno_prefix=GENO_PREFIX, qtl_regions=qtl_df, output_prefix="output/qtl_geno")
```

### 4. Haplotype analysis

```
calculate_qtl_haplo(geno_prefix=GENO_PREFIX, qtl_snp_list=peak_snps, output_prefix="output/haplo")
```

### 5. Allelic effect test

```
haplotype_test(pheno_df=pheno, geno_df=geno, snp_id=lead_snp, test_method="t-test")
plot_qtl_boxplot(pheno_file=PHENO_FILE, geno_prefix=GENO_PREFIX, qtl_df=qtl_regions)
```

### 6. Map QTL to genes

```
map_qtl_to_genes(qtl_df=qtl_regions, annotation_file=GTF_FILE)
```

### 7. Regional visualization

```
plot_qtl_region(gwas_file=GWAS_FILE, chr_val=chr, start=start, end=end, output_file="output/qtl_region.png")
```

### 8. Export results

```
export_qtl_bed(qtl_df=qtl_regions, output_file="output/qtl.bed")
qtl_summary(gwas_dir="output/gwas/", output_prefix="output/qtl_summary")
```

## Output

- QTL region table with lead SNPs and p-values
- Haplotype assignments per QTL
- Boxplots showing allelic effects
- BED file for downstream analyses
- Candidate gene list per QTL
