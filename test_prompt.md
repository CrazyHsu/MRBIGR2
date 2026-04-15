# MRBIGR2 Test Prompt Catalog

This document consolidates the testing prompts for the current MRBIGR MCP and CLI workflows. It is organized by module and written as natural-language requests that can be used during interactive testing.

---

## Genotype - 13 tools

```text
1. snp_qc: "Run SNP quality control on my genotype data with MAF > 0.05 and missing rate < 0.2."
2. subset_genotype: "Extract a genotype subset from chromosome 1, or keep a random 10% SNP subset with seed 42."
3. genotype_pca: "Calculate PCA for my genotype data and keep 10 principal components."
4. calculate_kinship: "Use GEMMA to calculate the kinship matrix from my genotype data."
5. calculate_ibd: "Calculate the identity-by-descent matrix."
6. snp_pruning: "Run LD-based SNP pruning with an r2 threshold of 0.5."
7. snp_clumping: "Run SNP clumping with r2 > 0.5."
8. snp_impute: "Impute missing genotype values using the mean-based method."
9. snp_stats: "Return basic SNP statistics for my genotype dataset."
10. convert_vcf: "Convert a VCF file to PLINK format."
11. convert_hapmap: "Convert a HapMap file to PLINK format."
12. plink_to_vcf: "Convert a PLINK dataset to VCF format."
13. build_tree: "Build a phylogenetic tree from my genotype dataset."
```

## Phenotype - 8 tools

```text
1. filter_abundance: "Filter phenotype features with abundance below 1."
2. filter_missing: "Filter samples with more than 20% missing values."
3. scale_phenotype: "Apply Z-score normalization to my phenotype matrix."
4. impute_phenotype: "Impute phenotype missing values using the mean."
5. remove_outliers: "Detect and mark phenotype outliers."
6. blup: "Run BLUP on a multi-environment phenotype matrix, or on long-format data using genotype, environment, and response columns."
7. blue: "Run BLUE on a multi-environment phenotype matrix, or on long-format data using genotype, environment, and response columns."
8. trait_correct: "Correct phenotype values using the genotype PCA results."
```

## GWAS - 9 tools

```text
1. gwas_lm: "Run GWAS using a linear model."
2. gwas_lmm: "Run GWAS using a linear mixed model."
3. gwas_plink: "Run GWAS using PLINK linear regression."
4. get_top_snps: "Extract the top 10 SNPs from my GWAS result file."
5. calculate_lambda: "Calculate the genomic inflation factor from my GWAS results."
6. ensure_rs_id: "Ensure that the SNP IDs in my PLINK BIM file are valid rs-style IDs."
7. add_rs_id_to_vcf: "Add rs-style IDs to a VCF file."
8. generate_clump: "Generate PLINK clumping input files from GWAS results."
9. plink_clump: "Run PLINK clumping on GWAS results."
```

## Visualization - 8 tools

```text
1. plot_manhattan: "Draw a Manhattan plot from my GWAS result file."
2. plot_qq: "Draw a QQ plot to assess GWAS inflation."
3. plot_pca: "Draw a PCA scatter plot."
4. plot_ld_heatmap: "Draw an LD heatmap for SNPs in a target region."
5. plot_pheno_hist: "Draw a phenotype histogram."
6. plot_pheno_boxplot: "Draw a phenotype boxplot."
7. plot_pheno_correlation: "Calculate and visualize phenotype correlations."
8. gwas_summary: "Generate both Manhattan and QQ plots for the same GWAS result."
```

## Annotation - 6 tools

```text
1. parse_gtf: "Parse this GTF annotation file and return gene_id, gene_name, transcript_id, and gene_type fields."
2. annotate_snps: "Annotate these SNPs functionally. The input may be a chr:pos list, a SNP table, or a structured SNP list."
3. get_genes_in_region: "Return genes that overlap or are fully contained within a target genomic interval. Test both overlap and contain modes."
4. qtl_annotation: "Annotate QTL points or QTL intervals with overlapping genes."
5. create_annotation_db: "Create a local annotation database from this GTF file using the given output prefix."
6. predict_variant_effect: "Predict variant effects in either simple mode or full GTF-based mode."
```

## QTL - 11 tools

```text
1. detect_qtl_regions: "Detect QTL regions from a GWAS result file."
2. get_lead_snp: "Return the lead SNP for each QTL region."
3. identify_peak_snps: "Identify peak SNPs from a GWAS result file."
4. extract_qtl_genotypes: "Extract genotype data for the detected QTL regions."
5. calculate_qtl_haplo: "Calculate haplotype summaries from SNPs inside QTL regions."
6. plot_qtl_region: "Visualize GWAS results for a QTL region."
7. qtl_summary: "Generate a QTL analysis summary report."
8. map_qtl_to_genes: "Map QTL intervals to candidate genes."
9. qtl_enrichment_test: "Test whether QTL candidate genes are enriched for a given gene set."
10. get_qtl_overlap: "Compare two QTL result sets and find overlapping intervals."
11. export_qtl_bed: "Export QTL regions in BED format."
```

## MR - 6 tools

```text
1. run_mr_analysis: "Run Mendelian randomization using exposure and outcome summary statistics. Test both IVW and MR-Egger."
2. calculate_mr_causal_estimate: "Calculate a causal estimate from exposure and outcome summary statistics."
3. test_mr_pleiotropy: "Test horizontal pleiotropy in the MR analysis."
4. test_mr_heterogeneity: "Test heterogeneity across MR instruments."
5. run_qtl_target_analysis: "Analyze candidate regulatory links using QTL results, TF genes, and target genes."
6. format_qtl_for_mr: "Format a QTL table for downstream MR analysis."
```

## GO/KEGG - 9 tools

```text
1. run_go_enrichment: "Run GO enrichment on a candidate gene list. Test online, local, and auto modes."
2. run_gsea_analysis: "Run GSEA using either an online library or a local gene-set file."
3. plot_go_enrichment: "Draw GO enrichment barplots or dotplots. Test split-ontology mode and combined mode."
4. plot_gsea_results: "Draw GSEA results as a standard enrichment curve, a trace plot, or an NES summary barplot. The default should be the standard enrichment curve."
5. run_kegg_enrichment: "Run KEGG enrichment on a candidate gene list."
6. get_enrichr_libraries: "List the currently available Enrichr gene-set libraries."
7. extract_go_from_gtf: "Extract gene-to-GO mappings from a GTF annotation file."
8. simplify_go_results: "Simplify a GO enrichment result by removing redundant terms."
9. export_go_report: "Export GO enrichment results as csv, tsv, or Excel."
```

## Network - 2 tools

```text
1. module_identify: "Identify network modules from an edge-weight table using ClusterONE when Java is available."
2. hub_identify: "Identify hub genes within the detected network modules."
```

## Utility - 2 tools

```text
1. parallel_run_commands: "Run multiple shell commands in parallel and return their execution status."
2. get_optimal_threads: "Return the recommended number of worker threads for the current environment."
```

## Example Requests

```text
1. "Run SNP quality control on my PLINK dataset, then calculate genotype PCA with 10 components."
2. "Use the LMM method to run GWAS on my phenotype matrix, then generate Manhattan and QQ plots."
3. "Detect QTL regions from this GWAS result and export the intervals as a BED file."
4. "Annotate these significant SNPs and return overlapping genes from the GTF file."
5. "Run GO enrichment on this maize gene list using local mode and the local go-basic.obo file."
6. "Run MR analysis using these exposure and outcome summary files, then test pleiotropy and heterogeneity."
7. "Identify network modules from this edge-weight file and then report the hub genes for each module."
```
