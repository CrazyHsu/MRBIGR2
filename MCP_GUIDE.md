# MRBIGR2 MCP Guide

## Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Tool Organization](#tool-organization)
4. [Workflow Patterns](#workflow-patterns)
5. [Prompting Notes](#prompting-notes)
6. [Interface Boundaries](#interface-boundaries)

---

## Overview

MRBIGR2 exposes the MRBIGR analysis stack through an MCP server implemented in `src/server.py`.

The current server exports **81 MCP tools** covering:

- genotype processing
- phenotype preprocessing
- GWAS
- visualization
- annotation
- QTL
- peak and haplotype analysis
- Mendelian randomization
- GO/KEGG/GSEA
- network analysis
- utility helpers

The MCP layer is intended for structured, agent-driven orchestration. The CLI layer in `src/mrbigr_cli.py` remains available for direct file-based shell usage.

Repository: `https://github.com/CrazyHsu/MRBIGR2`

Original MRBIGR repository: `https://gitee.com/crazyhsu/MRBIGR`

Original MRBIGR tutorial: `https://mrbigr.github.io/`

## Quick Start

### 1. Install dependencies

```bash
./install.sh
```

By default this creates or updates a conda environment named `mrbigr2`.

If you want strict ClusterONE execution for network module detection, install Java as well:

```bash
./install.sh --with-java
```

### 2. Start the MCP server

```bash
python src/server.py
```

### 3. Example MCP configuration

```yaml
mcps:
  mrbigr2:
    runtime: python
    path: .
    server_command: python
    server_args: [src/server.py]
    env_vars:
      MRBIGR_ROOT: .
```

## Tool Organization

### Genotype

Representative tools:

- `run_snp_qc`
- `run_subset_genotype`
- `run_genotype_pca`
- `calculate_kinship`
- `calculate_ibd`
- `run_snp_pruning`
- `run_snp_clumping`
- `run_snp_imputation`
- `get_snp_stats`
- `convert_vcf_to_plink`
- `convert_hapmap_to_plink`
- `convert_plink_to_vcf`
- `build_phylogenetic_tree`

### Phenotype

Representative tools:

- `filter_abundance`
- `filter_missing`
- `scale_phenotype`
- `impute_phenotype`
- `remove_outliers`
- `run_blup`
- `run_blue`
- `trait_correction`

### GWAS

Representative tools:

- `run_gwas_lm`
- `run_gwas_lmm`
- `run_gwas_plink`
- `get_top_snps`
- `calculate_lambda`
- `ensure_rs_id`
- `add_rs_id_to_vcf`
- `generate_clump_input`
- `run_plink_clump`

### Visualization

Representative tools:

- `plot_manhattan`
- `plot_qq`
- `plot_pca`
- `plot_ld_heatmap`
- `plot_pheno_hist`
- `plot_pheno_boxplot`
- `plot_pheno_correlation`
- `gwas_summary`

### Annotation

Representative tools:

- `parse_gtf`
- `annotate_snps`
- `qtl_annotation`
- `create_annotation_db`
- `predict_variant_effect`
- `get_genes_in_region`

### QTL

Representative tools:

- `detect_qtl_regions`
- `get_lead_snp`
- `identify_peak_snps`
- `extract_qtl_genotypes`
- `calculate_qtl_haplo`
- `plot_qtl_region`
- `qtl_summary`
- `map_qtl_to_genes`
- `qtl_enrichment_test`
- `get_qtl_overlap`
- `export_qtl_bed`

### Peak

Representative tools:

- `plot_qtl_boxplot`
- `haplotype_test`
- `gene_haplotype_test`
- `extract_haplotype_snps`

### Mendelian Randomization

Representative tools:

- `run_mr_analysis`
- `calculate_mr_causal_estimate`
- `test_mr_pleiotropy`
- `test_mr_heterogeneity`
- `run_qtl_target_analysis`
- `format_qtl_for_mr`

### GO / KEGG / GSEA

Representative tools:

- `run_go_enrichment`
- `run_gsea_analysis`
- `plot_go_enrichment`
- `plot_gsea_results`
- `run_kegg_enrichment`
- `get_enrichr_libraries`
- `extract_go_from_gtf`
- `simplify_go_results`
- `export_go_report`

### Network

Representative tools:

- `module_identify`
- `hub_identify`

### Utility

Representative tools:

- `parallel_run_commands`
- `get_optimal_threads`

## Workflow Patterns

The following workflow patterns are intentionally documented as placeholders for future skill-style packaging. They are currently documentation patterns, not installed Codex skills.

### Pattern 1: GWAS pipeline

1. Prepare or filter the phenotype matrix.
2. Run `run_gwas_lmm` by default unless a different method is explicitly requested.
3. Generate Manhattan and QQ plots.
4. Summarize top SNPs and lambda.
5. Optionally generate clumping inputs and run PLINK clumping.

### Pattern 2: QTL pipeline

1. Start from a GWAS result file.
2. Detect QTL regions.
3. Identify lead SNPs.
4. Plot target regions.
5. Map QTL intervals to candidate genes.
6. Export BED or summary tables as needed.

### Pattern 3: MR pipeline

1. Prepare exposure and outcome summary statistics.
2. Run MR estimation.
3. Test pleiotropy.
4. Test heterogeneity.
5. Interpret candidate regulatory relationships.

### Pattern 4: GO pipeline

1. Build or load a candidate gene list.
2. Run GO enrichment in local mode by default unless online mode is explicitly requested.
3. Plot enrichment results.
4. Simplify redundant terms if necessary.
5. Export a structured report.

### Pattern 5: Network pipeline

1. Build or load an edge-weight table.
2. Run module detection.
3. Identify hubs per module.
4. Review whether Java-backed ClusterONE or the NetworkX fallback was used.

## Prompting Notes

Recommended prompting style:

- specify the exact input files
- specify the desired output directory
- specify the method when multiple methods are available
- state whether intermediate plots are expected
- state whether online services are allowed

Examples:

```text
Run LMM-based GWAS on this phenotype file and write all outputs to output/gwas_run_01.
```

```text
Annotate these SNPs with this GTF file and return overlapping genes for the interval chr1:100000-200000.
```

```text
Run local GO enrichment on this maize gene list using the local go-basic.obo file and save both csv output and plots.
```

## Interface Boundaries

### MCP vs CLI

- The MCP layer is the primary structured interface for interactive workflows.
- The CLI layer is a direct file-based interface and currently exposes **55 CLI tools**.
- The MCP and CLI layers overlap substantially, but they are not identical in all modules.

### Java-dependent network behavior

- If Java is available, `module_identify` can use the bundled ClusterONE jar.
- If Java is not available, the same entry point falls back to a NetworkX-based implementation.
- The fallback path is useful for testing and basic workflows, but it should not be assumed to be numerically identical to ClusterONE.

### GO enrichment modes

- `online` uses remote Enrichr-style services.
- `local` uses local gene-set resources such as `maize.genes2go.txt` and `go-basic.obo`.
- `auto` selects the local path when a local gene-set source is provided.

## Author

- **CrazyHsu**
- **Email:** `crazyhsu9527@gmail.com`
