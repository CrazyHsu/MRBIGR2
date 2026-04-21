# MRBIGR2 MCP Guide

## Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Architecture](#architecture)
4. [Tool Organization](#tool-organization)
5. [Workflow Skills](#workflow-skills)
6. [Prompting Notes](#prompting-notes)
7. [Interface Boundaries](#interface-boundaries)

---

## Overview

MRBIGR2 exposes the MRBIGR analysis stack as **10 composable MCP servers** totalling **79 tools**, managed by a meta-orchestrator CLI. A single-server compatibility mode (`src/server.py`) is also available.

The MCP layer is intended for structured, agent-driven orchestration. The CLI layer in `src/mrbigr_cli.py` (55 tools) remains available for direct file-based shell usage.

## Quick Start

### Option A: Per-domain MCPs (recommended)

```bash
pip install -e .

# Register individual MCPs with a supported MCP client
mrbigr install geno_mcp
mrbigr install gwas_mcp --client codex
mrbigr install gwas_mcp --client gemini

# Or register all 10 at once
mrbigr install-all --client all

# Or export config for OpenCode / generic MCP clients
mrbigr export-config --format opencode
mrbigr export-config --format gemini
mrbigr export-config --format mcpservers

# Check status
mrbigr list --client all
```

### Option B: Single aggregated server

```bash
python src/server.py
```

### Option C: Client configuration

Use `mrbigr export-config` for clients that consume JSON configuration.
Claude Desktop can also use `claude_desktop_config.example.json` after you
update the paths.

### MCP configuration

The authoritative MCP registry is `mcps.yaml`. The meta-orchestrator reads it
to drive `mrbigr install/install-all/list/status/uninstall/export-config`.
For Claude Desktop, see `claude_desktop_config.example.json`.

## Architecture

```
mcps.yaml (registry)
    │
    ├── mrbigr install <name> --client claude|codex|gemini
    ├── mrbigr list --client all
    ├── mrbigr uninstall <name> --client claude|codex|gemini
    └── mrbigr export-config --format mcpservers  (or gemini/opencode)

tool-mcps/
    ├── geno_mcp/src/server.py   → standalone FastMCP (12 tools)
    ├── pheno_mcp/src/server.py  → standalone FastMCP (8 tools)
    ├── gwas_mcp/src/server.py   → standalone FastMCP (12 tools)
    └── ...                      → 7 more MCPs

src/server.py                    → aggregator: imports all register() functions
src/mrbigr/core/                 → shared domain modules
```

Each tool-mcp can run standalone (`python tool-mcps/<name>/src/server.py`) or be registered via the meta-orchestrator (`mrbigr install <name>`).

## Tool Organization

### geno_mcp (12 tools)

`run_snp_qc`, `subset_genotype`, `run_genotype_pca`, `run_calculate_ibd`, `convert_vcf`, `convert_hapmap`, `run_plink_to_vcf`, `calculate_kinship`, `impute_genotype`, `run_snp_pruning`, `run_snp_clumping`, `get_snp_statistics`

### pheno_mcp (8 tools)

`filter_abundance`, `filter_missing`, `scale_phenotype`, `impute_phenotype`, `remove_outliers`, `run_blup`, `run_blue`, `correct_trait`

### gwas_mcp (12 tools)

`run_gwas_lm`, `run_gwas_lmm`, `run_gwas_plink`, `add_rs_id_to_vcf`, `ensure_rs_id`, `run_simple_gwas`, `generate_clump`, `run_gwas_clump`, `get_top_snps`, `calculate_lambda`, `qq_plot_data`, `manhattan_plot_data`

### vis_mcp (9 tools)

`plot_manhattan`, `plot_qq`, `plot_pca`, `plot_tsne`, `plot_ld_heatmap`, `plot_phenotype_hist`, `plot_phenotype_boxplot`, `plot_phenotype_correlation`, `gwas_summary`

### anno_mcp (6 tools)

`parse_gtf`, `annotate_snps`, `qtl_annotation`, `create_annotation_db`, `predict_variant_effect`, `get_genes_in_region`

### qtl_mcp (11 tools)

`detect_qtl_regions`, `get_lead_snp`, `identify_peak_snps`, `extract_qtl_genotypes`, `calculate_qtl_haplo`, `plot_qtl_region`, `qtl_summary`, `map_qtl_to_genes`, `qtl_enrichment_test`, `get_qtl_overlap`, `export_qtl_bed`

### mr_mcp (6 tools)

`run_mr_analysis`, `calculate_mr_causal_estimate`, `test_mr_pleiotropy`, `test_mr_heterogeneity`, `run_qtl_target_analysis`, `format_qtl_for_mr`

### go_mcp (9 tools)

`run_go_enrichment`, `run_gsea_analysis`, `plot_go_enrichment`, `plot_gsea_results`, `run_kegg_enrichment`, `get_enrichr_libraries`, `extract_go_from_gtf`, `simplify_go_results`, `export_go_report`

### net_mcp (2 tools)

`module_identify`, `hub_identify`

### peak_mcp (4 tools)

`plot_qtl_boxplot`, `plot_grouped_boxplot`, `haplotype_test`, `multi_trait_qtl_plot`

## Workflow Skills

Five dialog-driven workflow skills are available under `skills/`. They are
the primary interface for agent-driven analysis: each skill defines a fixed
MCP tool sequence, the parameters the agent must ask the user about
(interaction checkpoints), and the failure conditions that should pause the
workflow rather than be auto-recovered.

`./install.sh` installs the active skills globally into all built-in
Agent Skills-compatible targets (`claude`, `codex`, `gemini`, `opencode`, and `agents`).
You can also manage them manually:

```bash
mrbigr-skill install-all --target all
mrbigr-skill install-all --target claude
mrbigr-skill install-all --target codex
mrbigr-skill install-all --target gemini
mrbigr-skill install-all --target opencode
mrbigr-skill install-all --target-dir ~/.someagent/skills
mrbigr-skill install gwas_pipeline
mrbigr-skill install qtl_to_target
mrbigr-skill install causal_network
mrbigr-skill install functional_enrichment
mrbigr-skill install reproduce_v1_case
mrbigr-skill list
```

### gwas-pipeline

Genotype QC → kinship → PCA → phenotype prep (filter / outlier / scale) →
GWAS (LMM by default) → λ + Manhattan + QQ → QTL region detection → lead
SNP → gene mapping. Six checkpoints (QC thresholds, normalization, PCA k,
model, significance threshold, flanking window).

### qtl-to-target

Post-GWAS candidate-gene prioritization for one QTL region: gene listing →
genotype extraction → SNP / variant-effect annotation → per-gene
expression → trait MR → haplotype effect test → composite ranking. Three
hard interaction gates (region confirmation, candidate confirmation,
deliverable scope).

### causal-network

MR + network construction: instrument clumping → IVW + Egger + weighted
median across exposure–outcome pairs → pleiotropy + heterogeneity filters
→ ClusterONE / NetworkX module identification → hub identification →
forest + network plots. Five interaction checkpoints.

### functional-enrichment

GO / KEGG / GSEA on a gene list or a ranked list (auto-detected by the
agent), with redundancy simplification and an HTML / markdown report.
Single interaction point: gene-list source confirmation.

### reproduce-v1-case

Meta-skill that chains `gwas-pipeline → qtl-to-target → causal-network →
functional-enrichment` end-to-end with three hard user-confirmation gates
between stages. Designed for case-study reproduction (e.g. the V1
P1–flavonoid analysis); emits `chat_log.md` + `tool_trace.jsonl` as the
primary deliverables.

The legacy skill files `qtl_mapping.md` and `mr_analysis.md` are retained
as redirect stubs and should not be used directly.

## Prompting Notes

Recommended prompting style:

- Specify the exact input files
- Specify the desired output directory
- Specify the method when multiple methods are available
- State whether intermediate plots are expected
- State whether online services are allowed

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

- The MCP layer is the primary structured interface for interactive workflows (79 tools across 10 MCPs).
- The CLI layer is a direct file-based interface (55 tools).
- Both layers share the same `mrbigr.core` modules.

### Java-dependent network behavior

- If Java is available, `module_identify` can use the bundled ClusterONE jar.
- If Java is not available, the same entry point falls back to a NetworkX-based implementation.

### GO enrichment modes

- `online` uses remote Enrichr-style services.
- `local` uses local gene-set resources such as `maize.genes2go.txt` and `go-basic.obo`.
- `auto` selects the local path when a local gene-set source is provided.

## Author

- **CrazyHsu**
- **Email:** `crazyhsu9527@gmail.com`
