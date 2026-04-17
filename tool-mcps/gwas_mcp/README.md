# gwas_mcp

Association analysis MCP — 12 tools. Wraps `mrbigr.core.gwas` and
`mrbigr.core.vis` (for `auto_plot=True` Manhattan + QQ output).

## Tools

GEMMA-based: `run_gwas_lm`, `run_gwas_lmm`
PLINK-based: `run_gwas_plink`, `run_gwas_clump`
Helpers: `add_rs_id_to_vcf`, `ensure_rs_id`, `run_simple_gwas`,
`generate_clump`, `get_top_snps`, `calculate_lambda`,
`qq_plot_data`, `manhattan_plot_data`.

Standalone: `python tool-mcps/gwas_mcp/src/server.py`
Register:   `mrbigr install gwas_mcp`
