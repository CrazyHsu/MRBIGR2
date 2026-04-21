# qtl_mcp

QTL analysis MCP — 11 tools backed by `mrbigr.core.qtl`.

## Tools

`detect_qtl_regions`, `get_lead_snp`, `identify_peak_snps`,
`extract_qtl_genotypes`, `calculate_qtl_haplo`, `plot_qtl_region`,
`qtl_summary`, `map_qtl_to_genes`, `qtl_enrichment_test`,
`get_qtl_overlap`, `export_qtl_bed`.

Standalone: `python tool-mcps/qtl_mcp/src/server.py`
Register:   `mrbigr install qtl_mcp` (Claude default; add `--client codex` or `--client gemini`)
Export:     `mrbigr export-config qtl_mcp --format mcpservers` (or `gemini` / `opencode`)
