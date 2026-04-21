# vis_mcp

Visualization MCP — 9 tools backed by `mrbigr.core.vis`.

## Tools

`plot_manhattan`, `plot_qq`, `plot_pca`, `plot_tsne`, `plot_ld_heatmap`,
`plot_phenotype_hist`, `plot_phenotype_boxplot`,
`plot_phenotype_correlation`, `gwas_summary`.

Standalone: `python tool-mcps/vis_mcp/src/server.py`
Register:   `mrbigr install vis_mcp` (Claude default; add `--client codex` or `--client gemini`)
Export:     `mrbigr export-config vis_mcp --format mcpservers` (or `gemini` / `opencode`)
