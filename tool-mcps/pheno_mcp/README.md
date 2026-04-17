# pheno_mcp

Phenotype preprocessing MCP — 8 tools built on `mrbigr.core.pheno`.

## Tools

| Tool | Purpose |
| --- | --- |
| `filter_abundance` | drop features below a minimum abundance |
| `filter_missing` | drop features above a missing-ratio threshold |
| `scale_phenotype` | log2 / log10 / zscore / minmax / robust / boxcox |
| `impute_phenotype` | fill missing values (mean/median/etc.) |
| `remove_outliers` | flag outliers by zscore / IQR |
| `run_blup` | BLUP via REML mixed model |
| `run_blue` | BLUE via OLS fixed-effects model |
| `correct_trait` | regress out principal components |

Standalone: `python tool-mcps/pheno_mcp/src/server.py`
Register:   `mrbigr install pheno_mcp`
