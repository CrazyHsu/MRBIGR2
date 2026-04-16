# geno_mcp

Genotype processing MCP — 12 tools built on top of `mrbigr.core.geno` and the
bundled PLINK / GEMMA / FastTree binaries under `utils/`.

## Tools

| Tool | Purpose |
| --- | --- |
| `run_snp_qc` | PLINK-based SNP QC (MAF, missingness, mind) |
| `subset_genotype` | subset PLINK data by chromosome or random proportion |
| `run_genotype_pca` | PCA on the genotype matrix |
| `run_calculate_ibd` | Identity-by-descent matrix |
| `convert_vcf` | VCF → PLINK |
| `convert_hapmap` | HapMap → PLINK |
| `run_plink_to_vcf` | PLINK → VCF |
| `calculate_kinship` | GEMMA kinship matrix |
| `impute_genotype` | fill missing genotype values |
| `run_snp_pruning` | LD-based PLINK pruning |
| `run_snp_clumping` | greedy LD clumping |
| `get_snp_statistics` | summary stats |

## Standalone run

```bash
python tool-mcps/geno_mcp/src/server.py
```

## Register with Claude Code

```bash
mrbigr install geno_mcp
```
