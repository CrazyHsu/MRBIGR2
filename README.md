# MRBIGR2

**MRBIGR2: An agentic toolbox for multi-omics association analysis and genetic regulation inference**

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview

MRBIGR2 is an agentic toolbox for multi-omics association analysis and genetic regulation inference. It extends the MRBIGR workflow into a pure Python environment and exposes the project through three complementary interfaces:

- **10 composable MCP servers** for agent-driven workflows (79 tools)
- A **meta-orchestrator CLI** (`mrbigr install/list/status`) for one-command MCP registration
- A **file-oriented CLI** (`mrbigr_cli.py`, 55 tools) for direct shell usage

Repository: `https://github.com/CrazyHsu/MRBIGR2`

Original MRBIGR: [Paper (Plant Communications 2025)](https://doi.org/10.1016/j.xplc.2024.101197) | [Repository](https://gitee.com/crazyhsu/MRBIGR) | [Tutorial](https://mrbigr.github.io/)

## Architecture

```
MRBIGR2/
├── src/mrbigr/                      # shared core + meta-orchestrator
│   ├── core/                        # 11 domain modules (geno, pheno, gwas, ...)
│   ├── mcp/                         # MCP dataclass + manager (reads mcps.yaml)
│   ├── skill/                       # skill deployer → ~/.claude/skills/
│   ├── mcp_cli.py                   # `mrbigr install|list|status|uninstall`
│   └── skill_cli.py                 # `mrbigr-skill install|list`
├── tool-mcps/                       # 10 per-domain MCP servers
│   ├── geno_mcp/   (12 tools)      # genotype QC, PCA, IBD, kinship, conversion
│   ├── pheno_mcp/  (8 tools)       # filtering, scaling, BLUP/BLUE, imputation
│   ├── gwas_mcp/   (12 tools)      # LM/LMM/PLINK GWAS, clumping, top SNPs
│   ├── vis_mcp/    (9 tools)       # Manhattan, QQ, PCA, LD, phenotype plots
│   ├── anno_mcp/   (6 tools)       # GTF parsing, SNP annotation, variant effect
│   ├── qtl_mcp/    (11 tools)      # QTL detection, haplotype, gene mapping
│   ├── mr_mcp/     (6 tools)       # IVW/MR-Egger, pleiotropy, heterogeneity
│   ├── go_mcp/     (9 tools)       # GO/KEGG/GSEA enrichment + plotting
│   ├── net_mcp/    (2 tools)       # module + hub identification
│   └── peak_mcp/   (4 tools)       # QTL boxplot, haplotype test
├── skills/                          # Claude Code workflow skills (.md)
├── notebooks/                       # Jupyter example workflows
├── mcps.yaml                        # authoritative MCP registry
├── src/server.py                    # compatibility aggregator (all 79 tools)
├── src/mrbigr_cli.py                # file-based CLI (55 tools)
└── utils/                           # bundled binaries (plink, gemma, ...)
```

## Installation

### Clone and install

```bash
git clone https://github.com/CrazyHsu/MRBIGR2.git
cd MRBIGR2
pip install -e .
```

Or use the full conda setup:

```bash
./install.sh              # creates mrbigr2 conda env
./install.sh --with-java  # + Java for ClusterONE
```

## Quick Start

### Mode 1: Per-domain MCPs (recommended)

Register individual MCPs with Claude Code:

```bash
mrbigr install geno_mcp      # register genotype tools
mrbigr install gwas_mcp      # register GWAS tools
mrbigr install-all            # register all 10 MCPs
mrbigr list                   # show status of all MCPs
```

Each MCP also runs standalone:

```bash
python tool-mcps/geno_mcp/src/server.py
```

### Mode 2: Single aggregated server (compatibility)

All 79 tools on one server — useful for Claude Desktop:

```bash
python src/server.py
```

See `claude_desktop_config.example.json` for Claude Desktop configuration.

### Mode 3: File-based CLI

```bash
python src/mrbigr_cli.py list
python src/mrbigr_cli.py gwas_lmm --phe data/pheno.csv --geno data/chr_HAMP
```

### Mode 4: Python API

```python
from mrbigr.core import geno, gwas, vis

pc_df, var_ratio = geno.calculate_pca("data/chr_HAMP", n_components=5)
vis.manhattan_plot("output/gwas_result.assoc.txt", output_file="output/manhattan.png")
```

## Skills

Install Claude Code workflow skills (deployed to `~/.claude/skills/`):

```bash
mrbigr-skill install gwas_pipeline           # full GWAS from raw genotype + phenotype
mrbigr-skill install qtl_to_target           # post-GWAS candidate-gene prioritization
mrbigr-skill install causal_network          # MR + module/hub network construction
mrbigr-skill install functional_enrichment   # GO / KEGG / GSEA + report
mrbigr-skill install reproduce_v1_case       # meta-skill chaining the four above
mrbigr-skill list
```

Each skill is a markdown file under `skills/` describing trigger phrases,
expected inputs, the fixed MCP tool sequence, interaction checkpoints, and
failure modes. Once installed, a Claude Code session in any directory will
match user intent against these skills automatically.

`qtl_mapping` and `mr_analysis` are kept as deprecated redirect stubs
pointing at `qtl_to_target` and `causal_network` respectively.

## Domain Modules

| Module | Tools | Purpose |
| --- | --- | --- |
| `geno_mcp` | 12 | Genotype QC, PCA, IBD, kinship, format conversion, pruning, clumping |
| `pheno_mcp` | 8 | Filtering, scaling, imputation, outlier removal, BLUP, BLUE |
| `gwas_mcp` | 12 | LM/LMM/PLINK GWAS, clumping, top SNPs, lambda, QQ/Manhattan data |
| `vis_mcp` | 9 | Manhattan, QQ, PCA, t-SNE, LD heatmap, phenotype plots |
| `anno_mcp` | 6 | GTF parsing, SNP annotation, variant effect prediction |
| `qtl_mcp` | 11 | QTL detection, peak SNPs, haplotype, gene mapping, BED export |
| `mr_mcp` | 6 | IVW/MR-Egger causal estimates, pleiotropy, heterogeneity |
| `go_mcp` | 9 | GO/KEGG/GSEA enrichment, plotting, simplification |
| `net_mcp` | 2 | Network module + hub identification |
| `peak_mcp` | 4 | QTL boxplot, grouped boxplot, haplotype test |

## External Tools

Bundled in `utils/`:

- **PLINK** for genotype processing
- **GEMMA** for LMM GWAS and kinship
- **FastTree** for tree construction
- **ClusterONE** for network module detection (requires Java)
- **ANNOVAR helpers** for annotation (Perl)

## Documentation

- [TUTORIAL.md](./TUTORIAL.md): module-by-module usage tutorial
- [MCP_GUIDE.md](./MCP_GUIDE.md): MCP usage guide and workflow patterns
- [CHANGELOG.md](./CHANGELOG.md): project change history
- `notebooks/`: Jupyter workflow examples (GWAS, MR)
- `skills/`: Claude Code skill definitions
- [Original MRBIGR Tutorial](https://mrbigr.github.io/)

## Testing

```bash
# Import check
python -c "from mrbigr.core import geno, pheno, gwas, vis, anno, qtl, mr, go, net, peak; print('OK')"

# CLI listing
python src/mrbigr_cli.py list

# MCP status
mrbigr list

# Run smoke tests
pytest tests/smoke_test.py
```

## Author

- **CrazyHsu**
- **Email:** `crazyhsu9527@gmail.com`

## License

MIT License

## Citation

```bibtex
@software{MRBIGR22026,
  title = {MRBIGR2: An agentic toolbox for multi-omics association analysis and genetic regulation inference},
  author = {CrazyHsu},
  year = {2026},
  url = {https://github.com/CrazyHsu/MRBIGR2}
}

@article{MRBIGR2025,
  title = {MRBIGR: A versatile toolbox for genetic regulation inference from population-scale multi-omics data},
  journal = {Plant Communications},
  year = {2025},
  doi = {10.1016/j.xplc.2024.101197}
}
```
