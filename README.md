# MRBIGR2

**MRBIGR2: An agentic toolbox for multi-omics association analysis and genetic regulation inference**

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview

MRBIGR2 is an agentic toolbox for multi-omics association analysis and genetic regulation inference. It extends the MRBIGR workflow into a pure Python environment, removes the original R dependency chain from the active pipeline, and exposes the project through both:

- an MCP server for agent-driven workflows
- a file-oriented CLI for direct shell usage

The project covers genotype processing, phenotype preprocessing, GWAS, QTL analysis, Mendelian randomization, annotation, GO/KEGG analysis, network analysis, visualization, and utility helpers.

Repository: `https://github.com/CrazyHsu/MRBIGR2`

Original MRBIGR repository: `https://gitee.com/crazyhsu/MRBIGR`

Original MRBIGR tutorial: `https://mrbigr.github.io/`

## Main Modules

| Module | Purpose | Typical Functions |
| --- | --- | --- |
| `geno` | Genotype preprocessing and conversion | QC, subset, PCA, kinship, IBD, pruning, clumping, imputation |
| `pheno` | Phenotype preprocessing and mixed-model summaries | scaling, imputation, outlier handling, BLUP, BLUE |
| `gwas` | Association analysis | LM, LMM, PLINK-based GWAS, top SNP extraction, lambda, clumping prep |
| `vis` | Visualization | Manhattan, QQ, PCA, LD heatmap, phenotype distribution plots |
| `anno` | Annotation | GTF parsing, SNP annotation, QTL annotation, variant effect prediction |
| `qtl` | QTL analysis | region detection, lead SNPs, genotype extraction, gene mapping, BED export |
| `peak` | Peak and haplotype analysis | haplotype tests, region plots, gene haplotype workflows |
| `mr` | Mendelian randomization | causal effect estimation, heterogeneity, pleiotropy, QTL target analysis |
| `go` | GO/KEGG/GSEA analysis | local and online enrichment, plotting, simplification, report export |
| `net` | Network analysis | edge weighting, module detection, hub identification, network statistics |
| `multi` | Parallel and utility helpers | shell execution, parallel map, safe map, thread recommendation |

## Installation

### Clone the repository

```bash
git clone https://github.com/CrazyHsu/MRBIGR2.git
cd MRBIGR2
```

### Recommended installation

```bash
./install.sh
```

This creates or updates a conda environment named `mrbigr2` by default.

### Optional Java support for ClusterONE

Java is **not** installed by default.

Without Java, `net.module_identify` still works, but it uses the NetworkX fallback implementation instead of the bundled `cluster_one-1.0.jar`.

If you want the stricter ClusterONE execution path, install Java into the conda environment explicitly:

```bash
./install.sh --with-java
```

### Editable installation

```bash
pip install -e .
```

## Quick Start

### Example 1: Genotype PCA

```python
import sys
sys.path.insert(0, "src")

from geno import calculate_pca

pc_df, variance_ratio = calculate_pca("data/geno_output_qc", n_components=10)
print(pc_df.shape)
print(variance_ratio[:5])
```

### Example 2: GWAS result visualization

```python
import sys
sys.path.insert(0, "src")

from vis import manhattan_plot, qq_plot

manhattan_plot("output/example.assoc.txt", output_file="output/example_manhattan.png")
qq_plot("output/example.assoc.txt", output_file="output/example_qq.png")
```

### Example 3: Phenotype preprocessing

```python
import sys
sys.path.insert(0, "src")

import pandas as pd
from pheno import zscore_scale

pheno_df = pd.read_csv("data/phenotype.csv", index_col=0)
scaled = zscore_scale(pheno_df)
print(scaled.head())
```

### Example 4: Start the MCP server

```bash
python src/server.py
```

## Interfaces

### MCP server

`src/server.py` currently exposes **81 MCP tools**.

These tools cover genotype, phenotype, GWAS, visualization, annotation, QTL, MR, GO/KEGG, network, and utility workflows.

### CLI

`src/mrbigr_cli.py` currently exposes **55 CLI tools** for file-based shell workflows.

The CLI is useful for direct scripting, while the MCP layer is better suited for interactive or agentic orchestration.

## Project Layout

```text
MRBIGR2/
├── src/
│   ├── anno.py
│   ├── geno.py
│   ├── go.py
│   ├── gwas.py
│   ├── mr.py
│   ├── multi.py
│   ├── mrbigr_cli.py
│   ├── net.py
│   ├── peak.py
│   ├── pheno.py
│   ├── qtl.py
│   ├── server.py
│   └── vis.py
├── utils/
│   ├── FastTree
│   ├── cluster_one-1.0.jar
│   ├── gemma.linux
│   ├── plink
│   └── *.pl
├── README.md
├── MCP_GUIDE.md
├── TUTORIAL.md
├── CHANGELOG.md
├── mcp_config.yaml
├── pyproject.toml
├── install.sh
└── quick_setup.sh
```

## External Tools

The repository bundles or expects the following external tools:

- **PLINK** for genotype processing
- **GEMMA** for LMM GWAS and kinship-related workflows
- **FastTree** for tree construction
- **ClusterONE** for network module detection when Java is available
- **ANNOVAR helper scripts** for annotation-related workflows that depend on Perl

## Documentation

- [TUTORIAL.md](./TUTORIAL.md): module-by-module usage tutorial
- [MCP_GUIDE.md](./MCP_GUIDE.md): MCP usage guide and workflow patterns
- [CHANGELOG.md](./CHANGELOG.md): project change history
- [test_prompt.all_55_tools.md](./test_prompt.all_55_tools.md): consolidated testing prompt catalog
- [Original MRBIGR Tutorial](https://mrbigr.github.io/): upstream workflow documentation
- [Original MRBIGR Repository](https://gitee.com/crazyhsu/MRBIGR): upstream source repository

## Testing

### Minimal import check

```bash
python - <<'PY'
import sys
sys.path.insert(0, 'src')

import anno, geno, go, gwas, mr, multi, net, peak, pheno, qtl, vis
print("All core modules imported successfully.")
PY
```

### CLI listing

```bash
python src/mrbigr_cli.py list
```

### MCP server startup

```bash
python src/server.py
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
