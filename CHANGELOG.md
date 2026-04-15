# Changelog

## Version 1.0.0

**Repository:** `https://github.com/CrazyHsu/MRBIGR2`  
**Author:** `CrazyHsu <crazyhsu9527@gmail.com>`

---

## Summary

MRBIGR2 is an agentic, pure Python implementation of the MRBIGR workflow with MCP support, CLI exposure, and expanded analysis coverage beyond the original R-based layout.

## Major Structural Changes

### Python-first implementation

- removed the original R dependency chain from the active workflow
- consolidated analysis logic in Python modules
- exposed the project through both MCP and CLI interfaces

### Expanded module layout

The current repository includes:

- `geno.py`
- `pheno.py`
- `gwas.py`
- `anno.py`
- `qtl.py`
- `peak.py`
- `mr.py`
- `go.py`
- `net.py`
- `vis.py`
- `multi.py`

### Interface counts

- `src/server.py` currently exposes **81 MCP tools**
- `src/mrbigr_cli.py` currently exposes **55 CLI tools**

## Current Functional Highlights

### Annotation

The annotation stack now includes:

- GTF parsing
- SNP functional annotation
- interval-based gene lookup
- QTL annotation
- local annotation database generation
- variant effect prediction

### QTL

The QTL stack now includes:

- region detection from GWAS output
- lead SNP identification
- peak SNP identification
- genotype extraction for target regions
- haplotype summaries
- regional plotting
- QTL summaries
- QTL-to-gene mapping
- enrichment testing
- overlap comparison
- BED export

### MR

The MR stack includes:

- causal effect estimation
- heterogeneity testing
- pleiotropy testing
- QTL target analysis
- QTL formatting for MR workflows

### GO / KEGG / GSEA

The GO stack supports:

- GO enrichment
- KEGG enrichment
- GSEA
- local and online modes
- GO term simplification
- report export
- configurable GO and GSEA plotting

### Network

The network stack supports:

- edge weighting
- module detection
- hub identification
- optional Java-backed ClusterONE execution
- a NetworkX fallback when Java is not available

## Documentation Cleanup in v4

The documentation set has been aligned for the current repository state:

- repository URL standardized to `https://github.com/CrazyHsu/MRBIGR2`
- author standardized to `CrazyHsu`
- contact email standardized to `crazyhsu9527@gmail.com`
- MCP tool count standardized to `81`
- default install environment standardized to `mrbigr2`
- prompt documentation consolidated into a single file
- legacy Chinese-only prompt duplication removed

## Recent Behavioral Clarifications

### GO plotting

- local GO terms that cannot be resolved through `go-basic.obo` are now intended to be filtered out rather than shown as synthetic `LOCAL` ontology plots
- GO plotting supports both split-ontology and combined-ontology outputs
- GSEA plotting can produce standard enrichment curves, trace plots, or NES summary barplots

### GWAS

- GWAS plotting remains part of the reporting workflow, even when analysis and plotting are triggered separately
- LMM should generally be preferred for structured maize populations

### Network

- ClusterONE and the NetworkX fallback should not be assumed to return identical module results
- Java remains optional at install time and can be added with `./install.sh --with-java`

## Compatibility Notes

- The MCP layer is the primary structured interface for interactive agent workflows
- The CLI layer remains useful for file-oriented batch usage
- The two interfaces overlap substantially, but they do not expose exactly the same surface area in every module

## Author

- **CrazyHsu**
- **Email:** `crazyhsu9527@gmail.com`
