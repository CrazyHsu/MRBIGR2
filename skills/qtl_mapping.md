---
name: qtl-mapping-deprecated
description: Deprecated — superseded by qtl-to-target (qtl_to_target.md).
tags: [deprecated]
---

# (deprecated) QTL Mapping

This skill has been replaced by **`qtl-to-target`**, which adds
variant-effect prediction, gene-to-trait MR, composite candidate-gene
ranking, and three hard interaction checkpoints on top of the original
QTL detection + haplotype flow.

See [`qtl_to_target.md`](./qtl_to_target.md).

Also note: QTL region *discovery* (detect_qtl_regions, identify_peak_snps,
etc.) is now part of the unified **`gwas-pipeline`** skill —
[`gwas_pipeline.md`](./gwas_pipeline.md) — rather than a separate
post-GWAS skill.

Kept as a thin redirect to avoid breaking references in old notebooks / docs.
