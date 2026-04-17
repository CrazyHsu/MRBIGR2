---
name: qtl-to-target
description: >
  Rank candidate causal genes under a QTL region using variant-effect
  prediction, gene–trait MR, and haplotype testing. Triggers on intent like
  "prioritize candidate genes for this QTL", "map the lead SNP to function",
  "which gene under this peak actually drives the phenotype?".
tags: [qtl, annotation, haplotype, prioritization, interactive]
mcps: [qtl_mcp, anno_mcp, geno_mcp, peak_mcp, mr_mcp, vis_mcp]
---

# qtl-to-target

Post-GWAS prioritization skill. Given a QTL region (bed) or a lead SNP and
flanking window, walks through: genes in region → variant-effect prediction
on the lead SNP(s) → gene-to-trait MR via expression → haplotype effect on
trait → locus visualization. Produces a ranked candidate table.

## When to use

Trigger on:
- "Which gene under the chr1 48 Mb peak is the causal one?"
- "Prioritize candidate genes for this QTL."
- "Annotate the lead SNP and test haplotype effects."

Do **not** use this skill:
- To discover peaks in the first place — use `gwas-pipeline`.
- For genome-wide MR sweeps — use `causal-network`.

## Inputs

| Key | Required | Default | Note |
|---|---|---|---|
| `region` | yes | — | either `qtl_bed` path OR (`chrom, start, end`) OR `lead_snp` |
| `flanking_kb` | no | 100 | used when only a lead SNP is given |
| `trait` | yes | — | phenotype / metabolite / trait column the region was called on |
| `expression_matrix` | no | — | gene expression (needed for gene-expr → trait MR) |
| `genotype_prefix` | yes | — | for haplotype extraction + VEP context |
| `gtf_file` | yes | — | gene coordinates |
| `anno_db` | no | auto | variant-effect annotation DB |
| `output_dir` | no | `output/qtl_target_<timestamp>` | |

## Interaction checkpoints

1. **Region confirmation.** "Using region chr<X>:<start>-<end> (N genes
   inside). Keep flanking ±100 kb? Or tighten?"
2. **VEP scope.** "Include all coding and non-coding consequence classes,
   or only missense / splice / start-lost / stop-gain?"
3. **MR direction.** "Test QTL → gene expression → trait (mediation-style),
   or only cis-eQTL → trait (SMR-style)?"
4. **Haplotype test.** "Minimum N per haplotype group (default 20) — drop
   rare haplotypes below this?"
5. **Ranking weights.** "Rank by composite score of (VEP severity × MR P
   × haplotype P), or let you inspect each column and decide yourself?"
   (default: ranked + full columns shown)

## Steps

1. **Resolve region**
   - If `qtl_bed`: read bed. Else build a ±flanking window around
     `lead_snp` via `get_lead_snp` if needed.
   - `get_genes_in_region(chrom, start, end, gtf)` → gene list (ID, name,
     biotype, TSS dist to lead SNP).

2. **Genotypes + annotation DB**
   - `extract_qtl_genotypes(geno_prefix, region, output_prefix)` —
     PLINK-format subset.
   - `create_annotation_db(gtf_file, output_db)` once (cached per GTF).
   - `qtl_annotation(qtl_df, anno_db)` — overlap QTL with gene features.

3. **Variant-effect prediction**
   - `predict_variant_effect(vcf_or_plink, anno_db, snp_ids=[lead_snp, ...])`
     — per-SNP consequence + impact class.
   - Also run `annotate_snps` on the top-10 LD partners of the lead SNP
     to widen the functional footprint.

4. **Gene–trait MR (per candidate gene)**
   - For each gene in the region, if `expression_matrix` provided:
     - Build a gene-specific instrument (cis-eQTL within region) from
       `gwas_mcp` outputs if available, else use `run_simple_gwas` on the
       gene expression column with the region's genotypes.
     - `run_qtl_target_analysis(qtl_df, tf_genes=[gene], target_genes=[trait])`
       — returns β, P for the gene → trait causal estimate.

5. **Haplotype effect**
   - `calculate_qtl_haplo(geno_prefix, qtl_snp_list)` — haplotype groups.
   - `haplotype_test(pheno_df, geno_df, snp_id=lead_snp, test_method=t-test)`
     (peak_mcp) — effect of each haplotype on trait.
   - `plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df)` per top haplotype.

6. **Visualization + ranking**
   - `plot_qtl_region(gwas_file, chr_val, start, end, output_file)` —
     regional Manhattan with genes.
   - `plot_ld_heatmap(geno_prefix, region)` — LD structure.
   - Assemble ranking table: gene, distance_to_lead, VEP_impact, MR_β,
     MR_P, haplotype_P, composite_score.

## Outputs

```
<output_dir>/
├── genes_in_region.csv
├── vep.csv                      # consequence per SNP
├── gene_trait_mr.csv            # β, P per candidate gene
├── haplotype_test.csv
├── ranking.csv                  # final prioritized list
├── figures/
│   ├── qtl_region.png
│   ├── ld_heatmap.png
│   └── haplotype_box_<gene>.png
└── run_summary.json
```

## Failure modes

| Condition | Agent behaviour |
|---|---|
| Region has > 200 genes | Ask user to tighten window before proceeding. |
| Lead SNP has no annotation in GTF range | Expand window by 50 kb once; if still empty, report and stop. |
| Expression matrix missing for the gene-trait MR step | Skip MR column; report "no eQTL evidence"; keep VEP + haplotype. |
| VEP returns all `intergenic` | Suggest widening flank; flag the region as non-coding-driven. |
| Haplotype test groups each < 20 samples | Merge groups or skip; report sample sizes. |

## Example dialog

```
USER:  For the chr1:47.5–49.5 Mb peak from the flavonoid GWAS,
       which gene is the most likely causal driver?

AGENT: 52 protein-coding genes in region. I'll:
       1. VEP on the lead SNP (chr1.s_48424403) and top 10 LD partners
       2. For each gene with a cis-eQTL within the region, test
          gene-expression → flavonoid-composite via MR
       3. Haplotype test on the lead SNP vs trait
       4. Rank by composite (VEP impact × MR P × haplotype P).
       OK to run?

USER:  Yes, but only missense / splice / regulatory for VEP.

AGENT: [get_genes_in_region → predict_variant_effect (scope=filtered)
       → run_qtl_target_analysis ×52 → calculate_qtl_haplo →
       haplotype_test → plot_qtl_region]
       Top candidate: Zm00001d028854 (P1 / MYB12).
         - lead SNP: upstream (165 kb); no missense nearby
         - gene-expr → trait MR β = 1.8, P = 3.2e-12
         - haplotype test: AA vs GG, P = 4.4e-18
       Second: Zm00001d028851 (distance 180 kb; MR P = 3e-4).
       Rank CSV and 3 plots saved to output/.
```
