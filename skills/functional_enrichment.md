---
name: functional-enrichment
description: >
  Run GO / KEGG over-representation and GSEA on a gene list or ranked list,
  simplify redundant terms, and produce a report. Triggers on intent like
  "what pathways are these genes enriched in?", "run GSEA on this MR ranking",
  "GO enrichment of my candidate genes".
tags: [go, kegg, gsea, enrichment, interactive]
mcps: [go_mcp, anno_mcp]
---

# functional-enrichment

Takes a gene list (foreground + optional background) or a ranked list
(e.g., by signed MR statistic) and produces GO / KEGG / GSEA results with
redundancy-simplified terms, enrichment plots, and an HTML report.

## When to use

Trigger on:
- "What pathways are these genes enriched in?"
- "GO enrichment on this list."
- "Run GSEA with the MR-ranked genes."
- "Find the biological processes over-represented here."

Do **not** use this skill to:
- Discover genes in the first place (use `gwas-pipeline` or `qtl-to-target`).
- Build a causal network (use `causal-network`).

## Inputs

| Key | Required | Default | Note |
|---|---|---|---|
| `gene_list` | yes (one of) | — | unranked list (CSV or .txt, one gene per line) |
| `ranked_list` | yes (one of) | — | CSV with gene + numeric ranking column |
| `ranking_column` | if ranked | auto-detect | column name used for GSEA ranking |
| `background` | no | all annotated genes | foreground universe for over-representation |
| `organism` | no | `maize` | used to pick GO DB / KEGG code |
| `gtf_file` | no | — | needed to build a GO DB if none exists locally |
| `go_db` | no | auto | pre-built GO mapping; if missing, trigger `extract_go_from_gtf` |
| `p_threshold` | no | 0.05 | adjusted P (BH) |
| `output_dir` | no | `output/enrich_<timestamp>` | |

## Interaction checkpoints

1. **Mode.** If the input is a CSV with a numeric column, ask:
   "Treat this as a **ranked** list for GSEA (ranking by <column>), or as an
   **unranked** foreground for over-representation? Default: ranked if a
   numeric column is present."
2. **Organism and ID system.** "Using maize (Zm_00001d*) — correct? Or is
   your list in Ensembl / RefSeq IDs?" If ID mapping hit-rate is low
   (< 30 %), stop and confirm.
3. **Databases.** "Run GO (BP / MF / CC), KEGG, and GSEA — all three? Or
   only a subset?"
4. **Simplification.** "Simplify redundant GO terms (similarity > 0.7)?"
   (default: yes)
5. **Background.** "Use all protein-coding genes in the GTF as background,
   or a custom expressed-gene background?"

## Steps

1. **DB bootstrap**
   - If `go_db` not found: `extract_go_from_gtf(gtf_file)` once, cached.
   - `get_enrichr_libraries` — list available Enrichr libraries for the
     organism so the user can pick if they want Enrichr-backed enrichment.

2. **Gene-ID sanity check**
   - Map the user's list to the GO DB's gene-ID space; report hit rate.
   - If hit rate < 30 %, stop and ask.

3. **Unranked path** (gene_list)
   - `run_go_enrichment(gene_list, organism, mode=local|enrichr, p_threshold)`
     → GO BP / MF / CC tables.
   - `run_kegg_enrichment(gene_list, organism)` → KEGG pathway table.

4. **Ranked path** (ranked_list)
   - `run_gsea_analysis(ranked_list, ranking_column, organism, permutations)`
     → GSEA result table + enrichment curves.

5. **Simplification + plots**
   - `simplify_go_results(go_table, similarity_threshold=0.7)`.
   - `plot_go_enrichment` for the top-N bubble / bar chart.
   - `plot_gsea_results` for the leading-edge curve on the top GSEA set.

6. **Report**
   - `export_go_report(output_dir)` — single HTML index pointing at all
     tables and plots.
   - Short prose summary in the agent's reply: top 5 processes + top
     KEGG pathway + most significant GSEA set.

## Outputs

```
<output_dir>/
├── go_bp.csv
├── go_mf.csv
├── go_cc.csv
├── kegg.csv
├── gsea.csv                     # if ranked
├── go_simplified.csv
├── figures/
│   ├── go_bubble.png
│   ├── gsea_<top_term>.png
│   └── kegg_barplot.png
├── report.html
└── run_summary.json
```

## Failure modes

| Condition | Agent behaviour |
|---|---|
| ID-mapping hit rate < 30 % | Stop; ask user which ID system the list uses. |
| Gene list < 5 | Refuse over-representation; propose GSEA only if a ranking is available. |
| No significant terms at P_adj < 0.05 | Report honestly; offer to relax to nominal P ≤ 0.05 OR inspect top-10 terms without a cutoff. |
| GO DB missing and no GTF provided | Offer to download the organism's go-basic.obo + species-GAF from a path the user supplies; do not silently use a different organism. |
| Permutations < 100 for GSEA | Warn; default to 1000. |

## Example dialog

```
USER:  These are the 56 genes downstream of P1 (output/mr_edges.filtered.csv,
       "outcome" column) — what biological processes are they enriched in?

AGENT: 56 genes is a reasonable foreground. Options:
       (a) Unranked GO + KEGG over-representation against all annotated
           maize genes as background.
       (b) GSEA using the MR beta column as ranking (more powerful if you
           have the full 568-gene list, not just the 56).
       Which path? (Default (b) on the full table, or (a) on the filtered
       56 list?)

USER:  Do (b) with the full 568-gene list, ranked by signed MR beta.

AGENT: [extract_go_from_gtf (cached) → run_gsea_analysis(ranking=mr_beta,
       perm=1000) → simplify_go_results → plot_gsea_results → export_go_report]
       Top enriched set: flavonoid biosynthetic process
       (NES = 6.3, P = 1e-3, leading-edge 12 genes including C2, CGT1,
       UGT1, PR1, FNS1). Phenylpropanoid biosynthetic process and
       anthocyanin-containing compound biosynthesis also significant.
       Report: output/enrich_<ts>/report.html.
```
