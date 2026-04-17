---
name: causal-network
description: >
  Build a directed causal network via Mendelian randomization — from an
  exposure (SNP / QTL / gene expression) to a set of outcomes (metabolites,
  traits, or other genes). Triggers on user intent like "run MR",
  "is gene X causal for trait Y?", "build a causal network around this locus".
tags: [mr, causal, network, interactive]
mcps: [gwas_mcp, geno_mcp, mr_mcp, net_mcp, vis_mcp]
---

# causal-network

Multi-step MR workflow that starts from GWAS summary statistics (from
`gwas-pipeline` or user-supplied), clumps to independent instruments,
estimates causal effects with IVW + MR-Egger, tests for
heterogeneity / pleiotropy, and assembles the surviving edges into a
directed network with hub identification.

## When to use

Trigger on:
- "Is *P1* causal for these flavonoids?"
- "Run two-sample MR from this gene's expression to yield traits."
- "Build a causal network around this locus."
- "Which genes are downstream targets of locus X?"

Do **not** use this skill for:
- Finding the GWAS peaks in the first place (use `gwas-pipeline`).
- Prioritizing candidate genes under a QTL with VEP + haplotype
  (use `qtl-to-target`).

## Inputs

| Key | Required | Default | Note |
|---|---|---|---|
| `exposure` | yes | — | one of: GWAS summary file, gene expression column name, or QTL bed path |
| `outcome_set` | yes | — | matrix (CSV) or list of GWAS summaries for outcomes |
| `clump_params` | no | r²=0.1, kb=500 | for independent instruments |
| `mr_p_threshold` | no | 0.05 | FDR adjusted reported separately |
| `build_network` | no | true | if false, stop after per-pair MR |
| `output_dir` | no | `output/mr_<timestamp>` | |

## Interaction checkpoints

1. **Instrument choice.** "Use all GWAS-significant SNPs after clumping, or
   restrict to cis-eQTLs within N kb of the gene?" (default: clumped
   significant SNPs at the user's threshold)
2. **MR method.** "IVW as primary, MR-Egger + weighted median as sensitivity
   — OK?"
3. **Sensitivity handling.** "Drop pairs with Egger-intercept P < 0.05 (likely
   pleiotropy) from the network?" (default: yes, drop and log)
4. **Network building.** "Min causal P for an edge in the network (default
   1e-5)? Direction: exposure → outcome only, or also reverse MR?"
5. **Hub identification.** "Call a node a hub at degree ≥ K (default K = 5)?"

## Steps

1. **Instrument building**
   - If starting from a QTL bed: `format_qtl_for_mr(qtl_file, output_file)`.
   - Otherwise: `generate_clump` or `run_snp_clumping` (geno_mcp) on the
     exposure GWAS with user-confirmed r²/kb.
   - Report: # instruments, mean F-statistic.

2. **MR pair-wise estimation**
   - For each outcome: `run_mr_analysis(exposure_gwas, outcome_gwas, method=ivw)`.
   - Also `calculate_mr_causal_estimate` for MR-Egger and weighted median as
     sensitivity.
   - Record β, SE, P, method, N-SNPs per pair.

3. **Sensitivity tests**
   - `test_mr_pleiotropy(exposure_gwas, outcome_gwas)` — Egger intercept.
   - `test_mr_heterogeneity(exposure_gwas, outcome_gwas)` — Cochran's Q.
   - Mark pairs failing pleiotropy; by default drop them from the network.

4. **Edge filtering**
   - Keep pairs with IVW P < user-chosen threshold AND pleiotropy OK.
   - Write filtered edge table with sign(β) and |β|.

5. **Optional gene-scale sweep**
   - `run_qtl_target_analysis(qtl_df, tf_genes, target_genes, window)` —
     use this when the user wants a one-call screen over many targets.

6. **Network assembly + hub ID**
   - Load edges into NetworkX-backed network builder.
   - `module_identify(edges)` (net_mcp) — overlapping module detection.
   - `hub_identify(edges, degree_threshold)` — report hubs.

7. **Visualization + report**
   - Forest plot of top causal pairs (use `plot_qtl_boxplot` style helpers
     or a custom plot emitted alongside the edge CSV).
   - Network PNG with hubs highlighted.
   - Table: exposure, outcome, β_IVW, P, SNPs, heterogeneity Q-P,
     pleiotropy P, included_in_network.

## Outputs

```
<output_dir>/
├── instruments.csv
├── mr_pairs.csv                 # all pairs: β, SE, P, methods
├── mr_edges.filtered.csv        # kept pairs
├── figures/
│   ├── forest_top20.png
│   └── network.png
├── hubs.csv
└── run_summary.json
```

## Failure modes

| Condition | Agent behaviour |
|---|---|
| < 3 independent instruments after clumping | Switch to single-SNP Wald-ratio; clearly flag as low-confidence. |
| All outcomes non-significant | Report the null result; ask if user wants to relax threshold or review instrument strength. |
| Pleiotropy P < 0.05 on most pairs | Stop; ask user to review instrument selection (often upstream GWAS problem). |
| Outcome sample overlap with exposure > 50 % and two-sample MR assumed | Warn: one-sample-bias risk; suggest split-sample or LOCO. |
| Network has > 1000 edges | Ask user if they want to tighten P threshold before rendering. |

## Example dialog

```
USER:  Given the GWAS output from the last skill run on the flavonoid
       metabolites, is P1 expression causally linked to those metabolites?
       Also check gene-to-gene MR downstream of P1.

AGENT: Using data/gwas/flavonoid_*.assoc.txt as exposure-side GWAS and
       data/expression/AMP_kernel_fpkm.csv for P1 expression. Clumping
       at r²=0.1, 500 kb. MR will be IVW primary + MR-Egger + weighted
       median as sensitivity; edges with Egger-intercept P < 0.05 will
       be dropped. OK?

USER:  Yes.

AGENT: [run_snp_clumping → format_qtl_for_mr → run_mr_analysis ×102
       outcomes → test_mr_pleiotropy → test_mr_heterogeneity → filter]
       P1 exposure: 30/102 metabolites significant at P < 0.01, 18/102
       at P < 1e-5. Gene-to-gene sweep: 568 candidates nominal, 56 at
       P < 1e-5. Network: 67 nodes, 94 edges, P1 is the top hub
       (degree = 23). Network PNG saved. Want me to enrich the
       56 downstream genes (functional-enrichment)?
```
