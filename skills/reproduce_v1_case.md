---
name: reproduce-v1-case
description: >
  Meta-skill that chains gwas-pipeline → qtl-to-target → causal-network →
  functional-enrichment on a user-supplied multi-omics path, for one-shot
  reproduction of published MRBIGR V1 case-study flows. Triggers on intent
  like "reproduce the P1–flavonoid case", "run the V1 pipeline end-to-end",
  "重跑一下 V1 那个 P1 案例".
tags: [meta, reproduction, case-study, interactive]
mcps: [geno_mcp, pheno_mcp, gwas_mcp, mr_mcp, go_mcp, qtl_mcp, anno_mcp, net_mcp, peak_mcp, vis_mcp]
---

# reproduce-v1-case

A meta-skill that orchestrates the four atomic skills to reproduce a complete
V1-style case study (discover peak → prioritize gene → build causal network
→ enrich downstream pathways) in one conversation.

Unlike the atomic skills, this one is **deliberately noisy** — it reports
every intermediate finding back to the user and stops for confirmation at
three hard gates. Designed for manuscript case studies where we want the
conversation itself to be the artifact.

## When to use

Trigger on:
- "Reproduce the P1–flavonoid case on the Maizego data at <path>."
- "Run the V1 MRBIGR analysis end-to-end on my data."
- "从 GWAS 走到通路富集，重现 V1 的 P1 分析。"

Use the atomic skills (`gwas-pipeline` / `causal-network` / …) instead when:
- The user wants one stage only.
- The user has already run GWAS elsewhere.

## Inputs

| Key | Required | Default | Note |
|---|---|---|---|
| `data_root` | yes | — | folder with genotype, phenotype, expression, GTF |
| `genotype_prefix` | auto / yes | auto-discover from `data_root` | |
| `phenotype_file` | yes | — | primary trait CSV (metabolite / agronomic) |
| `expression_file` | no | — | enables gene-to-trait MR |
| `gtf_file` | yes | — | |
| `trait_subset` | no | all | a list of target columns |
| `output_dir` | no | `output/case_<timestamp>` | |

## Interaction checkpoints (three hard gates)

**Gate 1 — after GWAS.**
> Agent: "Here are the peaks I found: <table>. Which region do you want me
> to follow up on? (Default: the most significant shared across traits.)"

**Gate 2 — after region prioritization.**
> Agent: "Top candidate gene in the region is <gene> (composite score X).
> Want me to use this gene as exposure for the causal-network step, or pick
> a different gene from the ranking?"

**Gate 3 — before enrichment.**
> Agent: "Network built: N nodes, M edges, hubs = [...]. Should I enrich
> (a) only the hub's direct neighbours, (b) all outcomes at P < 1e-5, or
> (c) the full ranked list by MR beta via GSEA?" (default c)

At each gate, the agent waits and does not proceed without a reply.

## Steps

1. **Discover data layout** — list files in `data_root`, infer genotype
   prefix, phenotype matrix, expression matrix, GTF. Print the plan.

2. **Stage 1: `gwas-pipeline`** — inherits the atomic skill's defaults,
   running on `trait_subset`. Emits per-trait Manhattan + QQ + QTL bed.
   **Stop at Gate 1.**

3. **Stage 2: `qtl-to-target`** — run on the region confirmed at Gate 1.
   Emits candidate ranking and regional visualization. **Stop at Gate 2.**

4. **Stage 3: `causal-network`** — exposure = gene confirmed at Gate 2
   (expression column from the provided matrix). Outcomes = `trait_subset`
   and/or all outcomes in `phenotype_file` (ask the user if ambiguous).
   Emits edge table, forest, network, hubs. **Stop at Gate 3.**

5. **Stage 4: `functional-enrichment`** — gene list per Gate 3 choice.
   Emits GO / KEGG / GSEA + HTML report.

6. **Final report** — a `case_report.md` file that concatenates the
   decisions made at each gate, all tool calls, all figures, and a short
   prose summary in the user's language.

## Outputs

```
<output_dir>/
├── stage1_gwas/                 # same layout as gwas-pipeline
├── stage2_qtl_target/           # same layout as qtl-to-target
├── stage3_causal_network/       # same layout as causal-network
├── stage4_enrichment/           # same layout as functional-enrichment
├── case_report.md               # Gate-by-Gate narrative + tool log
├── chat_log.md                  # full conversation (appended live)
└── tool_trace.jsonl             # one line per MCP tool call
```

`chat_log.md` and `tool_trace.jsonl` are the **primary deliverables** for
the manuscript's Supplementary Notes and Supplementary Table S3.

## Failure modes

| Condition | Agent behaviour |
|---|---|
| Data-root discovery finds no genotype / phenotype | Stop; ask user to supply explicit paths. |
| Stage-1 yields 0 QTL at default thresholds | Relax once (p1=1e-6 → 1e-5), flag it; if still 0, stop and ask. |
| Gate response unclear (e.g., user says "都行" / "whatever") | Offer the default explicitly and wait for explicit yes. |
| Any atomic skill errors | Surface the error + its tool-call trace; do not swallow it; do not auto-retry with changed parameters without asking. |

## Example dialog

```
USER:  Reproduce the P1–flavonoid case. Data is at
       /mnt/e/工作项目/wanglab/MRBIGR_data/MRBIGR_data.20241019/ .
       Use the 14 flavonoid metabolites.

AGENT: Discovered:
       - genotype: chr_HAMP.{bed,bim,fam}
       - phenotype: 102_metabolites.E3_log2.normalized_phe.csv
       - expression: AMP_kernel_transcriptome_v4_FPKM_correct.qqnorm...csv
       - GTF: Zea_mays.B73_RefGen_v4.50.gtf.gz
       Plan: gwas-pipeline → qtl-to-target → causal-network →
       functional-enrichment, with stops at three gates. OK?

USER:  Go.

AGENT: [Stage 1 ... 14 traits LMM ... 28 min]
       GATE 1: One shared peak, chr1:47.5–49.5 Mb (lead chr1.s_48424403,
       P=9.9e-30). Proceed to prioritize this region?
USER:  Yes.

AGENT: [Stage 2 ...]
       GATE 2: Top candidate = Zm00001d028854 (P1 / MYB12). Use P1
       expression as exposure for Stage 3?
USER:  Yes.

AGENT: [Stage 3 ... 568 pairs tested, 56 at P<1e-5, P1 is hub]
       GATE 3: Enrichment on (a) direct neighbours (23 genes), (b) 56
       at P<1e-5, (c) GSEA on full 568-gene ranked list. Default (c).
USER:  (c).

AGENT: [Stage 4 ...]
       Flavonoid biosynthetic process enriched NES=6.3, P=1e-3.
       Full case_report.md + chat_log.md + tool_trace.jsonl saved to
       output/case_<ts>/.
```
