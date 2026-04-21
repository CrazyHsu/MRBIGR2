---
name: gwas-pipeline
description: >
  Run a GWAS end-to-end from genotype + phenotype to QTL regions and gene
  mapping, in an interactive dialog. Triggers on user intent like "run GWAS",
  "find peaks for this trait", "do association analysis on these metabolites".
tags: [gwas, qtl, interactive]
mcps: [geno_mcp, pheno_mcp, gwas_mcp, vis_mcp, qtl_mcp, anno_mcp]
---

# gwas-pipeline

End-to-end GWAS as a dialog. The agent drives QC → kinship → PCA →
association → diagnostics → QTL detection → gene mapping, and **pauses at
defined checkpoints to confirm key parameters with the user**.

## When to use

Trigger on user intent such as:
- "run GWAS on these flavonoid metabolites"
- "find significant loci for 100-grain weight"
- "从基因型和表型开始做关联分析"
- "map QTL for this trait"

Do **not** use this skill when:
- The user already has GWAS summary statistics and only wants MR / enrichment
  (use `causal-network` or `functional-enrichment`).
- The user only wants post-GWAS QTL prioritization on a known region
  (use `qtl-to-target`).

## Inputs

| Key | Required | Default | Note |
|---|---|---|---|
| `genotype_prefix` | yes | — | PLINK bed/bim/fam prefix or VCF path |
| `phenotype_file` | yes | — | CSV, first column = sample ID |
| `traits` | no | all columns | list of phenotype columns |
| `model` | no | `lmm` | `lmm` \| `lm` \| `plink` |
| `gtf_file` | no | — | needed for gene mapping at the last step |
| `output_dir` | no | `output/gwas_<timestamp>` | all artifacts land here |

If `genotype_prefix` is a VCF, the skill first converts it via `convert_vcf`.
If the VCF lacks rs IDs, `ensure_rs_id` / `add_rs_id_to_vcf` is run first.

## Interaction checkpoints

The agent **must** ask before running the heavy steps:

1. **QC thresholds.** "Use MAF ≥ 0.05 and SNP missing rate ≤ 0.1? Or should I
   change them?"
2. **Phenotype normalization.** "Do you want me to remove outliers (z-score
   |z| > 3) and scale to z-score before GWAS?"
3. **PCA components.** "Use top 5 PCs as covariates, or a different number?"
4. **Model.** "LMM (GEMMA, accounts for kinship) or PLINK logistic/linear?
   LMM is default." — skip if `model` already specified.
5. **Significance threshold for QTL detection.** "Bonferroni (P = 0.05 / N)
   or a fixed threshold like 1e-6?"
6. **Gene-mapping flanking window.** "Extend QTL regions by ±100 kb when
   mapping to genes?"

The agent announces the concrete numbers chosen before each run.

## Steps

Tool names below are the actual `@mcp.tool()` registrations. Inputs/outputs
are kept as-is from the MCP schema; do not rename.

**Hard-stop rule:** if any MCP tool call returns `is_error=true`, stop the
workflow and report the failing tool, exact input class/shape/path, and error
message to the user. Do **not** use Grep, Read, Bash, Python snippets, or source
inspection to self-diagnose MCP internals after a tool error.

1. **Genotype ingest + QC**
   - If VCF: `convert_vcf` → PLINK prefix.
   - `run_snp_qc(input_prefix, maf, missing_rate, output_prefix)`
   - `get_snp_statistics` for a quick summary (report to user).

2. **Optional LD pruning / imputation**
   - If user asked for imputation OR post-QC missing > 2 %: `impute_genotype`.
   - If the user wants independent SNPs (e.g., PCA only):
     `run_snp_pruning(window_kb, step, r2)`.

3. **Population structure**
   - `calculate_kinship` (for LMM).
   - `run_genotype_pca(n_components)` + `plot_pca` if user wants a quick view.

4. **Phenotype prep** (only if user agreed to normalize)
   - `filter_missing(missing_ratio)` → `remove_outliers(method=zscore)` →
     `scale_phenotype(method=zscore)`.
   - Keep the returned phenotype object in memory and pass that object to GWAS.
     Do not create a custom `pheno_clean.csv` with Bash/Python as a substitute
     for the phenotype MCP calls.
   - Optional quick diagnostics: `plot_phenotype_hist`,
     `plot_phenotype_correlation` on the final matrix.

5. **Association**
   - `run_gwas_lmm(phe, geno_prefix, output_dir, auto_plot=True)`, where `phe`
     is either the phenotype object returned from step 4 or the original
     phenotype CSV path when normalization was explicitly skipped.
     (or `run_gwas_lm` / `run_gwas_plink` based on model choice).
   - After each trait: `calculate_lambda` — if λ ∉ [0.95, 1.15], flag it.

6. **Diagnostics and summary**
   - `plot_manhattan` + `plot_qq` per trait (if `auto_plot=False` above).
   - `gwas_summary` → one-row-per-trait summary CSV for the final report.

7. **QTL detection**
   - `detect_qtl_regions(gwas_file, p1, p2, window)` with user-confirmed p1/p2.
   - `get_lead_snp` + `identify_peak_snps` for each detected region.
   - `export_qtl_bed` for downstream use.

8. **Gene mapping** (if `gtf_file` provided)
   - `parse_gtf(gtf_file)` once (cached).
   - `map_qtl_to_genes(qtl_df, annotation_file, flanking_kb)`.

9. **Report back to user**
   - Table of regions × (lead SNP, P, λ, # genes in window).
   - Manhattan + QQ links per trait.
   - Ask: "Proceed to causal-network / qtl-to-target on one of these regions?"

## Outputs

```
<output_dir>/
├── qc/                         # post-QC PLINK files
├── kinship/                    # GRM matrix
├── pca/                        # PCs CSV + plot_pca PNG
├── pheno_clean.csv             # optional MCP-derived export; not a substitute for phenotype MCP calls
├── gwas/
│   ├── <trait>.assoc.txt       # raw GEMMA / PLINK output
│   ├── <trait>_manhattan.png
│   └── <trait>_qq.png
├── qtl/
│   ├── <trait>.bed
│   └── <trait>_genes.csv       # if gtf provided
└── run_summary.json            # tool calls, lambdas, timings
```

## Failure modes

| Condition | Agent behaviour |
|---|---|
| Sample ID overlap between genotype and phenotype < 50 | Stop; report overlap; ask user to reconcile IDs before proceeding. |
| λ > 1.5 for a trait | Pause; report inflation; ask whether to re-normalize or change model. |
| Kinship singular / GEMMA failure | Fall back to `run_gwas_lm` and flag the change in the report. |
| MAF threshold removes > 80 % of SNPs | Stop; ask user to confirm threshold or reconsider dataset. |
| No QTL region detected at user's P threshold | Report "0 regions"; offer to relax thresholds (do NOT silently relax). |

## Example dialog

```
USER:  I have the Maizego PLINK files at data/chr_HAMP.* and a
       flavonoid metabolite matrix at data/102_metabolites.csv.
       Run GWAS on 14 of the flavonoid columns listed in
       data/flavonoid_cols.txt. I have the v4 GTF at
       data/Zea_mays.B73_RefGen_v4.50.gtf.gz.

AGENT: Here's what I'll run:
       - geno: data/chr_HAMP  (PLINK)
       - pheno: data/102_metabolites.csv, 14 flavonoid columns
       - model: LMM via GEMMA (kinship + 5 PCs)
       - QC: MAF >= 0.05, missing <= 0.1
       - QTL detection: p1=1e-7, p2=1e-5, window=500kb
       - Gene mapping: ±100kb flanking
       Confirm, or change any of these?

USER:  Looks good, but use p1=1e-6.

AGENT: Running ... [calls run_snp_qc, calculate_kinship, run_genotype_pca,
       run_gwas_lmm per trait, detect_qtl_regions, map_qtl_to_genes ...]
       14 traits: mean λ=0.98 (range 0.94–1.03). 1 shared peak at
       chr1:47.5–49.5 Mb (lead SNP chr1.s_48424403, P=9.9e-30).
       52 genes in flank; P1/Zm00001d028854 is the closest TF.
       Want me to dig into this region (qtl-to-target) or run MR
       from P1 expression to all 102 metabolites (causal-network)?
```
