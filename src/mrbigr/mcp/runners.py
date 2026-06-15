"""Adapter functions invoked by ``mrbigr.mcp.longjob_worker`` for tools that
need Python-side data prep (DataFrames, CSV inputs, ...).

Each adapter takes only JSON-serializable kwargs (paths, scalars). It loads
DataFrames from disk on the worker side and calls the underlying
``mrbigr.core.*`` function. Keeps JobSpec JSON-clean.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional


def _ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def _write_json(path: str | Path, payload: dict) -> str:
    _ensure_parent(path)
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    return str(path)


def _write_frame(frame, output_file: str, *, sep: str = ",") -> str:
    _ensure_parent(output_file)
    if hasattr(frame, "to_csv"):
        frame.to_csv(output_file, index=False, sep=sep)
    else:
        Path(output_file).write_text(json.dumps(frame, ensure_ascii=False, default=str), encoding="utf-8")
    return output_file


def _coerce_frame(value):
    import pandas as pd
    from ..core._argjson import maybe_json_loads

    value = maybe_json_loads(value)
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, (str, Path)):
        return pd.read_csv(value, sep=None, engine="python")
    if isinstance(value, dict):
        return pd.DataFrame(value)
    return pd.DataFrame(value)


def _coerce_gene_list(value):
    import pandas as pd
    from ..core._argjson import maybe_json_loads

    value = maybe_json_loads(value)
    if isinstance(value, list):
        return pd.DataFrame({"genes": value})
    if isinstance(value, dict):
        return pd.DataFrame(value)
    return value


def _coerce_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _read_phe(phe_csv: str):
    import pandas as pd

    return pd.read_csv(phe_csv, index_col=0)


def _sample_match_summary(phe, geno_prefix: str) -> dict:
    """Summarize phenotype/FAM overlap before GEMMA rewrites the FAM file."""
    import pandas as pd

    fam_file = Path(str(geno_prefix) + ".fam")
    fam = pd.read_csv(
        fam_file,
        sep=r"\s+",
        header=None,
        names=["fam", "id", "pat", "mat", "sex", "pheno"],
    )
    fam_ids = set(fam["id"].astype(str))
    phe_index = phe.index.astype(str)
    phe_ids = set(phe_index)

    per_trait = {}
    for trait in phe.columns:
        nonmissing = phe[trait].notna()
        trait_ids = set(phe_index[nonmissing])
        per_trait[str(trait)] = {
            "phenotype_nonmissing": int(nonmissing.sum()),
            "matched_nonmissing": int(len(trait_ids & fam_ids)),
        }

    not_in_fam = sorted(phe_ids - fam_ids)
    not_in_phe = sorted(fam_ids - phe_ids)
    return {
        "phenotype_rows": int(len(phe)),
        "fam_rows": int(len(fam)),
        "traits": per_trait,
        "phenotype_ids_not_in_fam_count": int(len(not_in_fam)),
        "phenotype_ids_not_in_fam_first20": not_in_fam[:20],
        "fam_ids_not_in_phenotype_count": int(len(not_in_phe)),
        "fam_ids_not_in_phenotype_first20": not_in_phe[:20],
    }


def _auto_plot(result_files):
    from mrbigr.core import vis

    plots = []
    for result_file in result_files or []:
        if not result_file:
            continue
        stem = str(result_file).replace(".assoc.txt", "").replace(".assoc.linear", "")
        manhattan = vis.manhattan_plot(result_file, output_file=f"{stem}_manhattan.png")
        qq = vis.qq_plot(result_file, output_file=f"{stem}_qq.png")
        plots.append({"file": result_file, "manhattan": manhattan, "qq": qq})
    return plots


def gwas_lm_runner(
    phe_csv: str,
    geno_prefix: str,
    output_name: Optional[str] = None,
    output_dir: Optional[str] = None,
    num_threads: Optional[int] = None,
    auto_plot: bool = False,
):
    from mrbigr.core import gwas

    phe = _read_phe(phe_csv)
    result_files = gwas.gwas_lm(
        phe, geno_prefix,
        output_name=output_name, output_dir=output_dir, num_threads=num_threads,
    )
    if _coerce_bool(auto_plot):
        return {"result_files": result_files, "plots": _auto_plot(result_files)}
    return result_files


def gwas_lmm_runner(
    phe_csv: str,
    geno_prefix: str,
    output_name: Optional[str] = None,
    output_dir: Optional[str] = None,
    num_threads: Optional[int] = None,
    cov: Any = None,
    kinship_file: Optional[str] = None,
    auto_plot: bool = False,
):
    from mrbigr.core import gwas

    phe = _read_phe(phe_csv)
    sample_summary = _sample_match_summary(phe, geno_prefix)
    if isinstance(cov, dict) and cov.get("__json_file__"):
        import json as _json
        cov = _json.loads(Path(cov["__json_file__"]).read_text(encoding="utf-8"))
    result_files = gwas.gwas_lmm(
        phe, geno_prefix,
        output_name=output_name, output_dir=output_dir, num_threads=num_threads,
        cov=cov, kinship_file=kinship_file,
    )
    if output_dir is not None:
        summary_stem = output_name or "gwas_lmm"
        sample_summary_file = Path(output_dir) / f"{summary_stem}.sample_summary.json"
        _write_json(sample_summary_file, sample_summary)
    else:
        sample_summary_file = None
    if _coerce_bool(auto_plot):
        return {
            "result_files": result_files,
            "plots": _auto_plot(result_files),
            "sample_summary": sample_summary,
            "sample_summary_file": str(sample_summary_file) if sample_summary_file else None,
        }
    return {
        "result_files": result_files,
        "sample_summary": sample_summary,
        "sample_summary_file": str(sample_summary_file) if sample_summary_file else None,
    }


def gwas_plink_runner(
    phe_csv: str,
    geno_prefix: str,
    pheno_col: Optional[str] = None,
    output_name: str = "gwas",
    auto_plot: bool = False,
):
    from mrbigr.core import gwas
    import pandas as pd

    phe = pd.read_csv(phe_csv, index_col=0)
    result_file = gwas.gwas_plink(
        phe, geno_prefix, pheno_col=pheno_col, output_name=output_name,
    )
    result_files = [result_file] if result_file else []
    if _coerce_bool(auto_plot):
        return {"result_files": result_files, "plots": _auto_plot(result_files)}
    return result_file


def gwas_clump_runner(
    geno_prefix: str,
    p1: float = 0.001,
    p2: float = 0.05,
    clump_input_dir: Optional[str] = None,
    result_dir: Optional[str] = None,
    clump_kb: int = 500,
    clump_r2: float = 0.1,
    result_json: Optional[str] = None,
):
    from mrbigr.core import gwas

    out_dir = gwas.gwas_clump(
        geno_prefix,
        p1=float(p1), p2=float(p2),
        clump_input_dir=clump_input_dir, result_dir=result_dir,
        clump_kb=int(clump_kb), clump_r2=float(clump_r2),
    )
    out_path = Path(out_dir).resolve()
    clumped = sorted(str(path) for path in out_path.glob("*.clumped") if path.is_file())
    payload = {
        "result_dir": str(out_path),
        "clumped_files": clumped,
        "n_clumped": len(clumped),
    }
    if result_json:
        return _write_json(result_json, payload)
    return payload


def mr_analysis_runner(
    exposure_gwas_csv: str,
    outcome_gwas_csv: str,
    output_csv: str,
    method: str = "ivw",
):
    """Run mr.mr_analysis with CSVs on disk; persist the result DataFrame to
    ``output_csv`` so completion can be detected by file_nonempty marker."""
    from mrbigr.core import mr
    import pandas as pd

    exp_df = pd.read_csv(exposure_gwas_csv, sep=None, engine="python")
    out_df = pd.read_csv(outcome_gwas_csv, sep=None, engine="python")
    result = mr.mr_analysis(exp_df, out_df, method=method)
    if result is None:
        return None
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    if hasattr(result, "to_csv"):
        result.to_csv(output_csv, index=False)
    else:
        # Allow dict / scalar returns; serialize as JSON next to the CSV.
        import json as _json
        Path(output_csv).write_text(_json.dumps(result, default=str), encoding="utf-8")
    return output_csv


def qtl_target_runner(
    qtl_csv: str,
    tf_genes: list,
    target_genes: list,
    output_csv: str,
    window: int = 500000,
):
    from mrbigr.core import mr
    import pandas as pd

    qtl_df = pd.read_csv(qtl_csv, sep=None, engine="python")
    result = mr.qtl_target_analysis(qtl_df, tf_genes, target_genes, window=int(window))
    if result is None:
        return None
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    if hasattr(result, "to_csv"):
        result.to_csv(output_csv, index=False)
    else:
        import json as _json
        Path(output_csv).write_text(_json.dumps(result, default=str), encoding="utf-8")
    return output_csv


def net_module_identify_runner(edge_weight_csv: str, output_csv: str, module_size: int = 3):
    from mrbigr.core import net

    result = net.module_identify(edge_weight_csv, module_size=int(module_size))
    if result is None:
        return None
    return _write_frame(result, output_csv)


def net_hub_identify_runner(
    edge_weight_csv: str,
    cluster_one_csv: str,
    output_cluster_csv: str,
    output_hub_csv: str,
):
    from mrbigr.core import net
    import pandas as pd

    edge_df = pd.read_csv(edge_weight_csv)
    cluster_df = pd.read_csv(cluster_one_csv)
    cluster_res, hub_res = net.hub_identify(edge_df, cluster_df)
    _ensure_parent(output_cluster_csv)
    if cluster_res is not None:
        cluster_res.to_csv(output_cluster_csv, index=False)
    if hub_res is not None and len(hub_res):
        hub_res.to_csv(output_hub_csv, index=False)
    return output_cluster_csv


def qtl_annotation_runner(qtl_file: str, annotation_file: str, output_csv: str):
    from mrbigr.core import anno

    result = anno.qtl_annotation(qtl_file, annotation_file)
    if result is None:
        return None
    return _write_frame(result, output_csv)


def qtl_extract_runner(geno_prefix: str, qtl_regions_csv: str, output_prefix: str, result_json: str):
    from mrbigr.core import qtl
    import pandas as pd

    qtl_regions = pd.read_csv(qtl_regions_csv)
    prefixes = qtl.extract_qtl_genotypes(geno_prefix, qtl_regions, output_prefix)
    files = []
    for prefix in prefixes:
        for suffix in (".bed", ".bim", ".fam", ".log", ".nosex"):
            path = Path(f"{prefix}{suffix}")
            if path.exists():
                files.append(str(path))
    return _write_json(result_json, {
        "output_prefix": output_prefix,
        "prefixes": prefixes,
        "files": files,
        "n_regions": int(len(qtl_regions)),
    })


def qtl_haplo_runner(geno_prefix: str, qtl_snp_list, output_prefix: str, output_csv: str):
    from mrbigr.core import qtl

    result = qtl.calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix)
    if result is None:
        return None
    return _write_frame(result, output_csv)


def parse_gtf_runner(gtf_file: str, output_file: str):
    from mrbigr.core import anno

    result = anno.parse_gtf(gtf_file)
    if result is None:
        return None
    return _write_frame(result, output_file)


def peak_boxplot_runner(
    pheno_file: str,
    geno_prefix: str,
    qtl_csv: str,
    output_dir: str,
    test_method: str = "t-test",
    result_json: Optional[str] = None,
):
    from mrbigr.core import peak

    outputs = peak.plot_qtl_boxplot(
        pheno_file,
        geno_prefix,
        qtl_csv,
        output_dir=output_dir,
        test_method=test_method,
    )
    out_dir = Path(output_dir).resolve()
    # The core writes qtl_boxplot_diagnostics.json with a human-readable reason
    # when it produces no plots. Surface it so an empty result is a *completed*
    # job with a reason, never an undiagnosable rc=1/empty-log failure (the
    # worker codes a None return as returncode=1/failed).
    reason = None
    diag_path = out_dir / "qtl_boxplot_diagnostics.json"
    if diag_path.is_file():
        try:
            reason = json.loads(diag_path.read_text(encoding="utf-8")).get("reason")
        except Exception:
            reason = None
    payload = {
        "output_dir": str(out_dir),
        "output_files": list(outputs or []),
        "n_plots": len(outputs or []),
        "summary_file": str(out_dir / "qtl_boxplot_summary.csv"),
        "reason": reason,
    }
    if result_json:
        return _write_json(result_json, payload)
    return payload


def go_enrich_runner(
    gene_list,
    output_file: str,
    organism: str = "Mouse",
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.05,
    mode: str = "local",
    gene_sets_file: Optional[str] = None,
    background_genes=None,
    go_obo_file: Optional[str] = None,
    obo_cache_dir: Optional[str] = None,
    auto_download_obo: bool = True,
    drop_unmapped_terms: bool = True,
):
    from mrbigr.core import go

    result = go.go_enrich(
        _coerce_gene_list(gene_list),
        organism=organism,
        pvalue_cutoff=float(pvalue_cutoff),
        qvalue_cutoff=float(qvalue_cutoff),
        mode=mode,
        gene_sets_file=gene_sets_file,
        background_genes=background_genes,
        go_obo_file=go_obo_file,
        obo_cache_dir=obo_cache_dir,
        auto_download_obo=_coerce_bool(auto_download_obo),
        drop_unmapped_terms=_coerce_bool(drop_unmapped_terms),
    )
    if result is None:
        return None
    return _write_frame(result, output_file)


def gsea_enrich_runner(
    gene_rank,
    output_file: str,
    gene_symbol_col: str = "gene",
    score_col: str = "score",
    organism: str = "Mouse",
    mode: str = "local",
    gene_sets_file: Optional[str] = None,
):
    from mrbigr.core import go

    rank_df = _coerce_frame(gene_rank)
    result = go.gsea_enrich(
        None,
        gene_rank=rank_df,
        gene_symbol_col=gene_symbol_col,
        score_col=score_col,
        organism=organism,
        mode=mode,
        gene_sets_file=gene_sets_file,
    )
    if result is None:
        return None
    return _write_frame(result, output_file)


def kegg_enrich_runner(gene_list, output_file: str, organism: str = "mmu", pvalue_cutoff: float = 0.05):
    from mrbigr.core import go

    result = go.kegg_enrich(
        _coerce_gene_list(gene_list),
        organism=organism,
        pvalue_cutoff=float(pvalue_cutoff),
    )
    if result is None:
        return None
    return _write_frame(result, output_file)


def core_call_runner(module: str, function: str, kwargs: dict):
    """Generic adapter: import ``module`` and call ``module.function(**kwargs)``.
    Used for plain ``mrbigr.core.geno`` calls (impute_genotype, snp_pruning,
    snp_clumping, ibd, vcf_to_plink, hapmap_to_plink, plink_to_vcf) where
    inputs and outputs are file paths and no DataFrame plumbing is needed.
    """
    import importlib

    mod = importlib.import_module(module)
    fn = getattr(mod, function)
    return fn(**(kwargs or {}))


__all__ = [
    "gwas_lm_runner",
    "gwas_lmm_runner",
    "gwas_plink_runner",
    "gwas_clump_runner",
    "mr_analysis_runner",
    "qtl_target_runner",
    "net_module_identify_runner",
    "net_hub_identify_runner",
    "qtl_annotation_runner",
    "qtl_extract_runner",
    "qtl_haplo_runner",
    "parse_gtf_runner",
    "peak_boxplot_runner",
    "go_enrich_runner",
    "gsea_enrich_runner",
    "kegg_enrich_runner",
    "core_call_runner",
]
