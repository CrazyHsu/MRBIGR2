#!/usr/bin/env python3
"""
GOMCP - GO Enrichment Analysis Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/go.py

Functions:
- GO enrichment analysis (BP, MF, CC)
- GSEA enrichment analysis
- GO visualization (barplot, dotplot, network)
- Gene set enrichment
"""
import pandas as pd
import numpy as np
import os
import ast
import gzip
import re
import shutil
import subprocess
import tempfile
import urllib.request
import warnings
from typing import List, Dict, Optional, Union
from scipy.stats import hypergeom

try:
    import gseapy as gp
    HAS_GSEAPY = True
except ImportError:
    HAS_GSEAPY = False
    warnings.warn("gseapy not installed. GO enrichment will be limited.")

# External tools
from .paths import repo_root as _repo_root
SCRIPT_DIR = str(_repo_root())
PLINK_BIN = os.path.join(SCRIPT_DIR, "utils", "plink")
GO_BASIC_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
GO_NAMESPACE_MAP = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}
_GO_OBO_CACHE = {}


def _ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _register_go_term(metadata: Dict[str, Dict[str, str]], current: Dict[str, Union[str, List[str]]]) -> None:
    if not current or not current.get("id") or current.get("is_obsolete"):
        return

    namespace = current.get("namespace")
    payload = {
        "name": current.get("name", current["id"]),
        "namespace": namespace,
        "ontology": GO_NAMESPACE_MAP.get(namespace),
    }
    if payload["ontology"] is None:
        return

    metadata[current["id"]] = payload
    for alt_id in current.get("alt_ids", []):
        metadata[alt_id] = payload


def _safe_filename(value: str) -> str:
    text = str(value).strip()
    text = re.sub(r"[^\w.-]+", "_", text)
    return text[:80] or "plot"


def _coerce_sequence(value):
    if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        return list(value)
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return []
        try:
            parsed = ast.literal_eval(text)
        except (SyntaxError, ValueError):
            return []
        if isinstance(parsed, (list, tuple, np.ndarray, pd.Series)):
            return list(parsed)
    return []


def _parse_overlap_ratio(value: str) -> Optional[float]:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if "/" not in text:
        return None
    left, right = text.split("/", 1)
    try:
        numerator = float(left)
        denominator = float(right)
    except ValueError:
        return None
    if denominator == 0:
        return None
    return numerator / denominator


def _normalize_enrichr_organism(organism: str) -> str:
    aliases = {
        'mmu': 'Mouse',
        'mouse': 'Mouse',
        'mus musculus': 'Mouse',
        'hsa': 'Human',
        'human': 'Human',
        'hsapiens': 'Human',
        'rno': 'Rat',
        'rat': 'Rat',
        'dme': 'Fly',
        'fly': 'Fly',
        'dre': 'Zebrafish',
        'zebrafish': 'Zebrafish',
        'sce': 'Yeast',
        'yeast': 'Yeast',
        'cel': 'Celegan',
        'celegans': 'Celegan',
    }
    if organism is None:
        return "Mouse"
    return aliases.get(str(organism).strip().lower(), organism)


def _coerce_gene_list(gene_list: Union[List[str], pd.DataFrame, pd.Series, np.ndarray, set, tuple, None], column_name: str = "genes") -> List[str]:
    if gene_list is None:
        return []
    if isinstance(gene_list, pd.DataFrame):
        if column_name in gene_list.columns:
            gene_list = gene_list[column_name].tolist()
        elif len(gene_list.columns) > 0:
            gene_list = gene_list.iloc[:, 0].tolist()
        else:
            gene_list = []
    elif isinstance(gene_list, pd.Series):
        gene_list = gene_list.tolist()
    elif isinstance(gene_list, np.ndarray):
        gene_list = gene_list.tolist()
    elif isinstance(gene_list, set):
        gene_list = list(gene_list)
    elif isinstance(gene_list, tuple):
        gene_list = list(gene_list)

    return sorted({str(g).strip() for g in gene_list if pd.notna(g) and str(g).strip()})


def _bh_adjust(pvalues: List[float]) -> List[float]:
    if not pvalues:
        return []

    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = np.empty_like(ranked)
    n = len(ranked)
    prev = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = min(prev, ranked[i] * n / rank)
        adjusted[i] = value
        prev = value

    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result.tolist()


def _download_go_obo(target_file: str, timeout: int = 120) -> str:
    _ensure_parent_dir(target_file)
    req = urllib.request.Request(GO_BASIC_OBO_URL)
    req.add_header("User-Agent", "Python-MRBIGR/1.0")
    with urllib.request.urlopen(req, timeout=timeout) as response, open(target_file, "wb") as handle:
        shutil.copyfileobj(response, handle)
    return target_file


def _resolve_go_obo_file(
    go_obo_file: Optional[str] = None,
    gene_sets_file: Optional[str] = None,
    obo_cache_dir: Optional[str] = None,
    auto_download_obo: bool = True,
) -> Optional[str]:
    candidates = []

    if go_obo_file:
        candidates.append(go_obo_file)
    if gene_sets_file:
        candidates.append(os.path.join(os.path.dirname(os.path.abspath(gene_sets_file)), "go-basic.obo"))
    if obo_cache_dir:
        candidates.append(os.path.join(obo_cache_dir, "go-basic.obo"))

    candidates.extend([
        os.path.join(os.getcwd(), "go-basic.obo"),
        os.path.join(SCRIPT_DIR, "data", "go-basic.obo"),
        os.path.join(tempfile.gettempdir(), "mrbigr_go_cache", "go-basic.obo"),
    ])

    seen = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        if os.path.exists(candidate):
            return candidate

    if not auto_download_obo:
        return None

    cache_dir = obo_cache_dir or os.path.join(tempfile.gettempdir(), "mrbigr_go_cache")
    target_file = os.path.join(cache_dir, "go-basic.obo")
    if os.path.exists(target_file):
        return target_file

    try:
        return _download_go_obo(target_file)
    except Exception as e:
        warnings.warn(f"Failed to download go-basic.obo: {e}")
        return None


def _load_go_obo_metadata(
    go_obo_file: Optional[str] = None,
    gene_sets_file: Optional[str] = None,
    obo_cache_dir: Optional[str] = None,
    auto_download_obo: bool = True,
) -> Dict[str, Dict[str, str]]:
    obo_path = _resolve_go_obo_file(
        go_obo_file=go_obo_file,
        gene_sets_file=gene_sets_file,
        obo_cache_dir=obo_cache_dir,
        auto_download_obo=auto_download_obo,
    )
    if obo_path is None:
        return {}
    if obo_path in _GO_OBO_CACHE:
        return _GO_OBO_CACHE[obo_path]

    metadata = {}
    current = None
    open_fn = gzip.open if str(obo_path).lower().endswith(".gz") else open

    with open_fn(obo_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if line == "[Term]":
                _register_go_term(metadata, current)
                current = {}
                continue

            if current is None:
                continue

            if not line:
                _register_go_term(metadata, current)
                current = None
                continue

            if line.startswith("id: GO:"):
                current["id"] = line.split("id: ", 1)[1]
            elif line.startswith("alt_id: GO:"):
                current.setdefault("alt_ids", []).append(line.split("alt_id: ", 1)[1])
            elif line.startswith("name: "):
                current["name"] = line.split("name: ", 1)[1]
            elif line.startswith("namespace: "):
                current["namespace"] = line.split("namespace: ", 1)[1]
            elif line == "is_obsolete: true":
                current["is_obsolete"] = True

    _register_go_term(metadata, current)

    _GO_OBO_CACHE[obo_path] = metadata
    return metadata


def _load_local_gene_sets(
    gene_sets_file: str,
    min_set_size: int = 3,
    max_set_size: int = 500,
) -> Dict[str, set]:
    gene_sets = {}
    open_fn = gzip.open if str(gene_sets_file).lower().endswith(".gz") else open

    with open_fn(gene_sets_file, "rt", encoding="utf-8", errors="replace") as handle:
        first_line = ""
        for line in handle:
            if line.strip():
                first_line = line.rstrip("\n")
                break

    if not first_line:
        return gene_sets

    first_parts = first_line.split("\t")
    is_gmt = len(first_parts) >= 3

    if is_gmt:
        with open_fn(gene_sets_file, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                term = parts[0].strip()
                genes = {g.strip() for g in parts[2:] if g.strip()}
                if min_set_size <= len(genes) <= max_set_size:
                    gene_sets[term] = genes
        return gene_sets

    mapping = pd.read_csv(
        gene_sets_file,
        sep=r"\s+|\t+",
        engine="python",
        header=None,
        names=["gene", "term"],
        usecols=[0, 1],
    ).dropna()
    mapping["gene"] = mapping["gene"].astype(str).str.strip()
    mapping["term"] = mapping["term"].astype(str).str.strip()

    for term, term_df in mapping.groupby("term"):
        genes = set(term_df["gene"])
        if min_set_size <= len(genes) <= max_set_size:
            gene_sets[term] = genes

    return gene_sets


def _normalize_gsea_results(results) -> pd.DataFrame:
    if isinstance(results, pd.DataFrame):
        return results
    if isinstance(results, dict):
        records = []
        for term, payload in results.items():
            if not isinstance(payload, dict):
                continue
            records.append({
                "Term": term,
                "ES": payload.get("es"),
                "NES": payload.get("nes"),
                "NOM p-val": payload.get("pval"),
                "FDR q-val": payload.get("fdr"),
                "FWER p-val": payload.get("fwerp"),
                "Tag %": payload.get("tag %"),
                "Gene %": payload.get("gene %"),
                "Lead_genes": payload.get("lead_genes"),
                "matched_genes": payload.get("matched_genes"),
                "hits": payload.get("hits"),
                "RES": payload.get("RES"),
            })
        return pd.DataFrame(records)
    return pd.DataFrame()


def _collect_prerank_results(pre_res) -> pd.DataFrame:
    frame = pd.DataFrame()
    if pre_res is not None and hasattr(pre_res, "res2d") and pre_res.res2d is not None:
        frame = pd.DataFrame(pre_res.res2d).copy()

    details = _normalize_gsea_results(getattr(pre_res, "results", None))
    if not details.empty:
        if frame.empty:
            frame = details.copy()
        else:
            detail_cols = [c for c in details.columns if c not in frame.columns or c == "Term"]
            frame = frame.merge(details[detail_cols], on="Term", how="left")

    if frame.empty:
        return frame

    ranking = getattr(pre_res, "ranking", None)
    if ranking is not None:
        try:
            if isinstance(ranking, pd.Series):
                frame.attrs["rank_metric"] = ranking.astype(float).tolist()
            else:
                frame.attrs["rank_metric"] = list(ranking)
        except TypeError:
            pass

    return frame


def _prepare_go_plot_frame(enrich_result: pd.DataFrame) -> pd.DataFrame:
    if enrich_result is None or len(enrich_result) == 0:
        return pd.DataFrame()

    df = enrich_result.copy()
    if "ontology" not in df.columns and "ONTOLOGY" in df.columns:
        df["ontology"] = df["ONTOLOGY"]
    if "Term" not in df.columns and "Description" in df.columns:
        df["Term"] = df["Description"]
    if "Term" not in df.columns and "GO_ID" in df.columns:
        df["Term"] = df["GO_ID"]

    if "Count" not in df.columns:
        if "Overlap" in df.columns:
            df["Count"] = df["Overlap"].astype(str).str.split("/").str[0]
        elif "Genes" in df.columns:
            df["Count"] = df["Genes"].astype(str).str.split(r"[;,]").map(
                lambda genes: len([g for g in genes if str(g).strip()])
            )
    if "Count" in df.columns:
        df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
    if "Query_size" in df.columns:
        df["Query_size"] = pd.to_numeric(df["Query_size"], errors="coerce")
    if "Gene_set_size" in df.columns:
        df["Gene_set_size"] = pd.to_numeric(df["Gene_set_size"], errors="coerce")

    if "GeneRatio" not in df.columns:
        if "Query_size" in df.columns and df["Query_size"].notna().any():
            denominator = df["Query_size"].replace(0, np.nan)
            df["GeneRatio"] = df["Count"] / denominator
        elif "Overlap" in df.columns:
            df["GeneRatio"] = df["Overlap"].map(_parse_overlap_ratio)

    if "BgRatio" not in df.columns and "Background_size" in df.columns and "Gene_set_size" in df.columns:
        bg = pd.to_numeric(df["Background_size"], errors="coerce").replace(0, np.nan)
        df["BgRatio"] = df["Gene_set_size"] / bg

    return df

# ========== GO Enrichment Functions ==========

def go_enrich(
    gene_list: Union[List[str], pd.DataFrame],
    organism: str = "Mouse",
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.05,
    min_set_size: int = 3,
    max_set_size: int = 500,
    ontologies: List[str] = ["BP", "MF", "CC"],
    mode: str = "local",
    gene_sets_file: Optional[str] = None,
    background_genes: Optional[Union[List[str], pd.DataFrame, pd.Series]] = None,
    go_obo_file: Optional[str] = None,
    obo_cache_dir: Optional[str] = None,
    auto_download_obo: bool = True,
    drop_unmapped_terms: bool = True,
) -> pd.DataFrame:
    """GO enrichment analysis using gseapy.
    
    Args:
        gene_list: List of gene IDs or DataFrame with 'genes' column
        organism: Organism name ('Human', 'Mouse', 'Yeast', etc.)
        pvalue_cutoff: P-value cutoff
        qvalue_cutoff: Q-value cutoff (FDR)
        min_set_size: Minimum gene set size
        max_set_size: Maximum gene set size
        ontologies: List of ontologies to test ('BP', 'MF', 'CC')
    
    Returns:
        DataFrame with GO enrichment results
    """
    mode = str(mode).strip().lower()
    if mode not in {"auto", "online", "local"}:
        raise ValueError("mode must be one of: auto, online, local")

    gene_list = _coerce_gene_list(gene_list)
    if len(gene_list) == 0:
        return pd.DataFrame()

    if mode in {"auto", "local"} and gene_sets_file:
        return local_go_enrich(
            gene_list=gene_list,
            gene_sets_file=gene_sets_file,
            pvalue_cutoff=pvalue_cutoff,
            qvalue_cutoff=qvalue_cutoff,
            min_set_size=min_set_size,
            max_set_size=max_set_size,
            ontologies=ontologies,
            background_genes=background_genes,
            go_obo_file=go_obo_file,
            obo_cache_dir=obo_cache_dir,
            auto_download_obo=auto_download_obo,
            drop_unmapped_terms=drop_unmapped_terms,
        )

    if mode == "local" and not gene_sets_file:
        warnings.warn("Local GO enrichment requires gene_sets_file.")
        return pd.DataFrame()

    if not HAS_GSEAPY:
        return pd.DataFrame()
    
    results = []
    
    for ont in ontologies:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets='GO_Biological_Process_2023' if ont == 'BP' else ('GO_Molecular_Function_2023' if ont == 'MF' else 'GO_Cellular_Component_2023'),
                organism=_normalize_enrichr_organism(organism).lower(),
                outdir=None,
                no_plot=True
            )
            
            if enr is not None and hasattr(enr, 'results'):
                df = enr.results
                df['ontology'] = ont
                results.append(df)
        except Exception as e:
            warnings.warn(f"GO enrichment failed for {ont}: {e}")
            continue
    
    if len(results) == 0:
        return pd.DataFrame()
    
    # Combine results
    combined = pd.concat(results, ignore_index=True)
    
    # Filter by p-value and q-value
    if 'P-value' in combined.columns:
        combined = combined[combined['P-value'] <= pvalue_cutoff]
    if 'Adjusted P-value' in combined.columns:
        combined = combined[combined['Adjusted P-value'] <= qvalue_cutoff]
    
    return combined


def local_go_enrich(
    gene_list: Union[List[str], pd.DataFrame, pd.Series],
    gene_sets_file: str,
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.05,
    min_set_size: int = 3,
    max_set_size: int = 500,
    ontologies: Optional[List[str]] = None,
    background_genes: Optional[Union[List[str], pd.DataFrame, pd.Series]] = None,
    go_obo_file: Optional[str] = None,
    obo_cache_dir: Optional[str] = None,
    auto_download_obo: bool = True,
    drop_unmapped_terms: bool = True,
) -> pd.DataFrame:
    """Offline GO enrichment from a local mapping or GMT file."""
    query_genes = _coerce_gene_list(gene_list)
    if not query_genes:
        return pd.DataFrame()

    gene_sets = _load_local_gene_sets(
        gene_sets_file=gene_sets_file,
        min_set_size=min_set_size,
        max_set_size=max_set_size,
    )
    if not gene_sets:
        return pd.DataFrame()

    if background_genes is not None:
        background = set(_coerce_gene_list(background_genes))
    else:
        background = set()
        for genes in gene_sets.values():
            background.update(genes)

    query = set(query_genes) & background
    if not query:
        return pd.DataFrame()

    go_metadata = _load_go_obo_metadata(
        go_obo_file=go_obo_file,
        gene_sets_file=gene_sets_file,
        obo_cache_dir=obo_cache_dir,
        auto_download_obo=auto_download_obo,
    )
    if drop_unmapped_terms and not go_metadata:
        warnings.warn("No GO metadata available; local GO enrichment cannot resolve GO terms.")
        return pd.DataFrame()
    if go_metadata:
        ontology_filter = {str(o).strip().upper() for o in (ontologies or []) if str(o).strip()}
    else:
        ontology_filter = set()

    M = len(background)
    N = len(query)
    records = []

    for term_id, genes in gene_sets.items():
        term_background = genes & background
        n = len(term_background)
        if n < min_set_size or n > max_set_size:
            continue

        meta = go_metadata.get(term_id, {})
        if drop_unmapped_terms and not meta:
            continue
        ontology = meta.get("ontology", None if drop_unmapped_terms else "LOCAL")
        if ontology_filter and ontology not in ontology_filter:
            continue

        overlap_genes = sorted(query & term_background)
        k = len(overlap_genes)
        if k == 0:
            continue

        pvalue = float(hypergeom.sf(k - 1, M, n, N))
        records.append({
            "GO_ID": term_id,
            "Term": meta.get("name", term_id if not drop_unmapped_terms else None),
            "namespace": meta.get("namespace"),
            "ontology": ontology,
            "Overlap": f"{k}/{n}",
            "P-value": pvalue,
            "Genes": ",".join(overlap_genes),
            "Count": k,
            "Gene_set_size": n,
            "Query_size": N,
            "Background_size": M,
        })

    if not records:
        return pd.DataFrame()

    result = pd.DataFrame(records).sort_values("P-value", ascending=True).reset_index(drop=True)
    result["Adjusted P-value"] = _bh_adjust(result["P-value"].tolist())
    result = result[(result["P-value"] <= pvalue_cutoff) & (result["Adjusted P-value"] <= qvalue_cutoff)]
    return result.reset_index(drop=True)


def gsea_enrich(
    gene_list: Union[List[str], pd.DataFrame],
    gene_rank: Optional[pd.DataFrame] = None,
    gene_symbol_col: str = "gene",
    score_col: str = "score",
    organism: str = "Mouse",
    pvalue_cutoff: float = 0.05,
    min_set_size: int = 10,
    max_set_size: int = 500,
    processes: int = 1,
    mode: str = "local",
    gene_sets_file: Optional[str] = None,
) -> pd.DataFrame:
    """Gene Set Enrichment Analysis (GSEA).
    
    Args:
        gene_list: List of genes or DataFrame with ranked genes
        gene_rank: DataFrame with gene ranks (if None, use gene_list)
        gene_symbol_col: Column name for gene symbols
        score_col: Column name for ranking scores
        organism: Organism name
        pvalue_cutoff: P-value cutoff
        min_set_size: Minimum gene set size
        max_set_size: Maximum gene set size
        processes: Number of processes
    
    Returns:
        DataFrame with GSEA results
    """
    mode = str(mode).strip().lower()
    if mode not in {"auto", "online", "local"}:
        raise ValueError("mode must be one of: auto, online, local")

    if not HAS_GSEAPY:
        return pd.DataFrame()
    
    # Handle input format
    if isinstance(gene_list, pd.DataFrame):
        if gene_rank is None:
            gene_rank = gene_list
    
    if gene_rank is not None:
        if gene_symbol_col not in gene_rank.columns or score_col not in gene_rank.columns:
            return pd.DataFrame()
        gene_rank = gene_rank.sort_values(score_col, ascending=False)
        gene_list = gene_rank[gene_symbol_col].tolist()
    
    # Remove duplicates and NaN
    gene_list = [str(g).strip() for g in gene_list if pd.notna(g)]
    gene_list = list(set(gene_list))
    
    try:
        rank_df = gene_rank[[gene_symbol_col, score_col]].copy()
        if mode in {"auto", "local"} and gene_sets_file:
            gene_sets_source = gene_sets_file
        elif mode == "local":
            warnings.warn("Local GSEA requires gene_sets_file.")
            return pd.DataFrame()
        else:
            gene_sets_source = 'MSigDB_Hallmark_2020'

        pre_res = gp.prerank(
            rnk=rank_df,
            gene_sets=gene_sets_source,
            outdir=None,
            no_plot=True,
            min_size=min_set_size,
            max_size=max_set_size,
            threads=processes
        )

        frame = _collect_prerank_results(pre_res)
        if not frame.empty:
            return frame
    except Exception as e:
        warnings.warn(f"GSEA failed: {e}")
    
    return pd.DataFrame()


def go_plot(
    enrich_result: pd.DataFrame,
    output_prefix: str,
    plot_types: List[str] = ["barplot", "dotplot"],
    top_n: int = 20,
    figsize: tuple = (8, 6),
    split_ontology: bool = True,
    label_col: str = "Term",
    bar_x: str = "GeneRatio",
    bar_color: str = "Adjusted P-value",
    dot_x: str = "GeneRatio",
    dot_color: str = "Adjusted P-value",
    dot_size: str = "Count",
    sort_by: Optional[str] = None,
    ascending: Optional[bool] = None,
) -> Dict[str, str]:
    """Generate GO enrichment plots.
    
    Args:
        enrich_result: GO enrichment result DataFrame
        output_prefix: Output file prefix
        plot_types: Types of plots ('barplot', 'dotplot', 'network')
        top_n: Number of top terms to plot
        figsize: Figure size
    
    Returns:
        Dictionary with output file paths
    """
    if enrich_result is None or len(enrich_result) == 0:
        return {}

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    def _pick_sort_column(df: pd.DataFrame, preferred: Optional[str], fallback_x: str) -> str:
        if preferred and preferred in df.columns:
            return preferred
        for candidate in ["Adjusted P-value", "P-value", "GeneRatio", "Count", fallback_x]:
            if candidate in df.columns:
                return candidate
        return df.columns[0]

    def _default_ascending(metric: str) -> bool:
        return metric in {"Adjusted P-value", "P-value", "FDR q-val", "NOM p-val"}

    def _resolve_color_mapping(df: pd.DataFrame, column: Optional[str], fallback: str, categorical_fallback: Optional[str] = None):
        if column and column in df.columns:
            raw = df[column]
        elif categorical_fallback and categorical_fallback in df.columns:
            column = categorical_fallback
            raw = df[column]
        else:
            column = fallback
            raw = df.get(column, pd.Series([fallback] * len(df), index=df.index))

        numeric = pd.to_numeric(raw, errors="coerce")
        if numeric.notna().all():
            plot_vals = numeric.astype(float)
            label = column
            cmap = plt.cm.viridis
            if column in {"Adjusted P-value", "P-value", "FDR q-val", "NOM p-val"}:
                plot_vals = -np.log10(np.clip(plot_vals, 1e-300, None))
                label = f"-log10({column})"
                cmap = plt.cm.viridis_r
            vmin = float(plot_vals.min()) if len(plot_vals) else 0.0
            vmax = float(plot_vals.max()) if len(plot_vals) else 1.0
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0
            return {
                "kind": "numeric",
                "column": column,
                "values": plot_vals,
                "label": label,
                "norm": Normalize(vmin=vmin, vmax=vmax),
                "cmap": cmap,
            }

        categories = raw.fillna("Unknown").astype(str)
        unique = list(dict.fromkeys(categories.tolist()))
        palette = plt.cm.tab10(np.linspace(0, 1, max(len(unique), 1)))
        color_map = {cat: palette[i] for i, cat in enumerate(unique)}
        return {
            "kind": "categorical",
            "column": column,
            "values": categories,
            "label": column,
            "color_map": color_map,
        }

    def _finalize_subset(df: pd.DataFrame, group_name: str, metric: str, is_ascending: bool) -> pd.DataFrame:
        subset = df.copy()
        subset = subset.dropna(subset=[label_col])
        if metric in subset.columns:
            subset = subset.sort_values(metric, ascending=is_ascending).head(top_n)
        else:
            subset = subset.head(top_n)
        return subset.copy()

    def _plot_barplot(df: pd.DataFrame, output_file: str, title: str) -> None:
        subset = df.dropna(subset=[bar_x])
        if subset.empty:
            return
        subset = subset.sort_values(bar_x, ascending=True)
        labels = subset[label_col].astype(str).tolist()
        color_info = _resolve_color_mapping(subset, bar_color, fallback="steelblue", categorical_fallback="ontology")
        fig_height = max(figsize[1], len(subset) * 0.38 + 1.8)
        fig, ax = plt.subplots(figsize=(figsize[0], fig_height))
        xvals = pd.to_numeric(subset[bar_x], errors="coerce").fillna(0.0)

        if color_info["kind"] == "numeric":
            colors = color_info["cmap"](color_info["norm"](color_info["values"]))
        else:
            colors = color_info["values"].map(color_info["color_map"]).tolist()

        ax.barh(labels, xvals, color=colors, edgecolor="grey", linewidth=0.5)
        ax.set_xlabel(bar_x)
        ax.set_ylabel("")
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        if color_info["kind"] == "numeric":
            sm = ScalarMappable(norm=color_info["norm"], cmap=color_info["cmap"])
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.7, pad=0.02)
            cbar.set_label(color_info["label"])
        else:
            handles = [
                Line2D([0], [0], marker="s", color="w", label=cat, markerfacecolor=color, markersize=8)
                for cat, color in color_info["color_map"].items()
            ]
            ax.legend(handles=handles, title=color_info["label"], bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

        plt.tight_layout()
        _ensure_parent_dir(output_file)
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)

    def _plot_dotplot(df: pd.DataFrame, output_file: str, title: str) -> None:
        subset = df.dropna(subset=[dot_x, dot_size])
        if subset.empty:
            return
        subset = subset.sort_values(dot_x, ascending=True)
        labels = subset[label_col].astype(str).tolist()
        y_positions = np.arange(len(subset))
        xvals = pd.to_numeric(subset[dot_x], errors="coerce").fillna(0.0)
        size_raw = pd.to_numeric(subset[dot_size], errors="coerce").fillna(1.0)
        if np.isclose(size_raw.min(), size_raw.max()):
            sizes = np.repeat(220.0, len(size_raw))
        else:
            sizes = np.interp(size_raw, (size_raw.min(), size_raw.max()), (80.0, 500.0))

        color_info = _resolve_color_mapping(subset, dot_color, fallback="steelblue", categorical_fallback="ontology")
        fig_height = max(figsize[1], len(subset) * 0.38 + 1.8)
        fig, ax = plt.subplots(figsize=(figsize[0], fig_height))

        if color_info["kind"] == "numeric":
            scatter = ax.scatter(
                xvals,
                y_positions,
                s=sizes,
                c=color_info["values"],
                cmap=color_info["cmap"],
                norm=color_info["norm"],
                alpha=0.9,
                edgecolors="black",
                linewidths=0.4,
            )
            cbar = plt.colorbar(scatter, ax=ax, shrink=0.7, pad=0.02)
            cbar.set_label(color_info["label"])
        else:
            face_colors = color_info["values"].map(color_info["color_map"]).tolist()
            ax.scatter(
                xvals,
                y_positions,
                s=sizes,
                c=face_colors,
                alpha=0.9,
                edgecolors="black",
                linewidths=0.4,
            )
            handles = [
                Line2D([0], [0], marker="o", color="w", label=cat, markerfacecolor=color, markersize=8)
                for cat, color in color_info["color_map"].items()
            ]
            ax.legend(handles=handles, title=color_info["label"], bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

        legend_values = np.unique(np.round(np.linspace(size_raw.min(), size_raw.max(), num=min(3, len(size_raw))), 4))
        legend_handles = []
        legend_labels = []
        for value in legend_values:
            if np.isclose(size_raw.min(), size_raw.max()):
                marker_size = np.sqrt(220.0)
            else:
                marker_size = np.sqrt(np.interp(value, (size_raw.min(), size_raw.max()), (80.0, 500.0)))
            legend_handles.append(
                Line2D([0], [0], marker="o", color="grey", linestyle="", markersize=marker_size / 2)
            )
            legend_labels.append(f"{value:g}")
        if legend_handles:
            size_legend = ax.legend(
                legend_handles,
                legend_labels,
                title=dot_size,
                bbox_to_anchor=(1.02, 0.55),
                loc="upper left",
                frameon=False,
            )
            ax.add_artist(size_legend)

        ax.set_xlabel(dot_x)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels)
        ax.set_title(title)
        ax.grid(axis="x", linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        plt.tight_layout()
        _ensure_parent_dir(output_file)
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)

    df = _prepare_go_plot_frame(enrich_result)
    if df.empty:
        return {}

    plot_types = [plot_types] if isinstance(plot_types, str) else list(plot_types)
    if label_col not in df.columns:
        raise ValueError(f"label_col '{label_col}' not found in enrichment result")
    required_columns = set()
    if "barplot" in plot_types:
        required_columns.add(bar_x)
    if "dotplot" in plot_types:
        required_columns.update({dot_x, dot_size})
    for required_col in required_columns:
        if required_col not in df.columns:
            raise ValueError(f"Required plotting column '{required_col}' not found in enrichment result")

    output_files = {}
    ont_col = "ontology" if "ontology" in df.columns else None
    group_items = []
    if split_ontology and ont_col:
        preferred_order = [ont for ont in ["BP", "CC", "MF"] if ont in set(df[ont_col].dropna().astype(str))]
        extras = [ont for ont in df[ont_col].dropna().astype(str).unique().tolist() if ont not in preferred_order]
        for ont in preferred_order + extras:
            group_items.append((ont, df[df[ont_col].astype(str) == ont]))
    else:
        group_items.append(("combined" if ont_col else "all", df))

    for group_name, group_df in group_items:
        metric = _pick_sort_column(group_df, sort_by, bar_x)
        subset = _finalize_subset(group_df, group_name, metric, _default_ascending(metric) if ascending is None else ascending)
        if subset.empty:
            continue

        for ptype in plot_types:
            if group_name in {"all", "combined"}:
                suffix = ptype
                result_key = f"{group_name}_{ptype}" if group_name == "combined" else ptype
            else:
                suffix = f"{group_name}_{ptype}"
                result_key = suffix
            output_file = f"{output_prefix}_{suffix}.png"
            title = f"GO Enrichment ({group_name})" if group_name not in {"all", "combined"} else "GO Enrichment"
            try:
                if ptype == "barplot":
                    _plot_barplot(subset, output_file, title)
                elif ptype == "dotplot":
                    _plot_dotplot(subset, output_file, title)
                else:
                    warnings.warn(f"Unsupported plot type: {ptype}")
                    continue
                output_files[result_key] = output_file
            except Exception as e:
                warnings.warn(f"{ptype} failed for {group_name}: {e}")

    return output_files


def gsea_plot(
    gsea_result: pd.DataFrame,
    output_prefix: str,
    top_n: int = 20,
    figsize: tuple = (8, 6),
    format: str = "png",
    plot_mode: Union[str, List[str]] = "curve",
    terms: Optional[Union[str, List[str]]] = None,
    curve_terms: int = 1,
    trace_terms: int = 3,
    rank_metric: Optional[Union[pd.Series, List[float], np.ndarray]] = None,
) -> Dict[str, str]:
    """Generate GSEA plots.

    Args:
        gsea_result: GSEA result DataFrame (from gsea_enrich)
        output_prefix: Output file prefix
        top_n: Number of top pathways to plot
        figsize: Figure size
        format: Output format ('png', 'pdf')
        plot_mode: One or more plot types ('curve', 'trace', 'nes_barplot')
        terms: Optional explicit term or list of terms to plot
        curve_terms: Number of leading terms for curve plots when terms is None
        trace_terms: Number of leading terms for trace plots when terms is None

    Returns:
        Dictionary with output file paths
    """
    if gsea_result is None or len(gsea_result) == 0:
        return {}

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    output_files = {}
    plot_modes = [plot_mode] if isinstance(plot_mode, str) else list(plot_mode)
    df = gsea_result.copy()
    if "Term" not in df.columns and "Name" in df.columns:
        df["Term"] = df["Name"]

    nes_col = "NES" if "NES" in df.columns else None
    term_col = "Term" if "Term" in df.columns else None
    fdr_col = "FDR q-val" if "FDR q-val" in df.columns else ("Adjusted P-value" if "Adjusted P-value" in df.columns else None)
    pval_col = "NOM p-val" if "NOM p-val" in df.columns else ("P-value" if "P-value" in df.columns else None)
    if nes_col is None or term_col is None:
        warnings.warn("GSEA result missing NES or Term column for plotting.")
        return output_files

    if terms is None:
        selected_terms = None
    elif isinstance(terms, str):
        selected_terms = [terms]
    else:
        selected_terms = [str(term) for term in terms]

    rank_metric_values = rank_metric
    if rank_metric_values is None:
        rank_metric_values = df.attrs.get("rank_metric")
    if isinstance(rank_metric_values, pd.Series):
        rank_metric_values = rank_metric_values.tolist()
    elif isinstance(rank_metric_values, np.ndarray):
        rank_metric_values = rank_metric_values.tolist()

    sort_cols = []
    if fdr_col:
        sort_cols.append((fdr_col, True))
    if nes_col:
        df["__abs_nes"] = pd.to_numeric(df[nes_col], errors="coerce").abs()
        sort_cols.append(("__abs_nes", False))
    for col_name, asc in reversed(sort_cols):
        df = df.sort_values(col_name, ascending=asc, na_position="last")

    def _pick_rows(limit: int) -> pd.DataFrame:
        if selected_terms:
            return df[df[term_col].astype(str).isin(selected_terms)].copy()
        return df.head(limit).copy()

    if "curve" in plot_modes or "trace" in plot_modes:
        if "hits" not in df.columns or "RES" not in df.columns:
            warnings.warn("GSEA curve plotting requires 'hits' and 'RES' columns.")
        else:
            if "curve" in plot_modes:
                curve_df = _pick_rows(max(curve_terms, 1))
                for _, row in curve_df.iterrows():
                    term_name = str(row[term_col])
                    hits = [int(v) for v in _coerce_sequence(row.get("hits"))]
                    res = [float(v) for v in _coerce_sequence(row.get("RES"))]
                    if not hits or not res:
                        warnings.warn(f"Skipping GSEA curve plot for {term_name}: missing hits/RES data")
                        continue
                    output_file = f"{output_prefix}_{_safe_filename(term_name)}_curve.{format}"
                    _ensure_parent_dir(output_file)
                    gp.gseaplot(
                        term=term_name,
                        hits=hits,
                        nes=float(row.get(nes_col, np.nan)),
                        pval=float(row.get(pval_col, 1.0)) if pval_col else 1.0,
                        fdr=float(row.get(fdr_col, 1.0)) if fdr_col else 1.0,
                        RES=res,
                        rank_metric=rank_metric_values,
                        figsize=figsize,
                        ofname=output_file,
                    )
                    output_files[f"curve_{_safe_filename(term_name)}"] = output_file

            if "trace" in plot_modes:
                trace_df = _pick_rows(max(trace_terms, 1))
                trace_terms_list = []
                trace_hits = []
                trace_res = []
                for _, row in trace_df.iterrows():
                    hits = [int(v) for v in _coerce_sequence(row.get("hits"))]
                    res = [float(v) for v in _coerce_sequence(row.get("RES"))]
                    if not hits or not res:
                        continue
                    trace_terms_list.append(str(row[term_col]))
                    trace_hits.append(hits)
                    trace_res.append(res)
                if trace_terms_list:
                    output_file = f"{output_prefix}_trace.{format}"
                    _ensure_parent_dir(output_file)
                    gp.gseaplot2(
                        terms=trace_terms_list,
                        hits=trace_hits,
                        RESs=trace_res,
                        rank_metric=rank_metric_values,
                        figsize=figsize,
                        ofname=output_file,
                    )
                    output_files["trace"] = output_file

    if "nes_barplot" in plot_modes:
        plot_df = df.head(top_n).copy().sort_values(nes_col, ascending=True)
        terms_list = plot_df[term_col].astype(str).tolist()
        nes_vals = pd.to_numeric(plot_df[nes_col], errors="coerce").fillna(0.0).tolist()
        fdr_vals = pd.to_numeric(plot_df[fdr_col], errors="coerce").fillna(0.05).tolist() if fdr_col else [0.05] * len(nes_vals)

        labels = [term[:60] for term in terms_list]
        n = len(labels)
        bar_height = max(figsize[1], n * 0.35 + 1.5)
        fig, ax = plt.subplots(figsize=(figsize[0], bar_height))

        neg_log_fdr = [-np.log10(max(f, 1e-10)) for f in fdr_vals]
        norm = Normalize(vmin=min(neg_log_fdr), vmax=max(neg_log_fdr)) if neg_log_fdr else Normalize(0, 1)
        cmap = plt.cm.RdYlBu_r
        colors = [cmap(norm(v)) for v in neg_log_fdr]

        ax.barh(range(n), nes_vals, color=colors, edgecolor="grey", linewidth=0.5)
        ax.set_yticks(range(n))
        ax.set_yticklabels(labels, fontsize=max(6, 10 - n // 10))
        ax.set_xlabel("Normalized Enrichment Score (NES)")
        ax.set_title("GSEA Summary")
        ax.axvline(0, color="black", linewidth=0.8)

        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("-log10(FDR)")

        plt.tight_layout()
        output_file = f"{output_prefix}_nes_barplot.{format}"
        _ensure_parent_dir(output_file)
        fig.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close(fig)
        output_files["nes_barplot"] = output_file

    if "__abs_nes" in df.columns:
        df = df.drop(columns=["__abs_nes"])
    return output_files


def enrichr_library_list(organism: str = "Mouse") -> List[str]:
    """Get available gene set libraries for an organism.
    
    Args:
        organism: Organism name
    
    Returns:
        List of available gene set libraries
    """
    if not HAS_GSEAPY:
        return []
    
    try:
        libraries = gp.get_library_name(organism=organism)
        return libraries
    except Exception as e:
        warnings.warn(f"Failed to get library list: {e}")
        return []


def cheak_orgdb(orgdb_name: str) -> bool:
    """Check if organism database is available.
    
    Args:
        orgdb_name: Organism database name
    
    Returns:
        True if available
    """
    # For gseapy, check if organism is supported
    supported_organisms = ['Human', 'Mouse', 'Rat', 'Fly', 'Yeast', 'Zebrafish', 'Celegan']
    return orgdb_name in supported_organisms


def build_gene_set(
    gene_list: Union[List[str], pd.DataFrame],
    background_genes: Optional[List[str]] = None,
    min_size: int = 3,
    max_size: int = 500
) -> pd.DataFrame:
    """Build gene set from gene list with background.
    
    Args:
        gene_list: List of query genes
        background_genes: Background gene universe
        min_size: Minimum gene set size
        max_size: Maximum gene set size
    
    Returns:
        DataFrame with gene set info
    """
    if isinstance(gene_list, pd.DataFrame):
        if 'genes' in gene_list.columns:
            gene_list = gene_list['genes'].tolist()
    
    gene_list = [str(g).strip() for g in gene_list if pd.notna(g)]
    gene_set = list(set(gene_list))
    
    if background_genes is not None:
        background = set([str(g).strip() for g in background_genes if pd.notna(g)])
        gene_set = [g for g in gene_set if g in background]
    
    return pd.DataFrame({
        'gene': gene_set,
        'in_set': [True] * len(gene_set)
    })


def go_annotation_from_gtf(
    gtf_file: str,
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """Extract GO annotations from GTF file.
    
    Args:
        gtf_file: GTF annotation file
        output_file: Output file for annotations
    
    Returns:
        DataFrame with gene-GO mappings
    """
    go_mapping = []
    
    try:
        open_fn = gzip.open if str(gtf_file).lower().endswith('.gz') else open
        with open_fn(gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                # Check for GO annotation in attributes
                attrs = parts[8]
                gene_id = None
                gene_name = None
                go_terms = []
                
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                    elif 'GO:' in attr:
                        raw_terms = attr.split('"')[1] if '"' in attr else attr
                        go_terms.extend(re.findall(r'GO:\d{7}', raw_terms))
                
                gene_value = gene_id or gene_name
                for go in go_terms:
                    go_mapping.append({
                        'gene': gene_value,
                        'go_id': go
                    })
        
        result = pd.DataFrame(go_mapping)
        
        if output_file is not None:
            _ensure_parent_dir(output_file)
            result.to_csv(output_file, sep='\t', index=False)
        
        return result
    
    except Exception as e:
        warnings.warn(f"Failed to parse GTF: {e}")
        return pd.DataFrame()


def simplify_go_terms(
    enrich_result: pd.DataFrame,
    similarity_cutoff: float = 0.7
) -> pd.DataFrame:
    """Simplify GO terms by removing redundancy.
    
    Args:
        enrich_result: GO enrichment result
        similarity_cutoff: Similarity cutoff
    
    Returns:
        Simplified GO enrichment result
    """
    if len(enrich_result) == 0:
        return enrich_result
    
    # For now, return top terms per ontology
    ont_col = 'ontology' if 'ontology' in enrich_result.columns else 'ONTOLOGY'
    
    # Remove terms with high overlap
    simplified = []
    
    for ont in enrich_result[ont_col].unique():
        ont_data = enrich_result[enrich_result[ont_col] == ont].copy()
        pval_col = next((c for c in ['P-value', 'pvalue', 'Adjusted P-value', 'p_value'] if c in ont_data.columns), None)
        if pval_col:
            ont_data = ont_data.sort_values(pval_col)
        
        selected = []
        for _, row in ont_data.iterrows():
            genes_text = row.get('Genes', '')
            
            # Check similarity with selected terms
            is_redundant = False
            for sel in selected:
                sel_genes = set(str(sel.get('Genes', '')).split(','))
                curr_genes = set(str(genes_text).split(','))
                
                if len(curr_genes) > 0:
                    overlap = len(sel_genes & curr_genes) / len(curr_genes)
                    if overlap > similarity_cutoff:
                        is_redundant = True
                        break
            
            if not is_redundant:
                selected.append(row)
        
        simplified.extend(selected)
    
    return pd.DataFrame(simplified)


def export_go_report(
    enrich_result: pd.DataFrame,
    output_file: str,
    format: str = "excel"
) -> str:
    """Export GO enrichment results to file.
    
    Args:
        enrich_result: GO enrichment result
        output_file: Output file path
        format: Output format ('excel', 'csv', 'tsv')
    
    Returns:
        Output file path
    """
    _ensure_parent_dir(output_file)
    if format == "excel":
        enrich_result.to_excel(output_file, index=False)
    elif format == "csv":
        enrich_result.to_csv(output_file, index=False)
    else:
        enrich_result.to_csv(output_file, sep='\t', index=False)
    
    return output_file


# ========== KEGG Enrichment (bonus) ==========

def kegg_enrich(
    gene_list: Union[List[str], pd.DataFrame],
    organism: str = "mmu",
    pvalue_cutoff: float = 0.05,
    qvalue_cutoff: float = 0.05
) -> pd.DataFrame:
    """KEGG pathway enrichment analysis.
    
    Args:
        gene_list: List of gene IDs
        organism: KEGG organism code (e.g., 'mmu' for mouse, 'hsa' for human)
        pvalue_cutoff: P-value cutoff
        qvalue_cutoff: Q-value cutoff
    
    Returns:
        DataFrame with KEGG enrichment results
    """
    if not HAS_GSEAPY:
        return pd.DataFrame()
    
    # Handle gene_list format
    if isinstance(gene_list, pd.DataFrame):
        if 'genes' in gene_list.columns:
            gene_list = gene_list['genes'].tolist()
        elif len(gene_list) > 0:
            gene_list = gene_list.iloc[:, 0].tolist()
    
    gene_list = [str(g).strip() for g in gene_list if pd.notna(g)]
    gene_list = list(set(gene_list))
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets='KEGG_2019_Mouse',
            organism=_normalize_enrichr_organism(organism).lower(),
            outdir=None,
            no_plot=True
        )
        
        if enr is not None and hasattr(enr, 'results'):
            df = enr.results
            if 'P-value' in df.columns:
                df = df[df['P-value'] <= pvalue_cutoff]
            return df
    except Exception as e:
        warnings.warn(f"KEGG enrichment failed: {e}")
    
    return pd.DataFrame()


__all__ = [
    'go_enrich', 'local_go_enrich', 'gsea_enrich', 'go_plot', 'gsea_plot',
    'enrichr_library_list', 'cheak_orgdb', 'build_gene_set',
    'go_annotation_from_gtf', 'simplify_go_terms', 'export_go_report',
    'kegg_enrich'
]
