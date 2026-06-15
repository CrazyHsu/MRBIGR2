"""pheno_mcp — phenotype preprocessing MCP (9 tools)."""
from __future__ import annotations

import json
import sys
from pathlib import Path


def _bootstrap_syspath() -> None:
    here = Path(__file__).resolve()
    for candidate in here.parents:
        src_dir = candidate / "src" / "mrbigr"
        if src_dir.is_dir():
            src_parent = str(src_dir.parent)
            if src_parent not in sys.path:
                sys.path.insert(0, src_parent)
            return


_bootstrap_syspath()

from mrbigr.core import pheno  # noqa: E402


def _phenotype_frame(d):
    """Return a phenotype DataFrame from a CSV path or MCP JSON object."""
    import pandas as pd

    if isinstance(d, pd.DataFrame):
        df = d.copy()
    elif isinstance(d, (str, Path)):
        path = Path(d).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"phenotype CSV not found: {path}")
        df = pd.read_csv(path, compression="infer", sep=None, engine="python", index_col=0)
    else:
        df = pd.DataFrame(d)
        sample_col = _sample_id_column(df)
        if sample_col is not None:
            df = df.set_index(sample_col)

    if df.empty:
        raise ValueError("phenotype data is empty")
    df.index = df.index.astype(str)
    return df


def _derived_output_file(source, step_name):
    if not isinstance(source, (str, Path)):
        return None
    path = Path(source).expanduser()
    if not path.is_file():
        return None
    return path.with_name(f"{path.stem}.{step_name}.csv")


def _return_table_or_file(source, df, step_name):
    output_path = _derived_output_file(source, step_name)
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path)
        return {
            "output_file": str(output_path),
            "rows": int(df.shape[0]),
            "columns": [str(col) for col in df.columns],
        }
    return df.to_dict()


def _sample_id_column(df):
    id_names = {
        "id",
        "iid",
        "sample",
        "sample_id",
        "sampleid",
        "genotype",
        "accession",
    }
    for col in df.columns:
        col_text = str(col).strip()
        if col_text.lower() in id_names or col_text.startswith("Unnamed:"):
            return col
    return None


def _coerce_column_list(columns):
    if columns is None:
        raise ValueError("columns is required")
    if isinstance(columns, str):
        text = columns.strip()
        if not text:
            raise ValueError("columns is empty")
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            parsed = None
        if isinstance(parsed, list):
            values = parsed
        elif "," in text:
            values = [part.strip() for part in text.split(",")]
        else:
            values = [text]
    elif isinstance(columns, (list, tuple, set)):
        values = list(columns)
    else:
        values = [columns]

    result = [str(value).strip() for value in values if str(value).strip()]
    if not result:
        raise ValueError("columns is empty")
    return result


def _normalize_column_name(value):
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def _resolve_columns(df, columns):
    requested = _coerce_column_list(columns)
    exact = {str(col): col for col in df.columns}
    normalized = {}
    for col in df.columns:
        normalized.setdefault(_normalize_column_name(col), col)

    aliases = {
        "100gw": ["100grainweight", "100grainwt", "hundredgrainweight"],
    }

    resolved = []
    for name in requested:
        if name in exact:
            col = exact[name]
        else:
            norm = _normalize_column_name(name)
            col = normalized.get(norm)
            if col is None:
                for alias in aliases.get(norm, []):
                    col = normalized.get(alias)
                    if col is not None:
                        break
        if col is None:
            available = ", ".join(str(col) for col in df.columns)
            raise KeyError(f"phenotype column not found: {name}; available columns: {available}")
        if col not in resolved:
            resolved.append(col)
    return resolved


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def select_phenotype_columns(d, columns, output_file=None):
        """Select phenotype columns by exact name, comma-separated names, or JSON list. `100GW` is resolved to `100grainweight` when present."""
        df = _phenotype_frame(d)
        selected = df.loc[:, _resolve_columns(df, columns)]
        if output_file:
            output_path = Path(output_file).expanduser()
            output_path.parent.mkdir(parents=True, exist_ok=True)
            selected.to_csv(output_path)
            return {
                "output_file": str(output_path),
                "rows": int(selected.shape[0]),
                "columns": [str(col) for col in selected.columns],
            }
        return selected.to_dict()

    @mcp.tool()
    def filter_abundance(d, abundance):
        """Filter features by minimum abundance threshold. `d` may be a CSV path or phenotype object."""
        result = pheno.abundance_filter(_phenotype_frame(d), float(abundance))
        return _return_table_or_file(d, result, f"abundance_{abundance}")

    @mcp.tool()
    def filter_missing(d, missing_ratio):
        """Filter features by missing ratio threshold. `d` may be a CSV path or phenotype object."""
        result = pheno.missing_filter(_phenotype_frame(d), float(missing_ratio))
        return _return_table_or_file(d, result, f"missing_{missing_ratio}")

    @mcp.tool()
    def scale_phenotype(d, method='zscore'):
        """Scale phenotype data (log2, log10, zscore, minmax, robust, boxcox). `d` may be a CSV path or phenotype object."""
        result = pheno.scale_wrapper(_phenotype_frame(d), method)
        return _return_table_or_file(d, result, f"scaled_{method}")

    @mcp.tool()
    def impute_phenotype(d, method='mean'):
        """Impute missing phenotype values. `d` may be a CSV path or phenotype object."""
        result = pheno.pheno_imputer(_phenotype_frame(d), method=method)
        return _return_table_or_file(d, result, f"imputed_{method}")

    @mcp.tool()
    def remove_outliers(d, method='zscore'):
        """Remove outliers from phenotype data. `d` may be a CSV path or phenotype object."""
        result = pheno.outlier(_phenotype_frame(d), method=method)
        return _return_table_or_file(d, result, f"outliers_{method}")

    @mcp.tool()
    def run_blup(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUP via REML mixed model. `d` may be a CSV path or phenotype object."""
        return pheno.blup(_phenotype_frame(d), method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def run_blue(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUE via OLS fixed-effects model. `d` may be a CSV path or phenotype object."""
        return pheno.blue(_phenotype_frame(d), method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def correct_trait(pc, y):
        """Correct phenotype using principal components."""
        import pandas as pd
        import numpy as np
        from mrbigr.core._argjson import maybe_json_loads
        pc = maybe_json_loads(pc)
        y = maybe_json_loads(y)
        return pheno.trait_correct(pd.DataFrame(pc), np.array(y)).tolist()


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("pheno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
