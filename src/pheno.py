#!/usr/bin/env python3
"""
PhenoMCP - Phenotype Analysis Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/pheno.py
"""
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import norm
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
import warnings

warnings.filterwarnings("ignore")


def mean(d):
    """Calculate row-wise mean."""
    return d.mean(axis=1, skipna=True)


def transpose(d):
    """Transpose dataframe."""
    return d.T


def abundance_filter(d, abundance):
    """Filter features by minimum abundance threshold."""
    return d.loc[:, d.mean(skipna=True) >= abundance]


def isDigit(x):
    """Check if value is a digit."""
    try:
        float(x)
        return True
    except ValueError:
        return False


def missing_filter(d, missing_ratio):
    """Filter features by missing ratio threshold."""
    if d.map(np.isreal).all().sum() == d.shape[1]:
        d = d.loc[:, d.map(np.isnan).sum() <= d.shape[0] * missing_ratio]
        return d
    else:
        digit_mask = d.map(lambda x: isDigit(x))
        d = d.loc[:, digit_mask.sum() >= d.shape[0] * (1 - missing_ratio)]
        d_array = d.values.copy().astype(object)
        d_array[~digit_mask[d.columns].values] = np.nan
        d = pd.DataFrame(d_array, index=d.index, columns=d.columns)
        d = d.astype(float)
        return d


def log2_scale(d):
    """Log2 transformation."""
    return np.log2(d + 1)


def ln_scale(d):
    """Natural log transformation."""
    return np.log(d + 1)


def log10_scale(d):
    """Log10 transformation."""
    return np.log10(d + 1)


def boxcox_scale(d):
    """Box-Cox transformation."""
    result = d.copy()
    for col in result.columns:
        try:
            valid_mask = d[col].notna()
            transformed, _ = stats.boxcox(d[col].dropna().values + 1)
            result[col] = np.nan
            result.loc[valid_mask, col] = transformed
        except:
            pass
    return result


def minmax_scale(d):
    """Min-max normalization."""
    d_copy = d.copy()
    d_copy.loc[:, :] = MinMaxScaler().fit_transform(d_copy.values)
    return d_copy


def zscore_scale(d):
    """Z-score standardization."""
    d_copy = d.copy()
    d_copy.loc[:, :] = StandardScaler().fit_transform(d_copy.values)
    return d_copy


def robust_scale(d):
    """Robust scaling."""
    d_copy = d.copy()
    d_copy.loc[:, :] = RobustScaler().fit_transform(d_copy.values)
    return d_copy


def ppoints(n, a=None):
    """Calculate plotting positions for QQ plot."""
    try:
        n = float(len(n))
    except TypeError:
        n = float(n)
    if a is None:
        a = 3.0 / 8 if (n <= 10) else 1.0 / 2
    return (np.arange(n) + 1 - a) / (n + 1 - 2 * a)


def qqnorm(y):
    """Normal QQ transformation."""
    ina = np.isnan(y)
    if ina.sum() > 0:
        yN = y
        y = y[~ina]
    n = y.shape[0]
    if n == 0:
        print('y is empty or has only NAs')
        return np.array([])
    x = np.around(norm.ppf(ppoints(n)[np.argsort(np.argsort(y))]), decimals=15)
    if ina.sum() > 0:
        y = x
        x = yN
        x[~ina] = y
    return x


def pheno_imputer(d, method='mean'):
    """Impute missing values."""
    imputer = SimpleImputer(missing_values=np.nan, strategy=method)
    d_copy = d.copy()
    d_copy.loc[:, :] = imputer.fit_transform(d_copy.values)
    return d_copy


def trait_correct(pc, y):
    """Correct phenotype using principal components.
    
    Args:
        pc: Principal components DataFrame
        y: Phenotype values
    
    Returns:
        Corrected phenotype values
    """
    # Add intercept
    ones = pd.DataFrame(np.ones((y.shape[0], 1)), index=pc.index, columns=['intercept'])
    pc1 = pd.concat([ones, pc], axis=1)
    
    # Calculate coefficients using pseudoinverse
    vhat = np.dot(np.linalg.pinv(np.dot(pc1.T, pc1)), np.dot(pc1.T, y))
    
    if len(vhat.shape) == 1:
        y_corr = y - np.dot(pc, vhat[1:])
    else:
        y_corr = y - np.dot(pc, vhat[1:, :])
    
    return y_corr


def _coerce_numeric_frame(d):
    """Coerce a phenotype frame to numeric values where possible."""
    result = d.copy()
    for col in result.columns:
        result[col] = pd.to_numeric(result[col], errors='coerce')
    return result


def _matrix_to_long(d):
    """Convert a wide phenotype matrix (rows=genotype, cols=env/rep) to long format."""
    df = _coerce_numeric_frame(d)
    idx_name = d.index.name or 'genotype'
    long = (
        df.reset_index()
        .melt(id_vars=idx_name, var_name='env', value_name='y')
        .rename(columns={idx_name: 'genotype'})
    )
    long = long.dropna(subset=['y']).copy()
    long['genotype'] = long['genotype'].astype(str)
    long['env'] = long['env'].astype(str)
    long['y'] = pd.to_numeric(long['y'], errors='coerce')
    long = long.dropna(subset=['y'])
    return long


def _fit_mixedlm(model):
    """Fit MixedLM with a small fallback for convergence robustness."""
    try:
        return model.fit(reml=True, method='lbfgs', disp=False)
    except Exception:
        return model.fit(reml=True, disp=False)


def _blup_legacy(d):
    """Legacy pseudo-BLUP kept only for backward compatibility."""
    result = _coerce_numeric_frame(d)

    for col in result.columns:
        y = result[col].values.copy()
        valid = ~np.isnan(y)
        n = valid.sum()

        if n > 1:
            y_valid = y[valid]
            grand_mean = np.mean(y_valid)
            total_var = np.var(y_valid, ddof=1)

            sigma_g = total_var * 0.5
            sigma_e = total_var * 0.5
            h2 = sigma_g / (sigma_g + sigma_e) if (sigma_g + sigma_e) > 0 else 0

            blup_vals = y.copy()
            blup_vals[valid] = grand_mean + h2 * (y_valid - grand_mean)
            result[col] = blup_vals
        elif n == 1:
            result[col] = y

    return result


def _blue_legacy(d):
    """Legacy pseudo-BLUE kept only for backward compatibility."""
    result = _coerce_numeric_frame(d)

    for col in result.columns:
        y = result[col].values.copy()
        valid = ~np.isnan(y)

        if valid.sum() > 0:
            grand_mean = np.mean(y[valid])
            blue_vals = y.copy()
            blue_vals[~valid] = grand_mean
            result[col] = blue_vals

    return result


def blup(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
    """Best Linear Unbiased Prediction (BLUP).

    method='matrix':
        Input: rows=genotypes, cols=environments/replicates.
        Model: y_ij = mu + env_j(fixed) + g_i(random) + e_ij
        Output: DataFrame indexed by genotype with one column: 'blup'

    method='long':
        Input: long-format table with genotype, env, y columns.
        Model: y ~ C(env) [+ C(rep):C(env)] + (1|genotype)
        Output: DataFrame with columns [geno_col, 'blup']

    method='legacy':
        Keep the old shrinkage behavior for backward compatibility.
    """
    if method == 'matrix':
        return _blup_matrix(d)
    if method == 'long':
        return _blup_long(d, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col)
    if method == 'legacy':
        return _blup_legacy(d)
    raise ValueError(f"method must be 'matrix', 'long', or 'legacy', got '{method}'")


def _blup_matrix(d):
    import statsmodels.api as sm

    long = _matrix_to_long(d)

    if long['genotype'].nunique() < 2 or long['env'].nunique() < 2:
        raise ValueError("Need >= 2 genotypes and >= 2 environments for BLUP matrix mode.")

    model = sm.MixedLM.from_formula("y ~ C(env)", groups=long['genotype'], re_formula="1", data=long)
    fit = _fit_mixedlm(model)

    long = long.copy()
    long['fitted'] = fit.fittedvalues
    blup_means = long.groupby('genotype', sort=False)['fitted'].mean()

    out = pd.DataFrame(index=d.index.copy())
    out.index.name = d.index.name
    out['blup'] = [blup_means.get(str(g), np.nan) for g in d.index]
    return out


def _blup_long(d, y_col='y', geno_col='genotype', env_col='env', rep_col=None):
    import statsmodels.api as sm

    required = {y_col, geno_col, env_col}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Available: {list(d.columns)}")

    df = d.copy()
    df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
    df = df.dropna(subset=[y_col, geno_col, env_col]).copy()
    df[geno_col] = df[geno_col].astype(str)
    df[env_col] = df[env_col].astype(str)

    formula = f"{y_col} ~ C({env_col})"
    if rep_col and rep_col in df.columns:
        df[rep_col] = df[rep_col].astype(str)
        formula += f" + C({rep_col}):C({env_col})"

    if df[geno_col].nunique() < 2 or df[env_col].nunique() < 2:
        raise ValueError("Need >= 2 genotypes and >= 2 environments for BLUP long mode.")

    model = sm.MixedLM.from_formula(formula, groups=df[geno_col], re_formula="1", data=df)
    fit = _fit_mixedlm(model)

    df = df.copy()
    df['fitted'] = fit.fittedvalues
    blup_means = df.groupby(geno_col, sort=False)['fitted'].mean()

    return pd.DataFrame({geno_col: blup_means.index.tolist(), 'blup': blup_means.values.tolist()})


def blue(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
    """Best Linear Unbiased Estimation (BLUE).

    method='matrix':
        Input: rows=genotypes, cols=environments/replicates.
        Model: y_ij = mu + env_j + tau_i + e_ij
        Output: DataFrame indexed by genotype with one column: 'blue'

    method='long':
        Input: long-format table with genotype, env, y columns.
        Model: y ~ C(env) + C(genotype) [+ C(rep):C(env)]
        Output: DataFrame with columns [geno_col, 'blue']

    method='legacy':
        Keep the old mean-fill behavior for backward compatibility.
    """
    if method == 'matrix':
        return _blue_matrix(d)
    if method == 'long':
        return _blue_long(d, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col)
    if method == 'legacy':
        return _blue_legacy(d)
    raise ValueError(f"method must be 'matrix', 'long', or 'legacy', got '{method}'")


def _blue_matrix(d):
    import statsmodels.formula.api as smf

    long = _matrix_to_long(d)

    if long['genotype'].nunique() < 2 or long['env'].nunique() < 2:
        raise ValueError("Need >= 2 genotypes and >= 2 environments for BLUE matrix mode.")

    fit = smf.ols("y ~ C(env) + C(genotype)", data=long).fit()

    long = long.copy()
    long['fitted'] = fit.fittedvalues
    blue_means = long.groupby('genotype', sort=False)['fitted'].mean()

    out = pd.DataFrame(index=d.index.copy())
    out.index.name = d.index.name
    out['blue'] = [blue_means.get(str(g), np.nan) for g in d.index]
    return out


def _blue_long(d, y_col='y', geno_col='genotype', env_col='env', rep_col=None):
    import statsmodels.formula.api as smf

    required = {y_col, geno_col, env_col}
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Available: {list(d.columns)}")

    df = d.copy()
    df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
    df = df.dropna(subset=[y_col, geno_col, env_col]).copy()
    df[geno_col] = df[geno_col].astype(str)
    df[env_col] = df[env_col].astype(str)

    formula = f"{y_col} ~ C({env_col}) + C({geno_col})"
    if rep_col and rep_col in df.columns:
        df[rep_col] = df[rep_col].astype(str)
        formula += f" + C({rep_col}):C({env_col})"

    if df[geno_col].nunique() < 2 or df[env_col].nunique() < 2:
        raise ValueError("Need >= 2 genotypes and >= 2 environments for BLUE long mode.")

    fit = smf.ols(formula, data=df).fit()

    df = df.copy()
    df['fitted'] = fit.fittedvalues
    blue_means = df.groupby(geno_col, sort=False)['fitted'].mean()

    return pd.DataFrame({geno_col: blue_means.index.tolist(), 'blue': blue_means.values.tolist()})


def outlier(d, method='zscore'):
    """Remove outliers from phenotype data.
    
    Args:
        d: Phenotype DataFrame
        method: 'zscore' or 'boxplot'
    
    Returns:
        DataFrame with outliers removed
    """
    result = d.copy()
    
    if method == 'zscore':
        for col in result.columns:
            y = result[col].values
            valid = ~np.isnan(y)

            if valid.sum() > 1:
                mean = np.mean(y[valid])
                std = np.std(y[valid], ddof=1)

                if std > 0:
                    z_scores = np.abs((y - mean) / std)
                    result.loc[z_scores > 3, col] = np.nan
    
    elif method == 'boxplot':
        for col in result.columns:
            y = result[col].values
            valid = ~np.isnan(y)
            
            if valid.sum() > 0:
                q1 = np.percentile(y[valid], 25)
                q3 = np.percentile(y[valid], 75)
                iqr = q3 - q1
                
                lower = q1 - 1.5 * iqr
                upper = q3 + 1.5 * iqr
                
                # Set outliers to NaN
                result.loc[(y < lower) | (y > upper), col] = np.nan
    
    return result


def scale_wrapper(d, method):
    """Wrapper for different scaling methods."""
    if method == 'log2':
        return log2_scale(d)
    elif method == 'ln':
        return ln_scale(d)
    elif method == 'log10':
        return log10_scale(d)
    elif method == 'boxcox':
        return boxcox_scale(d)
    elif method == 'minmax':
        return minmax_scale(d)
    elif method == 'zscore':
        return zscore_scale(d)
    elif method == 'robust':
        return robust_scale(d)
    else:
        return d


__all__ = [
    'mean', 'transpose', 'abundance_filter', 'missing_filter',
    'log2_scale', 'ln_scale', 'log10_scale', 'boxcox_scale',
    'minmax_scale', 'zscore_scale', 'robust_scale',
    'qqnorm', 'pheno_imputer', 'trait_correct',
    'blup', 'blue', 'outlier', 'scale_wrapper'
]
