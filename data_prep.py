"""
One-time data preparation script.

Reads the raw GDSC files in `data/raw/` and produces a small, GitHub-friendly
processed CSV at `data/processed/erlotinib_dataset.csv` containing:

    - One row per cell line that has both an Erlotinib IC50 and an expression profile.
    - Metadata columns: COSMIC_ID, cell_line_name, ln_IC50.
    - Expression columns: the top ~1,500 most variable genes plus a curated set
      of canonical EGFR-pathway genes (so the notebook can show biological validation
      even if those genes aren't in the top-variance set).

Run once after downloading the raw files:

    python data_prep.py
"""
from pathlib import Path

import pandas as pd

# --- Paths ---
RAW_DIR        = Path('data/raw')
PROCESSED_DIR  = Path('data/processed')
DRUG_PATH      = RAW_DIR / 'GDSC2_fitted_dose_response.csv'
EXPR_GZ_PATH   = RAW_DIR / 'sanger_expression.txt.gz'
OUT_PATH       = PROCESSED_DIR / 'erlotinib_dataset.csv'

# --- Config ---
TARGET_DRUG     = 'Erlotinib'
TOP_VARIABLE_N  = 1500

# Canonical EGFR-pathway genes we want to retain regardless of variance,
# so the notebook can validate against known biology.
EGFR_PATHWAY = {
    'ENSG00000146648': 'EGFR',
    'ENSG00000133703': 'KRAS',
    'ENSG00000157764': 'BRAF',
    'ENSG00000141736': 'ERBB2',
    'ENSG00000065361': 'ERBB3',
    'ENSG00000105976': 'MET',
    'ENSG00000121879': 'PIK3CA',
    'ENSG00000171862': 'PTEN',
}


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print(f'Loading {DRUG_PATH.name} ...')
    drug_df = pd.read_csv(DRUG_PATH)
    erl = drug_df[drug_df['DRUG_NAME'] == TARGET_DRUG]
    erl = (
        erl.groupby(['COSMIC_ID', 'CELL_LINE_NAME'], as_index=False)['LN_IC50']
           .mean()
           .rename(columns={'LN_IC50': 'ln_IC50', 'CELL_LINE_NAME': 'cell_line_name'})
    )
    erl['COSMIC_ID'] = erl['COSMIC_ID'].astype(int)
    print(f'  Cell lines screened against {TARGET_DRUG}: {len(erl)}')

    print(f'Loading {EXPR_GZ_PATH.name} ...')
    # Tab-separated, genes in rows, cell lines in columns (prefixed 'DATA.<COSMIC_ID>')
    expr = pd.read_csv(EXPR_GZ_PATH, sep='\t', index_col=0)
    expr.columns = [c.replace('DATA.', '') for c in expr.columns]
    expr = expr.loc[:, expr.columns.str.isdigit()]
    expr.columns = expr.columns.astype(int)
    print(f'  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} cell lines')

    print('Selecting top variable genes + EGFR pathway genes ...')
    variances     = expr.var(axis=1)
    top_var_genes = variances.nlargest(TOP_VARIABLE_N).index.tolist()
    pathway_in    = [g for g in EGFR_PATHWAY if g in expr.index]
    selected      = sorted(set(top_var_genes) | set(pathway_in))
    print(f'  Top variable: {len(top_var_genes)} | EGFR pathway retained: {len(pathway_in)} | union: {len(selected)}')

    expr_subset = expr.loc[selected].T  # rows = cell lines, cols = genes
    expr_subset.index.name = 'COSMIC_ID'
    expr_subset = expr_subset.round(4)  # smaller CSV

    print('Merging drug response with expression ...')
    merged = erl.merge(expr_subset, left_on='COSMIC_ID', right_index=True, how='inner')
    cols = ['COSMIC_ID', 'cell_line_name', 'ln_IC50'] + [c for c in merged.columns if c.startswith('ENSG')]
    merged = merged[cols]
    print(f'  Final dataset: {merged.shape[0]} cell lines x {merged.shape[1]} columns ({merged.shape[1] - 3} genes)')

    merged.to_csv(OUT_PATH, index=False)
    size_mb = OUT_PATH.stat().st_size / 1e6
    print(f'\nWrote {OUT_PATH} ({size_mb:.1f} MB)')


if __name__ == '__main__':
    main()
