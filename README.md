# Biomarker Discovery & Drug Sensitivity Prediction: A Machine Learning Approach

An end-to-end computational pipeline that predicts cancer cell line sensitivity to **Erlotinib** (an EGFR-targeted therapy) from gene expression data, using real public datasets from the **Genomics of Drug Sensitivity in Cancer (GDSC)** project at the Wellcome Sanger Institute.

The pipeline takes ~1,500 gene expression measurements per cell line, isolates a small set of statistically informative features, trains a machine learning model to predict drug response (IC50), and validates the biological plausibility of the model's top biomarkers against the canonical EGFR pathway.

---

## Objective

A central challenge in **computational precision oncology** is the "needle-in-a-haystack" problem: tumors are profiled across tens of thousands of genes, but only a small subset actively drive therapeutic response or resistance. This high noise-to-signal ratio makes naive modeling unreliable.

This project tackles that problem on real cell line data:

- Load gene expression and drug response profiles for hundreds of cancer cell lines.
- Apply **statistical feature selection** to compress thousands of candidate genes into a tractable set.
- Train a **regression model** to predict log-IC50 from the selected features.
- Inspect the top biomarkers and check whether they recover **known EGFR pathway biology** (EGFR, KRAS, BRAF, MET, ERBB2 — the canonical drivers of Erlotinib response and resistance).

---

## Repository Layout

```
.
├── README.md
├── data_prep.py                       # one-time preparation: raw GDSC files → processed CSV
├── drug_resistance_pipeline.ipynb     # the analysis notebook (loads the processed CSV)
├── .gitignore
└── data/
    ├── raw/                           # raw GDSC downloads (gitignored, see "Data Sources" below)
    │   ├── GDSC2_fitted_dose_response.csv
    │   └── sanger_expression.txt.gz
    └── processed/
        └── erlotinib_dataset.csv      # the ~8 MB committed dataset the notebook reads
```

The notebook is fully self-contained at run time — it only reads `data/processed/erlotinib_dataset.csv`, which is committed to the repository. The raw files (~175 MB) are not committed; they're recreated by re-downloading from the Sanger Institute (see below) if anyone wants to regenerate the processed CSV from scratch.

---

## Data Sources

| Dataset | Description | Source | Local path |
|---|---|---|---|
| GDSC2 Fitted Dose Response | ln(IC50) values for ~1,000 cell lines × ~250 drugs (we filter to Erlotinib) | [Sanger FTP — current release](https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_24Jul22.csv) | `data/raw/GDSC2_fitted_dose_response.csv` (38 MB) |
| GDSC RMA-normalized Expression | Basal gene expression matrix (~17,000 Ensembl genes × ~1,000 cell lines) | [Sanger FTP — release 6.0](https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/sanger1018_brainarray_ensemblgene_rma.txt.gz) | `data/raw/sanger_expression.txt.gz` (135 MB) |

To rebuild the processed CSV from scratch:

```bash
# 1. download the raw files into data/raw/
mkdir -p data/raw && cd data/raw
curl -O https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_24Jul22.csv
mv GDSC2_fitted_dose_response_24Jul22.csv GDSC2_fitted_dose_response.csv
curl -O https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/sanger1018_brainarray_ensemblgene_rma.txt.gz
mv sanger1018_brainarray_ensemblgene_rma.txt.gz sanger_expression.txt.gz

# 2. run the prep script
cd ../.. && python data_prep.py
```

---

## Tech Stack & Libraries

| Category | Tools |
|---|---|
| Language | Python 3.10+ |
| Data Manipulation | `pandas`, `numpy` |
| Machine Learning | `scikit-learn` (SelectKBest, f_regression, RandomForestRegressor) |
| Visualization | `matplotlib`, `seaborn` |
| Environment | Jupyter Notebook |

```bash
pip install pandas numpy scikit-learn matplotlib seaborn jupyter
```

---

## Workflow

**One-time data prep** (`data_prep.py`):
1. Filter raw GDSC2 drug response to Erlotinib only.
2. Load the RMA-normalized expression matrix.
3. Select the top 1,500 most variable genes plus a curated set of EGFR-pathway genes (so biological validation is always possible).
4. Merge expression with Erlotinib IC50 by COSMIC cell line ID; write the result to `data/processed/erlotinib_dataset.csv`.

**Analysis notebook** (`drug_resistance_pipeline.ipynb`):
1. Load the processed CSV.
2. Inspect the IC50 distribution across the cell-line panel.
3. **Feature selection** — `SelectKBest` with F-regression to compress ~1,500 genes into the top 50 associated with Erlotinib response.
4. **Model training** — `RandomForestRegressor` on a standard train/test split.
5. **Evaluation** — MSE, R², and a predicted-vs-actual scatter plot.
6. **Biological validation** — rank gene importance and cross-reference against known EGFR pathway biology.

---

## Expected Results & Honest Limitations

Unlike toy synthetic data, **real drug response is noisy** — many factors beyond gene expression (mutations, copy number, tissue lineage, microenvironment) drive Erlotinib response. R² in the 0.10–0.30 range is realistic and consistent with published GDSC modeling work. The value of this project is *not* a high R² — it is **demonstrating the full pipeline** and showing that the top selected biomarkers are biologically interpretable.

---

## Key Resume Takeaway

This project demonstrates the **end-to-end computational thinking** required for a precision oncology research team:

- **Working with real public cancer data** — the same GDSC datasets used in dozens of published precision-oncology studies.
- **Handling high-dimensional, noisy biological data** — the dominant challenge in genomic and single-cell tumor analysis.
- **Building reproducible pipelines** — a clear `raw → processed → analysis` flow with a documented prep step.
- **Connecting model interpretability to biology** — feature importances are not just metrics; they are *candidate biomarkers* that can inform therapy selection and reveal mechanisms of resistance.
- **Setting honest expectations** — understanding that real biology imposes ceilings on predictive accuracy and that biological plausibility matters as much as raw model performance.

The workflow mirrors how research labs identify drug-resistance signatures in patient tumors — making this a realistic showcase of skills relevant to multimodal tumor atlas building, spatial biomarker discovery, and ML-driven therapeutic response modeling.
