#!/usr/bin/env python3
"""Patient-level paired tumor-normal validation of frozen macrophage programs."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import wilcoxon
import statsmodels.api as sm
import yaml


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GSE114725_rna_raw.csv.gz"
MARKER_SCORES = ROOT / "metadata" / "GSE114725_marker_program_scores_raw.tsv.gz"
LOCK_FILE = ROOT / "config" / "locked_macrophage_programs_v1.yml"
DISCOVERY = ROOT / "results" / "tables" / "GSE176078_macrophage_nmf_program_genes_pilot.tsv"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
METADATA = ROOT / "metadata"
INTERIM = ROOT / "data" / "interim"
LOGS = ROOT / "logs"

META = ["patient", "tissue", "replicate", "cluster", "cellid"]
LINEAGE_THRESHOLD = 0.4


def score_genes(matrix, var_names: pd.Index, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    lookup = {gene: index for index, gene in enumerate(var_names.astype(str))}
    present = [gene for gene in genes if gene in lookup]
    if not present:
        raise ValueError(f"No genes present from requested module: {genes}")
    columns = [lookup[gene] for gene in present]
    return np.asarray(matrix[:, columns].mean(axis=1)).ravel(), present


def fit_patient_level_model(frame: pd.DataFrame) -> pd.DataFrame:
    model_frame = frame.loc[frame["tissue"].isin(["TUMOR", "NORMAL"])].copy()
    model_frame["tumor"] = model_frame["tissue"].eq("TUMOR").astype(int)
    predictors = ["tumor", "hypoxia", "macrophage_identity", "median_log10_total_counts"]
    x = sm.add_constant(model_frame[predictors].astype(float))
    model = sm.OLS(model_frame["P4_primary"].astype(float), x).fit()
    rows = []
    for term in model.params.index:
        rows.append(
            {
                "term": term,
                "estimate": model.params[term],
                "standard_error": model.bse[term],
                "t_value": model.tvalues[term],
                "p_value": model.pvalues[term],
                "n_patient_tissues": model.nobs,
                "r_squared": model.rsquared,
                "note": "patient-tissue aggregate exploratory OLS; small sample",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    with LOCK_FILE.open() as handle:
        lock = yaml.safe_load(handle)
    discovery = pd.read_csv(DISCOVERY, sep="\t")
    discovery_controls = {
        "P1_inflammatory": discovery.loc[
            discovery["program"].eq("P1") & discovery["rank"].le(30), "gene"
        ].astype(str).tolist(),
        "P6_interferon": discovery.loc[
            discovery["program"].eq("P6") & discovery["rank"].le(30), "gene"
        ].astype(str).tolist(),
    }
    modules = {
        "P4_primary": lock["primary_program"]["genes"],
        "P3_secondary": lock["secondary_program"]["genes"],
        **discovery_controls,
        "hypoxia": [
            "HIF1A", "EGLN1", "CA9", "VEGFA", "LDHA", "PDK1", "SLC2A1",
            "BNIP3", "NDRG1", "ENO1", "ALDOA", "PGK1", "HK2", "PFKP",
        ],
        "macrophage_identity": ["C1QA", "C1QB", "C1QC", "CD68", "MSR1"],
        "phagolysosome": ["FCER1G", "TYROBP", "LST1", "AIF1", "CTSB", "LAMP1"],
    }

    markers = pd.read_csv(MARKER_SCORES, sep="\t")
    eligible = (
        markers["tissue"].isin(["TUMOR", "NORMAL"])
        & markers["macrophage_lineage_detected_fraction"].ge(LINEAGE_THRESHOLD)
        & markers["t_cell_detected_fraction"].lt(0.3)
        & markers["b_cell_detected_fraction"].lt(0.5)
    )
    selected_rows = set(markers.index[eligible].astype(int))
    raw = pd.read_csv(
        RAW,
        skiprows=lambda row_number: row_number > 0 and (row_number - 1) not in selected_rows,
        low_memory=False,
    )
    gene_columns = [column for column in raw.columns if column not in META]
    counts = sparse.csr_matrix(raw[gene_columns].to_numpy(dtype=np.float32))
    obs = raw[META].copy()
    obs.index = [f"GSE114725_pair_{i}" for i in range(obs.shape[0])]
    var = pd.DataFrame(index=pd.Index(gene_columns, name="gene"))
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    adata.obs["detected_genes"] = np.asarray((counts > 0).sum(axis=1)).ravel()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)

    coverage_rows = []
    for name, genes in modules.items():
        values, present = score_genes(adata.X, adata.var_names, genes)
        adata.obs[name] = values
        coverage_rows.append(
            {
                "module": name,
                "requested_genes": len(genes),
                "present_genes": len(present),
                "coverage": len(present) / len(genes),
                "genes_present": ";".join(present),
            }
        )
    pd.DataFrame(coverage_rows).to_csv(
        TABLES / "GSE114725_locked_program_module_coverage.tsv", sep="\t", index=False
    )

    score_columns = list(modules)
    cell_output = adata.obs[META + ["total_counts", "detected_genes"] + score_columns].copy()
    cell_output.to_csv(
        METADATA / "GSE114725_locked_program_cell_scores.tsv.gz",
        sep="\t",
        compression="gzip",
        index=True,
    )

    patient_tissue = (
        adata.obs.groupby(["patient", "tissue"], observed=True)[score_columns]
        .median()
        .reset_index()
    )
    qc = (
        adata.obs.assign(log10_total_counts=np.log10(adata.obs["total_counts"].clip(lower=1)))
        .groupby(["patient", "tissue"], observed=True)
        .agg(
            macrophage_cells=("cellid", "size"),
            median_total_counts=("total_counts", "median"),
            median_detected_genes=("detected_genes", "median"),
            median_log10_total_counts=("log10_total_counts", "median"),
        )
        .reset_index()
    )
    patient_tissue = patient_tissue.merge(qc, on=["patient", "tissue"], how="left")
    patient_tissue.to_csv(
        TABLES / "GSE114725_locked_program_scores_by_patient_tissue.tsv", sep="\t", index=False
    )

    paired_patients = sorted(
        set(patient_tissue.loc[patient_tissue["tissue"].eq("TUMOR"), "patient"])
        & set(patient_tissue.loc[patient_tissue["tissue"].eq("NORMAL"), "patient"])
    )
    paired = patient_tissue.loc[
        patient_tissue["patient"].isin(paired_patients)
        & patient_tissue["tissue"].isin(["TUMOR", "NORMAL"])
    ]
    test_rows = []
    for module in score_columns:
        wide = paired.pivot(index="patient", columns="tissue", values=module).dropna()
        differences = wide["TUMOR"] - wide["NORMAL"]
        test = wilcoxon(
            wide["TUMOR"], wide["NORMAL"], alternative="two-sided", method="exact"
        )
        test_rows.append(
            {
                "module": module,
                "paired_patients": wide.shape[0],
                "patients_tumor_greater": int((differences > 0).sum()),
                "median_tumor": wide["TUMOR"].median(),
                "median_normal": wide["NORMAL"].median(),
                "median_paired_difference": differences.median(),
                "mean_paired_difference": differences.mean(),
                "wilcoxon_statistic": test.statistic,
                "wilcoxon_exact_p": test.pvalue,
            }
        )
    paired_tests = pd.DataFrame(test_rows)
    paired_tests.to_csv(
        TABLES / "GSE114725_locked_program_paired_tumor_normal_tests.tsv", sep="\t", index=False
    )

    # Frozen primary program leave-one-gene-out sensitivity.
    p4_genes = lock["primary_program"]["genes"]
    loo_rows = []
    for omitted in [None] + p4_genes:
        retained = [gene for gene in p4_genes if gene != omitted]
        values, present = score_genes(adata.X, adata.var_names, retained)
        temp = adata.obs[["patient", "tissue"]].copy()
        temp["score"] = values
        aggregate = temp.groupby(["patient", "tissue"], observed=True)["score"].median().reset_index()
        wide = aggregate.loc[aggregate["patient"].isin(paired_patients)].pivot(
            index="patient", columns="tissue", values="score"
        ).dropna()
        diff = wide["TUMOR"] - wide["NORMAL"]
        loo_rows.append(
            {
                "omitted_gene": "none" if omitted is None else omitted,
                "retained_genes": len(present),
                "patients_tumor_greater": int((diff > 0).sum()),
                "median_paired_difference": diff.median(),
                "minimum_patient_difference": diff.min(),
                "maximum_patient_difference": diff.max(),
            }
        )
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(
        TABLES / "GSE114725_P4_leave_one_gene_out.tsv", sep="\t", index=False
    )

    regression = fit_patient_level_model(patient_tissue)
    regression.to_csv(
        TABLES / "GSE114725_P4_patient_tissue_adjusted_model.tsv", sep="\t", index=False
    )

    p4_plot = paired.pivot(index="patient", columns="tissue", values="P4_primary").dropna()
    fig, ax = plt.subplots(figsize=(4.2, 4.5))
    for patient, row in p4_plot.iterrows():
        ax.plot([0, 1], [row["NORMAL"], row["TUMOR"]], marker="o", linewidth=1.5, label=patient)
    ax.set_xticks([0, 1], ["Adjacent normal", "Tumor"])
    ax.set_ylabel("Frozen P4 score\nmedian log1p(CP10K)")
    ax.set_title("Patient-paired macrophage program")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=7, title="Patient")
    fig.tight_layout()
    fig.savefig(FIGURES / "GSE114725_P4_paired_tumor_normal.pdf")
    plt.close(fig)

    adata.write_h5ad(INTERIM / "GSE114725_tumor_normal_lineage_gated_macrophages.h5ad", compression="gzip")
    with (LOGS / "12_gse114725_paired_locked_programs.log").open("w") as handle:
        handle.write(f"lock_sha256\t{lock['source_table_sha256']}\n")
        handle.write(f"cells\t{adata.n_obs}\n")
        handle.write(f"patients\t{adata.obs['patient'].nunique()}\n")
        handle.write(f"paired_patients\t{len(paired_patients)}\n")
        handle.write("statistical_unit\tpatient\n")

    print(paired_tests.to_string(index=False))
    print("\nAdjusted patient-tissue model")
    print(regression.to_string(index=False))
    print("\nP4 leave-one-gene-out")
    print(loo.to_string(index=False))


if __name__ == "__main__":
    main()
