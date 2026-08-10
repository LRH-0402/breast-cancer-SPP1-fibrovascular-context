#!/usr/bin/env python3
"""Audit GSE246613 metadata without using clinical labels for feature discovery."""

from __future__ import annotations

import contextlib
from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
TABLES = ROOT / "results" / "tables"
METADATA = ROOT / "metadata"
LOGS = ROOT / "logs"
for folder in (TABLES, METADATA, LOGS):
    folder.mkdir(parents=True, exist_ok=True)


def count_table(obs: pd.DataFrame, rows: list[str], value_name: str = "cells") -> pd.DataFrame:
    result = obs.groupby(rows, observed=True).size().rename(value_name).reset_index()
    return result.sort_values(rows).reset_index(drop=True)


def main() -> None:
    adata = sc.read_h5ad(INPUT, backed="r")
    obs = adata.obs.copy()

    required = [
        "celltype", "subcluster", "myeloid_leiden_nbr30_res0.8", "cohort",
        "patient_treatment", "treatment", "pCR", "RCB", "batch",
        "n_counts", "n_genes_by_counts", "percent_mito", "scrublet",
    ]
    missing = sorted(set(required) - set(obs.columns))
    if missing:
        raise ValueError(f"Missing required metadata columns: {missing}")

    selected = obs[required].copy()
    selected.index.name = "cell_id"
    selected.to_csv(METADATA / "GSE246613_cell_metadata_selected.tsv.gz", sep="\t", compression="gzip")

    count_table(obs, ["celltype"]).to_csv(
        TABLES / "GSE246613_cells_by_major_type.tsv", sep="\t", index=False
    )
    count_table(obs, ["subcluster"]).to_csv(
        TABLES / "GSE246613_cells_by_subcluster.tsv", sep="\t", index=False
    )
    count_table(obs, ["cohort", "treatment", "pCR", "RCB", "celltype"]).to_csv(
        TABLES / "GSE246613_cells_by_patient_timepoint_outcome_type.tsv", sep="\t", index=False
    )

    myeloid = obs.loc[obs["celltype"].astype(str).eq("myeloid")].copy()
    count_table(myeloid, ["cohort", "treatment", "subcluster"]).to_csv(
        TABLES / "GSE246613_myeloid_cells_by_patient_timepoint_subcluster.tsv",
        sep="\t",
        index=False,
    )

    patient_timepoints = obs[["cohort", "treatment", "pCR", "RCB", "patient_treatment"]].drop_duplicates()
    patient_timepoints.to_csv(
        METADATA / "GSE246613_patient_timepoints.tsv", sep="\t", index=False
    )

    summary = pd.DataFrame(
        {
            "metric": [
                "cells", "genes", "patients", "patient_timepoints", "myeloid_cells",
                "patients_with_myeloid", "baseline_patients", "PD1_patients", "RTPD1_patients",
            ],
            "value": [
                adata.n_obs,
                adata.n_vars,
                obs["cohort"].nunique(),
                obs["patient_treatment"].nunique(),
                myeloid.shape[0],
                myeloid["cohort"].nunique(),
                obs.loc[obs["treatment"].astype(str).eq("Base"), "cohort"].nunique(),
                obs.loc[obs["treatment"].astype(str).eq("PD1"), "cohort"].nunique(),
                obs.loc[obs["treatment"].astype(str).eq("RTPD1"), "cohort"].nunique(),
            ],
        }
    )
    summary.to_csv(TABLES / "GSE246613_audit_summary.tsv", sep="\t", index=False)

    with (LOGS / "04_audit_gse246613_versions.log").open("w") as handle:
        with contextlib.redirect_stdout(handle):
            sc.logging.print_versions()

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
