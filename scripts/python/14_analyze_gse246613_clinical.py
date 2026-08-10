#!/usr/bin/env python3
"""Execute the frozen GSE246613 patient-level clinical analysis."""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2, mannwhitneyu, wilcoxon
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
import yaml


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "raw" / "GSE246613_PembroRT_immune_R100_final.h5ad"
SCORES = ROOT / "results" / "tables" / "GSE246613_locked_program_by_patient_timepoint_blinded.tsv"
SPEC = ROOT / "config" / "gse246613_clinical_analysis_spec_v1.yml"
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
LOGS = ROOT / "logs"
RNG_SEED = 20260721


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def zscore(values: pd.Series) -> pd.Series:
    sd = values.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardize {values.name}: non-positive SD")
    return (values - values.mean()) / sd


def fit_logistic_lrt(frame: pd.DataFrame, predictor: str, covariate: str | None = None) -> dict:
    y = frame["response"].astype(float)
    predictors = [predictor] + ([] if covariate is None else [covariate])
    x_full = sm.add_constant(frame[predictors].astype(float), has_constant="add")
    full = sm.GLM(y, x_full, family=sm.families.Binomial()).fit()
    reduced_predictors = [] if covariate is None else [covariate]
    x_reduced = sm.add_constant(frame[reduced_predictors].astype(float), has_constant="add")
    reduced = sm.GLM(y, x_reduced, family=sm.families.Binomial()).fit()
    lrt = 2 * (full.llf - reduced.llf)
    coefficient = full.params[predictor]
    standard_error = full.bse[predictor]
    return {
        "predictor": predictor,
        "covariate": "none" if covariate is None else covariate,
        "n": frame.shape[0],
        "responders": int(y.sum()),
        "nonresponders": int((1 - y).sum()),
        "coefficient": coefficient,
        "odds_ratio_per_sd": np.exp(coefficient),
        "ci95_low": np.exp(coefficient - 1.96 * standard_error),
        "ci95_high": np.exp(coefficient + 1.96 * standard_error),
        "wald_p": full.pvalues[predictor],
        "likelihood_ratio_chi2": lrt,
        "likelihood_ratio_p": chi2.sf(lrt, df=1),
        "aic": full.aic,
    }


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def auc_robustness(y: np.ndarray, score: np.ndarray) -> dict:
    rng = np.random.default_rng(RNG_SEED)
    observed = roc_auc_score(y, score)
    permuted = np.empty(10_000)
    for index in range(permuted.size):
        permuted[index] = roc_auc_score(rng.permutation(y), score)
    permutation_p = (1 + np.sum(np.abs(permuted - 0.5) >= abs(observed - 0.5))) / (
        permuted.size + 1
    )

    responder_indices = np.flatnonzero(y == 1)
    nonresponder_indices = np.flatnonzero(y == 0)
    bootstrapped = np.empty(5_000)
    for index in range(bootstrapped.size):
        sampled = np.concatenate(
            [
                rng.choice(responder_indices, size=len(responder_indices), replace=True),
                rng.choice(nonresponder_indices, size=len(nonresponder_indices), replace=True),
            ]
        )
        bootstrapped[index] = roc_auc_score(y[sampled], score[sampled])
    return {
        "auc_for_response": observed,
        "auc_bootstrap_ci95_low": np.quantile(bootstrapped, 0.025),
        "auc_bootstrap_ci95_high": np.quantile(bootstrapped, 0.975),
        "auc_two_sided_permutation_p": permutation_p,
        "permutations": permuted.size,
        "bootstrap_resamples": bootstrapped.size,
    }


def main() -> None:
    with SPEC.open() as handle:
        specification = yaml.safe_load(handle)
    actual_hash = file_sha256(SCORES)
    expected_hash = specification["blinded_score_table_sha256"]
    if actual_hash != expected_hash:
        raise RuntimeError(f"Blinded score hash mismatch: {actual_hash} != {expected_hash}")

    scores = pd.read_csv(SCORES, sep="\t")
    adata = sc.read_h5ad(INPUT, backed="r")
    outcomes = adata.obs[["cohort", "pCR", "RCB"]].copy()
    outcomes["pCR"] = outcomes["pCR"].astype(str)
    outcomes["RCB"] = outcomes["RCB"].astype(str)
    outcomes = outcomes.drop_duplicates()
    inconsistent = outcomes.groupby("cohort", observed=True).size()
    if (inconsistent > 1).any():
        raise ValueError(f"Inconsistent patient outcomes: {inconsistent[inconsistent > 1].to_dict()}")

    baseline = scores.loc[scores["treatment"].eq("Base")].merge(
        outcomes, on="cohort", how="left", validate="one_to_one"
    )
    baseline = baseline.loc[baseline["pCR"].isin(["R", "NR"])].copy()
    baseline["response"] = baseline["pCR"].map({"R": 1, "NR": 0}).astype(int)

    model_variables = [
        "P4_primary", "P3_secondary", "P1_inflammatory", "P6_interferon",
        "hypoxia", "macrophage_identity", "phagolysosome",
        "macrophage_like_fraction_of_immune", "median_n_counts",
    ] + [column for column in baseline if column.startswith("P4_loo_")]
    for variable in model_variables:
        baseline[f"z_{variable}"] = zscore(baseline[variable].astype(float))

    primary = fit_logistic_lrt(baseline, "z_P4_primary")
    auc_stats = auc_robustness(
        baseline["response"].to_numpy(), baseline["z_P4_primary"].to_numpy()
    )
    responders = baseline.loc[baseline["response"].eq(1), "z_P4_primary"]
    nonresponders = baseline.loc[baseline["response"].eq(0), "z_P4_primary"]
    rank_test = mannwhitneyu(responders, nonresponders, alternative="two-sided")
    primary.update(auc_stats)
    primary.update(
        {
            "responder_median_zscore": responders.median(),
            "nonresponder_median_zscore": nonresponders.median(),
            "mannwhitney_u": rank_test.statistic,
            "mannwhitney_p": rank_test.pvalue,
            "score_table_sha256": actual_hash,
        }
    )
    pd.DataFrame([primary]).to_csv(
        TABLES / "GSE246613_P4_primary_pCR_analysis.tsv", sep="\t", index=False
    )

    adjustment_rows = []
    for covariate in [
        "macrophage_identity", "hypoxia", "macrophage_like_fraction_of_immune", "median_n_counts"
    ]:
        adjustment_rows.append(
            fit_logistic_lrt(baseline, "z_P4_primary", f"z_{covariate}")
        )
    adjustments = pd.DataFrame(adjustment_rows)
    adjustments.to_csv(
        TABLES / "GSE246613_P4_pCR_one_at_a_time_adjustments.tsv", sep="\t", index=False
    )

    negative_rows = []
    for variable in [
        "P1_inflammatory", "P6_interferon", "P3_secondary", "macrophage_identity",
        "phagolysosome", "hypoxia", "macrophage_like_fraction_of_immune",
    ]:
        result = fit_logistic_lrt(baseline, f"z_{variable}")
        result["negative_control"] = variable
        negative_rows.append(result)
    negative = pd.DataFrame(negative_rows)
    negative["likelihood_ratio_fdr"] = bh(negative["likelihood_ratio_p"])
    negative.to_csv(
        TABLES / "GSE246613_pCR_negative_control_models.tsv", sep="\t", index=False
    )

    loo_rows = []
    for column in [column for column in baseline if column.startswith("P4_loo_")]:
        result = fit_logistic_lrt(baseline, f"z_{column}")
        result["omitted_gene"] = column.removeprefix("P4_loo_")
        loo_rows.append(result)
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(TABLES / "GSE246613_P4_pCR_leave_one_gene_out.tsv", sep="\t", index=False)

    joined = scores.merge(outcomes, on="cohort", how="left", validate="many_to_one")
    joined.to_csv(TABLES / "GSE246613_locked_program_clinical_join.tsv", sep="\t", index=False)
    wide = joined.pivot(index="cohort", columns="treatment", values="P4_primary")
    response_map = outcomes.set_index("cohort")["pCR"]
    longitudinal_rows = []
    for start, end in [("Base", "PD1"), ("PD1", "RTPD1"), ("Base", "RTPD1")]:
        paired = wide[[start, end]].dropna()
        change = paired[end] - paired[start]
        paired_test = wilcoxon(paired[end], paired[start], alternative="two-sided", method="auto")
        longitudinal_rows.append(
            {
                "contrast": f"{end}_minus_{start}",
                "test": "paired_wilcoxon_overall",
                "n": paired.shape[0],
                "median_change": change.median(),
                "statistic": paired_test.statistic,
                "p_value": paired_test.pvalue,
            }
        )
        groups = response_map.reindex(change.index)
        r_change = change.loc[groups.eq("R")]
        nr_change = change.loc[groups.eq("NR")]
        interaction = mannwhitneyu(r_change, nr_change, alternative="two-sided")
        longitudinal_rows.append(
            {
                "contrast": f"{end}_minus_{start}",
                "test": "mannwhitney_change_R_vs_NR",
                "n": change.shape[0],
                "median_change": r_change.median() - nr_change.median(),
                "statistic": interaction.statistic,
                "p_value": interaction.pvalue,
            }
        )
    longitudinal = pd.DataFrame(longitudinal_rows)
    longitudinal["fdr"] = bh(longitudinal["p_value"])
    longitudinal.to_csv(
        TABLES / "GSE246613_P4_longitudinal_response_analysis.tsv", sep="\t", index=False
    )

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    for x_position, label in enumerate(["NR", "R"]):
        values = baseline.loc[baseline["pCR"].eq(label), "P4_primary"]
        jitter = np.random.default_rng(RNG_SEED + x_position).normal(x_position, 0.045, len(values))
        axes[0].scatter(jitter, values, alpha=0.8, s=24)
        axes[0].plot([x_position - 0.18, x_position + 0.18], [values.median()] * 2, color="black")
    axes[0].set_xticks([0, 1], ["Non-response", "pCR"])
    axes[0].set_ylabel("Baseline frozen P4 score")
    axes[0].set_title("Patient-level baseline analysis")

    treatment_positions = {"Base": 0, "PD1": 1, "RTPD1": 2}
    for cohort, group in joined.groupby("cohort", observed=True):
        group = group.loc[group["treatment"].isin(treatment_positions)].copy()
        group["position"] = group["treatment"].map(treatment_positions)
        group = group.sort_values("position")
        color = "#2B6CB0" if group["pCR"].iloc[0] == "R" else "#C53030"
        axes[1].plot(group["position"], group["P4_primary"], color=color, alpha=0.28, linewidth=0.9)
    axes[1].set_xticks([0, 1, 2], ["Base", "PD1", "RT+PD1"])
    axes[1].set_ylabel("Frozen P4 score")
    axes[1].set_title("Longitudinal trajectories")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "GSE246613_P4_clinical_analysis.pdf")
    plt.close(fig)

    with (LOGS / "14_analyze_gse246613_clinical.log").open("w") as handle:
        handle.write(f"verified_blinded_score_sha256\t{actual_hash}\n")
        handle.write(f"baseline_patients_analyzed\t{baseline.shape[0]}\n")
        handle.write(f"responders\t{baseline['response'].sum()}\n")
        handle.write(f"nonresponders\t{(1 - baseline['response']).sum()}\n")
        handle.write("statistical_unit\tpatient\n")
        handle.write("feature_reselection_after_outcome_join\tfalse\n")

    print("Primary analysis")
    print(pd.DataFrame([primary]).to_string(index=False))
    print("\nOne-at-a-time adjustments")
    print(adjustments.to_string(index=False))
    print("\nNegative controls")
    print(negative.to_string(index=False))
    print("\nLongitudinal analyses")
    print(longitudinal.to_string(index=False))


if __name__ == "__main__":
    main()
