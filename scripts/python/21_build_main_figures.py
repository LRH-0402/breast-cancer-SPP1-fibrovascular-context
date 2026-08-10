#!/usr/bin/env python3
"""Build manuscript-ready multipanel Figures 1 and 2 from frozen tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr
import seaborn as sns
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures" / "main"
METADATA = ROOT / "metadata"
FIGSHARE = ROOT / "data" / "derived" / "figshare_21591429_scores"

# Okabe-Ito-derived palette: distinguishable in common color-vision deficiencies.
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
GREY = "#667085"
DIVERGING = "PuOr_r"


def panel_label(axis, label: str) -> None:
    axis.text(-0.16, 1.08, label, transform=axis.transAxes, fontsize=10,
              fontweight="bold", va="top")


def save(figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    # Exact 178-mm double-column artwork; PDF remains vector and PNG is 600 dpi.
    figure.savefig(FIGURES / f"{stem}.pdf", facecolor="white")
    figure.savefig(FIGURES / f"{stem}.png", dpi=600, facecolor="white")
    figure.savefig(FIGURES / f"{stem}.tiff", dpi=600, facecolor="white",
                   pil_kwargs={"compression": "tiff_lzw"})
    plt.close(figure)


def build_figure1() -> None:
    match = pd.read_csv(TABLES / "gate1_cross_cohort_program_matching_summary.tsv", sep="\t")
    core = pd.read_csv(TABLES / "locked_macrophage_program_core_genes_v1.tsv", sep="\t")
    exposure = pd.read_csv(TABLES / "GSE176078_macrophage_nmf_patient_scores_pilot.tsv", sep="\t")
    paired = pd.read_csv(TABLES / "GSE114725_locked_program_scores_by_patient_tissue.tsv", sep="\t")

    figure = plt.figure(figsize=(7.20, 5.25), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, width_ratios=[1.15, 1, 1.1], height_ratios=[1, 1.05])

    axis = figure.add_subplot(grid[0, 0])
    axis.axis("off")
    panel_label(axis, "A")
    boxes = [
        (0.03, 0.69, 0.94, 0.23, "Discovery", "GSE176078\n26 patients | 5,929 macrophages"),
        (0.03, 0.38, 0.94, 0.23, "Independent factorization", "GSE161529 | GSE246613 | GSE114725"),
        (0.03, 0.07, 0.94, 0.23, "Locked programs", "P4 primary: 13 genes\nP3 secondary: 8 genes"),
    ]
    for x, y, width, height, title, body in boxes:
        axis.add_patch(mpl.patches.FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.012", facecolor="#F4F7FA",
            edgecolor="#93A4B8", linewidth=1,
        ))
        axis.text(x + 0.04, y + height - 0.045, title, fontsize=7.2,
                  fontweight="bold", va="top")
        axis.text(x + 0.04, y + 0.035, body, fontsize=6.1, va="bottom")
    axis.annotate("", xy=(0.5, 0.61), xytext=(0.5, 0.68), arrowprops=dict(arrowstyle="-|>", color=GREY))
    axis.annotate("", xy=(0.5, 0.30), xytext=(0.5, 0.37), arrowprops=dict(arrowstyle="-|>", color=GREY))

    axis = figure.add_subplot(grid[0, 1])
    panel_label(axis, "B")
    selected = match.loc[match["discovery_program"].isin(["P3", "P4"])].copy()
    heat = selected.pivot(index="discovery_program", columns="cohort", values="hypergeom_fdr")
    heat = heat[["GSE161529", "GSE246613", "GSE114725"]]
    sns.heatmap(-np.log10(heat), annot=heat.map(lambda p: f"{p:.1e}"), fmt="",
                cmap="cividis", cbar_kws={"label": "−log10(FDR)"}, ax=axis,
                linewidths=0.8, linecolor="white")
    axis.set_xlabel("Independent cohort")
    axis.set_ylabel("Discovery program")
    axis.set_title("Cross-cohort factor matching", fontsize=8)
    axis.tick_params(axis="x", labelrotation=25)
    for tick in axis.get_xticklabels():
        tick.set_ha("right")

    axis = figure.add_subplot(grid[0, 2])
    panel_label(axis, "C")
    p4 = core.loc[core["program"].eq("P4")].sort_values(
        ["datasets_supporting_top50_membership", "gene"], ascending=[True, True]
    )
    colors = [ORANGE if value == 4 else BLUE for value in p4["datasets_supporting_top50_membership"]]
    axis.barh(p4["gene"], p4["datasets_supporting_top50_membership"], color=colors)
    axis.set_xlim(0, 4.3)
    axis.set_xticks([0, 1, 2, 3, 4])
    axis.set_xlabel("Datasets supporting each gene")
    axis.set_title("Locked P4 core", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[1, :2])
    panel_label(axis, "D")
    subtype_order = ["ER+", "HER2+", "TNBC"]
    sns.stripplot(data=exposure, x="subtype", y="P4", order=subtype_order, hue="subtype",
                  palette={"ER+": BLUE, "HER2+": GREEN, "TNBC": ORANGE}, size=6,
                  jitter=0.18, legend=False, ax=axis)
    sns.boxplot(data=exposure, x="subtype", y="P4", order=subtype_order, width=0.45,
                showfliers=False, boxprops={"facecolor": "none"}, ax=axis)
    axis.set_xlabel("")
    axis.set_ylabel("Patient-balanced P4 exposure")
    axis.set_title("P4 is continuous and patient-variable", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[1, 2])
    panel_label(axis, "E")
    wide = paired.pivot(index="patient", columns="tissue", values="P4_primary").dropna()
    for patient, row in wide.iterrows():
        color = ORANGE if row["TUMOR"] > row["NORMAL"] else BLUE
        axis.plot([0, 1], [row["NORMAL"], row["TUMOR"]], color=color, alpha=0.85, marker="o")
        axis.text(1.03, row["TUMOR"], patient, fontsize=7, va="center")
    axis.set_xticks([0, 1], ["Normal", "Tumor"])
    axis.set_xlim(-0.15, 1.3)
    axis.set_ylabel("P4 score")
    axis.set_title("Tumor–normal boundary\n2/4 increased; exact P = 0.875", fontsize=8)
    axis.spines[["top", "right"]].set_visible(False)

    save(figure, "Figure1_cross_cohort_P4_discovery")


def prepare_figshare_p8() -> pd.DataFrame:
    scores = pd.read_csv(FIGSHARE / "P8_scores.tsv", sep="\t")
    x = sm.add_constant(pd.DataFrame({
        "macrophage_identity": scores["macrophage_identity"],
        "hypoxia": scores["hypoxia"],
        "log1p_umi": np.log1p(scores["nCount_RNA"]),
    }), has_constant="add")
    scores["P4_residual"] = sm.OLS(scores["P4_primary"], x).fit().resid
    scores["fibrovascular"] = (
        (scores["fibroblast_ecm"] - scores["fibroblast_ecm"].mean()) / scores["fibroblast_ecm"].std()
        + (scores["endothelial_angiogenic"] - scores["endothelial_angiogenic"].mean())
        / scores["endothelial_angiogenic"].std()
    ) / 2
    return scores


def build_figure2() -> None:
    wu_spots = pd.read_csv(METADATA / "Wu_Visium_frozen_program_spot_scores.tsv.gz", sep="\t")
    wu_spots = wu_spots.loc[wu_spots["patientid"].eq("1160920F")]
    p8 = prepare_figshare_p8()
    wu_effect = pd.read_csv(TABLES / "Wu_Visium_spatial_neighbor_effects_by_patient.tsv", sep="\t")
    fs_effect = pd.read_csv(TABLES / "Figshare_spatial_neighbor_effects_by_section.tsv", sep="\t")

    figure = plt.figure(figsize=(7.20, 6.35), constrained_layout=True)
    grid = figure.add_gridspec(3, 4, height_ratios=[1, 1, 0.9])

    for column_index, (column, title) in enumerate([
        ("P4_residual", "P4 residual"), ("fibrovascular", "Fibrovascular score")
    ]):
        axis = figure.add_subplot(grid[0, column_index])
        if column_index == 0:
            panel_label(axis, "A")
        plotted = axis.scatter(wu_spots["array_col"], -wu_spots["array_row"],
                               c=wu_spots[column], cmap=DIVERGING, s=3.5, linewidths=0)
        axis.set_aspect("equal")
        axis.axis("off")
        axis.set_title(f"Wu 1160920F\n{title}", fontsize=9)
        figure.colorbar(plotted, ax=axis, fraction=0.04, pad=0.01)

    axis = figure.add_subplot(grid[0, 2:])
    panel_label(axis, "B")
    primary = wu_effect.loc[wu_effect["target"].eq("fibrovascular")].copy()
    primary["excess_over_null"] = primary["observed_neighbor_effect"] - primary["null_median"]
    sns.stripplot(data=primary, x="excess_over_null", y="patient", hue="subtype",
                  palette={"TNBC": ORANGE, "ER": BLUE}, size=7, ax=axis)
    axis.axvline(0, color="black", linewidth=0.8, linestyle="--")
    axis.set_xlabel("Observed effect − null median")
    axis.set_ylabel("")
    axis.set_title("Discovery patients: 5/6 positive\nCohort permutation P = 0.0005", fontsize=10)
    axis.legend(title="Subtype", frameon=False, loc="lower right")
    axis.spines[["top", "right"]].set_visible(False)

    for column_index, (column, title) in enumerate([
        ("P4_residual", "P4 residual"), ("fibrovascular", "Fibrovascular score")
    ]):
        axis = figure.add_subplot(grid[1, column_index])
        if column_index == 0:
            panel_label(axis, "C")
        plotted = axis.scatter(p8["y"], p8["x"], c=p8[column], cmap=DIVERGING, s=3.2, linewidths=0)
        axis.invert_yaxis()
        axis.set_aspect("equal")
        axis.axis("off")
        axis.set_title(f"Validation P8\n{title}", fontsize=9)
        figure.colorbar(plotted, ax=axis, fraction=0.04, pad=0.01)

    axis = figure.add_subplot(grid[1, 2:])
    panel_label(axis, "D")
    primary = fs_effect.loc[fs_effect["target"].eq("fibrovascular")].copy()
    primary["excess_over_null"] = primary["observed_neighbor_effect"] - primary["null_median"]
    sns.stripplot(data=primary, x="excess_over_null", y="section", color=PURPLE, size=7, ax=axis)
    axis.axvline(0, color="black", linewidth=0.8, linestyle="--")
    axis.set_xlabel("Observed effect − null median")
    axis.set_ylabel("")
    axis.set_title("Independent validation: 7/8 positive\nCohort permutation P = 0.0005", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[2, :])
    panel_label(axis, "E")
    targets = ["fibroblast_ecm", "endothelial_angiogenic", "malignant_epithelial", "cd8_cytotoxic", "treg"]
    labels = ["Fibroblast/ECM", "Endothelial/angiogenic", "Malignant epithelial", "CD8/cytotoxic", "Treg"]
    combined = []
    for cohort, table, unit in [("Discovery", wu_effect, "patient"), ("Validation", fs_effect, "section")]:
        selected = table.loc[table["target"].isin(targets)].copy()
        selected["excess_over_null"] = selected["observed_neighbor_effect"] - selected["null_median"]
        selected["cohort"] = cohort
        selected["target"] = pd.Categorical(selected["target"], categories=targets, ordered=True)
        combined.append(selected)
    component = pd.concat(combined)
    sns.stripplot(data=component, x="target", y="excess_over_null", hue="cohort",
                  palette={"Discovery": BLUE, "Validation": PURPLE}, dodge=True,
                  alpha=0.68, size=4, ax=axis)
    sns.pointplot(data=component, x="target", y="excess_over_null", hue="cohort",
                  palette={"Discovery": BLUE, "Validation": PURPLE}, dodge=0.35,
                  errorbar=("ci", 95), markers="D", linestyles="none", ax=axis)
    handles, legends = axis.get_legend_handles_labels()
    axis.legend(handles[:2], legends[:2], title="", frameon=False, ncol=2, loc="upper right")
    axis.axhline(0, color="black", linewidth=0.8, linestyle="--")
    axis.set_xticks(range(len(labels)), labels, rotation=15, ha="right")
    axis.set_xlabel("")
    axis.set_ylabel("Neighbor effect − null median")
    axis.set_title("Component specificity: fibrovascular and CD8 effects reproduce", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    save(figure, "Figure2_spatial_fibrovascular_validation")


def build_figure3() -> None:
    clinical = pd.read_csv(TABLES / "GSE246613_locked_program_clinical_join.tsv", sep="\t")
    primary = pd.read_csv(TABLES / "GSE246613_P4_primary_pCR_analysis.tsv", sep="\t").iloc[0]
    specificity = pd.read_csv(TABLES / "GSE246613_longitudinal_module_specificity.tsv", sep="\t")
    survival = pd.read_csv(TABLES / "TNBC94_ecosystem_survival_models.tsv", sep="\t")
    survival = survival.loc[survival["model"].eq("adjusted")]
    loo = pd.read_csv(TABLES / "TNBC94_ecosystem_leave_one_P4_gene_out.tsv", sep="\t")

    figure = plt.figure(figsize=(7.20, 6.25), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, width_ratios=[1, 1.15, 1.1])

    axis = figure.add_subplot(grid[0, 0])
    panel_label(axis, "A")
    baseline = clinical.loc[clinical["treatment"].eq("Base")].copy()
    baseline["P4_z"] = (baseline["P4_primary"] - baseline["P4_primary"].mean()) / baseline["P4_primary"].std()
    baseline["response"] = baseline["pCR"].map({"R": "Responder", "NR": "Nonresponder"})
    order = ["Nonresponder", "Responder"]
    sns.violinplot(data=baseline, x="response", y="P4_z", hue="response", order=order,
                   hue_order=order, palette={"Nonresponder": GREY, "Responder": GREEN},
                   legend=False, inner=None, cut=0, alpha=0.35, ax=axis)
    sns.stripplot(data=baseline, x="response", y="P4_z", order=order,
                  color="black", size=4, jitter=0.16, ax=axis)
    axis.set_xlabel("")
    axis.set_xticks([0, 1], ["Non-pCR", "pCR"])
    axis.set_ylabel("Baseline P4 (z score)")
    axis.set_title("Pathological response", fontsize=8.2, loc="left", x=0.08)
    axis.text(0.02, 0.98,
              f"OR/SD {primary['odds_ratio_per_sd']:.2f} "
              f"({primary['ci95_low']:.2f}–{primary['ci95_high']:.2f}); "
              f"P = {primary['likelihood_ratio_p']:.3f}",
              transform=axis.transAxes, fontsize=6.4, va="top")
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[0, 1])
    panel_label(axis, "B")
    plot = clinical.copy()
    order_map = {"Base": 0, "PD1": 1, "RTPD1": 2}
    plot["time"] = plot["treatment"].map(order_map)
    for patient, values in plot.groupby("cohort", observed=True):
        values = values.sort_values("time")
        response = values["pCR"].dropna().iloc[0]
        axis.plot(values["time"], values["P4_primary"], color=GREEN if response == "R" else GREY,
                  alpha=0.20, linewidth=0.8)
    summary = plot.groupby("time", observed=True)["P4_primary"].agg(["mean", "sem"]).reset_index()
    axis.errorbar(summary["time"], summary["mean"], yerr=summary["sem"], color=ORANGE,
                  marker="o", linewidth=2, capsize=3, label="Mean ± SEM")
    axis.set_xticks([0, 1, 2], ["Baseline", "PD-1", "RT +\nPD-1"])
    axis.set_ylabel("P4 score")
    axis.set_title("Longitudinal remodeling", fontsize=8.2, loc="left", x=0.08)
    axis.text(0.02, 0.98, "RT + PD-1 vs PD-1; FDR = 0.0094",
              transform=axis.transAxes, fontsize=6.4, va="top")
    axis.legend(frameon=False, loc="best")
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[0, 2])
    panel_label(axis, "C")
    module_order = ["P4_primary", "P3_secondary", "P1_inflammatory", "P6_interferon",
                    "hypoxia", "macrophage_identity", "phagolysosome"]
    contrast_order = ["PD1_minus_Base", "RTPD1_minus_PD1", "RTPD1_minus_Base"]
    heat = specificity.pivot(index="module", columns="contrast", values="median_change_in_baseline_sd")
    heat = heat.reindex(index=module_order, columns=contrast_order)
    fdr = specificity.pivot(index="module", columns="contrast", values="fdr_across_modules_and_contrasts")
    fdr = fdr.reindex(index=module_order, columns=contrast_order)
    annotations = heat.copy().astype(str)
    for row in heat.index:
        for col in heat.columns:
            star = "*" if fdr.loc[row, col] < 0.05 else ""
            annotations.loc[row, col] = f"{heat.loc[row, col]:.2f}{star}"
    sns.heatmap(heat, annot=annotations, fmt="", cmap=DIVERGING, center=0,
                cbar_kws={"label": "Median change (baseline SD)"}, linewidths=0.5, ax=axis)
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.set_xticklabels(["PD-1 − base", "RT+PD-1 − PD-1", "RT+PD-1 − base"], rotation=35, ha="right")
    axis.set_title("Module specificity\n*FDR < 0.05 across 21 tests",
                   fontsize=7.3, loc="left", x=0.08)

    axis = figure.add_subplot(grid[1, 0])
    panel_label(axis, "D")
    endpoint_labels = {"DRFS": "Distant relapse-free", "iDFS": "Invasive disease-free", "OS": "Overall survival"}
    survival = survival.assign(label=survival["endpoint"].map(endpoint_labels)).iloc[::-1]
    y = np.arange(len(survival))
    axis.errorbar(survival["hazard_ratio"], y,
                  xerr=[survival["hazard_ratio"] - survival["ci_low"], survival["ci_high"] - survival["hazard_ratio"]],
                  fmt="o", color=BLUE, capsize=3)
    axis.axvline(1, color="black", linestyle="--", linewidth=0.8)
    axis.set_yticks(y, survival["label"])
    axis.set_xlabel("Adjusted hazard ratio per SD")
    axis.set_title("TNBC94 survival tests", fontsize=8.2, loc="left", x=0.08)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[1, 1])
    panel_label(axis, "E")
    loo = loo.copy()
    loo["gene"] = loo["exposure"].str.replace("ecosystem_without_", "", regex=False)
    loo = loo.sort_values("hazard_ratio")
    axis.errorbar(loo["hazard_ratio"], np.arange(len(loo)),
                  xerr=[loo["hazard_ratio"] - loo["ci_low"], loo["ci_high"] - loo["hazard_ratio"]],
                  fmt="o", color=PURPLE, markersize=4, alpha=0.9)
    axis.axvline(1, color="black", linestyle="--", linewidth=0.8)
    axis.set_yticks(np.arange(len(loo)), loo["gene"], fontsize=7)
    axis.set_xlabel("DRFS hazard ratio")
    axis.set_title("Leave-one-gene-out", fontsize=8.2, loc="left", x=0.08)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[1, 2])
    panel_label(axis, "F")
    axis.axis("off")
    evidence = [
        ("Cross-cohort program\nrecurrence", True),
        ("Fibrovascular spatial\nreplication", True),
        ("Local CD8 exclusion", False),
        ("Pathological-response\nprediction", False),
        ("Adverse DRFS association", False),
    ]
    for index, (label, supported) in enumerate(evidence):
        y_pos = 0.88 - index * 0.17
        axis.add_patch(mpl.patches.FancyBboxPatch(
            (0.03, y_pos - 0.065), 0.92, 0.12, boxstyle="round,pad=0.01",
            facecolor="#ECF8F2" if supported else "#F5F5F5",
            edgecolor=GREEN if supported else "#A3A3A3", linewidth=1,
        ))
        axis.text(0.09, y_pos, "+" if supported else "−", color=GREEN if supported else GREY,
                  fontsize=14, fontweight="bold", va="center")
        axis.text(0.19, y_pos, label, fontsize=5.5, va="center", linespacing=0.95)
    axis.set_title("Evidence summary", fontsize=8.2, loc="left", x=0.08)

    save(figure, "Figure3_clinical_boundaries")


def build_figure4() -> None:
    effects = pd.read_csv(TABLES / "TNBC94_spatial_neighbor_effects_by_patient.tsv", sep="\t")
    primary = effects.loc[effects["target"].eq("fibrovascular")].copy()
    primary["patient"] = primary["patient"].astype(str)
    scores = pd.read_csv(TABLES / "TNBC94_frozen_ecosystem_scores.tsv", sep="\t", dtype={"patient": str})
    clinical = pd.read_csv(TABLES / "TNBC94_clinical_ecology_key.tsv", sep="\t", dtype={"patient": str})
    joined = primary.merge(scores[["patient", "ecosystem_score"]], on="patient").merge(
        clinical[["patient", "spatial_archetype", "TIME_pseudobulk"]], on="patient"
    )
    median_effect = primary["excess_over_null"].median()
    representative_patient = primary.iloc[
        (primary["excess_over_null"] - median_effect).abs().argmin()
    ]["patient"]
    spots = pd.read_csv(
        METADATA / "TNBC94_annotated_frozen_program_spot_scores.tsv.gz", sep="\t",
        dtype={"patient": str},
    )
    representative = spots.loc[spots["patient"].eq(representative_patient)].copy()

    archetype_groups = [
        group["excess_over_null"].to_numpy()
        for _, group in joined.dropna(subset=["spatial_archetype"]).groupby("spatial_archetype")
    ]
    kw = kruskal(*archetype_groups)
    rho, rho_p = spearmanr(joined["ecosystem_score"], joined["excess_over_null"], nan_policy="omit")
    pd.DataFrame([
        {"analysis": "spatial_effect_by_published_archetype", "statistic": kw.statistic, "p_value": kw.pvalue},
        {"analysis": "spatial_effect_vs_pseudobulk_ecosystem", "statistic": rho, "p_value": rho_p},
    ]).to_csv(TABLES / "TNBC94_spatial_ecology_exploratory_tests.tsv", sep="\t", index=False)

    figure = plt.figure(figsize=(7.20, 7.15), constrained_layout=True)
    grid = figure.add_gridspec(3, 3, height_ratios=[1, 0.95, 1.05])
    for column_index, (column, title, cmap) in enumerate([
        ("P4_residual", "P4 residual", DIVERGING),
        ("fibrovascular", "Fibrovascular score", DIVERGING),
        ("Classification", "Dominant pathology", "tab20"),
    ]):
        axis = figure.add_subplot(grid[0, column_index])
        if column_index == 0:
            panel_label(axis, "A")
        if column == "Classification":
            categories = pd.Categorical(representative[column])
            color_values = categories.codes
            plotted = axis.scatter(representative["y"], -representative["x"], c=color_values,
                                   cmap=cmap, s=8, linewidths=0)
            top_labels = representative[column].value_counts().head(5).index
            legend_handles = []
            for label in top_labels:
                code = list(categories.categories).index(label)
                legend_handles.append(mpl.lines.Line2D([], [], marker="o", linestyle="",
                    color=mpl.colormaps[cmap](code / max(len(categories.categories) - 1, 1)), label=label))
            axis.legend(handles=legend_handles, fontsize=6, frameon=False, loc="upper left",
                        bbox_to_anchor=(1.0, 1.0))
        else:
            plotted = axis.scatter(representative["y"], -representative["x"], c=representative[column],
                                   cmap=cmap, s=8, linewidths=0)
            figure.colorbar(plotted, ax=axis, fraction=0.04, pad=0.01)
        axis.set_aspect("equal")
        axis.axis("off")
        axis.set_title(f"Median-effect patient {representative_patient}\n{title}", fontsize=9)

    axis = figure.add_subplot(grid[1, :2])
    panel_label(axis, "B")
    ordered = primary.sort_values("excess_over_null").reset_index(drop=True)
    colors = np.where(ordered["excess_over_null"] > 0, GREEN, GREY)
    axis.bar(np.arange(len(ordered)), ordered["excess_over_null"], color=colors, width=0.85)
    axis.axhline(0, color="black", linewidth=0.8)
    axis.set_xlabel("Patients ordered by excess effect")
    axis.set_ylabel("Fibrovascular effect − null median")
    axis.set_title(r"82/94 patients positive; exact sign-test $P = 5.62 \times 10^{-14}$", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[1, 2])
    panel_label(axis, "C")
    sns.boxplot(data=joined, x="spatial_archetype", y="excess_over_null", color="#DCE6F1",
                showfliers=False, ax=axis)
    sns.stripplot(data=joined, x="spatial_archetype", y="excess_over_null", color=BLUE,
                  size=3, alpha=0.65, jitter=0.18, ax=axis)
    axis.axhline(0, color="black", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Published TNBC spatial archetype")
    axis.set_ylabel("Excess effect")
    axis.set_title("Published spatial archetypes", fontsize=8.2, loc="left", x=0.10)
    axis.text(0.02, 0.98, f"Kruskal–Wallis P = {kw.pvalue:.3f}",
              transform=axis.transAxes, fontsize=6.4, va="top")
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[2, :2])
    panel_label(axis, "D")
    component = effects.copy()
    target_order = ["fibroblast_ecm", "endothelial_angiogenic", "malignant_epithelial", "cd8_cytotoxic", "treg"]
    component = component.loc[component["target"].isin(target_order)]
    component["target"] = pd.Categorical(component["target"], target_order, ordered=True)
    sns.violinplot(data=component, x="target", y="excess_over_null", color="#E4EAF0",
                   inner=None, cut=0, ax=axis)
    sns.boxplot(data=component, x="target", y="excess_over_null", width=0.22, showfliers=False,
                boxprops={"facecolor": "white"}, ax=axis)
    axis.axhline(0, color="black", linestyle="--", linewidth=0.8)
    axis.set_xticks(
        range(5),
        ["Fibroblast/ECM", "Endothelial", "Malignant", "CD8/cytotoxic", "Treg"],
        rotation=15,
        ha="right",
    )
    axis.set_xlabel("")
    axis.set_ylabel("Neighbor effect − null median")
    axis.set_title("Component effects across 94 patients", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)

    axis = figure.add_subplot(grid[2, 2])
    panel_label(axis, "E")
    sns.regplot(data=joined, x="ecosystem_score", y="excess_over_null", scatter_kws={"s": 20, "alpha": 0.7},
                line_kws={"color": ORANGE}, color=BLUE, ax=axis)
    axis.axhline(0, color="black", linestyle="--", linewidth=0.8)
    axis.set_xlabel("Pseudobulk ecosystem score")
    axis.set_ylabel("Spatial excess effect")
    axis.set_title("Bulk–spatial relationship", fontsize=8.2, loc="left", x=0.10)
    axis.text(0.02, 0.98, f"Spearman ρ = {rho:.2f}; P = {rho_p:.3f}",
              transform=axis.transAxes, fontsize=6.4, va="top")
    axis.spines[["top", "right"]].set_visible(False)

    save(figure, "Figure4_TNBC94_spatial_extension")


def main() -> None:
    sns.set_theme(style="ticks", context="paper", font_scale=1.0)
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.titlesize": 8.2,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 6.7,
        "ytick.labelsize": 6.7,
        "legend.fontsize": 6.7,
        "legend.title_fontsize": 6.7,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "savefig.transparent": False,
    })
    build_figure1()
    build_figure2()
    build_figure3()
    build_figure4()
    print(f"Wrote main figures to {FIGURES}")


if __name__ == "__main__":
    main()
