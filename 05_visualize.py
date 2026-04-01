#!/usr/bin/env python3
"""
05_visualize.py — Generate charts, 3D binding poses, and report.

Outputs:
  - results/charts/*.png — score distribution, top 20, selectivity, dual-activity, CNS
  - results/poses_3d/*.html — interactive 3D binding poses (py3Dmol)
  - results/report.md — full report with embedded charts and links
"""

import textwrap
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import py3Dmol
import seaborn as sns

RESULTS_DIR = Path("results")
CHARTS_DIR = RESULTS_DIR / "charts"
POSES_DIR = RESULTS_DIR / "poses_3d"
DATA_DIR = Path("data")

sns.set_theme(style="whitegrid", font_scale=1.1)
plt.rcParams["figure.dpi"] = 150


def load_data():
    """Load all result CSVs."""
    ranked = pd.read_csv(RESULTS_DIR / "repurposing_ranked.csv")
    top20 = pd.read_csv(RESULTS_DIR / "top20_hits.csv")
    controls = pd.read_csv(RESULTS_DIR / "controls.csv")
    benchmark_metrics_path = RESULTS_DIR / "benchmark_metrics.csv"
    benchmark_metrics = (
        pd.read_csv(benchmark_metrics_path) if benchmark_metrics_path.exists() else pd.DataFrame()
    )

    c1_path = RESULTS_DIR / "campaign1_scores.csv"
    c1 = pd.read_csv(c1_path) if c1_path.exists() else pd.DataFrame()

    return ranked, top20, controls, c1, benchmark_metrics


def chart_score_distribution(c1, controls):
    """Histogram of all drug scores with control TKI lines."""
    fig, ax = plt.subplots(figsize=(10, 6))

    scores = pd.to_numeric(c1["score"], errors="coerce").dropna()
    ax.hist(scores, bins=60, color="#4C72B0", alpha=0.7, edgecolor="white")

    # Vertical lines for each control TKI
    colors = plt.cm.Set1(np.linspace(0, 1, len(controls)))
    for i, (_, row) in enumerate(controls.iterrows()):
        score = row.get("g2032r_score")
        if pd.notna(score):
            label = row["drug_name"]
            if len(label) > 20:
                label = label[:18] + "..."
            ax.axvline(score, color=colors[i], linestyle="--", linewidth=2, label=f"{label} ({score:.1f})")

    ax.set_xlabel("Binding Score (kcal/mol)")
    ax.set_ylabel("Count")
    ax.set_title("Drug Score Distribution vs Known ROS1 TKIs")
    ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "score_distribution.png")
    plt.close()
    print("  Saved score_distribution.png")


def chart_top20_bars(top20, controls):
    """Horizontal bar chart of top 20 hits + controls."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Combine top 20 + controls
    ctrl_df = controls[["drug_name", "g2032r_score", "cns_mpo"]].copy()
    ctrl_df["type"] = "Known TKI"
    hit_df = top20[["drug_name", "g2032r_score", "cns_mpo"]].copy()
    hit_df["type"] = "Hit"

    combined = pd.concat([hit_df, ctrl_df], ignore_index=True)
    combined = combined.sort_values("g2032r_score", ascending=True)

    # Color by CNS MPO
    colors = []
    for _, row in combined.iterrows():
        cns = row.get("cns_mpo", 0) or 0
        if row["type"] == "Known TKI":
            colors.append("#888888")  # gray for controls
        elif cns >= 4.0:
            colors.append("#2ca02c")  # green
        elif cns >= 3.0:
            colors.append("#ff7f0e")  # orange
        else:
            colors.append("#d62728")  # red

    # Truncate long names
    names = [n[:30] + "..." if len(n) > 30 else n for n in combined["drug_name"]]

    bars = ax.barh(range(len(combined)), combined["g2032r_score"], color=colors, edgecolor="white")
    ax.set_yticks(range(len(combined)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("G2032R Binding Score (kcal/mol)")
    ax.set_title("Top 20 Hits vs Known ROS1 TKIs")

    # Legend
    patches = [
        mpatches.Patch(color="#2ca02c", label="CNS MPO >= 4.0"),
        mpatches.Patch(color="#ff7f0e", label="CNS MPO 3.0-4.0"),
        mpatches.Patch(color="#d62728", label="CNS MPO < 3.0"),
        mpatches.Patch(color="#888888", label="Known TKI"),
    ]
    ax.legend(handles=patches, loc="lower right")
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "top20_scores.png")
    plt.close()
    print("  Saved top20_scores.png")


def chart_selectivity(top20, controls):
    """Scatter: G2032R vs WT scores."""
    # Need campaign 3 data
    c3_path = RESULTS_DIR / "campaign3_scores.csv"
    if not c3_path.exists():
        print("  Skipping selectivity chart (no campaign 3 data)")
        return

    c3 = pd.read_csv(c3_path)
    if "wt_score" not in c3.columns:
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot controls
    ctrl_in_c3 = c3[c3["drug_name"].isin(controls["drug_name"])]
    hits_in_c3 = c3[~c3["drug_name"].isin(controls["drug_name"])]

    if len(ctrl_in_c3) > 0:
        ax.scatter(ctrl_in_c3["wt_score"], ctrl_in_c3["g2032r_score"],
                   c="#888888", s=80, marker="D", label="Known TKI", zorder=3)
        for _, row in ctrl_in_c3.iterrows():
            if pd.notna(row["wt_score"]) and pd.notna(row["g2032r_score"]):
                name = row["drug_name"][:15]
                ax.annotate(name, (row["wt_score"], row["g2032r_score"]),
                            fontsize=7, ha="left", va="bottom")

    if len(hits_in_c3) > 0:
        ax.scatter(hits_in_c3["wt_score"], hits_in_c3["g2032r_score"],
                   c="#d62728", s=80, marker="o", label="Hit", zorder=3)
        for _, row in hits_in_c3.iterrows():
            if pd.notna(row["wt_score"]) and pd.notna(row["g2032r_score"]):
                name = row["drug_name"][:15]
                ax.annotate(name, (row["wt_score"], row["g2032r_score"]),
                            fontsize=7, ha="left", va="bottom")

    # Diagonal line (equal selectivity)
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "k--", alpha=0.3, label="Equal selectivity")
    ax.fill_between(lims, lims, [lims[0]] * 2, alpha=0.05, color="green")

    ax.set_xlabel("WT ROS1 Score (kcal/mol)")
    ax.set_ylabel("G2032R Score (kcal/mol)")
    ax.set_title("Mutant vs Wild-type Selectivity\n(below diagonal = prefers mutant)")
    ax.legend()
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "selectivity_g2032r_vs_wt.png")
    plt.close()
    print("  Saved selectivity_g2032r_vs_wt.png")


def chart_dual_activity(top20, controls):
    """Scatter: G2032R vs MET scores."""
    c3_path = RESULTS_DIR / "campaign3_scores.csv"
    if not c3_path.exists():
        print("  Skipping dual-activity chart (no campaign 3 data)")
        return

    c3 = pd.read_csv(c3_path)
    if "met_score" not in c3.columns:
        return

    fig, ax = plt.subplots(figsize=(8, 8))

    ctrl_in_c3 = c3[c3["drug_name"].isin(controls["drug_name"])]
    hits_in_c3 = c3[~c3["drug_name"].isin(controls["drug_name"])]

    if len(ctrl_in_c3) > 0:
        ax.scatter(ctrl_in_c3["met_score"], ctrl_in_c3["g2032r_score"],
                   c="#888888", s=80, marker="D", label="Known TKI", zorder=3)
        for _, row in ctrl_in_c3.iterrows():
            if pd.notna(row["met_score"]) and pd.notna(row["g2032r_score"]):
                ax.annotate(row["drug_name"][:15],
                            (row["met_score"], row["g2032r_score"]),
                            fontsize=7, ha="left", va="bottom")

    if len(hits_in_c3) > 0:
        ax.scatter(hits_in_c3["met_score"], hits_in_c3["g2032r_score"],
                   c="#d62728", s=80, marker="o", label="Hit", zorder=3)
        for _, row in hits_in_c3.iterrows():
            if pd.notna(row["met_score"]) and pd.notna(row["g2032r_score"]):
                ax.annotate(row["drug_name"][:15],
                            (row["met_score"], row["g2032r_score"]),
                            fontsize=7, ha="left", va="bottom")

    # Highlight sweet spot quadrant
    ax.axhline(-8.5, color="green", linestyle=":", alpha=0.4)
    ax.axvline(-8.5, color="green", linestyle=":", alpha=0.4)
    ax.text(ax.get_xlim()[0] + 0.2, -8.7, "Strong ROS1 + Strong MET",
            fontsize=8, color="green", alpha=0.7)

    ax.set_xlabel("MET Score (kcal/mol)")
    ax.set_ylabel("G2032R Score (kcal/mol)")
    ax.set_title("Dual ROS1/MET Activity")
    ax.legend()
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "dual_ros1_met.png")
    plt.close()
    print("  Saved dual_ros1_met.png")


def chart_cns_vs_binding(ranked, top20, controls):
    """CNS MPO vs binding score scatter for all drugs."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # All drugs as light dots
    all_drugs = ranked[ranked["g2032r_score"].notna()].copy()
    ax.scatter(all_drugs["g2032r_score"], all_drugs["cns_mpo"],
               c="#cccccc", s=10, alpha=0.5, label="All drugs")

    # Top 20 highlighted
    if len(top20) > 0 and "g2032r_score" in top20.columns:
        ax.scatter(top20["g2032r_score"], top20["cns_mpo"],
                   c="#d62728", s=60, zorder=3, label="Top 20 hits")
        for _, row in top20.iterrows():
            if pd.notna(row["g2032r_score"]):
                name = row["drug_name"][:15]
                ax.annotate(name, (row["g2032r_score"], row["cns_mpo"]),
                            fontsize=6, ha="left", va="bottom")

    # Controls
    if len(controls) > 0 and "g2032r_score" in controls.columns:
        ctrl_with_cns = controls[controls["cns_mpo"].notna()]
        if len(ctrl_with_cns) > 0:
            ax.scatter(ctrl_with_cns["g2032r_score"], ctrl_with_cns["cns_mpo"],
                       c="#888888", s=60, marker="D", zorder=3, label="Known TKI")

    # Quadrant lines
    ax.axhline(4.0, color="green", linestyle="--", alpha=0.4)
    ax.axvline(-8.5, color="green", linestyle="--", alpha=0.4)

    # Highlight sweet spot
    xlim = ax.get_xlim()
    ax.fill_between([xlim[0], -8.5], 4.0, 6.5, alpha=0.08, color="green")
    ax.text(xlim[0] + 0.2, 6.2, "Sweet spot:\nStrong binding + CNS-penetrant",
            fontsize=9, color="green")

    ax.set_xlabel("G2032R Binding Score (kcal/mol)")
    ax.set_ylabel("CNS MPO Score")
    ax.set_title("CNS Penetration vs Binding Affinity")
    ax.legend(loc="lower left")
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "cns_vs_binding.png")
    plt.close()
    print("  Saved cns_vs_binding.png")


def generate_3d_pose(drug_name, receptor_pdb, pose_pdbqt, output_html,
                     key_residues=None):
    """Generate interactive 3D view with py3Dmol."""
    if key_residues is None:
        key_residues = [1980, 2027, 2102, 2032]

    try:
        receptor_str = Path(receptor_pdb).read_text()
        pose_str = Path(pose_pdbqt).read_text()
    except FileNotFoundError:
        return False

    # Extract first pose from PDBQT (up to first ENDMDL)
    pose_lines = []
    for line in pose_str.split("\n"):
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            pose_lines.append(line)
    first_pose = "\n".join(pose_lines)

    view = py3Dmol.view(width=800, height=600)

    # Protein: cartoon (light gray) with surface near binding site
    view.addModel(receptor_str, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "#dddddd", "opacity": 0.8}})

    # Key residues as sticks with labels
    for resnum in key_residues:
        view.addStyle(
            {"model": 0, "resi": resnum},
            {"stick": {"colorscheme": "greenCarbon", "radius": 0.15}},
        )
        view.addResLabels(
            {"model": 0, "resi": resnum},
            {"fontSize": 10, "backgroundColor": "white", "backgroundOpacity": 0.7},
        )

    # Ligand as colored sticks
    view.addModel(first_pose, "pdb")
    view.setStyle(
        {"model": 1},
        {"stick": {"colorscheme": "cyanCarbon", "radius": 0.2}},
    )

    view.zoomTo({"model": 1})
    view.zoom(0.7)

    # Write HTML
    html = f"""<!DOCTYPE html>
<html>
<head><title>{drug_name} - ROS1 G2032R Binding Pose</title></head>
<body style="margin:0; font-family:sans-serif;">
<h3 style="padding:10px; margin:0;">{drug_name} docked into ROS1 G2032R</h3>
<div style="width:800px; height:600px;">
{view._make_html()}
</div>
<p style="padding:10px; color:#666;">
Green sticks: K1980, E2027, D2102, R2032 (catalytic triad + mutation site).
Cyan sticks: docked ligand.
</p>
</body>
</html>"""

    Path(output_html).write_text(html)
    return True


def generate_3d_poses(top20, controls):
    """Generate 3D binding poses for the current top 20 hits and reference."""
    POSES_DIR.mkdir(parents=True, exist_ok=True)

    receptor_pdb = str(DATA_DIR / "receptor_G2032R.pdb")
    expected_outputs = set()

    # Top 20 hits — look in campaign2 poses first, then campaign1
    generated = 0
    for _, row in top20.iterrows():
        name = row["drug_name"]
        safe_name = name.replace("/", "_").replace(" ", "_")

        # Find pose file
        pose_path = None
        for campaign in ["campaign2", "campaign1"]:
            candidate = RESULTS_DIR / "poses" / campaign / f"{safe_name}.pdbqt"
            if candidate.exists():
                pose_path = str(candidate)
                break

        if pose_path is None:
            continue

        output = str(POSES_DIR / f"{safe_name}.html")
        if generate_3d_pose(name, receptor_pdb, pose_path, output):
            generated += 1
            expected_outputs.add(Path(output))

    # Reference: zidesamtinib
    zid_pose = RESULTS_DIR / "poses" / "campaign2" / "zidesamtinib.pdbqt"
    if not zid_pose.exists():
        zid_pose = RESULTS_DIR / "poses" / "campaign1" / "zidesamtinib.pdbqt"
    if zid_pose.exists():
        generate_3d_pose(
            "zidesamtinib (reference)",
            receptor_pdb,
            str(zid_pose),
            str(POSES_DIR / "zidesamtinib_reference.html"),
        )
        generated += 1
        expected_outputs.add(POSES_DIR / "zidesamtinib_reference.html")

    if expected_outputs:
        for html_file in POSES_DIR.glob("*.html"):
            if html_file not in expected_outputs:
                html_file.unlink()

    print(f"  Generated {generated} interactive 3D pose views")
    return sorted(expected_outputs)


def generate_report(top20, controls, ranked, pose_files, benchmark_metrics):
    """Generate results/report.md."""
    n_scored = ranked["g2032r_score"].notna().sum()
    n_candidates = ranked["is_candidate"].sum() if "is_candidate" in ranked.columns else "N/A"
    best_tki = controls["g2032r_score"].min() if len(controls) > 0 else "N/A"

    # Validation RMSD
    rmsd_path = RESULTS_DIR / "validation_rmsd.txt"
    rmsd_text = rmsd_path.read_text().strip() if rmsd_path.exists() else "Not available"

    # Controls table
    ctrl_md = controls.to_markdown(index=False) if len(controls) > 0 else "No controls found"

    # Top 20 table
    display_cols = ["drug_name", "g2032r_score", "composite_score", "cns_mpo",
                    "cns_penetrant", "mw", "clogp", "target"]
    available = [c for c in display_cols if c in top20.columns]
    top20_md = top20[available].to_markdown(index=False) if len(top20) > 0 else "No hits found"

    # 3D pose links
    pose_links = ""
    for html_file in pose_files:
        pose_links += f"- [{html_file.stem}](poses_3d/{html_file.name})\n"

    has_admet = {"lipinski_pass", "bbb_permeant", "herg_liability"}.issubset(top20.columns)
    multiconf_path = RESULTS_DIR / "campaign2_scores_single_conf.csv"
    multiconf_pose_dir = RESULTS_DIR / "poses" / "campaign2_multiconf"
    has_multiconf = multiconf_path.exists() or (
        multiconf_pose_dir.exists() and any(multiconf_pose_dir.glob("*.pdbqt"))
    )
    validation_status = "Not available"
    if rmsd_path.exists():
        if "Aligned RMSD" in rmsd_text:
            validation_status = "Alignment-aware RMSD analysis completed in 06_improve.py."
        elif "deferred to 06_improve.py" in rmsd_text:
            validation_status = "Only preliminary re-docking score is available; detailed RMSD has not been rerun yet."
        else:
            validation_status = "Validation file present, but status could not be classified automatically."

    improvements = []
    improvements.append(f"- **Validation status** — {validation_status}")
    improvements.append(
        "- **Multi-conformer docking** — "
        + (
            "Campaign 2 was refreshed from multi-conformer docking outputs."
            if has_multiconf
            else "No multi-conformer evidence detected in the current outputs."
        )
    )
    improvements.append(
        "- **ADMET annotation** — "
        + (
            "Top-hit ADMET columns are present in `top20_hits.csv`."
            if has_admet
            else "ADMET columns are not present in the current `top20_hits.csv`."
        )
    )

    benchmark_section = "Benchmark outputs not generated yet."
    benchmark_chart = ""
    if not benchmark_metrics.empty:
        metric_cols = [
            "benchmark_id",
            "benchmark_mode",
            "score_column",
            "roc_auc",
            "average_precision",
            "ef1",
            "nef1",
            "ef5",
            "nef5",
            "ef10",
            "nef10",
            "median_active_rank",
            "hits_top_10",
            "hits_top_25",
            "hits_top_50",
        ]
        available_metric_cols = [c for c in metric_cols if c in benchmark_metrics.columns]
        benchmark_section = benchmark_metrics[available_metric_cols].to_markdown(index=False)
        benchmark_chart_path = CHARTS_DIR / "benchmark_enrichment.png"
        calibration_chart_path = CHARTS_DIR / "benchmark_calibration.png"
        if benchmark_chart_path.exists():
            benchmark_chart = "\n## Benchmark Enrichment\n\n![Benchmark enrichment](charts/benchmark_enrichment.png)\n"
        if calibration_chart_path.exists():
            benchmark_chart += "\n## Benchmark Calibration\n\n![Benchmark calibration](charts/benchmark_calibration.png)\n"

    report = textwrap.dedent(f"""\
# ROS1 Drug Repurposing Screen — Results

**Date:** {date.today().isoformat()}

**Patient context:** 39yo, lung adenocarcinoma.
EZR::ROS1 fusion (exon 34), MET IHC 3+ 50%, PD-L1 TPS 95%.
High risk of CNS metastasis.

## Methods

- **Targets:** ROS1 G2032R (PDB 9QEK), ROS1 WT (PDB 7Z5X), MET (PDB 2WGJ)
- **Grid box:** 25x25x25 A centered on ATP-binding pocket catalytic triad
- **Docking:** AutoDock Vina
  - Campaign 1: All drugs, exhaustiveness=8
  - Campaign 2: Top 50 + controls, exhaustiveness=32, 9 poses
  - Campaign 3: Top 20 + controls vs WT and MET, exhaustiveness=32
- **Library:** {n_scored} FDA-approved drugs (Selleckchem L1300 + L8000)

## Re-docking Validation

```
{rmsd_text}
```

## Known TKI Control Scores

{ctrl_md}

## Top 20 Repurposing Candidates

{top20_md}

## Score Distribution

![Score distribution](charts/score_distribution.png)

## Top 20 vs Controls

![Top 20 bar chart](charts/top20_scores.png)

## CNS Penetration vs Binding

![CNS vs binding](charts/cns_vs_binding.png)

## G2032R vs WT Selectivity

![Selectivity scatter](charts/selectivity_g2032r_vs_wt.png)

## Dual ROS1/MET Activity

![Dual activity](charts/dual_ros1_met.png)

## Benchmark Metrics

{benchmark_section}
{benchmark_chart}

## Interactive 3D Binding Poses

{pose_links if pose_links else "No poses generated yet."}

## Summary

- **Drugs screened:** {n_scored}
- **Best TKI score:** {f"{best_tki:.1f}" if isinstance(best_tki, (int, float)) else best_tki} kcal/mol
- **Repurposing candidates** (within 2 kcal/mol of best TKI): {n_candidates}
- **Composite scoring:** G2032R binding + CNS MPO bonus + MET dual-activity bonus + mutant selectivity bonus

## Limitations

- **Docking != binding.** Vina scores are approximations. Experimental validation is required.
- **Scoring function limitations.** Vina's empirical scoring may miss important interactions
  (e.g., cation-pi, halogen bonds, water-mediated contacts).
- **Static receptor.** No induced fit or protein flexibility modeled.
- **Off-target effects.** Predicted binding does not account for selectivity across the kinome.

## Current Improvement Status

{chr(10).join(improvements)}

## Next Steps

1. **Discuss with oncologist** — review candidates for clinical plausibility
2. **Experimental validation** — CRO testing (e.g., Eurofins, Reaction Biology)
   for top candidates: biochemical kinase assay, cell-based ROS1 activity
3. **Community** — share findings with ROS1ders patient network
4. **Molecular dynamics** — run MD simulations on top 3-5 hits for binding stability
""")

    (RESULTS_DIR / "report.md").write_text(report)
    print("  Saved results/report.md")


def main():
    CHARTS_DIR.mkdir(parents=True, exist_ok=True)
    POSES_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Generating Visualizations ===\n")
    ranked, top20, controls, c1, benchmark_metrics = load_data()

    # Charts
    print("Generating charts...")
    if len(c1) > 0:
        chart_score_distribution(c1, controls)
    chart_top20_bars(top20, controls)
    chart_selectivity(top20, controls)
    chart_dual_activity(top20, controls)
    chart_cns_vs_binding(ranked, top20, controls)

    # 3D poses
    print("\nGenerating 3D binding poses...")
    pose_files = generate_3d_poses(top20, controls)

    # Report
    print("\nGenerating report...")
    generate_report(top20, controls, ranked, pose_files, benchmark_metrics)

    print("\n=== Visualization complete ===")
    print(f"Charts: {CHARTS_DIR}")
    print(f"3D poses: {POSES_DIR}")
    print(f"Report: {RESULTS_DIR / 'report.md'}")


if __name__ == "__main__":
    main()
