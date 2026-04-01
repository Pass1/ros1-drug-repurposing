#!/usr/bin/env python3
"""
07_benchmark.py — Evaluate enrichment and active recovery for the ROS1 screen.

The benchmark layer is schema-driven:
- `data/benchmarks/benchmark_sets.csv` declares benchmark sets and modes
- `data/benchmarks/benchmark_compounds.csv` declares labeled compounds/entities

Two benchmark modes are supported:
- `active_vs_background`: labeled actives vs the rest of the scored library
- `labeled_only`: labeled actives vs labeled inactive/decoy compounds only

This keeps the current provisional ROS1-active recovery benchmark usable now
while making room for stronger literature-curated active/inactive or
active/decoy panels before the next full rerun.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA_DIR = Path("data")
BENCHMARK_DIR = DATA_DIR / "benchmarks"
RESULTS_DIR = Path("results")
CHARTS_DIR = RESULTS_DIR / "charts"

BENCHMARK_SETS_FILE = BENCHMARK_DIR / "benchmark_sets.csv"
BENCHMARK_COMPOUNDS_FILE = BENCHMARK_DIR / "benchmark_compounds.csv"
RANKED_FILE = RESULTS_DIR / "repurposing_ranked.csv"
LIBRARY_FILE = DATA_DIR / "drug_library.csv"

DEFAULT_SCORE_COLUMNS = ("g2032r_score", "composite_score")
TOP_K_VALUES = (10, 25, 50, 100)
EF_FRACTIONS = (0.01, 0.05, 0.10)
CALIBRATION_QUANTILES = 5

SET_REQUIRED_COLUMNS = {
    "benchmark_id",
    "description",
    "benchmark_mode",
    "score_columns",
    "target",
    "mutation",
}
COMPOUND_REQUIRED_COLUMNS = {
    "benchmark_id",
    "entity_id",
    "drug_name",
    "label",
}


def _require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise SystemExit(f"Missing columns in {path}: {', '.join(missing)}")


def _parse_score_columns(value: str) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return list(DEFAULT_SCORE_COLUMNS)
    return [part.strip() for part in str(value).split("|") if part.strip()]


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load benchmark declarations plus ranked results with library metadata."""
    benchmark_sets = pd.read_csv(BENCHMARK_SETS_FILE)
    benchmark_compounds = pd.read_csv(BENCHMARK_COMPOUNDS_FILE)
    ranked = pd.read_csv(RANKED_FILE)
    library = pd.read_csv(LIBRARY_FILE)

    _require_columns(benchmark_sets, SET_REQUIRED_COLUMNS, BENCHMARK_SETS_FILE)
    _require_columns(benchmark_compounds, COMPOUND_REQUIRED_COLUMNS, BENCHMARK_COMPOUNDS_FILE)

    extra_cols = [
        c for c in ["drug_name", "inchikey", "is_known_ros1_tki", "source", "target", "cas"]
        if c in library.columns
    ]
    ranked = ranked.merge(
        library[extra_cols].drop_duplicates(subset=["drug_name"]),
        on="drug_name",
        how="left",
        suffixes=("", "_library"),
    )
    return benchmark_sets, benchmark_compounds, ranked


def _roc_auc_from_scores(y_true: np.ndarray, scores: np.ndarray) -> float:
    """ROC AUC using the rank-sum formula. Higher score means better."""
    n_pos = int(y_true.sum())
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")

    order = np.argsort(scores)
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg_rank = (i + j + 2) / 2.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1

    pos_ranks = ranks[y_true == 1].sum()
    return float((pos_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def _average_precision(y_true: np.ndarray, scores: np.ndarray) -> float:
    """Average precision with higher scores meaning better."""
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return float("nan")

    order = np.argsort(-scores, kind="mergesort")
    y_sorted = y_true[order]
    cum_pos = np.cumsum(y_sorted)
    precision = cum_pos / np.arange(1, len(y_sorted) + 1)
    return float((precision * y_sorted).sum() / n_pos)


def _enrichment_factor(sorted_labels: np.ndarray, fraction: float) -> tuple[float, float, int, int]:
    """Compute EF and normalized EF at a given top-library fraction."""
    n_total = len(sorted_labels)
    n_pos = int(sorted_labels.sum())
    if n_total == 0 or n_pos == 0:
        return float("nan"), float("nan"), 0, 0

    top_k = max(1, math.ceil(n_total * fraction))
    hits = int(sorted_labels[:top_k].sum())
    expected = n_pos * (top_k / n_total)
    ef = hits / expected if expected else float("nan")
    max_hits = min(n_pos, top_k)
    max_ef = max_hits / expected if expected else float("nan")
    nef = ef / max_ef if max_ef and not math.isnan(max_ef) else float("nan")
    return float(ef), float(nef), hits, top_k


def _topk_recovery(sorted_labels: np.ndarray, top_k: int) -> tuple[int, float]:
    n_pos = int(sorted_labels.sum())
    if n_pos == 0:
        return 0, float("nan")
    k = min(top_k, len(sorted_labels))
    hits = int(sorted_labels[:k].sum())
    return hits, float(hits / n_pos)


def _entity_view(membership: pd.DataFrame, score_column: str) -> pd.DataFrame:
    """Collapse rows to one ranked entity per benchmark entity."""
    df = membership[membership[score_column].notna()].copy()
    if df.empty:
        return df

    entity_label_counts = (
        df.groupby("entity_id")["label"]
        .nunique()
        .reset_index(name="n_labels")
    )
    ambiguous = entity_label_counts[entity_label_counts["n_labels"] > 1]
    if not ambiguous.empty:
        ids = ", ".join(ambiguous["entity_id"].astype(str).tolist()[:10])
        raise SystemExit(
            f"Benchmark contains entity_ids with conflicting labels for {score_column}: {ids}"
        )

    aggregated = (
        df.sort_values([score_column, "drug_name"], ascending=[True, True])
        .groupby("entity_id", as_index=False)
        .agg(
            benchmark_id=("benchmark_id", "first"),
            benchmark_mode=("benchmark_mode", "first"),
            benchmark_name=("benchmark_name", "first"),
            label=("label", "first"),
            representative_name=("drug_name", "first"),
            score=(score_column, "min"),
            evidence_tier=("evidence_tier", "first"),
            source_ref=("source_ref", "first"),
            source_kind=("source_kind", "first"),
        )
        .sort_values(["score", "entity_id"], ascending=[True, True])
        .reset_index(drop=True)
    )

    aggregated["is_active"] = aggregated["label"].eq("active")
    aggregated["rank"] = np.arange(1, len(aggregated) + 1)
    aggregated["screen_fraction"] = aggregated["rank"] / len(aggregated)
    aggregated["cumulative_actives"] = aggregated["is_active"].cumsum()
    total_actives = int(aggregated["is_active"].sum())
    aggregated["active_recall"] = (
        aggregated["cumulative_actives"] / total_actives if total_actives else np.nan
    )
    aggregated["score_column"] = score_column
    return aggregated


def build_membership_for_benchmark(
    benchmark_set: pd.Series,
    benchmark_compounds: pd.DataFrame,
    ranked: pd.DataFrame,
) -> pd.DataFrame:
    """Build row-level membership table for one benchmark definition."""
    benchmark_id = benchmark_set["benchmark_id"]
    benchmark_mode = benchmark_set["benchmark_mode"]
    benchmark_name = benchmark_set.get("name", benchmark_id)

    compounds = benchmark_compounds[benchmark_compounds["benchmark_id"] == benchmark_id].copy()
    if compounds.empty:
        return pd.DataFrame()

    compounds["label"] = compounds["label"].astype(str).str.lower().str.strip()
    valid_labels = {"active", "inactive", "decoy"}
    invalid_labels = sorted(set(compounds["label"]) - valid_labels)
    if invalid_labels:
        raise SystemExit(f"{benchmark_id}: invalid labels: {', '.join(invalid_labels)}")

    labeled = ranked.merge(compounds, on="drug_name", how="left")
    labeled["benchmark_id"] = benchmark_id
    labeled["benchmark_name"] = benchmark_name
    labeled["benchmark_mode"] = benchmark_mode
    labeled["target_benchmark"] = benchmark_set.get("target", "")
    labeled["mutation_benchmark"] = benchmark_set.get("mutation", "")

    if benchmark_mode == "active_vs_background":
        labeled["label"] = labeled["label"].fillna("background")
        labeled["entity_id"] = labeled["entity_id"].fillna(labeled["drug_name"])
        labeled["notes"] = labeled["notes"].fillna("Unlabeled scored row treated as presumed background.")
        labeled["evidence_tier"] = labeled["evidence_tier"].fillna("background")
        labeled["source_kind"] = labeled["source_kind"].fillna("screen_background")
        labeled["source_ref"] = labeled["source_ref"].fillna("")
        return labeled

    if benchmark_mode == "labeled_only":
        labeled_only = labeled[labeled["label"].isin(["active", "inactive", "decoy"])].copy()
        return labeled_only

    raise SystemExit(f"{benchmark_id}: unsupported benchmark_mode `{benchmark_mode}`")


def compute_metrics_for_membership(
    membership: pd.DataFrame,
    benchmark_set: pd.Series,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], list[pd.DataFrame]]:
    """Compute metrics, entity tables, and calibration tables for one benchmark."""
    metrics_rows = []
    entity_tables: dict[str, pd.DataFrame] = {}
    calibration_tables: list[pd.DataFrame] = []

    score_columns = [
        c for c in _parse_score_columns(benchmark_set["score_columns"])
        if c in membership.columns
    ]

    for score_column in score_columns:
        entity_df = _entity_view(membership, score_column)
        if entity_df.empty:
            continue

        entity_tables[score_column] = entity_df
        y_true = entity_df["is_active"].astype(int).to_numpy()
        scores = -entity_df["score"].to_numpy()

        n_total = len(entity_df)
        n_actives = int(y_true.sum())
        n_negatives = n_total - n_actives
        explicit_negative_mask = entity_df["label"].isin(["inactive", "decoy"])
        background_mask = entity_df["label"].eq("background")
        active_ranks = entity_df.loc[entity_df["is_active"], "rank"].to_numpy()
        active_scores = entity_df.loc[entity_df["is_active"], "score"].to_numpy()
        negative_scores = entity_df.loc[~entity_df["is_active"], "score"].to_numpy()

        row: dict[str, object] = {
            "benchmark_id": benchmark_set["benchmark_id"],
            "benchmark_name": benchmark_set.get("name", benchmark_set["benchmark_id"]),
            "benchmark_mode": benchmark_set["benchmark_mode"],
            "target": benchmark_set.get("target", ""),
            "mutation": benchmark_set.get("mutation", ""),
            "score_column": score_column,
            "n_entities": n_total,
            "n_actives": n_actives,
            "n_negatives": n_negatives,
            "n_explicit_negatives": int(explicit_negative_mask.sum()),
            "n_background": int(background_mask.sum()),
            "roc_auc": _roc_auc_from_scores(y_true, scores),
            "average_precision": _average_precision(y_true, scores),
            "median_active_rank": float(np.median(active_ranks)) if len(active_ranks) else np.nan,
            "median_active_percentile": float(np.median(active_ranks / n_total)) if len(active_ranks) else np.nan,
            "best_active_rank": int(active_ranks.min()) if len(active_ranks) else np.nan,
            "worst_active_rank": int(active_ranks.max()) if len(active_ranks) else np.nan,
            "active_score_median": float(np.median(active_scores)) if len(active_scores) else np.nan,
            "negative_score_median": float(np.median(negative_scores)) if len(negative_scores) else np.nan,
            "score_gap_vs_negative_median": (
                float(np.median(negative_scores) - np.median(active_scores))
                if len(active_scores) and len(negative_scores)
                else np.nan
            ),
        }

        sorted_labels = y_true
        for fraction in EF_FRACTIONS:
            ef, nef, hits, top_k = _enrichment_factor(sorted_labels, fraction)
            label = int(fraction * 100)
            row[f"ef{label}"] = ef
            row[f"nef{label}"] = nef
            row[f"hits_top_{label}pct"] = hits
            row[f"top_{label}pct_size"] = top_k

        for top_k in TOP_K_VALUES:
            hits, recall = _topk_recovery(sorted_labels, top_k)
            row[f"hits_top_{top_k}"] = hits
            row[f"recall_top_{top_k}"] = recall

        metrics_rows.append(row)

        calibration = entity_df[["benchmark_id", "benchmark_name", "benchmark_mode", "score_column", "score", "is_active"]].copy()
        calibration = calibration.sort_values("score", ascending=True).reset_index(drop=True)
        calibration["bin"] = pd.qcut(
            calibration.index,
            q=min(CALIBRATION_QUANTILES, len(calibration)),
            labels=False,
            duplicates="drop",
        )
        calibration["bin"] = calibration["bin"].astype(int) + 1
        calibration_summary = (
            calibration.groupby("bin", as_index=False)
            .agg(
                benchmark_id=("benchmark_id", "first"),
                benchmark_name=("benchmark_name", "first"),
                benchmark_mode=("benchmark_mode", "first"),
                score_column=("score_column", "first"),
                n_entities=("score", "size"),
                score_min=("score", "min"),
                score_max=("score", "max"),
                n_actives=("is_active", "sum"),
                active_rate=("is_active", "mean"),
            )
            .sort_values("bin")
        )
        calibration_summary["cumulative_entities"] = calibration_summary["n_entities"].cumsum()
        calibration_summary["cumulative_actives"] = calibration_summary["n_actives"].cumsum()
        calibration_summary["cumulative_active_rate"] = (
            calibration_summary["cumulative_actives"] / calibration_summary["cumulative_entities"]
        )
        calibration_tables.append(calibration_summary)

    return pd.DataFrame(metrics_rows), entity_tables, calibration_tables


def plot_enrichment(all_entity_tables: dict[tuple[str, str], pd.DataFrame]) -> None:
    if not all_entity_tables:
        return

    CHARTS_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6))

    for (benchmark_id, score_column), entity_df in all_entity_tables.items():
        x = np.concatenate([[0.0], entity_df["screen_fraction"].to_numpy()])
        y = np.concatenate([[0.0], entity_df["active_recall"].to_numpy()])
        label = f"{benchmark_id} | {score_column}"
        ax.plot(x, y, linewidth=2.0, label=label)

    ax.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1.2, label="random")
    ax.set_xlabel("Fraction of ranked library screened")
    ax.set_ylabel("Fraction of benchmark actives recovered")
    ax.set_title("Benchmark Active Recovery")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "benchmark_enrichment.png")
    plt.close(fig)


def plot_calibration(calibration_df: pd.DataFrame) -> None:
    if calibration_df.empty:
        return

    CHARTS_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 6))

    for (benchmark_id, score_column), subdf in calibration_df.groupby(["benchmark_id", "score_column"]):
        label = f"{benchmark_id} | {score_column}"
        ax.plot(
            subdf["bin"],
            subdf["active_rate"],
            marker="o",
            linewidth=2.0,
            label=label,
        )

    ax.set_xlabel("Score quantile bin (1 = strongest scores)")
    ax.set_ylabel("Observed active fraction")
    ax.set_title("Benchmark Score Calibration")
    ax.legend(fontsize=8)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    fig.savefig(CHARTS_DIR / "benchmark_calibration.png")
    plt.close(fig)


def write_summary(
    benchmark_sets: pd.DataFrame,
    benchmark_compounds: pd.DataFrame,
    metrics_df: pd.DataFrame,
    active_ranks_df: pd.DataFrame,
) -> None:
    lines = [
        "# ROS1 Benchmark Summary",
        "",
        "This summary is generated from the benchmark registry in `data/benchmarks/`.",
        "",
    ]

    for _, bench in benchmark_sets.iterrows():
        benchmark_id = bench["benchmark_id"]
        set_metrics = metrics_df[metrics_df["benchmark_id"] == benchmark_id].copy()
        set_compounds = benchmark_compounds[benchmark_compounds["benchmark_id"] == benchmark_id].copy()
        set_active_ranks = active_ranks_df[active_ranks_df["benchmark_id"] == benchmark_id].copy()

        lines.extend([
            f"## `{benchmark_id}`",
            "",
            f"- Description: {bench['description']}",
            f"- Mode: `{bench['benchmark_mode']}`",
            f"- Target: {bench.get('target', '')}",
            f"- Mutation/state: {bench.get('mutation', '')}",
            f"- Declared labeled rows: {len(set_compounds)}",
        ])
        if "limitation_note" in bench and pd.notna(bench["limitation_note"]) and str(bench["limitation_note"]).strip():
            lines.append(f"- Limitation note: {bench['limitation_note']}")
        lines.append("")

        if set_metrics.empty:
            lines.extend(["No metrics were computed for this benchmark.", ""])
            continue

        metric_cols = [
            "score_column",
            "n_entities",
            "n_actives",
            "n_negatives",
            "n_explicit_negatives",
            "n_background",
            "roc_auc",
            "average_precision",
            "ef1",
            "nef1",
            "ef5",
            "nef5",
            "ef10",
            "nef10",
            "median_active_rank",
            "median_active_percentile",
            "hits_top_10",
            "hits_top_25",
            "hits_top_50",
            "hits_top_100",
        ]
        metrics_table = set_metrics[metric_cols].copy()
        for col in metrics_table.columns:
            if metrics_table[col].dtype.kind in {"f"}:
                metrics_table[col] = metrics_table[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
        lines.extend(["### Metrics", "", metrics_table.to_markdown(index=False), ""])

        if not set_active_ranks.empty:
            for score_column, subdf in set_active_ranks.groupby("score_column"):
                active_table = subdf[["representative_name", "score", "rank", "screen_fraction"]].copy()
                active_table["score"] = active_table["score"].map(lambda x: f"{x:.3f}")
                active_table["screen_fraction"] = active_table["screen_fraction"].map(lambda x: f"{100 * x:.2f}%")
                lines.extend([
                    f"### Active Ranks: `{score_column}`",
                    "",
                    active_table.to_markdown(index=False),
                    "",
                ])

    (RESULTS_DIR / "benchmark_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    print("=== ROS1 Benchmark / Enrichment Evaluation ===")

    for path in [BENCHMARK_SETS_FILE, BENCHMARK_COMPOUNDS_FILE, RANKED_FILE]:
        if not path.exists():
            raise SystemExit(f"Missing required input: {path}")

    benchmark_sets, benchmark_compounds, ranked = load_inputs()

    all_memberships = []
    all_metrics = []
    all_calibration = []
    all_active_ranks = []
    all_entity_tables: dict[tuple[str, str], pd.DataFrame] = {}
    registry_rows = []

    for _, benchmark_set in benchmark_sets.iterrows():
        benchmark_id = benchmark_set["benchmark_id"]
        membership = build_membership_for_benchmark(benchmark_set, benchmark_compounds, ranked)
        if membership.empty:
            print(f"{benchmark_id}: no matched rows, skipping")
            continue

        metrics_df, entity_tables, calibration_tables = compute_metrics_for_membership(membership, benchmark_set)
        if metrics_df.empty:
            print(f"{benchmark_id}: no metrics computed, skipping")
            continue

        all_memberships.append(membership)
        all_metrics.append(metrics_df)
        all_calibration.extend(calibration_tables)

        for score_column, entity_df in entity_tables.items():
            all_entity_tables[(benchmark_id, score_column)] = entity_df
            active_rows = entity_df[entity_df["is_active"]].copy()
            if not active_rows.empty:
                active_rows["benchmark_id"] = benchmark_id
                active_rows["score_column"] = score_column
                all_active_ranks.append(active_rows)

        registry_rows.append({
            "benchmark_id": benchmark_id,
            "benchmark_name": benchmark_set.get("name", benchmark_id),
            "benchmark_mode": benchmark_set["benchmark_mode"],
            "target": benchmark_set.get("target", ""),
            "mutation": benchmark_set.get("mutation", ""),
            "matched_rows": int(len(membership)),
            "matched_labeled_rows": int(membership["label"].isin(["active", "inactive", "decoy"]).sum()),
            "matched_actives": int(membership["label"].eq("active").sum()),
        })

    if not all_metrics:
        raise SystemExit("No benchmark metrics could be computed from current outputs.")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    CHARTS_DIR.mkdir(parents=True, exist_ok=True)

    membership_df = pd.concat(all_memberships, ignore_index=True)
    metrics_df = pd.concat(all_metrics, ignore_index=True)
    calibration_df = pd.concat(all_calibration, ignore_index=True) if all_calibration else pd.DataFrame()
    active_ranks_df = pd.concat(all_active_ranks, ignore_index=True) if all_active_ranks else pd.DataFrame()
    registry_df = pd.DataFrame(registry_rows)

    membership_cols = [
        "benchmark_id",
        "benchmark_name",
        "benchmark_mode",
        "drug_name",
        "entity_id",
        "label",
        "evidence_tier",
        "source_kind",
        "source_ref",
        "target_benchmark",
        "mutation_benchmark",
        "g2032r_score",
        "composite_score",
        "is_known_ros1_tki",
        "source",
        "target",
        "inchikey",
        "notes",
    ]
    membership_df[[c for c in membership_cols if c in membership_df.columns]].to_csv(
        RESULTS_DIR / "benchmark_membership.csv",
        index=False,
    )
    metrics_df.to_csv(RESULTS_DIR / "benchmark_metrics.csv", index=False)
    registry_df.to_csv(RESULTS_DIR / "benchmark_registry.csv", index=False)
    if not calibration_df.empty:
        calibration_df.to_csv(RESULTS_DIR / "benchmark_calibration.csv", index=False)
    if not active_ranks_df.empty:
        active_ranks_df.to_csv(RESULTS_DIR / "benchmark_active_ranks.csv", index=False)

    for (benchmark_id, score_column), entity_df in all_entity_tables.items():
        entity_df.to_csv(RESULTS_DIR / f"benchmark_entities_{benchmark_id}_{score_column}.csv", index=False)

    plot_enrichment(all_entity_tables)
    plot_calibration(calibration_df)
    write_summary(benchmark_sets, benchmark_compounds, metrics_df, active_ranks_df)

    print(metrics_df.to_string(index=False))
    print("Wrote:")
    print(f"  {RESULTS_DIR / 'benchmark_registry.csv'}")
    print(f"  {RESULTS_DIR / 'benchmark_membership.csv'}")
    print(f"  {RESULTS_DIR / 'benchmark_metrics.csv'}")
    print(f"  {RESULTS_DIR / 'benchmark_calibration.csv'}")
    print(f"  {RESULTS_DIR / 'benchmark_active_ranks.csv'}")
    print(f"  {RESULTS_DIR / 'benchmark_summary.md'}")
    print(f"  {CHARTS_DIR / 'benchmark_enrichment.png'}")
    print(f"  {CHARTS_DIR / 'benchmark_calibration.png'}")


if __name__ == "__main__":
    main()
