#!/usr/bin/env python3
"""
04_analyze.py — Merge docking scores with drug properties, rank, and annotate.

Inputs:
  - data/drug_library.csv
  - results/campaign1_scores.csv
  - results/campaign2_scores.csv
  - results/campaign3_scores.csv

Outputs:
  - results/repurposing_ranked.csv — all drugs ranked by composite score
  - results/top20_hits.csv — best 20 non-TKI hits, fully annotated
  - results/controls.csv — known TKI benchmark scores
"""

import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")


def load_and_merge():
    """Load all data sources and merge into a single DataFrame."""
    library = pd.read_csv(DATA_DIR / "drug_library.csv")
    c1 = pd.read_csv(RESULTS_DIR / "campaign1_scores.csv")
    c1["score"] = pd.to_numeric(c1["score"], errors="coerce")

    # Merge campaign 1 scores
    merged = library.merge(
        c1[["drug_name", "score"]].rename(columns={"score": "c1_score"}),
        on="drug_name",
        how="left",
    )

    # Merge campaign 2 scores (higher accuracy for top drugs)
    c2_path = RESULTS_DIR / "campaign2_scores.csv"
    if c2_path.exists():
        c2 = pd.read_csv(c2_path)
        merged = merged.merge(
            c2[["drug_name", "g2032r_score"]].rename(columns={"g2032r_score": "c2_score"}),
            on="drug_name",
            how="left",
        )
    else:
        merged["c2_score"] = np.nan

    # Use campaign 2 score when available, fall back to campaign 1
    merged["g2032r_score"] = merged["c2_score"].fillna(merged["c1_score"])

    # Merge campaign 3 scores (selectivity)
    c3_path = RESULTS_DIR / "campaign3_scores.csv"
    if c3_path.exists():
        c3 = pd.read_csv(c3_path)
        merged = merged.merge(
            c3[["drug_name", "wt_score", "met_score"]],
            on="drug_name",
            how="left",
        )
    else:
        merged["wt_score"] = np.nan
        merged["met_score"] = np.nan

    return merged


def establish_benchmark(df):
    """Extract known TKI scores as benchmark."""
    controls = df[df["is_known_ros1_tki"] == True].copy()

    if len(controls) == 0:
        print("WARNING: No known TKIs found in results!")
        return controls, -8.0  # fallback threshold

    controls = controls.sort_values("g2032r_score")

    print("\n=== Known TKI Benchmark Scores ===")
    cols = ["drug_name", "g2032r_score", "wt_score", "met_score", "cns_mpo", "cns_penetrant"]
    available = [c for c in cols if c in controls.columns]
    print(controls[available].to_string(index=False))

    best_tki_score = controls["g2032r_score"].min()
    print(f"\nBest TKI score: {best_tki_score:.1f} kcal/mol")

    # Save controls CSV
    controls[available].to_csv(RESULTS_DIR / "controls.csv", index=False)
    print(f"Saved results/controls.csv")

    return controls, best_tki_score


def compute_composite_score(df, best_tki_score):
    """Compute composite score for ranking."""
    # Start with raw G2032R binding score (more negative = better)
    df["composite_score"] = df["g2032r_score"].copy()

    # Bonus -0.5 if CNS MPO >= 4.0
    df.loc[df["cns_mpo"] >= 4.0, "composite_score"] -= 0.5

    # Bonus -0.3 if MET score within 2 kcal/mol of best MET control
    met_controls = df[df["is_known_ros1_tki"] == True]["met_score"].dropna()
    if len(met_controls) > 0:
        best_met = met_controls.min()
        df.loc[
            df["met_score"].notna() & (df["met_score"] <= best_met + 2.0),
            "composite_score",
        ] -= 0.3

    # Bonus -0.2 if G2032R score is better than WT score (mutant-selective)
    df.loc[
        df["wt_score"].notna() & (df["g2032r_score"] < df["wt_score"]),
        "composite_score",
    ] -= 0.2

    return df


def flag_candidates(df, best_tki_score):
    """Flag repurposing candidates."""
    df["is_candidate"] = (
        (df["g2032r_score"] <= best_tki_score + 2.0)
        & (df["is_known_ros1_tki"] == False)
        & df["g2032r_score"].notna()
    )

    n_candidates = df["is_candidate"].sum()
    print(f"\nRepurposing candidates (within 2 kcal/mol of best TKI): {n_candidates}")
    return df


def annotate_top20(df):
    """Annotate top 20 non-TKI hits."""
    # Filter out known TKIs and sort by composite score
    non_tki = df[df["is_known_ros1_tki"] == False].copy()
    top20 = non_tki.nsmallest(20, "composite_score")

    cols = [
        "drug_name", "catalog_id", "smiles", "g2032r_score", "wt_score", "met_score",
        "composite_score", "mw", "tpsa", "clogp", "hbd", "hba", "rotatable_bonds",
        "cns_mpo", "cns_penetrant", "target", "cas", "source", "is_candidate",
    ]
    available = [c for c in cols if c in top20.columns]
    top20 = top20[available]

    top20.to_csv(RESULTS_DIR / "top20_hits.csv", index=False)
    print(f"\nSaved results/top20_hits.csv")
    return top20


def main():
    print("=== ROS1 Drug Repurposing Analysis ===\n")

    # Load and merge all data
    print("Loading and merging data...")
    df = load_and_merge()
    n_scored = df["g2032r_score"].notna().sum()
    print(f"  {len(df)} drugs in library, {n_scored} with docking scores")

    # Establish benchmark from known TKIs
    controls, best_tki_score = establish_benchmark(df)

    # Compute composite score
    df = compute_composite_score(df, best_tki_score)

    # Flag candidates
    df = flag_candidates(df, best_tki_score)

    # Rank all drugs
    ranked = df.sort_values("composite_score").copy()
    rank_cols = [
        "drug_name", "g2032r_score", "wt_score", "met_score", "composite_score",
        "mw", "tpsa", "clogp", "cns_mpo", "cns_penetrant", "target",
        "is_known_ros1_tki", "is_candidate", "source",
    ]
    available = [c for c in rank_cols if c in ranked.columns]
    ranked[available].to_csv(RESULTS_DIR / "repurposing_ranked.csv", index=False)
    print(f"\nSaved results/repurposing_ranked.csv ({len(ranked)} drugs)")

    # Annotate top 20
    top20 = annotate_top20(df)

    # Print summary
    print("\n" + "=" * 60)
    print("TOP 20 REPURPOSING CANDIDATES")
    print("=" * 60)
    display_cols = ["drug_name", "g2032r_score", "composite_score", "cns_mpo", "target"]
    available = [c for c in display_cols if c in top20.columns]
    print(top20[available].to_string(index=False))

    print(f"\n{'=' * 60}")
    print(f"Total drugs screened: {n_scored}")
    print(f"Known TKI benchmark: {best_tki_score:.1f} kcal/mol")
    print(f"Candidates within 2 kcal/mol: {df['is_candidate'].sum()}")
    cns_candidates = df[df["is_candidate"] & df["cns_penetrant"]].shape[0]
    print(f"CNS-penetrant candidates: {cns_candidates}")


if __name__ == "__main__":
    main()
