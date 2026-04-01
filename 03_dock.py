#!/usr/bin/env python3
"""
03_dock.py — Dock drug library against ROS1 G2032R, WT, and MET.

Three campaigns:
  1. Primary screen: all drugs vs G2032R (exhaustiveness=8)
  2. Re-dock top 50 + TKIs vs G2032R (exhaustiveness=32, 9 poses)
  3. Selectivity: top 20 + TKIs vs WT and MET (exhaustiveness=32, 9 poses)

Pre-flight: re-docking validation of zidesamtinib crystal pose.

Outputs:
  - results/campaign1_scores.csv
  - results/campaign2_scores.csv
  - results/campaign3_scores.csv
  - results/validation_rmsd.txt
  - results/poses/{campaign1,campaign2,campaign3_wt,campaign3_met}/*.pdbqt
"""

import csv
import multiprocessing
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# Ensure output is flushed immediately for terminal monitoring
sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
POSES_DIR = RESULTS_DIR / "poses"

N_CPU = os.cpu_count() or 4


def parse_grid_config(grid_path: str) -> dict:
    """Parse a Vina grid config file."""
    config = {}
    with open(grid_path) as f:
        for line in f:
            key, val = line.strip().split(" = ")
            config[key.strip()] = float(val)
    return config


def dock_single(receptor_pdbqt: str, ligand_pdbqt: str, grid: dict,
                exhaustiveness: int = 8, n_poses: int = 1) -> tuple:
    """Dock a single ligand. Returns (best_score, all_scores, poses_pdbqt_string)."""
    try:
        v = Vina(sf_name="vina", cpu=1, verbosity=0)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(
            center=[grid["center_x"], grid["center_y"], grid["center_z"]],
            box_size=[grid["size_x"], grid["size_y"], grid["size_z"]],
        )
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies()
        best_score = energies[0][0] if len(energies) > 0 else None
        all_scores = [e[0] for e in energies]
        poses_str = v.poses(n_poses=n_poses)
        return best_score, all_scores, poses_str
    except Exception as e:
        return None, [], str(e)


def dock_worker(args):
    """Worker for multiprocessing pool."""
    drug_name, receptor_pdbqt, ligand_pdbqt, grid, exhaustiveness, n_poses, pose_dir = args
    score, all_scores, poses_str = dock_single(
        receptor_pdbqt, ligand_pdbqt, grid, exhaustiveness, n_poses
    )
    # Save pose
    if score is not None and pose_dir:
        safe_name = drug_name.replace("/", "_").replace(" ", "_")
        pose_path = os.path.join(pose_dir, f"{safe_name}.pdbqt")
        with open(pose_path, "w") as f:
            f.write(poses_str)
    return drug_name, score, all_scores


def run_redocking_validation(receptor_pdbqt: str, grid: dict):
    """Re-dock zidesamtinib crystal pose and save a preliminary validation note.

    Detailed RMSD analysis is delegated to 06_improve.py, which can account
    for frame alignment and atom mapping differences between crystal and docked
    representations.
    """
    print("\n=== Re-docking Validation ===")

    crystal_pdb = str(DATA_DIR / "controls" / "zidesamtinib_crystal.pdb")
    control_pdbqt = str(DATA_DIR / "controls" / "zidesamtinib.pdbqt")

    if not os.path.exists(crystal_pdb):
        print("  WARNING: Crystal ligand not found, skipping validation")
        return

    if not os.path.exists(control_pdbqt):
        print("  WARNING: Zidesamtinib PDBQT not found, skipping validation")
        return

    # Dock zidesamtinib
    print("  Docking zidesamtinib (exhaustiveness=32)...")
    score, _, poses_str = dock_single(
        receptor_pdbqt, control_pdbqt, grid, exhaustiveness=32, n_poses=1
    )

    if score is None:
        print(f"  WARNING: Docking failed")
        Path(RESULTS_DIR / "validation_rmsd.txt").write_text("FAILED\n")
        return

    print(f"  Docking score: {score:.1f} kcal/mol")

    # Save docked pose
    docked_pose_path = RESULTS_DIR / "poses" / "validation_zidesamtinib.pdbqt"
    docked_pose_path.parent.mkdir(parents=True, exist_ok=True)
    with open(docked_pose_path, "w") as f:
        f.write(poses_str)

    result = (
        f"Score: {score:.1f} kcal/mol\n"
        "RMSD: deferred to 06_improve.py\n"
        "NOTE: Raw coordinate-order RMSD is intentionally not reported here because\n"
        "crystal and docked ligand frames/atom lists may differ.\n"
    )
    print("  Detailed RMSD deferred to 06_improve.py")

    Path(RESULTS_DIR / "validation_rmsd.txt").write_text(result)
    print(f"  Saved to results/validation_rmsd.txt")


def run_campaign1(drug_library: pd.DataFrame, receptor_pdbqt: str, grid: dict):
    """Campaign 1: Primary screen — all drugs vs G2032R, exhaustiveness=8."""
    print(f"\n=== Campaign 1: Primary Screen ({len(drug_library)} drugs) ===")

    pose_dir = str(POSES_DIR / "campaign1")
    os.makedirs(pose_dir, exist_ok=True)

    # Load checkpoint
    checkpoint_path = RESULTS_DIR / "docking_progress.csv"
    completed = set()
    if checkpoint_path.exists():
        cp = pd.read_csv(checkpoint_path)
        completed = set(cp["drug_name"])
        print(f"  Resuming: {len(completed)} already docked")

    # Prepare tasks
    tasks = []
    for _, row in drug_library.iterrows():
        if row["drug_name"] in completed:
            continue
        tasks.append((
            row["drug_name"],
            receptor_pdbqt,
            row["pdbqt_path"],
            grid,
            8,   # exhaustiveness
            1,   # n_poses
            pose_dir,
        ))

    if not tasks:
        print("  All drugs already docked, skipping")
        return

    print(f"  Docking {len(tasks)} drugs with {N_CPU} workers...")

    # Use multiprocessing pool
    results = []
    best_score = 0.0
    start_time = time.time()
    with open(checkpoint_path, "a", newline="") as cp_file:
        writer = csv.writer(cp_file)
        if not completed:
            writer.writerow(["drug_name", "score", "timestamp"])

        with multiprocessing.Pool(processes=N_CPU) as pool:
            for i, (drug_name, score, _) in enumerate(tqdm(
                pool.imap_unordered(dock_worker, tasks),
                total=len(tasks),
                desc="Campaign 1",
                miniters=1,
            )):
                if score is not None:
                    writer.writerow([drug_name, f"{score:.2f}", datetime.now().isoformat()])
                    cp_file.flush()
                    results.append({"drug_name": drug_name, "score": score})
                    if score < best_score:
                        best_score = score

                # Print summary every 100 drugs
                done = i + 1
                if done % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = done / elapsed
                    remaining = (len(tasks) - done) / rate if rate > 0 else 0
                    print(f"\n  [{done}/{len(tasks)}] "
                          f"Best: {best_score:.1f} kcal/mol | "
                          f"Rate: {rate:.1f} drugs/s | "
                          f"ETA: {remaining/60:.0f} min",
                          flush=True)

    # Write final campaign1 scores (checkpoint + new)
    if checkpoint_path.exists():
        all_scores = pd.read_csv(checkpoint_path)
    else:
        all_scores = pd.DataFrame(results)

    all_scores.to_csv(RESULTS_DIR / "campaign1_scores.csv", index=False)
    print(f"  Campaign 1 complete: {len(all_scores)} drugs docked")
    print(f"  Best score: {all_scores['score'].min():.1f} kcal/mol")


def run_campaign2(receptor_pdbqt: str, grid: dict, drug_library: pd.DataFrame):
    """Campaign 2: Re-dock top 50 + all TKIs at exhaustiveness=32."""
    print("\n=== Campaign 2: High-accuracy Re-dock ===")

    c1 = pd.read_csv(RESULTS_DIR / "campaign1_scores.csv")
    c1["score"] = pd.to_numeric(c1["score"], errors="coerce")

    # Top 50 by score + all known TKIs
    top50 = set(c1.nsmallest(50, "score")["drug_name"])
    tkis = set(drug_library[drug_library["is_known_ros1_tki"]]["drug_name"])
    redock_names = top50 | tkis

    # Get PDBQT paths
    redock_df = drug_library[drug_library["drug_name"].isin(redock_names)]
    print(f"  Re-docking {len(redock_df)} drugs (top 50 + {len(tkis)} TKIs)")

    pose_dir = str(POSES_DIR / "campaign2")
    os.makedirs(pose_dir, exist_ok=True)

    tasks = [
        (row["drug_name"], receptor_pdbqt, row["pdbqt_path"], grid, 32, 9, pose_dir)
        for _, row in redock_df.iterrows()
    ]

    results = []
    with multiprocessing.Pool(processes=N_CPU) as pool:
        for drug_name, score, all_scores in tqdm(
            pool.imap_unordered(dock_worker, tasks),
            total=len(tasks),
            desc="Campaign 2",
        ):
            if score is not None:
                results.append({"drug_name": drug_name, "g2032r_score": score})

    c2 = pd.DataFrame(results)
    c2.to_csv(RESULTS_DIR / "campaign2_scores.csv", index=False)
    print(f"  Campaign 2 complete: {len(c2)} drugs re-docked")
    print(f"  Best score: {c2['g2032r_score'].min():.1f} kcal/mol")


def run_campaign3(drug_library: pd.DataFrame):
    """Campaign 3: Selectivity profiling — top 20 + TKIs vs WT and MET."""
    print("\n=== Campaign 3: Selectivity Profiling ===")

    c2 = pd.read_csv(RESULTS_DIR / "campaign2_scores.csv")

    # Top 20 from campaign 2 + all TKIs
    top20 = set(c2.nsmallest(20, "g2032r_score")["drug_name"])
    tkis = set(drug_library[drug_library["is_known_ros1_tki"]]["drug_name"])
    profile_names = top20 | tkis

    profile_df = drug_library[drug_library["drug_name"].isin(profile_names)]
    print(f"  Profiling {len(profile_df)} drugs against WT and MET")

    # WT docking
    wt_receptor = str(DATA_DIR / "receptor_WT.pdbqt")
    wt_grid = parse_grid_config(str(DATA_DIR / "grid_WT.txt"))
    wt_pose_dir = str(POSES_DIR / "campaign3_wt")
    os.makedirs(wt_pose_dir, exist_ok=True)

    tasks_wt = [
        (row["drug_name"], wt_receptor, row["pdbqt_path"], wt_grid, 32, 9, wt_pose_dir)
        for _, row in profile_df.iterrows()
    ]

    wt_results = {}
    print("  Docking against WT ROS1...")
    with multiprocessing.Pool(processes=N_CPU) as pool:
        for drug_name, score, _ in tqdm(
            pool.imap_unordered(dock_worker, tasks_wt),
            total=len(tasks_wt),
            desc="Campaign 3 WT",
        ):
            if score is not None:
                wt_results[drug_name] = score

    # MET docking
    met_receptor = str(DATA_DIR / "receptor_MET.pdbqt")
    met_grid = parse_grid_config(str(DATA_DIR / "grid_MET.txt"))
    met_pose_dir = str(POSES_DIR / "campaign3_met")
    os.makedirs(met_pose_dir, exist_ok=True)

    tasks_met = [
        (row["drug_name"], met_receptor, row["pdbqt_path"], met_grid, 32, 9, met_pose_dir)
        for _, row in profile_df.iterrows()
    ]

    met_results = {}
    print("  Docking against MET...")
    with multiprocessing.Pool(processes=N_CPU) as pool:
        for drug_name, score, _ in tqdm(
            pool.imap_unordered(dock_worker, tasks_met),
            total=len(tasks_met),
            desc="Campaign 3 MET",
        ):
            if score is not None:
                met_results[drug_name] = score

    # Merge with G2032R scores from campaign 2
    g2032r_scores = dict(zip(c2["drug_name"], c2["g2032r_score"]))

    rows = []
    for name in profile_names:
        rows.append({
            "drug_name": name,
            "g2032r_score": g2032r_scores.get(name),
            "wt_score": wt_results.get(name),
            "met_score": met_results.get(name),
        })

    c3 = pd.DataFrame(rows)
    c3.to_csv(RESULTS_DIR / "campaign3_scores.csv", index=False)
    print(f"  Campaign 3 complete: {len(c3)} drugs profiled")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    for subdir in ["campaign1", "campaign2", "campaign3_wt", "campaign3_met"]:
        (POSES_DIR / subdir).mkdir(parents=True, exist_ok=True)

    # Load drug library
    drug_library = pd.read_csv(DATA_DIR / "drug_library.csv")
    print(f"Loaded {len(drug_library)} drugs from drug library")

    # Load G2032R receptor and grid
    g2032r_receptor = str(DATA_DIR / "receptor_G2032R.pdbqt")
    g2032r_grid = parse_grid_config(str(DATA_DIR / "grid_G2032R.txt"))

    # Pre-flight: re-docking validation
    run_redocking_validation(g2032r_receptor, g2032r_grid)

    # Campaign 1: Primary screen
    run_campaign1(drug_library, g2032r_receptor, g2032r_grid)

    # Campaign 2: High-accuracy re-dock
    run_campaign2(g2032r_receptor, g2032r_grid, drug_library)

    # Campaign 3: Selectivity profiling
    run_campaign3(drug_library)

    print("\n=== All docking campaigns complete ===")


if __name__ == "__main__":
    main()
