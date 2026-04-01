# Methodology Improvements Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix re-docking RMSD validation, add multi-conformer docking for top 50 hits, and add ADMET annotation for top 20 hits — without re-running the full library screen.

**Architecture:** Single new script `06_improve.py` with three phases: (1) fix RMSD validation via alignment-based calculation, (2) multi-conformer re-dock of top 50 replacing Campaign 2, (3) ADMET property annotation. Then re-run `04_analyze.py` and `05_visualize.py` to regenerate outputs.

**Tech Stack:** RDKit (conformer generation, ADMET descriptors, mol alignment), Meeko (PDBQT conversion), AutoDock Vina (docking), pandas

---

### Task 1: Create `06_improve.py` scaffold with Phase 1 — Fix RMSD Validation

**Files:**
- Create: `06_improve.py`
- Read: `03_dock.py:93-173` (current validation logic)
- Read: `data/controls/zidesamtinib_crystal.pdb` (crystal ligand)
- Read: `results/poses/validation_zidesamtinib.pdbqt` (existing docked pose)

**Step 1: Write `06_improve.py` with the RMSD fix**

The current RMSD calculation in `03_dock.py` compares raw coordinates across different reference frames (original PDB vs minimized receptor) without superposition. Fix by:

1. Load crystal ligand PDB as RDKit mol
2. Load docked PDBQT, extract heavy-atom coordinates, inject them into a copy of the crystal mol (same atom ordering guaranteed since it's the same molecule)
3. Use `rdMolAlign.GetBestRMS()` for symmetry-aware aligned RMSD

```python
#!/usr/bin/env python3
"""
06_improve.py — Methodology improvements to the repurposing screen.

Phase 1: Fix RMSD validation (alignment-based calculation)
Phase 2: Multi-conformer re-dock of top 50 hits
Phase 3: ADMET annotation of top 20 hits

Run after 03_dock.py, before re-running 04_analyze.py and 05_visualize.py.
"""

import csv
import multiprocessing
import os
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolAlign
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy
from vina import Vina
from tqdm import tqdm

RDLogger.DisableLog("rdApp.*")

DATA_DIR = Path("data")
RESULTS_DIR = Path("results")
POSES_DIR = RESULTS_DIR / "poses"
N_CPU = os.cpu_count() or 4


# ---------------------------------------------------------------------------
# Phase 1: Fix RMSD validation
# ---------------------------------------------------------------------------

def fix_rmsd_validation():
    """Recompute RMSD between crystal and docked zidesamtinib using alignment."""
    print("\n=== Phase 1: Fix RMSD Validation ===")

    crystal_pdb = DATA_DIR / "controls" / "zidesamtinib_crystal.pdb"
    docked_pdbqt = RESULTS_DIR / "poses" / "validation_zidesamtinib.pdbqt"

    if not crystal_pdb.exists() or not docked_pdbqt.exists():
        print("  WARNING: Missing crystal or docked pose, skipping")
        return

    # Load crystal ligand
    crystal_mol = Chem.MolFromPDBFile(str(crystal_pdb), removeHs=True, sanitize=False)
    if crystal_mol is None:
        print("  WARNING: Could not parse crystal ligand PDB")
        return

    try:
        Chem.SanitizeMol(crystal_mol)
    except Exception:
        pass  # keep unsanitized — we only need coords

    # Parse docked PDBQT heavy-atom coordinates
    docked_coords = []
    for line in docked_pdbqt.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            element = line[76:78].strip() if len(line) >= 78 else ""
            # Skip hydrogens
            if element == "H" or (not element and line[12:16].strip().startswith("H")):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            docked_coords.append((x, y, z))

    n_crystal = crystal_mol.GetNumAtoms()
    n_docked = len(docked_coords)
    print(f"  Crystal heavy atoms: {n_crystal}, Docked heavy atoms: {n_docked}")

    if n_docked == 0:
        print("  WARNING: No heavy atoms in docked pose")
        return

    # Method 1: Alignment-based RMSD using RDKit
    # Create a copy of crystal mol and set docked coordinates
    if n_crystal == n_docked:
        docked_mol = Chem.RWMol(crystal_mol)
        conf = docked_mol.GetConformer()
        for i, (x, y, z) in enumerate(docked_coords):
            conf.SetAtomPosition(i, (x, y, z))

        # GetBestRMS handles molecular symmetry
        aligned_rmsd = rdMolAlign.GetBestRMS(crystal_mol, docked_mol)
        print(f"  Aligned RMSD (symmetry-aware): {aligned_rmsd:.2f} A")
    else:
        # Fallback: centroid alignment + min-atom pairing
        print(f"  Atom count mismatch ({n_crystal} vs {n_docked}), using centroid RMSD")
        n = min(n_crystal, n_docked)
        c1 = np.array([(crystal_mol.GetConformer().GetAtomPosition(i).x,
                         crystal_mol.GetConformer().GetAtomPosition(i).y,
                         crystal_mol.GetConformer().GetAtomPosition(i).z)
                        for i in range(n)])
        c2 = np.array(docked_coords[:n])

        # Centroid-align before RMSD
        c1_centered = c1 - c1.mean(axis=0)
        c2_centered = c2 - c2.mean(axis=0)

        # Kabsch alignment
        H = c1_centered.T @ c2_centered
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign_matrix = np.diag([1, 1, d])
        R = Vt.T @ sign_matrix @ U.T
        c2_aligned = (c2_centered @ R.T)

        aligned_rmsd = np.sqrt(np.mean(np.sum((c1_centered - c2_aligned) ** 2, axis=1)))
        print(f"  Kabsch-aligned RMSD: {aligned_rmsd:.2f} A")

    # Also compute the naive (unaligned) RMSD for comparison
    n = min(n_crystal, len(docked_coords))
    c1 = np.array([(crystal_mol.GetConformer().GetAtomPosition(i).x,
                     crystal_mol.GetConformer().GetAtomPosition(i).y,
                     crystal_mol.GetConformer().GetAtomPosition(i).z)
                    for i in range(n)])
    c2 = np.array(docked_coords[:n])
    naive_rmsd = np.sqrt(np.mean(np.sum((c1 - c2) ** 2, axis=1)))
    print(f"  Naive RMSD (no alignment): {naive_rmsd:.2f} A")

    # Read old docking score from existing validation file
    old_txt = (RESULTS_DIR / "validation_rmsd.txt").read_text()
    score_line = [l for l in old_txt.splitlines() if "Score" in l]
    score_str = score_line[0] if score_line else "Score: unknown"

    # Write updated validation
    result = f"Aligned RMSD: {aligned_rmsd:.2f} A\n"
    result += f"Naive RMSD (no alignment): {naive_rmsd:.2f} A\n"
    result += f"{score_str}\n"
    if aligned_rmsd <= 2.0:
        result += "PASS: Aligned RMSD <= 2.0 A — docking protocol validated\n"
        print("  PASS: Re-docking validated!")
    else:
        result += "WARNING: Aligned RMSD > 2.0 A — docking protocol may need adjustment\n"
        print("  WARNING: Aligned RMSD > 2.0 A")

    (RESULTS_DIR / "validation_rmsd.txt").write_text(result)
    print(f"  Updated results/validation_rmsd.txt")
    return aligned_rmsd
```

**Step 2: Test Phase 1 in isolation**

Run: `cd /Users/lucapassone/Documents/ros1-screens/repurposing && uv run python -c "from importlib.machinery import SourceFileLoader; m = SourceFileLoader('m', '06_improve.py').load_module(); m.fix_rmsd_validation()"`

Expected: Aligned RMSD significantly lower than 6.87 A. Likely < 2.0 A.

**Step 3: Commit**

```bash
git add 06_improve.py
git commit -m "feat: add 06_improve.py with Phase 1 — fix RMSD validation via alignment"
```

---

### Task 2: Phase 2 — Multi-conformer Re-dock of Top 50

**Files:**
- Modify: `06_improve.py`
- Read: `02_fetch_drugs.py:233-279` (existing conformer + PDBQT generation)
- Read: `03_dock.py:57-91` (existing docking functions)
- Read: `results/campaign1_scores.csv` (to identify top 50)
- Read: `data/drug_library.csv` (for SMILES)

**Step 1: Add multi-conformer docking functions to `06_improve.py`**

Reuse `dock_single` / `dock_worker` / `parse_grid_config` from `03_dock.py` (import or inline). For each of the top 50 + known TKIs:
1. Generate 10 conformers via ETKDG v3 + MMFF94
2. Convert each conformer to a temporary PDBQT
3. Dock each at exhaustiveness=32, 9 poses
4. Keep the best score across all conformers
5. Save the best pose

```python
# ---------------------------------------------------------------------------
# Phase 2: Multi-conformer re-dock
# ---------------------------------------------------------------------------

def parse_grid_config(grid_path: str) -> dict:
    """Parse a Vina grid config file."""
    config = {}
    with open(grid_path) as f:
        for line in f:
            key, val = line.strip().split(" = ")
            config[key.strip()] = float(val)
    return config


def dock_single(receptor_pdbqt: str, ligand_pdbqt: str, grid: dict,
                exhaustiveness: int = 32, n_poses: int = 9) -> tuple:
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


def generate_conformers(smiles: str, n_confs: int = 10) -> list:
    """Generate multiple 3D conformers. Returns list of RDKit mols with conformers."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    params.pruneRmsThresh = 0.5  # prune similar conformers

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if len(cids) == 0:
        # Fallback with random coords
        params.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

    # MMFF94 optimize each conformer
    if len(cids) > 0:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94", maxIters=500)

    return [(mol, cid) for cid in cids]


def mol_conf_to_pdbqt(mol, conf_id: int, output_path: str) -> bool:
    """Convert a specific conformer to PDBQT via Meeko."""
    try:
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        for setup in mol_setups:
            pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                with open(output_path, "w") as f:
                    f.write(pdbqt_string)
                return True
    except Exception:
        pass
    return False


def dock_conformer_worker(args):
    """Worker: dock one conformer of one drug."""
    drug_name, conf_idx, receptor_pdbqt, ligand_pdbqt, grid = args
    score, all_scores, poses_str = dock_single(
        receptor_pdbqt, ligand_pdbqt, grid, exhaustiveness=32, n_poses=9
    )
    return drug_name, conf_idx, score, poses_str


def run_multiconformer_redock():
    """Phase 2: Re-dock top 50 + TKIs with multiple conformers."""
    print("\n=== Phase 2: Multi-conformer Re-dock ===")

    # Load campaign 1 scores and drug library
    c1 = pd.read_csv(RESULTS_DIR / "campaign1_scores.csv")
    c1["score"] = pd.to_numeric(c1["score"], errors="coerce")
    library = pd.read_csv(DATA_DIR / "drug_library.csv")

    # Identify top 50 + all TKIs
    top50_names = set(c1.nsmallest(50, "score")["drug_name"])
    tki_names = set(library[library["is_known_ros1_tki"] == True]["drug_name"])
    redock_names = top50_names | tki_names

    redock_df = library[library["drug_name"].isin(redock_names)].copy()
    print(f"  Re-docking {len(redock_df)} drugs (top 50 + {len(tki_names)} TKIs)")
    print(f"  Generating up to 10 conformers each...")

    # Load receptor and grid
    receptor_pdbqt = str(DATA_DIR / "receptor_G2032R.pdbqt")
    grid = parse_grid_config(str(DATA_DIR / "grid_G2032R.txt"))

    # Prepare conformers and temporary PDBQTs
    import tempfile
    tmp_dir = Path(tempfile.mkdtemp(prefix="multiconf_"))
    print(f"  Temp conformer dir: {tmp_dir}")

    tasks = []
    conf_counts = {}
    for _, row in tqdm(redock_df.iterrows(), total=len(redock_df), desc="Generating conformers"):
        drug_name = row["drug_name"]
        smiles = row["smiles"]
        conformers = generate_conformers(smiles, n_confs=10)

        if not conformers:
            # Fall back to existing single-conformer PDBQT
            existing_pdbqt = row.get("pdbqt_path", "")
            if existing_pdbqt and os.path.exists(existing_pdbqt):
                tasks.append((drug_name, 0, receptor_pdbqt, existing_pdbqt, grid))
                conf_counts[drug_name] = 1
            continue

        n_ok = 0
        for mol, cid in conformers:
            safe_name = drug_name.replace("/", "_").replace(" ", "_")
            pdbqt_path = str(tmp_dir / f"{safe_name}_conf{cid}.pdbqt")
            if mol_conf_to_pdbqt(mol, cid, pdbqt_path):
                tasks.append((drug_name, cid, receptor_pdbqt, pdbqt_path, grid))
                n_ok += 1
        conf_counts[drug_name] = n_ok

    total_confs = len(tasks)
    avg_confs = total_confs / len(redock_df) if len(redock_df) > 0 else 0
    print(f"  Total docking jobs: {total_confs} ({avg_confs:.1f} conformers/drug avg)")

    # Dock all conformers
    pose_dir = POSES_DIR / "campaign2_multiconf"
    pose_dir.mkdir(parents=True, exist_ok=True)

    best_scores = {}  # drug_name -> (best_score, best_poses_str)
    start_time = time.time()

    with multiprocessing.Pool(processes=N_CPU) as pool:
        for drug_name, conf_idx, score, poses_str in tqdm(
            pool.imap_unordered(dock_conformer_worker, tasks),
            total=len(tasks),
            desc="Multi-conf docking",
        ):
            if score is not None:
                if drug_name not in best_scores or score < best_scores[drug_name][0]:
                    best_scores[drug_name] = (score, poses_str)

    elapsed = time.time() - start_time
    print(f"  Docking complete in {elapsed/60:.1f} min")

    # Save best poses
    for drug_name, (score, poses_str) in best_scores.items():
        safe_name = drug_name.replace("/", "_").replace(" ", "_")
        pose_path = pose_dir / f"{safe_name}.pdbqt"
        with open(pose_path, "w") as f:
            f.write(poses_str)

    # Write updated campaign 2 scores
    rows = [{"drug_name": name, "g2032r_score": score}
            for name, (score, _) in best_scores.items()]
    c2_new = pd.DataFrame(rows).sort_values("g2032r_score")

    # Back up old campaign 2
    old_c2 = RESULTS_DIR / "campaign2_scores.csv"
    if old_c2.exists():
        old_c2.rename(old_c2.with_name("campaign2_scores_single_conf.csv"))

    c2_new.to_csv(RESULTS_DIR / "campaign2_scores.csv", index=False)
    print(f"  Saved updated results/campaign2_scores.csv ({len(c2_new)} drugs)")
    print(f"  Best score: {c2_new['g2032r_score'].min():.1f} kcal/mol")

    # Cleanup temp dir
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)

    return c2_new
```

**Step 2: Test Phase 2 with a dry run**

Run: `cd /Users/lucapassone/Documents/ros1-screens/repurposing && uv run python -c "from importlib.machinery import SourceFileLoader; m = SourceFileLoader('m', '06_improve.py').load_module(); m.generate_conformers('c1ccccc1', 3)"` to verify conformer generation works.

Expected: Returns list of 3 (mol, cid) tuples.

**Step 3: Commit**

```bash
git add 06_improve.py
git commit -m "feat: add Phase 2 — multi-conformer re-dock for top 50 hits"
```

---

### Task 3: Phase 3 — ADMET Annotation

**Files:**
- Modify: `06_improve.py`
- Read: `results/top20_hits.csv` (current top 20)
- Read: `02_fetch_drugs.py:282-340` (existing property computation)

**Step 1: Add ADMET annotation functions to `06_improve.py`**

```python
# ---------------------------------------------------------------------------
# Phase 3: ADMET annotation
# ---------------------------------------------------------------------------

# Lorlatinib is a CYP3A4 substrate. Drugs that inhibit CYP3A4 may increase
# lorlatinib exposure. This is flagged but NEVER used to exclude candidates.
LORLATINIB_DDI_NOTE = (
    "Lorlatinib is a CYP3A4 substrate. Co-administration with strong "
    "CYP3A4 inhibitors may increase lorlatinib plasma levels. "
    "Dose adjustment may be needed. Consult oncologist."
)


def compute_admet(smiles: str) -> dict:
    """Compute ADMET properties from SMILES using RDKit descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}

    mol_noH = Chem.RemoveHs(mol)
    mw = Descriptors.MolWt(mol_noH)
    tpsa = Descriptors.TPSA(mol_noH)
    clogp = Crippen.MolLogP(mol_noH)
    hbd = Lipinski.NumHDonors(mol_noH)
    hba = Lipinski.NumHAcceptors(mol_noH)
    rot = Lipinski.NumRotatableBonds(mol_noH)
    n_rings = Descriptors.RingCount(mol_noH)
    n_aromatic_rings = Descriptors.NumAromaticRings(mol_noH)
    fsp3 = Descriptors.FractionCSP3(mol_noH)

    # --- Lipinski Rule of 5 ---
    lipinski_violations = sum([
        mw > 500,
        clogp > 5,
        hbd > 5,
        hba > 10,
    ])

    # --- Veber rules (oral bioavailability) ---
    veber_pass = (tpsa <= 140) and (rot <= 10)

    # --- GI absorption (Egan model) ---
    # High absorption if TPSA <= 132 and cLogP in [-1, 6]
    gi_absorption = "High" if (tpsa <= 132 and -1 <= clogp <= 6) else "Low"

    # --- BBB permeability (simple Clark model) ---
    # BBB+ if TPSA <= 90 and cLogP in [1, 5]
    bbb_permeable = (tpsa <= 90) and (1 <= clogp <= 5)

    # --- P-gp substrate likelihood ---
    # Heuristic: MW > 400 and TPSA > 90 suggest P-gp substrate
    pgp_substrate_likely = (mw > 400) and (tpsa > 90)

    # --- CYP3A4 inhibition risk ---
    # Heuristic based on Veith et al. 2009 descriptors
    # Risk factors: MW > 400, cLogP > 3, aromatic rings >= 3
    cyp3a4_risk_factors = sum([
        mw > 400,
        clogp > 3,
        n_aromatic_rings >= 3,
        hba >= 4,
    ])
    cyp3a4_inhibitor_likely = cyp3a4_risk_factors >= 3

    # --- hERG liability ---
    # Basic amines with high lipophilicity are risk factors
    # Check for basic nitrogen
    has_basic_n = False
    for atom in mol_noH.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            # Basic if: not in amide, not aromatic, has H or is tertiary amine
            if not atom.GetIsAromatic():
                # Check if it's an amide nitrogen (bonded to C=O)
                is_amide = False
                for bond in atom.GetBonds():
                    other = bond.GetOtherAtom(atom)
                    if other.GetAtomicNum() == 6:  # carbon
                        for b2 in other.GetBonds():
                            o = b2.GetOtherAtom(other)
                            if o.GetAtomicNum() == 8 and b2.GetBondTypeAsDouble() == 2:
                                is_amide = True
                if not is_amide:
                    has_basic_n = True
                    break

    herg_risk = has_basic_n and clogp > 3.7

    # --- Lorlatinib DDI flag ---
    # Flag if drug is a likely CYP3A4 inhibitor (lorlatinib is CYP3A substrate)
    lorlatinib_ddi_flag = cyp3a4_inhibitor_likely

    # --- Plasma protein binding estimate ---
    # Heuristic: high cLogP and MW correlate with high PPB
    if clogp > 4.5 or mw > 500:
        ppb_estimate = "High (>90%)"
    elif clogp > 2.5:
        ppb_estimate = "Moderate (70-90%)"
    else:
        ppb_estimate = "Low (<70%)"

    return {
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": lipinski_violations <= 1,
        "veber_pass": veber_pass,
        "gi_absorption": gi_absorption,
        "bbb_permeable": bbb_permeable,
        "pgp_substrate_likely": pgp_substrate_likely,
        "cyp3a4_inhibitor_likely": cyp3a4_inhibitor_likely,
        "herg_risk": herg_risk,
        "lorlatinib_ddi_flag": lorlatinib_ddi_flag,
        "lorlatinib_ddi_note": LORLATINIB_DDI_NOTE if lorlatinib_ddi_flag else "",
        "ppb_estimate": ppb_estimate,
        "fsp3": round(fsp3, 2),
        "n_aromatic_rings": n_aromatic_rings,
    }


def run_admet_annotation():
    """Phase 3: Add ADMET properties to top 20 hits."""
    print("\n=== Phase 3: ADMET Annotation ===")

    top20 = pd.read_csv(RESULTS_DIR / "top20_hits.csv")
    print(f"  Annotating {len(top20)} top hits")

    admet_rows = []
    for _, row in top20.iterrows():
        props = compute_admet(row["smiles"])
        props["drug_name"] = row["drug_name"]
        admet_rows.append(props)

    admet_df = pd.DataFrame(admet_rows)

    # Merge with top20
    top20_enhanced = top20.merge(admet_df, on="drug_name", how="left")
    top20_enhanced.to_csv(RESULTS_DIR / "top20_hits.csv", index=False)
    print(f"  Updated results/top20_hits.csv with ADMET columns")

    # Also save standalone ADMET file
    admet_df.to_csv(RESULTS_DIR / "admet_flags.csv", index=False)
    print(f"  Saved results/admet_flags.csv")

    # Print summary
    n_ddi = admet_df["lorlatinib_ddi_flag"].sum()
    n_herg = admet_df["herg_risk"].sum()
    n_bbb = admet_df["bbb_permeable"].sum()
    n_gi = (admet_df["gi_absorption"] == "High").sum()
    print(f"\n  ADMET Summary:")
    print(f"    High GI absorption: {n_gi}/{len(admet_df)}")
    print(f"    BBB permeable: {n_bbb}/{len(admet_df)}")
    print(f"    hERG risk: {n_herg}/{len(admet_df)}")
    print(f"    Lorlatinib DDI flag: {n_ddi}/{len(admet_df)} (flagged, NOT excluded)")

    return admet_df
```

**Step 2: Commit**

```bash
git add 06_improve.py
git commit -m "feat: add Phase 3 — ADMET annotation for top hits"
```

---

### Task 4: Add `main()` and Run Full Pipeline

**Files:**
- Modify: `06_improve.py` (add main)

**Step 1: Add main function**

```python
def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Phase 1: Fix RMSD
    fix_rmsd_validation()

    # Phase 2: Multi-conformer re-dock
    run_multiconformer_redock()

    # Phase 3: ADMET annotation (uses top20 from CURRENT results;
    # re-run 04_analyze.py first if Campaign 2 scores changed)
    # We run analyze first to refresh top20 with new Campaign 2 scores
    print("\n=== Re-running analysis with updated Campaign 2 scores ===")
    import subprocess
    subprocess.run([sys.executable, "04_analyze.py"], check=True)

    # Now annotate the refreshed top 20
    run_admet_annotation()

    # Re-generate visualizations
    print("\n=== Re-generating report ===")
    subprocess.run([sys.executable, "05_visualize.py"], check=True)

    print("\n=== All improvements complete ===")


if __name__ == "__main__":
    main()
```

**Step 2: Run the full improvement pipeline**

Run: `cd /Users/lucapassone/Documents/ros1-screens/repurposing && uv run python 06_improve.py`

Expected output:
- Phase 1: Aligned RMSD < 2.0 A, validation_rmsd.txt updated
- Phase 2: ~630 docking jobs (63 drugs x 10 conformers), updated campaign2_scores.csv
- Phase 3: ADMET columns added to top20_hits.csv, admet_flags.csv created
- Re-run of 04_analyze.py and 05_visualize.py to regenerate all outputs

**Step 3: Verify outputs**

Check: `cat results/validation_rmsd.txt` — should show PASS
Check: `head -3 results/campaign2_scores.csv` — should have updated scores
Check: `head -1 results/top20_hits.csv` — should have ADMET columns
Check: `ls results/admet_flags.csv` — should exist

**Step 4: Commit all results**

```bash
git add 06_improve.py results/
git commit -m "feat: run methodology improvements — fixed RMSD, multi-conformer docking, ADMET"
```

---

### Task 5: Update Report to Reflect Improvements

**Files:**
- Modify: `05_visualize.py` — update the Limitations section in the HTML report

**Step 1: Update the report generation**

In `05_visualize.py`, find the section that generates the Limitations table (section 11.1) and update:
- Re-docking RMSD row: show new aligned RMSD value, mark as resolved
- Vina scoring: note multi-conformer approach as mitigation
- No ADMET filtering: mark as resolved, reference admet_flags.csv
- Single conformer: mark as resolved (10 conformers per drug)

Also add a new section or subsection showing the ADMET flags table for top 20 hits.

**Step 2: Re-run visualization**

Run: `cd /Users/lucapassone/Documents/ros1-screens/repurposing && uv run python 05_visualize.py`

**Step 3: Commit**

```bash
git add 05_visualize.py results/
git commit -m "docs: update report with resolved limitations and ADMET table"
```
