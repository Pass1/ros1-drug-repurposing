# ROS1 Drug Repurposing Screen — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Screen ~2,500 FDA-approved drugs against ROS1 G2032R (and WT/MET) to find repurposing candidates for a patient with EZR::ROS1 fusion NSCLC and MET 3+ overexpression.

**Architecture:** Five modular Python scripts with intermediate files and checkpointing. Each script reads from `data/` and writes to `results/`. Any step can re-run independently.

**Tech Stack:** Python 3.11+, uv, rdkit, meeko, vina, biopython, pdbfixer, openmm, numpy, pandas, matplotlib, seaborn, py3Dmol, tqdm, requests

**Design doc:** `docs/plans/2026-03-30-repurposing-screen-design.md`

---

### Task 0: Project Setup

**Files:**
- Create: `pyproject.toml`

**Step 1: Initialize uv project**

```bash
cd /Users/lucapassone/Documents/ros1-screens/repurposing
uv init --python 3.11
```

**Step 2: Add dependencies**

```bash
uv add rdkit-pypi meeko vina biopython pdbfixer openmm numpy pandas matplotlib seaborn py3Dmol tqdm requests
```

**Step 3: Create directory structure**

```bash
mkdir -p data/ligands data/controls results/poses results/charts results/poses_3d
```

**Step 4: Initialize git**

```bash
git init
```

Create `.gitignore`:

```
.venv/
data/ligands/*.pdbqt
data/*.pdb
data/*.pdbqt
results/poses/*.pdbqt
__pycache__/
*.pyc
```

**Step 5: Verify environment**

```bash
uv run python -c "from rdkit import Chem; from vina import Vina; from pdbfixer import PDBFixer; print('All imports OK')"
```

Expected: `All imports OK`

**Step 6: Commit**

```bash
git add -A
git commit -m "chore: initialize project with uv and dependencies"
```

---

### Task 1: Receptor Preparation Script

**Files:**
- Create: `01_prepare_receptor.py`

**Step 1: Write the script**

`01_prepare_receptor.py` must:

1. Download PDB files (9QEK, 7Z5X, 2WGJ) from RCSB using `requests`
2. For each structure:
   - Load with PDBFixer
   - Extract chain A only
   - Remove heterogens (ligands, water, ions) but keep any crystallographic metals in the binding site
   - Find and add missing residues/atoms
   - Add hydrogens at pH 7.4
   - Energy minimize with OpenMM (AMBER14 force field, backbone heavy atoms restrained with 10 kcal/mol/A^2 harmonic restraint, 1500 steps L-BFGS)
   - Save cleaned PDB to `data/receptor_<name>.pdb`
3. Convert each cleaned PDB to PDBQT using Meeko's `mk_prepare_receptor` or `MoleculePreparation`
4. Calculate grid box center for each target:
   - ROS1 (both): average CA coordinates of K1980, E2027, D2102
   - MET: average CA coordinates of K1110, M1160, D1222
5. Write grid config files `data/grid_<target>.txt` with format:
   ```
   center_x = X.XXX
   center_y = Y.YYY
   center_z = Z.ZZZ
   size_x = 25.0
   size_y = 25.0
   size_z = 25.0
   ```
6. Extract the co-crystallized zidesamtinib from 9QEK before cleanup, save as `data/controls/zidesamtinib_crystal.pdb` for later RMSD validation

**Step 2: Run it**

```bash
uv run python 01_prepare_receptor.py
```

Expected output files:
- `data/receptor_G2032R.pdb` and `data/receptor_G2032R.pdbqt`
- `data/receptor_WT.pdb` and `data/receptor_WT.pdbqt`
- `data/receptor_MET.pdb` and `data/receptor_MET.pdbqt`
- `data/grid_G2032R.txt`, `data/grid_WT.txt`, `data/grid_MET.txt`
- `data/controls/zidesamtinib_crystal.pdb`

**Step 3: Verify**

```bash
# Check all files exist and are non-empty
ls -la data/receptor_*.pdbqt data/grid_*.txt data/controls/zidesamtinib_crystal.pdb
```

**Step 4: Commit**

```bash
git add 01_prepare_receptor.py data/grid_*.txt
git commit -m "feat: add receptor preparation script (9QEK, 7Z5X, 2WGJ)"
```

---

### Task 2: Drug Library Preparation Script

**Files:**
- Create: `02_fetch_drugs.py`

**Step 1: Write the script**

`02_fetch_drugs.py` must:

1. **Download DrugBank Open Data** — the approved drug structures SDF from DrugBank's open data release. If download requires manual step (account), print instructions and check if `data/drugbank_approved.sdf` already exists.
2. **Fallback to ZINC20** — if DrugBank SDF not found, download FDA-approved SMILES from ZINC20.
3. **Parse and deduplicate**:
   - Read all molecules with RDKit
   - Compute InChIKey for each
   - Deduplicate by InChIKey
   - Log how many drugs from each source
4. **Filter**:
   - Remove MW > 1000 Da (biologics/peptides)
   - Strip salts (`Chem.SaltRemover`)
   - Remove molecules that fail sanitization
5. **Tag known ROS1 TKIs** by name matching: crizotinib, entrectinib, repotrectinib, taletrectinib, lorlatinib, zidesamtinib, ceritinib. Set `is_known_ros1_tki = True`.
6. **Add positive control SMILES** if not found in library — hardcode the SMILES from the design doc for crizotinib, repotrectinib, lorlatinib. Look up zidesamtinib SMILES from PubChem (CID 137321816).
7. **Generate 3D conformers**:
   - `AllChem.EmbedMolecule` with ETKDG v3
   - `AllChem.MMFFOptimizeMolecule` (MMFF94 force field)
   - Skip molecules that fail embedding (log them)
8. **Convert to PDBQT** using Meeko for each molecule. Save to `data/ligands/<drugbank_id_or_name>.pdbqt`. Save controls separately to `data/controls/`.
9. **Compute properties** for every drug:
   - MW, TPSA, ClogP (Crippen), HBD, HBA, rotatable bonds
   - pKa of most basic center (estimate via RDKit or set to NaN if unavailable)
   - Pfizer CNS MPO score: sum of 6 desirability functions (0-1 each) for MW, ClogP, ClogD, TPSA, HBD, pKa. Score 0-6, >= 4.0 = CNS-favorable.
   - Simplified CNS MPO (without ClogD/pKa if unavailable): flag as CNS-penetrant if TPSA < 79 and HBD <= 2 and MW < 450
10. **Save** `data/drug_library.csv` with columns: drug_name, drugbank_id, smiles, inchikey, mw, tpsa, clogp, hbd, hba, rotatable_bonds, cns_mpo, cns_penetrant, is_known_ros1_tki, pdbqt_path, source

**Step 2: Run it**

```bash
uv run python 02_fetch_drugs.py
```

Expected: ~2,000-2,500 drugs processed, CSV written, PDBQT files in `data/ligands/`.

**Step 3: Verify**

```bash
uv run python -c "import pandas as pd; df = pd.read_csv('data/drug_library.csv'); print(f'{len(df)} drugs, {df.is_known_ros1_tki.sum()} controls, {df.cns_penetrant.sum()} CNS-penetrant')"
```

Expected: ~2000+ drugs, 7 controls, several hundred CNS-penetrant.

**Step 4: Commit**

```bash
git add 02_fetch_drugs.py
git commit -m "feat: add drug library fetch and preparation script"
```

---

### Task 3: Docking Script

**Files:**
- Create: `03_dock.py`

**Step 1: Write the script**

`03_dock.py` must:

1. **Re-docking validation** (runs first):
   - Prepare zidesamtinib from crystal pose (`data/controls/zidesamtinib_crystal.pdb`) as PDBQT
   - Dock into G2032R receptor with exhaustiveness=32
   - Compute RMSD between predicted and crystal pose (heavy atoms)
   - Write RMSD to `results/validation_rmsd.txt`
   - Print warning if RMSD > 2.0 A but continue anyway (log the warning)

2. **Campaign 1 — Primary screen**:
   - Load G2032R receptor PDBQT and grid config
   - Read `data/drug_library.csv` to get list of all drugs
   - For each drug:
     - Check if already in `results/docking_progress.csv` (checkpoint) — skip if so
     - Load ligand PDBQT from `data/ligands/` or `data/controls/`
     - Dock with Vina: exhaustiveness=8, n_poses=1, cpu=N (auto-detect)
     - Save best pose to `results/poses/campaign1/<drug_name>.pdbqt`
     - Append row to `results/docking_progress.csv`: drug_name, score, pose_path, timestamp
   - After all done, copy checkpoint to `results/campaign1_scores.csv`

3. **Campaign 2 — Re-dock top 50**:
   - Read campaign 1 scores, take top 50 by score + all known TKIs
   - Dock each against G2032R with exhaustiveness=32, n_poses=9
   - Save all poses to `results/poses/campaign2/`
   - Write `results/campaign2_scores.csv`

4. **Campaign 3 — Selectivity profiling**:
   - Take top 20 from campaign 2 + all known TKIs
   - Dock each against WT receptor (7Z5X grid) with exhaustiveness=32, n_poses=9
   - Dock each against MET receptor (2WGJ grid) with exhaustiveness=32, n_poses=9
   - Save poses to `results/poses/campaign3_wt/` and `results/poses/campaign3_met/`
   - Write `results/campaign3_scores.csv` with columns: drug_name, g2032r_score, wt_score, met_score

5. **Progress reporting**: Print progress every 50 drugs with ETA using tqdm.

6. **Parallelism**: Use `multiprocessing.Pool` with `os.cpu_count()` workers. Each worker runs one Vina docking. The checkpoint file must be written with a lock to avoid corruption.

**Step 2: Run it**

```bash
uv run python 03_dock.py
```

Expected: ~4-6 hours for campaign 1, <30 min for campaigns 2-3. Progress bar visible.

**Step 3: Verify**

```bash
uv run python -c "
import pandas as pd
c1 = pd.read_csv('results/campaign1_scores.csv')
c2 = pd.read_csv('results/campaign2_scores.csv')
c3 = pd.read_csv('results/campaign3_scores.csv')
print(f'Campaign 1: {len(c1)} drugs docked')
print(f'Campaign 2: {len(c2)} drugs re-docked')
print(f'Campaign 3: {len(c3)} drugs profiled')
print(f'Best score: {c1.score.min():.1f} kcal/mol')
"
```

**Step 4: Commit**

```bash
git add 03_dock.py
git commit -m "feat: add docking script with 3 campaigns and checkpointing"
```

---

### Task 4: Analysis Script

**Files:**
- Create: `04_analyze.py`

**Step 1: Write the script**

`04_analyze.py` must:

1. **Load and merge**:
   - `data/drug_library.csv` (properties, CNS MPO, TKI flag)
   - `results/campaign2_scores.csv` (high-accuracy G2032R scores for top 50)
   - `results/campaign1_scores.csv` (all drugs, lower accuracy)
   - `results/campaign3_scores.csv` (WT and MET scores for top 20)
   - Use campaign 2 score when available, fall back to campaign 1

2. **Establish benchmark**:
   - Filter known TKIs, report their scores
   - Best TKI score = reference threshold
   - Write `results/controls.csv` with: drug_name, g2032r_score, wt_score, met_score, cns_mpo

3. **Flag candidates**:
   - `is_candidate = True` if: score within 2 kcal/mol of best TKI AND `is_known_ros1_tki == False`

4. **Compute composite score** for ranking:
   - Start with raw G2032R binding score (more negative = better)
   - Bonus -0.5 kcal/mol if CNS MPO >= 4.0
   - Bonus -0.3 kcal/mol if binds MET within 2 kcal/mol of MET controls
   - Bonus -0.2 kcal/mol if G2032R score is better than WT score (mutant-selective)

5. **Rank all drugs** by composite score, write `results/repurposing_ranked.csv`

6. **Annotate top 20**:
   - All scores (G2032R, WT, MET)
   - All properties (MW, TPSA, ClogP, CNS MPO)
   - DrugBank metadata: indication, mechanism, safety notes
   - Write `results/top20_hits.csv`

7. **Print summary** to console: top 20 table, control scores, number of candidates found.

**Step 2: Run it**

```bash
uv run python 04_analyze.py
```

**Step 3: Verify**

```bash
head -25 results/top20_hits.csv
cat results/controls.csv
```

**Step 4: Commit**

```bash
git add 04_analyze.py
git commit -m "feat: add analysis script with composite scoring and annotation"
```

---

### Task 5: Visualization Script

**Files:**
- Create: `05_visualize.py`

**Step 1: Write the script**

`05_visualize.py` must:

1. **Score distribution histogram** (`results/charts/score_distribution.png`):
   - Histogram of all drug scores from campaign 1
   - Vertical dashed lines for each control TKI (labeled)
   - X-axis: binding score (kcal/mol), Y-axis: count

2. **Top 20 bar chart** (`results/charts/top20_scores.png`):
   - Horizontal bars for top 20 hits + controls
   - Color-coded: green if CNS MPO >= 4.0, orange if 3.0-4.0, red if < 3.0
   - Known TKIs in a different shade/pattern for distinction

3. **Selectivity scatter** (`results/charts/selectivity_g2032r_vs_wt.png`):
   - X: WT score, Y: G2032R score for top 20 + controls
   - Diagonal line (equal selectivity)
   - Points below the line = prefer mutant
   - Label each point with drug name

4. **Dual-activity scatter** (`results/charts/dual_ros1_met.png`):
   - X: MET score, Y: G2032R score
   - Highlight "sweet spot" quadrant (strong on both)
   - Label each point

5. **CNS MPO vs binding** (`results/charts/cns_vs_binding.png`):
   - All drugs as light dots
   - Top 20 highlighted and labeled
   - Quadrant lines at score=-8.5 and CNS MPO=4.0
   - "Sweet spot" quadrant highlighted

6. **3D binding poses** (py3Dmol → HTML + PNG):
   - For each top 20 hit: load receptor PDB + docked ligand pose
   - Create py3Dmol view:
     - Protein as cartoon (light gray) with surface near binding site
     - Key residues (K1980, E2027, D2102, R2032) as labeled sticks
     - Ligand as colored sticks
     - H-bonds as dashed yellow lines (distance-based detection: donor-acceptor < 3.5 A)
   - Save as `results/poses_3d/<drug_name>.html` (interactive)
   - Also render a static PNG snapshot for the report
   - Create reference view with zidesamtinib for comparison

7. **Generate report** (`results/report.md`):
   - Title, date, patient context (EZR::ROS1, MET 3+, brain met risk)
   - Methods: structures used, grid box, exhaustiveness, library size
   - Re-docking validation: zidesamtinib RMSD result
   - Controls table
   - Top 20 hits table (embedded from CSV)
   - Charts embedded as `![](charts/filename.png)`
   - Links to interactive 3D HTML files
   - Limitations section (docking ≠ binding, scoring function approximations, need experimental validation)
   - Next steps (discuss with oncologist, CRO testing at Eurofins, ROS1ders network)

**Step 2: Run it**

```bash
uv run python 05_visualize.py
```

**Step 3: Verify**

```bash
ls results/charts/*.png results/poses_3d/*.html results/report.md
open results/report.md
```

**Step 4: Commit**

```bash
git add 05_visualize.py
git commit -m "feat: add visualization script with charts, 3D poses, and report"
```

---

### Task 6: Final Review and Cleanup

**Step 1: Run the full pipeline end-to-end**

```bash
uv run python 01_prepare_receptor.py
uv run python 02_fetch_drugs.py
uv run python 03_dock.py        # ~4-6 hours
uv run python 04_analyze.py
uv run python 05_visualize.py
```

**Step 2: Review outputs**

- Open `results/report.md` — does it make sense?
- Open a few `results/poses_3d/*.html` — do the poses look reasonable?
- Check `results/controls.csv` — do TKI scores match expectations (< -8 kcal/mol)?
- Check `results/validation_rmsd.txt` — is RMSD < 2.0 A?

**Step 3: Final commit**

```bash
git add -A
git commit -m "feat: complete ROS1 drug repurposing screen pipeline"
```
