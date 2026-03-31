# ROS1 Drug Repurposing Screen — Design Document

**Date**: 2026-03-30
**Patient context**: 39yo, lung adenocarcinoma (lower lobe).
EZR::ROS1 fusion (exon 34), MET IHC 3+ 50%, PD-L1 TPS 95%, ARID2
pathogenic variant. High risk of CNS metastasis.

## Goal

Screen ~2,500 FDA-approved drugs against ROS1 G2032R to find existing
medicines with unexpected ROS1 activity. G2032R is the most common
resistance mutation after first-line ROS1 TKI therapy. If a hit is found,
a doctor could prescribe it off-label.

Secondary goals:
- Prioritize CNS-penetrant candidates (brain mets are common in ROS1+ NSCLC)
- Identify drugs with dual ROS1/MET activity (tumor has MET 3+)
- Profile selectivity: G2032R mutant vs wild-type ROS1

## Architecture

Five modular scripts with intermediate files. Any step can be re-run
independently. Checkpointing in the docking step so crashes don't lose work.

```
01_prepare_receptor.py  -->  data/receptor_*.pdbqt, data/grid_*.txt
02_fetch_drugs.py       -->  data/ligands/*.pdbqt, data/drug_library.csv
03_dock.py              -->  results/campaign*_scores.csv, results/poses/
04_analyze.py           -->  results/repurposing_ranked.csv, results/top20_hits.csv
05_visualize.py         -->  results/charts/, results/poses_3d/, results/report.md
```

## Step 1: Receptor Preparation

Three receptor structures:

| Target | PDB | Resolution | Rationale |
|--------|-----|-----------|-----------|
| ROS1 G2032R | 9QEK | 2.2 A | Only crystal structure with G2032R mutation present. Co-crystallized with zidesamtinib. Eliminates need for in-silico mutation modeling. |
| ROS1 WT | 7Z5X | 2.0 A | Highest resolution WT ROS1 kinase domain. For selectivity comparison. |
| MET | 2WGJ | 2.0 A | Crizotinib-bound MET kinase. Standard for Type I inhibitor screening. For dual-activity profiling given MET 3+. |

For each structure:
1. Download PDB, extract chain A
2. Remove co-crystallized ligand, water, ions
3. Add hydrogens at pH 7.4 (PDBFixer)
4. Constrained energy minimization (OpenMM, AMBER14 force field, backbone restrained)
5. Convert to PDBQT (Meeko `mk_prepare_receptor.py`)

Grid box per target: 25x25x25 A, centered on ATP-binding pocket catalytic
triad. For ROS1: K1980/E2027/D2102. For MET: K1110/M1160/D1222.
Config saved to `data/grid_<target>.txt`.

## Step 2: Drug Library Preparation

**Primary source**: DrugBank Open Data (SDF, CC0 license, ~2,500 approved drugs).
**Fallback**: ZINC20 FDA-approved subset.

Processing:
1. Download and parse structures
2. Deduplicate by InChIKey (merge sources, prefer DrugBank metadata)
3. Filter: remove biologics/peptides (MW > 1,000 Da), strip salts
4. Tag known ROS1 TKIs (crizotinib, entrectinib, repotrectinib,
   taletrectinib, lorlatinib, zidesamtinib, ceritinib) — kept in library
   as benchmarks, flagged `is_known_ros1_tki = True`
5. Generate 3D conformers (RDKit ETKDG + MMFF94 minimization)
6. Convert to PDBQT (Meeko)
7. Pre-compute properties: MW, TPSA, ClogP, HBD, HBA, rotatable bonds,
   Pfizer CNS MPO score (6-parameter desirability function, calculated in RDKit)

Output:
- `data/ligands/*.pdbqt` — one per drug
- `data/drug_library.csv` — drug name, SMILES, DrugBank ID, properties, CNS MPO, source, is_known_ros1_tki flag

## Step 3: Docking

Three sequential campaigns using AutoDock Vina on Mac Mini M4 (~10 CPU cores):

| Campaign | Target | Library | Exhaustiveness | Purpose |
|----------|--------|---------|---------------|---------|
| 1 | ROS1 G2032R | All ~2,500 drugs | 8 | Primary screen |
| 2 | ROS1 G2032R | Top 50 + all known TKIs | 32 | High-accuracy re-dock |
| 3 | ROS1 WT + MET | Top 20 + all known TKIs | 32 | Selectivity profiling |

**Re-docking validation**: Before the full screen, re-dock zidesamtinib
into 9QEK. Compare predicted pose to crystal pose (RMSD). Must be < 2.0 A
to validate the protocol. If not, adjust grid box.

**Checkpointing**: Append each result to `results/docking_progress.csv`
after completion. On restart, skip drugs already in that file.

**Pose saving**: Best pose for all drugs in campaign 1. All 9 poses for
campaigns 2 and 3.

**Estimated runtime**: ~4-6 hours for campaign 1. Campaigns 2-3 are small
and fast (<30 min total).

Output:
- `results/campaign1_scores.csv`
- `results/campaign2_scores.csv`
- `results/campaign3_scores.csv`
- `results/poses/` — organized by campaign
- `results/validation_rmsd.txt`

## Step 4: Analysis

1. Merge docking scores with drug_library.csv
2. Establish benchmark from known ROS1 TKI control scores
3. Flag repurposing candidates: score within 1-2 kcal/mol of best TKI
   AND not a known ROS1 inhibitor
4. Rank by composite score weighting:
   - G2032R binding score (primary)
   - CNS MPO >= 4.0 (bonus for brain met risk)
   - G2032R vs WT selectivity (mutant-selective = more interesting)
   - MET binding (dual ROS1/MET activity = bonus for MET 3+)
5. Annotate top 20 with: FDA indication, mechanism of action, typical
   dosing, safety profile, CNS penetration status

Output:
- `results/repurposing_ranked.csv` — all drugs ranked
- `results/top20_hits.csv` — best 20, fully annotated
- `results/controls.csv` — known TKI benchmark scores

## Step 5: Visualization

**Charts** (matplotlib/seaborn):
- Score distribution histogram with control TKI scores as vertical lines
- Top 20 bar chart, color-coded by CNS MPO
- G2032R vs WT selectivity scatter (top 20)
- ROS1 G2032R vs MET dual-activity scatter (top 20)
- CNS MPO vs binding score (all drugs, sweet-spot quadrant highlighted)

**3D binding poses** (py3Dmol, exported as interactive HTML + PNG):
- Top 20 hits in ROS1 G2032R binding pocket
- Key residues labeled (K1980, E2027, D2102, R2032)
- Hydrogen bonds as dashed lines
- Zidesamtinib reference pose for comparison
- Side-by-side: best control vs top hit

**Report** (`results/report.md`):
- Patient context
- Methods summary
- Re-docking validation result
- Control TKI scores table
- Top 20 hits table with full annotation
- Embedded charts
- Links to interactive 3D HTML files
- Limitations and next steps

## Positive Controls

| Drug | SMILES | Expected score | CNS penetration |
|------|--------|---------------|-----------------|
| Crizotinib | `Clc1cc(OC2(N)CC2)c(Cl)cc1-c1ccc2[nH]nc(-c3ccc(F)cc3F)c2n1` | < -8 kcal/mol | Poor |
| Repotrectinib | `CC1(c2cc(-c3cncc4[nH]nc(N)c34)ccn2)C(F)(F)C1` | < -8 kcal/mol | Good |
| Lorlatinib | `O=C(Nc1cc(-c2cccc(F)c2)cn1C1CCC(NC(=O)C2(F)CC2)CC1)C1CC1` | < -8 kcal/mol | Best |
| Zidesamtinib | (retrieve from DrugBank/PubChem) | < -9 kcal/mol | Good |

## What a "hit" looks like

A repurposing candidate should:
- Score < -8.5 kcal/mol against ROS1 G2032R (comparable to known TKIs)
- Be a small molecule, not a known ROS1 inhibitor
- Ideally CNS-penetrant (CNS MPO >= 4.0, TPSA < 79 A^2)
- Have a tolerable safety profile at effective doses
- Bonus: dual ROS1/MET activity

## Infrastructure

- **Execution**: Mac Mini M4 (local), ~10 CPU cores
- **Disk**: ~1.5-2 GB total
- **Runtime**: ~4-6 hours for primary screen, minutes for everything else
- **Key packages**: rdkit, meeko, vina, biopython, pdbfixer, openmm, numpy,
  pandas, matplotlib, seaborn, py3Dmol, tqdm, requests

## Data sources

- PDB 9QEK: https://www.rcsb.org/structure/9QEK (ROS1 G2032R + zidesamtinib)
- PDB 7Z5X: https://www.rcsb.org/structure/7Z5X (ROS1 WT, highest resolution)
- PDB 2WGJ: https://www.rcsb.org/structure/2WGJ (MET + crizotinib)
- DrugBank Open Data: https://go.drugbank.com/releases/latest
- ZINC20 FDA subset: https://zinc20.docking.org/
