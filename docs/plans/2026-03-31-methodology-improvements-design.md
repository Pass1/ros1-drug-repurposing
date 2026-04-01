# Methodology Improvements Design

Date: 2026-03-31

## Goal

Address three methodological limitations from the repurposing screen report without re-running the full docking campaign.

## A) Fix Re-docking RMSD Validation

**Root cause**: The 6.87 A RMSD is a measurement bug. Crystal ligand coords are from original 9QEK; docked pose is in the minimized receptor frame. RMSD computed on raw coordinates without superposition.

**Fix**: Use RDKit `rdMolAlign.GetBestRMS()` for symmetry-aware, alignment-based RMSD. Parse docked PDBQT back into RDKit mol by injecting coordinates into original mol template. No re-docking needed.

## B) ADMET Annotation for Top Hits

Compute for top 20 hits using RDKit descriptors and heuristic models:

- Lipinski / Veber rule violations
- GI absorption (Egan model: TPSA + cLogP)
- BBB permeability (already have CNS MPO)
- P-gp substrate likelihood (MW > 400 + TPSA > 90)
- CYP3A4 inhibition flag (MW, cLogP, rotatable bonds heuristic)
- hERG liability (cLogP > 3.7 + basic amine)
- Lorlatinib DDI flag (CYP3A4 involvement) -- FLAG ONLY, NEVER EXCLUDE
- Plasma protein binding estimate

Output: expanded top20_hits.csv + results/admet_flags.csv.

## C) Multi-conformer Docking for Top Hits

For top 50 from Campaign 1:

1. Generate 10 conformers per drug (ETKDG v3 + MMFF94 minimization)
2. Convert each to PDBQT via Meeko
3. Dock all conformers vs G2032R at exhaustiveness=32
4. Keep best score across conformers -> new Campaign 2

Replaces existing single-conformer Campaign 2.

## Implementation

Single new script `06_improve.py` with three phases:
1. Fix RMSD validation
2. Multi-conformer re-dock top 50
3. ADMET annotation top 20

Then re-run 04_analyze.py and 05_visualize.py to regenerate outputs.

## What does NOT change

- Campaign 1 scores (full library screen) stay as-is
- Campaign 3 selectivity profiling stays as-is
- Receptor preparation stays as-is
- Drug library stays as-is
