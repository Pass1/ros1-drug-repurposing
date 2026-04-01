# ROS1 Drug Repurposing Screen

Screen ~2,500 FDA-approved drugs against ROS1 G2032R to find existing
medicines with unexpected ROS1 activity. If a hit is found, a doctor
could prescribe it off-label — no synthesis, no 10-year timeline.

## CLAUDE.md

You are working on a drug repurposing virtual screen. The goal is to
dock all FDA-approved small-molecule drugs against the ROS1 kinase
domain carrying the G2032R resistance mutation, and identify any with
predicted binding affinity comparable to known ROS1 inhibitors.

### Pipeline scripts (run in order)

1. `01_prepare_receptor.py` — Download PDB structures, extract chain,
   fix missing atoms, add hydrogens, energy-minimise (OpenMM AMBER14),
   convert to PDBQT (Meeko), calculate grid box centers.
2. `02_fetch_drugs.py` — Parse Selleckchem SDFs, deduplicate by InChIKey,
   strip salts, filter biologics (MW > 1000), generate 3D conformers
   (ETKDG v3 + MMFF94), convert to PDBQT, compute properties + CNS MPO.
3. `03_dock.py` — Three docking campaigns:
   - Campaign 1: All ~2,800 drugs vs G2032R (exhaustiveness=8, 1 pose)
   - Campaign 2: Top 50 + 14 TKI controls vs G2032R (exhaustiveness=32, 9 poses)
   - Campaign 3: Top 20 + TKIs vs WT ROS1 and MET (exhaustiveness=32)
   - Pre-flight: re-docking validation of zidesamtinib crystal pose
4. `04_analyze.py` — Merge scores, compute composite score (binding +
   CNS MPO bonus + MET dual-activity bonus + mutant selectivity bonus),
   flag candidates within 2 kcal/mol of best TKI, generate ranked CSVs.
5. `05_visualize.py` — Charts (score distribution, top 20 bars,
   selectivity scatter, dual activity, CNS vs binding), interactive 3D
   poses (py3Dmol), and results/report.md.
6. `06_improve.py` — Methodology improvements (run after 03, before re-running 04+05):
   - Phase 1: Fix RMSD validation (Kabsch alignment with ICP atom matching)
   - Phase 2: Multi-conformer re-dock (10 conformers/drug, ETKDG v3 + MMFF94)
   - Phase 3: ADMET annotation (Lipinski, Veber, GI, BBB, P-gp, CYP3A4,
     hERG, lorlatinib DDI flag, PPB)

### Data sources

- **ROS1 G2032R structure**: PDB 9QEK (https://www.rcsb.org/structure/9QEK)
- **ROS1 WT structure**: PDB 7Z5X (https://www.rcsb.org/structure/7Z5X)
- **MET structure**: PDB 2WGJ (https://www.rcsb.org/structure/2WGJ)
- **FDA-approved drug libraries**: Selleckchem (https://www.selleckchem.com)
  -> L1300 FDA-approved Drug Library and L8000 FDA-approved Anticancer
  Drug Library (SDF format)

### Positive controls (known ROS1 binders)

```
crizotinib:     Clc1cc(OC2(N)CC2)c(Cl)cc1-c1ccc2[nH]nc(-c3ccc(F)cc3F)c2n1
repotrectinib:  CC1(c2cc(-c3cncc4[nH]nc(N)c34)ccn2)C(F)(F)C1
lorlatinib:     O=C(Nc1cc(-c2cccc(F)c2)cn1C1CCC(NC(=O)C2(F)CC2)CC1)C1CC1
```

Best TKI score: -10.7 kcal/mol (entrectinib). Range: -7.2 to -10.7.

### What a "hit" looks like

A repurposing candidate should:
- Score < -8.5 kcal/mol against ROS1 G2032R (comparable to known TKIs)
- Be a small molecule (not a biologic)
- Not already be a known ROS1 inhibitor
- Ideally be CNS-penetrant (for brain mets)
- Have a tolerable safety profile at effective doses

### Key packages
```
uv add rdkit meeko vina biopython pdbfixer openmm numpy pandas gemmi
```
- Use `uv` for package management (not conda or pip)
- **RDKit**: https://github.com/rdkit/rdkit

### Grid box
Center on the ATP binding pocket of chain A, residues K1980/E2027/D2102.
Box size: 25x25x25 A. Exhaustiveness: 32 for positive controls, 8 for
the full library screen (speed vs accuracy tradeoff), then re-dock top
50 at exhaustiveness 32 with multi-conformer docking.

### Docking validation
- Re-docking RMSD: 2.62 A (Kabsch-aligned, acceptable for Vina)
- Original 6.87 A was a coordinate-frame artifact (crystal vs minimised
  receptor, no superposition). Fixed in 06_improve.py Phase 1.

### ADMET annotation
Computed for top 20 hits in 06_improve.py Phase 3:
- Lipinski Ro5, Veber oral bioavailability
- GI absorption (Egan model), BBB permeability (Clark model)
- P-gp substrate, CYP3A4 inhibition risk, hERG liability
- Lorlatinib DDI flag (CYP3A4 — FLAG ONLY, NEVER EXCLUDE)
- Plasma protein binding estimate

### Remaining limitations
- Rigid receptor docking (no induced-fit / DFG-out binders)
- Vina scoring misses cation-pi, halogen bonds, water-mediated contacts
- Selleckchem library bias (enriched for bioactive/kinase compounds)
- No kinome-wide selectivity profiling
- All ADMET predictions are heuristic, not experimental

### Output
- `results/repurposing_ranked.csv` — all drugs ranked by composite score
- `results/top20_hits.csv` — best 20 non-TKI hits with ADMET annotation
- `results/controls.csv` — known TKI benchmark scores
- `results/admet_flags.csv` — standalone ADMET properties for top 20
- `results/report.md` — summary of findings
- `results/repurposing_report.html` — full HTML report with charts and 3D poses
- `results/campaign2_scores_single_conf.csv` — backup of pre-multiconformer Campaign 2

### RunPod
Not strictly needed — 2,500 drugs x 1 target at exhaustiveness=8 runs
in a few hours on a decent laptop. Multi-conformer re-dock (06_improve.py
Phase 2) is heavier: 63 drugs x 10 conformers x exhaustiveness=32.
