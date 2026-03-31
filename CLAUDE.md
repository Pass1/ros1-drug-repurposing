# ROS1 Drug Repurposing Screen

Screen ~2,500 FDA-approved drugs against ROS1 G2032R to find existing
medicines with unexpected ROS1 activity. If a hit is found, a doctor
could prescribe it off-label — no synthesis, no 10-year timeline.

## CLAUDE.md

You are working on a drug repurposing virtual screen. The goal is to
dock all FDA-approved small-molecule drugs against the ROS1 kinase
domain carrying the G2032R resistance mutation, and identify any with
predicted binding affinity comparable to known ROS1 inhibitors.

### What to do

1. **Fetch the ROS1 G2032R structure** from PDB 3ZBF, model the G2032R
   mutation, and prepare it for docking.
2. **Get the drug library** — download FDA-approved drug SMILES from
   DrugBank Open Data (CC0 license, free) or the ZINC20 FDA-approved
   subset (~2,000 drugs).
3. **Dock everything** against ROS1 G2032R using AutoDock Vina.
4. **Also dock known ROS1 TKIs** as positive controls (crizotinib,
   repotrectinib, zidesamtinib) so we know what "good" looks like.
5. **Rank results**. Flag anything that scores within 1-2 kcal/mol of
   the known TKIs — those are repurposing candidates.
6. **For the top 20 hits**, also dock against wild-type ROS1, TRK-A
   (to check selectivity), and run ADMET checks.

### Data sources (all free)

- **ROS1 structure**: PDB 3ZBF (https://www.rcsb.org/structure/3ZBF)
- **TRK-A structure**: PDB 4AOJ
- **DrugBank Open Data**: https://go.drugbank.com/releases/latest
  → Download "Approved Drug Structures" (SDF or SMILES, CC0 license)
- **ZINC FDA subset**: https://zinc20.docking.org/ → filter by
  "FDA approved" → download SMILES
- **ChEMBL kinase inhibitors**: Search for all approved kinase
  inhibitors as a focused sub-library

### Positive controls (known ROS1 binders — expect scores < -8 kcal/mol)

```
crizotinib:     Clc1cc(OC2(N)CC2)c(Cl)cc1-c1ccc2[nH]nc(-c3ccc(F)cc3F)c2n1
repotrectinib:  CC1(c2cc(-c3cncc4[nH]nc(N)c34)ccn2)C(F)(F)C1
lorlatinib:     O=C(Nc1cc(-c2cccc(F)c2)cn1C1CCC(NC(=O)C2(F)CC2)CC1)C1CC1
```

### What a "hit" looks like

A repurposing candidate should:
- Score < -8.5 kcal/mol against ROS1 G2032R (comparable to known TKIs)
- Be a small molecule (not a biologic)
- Not already be a known ROS1 inhibitor
- Ideally be CNS-penetrant (for brain mets)
- Have a tolerable safety profile at effective doses

### Key packages
```
pip install rdkit-pypi meeko vina biopython pdbfixer openmm numpy pandas
```

### Grid box
Center on the ATP binding pocket of chain A, residues K1980/E2027/D2102.
Box size: 25×25×25 Å. Exhaustiveness: 32 for positive controls, 8 for
the full library screen (speed vs accuracy tradeoff), then re-dock top
50 at exhaustiveness 32.

### Output
- `results/repurposing_ranked.csv` — all drugs ranked by binding score
- `results/top20_hits.csv` — best hits with full annotation
- `results/controls.csv` — known TKI scores for comparison
- `results/report.md` — summary of findings

### RunPod
Not strictly needed — 2,500 drugs × 1 target at exhaustiveness=8 runs
in a few hours on a decent laptop. Use RunPod if you want it done in
30 minutes (A40, ~$1-2).
