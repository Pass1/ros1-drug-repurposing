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
   normalize/strip salts before deduplication, filter biologics (MW > 1000), generate 3D conformers
   (ETKDG v3 + MMFF94), convert to PDBQT, compute properties + CNS MPO.
3. `03_dock.py` — Three docking campaigns:
   - Campaign 1: All ~2,800 drugs vs G2032R (exhaustiveness=8, 1 pose)
   - Campaign 2: Top 50 + 14 TKI controls vs G2032R (exhaustiveness=32, 9 poses)
   - Campaign 3: Top 20 + TKIs vs WT ROS1 and MET (exhaustiveness=32)
   - Pre-flight: preliminary re-docking of zidesamtinib crystal pose
     (score only; detailed RMSD is deferred to `06_improve.py`)
4. `04_analyze.py` — Merge scores, compute composite score (binding +
   CNS MPO bonus + MET dual-activity bonus + mutant selectivity bonus),
   flag candidates within 2 kcal/mol of best TKI, generate ranked CSVs.
5. `05_visualize.py` — Charts (score distribution, top 20 bars,
   selectivity scatter, dual activity, CNS vs binding), interactive 3D
   poses (py3Dmol), and results/report.md.
6. `06_improve.py` — Methodology improvements (run after 03, before re-running 04+05):
   - Phase 1: Fix RMSD validation (RDKit GetBestRMS when atom counts match;
     otherwise element-constrained Hungarian matching + Kabsch alignment)
   - Phase 2: Multi-conformer re-dock (10 conformers/drug, ETKDG v3 + MMFF94)
   - Phase 3: ADMET annotation (Lipinski, Veber, GI, BBB, P-gp, CYP3A4,
     hERG, lorlatinib DDI flag, PPB)
7. `07_benchmark.py` — Benchmark / enrichment evaluation:
   - Load benchmark sets from `data/benchmarks/benchmark_sets.csv`
     and compound labels from `data/benchmarks/benchmark_compounds.csv`
   - Support both `active_vs_background` and `labeled_only` benchmark modes
   - Collapse aliases/salts to chemotype-level entities
   - Evaluate ROC AUC, average precision, EF/nEF at 1/5/10%, top-N recovery,
     cumulative enrichment, and score-calibration tables
   - Treat the current ROS1 benchmark as provisional until a stronger
     literature-curated active/inactive or active/decoy panel is added

### Data sources

- **ROS1 G2032R structure**: PDB 9QEK (https://www.rcsb.org/structure/9QEK)
- **ROS1 WT structure**: PDB 7Z5X (https://www.rcsb.org/structure/7Z5X)
- **MET structure**: PDB 2WGJ (https://www.rcsb.org/structure/2WGJ)
- **FDA-approved drug libraries**: Selleckchem (https://www.selleckchem.com)
  -> L1300 FDA-approved Drug Library and L8000 FDA-approved Anticancer
  Drug Library (SDF format)
- **Benchmark registry**:
  - `data/benchmarks/benchmark_sets.csv` — benchmark-set metadata
  - `data/benchmarks/benchmark_compounds.csv` — labeled compounds/entities
  - `data/benchmarks/README.md` — schema and curation notes

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

### Benchmark workflow
- `07_benchmark.py` labels and evaluates compounds already present in
  `results/repurposing_ranked.csv`; it does **not** modify the underlying
  screened drug lists or replace `data/drug_library.csv`.
- Current seed benchmark:
  - `ros1_known_actives_background`
  - Mode: `active_vs_background`
  - Positives: known ROS1 inhibitors present in the screened library
    (collapsed to chemotype/entity level)
  - Negatives: the remaining scored library as presumed background
- This is useful for tracking early-recognition performance of methodology
  changes, but it remains weaker than a curated active/inactive or
  active/decoy benchmark.
- The benchmark infrastructure is now ready for additional sets, including:
  - literature-curated ROS1 WT or G2032R active/inactive panels
  - justified ROS1 decoy sets
  - target-specific calibration subsets for alternate ranking strategies

### Current benchmark outputs
- `results/benchmark_registry.csv` — benchmark-set summary and matched counts
- `results/benchmark_membership.csv` — row-level labels joined to ranked outputs
- `results/benchmark_metrics.csv` — enrichment / recovery metrics by benchmark and score column
- `results/benchmark_calibration.csv` — active-rate tables across score bins
- `results/benchmark_active_ranks.csv` — active ranks per benchmark / score column
- `results/benchmark_summary.md` — human-readable benchmark summary
- `results/charts/benchmark_enrichment.png` — cumulative active-recovery plot
- `results/charts/benchmark_calibration.png` — score-calibration plot

### Current benchmark snapshot
- Current benchmark set: `ros1_known_actives_background`
- Current matched positive rows: 13
- Current unique active entities after alias/salt collapse: 7
- Current metrics on existing ranked outputs:
  - `g2032r_score`: ROC AUC 0.935, AP 0.064, EF1 14.3
  - `composite_score`: ROC AUC 0.968, AP 0.093, EF1 28.5
- Interpretation:
  - The workflow retrieves known ROS1 actives better than random.
  - `composite_score` currently outperforms raw `g2032r_score`.
  - Benchmark strength is still limited by the lack of curated negatives.

### Docking validation
- `03_dock.py` intentionally does **not** report RMSD anymore; it writes a
  preliminary note and docking score only.
- Current `06_improve.py` validation output: geometry-aware aligned RMSD
  3.08 A, versus 7.76 A naive unaligned RMSD.
- Original 6.87 A was a coordinate-frame artifact (crystal vs minimised
  receptor, no superposition).
- When atom counts differ between crystal and docked ligands, the RMSD is
  still an approximate geometry-aware estimate, not a chemistry-exact atom map.

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
- Current benchmark layer uses known ROS1 actives recovered from the
  screened library against presumed-background decoys; this is useful for
  early-recognition tracking but is weaker than a literature-curated
  actives/inactives or actives/decoys benchmark
- No curated ROS1 inactive/decoy panel has been added yet
- No mutant-specific benchmark labels have been curated yet for G2032R vs WT
- Benchmark infrastructure is now ready for stronger labeled panels, but
  the repository still needs those labels curated from literature/ChEMBL
  before the next full scientific pass
- Existing ranked outputs may still reflect pre-remediation library
  normalization until the full `02 -> 03 -> 04` pipeline is rerun

### Recent remediation
- Salt stripping now happens before deduplication in `02_fetch_drugs.py`,
  so salt forms should collapse before ranking on the next full rerun.
- Known ROS1 TKI tagging now uses either name/synonym matching or `target`
  fields mentioning ROS1.
- `05_visualize.py` now removes stale 3D pose pages and writes state-aware
  report text instead of claiming methodology improvements unconditionally.
- `07_benchmark.py` now provides a registry-driven enrichment layer so
  methodology changes can be judged on recovery of known ROS1 actives,
  not only on docking score anecdotes.
- The benchmark layer now emits registry, calibration, active-rank, and
  benchmark-specific entity outputs suitable for comparing future reruns.

### Output
- `results/repurposing_ranked.csv` — all drugs ranked by composite score
- `results/top20_hits.csv` — best 20 non-TKI hits with ADMET annotation
- `results/controls.csv` — known TKI benchmark scores
- `results/admet_flags.csv` — standalone ADMET properties for top 20
- `results/benchmark_membership.csv` — row-level benchmark labels joined to ranked outputs
- `results/benchmark_registry.csv` — benchmark-set summary and matched counts
- `results/benchmark_metrics.csv` — enrichment / recovery metrics for each ranking column
- `results/benchmark_calibration.csv` — active-rate tables across score bins
- `results/benchmark_active_ranks.csv` — active ranks per benchmark / score column
- `results/benchmark_summary.md` — human-readable benchmark summary with active ranks
- `results/benchmark_entities_<benchmark_id>_<score_column>.csv` — entity-level ranked tables used for metric calculations
- `results/report.md` — summary of findings
- `results/repurposing_report.html` — full HTML report with charts and 3D poses
- `results/campaign2_scores_single_conf.csv` — backup of pre-multiconformer Campaign 2

### RunPod
Not strictly needed — 2,500 drugs x 1 target at exhaustiveness=8 runs
in a few hours on a decent laptop. Multi-conformer re-dock (06_improve.py
Phase 2) is heavier: 63 drugs x 10 conformers x exhaustiveness=32.

### Validation commands
Use `uv` for validation in this sandbox. The current working pattern is:
```bash
UV_CACHE_DIR=/tmp/uv-cache MPLCONFIGDIR=/tmp/mpl uv run python 07_benchmark.py
UV_CACHE_DIR=/tmp/uv-cache MPLCONFIGDIR=/tmp/mpl uv run python 05_visualize.py
UV_CACHE_DIR=/tmp/uv-cache uv run python -m py_compile 07_benchmark.py 05_visualize.py
```
