# Benchmark Registry

This directory defines virtual-screening benchmarks for the repository.

Files:
- `benchmark_sets.csv`: one row per benchmark set
- `benchmark_compounds.csv`: labeled compounds/entities belonging to each set

`benchmark_sets.csv` columns:
- `benchmark_id`: stable identifier
- `name`: short display name
- `description`: what the benchmark measures
- `benchmark_mode`: `active_vs_background` or `labeled_only`
- `score_columns`: pipe-separated ranking columns from `results/repurposing_ranked.csv`
- `target`: target protein for the benchmark
- `mutation`: target state or mutation
- `limitation_note`: concise scientific caveat

`benchmark_compounds.csv` columns:
- `benchmark_id`: foreign key into `benchmark_sets.csv`
- `entity_id`: chemotype-level identifier used to collapse aliases/salts
- `drug_name`: exact name expected in the screened library
- `label`: `active`, `inactive`, or `decoy`
- `evidence_tier`: provenance strength
- `source_kind`: source category such as `literature`, `clinical`, `screen_background`
- `source_ref`: short citation or identifier
- `target`: assay target
- `mutation`: assay mutation or state
- `assay_type`: biochemical, cellular, clinical, etc.
- `activity_relation`: e.g. `=`, `<`, `<=`
- `activity_value`: optional numeric activity value
- `activity_unit`: e.g. `nM`, `uM`
- `notes`: free-text context

Current state:
- `ros1_known_actives_background` is a provisional recovery benchmark only.
- It is useful for tracking whether methodology changes improve recovery of known ROS1 actives.
- It is not a substitute for a literature-curated active/inactive or active/decoy benchmark.
