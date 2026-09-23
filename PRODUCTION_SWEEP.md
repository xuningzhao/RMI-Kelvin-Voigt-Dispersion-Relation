# Manifest-driven production sweeps

This optional workflow uses the float64 solver with configurable targeted
high-precision refinement. It is distinct from the full MP reference sweep in
README.md. The full production configuration is computationally expensive.
Use the small reference example for an initial installation check.

`config/production_sweep_v1.json` describes Group A (343 cases with 501 Ek
points) and three Group B families (49 maps per family, 501 by 501 points per
map). Group A uses atomic case payloads; Group B uses checksummed tiles.

## Build manifests

```bash
python scripts/build_production_manifests.py \
  --config config/production_sweep_v1.json \
  --campaign-root results/production_v1
```

This writes manifests under `results/production_v1/manifests/` and records
configuration and environment provenance. It does not launch solver jobs.

Use `python scripts/run_production_bundle.py --help` to select a manifest,
case-worker count, cases per worker, and optional work limit. The number of
manifest records must equal workers times cases per worker. Before running a
full campaign, select resources and a configuration appropriate to your
machine. Cluster-specific submission scripts are not part of this release.

## Validate completed output

```bash
python scripts/validate_production_results.py \
  --config config/production_sweep_v1.json \
  --manifest results/production_v1/manifests/a.jsonl \
  --campaign-root results/production_v1
python scripts/update_production_index.py results/production_v1
python scripts/build_retry_manifest.py \
  --manifest results/production_v1/manifests/a.jsonl \
  --campaign-root results/production_v1 \
  --output results/production_v1/manifests/a_retry.jsonl
```

Validation fails while selected cases are incomplete. Restart skips only
payloads whose configuration hash, checksum, schema, and internal array
structure validate. Failed or damaged output units are recomputed. Keep the
original configuration and manifest with each campaign.
