# Reproducible examples

Run all commands from the repository root after installing the dependencies.

## Air–hydrogel reference point

```bash
python scripts/analyze_mp_parameter_point.py \
  --Arho -0.9974 --Amu -0.999998 --AG -1 --Ek 0 --precision 80
```

`air_hydrogel_reference.json` records the output from the tested environment.
Compare the admissible roots and counts; last digits of residuals can vary
across environments. The acceptance calculation is at 80 decimal digits,
while exported root coordinates are binary64.

## Small Cartesian reference sweep

```bash
python scripts/run_mp_reference_sweep.py \
  --Arho -0.9974 --Amu -0.999998 --samples 3 --workers 1 \
  --precision 80 --output-dir results/reference_smoke
```

## Asymmetric float64 benchmark with targeted refinement

`benchmarks/verify_Arho_m0p95_Amu_m0p8.yaml` records this 31 by 31 study. The
CLI does not read YAML; these Bash commands construct the corresponding grid:

```bash
AG_VALUES=$(python -c 'import numpy as np; print(" ".join(f"{x:.17g}" for x in np.linspace(-0.8, 0.8, 31)))')
EK_VALUES=$(python -c 'import numpy as np; print(" ".join(f"{x:.17g}" for x in np.linspace(-0.8, 0.8, 31)))')
python scripts/run_parameter_sweep.py \
  --Arho -0.95 --Amu -0.8 --sweep AG $AG_VALUES --sweep Ek $EK_VALUES \
  --output results/verify_Arho_m0p95_Amu_m0p8/sweep_output.jsonl \
  --enable-high-precision-refinement
python scripts/summarize_parameter_sweep.py \
  results/verify_Arho_m0p95_Amu_m0p8/sweep_output.jsonl \
  --output-dir results/verify_Arho_m0p95_Amu_m0p8
python scripts/plot_parameter_sweep.py \
  results/verify_Arho_m0p95_Amu_m0p8 \
  --output-dir results/verify_Arho_m0p95_Amu_m0p8/figures --all --color-by Ek
```

This benchmark is separate from the tiny MP smoke test; it was not rerun as
part of standalone-folder preparation. It is included as a reproduction recipe.
