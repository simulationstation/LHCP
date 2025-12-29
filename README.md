# wdalpha

Production-grade pipeline for measuring fractional fine-structure constant variation (Δα/α) from white dwarf spectra using the many-multiplet framework. The initial configuration targets G191-B2B and is designed to scale to multiple stars.

## Requirements

- Python 3.11+
- Dependencies: numpy, scipy, pandas, matplotlib, astropy, astroquery, specutils (optional), lmfit or scipy.optimize (scipy used here), pydantic, PyYAML, rich, pytest

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Build atomic table

```bash
wdalpha build-atomic \
  --nistsource data/atomic/nist_g191.csv \
  --supp data/atomic/preval.csv \
  --out data/atomic/combined.csv
```

### Download HLSP coadd and run full pipeline

```bash
wdalpha run \
  --target G191-B2B \
  --data-root data \
  --atomic data/atomic/combined.csv \
  --out results
```

### Holistic WD run (E140H baseline)

```bash
wdalpha holistic-run \
  --target G191-B2B \
  --data-root ./LHCP/ \
  --out ./LHCP/results/holistic_run/ \
  --atomic ./LHCP/data/atomic/analysis_lines_clean.csv
```

### Piecewise workflow

```bash
wdalpha download --target G191-B2B --mode hlsp --out data/spectra/
wdalpha fit-lines \
  --spectrum data/spectra/<file>.fits \
  --atomic data/atomic/combined.csv \
  --out results/run_manual/line_measurements.csv
wdalpha infer \
  --lines results/run_manual/line_measurements.csv \
  --atomic data/atomic/combined.csv \
  --out results/run_manual/inferred_alpha.json
wdalpha report --run results/run_manual/
```

## Reproducibility

- Runs use fixed seeds (default: 42).
- Each run writes `config_used.yaml`, `line_measurements.csv`, `inferred_alpha.json`, and diagnostic plots under `results/run_<timestamp>/`.
- Assumptions:
  - Gravitational redshift/systemic velocity are absorbed into the global redshift nuisance parameter.
  - Laboratory wavelength uncertainties can be included as an additional variance term during inference.

## Configuration

See `default_config.yaml` for the default run configuration.
