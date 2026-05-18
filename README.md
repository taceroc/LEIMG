# LEIMG
Simulate Light Echoes and return a raw image in fits format

*The content of this README.md was generated with the assistance of an AI tool.*

*The Python Package Skeleton Template was generated with the assistance of an AI tool.*

## What this does

`LEIMG/main.py` runs `SimulateLEInfPlane` for one or more parameter sets from a YAML file, then writes outputs to an output directory.

## Requirements

- Python 3.10+
- Recommended packages:
  - `numpy`
  - `scipy`
  - `pandas`
  - `matplotlib`
  - `astropy`
  - `pyyaml`

Install packages:

```bash
pip install numpy scipy pandas matplotlib astropy pyyaml
```

## How to run

From one directory from the project root (the simulation is then saved outside of the LEIMG repo):

```bash
python LEIMG/main.py SimulateLEInfPlane -file_to_parameters <path-to-input-yaml> -outdir <output-folder>
```

Optional flag:

- `--bool_save` / `--no-bool_save`: enable or disable saving output files (default: enabled).

Example:

```bash
python LEIMG/main.py SimulateLEInfPlane -file_to_parameters runs.yml -outdir results
```

## Input YAML format

The YAML file must be a mapping of run IDs to parameter dictionaries.

Example `runs.yml`:

```yaml
run_001:
  dt0: 100
  d: 1000
  dz0: 0.02
  ct: 100
  plane_coefficients: [1, 0, 1, 1]
  angles: [-180, 180]
  wave: 0.7499
  dust_env: "mw"
  composition: "both"

run_002:
  dt0: 120
  d: 1200
  dz0: 0.03
  ct: 110
  plane_coefficients: [1, 0, 1, 1.2]
  angles: [-180, 180]
  wave: 0.7499
  dust_env: "lmc"
  composition: "S"
```

Required keys per run:

- `dt0`
- `d`
- `dz0`
- `ct`
- `plane_coefficients` (must be 4 values)
- `angles` (must be 2 values)
- `wave`
- `dust_env` (`mw` or `lmc`)
- `composition` (`both`, `S`, or `C`)

## Outputs

For each run ID, files are written under:

```text
<outdir>/<run_id>/
```

And a combined manifest is written to:

```text
<outdir>/manifest.yml
```

