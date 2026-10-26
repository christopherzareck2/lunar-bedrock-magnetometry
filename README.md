# lunar-bedrock-magnetometry

Code and data for Narrett, Zareck, Weiss & Hoffman, *Lunar Rover Magnetometry: Recovering Paleofield Intensity and
Direction From Magnetized Bedrock* (The Planetary Science Journal, in preparation).

We simulate a rover-mounted magnetometer surveying magnetized lunar bedrock and invert the synthetic data for the
natural remanent magnetization of the outcrop using [SimPEG](https://simpeg.xyz).

## Installation

```
git clone https://github.com/christopherzareck2/lunar-bedrock-magnetometry.git
cd lunar-bedrock-magnetometry
conda env create -f environment.yml
conda activate lunar-mag
```

## Usage

Open `notebook/reproduce_paper.ipynb` and run all cells. Part 1 builds and inverts the isolated-outcrop model.
Part 2 redraws Figures 1-8 and recomputes the numbers quoted in the paper from the stored results in `results/`.

To regenerate the stored results, set `RERUN = True` in the last section of the notebook (about 3.5 hours on a
laptop). All random draws are seeded.

`python tests/reproduce_check.py` reruns parts of the analysis from scratch (about a minute) and compares them with
the stored results.

## Contents

| Folder | Contents |
|---|---|
| `engine/` | forward model, inversion and noise budget |
| `analysis/` | scripts that produce the results in the paper |
| `figures/` | plotting scripts and figures |
| `results/` | stored results |
| `notebook/` | reproduction notebook |
| `tests/` | reproduction check |

| Paper | Script |
|---|---|
| Table 1 | `engine/noise_budget.py` |
| Figures 1-3 | `figures/scripts/fig01_*.py` to `fig03_*.py` |
| Table 2 | `engine/inversion_approach_comparison.py` |
| Isolated outcrop, Figures 4 and 5 | `analysis/a02_*.py`, `analysis/a04_text_number_checks.py` |
| Table 3, Figure 8 | `analysis/a06_isolated_robustness.py` |
| Skate park, Figure 6 | `engine/skatepark_1cm_dr_sweep.py`, `analysis/a08_skatepark_neighbor_directions.py` |
| Three-layer wall, Table 4, Figure 7 | `engine/three_layer_wall.py`, `engine/scripts/elevated_traverse_height_sweep.py`, `analysis/a03_*.py`, `analysis/a07_*.py` |
| Background magnetization and correlated noise | `analysis/a05_background_and_correlated_noise.py` |
| Systematic errors and error budget, Table 5 | `analysis/a09_systematics.py` |

Results were produced with Python 3.10 and the package versions in `environment.yml`.

## Citation

Please cite the paper if you use this code:

> Narrett, I. S., Zareck, C. A., Weiss, B. P., & Hoffman, J. A. (in preparation). Lunar Rover Magnetometry:
> Recovering Paleofield Intensity and Direction From Magnetized Bedrock. *The Planetary Science Journal*.

See `CITATION.cff` for the software citation.

## License

MIT (see `LICENSE`).
