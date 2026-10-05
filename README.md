# Lunar rover magnetometry of magnetized bedrock

## Summary

This repository reproduces Narrett, Zareck, Weiss & Hoffman, *Lunar Rover Magnetometry: Recovering Paleofield
Intensity and Direction From Magnetized Bedrock* (Planetary Science Journal, in preparation). We simulate a rover
magnetometer surveying magnetized lunar bedrock and invert the synthetic data for the natural remanent
magnetization of the outcrop with [SimPEG](https://simpeg.xyz).

## Contents

There is one notebook in this repository:

- [notebook/reproduce_paper.ipynb](notebook/reproduce_paper.ipynb)

Part 1 builds the isolated-outcrop forward model with SimPEG and inverts it. Part 2 redraws Figures 1-8 and
recomputes every number quoted in the paper (374 of them) from the stored results in `results/`, each printed next
to the sentence it appears in.

## Usage

Install Python through [Miniforge](https://github.com/conda-forge/miniforge) or Anaconda, then

```
git clone https://github.com/christopherzareck2/lunar-bedrock-magnetometry.git
cd lunar-bedrock-magnetometry
conda env create -f environment.yml
conda activate lunar-mag
jupyter lab notebook/reproduce_paper.ipynb
```

## Running the notebook

Part 1 runs a forward model, generates synthetic noise draws, and performs inversions, even when the final
`RERUN` flag is `False`. Part 2 redraws figures and computes paper numbers from the stored results.
The final `RERUN` flag controls only the scripts listed in the last section; setting it to `True` regenerates
the stored results (about 3.5 hours on a laptop). All random draws are seeded.

To check the quoted numbers using only existing results, run from the repository root:

```
python -c "import sys; sys.path.insert(0, 'notebook'); from paper_numbers import count; matched, total = count(); print(f'{matched} of {total} numbers match'); assert matched == total"
```

This checks consistency with the quoted values embedded in `paper_numbers.py`; it does not independently
validate the manuscript or rerun the simulations.

`python tests/reproduce_check.py` reruns two parts of the analysis from scratch in about a minute and checks them
against the stored results.

| Folder | Contents |
|---|---|
| `notebook/` | the notebook and `paper_numbers.py` (the paper's numbers, computed from `results/`) |
| `engine/` | forward model and inversion (`core.py`, `inversion_helpers.py`) and the base simulation scripts (the two in `engine/scripts/` run through `analysis/a05`) |
| `analysis/` | the scripts behind the paper's numbers (`a02`-`a09`) |
| `figures/` | one plotting script per figure (`scripts/`) and the figures as PDF and PNG (`output/`) |
| `results/` | stored results (`engine/` from the base scripts, `analysis/` from `analysis/`) |
| `tests/` | `reproduce_check.py`, a one-minute rerun of two parts of the analysis |

| Paper | Script |
|---|---|
| Table 2 | `engine/inversion_approach_comparison.py` |
| Section 4.2, Figures 4 and 5 | `analysis/a02_isolated_outcrop_paper_geometry.py`, `analysis/a02b_detection_sweep_N1000.py`, `analysis/a04_text_number_checks.py` |
| Table 3, Figure 8 | `analysis/a06_isolated_robustness.py` |
| Section 4.3, Figure 6 | `engine/skatepark_1cm_dr_sweep.py`, `analysis/a08_skatepark_neighbor_directions.py` |
| Section 4.4, Table 4, Figure 7 | `engine/three_layer_wall.py`, `engine/scripts/elevated_traverse_height_sweep.py`, `analysis/a03_three_layer_direction_vs_distance.py`, `analysis/a07_three_layer_wall_extent.py` |
| Sections 3.3, 5.1 and 5.5 (background, correlated noise) | `analysis/a05_background_and_correlated_noise.py` |
| Sections 4.2 and 5.1-5.5, Table 5 | `analysis/a09_systematics.py` |

The stored results were produced with Python 3.10 and the package versions in `environment.yml`. Other versions
may change the last digit of some numbers.

## Issues

Please open an issue on this repository if something does not run or a
number does not match the paper.

## Citation

If you use this code, please cite the paper:

> Narrett, I. S., Zareck, C. A., Weiss, B. P., & Hoffman, J. A. (in preparation). Lunar Rover Magnetometry: Recovering
> Paleofield Intensity and Direction From Magnetized Bedrock. *The Planetary Science Journal*.

`CITATION.cff` has the software citation. The code is released under the MIT license (`LICENSE`).

## Source provenance

This release snapshot derives from `czareck2/lunar-bedrock-magnetometry` on MIT GitHub,
commit `ee48005e0459ab5533b1d6ddd812ea3365a69587`. The original repository and its Git history are
preserved separately. The public snapshot starts a separate release history to avoid publishing optional
contact information in historical metadata. Scientific authorship and the original copyright and license
are retained in `CITATION.cff` and `LICENSE`.
