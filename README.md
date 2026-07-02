# Synthetic Ultrasound RF Channel Data for Time-Harmonic Elastography

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21142261.svg)](https://doi.org/10.5281/zenodo.21142261)

This repository contains MATLAB code for generating synthetic ultrasound RF
channel data for ultrasound time-harmonic elastography (THE). The framework
links numerically simulated shear-wave motion to pulse-echo ultrasound
simulation by encoding the motion into a time-varying acoustic scatterer
distribution.

The code accompanies the manuscript:

> Synthetic ultrasound channel data generation for time-harmonic shear-wave
> elastography

The examples reproduce the two case studies described in the manuscript:

- a phantom-mimicking model simulated with the k-Wave elastic solver;
- a realistic three-dimensional liver model represented by a tetrahedral FEM
  mesh and simulated motion stored in an HDF5 file.

## Overview

The simulation pipeline follows three main steps:

1. Compute or load a time-dependent shear-wave motion field.
2. Convert the motion field into time-varying scatterer coordinates.
3. Use Field II to generate synthetic pulse-echo RF channel data.

Optional sections in the main scripts compute ground-truth motion fields on the
ultrasound image grid and beamform the RF channel data with USTB.

## Requirements

The code is written for MATLAB. It has been checked in this workspace with
MATLAB R2024a.

Required external toolboxes/software:

- [Field II](https://field-ii.dk/) for RF channel-data simulation.
- [k-Wave](http://www.k-wave.org/) for the phantom shear-wave simulation.
- [USTB](https://ultrasoundtoolbox.com/) for the optional beamforming sections.

If USTB is not installed, comment out the optional beamforming section in the
main script after the RF data have been generated.

## Setup

1. Clone or download this repository.

2. Start MATLAB in the repository root.

3. Edit the dependency paths at the top of the main scripts if needed:

```matlab
addpath('/path/to/field_II/')
addpath('/path/to/ustb')
```

4. Edit the output path `SAVE_DIR` near the end of each main script so that it
points to a writable location on your machine.

The liver case expects `Liver-sim-data.h5` to be available in the repository
root. The included file stores the tetrahedral mesh, cell labels and simulated
FEM displacement data used by `sim_liver_main.m`.

## Running the Examples

### Phantom Model

Run:

```matlab
sim_phantom_main
```

### Liver Model

Run:

```matlab
sim_liver_main
```

## Main Outputs

The main scripts save variables such as:

- `RF`: simulated RF channel data from Field II;
- `disp_gt`: ground-truth motion projected onto the ultrasound image direction;
- `channel_data`: USTB channel-data object that contain all the required parameters for beamforming;


The RF array is arranged as:

```text
[samples, receive elements, transmits, frames]
```

## Notes on Reproducibility

- The phantom example performs a full k-Wave shear-wave simulation and can be
  computationally expensive.
- A large number of scatterers is required to establish fully developed
  ultrasound speckle statistics. Consequently, the Field II simulation is also
  computationally expensive; in the current examples, one RF frame can take
  roughly half an hour to compute, depending on hardware and scatterer count.
- Field II simulation is run inside `parfor` in `FIELD_calc_RF.m`.

## Citation

If you use this code, please cite the accompanying manuscript. A final citation
entry will be added when the paper details are available.

```bibtex
@article{han2026syntheticTHE,
  title   = {Synthetic ultrasound channel data generation for time-harmonic shear-wave elastography},
  author  = {Han, Chaoran and N{\ae}sholm, Sven Peter and Austeng, Andreas and Ziksari, Mahsa Sotoodeh and Karabiyik, Y{\"u}cel},
  journal = {IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control},
  year    = {2026},
  note    = {Manuscript under review}
}
```

## License

This project is released under the BSD 3-Clause License. See
[`LICENSE`](LICENSE) for details.
