# Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels

Code accompanying **Chapter 4** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Paper
Rezaei, M., Chappell, M., and Noel, A. "Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels." *IEEE Transactions on Molecular, Biological, and Multi-Scale Communications*, 2026.

*Preliminary version: Rezaei, M., Chappell, M., and Noel, A. "General Molecular Communication Model in Multi-Layered Spherical Channels," IEEE International Conference on Communications (ICC), 2025.*

## Overview
This repository implements a generalized analytical framework for diffusion-based molecular communication in multi-layered spherical environments, supporting an arbitrary number of layers and flexible transmitter–receiver positioning. The Green's function for the boundary value diffusion problem is derived for a layer containing a point source, with the channel impulse response derived for all other layers, coupled through interface continuity and flux boundary conditions.

The framework is applied to a three-layer tumour spheroid case study (necrotic core, hypoxic layer, proliferating outer layer) embedded in an infinite medium, modelling two drug delivery strategies:

1. **Direct delivery** — point source transmitter positioned externally or internally to the spheroid
2. **pH-triggered delivery** — drug carriers that remain stable during transport and release their payload upon reaching the acidic necrotic core, modelled as a first-order degradation process

All analytical results are validated against a custom particle-based simulation (PBS), which accounts for short inter-layer distances by applying multiple diffusion-coefficient updates within a single molecule's time-step trajectory.

## Repository structure
src/

├── main.m                          — entry point; sets layer/geometry parameters and runs analysis

├── greens_function.m               — Green's function for the layer containing the point source

├── impulse_response.m              — channel impulse response for source-free layers

├── boundary_conditions.m           — interface continuity & flux matching across layers

├── external_delivery.m             — direct delivery: external point source scenario

├── internal_delivery.m             — direct delivery: internal point source scenario

├── ph_triggered_delivery.m         — pH-triggered carrier release & drug concentration model

├── pbs_simulation.m                — particle-based simulation (validation)

└── plot_results.m                  — generates concentration profiles / colormap figures
figures/

├── concentration_profile_external.png   — Fig. 2 equivalent: external transmitter case

├── concentration_profile_internal.png   — Fig. 3 equivalent: internal transmitter case

├── porosity_effect.png                  — Figs. 4–7 equivalent: layer porosity impact

├── colormap_distribution.png            — Fig. 8 equivalent: cross-sectional molecular density

└── ph_triggered_comparison.png          — Figs. 9–13 equivalent: pH-triggered vs. direct release
data/

└── spheroid_parameters.mat         — layer radii, porosities (ε₁=0.2964, ε₂=0.1196, ε₃=0.1697), diffusion coefficients

## Requirements
- MATLAB R2021b or later
- Toolboxes: [list only if actually required, e.g. none / Symbolic Math Toolbox]

## Usage
1. Open MATLAB and navigate to `src/`
2. Run the main script:
```matlab
   main
```
3. Output: generates concentration-vs-time profiles at specified observation points, validates them against PBS, and produces cross-sectional colormap visualizations of molecular distribution, saved to `figures/`.

Key parameters that can be adjusted in `main.m`:
- Number of layers and layer widths (`L_i`)
- Layer porosities (`ε_i`) and resulting effective diffusion coefficients (`D_i = ε_i/τ_i · D`)
- Transmitter position (external or internal, at radial position `r₀`)
- Delivery scenario: direct release vs. pH-triggered carrier release
- Carrier degradation rate `k₁` (necrotic core sensitivity) for pH-triggered delivery

## Model summary
- **Geometry**: spherical structure with `N_L` finite concentric layers plus an unbounded outer layer, each layer `i` characterized by an effective diffusion coefficient `D_i = (ε_i/τ_i)·D`, where `ε_i` is porosity and `τ_i` is tortuosity.
- **Boundary conditions**: flux continuity and a permeability-weighted concentration jump condition at each layer interface.
- **Solution method**: Green's function (for the source layer) and channel impulse response (source-free layers), expanded in spherical harmonics and solved via spherical Bessel/Hankel functions, with unknown coefficients determined from interface and source conditions.
- **Validation**: particle-based simulation in MATLAB with Δt = 0.5 s, Gaussian-distributed displacements per layer, and displacement-vector rescaling at every layer-boundary crossing to account for changing diffusion coefficients — including multiple crossings within a single time step.

## Citation
```bibtex
@article{rezaei2026multilayerspherical,
  title={Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels},
  author={Rezaei, Mitra and Chappell, Michael and Noel, Adam},
  journal={IEEE Transactions on Molecular, Biological, and Multi-Scale Communications},
  year={2026}
}
```

## Author
Mitra Rezaei — [LinkedIn](https://www.linkedin.com/in/mitra-rezaei-834784159/) · [Thesis hub](https://github.com/Mitra74/phd-thesis-warwick)
