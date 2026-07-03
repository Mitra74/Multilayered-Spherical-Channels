# Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels

Code accompanying **Chapter 4** of the PhD thesis *"Diffusion-Based Molecular Communication in Discrete Heterogeneous Environments"* (Mitra Rezaei, University of Warwick).

## Paper
Rezaei, M., Chappell, M., and Noel, A. "Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels." *IEEE Transactions on Molecular, Biological, and Multi-Scale Communications*, 2026.

*Preliminary version: Rezaei, M., Chappell, M., and Noel, A. "General Molecular Communication Model in Multi-Layered Spherical Channels," IEEE International Conference on Communications (ICC), 2025.*

## Overview
The paper derives a generalized analytical framework for diffusion-based molecular communication in multi-layered spherical environments, supporting an arbitrary number of layers and flexible transmitter–receiver positioning (transmitters cannot be located at layer boundaries). The Green's function for the boundary value diffusion problem is derived for the layer containing a point source, with the channel impulse response derived for all other layers, coupled through interface continuity and flux boundary conditions.

This repository implements the case study from the paper: a three-layer tumour spheroid (necrotic core, hypoxic layer, proliferating outer layer) embedded in an infinite medium, modelling direct delivery from a point-source transmitter positioned either externally or internally to the spheroid. The general N-layer framework itself is not implemented here in its arbitrary-layer form — the two scripts below solve the specific, closed-form three-layer boundary-value system.

All analytical results are validated against a custom particle-based simulation (PBS), which accounts for short inter-layer distances by applying multiple diffusion-coefficient updates within a single molecule's time-step trajectory.

## Repository structure
```
src/
├── external_transmitter_3Layered_spheroid.m   — direct delivery: external point source 
└── internal_transmitter_3Layered_spheroid.m   — direct delivery: internal point source 
```

Each script is self-contained: it builds and solves the boundary-value diffusion problem for a three-layer spheroid, reconstructs the frequency-domain concentration at a chosen observation point via a spherical-harmonic expansion, applies an inverse FFT to obtain the time-domain concentration profile, and plots the result. Comments throughout reference the corresponding equation numbers in the paper (Eqs. 6–45) so the code can be read alongside the derivation.

## Requirements
- MATLAB R2021b or later
- No additional toolboxes required (uses base MATLAB `besselj`, `bessely`, `besselh`, `legendre`)

## Usage
1. Open MATLAB and navigate to `src/`
2. Run either script directly, e.g.:
```matlab
   external_transmitter_3Layered_spheroid
```
   or
```matlab
   internal_transmitter_3Layered_spheroid
```
3. Output: each script produces a concentration-vs-time plot at the selected observation point.

Key parameters that can be adjusted at the top of each script:
- `N` — number of spherical-harmonic degrees kept in the truncated series expansion
- Layer radii (`R1`, `R2`, `R3`) and porosities (`eps_L1`/`eps_R1`, `eps_L2`/`eps_R2`, `eps_L3`/`eps_R3`), which set the effective diffusion coefficients `D_i = eps_i^1.5 · D`
- Transmitter position `r_tx` (and `phi_tx`, `theta_tx`)
- Observation point, via `obs_idx` (selects among the pre-defined layer/outside/center radii)
- Bulk diffusion coefficient `D_drugs`
- Layer-specific degradation rates (`kd_*`), left at zero here for the direct-delivery scenario

## Model summary
- **Geometry**: spherical structure with `N_L` finite concentric layers plus an unbounded outer layer, each layer `i` characterized by an effective diffusion coefficient `D_i = (ε_i/τ_i)·D`, where `ε_i` is porosity and `τ_i` is tortuosity.
- **Boundary conditions**: flux continuity and a permeability-weighted concentration jump condition at each layer interface.
- **Solution method**: Green's function (for the source layer) and channel impulse response (source-free layers), expanded in spherical harmonics and solved via spherical Bessel/Hankel functions, with unknown coefficients determined from interface and source conditions, then converted to the time domain via inverse FFT over a frequency sweep.
- **Validation**: particle-based simulation in MATLAB with Δt = 0.5 s, Gaussian-distributed displacements per layer, and displacement-vector rescaling at every layer-boundary crossing to account for changing diffusion coefficients — including multiple crossings within a single time step.

## pH-triggered delivery: approach (not included in this repository)

The paper also analyzes a second delivery scenario — pH-sensitive drug carriers that remain stable while transiting the spheroid, then degrade and release their drug payload specifically in the acidic necrotic core (Section IV-B). This scenario is **not implemented in this repository**; it builds on the two direct-delivery models above in the following way:

1. **Carrier transport (outside → necrotic core).** Drug *carriers* are released from an external point source, exactly as in `external_transmitter_spheroid_diffusion.m`, except the diffusion coefficient used throughout is the *carrier's* diffusion coefficient rather than the free drug's. This gives the carrier concentration `U_c^1(r,t;r0)` inside the necrotic core (region 1) as a function of time.

2. **Degradation / release at the necrotic core.** Carriers are assumed to degrade — releasing their drug payload — only within the necrotic core, via a first-order process with rate `k1` (all other layers use zero degradation, since they are assumed to be at a higher, non-triggering pH). The probability flux of carrier degradation at time `t` is obtained by integrating the carrier concentration over the necrotic-core volume and multiplying by `k1`:
   `Ψ(r0,t) = k1 ∫ u_c^1(r,t;r0) dV`.

3. **Release function.** Discretizing time into intervals gives the probability that a carrier releases its payload within each interval — this is the **release function**, `ρ(t)`. It is essentially a histogram of "how much drug becomes available, and when," driven entirely by the carrier-transport solution from step 1 and the degradation rate `k1`. Ψ(r0,t) itself can also be considered as a release function without considering the probability. 

4. **Drug propagation from the release point.** Separately, compute the impulse response for a drug (not carrier) released from the *center* of the necrotic core — this is the internal-transmitter model (`internal_transmitter_spheroid_diffusion.m`), but with `r_tx → 0` and using the free *drug's* diffusion coefficient.

5. **Convolution.** The final pH-triggered drug concentration at any observation point is the convolution of the release function from step 3 with the point-release impulse response from step 4:
   `T_pH-released(r,t) = ρ(t) * T(r,t | r0 = 0)`.

   Intuitively: the release function tells you *when and how much* drug gets released at the core, and the point-source impulse response tells you *how that drug then spreads* through the spheroid — convolving the two gives the full spatiotemporal drug concentration.

If you need the implementation of this pH-triggered scenario, please get in touch and I'm happy to share/discuss that part of the code.

## Citation
```bibtex
@article{rezaei2026multilayerspherical,
  author={Rezaei, Mitra and Chappell, Michael J. and Noel, Adam},
  journal={IEEE Transactions on Molecular, Biological, and Multi-Scale Communications}, 
  title={Molecular Communication Model for Drug Delivery in Multi-Layered Spherical Channels}, 
  year={2026},
  volume={12},
  number={},
  pages={309-322},
  keywords={Drugs;Tumors;Transmitters;Molecular communication;Extracellular;Toxicology;Biological system modeling;Receivers;Numerical models;Geometry;Molecular communication;diffusion;channel model;multi-layered spherical structures;tumor spheroids;drug delivery systems},
  doi={10.1109/TMBMC.2026.3657747}}

```

## Author
Mitra Rezaei — [LinkedIn](https://www.linkedin.com/in/mitra-rezaei-834784159/) · [Thesis hub](https://github.com/Mitra74/phd-thesis-warwick)
