# PhysicsProjects

## Overview

This repository contains several physics-focused MATLAB projects and assignment solutions. The work includes frequency analysis of audio signals, normal mode analysis of a DNA model, a supercapacitor energy storage model, and a time-dependent Schrödinger equation solver (Crank–Nicolson).

## Requirements

- MATLAB (Live Editor recommended for `.mlx` files)
- Basic toolchain: no external libraries required (pure MATLAB code)

## Repository structure and brief descriptions

- `FrequencyAnalysis/`
	- Live scripts: `Frequency_analysis_a_b.mlx`, `Frequency_analysis_c_e.mlx` — contain the main analysis, visual outputs, and commentary.
	- `Frequency_analysis_*.html` — exported HTML versions for quick viewing.
	- `PlainCodes/` — plain `.m` implementations (e.g., `Frequency_analysis_a_b.m`, `Frequency_analysis_c_e.m`) including a manual DFT implementation and FFT usage.
		- These scripts read an audio file (e.g., `single_tone.ogg`) and compute spectra, PSD, and approximations. The DFT implementation is intentionally slow (educational) and includes timing comparisons with `fft`.
	- `ApproximatedAudioFiles/` and `Plots/` — exported results (audio and figures) from the analyses.

- `NormalModeAnalysisDNA/`
	- `DNA.m` — builds a model (from `xyzm_dna.txt`) and computes eigenmodes using Hessian/Mass matrices and the power method; produces eigenvalue/eigenvector plots and animation code (some animations exported as GIFs). Some variants and helper scripts are in `ANM/`.
	- `ANM/ANM.m` — another version of normal mode analysis / visualization code.
	- `xyzm_dna.txt`, `time_evolution.txt` — data files used by the scripts.
	- Assignment PDF(s) providing problem statements are included.

- `SupercapacitorEnergyStorage/`
	- `SupercapacitorEnergyStorage.m` — numerical solution of a boundary-value problem for ionic charge distribution and capacitance calculations. Includes LU decomposition helper, comparisons with analytical solutions, and plots of capacitance versus width, Debye length and temperature.
	- Several exported plot images and `assignment_1.pdf` are included.

- `TheTimeDependentSchrodingerEquation/`
	- `Crank_Nicolson.m` — implements a Crank–Nicolson time integrator for the 1D time-dependent Schrödinger equation, with handlers for free propagation and propagation through potential barriers. Several exported GIFs illustrate dynamics.