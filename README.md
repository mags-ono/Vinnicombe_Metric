# Code for paper "Data-Driven Estimation of the Vinnicombe Metric"

This repository contains MATLAB code associated with the paper:

  **"Data-Driven Estimation of the Vinnicombe Metric"**  
  Margarita A. Guerrero, Henrik Sandberg, and Cristian R. Rojas  
  📄  [arXiv link](https://arxiv.org/pdf/2603.17545)  
  ✅ [Submitted to the 65th IEEE Conference on Decision and Control]

---

## Overview

This code implements a **data-driven method** to estimate the **Vinnicombe metric** (also known as the **ν-gap metric**) between two discrete-time linear time-invariant systems.
The method relies on **time-domain input-output experiments** and does **not require an explicit parametric model of the true plant**. Instead, it assumes that a **nominal model** is available and that experiments can be performed on the unknown system.

The approach combines two computations:

- a **power-iteration-inspired input redesign** to estimate the ν-gap from experimental data, and
- a **Welch-inspired frequency-domain check** of the admissibility condition, implemented through the winding numbers of the auxiliary functions $f_1$ and $f_2$.

As a baseline, the script also compares the data-driven estimate with MATLAB's `gapmetric` function from the Robust Control Toolbox.

---

## Contents

The repository provides a single, **self-contained MATLAB script**:

```text
nu_gap_time_with_C_welch.m
```

This script includes:

- utility functions to convert transfer functions into the internal simulation format,
- a discrete-time plant simulator with additive output noise,
- Welch-style spectral accumulation across iterations,
- a data-driven check of the admissibility condition \((G_0,G)\in\mathcal{C}\),
- the main ν-gap estimation loop,
- plotting routines for convergence and input-distribution evolution,
- and a ground-truth winding-number check from exact transfer functions.

---

## Main helper functions

Below is a brief description of the main routines defined in the script.

### `tf_to_cell_struct(P)`
Converts a MATLAB transfer function object into the cell-based structure used by the simulator.
Each entry stores the numerator and denominator coefficients as fields `b` and `a`.

### `apply_P(P, u, sigma)`
Wrapper that simulates the response of a plant to an input signal `u` with output-noise level `sigma`.

### `simulate_plant_general_discrete(G0, u_time, sigma)`
Simulates the time-domain response of a discrete LTI plant represented in cell format.  
For each output channel, the function filters the corresponding input channels and adds Gaussian output noise.

### `specacc_init(N)`
Initializes the spectral accumulator used for the Welch-inspired admissibility check.

### `specacc_update(acc, U, Y1, Y2)`
Updates the accumulated spectra using the current FFTs of the input and of the two output signals.

### `specacc_get_frf(acc, eps0)`
Builds frequency-response estimates $\hat{G_1}$ and $\hat{G_2}$ from the accumulated spectra.

### `check_C_condition_from_acc(acc, tol_f, eps0, do_smooth, smooth_bins)`
Checks the admissibility condition $(G_0,G)\in\mathcal{C}$ in a data-driven way.
It computes the functions

- $f_1=1+|\hat G_1|^2$,
- $f_2=1+\hat G_2\overline{\hat G_1}\$,

and verifies:

- that their minimum magnitudes remain above a threshold, and
- that their winding numbers coincide.

### `winding_number(z)`
Computes the winding number of a sampled complex curve by summing unwrapped phase increments.

### `hann_kernel(L)`
Returns a normalized Hann kernel used for frequency smoothing.

### `frf_dt_tfdata_z(P, w)`
Evaluates the exact discrete-time frequency response of a transfer function `P` over a frequency grid `w`.
This is used at the end of the script for the ground-truth winding-number check.

---

## Example included: Nominal Model (Fig. 2)

The script includes the example corresponding to the **Nominal Model case** in **Fig. 2**.
In this example, the unknown plant and the nominal model are both derived from simplified heavy-duty gas-turbine fuel-system dynamics based on **Rowen-type modeling**, using two parameter sets taken from the literature (here referred to as the **Tavakoli** and **Khormali** parameterizations).

The continuous-time transfer functions are defined as

```matlab
Ts = 0.05;
s = tf('s');

% Parameters
Tf  = 0.26;     Kf  = 0;   Tcd = 0.16;    Kt = 1.158;   ecr = 0.005;   % Real Plant
Tf0 = 0.41;     Kf0 = 0;   Tcd0 = 0.0823; Kt0 = 0.8302; ecr0 = 0.005;  % Nominal Plant

% Real Model
Gf  = 1/(Tf*s + 1);
Gcd = 1/(Tcd*s + 1);
Pc  = Kt * Gf * Gcd;

% Nominal Model
Gf0  = 1/(Tf0*s + 1);
Gcd0 = 1/(Tcd0*s + 1);
Pc0  = Kt0 * Gf0 * Gcd0;

% Discretization
P1 = c2d(Pc,  Ts, 'zoh');
P2 = c2d(Pc0, Ts, 'zoh');
```

That is, the two plants are modeled as cascades of first-order fuel-system and compressor-discharge dynamics, scaled by an overall turbine gain. After zero-order-hold discretization, the script treats:

- `P1` as the **true plant**, and
- `P2` as the **nominal model**.

The algorithm then:

1. computes the MATLAB baseline using `gapmetric`,
2. runs the iterative ν-gap estimator over multiple Monte Carlo realizations,
3. checks the admissibility condition using the Welch-inspired estimator of $f_1$ and $f_2$,
4. plots the convergence of the estimate, and
5. reports the dominant excited frequency together with the ground-truth winding numbers.

---

## Main script settings used in the example

The provided Nominal Model example uses the following main parameters:

```matlab
M        = 150;    % number of iterations / experiments
N        = 9000;   % signal length
Nacc     = 10;     % number of experiments used for the admissibility check
Nmc      = 100;    % number of Monte Carlo runs
gamma    = 1;      % input-energy scaling
sigma_y1 = 0.01;   % output-noise variance on plant 1
sigma_y2 = 0.01;   % output-noise variance on plant 2
tol_f    = 1e-2;   % threshold for min|f_i|
smooth_bins = 31;  % Hann smoothing width
```

---

## Example output

A typical console output for the Nominal Model case is:

```matlab
MC=100 | C-check n=10: min|f1|=1.00e+00 min|f2|=1.00e+00 | wind1=0 wind2=0 | inC=1

===== PEAK ANALYSIS =====
Peak input-distribution magnitude = 1.992817e+06
Positive peak freq  = 1.361357e-01 rad/sample
Mirror peak freq    = 6.147050e+00 rad/sample
|p0(w_peak)|        = 0.541266
angle p0(w_peak)    = -1.102954 rad
=====GROUND TRUTH Winding Number===== 
tfdata f1: wind=0 
tfdata f2: wind=0
```

This output indicates that:

- the admissibility condition is satisfied (`inC=1`),
- both estimated winding numbers are zero,
- and the input redesign concentrates energy around the dominant frequency reported in the peak analysis.

---

## Notes

- The script is written for the **SISO** case.
- The admissibility condition is checked numerically from data using the auxiliary functions $f_1$ and $f_2$.
- Noise is included in the simulations, while the exact transfer functions are only used for baseline comparison and for the final ground-truth winding-number verification.
- Additional case studies can be activated directly in the script by uncommenting the corresponding blocks.

---

## Citation

If you use this code in your work, please cite:

```bibtex
@inproceedings{Guerrero_26,
  title     = {Data-Driven Estimation of the Vinnicombe Metric},
  author    = {Guerrero, Margarita A., Sandberg, Henrik, and Rojas, Cristian R.},
  year      = {2026},
  note      = {Submitted to the 65th IEEE Conference on Decision and Control (CDC). Also available at arXiv:2603.17545}
}
```
