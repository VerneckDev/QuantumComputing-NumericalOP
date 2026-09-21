# Numerical Pulse Optimization in NMR (SPINQ Gemini Lab)

MATLAB codes for the work **"Integrated Quantum Technology: Numerical Pulse Optimization"**, from the *Quantum Technologies 2025/2026* course unit (FCTUC – Department of Physics, University of Coimbra).

The work uses a 2-qubit NMR quantum computer (**¹H** and **³¹P**, dimethyl methylphosphonate sample) to compare three ways of implementing logical operations:

| Method                                | Target                                   | Files                          | Fidelity reported in the report                |
| ------------------------------------- | ---------------------------------------- | ------------------------------ | ---------------------------------------------- |
| Hard pulses                           | Unitary gate `U = U¹ ⊗ U²`               | `FidCal.m`                     | 90.04 % (see [Note 1](#1-hard-pulse-fidelity)) |
| **SMP** (*Strongly Modulated Pulses*) | Same unitary gate                        | `SMP_OP.m`, `SMP_FidCal.m`     | 96.82 %                                        |
| **GRAPE**                             | Bell state `\|β₀₀⟩ = (\|00⟩ + \|11⟩)/√2` | `GRAPE_OP.m`, `GRAPE_FidCal.m` | 92.86 %                                        |

> The SMP and GRAPE targets are different, so the fidelities **do not** form a direct ranking of the methods (as stated in the report's conclusion).

---

## Repository Structure

```text

├── Quantum_Technology.pdf   # Complete report
├── FidCal.m                 # Fidelity and plots: hard pulses (theoretical vs. experimental)
├── SMP_OP.m                 # SMP optimization -> generates JSON for the equipment
├── SMP_FidCal.m             # Fidelity and plot: experimental SMP result
├── GRAPE_OP.m               # GRAPE optimization (Bell state) -> generates JSON
├── GRAPE_FidCal.m           # Fidelity and plot: experimental GRAPE result
└── README.md
```

There are two types of scripts:

* **`*_OP.m`** – *optimization*: calculate the pulses and export a `.json` file to be loaded into the SPINQ Gemini Lab.
* **`*_FidCal.m`** – *post-processing*: compare the experimental density matrix exported from the equipment with the theoretical one, calculate the fidelity, and plot the 3D histograms.

## Workflow

```text
   *_OP.m  ──►  .json file  ──►  SPINQ Gemini Lab  ──►  experimental density matrix

  (optimization)    (pulses)          (execution, 5 repetitions)             │
                                                                             ▼
                                                     copy to *_FidCal.m  ──►  fidelity + plots
```

## Requirements

* MATLAB **R2021a or later** (`jsonencode` with `'PrettyPrint'`).
* **Optimization Toolbox** – `fminunc` (SMP) and `fmincon` (hybrid function of GRAPE).
* **Global Optimization Toolbox** – `particleswarm` (GRAPE).

The `*_FidCal.m` scripts only require base MATLAB (`sqrtm`, `bar3`).

---

## Physical Parameters Used

Defined at the top of the `*_OP.m` scripts, using the values measured in the report:

| Symbol   | Value    | Meaning                                                                     |
| -------- | -------- | --------------------------------------------------------------------------- |
| `J`      | 695.8 Hz | ¹H – ³¹P scalar coupling constant (Eq. 15)                                  |
| `w1_max` | 6250 Hz  | Maximum Rabi frequency (measured: 6257.66 Hz for ¹H and 6187.98 Hz for ³¹P) |

Hamiltonian in the rotating frame (zero offsets, exact resonance):

```text
H(t) = H₀ + H₁(t)

H₀   = 2π·J · (σz ⊗ σz) / 4 (coupling, = 2πJ·Iz·Iz)

H₁   = 2π·A(t)·w1_max · [cos φ · σx/2 + sin φ · σy/2] (on each channel)
```

Each segment `k` has a constant Hamiltonian, and the total propagator is the ordered product

`U = ∏ exp(-i·Hₖ·tₖ)`.

Control parameters per segment (two channels: **H** = channel 1, **P** = channel 2): amplitude `A ∈ [0, 1]` (fraction of `w1_max`) and phase `φ`.

---

## Script Description

### `FidCal.m` — Hard Pulses

Defines the target gate from the report (Eq. 39):

```text
U = R²₋z(π/4) · R²₋y(π/4) · R¹z(π/4) · R¹y(π/4)   →   U = U¹ ⊗ U²
```

In the code, `R = cos(π/8)·I ∓ i·sin(π/8)·σ`, corresponding to rotations of `θ = π/4`. Then:

1. Checks unitarity (`U†U = I`).
2. Calculates `|ψ_teo⟩ = U|00⟩` and `ρ_teo = |ψ⟩⟨ψ|`.
3. Reads the experimental density matrix (`rho_exp_r + rho_exp_img`, manually pasted).
4. Calculates the fidelity `F = Tr √(√ρ_exp · ρ_teo · √ρ_exp)` (Eq. 48).
5. Generates **two figures**: theoretical and experimental matrices (real part in red, imaginary part in blue, using the magnitude).

### `SMP_OP.m` — SMP Optimization

* **Method:** `fminunc` (unconstrained), with `N = 6` segments.
* **Variables (5·N = 30):** `[A_H, φ_H, A_P, φ_P, t]`, with durations `t` in µs (initialized between 20 and 30 µs). The amplitudes are limited to `[0, 1]` inside the cost function; the duration uses `|t|`.
* **Cost:** `Q = 1 − |Tr(U_alvo† · U_SMP)| / 4` (Eqs. 36–37, gate fidelity).
* **Initialization:** random. The number of trials is `num_trials` (set to `1` in the file; increase it to explore more starting points, as described in Section 3.4 of the report).
* **Output:** `SMP_pulse_no_detuning.json`.

### `SMP_FidCal.m`

Same as `FidCal.m`, but using the experimental matrix obtained with the SMP pulses and only the experimental plot (title *"Experimental Density Matrix - SMP"*).

### `GRAPE_OP.m` — GRAPE Optimization (Bell State)

* **Discretization:** `N = 300` segments with fixed duration, `T_total = 1000 µs` (`dt ≈ 3.33 µs`).
* **Variables (4·N = 1200):** `[A_H, φ_H, A_P, φ_P]` per segment.
* **Optimizer:** `particleswarm` with a hybrid `fmincon` function, 5 trials (`num_trials`), keeping the best result.
* **Cost:** `Q = 1 − |⟨β₀₀|U_total|00⟩|²` (**state-to-state** fidelity).
* **Output:** `GRAPE_Bell_pulse.json`.

### `GRAPE_FidCal.m`

Uses the unitary target matrix from Eq. 50 (first column = Bell state), calculates `ρ_teo` for `|00⟩`, compares it with the GRAPE `ρ_exp`, and plots the experimental matrix.

---

## JSON File Format

The `generate_json_file` (SMP) and `generate_json_grape` (GRAPE) functions write the format read by the SPINQ Gemini Lab:

```text
description   → TITLE, OWNER, DATE, FIDELITY, TOTALPULSEWIDTH, SLICES

parameters    → channel offsets (all 0)

pulse
 ├─ channel1_pulse{k} → detuning (0), phase [degrees, 0–360), amplitude [%], width [µs]
 └─ channel2_pulse{k} → idem
```

Conversions performed during export: amplitude `A × 100` (percentage), phase `mod(φ·180/π, 360)` (degrees), width in µs.

## How to Run

1. **Optimize:** open `SMP_OP.m` or `GRAPE_OP.m` in MATLAB and run it. The `.json` file is created in the current folder.
2. **Experiment:** load the `.json` into the SPINQ Gemini Lab software, run it (the report repeats the experiment 5 times), and export the density matrix.
3. **Analyze:** copy the real and imaginary parts of the matrix into `rho_exp_r` and `rho_exp_img` in the corresponding `*_FidCal.m` script and run it to obtain the fidelity and plots.

> `GRAPE_OP.m` is computationally intensive (swarm of 100 particles, 1200 variables, 300 matrix exponentials per evaluation, 5 trials). Reduce `num_trials`, `MaxIterations`, or `N` for quick tests.

---

## Notes and Points to Consider

When reviewing the codes against the report, I found the following points worth taking into account:

### 1. Hard-Pulse Fidelity

Running `FidCal.m` with the experimental matrix included in the file gives **`F ≈ 0.9504 (95.04 %)`**, rather than 90.04 % as stated in the report (Section 3.4.1 and Conclusion). The SMP (**0.9682**) and GRAPE (**0.9286**) values match the report. If 95.04 % is the correct value, the improvement of SMP over hard pulses is approximately 1.8 percentage points (rather than approximately 6.8). It is worth confirming whether the matrix pasted into `FidCal.m` is the same one used in the report.

### 2. Fidelity Definition

The `*_FidCal.m` scripts use `F = Tr √(√ρ_exp ρ_teo √ρ_exp)` (the "root" fidelity, without squaring). For a pure target state, this is equivalent to `√⟨ψ|ρ_exp|ψ⟩`. The GRAPE cost uses `|⟨ψ_target|ψ⟩|²` (with the square), while SMP uses the gate fidelity `|Tr(U†V)|/4`. Therefore, the optimization (simulation) values and experimental values are not directly comparable without taking this into account.

### 3. "GRAPE" vs. Implementation

The report describes GRAPE mainly through its advantage of efficient gradient calculation. In the code, however, the optimization uses `particleswarm` (gradient-free) with refinement using `fmincon`, without an analytical gradient being provided. In practice, it is a fixed-segment pulse optimization, but not the classical GRAPE algorithm. It is worth clarifying this in the report or implementing the analytical gradient.

### 4. SMP Parameters

* The report (Section 3.2) lists the **RF frequency** `f_rf` as an SMP parameter, but the code only optimizes amplitude, phase, and duration (`detuning = 0`, hence the name `SMP_pulse_no_detuning.json`).
* The duration `t(k)` is **shared** by both channels in each segment.
* `SMP_OP.m` is set to `num_trials = 1`, while the report mentions several starting points.
* The durations have no upper bound, so the total pulse duration is not controlled.

### 5. Simulated vs. Experimental Fidelities

The `*_OP.m` scripts store the **simulated** fidelity in the `FIDELITY` field of the JSON. This value does not include relaxation (T₁, T₂), magnetic-field inhomogeneities, or calibration errors, which explain the difference from the experimental values.

### 6. Experimental Data

The experimental density matrices are manually inserted into the `*_FidCal.m` scripts (values exported from the SPINQ Gemini Lab). If the experiment is repeated, replace these values.

---

## References

1. G. Feng *et al.*, "A commercial three-qubit desktop quantum computer", *IEEE Nanotechnology Magazine*, 16, 2022. [doi:10.1109/MNANO.2022.3175392](https://doi.org/10.1109/MNANO.2022.3175392)

2. T. S. Mahesh and D. Suter, "Quantum-information processing using strongly dipolar coupled nuclear spins", *Phys. Rev. A*, 74, 2006. [doi:10.1103/PhysRevA.74.062312](https://doi.org/10.1103/PhysRevA.74.062312)

3. SpinQ Technology, *SpinQ Gemini Lab Experiment Manual*, 2024.

4. X. Yang *et al.*, "Assessing three closed-loop learning algorithms by searching for high-quality quantum control pulses", *Phys. Rev. A*, 102, 2020. [doi:10.1103/PhysRevA.102.062605](https://doi.org/10.1103/PhysRevA.102.062605)

## Authorship

**João Pedro Verneck** — Supervision: Prof.ª Maria Helena Almeida Vieira Alberto

FCTUC – Department of Physics, University of Coimbra · 2025/2026
