# MCCCII Current-Mode Biquad Filter

**主動式電流模式二階濾波器設計與實現**  
*Realization of Current-Mode Biquad using a Single Multi-output Current Controlled Conveyor (MCCCII)*

> Undergraduate electrical engineering project portfolio focused on analog current-mode filter design, transfer-function analysis, sensitivity analysis, and experimental verification.

---

## Project Overview

This project studies a current-mode biquadratic filter architecture implemented with a **single Multi-output Current Controlled Conveyor (MCCCII)** and passive elements. The design supports three second-order filter functions from the same active-device framework:

- Low-Pass Filter (LPF)
- High-Pass Filter (HPF)
- Band-Pass Filter (BPF)

The work covers circuit architecture, transfer-function derivation, center-frequency and quality-factor analysis, non-ideal behavior, sensitivity analysis, and comparison between theoretical and experimental/simulation results.

## Key Topics

`Analog Circuit Design` · `MCCCII` · `Current-Mode Circuit` · `LPF` · `HPF` · `BPF` · `Transfer Function` · `Sensitivity Analysis` · `MATLAB`

---

## 1. MCCCII Architecture

The MCCCII is described by the ideal relations:

```text
Vx = Vy + IxRx
Iy = 0
Iz+ = Ix
Iz- = -Ix
```

A single MCCCII combined with passive admittances can be configured to realize multiple second-order current-mode filtering functions.

---

## 2. Filter Realizations

### Low-Pass Filter (LPF)

For one realization, the project uses:

```text
Y1 = sC1
Y2 = G2
Y3 = sC3
Y4 = G4
```

The transfer-function form is:

```text
Io/Iin = G2G4 /
         [s²C1C3 + s(C1G2 + C3G2) + G2G4]
```

with center frequency:

```text
ω0 = sqrt(G2G4 / (C1C3))
```

### High-Pass Filter (HPF)

The high-pass configuration is implemented by selecting a different set of passive admittances around the same MCCCII structure. Its gain and phase responses are used to verify the expected second-order high-pass behavior.

### Band-Pass Filter (BPF)

The band-pass realization uses the same current-mode design philosophy and is evaluated through center frequency, quality factor, gain response, and phase response.

---

## 3. Mathematical Analysis

The project derives and evaluates:

- Transfer function H(s)
- Center frequency ω0
- Quality factor Q
- Frequency response
- Gain response
- Phase response

The main engineering objective is to verify that the proposed filter configuration can realize useful second-order functions while keeping the circuit compact.

---

## 4. Sensitivity Analysis

To study non-ideal behavior, the MCCCII equations are extended with tracking factors:

```text
Vx = αVy + IxRx
Iz+ = βIx
```

The project then evaluates how circuit parameters influence:

- Center frequency ω0
- Quality factor Q
- Passive sensitivity
- Active-device tracking error

A major conclusion is that the proposed current-mode filters exhibit low passive sensitivity, while key filter parameters remain comparatively robust against MCCCII current-tracking error.

---

## 5. Experimental / Simulation Verification

The project compares theoretical predictions with measured or simulated frequency-response results for:

- LPF gain and phase
- HPF gain and phase
- BPF gain and phase

The reported results show only small deviations between theory and experiment/simulation, supporting the feasibility of the proposed topology.

---

## 6. Engineering Skills Demonstrated

- Analog circuit analysis
- Active current-mode filter design
- Second-order system modeling
- Transfer-function derivation
- Frequency-domain analysis
- Sensitivity analysis
- Parameter interpretation
- MATLAB-assisted signal / response analysis
- Technical literature review
- Engineering report writing

---

## 7. Project Structure

```text
MCCCII-Current-Mode-Biquad-Filter/
├── README.md
├── docs/
│   └── project-summary.md
└── references/
    └── references.md
```

---

## 8. Project Summary

The project demonstrates that a **single MCCCII** can be used to realize compact current-mode second-order low-pass, high-pass, and band-pass filters. Through transfer-function analysis, sensitivity analysis, and response verification, the design shows the practical value of current-mode active circuits for flexible analog signal-processing applications.

---

## Author

Electrical Engineering undergraduate project portfolio  
Fu Jen Catholic University