# Middle Permian Seawater ⁸⁷Sr/⁸⁶Sr Box Modeling: Comprehensive Implementation Plan

**Project**: Three Consecutive Declines Shaping the Middle Permian Seawater ⁸⁷Sr/⁸⁶Sr Trough
**Target Journal**: *Earth-Science Reviews* (following Wang et al., 2021 methodology)
**Study Section**: Dukou Section, Upper Yangtze Platform
**Date**: March 2026

---

## 1. Project Overview

### 1.1 Scientific Objective
Construct a stochastic oceanic Sr box model to quantitatively evaluate the coupled climatic-tectonic controls on the Middle Permian seawater ⁸⁷Sr/⁸⁶Sr evolution. Specifically, test the hypothesis that three distinct declining intervals (Roadian, Wordian, mid-Capitanian) were driven by different forcing mechanisms:

- **Decline 1 (Roadian)**: Late Paleozoic P3 glaciation → Reduced continental weathering
- **Decline 2 (Wordian)**: Neo-Tethys opening → Enhanced hydrothermal input
- **Decline 3 (mid-Capitanian)**: Emeishan LIP eruption → Basalt weathering pulse

### 1.2 Methodological Framework
Adopt the two-step Monte Carlo approach from **Wang et al. (2021, ESR)**:
1. **Step 1**: Uniform random sampling across broad parameter ranges (5,000 iterations)
2. **Step 2**: Gaussian refinement based on filtered results (2,000 iterations)
3. **Filtering**: Retain solutions matching observed ⁸⁷Sr/⁸⁶Sr within 95% confidence intervals

---

## 2. Mathematical Core

### 2.1 Sr Mass Balance Equation
Seawater ⁸⁷Sr/⁸⁶Sr at steady state:

$$R_{sw}(t) = \frac{\sum F_i(t) \cdot R_i(t)}{\sum F_i(t)}$$

Where $i$ = {riv, highT, lowT, dia}

### 2.2 Critical Innovation: Dual-Component Riverine Model (Decline 3)
For the Emeishan LIP interval, partition riverine flux into crustal and basaltic end-members:

$$F_{riv} \cdot R_{riv} = F_{crust} \cdot R_{crust} + F_{basalt}(t) \cdot R_{basalt}$$

$$F_{riv} = F_{crust} + F_{basalt}(t)$$

Where:
- $R_{crust}$ = 0.71107 (continental crust)
- $R_{basalt}$ = 0.7040 ± 0.001 (Emeishan basalt)
- $F_{basalt}(t)$ = Transient pulse function (Gaussian)

---

## 3. Scenario Specifications

### 3.1 SCENARIO A: Roadian Cooling (Interval 2)
**Time Window**: 273-268 Ma (*Jinogondolella nankingensis* to *J. postserrata* zones)
**Stratigraphic Range**: DK-75 to DK-53 (~144-196 m)

| Parameter | Baseline Value | Scenario A Range | Physical Justification |
|-----------|---------------|------------------|----------------------|
| $F_{riv}$ | 47.6×10⁹ mol/yr | **28-38×10⁹** | Glacial suppression of chemical weathering (30-40% reduction) |
| $R_{riv}$ | 0.71107 | 0.7110-0.7120 | Slight increase due to reduced shield rock exposure |
| $F_{highT}$ | 8.4×10⁹ | 6-10×10⁹ | Background level (pre-Neo-Tethys expansion) |
| $F_{lowT}$ | 10×10⁹ | 8-12×10⁹ | Stable low-temperature hydrothermal input |

**Target**: Decline from ~0.70765 to ~0.70700 (gradual, ~0.00013/Myr)

---

### 3.2 SCENARIO B: Wordian Rifting (Interval 3)
**Time Window**: 268-265 Ma (*J. altudaensis-shannoni* to *J. prexuanhanensis* zones)
**Stratigraphic Range**: DK-52 to DK-36 (~197-221 m)

| Parameter | Baseline Value | Scenario B Range | Physical Justification |
|-----------|---------------|------------------|----------------------|
| $F_{riv}$ | 40×10⁹ | 45-55×10⁹ | Recovery from glaciation, stable weathering |
| $R_{riv}$ | 0.71107 | **Fixed at 0.71107** | No change in weathered lithology |
| $F_{highT}$ | 8.4×10⁹ | **15-25×10⁹** | Neo-Tethys MOR expansion (Merdith et al., 2021) |
| $F_{lowT}$ | 10×10⁹ | **13-18×10⁹** | Enhanced off-axis hydrothermal circulation |

**Target**: Maintain low values ~0.70690-0.70700 (plateau)

---

### 3.3 SCENARIO C: Capitanian LIP (Interval 4)
**Time Window**: 265-259.1 Ma (*J. prexuanhanensis* to *J. xuanhanensis* zones, pre-GLB)
**Stratigraphic Range**: DK-35 to DK-30 (~223-232 m)

#### Special Parameterization:

**A. Basalt Weathering Pulse Function**
$$F_{basalt}(t) = F_{basalt}^{max} \cdot \exp\left(-\frac{(t-263.5)^2}{2 \cdot \sigma^2}\right) \cdot H(t-265.5)$$

Where:
- $t_{onset}$ = 265.5 Ma (Emeishan eruption start)
- $t_{peak}$ = 263.5 Ma (maximum weathering intensity)
- $\sigma$ = 1.5 Myr (pulse duration)
- $F_{basalt}^{max}$ = **15-30×10⁹ mol/yr** (30-60% of total riverine flux at peak)

**B. Enhanced Total Weathering**
$$F_{riv}(t) = [F_{crust} + F_{basalt}(t)] \cdot \alpha_{climate}$$
- $F_{crust}$ = 35×10⁹ mol/yr (background)
- $\alpha_{climate}$ = 1.2-1.8 (greenhouse warming enhancement)

**C. Isotopic Composition**
- $R_{basalt}$ = **0.7040 ± 0.0005** (unradiogenic end-member)
- $R_{riv}^{effective}$ = 0.7055-0.7080 (mixture)

**D. Hydrothermal Constraint**
- $F_{highT}$ = 6-10×10⁹ (declining MOR activity post-Wordian peak)

**Target**: Sharp decline to **0.70682** (Paleozoic minimum) followed by rapid rebound at 259.1 Ma

---

## 4. Input Data Requirements

### 4.1 Dukou Section Sr Isotope Data
**Source**: `Supplementary Table S1.xlsx`

| Column | Data Type | Processing Required |
|--------|-----------|---------------------|
| Sample ID | String | Extract formation info (DK prefix) |
| Position (m) | Float | Convert to Age (Ma) using bio-stratigraphy |
| ⁸⁷Sr/⁸⁶Sr | Float (6 decimals) | Filter by diagenetic criteria (Mn/Sr &lt; 1) |
| Formation | Categorical | Group by Chihsia/Maokou/Wuchiaping |

### 4.2 Age Model Construction
**Biostratigraphic Anchors** (from manuscript Fig. 1C, 4):

| Sample | Position (m) | Conodont Zone | Age (Ma) | Reference |
|--------|-------------|---------------|----------|-----------|
| DK-85 | 120 | Pre-Roadian | 283.0 | Base Maokou |
| DK-75 | 144 | *J. nankingensis* | 272.95 | Roadian base |
| DK-53 | 196 | *J. postserrata* | 268.0 | Wordian base |
| DK-36 | 221 | *J. prexuanhanensis* | 265.1 | Mid-Capitanian |
| DK-30 | 232 | *J. xuanhanensis* | 259.1 | GLB (top Maokou) |
| DK-18 | 258 | Wuchiapingian | 254.14 | Wuchiaping base |

**Interpolation Method**: Linear interpolation between biozone boundaries, assuming constant sedimentation rates within zones.

### 4.3 Target Curve Generation
1. Filter raw data: Mn/Sr &lt; 1.0 (from `Supplementary Table S2.xlsx`)
2. LOWESS smoothing (frac=0.3) to generate $R_{target}(t)$
3. Calculate 95% confidence intervals for model filtering

---

## 5. Implementation Roadmap

### Phase 1: Infrastructure (Week 1)
- [ ] Implement base `SrBoxModel` class with mass balance solver
- [ ] Implement `ParameterGenerator` with scenario-specific logic
- [ ] Load and clean Dukou section data (depth-to-age conversion)
- [ ] Create `LOWESSSmoother` for target curve generation

### Phase 2: Individual Scenarios (Week 2)
**Priority Order**: C → B → A (Start with most complex Emeishan mechanism)

- [ ] **Scenario C (Emeishan)**:
  - Implement Gaussian pulse function for $F_{basalt}(t)$
  - Implement dual-endmember mixing ($R_{riv}$ calculation)
  - Run 5,000 iterations, filter to match 0.70682 minimum
  - Generate density heatmaps ($F_{basalt}$ vs $R_{basalt}$)

- [ ] **Scenario B (Neo-Tethys)**:
  - Fix $R_{riv}$ = 0.71107, vary $F_{highT}$ and $F_{lowT}$
  - Test if hydrothermal enhancement alone maintains low ⁸⁷Sr/⁸⁶Sr plateau

- [ ] **Scenario A (Glaciation)**:
  - Implement weathering suppression function ($F_{riv}$ reduction)
  - Couple with inverse relationship between $F_{riv}$ and $R_{riv}$

### Phase 3: Integration & Validation (Week 3)
- [ ] **Sequential Run**: Chain scenarios A→B→C (time-continuous boundary conditions)
- [ ] **Sensitivity Analysis**:
  - Test Decline 3 without basalt end-member (fail test)
  - Test Decline 3 without flux enhancement (fail test)
  - Prove both $R_{riv}$↓ and $F_{riv}$↑ are required
- [ ] **Bayesian Comparison**: BIC scores for single-mechanism vs. coupled-mechanism models

### Phase 4: Visualization (Week 4)
- [ ] **Fig. 7 Equivalent**: Multi-panel proxy correlation
  - Panel A: Global Sr curve context
  - Panel B: Dukou data with 3 declines highlighted
  - Panel C: Model $F_{riv}$ vs CO₂ curve (Foster et al., 2017)
  - Panel D: Model $F_{highT}$ vs MOR length (Merdith et al., 2021)
  - Panel E: $F_{basalt}$ pulse timing vs Emeishan stratigraphy
- [ ] **Fig. 8 Equivalent**: Density heatmaps (Wang et al. 2021 style)
- [ ] **Supplementary**: All 5,000 iteration traces with confidence envelopes

---

## 6. Code Architecture

```python
# Core Classes Structure

class MiddlePermianSrModel:
    ├── __init__(age_array, target_sr)
    ├── scenarios/
    │   ├── ScenarioA_Glaciation()      # Roadian cooling
    │   ├── ScenarioB_Rifting()         # Neo-Tethys opening
    │   └── ScenarioC_EmeishanLIP()     # Basalt weathering pulse (PRIORITY)
    ├── physics/
    │   ├── mass_balance_solver()       # Eq. 1 implementation
    │   ├── basalt_pulse_function()     # Gaussian pulse for Decline 3
    │   └── weathering_feedback()       # Climate coupling
    ├── montecarlo/
    │   ├── step1_uniform_sampling()    # 5,000 iterations
    │   ├── step2_gaussian_refinement() # 2,000 iterations
    │   └── filter_solutions()          # Chi-square test
    └── visualization/
        ├── plot_three_declines()       # Fig. 7 equivalent
        ├── plot_density_heatmaps()     # Wang 2021 style
        └── plot_flux_evolution()       # Time-series of $F_{riv}$, $F_{highT}$, $F_{basalt}$

7. Success Criteria

7.1 Quantitative Metrics

Goodness-of-fit: Reduced χ² < 2.0 for all three intervals
Minimum value: Model must reach ⁸⁷Sr/⁸⁶Sr ≤ 0.70685 in Interval 4
Rate matching:
Decline 1 slope: ~0.00013/Myr
Decline 3 slope: ~0.000045/Myr (accelerated)
7.2 Qualitative Validation

Scenario A: Must show F
riv
​
   reduction correlates with CO₂ decline (glaciation)
Scenario B: Must show F
highT
​
   increase correlates with MOR expansion (independent of sea-level fall)
Scenario C: Must show f
basalt
​
   peak at ~263.5 Ma, trailing Emeishan eruption by ~2 Myr (weathering lag)

8. Key References for Parameterization

| Parameter         | Modern Value    | Permian Range     | Source                            |
| ----------------- | --------------- | ----------------- | --------------------------------- |
| $F_{riv}$         | 47.6×10⁹ mol/yr | 10-190×10⁹        | Peucker-Ehrenbrink & Fiske (2019) |
| $R_{riv}$ (crust) | 0.71107         | 0.710-0.712       | Wang et al. (2021)                |
| $R_{basalt}$      | N/A             | **0.7038-0.7050** | Huang et al. (2019, IGR)          |
| $F_{highT}$       | 8.4×10⁹         | 2-35×10⁹          | Li & Elderfield (2013)            |
| Emeishan onset    | N/A             | **265.5±0.5 Ma**  | Huang et al. (2022, Geology)      |
| MOR length peak   | N/A             | **Wordian**       | Merdith et al. (2021)             |

9. Risk Mitigation

| Risk                              | Mitigation Strategy                                                             |
| --------------------------------- | ------------------------------------------------------------------------------- |
| **Age model uncertainty**         | Test sensitivity to ±1 Myr shifts in biozone boundaries                         |
| **Basalt end-member uncertainty** | Run ensemble with $R_{basalt}$ = 0.7035, 0.7040, 0.7045                         |
| **Overfitting**                   | Use BIC to penalize excessive parameter freedom                                 |
| **Non-uniqueness**                | Explicitly test alternative hypotheses (e.g., Decline 3 by MOR only) and reject |
ß