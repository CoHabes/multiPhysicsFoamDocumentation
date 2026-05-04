---
title: "NTGK"
mathjax: true
description: "Semi-empirical electrochemical and thermal model for lithium-ion battery cells, including a thermal abuse sub-model."
tags: ["region-type", "electrochemistry", "battery", "thermal"]
weight: 10
---

## Overview

`NTGK` implements a semi-empirical electrochemical and thermal model for
lithium-ion (Li-ion) battery cells. The acronym stands for the names of the
fitting functions used for the electrochemical conductance (NTGK) and open
circuit potential, which are polynomial functions of the depth of discharge
(DOD) with temperature corrections. The model is based on the formulations of
Kwon et al. (2006) and Kim et al. (2007).

The region type solves:

- Electric potential transport in the positive and negative electrodes
- Volumetric heat generation inside the cell (electrochemical and ohmic sources)
- Cell energy equation (temperature)
- Optional thermal abuse reactions (SEI decomposition, negative/positive solvent
  reactions, electrolyte decomposition) for over-temperature safety analysis

It is intended for use as the jelly-roll region of a battery cell, and can be
coupled to current-collector tab regions (`conductPosPotentialTemperature`,
`conductNegPotentialTemperature`) via `defaultInterface` for realistic battery
geometries.

## Governing Equations

### Electrode potential (positive electrode)

$$
\nabla \cdot \left( \sigma^+ \nabla \phi^+ \right) = \frac{1}{D^+} j
$$

### Electrode potential (negative electrode)

$$
\nabla \cdot \left( \sigma^- \nabla \phi^- \right) = -\frac{1}{D^-} j
$$

Both equations are assembled with a linearised source term (semi-implicit) to
improve convergence:

$$
\frac{j}{D^\pm} = \frac{Y}{D^\pm}\left(\phi^\pm \mp \phi^\mp \mp U_\text{OCV}\right)
$$

### Volumetric exchange current density (NTGK model)

$$
j = Y \left( \phi^+ - \phi^- - U_\text{OCV} \right)
$$

The temperature- and DOD-dependent electrochemical conductance $Y$ and open
circuit voltage $U_\text{OCV}$ are computed as 5th-order polynomial fits in DOD
with Arrhenius/linear temperature corrections:

$$
Y_0(\text{DOD}) = \sum_{i=0}^{5} a_i \,\text{DOD}^i
\qquad
Y = Y_0 \cdot \exp\!\left(-C_1\left(\frac{1}{T} - \frac{1}{T_\text{ref}}\right)\right)
$$

$$
U_0(\text{DOD}) = \sum_{i=0}^{5} b_i \,\text{DOD}^i
\qquad
U_\text{OCV} = U_0 + C_2 \left(T - T_\text{ref}\right)
$$

### Depth of discharge

DOD is updated globally (cell-average) each time step:

$$
\text{DOD}^{n+1} = \text{DOD}^n - \frac{\Delta t}{Q_\text{bat}} \sum_\text{cells} j \, V_\text{cell} / D
$$

### Energy equation

$$
\rho c_p \frac{\partial T}{\partial t} - \nabla \cdot \left( k \nabla T \right) = \dot{Q}_\text{ech} + \dot{Q}_\text{ohm} + \dot{Q}_\text{abuse}
$$

The electrochemical heat source and ohmic dissipation are:

$$
\dot{Q}_\text{ech} = \frac{j}{D}\left(\phi^+ - \phi^- - U_\text{OCV} + C_2 T\right)
$$

$$
\dot{Q}_\text{ohm} = \sigma^+ \left|\nabla\phi^+\right|^2 + \sigma^- \left|\nabla\phi^-\right|^2
$$

### Thermal abuse reactions

Four exothermic decomposition reactions are activated progressively as the cell
temperature exceeds threshold values:

| Reaction | Onset temperature | Rate law |
|----------|-------------------|----------|
| SEI decomposition | $T > 363.15\,\text{K}$ | $R_\text{SEI} = A_\text{SEI} \exp(-E_{A,\text{SEI}}/RT)\, c_\text{SEI}^{m_\text{SEI}}$ |
| Negative-electrode solvent | $T > 393.15\,\text{K}$ | $R_\text{NE} = A_\text{NE} \exp(-t_\text{SEI}/t_{\text{SEI},\text{ref}})\, c_\text{NE}^{m_\text{NE}} \exp(-E_{A,\text{NE}}/RT)$ |
| Positive-electrode solvent | $T > 393.15\,\text{K}$ | $R_\text{PE} = A_\text{PE}\, \alpha^{m_{PE,1}} (1-\alpha)^{m_{PE,2}} \exp(-E_{A,\text{PE}}/RT)$ |
| Electrolyte decomposition | $T > 473.15\,\text{K}$ | $R_\text{ELE} = A_\text{ELE} \exp(-E_{A,\text{ELE}}/RT)\, c_\text{ELE}^{m_\text{ELE}}$ |

Each reaction contributes a heat release $\dot{Q}_k = H_k W_k R_k$ to
$\dot{Q}_\text{abuse}$.

The concentrations governing these reactions are tracked by separate transport
equations (no spatial diffusion):

$$
\frac{\partial c_\text{SEI}}{\partial t} = -R_\text{SEI}, \quad
\frac{\partial c_\text{NE}}{\partial t} = -R_\text{NE}, \quad
\frac{\partial t_\text{SEI}}{\partial t} = R_\text{NE}, \quad
\frac{\partial \alpha}{\partial t} = R_\text{PE}, \quad
\frac{\partial c_\text{ELE}}{\partial t} = -R_\text{ELE}
$$

### Variable table

| Symbol | Field name | Description | Unit |
|--------|-----------|-------------|------|
| $\phi^+$ | `phiPos` | Positive electrode electric potential | V |
| $\phi^-$ | `phiNeg` | Negative electrode electric potential | V |
| $T$ | `T` | Temperature | K |
| $j$ | `j` | Volumetric exchange current density | A/m² |
| $Y$ | `Y_ech` | Electrochemical conductance | S/m² |
| $U_\text{OCV}$ | `U` | Open circuit voltage | V |
| $\text{DOD}$ | `DOD` | Depth of discharge | — |
| $\dot{Q}$ | `ST` | Total volumetric heat source | W/m³ |
| $c_\text{SEI}$ | `cSEI` | Dimensionless SEI reactant concentration | — |
| $c_\text{NE}$ | `cNE` | Dimensionless negative-electrode reactant | — |
| $t_\text{SEI}$ | `tSEI` | Dimensionless SEI layer thickness | — |
| $\alpha$ | `alpha` | Degree of conversion (positive electrode) | — |
| $c_\text{ELE}$ | `cELE` | Dimensionless electrolyte concentration | — |

## Numerical Algorithm

The solution proceeds as follows each time step:

1. **`correct()`** (called at initialisation and after each solve):
   - Update DOD from the current exchange current density $j$.
   - Evaluate $Y(\text{DOD}, T)$ and $U_\text{OCV}(\text{DOD}, T)$.
   - Update $j = Y(\phi^+ - \phi^- - U_\text{OCV})$.
   - Compute heat sources $\dot{Q}_\text{ech}$ and $\dot{Q}_\text{ohm}$.
   - Compute thermal abuse reaction rates and heat sources.

2. **`solveRegion()`**: Integrate the five abuse-species ODEs explicitly
   ($c_\text{SEI}$, $c_\text{NE}$, $t_\text{SEI}$, $\alpha$, $c_\text{ELE}$).

3. **`setCoupledEqns()`**: Assemble the FVM matrices for $\phi^+$, $\phi^-$,
   and $T$ into the global block system for monolithic coupling with adjacent
   tab regions.

4. **`postSolve()`**: Calls `correct()` again to update all derived fields
   after the potentials and temperature have been updated by the coupled solve.

## Required Fields

| Field | Type | Description | Unit |
|-------|------|-------------|------|
| `phiPos` | `volScalarField` | Positive electrode potential | V |
| `phiNeg` | `volScalarField` | Negative electrode potential | V |
| `T` | `volScalarField` | Temperature | K |

## Required Dictionaries

```
constant/<region>/transportProperties      – rho, cp, k, sigmaPos, sigmaNeg
constant/<region>/electrochemicalProperties – geometry (D, Dp, Dn, Dech),
                                              fitting coefficients (a0–a5, b0–b5,
                                              C1, C2, TRef), QBat, DOD
constant/<region>/thermalAbuseProperties   – Arrhenius parameters (A*, EA*,
                                              m*, H*, W*) and initial species
                                              concentrations (cSEI, cNE, tSEI,
                                              alpha, cELE)
```

## Case Setup

### `multiRegionProperties` entry (standalone battery cell)

```
regions
(
    ( batteryRegion (NTGK) )
);
```

### `multiRegionProperties` entry (battery with current-collector tabs)

```
regions
(
    ( batteryRegion (NTGK) )
    ( negTabRegion  (conductNegPotentialTemperature) )
    ( posTabRegion  (conductPosPotentialTemperature) )
);

DNA
{
    phiPos { maxCoupleIter 50; residualControl { maxJumpRes 1e-5; maxFluxRes 1e-5; } }
    phiNeg { maxCoupleIter 25; residualControl { maxJumpRes 1e-5; maxFluxRes 1e-5; } }
    T      { maxCoupleIter 10; residualControl { maxJumpRes 1e-5; maxFluxRes 1e-5; } }
}
```

### `controlDict` (singleRegionMode)

```
singleRegionMode    yes;
regionType          NTGK;

endTime             3700;   // typical: one-hour discharge (3600 s) + margin
deltaT              1;      // 1 s time step is typical
```

## Coupling

| Interface type | Coupled region types | Exchanged fields | Coupling mode |
|----------------|---------------------|-----------------|---------------|
| `defaultInterface` | `conductPosPotentialTemperature` | `phiPos`, `T` | partitioned (DNA) |
| `defaultInterface` | `conductNegPotentialTemperature` | `phiNeg`, `T` | partitioned (DNA) |

The electrode potential and temperature boundary conditions at the tab interfaces
are enforced via `genericRegionCoupledJump` / `genericRegionCoupledFlux` patch
fields on the battery side, matching the normal current flux and temperature
continuity across the jelly-roll / tab interface.

## Known Limitations

- **No spatial charge transport in electrolyte** — the model uses a
  volume-averaged (lumped) DOD rather than a spatially resolved lithium
  concentration in the electrolyte. Internal concentration gradients within the
  cell thickness are not captured.
- **Homogeneous jelly-roll** — the electrode stack is treated as a single
  homogeneous material with effective properties ($\rho$, $c_p$, $k$). Layered
  electrode/separator geometry is not resolved.
- **Polynomial fit validity** — the $Y$ and $U_\text{OCV}$ fits are only valid
  over the DOD and temperature range used for calibration. Extrapolation outside
  the calibration window will produce unphysical results.
- **Thermal abuse thresholds hard-coded** — the onset temperatures for abuse
  reactions ($363.15\,\text{K}$, $393.15\,\text{K}$, $473.15\,\text{K}$) are
  fixed constants in the source code and cannot be changed via dictionary input.
- **Global (cell-average) DOD** — the depth of discharge is tracked as a single
  scalar value integrated over the whole region volume. Local DOD variations
  within the cell are not resolved.
- **No electrochemical–mechanical coupling** — swelling, lithiation-induced
  stress, and mechanical degradation are not modelled.

## Tutorials

| Tutorial path | Description |
|---------------|-------------|
| `tutorials/singleRegion/NTGK/cylindricalCell18650` | Standalone 18650 cylindrical cell discharge (singleRegionMode) |
| `tutorials/electrochemistry/NTGKBatteryWithTabs` | Prismatic battery cell with positive and negative current-collector tabs coupled via `defaultInterface` |

## References

- Kwon, K. H. et al., *A two-dimensional modeling of a lithium-polymer battery*,
  Journal of Power Sources 163(1), 151–157 (2006).
  DOI: [10.1016/j.jpowsour.2005.12.051](https://doi.org/10.1016/j.jpowsour.2005.12.051)

- Kim, G.-H., Pesaran, A., Spotnitz, R., *A three-dimensional thermal abuse model
  for lithium-ion cells*, Journal of Power Sources 170(2), 476–489 (2007).
  DOI: [10.1016/j.jpowsour.2007.04.018](https://doi.org/10.1016/j.jpowsour.2007.04.018)

- Alkfri, A., Habes, C., Marschall, H., *multiRegionFoam: A Unified Multiphysics
  Framework for Multi-Region Coupled Continuum-Physical Problems*,
  Engineering with Computers (2023).
  DOI: [10.1007/s00366-024-01974-4](https://doi.org/10.1007/s00366-024-01974-4)
