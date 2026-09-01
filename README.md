# Vapor Kinetic Energy (VKE) and Integrated Vapor Kinetic Energy (IVKE) Budget

This repository provides Python code for calculating **Vapor Kinetic Energy (VKE)**, **Integrated Vapor Kinetic Energy (IVKE)**, **Integrated Vapor Transport (IVT)**, and selected VKE/IVKE tendency terms from pressure-level atmospheric data.

The VKE/IVKE budget implemented in this code follows **Eq. (3) of Lubis et al. (2026)** and is based on the Vapor Kinetic Energy framework introduced by **Ong and Yang (2024)**.

The code explicitly calculates the dynamical terms associated with:

- horizontal advection of VKE,
- vertical advection of VKE,
- potential energy (PE) conversion to kinetic energy (KE).

The remaining physical tendency terms can be diagnosed as a **residual**, or calculated explicitly if the corresponding model/reanalysis tendency variables are available, following **Eq. (10) of Ong and Yang (2024)**.

---

## References

### Lubis et al. (2026)

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**  
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*  
**Geophysical Research Letters.**

The VKE/IVKE budget used in this repository corresponds to **Eq. (3) of Lubis et al. (2026)**.

### Ong and Yang (2024)

**Ong, H., and D. Yang (2024).**  
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*  
**Nature Communications, 15**, 9428.  
https://doi.org/10.1038/s41467-024-53369-0

The original VKE framework and its complete physical tendency decomposition are described in **Eqs. (3), (8), and (10) of Ong and Yang (2024)**.

---

# 1. Definition of Vapor Kinetic Energy

Following Ong and Yang (2024), the horizontal velocity vector is

```math
\mathbf{u}=(u,v)
```

and the horizontal kinetic energy per unit mass is

```math
K
=
\frac{1}{2}
\left(
u^2+v^2
\right)
=
\frac{1}{2}|\mathbf{u}|^2
```

where

- `u` is zonal velocity,
- `v` is meridional velocity.

Vapor kinetic energy is defined as

```math
\mathrm{VKE}
=
q^2K
```

where `q` is specific humidity.

Following Ong and Yang (2024), vertically integrated vapor kinetic energy is

```math
\mathrm{IVKE}
\equiv
-\frac{1}{g}
\int_{p_B}^{p_T}
q^2K\,dp
```

where

- `g` is gravitational acceleration,
- `p_B` is the lower pressure boundary,
- `p_T` is the upper pressure boundary.

Because `p_B > p_T`, the equivalent form used computationally is

```math
\mathrm{IVKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}
q^2K\,dp
```

The default integration range in this code is **1000–200 hPa**.

---

# 2. VKE Tendency Equation

Following the notation of Ong and Yang (2024), the simplified VKE prognostic equation is

```math
\frac{\partial q^2K}{\partial t}
=
-\mathbf{u}\cdot\boldsymbol{\nabla}_p(q^2K)
-\omega\frac{\partial q^2K}{\partial p}
-q^2\mathbf{u}\cdot\boldsymbol{\nabla}_p\Phi
+q^2\mathbf{u}\cdot\mathbf{F}_T
+2KqS_M
+\mathrm{Other}
```

where

- `q` is specific humidity,
- `K` is horizontal kinetic energy,
- `u` is the horizontal velocity vector,
- `∇p` is the horizontal gradient operator on a constant-pressure surface,
- `ω = Dp/Dt` is pressure vertical velocity,
- `Φ` is geopotential,
- `F_T` is the turbulent momentum tendency,
- `S_M` is the apparent moisture source/sink associated with moist physics.

The present code explicitly calculates the **first three terms on the RHS**.

The most important correspondence is:

```text
Eq. (3), 1st RHS term = HAKE + HAV
Eq. (3), 2nd RHS term = VAKE + VAV
Eq. (3), 3rd RHS term = PEKE
```

For the vertically integrated IVKE budget:

```text
Eq. (3), 1st RHS term = IHAKE + IHAV
Eq. (3), 2nd RHS term = IVAKE + IVAV
Eq. (3), 3rd RHS term = IPEKE
```

---

# 3. First RHS Term: Horizontal Advection of VKE

The first RHS term in Eq. (3) is

```math
-\mathbf{u}\cdot\boldsymbol{\nabla}_p(q^2K)
```

Using the product rule,

```math
-\mathbf{u}\cdot\boldsymbol{\nabla}_p(q^2K)
=
-q^2\mathbf{u}\cdot\boldsymbol{\nabla}_pK
-
2Kq\mathbf{u}\cdot\boldsymbol{\nabla}_pq
```

Following the terminology of Ong and Yang (2024), these two components are **HAKE** and **HAV**.

## 3.1 HAKE: Horizontal Advection of Kinetic Energy

```math
\mathrm{HAKE}
=
-q^2\mathbf{u}\cdot\boldsymbol{\nabla}_pK
```

In component form,

```math
\mathrm{HAKE}
=
-q^2
\left(
u\frac{\partial K}{\partial x}
+
v\frac{\partial K}{\partial y}
\right)
```

## 3.2 HAV: Horizontal Advection of Vapor

```math
\mathrm{HAV}
=
-2Kq\mathbf{u}\cdot\boldsymbol{\nabla}_pq
```

In component form,

```math
\mathrm{HAV}
=
-2Kq
\left(
u\frac{\partial q}{\partial x}
+
v\frac{\partial q}{\partial y}
\right)
```

Therefore,

```math
-\mathbf{u}\cdot\boldsymbol{\nabla}_p(q^2K)
=
\mathrm{HAKE}
+
\mathrm{HAV}
```

Thus:

```text
Horizontal VKE advection = HAKE + HAV
```

and after vertical integration,

```text
Horizontal IVKE advection = IHAKE + IHAV
```

---

# 4. Second RHS Term: Vertical Advection of VKE

The second RHS term in Eq. (3) is

```math
-\omega
\frac{\partial q^2K}{\partial p}
```

Using the product rule,

```math
-\omega
\frac{\partial q^2K}{\partial p}
=
-q^2\omega\frac{\partial K}{\partial p}
-
2Kq\omega\frac{\partial q}{\partial p}
```

These two components correspond to **VAKE** and **VAV**.

## 4.1 VAKE: Vertical Advection of Kinetic Energy

```math
\mathrm{VAKE}
=
-q^2\omega
\frac{\partial K}{\partial p}
```

Since

```math
K
=
\frac{1}{2}(u^2+v^2)
```

then

```math
\frac{\partial K}{\partial p}
=
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}
```

and therefore

```math
\mathrm{VAKE}
=
-q^2\omega
\left(
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}
\right)
```

## 4.2 VAV: Vertical Advection of Vapor

```math
\mathrm{VAV}
=
-2Kq\omega
\frac{\partial q}{\partial p}
```

Therefore,

```math
-\omega
\frac{\partial q^2K}{\partial p}
=
\mathrm{VAKE}
+
\mathrm{VAV}
```

Thus:

```text
Vertical VKE advection = VAKE + VAV
```

and after vertical integration,

```text
Vertical IVKE advection = IVAKE + IVAV
```

---

# 5. Third RHS Term: Potential Energy Conversion to Kinetic Energy

The third RHS term in Eq. (3) is

```math
-q^2
\mathbf{u}\cdot
\boldsymbol{\nabla}_p\Phi
```

Following Ong and Yang (2024), this is the **potential energy conversion to kinetic energy** term.

The code outputs this term as **PEKE**:

```math
\mathrm{PEKE}
=
-q^2
\mathbf{u}\cdot
\boldsymbol{\nabla}_p\Phi
```

In component form,

```math
\mathrm{PEKE}
=
-q^2
\left(
u\frac{\partial\Phi}{\partial x}
+
v\frac{\partial\Phi}{\partial y}
\right)
```

Therefore:

```text
PE → KE conversion = PEKE
```

and after vertical integration:

```text
PE → KE conversion = IPEKE
```

---

# 6. Mapping of Eq. (3) to Code Variables

The relationship between **Eq. (3) of Lubis et al. (2026)** and the variables produced by this code is summarized below.

| Eq. (3) process | Pressure-level output | Vertically integrated output |
|---|---|---|
| Horizontal advection of VKE | `HAKE + HAV` | `IHAKE + IHAV` |
| Vertical advection of VKE | `VAKE + VAV` | `IVAKE + IVAV` |
| PE conversion to KE | `PEKE` | `IPEKE` |
| Turbulent/frictional KE tendency | Not currently calculated | Not currently calculated |
| Moist-physics/condensation tendency | Not currently calculated | Not currently calculated |
| Other physical processes | Not currently calculated | Not currently calculated |

In compact form:

```text
1st RHS term = HAKE + HAV
2nd RHS term = VAKE + VAV
3rd RHS term = PEKE
```

For IVKE:

```text
1st RHS term = IHAKE + IHAV
2nd RHS term = IVAKE + IVAV
3rd RHS term = IPEKE
```

The sum of the dynamical VKE tendency terms calculated by this code is

```math
\mathrm{VKE}_{\mathrm{dyn}}
=
\mathrm{HAKE}
+
\mathrm{HAV}
+
\mathrm{VAKE}
+
\mathrm{VAV}
+
\mathrm{PEKE}
```

For IVKE,

```math
\mathrm{IVKE}_{\mathrm{dyn}}
=
\mathrm{IHAKE}
+
\mathrm{IHAV}
+
\mathrm{IVAKE}
+
\mathrm{IVAV}
+
\mathrm{IPEKE}
```

---

# 7. Complete VKE Decomposition

The present code calculates only the resolved dynamical components required for the application in Lubis et al. (2026).

For reference, the complete VKE prognostic equation derived by Ong and Yang (2024) can be written as

```math
\frac{\partial q^2K}{\partial t}
=
-q^2\mathbf{u}\cdot\boldsymbol{\nabla}_pK
-q^2\omega\frac{\partial K}{\partial p}
-q^2\mathbf{u}\cdot\boldsymbol{\nabla}_p\Phi
+q^2\mathbf{u}\cdot\mathbf{F}_M
+q^2\mathbf{u}\cdot\mathbf{F}_T
+q^2\mathbf{u}\cdot\mathbf{F}_G
-2Kq\mathbf{u}\cdot\boldsymbol{\nabla}_pq
-2Kq\omega\frac{\partial q}{\partial p}
+2KqS_M
+2KqS_T
+2KqS_C
```

where the additional momentum tendencies are

- `F_M`: apparent momentum tendency from subgrid-scale moist convection,
- `F_T`: turbulent momentum tendency,
- `F_G`: gravity-wave-drag momentum tendency,

and the moisture tendencies are

- `S_M`: apparent moisture source/sink from moist physics,
- `S_T`: turbulent moisture tendency,
- `S_C`: chemistry-related moisture tendency.

---

# 8. Complete IVKE Terms Following Eq. (10) of Ong and Yang (2024)

Ong and Yang (2024) define the pressure integration operator as

```math
\left\langle \; \right\rangle
\equiv
-\frac{1}{g}
\int_{p_B}^{p_T}
(\;)\,dp
```

Using this notation, their complete IVKE budget contains the following terms.

## HAKE

Horizontal advection of kinetic energy:

```math
\mathrm{HAKE}
=
\left\langle
-q^2
\mathbf{u}\cdot
\boldsymbol{\nabla}_pK
\right\rangle
```

## VAKE

Vertical advection of kinetic energy:

```math
\mathrm{VAKE}
=
\left\langle
-q^2\omega
\frac{\partial K}{\partial p}
\right\rangle
```

## PEKE

Potential energy conversion to kinetic energy:

```math
\mathrm{PEKE}
=
\left\langle
-q^2
\mathbf{u}\cdot
\boldsymbol{\nabla}_p\Phi
\right\rangle
```

## MOKE

Moist-convection tendency of kinetic energy:

```math
\mathrm{MOKE}
=
\left\langle
q^2
\mathbf{u}\cdot
\mathbf{F}_M
\right\rangle
```

## TOKE

Turbulent tendency/dissipation of kinetic energy:

```math
\mathrm{TOKE}
=
\left\langle
q^2
\mathbf{u}\cdot
\mathbf{F}_T
\right\rangle
```

## GOKE

Gravity-wave-drag tendency of kinetic energy:

```math
\mathrm{GOKE}
=
\left\langle
q^2
\mathbf{u}\cdot
\mathbf{F}_G
\right\rangle
```

## HAV

Horizontal advection of vapor:

```math
\mathrm{HAV}
=
\left\langle
-2Kq
\mathbf{u}\cdot
\boldsymbol{\nabla}_pq
\right\rangle
```

## VAV

Vertical advection of vapor:

```math
\mathrm{VAV}
=
\left\langle
-2Kq\omega
\frac{\partial q}{\partial p}
\right\rangle
```

## COV

Condensation/moist-physics tendency of vapor:

```math
\mathrm{COV}
=
\left\langle
2KqS_M
\right\rangle
```

## TOV

Turbulent tendency of vapor:

```math
\mathrm{TOV}
=
\left\langle
2KqS_T
\right\rangle
```

## CMOV

Chemistry tendency of vapor:

```math
\mathrm{CMOV}
=
\left\langle
2KqS_C
\right\rangle
```

## SPTE

Surface-pressure tendency effect:

```math
\mathrm{SPTE}
=
\left[
\frac{q^2K}{g}
\right]_{p_B}
\frac{\partial p_B}{\partial t}
```

These additional terms are **not currently calculated by this code**.

However, if the corresponding model or reanalysis tendency variables are available, they can be calculated explicitly following **Eq. (10) of Ong and Yang (2024)**.

---

# 9. Residual Terms

If the total VKE tendency is available, all physical processes not explicitly calculated by the present code can be grouped into a residual.

For pressure-level VKE:

```math
R_{\mathrm{VKE}}
=
\frac{\partial q^2K}{\partial t}
-
\left(
\mathrm{HAKE}
+
\mathrm{HAV}
+
\mathrm{VAKE}
+
\mathrm{VAV}
+
\mathrm{PEKE}
\right)
```

For vertically integrated VKE:

```math
R_{\mathrm{IVKE}}
=
\frac{\partial\mathrm{IVKE}}{\partial t}
-
\left(
\mathrm{IHAKE}
+
\mathrm{IHAV}
+
\mathrm{IVAKE}
+
\mathrm{IVAV}
+
\mathrm{IPEKE}
\right)
```

The residual may contain contributions from:

- moist-convective momentum tendencies,
- turbulent/frictional momentum tendencies,
- gravity-wave drag,
- condensation and evaporation,
- subgrid-scale moist physics,
- turbulent moisture tendencies,
- chemistry or other moisture tendencies,
- surface-pressure tendency effects,
- analysis increments in reanalysis products,
- numerical discretization and budget-closure errors.

Therefore, the residual should **not be interpreted as a single physical process**.

If variables for friction, turbulence, convection, gravity-wave drag, condensation, subgrid-scale processes, or other parameterized tendencies are available, these contributions can instead be calculated explicitly following **Eq. (10) of Ong and Yang (2024)**.

---

# 10. Required Input Variables

The following variables are required by the current code.

| Variable | Description | Expected units |
|---|---|---|
| `q` | Specific humidity | kg kg⁻¹ |
| `u` | Zonal wind | m s⁻¹ |
| `v` | Meridional wind | m s⁻¹ |
| `omega` | Pressure vertical velocity | Pa s⁻¹ |
| `phi` | Geopotential | m² s⁻² |
| `ps` | Surface pressure | Pa |

The corresponding variable names in the input NetCDF files can be changed in the `VARIABLES` dictionary.

For example:

```python
VARIABLES = {
    "q":     "var133",
    "u":     "var131",
    "v":     "var132",
    "omega": "var135",
    "phi":   "var129",
    "ps":    "var151",
}
```

---

# 11. Geopotential

The variable `phi` must represent **geopotential**, denoted by `Φ`, with units

```text
m² s⁻²
```

If the available variable is geopotential height `Z` in meters, convert it first:

```math
\Phi=gZ
```

For example:

```python
phi = 9.81 * geopotential_height
```

Do **not** use geopotential height in meters directly in the PEKE calculation.

---

# 12. Pressure Coordinate

Pressure may be supplied in either **Pa** or **hPa**.

The script automatically converts pressure to Pa when necessary.

Pressure levels may also be supplied in either ascending or descending order.

For example, both

```text
1000, 925, 850, 700, ..., 300, 250, 200 hPa
```

and

```text
200, 250, 300, ..., 700, 850, 925, 1000 hPa
```

are acceptable.

Internally, the code arranges pressure from **high pressure to low pressure**:

```text
1000 → 925 → 850 → ... → 300 → 250 → 200 hPa
```

The default pressure settings are

```python
P_BOTTOM = 100000.0   # 1000 hPa
P_TOP    =  20000.0   #  200 hPa
P_EXTRA  =  15000.0   #  150 hPa
```

`P_EXTRA` is used only for calculating vertical derivatives near the upper boundary.

It is **not** included in the final pressure-level output or the 1000–200 hPa vertical integration.

---

# 13. Integrated Vapor Transport

The script also calculates Integrated Vapor Transport (IVT).

Following pressure-coordinate notation, the zonal vapor transport is

```math
Q_u
=
-\frac{1}{g}
\int_{p_B}^{p_T}
qu\,dp
```

and the meridional vapor transport is

```math
Q_v
=
-\frac{1}{g}
\int_{p_B}^{p_T}
qv\,dp
```

The IVT magnitude is then

```math
\mathrm{IVT}
=
\sqrt{
Q_u^2+Q_v^2
}
```

with units

```text
kg m⁻¹ s⁻¹
```

---

# 14. Vertical Integration

For any pressure-level quantity `X`, the vertically integrated quantity is calculated as

```math
\left\langle X \right\rangle
=
-\frac{1}{g}
\int_{p_B}^{p_T}
X\,dp
```

or equivalently,

```math
\left\langle X \right\rangle
=
\frac{1}{g}
\int_{p_T}^{p_B}
X\,dp
```

For example,

```math
\mathrm{IHAKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}
\mathrm{HAKE}\,dp
```

```math
\mathrm{IHAV}
=
\frac{1}{g}
\int_{p_T}^{p_B}
\mathrm{HAV}\,dp
```

```math
\mathrm{IVAKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}
\mathrm{VAKE}\,dp
```

```math
\mathrm{IVAV}
=
\frac{1}{g}
\int_{p_T}^{p_B}
\mathrm{VAV}\,dp
```

and

```math
\mathrm{IPEKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}
\mathrm{PEKE}\,dp
```

Surface pressure is used to account for terrain and exclude portions of pressure layers that lie below the local surface.

---

# 15. Output Modes

The code supports three output modes.

## Vertically Integrated Output Only

```python
OUTPUT_MODE = "2D"
```

This saves only vertically integrated quantities.

## Pressure-Level Output Only

```python
OUTPUT_MODE = "3D"
```

This saves only pressure-level quantities.

## Both

```python
OUTPUT_MODE = "both"
```

This saves both pressure-level and vertically integrated quantities.

---

# 16. Pressure-Level Output Variables

Dimensions:

```text
time, plev, lat, lon
```

| Variable | Description | Units |
|---|---|---|
| `VKE` | Vapor kinetic energy | m² s⁻² |
| `HAKE` | Horizontal advection of KE contribution | m² s⁻³ |
| `VAKE` | Vertical advection of KE contribution | m² s⁻³ |
| `PEKE` | PE conversion to KE contribution | m² s⁻³ |
| `HAV` | Horizontal advection of vapor contribution | m² s⁻³ |
| `VAV` | Vertical advection of vapor contribution | m² s⁻³ |

The terms corresponding to Eq. (3) are:

```text
Horizontal VKE advection = HAKE + HAV

Vertical VKE advection   = VAKE + VAV

PE → KE conversion       = PEKE
```

---

# 17. Vertically Integrated Output Variables

Dimensions:

```text
time, lat, lon
```

| Variable | Description | Units |
|---|---|---|
| `IVT` | Integrated vapor transport | kg m⁻¹ s⁻¹ |
| `IVKE` | Integrated vapor kinetic energy | kg s⁻² |
| `IHAKE` | Integrated horizontal advection of KE | kg s⁻³ |
| `IVAKE` | Integrated vertical advection of KE | kg s⁻³ |
| `IPEKE` | Integrated PE-to-KE conversion | kg s⁻³ |
| `IHAV` | Integrated horizontal advection of vapor | kg s⁻³ |
| `IVAV` | Integrated vertical advection of vapor | kg s⁻³ |

The terms corresponding to Eq. (3) are:

```text
Horizontal IVKE advection = IHAKE + IHAV

Vertical IVKE advection   = IVAKE + IVAV

PE → KE conversion        = IPEKE
```

---

# 18. Quick Summary

The relationship between the individual output variables and the VKE budget is:

```text
                         VKE tendency
                              |
         ---------------------------------------------
         |                    |                      |
     Horizontal            Vertical               PE → KE
     advection             advection             conversion
         |                    |                      |
     HAKE + HAV           VAKE + VAV                PEKE
```

For vertically integrated VKE:

```text
                        IVKE tendency
                              |
         ---------------------------------------------
         |                    |                      |
     Horizontal            Vertical               PE → KE
     advection             advection             conversion
         |                    |                      |
   IHAKE + IHAV         IVAKE + IVAV               IPEKE
```

Therefore, the key correspondence with **Eq. (3) of Lubis et al. (2026)** is:

```text
First RHS term  = IHAKE + IHAV
Second RHS term = IVAKE + IVAV
Third RHS term  = IPEKE
```

The remaining terms can be diagnosed as a residual or calculated explicitly, if the required physical tendency variables are available, following **Eq. (10) of Ong and Yang (2024)**.

---

# 19. Running the Code

Modify the user settings near the beginning of the Python script.

For example:

```python
START_YEAR = 1979
END_YEAR   = 2025

OUTPUT_MODE = "2D"

P_BOTTOM = 100000.0
P_TOP    =  20000.0
P_EXTRA  =  15000.0
```

Set the input file locations:

```python
FILES = {
    "q":     "/path/to/q.{year}.nc",
    "u":     "/path/to/u.{year}.nc",
    "v":     "/path/to/v.{year}.nc",
    "omega": "/path/to/omega.{year}.nc",
    "phi":   "/path/to/phi.{year}.nc",
    "ps":    "/path/to/ps.{year}.nc",
}
```

Then run:

```bash
python VKE_budget.py
```

Output files are written to the directory specified by

```python
OUTPUT_DIR = "./output"
```

with filenames following

```text
VKE.budget.{mode}.{year}.nc
```

---

# 20. Notes

1. The pressure coordinate can be in either Pa or hPa.

2. Pressure levels can be ascending or descending; the script automatically rearranges them internally.

3. `phi` must be geopotential in m² s⁻², not geopotential height in meters.

4. Actual surface pressure should be used whenever possible for correct terrain masking.

5. `P_EXTRA` is used only to improve the pressure derivative near the upper integration boundary.

6. The currently calculated dynamical terms do **not** represent the complete VKE budget.

7. If friction, turbulence, moist convection, gravity-wave drag, condensation, subgrid-scale tendencies, or other physical tendency variables are available, additional terms can be calculated explicitly following Eq. (10) of Ong and Yang (2024).

8. Otherwise, the contribution of unresolved processes can be estimated from the residual of the total VKE or IVKE tendency.

---

# 21. Citation

If you use or adapt this code, please cite both the application paper and the original VKE framework.

## Lubis et al. (2026)

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**  
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*  
**Geophysical Research Letters.**

The VKE/IVKE budget implemented in this repository corresponds to **Eq. (3) of Lubis et al. (2026)**.

## Ong and Yang (2024)

**Ong, H., and D. Yang (2024).**  
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*  
**Nature Communications, 15**, 9428.  
https://doi.org/10.1038/s41467-024-53369-0

The original VKE framework and the complete physical tendency decomposition are described in **Eqs. (3), (8), and (10)** of Ong and Yang (2024).

---

# 22. Author

**Sandro W. Lubis**  
Atmospheric Sciences and Global Change Division  
Pacific Northwest National Laboratory (PNNL)  
Richland, Washington, USA
