# Vapor Kinetic Energy (VKE) and Integrated Vapor Kinetic Energy (IVKE) Budget

This repository provides Python code for calculating **Vapor Kinetic Energy (VKE)**, **Integrated Vapor Kinetic Energy (IVKE)**, **Integrated Vapor Transport (IVT)**, and individual VKE/IVKE tendency terms from pressure-level atmospheric data.

The VKE budget implemented in this code follows **Eq. (3) of Lubis et al. (2026)** and is based on the VKE framework introduced by **Ong and Yang (2024)**.

---

## References

### Application of the VKE/IVKE budget

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*
**Geophysical Research Letters.**

The VKE/IVKE budget implemented in this repository corresponds to **Eq. (3) of Lubis et al. (2026)**.

### Original VKE framework

**Ong, H., and D. Yang (2024).**
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*
**Nature Communications, 15**, 9428.
doi:10.1038/s41467-024-53369-0

The complete VKE tendency decomposition, including additional parameterized physical processes, is given in **Eq. (10) of Ong and Yang (2024)**.

---

# 1. Definition of VKE

Following Ong and Yang (2024), the horizontal wind vector is

$$
\mathbf{u}=(u,v),
$$

and the horizontal kinetic energy per unit mass is

$$
K=\frac{1}{2}|\mathbf{u}|^2
=\frac{1}{2}(u^2+v^2).
$$

Vapor kinetic energy is defined as

$$
\mathrm{VKE}=q^2K,
$$

where $q$ is specific humidity.

The vertically integrated vapor kinetic energy is

$$
\mathrm{IVKE}
\equiv
-\frac{1}{g}
\int_{p_B}^{p_T}
q^2K\,dp,
$$

where

* $p_B$ is the lower pressure boundary,
* $p_T$ is the upper pressure boundary,
* $g$ is gravitational acceleration.

Because $p_B>p_T$, this is equivalent to

$$
\mathrm{IVKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}
q^2K\,dp.
$$

The default integration range used in this code is **1000–200 hPa**.

---

# 2. VKE Tendency Equation

Following the notation of Ong and Yang (2024), the simplified VKE tendency equation can be written as

$$
\frac{\partial q^2K}{\partial t}
=
-\mathbf{u}\cdot\nabla_p(q^2K)
-\omega\frac{\partial q^2K}{\partial p}
-q^2\mathbf{u}\cdot\nabla_p\Phi
+q^2\mathbf{u}\cdot\mathbf{F}_T
+2KqS_M
+\mathrm{Other}.
$$

Here,

* $\mathbf{u}$ is the horizontal wind vector,
* $\nabla_p$ is the horizontal gradient operator on a constant-pressure surface,
* $\omega=Dp/Dt$ is pressure vertical velocity,
* $\Phi$ is geopotential,
* $\mathbf{F}_T$ represents the turbulent/frictional momentum tendency,
* $S_M$ represents the apparent moisture source or sink associated with moist processes.

The present code explicitly calculates the first three dynamical terms on the right-hand side:

1. horizontal advection of VKE,
2. vertical advection of VKE,
3. potential-energy-to-kinetic-energy conversion.

---

# 3. Horizontal Advection of VKE

The first term on the RHS of **Eq. (3) of Lubis et al. (2026)** is

$$
-\mathbf{u}\cdot\nabla_p(q^2K).
$$

Applying the product rule,

$$
-\mathbf{u}\cdot\nabla_p(q^2K)
=
-q^2\mathbf{u}\cdot\nabla_pK
-
2Kq\,\mathbf{u}\cdot\nabla_pq.
$$

The code calculates these two components separately.

## HAKE

The kinetic-energy component is

$$
\mathrm{HAKE}
=
-q^2\mathbf{u}\cdot\nabla_pK.
$$

In component form,

$$
\mathrm{HAKE}
=
-q^2
\left(
u\frac{\partial K}{\partial x}
+
v\frac{\partial K}{\partial y}
\right).
$$

## HAV

The water-vapor component is

$$
\mathrm{HAV}
=
-2Kq\,\mathbf{u}\cdot\nabla_pq.
$$

In component form,

$$
\mathrm{HAV}
=
-2Kq
\left(
u\frac{\partial q}{\partial x}
+
v\frac{\partial q}{\partial y}
\right).
$$

Therefore,

$$
-\mathbf{u}\cdot\nabla_p(q^2K)
=
\mathrm{HAKE}+\mathrm{HAV}.
$$

Thus, the **first RHS term of Eq. (3)** is

```text
Pressure level:
    HAKE + HAV

Vertically integrated:
    IHAKE + IHAV
```

For IVKE,

$$
-\frac{1}{g}
\int_{p_B}^{p_T}
\left[
-\mathbf{u}\cdot\nabla_p(q^2K)
\right]dp
=
\mathrm{IHAKE}+\mathrm{IHAV}.
$$

---

# 4. Vertical Advection of VKE

The second term on the RHS of **Eq. (3) of Lubis et al. (2026)** is

$$
-\omega
\frac{\partial q^2K}{\partial p}.
$$

Applying the product rule,

$$
-\omega
\frac{\partial q^2K}{\partial p}
=
-q^2\omega\frac{\partial K}{\partial p}
-
2Kq\omega\frac{\partial q}{\partial p}.
$$

The code again calculates the two contributions separately.

## VAKE

The kinetic-energy component is

$$
\mathrm{VAKE}
=
-q^2\omega\frac{\partial K}{\partial p}.
$$

Since

$$
K=\frac{1}{2}(u^2+v^2),
$$

then

$$
\frac{\partial K}{\partial p}
=
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}.
$$

Therefore,

$$
\mathrm{VAKE}
=
-q^2\omega
\left(
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}
\right).
$$

## VAV

The water-vapor component is

$$
\mathrm{VAV}
=
-2Kq\omega
\frac{\partial q}{\partial p}.
$$

Therefore,

$$
-\omega
\frac{\partial q^2K}{\partial p}
=
\mathrm{VAKE}+\mathrm{VAV}.
$$

Thus, the **second RHS term of Eq. (3)** is

```text
Pressure level:
    VAKE + VAV

Vertically integrated:
    IVAKE + IVAV
```

For IVKE,

$$
-\frac{1}{g}
\int_{p_B}^{p_T}
\left[
-\omega
\frac{\partial q^2K}{\partial p}
\right]dp
=
\mathrm{IVAKE}+\mathrm{IVAV}.
$$

---

# 5. Potential Energy to Kinetic Energy Conversion

The third term on the RHS of **Eq. (3) of Lubis et al. (2026)** is

$$
-q^2\mathbf{u}\cdot\nabla_p\Phi.
$$

The code calculates this term as

$$
\mathrm{PEKE}
=
-q^2\mathbf{u}\cdot\nabla_p\Phi.
$$

In component form,

$$
\mathrm{PEKE}
=
-q^2
\left(
u\frac{\partial\Phi}{\partial x}
+
v\frac{\partial\Phi}{\partial y}
\right).
$$

This represents the contribution from **potential-energy-to-kinetic-energy conversion** to the VKE tendency.

Thus, the **third RHS term of Eq. (3)** is

```text
Pressure level:
    PEKE

Vertically integrated:
    IPEKE
```

For IVKE,

$$
\mathrm{IPEKE}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\left(
-q^2\mathbf{u}\cdot\nabla_p\Phi
\right)dp.
$$

---

# 6. Mapping of Eq. (3) to the Code Output

The relationship between the terms in Eq. (3) and the variables produced by the code is

| Physical process         | VKE equation                        | Pressure-level output | Integrated output |
| ------------------------ | ----------------------------------- | --------------------- | ----------------- |
| Horizontal VKE advection | $-\mathbf{u}\cdot\nabla_p(q^2K)$    | `HAKE + HAV`          | `IHAKE + IHAV`    |
| Vertical VKE advection   | $-\omega,\partial(q^2K)/\partial p$ | `VAKE + VAV`          | `IVAKE + IVAV`    |
| PE-to-KE conversion      | $-q^2\mathbf{u}\cdot\nabla_p\Phi$   | `PEKE`                | `IPEKE`           |

Therefore,

```text
Eq. (3), first RHS term  = HAKE + HAV
Eq. (3), second RHS term = VAKE + VAV
Eq. (3), third RHS term  = PEKE
```

and for the vertically integrated budget,

```text
Eq. (3), first RHS term  = IHAKE + IHAV
Eq. (3), second RHS term = IVAKE + IVAV
Eq. (3), third RHS term  = IPEKE
```

The sum of the resolved dynamical VKE tendency terms calculated by this code is

$$
\mathrm{VKE}_{\mathrm{dyn}}
=
\mathrm{HAKE}
+\mathrm{HAV}
+\mathrm{VAKE}
+\mathrm{VAV}
+\mathrm{PEKE}.
$$

For IVKE,

$$
\mathrm{IVKE}_{\mathrm{dyn}}
=
\mathrm{IHAKE}
+\mathrm{IHAV}
+\mathrm{IVAKE}
+\mathrm{IVAV}
+\mathrm{IPEKE}.
$$

---

# 7. Other Terms and Residual

The full VKE budget contains additional physical processes that are not explicitly calculated by the present version of the code.

In the notation of Ong and Yang (2024), these include terms such as

$$
q^2\mathbf{u}\cdot\mathbf{F}_T
$$

associated with turbulent/frictional momentum tendencies, and

$$
2KqS_M
$$

associated with moist-physics sources and sinks of water vapor.

The more complete formulation in **Eq. (10) of Ong and Yang (2024)** further separates momentum and moisture tendency contributions associated with parameterized physical processes.

For example, the momentum tendencies can be represented schematically by

$$
\mathbf{F}_M,\qquad
\mathbf{F}_T,\qquad
\mathbf{F}_G,
$$

where the subscripts denote tendencies associated with processes such as moist convection, turbulence, and gravity-wave drag.

Likewise, moisture sources and sinks can be represented by

$$
S_M,\qquad
S_T,\qquad
S_C,
$$

for moist-physics, turbulence, and other moisture tendencies.

If these model or reanalysis tendency variables are available, their VKE contributions can be calculated explicitly following **Eq. (10) of Ong and Yang (2024)**.

Otherwise, the uncalculated terms may be estimated as a residual.

If the total VKE tendency is available,

$$
\frac{\partial q^2K}{\partial t},
$$

then

$$
R_{\mathrm{VKE}}
=
\frac{\partial q^2K}{\partial t}
-
\left[
\mathrm{HAKE}
+\mathrm{HAV}
+\mathrm{VAKE}
+\mathrm{VAV}
+\mathrm{PEKE}
\right].
$$

Similarly,

$$
R_{\mathrm{IVKE}}
=
\frac{\partial\mathrm{IVKE}}{\partial t}
-
\left[
\mathrm{IHAKE}
+\mathrm{IHAV}
+\mathrm{IVAKE}
+\mathrm{IVAV}
+\mathrm{IPEKE}
\right].
$$

The residual should not automatically be interpreted as one specific physical process. It may include contributions from

* friction,
* turbulence,
* moist convection,
* condensation and evaporation,
* gravity-wave drag,
* other subgrid-scale processes,
* surface-pressure tendency effects,
* analysis increments in reanalysis products,
* and numerical budget-closure errors.

Whenever the corresponding tendency variables are available, these processes should preferably be calculated explicitly rather than included in the residual.

---

# 8. Output Modes

The script provides three output options.

### Vertically integrated fields only

```python
OUTPUT_MODE = "2D"
```

### Pressure-level fields only

```python
OUTPUT_MODE = "3D"
```

### Both

```python
OUTPUT_MODE = "both"
```

---

# 9. Pressure-Level Output

Dimensions:

```text
time, plev, lat, lon
```

| Variable | Definition                         | Units  |
| -------- | ---------------------------------- | ------ |
| `VKE`    | $q^2K$                             | m² s⁻² |
| `HAKE`   | $-q^2\mathbf{u}\cdot\nabla_pK$     | m² s⁻³ |
| `HAV`    | $-2Kq,\mathbf{u}\cdot\nabla_pq$    | m² s⁻³ |
| `VAKE`   | $-q^2\omega,\partial K/\partial p$ | m² s⁻³ |
| `VAV`    | $-2Kq\omega,\partial q/\partial p$ | m² s⁻³ |
| `PEKE`   | $-q^2\mathbf{u}\cdot\nabla_p\Phi$  | m² s⁻³ |

The complete advection terms are therefore

$$
\mathrm{Horizontal\ VKE\ Advection}
=
\mathrm{HAKE}+\mathrm{HAV},
$$

and

$$
\mathrm{Vertical\ VKE\ Advection}
=
\mathrm{VAKE}+\mathrm{VAV}.
$$

---

# 10. Vertically Integrated Output

Dimensions:

```text
time, lat, lon
```

| Variable | Description                     | Units      |
| -------- | ------------------------------- | ---------- |
| `IVT`    | Integrated vapor transport      | kg m⁻¹ s⁻¹ |
| `IVKE`   | Integrated vapor kinetic energy | kg s⁻²     |
| `IHAKE`  | Integrated HAKE                 | kg s⁻³     |
| `IHAV`   | Integrated HAV                  | kg s⁻³     |
| `IVAKE`  | Integrated VAKE                 | kg s⁻³     |
| `IVAV`   | Integrated VAV                  | kg s⁻³     |
| `IPEKE`  | Integrated PEKE                 | kg s⁻³     |

Thus,

$$
\mathrm{Horizontal\ IVKE\ Advection}
=
\mathrm{IHAKE}+\mathrm{IHAV},
$$

$$
\mathrm{Vertical\ IVKE\ Advection}
=
\mathrm{IVAKE}+\mathrm{IVAV},
$$

and

$$
\mathrm{PE\rightarrow KE\ Conversion}
=
\mathrm{IPEKE}.
$$

---

# 11. Required Input Variables

| Variable | Description                                | Units   |
| -------- | ------------------------------------------ | ------- |
| `q`      | Specific humidity                          | kg kg⁻¹ |
| `u`      | Zonal wind                                 | m s⁻¹   |
| `v`      | Meridional wind                            | m s⁻¹   |
| `omega`  | Pressure vertical velocity, $\omega=Dp/Dt$ | Pa s⁻¹  |
| `phi`    | Geopotential, $\Phi$                       | m² s⁻²  |
| `ps`     | Surface pressure                           | Pa      |

The variable names in the input files can be changed through the `VARIABLES` dictionary in the script.

---

# 12. Pressure Coordinate

Pressure can be supplied in either **Pa** or **hPa**.

The script automatically converts hPa to Pa when necessary.

Pressure levels may be supplied in either ascending or descending order. For example, both

```text
1000, 925, 850, ..., 300, 250, 200 hPa
```

and

```text
200, 250, 300, ..., 850, 925, 1000 hPa
```

are acceptable.

Internally, pressure is arranged from high pressure to low pressure:

```text
1000 → 925 → 850 → ... → 300 → 250 → 200 hPa
```

The default pressure settings are

```python
P_BOTTOM = 100000.0   # 1000 hPa
P_TOP    =  20000.0   #  200 hPa
P_EXTRA  =  15000.0   #  150 hPa
```

`P_EXTRA` is used only to improve vertical derivatives near the upper integration boundary and is not included in the final vertical integration.

---

# 13. Geopotential

The variable $\Phi$ must be **geopotential** with units

```text
m² s⁻²
```

If the available input is geopotential height $Z$ in meters, convert it using

$$
\Phi=gZ.
$$

For example,

```python
phi = 9.81 * geopotential_height
```

Geopotential height in meters should not be used directly in the PE-to-KE conversion term.

---

# 14. Integrated Vapor Transport

The code also calculates IVT.

Following pressure-coordinate notation, the vertically integrated zonal and meridional vapor transports are

$$
Q_u
=
-\frac{1}{g}
\int_{p_B}^{p_T}
qu\,dp,
$$

and

$$
Q_v
=
-\frac{1}{g}
\int_{p_B}^{p_T}
qv\,dp.
$$

The IVT magnitude is

$$
\mathrm{IVT}
=
\sqrt{Q_u^2+Q_v^2}.
$$

Its units are

```text
kg m⁻¹ s⁻¹
```

---

# 15. Vertical Integration

Following the pressure-coordinate convention of Ong and Yang (2024), the vertically integrated form of a pressure-level quantity $X$ is

$$
IX
\equiv
-\frac{1}{g}
\int_{p_B}^{p_T}
X\,dp.
$$

Equivalently,

$$
IX
=
\frac{1}{g}
\int_{p_T}^{p_B}
X\,dp.
$$

For example,

$$
\mathrm{IHAKE}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\mathrm{HAKE}\,dp,
$$

$$
\mathrm{IHAV}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\mathrm{HAV}\,dp,
$$

$$
\mathrm{IVAKE}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\mathrm{VAKE}\,dp,
$$

$$
\mathrm{IVAV}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\mathrm{VAV}\,dp,
$$

and

$$
\mathrm{IPEKE}
=
-\frac{1}{g}
\int_{p_B}^{p_T}
\mathrm{PEKE}\,dp.
$$

Surface pressure is used to account for terrain and to exclude portions of pressure layers that lie below the local surface.

---

# 16. Quick Interpretation

The correspondence between **Eq. (3) of Lubis et al. (2026)** and the code output is

```text
                 VKE tendency
                      |
       --------------------------------
       |               |              |
   Horizontal       Vertical        PE → KE
   advection        advection      conversion
       |               |              |
  HAKE + HAV      VAKE + VAV         PEKE
```

For the vertically integrated budget:

```text
                IVKE tendency
                      |
       --------------------------------
       |               |              |
   Horizontal       Vertical        PE → KE
   advection        advection      conversion
       |               |              |
 IHAKE + IHAV    IVAKE + IVAV        IPEKE
```

The key correspondence is therefore

```text
Eq. (3), 1st RHS = IHAKE + IHAV
Eq. (3), 2nd RHS = IVAKE + IVAV
Eq. (3), 3rd RHS = IPEKE
```

---

# 17. Citation

If you use or adapt this code, please cite both the application paper and the original VKE framework.

### Lubis et al. (2026)

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*
**Geophysical Research Letters.**

The VKE/IVKE budget used in this repository corresponds to **Eq. (3) of Lubis et al. (2026)**.

### Ong and Yang (2024)

**Ong, H., and D. Yang (2024).**
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*
**Nature Communications, 15**, 9428.
doi:10.1038/s41467-024-53369-0

The VKE framework and the complete physical tendency decomposition are described in that paper, particularly **Eqs. (3) and (10)**.

---

# 18. Author

**Sandro W. Lubis**
Atmospheric Sciences and Global Change Division
Pacific Northwest National Laboratory (PNNL)
Richland, Washington, USA
