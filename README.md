

# Vapor Kinetic Energy (VKE) and Integrated VKE (IVKE) Budget

**Sandro W. Lubis, Ph.D.**
Pacific Northwest National Laboratory (PNNL)

This repository provides Python code for calculating **Vapor Kinetic Energy (VKE)**, **Integrated Vapor Kinetic Energy (IVKE)**, **Integrated Vapor Transport (IVT)**, and individual VKE/IVKE tendency terms from pressure-level atmospheric data.

The budget formulation implemented in this code follows **Eq. (3) of Lubis et al. (2026)** and is based on the Vapor Kinetic Energy framework introduced by **Ong and Yang (2024)**.

### Primary application

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*
**Geophysical Research Letters.**

The VKE/IVKE budget used in that study is given in **Eq. (3) of Lubis et al. (2026)**.

### Original VKE framework

**Ong, H., & Yang, D. (2024).**
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*
**Nature Communications, 15**, 9428.
https://doi.org/10.1038/s41467-024-53369-0

The decomposition of additional physical tendency terms follows the more complete formulation given in **Eq. (10) of Ong and Yang (2024)**.

---

## VKE Definition

Horizontal kinetic energy is defined as

$$
K=\frac{1}{2}(u^2+v^2),
$$

and Vapor Kinetic Energy is

$$
\mathrm{VKE}=q^2K,
$$

where

* \(q\) = specific humidity,
* \(u\) = zonal wind,
* \(v\) = meridional wind.

The vertically integrated quantity is

$$
\mathrm{IVKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}q^2K\,dp.
$$

The default vertical integration range in this code is **1000–200 hPa**.

---

# VKE Budget

The budget implemented here corresponds to **Eq. (3) of Lubis et al. (2026)** and the corresponding VKE formulation of Ong and Yang (2024).

The resolved dynamical part of the VKE tendency can be written schematically as

$$
\frac{\partial (q^2K)}{\partial t}
=
-\mathbf{V}\cdot\nabla_p(q^2K)
-\omega\frac{\partial(q^2K)}{\partial p}
-q^2\mathbf{V}\cdot\nabla_p\Phi
+\mathrm{Other},
$$

where \(\mathbf{V}=(u,v)\).

The first three terms on the right-hand side are explicitly calculated by this code.

---

## 1. Horizontal Advection of VKE

The **first term on the RHS of Eq. (3) of Lubis et al. (2026)** is

$$
-\mathbf{V}\cdot\nabla_p(q^2K).
$$

Using the product rule,

$$
-\mathbf{V}\cdot\nabla_p(q^2K)
=
-q^2\mathbf{V}\cdot\nabla_pK
-
2qK\mathbf{V}\cdot\nabla_pq.
$$

In this code, these two components are

$$
\mathrm{HAKE}
=
-q^2
\left(
u\frac{\partial K}{\partial x}
+
v\frac{\partial K}{\partial y}
\right)
$$

and

$$
\mathrm{HAV}
=
-2qK
\left(
u\frac{\partial q}{\partial x}
+
v\frac{\partial q}{\partial y}
\right).
$$

Therefore,

$$
\boxed{
\mathrm{Horizontal\ VKE\ Advection}
=
\mathrm{HAKE}+\mathrm{HAV}
}
$$

and for the vertically integrated budget,

$$
\boxed{
\mathrm{Horizontal\ IVKE\ Advection}
=
\mathrm{IHAKE}+\mathrm{IHAV}.
}
$$

---

## 2. Vertical Advection of VKE

The **second term on the RHS of Eq. (3) of Lubis et al. (2026)** is

$$
-\omega
\frac{\partial(q^2K)}{\partial p}.
$$

Applying the product rule,

$$
-\omega
\frac{\partial(q^2K)}{\partial p}
=
-q^2\omega\frac{\partial K}{\partial p}
-
2qK\omega\frac{\partial q}{\partial p}.
$$

The code separates this into

$$
\mathrm{VAKE}
=
-q^2\omega
\left(
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}
\right)
$$

and

$$
\mathrm{VAV}
=
-2qK\omega
\frac{\partial q}{\partial p}.
$$

Therefore,

$$
\boxed{
\mathrm{Vertical\ VKE\ Advection}
=
\mathrm{VAKE}+\mathrm{VAV}
}
$$

and for IVKE,

$$
\boxed{
\mathrm{Vertical\ IVKE\ Advection}
=
\mathrm{IVAKE}+\mathrm{IVAV}.
}
$$

---

## 3. Potential Energy to Kinetic Energy Conversion

The **third term on the RHS of Eq. (3) of Lubis et al. (2026)** is

$$
-q^2\mathbf{V}\cdot\nabla_p\Phi.
$$

This term is directly output by the code as

$$
\boxed{\mathrm{PEKE}}
$$

where

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

The vertically integrated counterpart is

$$
\boxed{\mathrm{IPEKE}}.
$$

This term represents the **potential-energy-to-kinetic-energy conversion contribution to VKE**.

---

# Mapping to Eq. (3) of Lubis et al. (2026)

For the pressure-level VKE budget:

| Eq. (3) term                      | Code variable                                                    |
| --------------------------------- | ---------------------------------------------------------------- |
| 1st RHS: Horizontal VKE advection | `HAKE + HAV`                                                     |
| 2nd RHS: Vertical VKE advection   | `VAKE + VAV`                                                     |
| 3rd RHS: PE → KE conversion       | `PEKE`                                                           |
| Remaining physical processes      | Residual or calculated explicitly if tendency data are available |

For the vertically integrated IVKE budget:

| Eq. (3) term                       | Code variable                                                    |
| ---------------------------------- | ---------------------------------------------------------------- |
| 1st RHS: Horizontal IVKE advection | `IHAKE + IHAV`                                                   |
| 2nd RHS: Vertical IVKE advection   | `IVAKE + IVAV`                                                   |
| 3rd RHS: PE → KE conversion        | `IPEKE`                                                          |
| Remaining physical processes       | Residual or calculated explicitly if tendency data are available |

Therefore,

$$
\boxed{
\mathrm{VKE}_{dyn}
=
\mathrm{HAKE}
+\mathrm{HAV}
+\mathrm{VAKE}
+\mathrm{VAV}
+\mathrm{PEKE}
}
$$

and

$$
\boxed{
\mathrm{IVKE}_{dyn}
=
\mathrm{IHAKE}
+\mathrm{IHAV}
+\mathrm{IVAKE}
+\mathrm{IVAV}
+\mathrm{IPEKE}.
}
$$

---

# Residual and Additional Physical Processes

If the total VKE tendency is calculated,

$$
\frac{\partial \mathrm{VKE}}{\partial t},
$$

the remaining contribution can be diagnosed as

$$
R_{\mathrm{VKE}}
=
\frac{\partial \mathrm{VKE}}{\partial t}
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
\frac{\partial \mathrm{IVKE}}{\partial t}
-
\left[
\mathrm{IHAKE}
+\mathrm{IHAV}
+\mathrm{IVAKE}
+\mathrm{IVAV}
+\mathrm{IPEKE}
\right].
$$

These residuals can contain contributions from physical processes that are not explicitly calculated by the present code, as well as numerical and budget-closure errors.

If model or reanalysis tendency variables are available, additional terms can instead be calculated explicitly following **Eq. (10) of Ong and Yang (2024)**.

These may include contributions from:

* friction and turbulent momentum tendencies,
* moist convection,
* condensation and other moist-physics processes,
* gravity-wave drag,
* subgrid-scale processes,
* turbulent moisture tendencies,
* surface-pressure tendency effects,
* and other model-physics tendencies.

Thus, these processes do not necessarily need to be treated as a residual if the corresponding tendency variables are available.

---

# Output Variables

## Pressure-Level Output

| Variable | Description                                         | Units  |
| -------- | --------------------------------------------------- | ------ |
| `VKE`    | Vapor kinetic energy                                | m² s⁻² |
| `HAKE`   | Horizontal advection of kinetic energy contribution | m² s⁻³ |
| `HAV`    | Horizontal advection of water-vapor contribution    | m² s⁻³ |
| `VAKE`   | Vertical advection of kinetic energy contribution   | m² s⁻³ |
| `VAV`    | Vertical advection of water-vapor contribution      | m² s⁻³ |
| `PEKE`   | Potential energy → kinetic energy conversion        | m² s⁻³ |

Hence,

```text
Horizontal VKE advection = HAKE + HAV
Vertical VKE advection   = VAKE + VAV
PE → KE conversion       = PEKE
```

## Vertically Integrated Output

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

```text
Eq. (3), 1st RHS = IHAKE + IHAV
Eq. (3), 2nd RHS = IVAKE + IVAV
Eq. (3), 3rd RHS = IPEKE
```

---

# Required Input Variables

| Variable | Description                | Units   |
| -------- | -------------------------- | ------- |
| `q`      | Specific humidity          | kg kg⁻¹ |
| `u`      | Zonal wind                 | m s⁻¹   |
| `v`      | Meridional wind            | m s⁻¹   |
| `omega`  | Pressure vertical velocity | Pa s⁻¹  |
| `phi`    | Geopotential               | m² s⁻²  |
| `ps`     | Surface pressure           | Pa      |

Pressure may be supplied in **Pa or hPa**.

Pressure levels may also be either ascending or descending. The code automatically rearranges them internally into

```text
1000 → 925 → 850 → ... → 300 → 250 → 200 hPa
```

for the calculation.

If geopotential height \(Z\) in meters is used instead of geopotential, convert it first:

$$
\Phi=gZ.
$$

---

# Citation

If you use or adapt this code, please cite **both** papers.

### Application of the VKE/IVKE budget

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*
**Geophysical Research Letters.**

The VKE/IVKE budget implemented in this repository corresponds to **Eq. (3) of Lubis et al. (2026)**.

### Original VKE framework

**Ong, H., & Yang, D. (2024).**
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*
**Nature Communications, 15**, 9428.
https://doi.org/10.1038/s41467-024-53369-0

The additional physical tendency decomposition discussed in this repository follows **Eq. (10) of Ong and Yang (2024)**.

---

# Author

**Sandro W. Lubis**
Atmospheric Sciences and Global Change Division
Pacific Northwest National Laboratory (PNNL)
Richland, Washington, USA
