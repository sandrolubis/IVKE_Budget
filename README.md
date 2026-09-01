

# Vapor Kinetic Energy (VKE) and Integrated VKE (IVKE) Budget

**Sandro W. Lubis, Ph.D.**
Pacific Northwest National Laboratory (PNNL)

This repository provides Python code for calculating Vapor Kinetic Energy (VKE), Integrated Vapor Kinetic Energy (IVKE), and individual VKE/IVKE tendency terms from pressure-level atmospheric data.

The budget formulation implemented in this code follows **Eq. (3) of Lubis et al. (2026)** and is based on the Vapor Kinetic Energy framework introduced by **Ong and Yang (2024)**.

### Application

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*
**Geophysical Research Letters.**

The VKE/IVKE budget used in that study is given in **Eq. (3) of Lubis et al. (2026)**.

### Original VKE framework

Eq. (3) of **Ong, H., & Yang, D. (2024).**
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*
**Nature Communications, 15**, 9428.
https://doi.org/10.1038/s41467-024-53369-0

The decomposition of additional physical tendency terms follows the more complete formulation given in **Eq. (10) of Ong and Yang (2024)**.

---

1. Vapor Kinetic Energy

The horizontal kinetic energy per unit mass is

$$
K = \frac{1}{2}\left(u^2+v^2\right),
$$

where

$u$ = zonal wind,
$v$ = meridional wind.

Vapor Kinetic Energy is defined as

$$
\mathrm{VKE} = q^2 K,
$$

where $q$ is specific humidity.

The vertically integrated Vapor Kinetic Energy is

\frac{1}{g}
\int_{p_T}^{p_B}
q^2 K,dp,
$$

where

$g$ = gravitational acceleration,
$p_B$ = lower pressure boundary,
$p_T$ = upper pressure boundary.

The default vertical integration range in this code is 1000–200 hPa.

2. VKE Budget

The budget implemented here corresponds to Eq. (3) of Lubis et al. (2026) and follows the VKE framework of Ong and Yang (2024).

The resolved dynamical part of the VKE tendency can be written schematically as

-\mathbf{V}\cdot\nabla_p(q^2K)
-\omega\frac{\partial(q^2K)}{\partial p}
-q^2\mathbf{V}\cdot\nabla_p\Phi
+\mathrm{Other},
$$

where

$$
\mathbf{V}=(u,v).
$$

The present code explicitly calculates the first three terms on the right-hand side.

These correspond to:

Horizontal advection of VKE
Vertical advection of VKE
Potential-energy-to-kinetic-energy conversion
3. First RHS Term: Horizontal Advection of VKE

The first term on the RHS of Eq. (3) of Lubis et al. (2026) is

$$
-\mathbf{V}\cdot\nabla_p(q^2K).
$$

Using the product rule,

-q^2\mathbf{V}\cdot\nabla_pK

2qK\mathbf{V}\cdot\nabla_pq.
$$

The code separates this term into two components.

HAKE: Horizontal Advection of Kinetic Energy

-q^2
\left(
u\frac{\partial K}{\partial x}
+
v\frac{\partial K}{\partial y}
\right).
$$

HAV: Horizontal Advection of Water Vapor

-2qK
\left(
u\frac{\partial q}{\partial x}
+
v\frac{\partial q}{\partial y}
\right).
$$

Therefore, the complete horizontal VKE-advection term is

\mathrm{HAKE}+\mathrm{HAV}
}
$$

For the vertically integrated budget,

\mathrm{IHAKE}+\mathrm{IHAV}
}
$$

Thus,

Eq. (3), first RHS term:

Pressure level:
    HAKE + HAV

Vertically integrated:
    IHAKE + IHAV
4. Second RHS Term: Vertical Advection of VKE

The second term on the RHS of Eq. (3) of Lubis et al. (2026) is

$$
-\omega
\frac{\partial(q^2K)}{\partial p}.
$$

Using the product rule,

-q^2\omega\frac{\partial K}{\partial p}

2qK\omega\frac{\partial q}{\partial p}.
$$

Because

$$
K=\frac{1}{2}(u^2+v^2),
$$

we have

u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}.
$$

The code therefore separates the vertical VKE-advection term into the following two components.

VAKE: Vertical Advection of Kinetic Energy

-q^2\omega
\left(
u\frac{\partial u}{\partial p}
+
v\frac{\partial v}{\partial p}
\right).
$$

VAV: Vertical Advection of Water Vapor

-2qK\omega
\frac{\partial q}{\partial p}.
$$

Therefore,

\mathrm{VAKE}+\mathrm{VAV}
}
$$

For the vertically integrated budget,

\mathrm{IVAKE}+\mathrm{IVAV}
}
$$

Thus,

Eq. (3), second RHS term:

Pressure level:
    VAKE + VAV

Vertically integrated:
    IVAKE + IVAV
5. Third RHS Term: Potential Energy to Kinetic Energy Conversion

The third term on the RHS of Eq. (3) of Lubis et al. (2026) is

$$
-q^2\mathbf{V}\cdot\nabla_p\Phi.
$$

The code calculates this term directly as PEKE:

-q^2
\left(
u\frac{\partial\Phi}{\partial x}
+
v\frac{\partial\Phi}{\partial y}
\right).
$$

Therefore,

\mathrm{PEKE}
}
$$

For the vertically integrated budget,

\mathrm{IPEKE}
}
$$

This term represents the contribution associated with potential-energy-to-kinetic-energy conversion, weighted by $q^2$.

Thus,

Eq. (3), third RHS term:

Pressure level:
    PEKE

Vertically integrated:
    IPEKE



    
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
