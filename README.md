# Vapor Kinetic Energy (VKE) Budget

**Sandro W. Lubis** Pacific Northwest National Laboratory (PNNL) (sandro.lubis@pnnl.gov)

Python code for calculating Vapor Kinetic Energy (VKE), Integrated Vapor Kinetic Energy (IVKE), Integrated Vapor Transport (IVT), and selected VKE/IVKE tendency terms from pressure-level atmospheric data.

The implementation follows **Eq. (3) of Lubis et al. (2026)** and is based on the VKE framework of **Ong and Yang (2024)** (see **Eqs. 3 and 10**).
The attached code uses ERA5 as the input dataset, but it can be easily adapted for other reanalysis products or model output.

## Citation
If you use or adapt this code, please cite:

**Lubis, S. W., L. R. Leung, and M. Battalio (2026).**  
*More Frequent Atmospheric Rivers and Associated Precipitation Extremes Induced by the Baroclinic Annular Mode.*  
**Geophysical Research Letters.**

The VKE/IVKE budget implemented here corresponds to **Eq. (3) of Lubis et al. (2026)**.

**Ong, H., and D. Yang (2024).**  
*Vapor kinetic energy for the detection and understanding of atmospheric rivers.*  
**Nature Communications, 15**, 9428.  
https://doi.org/10.1038/s41467-024-53369-0

The original VKE framework and additional physical-process decomposition are given in **Ong and Yang (2024)**, including Eq. (10).


Please let me know if you find any bugs or issues with the code.

---

## VKE Definition

Following Ong and Yang (2024),

```math
K=\frac{1}{2}\left\lVert\vec{u}\right\rVert^2
```

and

```math
\mathrm{VKE}=q^2K
```

where $\mathbf{u}=(u,v)$ is the horizontal wind vector and $q$ is specific humidity.

The vertically integrated VKE is

```math
\mathrm{IVKE}
=
\frac{1}{g}
\int_{p_T}^{p_B}q^2K\,dp
```

The default integration range is **1000–200 hPa**.

---

## VKE Budget

The resolved part of **Eq. (3) of Lubis et al. (2026)** is

```math
\frac{\partial(q^2K)}{\partial t}
=
-\mathbf{u}\cdot\nabla_p(q^2K)
-\omega\frac{\partial(q^2K)}{\partial p}
-q^2\mathbf{u}\cdot\nabla_p\Phi
+\mathrm{Other}
```

The first three RHS terms are calculated by this code.

### 1. Horizontal VKE Advection

```math
-\mathbf{u}\cdot\nabla_p(q^2K)
=
-q^2\mathbf{u}\cdot\nabla_pK
-
2Kq\,\mathbf{u}\cdot\nabla_pq
```

Therefore,

```text
Horizontal VKE advection  = HAKE + HAV
Horizontal IVKE advection = IHAKE + IHAV
```

where

```math
\mathrm{HAKE}
=
-q^2\mathbf{u}\cdot\nabla_pK
```

and

```math
\mathrm{HAV}
=
-2Kq\,\mathbf{u}\cdot\nabla_pq
```

### 2. Vertical VKE Advection

```math
-\omega\frac{\partial(q^2K)}{\partial p}
=
-q^2\omega\frac{\partial K}{\partial p}
-
2Kq\omega\frac{\partial q}{\partial p}
```

Therefore,

```text
Vertical VKE advection  = VAKE + VAV
Vertical IVKE advection = IVAKE + IVAV
```

where

```math
\mathrm{VAKE}
=
-q^2\omega\frac{\partial K}{\partial p}
```

and

```math
\mathrm{VAV}
=
-2Kq\omega\frac{\partial q}{\partial p}
```

### 3. Potential Energy to Kinetic Energy Conversion

```math
\mathrm{PEKE}
=
-q^2\mathbf{u}\cdot\nabla_p\Phi
```

Therefore,

```text
Pressure level:       PEKE
Vertically integrated: IPEKE
```

---

### 4. Other Terms / Residual

The present code does not explicitly calculate all parameterized physical-process terms.

If the total VKE tendency is available, the remaining contribution can be estimated as

```math
R_{\mathrm{VKE}}
=
\frac{\partial(q^2K)}{\partial t}
-
\left(
\mathrm{HAKE}
+\mathrm{HAV}
+\mathrm{VAKE}
+\mathrm{VAV}
+\mathrm{PEKE}
\right)
```

The residual may contain contributions from friction, turbulence, moist convection, condensation, gravity-wave drag, subgrid-scale processes, surface-pressure effects, and numerical budget-closure errors.

If the corresponding model or reanalysis tendency variables are available, these terms can instead be calculated explicitly following **Eq. (10) of Ong and Yang (2024)**.

***NOTE***:

**If you wish to calculate the “other” dissipation and physical tendency terms explicitly, see the example provided in the test_MERRA2 directory. This example uses MERRA-2 tendency fields to diagnose additional VKE/IVKE budget terms associated with processes such as turbulence, moist convection, gravity-wave drag, and other parameterized tendencies.**

---

## Required Inputs

| Variable | Description | Units |
|---|---|---|
| `q` | Specific humidity | kg kg⁻¹ |
| `u` | Zonal wind | m s⁻¹ |
| `v` | Meridional wind | m s⁻¹ |
| `omega` | Pressure vertical velocity | Pa s⁻¹ |
| `phi` | Geopotential, $\Phi$ | m² s⁻² |
| `ps` | Surface pressure | Pa |

Pressure levels may be supplied in **Pa or hPa** and in either ascending or descending order.

`phi` must be geopotential. If geopotential height $Z$ is provided instead,

```math
\Phi=gZ
```

---

## Output

Set

```python
OUTPUT_MODE = "2D"    # vertically integrated
OUTPUT_MODE = "3D"    # pressure-level
OUTPUT_MODE = "both"  # both
```

### Pressure-level output

`VKE`, `HAKE`, `HAV`, `VAKE`, `VAV`, `PEKE`

### Vertically integrated output

`IVT`, `IVKE`, `IHAKE`, `IHAV`, `IVAKE`, `IVAV`, `IPEKE`






