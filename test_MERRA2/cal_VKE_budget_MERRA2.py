#!/usr/bin/env python3
# =============================================================================
# MERRA-2 IVT / VKE / IVKE DIAGNOSTICS AND BUDGET
# =============================================================================
#
# Author:
#   Sandro W. Lubis
#
# Affiliation:
#   Pacific Northwest National Laboratory (PNNL)
#   Atmospheric Sciences and Global Change Division
#   Richland, Washington, USA
#
# Usage:
#   python cal_VKE_budget_MERRA2.py 20190101 \
#       --input-dir /path/to/MERRA2 \
#       --output-dir ./output
#
# =============================================================================

"""
MERRA-2 IVT / VKE / IVKE diagnostics and budget terms.

Description
-----------
Calculates:
  - instantaneous IVT and IVKE
  - time-mean VKE
  - IVKE tendency and surface-pressure tendency effect
  - HAKE, VAKE, PEKE
  - MOKE, TOKE, GOKE
  - DOKE/AOKE
  - HAV, VAV, HCV, VCV
  - POV, TOV, COV
  - DOV/AOV and AOVKE
  - vertically integrated counterparts

Usage
-----
python cal_VKE_budget_MERRA2.py 20190101 \\
    --input-dir /path/to/MERRA2 \\
    --output-dir ./output

Notes
-----
1. All pressure calculations are performed in Pa internally.
2. The default integration layer is 1000--200 hPa, with 150 hPa retained
   as an extra level for vertical derivatives/interpolation.
3. The calculation uses spherical finite differences, pressure-level
   interpolation, and terrain-aware vertical integration.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import xarray as xr


# =============================================================================
# USER SETTINGS
# =============================================================================

G = 9.81
EARTH_RADIUS = 6.371e6

P_BOTTOM_HPA = 1000.0
P_TOP_HPA = 200.0
P_EXTRA_HPA = 150.0

P_BOTTOM = P_BOTTOM_HPA * 100.0
P_TOP = P_TOP_HPA * 100.0
P_EXTRA = P_EXTRA_HPA * 100.0

DT_SECONDS = 10800.0  # 3 hours

# Change these templates for your local MERRA-2 archive.
QDT_TEMPLATE = "MERRA2_400.tavg3_3d_qdt_Np.{date}.nc"
UDT_TEMPLATE = "MERRA2_400.tavg3_3d_udt_Np.{date}.nc"
ASM_TEMPLATE = "MERRA2_400.inst3_3d_asm_Np.{date}.nc"
ASMA_TEMPLATE = "MERRA2_400.tavg3_3d_asm_Nv.{date}.nc"

OUTPUT_DIR = "."

# MERRA-2 variable names
VARNAMES = {
    "ps": "PS",
    "q": "QV",
    "u": "U",
    "v": "V",
    "h": "H",
    "omega": "OMEGA",
    "pl": "PL",
    "dqdtana": "DQVDTANA",
    "dqdtdyn": "DQVDTDYN",
    "dqdtmst": "DQVDTMST",
    "dqdttrb": "DQVDTTRB",
    "dqdtchm": "DQVDTCHM",
    "dudtana": "DUDTANA",
    "dudtdyn": "DUDTDYN",
    "dudtmst": "DUDTMST",
    "dudttrb": "DUDTTRB",
    "dudtgwd": "DUDTGWD",
    "dvdtana": "DVDTANA",
    "dvdtdyn": "DVDTDYN",
    "dvdtmst": "DVDTMST",
    "dvdttrb": "DVDTTRB",
    "dvdtgwd": "DVDTGWD",
}


# =============================================================================
# BASIC HELPERS
# =============================================================================

def as_float32(a):
    return np.asarray(a, dtype=np.float32)


def pressure_to_pa(p):
    p = np.asarray(p, dtype=np.float64)
    if np.nanmax(np.abs(p)) < 2000.0:
        p = p * 100.0
    return p


def ensure_4d_tlevlatlon(da: xr.DataArray) -> np.ndarray:
    """
    Return a DataArray as numpy with dimension order:
        time, lev, lat, lon
    """
    dim_map = {}
    for d in da.dims:
        dl = d.lower()
        if dl == "time":
            dim_map["time"] = d
        elif dl in ("lev", "level", "plev"):
            dim_map["lev"] = d
        elif dl in ("lat", "latitude"):
            dim_map["lat"] = d
        elif dl in ("lon", "longitude"):
            dim_map["lon"] = d

    required = ("time", "lev", "lat", "lon")
    if not all(k in dim_map for k in required):
        raise ValueError(f"Cannot identify time/lev/lat/lon dimensions in {da.dims}")

    da = da.transpose(
        dim_map["time"], dim_map["lev"], dim_map["lat"], dim_map["lon"]
    )
    return as_float32(da.values)


def ensure_3d_tlatlon(da: xr.DataArray) -> np.ndarray:
    dim_map = {}
    for d in da.dims:
        dl = d.lower()
        if dl == "time":
            dim_map["time"] = d
        elif dl in ("lat", "latitude"):
            dim_map["lat"] = d
        elif dl in ("lon", "longitude"):
            dim_map["lon"] = d

    required = ("time", "lat", "lon")
    if not all(k in dim_map for k in required):
        raise ValueError(f"Cannot identify time/lat/lon dimensions in {da.dims}")

    da = da.transpose(dim_map["time"], dim_map["lat"], dim_map["lon"])
    return as_float32(da.values)


def select_pressure_range(da, lev_name, p_bottom_hpa, p_extra_hpa):
    """Select levels between p_bottom and p_extra, independent of ordering."""
    lev = np.asarray(da[lev_name].values)
    lev_hpa = lev / 100.0 if np.nanmax(np.abs(lev)) > 2000.0 else lev
    idx = np.where((lev_hpa <= p_bottom_hpa) & (lev_hpa >= p_extra_hpa))[0]
    if idx.size == 0:
        raise ValueError(
            f"No levels found between {p_bottom_hpa} and {p_extra_hpa} hPa."
        )
    return da.isel({lev_name: idx}), lev_hpa[idx]


def identify_lev_name(da):
    for d in da.dims:
        if d.lower() in ("lev", "level", "plev"):
            return d
    raise ValueError(f"No pressure-level dimension found in {da.dims}")


def identify_lat_lon(ds):
    lat_name = next((x for x in ds.coords if x.lower() in ("lat", "latitude")), None)
    lon_name = next((x for x in ds.coords if x.lower() in ("lon", "longitude")), None)
    if lat_name is None or lon_name is None:
        raise ValueError("Could not identify latitude/longitude coordinates.")
    return lat_name, lon_name


# =============================================================================
# DIFFERENTIAL OPERATORS
# =============================================================================

def pressure_derivative_variable_p(field, pressure):
    """
    d(field)/dp along level dimension for a pressure coordinate that may vary
    with time/latitude/longitude.

    field, pressure: (time, lev, lat, lon)
    """
    f = np.asarray(field, dtype=np.float64)
    p = np.asarray(pressure, dtype=np.float64)
    out = np.full_like(f, np.nan, dtype=np.float64)

    # Interior centered differences
    denom = p[:, 2:, :, :] - p[:, :-2, :, :]
    out[:, 1:-1, :, :] = (
        f[:, 2:, :, :] - f[:, :-2, :, :]
    ) / denom

    # One-sided endpoints
    out[:, 0, :, :] = (
        f[:, 1, :, :] - f[:, 0, :, :]
    ) / (
        p[:, 1, :, :] - p[:, 0, :, :]
    )

    out[:, -1, :, :] = (
        f[:, -1, :, :] - f[:, -2, :, :]
    ) / (
        p[:, -1, :, :] - p[:, -2, :, :]
    )

    return out.astype(np.float32)


def horizontal_gradient(field, lat, lon):
    """
    Horizontal gradient on a spherical lat-lon grid.

    Returns
    -------
    dfdx, dfdy
        Derivatives with respect to physical distance [per meter].

    field shape: (..., lat, lon)
    """
    f = np.asarray(field, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)

    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    dfdx = np.full_like(f, np.nan, dtype=np.float64)
    dfdy = np.full_like(f, np.nan, dtype=np.float64)

    # Zonal derivative: periodic central difference
    if len(lon_rad) < 3:
        raise ValueError("At least three longitude points are required.")

    # Allow regular longitude spacing, which MERRA-2 uses.
    dlon = float(np.nanmedian(np.diff(lon_rad)))
    roll_p = np.roll(f, -1, axis=-1)
    roll_m = np.roll(f, 1, axis=-1)

    coslat = np.cos(lat_rad)
    shape = [1] * f.ndim
    shape[-2] = len(lat)
    coslat_b = coslat.reshape(shape)

    with np.errstate(divide="ignore", invalid="ignore"):
        dfdx = (roll_p - roll_m) / (
            2.0 * EARTH_RADIUS * coslat_b * dlon
        )

    # Meridional derivative
    if len(lat_rad) < 3:
        raise ValueError("At least three latitude points are required.")

    # Interior centered differences
    dphi = lat_rad[2:] - lat_rad[:-2]
    shp = [1] * f.ndim
    shp[-2] = len(dphi)
    dphi_b = dphi.reshape(shp)

    dfdy[..., 1:-1, :] = (
        f[..., 2:, :] - f[..., :-2, :]
    ) / (EARTH_RADIUS * dphi_b)

    # One-sided latitude endpoints
    dphi0 = lat_rad[1] - lat_rad[0]
    dphi1 = lat_rad[-1] - lat_rad[-2]
    dfdy[..., 0, :] = (
        f[..., 1, :] - f[..., 0, :]
    ) / (EARTH_RADIUS * dphi0)
    dfdy[..., -1, :] = (
        f[..., -1, :] - f[..., -2, :]
    ) / (EARTH_RADIUS * dphi1)

    # At exact poles d/dx is singular; leave as NaN.
    pole = np.abs(coslat) < 1e-10
    if np.any(pole):
        dfdx[..., pole, :] = np.nan

    return dfdx.astype(np.float32), dfdy.astype(np.float32)


def spherical_divergence(u, v, lat, lon):
    """
    Divergence of horizontal vector (u,v) on a sphere:

      div = 1/(a cos phi) d(u)/dlambda
          + 1/(a cos phi) d(v cos phi)/dphi
    """
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)

    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)
    coslat = np.cos(lat_rad)

    dlon = float(np.nanmedian(np.diff(lon_rad)))

    shape = [1] * u.ndim
    shape[-2] = len(lat)
    cos_b = coslat.reshape(shape)

    du_dlambda = (
        np.roll(u, -1, axis=-1) - np.roll(u, 1, axis=-1)
    ) / (2.0 * dlon)

    vcos = v * cos_b
    d_vcos_dphi = np.full_like(vcos, np.nan)

    dphi = lat_rad[2:] - lat_rad[:-2]
    shp = [1] * u.ndim
    shp[-2] = len(dphi)
    dphi_b = dphi.reshape(shp)

    d_vcos_dphi[..., 1:-1, :] = (
        vcos[..., 2:, :] - vcos[..., :-2, :]
    ) / dphi_b

    d_vcos_dphi[..., 0, :] = (
        vcos[..., 1, :] - vcos[..., 0, :]
    ) / (lat_rad[1] - lat_rad[0])

    d_vcos_dphi[..., -1, :] = (
        vcos[..., -1, :] - vcos[..., -2, :]
    ) / (lat_rad[-1] - lat_rad[-2])

    with np.errstate(divide="ignore", invalid="ignore"):
        div = (
            du_dlambda + d_vcos_dphi
        ) / (EARTH_RADIUS * cos_b)

    pole = np.abs(coslat) < 1e-10
    if np.any(pole):
        div[..., pole, :] = np.nan

    return div.astype(np.float32)


# =============================================================================
# VERTICAL INTERPOLATION FROM MERRA-2 MODEL LEVELS TO PRESSURE LEVELS
# =============================================================================

def interpolate_variable_pressure_to_levels(field, pressure, target_p):
    """
    Linear interpolation from model-level pressure PL(time,lev,lat,lon)
    to fixed pressure levels.

    Parameters
    ----------
    field, pressure : ndarray
        (time, lev, lat, lon)
    target_p : 1D ndarray [Pa]

    Returns
    -------
    ndarray
        (time, target_lev, lat, lon)
    """
    f = np.asarray(field, dtype=np.float64)
    pl = np.asarray(pressure, dtype=np.float64)
    target_p = np.asarray(target_p, dtype=np.float64)

    if f.shape != pl.shape:
        raise ValueError("field and pressure must have identical shapes.")

    nt, nz, ny, nx = f.shape

    # Flatten all columns.
    fcol = np.moveaxis(f, 1, -1).reshape(-1, nz)
    pcol = np.moveaxis(pl, 1, -1).reshape(-1, nz)

    # Ensure pressure decreases with level in each column.
    desc = pcol[:, 0] > pcol[:, -1]
    if not np.all(desc):
        # Reverse columns that are ascending.
        asc_idx = np.where(~desc)[0]
        pcol[asc_idx] = pcol[asc_idx, ::-1]
        fcol[asc_idx] = fcol[asc_idx, ::-1]

    ncol = fcol.shape[0]
    out = np.full((ncol, len(target_p)), np.nan, dtype=np.float64)
    rows = np.arange(ncol)

    for kk, pt in enumerate(target_p):
        valid_col = (
            np.isfinite(pcol[:, 0])
            & np.isfinite(pcol[:, -1])
            & (pt <= pcol[:, 0])
            & (pt >= pcol[:, -1])
        )
        if not np.any(valid_col):
            continue

        pc = pcol[valid_col]
        fc = fcol[valid_col]
        r = rows[valid_col]

        # Number of model levels strictly greater than target pressure.
        j = np.sum(pc > pt, axis=1)
        j = np.clip(j, 1, nz - 1)

        rr = np.arange(pc.shape[0])
        j0 = j - 1
        j1 = j

        p0 = pc[rr, j0]
        p1 = pc[rr, j1]
        f0 = fc[rr, j0]
        f1 = fc[rr, j1]

        denom = p1 - p0
        with np.errstate(divide="ignore", invalid="ignore"):
            w = (pt - p0) / denom
            val = f0 + w * (f1 - f0)

        bad = (
            ~np.isfinite(p0)
            | ~np.isfinite(p1)
            | ~np.isfinite(f0)
            | ~np.isfinite(f1)
            | (denom == 0)
        )
        val[bad] = np.nan
        out[r, kk] = val

    out = out.reshape(nt, ny, nx, len(target_p))
    out = np.moveaxis(out, -1, 1)
    return out.astype(np.float32)


# =============================================================================
# TERRAIN-AWARE VERTICAL INTEGRATION
# =============================================================================

def vertical_integral_beta(field, pressure, surface_pressure,
                           p_bottom=P_BOTTOM, p_top=P_TOP):
    """
    Terrain-aware pressure integral:

        (1/g) integral(field dp)

    between p_top and min(surface_pressure, p_bottom).

    field : (time, lev, lat, lon)
    pressure : 1D fixed pressure levels [Pa]
    surface_pressure : (time, lat, lon) [Pa]
    """
    f = np.asarray(field, dtype=np.float64)
    p = np.asarray(pressure, dtype=np.float64)
    ps = np.asarray(surface_pressure, dtype=np.float64)

    # Put levels in descending pressure order.
    if p[0] < p[-1]:
        p = p[::-1]
        f = f[:, ::-1, :, :]

    use = (p <= p_bottom) & (p >= p_top)
    p = p[use]
    f = f[:, use, :, :]

    if len(p) < 1:
        raise ValueError("No pressure levels inside integration range.")

    ps_eff = np.minimum(ps, p_bottom)
    out = np.zeros_like(ps_eff, dtype=np.float64)

    for k in range(len(p)):
        if k == 0:
            p_lower = p_bottom
        else:
            p_lower = 0.5 * (p[k - 1] + p[k])

        if k == len(p) - 1:
            p_upper = p_top
        else:
            p_upper = 0.5 * (p[k] + p[k + 1])

        dp = np.minimum(ps_eff, p_lower) - p_upper
        dp = np.maximum(dp, 0.0)

        layer = f[:, k, :, :]
        out += np.where((dp > 0.0) & np.isfinite(layer), layer * dp, 0.0)

    return (out / G).astype(np.float32)


def fill_vertical_missing(field):
    """
    Fill missing values along pressure-level dimension by 1-D linear
    interpolation/extrapolation using nearest valid endpoint.

    Intended mainly to reproduce the use of linmsg_n before SPTE.
    """
    x = np.asarray(field, dtype=np.float32).copy()
    nt, nz, ny, nx = x.shape

    cols = np.moveaxis(x, 1, -1).reshape(-1, nz)
    z = np.arange(nz)

    for i in range(cols.shape[0]):
        c = cols[i]
        good = np.isfinite(c)
        if good.sum() == 0:
            continue
        if good.sum() == 1:
            c[~good] = c[good][0]
        elif not np.all(good):
            c[~good] = np.interp(z[~good], z[good], c[good])
        cols[i] = c

    out = cols.reshape(nt, ny, nx, nz)
    out = np.moveaxis(out, -1, 1)
    return out


# =============================================================================
# INPUT READING
# =============================================================================

def open_dataset(path):
    return xr.open_dataset(path, decode_times=True, mask_and_scale=True)


def read_fixed_pressure_fields(ds, ntime=None):
    lev_name = identify_lev_name(ds[VARNAMES["q"]])
    q_da, lev_hpa = select_pressure_range(
        ds[VARNAMES["q"]], lev_name, P_BOTTOM_HPA, P_EXTRA_HPA
    )
    idx = q_da[lev_name]

    def sel4(name):
        da = ds[name].sel({lev_name: idx})
        if ntime is not None:
            da = da.isel(time=slice(0, ntime))
        return ensure_4d_tlevlatlon(da)

    q = sel4(VARNAMES["q"])
    u = sel4(VARNAMES["u"])
    v = sel4(VARNAMES["v"])
    h = sel4(VARNAMES["h"])
    omega = sel4(VARNAMES["omega"])

    ps_da = ds[VARNAMES["ps"]]
    if ntime is not None:
        ps_da = ps_da.isel(time=slice(0, ntime))
    ps = ensure_3d_tlatlon(ps_da)
    ps = pressure_to_pa(ps).astype(np.float32)

    lev_pa = pressure_to_pa(lev_hpa)

    # Enforce descending pressure.
    if lev_pa[0] < lev_pa[-1]:
        lev_pa = lev_pa[::-1]
        q = q[:, ::-1]
        u = u[:, ::-1]
        v = v[:, ::-1]
        h = h[:, ::-1]
        omega = omega[:, ::-1]

    return ps, q, u, v, h, omega, lev_pa.astype(np.float32)


def read_tendency_fields(qdt, udt, ntime, fixed_lev_da):
    """
    Read pressure-level MERRA-2 tendency fields on the same fixed levels used
    for the budget.
    """
    lev_name_q = identify_lev_name(qdt[VARNAMES["dqdtana"]])
    lev_name_u = identify_lev_name(udt[VARNAMES["dudtana"]])

    # Select by pressure values from the fixed-pressure coordinate.
    target_hpa = np.asarray(fixed_lev_da, dtype=np.float64) / 100.0

    def select_by_range(ds, var, lev_name):
        da, lev_hpa = select_pressure_range(
            ds[var], lev_name, P_BOTTOM_HPA, P_EXTRA_HPA
        )
        da = da.isel(time=slice(0, ntime))
        arr = ensure_4d_tlevlatlon(da)
        lev_pa = pressure_to_pa(lev_hpa)
        if lev_pa[0] < lev_pa[-1]:
            arr = arr[:, ::-1]
            lev_pa = lev_pa[::-1]
        if len(lev_pa) != len(target_hpa) or not np.allclose(
            lev_pa, fixed_lev_da, atol=1.0
        ):
            raise ValueError(
                f"{var}: tendency levels do not match fixed pressure levels."
            )
        return arr

    out = {}
    for key in ("dqdtana", "dqdtdyn", "dqdtmst", "dqdttrb", "dqdtchm"):
        out[key] = select_by_range(qdt, VARNAMES[key], lev_name_q)

    for key in (
        "dudtana", "dudtdyn", "dudtmst", "dudttrb", "dudtgwd",
        "dvdtana", "dvdtdyn", "dvdtmst", "dvdttrb", "dvdtgwd",
    ):
        out[key] = select_by_range(udt, VARNAMES[key], lev_name_u)

    return out


# =============================================================================
# MAIN CALCULATION
# =============================================================================

def process_date(date, input_dir=".", output_dir=OUTPUT_DIR):
    date = str(date)
    next_date = (
        np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:8]}")
        + np.timedelta64(1, "D")
    )
    next_date = str(next_date).replace("-", "")

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    qdt_path = input_dir / QDT_TEMPLATE.format(date=date)
    udt_path = input_dir / UDT_TEMPLATE.format(date=date)
    asm_path = input_dir / ASM_TEMPLATE.format(date=date)
    asm1_path = input_dir / ASM_TEMPLATE.format(date=next_date)
    asma_path = input_dir / ASMA_TEMPLATE.format(date=date)

    paths = [qdt_path, udt_path, asm_path, asm1_path, asma_path]
    for pth in paths:
        if not pth.exists():
            raise FileNotFoundError(f"Missing input file: {pth}")

    print("=" * 72)
    print(f"Processing MERRA-2 date: {date}")
    print("=" * 72)
    print("Reading input files...")

    with (
        open_dataset(qdt_path) as f_qdt,
        open_dataset(udt_path) as f_udt,
        open_dataset(asm_path) as f_asm,
        open_dataset(asm1_path) as f_asm1,
        open_dataset(asma_path) as f_asma,
    ):
        time_avg = f_qdt["time"].values
        n_tim_avg = len(time_avg)

        lat_name, lon_name = identify_lat_lon(f_asm)
        lat = as_float32(f_asm[lat_name].values)
        lon = as_float32(f_asm[lon_name].values)

        # ---------------------------------------------------------------------
        # Instantaneous fixed-pressure fields: current day + next-day records
        # ---------------------------------------------------------------------
        ps0, q0, u0, v0, h0, omega0, p = read_fixed_pressure_fields(
            f_asm, ntime=n_tim_avg
        )

        # Need at least one next-day record for forward 3-h tendency.
        ps1, q1, u1, v1, h1, omega1, p1 = read_fixed_pressure_fields(
            f_asm1, ntime=2
        )

        if not np.allclose(p, p1):
            raise ValueError("Current and next-day fixed pressure levels differ.")

        ps = np.concatenate([ps0, ps1], axis=0)
        q = np.concatenate([q0, q1], axis=0)
        u = np.concatenate([u0, u1], axis=0)
        v = np.concatenate([v0, v1], axis=0)
        h = np.concatenate([h0, h1], axis=0)
        omega = np.concatenate([omega0, omega1], axis=0)

        n_tim = n_tim_avg

        # ---------------------------------------------------------------------
        # Time-averaged model-level fields
        # ---------------------------------------------------------------------
        PL = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["pl"]].isel(time=slice(0, n_tim))
        )
        PL = pressure_to_pa(PL).astype(np.float32)

        PS = ensure_3d_tlatlon(
            f_asma[VARNAMES["ps"]].isel(time=slice(0, n_tim))
        )
        PS = pressure_to_pa(PS).astype(np.float32)

        Q0m = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["q"]].isel(time=slice(0, n_tim))
        )
        U0m = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["u"]].isel(time=slice(0, n_tim))
        )
        V0m = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["v"]].isel(time=slice(0, n_tim))
        )
        H0m = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["h"]].isel(time=slice(0, n_tim))
        )
        OMEGA0m = ensure_4d_tlevlatlon(
            f_asma[VARNAMES["omega"]].isel(time=slice(0, n_tim))
        )

        # Reverse the model-level ordering to match the expected pressure orientation.
        PL = PL[:, ::-1, :, :]
        Q0m = Q0m[:, ::-1, :, :]
        U0m = U0m[:, ::-1, :, :]
        V0m = V0m[:, ::-1, :, :]
        H0m = H0m[:, ::-1, :, :]
        OMEGA0m = OMEGA0m[:, ::-1, :, :]

        # If reversal made PL ascending, reverse back.  We require descending
        # pressure internally; the operation is applied consistently to all.
        sample = np.nanmedian(PL[:, :, PL.shape[2] // 2, PL.shape[3] // 2], axis=0)
        if sample[0] < sample[-1]:
            PL = PL[:, ::-1]
            Q0m = Q0m[:, ::-1]
            U0m = U0m[:, ::-1]
            V0m = V0m[:, ::-1]
            H0m = H0m[:, ::-1]
            OMEGA0m = OMEGA0m[:, ::-1]

        tendencies = read_tendency_fields(f_qdt, f_udt, n_tim, p)

        print("Input reading complete.")
        print(f"Fixed pressure levels [hPa]: {p / 100.0}")
        print(f"Grid: time={n_tim}, lev={len(p)}, lat={len(lat)}, lon={len(lon)}")

        # =====================================================================
        # INSTANTANEOUS IVT / IVKE
        # =====================================================================

        print("Calculating instantaneous IVT/IVKE...")

        q_inst = q[: n_tim + 1]
        u_inst = u[: n_tim + 1]
        v_inst = v[: n_tim + 1]
        ps_inst = ps[: n_tim + 1]

        q2_inst = q_inst * q_inst
        k_inst = 0.5 * (u_inst * u_inst + v_inst * v_inst)
        vke_inst = q2_inst * k_inst

        iqu = vertical_integral_beta(q_inst * u_inst, p, ps_inst)
        iqv = vertical_integral_beta(q_inst * v_inst, p, ps_inst)
        IVKE_inst = vertical_integral_beta(vke_inst, p, ps_inst)
        IVT_inst = np.sqrt(iqu * iqu + iqv * iqv).astype(np.float32)

        # =====================================================================
        # MODEL-LEVEL DERIVATIVES AND HORIZONTAL OPERATORS
        # =====================================================================

        print("Calculating model-level derivatives...")

        dUdp0 = pressure_derivative_variable_p(U0m, PL)
        dVdp0 = pressure_derivative_variable_p(V0m, PL)
        dQdp0 = pressure_derivative_variable_p(Q0m, PL)
        dQWdp0 = pressure_derivative_variable_p(Q0m * OMEGA0m, PL)

        WdUdp0 = OMEGA0m * dUdp0
        WdVdp0 = OMEGA0m * dVdp0
        WdQdp0 = OMEGA0m * dQdp0

        K0m = 0.5 * (U0m * U0m + V0m * V0m)

        gradKx0, gradKy0 = horizontal_gradient(K0m, lat, lon)
        gradQx0, gradQy0 = horizontal_gradient(Q0m, lat, lon)

        UgradQx0 = U0m * gradQx0
        VgradQy0 = V0m * gradQy0

        divQ0 = spherical_divergence(Q0m * U0m, Q0m * V0m, lat, lon)

        # =====================================================================
        # INTERPOLATE TIME-AVERAGED MODEL-LEVEL FIELDS TO FIXED PRESSURE
        # =====================================================================

        print("Interpolating model-level fields to pressure levels...")

        Q = interpolate_variable_pressure_to_levels(Q0m, PL, p)
        U = interpolate_variable_pressure_to_levels(U0m, PL, p)
        V = interpolate_variable_pressure_to_levels(V0m, PL, p)
        H = interpolate_variable_pressure_to_levels(H0m, PL, p)

        WdUdp = interpolate_variable_pressure_to_levels(WdUdp0, PL, p)
        WdVdp = interpolate_variable_pressure_to_levels(WdVdp0, PL, p)
        WdQdp = interpolate_variable_pressure_to_levels(WdQdp0, PL, p)
        dQWdp = interpolate_variable_pressure_to_levels(dQWdp0, PL, p)

        gradKx = interpolate_variable_pressure_to_levels(gradKx0, PL, p)
        gradKy = interpolate_variable_pressure_to_levels(gradKy0, PL, p)
        UgradQx = interpolate_variable_pressure_to_levels(UgradQx0, PL, p)
        VgradQy = interpolate_variable_pressure_to_levels(VgradQy0, PL, p)
        divQ = interpolate_variable_pressure_to_levels(divQ0, PL, p)

        # =====================================================================
        # PRESSURE-LEVEL VKE BUDGET
        # =====================================================================

        print("Calculating VKE tendency terms...")

        Q2 = Q * Q
        K = 0.5 * (U * U + V * V)
        VKE = Q2 * K

        Q2U = Q2 * U
        Q2V = Q2 * V
        KQ2 = 2.0 * K * Q

        VAKE = -Q2U * WdUdp - Q2V * WdVdp
        VAV = -KQ2 * WdQdp
        VCV = -KQ2 * dQWdp

        PHI = H * G
        gradPHIx, gradPHIy = horizontal_gradient(PHI, lat, lon)

        PEKE = -Q2U * gradPHIx - Q2V * gradPHIy
        HAKE = -Q2U * gradKx - Q2V * gradKy
        HAV = -KQ2 * (UgradQx + VgradQy)
        HCV = -KQ2 * divQ

        # =====================================================================
        # INSTANTANEOUS TENDENCIES
        # =====================================================================

        DQDT = (q_inst[1:n_tim + 1] - q_inst[0:n_tim]) / DT_SECONDS
        DUDT = (u_inst[1:n_tim + 1] - u_inst[0:n_tim]) / DT_SECONDS
        DVDT = (v_inst[1:n_tim + 1] - v_inst[0:n_tim]) / DT_SECONDS

        VKETEND = (
            KQ2 * DQDT
            + Q2U * DUDT
            + Q2V * DVDT
        ).astype(np.float32)

        IVKETEND = (
            IVKE_inst[1:n_tim + 1] - IVKE_inst[0:n_tim]
        ) / DT_SECONDS

        pstend = (
            ps_inst[1:n_tim + 1] - ps_inst[0:n_tim]
        ) / DT_SECONDS

        vke_missing_filled = fill_vertical_missing(VKE)
        SPTE = (
            pstend * vke_missing_filled[:, 0, :, :] / G
        ).astype(np.float32)

        # =====================================================================
        # MERRA-2 TENDENCY TERMS
        # =====================================================================

        DQVDTANA = tendencies["dqdtana"]
        DQVDTDYN = tendencies["dqdtdyn"]
        DQVDTMST = tendencies["dqdtmst"]
        DQVDTTRB = tendencies["dqdttrb"]
        DQVDTCHM = tendencies["dqdtchm"]

        DUDTANA = tendencies["dudtana"]
        DUDTDYN = tendencies["dudtdyn"]
        DUDTMST = tendencies["dudtmst"]
        DUDTTRB = tendencies["dudttrb"]
        DUDTGWD = tendencies["dudtgwd"]

        DVDTANA = tendencies["dvdtana"]
        DVDTDYN = tendencies["dvdtdyn"]
        DVDTMST = tendencies["dvdtmst"]
        DVDTTRB = tendencies["dvdttrb"]
        DVDTGWD = tendencies["dvdtgwd"]

        MOKE = Q2U * DUDTMST + Q2V * DVDTMST
        TOKE = Q2U * DUDTTRB + Q2V * DVDTTRB
        GOKE = Q2U * DUDTGWD + Q2V * DVDTGWD
        DOKE = Q2U * DUDTDYN + Q2V * DVDTDYN
        AOKE = Q2U * DUDTANA + Q2V * DVDTANA

        POV = KQ2 * DQVDTMST
        TOV = KQ2 * DQVDTTRB
        COV = KQ2 * DQVDTCHM
        DOV = KQ2 * DQVDTDYN
        AOV = KQ2 * DQVDTANA
        AOVKE = AOKE + AOV

        # =====================================================================
        # VERTICAL INTEGRALS
        # =====================================================================

        print("Vertically integrating VKE terms...")

        IVKEAVG = vertical_integral_beta(VKE, p, PS)
        IPEKE = vertical_integral_beta(PEKE, p, PS)

        IDOKE = vertical_integral_beta(DOKE, p, PS)
        IHAKE = vertical_integral_beta(HAKE, p, PS)
        IVAKE = vertical_integral_beta(VAKE, p, PS)

        IMOKE = vertical_integral_beta(MOKE, p, PS)
        ITOKE = vertical_integral_beta(TOKE, p, PS)
        IGOKE = vertical_integral_beta(GOKE, p, PS)

        IDOV = vertical_integral_beta(DOV, p, PS)
        IHAV = vertical_integral_beta(HAV, p, PS)
        IVAV = vertical_integral_beta(VAV, p, PS)
        IHCV = vertical_integral_beta(HCV, p, PS)
        IVCV = vertical_integral_beta(VCV, p, PS)

        IPOV = vertical_integral_beta(POV, p, PS)
        ITOV = vertical_integral_beta(TOV, p, PS)
        ICOV = vertical_integral_beta(COV, p, PS)
        IAOVKE = vertical_integral_beta(AOVKE, p, PS)
        I_VKETEND = vertical_integral_beta(VKETEND, p, PS)

        # =====================================================================
        # WRITE OUTPUT
        # =====================================================================

        print("Writing output files...")

        coords_2d_inst = {
            "time": f_asm["time"].values[:n_tim],
            "lat": lat,
            "lon": lon,
        }

        ds_inst = xr.Dataset(
            data_vars={
                "IVKE": (
                    ("time", "lat", "lon"),
                    IVKE_inst[:n_tim],
                    {
                        "long_name": "Integrated vapor kinetic energy",
                        "units": "kg s-2",
                    },
                ),
                "IVT": (
                    ("time", "lat", "lon"),
                    IVT_inst[:n_tim],
                    {
                        "long_name": "Integrated vapor transport",
                        "units": "kg m-1 s-1",
                    },
                ),
            },
            coords=coords_2d_inst,
        )

        inst_file = output_dir / f"analysis_AR_inst{date}.nc"
        ds_inst.to_netcdf(inst_file, encoding={
            "IVKE": {"zlib": True, "complevel": 2, "dtype": "float32"},
            "IVT": {"zlib": True, "complevel": 2, "dtype": "float32"},
        })

        coords = {
            "time": time_avg,
            "lev": p / 100.0,  # retain MERRA-style hPa coordinate
            "lat": lat,
            "lon": lon,
        }

        ds_out = xr.Dataset(
            coords=coords,
            attrs={
                "title": "MERRA-2 vapor kinetic energy and IVKE budget",
                "source": "Python implementation of the MERRA-2 VKE/IVKE budget",
                "integration_range": f"{P_BOTTOM_HPA:.0f}-{P_TOP_HPA:.0f} hPa",
            },
        )

        def add_2d(name, data, long_name):
            ds_out[name] = xr.DataArray(
                np.asarray(data, dtype=np.float32),
                dims=("time", "lat", "lon"),
                attrs={"long_name": long_name, "units": "kg s-3"},
            )

        def add_3d(name, data, long_name, units="m2 s-3"):
            ds_out[name] = xr.DataArray(
                np.asarray(data, dtype=np.float32),
                dims=("time", "lev", "lat", "lon"),
                attrs={"long_name": long_name, "units": units},
            )

        add_3d("VKE", VKE, "Vapor kinetic energy", "m2 s-2")
        ds_out["IVKEAVG"] = xr.DataArray(
            IVKEAVG.astype(np.float32),
            dims=("time", "lat", "lon"),
            attrs={
                "long_name": "Integrated vapor kinetic energy",
                "units": "kg s-2",
            },
        )

        add_2d("IVKETEND", IVKETEND, "IVKE tendency")
        add_2d("SPTE", SPTE, "Surface pressure tendency effect")
        add_2d("I_VKETEND", I_VKETEND, "Integral of VKE tendency")

        add_2d("IDOKE", IDOKE, "IVKE tendency due to dynamics on KE")
        add_2d("IHAKE", IHAKE, "IVKE tendency due to horizontal advection of KE")
        add_2d("IVAKE", IVAKE, "IVKE tendency due to vertical advection of KE")
        add_2d("IPEKE", IPEKE, "IVKE tendency due to potential energy conversion to KE")
        add_2d("IMOKE", IMOKE, "IVKE tendency due to moist convection on KE")
        add_2d("ITOKE", ITOKE, "IVKE tendency due to turbulence on KE")
        add_2d("IGOKE", IGOKE, "IVKE tendency due to gravity wave drag on KE")

        add_2d("IDOV", IDOV, "IVKE tendency due to dynamics on vapor")
        add_2d("IHAV", IHAV, "IVKE tendency due to horizontal advection of vapor")
        add_2d("IVAV", IVAV, "IVKE tendency due to vertical advection of vapor")
        add_2d("IHCV", IHCV, "IVKE tendency due to horizontal convergence of vapor transport")
        add_2d("IVCV", IVCV, "IVKE tendency due to vertical convergence of vapor transport")
        add_2d("IPOV", IPOV, "IVKE tendency due to precipitation on vapor")
        add_2d("ITOV", ITOV, "IVKE tendency due to turbulence on vapor")
        add_2d("ICOV", ICOV, "IVKE tendency due to chemistry on vapor")
        add_2d("IAOVKE", IAOVKE, "IVKE tendency due to analysis on VKE")

        add_3d("VKETEND", VKETEND, "Vapor kinetic energy tendency")
        add_3d("DOKE", DOKE, "VKE tendency due to dynamics on KE")
        add_3d("HAKE", HAKE, "VKE tendency due to horizontal advection of KE")
        add_3d("VAKE", VAKE, "VKE tendency due to vertical advection of KE")
        add_3d("PEKE", PEKE, "VKE tendency due to potential energy conversion to KE")
        add_3d("MOKE", MOKE, "VKE tendency due to moist convection on KE")
        add_3d("TOKE", TOKE, "VKE tendency due to turbulence on KE")
        add_3d("GOKE", GOKE, "VKE tendency due to gravity wave drag on KE")

        add_3d("DOV", DOV, "VKE tendency due to dynamics on vapor")
        add_3d("HAV", HAV, "VKE tendency due to horizontal advection of vapor")
        add_3d("VAV", VAV, "VKE tendency due to vertical advection of vapor")
        add_3d("HCV", HCV, "VKE tendency due to horizontal convergence of vapor transport")
        add_3d("VCV", VCV, "VKE tendency due to vertical convergence of vapor transport")
        add_3d("POV", POV, "VKE tendency due to precipitation on vapor")
        add_3d("TOV", TOV, "VKE tendency due to turbulence on vapor")
        add_3d("COV", COV, "VKE tendency due to chemistry on vapor")
        add_3d("AOVKE", AOVKE, "VKE tendency due to analysis")

        tavg_file = output_dir / f"analysis_AR_tavg{date}.nc"

        encoding = {}
        for name, da in ds_out.data_vars.items():
            encoding[name] = {
                "zlib": True,
                "complevel": 2,
                "shuffle": True,
                "dtype": "float32",
            }

        ds_out.to_netcdf(tavg_file, encoding=encoding)

        print("Finished.")
        print(f"Instantaneous output: {inst_file}")
        print(f"Budget output       : {tavg_file}")


# =============================================================================
# COMMAND LINE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Calculate MERRA-2 IVT/VKE/IVKE diagnostics and budget terms."
    )
    parser.add_argument(
        "date",
        help="Date in YYYYMMDD format, e.g. 20190101",
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Directory containing the MERRA-2 input files.",
    )
    parser.add_argument(
        "--output-dir",
        default=OUTPUT_DIR,
        help="Directory for output NetCDF files.",
    )
    args = parser.parse_args()

    process_date(
        date=args.date,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
