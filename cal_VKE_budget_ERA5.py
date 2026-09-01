#!/usr/bin/env python3
# =============================================================================
# VAPOR KINETIC ENERGY (VKE) AND INTEGRATED VKE (IVKE) BUDGET
# =============================================================================
#
# Author:
#   Sandro W. Lubis (Oct 2025)
#
# Affiliation:
#   Pacific Northwest National Laboratory (PNNL)
#   Atmospheric Sciences and Global Change Division
#   Richland, Washington, USA
#
# =============================================================================
# REFERENCE
# =============================================================================
#
# Lubis, S. W., L. R. Leung, and M. Battalio (2026):
# More Frequent Atmospheric Rivers and Associated Precipitation Extremes
# Induced by the Baroclinic Annular Mode.
# Geophysical Research Letters.
#
# Please cite the paper above when using or adapting this code.
#
# =============================================================================
# DESCRIPTION
# =============================================================================
#
# This script calculates pressure-level vapor kinetic energy (VKE),
# vertically integrated vapor kinetic energy (IVKE), integrated vapor
# transport (IVT), and individual VKE/IVKE tendency terms.
#
# K = 1/2 (u^2 + v^2)
#
# VKE = q^2 K
#
# VKE tendency terms:
#
#   HAKE = -q^2 [u dK/dx + v dK/dy]
#
#   VAKE = -q^2 omega [u du/dp + v dv/dp]
#
#   PEKE = -q^2 [u dPhi/dx + v dPhi/dy]
#
#   HAV  = -2 q K [u dq/dx + v dq/dy]
#
#   VAV  = -2 q K omega dq/dp
#
# Vertically integrated terms:
#
#              1
#       IX = ----- integral X dp
#              g
#
# =============================================================================
# EXPECTED INPUT VARIABLES AND UNITS
# =============================================================================
#
# Variable      Description                         Unit
# -----------------------------------------------------------------------------
# q             Specific humidity                   kg kg-1
# u             Zonal wind                          m s-1
# v             Meridional wind                     m s-1
# omega         Pressure vertical velocity          Pa s-1
# phi           Geopotential                        m2 s-2
# ps            Surface pressure                    Pa
#
# Pressure coordinate:
#
# plev          Pressure                            Pa or hPa
#
# Pressure coordinates are automatically converted to Pa when needed.
# Input pressure levels may be ascending or descending.
#
# Internally, pressure is always arranged from high pressure to low pressure:
#
#       1000 -> 925 -> 850 -> ... -> 200 -> 150 hPa
#
# IMPORTANT:
#
#   "phi" must be geopotential [m2 s-2].
#
#   If the input variable is geopotential height [m], convert it first:
#
#       phi = g * geopotential_height
#
# =============================================================================
# OUTPUT VARIABLES
# =============================================================================
#
# 3D pressure-level output:
#
# Variable      Dimensions                         Unit
# -----------------------------------------------------------------------------
# VKE           time, plev, lat, lon               m2 s-2
# HAKE          time, plev, lat, lon               m2 s-3
# VAKE          time, plev, lat, lon               m2 s-3
# PEKE          time, plev, lat, lon               m2 s-3
# HAV           time, plev, lat, lon               m2 s-3
# VAV           time, plev, lat, lon               m2 s-3
#
# 2D vertically integrated output:
#
# Variable      Dimensions                         Unit
# -----------------------------------------------------------------------------
# IVT           time, lat, lon                     kg m-1 s-1
# IVKE          time, lat, lon                     kg s-2
# IHAKE         time, lat, lon                     kg s-3
# IVAKE         time, lat, lon                     kg s-3
# IPEKE         time, lat, lon                     kg s-3
# IHAV          time, lat, lon                     kg s-3
# IVAV          time, lat, lon                     kg s-3
#
# =============================================================================

import os
import numpy as np
from netCDF4 import Dataset


# =============================================================================
# USER SETTINGS
# =============================================================================

# -----------------------------------------------------------------------------
# Analysis years
# -----------------------------------------------------------------------------

START_YEAR = 1979
END_YEAR   = 1980


# -----------------------------------------------------------------------------
# Output mode
#
# "2D"   = vertically integrated fields only
# "3D"   = pressure-level fields only
# "both" = save both
# -----------------------------------------------------------------------------

OUTPUT_MODE = "2D"


# -----------------------------------------------------------------------------
# Pressure settings
#
# All values below are in Pa.
#
# P_EXTRA is used only for vertical derivatives near P_TOP.
# It is NOT included in the layerwise output or vertical integration.
# -----------------------------------------------------------------------------

P_BOTTOM = 100000.0       # 1000 hPa
P_TOP    =  20000.0       #  200 hPa
P_EXTRA  =  15000.0       #  150 hPa


# -----------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------------------------

GRAVITY      = 9.81        # m s-2
EARTH_RADIUS = 6.371e6     # m


# -----------------------------------------------------------------------------
# Number of time steps processed simultaneously
#
# Reduce this value if memory usage is too high.
# -----------------------------------------------------------------------------

TIME_BLOCK = 4


# -----------------------------------------------------------------------------
# Coordinate names
# -----------------------------------------------------------------------------

TIME_NAME = "time"
LEV_NAME  = "plev"
LAT_NAME  = "lat"
LON_NAME  = "lon"


# -----------------------------------------------------------------------------
# Input variable names
# -----------------------------------------------------------------------------

VARIABLES = {
    "q":     "var133",
    "u":     "var131",
    "v":     "var132",
    "omega": "var135",
    "phi":   "var129",
    "ps":    "var151",
}


# -----------------------------------------------------------------------------
# Input files
#
# {year} is replaced automatically.
#
# IMPORTANT:
# The "ps" variable should ideally be actual surface pressure for proper
# terrain masking. The current path/name below follows the original dataset
# configuration. Replace it if actual surface pressure is available.
# -----------------------------------------------------------------------------

FILES = {

    "q":
        "/pscratch/sd/s/slubis/ERA5_Data/SH/daily/q.{year}.nc",

    "u":
        "/pscratch/sd/s/slubis/ERA5_Data/U/daily/u.{year}.nc",

    "v":
        "/pscratch/sd/s/slubis/ERA5_Data/V/daily/v.{year}.nc",

    "omega":
        "/pscratch/sd/s/slubis/ERA5_Data/W/daily/w.{year}.nc",

    "phi":
        "/pscratch/sd/s/slubis/ERA5_Data/Z/daily/z3.{year}.nc",

    "ps":
        "/pscratch/sd/s/slubis/ERA5_Data/SLP/daily/mslp.{year}.nc",
}


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

OUTPUT_DIR = "./output"

OUTPUT_TEMPLATE = "VKE.budget.{mode}.{year}.nc"


# =============================================================================
# OUTPUT MODE
# =============================================================================

OUTPUT_MODE = OUTPUT_MODE.upper()

if OUTPUT_MODE not in ("2D", "3D", "BOTH"):
    raise ValueError(
        'OUTPUT_MODE must be "2D", "3D", or "both".'
    )

SAVE_2D = OUTPUT_MODE in ("2D", "BOTH")
SAVE_3D = OUTPUT_MODE in ("3D", "BOTH")


# =============================================================================
# OUTPUT METADATA
# =============================================================================

LAYER_VARIABLES = {

    "VKE": (
        "Vapor kinetic energy",
        "m2 s-2"
    ),

    "HAKE": (
        "VKE tendency due to horizontal advection of kinetic energy",
        "m2 s-3"
    ),

    "VAKE": (
        "VKE tendency due to vertical advection of kinetic energy",
        "m2 s-3"
    ),

    "PEKE": (
        "VKE tendency due to potential energy conversion",
        "m2 s-3"
    ),

    "HAV": (
        "VKE tendency due to horizontal advection of water vapor",
        "m2 s-3"
    ),

    "VAV": (
        "VKE tendency due to vertical advection of water vapor",
        "m2 s-3"
    ),
}


INTEGRATED_VARIABLES = {

    "IVKE": (
        "Integrated vapor kinetic energy",
        "kg s-2"
    ),

    "IVT": (
        "Integrated vapor transport",
        "kg m-1 s-1"
    ),

    "IHAKE": (
        "IVKE tendency due to horizontal advection of kinetic energy",
        "kg s-3"
    ),

    "IVAKE": (
        "IVKE tendency due to vertical advection of kinetic energy",
        "kg s-3"
    ),

    "IPEKE": (
        "IVKE tendency due to potential energy conversion",
        "kg s-3"
    ),

    "IHAV": (
        "IVKE tendency due to horizontal advection of water vapor",
        "kg s-3"
    ),

    "IVAV": (
        "IVKE tendency due to vertical advection of water vapor",
        "kg s-3"
    ),
}


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def to_float32(data):
    """
    Convert data to float32 and replace masked values with NaN.
    """

    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)

    return np.asarray(data, dtype=np.float32)


def pressure_to_pa(pressure):
    """
    Convert pressure from hPa to Pa when necessary.
    """

    pressure = to_float32(pressure)

    if np.nanmax(pressure) < 2000.0:
        pressure = pressure * 100.0

    return pressure


# =============================================================================
# PRESSURE DERIVATIVE
# =============================================================================

def pressure_derivative(field, pressure):
    """
    Calculate d(field)/dp.

    Parameters
    ----------
    field : ndarray
        Dimensions:
        (time, plev, lat, lon)

    pressure : ndarray
        Pressure levels in Pa, ordered from high to low pressure.

    Returns
    -------
    derivative : ndarray
        Same dimensions as field.
    """

    derivative = np.empty_like(
        field,
        dtype=np.float32
    )


    # -------------------------------------------------------------------------
    # Interior levels: centered difference
    # -------------------------------------------------------------------------

    dp = (
        pressure[2:]
        - pressure[:-2]
    )

    derivative[:, 1:-1, :, :] = (
        field[:, 2:, :, :]
        - field[:, :-2, :, :]
    ) / dp[None, :, None, None]


    # -------------------------------------------------------------------------
    # Bottom boundary: one-sided difference
    # -------------------------------------------------------------------------

    derivative[:, 0, :, :] = (
        field[:, 1, :, :]
        - field[:, 0, :, :]
    ) / (
        pressure[1]
        - pressure[0]
    )


    # -------------------------------------------------------------------------
    # Top boundary: one-sided difference
    # -------------------------------------------------------------------------

    derivative[:, -1, :, :] = (
        field[:, -1, :, :]
        - field[:, -2, :, :]
    ) / (
        pressure[-1]
        - pressure[-2]
    )


    return derivative


# =============================================================================
# HORIZONTAL ADVECTION
# =============================================================================

def horizontal_advection(field, u, v, lat, lon):
    """
    Calculate:

        u d(field)/dx + v d(field)/dy

    on a spherical latitude-longitude grid.

    Longitude is treated as periodic.

    The first and last latitude rows are returned as NaN because
    centered meridional differences cannot be calculated there.
    """

    lat_rad = np.deg2rad(lat).astype(np.float32)
    lon_rad = np.deg2rad(lon).astype(np.float32)

    advection = np.full(
        field.shape,
        np.nan,
        dtype=np.float32
    )


    # -------------------------------------------------------------------------
    # Meridional gradient
    #
    # dF/dy = 1/a dF/dphi
    # -------------------------------------------------------------------------

    dlat = (
        lat_rad[2:]
        - lat_rad[:-2]
    )

    dfdy = (
        field[:, :, 2:, :]
        - field[:, :, :-2, :]
    )

    dfdy /= (
        EARTH_RADIUS
        * dlat[None, None, :, None]
    )


    # -------------------------------------------------------------------------
    # Zonal gradient
    #
    # dF/dx = 1/(a cos(phi)) dF/dlambda
    # -------------------------------------------------------------------------

    dlon = lon_rad[1] - lon_rad[0]

    field_mid = field[:, :, 1:-1, :]

    dfdx = (
        np.roll(field_mid, -1, axis=-1)
        - np.roll(field_mid, 1, axis=-1)
    )

    coslat = np.cos(
        lat_rad[1:-1]
    )

    dfdx /= (
        2.0
        * EARTH_RADIUS
        * coslat[None, None, :, None]
        * dlon
    )


    # -------------------------------------------------------------------------
    # V . grad(field)
    # -------------------------------------------------------------------------

    advection[:, :, 1:-1, :] = (
        u[:, :, 1:-1, :] * dfdx
        +
        v[:, :, 1:-1, :] * dfdy
    )


    return advection


# =============================================================================
# VERTICAL INTEGRATION
# =============================================================================

def vertical_integral(
    field,
    pressure,
    surface_pressure
):
    """
    Calculate:

             1
        ----------- integral field dp
             g

    between P_BOTTOM and P_TOP.

    Surface pressure is used to account for terrain and to exclude
    portions of pressure layers located below the surface.

    Parameters
    ----------
    field : ndarray
        Dimensions:
        (time, plev, lat, lon)

    pressure : ndarray
        Pressure levels in Pa, ordered high -> low.

    surface_pressure : ndarray
        Dimensions:
        (time, lat, lon)

    Returns
    -------
    integrated : ndarray
        Dimensions:
        (time, lat, lon)
    """

    use = (
        (pressure <= P_BOTTOM)
        &
        (pressure >= P_TOP)
    )

    p = pressure[use]
    f = field[:, use, :, :]

    nt, _, ny, nx = f.shape

    integrated = np.zeros(
        (nt, ny, nx),
        dtype=np.float32
    )


    # Surface pressure cannot exceed the lower integration boundary
    ps = np.minimum(
        surface_pressure,
        P_BOTTOM
    )


    # -------------------------------------------------------------------------
    # Integrate one pressure layer at a time
    # -------------------------------------------------------------------------

    for k in range(len(p)):


        # ---------------------------------------------------------------------
        # Lower layer boundary
        # ---------------------------------------------------------------------

        if k == 0:

            p_lower = P_BOTTOM

        else:

            p_lower = 0.5 * (
                p[k - 1]
                + p[k]
            )


        # ---------------------------------------------------------------------
        # Upper layer boundary
        # ---------------------------------------------------------------------

        if k == len(p) - 1:

            p_upper = P_TOP

        else:

            p_upper = 0.5 * (
                p[k]
                + p[k + 1]
            )


        # ---------------------------------------------------------------------
        # Actual pressure thickness above the surface
        # ---------------------------------------------------------------------

        dp = (
            np.minimum(ps, p_lower)
            - p_upper
        )

        dp = np.maximum(
            dp,
            0.0
        )


        # Avoid NaN * 0 below the surface
        integrated += np.where(
            dp > 0.0,
            f[:, k, :, :] * dp,
            0.0
        )


    return integrated / GRAVITY


# =============================================================================
# READ PRESSURE-LEVEL VARIABLE
# =============================================================================

def read_pressure_variable(
    variable,
    t0,
    t1,
    lev_indices,
    reverse_pressure
):
    """
    Read one time block and the requested pressure levels.
    """

    data = variable[
        t0:t1,
        lev_indices,
        :,
        :
    ]

    data = to_float32(data)


    # Put data in high -> low pressure order
    if reverse_pressure:
        data = data[:, ::-1, :, :]


    return data


# =============================================================================
# MASK PRESSURE-LEVEL OUTPUT BELOW SURFACE
# =============================================================================

def mask_below_surface(
    field,
    pressure,
    surface_pressure
):
    """
    Mask pressure-level values located below the local surface.

    field:
        (time, plev, lat, lon)

    pressure:
        (plev)

    surface_pressure:
        (time, lat, lon)
    """

    valid = (
        pressure[None, :, None, None]
        <= surface_pressure[:, None, :, :]
    )

    return np.where(
        valid,
        field,
        np.nan
    )


# =============================================================================
# WRITE PRESSURE-LEVEL VARIABLE
# =============================================================================

def write_layer(
    out,
    name,
    field,
    output_levels,
    pressure_output,
    surface_pressure,
    t0,
    t1
):
    """
    Select the analysis pressure range, mask below-ground values,
    and write one 3D pressure-level variable.
    """

    layer = field[
        :,
        output_levels,
        :,
        :
    ]

    layer = mask_below_surface(
        layer,
        pressure_output,
        surface_pressure
    )

    out.variables[name][
        t0:t1,
        :,
        :,
        :
    ] = np.ma.masked_invalid(
        layer.astype(np.float32)
    )


# =============================================================================
# CREATE OUTPUT FILE
# =============================================================================

def create_output_file(
    filename,
    time,
    pressure,
    lat,
    lon,
    time_source
):
    """
    Create output NetCDF file.
    """

    out = Dataset(
        filename,
        "w",
        format="NETCDF4"
    )


    # -------------------------------------------------------------------------
    # Dimensions
    # -------------------------------------------------------------------------

    out.createDimension(
        TIME_NAME,
        None
    )

    out.createDimension(
        LAT_NAME,
        len(lat)
    )

    out.createDimension(
        LON_NAME,
        len(lon)
    )


    if SAVE_3D:

        out.createDimension(
            LEV_NAME,
            len(pressure)
        )


    # -------------------------------------------------------------------------
    # Time
    # -------------------------------------------------------------------------

    tvar = out.createVariable(
        TIME_NAME,
        time_source.dtype,
        (TIME_NAME,)
    )


    for attribute in time_source.ncattrs():

        if attribute != "_FillValue":

            tvar.setncattr(
                attribute,
                time_source.getncattr(attribute)
            )


    tvar[:] = time


    # -------------------------------------------------------------------------
    # Latitude
    # -------------------------------------------------------------------------

    yvar = out.createVariable(
        LAT_NAME,
        "f4",
        (LAT_NAME,)
    )

    yvar[:] = lat

    yvar.long_name = "latitude"
    yvar.standard_name = "latitude"
    yvar.units = "degrees_north"


    # -------------------------------------------------------------------------
    # Longitude
    # -------------------------------------------------------------------------

    xvar = out.createVariable(
        LON_NAME,
        "f4",
        (LON_NAME,)
    )

    xvar[:] = lon

    xvar.long_name = "longitude"
    xvar.standard_name = "longitude"
    xvar.units = "degrees_east"


    # -------------------------------------------------------------------------
    # Pressure
    # -------------------------------------------------------------------------

    if SAVE_3D:

        pvar = out.createVariable(
            LEV_NAME,
            "f4",
            (LEV_NAME,)
        )

        pvar[:] = pressure

        pvar.long_name = "pressure"
        pvar.standard_name = "air_pressure"
        pvar.units = "Pa"
        pvar.positive = "down"
        pvar.axis = "Z"


    # -------------------------------------------------------------------------
    # 3D pressure-level variables
    # -------------------------------------------------------------------------

    if SAVE_3D:

        chunk_3d = (
            1,
            1,
            min(len(lat), 181),
            min(len(lon), 360)
        )


        for name, (long_name, units) in LAYER_VARIABLES.items():

            var = out.createVariable(
                name,
                "f4",
                (
                    TIME_NAME,
                    LEV_NAME,
                    LAT_NAME,
                    LON_NAME
                ),
                zlib=True,
                complevel=2,
                shuffle=True,
                chunksizes=chunk_3d,
                fill_value=np.float32(9.96921e36)
            )

            var.long_name = long_name
            var.units = units


    # -------------------------------------------------------------------------
    # 2D vertically integrated variables
    # -------------------------------------------------------------------------

    if SAVE_2D:

        chunk_2d = (
            1,
            min(len(lat), 181),
            min(len(lon), 360)
        )


        for name, (long_name, units) in INTEGRATED_VARIABLES.items():

            var = out.createVariable(
                name,
                "f4",
                (
                    TIME_NAME,
                    LAT_NAME,
                    LON_NAME
                ),
                zlib=True,
                complevel=2,
                shuffle=True,
                chunksizes=chunk_2d,
                fill_value=np.float32(9.96921e36)
            )

            var.long_name = long_name
            var.units = units


    # -------------------------------------------------------------------------
    # Global attributes
    # -------------------------------------------------------------------------

    out.title = (
        "Vapor kinetic energy and integrated vapor kinetic energy budget"
    )

    out.author = "Sandro W. Lubis"

    out.affiliation = (
        "Pacific Northwest National Laboratory (PNNL), "
        "Atmospheric Sciences and Global Change Division"
    )

    out.reference = (
        "Lubis, S. W., L. R. Leung, and M. Battalio (2026): "
        "More Frequent Atmospheric Rivers and Associated Precipitation "
        "Extremes Induced by the Baroclinic Annular Mode. "
        "Geophysical Research Letters."
    )

    out.output_mode = OUTPUT_MODE

    out.integration_range = (
        f"{P_BOTTOM:.0f} to {P_TOP:.0f} Pa"
    )

    out.pressure_derivative_extra_level = (
        f"{P_EXTRA:.0f} Pa"
    )


    return out


# =============================================================================
# PROCESS ONE YEAR
# =============================================================================

def process_year(year):
    """
    Calculate VKE/IVKE diagnostics for one year.
    """

    print()
    print("=" * 70)
    print(f"Processing year {year}")
    print("=" * 70)


    # -------------------------------------------------------------------------
    # Construct filenames
    # -------------------------------------------------------------------------

    filenames = {
        key: template.format(year=year)
        for key, template in FILES.items()
    }


    # -------------------------------------------------------------------------
    # Check input files
    # -------------------------------------------------------------------------

    for filename in filenames.values():

        if not os.path.exists(filename):

            raise FileNotFoundError(
                f"Missing input file:\n{filename}"
            )


    # -------------------------------------------------------------------------
    # Open input files
    # -------------------------------------------------------------------------

    datasets = {
        key: Dataset(filename, "r")
        for key, filename in filenames.items()
    }


    try:

        fq   = datasets["q"]
        fu   = datasets["u"]
        fv   = datasets["v"]
        fw   = datasets["omega"]
        fphi = datasets["phi"]
        fps  = datasets["ps"]


        # =====================================================================
        # COORDINATES
        # =====================================================================

        time = fq.variables[
            TIME_NAME
        ][:]

        lat = to_float32(
            fq.variables[LAT_NAME][:]
        )

        lon = to_float32(
            fq.variables[LON_NAME][:]
        )

        pressure_all = pressure_to_pa(
            fq.variables[LEV_NAME][:]
        )


        # =====================================================================
        # SELECT PRESSURE LEVELS
        #
        # Keep P_BOTTOM through P_EXTRA for calculation.
        # =====================================================================

        lev_indices = np.where(
            (pressure_all <= P_BOTTOM)
            &
            (pressure_all >= P_EXTRA)
        )[0]


        if len(lev_indices) < 2:

            raise ValueError(
                "Not enough pressure levels found between "
                f"{P_BOTTOM} and {P_EXTRA} Pa."
            )


        pressure = pressure_all[
            lev_indices
        ]


        # =====================================================================
        # AUTOMATIC PRESSURE ORDERING
        #
        # Internal ordering:
        #
        # high pressure -> low pressure
        # =====================================================================

        reverse_pressure = (
            pressure[0]
            < pressure[-1]
        )


        if reverse_pressure:
            pressure = pressure[::-1]


        pressure = pressure.astype(
            np.float32
        )


        # Safety check
        if not np.all(
            np.diff(pressure) < 0
        ):

            raise ValueError(
                "Pressure levels must be monotonic."
            )


        # =====================================================================
        # PRESSURE LEVELS SAVED TO OUTPUT
        #
        # Exclude P_EXTRA.
        # =====================================================================

        output_levels = (
            (pressure <= P_BOTTOM)
            &
            (pressure >= P_TOP)
        )


        pressure_output = pressure[
            output_levels
        ]


        ntime = len(time)


        print(
            "Calculation levels [hPa]:",
            pressure / 100.0
        )

        print(
            "Output levels [hPa]:",
            pressure_output / 100.0
        )

        print(
            f"Grid: time={ntime}, "
            f"level={len(pressure)}, "
            f"lat={len(lat)}, "
            f"lon={len(lon)}"
        )

        print(
            f"Output mode: {OUTPUT_MODE}"
        )


        # =====================================================================
        # CREATE OUTPUT FILE
        # =====================================================================

        os.makedirs(
            OUTPUT_DIR,
            exist_ok=True
        )


        mode_name = OUTPUT_MODE.lower()


        output_file = os.path.join(
            OUTPUT_DIR,
            OUTPUT_TEMPLATE.format(
                mode=mode_name,
                year=year
            )
        )


        out = create_output_file(
            output_file,
            time,
            pressure_output,
            lat,
            lon,
            fq.variables[TIME_NAME]
        )


        try:

            # =================================================================
            # TIME-BLOCK LOOP
            # =================================================================

            for t0 in range(
                0,
                ntime,
                TIME_BLOCK
            ):

                t1 = min(
                    t0 + TIME_BLOCK,
                    ntime
                )


                print(
                    f"  Processing time steps "
                    f"{t0}:{t1} / {ntime}"
                )


                # =============================================================
                # READ INPUT VARIABLES
                # =============================================================

                q = read_pressure_variable(
                    fq.variables[
                        VARIABLES["q"]
                    ],
                    t0,
                    t1,
                    lev_indices,
                    reverse_pressure
                )


                u = read_pressure_variable(
                    fu.variables[
                        VARIABLES["u"]
                    ],
                    t0,
                    t1,
                    lev_indices,
                    reverse_pressure
                )


                v = read_pressure_variable(
                    fv.variables[
                        VARIABLES["v"]
                    ],
                    t0,
                    t1,
                    lev_indices,
                    reverse_pressure
                )


                omega = read_pressure_variable(
                    fw.variables[
                        VARIABLES["omega"]
                    ],
                    t0,
                    t1,
                    lev_indices,
                    reverse_pressure
                )


                phi = read_pressure_variable(
                    fphi.variables[
                        VARIABLES["phi"]
                    ],
                    t0,
                    t1,
                    lev_indices,
                    reverse_pressure
                )


                surface_pressure = pressure_to_pa(
                    fps.variables[
                        VARIABLES["ps"]
                    ][
                        t0:t1,
                        :,
                        :
                    ]
                )


                # =============================================================
                # BASIC VARIABLES
                # =============================================================

                q2 = q * q

                K = 0.5 * (
                    u * u
                    +
                    v * v
                )


                # =============================================================
                # VKE
                #
                # VKE = q^2 K
                # =============================================================

                VKE = q2 * K


                if SAVE_3D:

                    write_layer(
                        out,
                        "VKE",
                        VKE,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IVKE = vertical_integral(
                        VKE,
                        pressure,
                        surface_pressure
                    )


                del VKE


                # =============================================================
                # IVT
                #
                # IVT = sqrt(IQU^2 + IQV^2)
                # =============================================================

                if SAVE_2D:

                    IQU = vertical_integral(
                        q * u,
                        pressure,
                        surface_pressure
                    )

                    IQV = vertical_integral(
                        q * v,
                        pressure,
                        surface_pressure
                    )

                    IVT = np.sqrt(
                        IQU * IQU
                        +
                        IQV * IQV
                    )

                    del IQU, IQV


                # =============================================================
                # VERTICAL ADVECTION OF KINETIC ENERGY
                #
                # VAKE =
                # -q^2 omega (u du/dp + v dv/dp)
                # =============================================================

                dudp = pressure_derivative(
                    u,
                    pressure
                )

                dvdp = pressure_derivative(
                    v,
                    pressure
                )


                VAKE = (
                    -q2
                    * omega
                    * (
                        u * dudp
                        +
                        v * dvdp
                    )
                )


                if SAVE_3D:

                    write_layer(
                        out,
                        "VAKE",
                        VAKE,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IVAKE = vertical_integral(
                        VAKE,
                        pressure,
                        surface_pressure
                    )


                del dudp, dvdp, VAKE


                # =============================================================
                # VERTICAL ADVECTION OF WATER VAPOR
                #
                # VAV = -2 q K omega dq/dp
                # =============================================================

                dqdp = pressure_derivative(
                    q,
                    pressure
                )


                VAV = (
                    -2.0
                    * q
                    * K
                    * omega
                    * dqdp
                )


                if SAVE_3D:

                    write_layer(
                        out,
                        "VAV",
                        VAV,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IVAV = vertical_integral(
                        VAV,
                        pressure,
                        surface_pressure
                    )


                del dqdp, VAV


                # =============================================================
                # HORIZONTAL ADVECTION OF KINETIC ENERGY
                #
                # HAKE = -q^2 V . grad(K)
                # =============================================================

                adv_K = horizontal_advection(
                    K,
                    u,
                    v,
                    lat,
                    lon
                )


                HAKE = (
                    -q2
                    * adv_K
                )


                if SAVE_3D:

                    write_layer(
                        out,
                        "HAKE",
                        HAKE,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IHAKE = vertical_integral(
                        HAKE,
                        pressure,
                        surface_pressure
                    )


                del adv_K, HAKE


                # =============================================================
                # POTENTIAL ENERGY CONVERSION
                #
                # PEKE = -q^2 V . grad(phi)
                # =============================================================

                adv_phi = horizontal_advection(
                    phi,
                    u,
                    v,
                    lat,
                    lon
                )


                PEKE = (
                    -q2
                    * adv_phi
                )


                if SAVE_3D:

                    write_layer(
                        out,
                        "PEKE",
                        PEKE,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IPEKE = vertical_integral(
                        PEKE,
                        pressure,
                        surface_pressure
                    )


                del adv_phi, PEKE


                # =============================================================
                # HORIZONTAL ADVECTION OF WATER VAPOR
                #
                # HAV = -2 q K V . grad(q)
                # =============================================================

                adv_q = horizontal_advection(
                    q,
                    u,
                    v,
                    lat,
                    lon
                )


                HAV = (
                    -2.0
                    * q
                    * K
                    * adv_q
                )


                if SAVE_3D:

                    write_layer(
                        out,
                        "HAV",
                        HAV,
                        output_levels,
                        pressure_output,
                        surface_pressure,
                        t0,
                        t1
                    )


                if SAVE_2D:

                    IHAV = vertical_integral(
                        HAV,
                        pressure,
                        surface_pressure
                    )


                del adv_q, HAV


                # =============================================================
                # WRITE 2D INTEGRATED VARIABLES
                # =============================================================

                if SAVE_2D:

                    out.variables["IVKE"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IVKE.astype(np.float32)
                    )


                    out.variables["IVT"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IVT.astype(np.float32)
                    )


                    out.variables["IHAKE"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IHAKE.astype(np.float32)
                    )


                    out.variables["IVAKE"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IVAKE.astype(np.float32)
                    )


                    out.variables["IPEKE"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IPEKE.astype(np.float32)
                    )


                    out.variables["IHAV"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IHAV.astype(np.float32)
                    )


                    out.variables["IVAV"][
                        t0:t1
                    ] = np.ma.masked_invalid(
                        IVAV.astype(np.float32)
                    )


                out.sync()


                # =============================================================
                # CLEAN MEMORY
                # =============================================================

                del (
                    q,
                    q2,
                    u,
                    v,
                    omega,
                    phi,
                    K,
                    surface_pressure
                )


                if SAVE_2D:

                    del (
                        IVKE,
                        IVT,
                        IHAKE,
                        IVAKE,
                        IPEKE,
                        IHAV,
                        IVAV
                    )


            print(
                f"Finished year {year}"
            )

            print(
                f"Output: {output_file}"
            )


        finally:

            out.close()


    finally:

        for dataset in datasets.values():
            dataset.close()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    for year in range(
        START_YEAR,
        END_YEAR + 1
    ):

        process_year(year)
