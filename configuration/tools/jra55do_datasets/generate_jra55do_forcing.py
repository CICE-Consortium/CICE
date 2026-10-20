#! /usr/bin/env python3

# =========================================================================================
# Regrid one year of JRA55-do (input4MIPs) atmospheric forcing onto a CICE standalone grid,
# renaming/combining fields to match the variable names expected by CICE's `ice_forcing.F90`
# JRA55(do) forcing option.
#
# To run:
#   python generate_jra55do_forcing.py --year=<year> --hgrid-filename=<path-to-supergrid-file>
#     --output-filename=<path-to-output-file>
# these extra arguments are available if required:
#     --jra55do-dir=<path-to-jra55do-source-root>
#
# For more information, run `python generate_jra55do_forcing.py -h`
#
# Method
# ------
# CICE's JRA55(do) forcing option (ice_forcing.F90:JRA55_data) reads 7 fields directly by
# name from a single annual, 3-hourly (8x daily) NetCDF file on the model's own tracer grid:
#   airtmp, wndewd, wndnwd, spchmd  - instantaneous, 1st record 00z Jan 1
#   glbrad, dlwsfc, ttlpcp          - 3-hour averages, 1st record 00z-03z Jan 1
# This script builds that file directly from the raw annual JRA55-do source files (no
# intermediate RYF step is needed for a single real calendar year - the source files already
# have the correct record count/order for the requested year).
#
# Precipitation note: CICE's JRA55(do) reader has no rain/snow partitioning - whatever is read
# into `ttlpcp` is added straight to `fsnow` (see ice_forcing.F90:JRA55_data). Using JRA55-do's
# `prsn` (snowfall) alone for `ttlpcp`, rather than `prsn + prra` (total precipitation), is
# therefore the behaviour consistent with how this field is actually used.
#
# The run command and git commit of the script (if committed) is added to the metadata of the
# generated file, to help identify how a given forcing file was produced.
#
# See also: README.md in this directory, and the User/Developer Guide sections on forcing
# (docs/source/developer_guide/dg_forcing.rst) for the file naming convention CICE expects.
#
# Contact:
#   Anton Steketee <anton.steketee@anu.edu.au>
#
# Dependencies: see environment.yml in this directory (argparse, numpy, xarray, xesmf)
# =========================================================================================

import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import xarray as xr
import xesmf as xe

DEFAULT_JRA55DO_DIR = (
    "/g/data/qv56/replicas/input4MIPs/CMIP6Plus/OMIP/MRI/MRI-JRA55-do-1-6-0"
)

# variable_id -> (atmos frequency subfolder, CICE field name, instantaneous?)
VARIABLE_MAP = {
    "tas": ("3hrPt", "airtmp", True),
    "uas": ("3hrPt", "wndewd", True),
    "vas": ("3hrPt", "wndnwd", True),
    "huss": ("3hrPt", "spchmd", True),
    "rsds": ("3hr", "glbrad", False),
    "rlds": ("3hr", "dlwsfc", False),
    "prsn": ("3hr", "ttlpcp", False),
}


def find_source_file(jra55do_dir, variable_id, frequency, year):
    pattern = os.path.join(
        jra55do_dir,
        "atmos",
        frequency,
        variable_id,
        "gr",
        "v*",
        f"{variable_id}_input4MIPs_atmosphericState_OMIP_*_gr_{year}*.nc",
    )
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one source file for {variable_id} {year}, found "
            f"{len(matches)}: {matches} (pattern: {pattern})"
        )
    return matches[0]


def tracer_cell_centres(hgrid_filename):
    """
    Get T-cell centre lon/lat from a MOM supergrid file.
    """
    hgrid = xr.open_dataset(hgrid_filename)
    lon = hgrid.x[1::2, 1::2]
    lat = hgrid.y[1::2, 1::2]
    return lon.rename({"nyp": "ny", "nxp": "nx"}), lat.rename(
        {"nyp": "ny", "nxp": "nx"}
    )


def build_regridder(source_ds, dest_lon, dest_lat):
    grid_src = {"lon": source_ds["lon"], "lat": source_ds["lat"]}
    grid_dest = {"lon": dest_lon, "lat": dest_lat}
    return xe.Regridder(
        grid_src,
        grid_dest,
        method="bilinear",
        extrap_method="nearest_s2d",
        periodic=True,
    )


def regrid_variable(
    variable_id, frequency, jra55do_dir, year, regridder, dest_lon, dest_lat
):
    source_filename = find_source_file(jra55do_dir, variable_id, frequency, year)
    source_ds = xr.open_dataset(source_filename)

    regridded = regridder(source_ds[variable_id], keep_attrs=True)
    regridded = regridded.assign_coords(lon=dest_lon, lat=dest_lat)

    return source_filename, source_ds["time"], regridded


def get_provenance(input_files, runcmd):
    """
    Return a dict of global attributes documenting how this file was generated: the
    command used, the script's git commit (if this file is part of a git checkout), and
    the list of input files used.
    """
    script_path = os.path.abspath(__file__)
    history = f"Generated on {os.uname().nodename} using: {runcmd}"

    try:
        commit = (
            subprocess.check_output(
                ["git", "-C", os.path.dirname(script_path), "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        history += f" (script at git commit {commit})"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    return {
        "history": history,
        "inputFile": ", ".join(os.path.abspath(f) for f in input_files),
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Regrid one year of JRA55-do forcing onto a CICE standalone grid, in the "
            "format expected by CICE's ice_forcing.F90 JRA55(do) option."
        )
    )
    parser.add_argument(
        "--year", type=int, required=True, help="Forcing year to generate, e.g. 2005."
    )
    parser.add_argument(
        "--hgrid-filename",
        type=str,
        required=True,
        help="Path to the MOM supergrid file for the destination CICE grid.",
    )
    parser.add_argument(
        "--output-filename",
        type=str,
        required=True,
        help="Path to the file to be written.",
    )
    parser.add_argument(
        "--jra55do-dir",
        type=str,
        default=DEFAULT_JRA55DO_DIR,
        help=f"Root directory of the JRA55-do v1.6.0 (input4MIPs) source data. Default: {DEFAULT_JRA55DO_DIR}",
    )
    args = parser.parse_args()

    hgrid_filename = os.path.abspath(args.hgrid_filename)
    output_filename = os.path.abspath(args.output_filename)

    dest_lon, dest_lat = tracer_cell_centres(hgrid_filename)

    # All JRA55-do variables share the same source grid (input4MIPs "gr" TL319 gaussian
    # grid), regardless of frequency, so one regridder can be reused for every variable.
    first_variable_id, (first_frequency, _, _) = next(iter(VARIABLE_MAP.items()))
    first_source_filename = find_source_file(
        args.jra55do_dir, first_variable_id, first_frequency, args.year
    )
    regridder = build_regridder(
        xr.open_dataset(first_source_filename), dest_lon, dest_lat
    )

    data_vars = {}
    input_files = [hgrid_filename]
    reference_time = None

    for variable_id, (frequency, cice_name, _instantaneous) in VARIABLE_MAP.items():
        print(f"Regridding {variable_id} -> {cice_name} ({frequency}) ...", flush=True)
        source_filename, time, regridded = regrid_variable(
            variable_id,
            frequency,
            args.jra55do_dir,
            args.year,
            regridder,
            dest_lon,
            dest_lat,
        )
        input_files.append(source_filename)

        if reference_time is None:
            reference_time = time
        elif reference_time.sizes["time"] != time.sizes["time"]:
            raise ValueError(
                f"{variable_id} has {time.sizes['time']} timesteps, expected "
                f"{reference_time.sizes['time']} (from an earlier variable)"
            )

        regridded = regridded.rename(cice_name).astype("float32")
        regridded = regridded.reset_coords(drop=True)
        # Instantaneous and 3hr-average variables carry different literal time
        # coordinate values (e.g. 00:00 vs 01:30). Drop them so xr.Dataset(...) combines
        # the 7 variables positionally instead of outer-joining on disjoint timestamps.
        regridded = regridded.drop_vars("time")
        data_vars[cice_name] = regridded

    forcing = xr.Dataset(data_vars)
    forcing = forcing.assign_coords(
        time=reference_time.values,
        ny=("ny", np.arange(forcing.sizes["ny"], dtype="float32")),
        nx=("nx", np.arange(forcing.sizes["nx"], dtype="float32")),
        lon=(("ny", "nx"), dest_lon.values),
        lat=(("ny", "nx"), dest_lat.values),
    )
    forcing["lon"].attrs = dict(
        long_name="Longitude of T-cell center",
        standard_name="longitude",
        units="degrees_east",
    )
    forcing["lat"].attrs = dict(
        long_name="Latitude of T-cell center",
        standard_name="latitude",
        units="degree_north",
    )
    forcing["ny"].attrs = dict(axis="Y", cartesian_axis="Y")
    forcing["nx"].attrs = dict(axis="X", cartesian_axis="X")

    forcing.attrs = {
        "title": f"CICE standalone JRA55-do forcing for {args.year}",
        "source": "MRI JRA55-do 1.6.0, regridded for use with CICE's JRA55do standalone forcing option",
        "comment": "ttlpcp is JRA55-do snowfall (prsn) only.",
        "regrid_method": "bilinear",
    }
    forcing.attrs.update(
        get_provenance(
            input_files=input_files,
            runcmd=f"{sys.executable} {' '.join(sys.argv)}",
        )
    )

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    var_encoding = dict(zlib=True, complevel=4, _FillValue=None)
    encoding = {var: var_encoding for var in forcing.data_vars}
    encoding["time"] = {"dtype": "double"}

    print(f"Writing {output_filename} ...", flush=True)
    forcing.to_netcdf(output_filename, unlimited_dims=["time"], encoding=encoding)
    print("Done", flush=True)


if __name__ == "__main__":
    main()
