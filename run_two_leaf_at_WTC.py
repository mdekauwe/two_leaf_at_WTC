#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Apply the two-leaf model to the WTC experiments.

"""
import os
import sys
import numpy as np
import math
import pandas as pd

import constants as c
from radiation import calculate_absorbed_radiation
from two_leaf import CoupledModel as TwoLeaf

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.12.2018)"
__email__   = "mdekauwe@gmail.com"


def run_treatment(T, df, footprint):

    wind = 5.0
    pressure = 101325.0
    Ca = 400.0

    days = df.doy
    hod = df.hod
    ndays = int(len(days) / 24.)

    out = setup_output_dataframe(ndays)

    i = 0
    j = 0
    while i < len(df):
        year = df.index.year[i]
        doy = df.doy[i]
        hod = 0
        for k in range(24):

            (An, et, Tcan,
             apar, lai_leaf) = T.main(df.tair[i], df.par[i], df.vpd[i],
                                       wind, pressure, Ca, doy, hod, lat,
                                       lon, df.LAI[i])

            out = update_output_hourly(i, j, An, et, Tcan, apar, lai_leaf, df,
                                       footprint, out)

            hod += 1
            i += 1

        out = update_output_daily(j, year, doy, out)
        j += 1

    return (out)

def setup_output_dataframe(ndays):

    zero = np.zeros(ndays)
    out = pd.DataFrame({'year':zero, 'doy':zero,
                        'An_obs':zero, 'E_obs':zero,
                        'An_can':zero, 'An_sun':zero, 'An_sha':zero,
                        'E_can':zero, 'E_sun':zero, 'E_sha':zero,
                        'T_can':zero, 'T_sun':zero, 'T_sha':zero,
                        'APAR_can':zero, 'APAR_sun':zero, 'APAR_sha':zero,
                        'LAI_can':zero, 'LAI_sun':zero, 'LAI_sha':zero})
    return (out)

def update_output_hourly(i, j, An, et, Tcan, apar, lai_leaf, df, footprint, out):

    an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * c.SEC_TO_HR
    et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG * c.SEC_TO_HR
    sun_frac = lai_leaf[c.SUNLIT] / np.sum(lai_leaf)
    sha_frac = lai_leaf[c.SHADED] / np.sum(lai_leaf)
    out.An_can[j] += np.sum(An) * an_conv
    out.An_sun[j] += An[c.SUNLIT] * an_conv
    out.An_sha[j] += An[c.SHADED] * an_conv
    out.E_can[j] += np.sum(et) * et_conv
    out.E_sun[j] += et[c.SUNLIT] * et_conv
    out.E_sha[j] += et[c.SHADED] * et_conv
    out.T_can[j] += (Tcan[c.SUNLIT] * sun_frac) + (Tcan[c.SHADED] * sha_frac)
    out.T_sun[j] += Tcan[c.SUNLIT]
    out.T_sha[j] += Tcan[c.SHADED]
    out.APAR_can[j] += np.sum(apar)
    out.APAR_sun[j] += apar[c.SUNLIT]
    out.APAR_sha[j] += apar[c.SHADED]
    out.LAI_can[j] += np.sum(lai_leaf)
    out.LAI_sun[j] += lai_leaf[c.SUNLIT]
    out.LAI_sha[j] += lai_leaf[c.SHADED]

    # Convert from per tree to m-2
    out.An_obs[j] += df.FluxCO2[i] * c.MMOL_2_UMOL * an_conv / footprint
    out.E_obs[j] += df.FluxH2O[i] * et_conv / footprint

    return out

def update_output_daily(j, year, doy, out):

    out.year[j] = year
    out.doy[j] = doy
    out.T_can[j] /= 24.
    out.T_sun[j] /= 24.
    out.T_sha[j] /= 24.
    out.LAI_can[j] /= 24.
    out.LAI_sun[j] /= 24.
    out.LAI_sha[j] /= 24.

    return out

if __name__ == "__main__":

    output_dir = "outputs"
    fpath = "/Users/mdekauwe/Downloads/"
    fname = "met_data_gap_fixed.csv"
    fn = os.path.join(fpath, fname)
    df = pd.read_csv(fn)
    #df = df.drop(df.columns[0], axis=1)
    df.index = pd.to_datetime(df.DateTime)

    #
    ## Parameters - Dushan to set these ...
    #
    lat = -33.617778 # Ellsworth 2017, NCC
    lon = 150.740278
    g0 = 1E-09
    g1 = 3.8
    D0 = 1.5 # kpa # Not used so ignore ...
    Vcmax25 = 81.706
    Jmax25 = Vcmax25 * 1.67
    Rd25 = 2.0  # Need to discuss what you want here, "None" -> Vcmax = 0.015 Rd
    Eaj = 30000.0
    Eav = 60000.0
    deltaSj = 650.0
    deltaSv = 650.0
    Hdv = 200000.0
    Hdj = 200000.0
    Q10 = 2.0
    gamma = 0.0
    leaf_width = 0.02
    SW_abs = 0.8 # use canopy absorptance of solar radiation, not used anyway...

    diameter = 3.25 # chamber
    footprint = np.pi * (diameter / 2.)**2 # to convert from tree to m2

    # Add an LAI field, i.e. converting from per tree to m2 m-2
    df = df.assign(LAI = lambda x: x.leafArea/footprint)

    T = TwoLeaf(g0, g1, D0, gamma, Vcmax25, Jmax25, Rd25, Eaj, Eav, deltaSj,
                deltaSv, Hdv, Hdj, Q10, leaf_width, SW_abs, gs_model="medlyn")

    chambers = np.unique(df.chamber)
    chambers = ["C01"]
    for chamber in chambers:
        print(chamber)
        dfx = df[(df.T_treatment == "ambient") &
                 (df.Water_treatment == "control") &
                 (df.chamber == chamber)].copy()

        (out) = run_treatment(T, dfx, footprint)

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        ofname = os.path.join(output_dir, "wtc_two_leaf_%s.csv" % (chamber))
        if os.path.isfile(ofname):
            os.remove(ofname)
        out.to_csv(ofname, index=False)

    
