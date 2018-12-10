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
    ofname = os.path.join(output_dir, "wtc_two_leaf.csv")
    fpath = "/Users/mdekauwe/Downloads/"
    fname = "met_data.csv"
    fn = os.path.join(fpath, fname)
    df = pd.read_csv(fn)
    df = df.drop(df.columns[0], axis=1)
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

    # Not sure which treatments Dushan wants to run, so will just use this one
    # Easy to edit as I'm passing to a func
    dfx = df[(df.T_treatment == "ambient") &
             (df.Water_treatment == "control") &
             (df.chamber == "C01")]

    (out) = run_treatment(T, dfx, footprint)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if os.path.isfile(ofname):
        os.remove(self.out_fname)
    out.to_csv(ofname, index=False)

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(16,4))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot(out.An_can, label="Model")
    ax1.plot(out.An_obs, label="Observations")
    ax1.set_ylabel("GPP (g C m$^{-2}$ d$^{-1}$)")
    ax1.set_xlabel("Days", position=(1.1, 0.5))
    ax1.legend(numpoints=1, loc="best")

    ax2.plot(out.E_can, label="Model")
    ax2.plot(out.E_obs, label="Observations")
    ax2.set_ylabel("E (mm d$^{-1}$)")

    ax1.locator_params(nbins=6, axis="y")
    ax2.locator_params(nbins=6, axis="y")

    plt.show()
