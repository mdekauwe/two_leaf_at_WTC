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
import parameters as p
from radiation import calculate_absorbed_radiation
from big_leaf import Canopy as BigLeaf

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.12.2018)"
__email__   = "mdekauwe@gmail.com"


def run_treatment(B, df, p, wind, pressure, Ca, vary_vj=False):

    days = df.doy
    hod = df.hod
    ndays = int(len(days) / 24.)
    nhours = len(df)

    out = setup_output_dataframe(nhours)

    i = 0
    j = 0
    while i < len(df):
        year = df.index.year[i]
        doy = df.doy[i]
        hod = df.hod[i]

        if vary_vj:
            (An, gsw, et, Tcan) = B.main(df.tair[i], df.par[i], df.vpd[i], wind,
                                         pressure, Ca, doy, hod, df.lai[i],
                                         Vcmax25=df.Vcmax25[i],
                                         Jmax25=df.Jmax25[i])
        else:
            (An, gsw, et, Tcan) = B.main(df.tair[i], df.par[i], df.vpd[i], wind,
                                         pressure, Ca, doy, hod, df.lai[i],
                                         Vcmax25=p.Vcmax25, Jmax25=p.Jmax25)

        out = update_output_hourly(doy, i, An, et, Tcan, df, p.footprint, out)

        i += 1

    return (out)

def setup_output_dataframe(ndays):

    zero = np.zeros(ndays)
    out = pd.DataFrame({'year':zero, 'doy':zero,
                        'An_obs':zero, 'E_obs':zero,
                        'An_can':zero, 'E_can':zero,
                        'T_can':zero})
    return (out)

def update_output_hourly(doy, j, An, et, Tcan, df, footprint, out):

    out.An_can[j] = np.sum(An)
    out.E_can[j] = np.sum(et)
    out.T_can[j] = Tcan

    # Convert from per tree to m-2
    out.An_obs[j] = df.FluxCO2[j] * c.MMOL_2_UMOL / footprint
    out.E_obs[j] = df.FluxH2O[j] / footprint

    return out

if __name__ == "__main__":

    output_dir = "outputs"
    fpath = "/Users/mdekauwe/Downloads/"
    fname = "met_data_gap_fixed_V1.csv"
    fn = os.path.join(fpath, fname)
    df = pd.read_csv(fn)
    #df = df.drop(df.columns[0], axis=1)
    df.index = pd.to_datetime(df.DateTime)

    # Add an LAI field, i.e. converting from per tree to m2 m-2
    df = df.assign(lai = lambda x: x.leafArea / p.footprint)

    ##  Fixed met stuff
    #
    wind = 2.5
    pressure = 101325.0
    Ca = 400.0

    B = BigLeaf(p, gs_model="medlyn")

    chambers = np.unique(df.chamber)
    for chamber in chambers:
        print(chamber)
        dfx = df[(df.T_treatment == "ambient") &
                 (df.Water_treatment == "control") &
                 (df.chamber == chamber)].copy()

        (out) = run_treatment(B, dfx, p, wind, pressure, Ca, vary_vj=False)

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        ofname = os.path.join(output_dir, "wtc_big_leaf_%s.csv" % (chamber))
        if os.path.isfile(ofname):
            os.remove(ofname)
        out.to_csv(ofname, index=False)
