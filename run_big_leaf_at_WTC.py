#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Apply the big-leaf model to the WTC experiments.

"""
import os
import argparse
import numpy as np
import pandas as pd

import constants as c
import parameters as p
from big_leaf import Canopy as BigLeaf
from run_two_leaf_at_WTC import load_met

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.12.2018)"
__email__   = "mdekauwe@gmail.com"


def run_treatment(B, df, p, wind, pressure, Ca, vary_vj=False):

    rows = []
    for (dt, r) in df.iterrows():

        if vary_vj:
            (Vcmax25, Jmax25) = (r.Vcmax25, r.Jmax25)
        else:
            (Vcmax25, Jmax25) = (p.Vcmax25, p.Jmax25)

        (An, gsw, et, Tcan) = B.main(r.tair, r.par, r.vpd, wind, pressure, Ca,
                                     r.doy, r.hod, r.lai, Vcmax25=Vcmax25,
                                     Jmax25=Jmax25)

        rows.append({'year':dt.year, 'doy':r.doy, 'hod':r.hod,
                     # Convert from per tree to m-2
                     'An_obs':r.FluxCO2 * c.MMOL_2_UMOL / p.footprint,
                     'E_obs':r.FluxH2O / p.footprint,
                     'An_can':An, 'E_can':et, 'T_can':Tcan})

    return pd.DataFrame(rows, index=df.index)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("met_fname", nargs="?",
                        default="/Users/mdekauwe/Downloads/"
                                "met_data_gap_fixed_V1.csv")
    parser.add_argument("-o", "--output_dir", default="outputs")
    args = parser.parse_args()

    df = load_met(args.met_fname)

    ##  Fixed met stuff
    #
    wind = 2.5
    pressure = 101325.0
    Ca = 400.0

    B = BigLeaf(p, gs_model="medlyn")

    os.makedirs(args.output_dir, exist_ok=True)

    chambers = np.unique(df.chamber)
    for chamber in chambers:
        print(chamber)
        dfx = df[(df.T_treatment == "ambient") &
                 (df.Water_treatment == "control") &
                 (df.chamber == chamber)].copy()

        out = run_treatment(B, dfx, p, wind, pressure, Ca, vary_vj=False)

        ofname = os.path.join(args.output_dir, "wtc_big_leaf_%s.csv" % (chamber))
        out.to_csv(ofname, index=False)
