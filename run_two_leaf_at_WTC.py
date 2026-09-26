#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Apply the two-leaf model to the WTC experiments, writing hourly and daily
outputs for each chamber.

"""
import os
import argparse
import numpy as np
import pandas as pd

import constants as c
import parameters as p
from two_leaf import Canopy as TwoLeaf

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.12.2018)"
__email__   = "mdekauwe@gmail.com"


def load_met(fn):
    df = pd.read_csv(fn)
    df.index = pd.to_datetime(df.DateTime)

    # Add an LAI field, i.e. converting from per tree to m2 m-2
    df = df.assign(lai=lambda x: x.leafArea / p.footprint)

    return df

def run_treatment(T, df, p, wind, pressure, Ca, vary_vj=False):

    rows = []
    for (dt, r) in df.iterrows():

        if vary_vj:
            (Vcmax25, Jmax25) = (r.Vcmax25, r.Jmax25)
        else:
            (Vcmax25, Jmax25) = (p.Vcmax25, p.Jmax25)

        (An, et, Tcan,
         apar, lai_leaf) = T.main(r.tair, r.par, r.vpd, wind, pressure, Ca,
                                  r.doy, r.hod, r.lai, Vcmax25=Vcmax25,
                                  Jmax25=Jmax25)

        rows.append(hourly_output(dt, r, An, et, Tcan, apar, lai_leaf,
                                  p.footprint))

    out = pd.DataFrame(rows, index=df.index)
    out_day = daily_output(out)

    return (out, out_day)

def hourly_output(dt, r, An, et, Tcan, apar, lai_leaf, footprint):

    lai_can = np.sum(lai_leaf)
    if lai_can > 0.0:
        T_can = np.sum(Tcan * lai_leaf) / lai_can
    else:
        T_can = r.tair

    return {'year':dt.year, 'doy':r.doy, 'hod':r.hod,
            # Convert from per tree to m-2
            'An_obs':r.FluxCO2 * c.MMOL_2_UMOL / footprint,
            'E_obs':r.FluxH2O / footprint,
            'An_can':np.sum(An), 'An_sun':An[c.SUNLIT], 'An_sha':An[c.SHADED],
            'E_can':np.sum(et), 'E_sun':et[c.SUNLIT], 'E_sha':et[c.SHADED],
            'T_can':T_can, 'T_sun':Tcan[c.SUNLIT], 'T_sha':Tcan[c.SHADED],
            'APAR_can':np.sum(apar), 'APAR_sun':apar[c.SUNLIT],
            'APAR_sha':apar[c.SHADED],
            'LAI_can':lai_can, 'LAI_sun':lai_leaf[c.SUNLIT],
            'LAI_sha':lai_leaf[c.SHADED]}

def daily_output(out):
    """
    Daily totals of An (g C m-2 d-1) and E (mm d-1) and mean temperatures,
    keeping only complete days.
    """
    tstep = out.index.to_series().diff().median().total_seconds()
    steps_per_day = int(round(c.SEC_TO_DAY / tstep))

    an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * tstep
    et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG * tstep

    an_cols = ['An_can', 'An_sun', 'An_sha']
    et_cols = ['E_can', 'E_sun', 'E_sha']
    t_cols = ['T_can', 'T_sun', 'T_sha']

    grp = out.groupby(out.index.normalize())
    out_day = pd.concat([grp[['year', 'doy']].first(),
                         grp[an_cols].sum() * an_conv,
                         grp[et_cols].sum() * et_conv,
                         grp[t_cols].mean()], axis=1)

    return out_day[grp.size() == steps_per_day]


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

    T = TwoLeaf(p, gs_model="medlyn")

    os.makedirs(args.output_dir, exist_ok=True)

    chambers = np.unique(df.chamber)
    for chamber in chambers:
        print(chamber)
        dfx = df[(df.T_treatment == "ambient") &
                 (df.Water_treatment == "control") &
                 (df.chamber == chamber)].copy()

        (out,
         out_day) = run_treatment(T, dfx, p, wind, pressure, Ca, vary_vj=False)

        ofname = os.path.join(args.output_dir, "wtc_two_leaf_%s.csv" % (chamber))
        out.to_csv(ofname, index=False)

        ofname = os.path.join(args.output_dir,
                              "wtc_two_leaf_day_%s.csv" % (chamber))
        out_day.to_csv(ofname, index=False)
