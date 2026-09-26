#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare the big-leaf and two-leaf models over an idealised clear-sky summer
day at the WTC site, for a sparse and a dense canopy.

"""
import os
import argparse
import numpy as np

import matplotlib.pyplot as plt

import constants as c
import parameters as p
from big_leaf import Canopy as BigLeaf
from two_leaf import Canopy as TwoLeaf
from radiation import calculate_cos_zenith

__author__  = "Martin De Kauwe"
__version__ = "1.0 (26.09.2026)"
__email__   = "mdekauwe@gmail.com"

TWO_LEAF_COLOUR = "#2a78d6"
BIG_LEAF_COLOUR = "#eb6834"
INK = "#0b0b0b"
INK_MUTED = "#8a8985"
GRID = "#e6e5e1"

def idealised_day(doy, tmin=18.0, tmax=32.0, vpd_min=0.8, vpd_max=3.5,
                  par_max=2000.0):
    """ Half-hourly clear-sky forcing; temperature & VPD peak mid-afternoon """
    hod = np.arange(48) / 2.0 + 0.25
    cos_zen = np.array([calculate_cos_zenith(doy, p.lat, h)[0] for h in hod])
    par = np.where(cos_zen > 0.01, par_max * cos_zen, 0.0)
    shape = 0.5 * (1.0 + np.cos(2.0 * np.pi * (hod - 15.0) / 24.0))
    tair = tmin + (tmax - tmin) * shape
    vpd = vpd_min + (vpd_max - vpd_min) * shape

    return (hod, par, tair, vpd)

def run(doy, lai, wind=2.5, pressure=101325.0, Ca=400.0):

    (hod, par, tair, vpd) = idealised_day(doy)
    T = TwoLeaf(p, gs_model="medlyn")
    B = BigLeaf(p, gs_model="medlyn")

    out = {k:np.zeros(len(hod)) for k in ["An_2l", "E_2l", "T_2l",
                                          "An_bl", "E_bl", "T_bl"]}
    for i in range(len(hod)):
        (An, et, Tcan,
         apar, lai_leaf) = T.main(tair[i], par[i], vpd[i], wind, pressure, Ca,
                                  doy, hod[i], lai)
        out["An_2l"][i] = np.sum(An)
        out["E_2l"][i] = np.sum(et)
        out["T_2l"][i] = np.sum(Tcan * lai_leaf) / np.sum(lai_leaf)

        (An, gsw, et, Tcan) = B.main(tair[i], par[i], vpd[i], wind, pressure,
                                     Ca, doy, hod[i], lai)
        out["An_bl"][i] = An
        out["E_bl"][i] = et
        out["T_bl"][i] = Tcan

    return (hod, tair, out)

def daily_total(x, conv):
    # half-hourly flux -> daily total
    return np.sum(x) * conv * 1800.0

def set_style():
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.edgecolor'] = INK_MUTED
    plt.rcParams['axes.labelcolor'] = INK
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['xtick.color'] = INK_MUTED
    plt.rcParams['ytick.color'] = INK_MUTED
    plt.rcParams['xtick.labelcolor'] = INK
    plt.rcParams['ytick.labelcolor'] = INK
    plt.rcParams['text.color'] = INK

def main(doy, lais, ofname):

    set_style()
    fig, axes = plt.subplots(len(lais), 3, figsize=(13, 3.6 * len(lais)),
                             sharex=True, gridspec_kw={'hspace':0.4,
                                                       'wspace':0.3})
    an_conv = c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C
    et_conv = c.MOL_WATER_2_G_WATER * c.G_TO_KG

    for row, lai in enumerate(lais):
        (hod, tair, o) = run(doy, lai)

        panels = [("An", 1.0, "A$_{n}$ (μmol m$^{-2}$ s$^{-1}$)",
                   "%.1f g C m$^{-2}$ d$^{-1}$", an_conv),
                  ("E", c.MOL_TO_MMOL, "E (mmol m$^{-2}$ s$^{-1}$)",
                   "%.1f mm d$^{-1}$", et_conv),
                  ("T", 1.0, "T$_{canopy}$ − T$_{air}$ (°C)", None,
                   None)]

        for col, (var, scale, ylabel, total_fmt, conv) in enumerate(panels):
            ax = axes[row, col]
            y2l = o["%s_2l" % var]
            ybl = o["%s_bl" % var]
            if var == "T":
                (y2l, ybl) = (y2l - tair, ybl - tair)
                ax.axhline(0.0, color=INK_MUTED, lw=0.8)

            ax.plot(hod, y2l * scale, color=TWO_LEAF_COLOUR, lw=2,
                    label="Two-leaf")
            ax.plot(hod, ybl * scale, color=BIG_LEAF_COLOUR, lw=2,
                    label="Big-leaf")
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", color=GRID, lw=0.8)
            ax.set_axisbelow(True)
            ax.set_xlim(0, 24)
            ax.set_xticks([0, 6, 12, 18, 24])

            if total_fmt:
                txt = "Daily total: two-leaf " + total_fmt % \
                        daily_total(o["%s_2l" % var], conv) + \
                      "\n                   big-leaf " + total_fmt % \
                        daily_total(o["%s_bl" % var], conv)
                ax.set_title(txt, loc="left", fontsize=9, color=INK)

            if col == 0:
                ax.text(-0.3, 0.5, "LAI = %.2f" % (lai),
                        transform=ax.transAxes, rotation=90, va="center",
                        fontsize=12, fontweight="bold")
            if row == len(lais) - 1:
                ax.set_xlabel("Hour of day")

    axes[0, 2].legend(loc="upper right", frameon=False, fontsize=10)
    fig.suptitle("Big-leaf vs two-leaf, idealised clear summer day at the "
                 "WTC (day %d)" % (doy), x=0.07, y=1.04, ha="left", fontsize=14)
    fig.text(0.07, 0.985, "T$_{air}$ 18–32 °C, VPD 0.8–3.5 kPa, "
             "peak PAR 2000 μmol m$^{-2}$ s$^{-1}$, "
             "C$_{a}$ 400 μmol mol$^{-1}$, wind 2.5 m s$^{-1}$",
             ha="left", fontsize=10, color="#52514e")

    fig.savefig(ofname, dpi=150, bbox_inches="tight", facecolor="white")
    print("Wrote %s" % (ofname))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--doy", type=int, default=30)
    parser.add_argument("--lai", type=float, nargs="+", default=[0.5, 3.0])
    parser.add_argument("-f", "--ofname",
                        default=os.path.join("outputs",
                                             "big_leaf_vs_two_leaf.png"))
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.ofname) or ".", exist_ok=True)
    main(args.doy, args.lai, args.ofname)
