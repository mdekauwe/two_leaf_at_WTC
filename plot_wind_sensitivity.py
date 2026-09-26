#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sensitivity of the two-leaf model to the (unmeasured) chamber air speed and
within-canopy wind extinction coefficient, across the range of WTC LAI, for
the idealised clear-sky summer day used in compare_big_leaf_two_leaf.py.

"""
import os
import argparse
import numpy as np

import matplotlib.pyplot as plt

import constants as c
import parameters as p
from compare_big_leaf_two_leaf import run, set_style, INK_MUTED, GRID

__author__  = "Martin De Kauwe"
__version__ = "1.0 (26.09.2026)"
__email__   = "mdekauwe@gmail.com"

# light -> dark = low -> high
SHADES = ["#8fbcef", "#2a78d6", "#123f7a"]

def summarise(doy, lai, wind, a):

    p.wind_extinction = a
    (hod, tair, o) = run(doy, lai, wind=wind)
    midday = (hod >= 11.0) & (hod <= 14.0)
    dT = np.mean(o["T_2l"][midday] - tair[midday])
    an = np.sum(o["An_2l"]) * c.UMOL_TO_MOL * c.MOL_C_TO_GRAMS_C * 1800.0

    return (dT, an)

def main(doy, ofname):

    set_style()
    lais = np.linspace(0.25, 3.5, 14)
    winds = [0.5, 1.0, 2.5]
    extinctions = [0.0, 1.0, 2.0]
    default_a = p.wind_extinction

    cases = [("Chamber air speed (a = 0)",
              [("%.1f m s$^{-1}$" % w, w, 0.0) for w in winds]),
             ("Wind extinction a (air speed 2.5 m s$^{-1}$)",
              [("a = %.0f" % a, 2.5, a) for a in extinctions])]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True,
                             gridspec_kw={'hspace':0.2, 'wspace':0.25})

    for col, (title, series) in enumerate(cases):
        for (label, wind, a), colour in zip(series, SHADES):
            res = np.array([summarise(doy, lai, wind, a) for lai in lais])
            for row in range(2):
                ax = axes[row, col]
                ax.plot(lais, res[:, row], color=colour, lw=2)
                ax.plot(lais[-1], res[-1, row], "o", color=colour, ms=5)
                # label the top row only, lines keep the same colours below
                if row == 0:
                    ax.annotate(label, (lais[-1], res[-1, row]),
                                xytext=(6, 0), textcoords="offset points",
                                va="center", fontsize=9)
        axes[0, col].set_title(title, loc="left", fontsize=11)

    p.wind_extinction = default_a

    for ax in axes.flat:
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 4.4)
    for ax in axes[1]:
        ax.set_xlabel("LAI (m$^{2}$ m$^{-2}$)")
    for ax in axes[0]:
        ax.set_ylim(bottom=0)
    axes[0, 0].set_ylabel("Midday T$_{canopy}$ − T$_{air}$ (°C)")
    axes[1, 0].set_ylabel("Daily A$_{n}$ (g C m$^{-2}$ d$^{-1}$)")

    fig.suptitle("Two-leaf sensitivity to chamber air flow, idealised clear "
                 "summer day (day %d)" % (doy), x=0.07, y=0.99, ha="left",
                 fontsize=14)
    fig.text(0.07, 0.945, "Midday = 11:00–14:00 mean. T$_{air}$ "
             "18–32 °C, VPD 0.8–3.5 kPa, peak PAR 2000 "
             "μmol m$^{-2}$ s$^{-1}$", ha="left", fontsize=10,
             color="#52514e")

    fig.savefig(ofname, dpi=150, bbox_inches="tight", facecolor="white")
    print("Wrote %s" % (ofname))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--doy", type=int, default=30)
    parser.add_argument("-f", "--ofname",
                        default=os.path.join("outputs", "wind_sensitivity.png"))
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.ofname) or ".", exist_ok=True)
    main(args.doy, args.ofname)
