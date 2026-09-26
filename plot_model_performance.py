#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot model performance against the WTC flux observations: observed vs
modelled An and E (pooled across chambers), per-chamber bias and RMSE and,
where the outputs carry an hour-of-day column, the mean diurnal cycle.

"""
import os
import glob
import argparse
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import constants as c

__author__  = "Martin De Kauwe"
__version__ = "1.0 (26.09.2026)"
__email__   = "mdekauwe@gmail.com"

MODEL_COLOUR = "#2a78d6"
OBS_COLOUR = "#52514e"
INK = "#0b0b0b"
INK_MUTED = "#8a8985"
GRID = "#e6e5e1"

VARS = [("An", "A$_{n}$", "μmol m$^{-2}$ s$^{-1}$", 1.0),
        ("E", "E", "mmol m$^{-2}$ s$^{-1}$", c.MOL_TO_MMOL)]

def set_style():
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.labelsize'] = 11
    plt.rcParams['axes.edgecolor'] = INK_MUTED
    plt.rcParams['axes.labelcolor'] = INK
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['xtick.color'] = INK_MUTED
    plt.rcParams['ytick.color'] = INK_MUTED
    plt.rcParams['xtick.labelcolor'] = INK
    plt.rcParams['ytick.labelcolor'] = INK
    plt.rcParams['text.color'] = INK

def read_outputs(output_dir, model):

    frames = []
    for fn in sorted(glob.glob(os.path.join(output_dir,
                                            "wtc_%s_C??.csv" % (model)))):
        df = pd.read_csv(fn)
        if len(df) == 0:
            continue
        if model == "two_leaf" and "An_sun" not in df:
            print("Skipping %s: not two-leaf output" % (fn))
            continue
        df["chamber"] = os.path.basename(fn)[-7:-4]
        frames.append(df)

    return pd.concat(frames, ignore_index=True)

def daytime(df):
    # the model only runs when the sun is up, so compare those hours
    if "APAR_can" in df:
        return df[df.APAR_can.fillna(0.0) > 0.0]
    return df[df.An_can != 0.0]

def stats(obs, mod):
    ok = np.isfinite(obs) & np.isfinite(mod)
    (obs, mod) = (obs[ok], mod[ok])
    r = np.corrcoef(obs, mod)[0, 1]
    return {'n':len(obs), 'r2':r**2, 'bias':np.mean(mod - obs),
            'rmse':np.sqrt(np.mean((mod - obs)**2))}

def plot_obs_vs_model(ax, obs, mod, label, units):

    lo = min(np.nanmin(obs), np.nanmin(mod))
    hi = max(np.nanmax(obs), np.nanmax(mod))
    pad = 0.03 * (hi - lo)
    (lo, hi) = (lo - pad, hi + pad)
    hb = ax.hexbin(obs, mod, gridsize=45, cmap="Blues", mincnt=1, bins="log",
                   extent=(lo, hi, lo, hi), linewidths=0.2)
    ax.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=1, ls="--", zorder=3)
    ax.text(0.97*hi + 0.03*lo, 0.97*hi + 0.03*lo, "1:1", color=INK_MUTED,
            ha="right", va="top", fontsize=9)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("Observed %s (%s)" % (label, units))
    ax.set_ylabel("Modelled %s (%s)" % (label, units))

    s = stats(obs, mod)
    ax.text(0.04, 0.96,
            "n = %d\nR$^2$ = %.2f\nRMSE = %.2f\nbias = %+.2f" % \
            (s['n'], s['r2'], s['rmse'], s['bias']),
            transform=ax.transAxes, va="top", fontsize=10, color=INK)

    return hb

def plot_by_chamber(ax, day, var, label, units, scale):

    rows = []
    for chamber, g in day.groupby("chamber"):
        s = stats(g["%s_obs" % var].values * scale,
                  g["%s_can" % var].values * scale)
        rows.append((chamber, s['bias'], s['rmse']))
    (chambers, bias, rmse) = zip(*rows)
    y = np.arange(len(chambers))

    ax.axvline(0.0, color=INK_MUTED, lw=1)
    ax.hlines(y, 0.0, bias, color=GRID, lw=2, zorder=1)
    ax.scatter(bias, y, s=40, color=MODEL_COLOUR, zorder=3, label="Bias",
               edgecolor="white", linewidth=1.5)
    ax.scatter(rmse, y, s=40, marker="D", facecolor="white",
               edgecolor=OBS_COLOUR, linewidth=1.5, zorder=3, label="RMSE")
    ax.set_yticks(y)
    ax.set_yticklabels(chambers)
    ax.set_xlabel("%s error, model − obs (%s)" % (label, units))
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_ylim(len(chambers) - 0.5, -0.5)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2,
              frameon=False, fontsize=9, borderaxespad=0.2)

def plot_diurnal(ax, df, var, label, units, scale):

    g = df.groupby(df.hod.astype(int))
    obs = g["%s_obs" % var]
    mod = g["%s_can" % var]
    hrs = np.array(list(g.groups.keys()))

    ax.fill_between(hrs, obs.quantile(0.25) * scale,
                    obs.quantile(0.75) * scale, color=OBS_COLOUR, alpha=0.15,
                    lw=0)
    ax.plot(hrs, obs.mean() * scale, color=OBS_COLOUR, lw=2,
            label="Observed (mean, IQR)")
    ax.plot(hrs, mod.mean() * scale, color=MODEL_COLOUR, lw=2,
            label="Two-leaf model")
    ax.axhline(0.0, color=INK_MUTED, lw=0.8)
    ax.set_xlim(0, 23)
    ax.set_xticks([0, 6, 12, 18])
    ax.set_xlabel("Hour of day")
    ax.set_ylabel("%s (%s)" % (label, units))
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", frameon=False, fontsize=9)


def main(output_dir, model, ofname, title_note):

    set_style()
    df = read_outputs(output_dir, model)
    day = daytime(df)
    has_hod = "hod" in df

    nrows = 3 if has_hod else 2
    fig, axes = plt.subplots(nrows, 2, figsize=(10, 4.2 * nrows),
                             gridspec_kw={'hspace':0.5, 'wspace':0.35,
                                          'height_ratios':[1.3] + \
                                                          [1] * (nrows - 1)})
    fig.subplots_adjust(top=1.0 - 0.3 / nrows)

    for col, (var, label, units, scale) in enumerate(VARS):
        hb = plot_obs_vs_model(axes[0, col], day["%s_obs" % var].values * scale,
                               day["%s_can" % var].values * scale, label, units)
        cb = fig.colorbar(hb, ax=axes[0, col], shrink=0.8, pad=0.03)
        cb.set_label("Hours per bin")
        cb.outline.set_visible(False)

        plot_by_chamber(axes[1, col], day, var, label, units, scale)

        if has_hod:
            plot_diurnal(axes[2, col], df, var, label, units, scale)

    title = "%s model vs WTC observations, daytime hours" % \
                ("Two-leaf" if model == "two_leaf" else "Big-leaf")
    top = axes[0, 0].get_position().y1
    fig.text(0.06, top + 0.1 / nrows, title, ha="left", fontsize=14)
    if title_note:
        fig.text(0.06, top + 0.065 / nrows, title_note, ha="left",
                 fontsize=10, color=OBS_COLOUR)

    fig.savefig(ofname, dpi=150, bbox_inches="tight", facecolor="white")
    print("Wrote %s" % (ofname))

    for var, label, units, scale in VARS:
        print("\n%s (%s)" % (var, units.replace("$", "")))
        for chamber, g in day.groupby("chamber"):
            s = stats(g["%s_obs" % var].values * scale,
                      g["%s_can" % var].values * scale)
            print("  %s  n=%5d  R2=%.2f  RMSE=%.2f  bias=%+.2f" % \
                  (chamber, s['n'], s['r2'], s['rmse'], s['bias']))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output_dir", default="outputs")
    parser.add_argument("-m", "--model", default="two_leaf",
                        choices=["two_leaf", "big_leaf"])
    parser.add_argument("-f", "--ofname", default=None)
    parser.add_argument("-n", "--note", default="")
    args = parser.parse_args()

    ofname = args.ofname or os.path.join(args.output_dir,
                                         "performance_%s.png" % (args.model))
    main(args.output_dir, args.model, ofname, args.note)
