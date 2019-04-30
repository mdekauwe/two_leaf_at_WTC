#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot the A as a function of Tcan

"""
import os
import sys
import numpy as np
import math
import pandas as pd

import matplotlib.pyplot as plt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (30.04.2019)"
__email__   = "mdekauwe@gmail.com"

def plot_all_chambers(output_dir):

    fig = plt.figure(figsize=(14,9))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    for i, chamber in enumerate(chambers):


        fn = os.path.join(output_dir, "wtc_two_leaf_%s.csv" % (chamber))
        df_2l = pd.read_csv(fn)

        fn = os.path.join(output_dir, "wtc_big_leaf_%s.csv" % (chamber))
        df_bl = pd.read_csv(fn)

        ax = fig.add_subplot(4,3,1+i)

        #ax.plot(df_2l.T_can, df_2l.An_can, ".", label="2-leaf")
        #ax.plot(df_bl.T_can, df_bl.An_can, ".", label="Big leaf")

        ax.plot(df_2l.T_sun, df_2l.An_sun, ".", label="Sun")
        ax.plot(df_2l.T_sha, df_2l.An_sha, ".", label="Shade")


        if i == 10:
            ax.set_xlabel("T$_{canopy}$ ($^{\circ}\mathrm{C}$)")

        if i == 0:
            ax.legend(numpoints=1, loc="best")

        if i < 9:
            plt.setp(ax.get_xticklabels(), visible=False)

        if i != 0 and i != 3 and i != 6 and i != 9:
            plt.setp(ax.get_yticklabels(), visible=False)

        if i == 6:
            ax.set_ylabel("An (\u03BCmol m$^{-2}$s$^{-1}$)",
                          position=(2.5, 1.2))

    plt.show()

def plot_one_chamber(output_dir, chamber):

    fn = os.path.join(output_dir, "wtc_two_leaf_%s.csv" % (chamber))
    df_2l = pd.read_csv(fn)

    fn = os.path.join(output_dir, "wtc_big_leaf_%s.csv" % (chamber))
    df_bl = pd.read_csv(fn)

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    ax = fig.add_subplot(111)
    ax.plot(df_bl.T_can, df_bl.An_can, ".", label="Big leaf")
    ax.plot(df_2l.T_can, df_2l.An_can, ".", label="2-leaf")
    ax.plot(df_2l.T_sun, df_2l.An_sun, ".", label="Sun")
    ax.plot(df_2l.T_sha, df_2l.An_sha, ".", label="Shade")



    ax.set_xlabel("T$_{canopy}$ ($^{\circ}\mathrm{C}$)")
    ax.legend(numpoints=1, loc="best")
    ax.set_ylabel("An (\u03BCmol m$^{-2}$s$^{-1}$)")

    plt.show()


if __name__ == "__main__":

    output_dir = "outputs"
    chambers = ["C01", "C02", "C03", "C04", "C05", "C06", \
                "C07", "C08", "C09", "C10", "C11", "C12"]

    #plot_all_chambers(output_dir)
    plot_one_chamber(output_dir, chambers[0])
