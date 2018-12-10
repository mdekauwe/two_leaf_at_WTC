#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot the model vs obs.

"""
import os
import sys
import numpy as np
import math
import pandas as pd

import matplotlib.pyplot as plt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (07.12.2018)"
__email__   = "mdekauwe@gmail.com"

chamber = "C01"

output_dir = "outputs"
fn = os.path.join(output_dir, "wtc_two_leaf_%s.csv" % (chamber))
df = pd.read_csv(fn)

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

ax1.plot(df.An_can, label="Model")
ax1.plot(df.An_obs, label="Observations")
ax1.set_ylabel("GPP (g C m$^{-2}$ d$^{-1}$)")
ax1.set_xlabel("Days", position=(1.1, 0.5))
ax1.legend(numpoints=1, loc="best")

ax2.plot(df.E_can, label="Model")
ax2.plot(df.E_obs, label="Observations")
ax2.set_ylabel("E (mm d$^{-1}$)")

ax1.locator_params(nbins=6, axis="y")
ax2.locator_params(nbins=6, axis="y")

plt.show()
