#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Solve 30-minute coupled A-gs(E) using a two-leaf approximation roughly following
Wang and Leuning.

References:
----------
* Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.
* Dai et al. (2004) Journal of Climate, 17, 2281-2299.
* De Pury & Farquhar (1997) PCE, 20, 537-557.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, sin, exp, sqrt, acos, asin
import random
import math

import constants as c
from farq import FarquharC3
from penman_monteith_leaf import PenmanMonteith
from radiation import spitters
from radiation import calculate_absorbed_radiation
from radiation import calculate_cos_zenith, calc_leaf_to_canopy_scalar

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.11.2018)"
__email__   = "mdekauwe@gmail.com"


class Canopy(object):
    """Iteratively solve leaf temp, Ci, gs and An."""

    def __init__(self, p, peaked_Jmax=True, peaked_Vcmax=True, model_Q10=True,
                 gs_model=None, iter_max=100):

        self.p = p

        self.peaked_Jmax = peaked_Jmax
        self.peaked_Vcmax = peaked_Vcmax
        self.model_Q10 = model_Q10
        self.gs_model = gs_model
        self.iter_max = iter_max


    def main(self, tair, par, vpd, wind, pressure, Ca, doy, hod,
             lai, rnet=None, Vcmax25=None, Jmax25=None):
        """
        Parameters:
        ----------
        tair : float
            air temperature (deg C)
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        vpd : float
            Vapour pressure deficit (kPa, needs to be in Pa, see conversion
            below)
        wind : float
            wind speed (m s-1)
        pressure : float
            air pressure (using constant) (Pa)
        Ca : float
            ambient CO2 concentration
        doy : float
            day of day
        hod : float
            hour of day
        lat : float
            latitude
        lon : float
            longitude
        lai : floar
            leaf area index

        Returns:
        --------
        An : float
            net leaf assimilation (umol m-2 s-1)
        gs : float
            stomatal conductance (mol m-2 s-1)
        et : float
            transpiration (mol H2O m-2 s-1)
        """

        F = FarquharC3(peaked_Jmax=self.peaked_Jmax,
                       peaked_Vcmax=self.peaked_Vcmax,
                       model_Q10=self.model_Q10, gs_model=self.gs_model)
        PM = PenmanMonteith()

        An = np.zeros(2) # sunlit, shaded
        gsc = np.zeros(2)  # sunlit, shaded
        et = np.zeros(2) # sunlit, shaded
        Tcan = np.zeros(2) # sunlit, shaded
        lai_leaf = np.zeros(2)
        sw_rad = np.zeros(2) # VIS, NIR
        tcanopy = np.zeros(2)

        (cos_zenith, elevation) = calculate_cos_zenith(doy, self.p.lat, hod)

        sw_rad[c.VIS] = 0.5 * (par * c.PAR_2_SW) # W m-2
        sw_rad[c.NIR] = 0.5 * (par * c.PAR_2_SW) # W m-2

        # get diffuse/beam frac, just use VIS as the answer is the same for NIR
        (diffuse_frac, direct_frac) = spitters(doy, sw_rad[0], cos_zenith)

        (qcan, apar,
         lai_leaf, kb) = calculate_absorbed_radiation(self.p, par, cos_zenith,
                                                      lai, direct_frac,
                                                      diffuse_frac, doy, sw_rad,
                                                      tair)

        # Calculate scaling term to go from a single leaf to canopy,
        # see Wang & Leuning 1998 appendix C
        scalex = calc_leaf_to_canopy_scalar(lai, kb=kb, kn=self.p.kn)

        if lai_leaf[0] < 1.e-3: # to match line 336 of CABLE radiation
            scalex[0] = 0.

        # Is the sun up?
        if elevation > 0.0 and par > 50.:

            # sunlit / shaded loop
            for ileaf in range(2):

                # initialise values of Tleaf, Cs, dleaf at the leaf surface
                dleaf = vpd
                Cs = Ca
                Tleaf = tair
                Tleaf_K = Tleaf + c.DEG_2_KELVIN


                iter = 0
                while True:

                    if scalex[ileaf] > 0.:
                        (An[ileaf],
                         gsc[ileaf]) = F.photosynthesis(self.p, Cs=Cs,
                                                        Tleaf=Tleaf_K,
                                                        Par=apar[ileaf],
                                                        vpd=dleaf,
                                                        scalex=scalex[ileaf],
                                                        Vcmax25=Vcmax25,
                                                        Jmax25=Jmax25)
                    else:
                        An[ileaf], gsc[ileaf] = 0., 0.

                    # Calculate new Tleaf, dleaf, Cs
                    (new_tleaf, et[ileaf],
                     le_et, gbH, gw) = self.calc_leaf_temp(self.p, PM, Tleaf,
                                                           tair, gsc[ileaf],
                                                           None, vpd,
                                                           pressure, wind,
                                                           rnet=qcan[ileaf],
                                                           lai=lai_leaf[ileaf])

                    gbc = gbH * c.GBH_2_GBC
                    if gbc > 0.0 and An[ileaf] > 0.0:
                        Cs = Ca - An[ileaf] / gbc # boundary layer of leaf
                    else:
                        Cs = Ca

                    if np.isclose(et[ileaf], 0.0) or np.isclose(gw, 0.0):
                        dleaf = vpd
                    else:
                        dleaf = (et[ileaf] * pressure / gw) * c.PA_2_KPA # kPa

                    # Check for convergence...?
                    if math.fabs(Tleaf - new_tleaf) < 0.02:
                        Tcan[ileaf] = Tleaf
                        break

                    if iter > self.iter_max:
                        #raise Exception('No convergence: %d' % (iter))
                        An[ileaf] = 0.0
                        gsc[ileaf] = 0.0
                        et[ileaf] = 0.0
                        break

                    # Update temperature & do another iteration
                    Tleaf = new_tleaf
                    Tleaf_K = Tleaf + c.DEG_2_KELVIN
                    Tcan[ileaf] = Tleaf

                    iter += 1

        return (An, et, Tcan, apar, lai_leaf)

    def calc_leaf_temp(self, p, PM=None, tleaf=None, tair=None, gsc=None,
                       par=None, vpd=None, pressure=None, wind=None, rnet=None,
                       lai=None):
        """
        Resolve leaf temp

        Parameters:
        ----------
        P : object
            Penman-Montheith class instance
        tleaf : float
            leaf temperature (deg C)
        tair : float
            air temperature (deg C)
        gs : float
            stomatal conductance (mol m-2 s-1)
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        vpd : float
            Vapour pressure deficit (kPa, needs to be in Pa, see conversion
            below)
        pressure : float
            air pressure (using constant) (Pa)
        wind : float
            wind speed (m s-1)

        Returns:
        --------
        new_Tleaf : float
            new leaf temperature (deg C)
        et : float
            transpiration (mol H2O m-2 s-1)
        gbH : float
            total boundary layer conductance to heat for one side of the leaf
        gw : float
            total leaf conductance to water vapour (mol m-2 s-1)
        """
        tleaf_k = tleaf + c.DEG_2_KELVIN
        tair_k = tair + c.DEG_2_KELVIN

        air_density = pressure / (c.RSPECIFC_DRY_AIR * tair_k)

        # convert from mm s-1 to mol m-2 s-1
        cmolar = pressure / (c.RGAS * tair_k)

        (grn, gh, gbH, gw) = PM.calc_conductances(p, tair_k, tleaf, tair,
                                                  wind, gsc, cmolar, lai)

        if np.isclose(gsc, 0.0):
            et = 0.0
            le_et = 0.0
        else:
            (et, le_et) = PM.calc_et(tleaf, tair, vpd, pressure, wind, par,
                                     gh, gw, rnet)

        # D6 in Leuning. NB I'm doubling conductances, see note below E5.
        # Leuning isn't explicit about grn but I think this is right
        # NB the units or grn and gbH are mol m-2 s-1 and not m s-1, but it
        # cancels.
        Y = 1.0 / (1.0 + (2.0 * grn) / (2.0 * gbH))

        # sensible heat exchanged between leaf and surroundings
        H = Y * (rnet - le_et)

        # leaf-air temperature difference recalculated from energy balance.
        # NB. I'm using gh here to include grn and the doubling of conductances
        new_Tleaf = tair + H / (c.CP * air_density * (gh / cmolar))

        # Update leaf temperature:
        new_tleaf_k = tleaf_k + (new_Tleaf + c.DEG_2_KELVIN)

        # Update net radiation for canopy
        rnet -= c.CP * c.AIR_MASS * (new_tleaf_k - tair_k) * grn

        return (new_Tleaf, et, le_et, gbH, gw)
