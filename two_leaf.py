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

import math
import numpy as np

import constants as c
from farq import FarquharC3
from penman_monteith_leaf import PenmanMonteith, calc_leaf_temp
from radiation import spitters
from radiation import calculate_absorbed_radiation
from radiation import calculate_cos_zenith, calc_leaf_to_canopy_scalar
from radiation import calc_conductance_scalars

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.11.2018)"
__email__   = "mdekauwe@gmail.com"


class Canopy(object):
    """Iteratively solve leaf temp, Ci, gs and An."""

    def __init__(self, p, peaked_Jmax=True, peaked_Vcmax=True, model_Q10=True,
                 gs_model=None, iter_max=100):

        self.p = p
        self.iter_max = iter_max
        self.F = FarquharC3(peaked_Jmax=peaked_Jmax, peaked_Vcmax=peaked_Vcmax,
                            model_Q10=model_Q10, gs_model=gs_model)
        self.PM = PenmanMonteith()

    def main(self, tair, par, vpd, wind, pressure, Ca, doy, hod, lai,
             Vcmax25=None, Jmax25=None):
        """
        Parameters:
        ----------
        tair : float
            air temperature (deg C)
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        vpd : float
            Vapour pressure deficit (kPa)
        wind : float
            wind speed (m s-1)
        pressure : float
            air pressure (using constant) (Pa)
        Ca : float
            ambient CO2 concentration
        doy : float
            day of year
        hod : float
            hour of day
        lai : float
            leaf area index

        Returns:
        --------
        An : array
            net assimilation, sunlit/shaded (umol m-2 s-1)
        et : array
            transpiration, sunlit/shaded (mol H2O m-2 s-1)
        Tcan : array
            leaf temperature, sunlit/shaded (deg C)
        apar : array
            absorbed PAR, sunlit/shaded (umol m-2 s-1)
        lai_leaf : array
            leaf area, sunlit/shaded (m2 m-2)
        """
        p = self.p
        if Vcmax25 is None:
            Vcmax25 = p.Vcmax25
        if Jmax25 is None:
            Jmax25 = p.Jmax25

        An = np.zeros(2)            # sunlit, shaded
        gsc = np.zeros(2)           # sunlit, shaded
        et = np.zeros(2)            # sunlit, shaded
        Tcan = np.full(2, tair)     # sunlit, shaded; leaves at tair at night
        sw_rad = np.zeros(2)        # VIS, NIR

        (cos_zenith, elevation) = calculate_cos_zenith(doy, p.lat, hod)

        sw_rad[c.VIS] = 0.5 * (par * c.PAR_2_SW) # W m-2
        sw_rad[c.NIR] = 0.5 * (par * c.PAR_2_SW) # W m-2

        # get diffuse/beam frac, just use VIS as the answer is the same for NIR
        (diffuse_frac, direct_frac) = spitters(doy, sw_rad[c.VIS], cos_zenith)

        (qcan, apar,
         lai_leaf, kb, kd) = calculate_absorbed_radiation(p, par, cos_zenith,
                                                      lai, direct_frac,
                                                      diffuse_frac, doy, sw_rad,
                                                      tair)

        # Calculate scaling term to go from a single leaf to canopy,
        # see Wang & Leuning 1998 appendix C
        scalex = calc_leaf_to_canopy_scalar(lai, kb=kb, kn=p.kn)

        if lai_leaf[c.SUNLIT] < 1.e-3: # to match line 336 of CABLE radiation
            scalex[c.SUNLIT] = 0.

        # Scale boundary layer & radiative conductances from a single leaf to
        # the sunlit/shaded canopy, to match qcan and gsc
        (f_forced, f_free,
         f_rad) = calc_conductance_scalars(lai, p.wind_extinction, kd, kb,
                                           lai_leaf)

        # Is the sun up?
        if elevation > 0.0 and par > 50.:

            # sunlit / shaded loop
            for ileaf in range(2):

                # no leaves in this fraction, e.g. no sunlit leaves
                if lai_leaf[ileaf] < 1.e-3:
                    continue

                # initialise values of Tleaf, Cs, dleaf at the leaf surface
                dleaf = vpd
                Cs = Ca
                Tleaf = tair

                for _ in range(self.iter_max + 1):

                    if scalex[ileaf] > 0.:
                        Tleaf_K = Tleaf + c.DEG_2_KELVIN
                        (An[ileaf],
                         gsc[ileaf]) = self.F.photosynthesis(p, Cs=Cs,
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
                     le_et, gbH, gw) = calc_leaf_temp(p, self.PM, Tleaf, tair,
                                                      gsc[ileaf], vpd,
                                                      pressure, wind,
                                                      rnet=qcan[ileaf],
                                                      scalars=(f_forced[ileaf],
                                                               f_free[ileaf],
                                                               f_rad[ileaf]))

                    gbc = gbH * c.GBH_2_GBC
                    if gbc > 0.0 and An[ileaf] > 0.0:
                        Cs = Ca - An[ileaf] / gbc # boundary layer of leaf
                    else:
                        Cs = Ca

                    if np.isclose(et[ileaf], 0.0) or np.isclose(gw, 0.0):
                        dleaf = vpd
                    else:
                        dleaf = (et[ileaf] * pressure / gw) * c.PA_2_KPA # kPa

                    converged = math.fabs(Tleaf - new_tleaf) < 0.02

                    # Update temperature & do another iteration
                    Tleaf = new_tleaf
                    Tcan[ileaf] = Tleaf

                    if converged:
                        break
                else:
                    # No convergence
                    An[ileaf] = 0.0
                    gsc[ileaf] = 0.0
                    et[ileaf] = 0.0

        return (An, et, Tcan, apar, lai_leaf)
