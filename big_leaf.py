#!/usr/bin/env python
"""
Solve 30-minute coupled A-gs(E) using a big-leaf approximation, i.e. a single
leaf whose Vcmax, Jmax and Rd are scaled to the canopy assuming capacity
declines through the canopy in proportion to the light (Beer's law).
"""

import math
import numpy as np

import constants as c
from farq import FarquharC3
from penman_monteith_leaf import PenmanMonteith, calc_leaf_temp
from radiation import calculate_cos_zenith, calc_leaf_to_canopy_scalar
from radiation import calc_kd, calc_conductance_scalars

__author__  = "Martin De Kauwe"
__version__ = "1.0 (09.11.2018)"
__email__   = "mdekauwe@gmail.com"

class Canopy(object):
    """
    Iteratively solve leaf temp, Ci, gs and An using a big-leaf approach
    """

    def __init__(self, p, peaked_Jmax=True, peaked_Vcmax=True, model_Q10=True,
                 gs_model=None, iter_max=100):

        self.p = p
        self.iter_max = iter_max
        self.F = FarquharC3(peaked_Jmax=peaked_Jmax, peaked_Vcmax=peaked_Vcmax,
                            model_Q10=model_Q10, gs_model=gs_model)
        self.PM = PenmanMonteith()

    def main(self, tair, par, vpd, wind, pressure, Ca, doy, hod, lai,
             rnet=None, Vcmax25=None, Jmax25=None):
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
        rnet : float
            isothermal net radiation (W m-2); estimated from PAR if None

        Returns:
        --------
        an_canopy : float
            net canopy assimilation (umol m-2 s-1)
        gsw_canopy : float
            canopy stomatal conductance to water (mol m-2 s-1)
        et_canopy : float
            transpiration (mol H2O m-2 s-1)
        tcanopy : float
            canopy temperature (deg C)
        """
        p = self.p
        if Vcmax25 is None:
            Vcmax25 = p.Vcmax25
        if Jmax25 is None:
            Jmax25 = p.Jmax25

        # set initial values
        dleaf = vpd
        Cs = Ca
        Tleaf = tair

        (cos_zenith, elevation) = calculate_cos_zenith(doy, p.lat, hod)

        # Calculate big-leaf scaling term to go from a single leaf to canopy
        scalex = calc_leaf_to_canopy_scalar(lai, k=p.k, big_leaf=True)

        # Fraction of radiation absorbed by the canopy (Beer's law)
        fabs = 1.0 - np.exp(-p.k * lai)

        # PAR absorbed by the canopy (umol m-2 ground s-1). Together with
        # Vcmax, Jmax & Rd scaled by scalex this gives canopy-scale An & gsc,
        # consistent with the canopy-scale rnet used in the energy balance.
        apar = par * fabs

        # Scale boundary layer & radiative conductances from a single leaf to
        # the canopy, to match rnet and gsc
        (kd, _) = calc_kd(p, lai)
        scalars = calc_conductance_scalars(lai, p.wind_extinction, kd)

        if rnet is None:
            # isothermal rnet for a closed canopy, reduced for sparse canopies
            tair_k = tair + c.DEG_2_KELVIN
            rnet = self.PM.calc_rnet(p, par, tair, tair_k, vpd) * fabs

        # Is the sun up?
        if elevation > 0.0 and par > 50.0:

            for _ in range(self.iter_max + 1):
                Tleaf_K = Tleaf + c.DEG_2_KELVIN
                (An, gsc) = self.F.photosynthesis(p, Cs=Cs, Tleaf=Tleaf_K,
                                                  Par=apar, vpd=dleaf,
                                                  scalex=scalex,
                                                  Vcmax25=Vcmax25,
                                                  Jmax25=Jmax25)

                # Calculate new Tleaf, dleaf, Cs
                (new_tleaf, et,
                 le_et, gbH, gw) = calc_leaf_temp(p, self.PM, Tleaf, tair,
                                                  gsc, vpd, pressure, wind,
                                                  rnet=rnet, scalars=scalars)

                gbc = gbH * c.GBH_2_GBC
                if gbc > 0.0 and An > 0.0:
                    Cs = Ca - An / gbc # boundary layer of leaf
                else:
                    Cs = Ca

                if np.isclose(et, 0.0) or np.isclose(gw, 0.0):
                    dleaf = vpd
                else:
                    dleaf = (et * pressure / gw) * c.PA_2_KPA # kPa

                converged = math.fabs(Tleaf - new_tleaf) < 0.02

                # Update temperature & do another iteration
                Tleaf = new_tleaf

                if converged:
                    break
            else:
                # No convergence
                An = 0.0
                gsc = 0.0
                et = 0.0

            an_canopy = An
            gsw_canopy = gsc * c.GSC_2_GSW
            et_canopy = et
            tcanopy = Tleaf
        else:
            an_canopy = 0.0
            gsw_canopy = 0.0
            et_canopy = 0.0
            tcanopy = tair

        return (an_canopy, gsw_canopy, et_canopy, tcanopy)
