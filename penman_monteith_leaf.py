#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Isothermal Penman-Monteith

Reference:
==========
* Leuning et al. (1995) Leaf nitrogen, photosynthesis, conductance and
  transpiration: scaling from leaves to canopies
* Wang and Leuning (1998) A two-leaf model for canopy conductance,
  photosynthesis and partitioning of available energy I: Model description and
  comparison with a multi-layered model
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.07.2015)"
__email__ = "mdekauwe@gmail.com"


import math
import numpy as np

from utils import calc_esat
import constants as c


class PenmanMonteith(object):

    def calc_et(self, tair, vpd, pressure, gh, gw, rnet):
        """
        Calculate transpiration following Penman-Monteith at the leaf level
        accounting for effects of leaf temperature and feedback on evaporation.
        For example, if leaf temperature is above the leaf temp, it can increase
        vpd, but it also reduces the lw and thus the net rad available for
        evaporation.

        Parameters:
        ----------
        tair : float
            air temperature (deg C)
        vpd : float
            Vapour pressure deficit (kPa)
        pressure : float
            air pressure (using constant) (Pa)
        gh : float
            boundary layer conductance to heat - free & forced & radiative
            components (mol m-2 s-1)
        gw :float
            conductance to water vapour - stomatal & bdry layer components
            (mol m-2 s-1)
        rnet : float
            Isothermal net radiation (J m-2 s-1 = W m-2)

        Returns:
        --------
        et : float
            transpiration (mol H2O m-2 s-1)
        LE : float
            latent heat flux (W m-2)
        """
        # latent heat of water vapour at air temperature (J mol-1)
        lhv = (c.H2OLV0 - 2.365E3 * tair) * c.H2OMW

        # slope of sat. water vapour pressure (e_sat) to temperature curve
        # (Pa K-1)
        slope = (calc_esat(tair + 0.1) - calc_esat(tair)) / 0.1

        # psychrometric constant (Pa K-1)
        gamma = c.CP * c.AIR_MASS * pressure / lhv

        # Y cancels in eqn 10
        arg1 = slope * rnet + (vpd * c.KPA_2_PA) * gh * c.CP * c.AIR_MASS
        arg2 = slope + gamma * gh / gw

        # W m-2
        LE = max(0.0, arg1 / arg2)

        # transpiration, mol H20 m-2 s-1
        et = LE / lhv

        return (et, LE)

    def calc_conductances(self, p, tair_k, tleaf, tair, wind, gsc, cmolar,
                          scalars=None):
        """
        Both forced and free convection contribute to exchange of heat and mass
        through leaf boundary layers at the wind speeds typically encountered
        within plant canopies (<0-5ms~'). It is particularly imponant to
        includethe contribution of buoyancy forces to the boundary
        conductance for sunlit leaves deep within the canopy where wind speeds
        are low, for without this mechanism computed leaf temperatures become
        excessively high.

        Parameters:
        ----------
        tair_k : float
            air temperature (K)
        tleaf : float
            leaf temperature (deg C)
        tair : float
            air temperature (deg C)
        wind : float
            wind speed (m s-1)
        gsc : float
            stomatal conductance to CO2 (mol m-2 s-1)
        cmolar : float
            Conversion from m s-1 to mol m-2 s-1
        scalars : tuple
            (f_forced, f_free, f_rad) scalars from single leaf to canopy, see
            radiation.calc_conductance_scalars. If None, single leaf values
            are returned, i.e. (1, 1, 2).

        Returns:
        --------
        grn : float
            radiation conductance, both sides (mol m-2 s-1)
        gh : float
            total conductance to heat (mol m-2 s-1), *note* two sided.
        gbH : float
            total boundary layer conductance to heat for one side of the
            leaves
        gw : float
            total leaf conductance to water vapour (mol m-2 s-1)

        References
        ----------
        * Leuning 1995, appendix E
        * Medlyn et al. 2007 appendix, for need for cmolar
        """

        if scalars is None:
            (f_forced, f_free, f_rad) = (1.0, 1.0, 2.0)
        else:
            (f_forced, f_free, f_rad) = scalars

        # radiation conductance, single side of a leaf, Wang and Leuning
        # (1998) just below eqn 9 (NB. units already in mol m-2 s-1). Scaled
        # to both sides/the canopy by f_rad.
        grn = (4.0 * c.SIGMA * tair_k**3 * \
               p.emissivity_leaf) / (c.CP * c.AIR_MASS) * f_rad

        # boundary layer conductance for heat: single sided, forced convection
        # (mol m-2 s-1), wind is at the top of the canopy
        gbHw = 0.003 * math.sqrt(wind / p.leaf_width) * cmolar * f_forced

        if np.isclose(tleaf - tair, 0.0):
            gbHf = 0.0
        else:
            # grashof number
            grashof_num = max(1.e-06, 1.6E8 * math.fabs(tleaf - tair) * \
                                      p.leaf_width**3)

            # boundary layer conductance for heat: single sided, free convection
            # (mol m-2 s-1)
            gbHf = 0.5 * c.DHEAT * grashof_num**0.25 / p.leaf_width * cmolar
            gbHf = max(1.e-06, gbHf * f_free)

        # total boundary layer conductance for heat
        gbH = gbHw + gbHf

        # total conductance for heat (mol m-2 s-1) - two sided
        gh = 2.0 * gbH + grn

        # total leaf conductance for water vapour (mol m-2 s-1)
        gbw = gbH * c.GBH_2_GBW
        gsw = gsc * c.GSC_2_GSW
        gw = (gbw * gsw) / (gbw + gsw)

        return (grn, gh, gbH, gw)

    def calc_rnet(self, p, par, tair, tair_k, vpd):
        """
        Net isothermal radaiation (Rnet, W m-2), i.e. the net radiation that
        would be recieved if leaf and air temperature were the same.

        References:
        ----------
        Jarvis and McNaughton (1986)

        Parameters:
        ----------
        par : float
            Photosynthetically active radiation (umol m-2 s-1)
        tair : float
            air temperature (deg C)
        tair_k : float
            air temperature (K)
        vpd : float
            Vapour pressure deficit (kPa)

        Returns:
        --------
        rnet : float
            Net radiation (J m-2 s-1 = W m-2)

        """

        # Short wave radiation (W m-2), i.e. VIS + NIR. NB. par / 4.6 would
        # only be the VIS half of the SW.
        sw_rad = par * c.PAR_2_SW

        # atmospheric water vapour pressure (Pa)
        ea = max(0.0, calc_esat(tair) - (vpd * c.KPA_2_PA))

        # apparent emissivity for a hemisphere radiating at air temperature
        # eqn D4
        emissivity_atm = 0.642 * (ea / tair_k)**(1.0 / 7.0)

        # isothermal net LW radiaiton at top of canopy, assuming emissivity of
        # the canopy is 1
        net_lw_rad = (1.0 - emissivity_atm) * c.SIGMA * tair_k**4

        # isothermal net radiation (W m-2)
        rnet = p.SW_abs * sw_rad - net_lw_rad

        return rnet


def calc_leaf_temp(p, PM, tleaf, tair, gsc, vpd, pressure, wind, rnet,
                   scalars=None):
    """
    Resolve leaf temp

    Parameters:
    ----------
    p : module
        model parameters
    PM : object
        Penman-Montheith class instance
    tleaf : float
        leaf temperature (deg C)
    tair : float
        air temperature (deg C)
    gsc : float
        stomatal conductance to CO2 (mol m-2 s-1)
    vpd : float
        Vapour pressure deficit (kPa)
    pressure : float
        air pressure (using constant) (Pa)
    wind : float
        wind speed (m s-1)
    rnet : float
        isothermal net radiation (W m-2)
    scalars : tuple
        (f_forced, f_free, f_rad) leaf to canopy conductance scalars

    Returns:
    --------
    new_tleaf : float
        new leaf temperature (deg C)
    et : float
        transpiration (mol H2O m-2 s-1)
    le_et : float
        latent heat flux (W m-2)
    gbH : float
        total boundary layer conductance to heat for one side of the leaf
    gw : float
        total leaf conductance to water vapour (mol m-2 s-1)
    """
    tair_k = tair + c.DEG_2_KELVIN

    air_density = pressure / (c.RSPECIFC_DRY_AIR * tair_k)

    # convert from m s-1 to mol m-2 s-1
    cmolar = pressure / (c.RGAS * tair_k)

    (grn, gh, gbH, gw) = PM.calc_conductances(p, tair_k, tleaf, tair,
                                              wind, gsc, cmolar, scalars)

    if np.isclose(gsc, 0.0):
        et = 0.0
        le_et = 0.0
    else:
        (et, le_et) = PM.calc_et(tair, vpd, pressure, gh, gw, rnet)

    # Leaf-air temperature difference from the energy balance using the
    # isothermal net radiation, Leuning 1995, appendix D/E:
    #   dT = (Rn_iso - LE) / (cp * gh), where gh = 2 * gbH + grn
    # gh already includes the radiative conductance, so the Y factor
    # (gbH / (gbH + grn)) must not also be applied here - it only partitions
    # (Rn_iso - LE) into sensible heat, H = Y * (Rn_iso - LE). Dividing H by
    # gh double counts grn and underestimates dT.
    new_tleaf = tair + (rnet - le_et) / (c.CP * air_density * (gh / cmolar))

    return (new_tleaf, et, le_et, gbH, gw)
