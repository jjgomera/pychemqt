#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from unittest import TestCase

# Noble Gases
from lib.mEoS.He import He
from lib.mEoS.Ne import Ne
from lib.mEoS.Ar import Ar
from lib.mEoS.Kr import Kr
from lib.mEoS.Xe import Xe

# Gases
from lib.mEoS.H2 import H2
from lib.mEoS.D2 import D2, pD2, oD2
from lib.mEoS.pH2 import pH2
from lib.mEoS.oH2 import oH2
from lib.mEoS.N2 import N2
from lib.mEoS.O2 import O2
from lib.mEoS.F2 import F2
from lib.mEoS.H2O import H2O
from lib.mEoS.D2O import D2O
from lib.mEoS.CO2 import CO2
from lib.mEoS.CO import CO
from lib.mEoS.N2O import N2O
from lib.mEoS.SO2 import SO2
from lib.mEoS.COS import COS
from lib.mEoS.NH3 import NH3
from lib.mEoS.H2S import H2S

# Alkanes
from lib.mEoS.CH4 import CH4
from lib.mEoS.C2 import C2
from lib.mEoS.C3 import C3
from lib.mEoS.nC4 import nC4
from lib.mEoS.iC4 import iC4
from lib.mEoS.nC5 import nC5
from lib.mEoS.neoC5 import neoC5
from lib.mEoS.iC5 import iC5
from lib.mEoS.nC6 import nC6
from lib.mEoS.iC6 import iC6
from lib.mEoS.nC7 import nC7
from lib.mEoS.nC8 import nC8
from lib.mEoS.iC8 import iC8
from lib.mEoS.nC9 import nC9
from lib.mEoS.nC10 import nC10
from lib.mEoS.nC11 import nC11
from lib.mEoS.nC12 import nC12
from lib.mEoS.nC16 import nC16
from lib.mEoS.nC22 import nC22
from lib.mEoS.C3_pentane import C3_pentane
from lib.mEoS.C22_butane import C22_butane
from lib.mEoS.C23_butane import C23_butane

# Naphthenes
from lib.mEoS.Cyclopropane import Cyclopropane
from lib.mEoS.Cyclopentane import Cyclopentane
from lib.mEoS.Cyclohexane import Cyclohexane
from lib.mEoS.C1Cyclohexane import C1Cyclohexane
from lib.mEoS.C3Cyclohexane import C3Cyclohexane

# Alkenes
from lib.mEoS.Benzene import Benzene
from lib.mEoS.Toluene import Toluene
from lib.mEoS.oXylene import oXylene
from lib.mEoS.mXylene import mXylene
from lib.mEoS.pXylene import pXylene
from lib.mEoS.EthylBenzene import EthylBenzene
from lib.mEoS.Ethylene import Ethylene
from lib.mEoS.Propylene import Propylene
from lib.mEoS.Butene_1 import Butene_1
from lib.mEoS.iButene import iButene
from lib.mEoS.Cis_2_butene import Cis_2_butene
from lib.mEoS.Trans_2_butene import Trans_2_butene
from lib.mEoS.Propyne import Propyne
from lib.mEoS.C1Oleate import C1Oleate
from lib.mEoS.C1Linolenate import C1Linolenate
from lib.mEoS.C1Linoleate import C1Linoleate
from lib.mEoS.C1Palmitate import C1Palmitate
from lib.mEoS.C1Stearate import C1Stearate

# Heteroatom
from lib.mEoS.Methanol import Methanol
from lib.mEoS.Ethanol import Ethanol
from lib.mEoS.Acetone import Acetone
from lib.mEoS.EthyOxide import EthyOxide
from lib.mEoS.PropylenGlycol import PropylenGlycol
from lib.mEoS.DME import DME
from lib.mEoS.DEE import DEE
from lib.mEoS.DMC import DMC
from lib.mEoS.NF3 import NF3
from lib.mEoS.SF6 import SF6
from lib.mEoS.HCl import HCl

# CFCs
from lib.mEoS.R13I1 import R13I1
from lib.mEoS.R11 import R11
from lib.mEoS.R12 import R12
from lib.mEoS.R13 import R13
from lib.mEoS.R14 import R14
from lib.mEoS.R21 import R21
from lib.mEoS.R22 import R22
from lib.mEoS.R23 import R23
from lib.mEoS.R32 import R32
from lib.mEoS.R40 import R40
from lib.mEoS.R41 import R41
from lib.mEoS.R113 import R113
from lib.mEoS.R114 import R114
from lib.mEoS.R115 import R115
from lib.mEoS.R116 import R116
from lib.mEoS.R123 import R123
from lib.mEoS.R124 import R124
from lib.mEoS.R125 import R125
from lib.mEoS.R134a import R134a
from lib.mEoS.R141b import R141b
from lib.mEoS.R142b import R142b
from lib.mEoS.R143a import R143a
from lib.mEoS.R152a import R152a
from lib.mEoS.R161 import R161
from lib.mEoS.R218 import R218
from lib.mEoS.R227ea import R227ea
from lib.mEoS.R236ea import R236ea
from lib.mEoS.R236fa import R236fa
from lib.mEoS.R245ca import R245ca
from lib.mEoS.R245fa import R245fa
from lib.mEoS.R365mfc import R365mfc
from lib.mEoS.RC318 import RC318
from lib.mEoS.R1123 import R1123
from lib.mEoS.R1234yf import R1234yf
from lib.mEoS.R1234zeE import R1234zeE
from lib.mEoS.R1234zeZ import R1234zeZ
from lib.mEoS.R1243zf import R1243zf
from lib.mEoS.R1216 import R1216
from lib.mEoS.R1233zd import R1233zd
from lib.mEoS.R1336mzzZ import R1336mzzZ
from lib.mEoS.RE143a import RE143a
from lib.mEoS.RE245cb2 import RE245cb2
from lib.mEoS.RE245fa2 import RE245fa2
from lib.mEoS.RE347mcc import RE347mcc
from lib.mEoS.Novec649 import Novec649

# Siloxanes
from lib.mEoS.D4 import D4
from lib.mEoS.D5 import D5
from lib.mEoS.D6 import D6
from lib.mEoS.MDM import MDM
from lib.mEoS.MD2M import MD2M
from lib.mEoS.MD3M import MD3M
from lib.mEoS.MD4M import MD4M
from lib.mEoS.MM import MM

# PseudoCompound
from lib.mEoS.Air import Air
from lib.mEoS.R404a import R404a
from lib.mEoS.R407c import R407c
from lib.mEoS.R410a import R410a
from lib.mEoS.R507a import R507a

# Others
from lib.mEoS.LJ import LJ


# Component grouping by chemical class
Nobles = [He, Ne, Ar, Kr, Xe]
Gases = [H2, D2, pD2, oD2, pH2, oH2, N2, O2, F2, H2O, D2O, CO2, CO, N2O, SO2,
         COS, NH3, H2S]
Alkanes = [CH4, C2, C3, nC4, iC4, nC5, neoC5, iC5, nC6, iC6, nC7, nC8, iC8,
           nC9, nC10, nC11, nC12, nC16, nC22, C3_pentane, C22_butane,
           C23_butane]
Naphthenes = [Cyclopropane, Cyclopentane, Cyclohexane, C1Cyclohexane,
              C3Cyclohexane]
Alkenes = [Benzene, Toluene, oXylene, mXylene, pXylene, EthylBenzene,
           Ethylene, Propylene, Butene_1, iButene, Cis_2_butene,
           Trans_2_butene, Propyne, C1Oleate, C1Linolenate, C1Linoleate,
           C1Palmitate, C1Stearate]
Heteroatom = [Methanol, Ethanol, Acetone, EthyOxide, PropylenGlycol, DME, DEE,
              DMC, NF3, SF6, HCl]
CFCs = [R13I1, R11, R12, R13, R14, R21, R22, R23, R32, R40, R41, R113, R114,
        R115, R116, R123, R124, R125, R134a, R141b, R142b, R143a, R152a, R161,
        R218, R227ea, R236ea, R236fa, R245ca, R245fa, R365mfc, RC318, R1123,
        R1234yf, R1234zeE, R1234zeZ, R1243zf, R1216, R1233zd, R1336mzzZ,
        RE143a, RE245cb2, RE245fa2, RE347mcc, Novec649]
Siloxanes = [D4, D5, D6, MDM, MD2M, MD3M, MD4M, MM]
PseudoCompounds = [Air, R404a, R407c, R410a, R507a]
Others = [LJ]

__all__ = Nobles + Gases + Alkanes + Naphthenes + Alkenes + Heteroatom + \
    CFCs + Siloxanes + PseudoCompounds + Others


# Id of compound supported for meos library
id_mEoS = [i.id for i in __all__]


# Add references from equation hardcoded in __doi__ property
__doi__ = {}
for obj in __all__:
    subdict = {}
    for prop in ["eq", "_viscosity", "_thermal"]:
        if prop not in obj.__dict__ or not obj.__dict__[prop]:
            continue
        for i, eq in enumerate(obj.__dict__[prop]):
            if eq and "__doi__" in eq:
                key = "%s_%i" % (prop.replace("_", ""), i)
                subdict[key] = eq["__doi__"]
    if obj._surface and "__doi__" in obj._surface:
        subdict["surface"] = obj._surface["__doi__"]
    if obj._dielectric and "__doi__" in obj._dielectric:
        subdict["dielectric"] = obj._dielectric["__doi__"]
    if obj._melting and "__doi__" in obj._melting:
        subdict["melting"] = obj._melting["__doi__"]
    if obj._sublimation and "__doi__" in obj._sublimation:
        subdict["sublimation"] = obj._sublimation["__doi__"]

    __doi__[obj.__name__] = subdict


# TODO: Add 1-propanol from 10.1016_j.fluid.2004.06.028
# TODO: Add 2-propanol from 10.1063/1.3112608
# TODO: Add 1,2-dichloroethane from Thol thesis


class Test(TestCase):
    def test_meos(self):
        """Cycle input parameter from selected point to check iteration"""
        # The input pair T-h, P-s, h-u has inconsistency, several point has
        # equal values so are not good as input definition, they need another
        # input like saturation state

        # Subcooled liquid region
        P = 5e7
        T = 470
        f_tp = H2O(T=T, P=P)
        f_trho = H2O(T=f_tp.T, rho=f_tp.rho)
        f_ts = H2O(T=f_trho.T, s=f_trho.s)
        f_tu = H2O(T=f_ts.T, u=f_ts.u)
        f_prho = H2O(P=f_tu.P, rho=f_tu.rho)
        f_ph = H2O(P=f_prho.P, h=f_prho.h)
        f_pu = H2O(P=f_ph.P, u=f_ph.u)
        f_rhoh = H2O(rho=f_pu.rho, h=f_pu.h)
        f_rhos = H2O(rho=f_rhoh.rho, s=f_rhoh.s)
        f_rhou = H2O(rho=f_rhos.rho, u=f_rhos.u)
        f_hs = H2O(h=f_rhou.h, s=f_rhou.s)
        f_su = H2O(s=f_hs.s, u=f_hs.u)
        self.assertEqual(round(f_su.P-P, 1), 0)
        self.assertEqual(round(f_su.T-T, 5), 0)

        # Two phases region
        T = 470
        x = 0.3
        f_tx = H2O(T=T, x=x)
        f_px = H2O(P=f_tx.P, x=f_tx.x)
        f_trho = H2O(T=f_px.T, rho=f_px.rho)
        f_ts = H2O(T=f_trho.T, s=f_trho.s)
        f_tu = H2O(T=f_ts.T, u=f_ts.u)
        f_prho = H2O(P=f_tu.P, rho=f_tu.rho)
        f_ph = H2O(P=f_prho.P, h=f_prho.h)
        f_pu = H2O(P=f_ph.P, u=f_ph.u)
        f_rhoh = H2O(rho=f_pu.rho, h=f_pu.h)
        f_rhos = H2O(rho=f_rhoh.rho, s=f_rhoh.s)
        f_rhou = H2O(rho=f_rhos.rho, u=f_rhos.u)
        f_hs = H2O(h=f_rhou.h, s=f_rhou.s)
        f_su = H2O(s=f_hs.s, u=f_hs.u)
        self.assertEqual(round(f_su.x-x, 5), 0)
        self.assertEqual(round(f_su.T-T, 5), 0)
        self.assertEqual(round(f_su.P-f_tx.P, 1), 0)

        # Superheated gas region
        P = 5e5
        T = 470
        f_tp = H2O(T=T, P=P)
        f_trho = H2O(T=f_tp.T, rho=f_tp.rho)
        f_ts = H2O(T=f_trho.T, s=f_trho.s)
        f_tu = H2O(T=f_ts.T, u=f_ts.u)
        f_prho = H2O(P=f_tu.P, rho=f_tu.rho)
        f_ph = H2O(P=f_prho.P, h=f_prho.h)
        f_pu = H2O(P=f_ph.P, u=f_ph.u)
        f_rhoh = H2O(rho=f_pu.rho, h=f_pu.h)
        f_rhos = H2O(rho=f_rhoh.rho, s=f_rhoh.s)
        f_rhou = H2O(rho=f_rhos.rho, u=f_rhos.u)
        f_hs = H2O(h=f_rhou.h, s=f_rhou.s)
        f_su = H2O(s=f_hs.s, u=f_hs.u)
        self.assertEqual(round(f_su.P-P, 1), 0)
        self.assertEqual(round(f_su.T-T, 5), 0)

        # Supercritical region
        P = 8e7
        T = 800
        f_tp = H2O(T=T, P=P)
        f_trho = H2O(T=f_tp.T, rho=f_tp.rho)
        f_ts = H2O(T=f_trho.T, s=f_trho.s)
        f_tu = H2O(T=f_ts.T, u=f_ts.u)
        f_prho = H2O(P=f_tu.P, rho=f_tu.rho)
        f_ph = H2O(P=f_prho.P, h=f_prho.h)
        f_pu = H2O(P=f_ph.P, u=f_ph.u)
        f_rhoh = H2O(rho=f_pu.rho, h=f_pu.h)
        f_rhos = H2O(rho=f_rhoh.rho, s=f_rhoh.s)
        f_rhou = H2O(rho=f_rhos.rho, u=f_rhos.u)
        f_hs = H2O(h=f_rhou.h, s=f_rhou.s)
        f_su = H2O(s=f_hs.s, u=f_hs.u)
        self.assertEqual(round(f_su.P-P, 1), 0)
        self.assertEqual(round(f_su.T-T, 5), 0)
