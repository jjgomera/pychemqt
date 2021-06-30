#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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


from math import exp
from unittest import TestCase

from scipy.constants import Boltzmann, Avogadro, pi

from lib import unidades
from lib.meos import MEoS


class Ethylene(MEoS):
    """Multiparameter equation of state for ehylene"""
    name = "ethylene"
    CASNumber = "74-85-1"
    formula = "CH2=CH2"
    synonym = "R-1150"
    _refPropName = "ETHYLENE"
    _coolPropName = "Ethylene"
    rhoc = unidades.Density(214.24)
    Tc = unidades.Temperature(282.35)
    Pc = unidades.Pressure(5041.8, "kPa")
    M = 28.05376  # g/mol
    Tt = unidades.Temperature(103.989)
    Tb = unidades.Temperature(169.379)
    f_acent = 0.0866
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 22
    _Tr = unidades.Temperature(273.316763)
    _rhor = unidades.Density(216.108926)
    _w = 0.085703183

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           # Using a custom integration parameters for reference state
           # "ao_pow": [8.68815523, -4.47960564],
           "ao_pow": [8.67604219, -4.46075924323],
           "ao_exp": [2.49395851, 3.0027152, 2.5126584, 3.99064217],
           "titao": [4.43266896, 5.74840149, 7.8027825, 15.5851154]}

    CP1 = {"ao": 4.,
           "ao_exp": [1]*12,
           "exp": [4353.907145, 2335.2251475, 1930.913215, 1471.9256475,
                   4464.6972475, 1778.39697, 1365.4520425, 1356.8190475,
                   4469.013745, 1188.475645, 4300.6703425, 2077.67413]}

    CP2 = {"ao": 3.554495281,
           "an": [0.5603615762e6, -0.2141069802e5, 0.2532008897e3,
                  -0.9951927478e-2, 0.5108931070e-4, -0.1928667482e-7],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-0.2061703241e2], "exp": [3000]}

    smukala = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Smukala et "
                    "al. (2000)",
        "__doi__": {"autor": "Smukala, J., Span, R., Wagner, W.",
                    "title": "New equation of state for ethylene covering the "
                             "fluid region from the melting line to 450 K at "
                             "pressures up to 300 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29(5) (2000) 1053-1121",
                    "doi": "10.1063/1.1329318"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 300000.0, "rhomax": 27.03,

        "nr1": [1.861742910067, -3.0913708460844, -0.17384817095516,
                0.08037098569284, 0.23682707317354, 0.021922786610247],
        "d1": [1, 1, 1, 2, 2, 4],
        "t1": [0.5, 1, 2.5, 0, 2, 0.5],

        "nr2": [0.11827885813193, -0.021736384396776, 0.044007990661139,
                0.12554058863881, -0.13167945577241, -0.0052116984575897,
                0.00015236081265419, -2.4505335342756e-05, 0.28970524924022,
                -0.18075836674288, 0.15057272878461, -0.14093151754458,
                0.022755109070253, 0.014026070529061, 0.0061697454296214,
                -0.00041286083451333, 0.012885388714785, -0.069128692157093,
                0.10936225568483, -0.0081818875271794, -0.05641847211717,
                0.0016517867750633, 0.0095904006517001, -0.0026236572984886],
        "d2": [1, 1, 3, 4, 5, 7, 10, 11, 1, 1, 2, 2, 4, 4, 6, 7, 4, 5, 6, 6,
               7, 8, 9, 10],
        "t2": [1., 4., 1.25, 2.75, 2.25, 1., 0.75, 0.5, 2.5, 3.5, 4., 6., 1.5,
               5., 4.5, 15., 20., 23., 22., 29., 19., 15., 13., 10.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4,
               4, 4, 4],
        "gamma2": [1]*24,

        "nr3": [-50.242414011355, 7484.6420119299, -6873.4299232625,
                -935.77982814338, 941.33024786113],
        "d3": [2, 2, 2, 3, 3],
        "t3": [1., 0., 1., 2., 3.],
        "alfa3": [25.]*5,
        "beta3": [325, 300, 300, 300, 300],
        "gamma3": [1.16, 1.19, 1.19, 1.19, 1.19],
        "epsilon3": [1.]*5}

    jahangiri = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Jahangiri "
                    "(1986)",
        "__doi__": {"autor": "Jahangiri, M., Jacobsen, R.T, Stewart, R.B., "
                             "McCarty, R.D.",
                    "title": "Thermodynamic properties of ethylene from the "
                             "freezing line to 450 K at pressures to 260 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 15(2) (1986) 293-734",
                    "doi": "10.1063/1.555753"},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 29645.9, "so": 219.39},
        "M": 28.054, "Tc": 282.3452, "rhoc": 7.634,

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 260000.0, "rhomax": 26.67,

        "nr1": [3.248937034, -10.17278862, 7.386604053, -1.568916359,
                -0.08884514287, 0.06021068143, 0.1078324588, -0.02004025211,
                0.001950491412, 0.06718006403, -0.04200451469, -0.001620507626,
                0.0005555156795, 0.0007583671146, -0.0002878544021],
        "d1": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 6, 6, 6],
        "t1": [0.5, 1, 1.25, 1.75, 4, 2, 4, 5, 6, 0.25, 3, 0.25, 0.5, 2.5, 3],

        "nr2": [0.06258987063, -0.06418431160, -0.1368693752, 0.5179207660,
                -0.3026331319, 0.7757213872, -2.639890864, 2.927563554,
                -1.066267599, -0.05380471540, 0.1277921080, -0.07450152310,
                -0.01624304356, 0.1476032429, -0.2003910489, 0.2926905618,
                -0.1389040901, 5.913513541, -38.00370130, 96.91940570,
                -122.6256839, 77.02379476, -19.22684672, -0.003800045701,
                0.01118003813, 0.002945841426],
        "d2": [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4,
               4, 4, 8, 8, 8],
        "t2": [0.5, 1, 0.5, 2, 4, 3, 4, 5, 6, 2, 3, 4, 1.5, 0.5, 1.5, 4, 5, 1,
               2, 3, 4, 5, 6, 0.5, 1, 5],
        "c2": [3, 3, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 3, 2, 2, 2, 2, 4, 4, 4, 4,
               4, 4, 2, 2, 2],
        "gamma2": [1]*26}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ethylene of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 27.03,
        "M": 28.054,

        "nr1": [0.9096223, -0.24641015e1, 0.56175311, -0.19688013e-1,
                0.78831145e-1, 0.21478776e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.23151337, -0.37804454e-1, -0.20122739, -0.44960157e-1,
                -0.2834296e-1, 0.12652824e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    mccarty = {
        # Also referenced in
        # Younglove, B.A.",
        # Thermophysical Properties of Fluids. I. Argon, Ethylene,
        # Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen
        # J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)
        # Younglove use a different rhoc so the tabulated values may differ
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethylene of McCarty (1981)",
        "__doi__": {"autor": "McCarty, R.D., Jacobsen, R.T.",
                    "title": "An Equation of State for Fluid Ethylene",
                    "ref": "Natl. Bur. Stand., Tech. Note 1045, 1981.",
                    "doi": ""},

        "R": 8.31434, "M": 28.054,
        "Tt": 103.986, "Tc": 282.3428, "Pc": 5.0403, "rhoc": 7.634,
        "cp": CP2,
        "ref": {"Tref": 300, "Pref": 101.325, "so": 219.4, "ho": 29646.46},

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 40000.0, "rhomax": 23.343,

        "gamma": -0.0172,
        "b": [None, -0.2146684366683e-1, 0.1791433722534e1, -0.3675315603930e2,
              0.3707178934669e4, -0.3198282566709e6, 0.5809379774732e-3,
              -0.7895570824899, 0.1148620375835e3, 0.2713774629193e6,
              -0.8647124319107e-4, 0.1617727266385, -0.2731527496271e2,
              -0.2672283641459e-2, -0.4752381331990e-1, -0.6255637346217e2,
              0.4576234964434e-2, -0.7534839269320e-4, 0.1638171982209,
              -0.3563090740740e-2, -0.1833000783170e6, -0.1805074209985e8,
              -0.4794587918874e4, 0.3531948274957e8, -0.2562571039155e2,
              0.1044308253292e4, -0.1695303363659, -0.1710334224958e4,
              -0.2054114462372e-3, 0.6727558766661e-1, -0.1557168403328e-5,
              -0.1229814736077e-3, 0.4234325938573e-3]}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 29610, "so": 219.225},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [8.42278605e-1, 8.65139678e-1, -2.79801027, 6.74520156e-2,
                2.42445468e-4, -2.74767618e-3],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.48602227e-2, 1.29307481e-1, 3.74759088e-1, -1.25336440e-2,
                -2.33507187e-1, 1.38862785e-2, -4.88033330e-2, -2.38141707e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = smukala, mccarty, jahangiri, shortSpan, sun

    # Unimplemented mBWR equation
    # Bender, E.
    # Equations of state for Ethylene and propylene
    # Cryogenics (1975) 667-673
    # 10.1016/0011-2275(75)90100-9

    _PR = [-0.2301, -15.7070]

    _surface = {"sigma": [0.0477], "exp": [1.17]}
    _dielectric = {
        "eq": 1,
        "a": [10.725, 0], "b": [55.19, 49.5], "c": [-2045., -1154.],
        "Au": 0, "D": 1.9}

    _melting = {
        "__doi__": smukala["__doi__"],
        "Tmin": Tt, "Tmax": 550.0}

    @classmethod
    def _Melting_Pressure(cls, T):
        """Calculate the melting pressure using the fractional method of
        Smukala"""
        Tt2 = 110.369
        if T < Tt2:
            a = 2947001.84
            t = 2.045
            Tt = cls.Tt
            Pt = 122.65
        else:
            a = 6.82693421
            t = 1.089
            Tt = Tt2
            Pt = 46.8e6

        Pr = 1 + a*((T/Tt)**t-1)
        return unidades.Pressure(Pr*Pt)

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.3905741, 1.4060338, -1.6589923, 1.0278028, -2.5071716],
        "t": [1.0, 1.5, 2.5, 3.0, 4.5]}
    _liquid_Density = {
        "eq": 2,
        "n": [1.8673079, -0.61533892, -0.058973772, 0.10744720],
        "t": [0.343, 3/6, 8/6, 12/6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.9034556, -0.75837929, -3.7717969, -8.7478586, -23.885296,
              -54.197979],
        "t": [0.349, 4/6, 1, 14/6, 29/6, 56/6]}

    visco0 = {"__name__": "Holland (1983)",
              "__doi__": {
                  "autor": "Holland, P.M., Eaton, B.E., Hanley, H.J.M.",
                  "title": "A Correlation of the Viscosity and Thermal "
                           "Conductivity Data of Gaseous and Liquid Ethylene",
                  "ref": "J. Phys. Chem. Ref. Data 12(4) (1983) 917-932",
                  "doi": "10.1063/1.555701"},

              "eq": 1, "omega": 0,

              "no": [-3.5098225018e5, 2.5008406184e5, -5.8365540744e4,
                     4.5549146583e2, 2.2881683403e3, -4.7318682077e2,
                     4.5022249258e1, -2.1490688088, 4.1649263233e-2],
              "to": [-1, -2/3, -1/3, 0, 1/3, 2/3, 1, 4/3, 5/3],

              "special": "_mur"}

    def _mur(self, rho, T, fase):
        """Density correction for viscosity correlation"""
        # η1 in Eq 3 is always 0

        # Eq 4
        tita = (rho-221)/221
        j = [-4.8544486732, 1.3033585236e1, 2.7808928908e4, -1.8241971308e3,
             1.5913024509, -2.0513573927e2, -3.9478454708e4]
        mu2 = exp(j[0]+j[3]/T) * (exp(rho.gcc**0.1*(j[1]+j[2]/T**1.5) +
                                  tita*rho.gcc**0.5*(j[4]+j[5]/T+j[6]/T**2))-1)

        # The reurned values is in microP, convert to μPas
        return mu2/10

    _viscosity = visco0,

    thermo0 = {"__name__": "Assael (2016)",
               "__doi__": {
                   "autor": "Assael, M.J., Koutian, A., Huber, M.L., Perkins, "
                            "R.A.",
                   "title": "Reference Correlations of the Thermal "
                            "Conductivity of Ethene and Propene",
                   "ref": "J. Phys. Chem. Ref. Data 45(3) (2016) 033104",
                   "doi": "10.1063/1.4958984"},

               "eq": 1,

               "Toref": 282.35, "koref": 1e-3,
               "no_num": [-54.1761, 541.904, -656.108, 667.048, -109.992,
                          60.6511, -1.01377],
               "to_num": [0, 1, 2, 3, 4, 5, 6],
               "no_den": [26.5363, -20.1401, 19.4152, -2.92695, 1],
               "to_den": [0, 1, 2, 3, 4],

               "Tref_res": 282.35, "rhoref_res": 214.24, "kref_res": 1e-3,
               "nr": [0.261453e2, -0.218619e2, 0.362068e2, -0.136642e2,
                      0.184752e1, -0.113225e2, 0.269282e2, -0.223164e2,
                      0.390241e1, 0.668286],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.181e-9,
               "gam0": 0.058, "qd": 0.49e-9, "Tcref": 423.53}

    thermo1 = {"__name__": "Holland (1983)",
               "__doi__": {
                   "autor": "Holland, P.M., Eaton, B.E., Hanley, H.J.M.",
                   "title": "A Correlation of the Viscosity and Thermal "
                            "Conductivity Data of Gaseous and Liquid Ethylene",
                   "ref": "J. Phys. Chem. Ref. Data 12(4) (1983) 917-932",
                   "doi": "10.1063/1.555701"},

               "eq": 0,
               "method": "_thermo0"}

    def _thermo0(self, rho, T, fase):
        # λ1 in Eq 3 is always 0

        GT = [-2.903423528e5, 4.680624952e5, -1.8954783215e5, -4.8262235392e3,
              2.243409372e4, -6.6206354818e3, 8.9937717078e2, -6.0559143718e1,
              1.6370306422]
        lo = 0
        for i in range(-3, 6):
            lo += GT[i+3]*T**(i/3.)

        l2, lc = 0, 0
        if rho:
            tita = (rho-221)/221
            k = [-1.304503323e1, 1.8214616599e1, -9.903022496e3, 7.420521631e2,
                 -3.0083271933e-1, 9.6456068829e1, 1.350256962e4]
            l2 = exp(k[0]+k[3]/T) * (
                exp(rho.gcc**0.1*(k[1]+k[2]/T**1.5) +
                    tita*rho.gcc**0.5*(k[4]+k[5]/T+k[6]/T**2))-1)

            # Critical enhancement
            deltarho = (rho-221)/221
            deltaT = (T-282.34)/282.34

            xt = rho**2*fase.kappa*5.039/221**2
            B = abs(deltarho)/abs(deltaT)**1.19                         # Eq 11
            Gamma = xt*abs(deltaT)**1.19                                # Eq 12
            xi = 0.69/(B**2*5.039/Gamma/Boltzmann/282.34)**0.5          # Eq 14

            # Eq 19
            F = exp(-18.66*deltaT**2) * exp(-4.25*deltarho**4)

            # Eq 18
            c = (self.M/rho.gcc/Avogadro/Boltzmann/T)**0.5
            d = Boltzmann*T**2/6/pi/fase.mu.muPas/xi
            lc = c*d*fase.dpdT_rho**2*fase.kappa**0.5*F

        return unidades.ThermalConductivity(lo+l2+lc, "mWmK")

    _thermal = thermo0, thermo1


class Test(TestCase):

    def test_smukala(self):
        # Zero enthalpy-entropy reference state
        st = Ethylene(T=298.15, P=101325)
        self.assertEqual(round(st.h.kJkg, 2), 0)
        self.assertEqual(round(st.s.kJkgK, 3), 0)

        # Selected values from Table 32, Pag 1093, saturation state
        # Using custom parametr for reference state, the enthalpy and entropy
        # values are diferent to table
        st = Ethylene(T=103.989, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000122)
        self.assertEqual(round(st.Liquido.rho, 2), 654.60)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.6220)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.4294)
        self.assertEqual(round(st.Liquido.w, 1), 1766.6)
        self.assertEqual(round(st.Gas.rho, 5), 0.00396)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.89014)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.1868)
        self.assertEqual(round(st.Gas.w, 2), 202.67)

        st = Ethylene(T=150, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.027377)
        self.assertEqual(round(st.Liquido.rho, 2), 594.60)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.4275)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.4038)
        self.assertEqual(round(st.Liquido.w, 1), 1449.4)
        self.assertEqual(round(st.Gas.rho, 5), 0.62385)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.91795)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.2320)
        self.assertEqual(round(st.Gas.w, 2), 241.10)

        st = Ethylene(T=200, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.45548)
        self.assertEqual(round(st.Liquido.rho, 2), 521.22)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.3214)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 2.5287)
        self.assertEqual(round(st.Liquido.w, 1), 1069.9)
        self.assertEqual(round(st.Gas.rho, 4), 8.4936)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.0431)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.4920)
        self.assertEqual(round(st.Gas.w, 2), 261.94)

        st = Ethylene(T=250, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 2.3295)
        self.assertEqual(round(st.Liquido.rho, 2), 422.02)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.3680)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.3629)
        self.assertEqual(round(st.Liquido.w, 1), 628.10)
        self.assertEqual(round(st.Gas.rho, 3), 44.970)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 1.3344)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.6609)
        self.assertEqual(round(st.Gas.w, 2), 248.80)

        st = Ethylene(T=282, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 5.0022)
        self.assertEqual(round(st.Liquido.rho, 2), 253.12)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.2089)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 146.97)
        self.assertEqual(round(st.Liquido.w, 2), 188.89)
        self.assertEqual(round(st.Gas.rho, 2), 175.80)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 2.3981)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 225.24)
        self.assertEqual(round(st.Gas.w, 2), 191.32)

        # Selected values from Table 33, Pag 1097, single phase region
        st = Ethylene(T=104.003, P=1e5)
        self.assertEqual(round(st.rho, 2), 654.63)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6219)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4293)
        self.assertEqual(round(st.w, 1), 1767.0)

        st = Ethylene(T=205, P=5e5)
        self.assertEqual(round(st.rho, 4), 9.1140)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0511)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.5024)
        self.assertEqual(round(st.w, 2), 264.57)

        st = Ethylene(T=450, P=1e6)
        self.assertEqual(round(st.rho, 4), 7.6018)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.7717)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.0933)
        self.assertEqual(round(st.w, 2), 391.57)

        st = Ethylene(T=180, P=2e6)
        self.assertEqual(round(st.rho, 2), 554.35)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.3484)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4263)
        self.assertEqual(round(st.w, 1), 1244.3)

        st = Ethylene(T=104.690, P=5e6)
        self.assertEqual(round(st.rho, 2), 656.09)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6202)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4226)
        self.assertEqual(round(st.w, 1), 1789.1)

        st = Ethylene(T=255, P=1e7)
        self.assertEqual(round(st.rho, 2), 444.77)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.3569)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.8116)
        self.assertEqual(round(st.w, 2), 773.66)

        st = Ethylene(T=450, P=2e7)
        self.assertEqual(round(st.rho, 2), 176.61)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.8596)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.7794)
        self.assertEqual(round(st.w, 2), 421.74)

        st = Ethylene(T=340, P=5e7)
        self.assertEqual(round(st.rho, 2), 430.74)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.5848)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4990)
        self.assertEqual(round(st.w, 2), 889.85)

        st = Ethylene(T=127.136, P=1e8)
        self.assertEqual(round(st.rho, 2), 670.20)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6349)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3202)
        self.assertEqual(round(st.w, 1), 2003.4)

        st = Ethylene(T=320, P=2e8)
        self.assertEqual(round(st.rho, 2), 576.41)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6686)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.2550)
        self.assertEqual(round(st.w, 1), 1668.7)

    def test_mccarty(self):
        # Selected values from Table B-1, Pag 74, saturation state
        st = Ethylene(T=104, x=0.5, eq="mccarty")
        self.assertEqual(round(st.P.MPa, 5), 0.00012)
        self.assertEqual(round(st.Liquido.rhoM, 3), 23.389)
        self.assertEqual(round(st.Liquido.dpdrho_T.MPakgm3*st.M, 3), 46.027)
        self.assertEqual(round(st.Liquido.dpdT_rho.MPaK, 3), 2.907)
        self.assertEqual(round(st.Liquido.uM.Jmol, 2), 6615.72)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 6615.73)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 84.4)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 39.49)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 74.39)
        self.assertEqual(round(st.Liquido.w, 0), 1758)
        self.assertEqual(round(st.Gas.rhoM, 6), 0.000141)
        self.assertEqual(round(st.Gas.dpdrho_T.MPakgm3*st.M, 3), 0.864)
        self.assertEqual(round(st.Gas.dpdT_rho.MPaK, 3), 0.000)
        self.assertEqual(round(st.Gas.uM.Jmol, 2), 21685.99)
        self.assertEqual(round(st.Gas.hM.Jmol, 2), 22550.47)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 237.6)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 24.99)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 33.31)
        self.assertEqual(round(st.Gas.w, 0), 203)

        st = Ethylene(T=125, x=0.5, eq="mccarty")
        self.assertEqual(round(st.P.MPa, 5), 0.00252)
        self.assertEqual(round(st.Liquido.rhoM, 3), 22.345)
        self.assertEqual(round(st.Liquido.dpdrho_T.MPakgm3*st.M, 3), 47.453)
        self.assertEqual(round(st.Liquido.dpdT_rho.MPaK, 3), 2.136)
        self.assertEqual(round(st.Liquido.uM.Jmol, 2), 8067.58)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 8067.69)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 97.1)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 43.19)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 67.25)
        self.assertEqual(round(st.Liquido.w, 0), 1623)
        self.assertEqual(round(st.Gas.rhoM, 6), 0.002432)
        self.assertEqual(round(st.Gas.dpdrho_T.MPakgm3*st.M, 3), 1.034)
        self.assertEqual(round(st.Gas.dpdT_rho.MPaK, 3), 0.000)
        self.assertEqual(round(st.Gas.uM.Jmol, 2), 22204.76)
        self.assertEqual(round(st.Gas.hM.Jmol, 2), 23241.57)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 218.5)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 25.19)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 33.61)
        self.assertEqual(round(st.Gas.w, 0), 222)

        st = Ethylene(T=150, x=0.5, eq="mccarty")
        self.assertEqual(round(st.P.MPa, 5), 0.02735)
        self.assertEqual(round(st.Liquido.rhoM, 3), 21.202)
        self.assertEqual(round(st.Liquido.dpdrho_T.MPakgm3*st.M, 3), 35.511)
        self.assertEqual(round(st.Liquido.dpdT_rho.MPaK, 3), 1.687)
        self.assertEqual(round(st.Liquido.uM.Jmol, 2), 9738.07)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 9739.36)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 109.3)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 40.14)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 66.88)
        self.assertEqual(round(st.Liquido.w, 0), 1452)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0222)
        self.assertEqual(round(st.Gas.dpdrho_T.MPakgm3*st.M, 3), 1.213)
        self.assertEqual(round(st.Gas.dpdT_rho.MPaK, 3), 0.000)
        self.assertEqual(round(st.Gas.uM.Jmol, 2), 22798.71)
        self.assertEqual(round(st.Gas.hM.Jmol, 2), 24029.05)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 204.6)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 25.97)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 34.87)
        self.assertEqual(round(st.Gas.w, 0), 241)

        st = Ethylene(T=175, x=0.5, eq="mccarty")
        self.assertEqual(round(st.P.MPa, 5), 0.13956)
        self.assertEqual(round(st.Liquido.rhoM, 3), 19.955)
        self.assertEqual(round(st.Liquido.dpdrho_T.MPakgm3*st.M, 3), 24.970)
        self.assertEqual(round(st.Liquido.dpdT_rho.MPaK, 3), 1.317)
        self.assertEqual(round(st.Liquido.uM.Jmol, 2), 11423.46)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 11430.45)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 119.7)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 37.83)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 68.34)
        self.assertEqual(round(st.Liquido.w, 0), 1268)
        self.assertEqual(round(st.Gas.rhoM, 6), 0.100137)
        self.assertEqual(round(st.Gas.dpdrho_T.MPakgm3*st.M, 3), 1.332)
        self.assertEqual(round(st.Gas.dpdT_rho.MPaK, 3), 0.001)
        self.assertEqual(round(st.Gas.uM.Jmol, 2), 23341.06)
        self.assertEqual(round(st.Gas.hM.Jmol, 2), 24734.75)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 195.7)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 27.56)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 37.69)
        self.assertEqual(round(st.Gas.w, 0), 255)

        st = Ethylene(T=225, x=0.5, eq="mccarty")
        self.assertEqual(round(st.P.MPa, 4), 1.1279)
        self.assertEqual(round(st.Liquido.rhoM, 3), 17.006)
        self.assertEqual(round(st.Liquido.dpdrho_T.MPakgm3*st.M, 3), 10.002)
        self.assertEqual(round(st.Liquido.dpdT_rho.MPaK, 3), 0.721)
        self.assertEqual(round(st.Liquido.uM.Jmol, 2), 14966.73)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 15033.06)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 137.5)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 37.13)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 77.58)
        self.assertEqual(round(st.Liquido.w, 0), 863)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.73325)
        self.assertEqual(round(st.Gas.dpdrho_T.MPakgm3*st.M, 3), 1.208)
        self.assertEqual(round(st.Gas.dpdT_rho.MPaK, 3), 0.007)
        self.assertEqual(round(st.Gas.uM.Jmol, 2), 24130.08)
        self.assertEqual(round(st.Gas.hM.Jmol, 2), 25668.30)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 184.8)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 33.21)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 52.17)
        self.assertEqual(round(st.Gas.w, 0), 260)

        # Table B-2, pag 81, single phase region
        st = Ethylene(T=150, P=5e4, eq="mccarty")
        self.assertEqual(round(st.rhoM, 3), 21.203)
        self.assertEqual(round(st.dpdrho_T.MPakgm3*st.M, 3), 35.52)
        self.assertEqual(round(st.dpdT_rho.MPaK, 3), 1.687)
        self.assertEqual(round(st.uM.Jmol, 2), 9737.71)
        self.assertEqual(round(st.hM.Jmol, 2), 9740.07)
        self.assertEqual(round(st.sM.JmolK, 1), 109.3)
        self.assertEqual(round(st.cvM.JmolK, 2), 40.14)
        self.assertEqual(round(st.cpM.JmolK, 2), 66.87)
        self.assertEqual(round(st.w, 0), 1452)

        st = Ethylene(T=300, P=101325, eq="mccarty")
        self.assertEqual(round(st.rhoM, 6), 0.040847)
        self.assertEqual(round(st.dpdrho_T.MPakgm3*st.M, 3), 2.467)
        self.assertEqual(round(st.dpdT_rho.MPaK, 6), 0.000342)
        self.assertEqual(round(st.uM.Jmol, 2), 27165.87)
        self.assertEqual(round(st.hM.Jmol, 2), 29646.46)
        self.assertEqual(round(st.sM.JmolK, 1), 219.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.81)
        self.assertEqual(round(st.cpM.JmolK, 2), 43.32)
        self.assertEqual(round(st.w, 0), 331)

        st = Ethylene(T=450, P=1e6, eq="mccarty")
        self.assertEqual(round(st.rhoM, 5), 0.27072)
        self.assertEqual(round(st.dpdrho_T.MPakgm3*st.M, 3), 3.646)
        self.assertEqual(round(st.dpdT_rho.MPaK, 5), 0.00232)
        self.assertEqual(round(st.uM.Jmol, 2), 33363.38)
        self.assertEqual(round(st.hM.Jmol, 2), 37057.24)
        self.assertEqual(round(st.sM.JmolK, 1), 220.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 49.68)
        self.assertEqual(round(st.cpM.JmolK, 2), 58.74)
        self.assertEqual(round(st.w, 0), 392)

        st = Ethylene(T=110, P=1e7, eq="mccarty")
        self.assertEqual(round(st.rhoM, 3), 23.252)
        self.assertEqual(round(st.dpdrho_T.MPakgm3*st.M, 3), 49.352)
        self.assertEqual(round(st.dpdT_rho.MPaK, 2), 2.63)
        self.assertEqual(round(st.uM.Jmol, 2), 6939.63)
        self.assertEqual(round(st.hM.Jmol, 2), 7369.71)
        self.assertEqual(round(st.sM.JmolK, 1), 87.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 43.54)
        self.assertEqual(round(st.cpM.JmolK, 2), 72.05)
        self.assertEqual(round(st.w, 0), 1706)

        st = Ethylene(T=450, P=4e7, eq="mccarty")
        self.assertEqual(round(st.rhoM, 3), 10.569)
        self.assertEqual(round(st.dpdrho_T.MPakgm3*st.M, 3), 6.793)
        self.assertEqual(round(st.dpdT_rho.MPaK, 3), 0.202)
        self.assertEqual(round(st.uM.Jmol, 2), 28192.71)
        self.assertEqual(round(st.hM.Jmol, 2), 31977.53)
        self.assertEqual(round(st.sM.JmolK, 1), 180.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 53.35)
        self.assertEqual(round(st.cpM.JmolK, 2), 77.55)
        self.assertEqual(round(st.w, 0), 593)

    def test_jahangiri(self):
        # Table 24, Pag 635, Second virial coefficients
        st = Ethylene(T=200.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -310.248)
        st = Ethylene(T=250.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -199.921)
        st = Ethylene(T=300.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -138.087)
        st = Ethylene(T=350.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -98.356)
        st = Ethylene(T=400.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -70.599)
        st = Ethylene(T=450.15, P=101325, eq="jahangiri")
        self.assertEqual(round(st.virialB.ccg*st.M, 3), -50.099)

        # # Table 27, Pag 646, saturation states
        st = Ethylene(T=104, x=0.5, eq="jahangiri")
        self.assertEqual(round(st.P.MPa, 5), 0.00012)
        self.assertEqual(round(st.Liquido.rhoM, 3), 23.347)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), 6613.4)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 84.32)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 39.00)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 61.67)
        self.assertEqual(round(st.Liquido.w, 0), 1822)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.00014)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 22552)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 237.57)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 24.97)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 33.29)
        self.assertEqual(round(st.Gas.w, 0), 203)

        st = Ethylene(T=150, x=0.5, eq="jahangiri")
        self.assertEqual(round(st.P.MPa, 5), 0.02747)
        self.assertEqual(round(st.Liquido.rhoM, 3), 21.206)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), 9780.0)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 109.54)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 39.29)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 66.68)
        self.assertEqual(round(st.Liquido.w, 0), 1457)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.02231)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 24039)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 204.60)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 25.76)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 34.56)
        self.assertEqual(round(st.Gas.w, 0), 241)

        st = Ethylene(T=200, x=0.5, eq="jahangiri")
        self.assertEqual(round(st.P.MPa, 5), 0.45596)
        self.assertEqual(round(st.Liquido.rhoM, 3), 18.589)
        self.assertEqual(round(st.Liquido.hM.Jmol, 0), 13190)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 129.01)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 37.13)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 71.01)
        self.assertEqual(round(st.Liquido.w, 0), 1070)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.30297)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 25325)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 189.68)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 29.46)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 42.12)
        self.assertEqual(round(st.Gas.w, 0), 262)

        st = Ethylene(T=250, x=0.5, eq="jahangiri")
        self.assertEqual(round(st.P.MPa, 4), 2.3307)
        self.assertEqual(round(st.Liquido.rhoM, 3), 15.048)
        self.assertEqual(round(st.Liquido.hM.Jmol, 0), 17143)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 146.08)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 38.57)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 94.60)
        self.assertEqual(round(st.Liquido.w, 0), 627)
        self.assertEqual(round(st.Gas.rhoM, 4), 1.6042)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 25672)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 180.19)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 37.48)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 74.66)
        self.assertEqual(round(st.Gas.w, 0), 249)

        # Table 28, pag 656, single phase values
        st = Ethylene(T=105, P=1e4, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 23.302)
        self.assertEqual(round(st.uM.Jmol, 1), 6676.1)
        self.assertEqual(round(st.hM.Jmol, 1), 6676.5)
        self.assertEqual(round(st.sM.JmolK, 2), 84.92)
        self.assertEqual(round(st.cvM.JmolK, 2), 40.86)
        self.assertEqual(round(st.cpM.JmolK, 2), 63.76)
        self.assertEqual(round(st.w, 0), 1801)

        st = Ethylene(T=160, P=5e4, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.03831)
        self.assertEqual(round(st.uM.Jmol, 0), 23043)
        self.assertEqual(round(st.hM.Jmol, 0), 24348)
        self.assertEqual(round(st.sM.JmolK, 2), 201.69)
        self.assertEqual(round(st.cvM.JmolK, 2), 26.13)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.17)
        self.assertEqual(round(st.w, 0), 248)

        st = Ethylene(T=450, P=1e5, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.02676)
        self.assertEqual(round(st.uM.Jmol, 0), 33502)
        self.assertEqual(round(st.hM.Jmol, 0), 37238)
        self.assertEqual(round(st.sM.JmolK, 2), 239.82)
        self.assertEqual(round(st.cvM.JmolK, 2), 49.54)
        self.assertEqual(round(st.cpM.JmolK, 2), 57.93)
        self.assertEqual(round(st.w, 0), 394)

        st = Ethylene(T=300, P=101325, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.04085)
        self.assertEqual(round(st.uM.Jmol, 0), 27166)
        self.assertEqual(round(st.hM.Jmol, 0), 29646)
        self.assertEqual(round(st.sM.JmolK, 2), 219.39)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.77)
        self.assertEqual(round(st.cpM.JmolK, 2), 43.29)
        self.assertEqual(round(st.w, 0), 331)

        st = Ethylene(T=180, P=2e5, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 19.694)
        self.assertEqual(round(st.uM.Jmol, 0), 11786)
        self.assertEqual(round(st.hM.Jmol, 0), 11796)
        self.assertEqual(round(st.sM.JmolK, 2), 121.74)
        self.assertEqual(round(st.cvM.JmolK, 2), 37.35)
        self.assertEqual(round(st.cpM.JmolK, 2), 68.10)
        self.assertEqual(round(st.w, 0), 1231)

        st = Ethylene(T=195, P=3e5, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.19795)
        self.assertEqual(round(st.uM.Jmol, 0), 23794)
        self.assertEqual(round(st.hM.Jmol, 0), 25309)
        self.assertEqual(round(st.sM.JmolK, 2), 192.81)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.41)
        self.assertEqual(round(st.cpM.JmolK, 2), 39.43)
        self.assertEqual(round(st.w, 0), 264)

        st = Ethylene(T=160, P=5e5, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 20.731)
        self.assertEqual(round(st.uM.Jmol, 0), 10436)
        self.assertEqual(round(st.hM.Jmol, 0), 10460)
        self.assertEqual(round(st.sM.JmolK, 2), 113.78)
        self.assertEqual(round(st.cvM.JmolK, 2), 38.28)
        self.assertEqual(round(st.cpM.JmolK, 2), 66.60)
        self.assertEqual(round(st.w, 0), 1387)

        st = Ethylene(T=288, P=8e5, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.35246)
        self.assertEqual(round(st.uM.Jmol, 0), 26525)
        self.assertEqual(round(st.hM.Jmol, 0), 28794)
        self.assertEqual(round(st.sM.JmolK, 2), 199.66)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.21)
        self.assertEqual(round(st.cpM.JmolK, 2), 44.56)
        self.assertEqual(round(st.w, 0), 316)

        st = Ethylene(T=220, P=1e6, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 17.360)
        self.assertEqual(round(st.uM.Jmol, 0), 14606)
        self.assertEqual(round(st.hM.Jmol, 0), 14663)
        self.assertEqual(round(st.sM.JmolK, 2), 135.88)
        self.assertEqual(round(st.cvM.JmolK, 2), 37.32)
        self.assertEqual(round(st.cpM.JmolK, 2), 76.00)
        self.assertEqual(round(st.w, 0), 902)

        st = Ethylene(T=244, P=2e6, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 15.582)
        self.assertEqual(round(st.uM.Jmol, 0), 16479)
        self.assertEqual(round(st.hM.Jmol, 0), 16607)
        self.assertEqual(round(st.sM.JmolK, 2), 143.99)
        self.assertEqual(round(st.cvM.JmolK, 2), 38.18)
        self.assertEqual(round(st.cpM.JmolK, 2), 88.69)
        self.assertEqual(round(st.w, 0), 685)

        st = Ethylene(T=450, P=3e6, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 5), 0.83490)
        self.assertEqual(round(st.uM.Jmol, 0), 33027)
        self.assertEqual(round(st.hM.Jmol, 0), 36621)
        self.assertEqual(round(st.sM.JmolK, 2), 210.49)
        self.assertEqual(round(st.cvM.JmolK, 2), 49.92)
        self.assertEqual(round(st.cpM.JmolK, 2), 60.60)
        self.assertEqual(round(st.w, 0), 387)

        st = Ethylene(T=280, P=5e6, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 11.098)
        self.assertEqual(round(st.uM.Jmol, 0), 20068)
        self.assertEqual(round(st.hM.Jmol, 0), 20519)
        self.assertEqual(round(st.sM.JmolK, 2), 158.00)
        self.assertEqual(round(st.cvM.JmolK, 2), 43.83)
        self.assertEqual(round(st.cpM.JmolK, 1), 254.2)
        self.assertEqual(round(st.w, 0), 311)

        st = Ethylene(T=385, P=8e6, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 4), 3.1028)
        self.assertEqual(round(st.uM.Jmol, 0), 28666)
        self.assertEqual(round(st.hM.Jmol, 0), 31244)
        self.assertEqual(round(st.sM.JmolK, 2), 190.15)
        self.assertEqual(round(st.cvM.JmolK, 2), 45.30)
        self.assertEqual(round(st.cpM.JmolK, 2), 67.89)
        self.assertEqual(round(st.w, 0), 336)

        st = Ethylene(T=110, P=1e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 23.250)
        self.assertEqual(round(st.uM.Jmol, 1), 6922.6)
        self.assertEqual(round(st.hM.Jmol, 1), 7352.8)
        self.assertEqual(round(st.sM.JmolK, 2), 87.21)
        self.assertEqual(round(st.cvM.JmolK, 2), 46.18)
        self.assertEqual(round(st.cpM.JmolK, 2), 69.34)
        self.assertEqual(round(st.w, 0), 1779)

        st = Ethylene(T=244, P=2e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 17.396)
        self.assertEqual(round(st.uM.Jmol, 0), 15489)
        self.assertEqual(round(st.hM.Jmol, 0), 16638)
        self.assertEqual(round(st.sM.JmolK, 2), 139.68)
        self.assertEqual(round(st.cvM.JmolK, 2), 38.06)
        self.assertEqual(round(st.cpM.JmolK, 2), 69.84)
        self.assertEqual(round(st.w, 0), 992)

        st = Ethylene(T=390, P=3e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 11.067)
        self.assertEqual(round(st.uM.Jmol, 0), 24895)
        self.assertEqual(round(st.hM.Jmol, 0), 27606)
        self.assertEqual(round(st.sM.JmolK, 2), 172.52)
        self.assertEqual(round(st.cvM.JmolK, 2), 47.06)
        self.assertEqual(round(st.cpM.JmolK, 2), 79.47)
        self.assertEqual(round(st.w, 0), 565)

        st = Ethylene(T=130, P=4e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 22.926)
        self.assertEqual(round(st.uM.Jmol, 1), 8020.7)
        self.assertEqual(round(st.hM.Jmol, 1), 9765.5)
        self.assertEqual(round(st.sM.JmolK, 2), 96.52)
        self.assertEqual(round(st.cvM.JmolK, 2), 44.99)
        self.assertEqual(round(st.cpM.JmolK, 2), 67.44)
        self.assertEqual(round(st.w, 0), 1779)

        st = Ethylene(T=350, P=5e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 15.003)
        self.assertEqual(round(st.uM.Jmol, 0), 21092)
        self.assertEqual(round(st.hM.Jmol, 0), 24425)
        self.assertEqual(round(st.sM.JmolK, 2), 159.84)
        self.assertEqual(round(st.cvM.JmolK, 2), 44.46)
        self.assertEqual(round(st.cpM.JmolK, 2), 70.13)
        self.assertEqual(round(st.w, 0), 867)

        st = Ethylene(T=240, P=8e7, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 19.935)
        self.assertEqual(round(st.uM.Jmol, 0), 13969)
        self.assertEqual(round(st.hM.Jmol, 0), 17982)
        self.assertEqual(round(st.sM.JmolK, 2), 132.07)
        self.assertEqual(round(st.cvM.JmolK, 2), 40.20)
        self.assertEqual(round(st.cpM.JmolK, 2), 62.12)
        self.assertEqual(round(st.w, 0), 1434)

        st = Ethylene(T=350, P=1e8, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 17.490)
        self.assertEqual(round(st.uM.Jmol, 0), 19860)
        self.assertEqual(round(st.hM.Jmol, 0), 25578)
        self.assertEqual(round(st.sM.JmolK, 2), 154.41)
        self.assertEqual(round(st.cvM.JmolK, 2), 46.02)
        self.assertEqual(round(st.cpM.JmolK, 2), 65.80)
        self.assertEqual(round(st.w, 0), 1203)

        st = Ethylene(T=350, P=2.6e8, eq="jahangiri")
        self.assertEqual(round(st.rhoM, 3), 20.953)
        self.assertEqual(round(st.uM.Jmol, 0), 18436)
        self.assertEqual(round(st.hM.Jmol, 0), 30845)
        self.assertEqual(round(st.sM.JmolK, 2), 145.92)
        self.assertEqual(round(st.cvM.JmolK, 2), 48.98)
        self.assertEqual(round(st.cpM.JmolK, 2), 64.39)
        self.assertEqual(round(st.w, 0), 1792)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = Ethylene(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 2.7682)
        self.assertEqual(round(st.P.MPa, 3), 48.416)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.065)

        st2 = Ethylene(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 174.10)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.47681)

    def test_Assael(self):
        # Table 5, pag 8
        self.assertEqual(round(Ethylene(T=200, rho=0).k.mWmK, 2), 10.39)
        self.assertEqual(round(Ethylene(T=300, rho=0).k.mWmK, 2), 21.01)
        self.assertEqual(round(Ethylene(T=400, rho=0).k.mWmK, 2), 36.36)
        self.assertEqual(round(Ethylene(T=500, rho=0).k.mWmK, 2), 55.05)
        self.assertEqual(round(Ethylene(T=200, P=1e5).k.mWmK, 2), 10.54)
        self.assertEqual(round(Ethylene(T=300, P=1e5).k.mWmK, 2), 21.09)
        self.assertEqual(round(Ethylene(T=400, P=1e5).k.mWmK, 2), 36.40)
        self.assertEqual(round(Ethylene(T=500, P=1e5).k.mWmK, 2), 55.07)
        self.assertEqual(round(Ethylene(T=200, P=5e7).k.mWmK, 1), 190.4)
        self.assertEqual(round(Ethylene(T=300, P=5e7).k.mWmK, 1), 126.9)
        self.assertEqual(round(Ethylene(T=400, P=5e7).k.mWmK, 2), 98.08)
        self.assertEqual(round(Ethylene(T=500, P=5e7).k.mWmK, 2), 94.57)
        self.assertEqual(round(Ethylene(T=200, P=1e8).k.mWmK, 1), 223.5)
        self.assertEqual(round(Ethylene(T=300, P=1e8).k.mWmK, 1), 164.3)
        self.assertEqual(round(Ethylene(T=400, P=1e8).k.mWmK, 1), 132.0)
        self.assertEqual(round(Ethylene(T=500, P=1e8).k.mWmK, 1), 121.3)
        self.assertEqual(round(Ethylene(T=200, P=1.5e8).k.mWmK, 1), 252.9)
        self.assertEqual(round(Ethylene(T=300, P=1.5e8).k.mWmK, 1), 196.4)
        self.assertEqual(round(Ethylene(T=400, P=1.5e8).k.mWmK, 1), 161.7)
        self.assertEqual(round(Ethylene(T=500, P=1.5e8).k.mWmK, 1), 145.8)
        self.assertEqual(round(Ethylene(T=200, P=2e8).k.mWmK, 1), 280.5)
        self.assertEqual(round(Ethylene(T=300, P=2e8).k.mWmK, 1), 226.5)
        self.assertEqual(round(Ethylene(T=400, P=2e8).k.mWmK, 1), 190.3)
        self.assertEqual(round(Ethylene(T=500, P=2e8).k.mWmK, 1), 170.6)

        # Critical enhancement point, section 3.1.4, pag 8
        self.assertEqual(round(Ethylene(T=300, rho=300).k.mWmK, 2), 69.62)

    def test_Holland(self):
        # Single phase selected point
        # Viscosity, Table 5, pag 924
        # Thermal Conductivity, Table 6, pag 927
        st = Ethylene(T=110, P=1e5, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 5660.5)
        self.assertEqual(round(st.k.mWmK, 2), 261.77)

        st = Ethylene(T=140, P=1e6, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 2769.8)
        self.assertEqual(round(st.k.mWmK, 2), 223.14)

        st = Ethylene(T=200, P=5e6, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 1223.7)
        self.assertEqual(round(st.k.mWmK, 1), 158.5)

        st = Ethylene(T=300, P=1e5, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 103.8)
        self.assertEqual(round(st.k.mWmK, 2), 20.58)

        st = Ethylene(T=130, P=1e7, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 3278.5)
        self.assertEqual(round(st.k.mWmK, 2), 244.97)

        st = Ethylene(T=300, P=5e7, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 759.0)
        self.assertEqual(round(st.k.mWmK, 2), 129.40)

        st = Ethylene(T=500, P=1e5, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 165.1)
        self.assertEqual(round(st.k.mWmK, 2), 49.95)

        st = Ethylene(T=500, P=5e7, eq="mccarty", thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 394.1)
        self.assertEqual(round(st.k.mWmK, 2), 93.57)

        # Saturated liquid values, Table 7, pag 930
        # The density values differ so the calculated transport properties
        # difffer in that cases, I think the paper use the ancillary equation
        # for liquid saturated density
        st = Ethylene(T=200, x=0, eq="mccarty", thermal=1)
        self.assertEqual(round(st.rhoM, 3), 18.575)
        self.assertEqual(round(st.mu.microP, 1), 1204.9)
        self.assertEqual(round(st.k.mWmK, 1), 153.5)

        # Dilute gas values, Table 8, pag 930
        st = Ethylene(T=180, rho=0, thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 63.6)
        self.assertEqual(round(st.k.mWmK, 1), 10.0)

        st = Ethylene(T=250, rho=0, thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 86.6)
        self.assertEqual(round(st.k.mWmK, 1), 14.9)

        st = Ethylene(T=400, rho=0, thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 135.7)
        self.assertEqual(round(st.k.mWmK, 1), 34.6)

        st = Ethylene(T=500, rho=0, thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 164.8)
        self.assertEqual(round(st.k.mWmK, 1), 49.9)

        st = Ethylene(T=680, rho=0, thermal=1)
        self.assertEqual(round(st.mu.microP, 1), 211.6)
        self.assertEqual(round(st.k.mWmK, 1), 91.2)
