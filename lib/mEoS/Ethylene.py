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


from unittest import TestCase

from scipy import exp

from lib.meos import MEoS
from lib import unidades


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
           "an": [], "pow": [],
           "ao_exp": [1]*12,
           "exp": [4353.907145, 2335.2251475, 1930.913215, 1471.9256475,
                   4464.6972475, 1778.39697, 1365.4520425, 1356.8190475,
                   4469.013745, 1188.475645, 4300.6703425, 2077.67413],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.554495281,
           "an": [0.5603615762e6, -0.2141069802e5, 0.2532008897e3,
                  -0.9951927478e-2, 0.5108931070e-4, -0.1928667482e-7],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-0.2061703241e2],
           "exp": [3000],
           "ao_hyp": [], "hyp": []}

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
        "Pmin": 0.12265, "rhomin": 23.334,

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
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 29610, "so": 219.225},
        "M": 28.054, "Tc": 282.3452, "rhoc": 7.634,

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 260000.0, "rhomax": 26.67,
        "Pmin": 0.1225, "rhomin": 23.348,

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
        "Pmin": 0.12123, "rhomin": 23.34,
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

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethylene of McCarty and Jacobsen (1981)",
        "__doi__": {"autor": "McCarty, R.D., Jacobsen, R.T.",
                    "title": "An Equation of State for Fluid Ethylene",
                    "ref": "Natl. Bur. Stand., Tech. Note 1045, 1981.",
                    "doi": ""},
        "__test__":
            # Table, Pag 138
            """
            >>> st=Ethylene(T=110, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 23.2516 110 7369.71 87.438 43.54 72.05 1706.12
            >>> st=Ethylene(T=150, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 21.4687 150 10056.28 108.312 40.64 65.72 1507.50
            >>> st=Ethylene(T=200, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 19.0802 200 13387.64 127.463 37.33 68.11 1172.01
            >>> st=Ethylene(T=250, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 16.1946 250 16970.99 143.418 38.08 77.31 811.65
            >>> st=Ethylene(T=300, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 11.5343 300 21586.21 160.141 42.53 118.76 418.90
            >>> st=Ethylene(T=350, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 5.3640 350 27794.36 179.364 44.3 92.15 309.36
            >>> st=Ethylene(T=400, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 3.7626 400 31676.15 189.761 46.49 70.37 347.16
            """,
        "R": 8.31434,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 40000.0, "rhomax": 23.343,
        "Pmin": 0.1213, "rhomin": 23.343,

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
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 29610, "so": 219.225},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

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

    # eq = smukala, MBWR, jahangiri, shortSpan, sun
    eq = smukala, jahangiri, shortSpan, sun

    _surface = {"sigma": [0.0477], "exp": [1.17]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [10.725], "expt1": [0], "expd1": [1],
                   "a2": [55.19, 49.5, -2045, -1154.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.9, 2.9]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "Tmin": Tt, "Tmax": 450.0,
                "a1": [0.1225e-3, 0.357924e3, -0.357924e3],
                "exp1": [0, 0.20645e1, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
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

    visco1 = {"__name__": "NIST",
              "__doi__": {
                  "autor": "",
                  "title": "Coefficients are taken from NIST14, Version 9.08",
                  "ref": "",
                  "doi": ""},

              "eq": 2, "omega": 2,

              "ek": 224.7, "sigma": 0.4163,
              "n_chapman": 0.141374566253583/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-8.03553028329404, -439.8962514, 8.69536237617,
                    5773.08496161, .267589139152, -34.39391627, 66.4795135739],
              "rhoc": 7.63299886259}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Holland (1983)",
               "__doi__": {
                   "autor": "Holland, P.M., Eaton, B.E., Hanley, H.J.M.",
                   "title": "A Correlation of the Viscosity and Thermal "
                            "Conductivity Data of Gaseous and Liquid Ethylene",
                   "ref": "J. Phys. Chem. Ref. Data 12(4) (1983) 917-932",
                   "doi": "10.1063/1.555701"},

               "__test__":
                   # Table 6, pag 927
                   """
                   >>> st=Ethylene(T=110, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   261.77
                   >>> st=Ethylene(T=140, P=1e6)
                   >>> print "%0.2f" % st.k.mWmK
                   223.14
                   >>> st=Ethylene(T=200, P=5e6)
                   >>> print "%0.2f" % st.k.mWmK
                   158.50
                   >>> st=Ethylene(T=300, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   20.56
                   >>> st=Ethylene(T=130, P=1e7)
                   >>> print "%0.2f" % st.k.mWmK
                   244.97
                   >>> st=Ethylene(T=300, P=5e7)
                   >>> print "%0.2f" % st.k.mWmK
                   129.32
                   >>> st=Ethylene(T=500, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   49.95
                   >>> st=Ethylene(T=500, P=5e7)
                   >>> print "%0.2f" % st.k.mWmK
                   93.57
                   """}

    def _thermo0(self, rho, T, fase):
        GT = [-2.903423528e5, 4.680624952e5, -1.8954783215e5, -4.8262235392e3,
              2.243409372e4, -6.6206354818e3, 8.9937717078e2, -6.0559143718e1,
              1.6370306422]
        lo = 0
        for i in range(-3, 6):
            lo += GT[i+3]*T**(i/3.)

        tita = (rho.gcc-0.221)/0.221
        j = [0, -1.304503323e1, 1.8214616599e1, 9.903022496e3, 7.420521631e2,
             -3.0083271933e-1, 9.6456068829e1, 1.350256962e4]
        l1 = exp(j[1]+j[4]/T)*(exp(rho.gcc**0.1*(j[2]+j[3]/self.T**1.5)+tita*rho.gcc**0.5*(j[5]+j[6]/T+j[7]/T**2))-1.)

        lc = 0
        # FIXME: no sale
#        deltarho=(self.rho/self.M-0.221)/0.221
#        deltaT=(self.T-282.34)/282.34
#        xkt=(1.0/self.rho/self.M/self.derivative("P", "rho", "T")*1e3)**0.5
#        b=abs(deltarho)/abs(deltaT)**1.19
#        xts=(self.rho/self.M)**2*xkt*5.039/.221**2
#        g=xts*abs(deltaT)**1.19
#        xi=0.69/(b**2*5.039/g/Boltzmann/282.34)**0.5
#        f=exp(-18.66*deltaT**2-4.25*deltarho**4)
#        c=(self.M/self.rho.gcc/Avogadro/Boltzmann/self.T)**0.5
#        lc=c*Boltzmann*self.T**2/6.0/pi/self.mu.muPas/xi*self.dpdT**2*self.kappa**0.5*f
#        print lo, l1
        return unidades.ThermalConductivity(lo+l1+lc, "mWmK")


    thermo1 = {"eq": 1, "critical": 0,
               "__name__": "NIST14",
               "__doi__": {"autor": "",
                           "title": "Coefficients are taken from NIST14, Version 9.08",
                           "ref": "",
                           "doi": ""},

               "Tref": 224.7, "kref": 1e-3,
               "no": [1.35558587, -0.14207565869509, 1],
               "co": [0, -1, -96],

               "Trefb": 282.350007277, "rhorefb": 7.63299886259, "krefb": 1e-3,
               "nb": [15.3064493136, 25.0280721432, -15.4526955192,
                      0.8590418672, 3.32700049633, -0.333048907849],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6}

    _thermal = thermo0, thermo1

# TODO: Add MBWR equation of Younglove


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

    def test_jahangiri(self):
        # FIXME: The ideal properties from Cp expresion don't work
        # Work from _Cp0, but don't from phio
        pass

        # Selected values from Table 18, Pag 613, ideal properties
        # st = Ethylene(T=105, P=1e5, eq="jahangiri")
        # self.assertEqual(round(st.cpM0.JmolK, 3), 33.278)
        # self.assertEqual(round(st.hM0.Jmol, 3), -7023.816)
        # self.assertEqual(round(st.sM0.JmolK, 3), -37.155)

        # st = Ethylene(T=200, P=1e5, eq="jahangiri")
        # self.assertEqual(round(st.cpM0.JmolK, 3), 35.352)
        # self.assertEqual(round(st.hM0.Jmol, 3), -3801.885)
        # self.assertEqual(round(st.sM0.JmolK, 3), -15.367)

        # st = Ethylene(T=300, P=1e5, eq="jahangiri")
        # self.assertEqual(round(st.cpM0.JmolK, 3), 43.027)
        # self.assertEqual(round(st.hM0.Jmol, 3), 79.437)
        # self.assertEqual(round(st.sM0.JmolK, 3), 0.266)

        # st = Ethylene(T=400, P=1e5, eq="jahangiri")
        # self.assertEqual(round(st.cpM0.JmolK, 3), 52.989)
        # self.assertEqual(round(st.hM0.Jmol, 3), 4877.189)
        # self.assertEqual(round(st.sM0.JmolK, 3), 13.999)

        # st = Ethylene(T=500, P=1e5, eq="jahangiri")
        # self.assertEqual(round(st.cpM0.JmolK, 3), 62.411)
        # self.assertEqual(round(st.hM0.Jmol, 3), 10656.721)
        # self.assertEqual(round(st.sM0.JmolK, 3), 26.856)

        # # Table 24, Pag 635, Second virial coefficients
        # st = Ethylene(T=200.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -310.248)
        # st = Ethylene(T=250.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -199.921)
        # st = Ethylene(T=300.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -138.087)
        # st = Ethylene(T=350.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -98.356)
        # st = Ethylene(T=400.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -70.599)
        # st = Ethylene(T=450.15, x=0, eq="jahangiri")
        # self.assertEqual(round(st.virialB.ccg*st.M, 3), -50.099)

        # # Table 25, Pag 637, saturation states
        # st = Ethylene(T=238.18, x=0.5, eq="jahangiri")
        # self.assertEqual(round(st.P.MPa, 3), 1.681)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 16108.4)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 1), 25739.7)

        # st = Ethylene(T=258.15, x=0.5, eq="jahangiri")
        # self.assertEqual(round(st.P.MPa, 3), 2.869)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 17915.0)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 1), 25528.3)

        # st = Ethylene(T=278.15, x=0.5, eq="jahangiri")
        # self.assertEqual(round(st.P.MPa, 3), 4.589)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 20437.1)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 1), 24294.9)

        # st = Ethylene(T=282.15, x=0.5, eq="jahangiri")
        # self.assertEqual(round(st.P.MPa, 3), 5.018)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 21644.4)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 1), 23021.6)

        # st = Ethylene(T=282.25, x=0.5, eq="jahangiri")
        # self.assertEqual(round(st.P.MPa, 3), 5.029)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 21781.5)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 1), 22853.6)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = Ethylene(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 2.7682)
        self.assertEqual(round(st.P.MPa, 3), 48.416)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.065)

        st2 = Ethylene(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 174.10)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.47681)

    def test_Holland(self):
        # Table 5, pag 924

        # FIXME: The state point use MBWR mEoS
        return

        # st = Ethylene(T=110, P=1e5, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 5660.5)
        # st = Ethylene(T=140, P=1e6, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 2769.8)
        # st = Ethylene(T=200, P=5e6, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 1223.7)
        # st = Ethylene(T=300, P=1e5, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 103.8)
        # st = Ethylene(T=130, P=1e7, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 3278.5)
        # st = Ethylene(T=300, P=5e7, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 759.0)
        # st = Ethylene(T=500, P=1e5, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 165.1)
        # st = Ethylene(T=500, P=5e7, eq="MBWR")
        # self.assertEqual(round(st.mu.microP, 1), 394.1)
