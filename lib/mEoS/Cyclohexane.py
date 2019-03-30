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

from lib import unidades
from lib.meos import MEoS


class Cyclohexane(MEoS):
    """Multiparameter equation of state for cyclohexane"""
    name = "cyclohexane"
    CASNumber = "110-82-7"
    formula = "cyclo(CH2)6"
    synonym = ""
    _refPropName = "CYCLOHEX"
    _coolPropName = "CycloHexane"
    rhoc = unidades.Density(271.33016352)
    Tc = unidades.Temperature(553.6)
    Pc = unidades.Pressure(4080.5, "kPa")
    M = 84.15948  # g/mol
    Tt = unidades.Temperature(279.47)
    Tb = unidades.Temperature(353.865)
    f_acent = 0.2096
    momentoDipolar = unidades.DipoleMoment(0.3, "Debye")
    id = 38
    _Tr = unidades.Temperature(526.231121)
    _rhor = unidades.Density(274.647526)
    _w = 0.221837522

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [0.9891140602, 1.6359660572],
           "ao_exp": [0.83775, 16.036, 24.636, 7.1715],
           "titao": [773/Tc, 941/Tc, 2185/Tc, 4495/Tc],
           "ao_hyp": [], "hyp": []}

    CP1 = {"ao": 9.368327211,
           "an": [-56214088, 0.01526155409, -3.6352468e-6],
           "pow": [-3, 1, 2],
           "ao_exp": [23.766589], "exp": [2000],
           "ao_hyp": [], "hyp": []}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclohexane of Zhou et "
                    "al. (2014)",
        "__doi__": {"autor": "Zhou, Y., Liu, J., Penoncello, S.G., Lemmon, "
                             "E.W.",
                    "title": "An Equation of State for the Thermodynamic "
                             "Properties of Cyclohexane",
                    "ref": "J. Phys. Chem. Ref. Data 43 (2014) 043105",
                    "doi": "10.1063/1.4900538"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": 279.86, "Tmax": 700.0, "Pmax": 250000.0, "rhomax": 10.3,
        "Pmin": 5.2402, "rhomin": 9.403,

        "nr1": [0.05483581, 1.607734, -2.375928, -0.5137709, 0.1858417],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.37, 0.79, 1.075, 0.37],

        "nr2": [-0.9007515, -0.5628776, 0.2903717, -0.3279141, -0.03177644],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.4, 2.5, 0.5, 3, 1.06],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.8668676, -0.1962725, -0.1425992, 0.004197016, 0.1776584,
                -0.04433903, -0.03861246, 0.07399692, 0.02036006, 0.00272825],
        "d3": [1, 1, 3, 3, 2, 2, 3, 2, 3, 2],
        "t3": [1.6, 0.37, 1.33, 2.5, 0.9, 0.5, 0.73, 0.2, 1.5, 1.5],
        "alfa3": [0.99, 1.43, 0.97, 1.93, 0.92, 1.27, 0.87, 0.82, 1.4, 3],
        "beta3": [0.38, 4.2, 1.2, 0.9, 1.2, 2.6, 5.3, 4.4, 4.2, 25],
        "gamma3": [0.65, 0.63, 1.14, 0.09, 0.56, 0.4, 1.01, 0.45, 0.85, 0.86],
        "epsilon3": [0.73, 0.75, 0.48, 2.32, 0.2, 1.33, 0.68, 1.11, 1.47, .99]}

    penoncello = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclohexane of "
                    "Penoncello et al. (1995)",
        "__doi__": {"autor": "Penoncello, S.G., Jacobsen, R.T., Goodwin, "
                             "A.R.H.",
                    "title": "A Thermodynamic Property Formulation for "
                             "Cyclohexane",
                    "ref": "Int. J. Thermophys., 16(2) (1995) 519-531",
                    "doi": "10.1007/BF01441918"},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 279.47, "Pref": 101.325, "ho": 33884.8, "so": 96.612},
        "Tt": 279.47, "Tc": 553.64, "rhoc": 3.244, "M": 84.1608,

        "Tmin": 279.47, "Tmax": 700.0, "Pmax": 80000.0, "rhomax": 9.77,
        "Pmin": 5.2538, "rhomin": 9.4045,

        "nr1": [0.8425412659, -0.3138388327e1, 0.1679072631e1, -0.153819249,
                0.1984911143, -0.144532594, 0.3746346428e-3, 0.1861479616e-3,
                0.1745721652e-3],
        "d1": [1, 1, 1, 2, 3, 3, 7, 6, 6],
        "t1": [0, 1.5, 2.5, 1.5, 1, 2.5, 2, 0.5, 3],

        "nr2": [-0.6427428062, 0.2280757615, -0.1868116802e1, -0.1028243711e1,
                0.5821457418, -0.255891152, 0.1276844113e-1, -0.5158613166e-2,
                0.6334794755e-1, -0.6014686589e-1, 0.4439056828, -0.6264920642,
                0.2132589969e1, -0.3620300991e-2, 0.2534453992,
                0.1669144715e-1, 0.3985052291e-2],
        "d2": [1, 1, 2, 3, 3, 5, 8, 10, 3, 4, 1, 1, 2, 2, 4, 4, 8],
        "t2": [5, 6, 5.5, 3, 7, 6, 6.5, 5.5, 11, 11, 0.5, 1, 4, 4, 1.5, 2, .5],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 3, 2, 6, 2, 4, 2],
        "gamma2": [1]*17}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for cyclohexane of "
                    "Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": {"Tref": 279.47, "Pref": 101.325, "ho": 33884.8, "so": 96.612},
        "Tc": 553.60, "rhoc": 273.02/84.161, "M": 84.161,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 9.77,
        "Pmin": 5.2428, "rhomin": 9.3999,

        "nr1": [0.10232354e1, -0.29204964e1, 0.10736630e1, -0.19573985,
                0.12228111, 0.28943321e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.27231767, -0.4483332e-1, -0.38253334, -0.89835333e-1,
                -0.24874965e-1, 0.10836132e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclohexane of Sun and "
                    "Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 279.47, "Pref": 101.325, "ho": 33884.8, "so": 96.612},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.27436292, 1.15372124, -3.86726473, 8.84627298e-2,
                2.76478090e-4, 7.26682313e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [7.10849914e-2, 4.46376742e-1, 7.64476190e-1, -4.23520282e-2,
                -0.396468623, -1.41250071e-2, -1.08371284e-1, -2.50082884e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = zhou, penoncello, shortSpan, sun

    _surface = {"sigma": [0.06485], "exp": [1.263]}
    _melting = {"eq": 1, "Tref": 1, "Pref": 700,
                "Tmin": Tt, "Tmax": 370.0,
                "a1": [0.1329969885, -374.255624], "exp1": [1.41, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.0342, 1.7311, -1.7572, -3.3406],
        "t": [1., 1.5, 2.3, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [5.5081, -14.486, 38.241, -64.589, 57.919, -20.55],
        "t": [0.51, 0.94, 1.4, 1.9, 2.4, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.69006, -41.4239, 220.914, -443.72, 491.49, -296.373],
        "t": [0.446, 1.98, 2.75, 3.3, 4.1, 4.8]}

    visco0 = {"__name__": "Tariq (2014)",
              "__doi__": {
                  "autor": "Tariq, U., Jusoh, A.R.B., Riesco, N., Vesovic, V.",
                  "title": "Reference Correlation of the Viscosity of "
                           "Cyclohexane from the Triple Point to 700K and up "
                           "to 110 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 43(3) (2014) 033101",
                  "doi": "10.1063/1.4891103"},

              "eq": 1, "omega": 3,
              "collision": [-1.5093, 364.87, -39537],

              "sigma": 1,
              "n_chapman": 0.19592/M**0.5,

              "Tref_res": 553.6, "rhoref_res": 3.224*M,
              "nr": [335.234, 7.8494803, -687.3976, 362.0868, -10.4793856,
                     2.5521774, 17.2734993, -5.9372242, -10.6186149, 4.3982781,
                     2.8894928, -1.3468174, -0.2938491, 0.1487134],
              "dr": [2.2, 2.5, 2.5, 2.8, 10, 10, 11, 11, 12, 12, 13, 13, 14,
                     14],
              "tr": [1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
              "gr": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              "cr": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              "special": "_vir"}

    def _vir(self, rho, T, fase):
        # The initial density dependence has a different expresion, without muo
        # and other normal method calculation so hardcoded here
        muB = 0
        if rho:
            for i, n in enumerate([5.09643, -3387.21, 337477]):
                muB += n/T**i
        return muB*rho/self.M

    _viscosity = visco0,


class Test(TestCase):

    def test_zhou(self):
        # Table 5, Pag 17

        # Possible Erratum reference Table 5, in enthalpy and entropy values,
        # only work the two phases region, and the reference state, but not
        # the other point.
        # My values are coherent with the values returned by CoolProp

        st = Cyclohexane(T=300, rhom=9.4)
        self.assertEqual(round(st.P.MPa, 6), 24.173705)
        self.assertEqual(round(st.cvM.JmolK, 5), 115.28600)
        self.assertEqual(round(st.cpM.JmolK, 5), 154.76956)
        self.assertEqual(round(st.w, 4), 1383.3878)
        # self.assertEqual(round(st.hM.Jmol, 4), -8400.0834)
        # self.assertEqual(round(st.sM.JmolK, 6), -28.889069)

        st = Cyclohexane(T=500, rhom=6.5)
        self.assertEqual(round(st.P.MPa, 7), 3.9246630)
        self.assertEqual(round(st.cvM.JmolK, 5), 192.52056)
        self.assertEqual(round(st.cpM.JmolK, 5), 255.57087)
        self.assertEqual(round(st.w, 5), 434.13064)
        # self.assertEqual(round(st.hM.Jmol, 3), 31070.127)
        # self.assertEqual(round(st.sM.JmolK, 6), 70.891447)

        st = Cyclohexane(T=500, rhom=0.7)
        self.assertEqual(round(st.P.MPa, 7), 1.9981172)
        self.assertEqual(round(st.cvM.JmolK, 5), 191.96446)
        self.assertEqual(round(st.cpM.JmolK, 5), 235.52281)
        self.assertEqual(round(st.w, 5), 155.34800)
        # self.assertEqual(round(st.hM.Jmol, 3), 52757.706)
        # self.assertEqual(round(st.sM.JmolK, 5), 122.92657)

        st = Cyclohexane(T=600, rhom=3.5)
        self.assertEqual(round(st.P.MPa, 7), 6.8225506)
        self.assertEqual(round(st.cvM.JmolK, 5), 232.79222)
        self.assertEqual(round(st.cpM.JmolK, 5), 388.55185)
        self.assertEqual(round(st.w, 5), 150.53318)
        # self.assertEqual(round(st.hM.Jmol, 3), 70150.132)
        # self.assertEqual(round(st.sM.JmolK, 5), 143.42323)

        st = Cyclohexane(T=553.6, rhom=3.3)
        self.assertEqual(round(st.P.MPa, 7), 4.0805433)
        self.assertEqual(round(st.cvM.JmolK, 5), 224.19555)
        self.assertEqual(round(st.cpM.JmolK, 2), 199224.62)
        self.assertEqual(round(st.w, 6), 87.913911)
        # self.assertEqual(round(st.hM.Jmol, 3), 58532.604)
        # self.assertEqual(round(st.sM.JmolK, 5), 123.59810)

        st = Cyclohexane(P=101325, x=0.5)
        self.assertEqual(round(st.T, 6), 353.864939)
        self.assertEqual(round(st.Liquido.rhoM, 7), 8.5487851)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 5), 134.61630)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 5), 179.07223)
        self.assertEqual(round(st.Liquido.w, 5), 994.05862)
        self.assertEqual(round(st.Liquido.hM.Jmol, 3), 0)
        self.assertEqual(round(st.Liquido.sM.JmolK, 6), 0)
        self.assertEqual(round(st.Gas.rhoM, 9), 0.035779032)
        self.assertEqual(round(st.Gas.cvM.JmolK, 5), 123.43050)
        self.assertEqual(round(st.Gas.cpM.JmolK, 5), 133.35895)
        self.assertEqual(round(st.Gas.w, 5), 186.91349)
        self.assertEqual(round(st.Gas.hM.Jmol, 3), 29991.286)
        self.assertEqual(round(st.Gas.sM.JmolK, 6), 84.753484)

        # Reference state
        st = Cyclohexane(T=300, P=1000)
        self.assertEqual(round(st.sM.JmolK, 4), 104.2926)
        self.assertEqual(round(st.hM.Jmol, 2), 23949.02)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = Cyclohexane(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.0278)
        self.assertEqual(round(st.P.MPa, 3), 9.007)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.5927)

        st2 = Cyclohexane(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 206.82)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.31449)

    def test_tariq(self):
        # Table 8, pag 11
        st = Cyclohexane(T=300, rhom=0)
        self.assertEqual(round(st.mu.muPas, 3), 7.058)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = Cyclohexane(T=300, rhom=0.0430)
        mu = st._Viscosity(0.0430*st.M, 300, None)
        self.assertEqual(round(mu.muPas, 3), 6.977)

        st = Cyclohexane(T=300, rhom=9.1756)
        self.assertEqual(round(st.mu.muPas, 2), 863.66)
        st = Cyclohexane(T=300, rhom=9.9508)
        self.assertEqual(round(st.mu.muPas, 2), 2850.18)
        st = Cyclohexane(T=500, rhom=0)
        self.assertEqual(round(st.mu.muPas, 3), 11.189)

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = Cyclohexane(T=500, rhom=6.0213)
        mu = st._Viscosity(6.0213*st.M, 500, None)
        self.assertEqual(round(mu.muPas, 3), 94.842)

        st = Cyclohexane(T=500, rhom=8.5915)
        self.assertEqual(round(st.mu.muPas, 2), 380.04)
        st = Cyclohexane(T=700, rhom=0)
        self.assertEqual(round(st.mu.muPas, 3), 15.093)
        st = Cyclohexane(T=700, rhom=7.4765)
        self.assertEqual(round(st.mu.muPas, 3), 176.749)
