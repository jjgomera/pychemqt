#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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

from lib import unidades
from lib.meos import MEoS


class nPropanol(MEoS):
    """Multiparameter equation of state for 1-propanol"""
    name = "1-propanol"
    CASNumber = "71-23-8"
    formula = "CH3-CH2-CH0H"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(270.90835016)
    Tc = unidades.Temperature(536.85)
    Pc = unidades.Pressure(5180.1, "kPa")
    M = 60.09502  # g/mol
    Tt = unidades.Temperature(148.764)
    Tb = unidades.Temperature(370.267)
    f_acent = 0.619
    momentoDipolar = unidades.DipoleMoment(1.68, "Debye")
    id = 146

    CP1 = {"ao": 4.8,
           "ao_pow": [5.8589, 15.295, 6.7218],
           "ao_exp": [583, 1851, 4963]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Gao (2023)",
        "__doi__": {
            "autor": "Cui, J.W., Gao, K.H., Wu, J.T.",
            "title": "Reference Correlation of the Viscosity of n-Propyl"
                     " Alcohol from 153 to 618 K and up to 118 MPa",
            "ref": "J. Phys. Chem. Ref. Data 54(3) (2025) 033106",
            "doi": "10.1063/5.0280205"},

        # Real reference
        # Gao, K.H., Wu, J.T., Lemmon, E.W.
        # Reference equation of state including new association terms form
        # thermodynamic properties of strong hydrogen bonding fluids:
        # n-Propyl alcohol
        # Coefficient from supplementary material of Cui viscosity correlation

        "R": 8.3144598,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650, "Pmax": 1200000, "rhomax": 20,

        "nr1": [0.8501919587338e-2, 0.1769572654173e1, -0.2601616361951e1,
                -0.8069586507251, 0.2699091898091, -0.1517099751574],
        "d1": [4, 1, 1, 2, 3, 2],
        "t1": [1, 0.383486245565, 1.00673506221, 1.000655973878,
               0.737864791347, 0.75644462067],

        "nr2": [-0.136826694139, 0.8186955171291, -0.1203832429174,
                -0.4787596045266, -0.8352292782826, -0.1183324061791e-1,
                0.1434056590817, -0.1207715298215e1, 0.1691887723885e-1],
        "d2": [3, 2, 3, 1, 1, 1, 2, 2, 2],
        "t2": [1.467225414428, 2.436829620603, 0.473745795953, 0.465572868745,
               1.890981559655, 12.649968615492, 1.264813959298, 4.289888248699,
               4.52793432034],
        "c2": [2, 2, 1, 1, 2, 3, 1, 2, 3],
        "gamma2": [1]*9,

        "nr3": [0.10839278117, 0.1146462676421e1, -0.1450150007612,
                0.6173688940994, -0.2525141509726, 0.2968112512824e-1,
                -0.9084481771871e1, -0.1141221110371e1],
        "d3": [1, 1, 3, 1, 2, 2, 1, 1],
        "t3": [14.73606008844, 0.559526227486, 1.617872490534, 2.115863103468,
               2.011307384187, 4.628909223213, 2.469668299598, 17.167683459526],
        "alfa3": [1499.87889513, 0.05730523, 0.0609921, 0.42166279, 0.33164463,
                  0.25315696, 0.97960817, 228.01781612],
        "beta3": [20.2675115, 0.6422969, 0.5898219, 0.1843193, 0.2186612,
                  3.2564005, 5.0181177, 20.0187472],
        "gamma3": [1.841707031, 1.061285623, 0.596951091, 0.833305507,
                   1.503334094, 0.565023908, 1.214102874, 1.809889063],
        "epsilon3": [0.015068023, 2.018661345, -0.589364213, 1.9423568,
                     1.085816589, 1.510740203, -0.122527829, 0.135063877],

        "nr_ass": [-0.5110634296449e-1, -0.1938787172895e-2,
                   -0.1865834030357e-1, 0.1194615747104e-2],
        "d_ass": [1, 1, 1, 1],
        "t_ass": [0.606487738911, 2.538481255171, 1.636947334141,
                  3.175214539812],
        "alfa_ass": [0.72884509, 1.09047435, 0.20356499, 0.76583431],
        "beta_ass": [2.4878069, 0.0890810, 0.0697761, 0.0583633],
        "gamma_ass": [1.127828698, 1.135904338, 1.688273219, 1.116235043],
        "epsilon_ass": [0.474063025, -0.67936823, -0.206270268, -0.454639836],
        "b_ass": [0.553635447, 0.11706642, 0.305768855, 0.105098566]}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-propanol of Sun and Ely"
                    " (2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650, "Pmax": 1200000, "rhomax": 20,

        "nr1": [-6.4846669, 6.3481226e-1, 5.34271316, 3.59156552e-2,
                3.91173758e-4, -4.4277807e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.33146361, 1.71475104, -1.20634979e-2, 2.02582101e-1,
                4.4959531e-2, -8.06185866e-1, -1.97404896e-2, 4.98309152e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = (gao, sun)

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Sanjuán, E.L.",
                    "title": "Surface Tension of Alcohols. Data Selection and "
                             "Recommended Correlations",
                    "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033104",
                    "doi": "10.1063/1.4927858"},
        "sigma": [-10.2192, 0.177644, 10.2509],
        "exp": [1.46876, 7.99999, 1.46348]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.33463, 0.765236, -17.9082, 21.7367, -10.1526],
        "t": [1, 1.5, 3.3, 4.4, 6]}
    _liquid_Density = {
        "eq": 1,
        "n": [4.6250, -8.1472, 16.761, -21.674, 11.590],
        "t": [0.46, 0.92, 1.4, 2.0, 2.5]}

    _vapor_Density = {
        "eq": 2,
        "n": [-3.42592, -8.524, -37.0732, -57.7727, -142.86, -2171.92],
        "t": [0.414, 1.25, 3.3, 7.1, 14, 38]}

    visco0 = {
        "__name__": "Cui (2025)",
        "__doi__": {
            "autor": "Cui, J.W., Gao, K.H., Wu, J.T.",
            "title": "Reference Correlation of the Viscosity of n-Propyl"
                     " Alcohol from 153 to 618 K and up to 118 MPa",
            "ref": "J. Phys. Chem. Ref. Data 54(3) (2025) 033106",
            "doi": "10.1063/5.0280205"},

        "eq": 1, "omega": 0,

        "Toref": 1000,
        "no_exp": [27.899, -23.467, -1.216],
        "to_exp": [0.54998, 0.62021, 0],
        "logo_exp": [0, 0, 0],

        "special": "_mur"}

    def _mur(self, rho, T, fase):
        """Special term of residual viscosity for Cui correlation"""
        Tr = T/self.Tc
        rhor = rho/270.91

        # Eq 9
        f1 = exp(-(2*rho/50)**3*(2*rho/50-1)**3)
        mugh = (-0.2715*rhor**3+14.966*rhor**2+5.6515*rhor-1.3907)*(1-f1) \
            - 0.1*f1

        # Eq 14
        c = [-1.3734e1, 1.7997, 4.1456e-1, -1.228e-1, 1.1279e-2]
        mu01 = 1e6*exp(sum(ci*(1000/T)**i for i, ci in enumerate(c)))

        # Eq 15
        b = [9.40527, -6.625e-1, 6.05537, 2.01883e-4, 9.80519e-1, -1.50583e1,
             -7.98126e-1, 7.74207e-1, 7.17048e-1]
        Dmu = rhor**(2/3) * Tr**0.5 * mu01 * (
            b[0] + b[1]/Tr + b[2]*Tr**2 + b[3]*rhor**8 + b[4]*Tr*rhor**2
            + b[5]*Tr + b[6]*rhor**2
            + b[7]*rhor**1.2/(Tr**-13.4+b[8]*rhor**3.45))

        kw = {"eq": 1, "omega": 0,
              "ek": 481.5, "sigma": 0.4738, "M": 60.095,

              "Toref": 1000,
              "no_exp": [27.899, -23.467, -1.216],
              "to_exp": [0.54998, 0.62021, 0],
              "logo_exp": [0, 0, 0],

              "Tref_virial": 481.5,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5]}
        mu1 = self._Viscosity(rho, T, None, kw, True)

        # Eq 10
        muh = mu1 + Dmu
        # Eq 16
        f2 = exp(-(2*rhor)**3*(2*rhor-1)**3)
        mur = f2*mugh + (1-f2)*muh

        return mur

    _viscosity = (visco0, )


class Test(TestCase):
    """Testing"""

    def test_cui(self):
        """Table 9, pag 12"""
        # The table include too density checking of Gao EoS
        st = nPropanol(T=160, P=0)
        self.assertEqual(round(st.rho, 5), 0)
        self.assertEqual(round(st.mu.muPas, 7), 4.1039044)

        st = nPropanol(T=250, P=0)
        self.assertEqual(round(st.rho, 5), 0)
        self.assertEqual(round(st.mu.muPas, 7), 6.3708690)

        st = nPropanol(T=450, P=0)
        self.assertEqual(round(st.rho, 5), 0)
        self.assertEqual(round(st.mu.muPas, 6), 11.671487)

        st = nPropanol(T=620, P=0)
        self.assertEqual(round(st.rho, 5), 0)
        self.assertEqual(round(st.mu.muPas, 6), 16.132564)

        st = nPropanol(T=160, P=1e6)
        self.assertEqual(round(st.rho, 5), 911.50799)
        self.assertEqual(round(st.mu.muPas, 1), 2570089.2)

        st = nPropanol(T=250, P=1e6)
        self.assertEqual(round(st.rho, 5), 838.00112)
        self.assertEqual(round(st.mu.muPas, 4), 7720.2906)

        st = nPropanol(T=450, P=1e6)
        self.assertEqual(round(st.rho, 6), 19.272004)
        self.assertEqual(round(st.mu.muPas, 6), 11.665344)

        st = nPropanol(T=620, P=1e6)
        self.assertEqual(round(st.rho, 6), 12.125910)
        self.assertEqual(round(st.mu.muPas, 6), 16.149306)

        st = nPropanol(T=160, P=1e7)
        self.assertEqual(round(st.rho, 5), 915.29782)
        self.assertEqual(round(st.mu.muPas, 1), 2778735.9)

        st = nPropanol(T=250, P=1e7)
        self.assertEqual(round(st.rho, 5), 843.48990)
        self.assertEqual(round(st.mu.muPas, 4), 8285.0350)

        st = nPropanol(T=450, P=1e7)
        self.assertEqual(round(st.rho, 5), 657.77167)
        self.assertEqual(round(st.mu.muPas, 5), 183.50757)

        st = nPropanol(T=620, P=1e7)
        self.assertEqual(round(st.rho, 5), 201.10060)
        # self.assertEqual(round(st.mu.muPas, 6), 34.668332)
        self.assertEqual(round(st.mu.muPas, 6), 34.668362)

        st = nPropanol(T=160, P=1.2e8)
        self.assertEqual(round(st.rho, 5), 953.69232)
        self.assertEqual(round(st.mu.muPas, 1), 5564764.3)

        st = nPropanol(T=250, P=1.2e8)
        self.assertEqual(round(st.rho, 5), 893.54372)
        self.assertEqual(round(st.mu.muPas, 3), 16004.243)

        st = nPropanol(T=450, P=1.2e8)
        self.assertEqual(round(st.rho, 5), 773.36561)
        self.assertEqual(round(st.mu.muPas, 5), 400.61496)

        st = nPropanol(T=620, P=1.2e8)
        self.assertEqual(round(st.rho, 5), 663.20347)
        self.assertEqual(round(st.mu.muPas, 5), 158.25547)
