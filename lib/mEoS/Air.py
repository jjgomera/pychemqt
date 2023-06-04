#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from math import log
from unittest import TestCase

from lib import unidades
from lib.meos import MEoSBlend


class Air(MEoSBlend):
    """Multiparameter equation of state for Air as pseudocomponent"""
    name = "air"
    CASNumber = "1"
    formula = "N2+Ar+O2"
    synonym = "R-729"
    _refPropName = "AIR"
    _coolPropName = "Air"
    rhoc = unidades.Density(342.60456)
    Tc = unidades.Temperature(132.6306)
    Pc = unidades.Pressure(3786.0, "kPa")
    M = 28.9586  # g/mol
    Tt = unidades.Temperature(59.75)
    Tb = unidades.Temperature(78.903)
    f_acent = 0.0335
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 475

    Fi1 = {"ao_log": [1, 2.490888032],
           "pow": [-3, -2, -1, 0, 1, 1.5],
           "ao_pow": [0.6057194e-7, -0.210274769e-4, -0.158860716e-3,
                      -13.841928076, 17.275266575, -0.19536342e-3],
           "ao_exp": [0.791309509, 0.212236768],
           "titao": [25.36365, 16.90741],
           "ao_exp2": [-0.197938904],
           "titao2": [87.31279],
           "sum2": [2./3]
           }

    Fi2 = {"ao_log": [1, 2.494156],
           "pow": [-3, -2, -1, 0, 1, 2, 3, 1.5],
           "ao_pow": [0.3893909e-6, -0.4437673e-4, 0.4532773e-3, -0.1382389e2,
                      -0.01204447, -0.8418768e-3, 0.4671333e-4, -0.1953634e-3],
           "tau*logtau": 0.3544246e-2,
           "ao_exp": [0.2122368, 0.7872444],
           "titao": [0.1690741e2, 0.2528369e2]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Lemmon et al. "
                    "(2000)",
        "__doi__": {
            "autor": "Lemmon, E.W., Jacobsen, R.T, Penoncello, S.G., Friend, "
                     "D.G.",
            "title": "Thermodynamic Properties of Air and Mixtures of "
                     "Nitrogen, Argon, and Oxygen From 60 to 2000 K at "
                     "Pressures to 2000 MPa",
            "ref": "J. Phys. Chem. Ref. Data 29, 331 (2000)",
            "doi": "10.1063/1.1285884"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 101.325, "ho": 0, "so": 0},

        "M": 28.9586, "Tc": 132.6312, "rhoc": 10.4477,

        "Tmin": Tt, "Tmax": 2000., "Pmax": 2000000.0, "rhomax": 53.73,

        "Tj": 132.6312, "Pj": 3.78502,
        "dew": {"i": [1, 2, 5, 8],
                "n": [-0.1567266, -5.539635, 0.7567212, -3.514322]},
        "bubble": {"i": [1, 2, 3, 4, 5, 6],
                   "n": [0.2260724, -7.080499, 5.700283, -12.44017, 17.81926,
                         -10.81364]},

        "nr1": [0.118160747229, 0.713116392079, -0.161824192067e1,
                0.714140178971e-1, -0.865421396646e-1, 0.134211176704,
                0.112626704218e-1, -0.420533228842e-1, 0.349008431982e-1,
                0.164957183186e-3],
        "d1": [1, 1, 1, 2, 3, 3, 4, 4, 4, 6],
        "t1": [0, 0.33, 1.01, 0, 0, 0.15, 0, 0.2, 0.35, 1.35],

        "nr2": [-0.101365037912, -0.173813690970, -0.472103183731e-1,
                -0.122523554253e-1, -0.146629609713, -0.316055879821e-1,
                0.233594806142e-3, 0.148287891978e-1, -0.938782884667e-2],
        "d2": [1, 3, 5, 6, 1, 3, 11, 1, 3],
        "t2": [1.6, 0.8, 0.95, 1.25, 3.6, 6, 3.25, 3.5, 15],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3, 3],
        "gamma2": [1]*9}

    jacobsen = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Jacobsen et al. "
                    "(1992)",
        "__doi__": {
            "autor": "Jacobsen, R.T, Penoncello, S.G., Beyerlein, S.W., "
                     "Clarke, W.P., Lemmon, E.W.",
            "title": "A Thermodynamic Property Formulation for Air",
            "ref": "Fluid Phase Equilibria, 79:113-124, 1992.",
            "doi": "10.1016/0378-3812(92)85124-Q"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 8649.34, "so": 194.},
        "M": 28.9586, "Tc": 132.6312, "rhoc": 10.4477,

        "Tmin": Tt, "Tmax": 870.0, "Pmax": 70000.0, "rhomax": 34.628,

        "Tj": 132.61738, "Pj": 3.78502,
        "dew": {"i": [1, 2, 10, 11, 13, 14],
                "n": [-0.1537763029, -5.544542064, 312.7182733, -895.9553274,
                      1834.176566, -1321.892808]},
        "bubble": {"i": [1, 2, 4, 5, 6, 7, 12],
                   "n": [0.2095592444, -6.654905539, 22.13718815, -84.14553609,
                         135.9753732, -83.66895082, 17.97856602]},

        "nr1": [0.206604930965, 0.367099749382, -0.943192015369,
                0.382519513142e-2, -0.865385542309e-1, 0.323019987452,
                0.608695449299e-2, 0.128352106296e-3, -0.400058181940e-5],
        "d1": [1, 1, 1, 1, 2, 2, 4, 6, 7],
        "t1": [0, 0.25, 1, 3.5, 0, 0.25, 0.5, 2, 3],

        "nr2": [-0.544697915817, -0.526471065792, -0.608529300347,
                -0.124174277875, -0.595578533411e-2, -0.157523548353,
                -0.346463251040e-2, 0.837023084176e-2, -0.316701981142e-1,
                -0.721856676857e-2, 0.276838040645e-3, 0.160877459321e-4,
                0.409235806738e-1, 0.652776125216e-3, -0.952903961290e-2,
                -0.100337820004e-1, 0.701111041628e-2, -0.472754336912e-2,
                0.399257638569e-2, 0.968453675994e-2, -0.106826283630e-1,
                -0.489679885832e-2],
        "d2": [1, 2, 3, 5, 6, 1, 1, 2, 2, 3, 11, 11, 1, 1, 2, 3, 7, 8, 2, 4, 5,
               2],
        "t2": [1.5, 1, 1, 1, 2, 3, 8, 0.5, 5.5, 9, 3, 6, 3, 9, 2, 13, 11, 11,
               8, 22, 23, 11],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4,
               5],
        "gamma2": [1]*22}

    panasiti = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Panasiti et al. "
                    "(1999)",
        "__doi__": {
            "autor": "Panasiti, M.D., Lemmon, E.W., Penoncello, S.G., "
                     "Jacobsen, R.T, Friend, D.G.",
            "title": "Thermodynamic Properties of Air from 60 to 2000 K at "
                     "Pressures up to 2000 MPa",
            "ref": "Int. J. Themophys. 20(1) (1999) 217-228",
            "doi": "10.1023/a:1021450818807"},

        "R": 8.31451,
        "cp": Fi2,
        "ref": {"Tref": 273.15, "Pref": 101.325, "ho": 0, "so": 0},

        "M": 28.9585, "Tc": 132.6312, "rhoc": 10.4477,

        "Tmin": Tt, "Tmax": 2000., "Pmax": 2000000.0, "rhomax": 53.73,

        "Tj": 132.6312, "Pj": 3.78502,
        "dew": {"i": [1, 2, 5, 8],
                "n": [-0.1567266, -5.539635, 0.7567212, -3.514322]},
        "bubble": {"i": [1, 2, 3, 4, 5, 6],
                   "n": [0.2260724, -7.080499, 5.700283, -12.44017, 17.81926,
                         -10.81364]},

        "nr1": [0.653047013811e-1, 0.678835009195, -0.143490786106e1,
                -0.321903633516e-2, 0.583911230624e-1, 0.239601980951e-1,
                -0.584000786794e-1, 0.137011320694, 0.219932875035e-3,
                0.837104476841e-2, -0.306746271489e-1, 0.338386156676e-4,
                0.164630578408e-5],
        "d1": [1, 1, 1, 1, 2, 4, 3, 3, 3, 4, 4, 6, 8],
        "t1": [0, 0.25, 1, 3.5, 0, 0.5, 0, 0.25, 3.5, 0, 0.25, 0.5, 2],

        "nr2": [-0.176893815222, -0.353933053336, -0.566621330528e-1,
                -0.190100444223, -0.215976092083e-2, -0.308948799371e-2,
                0.293942842065e-3, 0.293832797468e-1, -0.321576125494e-1,
                0.554138225907e-2, -0.542472707547e-2, 0.205042936393e-2,
                -0.219851045253e-1, 0.128752208849e-1, -0.219001412602e-1,
                0.390910141344e-5, -0.574114394530e-2, 0.182974806200e-1,
                -0.142919356764e-2, -0.225716528368e-2],
        "d2": [1, 3, 5, 1, 1, 3, 11, 1, 3, 4, 5, 2, 6, 3, 3, 12, 3, 4, 9, 3],
        "t2": [1.5, 1, 1, 3, 8, 9, 3, 3, 13, 22, 23, 11, 1, 0.5, 5.5, 6, 2, 13,
               11, 11],
        "c2": [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 1, 2, 2, 2, 3, 3, 3, 5],
        "gamma2": [1]*20}

    eq = lemmon, jacobsen, panasiti

    _surface = {"sigma": [0.03046], "exp": [1.28]}

    _melting = {
        "eq": 1,
        "__doi__": lemmon["__doi__"],
        "Tmin": 59.75, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 5265,

        "a0": 1,
        "a2": [35493.5], "exp2": [1.78963]}

    # Dew-point density
    _liquid_Density = {
        "eq": 2,
        "n": [-2.0466, -4.752, -13.259, -47.652],
        "t": [0.41, 1, 2.8, 6.5]}

    def _Vapor_Density(self, T):
        """Bubble-point density ancillary equation, Eq 2 in Lemmon reference"""
        if T < self.Tt:
            T = self.Tt
        elif T > 132.6312:
            T = 132.6312

        Tita = 1-T/132.6312

        N = [44.3413, -240.073, 285.139, -88.3366, -0.892181]
        rhor = 1 + N[0]*Tita**0.65 + N[1]*Tita**0.85 + N[2]*Tita**0.95 + \
            N[3]*Tita**1.1 + N[4]*log(T/132.6312)

        return unidades.Density(rhor*10.4477*self.M)

    visco0 = {"__name__": "Lemmon (2004)",
              "__doi__": {
                  "autor": "Lemmon, E.W., Jacobsen, R.T.",
                  "title": "Viscosity and Thermal Conductivity Equations for "
                           "Nitrogen, Oxygen, Argon, and Air",
                  "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                  "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

              "eq": 1,
              "omega": 1,
              "ek": 103.3, "sigma": 0.36,

              "Tref_res": 132.6312, "rhoref_res": 10.4477*M,
              "nr": [10.72, 1.122, 0.002019, -8.876, -0.02916],
              "tr": [.2, .05, 2.4, .6, 3.6],
              "dr": [1, 4, 9, 1, 8],
              "gr": [0, 0, 0, 1, 1],
              "cr": [0, 0, 0, 1, 1]}

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Lemmon (2004)",
               "__doi__": {
                   "autor": "Lemmon, E.W., Jacobsen, R.T.",
                   "title": "Viscosity and Thermal Conductivity Equations for "
                            "Nitrogen, Oxygen, Argon, and Air",
                   "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                   "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

               "eq": 1,
               "Pc": 3.78502e6,

               "Toref": 132.6312, "koref": 1e-3,
               "no_visco": 1.308,
               "no": [1.405, -1.036],
               "to": [1.1, 0.3],

               "Tref_res": 132.6312, "rhoref_res": 10.4477*M, "kref_res": 1e-3,
               "nr": [8.743, 14.76, -16.62, 3.793, -6.142, -0.3778],
               "tr": [0.1, 0, 0.5, 2.7, 0.3, 1.3],
               "dr": [1, 2, 3, 7, 7, 11],
               "cr": [0, 0, 2, 2, 2, 2],
               "gr": [0, 0, 1, 1, 1, 1],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.11e-9, "gam0": 0.055, "qd": 0.31e-9, "Tcref": 265.262}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_lemmon(self):
        """Selected point from Table A1, Pag 363"""
        self.assertEqual(round(Air._bubbleP(59.75).MPa, 6), 0.005265)
        self.assertEqual(round(Air._dewP(59.75).MPa, 5), 0.00243)
        self.assertEqual(round(Air._bubbleP(70).MPa, 5), 0.03191)
        self.assertEqual(round(Air._dewP(70).MPa, 5), 0.01943)
        self.assertEqual(round(Air._bubbleP(80).MPa, 5), 0.11462)
        self.assertEqual(round(Air._dewP(80).MPa, 5), 0.08232)
        self.assertEqual(round(Air._bubbleP(100).MPa, 5), 0.66313)
        self.assertEqual(round(Air._dewP(100).MPa, 5), 0.56742)
        self.assertEqual(round(Air._bubbleP(120).MPa, 5), 2.15573)
        self.assertEqual(round(Air._dewP(120).MPa, 5), 2.00674)
        self.assertEqual(round(Air._bubbleP(130).MPa, 5), 3.42947)
        self.assertEqual(round(Air._dewP(130).MPa, 5), 3.30835)

        # Selected point from Table A2, Pag 366
        st = Air(T=100, P=101325, rho0=1)
        self.assertEqual(round(st.rhoM, 5), 0.12449)
        self.assertEqual(round(st.uM.kJkmol, 1), 2028.2)
        self.assertEqual(round(st.hM.kJkmol, 1), 2842.1)
        self.assertEqual(round(st.sM.kJkmolK, 2), 166.61)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 21.09)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 30.13)
        self.assertEqual(round(st.w, 1), 198.2)

        st = Air(T=500, P=2e5)
        self.assertEqual(round(st.rhoM, 6), 0.048077)
        self.assertEqual(round(st.uM.kJkmol, 0), 10418.0)
        self.assertEqual(round(st.hM.kJkmol, 0), 14578.0)
        self.assertEqual(round(st.sM.kJkmolK, 2), 208.20)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 21.51)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 29.84)
        self.assertEqual(round(st.w, 1), 446.6)

        st = Air(T=130, P=1e6)
        self.assertEqual(round(st.rhoM, 4), 1.0295)
        self.assertEqual(round(st.uM.kJkmol, 1), 2461.1)
        self.assertEqual(round(st.hM.kJkmol, 1), 3432.5)
        self.assertEqual(round(st.sM.kJkmolK, 2), 153.79)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 22.05)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 34.69)
        self.assertEqual(round(st.w, 1), 216.8)

        st = Air(T=2000, P=10e6)
        self.assertEqual(round(st.rhoM, 5), 0.59094)
        self.assertEqual(round(st.uM.kJkmol, 0), 48600)
        self.assertEqual(round(st.hM.kJkmol, 0), 65522)
        self.assertEqual(round(st.sM.kJkmolK, 2), 221.44)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 27.93)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 36.25)
        self.assertEqual(round(st.w, 1), 878.6)

        st = Air(T=2000, P=500e6)
        self.assertEqual(round(st.rhoM, 2), 16.48)
        self.assertEqual(round(st.uM.kJkmol, 0), 48857)
        self.assertEqual(round(st.hM.kJkmol, 0), 79198)
        self.assertEqual(round(st.sM.kJkmolK, 2), 188.66)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 29.07)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 37.27)
        self.assertEqual(round(st.w, 1), 1497.6)

        # Custom cycle
        P = 50   # MPa
        T = 470  # K
        f_pt = Air(P=P, T=T)
        f_prho = Air(P=f_pt.P, rho=f_pt.rho)
        self.assertEqual(round(f_prho.P-P, 6), 0)
        self.assertEqual(round(f_prho.T-T, 6), 0)

    def test_LemmonTransport(self):
        """Table V, pag 28"""
        # Viscosity
        self.assertEqual(round(Air(T=100, rhom=0).mu.muPas, 5), 7.09559)
        self.assertEqual(round(Air(T=300, rhom=0).mu.muPas, 4), 18.5230)
        self.assertEqual(round(Air(T=100, rhom=28).mu.muPas, 3), 107.923)
        self.assertEqual(round(Air(T=200, rhom=10).mu.muPas, 4), 21.1392)
        self.assertEqual(round(Air(T=300, rhom=5).mu.muPas, 4), 21.3241)
        self.assertEqual(round(Air(T=132.64, rhom=10.4).mu.muPas, 4), 17.7623)

        # Thermal Conductivity
        self.assertEqual(round(Air(rhom=0, T=100).k.mWmK, 5), 9.35902)
        self.assertEqual(round(Air(rhom=0, T=300).k.mWmK, 4), 26.3529)
        self.assertEqual(round(Air(rhom=28, T=100).k.mWmK, 3), 119.221)
        self.assertEqual(round(Air(rhom=10, T=200).k.mWmK, 4), 35.3185)
        self.assertEqual(round(Air(rhom=5, T=300).k.mWmK, 4), 32.6062)
        self.assertEqual(round(Air(rhom=10.4, T=132.64).k.mWmK, 4), 75.6231)
