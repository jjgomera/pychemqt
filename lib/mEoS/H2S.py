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


from math import exp
from unittest import TestCase

from lib import unidades
from lib.meos import MEoS
from lib.mEoS.C3 import C3


class H2S(MEoS):
    """Multiparameter equation of state for hydrogen sulfide"""
    name = "hydrogen sulfide"
    CASNumber = "7783-06-4"
    formula = "H2S"
    synonym = ""
    _refPropName = "H2S"
    _coolPropName = "HydrogenSulfide"
    rhoc = unidades.Density(347.2841672)
    Tc = unidades.Temperature(373.1)
    Pc = unidades.Pressure(9000.0, "kPa")
    M = 34.08088  # g/mol
    Tt = unidades.Temperature(187.7)
    Tb = unidades.Temperature(212.85)
    f_acent = 0.1005
    momentoDipolar = unidades.DipoleMoment(0.97, "Debye")
    id = 50

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1.5],
           "ao_pow": [-4.0740770957, 3.7632137341, -0.002753352822675789],
           "ao_exp": [1.1364, 1.9721],
           "titao": [1823/Tc, 3965/Tc]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [9.336197742, -16.266508995],
           "ao_sinh": [3.11942], "sinh": [4.914580541],
           "ao_cosh": [1.00243], "cosh": [2.27065398]}

    Fi3 = {"ao_log": [1, 3.],
           "pow": [0, 1], "ao_pow": [7.881037, -3.20986],
           "ao_exp": [0.9767422, 2.151898], "titao": [4.506266, 10.15526]}

    CP3 = {"ao": 4.1012105,
           "an": [-1.6720073e-3, 7.5303152e-6, -0.62421053e-8, 0.18098453e-11],
           "pow": [1, 2, 3, 4]}

    f = 8.3159524/4.184
    CP4 = {"ao": 7.9468*f,
           "ao_sinh": [-1.5769761e4*f, 1.3861204e7*f],
           "sinh": [4.33801e2, 1.48143e3],
           "ao_cosh": [2.0329947e6*f, -3.5044957e6*f],
           "cosh": [8.43792e2, 1.10223e3]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for hydrogen sulfide of"
                    " Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 32,

        "nr1": [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.00019076],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.0045348, -0.22230, -0.034714, -0.014885, .0074154],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    sakoda = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of "
                    "Sakoda and Uematsu (2004)",
        "__doi__": {"autor": "Sakoda, N., Uematsu, M.",
                    "title": "A Thermodynamic Property Model for Fluid Phase "
                             "Hydrogen Sulfide",
                    "ref": "Int. J. Thermophys., 25(3) (2004) 709-737",
                    "doi": "10.1023/B:IJOT.0000034234.06341.8a"},

        "R": 8.314472,
        "cp": Fi3,
        "ref": {"Tref": 298.15, "Pref": 100., "ho": 0, "so": 0},

        "Tmin": 187.67, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.13,

        "nr1": [0.1545780, -0.1717693e1, -0.1595211e1, 0.2046589e1,
                -0.1690358e1, 0.9483623, -0.6800772e-1, 0.4372273e-2,
                0.3788552e-4, -0.3680980e-4, 0.8710726e-5],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 8, 9, 10],
        "t1": [.241, .705, 1, .626, 1.15, 1.63, 0.21, 3.08, 0.827, 3.05, 3.05],

        "nr2": [0.6886876, 0.2751922e1, -0.1492558e1, 0.9202832, -0.2103469,
                0.1084359e-2, 0.3754723e-1, -0.5885793e-1, -0.2329265e-1,
                -0.1272600e-3, -0.1336824e-1, 0.1053057e-1],
        "d2": [1, 1, 1, 2, 5, 1, 4, 4, 3, 8, 2, 3],
        "t2": [0.11, 1.07, 1.95, 0.142, 2.13, 4.92, 1.75, 3.97, 11.8, 10,
               9.83, 14.2],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4],
        "gamma2": [1]*12}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Polt "
                    "et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 142000.0, "rhomax": 29.1,

        "nr1": [0.135782366339e1, -0.153224981014e1, 0.329107661253,
                0.195802782279e1, -0.301125182071e1, -0.126614059078e1,
                0.129960331548e1, -0.185645977138, -0.160919744092e1,
                0.234395817019e1, -0.378573094883, 0.758423219040,
                -0.973372615169, -0.120786235447, 0.209004959689,
                -0.919656385346e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.135782366339e1, 0.153224981014e1, -0.329107661253,
                0.891427552242, -0.204776100441e1, 0.101366381241e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9873538]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of "
                    "Starling (1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 29.578,

        "nr1": [0.110928333109e1, 0.188834546108, -0.930906931583,
                -0.411249591635, 0.140676923412e-1, -0.169077883177e-4,
                0.510265859853, -0.572402742986, -0.828859606622e-3,
                0.971664064871e-2, 0.140700425434e-4],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.110928333109e1, -0.26913741657],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.48524558]*2}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.12,

        "nr1": [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.19076e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.45348e-2, -0.2223, -0.034714, -.014885, .74154e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon, sakoda, polt, starling, GERG

    _surface = {"sigma": [0.078557], "exp": [1.2074]}
    _dielectric = {
        "__doi__": {
            "autor": "Harvey, A.H., Mountain, R.D.",
            "title": "Correlations for the Dielectric Constants of H2S, SO2 "
                     "and SF6",
            "ref": "Int. J. Thermophys. 38 (2017) 147",
            "doi": "10.1007/s10765-017-2279-6"},
        "eq": 3,

        "alfa": 3.66e-24, "mu": 0.978325, "Tc": 380, "rhoc": 0.029,
        "cu": 0.241, "au": 2.83, "nu": 0.5, "cg": 1.18,
        "a1": 2.546, "n1": 0.90, "a2": 1.883, "n2": 3.5}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.5884, 2.1582, -1.6054, -2.3870],
        "t": [1, 1.5, 2, 4.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [11.833, -17.019, 7.8047],
        "t": [1.63/3, 1.95/3, 2.3/3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.9156, -7.7093, -19.543, -49.418],
        "t": [0.49, 1.69, 4, 8]}

    visco0 = {"__name__": "Quiñones-Cisneros (2012)",
              "__doi__": {
                  "autor": "Quinones-Cisneros, S.E., Schmidt, K.A.G., Giri, B."
                           "R, Blais, P., Marriott, R.A.,",
                  "title": "Reference Correlation for the Viscosity Surface "
                           "of Hydrogen Sulfide",
                  "ref": "J. Chem. Eng. Data 57(11) (2012) 3014-3018",
                  "doi": "10.1021/je300601h"},

              "eq": 1, "omega": 4,

              "ek": 276, "sigma": 0.3565,
              "n_chapman": 0.87721/M**0.5*0.3565**2,
              "collision": [0.53242, 0.93715, -0.69339, 1.16432, -0.84306,
                            0.20534],

              "Tref_virial": 355.8,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "special": "_mur"}

    def _mur(self, rho, T, fase):
        # Modified friction theory for residual viscosity contribution

        Gamma = self.Tc/T
        psi1 = exp(Gamma)                                               # Eq 15
        psi2 = exp(Gamma**2)                                            # Eq 16

        a = [68.9659e-6, -22.0494e-6, -42.6126e-6]
        b = [153.406e-6, 8.45198e-6, -113.967e-6]
        A = [0.78238e-9, -0.64717e-9, 1.39066e-9]
        B = [-9.75792e-9, -3.19303e-9, 12.4263e-9]

        ka = (a[0] + a[1]*psi1 + a[2]*psi2) * Gamma                     # Eq 11
        kr = (b[0] + b[1]*psi1 + b[2]*psi2) * Gamma                     # Eq 12
        kaa = (A[0] + A[1]*psi1 + A[2]*psi2) * Gamma                    # Eq 13
        krr = (B[0] + B[1]*psi1 + B[2]*psi2) * Gamma                    # Eq 14

        # All parameteres has pressure units of bar
        Patt = -fase.IntP.bar
        Prep = T*fase.dpdT_rho.barK
        Pid = rho*self.R*self.T/1e5
        delPr = Prep-Pid

        # Eq 5
        mur = kr*delPr + ka*Patt + krr*Prep**2 + kaa*Patt**2
        return mur*1e3

    visco1 = {"__name__": "Schmidt (2007)",
              "__doi__": {
                  "autor": "Schmidt, K.A.G., Carroll, J.J., Quinones-Cisneros,"
                           " S.E., Kvamme, B.",
                  "title": "Hydrogen Sulfide Viscosity Modeling",
                  "ref": "Energy & Fuels 22(5) (2008) 3424-3434",
                  "doi": "10.1021/ef700701h"},

              "eq": 4, "omega": 0,

              "Toref": 373.1,
              "no": [43.6694, -121.530, 93.5279],
              "to": [0, 0.25, 0.5],

              "a": [5.46919e-5, -7.32295e-6, -7.35622e-6],
              "b": [4.56159e-5, -1.82572e-5, -6.59654e-6],
              "c": [-4.33882e-6, 6.13716e-6, 0.0],
              "A": [6.67324e-9, -2.16365e-9, 0.0],
              "B": [-1.53973e-9, 2.17652e-9, 0.0],
              "C": [3.54228e-7, -4.76258e-8, 0.0]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 296.28, "sigma": 0.3732, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1], "psi_d": [0],
              "fint": [1.5603e-4, 1.78874e-6, -6.75136e-10],
              "fint_t": [0, 1, 2],
              "chi": [1], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.164e-9, "gam0": 0.058, "qd": 0.447e-9, "Tcref": 1.5*Tc}

    _viscosity = visco0, visco1, trnECS
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = H2S(T=375, rhom=10)
        self.assertEqual(round(st.P.kPa, 3), 9289.914)
        self.assertEqual(round(st.hM.kJkmol, 3), 15571.227)
        self.assertEqual(round(st.sM.kJkmolK, 3), 49.880)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 42.574)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2645.392)
        self.assertEqual(round(st.w, 3), 248.512)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = H2S(T=335.8, rhom=20.154, visco=2)
        # self.assertEqual(round(st.mu.muPas, 5), 97.16912)
        self.assertEqual(round(st.mu.muPas, 5), 97.16763)
        self.assertEqual(round(st.k.mWmK, 4), 121.2579)

    def test_dielectric(self):
        self.assertEqual(round(H2S(T=400, rhom=30).epsilon, 4), 5.0239)
