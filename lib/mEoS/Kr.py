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


from math import exp, log
from unittest import TestCase

from scipy.constants import Avogadro, Boltzmann
from lib import unidades
from lib.meos import MEoS
from lib.mEoS.N2 import N2


class Kr(MEoS):
    """Multiparameter equation of state for krypton"""
    name = "krypton"
    CASNumber = "7439-90-9"
    formula = "Kr"
    synonym = "R-784"
    _refPropName = "KRYPTON"
    _coolPropName = "Krypton"
    rhoc = unidades.Density(909.2083)
    Tc = unidades.Temperature(209.48)
    Pc = unidades.Pressure(5525.0, "kPa")
    M = 83.798  # g/mol
    Tt = unidades.Temperature(115.775)
    Tb = unidades.Temperature(119.73)
    f_acent = -0.000894
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 971

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-3.7506412806, 3.7798018435]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for krypton of "
                    "Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data 51(3) (2006) 785-850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 200000.0, "rhomax": 33.42,

        "nr1": [0.83561, -2.3725, 0.54567, 0.014361, 0.066502, 0.00019310],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.16818, -0.033133, -0.15008, -0.022897, -0.021454, 0.0069397],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for krypton of Polt (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 - 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 780.0, "Pmax": 375000.0, "rhomax": 33.55,

        "nr1": [-0.402218741560, 0.679250544381, -0.1878869802860,
                0.603399982935, -0.177297564389e1, 0.581208430222,
                -0.733585469788, 0.164651929067, -0.319923148922e-1,
                0.333278228743, 0.219652478083e-1, 0.751994891628e-1,
                -0.212109737251, -0.645185506524e-2, 0.409175610200e-1,
                0.169416098754e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.402218741560, -0.679250544381, 0.187886980286,
                0.108265263587, -0.137102675805, -0.110549803007],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2, 2, 2, 2, 2, 2],
        "gamma2": [1]*6}

    eq = lemmon, polt

    _surface = {"sigma": [0.0447], "exp": [1.245]}
    _dielectric = {
        "eq": 1,
        "a": [6.273, 0], "b": [6.485, 13.48], "c": [-82.51, -170.4],
        "Au": 0, "D": 1.7}

    _melting = {
        "eq": 2,
        "__doi__": {"autor": "Michels, A., Prins, C.",
                    "title": "The Melting Lines of Argon, Krypton and Xenon "
                             "up to 1500 Atm; Representation of the Results "
                             "by a Law of Corresponding States",
                    "ref": "Physica 28 (1962) 101-116",
                    "doi": "10.1016/0031-8914(62)90096-4"},

        "Tmin": Tt, "Tmax": 800.0,
        "Tref": 1, "Pref": -2345*101325,
        "a1": [1.08047668519*101325], "exp1": [1.6169841]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.59697e1, 0.12673e1, -0.95609, -0.35630e2, 0.56884e2],
        "t": [1.0, 1.5, 2.95, 9.3, 10.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.20593e2, -0.65490e2, 0.94407e2, -0.69678e2, 0.22810e2],
        "t": [0.62, 0.84, 1.07, 1.34, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.64163e1, .89956e1, -.10216e2, -.13477e2, -.21152e3, .21375e3],
        "t": [0.525, 0.77, 1.04, 3.2, 8.3, 9.0]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": N2,
              "visco": "visco1",

              "ek": 178.9, "sigma": 0.3655, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1.008291,

              "psi": [1], "psi_d": [0],
              "fint": [0.00132], "fint_t": [0],
              "chi": [0.962573, -0.0118156], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.168e-9, "gam0": 0.058, "qd": 0.437e-9, "Tcref": 1.5*Tc}

    visco1 = {"__name__": "Polychroniadou (2022)",

              "__doi__": {
                  "autor": "Polychroniadou, S., Antoniadis, K.D., Assael, "
                           "M.J., Bell, I.H.",
                  "title": "A Reference Correlation for the Viscosity of "
                           "Krypton From Entropy Scaling",
                  "ref": "Int. J. Thermophys. 43 (2022) 6",
                  "doi": "10.1007/s10765-021-02927-5"},

              "eq": 1,
              "special": "_mur"}

    def _mu0(self, T):
        """Special term for zero-density viscosity for Polychroniadou
        correlation"""
        # Parameters get this reference
        # X. Xiao, X., Rowland, D., Ghafri, Al Ghafri, S.Z.S., May, E.F.
        # Wide-Ranging Reference Correlations for Dilute Gas Transport
        # Properties Based on Ab Initio Calculations and Viscosity Ratio
        # Measurements
        # J. Phys. Chem. Ref. Data 49(1) (2020) 013101
        # https://doi.org/10.1063/1.5125100
        a = [9.129712e-1, -1.001470e-1, -2.454742e-2, 3.145009e-2,
             -4.456257e-3, -4.511243e-3, 2.237544e-3, -1.455422e-4,
             -2.006385e-4, 8.341288e-5, -1.520236e-5, 1.159085e-6]

        # Eq 8
        muo = exp(sum(ai*log(T/298.15)**(i+1) for i, ai in enumerate(a)))
        return 25.3062*muo

    def _mur(self, rho, T, fase=None):
        """Hardcoded viscosity correlation from Polychroniadou"""
        muo = self._mu0(T)

        if rho and fase:
            B = fase.virialB
            dBT = fase.dBt
            teta = (B+T*dBT)/Avogadro*self.M/1000
            mu_dl = muo*1e-6/(self.M*1e-3/Avogadro*Boltzmann*T)**0.5*teta**(2/3)

            s = self.R.kJkgK*(self.Tc/T*fase.firt-fase.fir)*self.M
            splus = -s/self.R.kJkgK/self.M
            clj = [0.125364, 0.220795, -0.0313726, 0.00313907]
            lj = exp(sum(ci*splus**(i+1) for i, ci in enumerate(clj))) - 1
            mu_res = lj * 1.05

            muPlus = mu_res + mu_dl
            rhon = rho/self.M*1000*Avogadro
            return 1e6*rhon**(2/3)*(self.M/1000/Avogadro*Boltzmann*T)**0.5*muPlus/splus**(2/3)

        else:
            return muo

    _viscosity = (trnECS, visco1)
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = Kr(T=211, rhom=10)
        self.assertEqual(round(st.P.kPa, 3), 5741.445)
        self.assertEqual(round(st.hM.kJkmol, 3), 6700.326)
        self.assertEqual(round(st.sM.kJkmolK, 3), 36.936)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 27.390)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1667.678)
        self.assertEqual(round(st.w, 3), 137.838)

    def test_Michels(self):
        """Table I, pag 105"""
        self.assertEqual(round(Kr._Melting_Pressure(115.893).atm, 2), 5.52)
        self.assertEqual(round(Kr._Melting_Pressure(119.069).atm, 2), 110.55)
        self.assertEqual(round(Kr._Melting_Pressure(138.598).atm, 2), 794.08)
        self.assertEqual(round(Kr._Melting_Pressure(157.033).atm, 2), 1496.47)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Kr(T=188.5, rhom=21.274)
        self.assertEqual(round(st.mu.muPas, 4), 115.6593)
        self.assertEqual(round(st.k.mWmK, 4), 48.3068)

    def test_Polychroniadou(self):
        st = Kr(T=200, rhom=1e-9, visco=1)
        self.assertEqual(round(st.mu.muPas, 8), 17.33865199)

        # This point isn't real, it's in two phases region so need force
        # calculation
        # st = Kr(T=200, rhom=13.020, visco=1)
        # self.assertEqual(round(mu.muPas, 15), 56.44764255291952)

        st = Kr(T=298.15, rhom=1e-9, visco=1)
        self.assertEqual(round(st.mu.muPas, 11), 25.30620000081)

        st = Kr(T=400, rhom=1e-9, visco=1)
        self.assertEqual(round(st.mu.muPas, 12), 32.79555841859)

        # Difference wieh paper by Avogadro value
        st = Kr(T=400, rhom=13.020, visco=1)
        self.assertEqual(round(st.mu.muPas, 5), 64.80147)
