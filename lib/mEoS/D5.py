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


from unittest import TestCase

from lib import unidades
from lib.meos import MEoS
from lib.mEoS import N2


class D5(MEoS):
    """Multiparameter equation of state for decamethylcyclopentasiloxane"""
    name = "decamethylcyclopentasiloxane"
    CASNumber = "541-02-6"
    formula = "C10H30O5Si5"
    synonym = "D5"
    _refPropName = "D5"
    _coolPropName = "D5"
    rhoc = unidades.Density(300.323457)
    Tc = unidades.Temperature(618.3)
    Pc = unidades.Pressure(1077.7, "kPa")
    M = 370.7697  # g/mol
    Tt = unidades.Temperature(224.65)
    Tb = unidades.Temperature(484.099)
    f_acent = 0.631
    momentoDipolar = unidades.DipoleMoment(1.349, "Debye")
    # id=1671

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [94.3892428631, -31.1102222402],
           "ao_exp": [51, 57.9, 35],
           "titao": [221/Tc, 1733/Tc, 4544/Tc]}

    f = 8.314472
    CP1 = {"ao": -34.898/f,
           "an": [1861.5e-3/f, -1403.4e-6/f, 500.0e-9/f],
           "pow": [1, 2, 3]}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for "
                    "decamethylcyclopentasiloxane of Thol (2019).",
        "__doi__": {
            "autor": "Thol, M., Javed, M.A., Baumhögger, E., Span, R., "
                     "Vrabec, J.",
            "title": "Thermodynamic Properties of Dodecamethylpentasiloxane, "
                     "Tetradecamethylhexasiloxane, and "
                     "Decamethylcyclopentasiloxane",
            "ref": "Ind. Eng. Chem. Res. 58(22) (2019) 9617-9635",
            "doi": "10.1021/acs.iecr.9b00608"},

        "R": 8.3144598,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 655, "Pmax": 125000,

        "nr1": [0.0177345, 4.3133088, -6.1586863, -1.4503945, 0.9519342],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.43, 0.754, 0.84, 0.72],

        "nr2": [-2.3848036, -1.4114529, 0.7255071, -2.9966803, -0.0902228],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.35, 2.58, 0.66, 1.71, 1.0163],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [6.3033323, -1.0592923, 0.79365281, -1.8982515, -0.01351964],
        "d3": [1, 3, 2, 2, 1],
        "t3": [1.114, 1.85, 0.9, 1.05, 1.09],
        "alfa3": [1.046, 0.993, 0.545, 1.128, 13.9],
        "beta3": [0.37, 0.11, 0.1, 0.37, 519],
        "gamma3": [1.626, 1.05, 1.11, 1.22, 1.083],
        "epsilon3": [0.787, 0.567, 0.685, 0.577, 0.936]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for "
                    "decamethylcyclopentasiloxane of Colonna (2006).",
        "__doi__": {
            "autor": "Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W.",
            "title": "Multiparameter Equations of State for Selected "
                     "Siloxanes",
            "ref": "Fluid Phase Equilibria, 244 (2006) 193-211",
            "doi": "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",
        "Tc": 619.15, "Pc": 1160, "rhoc": 1/1.216,

        "Tmin": 300, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.83,

        "nr1": [1.40844725, -2.29248044, 0.42851607, -0.73506382, 0.16103808,
                0.29643278e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.82412481, 0.15214274, -0.68495890, -0.55703624e-1,
                0.13055391e-1, -0.31853761e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = thol, colonna

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.04408], "exp": [1.357]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-9.256, 3.987, -11.02, -19.286, 16.524, -8.14],
        "t": [1, 1.5, 2.24, 3.48, 2.86, 11.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.0938, 5.254, -12.31, 19.364, -15.81, 5.983],
        "t": [0.25, 0.79, 1.33, 1.9, 2.52, 3.22]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.916, -5.911, -18.617, -74.29, -154.4, -284.1],
        "t": [0.23, 0.68, 2.24, 5.1, 10.7, 18.9]}

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

              "ek": 491, "sigma": 0.864, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [-2.49055, 4.63356, -1.89292, 0.247782],
              "psi_d": [0, 1, 2, 3],
              "fint": [0.00132], "fint_t": [0],
              "chi": [1.40287, 0.0940128], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.319e-9, "gam0": 0.064, "qd": 1.068e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_thol(self):
        """Table 6, Pag. 9621"""
        st = D5(T=290, rhom=2.7)
        self.assertEqual(round(st.P.MPa, 7), 36.3487297)
        self.assertEqual(round(st.hM.Jmol, 3), -122272.731)
        self.assertEqual(round(st.sM.JmolK, 6), -359.629958)
        self.assertEqual(round(st.w, 5), 1151.09861)
        self.assertEqual(round(st.aM.Jmol, 4), -31442.5359)

        st = D5(T=390, rhom=0.001)
        self.assertEqual(round(st.P.MPa, 10), 0.0032226439)
        self.assertEqual(round(st.hM.Jmol, 4), -14185.4999)
        self.assertEqual(round(st.sM.JmolK, 7), -14.0834572)
        self.assertEqual(round(st.w, 7), 93.6614237)
        self.assertEqual(round(st.aM.Jmol, 4), -11915.5955)

        st = D5(T=450, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 10), 0.0358844583)
        self.assertEqual(round(st.hM.Jmol, 4), 20404.0115)
        self.assertEqual(round(st.sM.JmolK, 7), 48.6842603)
        self.assertEqual(round(st.w, 7), 97.0959266)
        self.assertEqual(round(st.aM.Jmol, 5), -5092.35152)

        st = D5(T=450, rhom=2.5)
        self.assertEqual(round(st.P.MPa, 7), 77.0798056)
        self.assertEqual(round(st.hM.Jmol, 5), -4880.23864)
        self.assertEqual(round(st.sM.JmolK, 7), -81.6230026)
        self.assertEqual(round(st.w, 5), 1044.97883)
        self.assertEqual(round(st.aM.Jmol, 5), 1018.19028)

        st = D5(T=650, rhom=1.8)
        self.assertEqual(round(st.P.MPa, 7), 14.8882334)
        self.assertEqual(round(st.hM.Jmol, 3), 129408.704)
        self.assertEqual(round(st.sM.JmolK, 6), 215.447596)
        self.assertEqual(round(st.w, 6), 415.207142)
        self.assertEqual(round(st.aM.Jmol, 4), -18903.4744)

    # def test_Huber(self):
    #     """Table 7, pag 266"""
    #     st = D5(T=556.5, rhom=1.718)
    #     self.assertEqual(round(st.mu.muPas, 4), 184.5268)
    #     self.assertEqual(round(st.k.mWmK, 4), 59.4893)
