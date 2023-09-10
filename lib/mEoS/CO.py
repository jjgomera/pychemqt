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

from lib import unidades
from lib.meos import MEoS
from lib.mEoS.N2 import N2


class CO(MEoS):
    """Multiparameter equation of state for carbon monoxide"""
    name = "carbon monoxide"
    CASNumber = "630-08-0"
    formula = "CO"
    synonym = ""
    _refPropName = "CO"
    _coolPropName = "CarbonMonoxide"
    rhoc = unidades.Density(303.909585)
    Tc = unidades.Temperature(132.86)
    Pc = unidades.Pressure(3494.0, "kPa")
    M = 28.0101  # g/mol
    Tt = unidades.Temperature(68.16)
    Tb = unidades.Temperature(81.64)
    f_acent = 0.0497
    momentoDipolar = unidades.DipoleMoment(0.1, "Debye")
    id = 48

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1.5],
           "ao_pow": [-3.3728318564, 3.3683460039, -9.111274701235156e-5],
           "ao_exp": [1.0128],
           "titao": [3089/Tc]}

    Fi2 = {"ao_log": [1, 2.50055],
           "pow": [0, 1],
           "ao_pow": [10.813340744, -19.834733959],
           "ao_exp": [], "titao": [],
           "ao_sinh": [1.02865], "sinh": [11.6698028],
           "ao_cosh": [0.00493], "cosh": [5.302762306]}

    CP3 = {"ao": 0.36028218e1,
           "an": [-0.20871594e5, 0.89208708e3, -0.14157993e2, -0.34021345e-3,
                  0.44616091e-6, -0.15154703e-9],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [0.90426143], "exp": [30000]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbon monoxide of "
                    "Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500., "Pmax": 100000.0, "rhomax": 33.84,

        "nr1": [0.90554, -2.4515, 0.53149, 0.024173, 0.072156, 0.00018818],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.19405, -0.043268, -0.12778, -0.027896, -0.034154, 0.016329],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    mccarty = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for carbon monoxide of McCarty "
                    "(1989)",
        "__doi__": {"autor": "McCarty, R.D.",
                    "title": "Correlations for the Thermophysical Properties "
                             "of Carbon Monoxide",
                    "ref": "NIST, Boulder, CO, 1989",
                    "doi": ""},

        "R": 8.31434,
        "M": 28.011, "Tc": 132.8, "Pc": 3493.5, "rhoc": 10.85,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000., "Pmax": 30000.0, "rhomax": 30.25,

        "b": [None, 0.8845582109949e-2, -0.2236741566840, 0.1742275796442e1,
              -0.2169146998363e3, 0.1721504267082e4, -0.3990514770703e-4,
              0.1036880040451, -0.3376308165071e2, 0.2061895161095e5,
              0.2993711656350e-5, 0.1856003597097e-2, -0.2114419664527,
              -0.2436986935194e-5, -0.1858029609177e-2, -0.1734563867767e1,
              0.1509970839260e-3, -0.2282721433205e-5, 0.2202780295674e-2,
              -0.3313357789163e-4, -0.1473412120276e5, -0.3141136651147e6,
              -0.1451168999234e3, 0.6323441221817e5, -0.2203560539926,
              -0.2087738308480e2, -0.1508165207553e-2, 0.2740740634030e1,
              0.8687687989627e-6, -0.1451419251928e-3, -0.3040346241285e-8,
              0.4712050805815e-8, -0.2639772456566e-5]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon monoxide of Kunz "
                    "and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 500., "Pmax": 100000.0, "rhomax": 33.84,

        "nr1": [0.92310041400851, -0.248858452058e1, 0.58095213783396,
                0.028859164394654, 0.070256257276544, 0.21687043269488e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 0.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.13758331015182, -0.51501116343466e-1, -0.14865357483379,
                -0.03885710088681, -0.029100433948943, 0.14155684466279e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon, mccarty, GERG
    _PR = [-1.1087, -10.2779]

    _surface = {"sigma": [0.02843],
                "exp": [1.148]}

    _melting = {
        "eq": 1,
        "__doi__": {"autor": "Barreiros, S.F., Calado, J.C.G., Nunes da "
                             "Ponte, M.",
                    "title": "The melting curve of carbon monoxide",
                    "ref": "J. Chem. Thermodynamics 14 (1982) 1197-1198",
                    "doi": "10.1016/0021-9614(82)90044-1"},

        "Tmin": Tt, "Tmax": 1000.0,
        "Tref": 1, "Pref": 1e6,
        "a0": -142.941,
        "a1": [0.0195608], "exp1": [2.10747]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.61192e1, 0.10411e1, -0.62162e1, 0.10437e2, -0.76813e1],
        "t": [1.0, 1.5, 3.9, 4.6, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.29570e1, -0.42880e1, 0.87643e1, -0.84001e1, 0.36372e1],
        "t": [0.398, 0.735, 1.08, 1.5, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.25439e1, -0.55601e1, -0.85276e1, -0.51163e1, -0.17701e2,
              -0.29858e2],
        "t": [0.395, 1.21, 3.0, 3.5, 6.0, 8.0]}

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

              "ek": 103.697, "sigma": 0.3615, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.07369, -0.0283067], "psi_d": [0, 1],
              "fint": [3.29558e-4, 3.05976e-6, -3.13222e-9],
              "fint_t": [0, 1, 2],
              "chi": [1.00037, -0.0082682], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.164e-9, "gam0": 0.059, "qd": 0.437e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = CO(T=134, rhom=10)
        self.assertEqual(round(st.P.kPa, 3), 3668.867)
        self.assertEqual(round(st.hM.kJkmol, 3), 4838.507)
        self.assertEqual(round(st.sM.kJkmolK, 3), 41.601)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 38.702)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1642.142)
        self.assertEqual(round(st.w, 3), 168.632)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = CO(T=119.6, rhom=21.287)
        self.assertEqual(round(st.mu.muPas, 5), 56.70495)
        self.assertEqual(round(st.k.mWmK, 4), 75.7297)
