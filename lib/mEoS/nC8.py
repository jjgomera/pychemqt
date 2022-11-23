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


class nC8(MEoS):
    """Multiparameter equation of state for n-octane"""
    name = "octane"
    CASNumber = "111-65-9"
    formula = "CH3-(CH2)6-CH3"
    synonym = ""
    _refPropName = "OCTANE"
    _coolPropName = "n-Octane"
    rhoc = unidades.Density(234.9)
    Tc = unidades.Temperature(569.32)
    Pc = unidades.Pressure(2497.0, "kPa")
    M = 114.2285  # g/mol
    Tt = unidades.Temperature(216.37)
    Tb = unidades.Temperature(398.77)
    f_acent = 0.395
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 12
    _Tr = unidades.Temperature(565.427917)
    _rhor = unidades.Density(234.605116)
    _w = 0.402698435

    CP1 = {"ao": 4,
           "ao_sinh": [15.6865, 48.1731], "sinh": [158.9220, 1693.07],
           "ao_cosh": [33.8029], "cosh": [815.064]}

    Fi1 = {"ao_log": [1, 3.0],
           "ao_pow": [15.864687161, -97.370667555], "pow": [0, 1],
           "ao_sinh": [15.6865, 48.1731], "sinh": [158.9220/Tc, 1693.07/Tc],
           "ao_cosh": [33.8029], "cosh": [815.064/Tc]}

    CP3 = {"ao": 3.018753,
           "an": [0.07297005, -0.14171168e-4, -0.1225317e-7,  0.12912645e-11],
           "pow": [1, 2, 3, 4]}

    f = 8.3159524/4.184
    CP4 = {"ao": 34.0847*f,
           "ao_sinh": [2.603664e8*f], "sinh": [1.6115500e3],
           "ao_cosh": [4.1241363e7*f], "cosh": [7.6884700e2]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for octane of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",
        "M": 114.231, "Tc": 569.32, "rhoc": 234.9/114.231,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 6.69,

        "nr1": [0.10722545e1, -0.24632951e1, 0.65386674, -0.36324974,
                0.1271327, 0.30713573e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.52656857, 0.19362863e-1, -0.58939427, -0.14069964,
                -0.78966331e-2, 0.33036598e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-octane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 6.69,

        "nr1": [0.10722544875633e1, -0.24632951172003e1, 0.65386674054928,
                -0.36324974085628, 0.12713269626764, 0.30713572777930e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.5265685698754, 0.19362862857653e-1, -0.58939426849155,
                -0.14069963991934, -0.78966330500036e-2, 0.33036597968109e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octane of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "M": 114.233, "Tc": 569.35, "Pc": 2517, "rhoc": 2.0571989,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": 258.0, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 6.6355607,

        "nr1": [0.266117347782e1, -0.343810366899e1, 0.700476763325,
                0.573101545749e1, -0.411975339382e1, -0.771251551395e1,
                0.526137115388e1, -0.716144047789, -0.584632875151e1,
                0.736422551908e1, -0.100540027381e1, 0.158387242200e1,
                -0.153643650819e1, -0.142010818863, 0.333126039209e-1,
                0.271948869925e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.266117347782e1, 0.343810366899e1, -0.700476763325,
                0.443217980268e1, -0.123858312597e2, 0.803373487925e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9995725]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octane of Starling "
                    "(1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": 255.372, "Tmax": 644.0, "Pmax": 55000.0, "rhomax": 6.36203,

        "nr1": [0.253526486527e1, 0.616872653050, -0.941731168114,
                -0.109609729872e1, 0.849362892312e-1, -0.363538456997e-3,
                0.849748115039e-1, -0.961236603829e-1, -0.132591135067,
                0.269748328453e-2, 0.372085674947e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.253526486527e1, -0.447291258549],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.35285564]*2}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [1.57750154, 1.15745614, -3.54867092, 1.18030671e-1,
                3.02753897e-4, -2.63074957e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [2.55299486e-2, -1.26632996e-1, 4.48343319e-1, -9.46702997e-3,
                -0.443927529, -1.68224827e-2, -1.15864640e-1, -1.32417591e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = shortSpan, GERG, polt, starling, sun
    _PR = [0.0526, -21.0808]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.1708, -0.2337, 0.1194], "exp": [1.466, 1.81, 2.127]}
    _dielectric = {
        "eq": 1,
        "a": [39.74, 0.040], "b": [348.01, 494.18], "c": [-76838, -65772],
        "Au": 29.84, "D": 2}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.79713e1, 0.17915e1, -0.34540e1, -0.82509e1, 0.53357e1],
        "t": [1.0, 1.5, 2.64, 5.5, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.56814e1, 0.38908e2, -0.75923e2, 0.59548e2, -0.19651e2],
        "t": [0.1, 0.75, 0.9, 1.1, 1.25]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.16556, -5.9337, -18.915, -0.36484e3, 0.72686e3, -0.50392e3],
        "t": [0.09, 0.59, 2.4, 7.0, 8.0, 9.0]}

    visco0 = {"__name__": "Huber (2004)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A. Xiang, H.W.",
                  "title": "Viscosity correlations for minor constituent "
                           "fluids in natural gas: n-octane, n-nonane and "
                           "n-decane",
                  "ref": "Fluid Phase Equilibria 224 (2004) 263-270",
                  "doi": "10.1016/j.fluid.2004.07.012"},

              "eq": 1, "omega": 1,

              "ek": 452.09, "sigma": 0.63617,
              "n_chapman": 0.021357,
              "collision": [0.335103, -0.467898],

              "Tref_virial": 452.09,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 569.32, "rhoref_res": 234.9, "muref_res": 1000,
              "nr": [-0.103924, 0.113327e-1, 0.992302e-1, -0.322455e-1],
              "tr": [1, 1, 2, 2],
              "dr": [2, 3, 2, 3],

              "CPf": 606.122,
              "CPg1": 2.0651,
              "CPgi": [3.07843/2.0651, -0.879088/2.0651],
              "CPti": [-0.5, -1]}

    visco1 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 569.32,
              "no": [16.7562, -53.1705, 46.9105],
              "to": [0, 0.25, 0.5],

              "a": [8.68736e-5, 0.0, -2.69591e-5],
              "b": [1.46267e-4, 0.0, -5.44584e-5],
              "c": [1.28673e-4, -1.76442e-5, 0.0],
              "A": [-2.40884e-9, 5.20715e-11, 0.0],
              "B": [0.0, 6.62141e-9, 1.60012e-9],
              "C": [-9.50545e-7, 1.03767e-6, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"__name__": "Huber (2005)",
               "__doi__": {
                   "autor": "Huber, M.L., Perkins, R.A.",
                   "title": "Thermal conductivity correlations for minor "
                            "constituent fluids in natural gas: n-octane, "
                            "n-nonane and n-decane",
                   "ref": "Fluid Phase Equilibria 227 (2005) 47-55",
                   "doi": "10.1016/j.fluid.2004.10.031"},

               "eq": 1,

               "Toref": 569.32, "koref": 1,
               "no": [0.772930e-2, -0.371138e-1, 0.977580e-1, -0.288707e-1],
               "to": [0, 1, 2, 3],

               "Tref_res": 569.32, "rhoref_res": 234.9, "kref_res": 1,
               "nr": [0.285553e-1, -0.926155e-2, -0.171398e-1, 0.659971e-2,
                      0.153496e-2],
               "tr": [0, -1, 0, 0, -1],
               "dr": [1, 1, 2, 3, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 0.68628e-9, "Tcref": 853.98}

    _thermal = thermo0,


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 46
        st = nC8(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.1537)
        self.assertEqual(round(st.P.MPa, 3), 6.363)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.8007)

        st2 = nC8(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 211.79)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.31183)

    def test_viscoHuber(self):
        # Section 3.1 pag 265
        self.assertEqual(round(nC8(T=300, P=1e7).mu.muPas, 2), 553.60)

    def test_thermoHuber(self):
        # Section 3.1 pag 04
        self.assertEqual(round(nC8(T=300, P=1e7).k.mWmK, 2), 128.36)
