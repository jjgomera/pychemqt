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


class nC6(MEoS):
    """Multiparameter equation of state for n-hexane"""
    name = "hexane"
    CASNumber = "110-54-3"
    formula = "CH3-(CH2)4-CH3"
    synonym = ""
    _refPropName = "HEXANE"
    _coolPropName = "n-Hexane"
    rhoc = unidades.Density(233.1819)
    Tc = unidades.Temperature(507.82)
    Pc = unidades.Pressure(3034.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(177.83)
    Tb = unidades.Temperature(341.86)
    f_acent = 0.299
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 10
    _Tr = unidades.Temperature(487.762087)
    _rhor = unidades.Density(235.700888)
    _w = 0.298052404

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [14.345969349, -96.165722367],
           "ao_exp": [], "titao": [],
           "ao_sinh": [11.6977, 38.6164], "sinh": [182.326/Tc, 1826.59/Tc],
           "ao_cosh": [26.8142], "cosh": [859.207/Tc]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_sinh": [11.6977, 38.6164], "sinh": [182.326, 1826.59],
           "ao_cosh": [26.8142], "cosh": [859.207]}

    CP3 = {"ao": 2.5200507,
           "an": [0.05280653, -5.7861557e-6, -1.0899040e-8, -1.8988742e-13],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    CP4 = {"ao": 26.6225/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_sinh": [2.3738446e8/8.3159524*4.184], "sinh": [1.71849e3],
           "ao_cosh": [3.5806766e7/8.3159524*4.184], "cosh": [8.02069e2]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for hexane of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",
        "M": 86.177, "Tc": 507.82, "rhoc": 233.18/86.177,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 8.85,
        "Pmin": 0.001277, "rhomin": 8.8394,

        "nr1": [0.10553238013661e1, -0.26120615890629e1, 0.76613882967260,
                -0.29770320622459, 0.11879907733358, 0.27922861062617e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.46347589844105, 0.11433196980297e-1, -0.48256968738131,
                -0.093750558924659, -0.0067273247155994, -0.0051141583585428],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-hexane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 8.85,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.10553238013661e1, -0.26120615890629e1, 0.76613882967260,
                -0.29770320622459, 0.11879907733358, 0.27922861062617e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.46347589844105, 0.011433196980297, -0.48256968738131,
                -0.093750558924659, -0.0067273247155994, -0.0051141583585428],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": 223.0, "Tmax": 623.0, "Pmax": 510000.0, "rhomax": 8.726125,
        "Pmin": 0.001277, "rhomin": 8.8394,

        "nr1": [-0.157654494847e1, 0.178731485778e1, -0.341262936801,
                0.114919468260e1, -0.381451065649e1, 0.356688884337e1,
                -0.274863278063e1, 0.391987699726, 0.346062554746,
                -0.139140552239, 0.489013943543, -0.529751545354e-1,
                -0.149303737787, 0.455990262306e-1, -0.564866336099e-1,
                0.152437539639e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.157654494847e1, -0.178731485778e1, 0.341262936801,
                0.139479099785, 0.5076238131, -0.655600474113],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.00773692]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Starling "
                    "(1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": 222.04, "Tmax": 644.0, "Pmax": 55000.0, "rhomax": 8.6724844,
        "Pmin": 0.001277, "rhomin": 8.8394,

        "nr1": [0.261128818398e1, 0.451396780770, -0.783362300734,
                -0.108785843809e1, 0.124906986929, -0.155020819852e-1,
                0.42399441457, -0.636532521368, -0.524764104726e-1,
                0.120405133154e-1, 0.992632580157e-3],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.261128818398e1, -0.558196781075],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.42752599]*2}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexane of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [2.43433265, 1.18137185, -4.24411947, 1.08655334e-1,
                2.87828538e-4, -2.51781047e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [2.16096570e-2, -4.58052979e-1, 1.63940974e-1, -2.55034034e-2,
                -0.247418231, -8.05544799e-3, -7.78926202e-2, -2.69044742e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = shortSpan, GERG, polt, starling, sun

    _surface = {"sigma": [0.210952, -0.158485], "exp": [1.0962, 1.05893]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [30.18, 0.03], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [222.31, 232.62, -36872, -25733],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.74172e1, 0.12897e1, -0.32544e1, -0.14609e1, 0.81765e-1],
        "t": [1., 1.5, 3.1, 5.3, 5.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.14686e3, -0.26585e3, 0.12200e3],
        "t": [0.75, 0.81, 0.88]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.13309, -5.0653, -11.602, -28.530, -51.731, -134.82],
        "t": [0.107, 0.553, 2.006, 4.46, 8.0, 16.]}

    visco0 = {"__name__": "Michailidou (2013)",
              "__doi__": {
                  "autor": "Michailidou, E.K., Assael, M.J., Huber, M.L., "
                           "Perkins, R.A.",
                  "title": "Reference Correlation of the Viscosity of "
                           "n-Hexane from the Triple Point to 600 K and up "
                           "to 100 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 42(3) (2013) 033104",
                  "doi": "10.1063/1.4818980"},

              "eq": 1, "omega": 1,

              "ek": 378.4, "sigma": 0.6334,
              "n_chapman": 0.021357,
              "collision": [0.1876, -0.4843, 0.04477],

              "Tref_virial": 378.4,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "special": "_vir"}

    def _vir(self, rho, T, fase):
        """Special residual term for Michailidou viscosity correlation"""
        Tr = T/507.82
        rhor = rho/233.182

        # Eq 8
        vir = 2.53402335/Tr - \
            9.724061002/(0.469437316+Tr+158.5571631*rhor**2) + \
            72.42916856*(1+rhor)/(10.60751253+8.628373915*Tr-6.61346441*rhor +
                                  rhor**2-2.212724566*rhor*Tr)
        vir *= rhor**(2/3)*Tr**0.5
        return vir

    visco1 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 507.82,
              "no": [16.9975, -54.2985, 48.0065],
              "to": [0, 0.25, 0.5],

              "a": [-6.63501e-5, -2.14252e-5, 0],
              "b": [1.64280e-4, -1.34908e-4, 0],
              "c": [7.25571e-5, -3.12153e-6, 0.0],
              "A": [1.45984e-9, -8.15150e-10, 0.0],
              "B": [2.59524e-8, 1.69362e-9, 0.0],
              "C": [-2.29226e-6, 1.18011e-6, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"__name__": "Assael (2013)",
               "__doi__": {
                   "autor": "Assael, M.J., Mylona, S.K., Tsiglifisi, Ch.A., "
                            "Huber, M.L., Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal "
                            "Conductivity of n-Hexane from the Triple Point "
                            "to 600 K and up to 500 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 42(1) (2013) 013106",
                   "doi": "10.1063/1.4793335"},

               "eq": 1,

               "Toref": Tc, "koref": 1e-3,
               "no": [6.6742, -23.7619, 72.0155, -18.3714],
               "to": [0, 1, 2, 3],

               "Tref_res": Tc, "rhoref_res": 233.182, "kref_res": 1,
               "nr": [-3.01408e-2, 1.67975e-1, -1.29739e-1, 3.82833e-2,
                      -3.70294e-3, 2.18208e-2, -1.00833e-1, 7.7418e-2,
                      -2.15945e-2, 2.12487e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.2364e-9,
               "gam0": 0.05803, "qd": 0.737e-9, "Tcref": 761.7}

    _thermal = thermo0,


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 46
        st = nC6(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.1802)
        self.assertEqual(round(st.P.MPa, 3), 10.221)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.6535)

        st2 = nC6(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 213.09)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.33219)

    def test_Michailidou(self):
        # Table 6, Pag 10
        self.assertEqual(round(nC6(T=250, rho=0).mu.muPas, 4), 5.2584)
        self.assertEqual(round(nC6(T=400, rho=0).mu.muPas, 4), 8.4150)
        self.assertEqual(round(nC6(T=550, rho=0).mu.muPas, 3), 11.443)
        self.assertEqual(round(nC6(T=250, rho=700).mu.muPas, 2), 528.20)
        self.assertEqual(round(nC6(T=400, rho=600).mu.muPas, 2), 177.62)
        self.assertEqual(round(nC6(T=550, rho=500).mu.muPas, 3), 95.002)

    def test_Assael(self):
        # Table 4, Pag 8
        self.assertEqual(round(nC6(T=250, rho=700).k.mWmK, 2), 137.62)
        self.assertEqual(round(nC6(T=400, rho=2).k.mWmK, 3), 23.558)
        self.assertEqual(round(nC6(T=400, rho=650).k.mWmK, 2), 129.28)
        self.assertEqual(round(nC6(T=510, rho=2).k.mWmK, 3), 36.772)
