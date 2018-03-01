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

from lib.meos import MEoS
from lib import unidades


class Toluene(MEoS):
    """Multiparameter equation of state for toluene
    """
    name = "toluene"
    CASNumber = "108-88-3"
    formula = "C6H5-CH3"
    synonym = ""
    rhoc = unidades.Density(291.98665)
    Tc = unidades.Temperature(591.75)
    Pc = unidades.Pressure(4126.3, "kPa")
    M = 92.13842  # g/mol
    Tt = unidades.Temperature(178.0)
    Tb = unidades.Temperature(383.75)
    f_acent = 0.2657
    momentoDipolar = unidades.DipoleMoment(0.36, "Debye")
    id = 41

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [3.5241174832, 1.1360823464],
           "ao_exp": [1.6994, 8.0577, 17.059, 8.4567, 8.6423],
           "titao": [190/Tc, 797/Tc, 1619/Tc, 3072/Tc, 7915/Tc]}

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [1.6994, 8.0577, 17.059, 8.4567, 8.6423],
           "exp": [190, 797, 1619, 3072, 7915],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": -0.321892/8.3143*92.142,
           "an": [0.579338e-2/8.3143*92.142, -0.348446e-5/8.3143*92.142,
                  0.143577e-8/8.3143*92.142, -0.71935e-12/8.3143*92.142],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for toluene of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 500000.0, "rhomax": 10.581,
        "Pmin": 0.000039, "rhomin": 10.58,

        "nr1": [0.96464, -2.7855, 0.86712, -0.18860, 0.11804, 0.00025181],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.57196, -0.029287, -0.43351, -0.12540, -0.028207, 0.014076],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for toluene of Polt et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": 298.0, "Tmax": 673.0, "Pmax": 25000.0, "rhomax": 9.7242,
        "Pmin": 3.774, "rhomin": 9.3606,

        "nr1": [-0.343905499875, 0.737562743137, -0.158601557810,
                0.113243121503e1, -0.253681929563e1, 0.104584338973e1,
                -0.115732119380e1, 0.176205273278, -0.242942016719,
                0.398925293195, 0.193881828889, 0.199426230143, -0.306598708746,
                -0.114697533947e-1, 0.230068676459e-1, 0.658341220591e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.343905499875, -0.737562743137, 0.15860155781, 0.40707928397,
                -0.68140614165, 0.110425925004],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.841]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for toluene of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.34060172, 1.01624262, -3.27810202, 9.69209624e-2,
                2.61950176e-4, -1.58891991e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [6.28559812e-2, -8.42364946e-2, 4.49701117e-1, -1.08658876e-2,
                -3.83733669e-1, 2.21127543e-2, -9.54658223e-2, -1.77905259e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = lemmon, polt, sun

    _surface = {"sigma": [0.06897], "exp": [1.291]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.45201, 2.03681, -1.43777, -3.51652, -1.75818],
        "exp": [1.0, 1.5, 2.13, 4.0, 12.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [14.0531, -32.5072, 35.1091, -16.0694, 2.38699],
        "exp": [0.54, 0.72, 0.93, 1.2, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-2.97587, -5.34939, -19.1781, -24.0058, -32.4034, -140.645],
        "exp": [0.425, 1.06, 3.0, 6.3, 7.0, 15.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Lemmon (2010",
              "__doc__": """Lemmon, E.W. and Laesecke, A., 2010. Unpublished preliminary equation for the viscosity of toluene.""",
              "ek": 469.90, "sigma": 0.5507,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 9.598876,

              "Tref_res": 591.75, "rhoref_res": 3.169*M,
              "n_poly": [0.157560701809e2, 0.658234203776e2, -0.909162962259e2,
                         -0.806740654754e2, 0.395093273404e1, 0.867277691823e-1,
                         -0.928414042924e-2, 0.982264892850e-5,
                         -0.785434913708e-3, 0.169683455336e-7],
              "t_poly": [-0.2843, -2.4238, -2.7667, -3.0019, -3.2869, -6.0789,
                         -6.1564, -6.8541, -5.5123, -4.1175],
              "d_poly": [1, 2, 2, 4, 6, 9, 11, 12, 17, 19],
              "g_poly": [0, 0, 1, 1, 1, 1, 1, 0, 1, 0],
              "c_poly": [0, 0, 1, 1, 2, 1, 1, 0, 2, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Lemmon (2010)",
               "__doc__": """Lemmon, E.W. and Laesecke, A., 2010. Unpublished preliminary equation for the thermal conductivity of toluene.""",

               "Tref": 591.75, "kref": 1e-3,
               "no": [28.96745197, -167.24996945, 180.04690463],
               "co": [1.20532335, 1.58866032, 1.71267964],

               "Trefb": 591.75, "rhorefb": 3.169, "krefb": 1e-3,
               "nb": [-0.318905053658e1, 0.258544682121e2, -0.263059677817e2,
                      -0.691196173614, 0.542428651638e-1, -0.326501347819],
               "tb": [-0.53316, -0.27224, -0.09974, -5.53274, -6.84315, -0.39659],
               "db": [4, 3, 5, 7, 8, 3],
               "cb": [0, 0, 1, 2, 2, 2],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.33442441e-9, "gam0": 0.55e-1, "qd": 0.71763799e-9,
               "Tcref": 1183.5}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = Toluene(T=593, rhom=3)
        self.assertEqual(round(st.P.kPa, 3), 4186.620)
        self.assertEqual(round(st.hM.kJkmol, 3), 52937.550)
        self.assertEqual(round(st.sM.kJkmolK, 3), 105.422)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 214.488)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 7705.724)
        self.assertEqual(round(st.w, 3), 89.464)
