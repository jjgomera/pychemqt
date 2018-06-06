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


class nC10(MEoS):
    """Multiparameter equation of state for n-decane"""
    name = "decane"
    CASNumber = "124-18-5"
    formula = "CH3-(CH2)8-CH3"
    synonym = ""
    rhoc = unidades.Density(233.342)
    Tc = unidades.Temperature(617.7)
    Pc = unidades.Pressure(2103.0, "kPa")
    M = 142.28168  # g/mol
    Tt = unidades.Temperature(243.5)
    Tb = unidades.Temperature(447.27)
    f_acent = 0.4884
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 14

    Fi1 = {"ao_log": [1, 18.109],
           "pow": [0, 1],
           "ao_pow": [13.9361966549, -10.5265128286],
           "ao_exp": [25.685, 28.233, 12.417, 10.035],
           "titao": [1193/Tc, 2140/Tc, 4763/Tc, 10862/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.870791919, -108.858547525],
           "ao_exp": [], "titao": [],
           "ao_hyp": [21.0069, 43.4931, 58.3657, 0],
           "hyp": [0.267034159, 1.353835195, 2.833479035, 0]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for decane of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 800000.0, "rhomax": 5.41,
        "Pmin": 0.0014, "rhomin": 5.41,

        "nr1": [1.0461, -2.4807, 0.74372, -0.52579, 0.15315, 0.00032865],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.84178, 0.055424, -0.73555, -0.18507, -0.020775, 0.012335],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

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

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 800000.0, "rhomax": 5.41,
        "Pmin": 0.0014, "rhomin": 5.41,

        "nr1": [0.10461e1, -0.24807e1, 0.74372, -0.52579, 0.15315, 0.32865e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [.84178, .55424e-1, -.73555, -.18507, -.20775e-1, .12335e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon, GERG

    _surface = {"sigma": [0.05473], "exp": [1.29]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [49.32, 0.05], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [220.15, -316.3, -88358, 53511],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.87738e1, 0.40864e1, -0.40775e1, -0.64910e1, 0.15598e1],
        "exp": [1.0, 1.5, 1.93, 4.14, 4.7]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.92435e1, -0.16288e2, 0.20445e2, -0.17624e2, 0.73796e1],
        "exp": [0.535, 0.74, 1.0, 1.28, 1.57]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-5.0378, -3.4694, -.15906e2, -0.82894e2, 0.29336e2, -0.10985e3],
        "exp": [0.4985, 1.33, 2.43, 5.44, 5.8, 11.0]}

    visco0 = {"__name__": "Huber (2004)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A. Xiang, H.W.",
                  "title": "Viscosity correlations for minor constituent "
                           "fluids in natural gas: n-octane, n-nonane and "
                           "n-decane",
                  "ref": "Fluid Phase Equilibria 224 (2004) 263-270",
                  "doi": "10.1016/j.fluid.2004.07.012"},

              "eq": 1, "omega": 1,

              "ek": 490.51, "sigma": 0.686,
              "n_chapman": 0.021357,
              "collision": [0.343267, -0.460514],

              "Tref_virial": 490.51,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 617.7, "rhoref_res": 1.64*M, "muref_res": 1000,
              "nr": [-.0402094, 0.0404435, -0.0142063],
              "tr": [1, 2, 2],
              "dr": [2, 2, 3],

              "CPf": 453.387,
              "CPg1": 2.55105,
              "CPgi": [1.71465/2.55105],
              "CPti": [-0.5]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Huber (2005)",
               "__doi__": {"autor": "Huber, M.L. and Perkins, R.A.",
                           "title": "Thermal conductivity correlations for minor constituent fluids in natural gas: n-octane, n-nonane and n-decane",
                           "ref": "Fluid Phase Equilibria 227 (2005) 47-55",
                           "doi": "10.1016/j.fluid.2004.10.031"},
               "__test__": """
                   >>> st=nC10(T=300, rhom=5.1504)
                   >>> print "%0.2f" % st.k
                   132.80
                   """,  # Section 3.3 pag 54

               "Tref": 617.7, "kref": 1,
               "no": [0.105543e-1, -0.514530e-1, 0.118979, -0.372442e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 617.7, "rhorefb": 1.64, "krefb": 1,
               "nb": [-.294394e-1, .150509e-1, .499245e-1, 0., -.142700e-1,
                      -0.138857e-1, 0.150828e-2, 0.433326e-2, 0.0, 0.0],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 7.086368e-10, "Tcref": 926.55}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = nC10(T=619, rhom=1)
        self.assertEqual(round(st.P.kPa, 3), 2071.025)
        self.assertEqual(round(st.hM.kJkmol, 3), 89742.553)
        self.assertEqual(round(st.sM.kJkmolK, 3), 164.787)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 437.033)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 1043.328)
        self.assertEqual(round(st.w, 3), 74.576)

    def test_viscoHuber(self):
        # Section 3.3 pag 269
        self.assertEqual(round(nC10(T=300, rhom=5.1504).mu.muPas, 2), 926.44)
