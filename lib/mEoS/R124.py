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
from lib.mEoS import C3


class R124(MEoS):
    """Multiparameter equation of state for R124"""
    name = "1-chloro-1,2,2,2-tetrafluoroethane"
    CASNumber = "2837-89-0"
    formula = "CHClFCF3"
    synonym = "R124"
    _refPropName = "R124"
    _coolPropName = "R124"
    rhoc = unidades.Density(560.)
    Tc = unidades.Temperature(395.425)
    Pc = unidades.Pressure(3624.295, "kPa")
    M = 136.475  # g/mol
    Tt = unidades.Temperature(74.)
    Tb = unidades.Temperature(261.187)
    f_acent = 0.28810
    momentoDipolar = unidades.DipoleMoment(1.469, "Debye")
    # id = 1629

    CP1 = {"ao": 3.175638,
           "an": [14.77947/Tc, -5.2420986/Tc**2, 1.3381596/Tc**3],
           "pow": [1, 2, 3]}

    CP2 = {"ao": 3.20532538,
           "an": [13.4403357/395.62, -2.32192933/395.62**2,
                  -0.422826803/395.62**3],
           "pow": [1, 2, 3]}

    vries = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-124 of de Vries (1995)",
        "__doi__": {"autor": "de Vries, B., Tillner-Roth, R., Baehr, H.D.",
                    "title": "Thermodynamic Properties of HCFC 124,",
                    "ref": "19th International Congress of Refrigeration, The "
                           "Hague, The Netherlands, IIR, IVa:582-589, 1995",
                    "doi": ""},

        "R": 8.314471,
        "M": 136.475, "rhoc": 4.1033156,

        "cp": CP1,
        "ref": "IIR",

        "Tmin": 120.0, "Tmax": 470.0, "Pmax": 40000.0, "rhomax": 20,

        "nr1": [-0.1262962e-1, 0.2168373e1, -0.3330033e1, 0.1610361,
                -0.9666145e-4, 0.1191310e-1, -0.2880217e-2, 0.1681346e-2,
                0.1594968e-4],
        "d1": [1, 1, 1, 2, 2, 3, 5, 6, 8],
        "t1": [2, 0.5, 1, 0.5, 2.5, -1, 1, 0, -0.5],

        "nr2": [0.1289674, 0.1182213e-4, -0.4713997, -0.2412873, 0.6868066,
                -0.8621095e-1, 0.4728645e-5, 0.1487933e-1, -0.3001338e-1,
                0.1849606e-2, 0.4126073e-3],
        "d2": [2, 12, 1, 1, 1, 1, 15, 3, 3, 4, 9],
        "t2": [1.5, 1, 2.5, -0.25, 1, 5, 2, 15, 20, 15, 45],
        "c2": [1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4],
        "gamma2": [1]*11}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-124 of Younglove and "
                    "McLinden (1994)",
        "__doi__": {"autor": "McLinden, M.O., Younglove, B.A., Sandarusi, J.",
                    "title": "Measurement of the PVT properties and "
                             "formulation of an equation of state for "
                             "refrigerant 124 (1-chloro-1,2,2,2-"
                             "tetrafluoroethane)",
                    "ref": "1994. (unpublished manuscript)",
                    "doi": ""},

        "R": 8.314471,
        "M": 136.4762, "Tc": 395.62, "Pc": 3637., "rhoc": 4.101527,

        "cp": CP2,
        "ref": "IIR",

        "Tmin": 120.0, "Tmax": 475.0, "Pmax": 36000.0, "rhomax": 13.98,

        "b": [None, -0.195111839846e-1, 0.299978502039e1, -0.845849168162e2,
              0.146720754658e5, -0.232549336572e7, 0.938866046628e-3,
              -0.425069993257e1, 0.304859131600e4, 0.221314829910e7,
              -0.601971995213e-4, 0.100335188373e1, -0.468461812962e3,
              -0.927654315163e-2, -0.125426962519e2, -0.228534445089e4,
              0.168197835599e1, -0.537322295315e-1, 0.157915168095e2,
              -0.550297175283, -0.244349954189e7, -0.625153016263e8,
              -0.156149231820e6, 0.344268154495e10, -0.289212955106e4,
              0.108351996828e6, -0.404809912845e2, -0.220587292481e7,
              -0.564677367857, 0.175581172016e3, -0.762146322899e-3,
              -0.210617958917e1, 0.319236066221e2]}

    eq = vries, younglove

    _surface = {"sigma": [0.05175], "exp": [1.197]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.75146e1, 0.37481e1, -0.30124e1, -0.37808e1, -0.53114],
        "t": [1.0, 1.5, 1.68, 3.8, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.19127e1, 0.67778, -0.35129e-1, 0.30407, 0.69503e-1],
        "t": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.8551, -6.385, -17.616, -37.828, -23.785, -134.59],
        "t": [0.388, 1.17, 3.0, 6.0, 8.0, 15.0]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",

              "ek": 275.8, "sigma": 0.5501, "omega": 6,
              "n_chapman": 0.026692,

              "psi": [1.04253, 1.38528e-3], "psi_d": [0, 1],
              "fint": [1.1769e-3, 6.78397e-7], "fint_t": [0, 1],
              "chi": [1.0898, -1.54229e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_huber(self):
        """Table 4 from lib.coolprop reference, pag 2504"""
        st = R124(T=350, x=0)
        self.assertEqual(round(st.mu.muPas, 3), 138.059)
        self.assertEqual(round(st.k.mWmK, 3), 52.794)
