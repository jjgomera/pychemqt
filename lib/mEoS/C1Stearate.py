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


class C1Stearate(MEoS):
    """Multiparameter equation of state for methyl stearate"""
    name = "methyl stearate"
    CASNumber = "112-61-8"
    formula = "C19H38O2"
    synonym = ""
    _refPropName = "MSTEARAT"
    _coolPropName = "MethylStearate"
    rhoc = unidades.Density(237.101584226)
    Tc = unidades.Temperature(775.0)
    Pc = unidades.Pressure(1239.0, "kPa")
    M = 298.50382  # g/mol
    Tt = unidades.Temperature(311.84)
    Tb = unidades.Temperature(629.56)
    f_acent = 1.02
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")

    f = 8.314472
    CP1 = {"an": [247.115/f], "pow": [-0.0916606],
           "ao_exp": [276.94/f, 408.997/f, 472.702/f],
           "exp": [556.17, 1311.85, 2825.71]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl estearate of Huber"
                    " et al. (2009).",
        "__doi__": {"autor": "Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, "
                             "L.S., Bruno, T.J.",
                    "title": "Model for the Thermodynamic Properties of a "
                             "Biodiesel Fuel",
                    "ref": "Energy Fuels, 23 (7) (2009) 3790–3797",
                    "doi": "10.1021/ef900159g"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 2.86,

        "nr1": [0.3959635e-1, 2.466654, -3.89595, -0.1167375, 0.4127229e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 1.25, 1.65, 0.8],

        "nr2": [-1.403734, -0.6465264, 1.934675, -1.608124, -0.1113813e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [3.1, 3.4, 2.3, 3.8, 1.2],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.125325, -0.7772671, -0.4183684],
        "d3": [1, 1, 3],
        "t3": [3.2, 3.8, 3.8],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = (huber, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.02313, 0.04567], "exp": [3.242, 1.163]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.14597e2, 0.13836e2, -0.14484e2, -0.51856e1, -0.27076e1],
        "t": [1.0, 1.5, 2.12, 4.7, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.11202e2, 0.78636e2, -0.12554e3, 0.72942e2, -0.11524e2],
        "t": [0.439, 0.59, 0.73, 0.9, 1.2]}
    _vapor_Density = {
        "eq": 2,
        "n": [-18.187, 81.619, -90.21, -5.2888e2, 1.1270e3, -8.4453e2],
        "t": [0.71, 1.3, 1.5, 6.0, 7.0, 8.0]}

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

              "ek": 615.42, "sigma": 0.8735, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.46654, -0.260069, 0.0354629], "psi_d": [0, 1, 2]}

    _viscosity = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # C1Stearate(T=697.5, rhom=1.79).mu.muPas, 4), 151.1512)
            C1Stearate(T=697.5, rhom=1.79).mu.muPas, 4), 151.1649)
