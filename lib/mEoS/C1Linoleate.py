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
from lib.mEoS import C3


class C1Linoleate(MEoS):
    """Multiparameter equation of state for methyl linoleate"""
    name = "methyl linoleate"
    CASNumber = "112-63-0"
    formula = "C19H34O2"
    synonym = ""
    _refPropName = "MLINOLEA"
    _coolPropName = "MethylLinoleate"
    rhoc = unidades.Density(238.051213304)
    Tc = unidades.Temperature(799.0)
    Pc = unidades.Pressure(1341.0, "kPa")
    M = 294.47206  # g/mol
    Tt = unidades.Temperature(238.1)
    Tb = unidades.Temperature(628.84)
    f_acent = 0.805
    momentoDipolar = unidades.DipoleMoment(1.79, "Debye")

    f = 8.314472
    CP1 = {"an": [190.986/f],
           "pow": [0.020213],
           "ao_exp": [437.371/f, 287.222/f, 321.956/f],
           "exp": [3052.11, 746.631, 1624.33]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl linoleate of Huber"
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

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.16,

        "nr1": [0.3183187e-1, 0.1927286e1, -0.3685053e1, 0.8449312e-1],
        "d1": [4, 1, 1, 3],
        "t1": [1, 0.2, 1.2, 1.0],

        "nr2": [-0.9766643, -0.4323178, 2.00047, -1.75203, -0.01726895],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.2, 2.5, 1.8, 1.92, 1.47],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2116515e1, -0.7884271, -0.3811699],
        "d3": [1, 1, 3],
        "t3": [1.7, 2.3, 2.1],
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
        "sigma": [0.072487], "exp": [1.9014]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.10946e2, 0.48849e1, -0.46773e1, -0.80201e1, -0.89572e1],
        "t": [1.0, 1.5, 2.22, 3.6, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.22705e3, -0.66763e3, 0.72323e3, -0.49244e3, 0.21391e3],
        "t": [0.83, 0.98, 1.17, 1.5, 1.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-8.588, 14.766, -24.195, -374.74, 326.89, -191.25],
        "t": [0.568, 1.08, 1.4, 4.8, 5.0, 9.0]}

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

              "ek": 634.48, "sigma": 0.8684, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.49489, -0.25538, 0.0306593], "psi_d": [0, 1, 2]}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Perkins (2010)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivities of Biodiesel Constituent Fluids: "
                            "Methyl Oleate and Methyl Linoleate",
                   "ref": "Energy Fuels 25(5) (2011) 2383-2388",
                   "doi": "10.1021/ef200417x"},

               "eq": 1,

               "Toref": 799.0, "koref": 1,
               "no": [-1.09042e-4, 2.40543e-3, 0.0407364, -0.0105928],
               "to": [0, 1, 2, 3],

               "Tref_res": 799.0, "rhoref_res": 238.05, "kref_res": 1.,
               "nr": [-0.0713126, 0.0989415, 0.0466421, -0.065785, -0.00557406,
                      0.0128922],
               "tr": [0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1198.5}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # C1Linoleate(T=719.1, rhom=1.861).mu.muPas, 3), 148.324)
            C1Linoleate(T=719.1, rhom=1.861).mu.muPas, 3), 148.334)

    def test_Perkins(self):
        """Table 3, Pag 2386"""
        st = C1Linoleate(T=450, P=1e2)
        self.assertEqual(round(st.rho, 8), 0.00787223)
        self.assertEqual(round(st.k, 7), 0.0122743)

        st = C1Linoleate(T=450, P=1e6)
        self.assertEqual(round(st.rho, 3), 778.176)
        self.assertEqual(round(st.k, 6), 0.122742)

        st = C1Linoleate(T=450, P=2e7)
        self.assertEqual(round(st.rho, 3), 799.160)
        self.assertEqual(round(st.k, 6), 0.131867)
