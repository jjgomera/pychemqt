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


from lib.meos import MEoS
from lib import unidades


class RE245fa2(MEoS):
    """Multiparameter equation of state for RE245fa2"""
    name = "2,2,2-trifluoroethyl-difluoromethyl-ether"
    CASNumber = "1885-48-9"
    formula = "CHF2OCH2CF3"
    synonym = "HFE-245fa2"
    _refPropName = "RE245FA2"
    _coolPropName = ""
    rhoc = unidades.Density(515.001169364688)
    Tc = unidades.Temperature(444.88)
    Pc = unidades.Pressure(3433., "kPa")
    M = 150.047336  # g/mol
    Tt = unidades.Temperature(250)
    Tb = unidades.Temperature(302.4)
    f_acent = 0.387
    momentoDipolar = unidades.DipoleMoment(1.631, "Debye")
    # id = 1874

    CP1 = {"ao": 5.259865,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [12.12843, 13.25677, 0.521867, 0],
           "hyp": [486, 1762, 7631, 0]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE245fa2 of Zhou (2012)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., Mahmoud, A.M.",
                    "title": "Equations of state for RE245cb2, RE347mcc, "
                             "RE245fa2 and R1216",
                    "ref": "Preliminary equation",
                    "doi":  ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 400000.0, "rhomax": 10.02,
        "Pmin": 8.272, "rhomin": 10.,

        "nr1": [0.47771378e-1, 0.15745383e1, -0.24763491e1, -0.49414564,
                0.19380498],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.32, 0.91, 1.265, 0.4266],

        "nr2": [-0.97863158, -0.42660297, 0.85352583, -0.53380114,
                -0.29780036e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.24, 1.64, 1.65, 3.28, 0.855],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.97659111, -0.33121365, -0.14122591, -0.15312295e2],
        "d3": [1, 1, 3, 3],
        "t3": [1.227, 3.0, 4.3, 2.5],
        "alfa3": [1.005, 1.515, 1.156, 17.7],
        "beta3": [2, 3.42, 1.37, 471],
        "gamma3": [1.084, 0.72, 0.49, 1.152],
        "epsilon3": [0.723, 0.9488, 0.818, 0.891]}

    eq = zhou,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-8.9235, 10.527, -23.058, 30.291, -20.913, -26.745],
        "exp": [1, 1.5, 1.9, 2.4, 2.9, 3.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [1.2479, 5.5732, -12.26, 13.964, -6.0384],
        "exp": [0.34, 0.75, 1.2, 1.7, 2.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.667, -5.8238, -26.927, 21.574, -65.645],
        "exp": [0.28, 0.66, 2.6, 3.5, 5.2]}
