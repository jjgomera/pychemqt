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


class MM(MEoS):
    """Multiparameter equation of state for hexamethyldisiloxane"""
    name = "hexamethyldisiloxane"
    CASNumber = "107-46-0"
    formula = "C6H18OSi2"
    synonym = "MM"
    rhoc = unidades.Density(304.4043888253152)
    Tc = unidades.Temperature(518.69997204)
    Pc = unidades.Pressure(1939.39, "kPa")
    M = 162.37752  # g/mol
    Tt = unidades.Temperature(204.93)
    Tb = unidades.Temperature(373.401)
    f_acent = 0.418
    momentoDipolar = unidades.DipoleMoment(0.801, "Debye")
    id=1376

    Fi1 = {"ao_log": [1, 50.894],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [229.69732080664645, -86.53450336886623, -192.26652900000002,
                      18.654111840000002, -0.8140770995175003],
           "ao_exp": [], "titao": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of Colonna et al. (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},
        "__test__": """
            >>> st=MM(T=518.69997204, P=1939390)
            >>> print "%0.6f" % st.v
            0.003285
            """, # Table 16, Pag 202

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": 273.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 5.21,
        "Pmin": 0.00269, "rhomin": 5.2,

        "nr1": [1.01686012, -2.19713029, 0.75443188, -0.68003426, 0.19082162,
                0.10530133e-2],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.6284595, 0.30903042e-1, -0.83948727, -0.20262381,
                -0.35131597e-1, 0.25902341e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.86671e1, 0.11649e2, -0.11484e2, -0.53256e1],
        "exp": [1.0, 1.5, 1.65, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14533e2, -0.49804e2, 0.83748e2, -0.70321e2, 0.24283e2],
        "exp": [0.584, 0.8, 1.02, 1.26, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.35719e1, -0.14740e3, 0.40699e3, -0.69676e3, 0.12541e4, -0.91199e3],
        "exp": [0.373, 2.15, 2.6, 3.3, 4.2, 4.6]}
