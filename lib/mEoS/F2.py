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


from lib import unidades
from lib.meos import MEoS


class F2(MEoS):
    """Multiparameter equation of state for fluorine"""
    name = "fluorine"
    CASNumber = "7782-41-4"
    formula = "F2"
    synonym = ""
    _refPropName = "FLUORINE"
    _coolPropName = "Fluorine"
    rhoc = unidades.Density(592.864)
    Tc = unidades.Temperature(144.414)
    Pc = unidades.Pressure(5172.4, "kPa")
    M = 37.99681  # g/mol
    Tt = unidades.Temperature(53.4811)
    Tb = unidades.Temperature(85.0368)
    f_acent = 0.0449
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 208

    CP1 = {"ao": 3.5011231,
           "an": [-0.60936946e-4/144.414**4, 0.63196690e-3/144.414**3,
                  -0.74069617e-4/144.414**-2],
           "pow": [4, 3, -2],
           "ao_exp": [1.0127670], "exp": [1286.12]}

    CP2 = {"ao": 0.7593432/8.3143*37.997,
           "an": [0.2883653e-3/8.3143*37.997, -0.4192916e-5/8.3143*37.997,
                  0.2309778e-7/8.3143*37.997, -0.3291582e-10/8.3143*37.997],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    reuck = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for fluorine of de Reuck "
                    "(1990).",
        "__doi__": {"autor": "de Reuck, K.M.",
                    "title": "International thermodynamic tables of the fluid "
                             "state: Vol. 11 - fluorine",
                    "ref": "Pergamon Press, Oxford, 1990.",
                    "doi": ""},

        "R": 8.31448,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 300.0, "Pmax": 20000.0, "rhomax": 45.47,
        "Pmin": 0.23881, "rhomin": 44.917,

        "nr1": [0.151144749736e1, -0.298666288409e1, 0.329644905098e1,
                -0.298458624201e1, -0.228688966459e1, -0.109492193400e1,
                0.304775277572e1, 0.115689564208, -0.116100171627e1,
                0.295656394476, 0.711482542928e-1, -0.171363832155e-2,
                0.665317955515e-3],
        "d1": [1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 8, 9],
        "t1": [0, 0.5, 1.5, 2., 0.5, 1., 0.5, 2., 0.5, 1., 0., 0.5, 0],

        "nr2": [0.506026676251e1, -0.629268435440e1, 0.617784808739e1,
                -0.155366191788e1, -0.287170687343e1, 0.317214480494e1,
                -0.267969025215e1, 0.271865479252e1, -0.107191065039e1,
                0.126597342291e1, -0.706244695489, 0.268707888826,
                0.527251190274e-1, 0.544411481926e-1, 0.228949994105e-3,
                -0.547908264304e-9, -0.964273224950e-1, 0.368084486225e-3],
        "d2": [2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 7, 8, 12, 4, 6, 6],
        "t2": [1, 3, 4, 5, 1, 4, 5, 1, 3, 5, 4, 4, 1, 1, 5, 30, 20, 25],
        "c2": [2]*18,
        "gamma2": [1.07810258]*15+[2.15620515, 3.23430773, 3.23430773]}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for fluorine of Polt (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 300.0, "Pmax": 25000.0, "rhomax": 45.14,
        "Pmin": 0.25394, "rhomin": 44.89,

        "nr1": [0.862212325175e-2, 0.162286882091, -0.228707299586e-1,
                0.624951179331, -0.158918489879e1, 0.195171867807,
                -0.438453517535, 0.402200928405e-1, 0.319444405579e-1,
                0.161784325978e-1, 0.230132378392, 0.819206229044e-1,
                -0.173741828076, 0.137942204542e-1, -0.449971813506e-2,
                0.756554661780e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.862212325175e-2, -0.162286882091, 0.228707299586e-1,
                0.184612089745, -0.425779777811, 0.825656492996e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9225328]*6}

    eq = reuck, polt

    _surface = {"sigma": [0.03978], "exp": [1.218]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "Tmin": Tt, "Tmax": 300.0,
                "a1": [.000252, 249.975, -249.9750131], "exp1": [0, 2.1845, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.37061e1, -0.81517e2, 0.13743e3, -0.58617e2, -0.13528e1],
        "t": [1.0, 1.50, 1.61, 1.77, 7.3]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.21286, 0.44011e1, -0.53959e1, 0.41347e1, -0.97544],
        "t": [0.228, 0.58, 0.908, 1.24, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.6218, -55.363, 122.14, -230.92, -338.61, 432.18],
        "t": [0.454, 2.3, 2.9, 4.0, 6.0, 5.3]}
