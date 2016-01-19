#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


class Propylene(MEoS):
    """Multiparameter equation of state for propylene"""
    name = "propylene"
    CASNumber = "115-07-1"
    formula = "CH2=CH-CH3"
    synonym = "R-1270"
    rhoc = unidades.Density(230.08)
    Tc = unidades.Temperature(364.211)
    Pc = unidades.Pressure(4555.0, "kPa")
    M = 42.07974  # g/mol
    Tt = unidades.Temperature(87.953)
    Tb = unidades.Temperature(225.531)
    f_acent = 0.146
    momentoDipolar = unidades.DipoleMoment(0.366, "Debye")
    id = 23

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-5.1823279651, 4.3639902765],
           "ao_exp": [1.544, 4.013, 8.923, 6.020],
           "titao": [324/Tc, 973/Tc, 1932/Tc, 4317/Tc]}

    Fi2 = {"ao_log": [1, 3.07317535],
           "pow": [0, 1],
           "ao_pow": [9.48120502357782, -4.47976952867319],
           "ao_exp": [1.7018443, 3.61342025, 8.83689058, 6.27183616],
           "titao": [1.01164134251849, 2.75278088800174, 5.16557061703243,
                     11.68984352477]}

    CP1 = {"ao": 0.65591381,
           "an": [0.44359641e-1, -.36650786e-4, 0.16822223e-7, -.32651013e-11,
                  0.33747826e4],
           "pow": [1, 2, 3, 4, -2],
           "ao_exp": [-4.7032420], "exp": [615.8],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylene of Lemmon et al. (2013).",
        "__doi__": {"autor": "Lemmon, E.W., Overhoff, U., McLinden, M.O., Wagner, W.",
                    "title": "A reference equation of state for the thermodynamic properties of propene for temperatures from the melting line to 575 K and pressures up to 1000 MPa",
                    "ref": "to be submitted to J. Phys. Chem. Ref. Data",
                    "doi": ""},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 1000000.0, "rhomax": 23.1,
        "Pmin": 0.00000075, "rhomin": 18.255,

        "nr1": [0.4341002e-1, 0.1136592e1, -0.8528611, 0.5216669, -0.1382953e1,
                0.1214347],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0, 0.205, 0.56, 0.676, 1.0, 0.5],

        "nr2": [-0.5984662, -0.1391883e1, -0.1008434e1, 0.1961249, -0.3606930,
                -0.2407175e-2],
        "d2": [1, 1, 3, 2, 2, 8],
        "t2": [1.0, 1.94, 2.0, 1.0, 2.66, 0.83],
        "c2": [1, 2, 2, 1, 2, 1],
        "gamma2": [1]*6,

        "nr3": [0.7432121, 0.1475162, -0.2503391e-1, -0.2734409, 0.6378889e-2,
                0.1502940e-1, -0.3162971e-1, -0.4107194e-1, -0.1190241e1],
        "d3": [1, 1, 2, 3, 3, 2, 1, 2, 3],
        "t3": [1.6, 2.5, 3.0, 2.5, 2.72, 4.0, 4.0, 1.0, 4.0],
        "alfa3": [1.07, 0.66, 1.2, 1.12, 1.47, 1.93, 3.3, 15.4, 6],
        "beta3": [0.77, 0.83, 0.607, 0.4, 0.66, 0.07, 3.1, 387, 41],
        "gamma3": [1.21, 1.08, 0.83, 0.56, 1.22, 1.81, 1.54, 1.12, 1.4],
        "epsilon3": [0.78, 0.82, 1.94, 0.69, 1.96, 1.3, 0.38, 0.91, 0.7]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylene of Overhoff (2006).",
        "__doi__": {"autor": "Overhoff, U.",
                    "title": "Development of a new equation of state for the fluid region of propene for temperatures from the melting line to 575 K with pressures to 1000 MPa as well as software for the computation of thermodynamic properties of fluids",
                    "ref": "Ph.D. Dissertation, Ruhr University, Bochum, Germany, 2006.",
                    "doi": ""},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 1000000.0, "rhomax": 23.4,
        "Pmin": 0.00000074, "rhomin": 18.26,

        "nr1": [0.11167427541961e1, -0.76114879497376, -0.18654354344883e1,
                0.41500701892893e-1, 0.10706545719025e-1, 0.17481482892991e-1],
        "d1": [1, 1, 1, 3, 4, 4],
        "t1": [0.125, 0.625, 1.25, 0, 0.25, 1.25],

        "nr2": [0.56509607629258, 0.99156795771235, -0.16341922173416,
                -0.37037920319844e-1, -0.80058345775777e-1, 0.17004662808796,
                0.81351262137108e-1, -0.23817885171378, 0.12962562859214e-1,
                0.22577442976798e2, -0.43611886043491e2, 0.21944325628071e2,
                -0.66234078215924, -0.22258580712469e1, 0.29538388307646e1,
                -0.10257185828694e1, 0.20521625234481e-1, -0.36462809205891e-1,
                0.17625833164005e-1],
        "d2": [2, 3, 3, 3, 4, 4, 5, 5, 6, 1, 1, 1, 1, 2, 2, 2, 5, 6, 1],
        "t2": [2.25, 1.25, 2.125, 2.75, 0.125, 2, 1.125, 1.5, 1.375, 3.5, 3.75,
               4, 5, 3, 3.5, 4.5, 4.75, 3.25, 3, ],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3],
        "gamma2": [1]*19,

        "nr3": [0.31819374579431, -0.32648950998998, -0.37684374593786e2,
                0.72265437094447e2, -0.34814669335983e2, -0.39854778355193e1,
                0.37313453915501],
        "d3": [2, 2, 1, 1, 1, 2, 2],
        "t3": [3, 4, 2, 3, 4, 1, 1],
        "alfa3": [10, 10, 11, 11, 11, 25, 30],
        "beta3": [150, 150, 225, 225, 225, 300, 350],
        "gamma3": [1.13, 1.13, 1.19, 1.19, 1.19, 1.19, 1.19],
        "epsilon3": [0.85, 0.85, 1, 1, 1, 1, 1]}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylene of Angus et al. (1980)",
        "__doi__": {"autor": "Angus, S., Armstrong, B., and de Reuck, K.M.",
                    "title": "International Thermodynamic Tables of the Fluid State-7 Propylene",
                    "ref": "International Union of Pure and Applied Chemistry, Pergamon Press, Oxford, 1980.",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 100.0, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 9.73,
        "Pmin": 0.48475e-4, "rhomin": 17.938,

        "nr1": [0.631922681460, 0.102655250604, -0.70798923e-2, 0.18624829,
                -0.1292611017e1, -0.5410160974e-1, 0.5069017035, -0.10606146125e1,
                0.763136083, -0.850733053e-1, 0.438262575, 0.2316495716e-1,
                0.25503741325e-1, -0.57327581, -0.1141334722e-1, 0.2502895522,
                -0.468392547833e-1, 0.325228355714e-2],
        "d1": [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 6, 7],
        "t1": [3, 4, 5, 1, 2, 3, 0, 1, 2, 2, 3, 0, 1, 3, -1, 3, 3, 3],

        "nr2": [-0.63192268146, -0.102655250604, 0.70798923e-2, -0.63192268146,
                -0.102655250604, -0.11049992895, -0.31596134073, -0.51327625302e-1,
                -0.4918627871e-1, -0.17109208434e-1, -0.1492467645e-1,
                -0.42773021085e-2, -0.8554604217e-3, -0.14257673695e-3],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 6, 6, 8, 8, 10, 12],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 3, 4, 3, 4, 4, 4],
        "c2": [2]*14,
        "gamma2": [1]*14}

    eq = helmholtz1, helmholtz2, helmholtz3

    _surface = {"sigma": [0.05268], "exp": [1.186]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.48475e-4,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [-6593000000, 6593000001], "exp1": [0, 2.821],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.75625, 2.02700, -1.35883, -2.74671, -0.936445],
        "exp": [1.0, 1.5, 1.9, 4.3, 15.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.405430, 2.02481, 0.304022, 0.179159],
        "exp": [0.195, 0.47, 2.25, 8.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.59841, -4.73840, -10.8886, -31.0312, -56.9431, -143.544],
        "exp": [0.309, 0.853, 2.37, 5.2, 10., 20.]}
