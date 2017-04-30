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


class Kr(MEoS):
    """Multiparameter equation of state for krypton"""
    name = "krypton"
    CASNumber = "7439-90-9"
    formula = "Kr"
    synonym = "R-784"
    rhoc = unidades.Density(909.2083)
    Tc = unidades.Temperature(209.48)
    Pc = unidades.Pressure(5525.0, "kPa")
    M = 83.798  # g/mol
    Tt = unidades.Temperature(115.775)
    Tb = unidades.Temperature(119.73)
    f_acent = -0.000894
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 971

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-3.7506412806, 3.7798018435],
           "ao_exp": [], "titao": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for krypton of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},
        "__test__": """
            >>> st=Kr(T=211, rho=10*83.798)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            211 10 5741.445 6700.326 36.936 27.390 1667.678 137.838
            """, # Table 10, Pag 842

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 200000.0, "rhomax": 33.42,
        "Pmin": 73.5, "rhomin": 29.2,

        "nr1": [0.83561, -2.3725, 0.54567, 0.014361, 0.066502, 0.00019310],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.16818, -0.033133, -0.15008, -0.022897, -0.021454, 0.0069397],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for krypton of Polt et al. (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 780.0, "Pmax": 375000.0, "rhomax": 33.55,
        "Pmin": 73.476, "rhomin": 29.249,

        "nr1": [-0.402218741560, 0.679250544381, -0.1878869802860,
                0.603399982935, -0.177297564389e1, 0.581208430222,
                -0.733585469788, 0.164651929067, -0.319923148922e-1,
                0.333278228743, 0.219652478083e-1, 0.751994891628e-1,
                -0.212109737251, -0.645185506524e-2, 0.409175610200e-1,
                0.169416098754e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.402218741560, -0.679250544381, 0.187886980286,
                0.108265263587, -0.137102675805, -0.110549803007],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2, 2, 2, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.0447], "exp": [1.245]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [6.273], "expt1": [0], "expd1": [1],
                   "a2": [6.485, 13.48, -82.51, -170.4],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.7, 2.7]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 101.325,
                "Tmin": Tt, "Tmax": 800.0,
                "a1": [-2345.757, 1.080476685], "exp1": [0, 1.6169841],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 73.197,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-11.5616], "exp2": [1],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.59697e1, 0.12673e1, -0.95609, -0.35630e2, 0.56884e2],
        "exp": [1.0, 1.5, 2.95, 9.3, 10.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.20593e2, -0.65490e2, 0.94407e2, -0.69678e2, 0.22810e2],
        "exp": [0.62, 0.84, 1.07, 1.34, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.64163e1, 0.89956e1, -0.10216e2, -0.13477e2, -0.21152e3, 0.21375e3],
        "exp": [0.525, 0.77, 1.04, 3.2, 8.3, 9.0]}
