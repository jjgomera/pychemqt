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


class CO(MEoS):
    """Multiparameter equation of state for carbon monoxide"""
    name = "carbon monoxide"
    CASNumber = "630-08-0"
    formula = "CO"
    synonym = ""
    rhoc = unidades.Density(303.909585)
    Tc = unidades.Temperature(132.86)
    Pc = unidades.Pressure(3494.0, "kPa")
    M = 28.0101  # g/mol
    Tt = unidades.Temperature(68.16)
    Tb = unidades.Temperature(81.64)
    f_acent = 0.0497
    momentoDipolar = unidades.DipoleMoment(0.1, "Debye")
    id = 48

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1.5],
           "ao_pow": [-3.3728318564, 3.3683460039, 9.111274701235156e-5],
           "ao_exp": [1.0128],
           "titao": [3089/Tc]}

    Fi2 = {"ao_log": [1, 2.50055],
           "pow": [0, 1],
           "ao_pow": [10.813340744, -19.834733959],
           "ao_exp": [], "titao": [],
           "ao_hyp": [1.02865, 0.00493, 0, 0],
           "hyp": [11.6698028, 5.302762306, 0, 0]}

    CP3 = {"ao": 0.36028218e1,
           "an": [-0.20871594e5, 0.89208708e3, -0.14157993e2, -0.34021345e-3,
                  0.44616091e-6, -0.15154703e-9],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [0.90426143], "exp": [30000],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbon monoxide of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},
        "__test__": """
            >>> st=CO(T=134, rho=10*28.0101)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            134 10 3668.867 4838.507 41.601 38.702 1642.142 168.632
            """, # Table 10, Pag 842

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500., "Pmax": 100000.0, "rhomax": 33.84,
        "Pmin": 15.45, "rhomin": 30.33,

        "nr1":  [0.90554, -2.4515, 0.53149, 0.024173, 0.072156, 0.00018818],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.19405, -0.043268, -0.12778, -0.027896, -0.034154, 0.016329],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 15.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for carbon monoxide of McCarty (1989)",
        "__doi__": {"autor": "McCarty, R.D.",
                    "title": "Correlations for the Thermophysical Properties of Carbon Monoxide",
                    "ref": "National Institute of Standards and Technology, Boulder, CO, 1989.",
                    "doi":  ""},

        "R": 8.31434,
        "cp": CP3,

        "Tmin": Tt, "Tmax": 1000., "Pmax": 30000.0, "rhomax": 30.25,
        "Pmin": 15.423, "rhomin": 30.249,

        "b": [None, 0.8845582109949e-2, -0.2236741566840, 0.1742275796442e1,
              -0.2169146998363e3, 0.1721504267082e4, -0.3990514770703e-4,
              0.1036880040451, -0.3376308165071e2, 0.2061895161095e5,
              0.2993711656350e-5, 0.1856003597097e-2, -0.2114419664527,
              -0.2436986935194e-5, -0.1858029609177e-2, -0.1734563867767e1,
              0.1509970839260e-3, -0.2282721433205e-5, 0.2202780295674e-2,
              -0.3313357789163e-4, -0.1473412120276e5, -0.3141136651147e6,
              -0.1451168999234e3, 0.6323441221817e5, -0.2203560539926,
              -0.2087738308480e2, -0.1508165207553e-2, 0.2740740634030e1,
              0.8687687989627e-6, -0.1451419251928e-3, -0.3040346241285e-8,
              0.4712050805815e-8, -0.2639772456566e-5]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon monoxide of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004",
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"},
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 500., "Pmax": 100000.0, "rhomax": 33.84,
        "Pmin": 15.45, "rhomin": 30.33,

        "nr1":  [0.92310041400851, -0.248858452058e1, 0.58095213783396,
                 0.28859164394654e-1, 0.70256257276544e-1, 0.21687043269488e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 0.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.13758331015182, -0.51501116343466e-1, -0.14865357483379,
                -0.38857100886810e-1, -0.29100433948943e-1, 0.14155684466279e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1, MBWR, GERG

    _surface = {"sigma": [0.02843],
                "exp": [1.148]}
    _melting = {"eq": 1, "Tref": 1, "Pref": 1000,
                "Tmin": Tt, "Tmax": 1000.0,
                "a1": [-142.941, 0.0195608], "exp1": [0, 2.10747],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.61192e1, 0.10411e1, -0.62162e1, 0.10437e2, -0.76813e1],
        "exp": [1.0, 1.5, 3.9, 4.6, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.29570e1, -0.42880e1, 0.87643e1, -0.84001e1, 0.36372e1],
        "exp": [0.398, 0.735, 1.08, 1.5, 1.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.25439e1, -0.55601e1, -0.85276e1, -0.51163e1, -0.17701e2, -0.29858e2],
        "exp": [0.395, 1.21, 3.0, 3.5, 6.0, 8.0]}

    visco0 = {"eq": 2, "omega": 2,
              "__name__": "NIST",
              "__doi__": {"autor": "",
                          "title": "Coefficients are taken from NIST14, Version 9.08",
                          "ref": "",
                          "doi": ""},
              "ek": 91.7, "sigma": 0.369,
              "n_chapman": 0.141374566/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-8.89560281339404, -507.15174441, 9.03858480666,
                    5287.58110665, 0.322741446715, -49.2143622937, -23.7275041203],
              "rhoc": 10.85}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "NIST14",
               "__doi__": {"autor": "",
                           "title": "Coefficients are taken from NIST14, Version 9.08",
                           "ref": "",
                           "doi": ""},

               "Tref": 91.7, "kref": 1e-3,
               "no": [1.35558587, -0.16380500617, 1],
               "co": [0, -1, -96],

               "Trefb": 132.85, "rhorefb": 10.85, "krefb": 1e-3,
               "nb": [4.57815545028, 62.5180582967, -2.52594417152,
                      -65.0403809001, 4.06956197472, 18.0214260963],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.4449e-9, "Tcref": 199.29}

    _thermal = thermo0,
