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


class iC6(MEoS):
    """Multiparameter equation of state for isohexane"""
    name = "isohexane"
    CASNumber = "107-83-5"
    formula = "(CH3)2-CH-(CH2)2-CH3"
    synonym = ""
    rhoc = unidades.Density(233.966)
    Tc = unidades.Temperature(497.7)
    Pc = unidades.Pressure(3040.0, "kPa")
    M = 86.17536  # g/mol
    Tt = unidades.Temperature(119.6)
    Tb = unidades.Temperature(333.36)
    f_acent = 0.2797
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 52

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [6.9259123919, -0.3128629679],
           "ao_exp": [7.9127, 16.871, 19.257, 14.075],
           "titao": [325/Tc, 1150/Tc, 2397/Tc, 5893/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isohexane of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},
        "__test__": """
            >>> st=iC6(T=499, rho=2*86.17536)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            499 2 3058.917 48733.740 113.316 233.627 1129.816 90.210
            """, # Table 10, Pag 842

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 1000000.0, "rhomax": 9.38,
        "Pmin": 7.34e-9, "rhomin": 9.37,

        "nr1": [1.1027, -2.9699, 1.0295, -0.21238, 0.11897, 0.00027738],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40103, -0.034238, -0.43584, -0.11693, -0.019262, 0.0080783],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.05024], "exp": [1.194]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74130e1, 0.16267e1, -0.22311e1, -0.26040e1, -0.29490e1],
        "exp": [1.0, 1.5, 2.62, 4.56, 16.3]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18489e2, -0.43541e2, 0.43985e2, -0.16581e2, 0.64563],
        "exp": [0.59, 0.77, 0.96, 1.15, 3.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.41180e1, -0.61956e1, -0.21190e2, -0.58972e2, -0.15824e3],
        "exp": [0.4824, 1.418, 3.32, 7.1, 16.1]}

    visco0 = {"eq": 2, "omega": 3,
              "__name__": "NIST",
              "__doi__": {"autor": "",
                          "title": "Coefficients are taken from NIST14, Version 9.08",
                          "ref": "",
                          "doi": ""},

              "ek": 368.52, "sigma": 0.61222,
              "n_chapman": 0.2267237/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-13.294469653994, -466.41004563, 15.438316998,
                    -3363.2028894, -0.11398677788, 171.32077134, 2849.7100897],
              "rhoc": 2.727}

    _viscosity = visco0,

    thermo0 = {"eq": 1, "critical": 0,
               "__name__": "NIST14",
               "__doi__": {"autor": "",
                           "title": "Coefficients are taken from NIST14, Version 9.08",
                           "ref": "",
                           "doi": ""},

               "Tref": 368.52, "kref": 1e-3,
               "no": [1.35558587, -0.152808259573429, 1],
               "co": [0, -1, -96],

               "Trefb": 498.05, "rhorefb": 2.727, "krefb": 1e-3,
               "nb": [13.747515904, 10.1607102792, -7.75232868497,
                      0.627943006907, 1.9518640415, -0.293574041046],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6}

    _thermal = thermo0,
