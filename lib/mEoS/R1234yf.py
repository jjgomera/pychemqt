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


class R1234yf(MEoS):
    """Multiparameter equation of state for R1234yf"""
    name = "2,3,3,3-tetrafluoropropene"
    CASNumber = "754-12-1"
    formula = "CF3CF=CH2"
    synonym = "R-1234yf"
    rhoc = unidades.Density(475.553441976)
    Tc = unidades.Temperature(367.85)
    Pc = unidades.Pressure(3382.2, "kPa")
    M = 114.0415928  # g/mol
    Tt = unidades.Temperature(220.)
    Tb = unidades.Temperature(243.7)
    f_acent = 0.276
    momentoDipolar = unidades.DipoleMoment(2.48, "Debye")
#    id = 671

    Fi1 = {"ao_log": [1, 4.944],
           "pow": [0, 1],
           "ao_pow": [-12.837928, 8.042605],
           "ao_exp": [7.549, 1.537, 2.03, 7.455],
           "titao": [718/Tc, 877/Tc, 4465/Tc, 1755/Tc],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-11.412027, -52.9180363],
           "ao_exp": [5.2829, 6.96022, 7.04266],
           "titao": [354/Tc, 965/Tc, 1981/Tc],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234yf of Richter et al. (2011).",
        "__doi__": {"autor": "Richter, M., McLinden, M.O., and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf): Vapor Pressure and p-rho-T Measurements and an Equation of State",
                    "ref": "J. Chem. Eng. Data, 2011, 56 (7), pp 3254–3264",
                    "doi": "10.1021/je200369m"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 410.0, "Pmax": 30000.0, "rhomax": 11.64,
        "Pmin": 31.5, "rhomin": 11.63,

        "nr1": [0.4592563e-1, 0.1546958e1, -0.2355237e1, -0.4827835, 0.1758022],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.32, 0.929, 0.94, 0.38],

        "nr2": [-0.1210006e1, -0.6177084, 0.6805262, -0.6968555, -0.2695779e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.28, 1.76, 0.97, 2.44, 1.05],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*7,

        "nr3": [0.1389966e1, -0.4777136, -0.1975184, -0.1147646e1, 0.3428541e-3],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.4, 3.0, 3.5, 1.0, 3.5],
        "alfa3": [1.02, 1.336, 1.055, 5.84, 16.2],
        "beta3": [1.42, 2.31, 0.89, 80., 108.],
        "gamma3": [1.13, 0.67, 0.46, 1.28, 1.2],
        "epsilon3": [0.712, 0.910, 0.677, 0.718, 1.64],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1234yf of Akasaka (2011).",
        "__doi__": {"autor": "Akasaka, R.",
                    "title": "New Fundamental Equations of State with a Common Functional Form for 2,3,3,3-Tetrafluoropropene (R-1234yf) and trans-1,3,3,3-Tetrafluoropropene (R-1234ze(E))",
                    "ref": "Int J Thermophys (2011) 32:1125–1147",
                    "doi": "10.1007/s10765-011-0992-0"},
        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": 240., "Tmax": 400.0, "Pmax": 40000.0, "rhomax": 11.64,
        "Pmin": 31.5, "rhomin": 11.63,

        "nr1": [0.83266757e1, -0.92588001e1, -0.24906043, 0.14422208,
                0.11679917e-1],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.66886, 0.83392, 1.6982, 1.8030, 0.36657],

        "nr2": [-0.16465103, 0.10580795, 0.17135586e-1, -0.16764798e-1,
                -0.12781115e-1, 0.36440802, -0.28535370, -0.96835199e-1,
                0.88063705e-1, 0.18736343e-1, -0.16872191e-1, 0.70032274e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.8666, 1.0194, 0, 1.1655, 8.3101, 6.1459, 8.3495, 6.0422,
               7.444, 15.433, 21.543, 15.499],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.06274], "exp": [1.394]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74697e1, 0.27915e1, -0.21312e1, -0.29531e1],
        "exp": [1.0, 1.5, 1.8, 3.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.19083e1, -0.21383e1, 0.93653e1, -0.98659e1, 0.35859e1],
        "exp": [0.32, 0.56, 0.8, 1.0, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.23511e1, -0.11515e2, -0.53984e1, -0.37937e2],
        "exp": [0.355, 2.45, 1.0, 5.1]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2011)",
               "__doi__": {"autor": "Perkins, R.A. and Huber, M.L.",
                            "title": "Measurement and Correlation of the Thermal Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf) and trans-1,3,3,3-Tetrafluoropropene (R1234ze(E))",
                            "ref": "J. Chem. Eng. Data, 2011, 56 (12), pp 4868–4874",
                            "doi": "10.1021/je200811n"},
               "__test__": """
                    >>> st=R1234yf(T=250, P=5e4)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    2.80006 0.0098481
                    >>> st=R1234yf(T=300, P=1e5)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    4.671556 0.013996
                    >>> st=R1234yf(T=250, P=2e7)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    1299.5 0.088574
                    >>> st=R1234yf(T=300, P=2e7)
                    >>> print "%0.6g %0.5g" % (st.rho, st.k)
                    1182.05 0.075245
                    """, # Table 2, Pag 4872

               "Tref": 367.85, "kref": 1,
               "no": [-0.0102778, 0.0291098, 0.000860643],
               "co": [0, 1, 2],


               "Trefb": 367.85, "rhorefb": 4.17, "krefb": 1.,
               "nb": [-0.368219e-1, 0.883226e-1, -0.705909e-1, 0.259026e-1,
                      -0.32295e-2, 0.397166e-1, -0.772390e-1, 0.664707e-1,
                      -0.249071e-1, 0.336228e-2],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.835e-10, "Tcref": 551.775}

    _thermal = thermo0,
