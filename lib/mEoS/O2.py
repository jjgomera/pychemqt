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


from unittest import TestCase

from lib.meos import MEoS
from lib import unidades


class O2(MEoS):
    """Multiparameter equation of state for oxygen"""
    name = "oxygen"
    CASNumber = "7782-44-7"
    formula = "O2"
    synonym = "R-732"
    rhoc = unidades.Density(436.143644)
    Tc = unidades.Temperature(154.581)
    Pc = unidades.Pressure(5043.0, "kPa")
    M = 31.9988  # g/mol
    Tt = unidades.Temperature(54.361)
    Tb = unidades.Temperature(90.1878)
    f_acent = 0.0222
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 47
    _Tr = unidades.Temperature(150.875090)
    _rhor = unidades.Density(439.519141)
    _w = 0.023479051

    # Use the Cp0 expresion from alternate reference with standard term
    # Stewart, R.B., Jacobsen, R.T., Wagner, W.
    # Thermodynamic Properties of Oxygen from the Triple Point to 300 K with
    # Pressures to 80 MPa
    # J. Phys. Chem. Ref. Data 20(5) (1991) 917-1021
    # doi: 10.1063_1.555897
    CP1 = {"ao": 3.521876773671,
           "an": [-0.4981998537119e4, 0.2302477799952e3, -0.3455653235107e1,
                  -0.4354202160244e-4, 0.1346353450132e-7, .1620598259591e-10],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [1.031468515726], "exp": [2239.18105],
           "ao_hyp": [], "hyp": []}

    Fi1 = {"ao_log": [1, 2.50146],
           "pow": [0, 1],
           "ao_pow": [10.001843586, -14.996095135],
           "ao_exp": [], "titao": [],
           "ao_hyp": [1.07558, 1.01334, 0, 0],
           "hyp": [14.461722565, 7.223325463, 0, 0]}

    schmidt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Schmidt and "
                    "Wagner (1985)",
        "__doi__": {"autor": "Schmidt, R., Wagner, W.",
                    "title": "A New Form of the Equation of State for Pure "
                             "Substances and its Application to Oxygen",
                    "ref": "Fluid Phase Equuilibria. 19 (1985) 175-200",
                    "doi": "10.1016/0378-3812(85)87016-3"},

        "R": 8.31434,
        "cp": CP1,
        # "ref": "OTO",
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 8682, "so": 205.037},

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 82000.0, "rhomax": 43.348,
        "Pmin": 0.14628, "rhomin": 40.816,

        "nr1": [0.39837687490, -0.1846157454e1, 0.4183473197, 0.2370620711e-1,
                0.9771730573e-1, 0.3017891294e-1, 0.2273353212e-1,
                0.1357254086e-1, -0.40526989430e-1, 0.54546285150e-3,
                0.51131822770e-3, 0.29534668830e-6, -0.86876450720e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.2127082589, 0.8735941958e-1, 0.127550919, -0.9067701064e-1,
                -0.3540084206e-1, -0.3623278059e-1, 0.132769929e-1,
                -0.3254111865e-3, -0.8313582932e-2, 0.2124570559e-2,
                -0.8325206232e-3, -0.2626173276e-4, 0.2599581482e-2,
                0.9984649663e-2, 0.2199923153e-2, -0.2591350486e-1,
                -0.12596308480, 0.14783556370, -0.10112510780e-1],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*20}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 82000.0, "rhomax": 43.348,
        "Pmin": 73.5, "rhomin": 29.2,

        "nr1": [0.88878286369701, -0.24879433312148e1, 0.59750190775886,
                0.96501817061881e-2, 0.71970428712770e-1, 0.22337443000195e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18558686391474, -0.038129368035760, -0.15352245383006,
                -0.026726814910919, -0.025675298677127, 0.95714302123668e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.623, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*20,

        "nr3": [],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for oxygen of Younglove (1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, "
                             "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                             "Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                    "doi": ""},

        "R": 8.31411,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 40.820,
        "Pmin": 0.148, "rhomin": 40.820,

        "b": [None, -0.4365859650e-3, 0.2005820677, -0.4197909916e1,
              0.1878215317e3, -0.1287473398e5, 0.1556745888e-4, 0.001343639359,
              -2.228415518, 0.4767792275e4, 0.4790846641e-6, 0.002462611107,
              -0.192189168, -0.6978320847e-5, -0.6214145909e-3, -0.1860852567,
              0.2609791417e-4, -0.2447611408e-6, 0.1457743352e-3,
              -0.1726492873e-5, -0.238489252e4, -0.2301807796e6, -27.90303526,
              0.9400577575e5, -0.04169449637, 2.008497853, -0.125607652e-3,
              -0.6406362964, -0.2475580168e-7, 0.1346309703e-4,
              -0.116150247e-9, -0.1034699798e-7, 0.2365936964e-6]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for oxygen of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",
        "M": 31.999, "Tc": 154.595, "rhoc": 436.14/31.999,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 43.348,
        "Pmin": 0.14603, "rhomin": 40.885,

        "nr1": [0.88878286, -0.24879433e1, 0.59750191, 0.96501817e-2,
                0.71970429e-1, 0.22337443e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18558686, -0.38129368e-1, -0.15352245, -0.26726815e-1,
                -0.25675299e-1, 0.95714302e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = schmidt, MBWR, GERG, shortSpan
    _PR = -0.003157

    _surface = {"sigma": [0.03843], "exp": [1.225]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [3.9578, 0.0065], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.575, 1.028, -8.96, -5.15],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.5, 2.5]}
    _melting = {"eq": 2, "Tref": Tt, "Pref": 0.14633,
                "Tmin": Tt, "Tmax": 300.0,
                "a1": [], "exp1": [],
                "a2": [-0.32463539e2, 0.14278011e3, -0.14702341e3, 0.520012e2],
                "exp2": [0.0625, 0.125, 0.1875, 0.25],
                "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 0.14633,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-20.714], "exp2": [1.06],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.60595e1, 0.13050e1, -0.54178, -0.19410e1, 0.35514],
        "exp": [1., 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16622e1, 0.76846, -0.10041, 0.20480, 0.11551e-1],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.22695e1, -0.46578e1, -0.99480e1, -0.22845e2, -0.45190e2,
               -0.25101e2],
        "exp": [0.3785, 1.07, 2.7, 5.5, 10., 20.]}

    visco0 = {"__name__": "Lemmon (2004)",
              "__doi__": {
                  "autor": "Lemmon, E.W., Jacobsen, R.T.",
                  "title": "Viscosity and Thermal Conductivity Equations for "
                           "Nitrogen, Oxygen, Argon, and Air",
                  "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                  "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

              "eq": 1, "omega": 1,
              "ek": 118.5, "sigma": 0.3428,

              "Tref_res": 154.581, "rhoref_res": 13.63*M,
              "nr": [17.67, 0.4042, 0.0001077, 0.3510, -13.67],
              "tr": [0.05, 0, 2.1, 0, 0.5],
              "dr": [1, 5, 12, 8, 1],
              "gr": [0, 0, 0, 1, 1],
              "cr": [0, 0, 0, 1, 2]}

    visco1 = {"__name__": "Younglove (1982)",
              "__doi__": {
                  "autor": "Younglove, B.A.",
                  "title": "Thermophysical Properties of Fluids. I. Argon, "
                           "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                           "Trifluoride, and Oxygen",
                  "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                  "doi": ""},

              "eq": 2, "omega": 2,

              "ek": 113., "sigma": 0.3437,
              "n_chapman": 0.15099557923496,
              "t_chapman": 0.0,
              "collision": [-67.2093902106092, 277.148660965491, -399.192753863192,
                            166.828729537446, 143.163477478684, -191.767060368781,
                            98.4332230147836, -22.9410694301649, 2.12402264924749],

              "F": [1.39279625307e-2, -6.51536010579e-3, 1.4, 100],
              "E": [-14.45497211, 243.40689667, 12.9006761056004,
                    -1949.07966423848, -5.62078436742e-2,
                    21.3075467849104, 48.9965711691056],
              "rhoc": 13.5942597847419}

    visco2 = {"__name__": "Laesecke (1990)",
              "__doi__": {
                  "autor": "Laesecke, A., Krauss, R., Stephan, K., Wagner, W.",
                  "title": "Transport Properties of Fluid Oxygen",
                  "ref": "J. Phys. Chem. Ref. Data 19(5) (1990) 1089-1122",
                  "doi": "10.1063/1.555863"},

              "eq": 1, "omega": 1,

              "ek": 116.2, "sigma": 0.34318867,
              "collision": [0.46649, -0.57015, 0.19164, -0.03708, 0.00241],

              "Tref_res": 1, "rhoref_res": 13.63*M, "muref_res": 18.8928,
              "nr": [-1.7993647367, -0.397230772, 0.312536267, -0.0615559341],
              "tr": [0, 0, 0, 0],
              "dr": [0, 1, 2, 3],
              "nr_num": [-5.60288207],
              "tr_num": [0],
              "dr_num": [0],
              "nr_den": [1.0, -3.1138112],
              "tr_den": [0, 0],
              "dr_den": [1, 0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Lemmon (2004)",
               "__doi__": {
                  "autor": "Lemmon, E.W., Jacobsen, R.T.",
                  "title": "Viscosity and Thermal Conductivity Equations for "
                           "Nitrogen, Oxygen, Argon, and Air",
                  "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                  "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

               "eq": 1,

               "Toref": 154.581, "koref": 1e-3,
               "no_visco": 1.036,
               "no": [6.283, -4.262],
               "to": [0.9, 0.6],

               "Tref_res": 154.581, "rhoref_res": 13.63*M, "kref_res": 1e-3,
               "nr": [15.31, 8.898, -0.7336, 6.728, -4.374, -0.4747],
               "tr": [0, 0, 0.3, 4.3, 0.5, 1.8],
               "dr": [1, 3, 4, 5, 7, 10],
               "cr": [0, 0, 0, 2, 2, 2],
               "gr": [0, 0, 0, 1, 1, 1],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.24e-9, "gam0": 0.055, "qd": 0.51e-9, "Tcref": 309.162}

    thermo1 = {"eq": 3,
               "__name__": "Younglove (1982)",
               "__doi__": {"autor": "Younglove, B.A.",
                           "title": "Thermophysical Properties of Fluids. I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen",
                           "ref": "J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.",
                           "doi": ""},

               "ek": 113, "sigma": 0.3437,
               "Nchapman": 0.15099557923496,
               "tchapman": 0,
               "b": [-1.41202117453516, 8.06267523869911, -19.44147946395,
                     25.78193316324, -20.5167203343277, 10.0087040966906,
                     -2.90450673487991, 0.459605807669332, -3.01906029521e-2],
               "F": [0.00097916328, 0.00089116658, 1.12, 100],
               "E": [-21.520741137, 473.50508788, 11.9072051301147,
                     -2122.44247203833, 0, 0, 0],
               "rhoc": 31.251171918947,
               "ff": 2.21064,
               "rm": 0.000000038896}

    thermo2 = {"eq": 1,
               "__name__": "Laesecke (1990)",
               "__doi__": {
                  "autor": "Laesecke, A., Krauss, R., Stephan, K., Wagner, W.",
                  "title": "Transport Properties of Fluid Oxygen",
                  "ref": "J. Phys. Chem. Ref. Data 19(5) (1990) 1089-1122",
                  "doi": "10.1063/1.555863"},

               "__test__": """
                    >>> st=O2(T=70, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    181.3 5.981
                    >>> st=O2(T=80, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    167.2 7.105
                    >>> st=O2(T=90, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    153.1 8.258
                    >>> st=O2(T=100, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    139.0 9.477
                    >>> st=O2(T=110, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    124.8 10.82
                    >>> st=O2(T=120, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    110.3 12.41
                    >>> st=O2(T=130, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    95.44 14.44
                    >>> st=O2(T=140, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    79.48 17.41
                    >>> st=O2(T=150, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    60.11 23.16
                    >>> st=O2(T=154.581, x=0.5, thermo=2)
                    >>> print "%0.4g %0.4g" % (st.Liquido.k.mWmK, st.Gas.k.mWmK)
                    37.47 37.47
                    """ # Table 3, Pag 1100
                    """
                    >>> st=O2(T=70, P=1e5, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    181.0
                    >>> st=O2(T=180, P=1e5, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    16.80
                    >>> st=O2(T=110, P=1e6, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    125.1
                    >>> st=O2(T=70, P=1e7, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    184.2
                    >>> st=O2(T=120, P=2e7, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    129.2
                    >>> st=O2(T=50, P=1e8, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    201.6
                    >>> st=O2(T=400, P=1e5, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    33.24
                    >>> st=O2(T=200, P=1e8, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    120.9
                    >>> st=O2(T=500, P=1e6, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    40.29
                    >>> st=O2(T=800, P=1e5, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    59.15
                    >>> st=O2(T=500, P=1e8, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    70.64
                    >>> st=O2(T=1400, P=1e6, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    90.84
                    >>> st=O2(T=1000, P=1e7, thermo=2)
                    >>> print "%0.4g" % st.k.mWmK
                    71.59
                    """, # Table 5, Pag 1109

               "Tref": 1, "kref": 1e-3,
               "no": [0.5825413, 0.0321266],
               "co": [-97, -98],

               "Trefb": 1, "rhorefb": 13.63, "krefb": 4.909e-3,
               "nb": [2.32825085, 4.23024231, -3.60798307, 2.01675631, -0.289731736],
               "tb": [0, 0, 0, 0, 0],
               "db": [1, 2, 3, 4, 5],
               "cb": [0, 0, 0, 0, 0],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 1.6e-10, "gam0": 0.08391, "qd": 0.4167e-9, "Tcref": 309.162}

    _thermal = thermo0, thermo1, thermo2


class Test(TestCase):

    def test_schmidt(self):
        # Table 4 pag 191, Saturation densities
        st = O2(T=150, x=0.5)
        self.assertEqual(round(st.Liquido.rhoM, 4), 21.1096)
        self.assertEqual(round(st.Gas.rhoM, 5), 6.71701)

        st = O2(T=154.515, x=0.5)
        self.assertEqual(round(st.Liquido.rhoM, 4), 15.0637)
        self.assertEqual(round(st.Gas.rhoM, 4), 11.9252)

        st = O2(T=154.565, x=0.5)
        self.assertEqual(round(st.Liquido.rhoM, 4), 14.4607)
        self.assertEqual(round(st.Gas.rhoM, 4), 12.3876)

        # Selected values from tables in Stewart, Jacobsen and Wagner paper
        # Table 9 pag 948, ideal gas properties
        # FIXME: Fix 1/τ terms in cpo to phio conversion

        # st = O2()
        # self.assertEqual(round(st._Cp0(60), 2), 158.32)
        # st = O2(T=60, P=1e5)
        # self.assertEqual(round(st.sM0.JmolK, 2), 158.32)
        # self.assertEqual(round(st.hM0.Jmol, 1), 1738.2)
        # self.assertEqual(round(st.cv0M.JmolK, 3), 20.809)
        # self.assertEqual(round(st.cp0M.JmolK, 3), 29.123)

        # st = O2(T=200, P=1e5)
        # self.assertEqual(round(st.sM0.JmolK, 2), 193.37)
        # self.assertEqual(round(st.hM0.Jmol, 1), 5814.4)
        # self.assertEqual(round(st.cv0M.JmolK, 3), 20.812)
        # self.assertEqual(round(st.cp0M.JmolK, 3), 29.127)

        # st = O2(T=500, P=1e5)
        # self.assertEqual(round(st.sM0.JmolK, 2), 220.58)
        # self.assertEqual(round(st.hM0.Jmol, 1), 14766.5)
        # self.assertEqual(round(st.cv0M.JmolK, 3), 22.778)
        # self.assertEqual(round(st.cp0M.JmolK, 3), 31.093)

        # st = O2(T=1000, P=1e5)
        # self.assertEqual(round(st.sM0.JmolK, 2), 243.47)
        # self.assertEqual(round(st.hM0.Jmol, 1), 31387.4)
        # self.assertEqual(round(st.cv0M.JmolK, 3), 26.563)
        # self.assertEqual(round(st.cp0M.JmolK, 3), 34.877)

        # st = O2(T=2000, P=1e5)
        # self.assertEqual(round(st.sM0.JmolK, 2), 268.66)
        # self.assertEqual(round(st.hM0.Jmol, 1), 67882.8)
        # self.assertEqual(round(st.cv0M.JmolK, 3), 29.469)
        # self.assertEqual(round(st.cp0M.JmolK, 3), 37.783)

        # Table 10 pag 949, saturated states
        # st = O2(T=54.361, x=0.5)
        # self.assertEqual(round(st.P.MPa, 6), 0.000146)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 40.816)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), -6193.7)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 66.94)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 38.25)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 53.54)
        # self.assertEqual(round(st.Liquido.w, 0), 1123)
        # self.assertEqual(round(st.Gas.rhoM, 7), 0.0003237)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 1573.1)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 209.81)
        # self.assertEqual(round(st.Gas.w, 0), 140)

        # st = O2(T=80, x=0.5)
        # self.assertEqual(round(st.P.MPa, 5), 0.03012)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 37.203)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), -4817.7)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 87.66)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 31.03)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 53.81)
        # self.assertEqual(round(st.Liquido.w, 0), 987)
        # self.assertEqual(round(st.Gas.rhoM, 6), 0.045891)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 2295.9)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 176.58)
        # self.assertEqual(round(st.Gas.w, 0), 168)

        # st = O2(T=100, x=0.5)
        # self.assertEqual(round(st.P.MPa, 4), 0.2540)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 34.092)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), -3724.6)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 99.78)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 28.63)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 55.60)
        # self.assertEqual(round(st.Liquido.w, 0), 822)
        # self.assertEqual(round(st.Gas.rhoM, 5), 0.32579)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 2758.9)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 164.61)
        # self.assertEqual(round(st.Gas.cvM.JmolK, 2), 21.61)
        # self.assertEqual(round(st.Gas.cpM.JmolK, 2), 32.20)
        # self.assertEqual(round(st.Gas.w, 0), 184)

        # st = O2(T=120, x=0.5)
        # self.assertEqual(round(st.P.MPa, 4), 1.0223)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 30.434)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), -2554.8)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 110.21)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 26.98)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 61.67)
        # self.assertEqual(round(st.Liquido.w, 0), 641)
        # self.assertEqual(round(st.Gas.rhoM, 4), 1.2284)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 3002.0)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 156.52)
        # self.assertEqual(round(st.Gas.cvM.JmolK, 2), 23.73)
        # self.assertEqual(round(st.Gas.cpM.JmolK, 2), 40.84)
        # self.assertEqual(round(st.Gas.w, 0), 189)

        # st = O2(T=140, x=0.5)
        # self.assertEqual(round(st.P.MPa, 4), 2.7878)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 25.415)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), -1172.2)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 120.35)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 26.63)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 86.10)
        # self.assertEqual(round(st.Liquido.w, 0), 423)
        # self.assertEqual(round(st.Gas.rhoM, 4), 3.6487)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 2833.1)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 148.96)
        # self.assertEqual(round(st.Gas.cvM.JmolK, 2), 28.27)
        # self.assertEqual(round(st.Gas.cpM.JmolK, 2), 75.82)
        # self.assertEqual(round(st.Gas.w, 0), 182)

        # st = O2(T=154, x=0.5)
        # self.assertEqual(round(st.P.MPa, 4), 4.9307)
        # self.assertEqual(round(st.Liquido.rhoM, 3), 17.096)
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), 500.0)
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), 130.97)
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 36.13)
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 1192.05)
        # self.assertEqual(round(st.Liquido.w, 0), 178)
        # self.assertEqual(round(st.Gas.rhoM, 3), 10.213)
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), 1658.9)
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), 138.49)
        # self.assertEqual(round(st.Gas.cvM.JmolK, 2), 40.31)
        # self.assertEqual(round(st.Gas.cpM.JmolK, 2), 1633.93)
        # self.assertEqual(round(st.Gas.w, 0), 161)

        # st = O2(T=, x=0.5)
        # self.assertEqual(round(st.P.MPa, 2), )
        # self.assertEqual(round(st.Liquido.rhoM, 3), )
        # self.assertEqual(round(st.Liquido.hM.Jmol, 1), )
        # self.assertEqual(round(st.Liquido.sM.JmolK, 2), )
        # self.assertEqual(round(st.Liquido.cvM.JmolK, 2), )
        # self.assertEqual(round(st.Liquido.cpM.JmolK, 2), )
        # self.assertEqual(round(st.Liquido.w, 0), )
        # self.assertEqual(round(st.Gas.rhoM, 7), )
        # self.assertEqual(round(st.Gas.hM.Jmol, 1), )
        # self.assertEqual(round(st.Gas.sM.JmolK, 2), )
        # self.assertEqual(round(st.Gas.w, 0), )

    def test_shortSpan(self):
        # Table III, Pag 46
        st = O2(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.0308)
        self.assertEqual(round(st.P.MPa, 3), 41.051)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.09696)

        st2 = O2(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 50.89)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.26457)

    def test_LemmonTransport(self):
        # Table V, pag 28
        # Viscosity
        self.assertEqual(round(O2(T=100, rhom=0).mu.muPas, 5), 7.70243)
        self.assertEqual(round(O2(T=300, rhom=0).mu.muPas, 4), 20.6307)
        self.assertEqual(round(O2(T=100, rhom=35).mu.muPas, 3), 172.136)
        self.assertEqual(round(O2(T=200, rhom=10).mu.muPas, 4), 22.4445)
        self.assertEqual(round(O2(T=300, rhom=5).mu.muPas, 4), 23.7577)
        self.assertEqual(round(O2(T=154.6, rhom=13.6).mu.muPas, 4), 24.7898)

        # Thermal Conductivity
        self.assertEqual(round(O2(rhom=0, T=100).k.mWmK, 5), 8.94334)
        self.assertEqual(round(O2(rhom=0, T=300).k.mWmK, 4), 26.4403)
        self.assertEqual(round(O2(rhom=35, T=100).k.mWmK, 3), 146.044)
        self.assertEqual(round(O2(rhom=10, T=200).k.mWmK, 4), 34.6125)
        self.assertEqual(round(O2(rhom=5, T=300).k.mWmK, 4), 32.5491)
        self.assertEqual(round(O2(rhom=13.6, T=154.6).k.mWmK, 1), 377.5)

    def test_Laesecke(self):
        kw = {"visco": 2, "thermal": 2}
        # Table 3, Pag 1100, the saturation state is calculated with auxiliary
        # equation so the density used are not exact, so it can possible get
        # differences in returned values
        st = O2(T=100, x=0.5, **kw)
        self.assertEqual(round(st.Liquido.mu.muPas, 0), 146)
        self.assertEqual(round(st.Gas.mu.muPas, 3), 7.436)
        # self.assertEqual(round(st.Liquido.k.mWmK, 1), 139.0)
        # self.assertEqual(round(st.Gas.k.mWmK, 3), 9.477)

        st = O2(T=154, x=0.5, **kw)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 32.06)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 19.53)
        # self.assertEqual(round(st.Liquido.k.mWmK, 2), 46.61)
        # self.assertEqual(round(st.Gas.k.mWmK, 2), 30.02)

        # Table 4, Pag 1101, Viscosity single phase region
        self.assertEqual(round(O2(T=70, P=1e5, **kw).mu.muPas, 1), 352.1)
        self.assertEqual(round(O2(T=180, P=1e5, **kw).mu.muPas, 1), 13.4)
        self.assertEqual(round(O2(T=110, P=1e6, **kw).mu.muPas, 1), 118.0)
        self.assertEqual(round(O2(T=70, P=1e7, **kw).mu.muPas, 1), 401.1)
        self.assertEqual(round(O2(T=120, P=2e7, **kw).mu.muPas, 1), 124.3)
        self.assertEqual(round(O2(T=150, P=1e8, **kw).mu.muPas, 1), 166.9)
        self.assertEqual(round(O2(T=400, P=1e5, **kw).mu.muPas, 2), 25.91)
        self.assertEqual(round(O2(T=200, P=1e8, **kw).mu.muPas, 2), 105.35)
        self.assertEqual(round(O2(T=500, P=1e6, **kw).mu.muPas, 2), 30.56)
        self.assertEqual(round(O2(T=800, P=1e5, **kw).mu.muPas, 2), 42.27)
        self.assertEqual(round(O2(T=500, P=1e8, **kw).mu.muPas, 2), 49.74)
        self.assertEqual(round(O2(T=1400, P=1e6, **kw).mu.muPas, 2), 61.29)
        self.assertEqual(round(O2(T=1000, P=1e7, **kw).mu.muPas, 2), 49.44)
