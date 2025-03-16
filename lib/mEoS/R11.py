#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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
from lib.mEoS import R134a


class R11(MEoS):
    """Multiparameter equation of state for R11"""
    name = "trichlorofluoromethane"
    CASNumber = "75-69-4"
    formula = "CCl3F"
    synonym = "R11"
    _refPropName = "R11"
    _coolPropName = "R11"
    rhoc = unidades.Density(554.)
    Tc = unidades.Temperature(471.11)
    Pc = unidades.Pressure(4407.638, "kPa")
    M = 137.368  # g/mol
    Tt = unidades.Temperature(162.68)
    Tb = unidades.Temperature(296.858)
    f_acent = 0.18875
    momentoDipolar = unidades.DipoleMoment(0.450, "Debye")
    id = 217

    CP1 = {"ao": 4.00564923,
           "an": [2.228875e-4], "pow": [1],
           "ao_exp": [1, 2, 1, 2, 1, 2],
           "exp": [1561.076, 1218.647, 770.035, 572.634, 502.854, 346.746]}

    CP2 = {"ao": 4.0000024,
           "ao_exp": [0.32960961e1, 0.28401126e1, 0.40350474, 0.30739271e1],
           "exp": [381.63168, 1368.22648, 3435.66931, 689.55053]}

    jacobsen = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-11 of Jacobsen (1992)",
        "__doi__": {"autor": "Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.",
                    "title": "A Fundamental Equation for "
                             "Trichlorofluoromethane (R-11)",
                    "ref": "Fluid Phase Equilibria, 80 (1992) 45-56",
                    "doi": "10.1016/0378-3812(92)87054-Q"},

        "R": 8.31451,
        "cp": CP1,
        "ref": {"Tref": Tt, "Pref": 1.0, "ho": 53727.59, "so": 264.0369},

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 30000.0, "rhomax": 12.88,

        "nr1": [0.125993633881e1, -0.260818574641e1, 0.982122542463e-2,
                -0.106085385839e1, 0.122820363510e1, 0.118000776439,
                -0.698956926463e-3, -0.355428373358e-1, 0.197169579643e-2,
                -0.848363012252e-2, 0.417997567653e-2, -0.242772533848e-3,
                0.313371368974e-2, 0.396182646586e-5],
        "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 8],
        "t1": [0.5, 1.5, 5, 1, 1.5, 0, 5, 2, 3, 1, 2, 4, 1, 4],

        "nr2": [0.339736319502, -0.203010634531, -0.1060178599, 0.45156488259,
                -0.339265767612, 0.114338523359, 0.319537833995e-1,
                0.367908259780e-1, -0.961768948364e-5, 0.246717966418e-2,
                -0.167030256045e-2, 0.240710110806e-2, 0.156214678738e-2,
                -0.323352596704e-2],
        "d2": [1, 1, 2, 2, 2, 3, 4, 6, 10, 3, 5, 8, 9, 9],
        "t2": [5, 6, 3.5, 5.5, 7.5, 3, 2.5, 5, 1.5, 11, 9, 13, 5, 9],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 6, 6, 6, 6],
        "gamma2": [1]*14}

    marx = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-11 of Marx (1992)",
        "__doi__": {"autor": "Marx, V., Pruss, A., Wagner, W.",
                    "title": "Neue Zustandsgleichungen fuer R 12, R 22, R 11 "
                             "und R 113. Beschreibung des thermodynamishchen "
                             "Zustandsverhaltens bei Temperaturen bis 525 K "
                             "und Druecken bis 200 MPa",
                    "ref": "Düsseldorf: VDI Verlag, Series 19 "
                           "(Waermetechnik/Kaeltetechnik), No. 57, 1992.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 30000.0, "rhomax": 13.0,

        "nr1": [-0.219644325e1, 0.8562148696, 0.185864982e-1, 0.2807246052,
                -0.8526398864e-1, 0.1090334698e-1],
        "d1": [1, 1, 1, 2, 3, 5],
        "t1": [1.5, 2, 3, 0, 1.5, 1],

        "nr2": [0.4138515982, -0.3125498519, 0.1545749737, 0.1752299625,
                0.2295443969e-1, -0.2094422944e-2, -0.1267942875e-8,
                0.797272861e-2, -0.1520330549, 0.6448637628e-1,
                0.2046144277e-3, -0.4100829613e-4, -0.123188575e-1,
                0.6681486552e-2, -0.6742271171e-7],
        "d2": [1, 1, 2, 3, 5, 7, 14, 1, 2, 3, 11, 11, 4, 4, 10],
        "t2": [-0.5, 3.5, -0.5, 1, -0.5, 2, 4, 8, 8, 8, 4, 6, 18, 21, 33],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4],
        "gamma2": [1]*15}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-11 of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": CP2,
        "ref": "IIR",
        "M": 137.368, "Tc": 471.06, "rhoc": 565/137.368,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 13.0,

        "nr1": [0.10656383e1, -0.32495206e1, 0.87823894, 0.87611569e-1,
                0.29950049e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.42896949, 0.70828452, -0.17391823e-1, -0.37626521,
                0.11605284e-1, -0.89550567e-1, -0.30063991e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = jacobsen, marx, shortSpan
    _PR = [-0.2807, -15.0917]

    _surface = {"sigma": [0.06212], "exp": [1.247]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70742e1, 0.38118e1, -0.32850e1, -0.76340e1, 0.50598e1],
        "t": [1.0, 1.5, 1.73, 5.2, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.20368e1, 0.12850e2, -0.22521e2, 0.11340e2, -0.94375],
        "t": [0.357, 1.5, 1.7, 2.0, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0296, -6.0723, -15.890, -63.024, 87.167, -157.15],
        "t": [0.417, 1.25, 3.1, 6.8, 10.0, 12.0]}

    trnECS = {"__name__": "Klein (1997)",

              "__doi__": {
                  "autor": "Klein, S.A., McLinden, M.O., Laesecke, A.",
                  "title": "An improved extended corresponding states method "
                           "for estimation of viscosity of pure refrigerants "
                           "and mixtures",
                  "ref": "Int. J. Refrig. 20(3) (1997) 208-217",
                  "doi": "10.1016/s0140-7007(96)00073-4"},

              # Themal conductivity correlation from
              # McLinden, M.O., Klein, S.A., Perkins, R.A.
              # An Extended corresponding states model for the thermal
              # conductivity of refrigerants and refrigerant mixtures
              # Int. J. Refrigeration 23 (2000) 43-63
              # 10.1016/s0140-7007(99)00024-9

              # Parameters updated to upgrades thermal conductivity of
              # reference fluid R134a and reported in refprop
              # Lemmon, E.W., Huber, M.L., McLinden, M.O.
              # NIST Standard Reference Database 23:  Reference Fluid
              # Thermodynamic and Transport Properties-REFPROP, Version 9.1,
              # National Institute of Standards and Technology, Standard
              # Reference Data Program, Gaithersburg, 2013.

              "eq": "ecs",
              "ref": R134a,

              "ek": 363.609, "sigma": 0.5447, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.4e-3], "psi_d": [0],
              "fint": [1.0653851, -0.0250121], "fint_t": [0, 1],
              "chi": [1.0724, -2.2672e-2], "chi_d": [0, 1]}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Testing"""
    def test_shortSpan(self):
        """Table III, Pag 117"""
        st = R11(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.6879)
        self.assertEqual(round(st.P.MPa, 3), 6.077)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3618)

        st2 = R11(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 129.72)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.26917)
