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


class Hexene_1(MEoS):
    """Multiparameter equation of state for 1-hexene"""
    name = "1-hexene"
    CASNumber = "592-41-6"
    formula = "CH2=CH-(CH2)3-CH3"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(238.1713284)
    Tc = unidades.Temperature(504)
    Pc = unidades.Pressure(3062.97, "kPa")
    M = 84.15948  # g/mol
    Tt = unidades.Temperature(133.39)
    Tb = unidades.Temperature(336.61)
    f_acent = 0.27
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 35

    CP0 = {"ao": 4,
           "ao_exp": [8.65, 14.1, 21.9],
           "exp": [360, 3534, 1473]}

    betken = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for 1-hexene of Betken et al."
                    "(2023)",
        "__doi__": {"autor": "Betken, B., Beckmüller, R., Javed, M.A., "
                             "Baumhögger, E., Span, R. Vrabec, J., Thol, M.",
                    "title": "Thermodynamic properties fo 1-hexene - "
                             "Measurements and Modeling",
                    "ref": "J. Chem. Thermo., 176 (2023) 106881",
                    "doi": "10.1016/j.jct.2022.106881"},

        "R": 8.314462618,
        "cp": CP0,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 245000.0, "rhomax": 9.66,

        # Typo in last decimal in first parameters, used the value in FLD
        "nr1": [0.040441989, 1.8522012, -2.1391357, -0.77947556, 0.21159454],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.371, 0.855, 0.995, 0.553],

        "nr2": [-3.3264005, -1.0902532, 0.59957238, -1.2866639, -0.02127171],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.31, 2.3, 0.679, 1.45, 1.08],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [3.9185489, -0.45724185, -0.82698194, -1.0764178, -0.002580659],
        "d3": [1, 3, 2, 2, 1],
        "t3": [0.751, 0.863, 0.67, 0.638, 0.677],
        "alfa3": [0.862, 1.09, 0.802, 1.14, 6.47],
        "beta3": [0.766, 0.903, 0.714, 1.3, 212.14],
        "gamma3": [1.193, 1.297, 1.11, 0.879, 1.09],
        "epsilon3": [0.765, 0.746, 0.728, 0.498, 0.912]}

    eq = (betken, )

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.526, 2.924, -2.8748, -3.584, -3.229],
        "t": [1., 1.5, 2.02, 4.45, 17.8]}

    _liquid_Density = {
        "eq": 1,
        "n": [1.5987, 1.1019, 2.5409, -5.0722, 3.0936],
        "t": [0.29, 0.744, 5.68, 6.76, 7.66]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.3518, -5.7777, -16.194, -48.539, -81.266, -174.12],
        "t": [0.33, 1.04, 2.76, 5.98, 11.91, 21.4]}

    visco0 = {
        "__name__": "Sotiriadou (2023)",
        "__doi__": {
            "autor": "Sotiriadou, S., Ntonti, E., Assael, M.J., Huber, M.L.,",
            "title": "Reference Correlations of the Viscosity and Thermal "
                     "Conductivity of 1-Hexene from the Triple Point to High "
                     "Temmperatures and Pressures",
            "ref": "Int. J. Thermophysics 44 (2023) 108",
            "doi": "10.1007/s10765-023-03217-y"},

        "eq": 1, "omega": 0,
        "ek": 318, "sigma": 0.62,

        "Toref": Tc,
        "no_num": [0.0939564, 3.70702, -5.09947, 12.9844, -1.74165, 0.149722],
        "to_num": [0, 1, 2, 3, 4, 5],
        "no_den": [0.373772, -0.467423, 1.00],
        "to_den": [0, 1, 2],

        "Tref_virial": 318,
        "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717,
                     2491.6597, -787.26086, 14.085455, -0.34664158],
        "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

        "special": "_mur"}

    def _mur(self, rho, T, fase):
        """Special term of residual viscosity for Sotiriadou correlation"""
        Tr = T/self.Tc
        rhor = rho/238.1713

        # Eq 8
        mur = rhor**(2/3)*Tr**0.5 * (
            11.277217*rhor*(1+1/Tr) + (-12.898822*rhor + 1.0917037*rhor**5)
            / (10.742116+8.475744*Tr-6.74312*rhor+rhor**2-1.842042*Tr*rhor))

        return mur

    _viscosity = (visco0, )

    thermo0 = {
        "__name__": "Sotiriadou (2023)",
        "__doi__": {
            "autor": "Sotiriadou, S., Ntonti, E., Assael, M.J., Huber, M.L.,",
            "title": "Reference Correlations of the Viscosity and Thermal "
                     "Conductivity of 1-Hexene from the Triple Point to High "
                     "Temmperatures and Pressures",
            "ref": "Int. J. Thermophysics 44 (2023) 108",
            "doi": "10.1007/s10765-023-03217-y"},

        "eq": 1,

        "Toref": Tc,
        "no_num": [-4.02537e-4, 3.55864e-3, -3.98674e-4, -1.72629e-2,
                   5.23353e-2, -1.38719e-2, 1.30184e-3],
        "to_num": [0, 1, 2, 3, 4, 5, 6],
        "no_den": [0.238788, -0.526011, 1],
        "to_den": [0, 1, 2],

        "Tref_res": Tc, "rhoref_res": 238.1713, "kref_res": 1e-3,
        "nr": [17.4154, -5.92427, -53.4126, 25.7922, 59.5765, -20.8764,
               -19.3526, 3.91783, 2.09622, 0.158646],
        "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
        "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

        "critical": 3,
        "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
        "Xio": 0.235e-9, "gam0": 0.056, "qd": 0.698e-9, "Tcref": 1.5*Tc}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_betken(self):
        """Table 17, Pag 20"""
        st = Hexene_1(T=300, rhom=0.001)
        self.assertEqual(round(st.P.MPa, 9), 0.002490162)
        self.assertEqual(round(st.cpM.kJkmolK, 7), 130.2264437)
        self.assertEqual(round(st.w, 7), 177.6851621)
        self.assertEqual(round(st.hM.kJkmol, 5), 23687.70871)
        self.assertEqual(round(st.sM.kJkmolK, 8), 99.94829959)
        self.assertEqual(round(st.aM.kJkmol, 6), -8786.943313)

        st = Hexene_1(T=300, rhom=8)
        self.assertEqual(round(st.P.MPa, 9), 6.036182516)
        self.assertEqual(round(st.cpM.kJkmolK, 7), 182.8600112)
        self.assertEqual(round(st.w, 6), 1105.511473)
        self.assertEqual(round(st.hM.kJkmol, 5), -6538.47503)
        self.assertEqual(round(st.sM.kJkmolK, 7), -22.9515988)
        self.assertEqual(round(st.aM.kJkmol, 6), -407.518201)

        st = Hexene_1(T=450, rhom=5.8)
        self.assertEqual(round(st.P.MPa, 9), 1.450738906)
        self.assertEqual(round(st.cpM.kJkmolK, 7), 257.9030576)
        self.assertEqual(round(st.w, 7), 403.7924555)
        self.assertEqual(round(st.hM.kJkmol, 5), 25284.94461)
        self.assertEqual(round(st.sM.kJkmolK, 8), 63.84057209)
        self.assertEqual(round(st.aM.kJkmol, 6), -3693.440233)

        st = Hexene_1(T=450, rhom=0.07)
        self.assertEqual(round(st.P.MPa, 9), 0.250858298)
        self.assertEqual(round(st.cpM.kJkmolK, 7), 187.3573887)
        self.assertEqual(round(st.w, 7), 207.5147257)
        self.assertEqual(round(st.hM.kJkmol, 5), 46846.16148)
        self.assertEqual(round(st.sM.kJkmolK, 7), 124.0529255)
        self.assertEqual(round(st.aM.kJkmol, 5), -12561.34495)

        st = Hexene_1(T=600, rhom=3)
        self.assertEqual(round(st.P.MPa, 9), 8.033819707)
        self.assertEqual(round(st.cpM.kJkmolK, 7), 304.1341828)
        self.assertEqual(round(st.w, 7), 197.8662261)
        self.assertEqual(round(st.hM.kJkmol, 5), 66611.43496)
        self.assertEqual(round(st.sM.kJkmolK, 7), 140.0031906)
        self.assertEqual(round(st.aM.kJkmol, 5), -20068.41931)

    def test_Sotiriadou(self):
        """Section 4.2, Pag 22"""
        st = Hexene_1(T=300, rhom=0)
        self.assertEqual(round(st.mu.muPas, 4), 6.7237)
        self.assertEqual(round(st.k.mWmK, 3), 12.589)

        st = Hexene_1(T=300, rho=700)
        self.assertEqual(round(st.mu.muPas, 2), 364.37)
        # Critical enhancement don't reproduce paper, maybe a error in it
        # because critical enhancement at a point so far to critical point
        # self.assertEqual(round(st.k.mWmK, 3), 132.139)

        st = Hexene_1(T=500, x=0.5)
        self.assertEqual(round(st.Liquido.rho, 2), 339.09)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 53.46)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 64.80)
        self.assertEqual(round(st.Gas.rho, 2), 142.38)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 21.70)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 50.65)
