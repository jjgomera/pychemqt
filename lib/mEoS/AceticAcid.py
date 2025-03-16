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


class AceticAcid(MEoS):
    """Multiparameter equation of state for Acetic acid"""
    name = "acetic acid"
    CASNumber = "64-19-7"
    formula = "CH3COOH"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(351)
    Tc = unidades.Temperature(590.7)
    Pc = unidades.Pressure(5786, "kPa")
    M = 60.05196  # g/mol
    Tt = unidades.Temperature(289.8)
    Tb = unidades.Temperature(391.2)
    f_acent = 0.4665
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 130

    Fi1 = {"ao_log": [1, 3.6676653],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [-3.94616949, 5.48487930, -0.210687796, -0.781330239,
                      0.130979005],
           "titao": [2.09502491],
           "ao_exp": [6.28891793]}

    piazza = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for acetic acid of Piazza "
                    "(2011)",
        "__doi__": {
            "autor": "Piazza, L., Span, R.",
            "title": "An equation of state for acetic acid including the "
                     "association term of SAFT",
            "ref": "Fluid Phase Equilib. 303 (2011) 134-149",
            "doi": "10.1016/j.fluid.2011.01.008"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 540.0, "Pmax": 30000.0, "rhomax": 35.57,

        "nr1": [-1.5624834164583, -0.874703669570960, 4.6968858010355,
                0.0097367136204905, -0.0049055972708048],
        "d1": [1, 1, 2, 2, 6],
        "t1": [-1, 1.375, 1, 1.375, 0.75],

        "nr2": [24.499997808125, -31.443235067567, -1.3768156877983,
                1.4849435860881, 1.1374909453775, -2.6039791873344,
                -0.030484923493199, 5.3316386834696, -5.6733952193640,
                -0.126785566440530],
        "d2": [3, 3, 3, 4, 4, 5, 5, 5, 5, 2],
        "t2": [-0.25, 0, 2.25, 0.125, 2.125, 1.25, 2.25, 2.125, 2.375, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
        "gamma2": [1]*10,

        "type_ass": "1",
        "m_ass": 1.01871348,
        "v_ass": 0.0444215309,
        "k_ass": 0.109117041e-4,
        "e_ass": 12.2735737}

    eq = (piazza, )

    # TODO: Add ancillary equations
    # _vapor_Pressure = {
    #     "eq": 3,
    #     "n": [],
    #     "t": []}
    # _liquid_Density = {
    #     "eq": 1,
    #     "n": [],
    #     "t": []}
    # _vapor_Density = {
    #     "eq": 2,
    #     "n": [],
    #     "t": []}


class Test(TestCase):
    """Testing"""

    def test_Piazza(self):
        """Table 4, Pag 143"""

        # Tiny desviation in enthalpy and derived properties, probably last
        # decimal in ideal state integration parameters
        st = AceticAcid(T=290, rho=0)
        self.assertEqual(round(st.P.MPa, 8), 0)
        self.assertEqual(round(st.Z, 5), 1)
        # self.assertEqual(round(st.h.kJkg, 5), 669.46135)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.89579)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.03425)
        self.assertEqual(round(st.w, 5), 215.30851)

        st = AceticAcid(T=290, rho=0.01)
        self.assertEqual(round(st.P.MPa, 8), 0.23586e-3)
        self.assertEqual(round(st.Z, 5), 0.58743)
        self.assertEqual(round(st.h.kJkg, 5), 238.11545)
        self.assertEqual(round(st.s.kJkgK, 5), 1.05462)
        self.assertEqual(round(st.cv.kJkgK, 5), 4.64706)
        self.assertEqual(round(st.cp.kJkgK, 5), 5.35800)
        self.assertEqual(round(st.w, 5), 158.69277)

        st = AceticAcid(T=290, rho=0.025)
        self.assertEqual(round(st.P.MPa, 8), 0.55619e-3)
        self.assertEqual(round(st.Z, 5), 0.55408)
        self.assertEqual(round(st.h.kJkg, 5), 203.44937)
        self.assertEqual(round(st.s.kJkgK, 5), 0.86742)
        self.assertEqual(round(st.cv.kJkgK, 5), 3.46806)
        self.assertEqual(round(st.cp.kJkgK, 5), 3.91294)
        self.assertEqual(round(st.w, 5), 154.08506)

        st = AceticAcid(T=290, rho=1060)
        self.assertEqual(round(st.P.MPa, 4), 9.5458)
        self.assertEqual(round(st.Z, 5), 0.22429)
        self.assertEqual(round(st.h.kJkg, 5), -218.59533)
        self.assertEqual(round(st.s.kJkgK, 5), -0.66978)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.66179)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.99658)
        self.assertEqual(round(st.w, 5), 1221.38253)

        st = AceticAcid(T=290, rho=1070)
        self.assertEqual(round(st.P.MPa, 3), 22.580)
        self.assertEqual(round(st.Z, 5), 0.52557)
        self.assertEqual(round(st.h.kJkg, 5), -209.68926)
        self.assertEqual(round(st.s.kJkgK, 5), -0.68127)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.65651)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.98783)
        self.assertEqual(round(st.w, 5), 1278.95838)

        st = AceticAcid(T=450, rho=0)
        self.assertEqual(round(st.P.MPa, 8), 0)
        self.assertEqual(round(st.Z, 5), 1)
        self.assertEqual(round(st.h.kJkg, 5), 868.66332)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.31313)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.45159)
        self.assertEqual(round(st.w, 5), 262.43832)

        st = AceticAcid(T=450, rho=1)
        self.assertEqual(round(st.P.MPa, 6), 0.55569e-1)
        self.assertEqual(round(st.Z, 5), 0.89189)
        self.assertEqual(round(st.h.kJkg, 5), 760.34156)
        self.assertEqual(round(st.s.kJkgK, 5), 1.94166)
        self.assertEqual(round(st.cv.kJkgK, 5), 3.74492)
        self.assertEqual(round(st.cp.kJkgK, 5), 4.41190)
        self.assertEqual(round(st.w, 5), 244.91800)

        st = AceticAcid(T=450, rho=6)
        self.assertEqual(round(st.P.MPa, 5), 0.26514)
        self.assertEqual(round(st.Z, 5), 0.70925)
        self.assertEqual(round(st.h.kJkg, 5), 599.36356)
        self.assertEqual(round(st.s.kJkgK, 5), 1.40868)
        self.assertEqual(round(st.cv.kJkgK, 5), 4.43709)
        self.assertEqual(round(st.cp.kJkgK, 5), 5.53182)
        self.assertEqual(round(st.w, 5), 212.41680)

        st = AceticAcid(T=450, rho=870)
        self.assertEqual(round(st.P.MPa, 4), 2.7252)
        self.assertEqual(round(st.Z, 5), 0.05028)
        self.assertEqual(round(st.h.kJkg, 5), 152.00061)
        self.assertEqual(round(st.s.kJkgK, 5), 0.35592)
        self.assertEqual(round(st.cv.kJkgK, 5), 2.48256)
        self.assertEqual(round(st.cp.kJkgK, 5), 2.74525)
        self.assertEqual(round(st.w, 5), 536.47036)

        st = AceticAcid(T=450, rho=880)
        self.assertEqual(round(st.P.MPa, 4), 5.6865)
        self.assertEqual(round(st.Z, 5), 0.10372)
        self.assertEqual(round(st.h.kJkg, 5), 153.23344)
        self.assertEqual(round(st.s.kJkgK, 5), 0.35114)
        self.assertEqual(round(st.cv.kJkgK, 5), 2.46536)
        self.assertEqual(round(st.cp.kJkgK, 5), 2.73491)
        self.assertEqual(round(st.w, 5), 607.71388)
