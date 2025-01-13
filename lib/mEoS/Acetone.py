#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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
from lib.mEoS import C3


class Acetone(MEoS):
    """Multiparameter equation of state for Acetone"""
    name = "acetone"
    CASNumber = "67-64-1"
    formula = "CH3COCH3"
    synonym = ""
    _refPropName = "ACETONE"
    _coolPropName = "Acetone"
    rhoc = unidades.Density(272.971958)
    Tc = unidades.Temperature(508.1)
    Pc = unidades.Pressure(4700.0, "kPa")
    M = 58.07914  # g/mol
    Tt = unidades.Temperature(178.5)
    Tb = unidades.Temperature(329.22)
    f_acent = 0.3071
    momentoDipolar = unidades.DipoleMoment(2.88, "Debye")
    id = 140

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-9.4883659997, 7.1422719708],
           "ao_exp": [3.7072, 7.0675, 11.012],
           "titao": [310/Tc, 3480/Tc, 1576/Tc]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for acetone of Lemmon "
                    "and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 700000.0, "rhomax": 15.73,

        "nr1": [0.90041, -2.1267, -0.083409, 0.065683, 0.00016527],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.039663, 0.72085, 0.0092318, -0.17217, -0.14961, -0.076124,
                -0.018166],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = (lemmon, )

    _surface = {"sigma": [0.0633], "exp": [1.16]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.76214e1, 0.17441e1, -0.20514e1, -0.26644e1, -0.69437],
        "t": [1, 1.5, 2.57, 4.43, 15.]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.11118e2, -0.29507e2, 0.35255e2, -0.14712e2, 0.95560],
        "t": [0.456, 0.626, 0.8, 1.0, 2.47]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.25200e1, -.66065e1, -.25751e2, .78120e1, -.53778e2, -116.84],
        "t": [0.36, 1.05, 3.2, 4.0, 6.5, 14.0]}

    visco0 = {
        "__name__": "Sotiriadou (2025)",
        "__doi__": {
            "autor": "Sotiriadou, S.G., Ntonti, E., Assael, M.J., Huber, M.L.",
            "title": "Reference Correlation for the Viscosity and Thermal "
                     "Conductivity of Acetone from the Triple Point to High "
                     "Temperatures and Pressures",
            "ref": "Int. J. Thermophys. 46(1) (2025) 3",
            "doi": "10.1007/s10765-024-03465-6"},

        "eq": 1, "omega": 0,
        "sigma": 0.49,

        "Toref": Tc,
        "no_num": [0.931015, 13.4773, -6.84412, 3.30874, 4.78248, -1.45555,
                   0.149281],
        "to_num": [0, 1, 2, 3, 4, 5, 6],
        "no_den": [1.46335, -1.36059, 1],
        "to_den": [0, 1, 2],

        "Tref_virial": 432,
        "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                     -3375.1717, 2491.6597, -787.26086, 14.085455,
                     -0.34664158],
        "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

        "special": "_mur"}

    def _mur(self, rho, T, fase):
        """Special term of residual viscosity for Sotiriadou correlation"""
        Tr = T/self.Tc
        rhor = rho/self.rhoc

        # Eq 7
        mur = rhor**(2/3)*Tr**0.5 * (
            6.2435628350 * rhor
            + (0.16610522013 + 8.9088278828*rhor + 0.16610522013*rhor**5
               + 0.0069857927082*Tr**2*rhor**8)/(Tr -0.088521102246*rhor))
        return mur

    thermo0 = {
        "__name__": "Sotiriadou (2025)",
        "__doi__": {
            "autor": "Sotiriadou, S.G., Ntonti, E., Assael, M.J., Huber, M.L.",
            "title": "Reference Correlation for the Viscosity and Thermal "
                     "Conductivity of Acetone from the Triple Point to High "
                     "Temperatures and Pressures",
            "ref": "Int. J. Thermophys. 46(1) (2025) 3",
            "doi": "10.1007/s10765-024-03465-6"},

        "eq": 1,

        "Toref": Tc, "koref": 1e-3,
        "no_num": [-5.98797, 46.9565, -149.748, 241.207, -43.1278, 3.52029],
        "to_num": [0, 1, 2, 3, 4, 5],
        "no_den": [-0.614176, 2.57584, 1],
        "to_den": [0, 1, 2],

        "Tref_res": Tc, "rhoref_res": rhoc, "kref_res": 1e-3,
        "nr": [149.900120, -59.8846154, -223.191952, 117.823591, 130.528948,
               -70.7635055, -29.0922187, 14.9380192, 2.14538883, -0.506124251],
        "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
        "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

        "critical": 3,
        "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
        "Xio": 0.196e-9, "gam0": 0.052, "qd": 0.586e-9, "Tcref": 762.15}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 519, "sigma": 0.4669, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [1.25183, -0.239533, 0.0485815], "psi_d": [0, 1, 2],
              "fint": [0.954299e-3, 0.522303e-6], "fint_t": [0, 1],
              "chi": [1.08482, -0.0313081], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.196e-9, "gam0": 0.052, "qd": 0.586e-9, "Tcref": 1.5*Tc}

    _viscosity = (visco0, trnECS)
    _thermal = (thermo0, trnECS)


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = Acetone(T=510, rhom=4)
        self.assertEqual(round(st.P.kPa, 3), 4807.955)
        self.assertEqual(round(st.hM.kJkmol, 3), 51782.004)
        self.assertEqual(round(st.sM.kJkmolK, 3), 157.331)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 138.449)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3766.619)
        self.assertEqual(round(st.w, 3), 125.351)

    def test_Sotiriadou(self):
        """Checking values given in section 4.2"""
        st = Acetone(T=300, rho=0)
        self.assertEqual(round(st.mu.muPas, 4), 7.6011)
        self.assertEqual(round(st.k.mWmK, 3), 11.306)

        st = Acetone(T=300, rho=785)
        self.assertEqual(round(st.mu.muPas, 2), 309.65)
        self.assertEqual(round(st.k.mWmK, 2), 157.63)

    def test_Huber(self):
        """Table 7, pag 266"""
        st = Acetone(T=457.3, rhom=9.8, visco=1, thermal=1)
        # self.assertEqual(round(st.mu.muPas, 5), 99.01729)
        # self.assertEqual(round(st.k.mWmK, 4), 96.5053)
        self.assertEqual(round(st.mu.muPas, 5), 99.02148)
        self.assertEqual(round(st.k.mWmK, 4), 96.5060)
