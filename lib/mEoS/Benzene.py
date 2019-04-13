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

from lib import unidades
from lib.meos import MEoS


class Benzene(MEoS):
    """Multiparameter equation of state for benzene"""
    name = "benzene"
    CASNumber = "71-43-2"
    formula = "C6H6"
    synonym = ""
    _refPropName = "BENZENE"
    _coolPropName = "Benzene"
    rhoc = unidades.Density(304.7922436)
    Tc = unidades.Temperature(562.02)
    Pc = unidades.Pressure(4894, "kPa")
    M = 78.11184  # g/mol
    Tt = unidades.Temperature(278.674)
    Tb = unidades.Temperature(353.22)
    f_acent = 0.211
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 40

    Fi1 = {"ao_log": [1, 2.94645],
           "pow": [0, 1],
           "ao_pow": [-0.6740687105, 2.5560188958],
           "ao_exp": [7.36374, 18.6490, 4.01834],
           "titao": [4116/Tc, 1511/Tc, 630/Tc]}

    CP2 = {"ao": -0.478176/8.3143*78.108,
           "an": [0.618649e-2/8.3143*78.108, -0.380363e-5/8.3143*78.108,
                  0.699648e-9/8.3143*78.108, 0.42661e-13/8.3143*78.108],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    thol = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Thol (2013).",
        "__doi__": {"autor": "Thol, M., Lemmon, E.W., Span, R.",
                    "title": "Equation of State for Benzene for Temperatures "
                             "from the Melting Line up to 750 K with Pressures"
                             " up to 500 MPa",
                    "ref": "High Temperatures-High Pressures 41 (2012) 81-97",
                    "doi": ""},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 725, "Pmax": 500000.0, "rhomax": 11.45,

        "nr1": [0.03513062, 2.229707, -3.100459, -0.5763224, 0.2504179],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.744, 1.174, 0.68],

        "nr2": [-0.7049091, -0.1393433, 0.8319673, -0.3310741, -0.02793578],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.5, 3.67, 1.26, 2.6, 0.95],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1.]*5,

        "nr3": [0.7087408, -0.3723906, -0.06267414, -0.86295],
        "d3": [1, 1, 3, 3],
        "t3": [1, 2.47, 3.35, 0.75],
        "alfa3": [1.032, 1.423, 1.071, 14.35],
        "beta3": [1.867, 1.766, 1.824, 297.5],
        "gamma3": [1.1180, 0.6392, 0.6536, 1.1640],
        "epsilon3": [0.7289, 0.9074, 0.7655, 0.8711],

        "nr4": []}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Polt (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 278.7, "Tmax": 635.0, "Pmax": 78000.0, "rhomax": 11.45,

        "nr1": [-0.918572178424, 0.155357491575e1, -0.356149241161,
                0.817273664265, -0.331303917534e1, 0.335336626528e1,
                -0.256976312022e1, 0.427304812515, 0.406483484297,
                -0.329744378187, 0.208907540720, 0.777471199254e-1,
                -0.202621443063, -0.148580350700e-1, 0.503167715817e-1,
                0.293012717053e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.918572178424, -0.155357491575e1, 0.356149241161,
                -0.447029533153e-1, 0.957712367542, -0.114688433057e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.95481]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [1.76284970, 1.02610647, -3.74263321, 9.57682041e-2,
                2.59179321e-4, -1.03082188e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [1.07359246e-1, -1.12562310e-1, 3.18737987e-1, -3.07549016e-2,
                -3.25082386e-1, 2.28099159e-2, -7.07431076e-2, -1.96809158e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = thol, polt, sun

    _surface = {"sigma": [0.07298, -0.0007802, -0.0001756],
                "exp": [1.232, 0.8635, 0.3065]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.71661e1, 0.21551e1, -0.20297e1, -0.40668e1, 0.38092],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.18160e2, -0.56879e2, 0.87478e2, -0.64365e2, 0.18500e2],
        "t": [0.534, 0.686, 0.84, 1.0, 1.2]}

    _vapor_Density = {
        "eq": 2,
        "n": [-3.1147, -4.6689, -16.161, -146.5, 518.87, -827.72],
        "t": [0.419, 1.12, 2.8, 7.3, 10, 12]}

    visco0 = {"__name__": "Avgeri (2014)",
              "__doi__": {
                  "autor": "Avgeri, S., Assael, M.J., Huber, M.L., Perkins, "
                           "R.A",
                  "title": "Reference Correlation of the Viscosity of Benzene "
                           "from the Triple Point to 675K and up to 300MPa",
                  "ref": "J. Phys. Chem. Ref. Data 43 (2014) 033103",
                  "doi": "10.1063/1.4892935"},

              "eq": 1, "omega": 1,

              "n_chapman": 0.021357,
              "ek": 412,
              "sigma": 0.54,
              "collision": [0.234018, -0.476136, 0, -0.015269],

              "Tref_virial": 412,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "special": "_mur"}

    def _mur(self, rho, T, fase):
        """Special form for residual term of Avgeri viscosity correlation, eq 8
        in paper"""
        Tr = T/self.Tc
        rhor = rho/self.rhoc
        c = [-9.98945, 86.0626, 2.74872, 1.11130, -1, -134.133, -352.473,
             6.60989, 88.4174]
        F = c[0]*rhor**2 + c[1]*rhor/(c[2]+c[3]*Tr+c[4]*rhor) + \
            (c[5]*rhor+c[6]*rhor**2)/(c[7]+c[8]*rhor**2)
        mur = rhor**(2/3)*Tr**0.5*F
        return mur

    _viscosity = visco0,

    thermo0 = {"__name__": "Assael (2012)",
               "__doi__": {
                   "autor": "Assael, M.J., Mihailidou, E., Huber, M.L., "
                            "Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal "
                            "Conductivity of Benzene from the Triple Point to "
                            "725 K and up to 500 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 41(4) (2012) 043102",
                   "doi": "10.1063/1.4755781"},

               "eq": 1,

               "Toref": 562.02, "koref": 1e-3,
               "no_num": [101.404, -521.44, 868.266],
               "to_num": [0, 1, 2],
               "no_den": [1, 9.714, 1.467],
               "to_den": [0, 1, 2],

               "Tref_res": 562.02, "rhoref_res": 304.792, "kref_res": 1,
               "nr": [2.82489e-2, -7.73415e-2, 7.14001e-2, -2.36798e-2,
                      3.00875e-3, -1.19268e-2, 8.33389e-2, -8.98176e-2,
                      3.63025e-2, -4.90052e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 2.16e-10, "gam0": 0.0569, "qd": 6.2e-10, "Tcref": 843}

    _thermal = thermo0,


class Test(TestCase):

    def test_viscoAvgeri(self):
        # Table 8, pag 11
        self.assertEqual(round(Benzene(T=300, rho=0).mu.muPas, 3), 7.625)
        self.assertEqual(round(Benzene(T=400, rho=0).mu.muPas, 3), 10.102)
        self.assertEqual(round(Benzene(T=550, rho=0).mu.muPas, 3), 13.790)
        self.assertEqual(round(Benzene(T=300, rho=875).mu.muPas, 2), 608.53)
        self.assertEqual(round(Benzene(T=400, rho=760).mu.muPas, 2), 211.74)
        self.assertEqual(round(Benzene(T=550, rho=500).mu.muPas, 3), 60.511)

        # Table 5, Saturation line, give too liquid density values to check
        # Thol mEoS
        st = Benzene(T=280, x=0)
        self.assertEqual(round(st.Liquido.rho, 3), 892.701)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 795.3)

        st = Benzene(T=400, x=0)
        self.assertEqual(round(st.Liquido.rho, 3), 758.649)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 209.9)

        st = Benzene(T=540, x=0)
        self.assertEqual(round(st.Liquido.rho, 3), 508.838)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 62.26)

    def test_thermoAssael(self):
        # Table 4, pag 8
        # The viscosity correlation used in the paper is different to Avgeri
        # visco0 correlation, so tiny desviations are for that cause
        self.assertEqual(round(Benzene(T=290, rho=890).k.mWmK, 2), 147.65)
        self.assertEqual(round(Benzene(T=500, rho=2).k.mWmK, 3), 30.174)
        self.assertEqual(round(Benzene(T=500, rho=32).k.mWmK, 3), 32.209)
        self.assertEqual(round(Benzene(T=500, rho=800).k.mWmK, 2), 141.22)
        self.assertEqual(round(Benzene(T=570, rho=1.7).k.mWmK, 3), 37.763)
