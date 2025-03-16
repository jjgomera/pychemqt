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


from math import exp, log, log10
from unittest import TestCase

from lib import unidades
from lib.meos import MEoS


class H2(MEoS):
    """Multiparameter equation of state for hydrogen (normal)"""
    name = "hydrogen"
    CASNumber = "1333-74-0"
    formula = "H2"
    synonym = "R-702"
    _refPropName = "HYDROGEN"
    _coolPropName = "Hydrogen"
    rhoc = unidades.Density(31.26226704)
    Tc = unidades.Temperature(33.145)
    Pc = unidades.Pressure(1296.4, "kPa")
    M = 2.01588  # g/mol
    Tt = unidades.Temperature(13.957)
    Tb = unidades.Temperature(20.369)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4579856475, 1.888076782],
           "ao_exp": [1.616, -0.4117, -0.792, 0.758, 1.217],
           "titao": [16.0205159149, 22.6580178006, 60.0090511389,
                     74.9434303817, 206.9392065168]}

    Fi2 = {"ao_log": [1, 1.47906],
           "pow": [0, 1],
           "ao_pow": [13.796443393, -175.864487294],
           "ao_sinh": [0.95806, 1.56039], "sinh": [6.891654113, 49.76529075],
           "ao_cosh": [0.45444, -1.3756], "cosh": [9.84763483, 50.367279301]}

    CP1 = {"ao": 0.72480209e3,
           "an": [0.12155215e11, -0.36396763e10, 0.43375265e9, -0.23085817e8,
                  -0.38680927e4, 0.88240136e5, -0.78587085e4, -0.18426806e3,
                  0.21801550e2, -0.13051820e1, 0.21003175e-1, 0.23911604e-2,
                  -0.18240547e-3, 0.56149561e-5, -0.73803310e-7,
                  0.66357755e-11],
           "pow": [-7, -6, -5, -4, -3, -2, -1.001, 0.5, 1, 1.5, 2, 2.5, 3, 3.5,
                   4, 5]}

    leachman = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for normal hydrogen of "
                    "Leachman et al. (2009).",
        "__doi__": {"autor": "Leachman, J.W., Jacobsen, R.T, Penoncello, S.G.,"
                             " Lemmon, E.W.",
                    "title": "Fundamental equations of state for Parahydrogen,"
                             " Normal Hydrogen, and Orthohydrogen",
                    "ref": "J. Phys. Chem. Ref. Data, 38(3) (2009) 721-748",
                    "doi": "10.1063/1.3160306"},

        "R": 8.314472,
        "Tc": 33.145, "Pc": 1296.4, "rhoc": 15.508,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 102.0,

        "nr1": [-6.93643, 0.01, 2.1101, 4.52059, 0.732564, -1.34086, 0.130985],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.6844, 1., 0.989, 0.489, 0.803, 1.1444, 1.409],

        "nr2": [-0.777414, 0.351944],
        "d2": [1, 3],
        "t2": [1.754, 1.311],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.0211716, 0.0226312, 0.032187, -0.0231752, 0.0557346],
        "d3": [2, 1, 3, 1, 1],
        "t3": [4.187, 5.646, 0.791, 7.249, 2.986],
        "alfa3": [1.685, 0.489, 0.103, 2.506, 1.607],
        "beta3": [0.171, 0.2245, 0.1304, 0.2785, 0.3967],
        "gamma3": [0.7164, 1.3444, 1.4517, 0.7204, 1.5445],
        "epsilon3": [1.506, 0.156, 1.736, 0.67, 1.662]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 38.148,

        "nr1": [0.53579928451252e1, -0.62050252530595e1, 0.13830241327086,
                -0.71397954896129e-1, 0.15474053959733e-1],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.5, 0.625, 0.384, 0.625, 1.125],

        "nr2": [-0.14976806405771, -0.26368723988451e-1, 0.56681303156066e-1,
                -0.60063958030436e-1, -0.45043942027132, 0.42478840244500,
                -0.021997640827139, -0.01049952137453, -0.28955902866816e-2],
        "d2": [1, 5, 5, 5, 1, 1, 2, 5, 1],
        "t2": [2.625, 0, 0.25, 1.375, 4, 4.25, 5, 8, 8],
        "c2": [1, 1, 1, 1, 2, 2, 3, 3, 5],
        "gamma2": [1]*9}

    bender = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen of Bender "
                    "(1982).",
        "__doi__": {"autor": "Bender, E.",
                    "title": "Equation of state of normal hydrogen in the "
                             "range 18 to 700 K and 1 to 500 bar",
                    "ref": "VDI-Forschungsheft, no. 609, 1982, p. 15-20",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 18.0, "Tmax": 700.0, "Pmax": 50000.0, "rhomax": 38.74,

        "nr1": [0.133442326203e1, -0.104116843433e1, 0.227202245707,
                0.300374270906, -0.463984214813, -0.178010492282e1,
                0.100460103605e1, -0.187200622541, 0.980276957749e-2,
                0.543224866339e-1, -0.263496312610e-1, 0.315432315759e-1,
                -0.525788294155e-1, -0.685380627808e-2, 0.344540276656e-1,
                -0.555747275982e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.133442326203e1, 0.104116843433e1, -0.227202245707,
                -0.378598758038, 0.249888797892, -0.498847982876e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.711139834571]*6,

        "nr3": [],
        "nr4": []}

    eq = leachman, GERG, bender

    _surface = {"sigma": [-1.4165, 0.746383, 0.675625],
                "exp": [0.63882, 0.659804, 0.619149]}
    _dielectric = {
        "eq": 1,
        "a": [2.0306, 0.0056], "b": [0.181, 0.021], "c": [-7.4, 0],
        "Au": 0, "D": 2}

    _melting = {
        "eq": 1,
        "__doi__": {
            "autor": "Datchi, F., Loubeyre, P., LeToullec, R.",
            "title": "Extended and accuracy determination of the melting "
                     "curves of argon, helium, ice (H2O), and hydrogen (H2)",
            "ref": "Physical Review B, 61(10) (2000) 6535-6546",
            "doi": "10.1103/PhysRevB.61.6535"},

        # It used the Eq. 9 in paper, the Eq. 10 is only prefered at very high
        # pressure, solved only iteratively, out of interest range in pychemqt
        "Tmin": Tt, "Tmax": 1100,
        "Tref": 1, "Pref": 1e9,
        "a1": [1.63e-4], "exp1": [1.824]}

    _sublimation = {
        "__doi__": {
            "autor": "McCarty, R.D., Hord, J., Roder, H.M.",
            "title": "Selected Properties of Hydrogen (Engineering Design "
                     "Data)",
            "ref": "NBS Monograph 168, NBS 1981.",
            "doi": ""},
        "Tmin": Tt, "Tmax": Tt}

    @classmethod
    def _Sublimation_Pressure(cls, T):
        """Special sublimation pressure correlation"""
        # Use decimal logarithm
        P = 10**(-43.39/T+2.5*log10(T)+2.047)
        return unidades.Pressure(P, "mmHg")

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.489789e1, 0.988558, 0.349689, 0.499356],
        "t": [1.0, 1.5, 2.0, 2.85]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.15456e2, -0.41720e2, 0.50276e2, -0.27947e2, 0.56718e1],
        "t": [0.62, 0.83, 1.05, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.9962, -16.724, 15.819, -16.852, 34.586, -53.754],
        "t": [0.466, 2, 2.4, 4., 7., 8.]}

    visco0 = {"__name__": "Muzny (2013)",
              "__doi__": {
                  "autor": "Muzny, C.D., Huber, M.L., Kazakov, A.F.",
                  "title": "Correlation for the Viscosity of Normal Hydrogen "
                           "Obtained from Symbolic Regression",
                  "ref": "J. Chem. Eng. Data 58(4) (2013) 969-979",
                  "doi": "10.1021/je301273j"},

              "eq": 1, "omega": 1,
              "ek": 30.41, "sigma": 0.297,
              "n_chapman": 0.021357,
              "collision": [0.20963, -0.455274, 0.143602, -0.0335325,
                            0.00276981],

              "Tref_virial": 30.41,
              "n_virial": [-0.187, 2.4871, 3.7151, -11.0972, 9.0965, -3.8292,
                           0.5166],
              "t_virial": [0, -1, -2, -3, -4, -5, -6],

              "special": "_mur",
              }

    def _mur(self, rho, T, fase):
        # Simbolic Regression, Eq. 9

        # Update scale factor for density and add test from erratum paper
        # Muzny, C.D., Huber, M.L., Kazakov, A.F.
        # Erratum: Correlation for the Viscosity of Normal Hydrogen Obtained
        # from Symbolic Regression
        # J. Chem. Eng. Data 67(9) (2022) 2855
        # doi: 10.1021/acs.jced.2c00523

        rhor = rho/90.909090909
        Tr = T/self.Tc
        c = [6.43449673, 4.56334068e-2, 2.32797868e-1, 9.5832612e-1,
             1.27941189e-1, 3.63576595e-1]
        nr = c[0]*rhor**2*exp(c[1]*Tr + c[2]/Tr + c[3]*rhor**2/(c[4]+Tr)
                              + c[5]*rhor**6)

        return nr

    visco1 = {"__name__": "McCarty (1972)",
              "__doi__": {
                  "autor": "McCarty, R.D., Weber, L.A.",
                  "title": "Thermophysical Properties of Parahydrogen from "
                           "the Freezing Liquid Line to 5000 R for Pressures "
                           "to 10000 psia",
                  "ref": "NBS Technical Note 617",
                  "doi": ""},

              "eq": 0,
              "method": "_visco1"}

    def _visco1(self, rho, T, fase):
        def muo(T):
            b = [-0.1841091042788e2, 0.3185762039455e2, -0.2308233586574e2,
                 0.9129812714730e1, -0.2163626387630e1, 0.3175128582601,
                 -0.2773173035271e-1, 0.1347359367871e-2, -0.2775671778154e-4]
            suma = 0
            for i, b in enumerate(b):
                suma += b*T**((-3.+i)/3)
            return suma*100

        def mu1(rho, T):
            if not rho:
                return 0
            A = exp(5.7694 + log(rho.gcc) + 0.65e2*rho.gcc**1.5
                    - 6e-6*exp(127.2*rho.gcc))
            B = 10 + 7.2*((rho.gcc/0.07)**6-(rho.gcc/0.07)**1.5) - \
                17.63*exp(-58.75*(rho.gcc/0.07)**3)
            return A*exp(B/T)*0.1

        def mu2(rho, T):
            c = [-0.1324266117873e2, 0.1895048470537e2, 0.2184151514282e2,
                 0.9771827164811e5, -0.1157010275059e4, 0.1911147702539e3,
                 -0.3186427506942e4, 0.0705565000000]
            R2 = rho.gcc**0.5*(rho.gcc-c[7])/c[7]
            A = c[0] + c[1]*R2 + c[2]*rho.gcc**0.1 + c[3]*R2/T**2 + \
                c[4]*rho.gcc**0.1/T**1.5 + c[5]/T + c[6]*R2/T
            B = c[0]+c[5]/T
            return 0.1*(exp(A)-exp(B))

        def mur(rho1, T1, rho2, T2):
            return muo(T1) + mu2(rho1, T2) - muo(T2) - mu2(rho2, T2)

        if T > 100:
            mu = muo(100) + mu1(rho, 100) + mur(rho, T, rho, 100)
        else:
            mu = muo(T)+mu1(rho, T)
        return unidades.Viscosity(mu, "muPas")

    visco2 = {"__name__": "Vargaftik (1996)",
              "__doi__": {
                  "autor": "Vargaftik, N.B., Vinogradov, Y.K., Yargin, V.S.",
                  "title": "Handbook of Physical Properties of Liquids and "
                           "Gases",
                  "ref": "Begell House, New York, 1996",
                  "doi": ""},

              "eq": 1, "omega": 1,

              "ek": 59.7, "sigma": 0.2827,
              "Tref_virial": 32.938, "etaref_virial": 1.*M,
              "n_virial": [-2.1505e-1, 10.727e-1, -16.935e-1, 0.0, 22.702e-1,
                           2.2123e-1, 0.34163e-1, -0.043206e-1],
              "t_virial": [-1.5, -1, -0.5, 0, 0.5, 1.5, 2.],

              "Tref_res": 32.938, "rhoref_res": 15.556*M,
              "nr": [-0.922703, 6.41602, -5.98018, 0.289715, 2.36429,
                     -0.27887, -11.0595, 11.1582, 7.18928, -7.76971, -1.21827,
                     1.47193],
              "tr": [0, -1, -2, -3, 0, 0, -1, -2, -1, -2, -1, -2],
              "dr": [1, 1, 1, 1, 2, 3, 3, 3, 4, 4, 5, 5]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Assael (2011)",
               "__doi__": {
                   "autor": "Assael, M.J., Assael. J.-A.M., Huber, M.L., "
                            "Perkins, R.A., Takata, Y.",
                   "title": "Correlation of the Thermal Conductivity of "
                            "Normal and Parahydrogen from the Triple Point to "
                            "1000 K and up to 100 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 40(3) (2011) 033101",
                   "doi": "10.1063/1.3606499"},

               "eq": 1,
               "rhoc": 31.262,

               "Toref": 33.145, "koref": 1,
               "no_num": [-3.40976e-1, 4.58820, -1.45080, 3.26394e-1,
                          3.16939e-3, 1.90592e-4, -1.139e-6],
               "to_num": [0, 1, 2, 3, 4, 5, 6],
               "no_den": [1.38497e2, -2.21878e1, 4.57151, 1],
               "to_den": [0, 1, 2, 3],

               "Tref_res": 33.145, "rhoref_res": 31.262, "kref_res": 1.,
               "nr": [3.63081e-2, -2.07629e-2, 3.1481e-2, -1.43097e-2,
                      1.7498e-3, 1.8337e-3, -8.86716e-3, 1.5826e-2,
                      -1.06283e-2, 2.80673e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 0.4e-9, "Tcref": 49.7175}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""
    def test_leachman(self):
        """Selected point from Table 14, Pag 746, saturation states"""
        st = H2(T=13.957, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 7.358)
        self.assertEqual(round(st.Liquido.rho, 3), 77.004)
        self.assertEqual(round(st.Gas.rho, 5), 0.12985)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -53.926)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 399.83)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -3.0723)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 29.438)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.1616)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.2433)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 7.0212)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 10.564)
        self.assertEqual(round(st.Liquido.w, 1), 1269.2)
        self.assertEqual(round(st.Gas.w, 2), 307.14)

        st = H2(T=20, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 90.717)
        self.assertEqual(round(st.Liquido.rho, 3), 71.265)
        self.assertEqual(round(st.Gas.rho, 4), 1.2059)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), -3.6672)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 446.64)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.17429)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 22.341)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.6369)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.4343)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 9.5697)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 11.892)
        self.assertEqual(round(st.Liquido.w, 1), 1129.1)
        self.assertEqual(round(st.Gas.w, 2), 354.31)

        st = H2(T=30, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 804.32)
        self.assertEqual(round(st.Liquido.rho, 3), 54.538)
        self.assertEqual(round(st.Gas.rho, 3), 10.445)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 140.30)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 441.19)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 5.0661)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 15.096)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 6.3535)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 7.4945)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 25.284)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 30.425)
        self.assertEqual(round(st.Liquido.w, 2), 713.07)
        self.assertEqual(round(st.Gas.w, 2), 379.07)

        st = H2(T=33, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 1269.3)
        self.assertEqual(round(st.Liquido.rho, 3), 38.079)
        self.assertEqual(round(st.Gas.rho, 3), 24.637)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 255.69)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 343.40)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 8.3842)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 11.042)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 7.6982)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 8.5381)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 484.58)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 604.72)
        self.assertEqual(round(st.Liquido.w, 2), 423.53)
        self.assertEqual(round(st.Gas.w, 2), 373.93)

    def test_Assael(self):
        """Table 7, Pag 11"""
        kw = {"visco": 1}
        self.assertEqual(round(H2(T=298.15, rho=0, **kw).k.mWmK, 2), 185.67)
        self.assertEqual(round(
            H2(T=298.15, rho=0.80844, **kw).k.mWmK, 2), 186.97)
        self.assertEqual(round(
            H2(T=298.15, rho=14.4813, **kw).k.mWmK, 2), 201.35)
        self.assertEqual(round(H2(T=35, rho=0, **kw).k.mWmK, 3), 26.988)
        self.assertEqual(round(H2(T=35, rho=30, **kw).k.mWmK, 3), 75.595)
        self.assertEqual(round(H2(T=18, rho=0, **kw).k.mWmK, 3), 13.875)
        self.assertEqual(round(H2(T=18, rho=75, **kw).k.mWmK, 2), 104.48)

    def test_Muzny(self):
        """Table 1 from erratum article"""
        self.assertEqual(round(H2(T=40, rho=0).mu.muPas, 4), 1.9772)
        self.assertEqual(round(H2(T=40, rho=50).mu.muPas, 4), 5.9905)
        self.assertEqual(round(H2(T=40, rho=100).mu.muPas, 3), 49.034)
