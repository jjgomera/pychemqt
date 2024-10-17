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


from math import exp
from unittest import TestCase

from scipy.constants import pi, Avogadro

from lib import unidades
from lib.meos import MEoS


class Methanol(MEoS):
    """Multiparameter equation of state for methanol"""
    name = "Methanol"
    CASNumber = "67-56-1"
    formula = "CH3OH"
    synonym = ""
    _refPropName = "METHANOL"
    _coolPropName = "Methanol"
    rhoc = unidades.Density(275.5626)
    Tc = unidades.Temperature(512.6)
    Pc = unidades.Pressure(8103.5, "kPa")
    M = 32.04216  # g/mol
    Tt = unidades.Temperature(175.61)
    Tb = unidades.Temperature(337.632)
    f_acent = 0.5625
    momentoDipolar = unidades.DipoleMoment(1.7, "Debye")
    id = 117

    Fi1 = {"ao_log": [1, 3.1950423807804],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [13.9864114647, 3.2006369296e3, -1.14289818828912e-3,
                      -2.62687155181005e-7, 6.42610441977784e-11],
           "titao": [3.7664265756],
           "ao_exp": [4.70118076896145]}

    CP1 = {"ao": 3.9007912,
           "ao_exp": [0.10992677e2, 0.18336830e2, -0.16366004e2, -0.62332348e1,
                      0.28035363e1, 0.10778099e1, 0.96965697],
           "exp": [2115.01542, 1676.18569, 1935.16717, 1504.97016, 4222.83691,
                   5296.17127, 273.36934]}

    CP2 = {"ao": 0.964220/8.3143*32.,
           "an": [0.532325e-4/8.3143*32., 0.672819e-5/8.3143*32.,
                  -0.768411e-8/8.3143*32., 0.275220e-11/8.3143*32.],
           "pow": [1, 2, 3, 4]}

    reuck = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of de Reuck and "
                    "Craven (1993)",
        "__doi__": {
            "autor": "de Reuck, K.M., Craven, R.J.B.",
            "title": "Methanol, International Thermodynamic Tables of the "
                     "Fluid State - 12",
            "ref": "IUPAC, Blackwell Scientific Publications, London, 1993.",
            "doi": ""},

        "R": 8.31448,
        "cp": CP1,
        "ref": "NBP",

        "Tc": 513.38, "rhoc": 8.78517,
        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 35.57,

        "nr1": [-0.280062505988e1, 0.125636372418e2, -0.130310563173e2,
                0.326593134060e1, -0.411425343805e1, 0.346397741254e1,
                -0.836443967590e-1, -0.369240098923, 0.313180842152e-2,
                0.603201474111, -0.231158593638, 0.106114844945,
                -0.792228164995e-1, -0.422419150975e-4, 0.758196739214e-2,
                -0.244617434701e-4, 0.115080328802e-5],
        "d1": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 7],
        "t1": [0, 1, 2, 3, 1, 2, 3, 4, 6, 0, 3, 4, 0, 7, 1, 6, 7],

        "nr2": [-0.125099747447e2, 0.270392835391e2, -0.212070717086e2,
                0.632799472270e1, 0.143687921636e2, -0.287450766617e2,
                0.185397216068e2, -0.388720372879e1, -0.416602487963e1,
                0.529665875982e1, 0.509360272812, -0.330257604839e1,
                -0.311045210826, 0.273460830583, 0.518916583979,
                -0.227570803104e-2, 0.211658196182e-1, -0.114335123221e-1,
                0.249860798459e-2],
        "d2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 5, 5, 5, 6, 9, 6, 6, 4],
        "t2": [1, 2, 3, 4, 1, 2, 3, 5, 1, 2, 1, 2, 4, 5, 2, 5, 9, 14, 19],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 6],
        "gamma2": [1.01733510223052]*16 + [1.03497071023039]*2
        + [1.05291203329783],

        "nr3": [-0.819291988442e1, 0.478601004557, -0.444161392885,
                0.179621810410, -0.687602278259, 0.240459848295e1,
                -0.688463987466e1, 0.113992982501e1],
        "d3": [1, 1, 1, 1, 1, 3, 3, 3],
        "t3": [0]*8,
        "alfa3": [4.06934040892209, 8.20892015621185, 9.15601592007471,
                  83.8326275286616, 16.2773616356884, 27.705105527215,
                  16.2773616356884, 264.95250181898],
        "beta3": [-3.8940745646517, -3.8940745646517, -3.8940745646517,
                  -3.8940745646517, -3.8940745646517, -23.0649031906293,
                  -23.0649031906293, -23.0649031906293],
        "gamma3": [1.54080254509371, 1.54080254509371, 1.54080254509371,
                   1.54080254509371, 1.54080254509371, 1.08389789427588,
                   1.08389789427588, 1.08389789427588],
        "epsilon3": [0]*8,
        "exp1": [2, 3, 2, 4, 2, 3, 2, 4],
        "exp2": [1]*8}

    piazza = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for acetic acid of Piazza "
                    "(2013)",
        "__doi__": {
            "autor": "Piazza, L., Span, R.",
            "title": "An equation of state for methanol including the "
                     "association term of SAFT",
            "ref": "Fluid Phase Equilib. 349 (2013) 12-24",
            "doi": "10.1016/j.fluid.2013.03.024"},

        "R": 8.314472,
        # Error in Cp ideal gas parameters in paper
        # "cp": Fi1,
        "cp": CP1,
        "ref": "OTO",
        "Tc": 512.5, "rhoc": 273/32.04186, "M": 32.04186,

        "Tmin": Tt, "Tmax": 540.0, "Pmax": 30000.0, "rhomax": 35.57,

        "nr1": [0.096352729792779, -1.0848826325874, 0.029919647090261,
                -0.0017963419593895, 0.000047354317752015],
        "d1": [1, 1, 4, 5, 7],
        "t1": [-0.125, 1.5, 0, -0.875, 1.25],

        "nr2": [1.0013578850486, -1.2555691488591, 0.85469725717500,
                -0.058295570793694, 0.026935675584229, 0.11504892676606,
                -0.0051081766133636, 0.0019167368789348, -0.28618221186953,
                0.48168213019845, -0.33081091251828, 0.092842083313630,
                -0.035936470747247],
        "d2": [1, 1, 3, 4, 5, 1, 7, 9, 2, 3, 4, 6, 7],
        "t2": [0.25, 2, 1.75, 2.5, 2.375, 6.875, 5.875, 5, 18.5, 19, 17.5, 14, 12],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
        "gamma2": [1]*13,

        "type_ass": "2B",
        "m_ass": 0.977118832,
        "v_ass": 0.204481952,
        "k_ass": 0.148852832e-2,
        "e_ass": 5.46341463}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [-5.24578394, 1.39060027, 8.56114069e-1, -4.20843418e-2,
                3.63682442e-5, 7.05598662e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [3.70573369e-1, 2.46303468, 1.50253790, 7.47553687e-2,
                -3.06417876e-1, -7.48402758e-1, -1.01432849e-1, 8.06830693e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 298., "Tmax": 703.0, "Pmax": 63000.0, "rhomax": 26.0625,

        "nr1": [-4.12043979985, 5.41210456547, -0.974639417666,
                -0.909437999343, -0.143467597275, 5.57052459597,
                -6.97445416557, 0.860535902136, 2.44117735035, -4.49073510921,
                2.23855290012, -0.71733653794, 0.876135006507, 0.151777405466,
                -0.233178058896, 0.0140022534721],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.412043979985e1, -0.541210456547e1, 0.974639417666,
                -0.4642672133, 0.944015617353, -0.449348200461],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.591872]*6}

    eq = reuck, piazza, sun, polt

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Sanjuán, E.L.",
                    "title": "Surface Tension of Alcohols. Data Selection and "
                             "Recommended Correlations",
                    "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033104",
                    "doi": "10.1063/1.4927858"},
        "sigma": [0.0759, -2.449, 2.47], "exp": [1.134, 3.508, -3.58]}

    _melting = {
        "eq": 1,
        "__doi__": reuck["__doi__"],

        "Tmin": Tt, "Tmax": 650,
        "Tref": Tt, "Pref": 0.187,
        "a0": 1,
        "a2": [5.320770e9, 4.524780e9, 3.888861e10], "exp2": [1, 1.5, 4]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.87414e1, 0.15035e1, -0.28720e1, -0.51345],
        "t": [1., 1.5, 2.5, 5.]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.60230e-1, 0.18855e2, -0.27626e2, 0.11213e2, 0.69039],
        "t": [0.1, 0.65, 0.79, 0.95, 4.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.81104, -0.55661e1, -0.79326e3, 0.19234e4, -0.29219e4,
              0.29660e4, -0.13210e4],
        "t": [0.25, 0.6, 3.5, 4.0, 5.0, 6.0, 7.0]}

    visco0 = {"__name__": "Xiang (2006)",
              "__doi__": {
                  "autor": "Xiang, H.W., Laesecke, A., Huber, M.L.",
                  "title": "A New Reference Correlation for the Viscosity of "
                           "Methanol",
                  "ref": "J. Phys. Chem. Ref. Data 35(4) (2006) 1597",
                  "doi": "10.1063/1.2360605"},

              "eq": 0, "omega": 1,

              "ek": 577.87, "sigma": 0.3408,

              "method": "_visco0"}

    def _Omega(self, T, coef=None):
        """Custom Collision integral calculation procedure"""
        delta = 0.4575
        a = [1.16145, -0.14874, 0.52487, -0.7732, 2.16178, -2.43787, 9.5976e-4,
             0.10225, -0.97346, 0.10657, -0.34528, -0.44557, -2.58055]
        T = T/self._viscosity["ek"]

        OmegaLJ = a[0]*T**a[1] + a[2]*exp(a[3]*T) + a[4]*exp(a[5]*T)     # Eq 4
        OmegaD = a[7]*T**a[8] + a[9]*exp(a[10]*T) + a[11]*exp(a[12]*T)   # Eq 9
        OmegaSM = OmegaLJ*(1+delta**2*OmegaD/(1+a[6]*delta**6))          # Eq 8

        return OmegaSM

    def _visco0(self, rho, T, fase):
        # FIXME: Good values, but tiny desviation

        # Correlation parameters, Table 3
        rhoc = 273.
        sigma0 = 0.3408e-9
        sigmac = 0.7193422e-9

        rhom = rho/self.M*1000
        T_ = T/self._viscosity["ek"]
        Tr = T/self.Tc
        rhor = rho/rhoc

        # Zero-density limit viscosity, Eq 3
        muo = self._Visco0(T)

        # Second viscosity virial coefficient, Eq 13
        bi = [-19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717,
              2491.6597, -787.26086, 14.085455, -0.34664158]
        ti = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.5, 5.5]
        B_ = 0
        for t, b in zip(ti, bi):
            B_ += b/T_**t
        B = B_*Avogadro*sigma0**3

        # Third viscosity virial coefficient, Eq 14
        C = 1.86222085e-3*T_**3*exp(9.990338/T_**0.5)*(Avogadro*sigma0**3)**2

        # Gas phase viscosity, Eq 12
        mug = 1 + B*rhom + C*rhom**2

        # High-density region
        di = [-1.181909, 0.5031030, -0.6268461, 0.5169312, -0.2351349,
              0.053980235, -4.9069617e-3]
        ei = [0, 4.018368, -4.239180, 2.245110, -0.5750698, 2.3021026e-2,
              2.5696775e-2, -6.8372749e-3, 7.2707189e-4, -2.9255711e-5]
        suma = 0
        for i, d in enumerate(di):
            suma += d/Tr**i
        for j, e in enumerate(ei):
            suma += e*rhor**j

        Shs = sigmac*suma                                               # Eq 17
        b = 2*pi*Avogadro*Shs**3/3                        # Close-packed volume
        Xi = b*rhom/4                                        # Packing fraction
        g = (1-0.5*Xi)/(1-Xi)**3                 # Radial distribution function

        # Eq 15
        mue = 1/g + 0.8*b*rhom + 0.761*g*(b*rhom)**2

        # Eq 18
        f = 1/(1+exp(5*(rhor-1)))
        mu = muo*1e-6 * (f*mug + (1-f)*mue)
        return unidades.Viscosity(mu, "Pas")

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Sykioti (2013)",
               "__doi__": {
                   "autor": "Sykioti, E.A., Assael, M.J., Huber, M.L.,"
                            "Perkins, R.A.",
                   "title": "Reference Correlations of the Thermal "
                            "Conductivity of Methanol from the Triple Point "
                            "to 660 K and up to 245 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 42(4) (2013) 043101",
                   "doi": "10.1063/1.4829449"},

               "eq": 1,

               "Toref": Tc, "koref": 1e-3,
               "no_num": [-3.57796, 62.9638, -37.3047, -52.1182, 231.607,
                          44.1575],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [3.33313, -6.08398, 8.18739, -0.261074, 1],
               "to_den": [0, 1, 2, 3, 4],

               "Tref_res": Tc, "rhoref_res": 275.563, "kref_res": 1,
               "nr": [5.56918e-2, 1.04771e-2, 1.12174e-1, -7.45272e-2,
                      -8.43893e-2, 6.37569e-2, 1.97525e-2, -2.46826e-2,
                      -1.5253e-3, 4.34656e-3],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.1487e-9,
               "gam0": 0.05283, "qd": 7.2e-10, "Tcref": 768.9}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_xiang(self):
        """Table 5, Pag 15, saturation state
        Include too  a basic test for Reuck mEoS"""

        # Tiny desviation
        st = Methanol(T=175.63, x=0.5)
        self.assertEqual(round(st.P.MPa, 9), 1.87e-7)
        self.assertEqual(round(st.Gas.rho, 8), 4.10e-6)
        self.assertEqual(round(st.Gas.mu.mPas, 6), 0.005822)
        self.assertEqual(round(st.Liquido.rho, 2), 904.54)
        self.assertEqual(round(st.Liquido.mu.mPas, 2), 12.91)

        st = Methanol(T=200, x=0.5)
        self.assertEqual(round(st.P.MPa, 9), 0.000006096)
        self.assertEqual(round(st.Gas.rho, 8), 0.00011754)
        self.assertEqual(round(st.Gas.mu.mPas, 6), 0.006563)
        self.assertEqual(round(st.Liquido.rho, 2), 880.28)
        self.assertEqual(round(st.Liquido.mu.mPas, 3), 4.544)

        st = Methanol(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.018682)
        self.assertEqual(round(st.Gas.rho, 5), 0.24623)
        self.assertEqual(round(st.Gas.mu.mPas, 6), 0.009678)
        self.assertEqual(round(st.Liquido.rho, 2), 784.51)
        self.assertEqual(round(st.Liquido.mu.mPas, 4), 0.532)

        st = Methanol(T=400, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.77374)
        self.assertEqual(round(st.Gas.rho, 4), 8.7343)
        self.assertEqual(round(st.Gas.mu.mPas, 5), 0.01251)
        self.assertEqual(round(st.Liquido.rho, 2), 678.59)
        self.assertEqual(round(st.Liquido.mu.mPas, 4), 0.1720)

        st = Methanol(T=500, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 6.5250)
        self.assertEqual(round(st.Gas.rho, 2), 109.88)
        self.assertEqual(round(st.Gas.mu.mPas, 5), 0.01891)
        self.assertEqual(round(st.Liquido.rho, 2), 451.53)
        self.assertEqual(round(st.Liquido.mu.mPas, 5), 0.05757)

        # Table 6, Pag 16
        st = Methanol(T=180, P=1e4)
        self.assertEqual(round(st.rho, 2), 900.27)
        self.assertEqual(round(st.mu.mPas, 2), 10.54)

        st = Methanol(T=200, P=1.5e8)
        self.assertEqual(round(st.rho, 2), 944.08)
        self.assertEqual(round(st.mu.mPas, 2), 11.57)

        st = Methanol(T=260, P=1e5)
        self.assertEqual(round(st.rho, 2), 822.11)
        self.assertEqual(round(st.mu.mPas, 3), 1.026)

        st = Methanol(T=300, P=1e6)
        self.assertEqual(round(st.rho, 2), 785.49)
        self.assertEqual(round(st.mu.mPas, 4), 0.5352)

        st = Methanol(T=400, P=1e5)
        self.assertEqual(round(st.rho, 5), 0.97575)
        self.assertEqual(round(st.mu.mPas, 5), 0.01297)

        st = Methanol(T=500, P=1e8)
        self.assertEqual(round(st.rho, 2), 705.57)
        self.assertEqual(round(st.mu.mPas, 4), 0.1548)

        st = Methanol(T=600, P=1e4)
        self.assertEqual(round(st.rho, 6), 0.064244)
        self.assertEqual(round(st.mu.mPas, 5), 0.01981)

    def test_Sykioti(self):
        """Table 6, Pag 9"""
        self.assertEqual(round(Methanol(T=300, rho=850).k.mWmK, 2), 241.48)
        self.assertEqual(round(Methanol(T=400, rho=2).k.mWmK, 3), 25.803)
        self.assertEqual(round(Methanol(T=400, rho=690).k.mWmK, 2), 183.57)
        self.assertEqual(round(Methanol(T=500, rho=10).k.mWmK, 3), 40.495)

    def test_Piazza(self):
        """Table 4, Pag 22"""

        st = Methanol(T=180, rho=6e-6, eq="piazza")
        self.assertEqual(round(st.P.MPa, 10), 0.2802e-6)
        self.assertEqual(round(st.Z, 5), 0.99982)
        # self.assertEqual(round(st.h.kJkg, 5), -151.44350)
        # self.assertEqual(round(st.s.kJkgK, 5), 3.27504)
        # self.assertEqual(round(st.cv.kJkgK, 5), 0.96241)
        # self.assertEqual(round(st.cp.kJkgK, 5), 1.22336)
        # self.assertEqual(round(st.w, 5), 243.62052)

        st = Methanol(T=180, rho=910, eq="piazza")
        self.assertEqual(round(st.P.MPa, 3), 20.184)
        self.assertEqual(round(st.Z, 5), 0.47486)
        # self.assertEqual(round(st.h.kJkg, 5), -1442.54846)
        # self.assertEqual(round(st.s.kJkgK, 5), -4.69592)
        # self.assertEqual(round(st.cv.kJkgK, 5), 1.77392)
        # self.assertEqual(round(st.cp.kJkgK, 5), 2.17984)
        # self.assertEqual(round(st.w, 5), 1624.82189)

        st = Methanol(T=180, rho=1000, eq="piazza")
        self.assertEqual(round(st.P.MPa, 2), 299.59)
        self.assertEqual(round(st.Z, 5), 6.41420)
        # self.assertEqual(round(st.h.kJkg, 5), -1192.49438)
        # self.assertEqual(round(st.s.kJkgK, 5), -4.92516)
        # self.assertEqual(round(st.cv.kJkgK, 5), 1.81615)
        # self.assertEqual(round(st.cp.kJkgK, 5), 2.11901)
        # self.assertEqual(round(st.w, 5), 2200.22161)

        st = Methanol(T=180, rho=1040, eq="piazza")
        self.assertEqual(round(st.P.MPa, 2), 487.06)
        self.assertEqual(round(st.Z, 5), 10.02671)
        # self.assertEqual(round(st.h.kJkg, 5), -1028.09811)
        # self.assertEqual(round(st.s.kJkgK, 5), -5.03225)
        # self.assertEqual(round(st.cv.kJkgK, 5), 1.82635)
        # self.assertEqual(round(st.cp.kJkgK, 5), 2.09866)
        # self.assertEqual(round(st.w, 5), 2455.69417)

        st = Methanol(T=450, rho=10, eq="piazza")
        self.assertEqual(round(st.P.MPa, 4), 1.0641)
        self.assertEqual(round(st.Z, 5), 0.91132)
        # self.assertEqual(round(st.h.kJkg, 5), 184.85996)
        # self.assertEqual(round(st.s.kJkgK, 5), -0.068334)
        # self.assertEqual(round(st.cv.kJkgK, 5), 1.99022)
        # self.assertEqual(round(st.cp.kJkgK, 5), 2.53739)
        # self.assertEqual(round(st.w, 5), 348.93454)

        st = Methanol(T=450, rho=25, eq="piazza")
        self.assertEqual(round(st.P.MPa, 4), 2.2321)
        self.assertEqual(round(st.Z, 5), 0.76463)
        # self.assertEqual(round(st.h.kJkg, 5), 84.80645)
        # self.assertEqual(round(st.s.kJkgK, 5), -0.45469)
        # self.assertEqual(round(st.cv.kJkgK, 5), 3.71647)
        # self.assertEqual(round(st.cp.kJkgK, 5), 5.95074)
        # self.assertEqual(round(st.w, 5), 309.79673)

        st = Methanol(T=450, rho=610, eq="piazza")
        self.assertEqual(round(st.P.MPa, 4), 4.5142)
        self.assertEqual(round(st.Z, 5), 0.06337)
        # self.assertEqual(round(st.h.kJkg, 5), -697.62906)
        # self.assertEqual(round(st.s.kJkgK, 5), -2.22522)
        # self.assertEqual(round(st.cv.kJkgK, 5), 3.00081)
        # self.assertEqual(round(st.cp.kJkgK, 5), 4.23642)
        # self.assertEqual(round(st.w, 5), 622.75228)

        st = Methanol(T=450, rho=900, eq="piazza")
        self.assertEqual(round(st.P.MPa, 2), 456.15)
        self.assertEqual(round(st.Z, 5), 4.34045)
        # self.assertEqual(round(st.h.kJkg, 5), -364.62474)
        # self.assertEqual(round(st.s.kJkgK, 5), -2.73998)
        # self.assertEqual(round(st.cv.kJkgK, 5), 2.66908)
        # self.assertEqual(round(st.cp.kJkgK, 5), 3.10135)
        # self.assertEqual(round(st.w, 5), 2047.99620)
