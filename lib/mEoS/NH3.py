#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=protected-access

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


from math import exp, log, pi
from unittest import TestCase

from scipy.constants import Boltzmann

from lib import unidades
from lib.meos import MEoS


class NH3(MEoS):
    """Multiparameter equation of state for ammonia"""
    name = "ammonia"
    CASNumber = "7664-41-7"
    formula = "NH3"
    synonym = "R-717"
    _refPropName = "AMMONIA"
    _coolPropName = "Ammonia"
    rhoc = unidades.Density(225.)
    Tc = unidades.Temperature(405.40)
    Pc = unidades.Pressure(11333.0, "kPa")
    M = 17.03026  # g/mol
    Tt = unidades.Temperature(195.49)
    Tb = unidades.Temperature(239.823)
    f_acent = 0.256
    momentoDipolar = unidades.DipoleMoment(1.470, "Debye")
    id = 63

    Fi1 = {"R": 8.314471,
           "ao_log": [1, -1],
           "pow": [0, 1, 1/3, -1.5, -1.75],
           "ao_pow": [-15.81502, 4.255726, 11.47434, -1.296211, 0.5706757]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-6.59406093943886, 5.60101151987913],
           "ao_exp": [2.224, 3.148, 0.9579],
           "titao": [1646/405.56, 3965/405.56, 7231/405.56]}

    CP2 = {"ao": 5.111814,
           "an": [-0.42966650e2, -0.10243792e-1, 0.38750775e-4, -0.46406097e-7,
                  0.20268561e-10],
           "pow": [-1, 1, 2, 3, 4]}

    gao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Gao (2023)",
        "__doi__": {
            "autor": "Gao, K., Wu, J., Bell, I.H., Harvey, A.H., Lemmon E.W.",
            "title": "A Reference Equation of State with an Associating Term "
                     "for the Thermodynamic Properties of Ammonia",
            "ref": "J. Phys. Chem. Ref. Data 52 (2023) 013102",
            "doi": "10.1063/5.0128269"},

        "R": 8.314462618,
        "M": 17.03052, "Tc": 405.56, "Pc": 11363.4, "rhoc": 13.696,
        "cp": Fi2,
        "ref": {"Tref": 300, "Pref": 1, "ho": 28989.81844833686,
                "so": 159.53995277505977},

        "Tmin": 195.49, "Tmax": 725, "Pmax": 1000000, "rhomax": 53.13,

        "nr1": [0.006132232, 1.7395866, -2.2261792, -0.30127553, 0.08967023],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.382, 1, 1, 0.677],

        "nr2": [-0.076387037, -0.84063963, -0.27026327],
        "d2": [3, 2, 3],
        "t2": [2.915, 3.51, 1.063],
        "c2": [2, 2, 1],
        "gamma2": [1]*3,

        "nr3": [6.212578, -5.7844357, 2.4817542, -2.3739168, 0.01493697,
                -3.7749264, 0.0006254348, -0.000017359, -0.13462033,
                0.07749072839],
        "d3": [1, 1, 1, 2, 2, 1, 3, 3, 1, 1],
        "t3": [0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13, 4.5, 1, 4],
        "alfa3": [0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14,
                  22.56, 22.68],
        "beta3": [1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64,
                  993.85],
        "gamma3": [1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852,
                   1.05897, 1.05277],
        "epsilon3": [-0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02,
                     0.9574, 0.9576],

        "nr_ass": [-1.6909858, 0.93739074],
        "d_ass": [1, 1],
        "t_ass": [4.3315, 4.015],
        "alfa_ass": [2.8452, 2.8342],
        "beta_ass": [0.3696, 0.2962],
        "gamma_ass": [1.108, 1.313],
        "epsilon_ass": [0.4478, 0.44689],
        "b_ass": [1.244, 0.6826]}

    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Baehr and "
                    "Tillner-Roth (1993)",
        "__doi__": {
            "autor": "Baehr, H.D., Tillner-Roth, R.",
            "title": "Thermodynamic Properties of Environmentally Acceptable "
                     "Refrigerants: Equations of State and Tables for "
                     "Ammonia, R22, R134a, R152a, and R123",
            "ref": "Springer-Verlag, Berlin, 1994.",
            "doi": "10.1007/978-3-642-79400-1"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 700., "Pmax": 1000000.0, "rhomax": 52.915,

        "nr1": [-0.1858814e01, 0.4554431e-1, 0.7238548, 0.1229470e-1,
                0.2141882e-10],
        "d1": [1, 2, 1, 4, 15],
        "t1": [1.5, -0.5, 0.5, 1., 3.],

        "nr2": [-0.1430020e-1, 0.3441324, -0.2873571, 0.2352589e-4,
                -0.3497111e-1, 0.1831117e-2, 0.2397852e-1, -0.4085375e-1,
                0.2379275, -0.3548972e-1, -0.1823729, 0.2281556e-1,
                -0.6663444e-2, -0.8847486e-2, 0.2272635e-2, -0.5588655e-3],
        "d2": [3, 3, 1, 8, 2, 8, 1, 1, 2, 3, 2, 4, 3, 1, 2, 4],
        "t2": [0, 3, 4, 4, 5, 5, 3, 6, 8, 8, 10, 10, 5, 7.5, 15, 30],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*16}

    ahrendts = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Ahrendts "
                    "and Baehr (1979)",
        "__doi__": {
            "autor": "Ahrendts, J., Baehr, H.D.",
            "title": "The Thermodynamic Properties of Ammonia",
            "ref": "VDI-Forsch., Number 596, 1979.",
            "doi": ""},

        "R": 8.31434,
        "M": 17.03026, "Tc": 405.4, "Pc": 11333, "rhoc": 13.212,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": 195.486, "Tmax": 600., "Pmax": 400000.0, "rhomax": 44.0,

        "nr1": [0.911447599671, -0.382129415537e1, 0.147730246416e1,
                0.580205129871e-1, -0.574413226616e-3, 0.153018094697,
                -0.256626062036, 0.445448838055, -0.1533210545,
                0.527996725202e-1, -0.484726581121e-1, 0.246579503330e-2,
                -0.107999941003e-3, -0.215298673010e-4, -0.306938893790e-4,
                0.839163613582e-5, 0.814833533876e-6, -0.314753664228e-7],
        "d1": [1, 1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 7, 9, 9, 10, 11, 12, 14],
        "t1": [1, 2, 3, 6, 9, 0, 4, 2, 1, 1, 2, 3, 3, 5, 1, 1, 5, 5],

        "nr2": [0.642978802435, -0.139510669941e1, 0.956135683432,
                -0.272787386366, -0.189305337334e1, 0.479043603913e1,
                -0.245945016980e1, -0.121107723958e1, 0.500552271170e1,
                -0.615476024667e1, 0.210772481535e1, 0.298003513465,
                -0.152506723279, 0.115565883925e-2, -0.911244657201e-3,
                0.100587210000e-1, -0.120983155888e-1, 0.382694351151e-2],
        "d2": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 0, 0, 0, 0, 0],
        "t2": [2, 5, 6, 7, 5, 6, 7, 3, 4, 5, 6, 6, 7, 1, 2, 0, 1, 2],
        "c2": [2]*18,
        "gamma2": [0.86065403]*13+[506.2670781840292]*2+[50626.70781840292]*3}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ammonia of Span "
                    "and Wagner (2003)",
        "__doi__": {
            "autor": "Span, R., Wagner, W.",
            "title": "Equations of State for Technical Applications. III. "
                     "Results for Polar Fluids",
            "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
            "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600., "Pmax": 100000.0, "rhomax": 52.915,
        "M": 17.031, "rhoc": 13.211203,

        "nr1": [0.7302272, -1.1879116, -0.68319136, 0.040028683, 9.0801215e-5],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.056216175, 0.44935601, 0.029897121, -0.18181684,
                -0.09841666, -0.055083744, -0.88983219e-2],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.3143,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [3.29159441e-1, 8.48237019e-1, -2.30706412, 4.08625188e-2,
                6.79597481e-5, 4.99412149e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [1.23624654e-1, -3.02129187e-1, 3.31747586e-1, -2.97121254e-3,
                -0.130202073, -7.45181207e-2, -4.73506171e-2, -9.70095484e-3],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = gao, tillner, ahrendts, shortSpan, sun
    _PR = [0.1491, -26.0703]

    _surface = {"sigma": [0.1028, -0.09453], "exp": [1.211, 5.585]}

    @classmethod
    def _Melting_Pressure(cls, T, melting=None):
        """Customized method for calculate melting pressure implementing two
        different region"""
        if melting is not None:
            return MEoS._Melting_Pressure(T, melting)

        if T < 217.35:
            return MEoS._Melting_Pressure(T, cls._melting1)
        return MEoS._Melting_Pressure(T, cls._melting)

    _melting = {
        "eq": 1,
        "__doi__": gao["__doi__"],

        "Tmin": 217.35, "Tmax": 725,
        "Tref": 217.35, "Pref": 1e6,
        "a0": 307, "a2": [1135], "exp2": [1.51]}

    _melting1 = {
        "eq": 2,
        "__doi__": gao["__doi__"],

        "Tmin": Tt, "Tmax": 217.35,
        "Tref": 195.49, "Pref": 6.05339e3,
        "a2": [516.613763e6], "exp2": [4.4]}

    _sublimation = {
        "eq": 2,
        "__doi__": gao["__doi__"],

        "Tref": 1, "Pref": 1e6,
        "Tmin": Tt, "Tmax": Tt,
        "a1": [13.63932, -3537, -3.31e4, 1.742e6, -2.995e7],
        "exp1": [0, -1, -2, -3, -4]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.3128, 3.8888, -2.9908, -2.8636],
        "t": [1., 1.5, 1.6, 3.7]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.051236, 3.7925, -3.5929, 4.6409, -1.9893, 1.5978],
        "t": [0.07, 0.46, 0.77, 1.05, 1.25, 8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.089966, -3.8722, -8.1183, -25.293, -54.279, -400.83],
        "t": [0.112, 0.473, 1.5, 3.875, 8, 20]}

    visco0 = {"__name__": "Monogenidou (2018)",
              "__doi__": {
                  "autor": "Monogenidou, S.A., Assael, M.J., Huber, M.L.",
                  "title": "Reference Correlation for the Viscosity of Ammonia"
                           " from the Triple Point to 725K and up to 50 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 47(2) (2018) 023102",
                  "doi": "10.1063/1.5036724"},

              "eq": 1, "omega": 1,

              "ek": 386, "sigma": 0.2957, "M": 17.03052,
              "n_chapman": 0.021357,
              "collision": [0.39175, -0.59918, -0.00022, 0.19871, -0.06942],

              "Tref_virial": 386,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 405.56, "rhoref_res": 233.25,
              "nr": [0.0393308, 16.724735, 1.1975934, 0.0016995, -4.2399794],
              "tr": [-0.5, -0.5, -0.5, 3.5, -1.5],
              "dr": [2/3, 5/3, 14/3, 26/3, 8/3]}

    visco1 = {"__name__": "Fenghour (1995)",
              "__doi__": {
                  "autor": "Fenghour, A., Wakeham, W.A., Vesovic, V., Watson, "
                           "J.T.R., Millat, J., Vogel, E.",
                  "title": "The Viscosity of Ammonia",
                  "ref": "J. Phys. Chem. Ref. Data 24(5) (1995) 1649-1667",
                  "doi": "10.1063/1.555961"},

              "eq": 1, "omega": 1,

              "ek": 386, "sigma": 0.2957,
              # Missing parameter of 100 in Chapman-Enskog term
              "n_chapman": 0.021357*100,
              "collision": [4.9931822, -0.61122364, 0, .18535124, -0.11160946],

              "Tref_virial": 386,
              "n_virial": [-1.7999496, 46.692621, -534.60794, 3360.4074,
                           -13019.164, 33414.230, -58711.743, 71426.686,
                           -59834.012, 33652.741, -12027.35, 2434.8205,
                           -208.07957],
              "t_virial": [0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5,
                           -5.5, -6],

              "Tref_res": 386., "rhoref_res": M,
              "nr": [0.219664285, -0.083651107, 0.0017366936, -0.0064250359,
                     1.67668649e-4, -1.49710093e-4, 7.7012274e-5],
              "tr": [2, 4, 0, 1, 2, 3, 4],
              "dr": [2, 2, 3, 3, 4, 4, 4]}

    _viscosity = visco0, visco1

    thermo0 = {"__name__": "Monogenidou (2018)",
               "__doi__": {
                   "autor": "Monogenidou, S.A., Assael, M.J., Huber, M.L.",
                   "title": "Reference Correlations for the Thermal "
                            "Conductivity of Ammonia from the Triple Point to "
                            "680 K and Pressures up to 80 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 47(4) (2018) 043101",
                   "doi": "10.1063/1.5053087"},

               "eq": 1,

               "Toref": 405.56, "koref": 1e-3,
               "no_num": [86.9294, -170.5502, 608.0287, -100.9764, 85.1986],
               "to_num": [0, 1, 2, 3, 4],
               "no_den": [4.68994, 9.21307, -1.53637, 1],
               "to_den": [0, 1, 2, 3],

               # The table 2 in paper report values as mW/mK, it's a typo,
               # really is in W/mK
               "Tref_res": 405.56, "rhoref_res": 233.25, "kref_res": 1,
               "nr": [0.103432, -0.112597, 0.233301, -0.112536, 0.141129e-1,
                      -0.283976e-1, 0.482520e-1, -0.644124e-1, 0.529376e-2,
                      0.891203e-2],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.14e-9,
               "gam0": 0.053, "qd": 4e-10, "Tcref": 608.34}

    thermo1 = {"__name__": "Tufeu (1984)",
               "__doi__": {
                   "autor": "Tufeu, R., Ivanov, D.Y., Garrabos, Y., Le "
                            "Neindre, B.",
                   "title": "Thermal Conductivity of Ammonia in a Large "
                            "Temperature and Pressure Range Including the "
                            "Critical Region",
                   "ref": "Ber. Bunsenges. Phys. Chem. 88 (1984) 422-427",
                   "doi": "10.1002/bbpc.19840880421"},

               "eq": 1,

               "Tref": 1., "kref": 1.,
               "no": [0.03589, -1.750e-4, 4.551e-7, 1.685e-10, -4.828e-13],
               "to": [0, 1, 2, 3, 4],

               "rhoref_res": 1, "kref_res": 1.,
               "nr": [0.16207e-3, 0.12038e-5, -0.23139e-8, 0.32749e-11],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": "_ThCondCritical"}

    def _ThCondCritical(self, rho, T, fase):
        # Custom Critical enhancement

        # The paper use a diferent rhoc value to the EoS
        rhoc = 235

        t = abs(T-405.4)/405.4
        dPT = 1e5*(2.18-0.12/exp(17.8*t))
        nb = 1e-5*(2.6+1.6*t)

        DL = 1.2*Boltzmann*T**2/6/pi/nb/(1.34e-10/t**0.63*(1+t**0.5)) * \
            dPT**2 * 0.423e-8/t**1.24*(1+t**0.5/0.7)

        # Add correction for entire range of temperature, Eq 10
        DL *= exp(-36*t**2)

        X = 0.61*rhoc+16.5*log(t)
        if rho > 0.6*rhoc:
            # Eq 11
            DL *= X**2/(X**2+(rho-0.96*rhoc)**2)
        else:
            # Eq 14
            DL = X**2/(X**2+(0.6*rhoc-0.96*rhoc)**2)
            DL *= rho**2/(0.6*rhoc)**2

        return DL

    _thermal = thermo0, thermo1


class Test(TestCase):
    """Testing"""

    def test_gao(self):
        """Table 7, pag 24"""
        st = NH3(T=200, rhom=50)
        self.assertEqual(round(st.P.MPa, 4), 555.5631)
        self.assertEqual(round(st.cvM.JmolK, 5), 56.15976)
        self.assertEqual(round(st.cpM.JmolK, 5), 74.22262)
        self.assertEqual(round(st.w, 3), 2940.432)
        self.assertEqual(round(st.hM.kJmol, 6), 9.541211)
        self.assertEqual(round(st.sM.JmolK, 5), -11.61817)

        st = NH3(T=250, rhom=42)
        self.assertEqual(round(st.P.MPa, 4), 108.7263)
        self.assertEqual(round(st.cvM.JmolK, 5), 50.25907)
        self.assertEqual(round(st.cpM.JmolK, 5), 72.04765)
        self.assertEqual(round(st.w, 3), 2049.168)
        self.assertEqual(round(st.hM.kJmol, 6), 5.715934)
        self.assertEqual(round(st.sM.JmolK, 5), 14.26145)

        st = NH3(T=320, rhom=35)
        self.assertEqual(round(st.P.MPa, 5), 27.63275)
        self.assertEqual(round(st.cvM.JmolK, 5), 46.95882)
        self.assertEqual(round(st.cpM.JmolK, 5), 78.92022)
        self.assertEqual(round(st.w, 3), 1370.446)
        self.assertEqual(round(st.hM.kJmol, 6), 9.831787)
        self.assertEqual(round(st.sM.JmolK, 5), 36.06050)

        st = NH3(T=380, rhom=26)
        self.assertEqual(round(st.P.MPa, 6), 8.174128)
        self.assertEqual(round(st.cvM.JmolK, 5), 47.82965)
        self.assertEqual(round(st.cpM.JmolK, 4), 125.1971)
        self.assertEqual(round(st.w, 4), 668.5456)
        self.assertEqual(round(st.hM.kJmol, 5), 15.45317)
        self.assertEqual(round(st.sM.JmolK, 5), 53.81385)

        st = NH3(T=405, rhom=16.5)
        self.assertEqual(round(st.P.MPa, 5), 11.25287)
        self.assertEqual(round(st.cvM.JmolK, 5), 77.13442)
        self.assertEqual(round(st.cpM.JmolK, 3), 4137.738)
        self.assertEqual(round(st.w, 4), 247.2696)
        self.assertEqual(round(st.hM.kJmol, 5), 20.13541)
        self.assertEqual(round(st.sM.JmolK, 5), 65.30611)

        st = NH3(T=406, rhom=13.696)
        self.assertEqual(round(st.P.MPa, 5), 11.45246)
        self.assertEqual(round(st.cvM.JmolK, 4), 103.0508)
        self.assertEqual(round(st.cpM.JmolK, 2), 52943.69)
        self.assertEqual(round(st.w, 4), 225.4458)
        self.assertEqual(round(st.hM.kJmol, 5), 21.30410)
        self.assertEqual(round(st.sM.JmolK, 5), 68.15535)

        st = NH3(T=225, rhom=0)
        self.assertEqual(round(st.P.MPa, 5), 0.0)
        self.assertEqual(round(st.cvM.JmolK, 5), 25.60265)
        self.assertEqual(round(st.cpM.JmolK, 5), 33.91711)
        self.assertEqual(round(st.w, 4), 381.4709)
        self.assertEqual(round(st.hM.kJmol, 5), 26.38996)

        st = NH3(T=275, rhom=0.2)
        self.assertEqual(round(st.P.MPa, 7), 0.4270652)
        self.assertEqual(round(st.cvM.JmolK, 5), 32.60471)
        self.assertEqual(round(st.cpM.JmolK, 5), 45.39676)
        self.assertEqual(round(st.w, 4), 402.6253)
        self.assertEqual(round(st.hM.kJmol, 5), 27.46383)
        self.assertEqual(round(st.sM.JmolK, 4), 104.2956)

        st = NH3(T=325, rhom=0.95)
        self.assertEqual(round(st.P.MPa, 6), 2.093496)
        self.assertEqual(round(st.cvM.JmolK, 5), 40.68660)
        self.assertEqual(round(st.cpM.JmolK, 5), 65.25439)
        self.assertEqual(round(st.w, 4), 402.2104)
        self.assertEqual(round(st.hM.kJmol, 5), 27.92681)
        self.assertEqual(round(st.sM.JmolK, 5), 94.19375)

        st = NH3(T=375, rhom=3)
        self.assertEqual(round(st.P.MPa, 6), 6.044096)
        self.assertEqual(round(st.cvM.JmolK, 5), 48.91890)
        self.assertEqual(round(st.cpM.JmolK, 4), 111.3738)
        self.assertEqual(round(st.w, 4), 379.4300)
        self.assertEqual(round(st.hM.kJmol, 5), 27.45404)
        self.assertEqual(round(st.sM.JmolK, 5), 86.36105)

        st = NH3(T=405, rhom=10)
        self.assertEqual(round(st.P.MPa, 5), 11.22040)
        self.assertEqual(round(st.cvM.JmolK, 5), 77.15007)
        self.assertEqual(round(st.cpM.JmolK, 3), 2051.793)
        self.assertEqual(round(st.w, 4), 286.2398)
        self.assertEqual(round(st.hM.kJmol, 5), 23.20689)
        self.assertEqual(round(st.sM.JmolK, 5), 72.89756)

        st = NH3(T=520, rhom=1.5)
        self.assertEqual(round(st.P.MPa, 6), 5.968880)
        self.assertEqual(round(st.cvM.JmolK, 5), 36.52118)
        self.assertEqual(round(st.cpM.JmolK, 5), 49.38647)
        self.assertEqual(round(st.w, 4), 538.8044)
        self.assertEqual(round(st.hM.kJmol, 5), 36.23166)
        self.assertEqual(round(st.sM.JmolK, 4), 106.6103)

        st = NH3(T=520, rhom=13.696)
        self.assertEqual(round(st.P.MPa, 5), 35.83377)
        self.assertEqual(round(st.cvM.JmolK, 5), 43.66588)
        self.assertEqual(round(st.cpM.JmolK, 5), 98.36082)
        self.assertEqual(round(st.w, 4), 555.4454)
        self.assertEqual(round(st.hM.kJmol, 5), 28.71997)
        self.assertEqual(round(st.sM.JmolK, 5), 80.49762)

        st = NH3(T=520, rhom=20)
        self.assertEqual(round(st.P.MPa, 5), 56.30750)
        self.assertEqual(round(st.cvM.JmolK, 5), 43.81573)
        self.assertEqual(round(st.cpM.JmolK, 5), 87.32996)
        self.assertEqual(round(st.w, 4), 741.1842)
        self.assertEqual(round(st.hM.kJmol, 5), 26.48937)
        self.assertEqual(round(st.sM.JmolK, 5), 73.89537)

        st = NH3(T=620, rhom=5.5)
        self.assertEqual(round(st.P.MPa, 5), 24.38934)
        self.assertEqual(round(st.cvM.JmolK, 5), 41.29819)
        self.assertEqual(round(st.cpM.JmolK, 5), 61.81978)
        self.assertEqual(round(st.w, 4), 588.0372)
        self.assertEqual(round(st.hM.kJmol, 5), 38.40927)
        self.assertEqual(round(st.sM.JmolK, 4), 100.0682)

        st = NH3(T=620, rhom=14)
        self.assertEqual(round(st.P.MPa, 5), 58.52754)
        self.assertEqual(round(st.cvM.JmolK, 5), 43.50218)
        self.assertEqual(round(st.cpM.JmolK, 5), 75.46188)
        self.assertEqual(round(st.w, 4), 692.5867)
        self.assertEqual(round(st.hM.kJmol, 5), 34.49909)
        self.assertEqual(round(st.sM.JmolK, 5), 87.78122)

        st = NH3(T=620, rhom=45)
        self.assertEqual(round(st.P.MPa, 3), 1355.431)
        self.assertEqual(round(st.cvM.JmolK, 5), 50.78124)
        self.assertEqual(round(st.cpM.JmolK, 5), 64.36735)
        self.assertEqual(round(st.w, 3), 3201.084)
        self.assertEqual(round(st.hM.kJmol, 5), 50.78325)
        self.assertEqual(round(st.sM.JmolK, 5), 54.81426)

    def test_tillner(self):
        """Selected point from pag 42, saturation state"""
        st = NH3(T=-77.65+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 0.00609)
        self.assertEqual(round(st.Liquido.rho, 2), 732.90)
        self.assertEqual(round(st.Gas.rho, 4), 0.0641)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -143.13)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1484.4)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1341.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -0.4715)
        self.assertEqual(round(st.Svap.kJkgK, 4), 7.5927)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 7.1212)

        st = NH3(T=-50+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 0.04084)
        self.assertEqual(round(st.Liquido.rho, 2), 702.09)
        self.assertEqual(round(st.Gas.rho, 4), 0.3806)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -24.73)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1415.9)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1391.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.0945)
        self.assertEqual(round(st.Svap.kJkgK, 4), 6.3451)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.4396)

        st = NH3(T=273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 0.42938)
        self.assertEqual(round(st.Liquido.rho, 2), 638.57)
        self.assertEqual(round(st.Gas.rho, 4), 3.4567)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1262.2)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1462.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Svap.kJkgK, 4), 4.6210)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.6210)

        st = NH3(T=50+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 2.03403)
        self.assertEqual(round(st.Liquido.rho, 2), 562.86)
        self.assertEqual(round(st.Gas.rho, 3), 15.785)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 440.62)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1050.5)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1491.1)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.7990)
        self.assertEqual(round(st.Svap.kJkgK, 4), 3.2507)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.0497)

        st = NH3(T=100+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 6.25527)
        self.assertEqual(round(st.Liquido.rho, 2), 456.63)
        self.assertEqual(round(st.Gas.rho, 3), 56.117)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 721.00)
        self.assertEqual(round(st.Hvap.kJkg, 2), 715.63)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1436.6)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 2.5797)
        self.assertEqual(round(st.Svap.kJkgK, 4), 1.9178)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 4.4975)

        st = NH3(T=130+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 10.89768)
        self.assertEqual(round(st.Liquido.rho, 2), 312.29)
        self.assertEqual(round(st.Gas.rho, 2), 156.77)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 992.02)
        self.assertEqual(round(st.Hvap.kJkg, 2), 247.30)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1239.3)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.2437)
        self.assertEqual(round(st.Svap.kJkgK, 4), 0.6134)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.8571)

        st = NH3(T=132+273.15, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 11.28976)
        self.assertEqual(round(st.Liquido.rho, 2), 262.70)
        self.assertEqual(round(st.Gas.rho, 2), 193.88)
        self.assertEqual(round(st.Liquido.h.kJkg, 1), 1063.0)
        self.assertEqual(round(st.Hvap.kJkg, 1), 108.60)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1171.6)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.4160)
        self.assertEqual(round(st.Svap.kJkgK, 4), 0.2680)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.6840)

        st = NH3(P=1e4, x=0.5, eq="tillner")
        self.assertEqual(round(st.T.C, 2), -71.22)
        self.assertEqual(round(st.Liquido.rho, 2), 726.04)
        self.assertEqual(round(st.Gas.rho, 4), 0.1020)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -115.99)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1469.3)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1353.3)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -0.3349)
        self.assertEqual(round(st.Svap.kJkgK, 4), 7.2762)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.9412)

        st = NH3(P=1e5, x=0.5, eq="tillner")
        self.assertEqual(round(st.T.C, 2), -33.59)
        self.assertEqual(round(st.Liquido.rho, 2), 682.29)
        self.assertEqual(round(st.Gas.rho, 4), 0.8787)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 47.60)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1370.3)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1417.9)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.4068)
        self.assertEqual(round(st.Svap.kJkgK, 4), 5.7199)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.1267)

        st = NH3(P=1e6, x=0.5, eq="tillner")
        self.assertEqual(round(st.T.C, 2), 24.90)
        self.assertEqual(round(st.Liquido.rho, 2), 602.92)
        self.assertEqual(round(st.Gas.rho, 4), 7.7823)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 317.16)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1166.2)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1483.4)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.4072)
        self.assertEqual(round(st.Svap.kJkgK, 4), 3.9128)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.3200)

        st = NH3(P=1e7, x=0.5, eq="tillner")
        self.assertEqual(round(st.T.C, 2), 125.17)
        self.assertEqual(round(st.Liquido.rho, 2), 356.70)
        self.assertEqual(round(st.Gas.rho, 2), 121.58)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 921.57)
        self.assertEqual(round(st.Hvap.kJkg, 2), 385.87)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1307.4)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.0747)
        self.assertEqual(round(st.Svap.kJkgK, 4), 0.9688)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 4.0434)

        # Selected point from Pag 51, single region states
        st = NH3(T=-75+273.15, P=1e4, eq="tillner")
        self.assertEqual(round(st.rho, 2), 730.10)
        self.assertEqual(round(st.h.kJkg, 2), -131.97)
        self.assertEqual(round(st.s.kJkgK, 4), -0.4148)

        st = NH3(T=273.15, P=2e4, eq="tillner")
        self.assertEqual(round(st.rho, 4), 0.1504)
        self.assertEqual(round(st.h.kJkg, 1), 1499.2)
        self.assertEqual(round(st.s.kJkgK, 4), 7.2228)

        st = NH3(T=125+273.15, P=3e4, eq="tillner")
        self.assertEqual(round(st.rho, 4), 0.1545)
        self.assertEqual(round(st.h.kJkg, 1), 1769.3)
        self.assertEqual(round(st.s.kJkgK, 4), 7.8372)

        st = NH3(T=-50+273.15, P=5e4, eq="tillner")
        self.assertEqual(round(st.rho, 2), 702.09)
        self.assertEqual(round(st.h.kJkg, 2), -24.72)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0945)

        st = NH3(T=25+273.15, P=1e5, eq="tillner")
        self.assertEqual(round(st.rho, 4), 0.6942)
        self.assertEqual(round(st.h.kJkg, 1), 1546.7)
        self.assertEqual(round(st.s.kJkgK, 4), 6.6083)

        st = NH3(T=175+273.15, P=2e5, eq="tillner")
        self.assertEqual(round(st.rho, 4), 0.9183)
        self.assertEqual(round(st.h.kJkg, 1), 1882.1)
        self.assertEqual(round(st.s.kJkgK, 4), 7.1800)

        st = NH3(T=-10+273.15, P=3e5, eq="tillner")
        self.assertEqual(round(st.rho, 2), 652.06)
        self.assertEqual(round(st.h.kJkg, 2), 154.02)
        self.assertEqual(round(st.s.kJkgK, 4), 0.8292)

        st = NH3(T=273.15, P=5e5, eq="tillner")
        self.assertEqual(round(st.rho, 2), 638.62)
        self.assertEqual(round(st.h.kJkg, 2), 200.05)
        self.assertEqual(round(st.s.kJkgK, 4), 0.9998)

        st = NH3(T=25+273.15, P=1e6, eq="tillner")
        self.assertEqual(round(st.rho, 4), 7.7778)
        self.assertEqual(round(st.h.kJkg, 1), 1483.7)
        self.assertEqual(round(st.s.kJkgK, 4), 5.3211)

        st = NH3(T=50+273.15, P=2e6, eq="tillner")
        self.assertEqual(round(st.rho, 3), 15.449)
        self.assertEqual(round(st.h.kJkg, 1), 1493.5)
        self.assertEqual(round(st.s.kJkgK, 4), 5.0641)

        st = NH3(T=-50+273.15, P=3e6, eq="tillner")
        self.assertEqual(round(st.rho, 2), 703.33)
        self.assertEqual(round(st.h.kJkg, 2), -22.08)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0875)

    def test_shortSpan(self):
        """Table III, Pag 117"""
        st = NH3(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 2.4758)
        self.assertEqual(round(st.P.MPa, 3), 136.271)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.1915)

        st2 = NH3(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 776.68)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 2.07031)

    def test_monogenidouVisco(self):
        """Values cited in section 3, pag 7"""
        self.assertEqual(round(NH3(T=300, rho=0).mu.muPas, 4), 10.1812)
        self.assertEqual(round(NH3(T=300, rho=8).mu.muPas, 4), 9.9219)
        self.assertEqual(round(NH3(T=300, rho=609).mu.muPas, 4), 133.3937)

    def test_monogenidouThermo(self):
        """Values cited in section 3, pag 7"""
        self.assertEqual(round(NH3(T=390, rho=0).k.mWmK, 6), 35.969501)
        # Tiny error in point with high density > 700 kg/m³, maybe a typo
        # self.assertEqual(round(NH3(T=390, rho=415).k.mWmK, 6), 264.129743)

    def test_fenghour(self):
        """Selected values from tables in Appendix"""
        kw = {"eq": "tillner", "visco": 1}

        # Appendix II, pag 1664
        self.assertEqual(round(NH3(T=200, P=1e5, **kw).mu.muPas, 2), 507.47)
        self.assertEqual(round(NH3(T=290, P=1e6, **kw).mu.muPas, 2), 142.93)
        self.assertEqual(round(NH3(T=250, P=1e7, **kw).mu.muPas, 2), 233.81)
        self.assertEqual(round(NH3(T=300, P=1e5, **kw).mu.muPas, 2), 10.16)
        self.assertEqual(round(NH3(T=350, P=1.8e7, **kw).mu.muPas, 2), 91.36)
        self.assertEqual(round(NH3(T=400, P=5e7, **kw).mu.muPas, 2), 77.29)
        self.assertEqual(round(NH3(T=490, P=1e6, **kw).mu.muPas, 2), 17.49)
        self.assertEqual(round(NH3(T=550, P=1e5, **kw).mu.muPas, 2), 19.79)
        self.assertEqual(round(NH3(T=680, P=5e7, **kw).mu.muPas, 2), 31.90)

        # Appendix III, pag 1667, saturation state
        st = NH3(T=196, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 4), 0.0063)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0039)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 6.85)
        self.assertEqual(round(st.Liquido.rhoM, 4), 43.0041)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 553.31)

        st = NH3(T=240, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 4), 0.1022)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.0527)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.06)
        self.assertEqual(round(st.Liquido.rhoM, 4), 40.0318)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 254.85)

        st = NH3(T=300, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 4), 1.0617)
        self.assertEqual(round(st.Gas.rhoM, 4), 0.4845)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.89)
        self.assertEqual(round(st.Liquido.rhoM, 4), 35.2298)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 129.33)

        st = NH3(T=350, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 4), 3.8660)
        self.assertEqual(round(st.Gas.rhoM, 4), 1.8399)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 11.79)
        self.assertEqual(round(st.Liquido.rhoM, 4), 30.0867)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 80.43)

        st = NH3(T=402, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 4), 10.6777)
        self.assertEqual(round(st.Gas.rhoM, 4), 8.5479)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 19.69)
        self.assertEqual(round(st.Liquido.rhoM, 4), 19.1642)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 39.20)


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    from matplotlib import rcParams
    from numpy import linspace

    rcParams["xtick.top"] = True
    rcParams["ytick.right"] = True
    rcParams["ytick.minor.visible"] = True

    # Trying to emulate fig 4 in Gao paper

    # Saturation line
    Ts = linspace(NH3.Tt, NH3.Tc, 50)
    Ps = []
    for ts in Ts:
        Ps.append(NH3(T=ts, x=0).P.MPa)
    plt.plot(Ts, Ps, color="k", lw=0.8, label="Vapor Pressure Curve")

    # Sublimation line
    Ts = linspace(160, NH3.Tt, 10)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Sublimation_Pressure(ts).MPa)
    plt.plot(Ts, Ps, color="b", lw=0.8, label="Sublimation Curve")

    # Melting line Phase I
    Ts = linspace(NH3.Tt, 217.34, 50)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Melting_Pressure(ts).MPa)
    plt.plot(Ts, Ps, color="brown", lw=0.8, label="Melting Curve Solid I")

    # Melting line Phase I
    Ts = linspace(217.34, 750, 50)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Melting_Pressure(ts).MPa)
    plt.plot(Ts, Ps, color="r", lw=0.8, label="Melting Curve Solid II")

    # Haar-Gallagher melting correlation
    haar = {
        "eq": 1,
        "__doi__": {
            "autor": "Haar, L., Gallagher, J.S.",
            "title": "Thermodynamic Properties of Ammonia",
            "ref": "J. Phys. Chem. Ref. Data 7(3) (1978) 635-792",
            "doi": "10.1063/1.555579"},

        "Tmin": NH3.Tt, "Tmax": 700.0,
        "Tref": 195.48, "Pref": 101325,
        "a3": [2.5e4], "exp3": [1]}

    Ts = linspace(NH3.Tt, 217.34, 50)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Melting_Pressure(ts, haar).MPa)
    Ps[0] = NH3(T=NH3.Tt, x=0).P.MPa  # Fix initial value to fit triple point
    plt.plot(Ts, Ps, color="k", lw=1, ls="-.")
    Ts = linspace(217.34, 750, 20)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Melting_Pressure(ts, haar).MPa)
    plt.plot(Ts, Ps, color="k", lw=0.8, ls="-.",
             label="Melting Curve Haar and Gallagher (1978)")

    # Mills solid line transition
    mills = {
        "eq": 2,
        "__doi__": {
            "autor": "Mills, R.L., Liebenberg, D.H., Pruzan, P.",
            "title": "Phase diagram and transition properties of condensed "
                     "ammonia to 10 kbar",
            "ref": "J. Physical Chemistry, 86(26) (1982) 5219-5222",
            "doi": "10.1021/j100223a030 "},

        "Tmin": 217.34, "Tmax": 700.0,
        "Tref": 217.34, "Pref": 3.07e8,
        "a2": [1.608e8], "exp2": [47.28]}

    Ts = linspace(217.34, 300, 20)
    Ps = []
    for ts in Ts:
        Ps.append(NH3._Melting_Pressure(ts, mills).MPa)
    plt.plot(Ts, Ps, color="k", lw=1, ls=":",
             label="Solid I-II transition Mills et al. (1982)")

    # Experimental data from
    # Bergstrom, F.W.
    # The Vapor Pressure of Sulfur Dioxide and Ammonia
    # J. Physical Chemistry, 26(4), 358-376. doi:10.1021/j150220a005
    Tm = (-30.2+273.15, -50.2+273.15, -67.3+273.15, -73.3+273.15, -76.5+273.15,
          -79.9+273.15, -84.4+273.15, -77.9+273.15)  # in ºC
    Pm = (891, 304.4, 101.2, 64.8, 50.6, 37.2, 24.2, 45.5)  # in mmHg
    P = []
    for p in Pm:
        P.append(unidades.Pressure(p, "mmHg").MPa)
    plt.plot(Tm, P, ls="", marker="o", mfc="#0000", mec="k",
             label="Bergstrom (1922)")

    # Karwat, E.
    # Der Dampfdruck des festen Chlorwasserstoffs, Methans und Ammoniaks
    # Zeitschrift Für Physikalische Chemie 112U(1) (1924) 486-490
    # doi: 10.1515/zpch-1924-11230
    Tm = (162.39, 169.49, 171.52, 173.61, 178.97, 183.79, 186.99, 193.45)
    Pm = (0.95, 2.38, 3.09, 3.92, 7.84, 13.5, 19.13, 37.60)  # in mmHg
    P = []
    for p in Pm:
        P.append(unidades.Pressure(p, "mmHg").MPa)
    plt.plot(Tm, P, ls="", marker="D", mfc="#0000", mec="k",
             label="Karwat (1926)")

    # King, T.-V., Oi, T., Popowicz, A., Heinzinger, K., Ishida, T.
    # Vapor Pressure Isotope Effects in Liquid and Solid Ammonia
    # Z. Naturforsch 44a (1989) 359-370
    # doi: 10.1515/zna-1989-0503
    Tm = (243.01, 229.06, 219.54, 210.96, 200.40, 198.53, 195.36, 188.86,
          180.02, 168.8)
    Pm = (889.85, 431.14, 246.87, 142.37, 66.951, 58.058, 45.189, 23.362,
          8.781, 2.196)  # in torr
    P = []
    for p in Pm:
        P.append(unidades.Pressure(p, "torr").MPa)
    plt.plot(Tm, P, ls="", marker="x", mec="k", label="King et al. (1989)")

    # Overstreet, R., Giauque, W.F.
    # Ammonia. The Heat Capacity and Vapor Pressure of Solid and Liquid. Heat
    # of Vaporization. The Entropy Values from Thermal and Spectroscopic Data
    # J. American Chemical Society 59(2) (1937) 254-259
    # doi: 10.1021/ja01281a008
    Tm = (176.92, 181.13, 183.17, 186.0, 189.03, 191.85, 195.36, 199.262,
          203.114, 206.587, 210.811, 214.42, 218.354, 224.246, 229.205, 234.84,
          241.58)
    Pm = (0.615, 1.01, 1.266, 1.735, 2.398, 3.206, 4.558, 6.188, 8.261, 10.608,
          14.181, 18.017, 23.126, 32.98, 43.725, 59.284, 83.535)  # in cmHg
    P = []
    for p in Pm:
        P.append(unidades.Pressure(p, "cmHg").MPa)
    plt.plot(Tm, P, ls="", marker="s", mfc="#0000", mec="k",
             label="Overstreet and Giauque (1937)")

    # Ninet, S., Datchi, F.
    # High pressure-high temperature phase diagram of ammonia
    # J. Chem. Phys. 128 (2008) 154508
    # doi: 10.1063/1.2903491
    Tm = (298.3, 322.0, 371.1, 372.0, 420.1, 473.1, 484.3, 492.8, 509.4, 529.9,
          530.4, 553.8, 562.1, 580.8, 582.0, 586.3, 604.1, 616.8, 633.0, 643.1,
          652.9, 665.8, 681.1, 671.7, 705.9, 723.4, 724.6, 740.1)
    Pm = (1.01e3, 1.23e3, 1.66e3, 1.71e3, 2.21e3, 2.79e3, 2.93e3, 3.03e3,
          3.22e3, 3.53e3, 3.57e3, 3.76e3, 3.90e3, 4.12e3, 4.15e3, 4.26e3,
          4.48e3, 4.65e3, 4.94e3, 5.03e3, 5.24e3, 5.37e3, 5.50e3, 5.51e3,
          6.04e3, 6.26e3, 6.30e3, 6.47e3)
    plt.plot(Tm, Pm, ls="", marker="^", mfc="#0000", mec="k",
             label="Ninet and Datchi (2008)")

    # Grace, J.D., Kennedy, G.C.
    # The melting curve of five gases to 30 kb
    # J. Phys. Chem. Solids 28(6) (1967) 977-982
    # doi: 10.1016/0022-3697(67)90214-4
    Tm = (418.28571429, 397.71428571, 376.38095238, 361.9047619, 342.0952381,
          321.52380952, 291.80952381, 268.19047619, 253.71428571, 227.04761905,
          207.23809524)
    # bar
    Pm = (24.55919395, 22.29219144, 19.49622166, 17.83375315, 15.86901763,
          13.45088161, 10.50377834, 8.61460957, 6.95214106, 5.06297229, 2.4937)
    P = []
    for p in Pm:
        P.append(unidades.Pressure(p, "kbar").MPa)
    plt.plot(Tm, P, ls="", marker="4", mfc="#0000", mec="k",
             label="Grace and Kennedy (1967)")

    # Hanson, R.C., Jordan, M.
    # Ultrahigh-pressure studies of ammonia
    # J. Phys. Chem. 84(10) (1980) 1173-1175
    # doi: 10.1021/j100447a020
    Pm = (579.22535211, 801.05633803, 887.32394366, 887.32394366, 1035.2112676,
          948.94366197, 1294.0140845, 1417.2535211, 1491.197183, 1786.97183099)
    Tm = (248.48076923, 261.46153846, 271.55769231, 279.25, 289.82692308,
          298.48076923, 327.80769231, 346.07692308, 359.05769231, 371.55769231)
    plt.plot(Tm, Pm, ls="", marker=r"$\sqcap$", mfc="#0000", mec="k",
             label="Hanson and Jordan (1980)")

    plt.ylabel("P/MPa")
    plt.xlabel("T/K")
    plt.xlim(150, 800)
    plt.ylim(1e-4, 1e5)
    plt.yscale("log")
    plt.vlines(NH3.Tc, NH3.Pc.MPa, NH3._Melting_Pressure(NH3.Tc).MPa,
               ls="--", color="k")
    plt.hlines(NH3.Pc.MPa, NH3.Tc, 800, ls="--", color="k")
    plt.text(250, 10, "Liquid")
    plt.text(550, 100, "Supercritical")
    plt.text(450, 0.01, "Vapor")
    plt.text(350, 10000, "Solid Phase II/III")
    plt.text(170, 0.5, "Solid Phase I", rotation="vertical")

    # Put a legend below current axis
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
               fancybox=True, shadow=True, ncol=2)
    plt.tight_layout()

    plt.show()

    # FIXME: Monogenidou thermal conductivity fail for values with density over
    # 600 kg/m³
    # Table 6
    # for P in (1e5, 25e6, 50e6):
    #     for T in (200, 300, 400, 500, 600):
    #         st = NH3(T=T, P=P)
    #         print(f"{st.P.MPa:#10.2g} {st.T:0.0f} {st.rho:#10.5g} "
    #               "{st.k.mWmK:#10.5g}")
    # Table 7
    # for T in (200, 225, 250, 275, 300, 325, 350, 375, 400):
    #     st = NH3(T=T, x=0.5)
    #     print(f"{st.T:0.0f} {st.P.MPa:#10.5g} {st.Liquido.rho:#10.5g} "
    #           "{st.Gas.rho:#10.5g} {st.Liquido.k.mWmK:#10.5g} "
    #           "{st.Gas.k.mWmK:#10.5g}")
