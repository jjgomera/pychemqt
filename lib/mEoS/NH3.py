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


class NH3(MEoS):
    """Multiparameter equation of state for ammonia"""
    name = "ammonia"
    CASNumber = "7664-41-7"
    formula = "NH3"
    synonym = "R-717"
    rhoc = unidades.Density(225.)
    Tc = unidades.Temperature(405.40)
    Pc = unidades.Pressure(11333.0, "kPa")
    M = 17.03026  # g/mol
    Tt = unidades.Temperature(195.495)
    Tb = unidades.Temperature(239.823)
    f_acent = 0.25601
    momentoDipolar = unidades.DipoleMoment(1.470, "Debye")
    id = 63

    Fi1 = {"R": 8.314471,
           "ao_log": [1, -1],
           "pow": [0, 1, 1/3, -1.5, -1.75],
           "ao_pow": [-15.81502, 4.255726, 11.47434, -1.296211, 0.5706757],
           "ao_exp": [], "titao": [],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 5.111814,
           "an": [-0.42966650e2, -0.10243792e-1, 0.38750775e-4, -0.46406097e-7,
                  0.20268561e-10],
           "pow": [-1.001, 1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Baehr and "
                    "Tillner-Roth (1993)",
        "__doi__": {"autor": "Baehr, H.D., Tillner-Roth, R.",
                    "title": "Thermodynamic Properties of Environmentally "
                             "Acceptable Refrigerants: Equations of State and "
                             "Tables for Ammonia, R22, R134a, R152a, and R123",
                    "ref": "Springer-Verlag, Berlin, 1994.",
                    "doi": "10.1007/978-3-642-79400-1"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 700., "Pmax": 1000000.0, "rhomax": 52.915,
        "Pmin": 6.09, "rhomin": 43.035,

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
        "__doi__": {"autor": "Ahrendts, J. and Baehr, H.D.",
                    "title": "The Thermodynamic Properties of Ammonia",
                    "ref": "VDI-Forsch., Number 596, 1979.",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": 195.486, "Tmax": 600., "Pmax": 400000.0, "rhomax": 44.0,
        "Pmin": 6.0339, "rhomin": 43.137,

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
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600., "Pmax": 100000.0, "rhomax": 52.915,
        "Pmin": 6.0531, "rhomin": 43.158,
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
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.3143,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

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

    eq = tillner, ahrendts, shortSpan, sun

    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "Tmin": Tt, "Tmax": 700.0,
                "a1": [], "exp1": [], "a2": [], "exp2": [],
                "a3": [0.2533125e4], "exp3": [1]}
    _surface = {"sigma": [0.1028, -0.09453], "exp": [1.211, 5.585]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70993e1, -0.24330e1, 0.87591e1, -0.64091e1, -0.21185e1],
        "exp": [1., 1.5, 1.7, 1.95, 4.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.34488e2, -0.12849e3, 0.17382e3, -0.10699e3, 0.30339e2],
        "exp": [0.58, 0.75, 0.9, 1.1, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.38435, -4.0846, -6.6634, -31.881, 213.06, -246.48],
        "exp": [0.218, 0.55, 1.5, 3.7, 5.5, 5.8]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [4.99318220, -0.61122364, 0.0, 0.18535124, -0.11160946],
              "__name__": "Fenghour (1995)",
              "__doi__": {"autor": "Fenghour, A., Wakeham, W.A., Vesovic, V., Watson, J.T.R., Millat, J., and Vogel, E.",
                          "title": "The viscosity of ammonia",
                          "ref": "J. Phys. Chem. Ref. Data 24, 1649 (1995)",
                          "doi": "10.1063/1.555961"},
              "__test__":
                    # Appendix II, pag 1664
                    """
                    >>> st=NH3(T=200, P=1e5)
                    >>> print("%0.2f" % st.mu.muPas)
                    507.47
                    >>> st=NH3(T=290, P=1e6)
                    >>> print("%0.2f" % st.mu.muPas)
                    142.93
                    >>> st=NH3(T=250, P=1e7)
                    >>> print("%0.2f" % st.mu.muPas)
                    233.81
                    >>> st=NH3(T=300, P=1e5)
                    >>> print("%0.2f" % st.mu.muPas)
                    10.16
                    >>> st=NH3(T=350, P=1.8e7)
                    >>> print("%0.2f" % st.mu.muPas)
                    91.36
                    >>> st=NH3(T=400, P=5e7)
                    >>> print("%0.2f" % st.mu.muPas)
                    77.29
                    >>> st=NH3(T=490, P=1e6)
                    >>> print("%0.2f" % st.mu.muPas)
                    17.49
                    >>> st=NH3(T=550, P=1e5)
                    >>> print("%0.2f" % st.mu.muPas)
                    19.79
                    >>> st=NH3(T=680, P=5e7)
                    >>> print("%0.2f" % st.mu.muPas)
                    31.90
                    """
                    # Appendix III, pag 1667
                    """
                    >>> st=NH3(T=196, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    196 0.0063 0.0039 6.85 43.0041 553.31
                    >>> st=NH3(T=240, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    240 0.1022 0.0527 8.06 40.0318 254.85
                    >>> st=NH3(T=280, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    280 0.5509 0.2573 9.27 36.9389 158.12
                    >>> st=NH3(T=300, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    300 1.0617 0.4845 9.89 35.2298 129.33
                    >>> st=NH3(T=340, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    340 3.0803 1.4325 11.33 31.2641 88.55
                    >>> st=NH3(T=360, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    360 4.7929 2.3598 12.35 28.7879 65.49
                    >>> st=NH3(T=380, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    380 7.1403 3.9558 14.02 25.6059 58.31
                    >>> st=NH3(T=398, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    398 9.9436 7.0447 17.67 21.0667 43.95
                    >>> st=NH3(T=402, x=0.5)
                    >>> print "%0.0f %0.4f %0.4f %0.2f %0.4f %0.2f " % ( \
                        st.T, st.P.MPa, st.Gas.rhoM, st.Gas.mu.muPas, st.Liquido.rhoM, st.Liquido.mu.muPas)
                    402 10.6777 8.5479 19.69 19.0642 39.20
                    """
                    ,

              "ek": 386., "sigma": 0.2957,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 8.8135503/M**0.5,

              "n_virial": [-0.17999496e1, 0.46692621e2, -0.53460794e3,
                           0.33604074e4, -0.13019164e5, 0.33414230e5,
                           -0.58711743e5, 0.71426686e5, -0.59834012e5,
                           0.33652741e5, -0.1202735e5, 0.24348205e4, -0.20807957e3],
              "t_virial": [0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5,
                           -5.5, -6],
              "Tref_virial": 386., "etaref_virial": 0.015570557,

              "Tref_res": 386., "rhoref_res": 1.*M, "etaref_res": 1,
              "n_poly": [2.19664285e-1, -0.83651107e-1, 0.17366936e-2,
                         -0.64250359e-2, 1.67668649e-4, -1.49710093e-4, 0.77012274e-4],
              "t_poly": [-2, -4, -0, -1, -2, -3, -4],
              "d_poly": [2, 2, 3, 3, 4, 4, 4],
              "g_poly": [0, 0, 0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 0, 0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1, "critical": "NH3",
               "__name__": "Tufeu (1984)",
               "__doi__": {"autor": "Tufeu, R., Ivanov, D.Y., Garrabos, Y., and Le Neindre, B.",
                            "title": "Thermal conductivity of ammonia in a large temperature and pressure range including the critical region",
                            "ref": "Ber. Bunsenges. Phys. Chem., 88:422-427, 1984",
                            "doi": "10.1002/bbpc.19840880421"},

               "Tref": 1., "kref": 1.,
               "no": [0.3589e-1, -0.1750e-3, 0.4551e-6, 0.1685e-9, -0.4828e-12],
               "co": [0, 1, 2, 3, 4],

               "Trefb": 1., "rhorefb": 0.05871901, "krefb": 1.,
               "nb": [0.16207e-3, 0.12038e-5, -0.23139e-8, 0.32749e-11],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0]}

    _thermal = thermo0,


class Test(TestCase):

    def test_tillner(self):
        # Selected point from pag 42, saturation state
        st = NH3(T=-77.65+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00609)
        self.assertEqual(round(st.Liquido.rho, 2), 732.90)
        self.assertEqual(round(st.Gas.rho, 4), 0.0641)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -143.13)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1484.4)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1341.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -0.4715)
        self.assertEqual(round(st.Svap.kJkgK, 4), 7.5927)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 7.1212)

        st = NH3(T=-50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.04084)
        self.assertEqual(round(st.Liquido.rho, 2), 702.09)
        self.assertEqual(round(st.Gas.rho, 4), 0.3806)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -24.73)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1415.9)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1391.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.0945)
        self.assertEqual(round(st.Svap.kJkgK, 4), 6.3451)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.4396)

        st = NH3(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.42938)
        self.assertEqual(round(st.Liquido.rho, 2), 638.57)
        self.assertEqual(round(st.Gas.rho, 4), 3.4567)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1262.2)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1462.2)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Svap.kJkgK, 4), 4.6210)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.6210)

        st = NH3(T=50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 2.03403)
        self.assertEqual(round(st.Liquido.rho, 2), 562.86)
        self.assertEqual(round(st.Gas.rho, 3), 15.785)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 440.62)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1050.5)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1491.1)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.7990)
        self.assertEqual(round(st.Svap.kJkgK, 4), 3.2507)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.0497)

        st = NH3(T=100+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 6.25527)
        self.assertEqual(round(st.Liquido.rho, 2), 456.63)
        self.assertEqual(round(st.Gas.rho, 3), 56.117)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 721.00)
        self.assertEqual(round(st.Hvap.kJkg, 2), 715.63)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1436.6)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 2.5797)
        self.assertEqual(round(st.Svap.kJkgK, 4), 1.9178)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 4.4975)

        st = NH3(T=130+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 10.89768)
        self.assertEqual(round(st.Liquido.rho, 2), 312.29)
        self.assertEqual(round(st.Gas.rho, 2), 156.77)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 992.02)
        self.assertEqual(round(st.Hvap.kJkg, 2), 247.30)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1239.3)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.2437)
        self.assertEqual(round(st.Svap.kJkgK, 4), 0.6134)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.8571)

        st = NH3(T=132+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 11.28976)
        self.assertEqual(round(st.Liquido.rho, 2), 262.70)
        self.assertEqual(round(st.Gas.rho, 2), 193.88)
        self.assertEqual(round(st.Liquido.h.kJkg, 1), 1063.0)
        self.assertEqual(round(st.Hvap.kJkg, 1), 108.60)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1171.6)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.4160)
        self.assertEqual(round(st.Svap.kJkgK, 4), 0.2680)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 3.6840)

        st = NH3(P=1e4, x=0.5)
        self.assertEqual(round(st.T.C, 2), -71.22)
        self.assertEqual(round(st.Liquido.rho, 2), 726.04)
        self.assertEqual(round(st.Gas.rho, 4), 0.1020)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -115.99)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1469.3)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1353.3)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -0.3349)
        self.assertEqual(round(st.Svap.kJkgK, 4), 7.2762)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.9412)

        st = NH3(P=1e5, x=0.5)
        self.assertEqual(round(st.T.C, 2), -33.59)
        self.assertEqual(round(st.Liquido.rho, 2), 682.29)
        self.assertEqual(round(st.Gas.rho, 4), 0.8787)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 47.60)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1370.3)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1417.9)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.4068)
        self.assertEqual(round(st.Svap.kJkgK, 4), 5.7199)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.1267)

        st = NH3(P=1e6, x=0.5)
        self.assertEqual(round(st.T.C, 2), 24.90)
        self.assertEqual(round(st.Liquido.rho, 2), 602.92)
        self.assertEqual(round(st.Gas.rho, 4), 7.7823)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 317.16)
        self.assertEqual(round(st.Hvap.kJkg, 1), 1166.2)
        self.assertEqual(round(st.Gas.h.kJkg, 1), 1483.4)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.4072)
        self.assertEqual(round(st.Svap.kJkgK, 4), 3.9128)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.3200)

        st = NH3(P=1e7, x=0.5)
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
        st = NH3(T=-75+273.15, P=1e4)
        self.assertEqual(round(st.rho, 2), 730.10)
        self.assertEqual(round(st.h.kJkg, 2), -131.97)
        self.assertEqual(round(st.s.kJkgK, 4), -0.4148)

        st = NH3(T=273.15, P=2e4)
        self.assertEqual(round(st.rho, 4), 0.1504)
        self.assertEqual(round(st.h.kJkg, 1), 1499.2)
        self.assertEqual(round(st.s.kJkgK, 4), 7.2228)

        st = NH3(T=125+273.15, P=3e4)
        self.assertEqual(round(st.rho, 4), 0.1545)
        self.assertEqual(round(st.h.kJkg, 1), 1769.3)
        self.assertEqual(round(st.s.kJkgK, 4), 7.8372)

        st = NH3(T=-50+273.15, P=5e4)
        self.assertEqual(round(st.rho, 2), 702.09)
        self.assertEqual(round(st.h.kJkg, 2), -24.72)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0945)

        st = NH3(T=25+273.15, P=1e5)
        self.assertEqual(round(st.rho, 4), 0.6942)
        self.assertEqual(round(st.h.kJkg, 1), 1546.7)
        self.assertEqual(round(st.s.kJkgK, 4), 6.6083)

        st = NH3(T=175+273.15, P=2e5)
        self.assertEqual(round(st.rho, 4), 0.9183)
        self.assertEqual(round(st.h.kJkg, 1), 1882.1)
        self.assertEqual(round(st.s.kJkgK, 4), 7.1800)

        st = NH3(T=-10+273.15, P=3e5)
        self.assertEqual(round(st.rho, 2), 652.06)
        self.assertEqual(round(st.h.kJkg, 2), 154.02)
        self.assertEqual(round(st.s.kJkgK, 4), 0.8292)

        st = NH3(T=273.15, P=5e5)
        self.assertEqual(round(st.rho, 2), 638.62)
        self.assertEqual(round(st.h.kJkg, 2), 200.05)
        self.assertEqual(round(st.s.kJkgK, 4), 0.9998)

        st = NH3(T=25+273.15, P=1e6)
        self.assertEqual(round(st.rho, 4), 7.7778)
        self.assertEqual(round(st.h.kJkg, 1), 1483.7)
        self.assertEqual(round(st.s.kJkgK, 4), 5.3211)

        st = NH3(T=50+273.15, P=2e6)
        self.assertEqual(round(st.rho, 3), 15.449)
        self.assertEqual(round(st.h.kJkg, 1), 1493.5)
        self.assertEqual(round(st.s.kJkgK, 4), 5.0641)

        st = NH3(T=-50+273.15, P=3e6)
        self.assertEqual(round(st.rho, 2), 703.33)
        self.assertEqual(round(st.h.kJkg, 2), -22.08)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0875)

    def test_shortSpan(self):
        # Table III, Pag 117
        st = NH3(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 2.4758)
        self.assertEqual(round(st.P.MPa, 3), 136.271)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.1915)

        st2 = NH3(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 776.68)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 2.07031)
