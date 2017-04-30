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

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, 1./3, -1.5, -1.75],
           "ao_pow": [-15.81502, 4.255726, 11.47434, -1.296211, 0.5706757],
           "ao_exp": [], "titao": [],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 5.111814,
           "an": [-0.42966650e2, -0.10243792e-1, 0.38750775e-4, -0.46406097e-7,
                  0.20268561e-10],
           "pow": [-1.001, 1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Baehr and Tillner-Roth (1993)",
        "__doi__": {"autor": "Baehr, H.D. and Tillner-Roth, R.",
                    "title": "Thermodynamic Properties of Environmentally Acceptable Refrigerants; Equations of State and Tables for Ammonia, R22, R134a, R152a, and R123",
                    "ref": "Springer-Verlag, Berlin, 1994.",
                    "doi": ""},
        "__test__":
            # Table, Pag 42
            """
            >>> st=NH3(T=-77.65+273.15, x=0.5)
            >>> print "%0.2f %0.5f %0.5g %0.3g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            -77.65 0.00609 732.9 0.0641 -143.14 1484.4 1341.2 -0.4715 7.5928 7.1213
            >>> st=NH3(T=-50+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.4g %0.4g %0.5g %0.5g %0.3g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            -50 0.04084 702.09 0.3806 -24.73 1415.9 1391.2 0.0945 6.3451 6.4396
            >>> st=NH3(T=-25+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.4g %0.5g %0.5g %0.4g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            -25 0.15147 671.53 1.2959 86.01 1344.6 1430.7 0.5641 5.4187 5.9827
            >>> st=NH3(T=273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            0 0.42938 638.57 3.4567 200 1262.2 1462.2 1 4.621 5.621
            >>> st=NH3(T=25+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            25 1.00324 602.76 7.8069 317.67 1165.8 1483.4 1.4089 3.91 5.3188
            >>> st=NH3(T=50+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            50 2.03403 562.86 15.785 440.62 1050.5 1491.1 1.799 3.2507 5.0497
            >>> st=NH3(T=75+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            75 3.71045 516.23 29.923 572.38 907.35 1479.7 2.1823 2.6062 4.7885
            >>> st=NH3(T=100+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            100 6.25527 456.63 56.117 721 715.63 1436.6 2.5797 1.9178 4.4975
            >>> st=NH3(T=125+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            125 9.97022 357.8 120.73 919.68 389.44 1309.1 3.0702 0.9781 4.0483
            >>> st=NH3(T=132+273.15, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Hvap.kJkg, st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Svap.kJkgK, st.Gas.s.kJkgK)
            132 11.28976 262.7 193.88 1063 108.59 1171.6 3.416 0.268 3.684
            """
            # Table, Pag 51
            """
            >>> st=NH3(T=-75+273.15, P=1e4)
            >>> print "%0.0f %0.2f %0.5g %0.4g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -75 730.10 -131.97 -0.4148
            >>> st=NH3(T=-30+273.15, P=2e4)
            >>> print "%0.0f %0.4f %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -30 0.1693 1436.9 6.9812
            >>> st=NH3(T=273.15, P=3e4)
            >>> print "%0.0f %0.4f %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 0.2260 1498.3 7.0226
            >>> st=NH3(T=175+273.15, P=5e4)
            >>> print "%0.0f %0.4f %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            175 0.2288 1884.8 7.8612
            >>> st=NH3(T=-50+273.15, P=1e5)
            >>> print "%0.0f %0.2f %0.4g %0.3g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -50 702.11 -24.67 0.0944
            >>> st=NH3(T=-30+273.15, P=1e5)
            >>> print "%0.0f %0.4f %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -30 0.8641 1426.1 6.1607
            >>> st=NH3(T=273.15, P=2e5)
            >>> print "%0.0f %0.4f %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 1.5468 1483.8 6.0557
            >>> st=NH3(T=125+273.15, P=3e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            125 1.5597 1762.6 6.7009
            >>> st=NH3(T=-50+273.15, P=5e5)
            >>> print "%0.0f %0.5g %0.4g %0.3g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -50 702.28 -24.32 0.0934
            >>> st=NH3(T=273.15, P=1e6)
            >>> print "%0.0f %0.5g %0.5g %0.4g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 638.97 200.37 0.9981
            >>> st=NH3(T=273.15, P=2e6)
            >>> print "%0.0f %0.5g %0.5g %0.4g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 639.67 201.01 0.9947
            >>> st=NH3(T=-50+273.15, P=3e6)
            >>> print "%0.0f %0.5g %0.4g %0.3g" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -50 703.33 -22.08 0.0875
            """,

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

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Ahrendts and Baehr (1979)",
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

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ammonia of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"},
        "__test__": """
            >>> st=NH3(T=700, rho=200, eq=2)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            2.4759 136.271 4.1916
            >>> st2=NH3(T=750, rho=100, eq=2)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            776.69 2.07036
            """, # Table III, Pag 117

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600., "Pmax": 100000.0, "rhomax": 52.915,
        "Pmin": 6.0531, "rhomin": 43.158,

        "nr1": [0.7302272, -0.11879116e1, -0.68319136, 0.40028683e-1, 0.90801215e-4],
        "d1": [1, 1, 1, 3, 7],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.56216175e-1, 0.44935601, 0.29897121e-1, -0.18181684,
                -0.9841666e-1, -0.55083744e-1, -0.88983219e-2],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
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
                -1.30202073e-1, -7.45181207e-2, -4.73506171e-2, -9.70095484e-3],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2, helmholtz3, helmholtz4

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
        "ao": [-0.38435, -0.40846e1, -0.66634e1, -0.31881e2, 0.21306e3, -0.24648e3],
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
