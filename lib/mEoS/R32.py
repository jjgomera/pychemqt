#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


class R32(MEoS):
    """Multiparameter equation of state for R32"""
    name = "difluoromethane"
    CASNumber = "75-10-5"
    formula = "CH2F2"
    synonym = "R32"
    rhoc = unidades.Density(424.)
    Tc = unidades.Temperature(351.255)
    Pc = unidades.Pressure(5782., "kPa")
    M = 52.024  # g/mol
    Tt = unidades.Temperature(136.34)
    Tb = unidades.Temperature(221.499)
    f_acent = 0.2769
    momentoDipolar = unidades.DipoleMoment(1.978, "Debye")
    id = 645

    Fi1 = {"ao_log": [1, 3.004486],
           "pow": [0, 1],
           "ao_pow": [-8.258096, 6.353098],
           "ao_exp": [1.160761, 2.645151, 5.794987, 1.129475],
           "titao": [798/Tc, 4185/Tc, 1806/Tc, 11510/Tc]}

    Fi2 = {"ao_log": [1, 2.999660],
           "pow": [0, 1],
           "ao_pow": [-8.253834, 6.351918],
           "ao_exp": [3.12115, 0.9994221, 2.412721, 3.055435],
           "titao": [4.559777, 2.164788, 1.234687e1, 5.877902]}

    CP2 = {"ao": 36.79959/8.314471,
           "an": [-0.06304821/8.314471, 3.757936e-4/8.314471, -3.219812e-7/8.314471],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}


    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Tillner-Roth & Yokozeki (1997)",
        "__doi__": {"autor": "Tillner-Roth, R., Yokozeki, A.",
                    "title": "An international standard equation of state for difluoromethane (R-32) for temperatures from the triple point at 136.4 K to 435 K at pressures up to 70 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 26 (1997), 1273 – 1328.",
                    "doi": "10.1063/1.556002"},

        "__test__":
            #Table 12, Pag 1293
            """
            >>> st=R32(T=136.34, x=0.5)
            >>> print "%0.2f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.3f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -136.81 0.05 1429.3 0.0022 -19.07 463.38 444.31 -0.104 3.2937 1.593 0.660
            >>> st=R32(T=-120+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -120 0.48 1388.4 0.0195 7.52 447.81 455.33 0.079 3.003 1.573 0.674
            >>> st=R32(T=-100+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -100 3.81 1339 0.1385 38.83 429.48 468.31 0.2711 2.7515 1.560 0.703
            >>> st=R32(T=-80+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -80 18.65 1288.4 0.6129 90.02 410.71 480.73 0.4415 2.5679 1.561 0.754
            >>> st=R32(T=-60+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -60 64.96 1235.7 1.9690 101.38 390.73 492.11 0.5958 2.4290 1.576 0.833
            >>> st=R32(T=-40+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.4f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -40 177.41 1180.2 5.0651 133.23 368.79 502.02 0.7382 2.3200 1.608 0.940
            >>> st=R32(T=-20+273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.3f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            -20 405.75 1120.6 11.157 165.94 344.03 509.97 0.8720 2.2310 1.661 1.075
            >>> st=R32(T=273.15, x=0.5)
            >>> print "%0.0f %0.2f %0.5g %0.3f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            0 813.10 1055.3 22.091 200.00 315.30 515.30 1.0000 2.1543 1.745 1.251
            >>> st=R32(T=20+273.15, x=0.5)
            >>> print "%0.0f %0.1f %0.5g %0.3f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            20 1474.5 981.38 40.856 236.12 280.78 516.90 1.1253 2.0831 1.886 1.514
            >>> st=R32(T=40+273.15, x=0.5)
            >>> print "%0.0f %0.1f %0.5g %0.3f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            40 2478.3 893.04 73.268 275.61 237.09 512.71 1.2520 2.0091 2.163 2.001
            >>> st=R32(T=60+273.15, x=0.5)
            >>> print "%0.0f %0.1f %0.5g %0.2f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            60 3933.2 773.31 135.21 321.93 175.51 497.44 1.3898 1.9166 3.001 3.441
            >>> st=R32(T=74+273.15, x=0.5)
            >>> print "%0.0f %0.1f %0.5g %0.2f %0.2f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f" % (\
                st.T.C, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Hvap.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK)
            74 5304.5 624.57 240.12 367.53 98.88 466.41 1.5179 1.8028 8.051 12.094
            """
            # Table 14, Pag 1297
            """
            >>> st=R32(T=-85+273.15, P=1e4)
            >>> print "%0.0f %0.4f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -85 0.3357 478.09 2.6526
            >>> st=R32(T=-80+273.15, P=2e4)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -80 1288.4 70.02 0.4415
            >>> st=R32(T=273.15, P=3e4)
            >>> print "%0.0f %0.4f %0.1f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 0.6907 541.3 2.7544
            >>> st=R32(T=100+273.15, P=5e4)
            >>> print "%0.0f %0.4f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            100 0.6721 628.13 2.9782
            >>> st=R32(T=-55+273.15, P=1e5)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -55 1222.2 109.29 0.6324
            >>> st=R32(T=160+273.15, P=2e5)
            >>> print "%0.0f %0.4f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            160 2.9073 686.47 2.8669
            >>> st=R32(T=-30+273.15, P=3e5)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -30 1151.0 149.46 0.8059
            >>> st=R32(T=273.15, P=5e5)
            >>> print "%0.0f %0.3f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 12.594 526.88 2.2651
            >>> st=R32(T=-40+273.15, P=1e6)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -40 1181.9 133.53 0.7365
            >>> st=R32(T=30+273.15, P=2e6)
            >>> print "%0.0f %0.2f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            30 940.17 255.28 1.1877
            >>> st=R32(T=150+273.15, P=3e6)
            >>> print "%0.0f %0.3f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            150 49.577 650.52 2.3649
            >>> st=R32(T=-85+273.15, P=5e6)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            -85 1307.2 64.65 0.3932
            >>> st=R32(T=100+273.15, P=1e7)
            >>> print "%0.0f %0.2f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            100 566.69 416.83 1.6330
            >>> st=R32(T=273.15, P=2e7)
            >>> print "%0.0f %0.1f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            0 1112.9 204.05 0.9501
            >>> st=R32(T=160+273.15, P=4e7)
            >>> print "%0.0f %0.2f %0.2f %0.4f" % (st.T.C, st.rho, st.h.kJkg, st.s.kJkgK)
            160 760.12 475.86 1.6700
            """
            # Table 15, Pag 1319
            """
            >>> st=R32(T=-20+273.15, P=5e5)
            >>> print "%0.3f" % st.cp.kJkgK
            1.660
            >>> st=R32(T=155+273.15, P=2e4)
            >>> print "%0.3f" % st.cp.kJkgK
            1.033
            >>> st=R32(T=5+273.15, P=1e6)
            >>> print "%0.3f" % st.cp.kJkgK
            1.773
            >>> st=R32(T=80+273.15, P=6e6)
            >>> print "%0.3f" % st.cp.kJkgK
            61.802
            >>> st=R32(T=65+273.15, P=1e7)
            >>> print "%0.3f" % st.cp.kJkgK
            2.140
            """
            # Table 16, Pag 1322
            """
            >>> st=R32(T=-20+273.15, P=5e5)
            >>> print "%0.3f" % st.cv.kJkgK
            0.934
            >>> st=R32(T=155+273.15, P=2e4)
            >>> print "%0.3f" % st.cv.kJkgK
            0.873
            >>> st=R32(T=5+273.15, P=1e6)
            >>> print "%0.3f" % st.cv.kJkgK
            0.941
            >>> st=R32(T=80+273.15, P=6e6)
            >>> print "%0.3f" % st.cv.kJkgK
            1.411
            >>> st=R32(T=65+273.15, P=1e7)
            >>> print "%0.3f" % st.cv.kJkgK
            0.993
            """
            # Table 17, Pag 1325
            """
            >>> st=R32(T=-20+273.15, P=5e5)
            >>> print "%0.2f" % st.w
            805.61
            >>> st=R32(T=155+273.15, P=2e4)
            >>> print "%0.2f" % st.w
            284.42
            >>> st=R32(T=5+273.15, P=1e6)
            >>> print "%0.2f" % st.w
            669.69
            >>> st=R32(T=80+273.15, P=6e6)
            >>> print "%0.2f" % st.w
            150.04
            >>> st=R32(T=65+273.15, P=1e7)
            >>> print "%0.2f" % st.w
            439.83
            """,

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 435.0, "Pmax": 70000.0, "rhomax": 27.4734,
        "Pmin": 0.480e-1, "rhomin": 27.4734,

        "nr1": [0.1046634e1, -0.5451165, -0.2448595e-2, -0.4877002e-1,
                0.3520158e-1, 0.1622750e-2, 0.2377225e-4, 0.2914900e-1],
        "d1": [1, 2, 5, 1, 1, 3, 8, 4],
        "t1": [0.25, 1., -0.25, -1., 2., 2., 0.75, 0.25],

        "nr2": [0.3386203e-2, -0.4202444e-2, 0.4782025e-3, -0.5504323e-2,
                -0.2418396e-1, 0.4209034, -0.4616537, -0.1200513e1,
                -0.2591550e1, -0.1400145e1,  0.8263017],
        "d2": [4, 4, 8, 3, 5, 1, 1, 3, 1, 2, 3],
        "t2": [18., 26., -1., 25., 1.75, 4., 5., 1., 1.5, 1., 0.5],
        "c2": [4, 3, 1, 4, 1, 2, 2, 1, 1, 1, 1],
        "gamma2": [1]*11}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-32 of Span and Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"},
        "__test__": """
            >>> st=R32(T=700, rho=200, eq=2)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            1.1421 30.358 1.8392
            >>> st2=R32(T=750, rho=100, eq=2)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            235.85 0.59791
            """, # Table III, Pag 117

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 27.41,
        "Pmin": 0.047922, "rhomin": 27.41,

        "nr1": [.93080907, -.24777491e1, .41470439, .54859755e-1, .11475587e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.26225654, 0.41118822, 0.34970526e-2, -0.96790506e-1,
                -0.1172821, -0.4242838e-1, -0.12690083e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Astina and Sato (2003)",
        "__doi__": {"autor": "Astina, I.M. and Sato, H.",
                    "title": "A Rational Helmholtz Fundamental Equation of State for Difluoromethane with an Intermolecular Potential Background",
                    "ref": "Int. J. Thermophys., 34(4):963-990, 2003.",
                    "doi": "10.1023/A:1025096716493"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 72000.0, "rhomax": 27.48,
        "Pmin": 0.0485, "rhomin": 27.47,

        "nr1": [2.118688, -4.531096, 1.442456, 2.053906e-1, -1.311675e-1,
                1.022272e-2],
        "d1": [1, 1, 1, 3, 3, 4],
        "t1": [0.5, 1.125, 1.625, 0.875, 1.5, 1.75],

        "nr2": [4.873982e-1, -1.062213, -4.542051e-3, -6.933347e-4,
                -3.510307e-2, -5.606161e-2, 8.849625e-2, -1.850758e-2,
                7.878071e-3, -3.384115e-2, 1.641979e-4, -1.459172e-3],
        "d2": [1, 1, 5, 5, 6, 1, 2, 5, 6, 2, 2, 8],
        "t2": [1.75, 2.75, 0.25, 3.75, 1, 6.5, 2.5, 7.5, 7.5, 11, 16, 13],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-32 of Outcalt and McLinden (1995)",
        "__doi__": {"autor": "Outcalt, S.L. and McLinden, M.O.",
                    "title": "Equations of state for the thermodynamic properties of R32 (difluoromethane) and R125 (pentafluoroethane)",
                    "ref": "Int. J. Thermophysics, 16:79-89, 1995.",
                    "doi": "10.1007/BF01438959"},

        "R": 8.314471,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 27.48,
        "Pmin": 0.0477, "rhomin": 27.48,

        "b": [None, -0.131275405202e-3, 0.899927934911, -0.281400805178e2,
              0.436091182784e4, -0.837235280004e6, -0.782176408963e-6,
              -0.111226606825e1, 0.539331431878e3, 0.288600276863e6,
              -0.352264609289e-4, 0.189661830119, -0.686549003993e2,
              -0.349007064245e-2, -0.749983559476e-1, -0.321524283063e2,
              0.913057921906e-2, -0.171082181849e-3, 0.503986984347e-1,
              -0.830354867752e-3, -0.245522676708e6, -0.107859056038e8,
              -0.429514279646e4, 0.808724729567e8, -0.125945229993e2,
              -0.105735009761e4, -0.904064745354e-1, -0.183578733048e4,
              -0.169690612464e-3, 0.639250820631e-1, -0.204925767440e-6,
              -0.165629700870e-3, -0.932607493424e-2]}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-32 of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [2.75866232e-1, 9.26526641e-1, -2.44296579, 5.34289357e-2,
                1.06739638e-4, 3.46487335e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [9.07435007e-2, -1.93104843e-1, 5.11370826e-1, 3.09453923e-3,
                -1.53328967e-1, -1.03816916e-1, -3.8066998e-2, -1.16075825e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, helmholtz2, helmholtz3, helmholtz4
    _PR = 0.00585

    _surface = {"sigma": [0.07147], "exp": [1.246]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74883e1, 0.19697e1, -0.17496e1, -0.40224e1, 0.15209e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.12584e1, 0.46410e1, -0.54870e1, 0.33115e1, -0.61370],
        "exp": [0.27, 0.8, 1.1, 1.5, 1.8]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.22002e1, -.5972e1, -.14571e2, -.42598e2, .42686e1, -.73373e2],
        "exp": [0.336, 0.98, 2.7, 5.7, 6.5, 11.0]}

    thermo0 = {"eq": 1,
               "__name__": "Marsh (2002)",
               "__doc__": """Unpublished; however the fit uses the functional form found in: Marsh, K., Perkins, R., and Ramires, M.L.V., "Measurement and Correlation of the Thermal Conductivity of Propane from 86 to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data, 47(4):932-940, 2002""",

               "Tref": 351.255, "kref": 1,
               "no": [0.106548e-1, -0.194174e-1, 0.254295e-1],
               "co": [0, 1, 2],

               "Trefb": 351.255, "rhorefb": 8.1500846, "krefb": 1,
               "nb": [0.221878e-1, -0.215336e-1, 0.283523, -0.169164, -0.297237,
                      .191614, .105727, -.665397e-1, -.123172e-1, 0.766378e-2],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.582925e-10, "Tcref": 526.8825}

    _thermal = thermo0,
