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


class SF6(MEoS):
    """Multiparameter equation of state for sulfur hexafluoride"""
    name = "sulfur hexafluoride"
    CASNumber = "2551-62-4"
    formula = "SF6"
    synonym = ""
    rhoc = unidades.Density(742.3)
    Tc = unidades.Temperature(318.7232)
    Pc = unidades.Pressure(3754.983, "kPa")
    M = 146.0554192  # g/mol
    Tt = unidades.Temperature(223.555)
    Tb = unidades.Temperature(204.9)
    f_acent = 0.21
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 953
    id = 1
    _Tr = unidades.Temperature(304.013497)
    _rhor = unidades.Density(747.815849)
    _w = 0.181815238

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [11.638611086, -6.392241811],
           "ao_exp": [3.66118232, 7.87885103, 3.45981679],
           "titao": [1.617282065, 2.747115139, 4.232907175],
           "ao_hyp": [], "hyp": []}

    CP1 = {"ao": 3.9837756784,
           "an": [], "pow": [],
           "ao_exp": [2.2181851010, -1.0921337374e1, 3.3102497939,
                      17.5189671483, 2.8903523803],
           "exp": [1114.38, 925.64, 499.26, 884.9, 1363.93],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": -0.376915e-1/8.3143*146.05,
           "an": [0.305814e-2/8.3143*146.05, -0.237654e-5/8.3143*146.05],
           "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of Guder and Wagner (2009)",
        "__doi__": {"autor": "Guder, C. and Wagner, W.",
                    "title": "A Reference Equation of State for the Thermodynamic Properties of Sulfur Hexafluoride for Temperatures from the Melting Line to 625 K and Pressures up to 150 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 38, 33 (2009)",
                    "doi": "10.1063/1.3037344"},

        "__test__":
            # Table 27, Pag 57
            """
            >>> wt=SF6()
            >>> tau=wt.Tc/350
            >>> delta=436.9770888/wt.rhoc
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phi0(wt._constants["cp"], tau, delta)
            3.30559888 0.91277072 -14.4662979 1.69871606 -2.88563626 0
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phir(tau, delta)[:6]
            -0.496581463 -1.37926327 -1.37917096 -0.723171558 0.405086373 -2.09574715
            """
            # Table 28, Pag 71
            """
            >>> st=SF6(T=223.555, x=0.5)
            >>> print "%0.6g %0.6g %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            223.555 0.231425 1845.03 19.56 -158.14 -47.521 -0.72226 -0.22745 0.52753 0.48331 0.83712 0.56309 552.26 112.84
            >>> st=SF6(T=230, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            230 0.30133 1812.92 25.143 -152.66 -44.796 -0.69829 -0.2293 0.54069 0.49958 0.8573 0.58507 521.49 112.81
            >>> st=SF6(T=240, x=0.5)
            >>> print "%0.6g %0.6g %0.2f %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            240 0.440064 1760.80 36.21 -143.91 -40.63 -0.66135 -0.23103 0.56043 0.52495 0.89018 0.62193 474.67 112.27
            >>> st=SF6(T=250, x=0.5)
            >>> print "%0.6g %0.6g %0.2f %0.3f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            250 0.621964 1705.31 50.867 -134.81 -36.589 -0.62463 -0.23177 0.57964 0.55056 0.92638 0.66356 428.57 111.04
            >>> st=SF6(T=260, x=0.5)
            >>> print "%0.6g %0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            260 0.854649 1645.54 70.103 -125.33 -32.747 -0.58803 -0.23193 0.5986 0.57618 0.96814 0.71212 382.8 109.02
            >>> st=SF6(T=270, x=0.5)
            >>> print "%0.6g %0.7g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            270 1.146251 1580.17 95.322 -115.43 -29.213 -0.55134 -0.23202 0.61759 0.60262 1.0194 0.77387 336.89 106.11
            >>> st=SF6(T=280, x=0.5)
            >>> print "%0.6g %0.7g %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            280 1.505551 1507.05 128.67 -105.01 -26.139 -0.51431 -0.23262 0.63716 0.63271 1.0875 0.86231 290.16 102.15
            >>> st=SF6(T=290, x=0.5)
            >>> print "%0.6g %0.7g %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            290 1.942334 1422.51 173.85 -93.928 -23.766 -0.47647 -0.23453 0.6589 0.66504 1.1897 0.99972 241.62 96.919
            >>> st=SF6(T=300, x=0.5)
            >>> print "%0.6g %0.7g %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            300 2.468172 1319.00 238.43 -81.864 -22.596 -0.43688 -0.23932 0.68633 0.7029 1.3801 1.2722 189.8 90.145
            >>> st=SF6(T=310, x=0.5)
            >>> print "%0.6g %0.7g %0.2f %0.2f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            310 3.098401 1174.67 344.14 -67.924 -24.009 -0.39285 -0.25119 0.73084 0.76178 1.9613 2.1947 131.47 81.329
            >>> st=SF6(T=318, x=0.5)
            >>> print "%0.6g %0.7g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            318 3.695565 917.87 569.56 -50.87 -33.648 -0.34041 -0.28626 0.8803 0.93693 16.37 26.348 70.232 69.747
            >>> st=SF6(T=318.7, x=0.5)
            >>> print "%0.6g %0.7g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.6g %0.6g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            318.7 3.753053 787.75 696.67 -44.789 -40.391 -0.32152 -0.30773 1.0458 1.0721 1027.68 1273.3 60.578 62.816
            """
            # Table 29, Pag 73
            """
            >>> st=SF6(T=390, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            390 4.5244 44.587 66.69 0.19538 0.73041 0.78868 154.14
            >>> st=SF6(T=300, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 31.041 -18.062 -1.9544 -0.094212 0.61818 0.69206 130.24
            >>> st=SF6(T=225, P=1e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            225 1841.9 -157.3 -156.76 -0.71798 0.53068 0.83889 550.68
            >>> st=SF6(T=600, P=1.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            600 44.273 213.8 247.68 0.41151 0.87699 0.94036 189.91
            >>> st=SF6(T=400, P=2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 95.798 46.897 67.774 0.032259 0.74992 0.83896 146.21
            >>> st=SF6(T=305, P=3e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            305 1272.55 -78.218 -75.861 -0.41839 0.6977 1.4717 173.33
            >>> st=SF6(T=300, P=4e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 1382.35 -86.689 -83.795 -0.44709 0.6729 1.1916 234.55
            >>> st=SF6(T=625, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            625 142.37 231.78 266.9 0.37513 0.89076 0.9678 194.91
            >>> st=SF6(T=260, P=6e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            260 1693.94 -128.63 -125.08 -0.59892 0.59673 0.92285 436.07
            >>> st=SF6(T=300, P=1e7)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 1505.13 -92.926 -86.282 -0.46916 0.66145 1.0186 328.66
            >>> st=SF6(T=240, P=2e7)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            240 1862.82 -150.39 -139.65 -0.68852 0.56357 0.83019 603.24
            >>> st=SF6(T=400, P=5e7)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 1479.33 -19.467 14.333 -0.25764 0.77463 0.98607 437.03
            >>> st=SF6(T=500, P=1e8)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            500 1515.39 60.333 126.32 -0.082538 0.85386 1.0148 553.13
            >>> st=SF6(T=625, P=1.5e8)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            625 1505.25 171.72 271.37 0.11704 0.91497 1.0498 633.74
            """,

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 150000.0, "rhomax": 14.5,
        "Pmin": 231.429, "rhomin": 12.632,

        "nr1": [.54958259132835, -.87905033269396, -.84656969731452,
                .27692381593529, -.49864958372345e01, .48879127058055e01,
                .36917081634281e-1, .37030130305087e-3, .39389132911585e-1,
                .42477413690006e-3],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4, 6],
        "t1": [0.125, 1.25, 1.875, 0.125, 1.5, 1.625, 1.5, 5.625, 0.625, 0.25],

        "nr2": [-.24150013863890e-1,  .59447650642255e-1, -.38302880142267,
                .32606800951983, -.29955940562031e-1, -.86579186671173e-1,
                .41600684707562e01, -.41398128855814e01, -.55842159922714,
                .56531382776891,  .82612463415545e-2, -.10200995338080e-1],
        "d2": [1, 2, 2, 2, 3, 6, 2, 2, 4, 4, 2, 2],
        "t2": [6., 0.25, 4.75, 5.375, 5.875, 2., 5.875, 6., 5.625, 5.75, 0., 0.5],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*12,

        "nr3": [-.21662523861406e-1, .34650943893908e-1, -.28694281385812e-1,
                .84007238998053e-2, -.26969359922498, .90415215646344e01,
                -.37233103557977e01, -.27524670823704e04, .57711861697319e04,
                -.30234003119748e04, .22252778435360e07, -.23056065559032e07,
                .63918852944475e07, -.60792091415592e07],
        "d3": [1, 3, 4, 1, 1, 4, 3, 4, 4, 4, 1, 1, 3, 3],
        "t3": [4, 1, 3, 2, 4, 3, 4, 1, 2, 3, 3, 4, 3, 4],
        "alfa3": [10, 10, 10, 10, 11, 25, 30, 30, 30, 30, 30, 30, 30, 30],
        "beta3": [150, 150, 150, 150, 225, 300, 350, 350, 350, 350, 400, 400,
                  400, 400],
        "gamma3": [1.13, 1.13, 1.13, 1.16, 1.19, 1.19, 1.16, 1.16, 1.16, 1.16,
                   1.22, 1.22, 1.22, 1.22],
        "epsilon3": [0.85, 0.85, 0.85, 0.85, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of de Reuck et al. (1991)",
        "__doi__": {"autor": "de Reuck, K.M., Craven, R.J.B., and Cole, W.A.",
                    "title": "Report on the Development of an Equation of State for Sulphur Hexafluoride",
                    "ref": "IUPAC Thermodynamic Tables Project Centre, London, 1991.",
                    "doi": ""},

        "R": 8.31448,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 55000.0, "rhomax": 12.7,
        "Pmin": 224.36, "rhomin": 12.677,

        "nr1": [0.26945570453, -0.554046585076, -0.929624636454, 0.505661081063,
                -0.683495847809, 0.579161832426, -0.122636218956,
                -0.260339227668e-1, 0.222201648687e-1, -0.118992341472e-2,
                0.292000609763e-2, -0.243315775571e-2, 0.689778297550e-3],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 5, 10, 10, 10, 10],
        "t1": [0, 1.5, 2, 0, 1, 2, 0, 2, 0, 0.5, 1, 1.5, 2],

        "nr2": [-0.147585329235e1, 0.275952303526e1, -0.142721418498e1,
                0.598794196648e-1, 0.219991168025e-2, 0.746554473361e-2,
                0.345233637389e-2, -0.253226231963e-1, 0.433906886402e-1,
                -0.249349699078e-1, 0.338560952242e-2, 0.539985899700e-3],
        "d2": [2, 2, 2, 3, 7, 7, 9, 4, 4, 4, 6, 4],
        "t2": [3, 4, 5, 5, 1, 5, 1, 9, 14, 24, 24, 9],
        "c2": [2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 6],
        "gamma2": [1]*12}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for sulfur hexafluoride of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"},
        "__test__": """
            >>> st=SF6(T=700, rho=200, eq=2)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            0.9671 8.094 0.9958
            >>> st2=SF6(T=750, rho=100, eq=2)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            52.80 0.10913
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 12.65,
        "Pmin": 221.22, "rhomin": 12.645,

        "nr1": [0.12279403e1, -0.33035623e1, 0.12094019e1, -0.12316,
                0.11044657, 0.32952153e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.27017629, -0.62910351e-1, -0.3182889, -0.99557419e-1,
                -0.36909694e-1, 0.19136427e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*7}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of Polt et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 523.0, "Pmax": 40000.0, "rhomax": 13.133,
        "Pmin": 236.73, "rhomin": 12.712,

        "nr1": [0.131111896375, -0.792338803106, 0.580899809209,
                0.153233600406e1, -0.485096079094e1, 0.482411603806e1,
                -0.311285647219e1, 0.442141211276, 0.206313183222,
                -0.372305169645, 0.443536383059, -0.476354850910e-1,
                0.116313319336, 0.570240883234e-1, -0.152963195118,
                0.259842094503e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.131111896375, 0.792338803106, -0.580899809209,
                -0.744763581796, 0.204368923925e1, -0.129335324120e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.32678063]*6}

    eq = helmholtz1, helmholtz2, helmholtz3, helmholtz4

    _surface = {"sigma": [0.0538, -4.064e-5], "exp": [1.271, 0.2116]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.48475e-4,
                "Tmin": Tt, "Tmax": 800.0,
                "a1": [1., -30.0468473, 30.0468473, 359.771253, -359.771253],
                "exp1": [0, -20., 0, 3.25, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 2, "Tref": Tt, "Pref": 231.429,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [-11.6942141, 11.6942141], "exp1": [-1.07, 0],
                    "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 6,
        "ao": [-7.09634642, 1.676662, -2.3921599, 5.86078302, -9.02978735],
        "exp": [2.0, 3.0, 5.0, 8.0, 9.0]}
    _liquid_Density = {
        "eq": 6,
        "ao": [2.31174688, -1.12912486, -1.439347, 0.282489982],
        "exp": [1.065, 1.5, 4., 5.]}
    _vapor_Density = {
        "eq": 6,
        "ao": [23.68063442, 0.513062232, -24.4706238, -4.6715244, -1.7536843,
               -6.65585369],
        "exp": [1.044, 0.5, 1.0, 2.0, 8.0, 17.]}
