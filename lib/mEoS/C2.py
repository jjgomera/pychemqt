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


class C2(MEoS):
    """Multiparameter equation of state for ethane"""
    name = "ethane"
    CASNumber = "74-84-0"
    formula = "CH3CH3"
    synonym = "R-170"
    rhoc = unidades.Density(206.18)
    Tc = unidades.Temperature(305.322)
    Pc = unidades.Pressure(4872.2, "kPa")
    M = 30.06904  # g/mol
    Tt = unidades.Temperature(90.368)
    Tb = unidades.Temperature(184.569)
    f_acent = 0.0995
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 3
    _Tr = unidades.Temperature(295.159630)
    _rhor = unidades.Density(207.557649)
    _w = 0.095234716

    Fi1 = {"R": 8.314472,
           "ao_log": [1, 3.003039265],
           "pow": [0, 1],
           "ao_pow": [9.212802589, -4.68224855],
           "ao_exp": [1.117433359, 3.467773215, 6.941944640, 5.970850948],
           "titao": [1.4091052332, 4.0099170712, 6.5967098342, 13.9798102659]}

    Fi2 = {"R": 8.31451,
           "ao_log": [1, 3.00263],
           "pow": [0, 1],
           "ao_pow": [24.675437527, -77.42531376],
           "ao_exp": [], "titao": [],
           "ao_hyp": [4.33939, 1.23722, 13.1974, -6.01989],
           "hyp": [559.314/Tc, 223.284/Tc, 1031.38/Tc, 1071.29/Tc]}

    Fi3 = {"ao_log": [1, 3.8159476],
           "pow": [0, -1./3, -2./3, -1],
           "ao_pow": [-23.446765, 8.6021299, -3.3075735, -.55956678],
           "ao_exp": [5.0722267], "titao": [5.5074874],
           "ao_hyp": [], "hyp": []}

    CP5 = {"ao": 9.9507922459,
           "an": [-6.9341406909e5, 3.1534834135e4, -6.103375287e2,
                  -2.8657877948e-2, 9.0922897821e-5, -5.2750109915e-8],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-1.4243593411e1], "exp": [3000],
           "ao_hyp": [], "hyp": []}

    buecker = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethane of Buecker and "
                    "Wagner (2006)",
        "__doi__": {"autor": "Bücker, D., Wagner, W.",
                    "title": "A Reference Equation of State for the "
                             "Thermodynamic Properties of Ethane for "
                             "Temperatures from the Melting Line to 675 K and "
                             "Pressures up to 900 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 35(1) (2006) 205-266",
                    "doi": "10.1063/1.1859286"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 900000.0, "rhomax": 22.419,
        "Pmin": 0.00114, "rhomin": 21.668,

        "nr1": [0.83440745735241, -0.14287360607171e1, 0.34430242210927,
                -0.42096677920265, 0.12094500886549e-1],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.25, 1.00, 0.25, 0.75, 0.75],

        "nr2": [-0.57976201597341, -0.33127037870838e-1, -0.11751654894130,
                -0.11160957833067, 0.62181592654406e-1, 0.98481795434443e-1,
                -0.098268582682358, -0.23977831007049e-3, 0.69885663328821e-3,
                0.19665987803305e-4, -0.014586152207928, 0.46354100536781e-1,
                0.60764622180645e-2, -0.26447330147828e-2, -0.042931872689904,
                0.29987786517263e-2, 0.52919335175010e-2, -0.10383897798198e-2,
                -0.54260348214694e-1, -0.21959362918493, 0.35362456650354,
                -0.12477390173714, 0.18425693591517, -0.16192256436754,
                -0.82770876149064e-1, 0.50160758096437e-1, 0.93614326336655e-2,
                -0.27839186242864e-3, 0.23560274071481e-4, 0.39238329738527e-2,
                -0.76488325813618e-3, -0.49944304440730e-2,
                0.18593386407186e-2, -0.61404353331199e-3],
        "d2": [1, 1, 2, 2, 3, 6, 6, 7, 9, 10, 2, 4, 4, 5, 5, 6, 8, 9, 2, 3, 3,
               3, 4, 4, 5, 5, 6, 11, 14, 3, 3, 4, 8, 10],
        "t2": [2.00, 4.25, 0.75, 2.25, 3.00, 1.00, 1.25, 2.75, 1.00, 2.00,
               2.50, 5.50, 7.00, 0.50, 5.50, 2.50, 4.00, 2.00, 10.00, 16.00,
               18.00, 20.00, 14.00, 18.00, 12.00, 19.00, 7.00, 15.00, 9.00,
               26.00, 28.00, 28.00, 22.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3,
               3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4],
        "gamma2": [1]*34,

        "nr3": [-0.23312179367924e-2, 0.29301047908760e-2,
                -0.26912472842883e-3, 184.13834111814, -10.397127984854],
        "d3": [1, 1, 3, 3, 2],
        "t3": [0., 3., 3., 0., 3.],
        "alfa3": [15, 15, 15, 20, 20],
        "beta3": [150, 150, 150, 275, 400],
        "gamma3": [1.05, 1.05, 1.05, 1.22, 1.16],
        "epsilon3": [1]*5}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethane of Younglove and Ely "
                    "(1987)",
        "__doi__": {"autor": "Younglove, B.A. and Ely, J.F.",
                    "title": "Thermophysical Properties of Fluids. II. "
                             "Methane, Ethane, Propane, Isobutane, and Normal "
                             "Butane",
                    "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                    "doi": "10.1063/1.555785"},

        "R": 8.31434,
        "cp": CP5,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 11874.2, "so": 229.116},

        "Tmin": 90.348, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 21.68,
        "Pmin": 1.1308e-3, "rhomin": 21.68,

        "b": [None, -0.3204748852e-2, 0.6529792241, -0.1669704591e2,
              0.1147983381e4, -0.1854721998e6, 0.4994149431e-3, -0.4858871291,
              0.1225345776e3, 0.8622615988e5, -0.1081290283e-4, 0.06279096996,
              -17.16912675, -0.1640779401e-3, -0.4356516111e-1, -19.66649699,
              0.4026724698e-2, -0.6498241861e-4, 0.05111594139,
              -0.1113010349e-2, -0.7157747547e4, -0.1848571024e8,
              -0.2137365569e4, 0.6275079986e8, -0.9974911056e1, 0.1129115014e4,
              -0.1026469558, -0.5660525915e4, -0.4209846430e-3, 0.2374523553,
              -0.1289637823e-5, -0.5423801068e-3, 0.2239717230e-1]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethane of Kunz and Wagner"
                    " (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 900000.0, "rhomax": 22.419,
        # "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.63596780450714, -0.17377981785459e1, 0.28914060926272,
                -0.33714276845694, 0.22405964699561e-1, 0.15715424886913e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.11450634253745, 0.10612049379745e1, -0.12855224439423e1,
                0.39414630777652, 0.31390924682041, -0.21592277117247e-1,
                -0.21723666564905, -0.28999574439489, 0.42321173025732,
                0.46434100259260e-1, -0.13138398329741, 0.11492850364368e-1,
                -0.33387688429909e-1, 0.015183171583644, -0.47610805647657e-2,
                0.46917166277885e-1, -0.039401755804649, -0.32569956247611e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1]*18}

    friend = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethane of Friend (1991)",
        "__doi__": {"autor": "Friend, D.G., Ingham, H., and Ely, J.F.",
                    "title": "Thermophysical Properties of Ethane",
                    "ref": "J. Phys. Chem. Ref. Data 20, 275 (1991)",
                    "doi": "10.1063/1.555881"},
        "__test__":
            # Table A1, Pag 336
            """
            >>> st=C2(T=500, P=1e5, eq=3)
            >>> print "%0.6g %0.1f %0.3f %0.3f %0.3f %0.3f %0.2f" % (\
                st.T, st.aM0.kJkmol, st.hM0.kJkmol, st.sM0.kJkmolK, st.cpM0.kJkmolK)
            500 -110.311 25.059 262.43 77.987
            """
            # Table A2, Pag 337
            """
            >>> st=C2(T=92, x=0.5, eq=3)
            >>> print "%0.6g %0.1e %0.2f %0.2e %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            92 1.7e-06 21.61 2.27e-06 67.74 1987.2 1193.00 254.4
            >>> st=C2(T=100, x=0.5, eq=3)
            >>> print "%0.6g %0.1e %0.2f %0.2e %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            100 1.1e-05 21.32 1.33e-05 70.09 1937.6 876.96 248.1
            >>> st=C2(T=150, x=0.5, eq=3)
            >>> print "%0.6g %0.1e %0.2f %0.2e %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            150 9.7e-3 19.47 7.80e-3 70.27 1573.2 270.35 201.0
            >>> st=C2(T=200, x=0.5, eq=3)
            >>> print "%0.6g %0.3f %0.2f %0.3f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            200 0.217 17.42 0.139 74.86 1194.4 138.17 152.5
            >>> st=C2(T=250, x=0.5, eq=3)
            >>> print "%0.6g %0.3f %0.2f %0.3f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            250 1.301 14.89 0.787 87.29 794.6 78.06 109.1
            >>> st=C2(T=300, x=0.5, eq=3)
            >>> print "%0.6g %0.3f %0.2f %0.3f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            300 4.356 10.10 3.813 182.06 278.4 35.01 71.3
            >>> st=C2(T=302, x=0.5, eq=3)
            >>> print "%0.6g %0.3f %0.2f %0.3f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            302 4.543 9.59 4.262 223.66 246.4 32.44 72.0
            >>> st=C2(T=304, rhom=8.82, eq=3)
            >>> print "%0.6g %0.3f %0.2f %0.3f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.cpM.JmolK, \
                st.Liquido.w, st.Liquido.mu.muPas, st.Liquido.k.mWmK)
            304 4.738 8.82 4.969 354.78 209.4 28.97 79.0
            """
            # Table A3, Pag 339
            """
            >>> st=C2(T=130, P=1e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            130 1 20.24 -12.071 102.03 45.01 70.10 1726.9 392.40 221.3
            >>> st=C2(T=140, P=6e7, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            140 60 20.80 -9.131 102.52 46.34 67.67 1921.7 476.29 245.7
            >>> st=C2(T=160, P=2e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            160 2 19.13 -9.933 116.48 43.04 70.44 1511.1 235.10 192.5
            >>> st=C2(T=180, P=1e5, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            180 0.1 18.28 -8.571 125.09 42.65 72.41 1347.8 176.42 171.5
            >>> st=C2(T=200, P=1e7, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            200 10 17.79 -6.804 131.51 43.41 73.00 1281.7 151.38 161.5
            >>> st=C2(T=240, P=1e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            240 1 15.47 -3.894 147.18 44.93 85.36 878.8 87.70 117.4
            >>> st=C2(T=270, P=2e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            270 2 1.20 8.589 194.29 47.40 76.57 245.2 9.33 21.6
            >>> st=C2(T=280, P=5e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            280 5 13.26 -0.228 160.21 48.73 103.93 603.7 57.96 90.7
            >>> st=C2(T=300, P=1e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            300 1 0.43 11.364 209.01 45.59 57.20 296.8 9.65 22.2
            >>> st=C2(T=330, P=5e5, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            330 0.5 0.19 13.366 220.86 48.51 57.89 320.8 10.37 25.6
            >>> st=C2(T=360, P=2e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            360 2 0.73 14.5 213.23 53.11 65.46 319.6 11.65 31.1
            >>> st=C2(T=400, P=5e6, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            400 5 1.77 16.051 210.58 59.05 76.57 322.4 13.91 40.0
            >>> st=C2(T=430, P=2e7, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            430 20 7.42 14.158 197.14 64.79 101.22 409.8 27.52 66.5
            >>> st=C2(T=480, P=1e5, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            480 0.1 0.03 23.500 259.25 67.28 75.67 385.8 14.28 50.1
            >>> st=C2(T=500, P=6e7, eq=3)
            >>> print "%0.6g %0.2g %0.2f %0.3f %0.2f %0.2f %0.2f %0.1f %0.2f %0.1f" % (\
                st.T, st.P.MPa, st.rhoM, st.hM.kJmol, st.sM.JmolK, st.cvM.JmolK, \
                st.cpM.JmolK, st.w, st.mu.muPas, st.k.mWmK)
            500 60 11.21 19.385 199.38 73.24 95.28 752.5 48.34 101.4
            """,

        "R": 8.31451,
        "cp": Fi3,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 11874, "so": 229.12},
        "Tt": 90.352, "Tc": 305.33, "Pc": 4871.8, "rhoc": 6.87, "M": 30.07,

        "Tmin": 90.352, "Tmax": 625.0, "Pmax": 70000.0, "rhomax": 22.419,
        "Pmin": 1.130e-3, "rhomin": 21.665,

        "nr1": [0.4621543056, -1.9236936387, 0.39878604003, 0.16054532372e-1,
                0.12895242219, 0.35458320491e-1, 0.34927844540e-1,
                -0.11306183380e-1, -0.39809032779e-1, 0.83031936834e-3,
                0.45921575183e-3, 0.17530287917e-6, -0.70919516126e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.23436162249, 0.84574697645e-1, 0.1486105201, -0.10016857867,
                -0.59264824388e-1, -0.41263514217e-1, 0.21855161869e-1,
                -0.74552720958e-4, -0.98859085572e-2, 0.10208416499e-2,
                -0.52189655847e-3, 0.98592162030e-4, 0.46865140856e-1,
                -0.19558011646e-1, -0.46557161651e-1, 0.32877905376e-2,
                0.13572090185, -0.10846471455, -0.67502836903e-2],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*19}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ethane of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi2,
        "ref": "OTO",
        "M": 30.07, "Tc": 305.322, "Pc": 4872.0, "rhoc": 206.6/30.07,

        "Tmin": 90.352, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 22.419,
        "Pmin": 0.0010902, "rhomin": 21.721,

        "nr1": [0.97628068, -0.26905251e1, 0.73498222, -0.35366206e-1,
                0.84692031e-1, 0.24154594e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.23964954, -0.42780093e-1, -0.22308832, -0.51799954e-1,
                -0.27178426e-1, 0.11246305e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethane of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 900000.0, "rhomax": 22.419,
        "Pmin": 0.00114, "rhomin": 21.668,

        "nr1": [1.32031629, 9.47177394e-1, -3.21919278, 7.47287278e-2,
                2.74919584e-4, -6.33952115e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-5.17685674e-2, 3.65838926e-2, 2.57753669e-1, -1.34856586e-2,
                -0.221551776, -6.89219870e-4, -4.47904791e-2, -2.15665728e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = buecker, MBWR, GERG, friend, shortSpan, sun

    _surface = {"sigma": [0.07602, -0.02912], "exp": [1.32, 1.676]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [11.1552, 0.0112], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [36.759, 23.639, -808.03, -378.84],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.75, 2.75]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.0011421,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [1, 1.05262374e8, -1.05262374e8], "exp1": [0, 2.55, 0],
                "a2": [2.23626315e8], "exp2": [1], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.48647577, 1.47010078, -1.66261122, 3.57898378, -4.79105705],
        "exp": [1, 1.5, 2.5, 3.5, 4]}
    _liquid_Density = {
        "eq": 4,
        "ao": [1.56138026, -0.381552776, 0.078537204, 0.0370315089],
        "exp": [0.987, 2, 4, 9.5]}
    _vapor_Density = {
        "eq": 6,
        "ao": [-1.89879145, -3.65459262, 0.850562745, 0.363965487, -1.50005943,
               -2.26690389],
        "exp": [1.038, 2.5, 3, 6, 9, 15]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.17067154, -0.48879666, 0.039038856],
              "__name__": "Friend (1991)",
              "__doi__": {"autor": "Friend, D.G., Ingham, H., and Ely, J.F.",
                          "title": "Thermophysical Properties of Ethane",
                          "ref": "J. Phys. Chem. Ref. Data 20, 275 (1991)",
                          "doi": "10.1063/1.555881"},

              "ek": 245.0, "sigma": 0.43682,
              "Tref": 1, "rhoref": 1.*M,
              "n_chapman": 0.1463897/M**0.5,

              "Tref_res": 305.33, "rhoref_res": 6.87*M, "etaref_res": 15.977,
              "n_num": [0.47177003, -0.23950311, 0.39808301, -0.27343335,
                        0.35192260, -0.21101308, -0.00478579, 0.07378129,
                        -0.030425255],
              "t_num": [0, -1, 0, -1, -1.5, 0, -2, 0, -1],
              "d_num": [1, 1, 2, 2, 2, 3, 3, 4, 4],
              "g_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "c_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "n_den": [1., -0.30435286, 0.001215675],
              "t_den": [0, 0, -1],
              "d_den": [0, 1, 1],
              "g_den": [0, 0, 0],
              "c_den": [0, 0, 0]}

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "Younglove (1987)",
              "__doi__": {"autor": "Younglove, B.A. and Ely, J.F.",
                          "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, Propane, Isobutane, and Normal Butane ",
                          "ref": "J. Phys. Chem. Ref. Data 16, 577 (1987)",
                          "doi": "10.1063/1.555785"},

              "ek": 240.0, "sigma": 0.440110,
              "n_chapman": 0.146388493/M**0.5,
              "F": [0.2102436247e1, -0.1065920192e1, 1.4, 305.33],
              "E": [-0.1903481042e2, 0.1799260494e4, 0.1561316986e2,
                    -0.1497221136e5, 0.1130374601, -0.2186440756e2,
                    0.8235954037e4],
              "rhoc": 6.875}

    visco2 = {"eq": 4, "omega": 1,
              "__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {"autor": "S.E.Quiñones-Cisneros and U.K. Deiters",
                          "title": "Generalization of the Friction Theory for Viscosity Modeling",
                          "ref": "J. Phys. Chem. B, 2006, 110 (25), pp 12820–12834",
                          "doi": "10.1021/jp0618577"},

              "Tref": 305.322, "muref": 1.0,
              "ek": 240.0, "sigma": 0.440110, "n_chapman": 0,
              "n_ideal": [15.9252, -49.7734, 43.4368],
              "t_ideal": [0, 0.25, 0.5],

              "a": [-7.50685764546476e-6, -1.50327318940575e-6, 5.58090793793288e-15],
              "b": [6.72861662009487e-5, -4.36450942982638e-5, -7.97441663817752e-14],
              "c": [3.88039503242230e-5, -1.38523739665972e-5, -2.64094611051755e-15],
              "A": [7.68043111364307e-10, -1.32047872761278e-10, 0.0],
              "B": [9.15406537766279e-9, 4.13028199950288e-10, 0.0],
              "C": [-1.45842039761136e-7, 2.39764228120527e-7, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Friend (1991)",
               "__doi__": {"autor": "Friend, D.G., Ingham, H., and Ely, J.F.",
                           "title": "Thermophysical Properties of Ethane",
                           "ref": "J. Phys. Chem. Ref. Data 20, 275 (1991)",
                           "doi": "10.1063/1.555881"},

               "Tref": 245.0, "kref": 1e-3,
               "no": [1.7104147, -0.6936482, 0],
               "co": [0, -1, -96],

               "Trefb": 305.33, "rhorefb": 6.87, "krefb": 4.41786e-3,
               "nb": [0.96084322, 2.7500235, -0.026609289, -0.078146729,
                      0.21881339, 2.3849563, -0.75113971],
               "tb": [0, 0, 0, 0, 0, -1.5, -1],
               "db": [1, 2, 3, 4, 5, 1, 3],
               "cb": [0]*7,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.242, "R0": 1.01,
               "Xio": 0.19e-9, "gam0": 0.0563, "qd": -0.545e-9, "Tcref": 610.66}

    _thermal = thermo0,


class Test(TestCase):

    def test_buecker(self):
        # Selected point from Table 29, Pag 238, saturation states
        st = C2(T=90.368, x=0.5)
        self.assertEqual(round(st.P.MPa, 7), 0.0000011)
        self.assertEqual(round(st.Liquido.rho, 5), 651.52948)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -888.90)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -5.058)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.605)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.326)
        self.assertEqual(round(st.Liquido.w, 2), 2008.69)
        self.assertEqual(round(st.Gas.rho, 6), 0.000046)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -294.12)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 1.524)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.892)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.168)
        self.assertEqual(round(st.Gas.w, 2), 180.93)

        st = C2(T=100, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000011)
        self.assertEqual(round(st.Liquido.rho, 5), 640.94852)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -866.74)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -4.825)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.541)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.283)
        self.assertEqual(round(st.Liquido.w, 2), 1938.44)
        self.assertEqual(round(st.Gas.rho, 5), 0.00040)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -282.78)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 1.015)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.911)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.187)
        self.assertEqual(round(st.Gas.w, 2), 189.86)

        st = C2(T=130, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.001284)
        self.assertEqual(round(st.Liquido.rho, 5), 607.82999)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -798.36)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -4.227)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.462)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.293)
        self.assertEqual(round(st.Liquido.w, 2), 1722.03)
        self.assertEqual(round(st.Gas.rho, 5), 0.03576)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -246.43)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 0.019)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.977)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.256)
        self.assertEqual(round(st.Gas.w, 2), 214.69)

        st = C2(T=150, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.009638)
        self.assertEqual(round(st.Liquido.rho, 5), 585.16884)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -752.12)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -3.896)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.442)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.333)
        self.assertEqual(round(st.Liquido.w, 2), 1575.53)
        self.assertEqual(round(st.Gas.rho, 5), 0.23373)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -221.71)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.360)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.027)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.312)
        self.assertEqual(round(st.Gas.w, 2), 228.84)

        st = C2(T=180, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.078638)
        self.assertEqual(round(st.Liquido.rho, 5), 549.50874)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -680.84)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -3.464)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.434)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.421)
        self.assertEqual(round(st.Liquido.w, 2), 1350.47)
        self.assertEqual(round(st.Gas.rho, 5), 1.62533)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -185.53)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.712)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.098)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.409)
        self.assertEqual(round(st.Gas.w, 2), 245.54)

        st = C2(T=210, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.333796)
        self.assertEqual(round(st.Liquido.rho, 5), 510.45075)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -605.90)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -3.081)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.454)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.572)
        self.assertEqual(round(st.Liquido.w, 2), 1117.27)
        self.assertEqual(round(st.Gas.rho, 5), 6.23900)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -153.48)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.927)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.228)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.622)
        self.assertEqual(round(st.Gas.w, 2), 254.02)

        st = C2(T=240, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.966788)
        self.assertEqual(round(st.Liquido.rho, 5), 465.30887)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -524.72)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -2.726)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.507)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.847)
        self.assertEqual(round(st.Liquido.w, 2), 873.25)
        self.assertEqual(round(st.Gas.rho, 5), 17.43487)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -128.82)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -1.077)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.388)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.976)
        self.assertEqual(round(st.Gas.w, 2), 252.14)

        st = C2(T=270, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 2.209980)
        self.assertEqual(round(st.Liquido.rho, 5), 407.71776)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -432.13)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -2.375)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.605)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 3.491)
        self.assertEqual(round(st.Liquido.w, 2), 608.92)
        self.assertEqual(round(st.Gas.rho, 5), 42.08922)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -118.38)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -1.212)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.595)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 2.815)
        self.assertEqual(round(st.Gas.w, 2), 237.02)

        st = C2(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 4.357255)
        self.assertEqual(round(st.Liquido.rho, 5), 303.50879)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -305.32)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -1.952)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.912)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 10.022)
        self.assertEqual(round(st.Liquido.w, 2), 274.91)
        self.assertEqual(round(st.Gas.rho, 5), 114.50091)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -155.61)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -1.453)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 2.089)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 13.299)
        self.assertEqual(round(st.Gas.w, 2), 200.51)

        st = C2(T=305, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 4.839225)
        self.assertEqual(round(st.Liquido.rho, 5), 241.96149)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -255.73)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -1.794)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 2.470)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 164.093)
        self.assertEqual(round(st.Liquido.w, 2), 175.12)
        self.assertEqual(round(st.Gas.rho, 5), 170.75482)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -202.19)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -1.619)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 2.623)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 247.460)
        self.assertEqual(round(st.Gas.w, 2), 178.83)

        # Table 30, Pag 243, Single phase points
        st = C2(T=90.384, P=1e5)
        self.assertEqual(round(st.rho, 2), 651.55)
        self.assertEqual(round(st.u.kJkg, 2), -888.88)
        self.assertEqual(round(st.h.kJkg, 2), -888.73)
        self.assertEqual(round(st.s.kJkgK, 4), -5.0574)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6051)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3256)
        self.assertEqual(round(st.w, 2), 2008.97)

        st = C2(T=135, P=5e5)
        self.assertEqual(round(st.rho, 2), 602.50)
        self.assertEqual(round(st.u.kJkg, 2), -787.09)
        self.assertEqual(round(st.h.kJkg, 2), -786.26)
        self.assertEqual(round(st.s.kJkgK, 4), -4.1415)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.4563)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3009)
        self.assertEqual(round(st.w, 2), 1688.21)

        st = C2(T=220, P=1e6)
        self.assertEqual(round(st.rho, 2), 497.12)
        self.assertEqual(round(st.u.kJkg, 2), -581.36)
        self.assertEqual(round(st.h.kJkg, 2), -579.35)
        self.assertEqual(round(st.s.kJkgK, 4), -2.9641)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.4681)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.6365)
        self.assertEqual(round(st.w, 2), 1044.02)

        st = C2(T=110, P=1.5e6)
        self.assertEqual(round(st.rho, 2), 630.62)
        self.assertEqual(round(st.u.kJkg, 2), -844.43)
        self.assertEqual(round(st.h.kJkg, 2), -842.05)
        self.assertEqual(round(st.s.kJkgK, 4), -4.6118)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.5041)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.2713)
        self.assertEqual(round(st.w, 2), 1872.62)

        st = C2(T=675, P=2e6)
        self.assertEqual(round(st.rho, 3), 10.756)
        self.assertEqual(round(st.u.kJkg, 2), 754.73)
        self.assertEqual(round(st.h.kJkg, 2), 940.67)
        self.assertEqual(round(st.s.kJkgK, 4), 1.1385)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.9468)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.2442)
        self.assertEqual(round(st.w, 2), 451.69)

        st = C2(T=290, P=3e6)
        self.assertEqual(round(st.rho, 3), 55.401)
        self.assertEqual(round(st.u.kJkg, 2), -153.13)
        self.assertEqual(round(st.h.kJkg, 3), -98.979)
        self.assertEqual(round(st.s.kJkgK, 4), -1.2014)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6618)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.0580)
        self.assertEqual(round(st.w, 2), 239.18)

        st = C2(T=160, P=1e7)
        self.assertEqual(round(st.rho, 2), 580.45)
        self.assertEqual(round(st.u.kJkg, 2), -734.04)
        self.assertEqual(round(st.h.kJkg, 2), -716.81)
        self.assertEqual(round(st.s.kJkgK, 4), -3.7788)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.4493)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3263)
        self.assertEqual(round(st.w, 2), 1563.69)

        st = C2(T=500, P=2e7)
        self.assertEqual(round(st.rho, 2), 164.96)
        self.assertEqual(round(st.u.kJkg, 2), 184.25)
        self.assertEqual(round(st.h.kJkg, 2), 305.49)
        self.assertEqual(round(st.s.kJkgK, 5), -0.56870)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.3996)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.2172)
        self.assertEqual(round(st.w, 2), 416.34)

        st = C2(T=100, P=5e7)
        self.assertEqual(round(st.rho, 2), 658.54)
        self.assertEqual(round(st.u.kJkg, 2), -877.76)
        self.assertEqual(round(st.h.kJkg, 2), -801.84)
        self.assertEqual(round(st.s.kJkgK, 4), -4.9448)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6011)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.2516)
        self.assertEqual(round(st.w, 2), 2107.34)

        st = C2(T=450, P=1e8)
        self.assertEqual(round(st.rho, 2), 428.87)
        self.assertEqual(round(st.u.kJkg, 2), -108.47)
        self.assertEqual(round(st.h.kJkg, 2), 124.70)
        self.assertEqual(round(st.s.kJkgK, 4), -1.4710)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.2729)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.9465)
        self.assertEqual(round(st.w, 2), 1075.84)

        st = C2(T=675, P=9e8)
        self.assertEqual(round(st.rho, 2), 632.88)
        self.assertEqual(round(st.u.kJkg, 2), 443.09)
        self.assertEqual(round(st.h.kJkg, 2), 1865.16)
        self.assertEqual(round(st.s.kJkgK, 5), -0.95311)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.2264)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.6380)
        self.assertEqual(round(st.w, 2), 2628.58)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = C2(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 3), 3.299)
        self.assertEqual(round(st.P.MPa, 3), 44.781)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.6276)

        st2 = C2(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 209.07)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.50714)
