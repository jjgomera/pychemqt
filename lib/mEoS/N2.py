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


from math import exp
from unittest import TestCase

from scipy.constants import Boltzmann, Avogadro

from lib import unidades
from lib.meos import MEoS


class N2(MEoS):
    """Multiparamente equation of state for nitrogen"""
    name = "nitrogen"
    CASNumber = "7727-37-9"
    formula = "N2"
    synonym = "R-728"
    _refPropName = "NITROGEN"
    _coolPropName = "Nitrogen"
    rhoc = unidades.Density(313.299958972)
    Tc = unidades.Temperature(126.192)
    Pc = unidades.Pressure(3395.8, "kPa")
    M = 28.01348  # g/mol
    Tt = unidades.Temperature(63.151)
    Tb = unidades.Temperature(77.355)
    f_acent = 0.0372
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 46
    _Tr = unidades.Temperature(122.520245)
    _rhor = unidades.Density(316.134310)
    _w = 0.043553140

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [-12.76952708, -0.00784163, -1.934819e-4,
                      -1.247742e-5, 6.678326e-8],
           "ao_exp": [1.012941],
           "titao": [26.65788]}

    Fi2 = {"ao_log": [1, 2.50031],
           "pow": [0, 1],
           "ao_pow": [11.083407489, -22.202102428],
           "ao_sinh": [0.13732, 0.90066], "sinh": [5.25182262, 13.788988208],
           "ao_cosh": [-0.1466], "cosh": [-5.393067706]}

    CP1 = {"ao": 3.50404228308756,
           "an": [-0.735210401157252e3, 0.342239980411978e2, -0.55764828456762,
                  -1.73390185081005e-5, 1.74650849766463e-8,
                  -3.56892033544348e-12],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [1.00538722808834], "exp": [3353.4061]}

    CP2 = {"ao": 3.50418363823,
           "an": [-0.837079888737e3, 0.379147114487e2, -0.601737844275,
                  -0.874955653028e-5, 0.148958507239e-7, -0.256370354277e-11],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [1.00773735767], "exp": [3353.4061]}

    span = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Span (2000)",
        "__doi__": {"autor": "Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, "
                             "W., Yokozeki, A.",
                    "title": "A Reference Equation of State for the "
                             "Thermodynamic Properties of Nitrogen for "
                             "Temperatures from 63.151 to 1000 K and "
                             "Pressures to 2200 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29(6) (2000) 1361-1433",
                    "doi":  "10.1063/1.1349047"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101325., "ho": 8670, "so": 191.5},

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 2200000.0, "rhomax": 53.15,

        "nr1": [0.924803575275, -0.492448489428, 0.661883336938,
                -0.192902649201e1, -0.622469309629e-1, 0.349943957581],
        "d1": [1, 1, 2, 2, 3, 3],
        "t1": [0.25, 0.875, 0.5, 0.875, 0.375, 0.75],

        "nr2": [0.564857472498, -0.161720005987e1, -0.481395031883,
                0.421150636384, -0.161962230825e-1, 0.172100994165,
                0.735448924933e-2, 0.168077305479e-1, -0.107626664179e-2,
                -0.137318088513e-1, 0.635466899859e-3, 0.304432279419e-2,
                -0.435762336045e-1, -0.723174889316e-1, 0.389644315272e-1,
                -0.212201363910e-1, 0.408822981509e-2, -0.551990017984e-4,
                -0.462016716479e-1, -0.300311716011e-2, 0.368825891208e-1,
                -0.255856846220e-2, 0.896915264558e-2, -0.441513370350e-2,
                0.133722924858e-2, 0.264832491957e-3],
        "d2": [1, 1, 1, 3, 3, 4, 6, 6, 7, 7, 8, 8, 1, 2, 3, 4, 5, 8, 4, 5, 5,
               8, 3, 5, 6, 9],
        "t2": [0.5, 0.75, 2., 1.25, 3.5, 1., 0.5, 3., 0., 2.75, 0.75, 2.5, 4.,
               6., 6., 3., 3., 6., 16., 11., 15., 12., 12., 7., 4., 16.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3,
               3, 4, 4, 4, 4],
        "gamma2": [1]*26,

        "nr3": [0.196688194015e2, -0.209115600730e2, 0.167788306989e-1,
                0.262767566274e4],
        "d3": [1, 1, 3, 2],
        "t3": [0., 1., 2., 3.],
        "alfa3": [20, 20, 15, 25],
        "beta3": [325, 325, 300, 275],
        "gamma3": [1.16, 1.16, 1.13, 1.25],
        "epsilon3": [1]*4}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for nitrogen of Younglove (1982).",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, "
                             "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                             "Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                    "doi": ""},

        "R": 8.31434,
        "M": 28.013, "Tt": 63.15, "Tc": 126.26, "rhoc": 11.21, "Pc": 3399.08,
        "cp": CP1,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 8716.2, "so": 191.666},

        "Tmin": 63.15, "Tmax": 1900.0, "Pmax": 1013000.0, "rhomax": 32,

        "gamma": -0.0056,
        "b": [None, 0.1380297474657e-2, 0.1084506501349, -0.2471324064362e1,
              0.3455257980807e2, -0.4279707690666e4, 0.1064911566998e-3,
              -0.1140867079735e-1, 0.1444902497287e-3, 0.1871457567553e5,
              0.8218876886831e-7, 0.2360990493348e-2, -0.5144803081201,
              0.4914545013668e-4, -0.1151627162399e-2, -0.7168037246650,
              0.7616667619500e-4, -0.1130930066213e-5, 0.3736831166831e-3,
              -0.2039851507581e-5, -0.1719662008990e5, -0.1213055199748e6,
              -0.9881399141428e2, 0.5619886893511e5, -0.1823043964118,
              -0.2599826498477e1, -0.4191893423157e-3, -0.2596406670530,
              -0.1258683201921e-6, 0.1049286599400e-4, -0.5458369305152e-9,
              -0.7674511670597e-8, 0.5931232870994e-7]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 2200000.0, "rhomax": 53.15,

        "nr1": [0.59889711801201, -0.16941557480731e1, 0.24579736191718,
                -0.23722456755175, 0.17954918715141e-1, 0.14592875720215e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.10008065936206, 0.73157115385532, -0.88372272336366,
                0.31887660246708, 0.20766491728799, -0.19379315454158e-1,
                -0.16936641554983, 0.13546846041701, -0.33066712095307e-1,
                -0.60690817018557e-1, 0.12797548292871e-1, 0.58743664107299e-2,
                -0.018451951971969, 0.47226622042472e-2, -0.52024079680599e-2,
                0.043563505956635, -0.36251690750939e-1, -0.28974026866543e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1]*18,

        "nr3": [],
        "nr4": []}

    jacobsen = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Jacobsen et "
                    "al. (1986).",
        "__doi__": {"autor": "Jacobsen, R.T, Stewart, R.B., Jahangiri, M.",
                    "title": "Thermodynamic Properties of Nitrogen from the "
                             "Freezing Line to 2000K at Pressures to 1000MPa",
                    "ref": "J. Phys. Chem. Ref. Data, 15(2) (1986) 735-908",
                    "doi": "10.1007/BF00502385"},

        "R": 8.31434,
        "cp": CP2,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 8716.48, "so": 191.66},
        "Tc": 126.193, "Pc": 3397.8, "rhoc": 11.177, "M": 28.0134,
        "Tt": 63.148,

        "Tmin": 63.15, "Tmax": 2000.0, "Pmax": 1000000.0, "rhomax": 30.96,

        "nr1": [0.9499541827, 0.2481718513, -0.2046287122, -0.1748429008,
                0.6387017148, -0.5272986168, -0.2049741504e1, 0.5551383553e-1,
                -0.8191106396e-3, 0.2650110798, 0.7311459372e-1,
                -0.2813080718e-1, 0.1659823569e-2, -0.3785445194, 0.1895290433,
                -0.7001895093e-2],
        "d1": [1, 2, 3, 2, 3, 3, 1, 4, 6, 1, 2, 4, 6, 1, 2, 4],
        "t1": [.25, .25, .25, .5, .5, .75, 1, 1, 1, 1.5, 2, 2, 2, 3, 3, 3],

        "nr2": [-0.5032519699e-1, 0.6012817812e-1, -0.4927710927e-1,
                0.6512013679e-1, 0.1138121942, -0.955140963197e-1,
                0.211835414e-1, -0.1100721771e-1, 0.128443221e-1,
                -0.105447491e-1, -0.1484600538e-3, -0.5806483467e-2],
        "d2": [2, 2, 1, 4, 1, 2, 4, 2, 4, 4, 2, 3],
        "t2": [1, 2, 3, 4, 4, 5, 6, 8, 14, 18, 20, 22],
        "c2": [2, 2, 3, 2, 3, 2, 2, 4, 4, 4, 4, 3],
        "gamma2": [1]*12}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for nitrogen of Span "
                    "and Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101325., "ho": 8670, "so": 191.5},
        "M": 28.013, "Tc": 126.192, "rhoc": 313.3/28.013,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 53.15,

        "nr1": [0.92296567, -0.25575012e1, 0.64482463, 0.1083102e-1,
                0.73924167e-1, 0.23532962e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18024854, -0.45660299e-1, -0.1552106, -0.3811149e-1,
                -0.31962422e-1, 0.15513532e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101325., "ho": 8670, "so": 191.5},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [9.57664698e-1, 8.68692283e-1, -2.88536117, 6.12953165e-2,
                2.55919463e-4, 1.69423647e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-4.43639900e-2, 1.37987734e-1, 2.77148365e-1, -1.44381707e-2,
                -1.69955805e-1, 5.46894457e-3, -2.87747274e-2, -2.38630424e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = span, younglove, GERG, jacobsen, shortSpan, sun
    _PR = [-0.7468, -11.9697]

    _surface = {"sigma": [0.02898], "exp": [1.246]}
    _dielectric = {
        "eq": 1,
        "a": [4.3872, 0.00226], "b": [2.206, 1.135], "c": [-169.0, -35.83],
        "Au": 0, "D": 2.1}

    _melting = {
        "eq": 1,
        "__doi__": span["__doi__"],

        "Tmin": Tt, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 12523,
        "a0": 1,
        "a2": [12798.61], "exp2": [1.78963]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.12445284, 1.2632722, -0.765910082, -1.77570564],
        "t": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 2,
        "n": [1.48654237, -0.280476066, 0.894143085e-1, -0.119879866],
        "t": [0.3294, 2/3, 8/3, 35/6]}
    _vapor_Density = {
        "eq": 3,
        "n": [-1.70127164, -3.70402649, 1.29859383, -0.561424977, -2.68505381],
        "t": [0.34, 5/6, 7/6, 13/6, 14/3]}

    visco0 = {"__name__": "Lemmon (2004)",
              "__doi__": {
                  "autor": "Lemmon, E.W., Jacobsen, R.T.",
                  "title": "Viscosity and Thermal Conductivity Equations for "
                           "Nitrogen, Oxygen, Argon, and Air",
                  "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                  "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

              "eq": 1, "omega": 1,
              "ek": 98.94, "sigma": 0.3656,

              "Tref_res": 126.192, "rhoref_res": 11.1839*M,
              "nr": [10.72, 0.03989, 0.001208, -7.402, 4.62],
              "tr": [.1, .25, 3.2, .9, 0.3],
              "dr": [2, 10, 12, 2, 1],
              "gr": [0, 1, 1, 1, 1],
              "cr": [0, 1, 1, 2, 3]}

    visco1 = {"__name__": "Younglove (1982)",
              "__doi__": {
                  "autor": "Younglove, B.A.",
                  "title": "Thermophysical Properties of Fluids. I. Argon, "
                           "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                           "Trifluoride, and Oxygen",
                  "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                  "doi": ""},

              "eq": 2, "omega": 0,

              "no": [-0.1822424e5, 0.19915327374e5, -0.91542324494e4,
                     0.23255484059e4, -0.36307214228e3, 0.36457506811e2,
                     -0.22261880817e1, 0.78053904895e-1, -0.11894029104e-2],
              "to": [-1, -2/3, -1/3, 0, 1/3, 2/3, 1, 4/3, 5/3],

              "mod": True,
              "F": [-0.11217739623, 0.32912317244e-1, 1.4, 118],
              "E": [-12.128154129, 0.57156092139, 16.094611148, 3.6954086158e3,
                    -8.088980118e2, 68.46443564, -2.1241135912],
              "rhoc": 0.315*1000/28.013}

    visco2 = {"__name__": "Stephan (1987)",
              "__doi__": {
                  "autor": "Stephan, K., Krauss, R., Laesecke, A.",
                  "title": "Viscosity and Thermal Conductivity of Nitrogen "
                           "for a Wide Range of Fluid States",
                  "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 993-1023",
                  "doi": "10.1063/1.555798"},

              "eq": 1, "omega": 1,

              "ek": 100.01654, "sigma": 0.36502496,
              "n_chapman": 0.141290/M**0.5,
              "collision": [0.46649, -0.57015, 0.19164, -0.03708,  0.00241],

              "rhoref_res": 314, "muref_res": 14.,
              "nr": [-5.8470232, -1.4470051, -0.027766561, -0.21662362],
              "tr": [0, 0, 0, 0],
              "dr": [0, 1, 2, 3],
              "nr_num": [-20.09997],
              "tr_num": [0],
              "dr_num": [0],
              "nr_den": [1.0, -3.4376416],
              "tr_den": [0, 0],
              "dr_den": [1, 0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Lemmon (2004)",
               "__doi__": {
                   "autor": "Lemmon, E.W., Jacobsen, R.T.",
                   "title": "Viscosity and Thermal Conductivity Equations for "
                            "Nitrogen, Oxygen, Argon, and Air",
                   "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                   "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

               "eq": 1,

               "Toref": 126.192, "koref": 1e-3,
               "no_visco": 1.511,
               "no": [2.117, -3.332],
               "to": [1, 0.7],

               "Tref_res": 126.192, "rhoref_res": 11.1839*M, "kref_res": 1e-3,
               "nr": [8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672],
               "tr": [0, 0.03, 0.2, 0.8, 0.6, 1.9],
               "dr": [1, 2, 3, 4, 8, 10],
               "cr": [0, 0, 1, 2, 2, 2],
               "gr": [0, 0, 1, 1, 1, 1],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.17e-9, "gam0": 0.055, "qd": 0.40e-9, "Tcref": 252.384}

    thermo1 = {"__name__": "Younglove (1982)",
               "__doi__": {
                   "autor": "Younglove, B.A.",
                   "title": "Thermophysical Properties of Fluids. I. Argon, "
                            "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                            "Trifluoride, and Oxygen",
                   "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                   "doi": ""},

               "eq": 2,

               "no": [-0.20029573972e2, 0.49765746684e1, 0.80188959378e1,
                      -0.55022716888e1, 0.15363738965e1, -0.22974737257,
                      0.19360547346e-1, -0.85677385768e-3, 0.15564670935e-4],
               "to": [-1., -2/3, -1/3, 0, 1/3, 2/3, 1., 4/3, 5/3],

               "F": [0.53875666637e-1, 0.61027911104e-2, 0.12e1, 0.118e3],
               "E": [-0.38613291627e2, 0, 0.37201743333e2, 0, -0.39013509079e2,
                     -0.31826109485e2, 0, 1],

               "critical": 1,
               "rhoc": 0.3139, "Tc": 126.24,
               "ek": 118, "f": 1.67108, "gm": 3.933e-10}

    thermo2 = {"__name__": "Stephan (1987)",
               "__doi__": {
                   "autor": "Stephan, K., Krauss, R., Laesecke, A.",
                   "title": "Viscosity and Thermal Conductivity of Nitrogen "
                            "for a Wide Range of Fluid States",
                   "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 993-1023",
                   "doi": "10.1063/1.555798"},

               "eq": 1, "critical": 0,

               "special": "_thermo0",
               # "Toref": 1, "koref": 1e-3,
               # "no": [0.6950401, 0.03643102],
               # "to": [-97, -98],

               "rhoref_res": 314, "kref_res": 4.17e-3,
               "nr": [3.3373542, 0.37098251, 0.89913456, 0.16972505],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4]}

    def _thermo0(self, rho, T, fase):
        """Custom dilute-gas limit form of thermal conductivity"""
        X1 = 0.95185202
        X2 = 1.0205422
        M = 28.013

        # Used the ideal isochoric heat capacity of paper because differ in
        # about a 20% of values using the ideal correlation in meos
        def cv(T):
            ni = [-0.837079888737e3, 0.379147114487e2, -0.601737844275,
                  0.350418363823e1, -0.874955653028e-5, 0.148968607239e-7,
                  -0.256370354277e-11, 0.100773735767e1, 0.335340610e4]
            sum1 = 0
            for i, n in enumerate(ni[:-2]):
                sum1 += n*T**(i-3)

            u = ni[8]/T
            eu1 = exp(u)-1
            sum2 = (ni[7]*u**2*(eu1+1))/eu1**2
            cv = 8.31434*(sum1+sum2-1)
            return cv

        muo = self._Visco0(T, self.visco2)
        F = Boltzmann*Avogadro*muo/M                                     # Eq 9
        ltr = 2.5*F*(1.5-X1)                                             # Eq 7
        lint = F*X2*(cv(T)/Boltzmann/Avogadro+X1)                        # Eq 8

        return (ltr+lint)*1e-3

    _thermal = thermo0, thermo1, thermo2


class Test(TestCase):

    def test_Span(self):
        # Selected point of pag 1403, saturation state
        st = N2(T=63.151, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.01252)
        self.assertEqual(round(st.Liquido.rhoM, 3), 30.957)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -4222.6)
        self.assertEqual(round(st.Liquido.sM.JmolK, 3), 67.951)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 32.95)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 56.03)
        self.assertEqual(round(st.Liquido.w, 1), 995.3)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.02407)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 1814.7)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 163.55)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 21.01)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 29.65)
        self.assertEqual(round(st.Gas.w, 1), 161.1)

        st = N2(T=80, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.13687)
        self.assertEqual(round(st.Liquido.rhoM, 3), 28.341)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -3265.7)
        self.assertEqual(round(st.Liquido.sM.JmolK, 3), 81.317)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 29.95)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 57.58)
        self.assertEqual(round(st.Liquido.w, 1), 824.4)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.21737)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2215.8)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 149.84)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 21.78)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 32.07)
        self.assertEqual(round(st.Gas.w, 1), 176.7)

        st = N2(T=100, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.77827)
        self.assertEqual(round(st.Liquido.rhoM, 3), 24.608)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -2050.8)
        self.assertEqual(round(st.Liquido.sM.JmolK, 3), 94.576)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 27.54)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 64.93)
        self.assertEqual(round(st.Liquido.w, 1), 605.2)
        self.assertEqual(round(st.Gas.rhoM, 4), 1.1409)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2458.6)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 139.67)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 23.95)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 42.09)
        self.assertEqual(round(st.Gas.w, 1), 183.3)

        st = N2(T=120, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 2.51058)
        self.assertEqual(round(st.Liquido.rhoM, 3), 18.682)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), -500.60)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 107.89)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 28.31)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 1), 126.3)
        self.assertEqual(round(st.Liquido.w, 1), 317.3)
        self.assertEqual(round(st.Gas.rhoM, 4), 4.4653)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2077.8)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 129.38)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 30.77)
        self.assertEqual(round(st.Gas.cpM.JmolK, 1), 129.7)
        self.assertEqual(round(st.Gas.w, 1), 172.6)

        st = N2(T=126, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 3.36453)
        self.assertEqual(round(st.Liquido.rhoM, 3), 13.281)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), 492.37)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 115.50)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 43.40)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 0), 3138)
        self.assertEqual(round(st.Liquido.w, 1), 151.0)
        self.assertEqual(round(st.Gas.rhoM, 4), 9.1106)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 1194.9)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 121.08)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 47.44)
        self.assertEqual(round(st.Gas.cpM.JmolK, 0), 4521)
        self.assertEqual(round(st.Gas.w, 1), 148.4)

        st = N2(P=101325, x=0.5)
        self.assertEqual(round(st.T, 3), 77.355)
        self.assertEqual(round(st.Liquido.rhoM, 3), 28.775)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -3418.2)
        self.assertEqual(round(st.Liquido.sM.JmolK, 3), 79.395)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 30.37)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 57.19)
        self.assertEqual(round(st.Liquido.w, 1), 851.4)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.16464)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2161.5)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 151.53)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 21.61)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 31.49)
        self.assertEqual(round(st.Gas.w, 1), 174.8)

        # Selected point from pag 1410, single phase region
        st = N2(T=80, P=1e5)
        self.assertEqual(round(st.rhoM, 5), 0.15633)
        self.assertEqual(round(st.uM.Jmol, 1), 1605.7)
        self.assertEqual(round(st.hM.Jmol, 1), 2245.4)
        self.assertEqual(round(st.sM.JmolK, 2), 152.70)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.48)
        self.assertEqual(round(st.cpM.JmolK, 2), 31.15)
        self.assertEqual(round(st.w, 1), 178.3)

        st = N2(T=90, P=5e5)
        self.assertEqual(round(st.rhoM, 3), 26.615)
        self.assertEqual(round(st.uM.Jmol, 1), -2692.7)
        self.assertEqual(round(st.hM.Jmol, 1), -2673.9)
        self.assertEqual(round(st.sM.JmolK, 3), 88.130)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.57)
        self.assertEqual(round(st.cpM.JmolK, 2), 59.87)
        self.assertEqual(round(st.w, 1), 720.8)

        st = N2(T=63.368, P=1e6)
        self.assertEqual(round(st.rhoM, 3), 30.986)
        self.assertEqual(round(st.uM.Jmol, 1), -4220.4)
        self.assertEqual(round(st.hM.Jmol, 1), -4188.1)
        self.assertEqual(round(st.sM.JmolK, 3), 67.993)
        self.assertEqual(round(st.cvM.JmolK, 2), 32.98)
        self.assertEqual(round(st.cpM.JmolK, 2), 55.89)
        self.assertEqual(round(st.w, 1), 998.9)

        st = N2(T=1000, P=2e6)
        self.assertEqual(round(st.rhoM, 5), 0.23889)
        self.assertEqual(round(st.uM.Jmol, 0), 21800)
        self.assertEqual(round(st.hM.Jmol, 0), 30172)
        self.assertEqual(round(st.sM.JmolK, 2), 203.24)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.40)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.75)
        self.assertEqual(round(st.w, 1), 635.5)

        st = N2(T=100, P=5e6)
        self.assertEqual(round(st.rhoM, 3), 25.436)
        self.assertEqual(round(st.uM.Jmol, 1), -2217.6)
        self.assertEqual(round(st.hM.Jmol, 1), -2021.0)
        self.assertEqual(round(st.sM.JmolK, 3), 93.188)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.71)
        self.assertEqual(round(st.cpM.JmolK, 2), 59.87)
        self.assertEqual(round(st.w, 1), 673.2)

        st = N2(T=70, P=1e7)
        self.assertEqual(round(st.rhoM, 3), 30.605)
        self.assertEqual(round(st.uM.Jmol, 1), -3944.5)
        self.assertEqual(round(st.hM.Jmol, 1), -3617.7)
        self.assertEqual(round(st.sM.JmolK, 3), 72.168)
        self.assertEqual(round(st.cvM.JmolK, 2), 32.34)
        self.assertEqual(round(st.cpM.JmolK, 2), 54.75)
        self.assertEqual(round(st.w, 1), 989.6)

        st = N2(T=350, P=2e7)
        self.assertEqual(round(st.rhoM, 4), 6.3552)
        self.assertEqual(round(st.uM.Jmol, 1), 6432.4)
        self.assertEqual(round(st.hM.Jmol, 1), 9579.5)
        self.assertEqual(round(st.sM.JmolK, 2), 150.07)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.64)
        self.assertEqual(round(st.cpM.JmolK, 2), 34.19)
        self.assertEqual(round(st.w, 1), 449.7)

        st = N2(T=1000, P=5e7)
        self.assertEqual(round(st.rhoM, 4), 5.1051)
        self.assertEqual(round(st.uM.Jmol, 0), 21447)
        self.assertEqual(round(st.hM.Jmol, 0), 31241)
        self.assertEqual(round(st.sM.JmolK, 2), 176.14)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.79)
        self.assertEqual(round(st.cpM.JmolK, 2), 33.68)
        self.assertEqual(round(st.w, 1), 748.9)

        st = N2(T=100, P=1e8)
        self.assertEqual(round(st.rhoM, 3), 31.793)
        self.assertEqual(round(st.uM.Jmol, 1), -3136.1)
        self.assertEqual(round(st.hM.Jmol, 2), 9.22)
        self.assertEqual(round(st.sM.JmolK, 3), 81.024)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.89)
        self.assertEqual(round(st.cpM.JmolK, 2), 47.93)
        self.assertEqual(round(st.w, 0), 1208)

        st = N2(T=1000, P=1e9)
        self.assertEqual(round(st.rhoM, 3), 30.189)
        self.assertEqual(round(st.uM.Jmol, 0), 22346)
        self.assertEqual(round(st.hM.Jmol, 0), 55471)
        self.assertEqual(round(st.sM.JmolK, 2), 149.55)
        self.assertEqual(round(st.cvM.JmolK, 2), 29.21)
        self.assertEqual(round(st.cpM.JmolK, 2), 36.45)
        self.assertEqual(round(st.w, 0), 2002)

    def test_younglove(self):
        kw = {"eq": "younglove", "visco": 1, "thermal": 1}

        # Selected point from Appendix I, Pag 1-162, saturation states
        # The pressure and density in saturation is calculate in tables using
        # the ancillary equation used in paper so the calculated point differ
        # of implement eq
        st = N2(T=N2.Tt, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 5), 0.01268)
        self.assertEqual(round(st.Liquido.rho, 1), 867.8)
        self.assertEqual(round(st.Liquido.rhoM, 2), 30.98)
        self.assertEqual(round(st.Liquido.uM.Jmol, 0), -4212)
        self.assertEqual(round(st.Liquido.hM.Jmol, 0), -4212)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 68.01)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 26.69)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 54.68)
        self.assertEqual(round(st.Liquido.w, 0), 1326)
        self.assertEqual(round(st.Liquido.mu.muPas, 0), 283)
        self.assertEqual(round(st.Liquido.k.WmK, 3), 0.151)
        self.assertEqual(round(st.Gas.rho, 4), 0.6829)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.02438)
        self.assertEqual(round(st.Gas.uM.Jmol, 0), 1296)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 1816)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 163.5)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 21.00)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 29.65)
        self.assertEqual(round(st.Gas.w, 1), 161.1)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 4.22)
        self.assertEqual(round(st.Gas.k.WmK, 5), 0.00569)

        st = N2(T=110, x=0.5, **kw)
        self.assertEqual(round(st.P.MPa, 3), 1.467)
        self.assertEqual(round(st.Liquido.rho, 1), 620.0)
        self.assertEqual(round(st.Liquido.rhoM, 2), 22.13)
        self.assertEqual(round(st.Liquido.uM.Jmol, 0), -1412)
        self.assertEqual(round(st.Liquido.hM.Jmol, 0), -1346)
        self.assertEqual(round(st.Liquido.sM.JmolK, 1), 101.0)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 26.21)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 75.93)
        self.assertEqual(round(st.Liquido.w, 1), 483.6)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 53.1)
        self.assertEqual(round(st.Liquido.k.WmK, 4), 0.0810)
        self.assertEqual(round(st.Gas.rho, 2), 62.37)
        self.assertEqual(round(st.Gas.rhoM, 3), 2.226)
        self.assertEqual(round(st.Gas.uM.Jmol, 0), 1759)
        self.assertEqual(round(st.Gas.hM.Jmol, 0), 2418)
        self.assertEqual(round(st.Gas.sM.JmolK, 1), 135.2)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 25.12)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 55.27)
        self.assertEqual(round(st.Gas.w, 1), 181.4)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.39)
        self.assertEqual(round(st.Gas.k.WmK, 4), 0.0158)

        # Single phase region
        st = N2(T=65, P=1e4, **kw)
        self.assertEqual(round(st.rho, 4), 0.5219)
        self.assertEqual(round(st.rhoM, 5), 0.01863)
        self.assertEqual(round(st.uM.Jmol, 0), 1337)
        self.assertEqual(round(st.hM.Jmol, 0), 1874)
        self.assertEqual(round(st.sM.JmolK, 1), 166.3)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.94)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.48)
        self.assertEqual(round(st.w, 1), 163.7)
        self.assertEqual(round(st.mu.muPas, 2), 4.35)
        self.assertEqual(round(st.k.WmK, 5), 0.00592)

        st = N2(T=240, P=2e4, **kw)
        self.assertEqual(round(st.rho, 4), 0.2808)
        self.assertEqual(round(st.rhoM, 5), 0.01002)
        self.assertEqual(round(st.uM.Jmol, 0), 4980)
        self.assertEqual(round(st.hM.Jmol, 0), 6975)
        self.assertEqual(round(st.sM.JmolK, 1), 198.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.78)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.11)
        self.assertEqual(round(st.w, 1), 315.8)
        self.assertEqual(round(st.mu.muPas, 1), 15.0)
        self.assertEqual(round(st.k.WmK, 4), 0.0214)

        st = N2(T=75, P=5e4, **kw)
        self.assertEqual(round(st.rho, 3), 2.298)
        self.assertEqual(round(st.rhoM, 5), 0.08205)
        self.assertEqual(round(st.uM.Jmol, 0), 1525)
        self.assertEqual(round(st.hM.Jmol, 0), 2134)
        self.assertEqual(round(st.sM.JmolK, 1), 156.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.22)
        self.assertEqual(round(st.cpM.JmolK, 2), 30.31)
        self.assertEqual(round(st.w, 1), 174.2)
        self.assertEqual(round(st.mu.muPas, 2), 5.09)
        self.assertEqual(round(st.k.WmK, 5), 0.00719)

        st = N2(T=1900, P=1e5, **kw)
        self.assertEqual(round(st.rho, 4), 0.1773)
        self.assertEqual(round(st.rhoM, 6), 0.006329)
        self.assertEqual(round(st.uM.Jmol, -1), 45300)
        self.assertEqual(round(st.hM.Jmol, -1), 61100)
        self.assertEqual(round(st.sM.JmolK, 1), 250.1)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.35)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.67)
        self.assertEqual(round(st.w, 1), 857.7)
        self.assertEqual(round(st.mu.muPas, 1), 62.9)
        self.assertEqual(round(st.k.WmK, 3), 0.111)

        st = N2(T=300, P=101325, **kw)
        self.assertEqual(round(st.rho, 3), 1.138)
        self.assertEqual(round(st.rhoM, 5), 0.04063)
        self.assertEqual(round(st.uM.Jmol, 0), 6222)
        self.assertEqual(round(st.hM.Jmol, 0), 8716)
        self.assertEqual(round(st.sM.JmolK, 1), 191.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.80)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.15)
        self.assertEqual(round(st.w, 1), 353.2)
        self.assertEqual(round(st.mu.muPas, 1), 18.0)
        self.assertEqual(round(st.k.WmK, 4), 0.0258)

        st = N2(T=70, P=2e5, **kw)
        self.assertEqual(round(st.rho, 1), 841.0)
        self.assertEqual(round(st.rhoM, 2), 30.02)
        self.assertEqual(round(st.uM.Jmol, 0), -3829)
        self.assertEqual(round(st.hM.Jmol, 0), -3823)
        self.assertEqual(round(st.sM.JmolK, 2), 73.77)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.57)
        self.assertEqual(round(st.cpM.JmolK, 2), 57.14)
        self.assertEqual(round(st.w, 0), 1092)
        self.assertEqual(round(st.mu.muPas, 0), 204)
        self.assertEqual(round(st.k.WmK, 3), 0.144)

        st = N2(T=260, P=3e5, **kw)
        self.assertEqual(round(st.rho, 3), 3.895)
        self.assertEqual(round(st.rhoM, 4), 0.1390)
        self.assertEqual(round(st.uM.Jmol, 0), 5376)
        self.assertEqual(round(st.hM.Jmol, 0), 7533)
        self.assertEqual(round(st.sM.JmolK, 1), 178.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.81)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.29)
        self.assertEqual(round(st.w, 1), 328.9)
        self.assertEqual(round(st.mu.muPas, 1), 16.1)
        self.assertEqual(round(st.k.WmK, 4), 0.0231)

        st = N2(T=95, P=4e5, **kw)
        self.assertEqual(round(st.rho, 2), 15.73)
        self.assertEqual(round(st.rhoM, 4), 0.5615)
        self.assertEqual(round(st.uM.Jmol, 0), 1819)
        self.assertEqual(round(st.hM.Jmol, 0), 2531)
        self.assertEqual(round(st.sM.JmolK, 1), 145.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 22.30)
        self.assertEqual(round(st.cpM.JmolK, 2), 34.64)
        self.assertEqual(round(st.w, 1), 187.7)
        self.assertEqual(round(st.mu.muPas, 2), 6.64)
        self.assertEqual(round(st.k.WmK, 4), 0.0101)

        st = N2(T=900, P=6e5, **kw)
        self.assertEqual(round(st.rho, 3), 2.241)
        self.assertEqual(round(st.rhoM, 5), 0.08001)
        self.assertEqual(round(st.uM.Jmol, -1), 19380)
        self.assertEqual(round(st.hM.Jmol, -1), 26880)
        self.assertEqual(round(st.sM.JmolK, 1), 209.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 23.72)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.05)
        self.assertEqual(round(st.w, 1), 602.1)
        self.assertEqual(round(st.mu.muPas, 1), 38.8)
        self.assertEqual(round(st.k.WmK, 4), 0.0622)

        st = N2(T=350, P=8e5, **kw)
        self.assertEqual(round(st.rho, 3), 7.693)
        self.assertEqual(round(st.rhoM, 4), 0.2746)
        self.assertEqual(round(st.uM.Jmol, 0), 7231)
        self.assertEqual(round(st.hM.Jmol, -1), 10140)
        self.assertEqual(round(st.sM.JmolK, 1), 178.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.87)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.40)
        self.assertEqual(round(st.w, 1), 383.0)
        self.assertEqual(round(st.mu.muPas, 1), 20.3)
        self.assertEqual(round(st.k.WmK, 4), 0.0295)

        st = N2(T=100, P=1e6, **kw)
        self.assertEqual(round(st.rho, 1), 690.1)
        self.assertEqual(round(st.rhoM, 2), 24.63)
        self.assertEqual(round(st.uM.Jmol, 0), -2073)
        self.assertEqual(round(st.hM.Jmol, 0), -2033)
        self.assertEqual(round(st.sM.JmolK, 2), 94.65)
        self.assertEqual(round(st.cvM.JmolK, 2), 26.23)
        self.assertEqual(round(st.cpM.JmolK, 2), 64.36)
        self.assertEqual(round(st.w, 1), 616.8)
        self.assertEqual(round(st.mu.muPas, 1), 73.1)
        self.assertEqual(round(st.k.WmK, 3), 0.099)

        st = N2(T=120, P=2e6, **kw)
        self.assertEqual(round(st.rho, 2), 79.35)
        self.assertEqual(round(st.rhoM, 3), 2.833)
        self.assertEqual(round(st.uM.Jmol, 0), 1875)
        self.assertEqual(round(st.hM.Jmol, 0), 2581)
        self.assertEqual(round(st.sM.JmolK, 1), 134.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.98)
        self.assertEqual(round(st.cpM.JmolK, 2), 57.29)
        self.assertEqual(round(st.w, 1), 190.2)
        self.assertEqual(round(st.mu.muPas, 2), 9.39)
        self.assertEqual(round(st.k.WmK, 4), 0.0190)

        st = N2(T=1900, P=3e6, **kw)
        self.assertEqual(round(st.rho, 3), 5.289)
        self.assertEqual(round(st.rhoM, 4), 0.1888)
        self.assertEqual(round(st.uM.Jmol, -1), 45300)
        self.assertEqual(round(st.hM.Jmol, -1), 61190)
        self.assertEqual(round(st.sM.JmolK, 1), 221.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.36)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.68)
        self.assertEqual(round(st.w, 1), 862.6)
        self.assertEqual(round(st.mu.muPas, 1), 63.0)
        self.assertEqual(round(st.k.WmK, 3), 0.112)

        st = N2(T=127, P=4e6, **kw)
        self.assertEqual(round(st.rho, 1), 457.4)
        self.assertEqual(round(st.rhoM, 2), 16.33)
        self.assertEqual(round(st.uM.Jmol, 1), -109.4)
        self.assertEqual(round(st.hM.Jmol, 1), 135.5)
        self.assertEqual(round(st.sM.JmolK, 1), 112.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.52)
        self.assertEqual(round(st.cpM.JmolK, 1), 176.9)
        self.assertEqual(round(st.w, 1), 273.7)
        self.assertEqual(round(st.mu.muPas, 1), 28.8)
        self.assertEqual(round(st.k.WmK, 4), 0.0565)

        st = N2(T=800, P=5e6, **kw)
        self.assertEqual(round(st.rho, 2), 20.65)
        self.assertEqual(round(st.rhoM, 4), 0.7371)
        self.assertEqual(round(st.uM.Jmol, -1), 16980)
        self.assertEqual(round(st.hM.Jmol, -1), 23760)
        self.assertEqual(round(st.sM.JmolK, 1), 188.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 23.13)
        self.assertEqual(round(st.cpM.JmolK, 2), 31.60)
        self.assertEqual(round(st.w, 1), 580.9)
        self.assertEqual(round(st.mu.muPas, 1), 36.1)
        self.assertEqual(round(st.k.WmK, 4), 0.0575)

        st = N2(T=1900, P=6e6, **kw)
        self.assertEqual(round(st.rho, 2), 10.52)
        self.assertEqual(round(st.rhoM, 4), 0.3754)
        self.assertEqual(round(st.uM.Jmol, -1), 45300)
        self.assertEqual(round(st.hM.Jmol, -1), 61280)
        self.assertEqual(round(st.sM.JmolK, 1), 216.1)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.38)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.69)
        self.assertEqual(round(st.w, 1), 867.6)
        self.assertEqual(round(st.mu.muPas, 1), 63.0)
        self.assertEqual(round(st.k.WmK, 3), 0.112)

        st = N2(T=166, P=8e6, **kw)
        self.assertEqual(round(st.rho, 1), 237.2)
        self.assertEqual(round(st.rhoM, 3), 8.466)
        self.assertEqual(round(st.uM.Jmol, 0), 2040)
        self.assertEqual(round(st.hM.Jmol, 0), 2985)
        self.assertEqual(round(st.sM.JmolK, 1), 129.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.39)
        self.assertEqual(round(st.cpM.JmolK, 2), 64.41)
        self.assertEqual(round(st.w, 1), 265.1)
        self.assertEqual(round(st.mu.muPas, 1), 17.0)
        self.assertEqual(round(st.k.WmK, 4), 0.0318)

        st = N2(T=1900, P=1e7, **kw)
        self.assertEqual(round(st.rho, 2), 17.39)
        self.assertEqual(round(st.rhoM, 4), 0.6208)
        self.assertEqual(round(st.uM.Jmol, -1), 45300)
        self.assertEqual(round(st.hM.Jmol, -1), 61410)
        self.assertEqual(round(st.sM.JmolK, 1), 211.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.39)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.71)
        self.assertEqual(round(st.w, 1), 874.4)
        self.assertEqual(round(st.mu.muPas, 1), 63.1)
        self.assertEqual(round(st.k.WmK, 3), 0.113)

        st = N2(T=80, P=2e7, **kw)
        self.assertEqual(round(st.rho, 1), 837.7)
        self.assertEqual(round(st.rhoM, 2), 29.91)
        self.assertEqual(round(st.uM.Jmol, 0), -3500)
        self.assertEqual(round(st.hM.Jmol, 0), -2831)
        self.assertEqual(round(st.sM.JmolK, 2), 78.18)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.81)
        self.assertEqual(round(st.cpM.JmolK, 2), 52.70)
        self.assertEqual(round(st.w, 1), 988.1)
        self.assertEqual(round(st.mu.muPas, 0), 175)
        self.assertEqual(round(st.k.WmK, 3), 0.146)

        st = N2(T=370, P=4e7, **kw)
        self.assertEqual(round(st.rho, 1), 293.2)
        self.assertEqual(round(st.rhoM, 2), 10.47)
        self.assertEqual(round(st.uM.Jmol, 0), 6340)
        self.assertEqual(round(st.hM.Jmol, -1), 10160)
        self.assertEqual(round(st.sM.JmolK, 1), 145.1)
        self.assertEqual(round(st.cvM.JmolK, 2), 22.07)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.41)
        self.assertEqual(round(st.w, 1), 551.9)
        self.assertEqual(round(st.mu.muPas, 1), 29.7)
        self.assertEqual(round(st.k.WmK, 4), 0.0492)

        st = N2(T=1000, P=6e7, **kw)
        self.assertEqual(round(st.rho, 1), 166.3)
        self.assertEqual(round(st.rhoM, 3), 5.937)
        self.assertEqual(round(st.uM.Jmol, -1), 21370)
        self.assertEqual(round(st.hM.Jmol, -1), 31470)
        self.assertEqual(round(st.sM.JmolK, 1), 174.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.79)
        self.assertEqual(round(st.cpM.JmolK, 2), 33.76)
        self.assertEqual(round(st.w, 1), 774.3)
        self.assertEqual(round(st.mu.muPas, 1), 45.0)
        self.assertEqual(round(st.k.WmK, 3), 0.078)

        st = N2(T=500, P=8e7, **kw)
        self.assertEqual(round(st.rho, 1), 360.3)
        self.assertEqual(round(st.rhoM, 2), 12.86)
        self.assertEqual(round(st.uM.Jmol, 0), 8966)
        self.assertEqual(round(st.hM.Jmol, -1), 15180)
        self.assertEqual(round(st.sM.JmolK, 1), 148.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 22.63)
        self.assertEqual(round(st.cpM.JmolK, 2), 34.07)
        self.assertEqual(round(st.w, 1), 732.3)
        self.assertEqual(round(st.mu.muPas, 1), 38.5)
        self.assertEqual(round(st.k.WmK, 4), 0.0641)

        st = N2(T=100, P=1e8, **kw)
        self.assertEqual(round(st.rho, 1), 889.5)
        self.assertEqual(round(st.rhoM, 2), 31.75)
        self.assertEqual(round(st.uM.Jmol, 0), -3097)
        self.assertEqual(round(st.hM.Jmol, 1), 52.5)
        self.assertEqual(round(st.sM.JmolK, 2), 81.46)
        self.assertEqual(round(st.cvM.JmolK, 2), 30.39)
        self.assertEqual(round(st.cpM.JmolK, 2), 44.04)
        self.assertEqual(round(st.w, 0), 1221)
        self.assertEqual(round(st.mu.muPas, 0), 197)
        self.assertEqual(round(st.k.WmK, 3), 0.178)

        st = N2(T=1400, P=2e8, **kw)
        self.assertEqual(round(st.rho, 1), 315.5)
        self.assertEqual(round(st.rhoM, 2), 11.26)
        self.assertEqual(round(st.uM.Jmol, -1), 31490)
        self.assertEqual(round(st.hM.Jmol, -1), 49240)
        self.assertEqual(round(st.sM.JmolK, 1), 175.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 26.91)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.68)
        self.assertEqual(round(st.w, 0), 1119)
        self.assertEqual(round(st.mu.muPas, 1), 61.1)
        self.assertEqual(round(st.k.WmK, 3), 0.113)

        st = N2(T=140, P=4e8, **kw)
        self.assertEqual(round(st.rho, 0), 1012)
        self.assertEqual(round(st.rhoM, 2), 36.11)
        self.assertEqual(round(st.uM.Jmol, 0), -1996)
        self.assertEqual(round(st.hM.Jmol, 0), 9080)
        self.assertEqual(round(st.sM.JmolK, 2), 83.50)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.11)
        self.assertEqual(round(st.cpM.JmolK, 2), 37.92)
        self.assertEqual(round(st.w, 0), 1794)
        self.assertEqual(round(st.mu.muPas, 0), 291)
        self.assertEqual(round(st.k.WmK, 3), 0.272)

        st = N2(T=1800, P=1e9, **kw)
        self.assertEqual(round(st.rho, 1), 675.8)
        self.assertEqual(round(st.rhoM, 2), 24.13)
        self.assertEqual(round(st.uM.Jmol, -1), 43180)
        self.assertEqual(round(st.hM.Jmol, -1), 84630)
        self.assertEqual(round(st.sM.JmolK, 1), 170.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.79)
        self.assertEqual(round(st.cpM.JmolK, 2), 37.80)
        self.assertEqual(round(st.w, 0), 2053)
        self.assertEqual(round(st.mu.muPas, 0), 111)
        self.assertEqual(round(st.k.WmK, 3), 0.210)

    def test_jacobsen(self):
        # Selected point from Table 21, Pag 795, saturation states
        st = N2(T=63.15, x=0.5, eq="jacobsen")
        self.assertEqual(round(st.P.MPa, 5), 0.01252)
        self.assertEqual(round(st.Liquido.rhoM, 3), 31.046)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -4227.5)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 67.89)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 31.29)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 56.56)
        self.assertEqual(round(st.Liquido.w, 0), 1022)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.02410)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 1806.3)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 163.43)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 23.94)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 33.26)
        self.assertEqual(round(st.Gas.w, 0), 159)

        st = N2(T=80, x=0.5, eq="jacobsen")
        self.assertEqual(round(st.P.MPa, 5), 0.13698)
        self.assertEqual(round(st.Liquido.rhoM, 3), 28.351)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -3269.0)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 81.28)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 29.63)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 57.65)
        self.assertEqual(round(st.Liquido.w, 0), 821)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.21799)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2202.6)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 149.67)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 26.52)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 38.33)
        self.assertEqual(round(st.Gas.w, 0), 174)

        st = N2(T=100, x=0.5, eq="jacobsen")
        self.assertEqual(round(st.P.MPa, 5), 0.77917)
        self.assertEqual(round(st.Liquido.rhoM, 3), 24.584)
        self.assertEqual(round(st.Liquido.hM.Jmol, 1), -2050.6)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 94.57)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 27.84)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 2), 65.09)
        self.assertEqual(round(st.Liquido.w, 0), 601)
        self.assertEqual(round(st.Gas.rhoM, 4), 1.1442)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2450.6)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 139.58)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 27.43)
        self.assertEqual(round(st.Gas.cpM.JmolK, 2), 47.48)
        self.assertEqual(round(st.Gas.w, 0), 181)

        st = N2(T=120, x=0.5, eq="jacobsen")
        self.assertEqual(round(st.P.MPa, 4), 2.5130)
        self.assertEqual(round(st.Liquido.rhoM, 3), 18.644)
        self.assertEqual(round(st.Liquido.hM.Jmol, 2), -493.36)
        self.assertEqual(round(st.Liquido.sM.JmolK, 2), 107.95)
        self.assertEqual(round(st.Liquido.cvM.JmolK, 2), 29.06)
        self.assertEqual(round(st.Liquido.cpM.JmolK, 1), 128.9)
        self.assertEqual(round(st.Liquido.w, 0), 309)
        self.assertEqual(round(st.Gas.rhoM, 4), 4.4655)
        self.assertEqual(round(st.Gas.hM.Jmol, 1), 2081.3)
        self.assertEqual(round(st.Gas.sM.JmolK, 2), 129.40)
        self.assertEqual(round(st.Gas.cvM.JmolK, 2), 31.79)
        self.assertEqual(round(st.Gas.cpM.JmolK, 1), 131.4)
        self.assertEqual(round(st.Gas.w, 0), 171)

        # Selected values from Table 22, pag 800, single phase region
        st = N2(T=64, P=2e4, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 30.915)
        self.assertEqual(round(st.uM.Jmol, 1), -4179.9)
        self.assertEqual(round(st.hM.Jmol, 1), -4179.3)
        self.assertEqual(round(st.sM.JmolK, 1), 68.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.22)
        self.assertEqual(round(st.cpM.JmolK, 2), 56.49)
        self.assertEqual(round(st.w, 0), 1010)

        st = N2(T=72, P=4e4, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.06824)
        self.assertEqual(round(st.uM.Jmol, 1), 1459.9)
        self.assertEqual(round(st.hM.Jmol, 1), 2046.1)
        self.assertEqual(round(st.sM.JmolK, 2), 157.48)
        self.assertEqual(round(st.cvM.JmolK, 2), 23.32)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.87)
        self.assertEqual(round(st.w, 0), 170)

        st = N2(T=2000, P=6e4, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.00361)
        self.assertEqual(round(st.uM.Jmol, 0), 48178)
        self.assertEqual(round(st.hM.Jmol, 0), 64809)
        self.assertEqual(round(st.sM.JmolK, 2), 256.32)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.66)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.98)
        self.assertEqual(round(st.w, 0), 879)

        st = N2(T=100, P=8e4, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.09773)
        self.assertEqual(round(st.uM.Jmol, 1), 2046.2)
        self.assertEqual(round(st.hM.Jmol, 1), 2864.8)
        self.assertEqual(round(st.sM.JmolK, 2), 161.43)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.92)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.71)
        self.assertEqual(round(st.w, 0), 202)

        st = N2(T=76, P=1e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 29.017)
        self.assertEqual(round(st.uM.Jmol, 1), -3502.4)
        self.assertEqual(round(st.hM.Jmol, 1), -3498.9)
        self.assertEqual(round(st.sM.JmolK, 2), 78.34)
        self.assertEqual(round(st.cvM.JmolK, 2), 30.04)
        self.assertEqual(round(st.cpM.JmolK, 2), 57.02)
        self.assertEqual(round(st.w, 0), 865)

        st = N2(T=300, P=101325, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.04063)
        self.assertEqual(round(st.uM.Jmol, 1), 6222.7)
        self.assertEqual(round(st.hM.Jmol, 1), 8716.5)
        self.assertEqual(round(st.sM.JmolK, 2), 191.66)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.82)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.17)
        self.assertEqual(round(st.w, 0), 353)

        st = N2(T=84, P=2e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.30798)
        self.assertEqual(round(st.uM.Jmol, 1), 1635.9)
        self.assertEqual(round(st.hM.Jmol, 1), 2285.3)
        self.assertEqual(round(st.sM.JmolK, 2), 147.73)
        self.assertEqual(round(st.cvM.JmolK, 2), 26.17)
        self.assertEqual(round(st.cpM.JmolK, 2), 38.54)
        self.assertEqual(round(st.w, 0), 177)

        st = N2(T=120, P=3e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.31156)
        self.assertEqual(round(st.uM.Jmol, 1), 2418.1)
        self.assertEqual(round(st.hM.Jmol, 1), 3380.9)
        self.assertEqual(round(st.sM.JmolK, 2), 155.42)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.10)
        self.assertEqual(round(st.cpM.JmolK, 2), 30.62)
        self.assertEqual(round(st.w, 0), 219)

        st = N2(T=90, P=4e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 26.587)
        self.assertEqual(round(st.uM.Jmol, 1), -2692.1)
        self.assertEqual(round(st.hM.Jmol, 1), -2677.0)
        self.assertEqual(round(st.sM.JmolK, 2), 88.13)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.64)
        self.assertEqual(round(st.cpM.JmolK, 2), 60.16)
        self.assertEqual(round(st.w, 0), 714)

        st = N2(T=64, P=5e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 30.945)
        self.assertEqual(round(st.uM.Jmol, 1), -4184.9)
        self.assertEqual(round(st.hM.Jmol, 1), -4168.7)
        self.assertEqual(round(st.sM.JmolK, 2), 68.56)
        self.assertEqual(round(st.cvM.JmolK, 2), 31.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 56.43)
        self.assertEqual(round(st.w, 0), 1014)

        st = N2(T=1100, P=6e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.06548)
        self.assertEqual(round(st.uM.Jmol, 0), 24279)
        self.assertEqual(round(st.hM.Jmol, 0), 33442)
        self.assertEqual(round(st.sM.JmolK, 2), 216.41)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.93)
        self.assertEqual(round(st.cpM.JmolK, 2), 33.26)
        self.assertEqual(round(st.w, 0), 661)

        st = N2(T=102, P=8e5, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 1.1399)
        self.assertEqual(round(st.uM.Jmol, 1), 1823.9)
        self.assertEqual(round(st.hM.Jmol, 1), 2525.7)
        self.assertEqual(round(st.sM.JmolK, 2), 140.15)
        self.assertEqual(round(st.cvM.JmolK, 2), 25.66)
        self.assertEqual(round(st.cpM.JmolK, 2), 44.07)
        self.assertEqual(round(st.w, 0), 185)

        st = N2(T=2000, P=1e6, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 5), 0.06002)
        self.assertEqual(round(st.uM.Jmol, 0), 48178)
        self.assertEqual(round(st.hM.Jmol, 0), 64838)
        self.assertEqual(round(st.sM.JmolK, 2), 232.93)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.67)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.98)
        self.assertEqual(round(st.w, 0), 880)

        st = N2(T=116, P=2e6, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 3.1937)
        self.assertEqual(round(st.uM.Jmol, 1), 1684.5)
        self.assertEqual(round(st.hM.Jmol, 1), 2310.7)
        self.assertEqual(round(st.sM.JmolK, 2), 132.50)
        self.assertEqual(round(st.cvM.JmolK, 2), 29.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 80.10)
        self.assertEqual(round(st.w, 0), 177)

        st = N2(T=440, P=4e6, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 1.0773)
        self.assertEqual(round(st.uM.Jmol, 1), 9020.1)
        self.assertEqual(round(st.hM.Jmol, 0), 12733)
        self.assertEqual(round(st.sM.JmolK, 2), 172.00)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.16)
        self.assertEqual(round(st.cpM.JmolK, 2), 30.07)
        self.assertEqual(round(st.w, 0), 438)

        st = N2(T=300, P=8e6, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 3.2057)
        self.assertEqual(round(st.uM.Jmol, 1), 5777.2)
        self.assertEqual(round(st.hM.Jmol, 1), 8272.8)
        self.assertEqual(round(st.sM.JmolK, 2), 153.92)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.32)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.67)
        self.assertEqual(round(st.w, 0), 372)

        st = N2(T=950, P=1e7, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 1.2218)
        self.assertEqual(round(st.uM.Jmol, 0), 20508)
        self.assertEqual(round(st.hM.Jmol, 0), 28693)
        self.assertEqual(round(st.sM.JmolK, 2), 188.11)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.18)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.69)
        self.assertEqual(round(st.w, 0), 640)

        st = N2(T=112, P=2e7, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 25.665)
        self.assertEqual(round(st.uM.Jmol, 1), -1926.0)
        self.assertEqual(round(st.hM.Jmol, 1), -1146.7)
        self.assertEqual(round(st.sM.JmolK, 2), 95.89)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.59)
        self.assertEqual(round(st.cpM.JmolK, 2), 54.01)
        self.assertEqual(round(st.w, 0), 742)

        st = N2(T=370, P=4e7, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 10.462)
        self.assertEqual(round(st.uM.Jmol, 1), 6343.9)
        self.assertEqual(round(st.hM.Jmol, 0), 10167)
        self.assertEqual(round(st.sM.JmolK, 2), 145.07)
        self.assertEqual(round(st.cvM.JmolK, 2), 22.07)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.41)
        self.assertEqual(round(st.w, 0), 552)

        st = N2(T=740, P=6e7, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 7.6595)
        self.assertEqual(round(st.uM.Jmol, 0), 14972)
        self.assertEqual(round(st.hM.Jmol, 0), 22806)
        self.assertEqual(round(st.sM.JmolK, 2), 164.52)
        self.assertEqual(round(st.cvM.JmolK, 2), 23.39)
        self.assertEqual(round(st.cpM.JmolK, 2), 33.07)
        self.assertEqual(round(st.w, 0), 714)

        st = N2(T=1200, P=8e7, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 4), 6.4346)
        self.assertEqual(round(st.uM.Jmol, 0), 26438)
        self.assertEqual(round(st.hM.Jmol, 0), 38871)
        self.assertEqual(round(st.sM.JmolK, 2), 178.36)
        self.assertEqual(round(st.cvM.JmolK, 2), 25.92)
        self.assertEqual(round(st.cpM.JmolK, 2), 34.67)
        self.assertEqual(round(st.w, 0), 859)

        st = N2(T=360, P=2e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 24.003)
        self.assertEqual(round(st.uM.Jmol, 1), 4656.8)
        self.assertEqual(round(st.hM.Jmol, 0), 12989)
        self.assertEqual(round(st.sM.JmolK, 1), 128.60)
        self.assertEqual(round(st.cvM.JmolK, 2), 24.96)
        self.assertEqual(round(st.cpM.JmolK, 2), 36.48)
        self.assertEqual(round(st.w, 0), 1126)

        st = N2(T=260, P=4e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 32.343)
        self.assertEqual(round(st.uM.Jmol, 1), 1647.9)
        self.assertEqual(round(st.hM.Jmol, 0), 14016)
        self.assertEqual(round(st.sM.JmolK, 2), 109.04)
        self.assertEqual(round(st.cvM.JmolK, 2), 29.22)
        self.assertEqual(round(st.cpM.JmolK, 2), 40.15)
        self.assertEqual(round(st.w, 0), 1599)

        st = N2(T=140, P=5e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 37.663)
        self.assertEqual(round(st.uM.Jmol, 1), -1921.2)
        self.assertEqual(round(st.hM.Jmol, 0), 11355)
        self.assertEqual(round(st.sM.JmolK, 2), 80.59)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.39)
        self.assertEqual(round(st.cpM.JmolK, 2), 42.49)
        self.assertEqual(round(st.w, 0), 1901)

        st = N2(T=1200, P=6e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 23.051)
        self.assertEqual(round(st.uM.Jmol, 0), 26812)
        self.assertEqual(round(st.hM.Jmol, 0), 52841)
        self.assertEqual(round(st.sM.JmolK, 2), 160.83)
        self.assertEqual(round(st.cvM.JmolK, 2), 28.34)
        self.assertEqual(round(st.cpM.JmolK, 2), 36.02)
        self.assertEqual(round(st.w, 0), 1595)

        st = N2(T=175, P=8e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 39.839)
        self.assertEqual(round(st.uM.Jmol, 2), -375.97)
        self.assertEqual(round(st.hM.Jmol, 0), 19705)
        self.assertEqual(round(st.sM.JmolK, 2), 84.51)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.40)
        self.assertEqual(round(st.cpM.JmolK, 2), 41.01)
        self.assertEqual(round(st.w, 0), 2174)

        st = N2(T=1000, P=9e8, eq="jacobsen")
        self.assertEqual(round(st.rhoM, 3), 29.107)
        self.assertEqual(round(st.uM.Jmol, 0), 22222)
        self.assertEqual(round(st.hM.Jmol, 0), 53143)
        self.assertEqual(round(st.sM.JmolK, 2), 150.65)
        self.assertEqual(round(st.cvM.JmolK, 2), 29.15)
        self.assertEqual(round(st.cpM.JmolK, 2), 36.17)
        self.assertEqual(round(st.w, 0), 1899)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = N2(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.0979)
        self.assertEqual(round(st.P.MPa, 3), 51.268)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1719)

        st2 = N2(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 41.82)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.31053)

    def test_LemmonTransport(self):
        # Table V, pag 28
        # Viscosity
        self.assertEqual(round(N2(rhom=0, T=100).mu.muPas, 5), 6.90349)
        self.assertEqual(round(N2(rhom=0, T=300).mu.muPas, 4), 17.8771)
        self.assertEqual(round(N2(rhom=25, T=100).mu.muPas, 4), 79.7418)
        self.assertEqual(round(N2(rhom=10, T=200).mu.muPas, 4), 21.0810)
        self.assertEqual(round(N2(rhom=5, T=300).mu.muPas, 4), 20.7430)
        self.assertEqual(round(N2(rhom=11.18, T=126.195).mu.muPas, 4), 18.2978)

        # Thermal Conductivity
        self.assertEqual(round(N2(rhom=0, T=100).k.mWmK, 5), 9.27749)
        self.assertEqual(round(N2(rhom=0, T=300).k.mWmK, 4), 25.9361)
        self.assertEqual(round(N2(rhom=25, T=100).k.mWmK, 3), 103.834)
        self.assertEqual(round(N2(rhom=10, T=200).k.mWmK, 4), 36.0099)
        self.assertEqual(round(N2(rhom=5, T=300).k.mWmK, 4), 32.7694)
        self.assertEqual(round(N2(rhom=11.18, T=126.195).k.mWmK, 1), 675.8)

    def test_stephan(self):
        # The paper use a old Jacobsen EoS with MBWR format
        # Jacobsen, R.T., Stewart, R.B.
        # Thermodynamic Properties of Nitrogen Including Liquid and Vapor
        # Phases from 63K to 2000K with Pressures to 10000Bar
        # J. Phys. Chem. Ref. Data 2(4) (1973) 757-922
        # doi: 10.1063/1.3253132

        # So using TP as input parameter may differ, specially in region near
        # critical point

        kw = {"visco": 2, "thermal": 2}
        # Table A1, Pag 1013
        self.assertEqual(round(N2(T=80, P=1e5, **kw).mu.muPas, 2), 5.24)
        self.assertEqual(round(N2(T=300, P=1e6, **kw).mu.muPas, 2), 18.03)
        self.assertEqual(round(N2(T=1100, P=1e7, **kw).mu.muPas, 2), 44.67)
        self.assertEqual(round(N2(T=500, P=2e7, **kw).mu.muPas, 2), 28.37)
        self.assertEqual(round(N2(T=200, P=5e7, **kw).mu.muPas, 2), 49.40)
        self.assertEqual(round(N2(T=1100, P=1e8, **kw).mu.muPas, 2), 50.22)

        # Table A2, Pag 1016
        st = N2(T=100, x=0.5, **kw)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 7.08)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 72.88)

        # Table B1, Pag 1018
        self.assertEqual(round(N2(T=80, P=1e5, **kw).k.mWmK, 2), 7.73)
        self.assertEqual(round(N2(T=300, P=1e6, **kw).k.mWmK, 2), 26.51)
        self.assertEqual(round(N2(T=1100, P=1e7, **kw).k.mWmK, 2), 72.32)
        self.assertEqual(round(N2(T=500, P=2e7, **kw).k.mWmK, 2), 44.17)
        self.assertEqual(round(N2(T=200, P=5e7, **kw).k.mWmK, 2), 80.67)
        self.assertEqual(round(N2(T=1100, P=1e8, **kw).k.mWmK, 2), 83.72)

        # Table B2, Pag 1021
        self.assertEqual(round(st.Gas.k.mWmK, 2), 11.10)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 103.78)
