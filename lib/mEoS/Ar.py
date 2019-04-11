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


from math import exp, pi
from unittest import TestCase

from scipy.constants import Boltzmann

from lib import unidades
from lib.meos import MEoS


class Ar(MEoS):
    """Multiparamter equation of state for argon"""
    name = "argon"
    CASNumber = "7440-37-1"
    formula = "Ar"
    synonym = "R-740"
    _refPropName = "ARGON"
    _coolPropName = "Argon"
    rhoc = unidades.Density(535.6)
    Tc = unidades.Temperature(150.687)
    Pc = unidades.Pressure(4863, "kPa")
    M = 39.948  # g/mol
    Tt = unidades.Temperature(83.8058)
    Tb = unidades.Temperature(87.302)
    f_acent = -0.00219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 98
    _Tr = unidades.Temperature(147.707801)
    _rhor = unidades.Density(540.014968)
    _w = 0.000305675

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [8.31666243, -4.94651164],
           "ao_exp": [], "titao": []}

    CP1 = {"ao": 2.5,
           "an": [], "pow": [], "ao_exp": [], "exp": []}

    Fi2 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [8.3166315, -4.9465026],
           "ao_exp": [], "titao": []}

    tegeler = {
        "__type__": "Helmholtz",
        "__name__": "FEQ Helmholtz equation of state for argon of Tegeler et "
                    "al. (1999).",
        "__doi__": {"autor": "Tegeler, Ch., Span, R., Wagner, W.",
                    "title": "A New Equation of State for Argon Covering the "
                             "Fluid Region for Temperatures From the Melting "
                             "Line to 700 K at Pressures up to 1000 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 28, 779 (1999)",
                    "doi": "10.1063/1.556037"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 2000., "Pmax": 1000000.0, "rhomax": 50.65,
        "Pmin": 68.891, "rhomin": 35.465,

        "nr1": [0.887223049900e-1, 0.705148051673, -0.168201156541e1,
                -0.149090144315, -0.120248046009, -0.121649787986,
                0.400359336268, -0.271360626991, 0.242119245796,
                0.578895831856e-2, -0.410973356153e-1, 0.247107615416e-1],
        "d1": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4],
        "t1": [0., 0.25, 1., 2.75, 4.0, 0., 0.25, 0.75, 2.75, 0.0, 2.0, 0.75],

        "nr2": [-0.321813917507, 0.332300176958, 0.310199862873e-1,
                -0.307770860024e-1, 0.938911374196e-1, -0.906432106820e-1,
                -0.457783492767e-3, -0.826597290252e-4, 0.130134156031e-3,
                -0.113978400020e-1, -0.244551699605e-1, -0.643240671760e-1,
                0.588894710937e-1, -0.649335521130e-3, -0.138898621584e-1,
                0.404898392969, -0.386125195947, -0.188171423322,
                0.159776475965, 0.539855185139e-1, -0.289534179580e-1,
                -0.130254133814e-1, 0.289486967758e-2, -0.226471343048e-2,
                0.176164561964e-2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
               3, 3, 4, 4],
        "d2": [1, 1, 3, 4, 4, 5, 7, 10, 10, 2, 2, 4, 4, 8, 3, 5, 5, 6, 6, 7, 7,
               8, 9, 5, 6],
        "t2": [3., 3.5, 1., 2., 4., 3., 0., 0.5, 1., 1., 7., 5., 6., 6., 10.,
               13., 14., 11., 14., 8., 14., 6., 7., 24., 22.],
        "gamma2": [1]*25,

        "nr3": [0.585524544828e-2, -0.6925190827, 0.153154900305e1,
                -0.273804474498e-2],
        "d3": [2, 1, 2, 3],
        "t3": [3, 1, 0, 0],
        "alfa3": [20]*4,
        "beta3": [250, 375, 300, 225],
        "gamma3": [1.11, 1.14, 1.17, 1.11],
        "epsilon3": [1, 1, 1, 1],
        "nr4": []}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for argon of Kunz and Wagner "
                    "(2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": 83.8, "Tmax": 700., "Pmax": 1000000.0, "rhomax": 50.65,
        "Pmin": 68.891, "rhomin": 35.465,

        "nr1": [0.85095714803969, -0.24003222943480e1, 0.54127841476466,
                0.16919770692538e-1, 0.68825965019035e-1, 0.21428032815338e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.17429895321992, -0.033654495604194, -0.13526799857691,
                -0.016387350791552, -0.024987666851475, 0.0088769204815709],
        "c2": [1, 1, 2, 2, 3, 3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    stewart = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for argon of Stewart and "
                    "Jacobsen (1989).",
        "__doi__": {"autor": "Stewart, R.B., Jacobsen, R.T.",
                    "title": "Thermodynamic Properties of Argon from the "
                             "Triple Point to 1200 K at Pressures to 1000 MPa",
                    "ref": "J. Phys. Chem. Ref. Data, 18(2):639-798, 1989",
                    "doi": "10.1063/1.555829"},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 6227.9, "so": 154.84},
        "Tc": 150.6633, "Pc": 4860, "rhoc": 13.29, "Tt": 83.804,

        "Tmin": 83.804, "Tmax": 1200., "Pmax": 1000000.0, "rhomax": 45.814,
        "Pmin": 68.961, "rhomin": 35.475,

        "nr1": [0.7918675715, -1.633346151, -0.439530293, 0.1033899999,
                0.2061801664, -0.2888681776, 0.439801055, -0.08429550391,
                -0.2155658654, 0.4786509099, -0.3525884593, 0.03015073692,
                0.02987679059, -0.01522568583, 0.0007435785786],
        "d1": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 6],
        "t1": [0.25, 1, 3, 4, 0.25, 1, 2.5, 3.5, 0.75, 1, 1.5, 2.5, 1, 2, 2],

        "nr2": [0.07099541624, -0.02904237185, -0.06223078525, 0.0001410895187,
                -0.001481241783, 0.03023342784, -0.06126784685, 0.0270996709,
                0.09411034405, -0.007291645114, -0.001586314976,
                0.0009510948813, 0.0007786181844],
        "c2": [3, 3, 2, 4, 6, 3, 3, 3, 2, 2, 4, 2, 2],
        "d2": [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 8, 8],
        "t2": [5, 7, 5, 22, 16, 10, 14, 16, 4, 8, 10, 5, 6],
        "gamma2": [1]*13,

        "nr3": [],
        "nr4": []}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for argon of Span and "
                    "Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 750., "Pmax": 100000.0, "rhomax": 50.65,
        "Pmin": 69.026, "rhomin": 35.498,

        "nr1": [0.85095715, -0.24003223e1, 0.54127841, 0.16919771e-1,
                0.68825965e-1, 0.21428033e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.17429895, -0.33654496e-1, -0.135268, -0.16387351e-1,
                -0.24987667e-1, 0.88769205e-2],
        "c2": [1, 1, 2, 2, 3, 3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for argon of Younglove (1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, "
                             "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                             "Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 9715.8, "so": 156.65},
        "Tt": 83.8, "Tc": 150.86, "rhoc": 13.41,

        "Tmin": 83.80, "Tmax": 500., "Pmax": 101000.0, "rhomax": 50.65,
        "Pmin": 68.906, "rhomin": 35.4,

        "gamma": -0.0055542372,
        "b": [None, -0.6569731294e-3, 0.1822957801, -0.3649470141e1,
              0.1232012107e3, -0.8613578274e4, 0.7978579691e-4, -0.02911489110,
              0.7581821758e1, 0.8780488169e4, 0.1423145989e-6, 0.1674146131e-2,
              -0.3200447909, 0.2561766372e-4, -0.5475934941e-3, -0.4505032058,
              0.2013254653e-4, -0.1678941273e-6, 0.4207329271e-3,
              -0.5444212996e-5, -0.8004855011e4, -0.1319304201e6, -49.54923930,
              0.8092132177e5, -0.9870104061e-1, 2.020441562, -0.1637417205e-3,
              -0.7038944136, -0.1154324539e-6, 0.1555990117e-4,
              -0.1492178536e-9, -0.1001356071e-7, 0.2933963216e-6]}

    eq = tegeler, younglove, GERG, stewart, shortSpan
    _PR = -0.0034

    _surface = {"sigma": [0.037], "exp": [1.25]}
    _dielectric = {
        "eq": 1,
        "a": [4.1414], "b": [1.597, 0.262], "c": [-117.9],
        "Au": 0, "D": 2.1}

    _melting = {
        "eq": 1,
        "__doi__": tegeler["__doi__"],
        "Tmin": Tt, "Tmax": 700.0,
        "Tref": Tt, "Pref": 68891,

        "a0": 1,
        "a2": [-7476.2665, 9959.0613], "exp2": [1.05, 1.275]}

    _sublimation = {
        "eq": 3,
        "__doi__": tegeler["__doi__"],
        "Tmin": Tt, "Tmax": Tt,
        "Tref": Tt, "Pref": 68891,

        "a1": [], "exp1": [],
        "a2": [-11.391604, -0.39513431], "exp2": [1, 2.7],
        "a3": [], "exp3": []}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-5.9409785, 1.3553888, -0.4649761, -1.5399043],
        "t": [1., 1.5, 2., 4.5]}
    _liquid_Density = {
        "eq": 2,
        "n": [1.5004264, -0.3138129, 0.086461622, -0.041477525],
        "t": [0.334, 2./3, 7./3, 4]}
    _vapor_Density = {
        "eq": 3,
        "n": [-0.29182e1, 0.97930e-1, -0.13721e1, -0.22898e1],
        "t": [0.72, 1.25, 0.32, 4.34]}

    visco0 = {"__name__": "Lemmon (2004)",
              "__doi__": {
                   "autor": "Lemmon, E.W., Jacobsen, R.T.",
                   "title": "Viscosity and Thermal Conductivity Equations for "
                            "Nitrogen, Oxygen, Argon, and Air",
                   "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                   "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

              "eq": 1, "omega": 1,
              "ek": 143.2, "sigma": 0.335,

              "nr": [12.19, 13.99, 0.005027, -18.93, -6.698, -3.827],
              "tr": [0.42, 0.0, 0.95, 0.5, 0.9, 0.8],
              "dr": [1, 2, 10, 5, 1, 2],
              "gr": [0, 0, 0, 1, 1, 1],
              "cr": [0, 0, 0, 2, 4, 4]}

    visco1 = {"__name__": "Younglove-Hanley (1986)",
              "__doi__": {
                  "autor": "Younglove, B.A., Hanley, H.J.M.",
                  "title": "The Viscosity and Thermal Conductivity "
                           "Coefficients of Gaseous and Liquid Argon",
                  "ref": "J. Phys. Chem. Ref. Data 15(4) (1986) 1323-1337",
                  "doi": "10.1063/1.555765"},

              "eq": 1, "omega": 0,

              "Toref": 1,
              "no": [-0.8973188257e5, 0.8259113473e5, -0.2766475915e5,
                     0.3068539784e4, 0.4553103615e3, -0.1793443839e3,
                     0.2272225106e2, -0.1350672796e1, 0.3183693230e-1],
              "to": [-1., -2/3, -1/3, 0, 1/3, 2/3, 1., 4/3, 5/3],

              "Tref_res": 1, "rhoref_res": 1*M,
              "nr_num": [0.5927733783, -0.4251221169e2, -0.2698477165e-1,
                         0.3727762288e2, -0.3958508720e4, 0.3636730841e-2,
                         -0.2633471347e1, 0.2936563322e3, -0.3811869019e-4,
                         0.4451947464e-1, -0.5385874487e1],
              "tr_num": [0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2],
              "dr_num": [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
              "nr_den": [1.0, -0.1115054926e-1, -0.1328893444e1],
              "tr_den": [0, 0, 1],
              "dr_den": [0, 1, 1]}

    visco2 = {"__name__": "Younglove (1982)",
              "__doi__": {
                  "autor": "Younglove, B.A.",
                  "title": "Thermophysical Properties of Fluids. I. Argon, "
                           "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                           "Trifluoride, and Oxygen",
                  "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                  "doi": ""},

              "eq": 2, "omega": 0,
              "mod": True,

              "no": [0.61145472787e4, -0.10394390312e5, 0.67594614619e4,
                     -0.22536509380e4, 0.42593950138e3, -0.47252671093e2,
                     0.31795275425e1, -0.11629083780, 0.18043010592e-2],
              "to": [-1, -2/3, -1/3, 0, 1/3, 2/3, 1, 4/3, 5/3],

              "F": [0.14653652433, -0.77487424965e-1, 0.14e1, 0.1528e3],
              "E": [-0.12313579086e2, 0.20694685712, 0.16029145122e2,
                    0.11717461351e4, -0.5699589878e3, 0.40136071933e2,
                    0.39870122403e2],
              "rhoc": 0.537*1000/39.948}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Lemmon (2004)",
               "__doi__": {
                   "autor": "Lemmon, E.W., Jacobsen, R.T.",
                   "title": "Viscosity and Thermal Conductivity Equations for "
                            "Nitrogen, Oxygen, Argon, and Air",
                   "ref": "Int. J. Thermophys., 25(1) (2004) 21-69",
                   "doi": "10.1023/B:IJOT.0000022327.04529.f3"},

               "eq": 1,

               "Toref": 150.687, "koref": 1e-3,
               "no_visco": 0.8158,
               "no": [-0.432],
               "to": [0.77],

               "Tref_res": 150.687, "rhoref_res": 13.40743*M, "kref_res": 1e-3,
               "nr": [13.73, 10.07, 0.7375, -33.96, 20.47, -2.274, -3.973],
               "tr": [0.0, 0.0, 0.0, 0.8, 1.2, 0.8, 0.5],
               "dr": [1, 2, 4, 5, 6, 9, 1],
               "cr": [0, 0, 0, 2, 2, 2, 4],
               "gr": [0, 0, 0, 1, 1, 1, 1],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.13e-9, "gam0": 0.055, "qd": 0.32e-9, "Tcref": 301.374}

    thermo1 = {"__name__": "Younglove (1982)",
               "__doi__": {
                   "autor": "Younglove, B.A.",
                   "title": "Thermophysical Properties of Fluids. I. Argon, "
                            "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                            "Trifluoride, and Oxygen",
                   "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                   "doi": ""},

               "eq": 2,

               "no": [0.62777703742e1, -0.96096376637e1, 0.58887549191e1,
                      -0.1892092632e1, 0.34886571437, -0.38016786193e-1,
                      0.25207283167e-2, -0.91098744478e-4, 0.13990842942e-5],
               "to": [-1., -2/3, -1/3, 0, 1/3, 2/3, 1., 4/3, 5/3],

               "F": [0.2414210327e-1, 0.75696234255e-2, 1, 0.1528e3],
               "E": [-0.33327027332e2, 0, 0.30694859971e2, 0, 0.22956551674e4,
                     -0.35559415848e3, 0, 1],

               "critical": 1,
               "rhoc": 0.533, "Tc": 150.725,
               "ek": 52.8, "f": 1.7124, "gm": 3.669e-10}

    thermo2 = {"__name__": "Perkins (1991)",
               "__doi__": {
                   "autor": "Perkins, R.A., Friend, D.G., Roder, H.M., Nieto "
                            "de Castro, C.A.",
                   "title": "Thermal Conductivity Surface of Argon: A Fresh "
                            "Analysis",
                   "ref": "Int. J. Thermophys., 12(6) (1991) 965-984",
                   "doi": "10.1007/BF00503513"},

               "eq": 1,
               "Pc": 4.860e6, "rhoc": 13.29*M,

               # Dilute-gas term from Younglove and Hanley correlation
               "Toref": 1, "koref": 1e-3,
               "no": [-0.6700976182e5, 0.6152255283e5, -0.2049218286e5,
                      0.2216966254e4, 0.3579189325e3, -0.1364658914e3,
                      0.1718671649e2, -0.1018933154e1, 0.2397996932e-1],
               "to": [-1, -2./3, -1./3, 0, 1./3, 2./3, 1., 4./3, 5./3],

               "rhoref_res": M, "kref_res": 1.,
               "nr": [0.757894e-3, 0.612624e-4, -0.205353e-5,  0.745621e-7],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.16e-9, "gam0": 0.075, "qd": 0.2e-9, "Tcref": 376.65825}

    thermo3 = {"__name__": "Younglove-Hanley (1986)",
               "__doi__": {
                   "autor": "Younglove, B.A., Hanley, H.J.M.",
                   "title": "The Viscosity and Thermal Conductivity "
                            "Coefficients of Gaseous and Liquid Argon",
                   "ref": "J. Phys. Chem. Ref. Data 15(4) (1986) 1323-1337",
                   "doi": "10.1063/1.555765"},

               "eq": 1,

               "Toref": 1, "koref": 1e-3,
               "no": [-0.6700976182e5, 0.6152255283e5, -0.2049218286e5,
                      0.2216966254e4, 0.3579189325e3, -0.1364658914e3,
                      0.1718671649e2, -0.1018933154e1, 0.2397996932e-1],
               "to": [-1, -2./3, -1./3, 0, 1./3, 2./3, 1., 4./3, 5./3],

               "Tref_res": 1, "rhoref_res": M, "kref_res": 1e-3,
               "nr_num": [0.1536300190e1, -0.2332533199e3, -0.3027085824e-1,
                          0.1896279196e2, 0.1054230664e2, 0.2588139028e-4,
                          -0.4546798772, 0.4320206998e1, 0.1593643304e-4,
                          0.1262253904e-3, -0.2937213042e-2],
               "tr_num": [0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2],
               "dr_num": [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
               "nr_den": [1, -0.2262773007e-1, -0.1445619495],
               "tr_den": [0, 0, 1],
               "dr_den": [0, 1, 1],

               "critical": "_ThCritical"}

    def _ThCritical(self, rho, T, fase):
        """Younglove-Hanley thermal conductivity critical enhancement"""
        rhom = rho/self.M
        Tc = 150.86
        Pc = 4.9058e6
        rhocm = 13.41
        rhoc = rhocm*self.M

        # Eq 9
        T_ = T/Tc
        rho_ = rhom/rhocm
        DelT = (T-Tc)/Tc
        Delrho = (rhom-rhocm)/rhocm
        Xt = rho*fase.drhodP_T*Pc/rhoc**2

        # Eq 10
        tc = 1.02*Boltzmann*Pc/6.0795e-1/6/pi/fase.mu*(T_/rho_)**2 * \
            fase.dpdT_rho**2*Xt**0.46807*exp(-39.8*DelT**2-5.45*Delrho**4)

        return tc

    _thermal = thermo0, thermo1, thermo2, thermo3


class Test(TestCase):

    def test_Tegeler(self):
        # Selected point from Table 33, Pag 828, saturation states
        st = Ar(T=83.8058, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.068891)
        self.assertEqual(round(st.Liquido.rho, 2), 1416.77)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -276.56)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -2.5440)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.54960)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.1157)
        self.assertEqual(round(st.Liquido.w, 2), 862.43)
        self.assertEqual(round(st.Gas.rho, 4), 4.0546)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -112.85)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.59044)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.32471)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.55503)
        self.assertEqual(round(st.Gas.w, 2), 168.12)

        st = Ar(T=90, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.13351)
        self.assertEqual(round(st.Liquido.rho, 2), 1378.63)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -269.61)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -2.4645)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.52677)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.1212)
        self.assertEqual(round(st.Liquido.w, 2), 819.45)
        self.assertEqual(round(st.Gas.rho, 4), 7.4362)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -110.55)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.69718)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.33094)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.57569)
        self.assertEqual(round(st.Gas.w, 2), 172.83)

        st = Ar(T=120, x=0.5)
        self.assertEqual(round(st.P.MPa, 3), 1.2130)
        self.assertEqual(round(st.Liquido.rho, 2), 1162.82)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -233.48)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -2.1274)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.45763)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.3324)
        self.assertEqual(round(st.Liquido.w, 2), 584.19)
        self.assertEqual(round(st.Gas.rho, 3), 60.144)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -106.71)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.0710)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.38934)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.86265)
        self.assertEqual(round(st.Gas.w, 2), 185.09)

        st = Ar(T=150, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 4.7346)
        self.assertEqual(round(st.Liquido.rho, 2), 680.43)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -173.01)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -1.7145)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.70603)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 23.582)
        self.assertEqual(round(st.Liquido.w, 2), 174.74)
        self.assertEqual(round(st.Gas.rho, 2), 394.50)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -143.60)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.5185)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.82182)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 35.468)
        self.assertEqual(round(st.Gas.w, 2), 157.01)

        # Selected points from Table 34, Pag 830, Single phase points

        st = Ar(T=83.814, P=1e5)
        self.assertEqual(round(st.rho, 2), 1416.80)
        self.assertEqual(round(st.u.kJkg, 2), -276.61)
        self.assertEqual(round(st.h.kJkg, 2), -276.54)
        self.assertEqual(round(st.s.kJkgK, 4), -2.5440)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.54961)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1156)
        self.assertEqual(round(st.w, 2), 862.52)

        st = Ar(T=700, P=1e5)
        self.assertEqual(round(st.rho, 5), 0.68619)
        self.assertEqual(round(st.u.kJkg, 3), 63.355)
        self.assertEqual(round(st.h.kJkg, 2), 209.09)
        self.assertEqual(round(st.s.kJkgK, 5), 0.44677)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.31223)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.5205)
        self.assertEqual(round(st.w, 2), 492.95)

        st = Ar(T=150, P=5e5)
        self.assertEqual(round(st.rho, 3), 16.605)
        self.assertEqual(round(st.u.kJkg, 2), -110.45)
        self.assertEqual(round(st.h.kJkg, 3), -80.334)
        self.assertEqual(round(st.s.kJkgK, 5), -0.70404)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.32098)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.55987)
        self.assertEqual(round(st.w, 2), 224.97)

        st = Ar(T=170, P=1e6)
        self.assertEqual(round(st.rho, 3), 29.723)
        self.assertEqual(round(st.u.kJkg, 2), -105.62)
        self.assertEqual(round(st.h.kJkg, 3), -71.972)
        self.assertEqual(round(st.s.kJkgK, 5), -0.78987)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.32356)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.57801)
        self.assertEqual(round(st.w, 2), 238.88)

        st = Ar(T=125, P=2e6)
        self.assertEqual(round(st.rho, 2), 1122.34)
        self.assertEqual(round(st.u.kJkg, 2), -228.45)
        self.assertEqual(round(st.h.kJkg, 2), -226.66)
        self.assertEqual(round(st.s.kJkgK, 4), -2.0773)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.45179)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.4048)
        self.assertEqual(round(st.w, 2), 544.65)

        st = Ar(T=135, P=3e6)
        self.assertEqual(round(st.rho, 2), 1020.52)
        self.assertEqual(round(st.u.kJkg, 2), -214.63)
        self.assertEqual(round(st.h.kJkg, 2), -211.69)
        self.assertEqual(round(st.s.kJkgK, 4), -1.9694)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.44845)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.7159)
        self.assertEqual(round(st.w, 2), 445.83)

        st = Ar(T=150, P=4e6)
        self.assertEqual(round(st.rho, 2), 209.45)
        self.assertEqual(round(st.u.kJkg, 2), -134.47)
        self.assertEqual(round(st.h.kJkg, 2), -115.38)
        self.assertEqual(round(st.s.kJkgK, 4), -1.3116)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.46106)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.8982)
        self.assertEqual(round(st.w, 2), 193.39)

        st = Ar(T=675, P=5e6)
        self.assertEqual(round(st.rho, 3), 35.123)
        self.assertEqual(round(st.u.kJkg, 3), 53.149)
        self.assertEqual(round(st.h.kJkg, 2), 195.50)
        self.assertEqual(round(st.s.kJkgK, 5), -0.38989)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.31374)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.52898)
        self.assertEqual(round(st.w, 2), 493.26)

        st = Ar(T=400, P=1e7)
        self.assertEqual(round(st.rho, 2), 119.43)
        self.assertEqual(round(st.u.kJkg, 3), -40.267)
        self.assertEqual(round(st.h.kJkg, 3), 43.464)
        self.assertEqual(round(st.s.kJkgK, 5), -0.82699)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.32049)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.57843)
        self.assertEqual(round(st.w, 2), 391.51)

        st = Ar(T=110, P=1e8)
        self.assertEqual(round(st.rho, 2), 1494.28)
        self.assertEqual(round(st.u.kJkg, 2), -268.41)
        self.assertEqual(round(st.h.kJkg, 2), -201.49)
        self.assertEqual(round(st.s.kJkgK, 4), -2.4764)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.56101)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.94752)
        self.assertEqual(round(st.w, 1), 1074.3)

        st = Ar(T=330, P=1e9)
        self.assertEqual(round(st.rho, 2), 1764.39)
        self.assertEqual(round(st.u.kJkg, 2), -144.17)
        self.assertEqual(round(st.h.kJkg, 2), 422.60)
        self.assertEqual(round(st.s.kJkgK, 4), -2.0896)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.55397)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.75591)
        self.assertEqual(round(st.w, 1), 1851.9)

    def test_Younglove(self):
        # Selected point from Appendix F, Pag 1-12, saturation states
        # The pressure and density in saturation is calculate in tables using
        # the ancillary equation used in paper so the calculated point differ
        # of implement eq
        st = Ar(T=400, P=8e4, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 4), 0.9610)
        self.assertEqual(round(st.rhoM, 5), 0.02406)
        self.assertEqual(round(st.uM.Jmol, 0), 8473)
        self.assertEqual(round(st.hM.Jmol, -1), 11800)
        self.assertEqual(round(st.sM.JmolK, 1), 164.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.48)
        self.assertEqual(round(st.cpM.JmolK, 2), 20.81)
        self.assertEqual(round(st.w, 1), 372.6)
        self.assertEqual(round(st.mu.muPas, 1), 28.9)
        self.assertEqual(round(st.k.WmK, 4), 0.0226)

        st = Ar(T=100, P=1e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 3), 4.915)
        self.assertEqual(round(st.rhoM, 4), 0.1230)
        self.assertEqual(round(st.uM.Jmol, 0), 4700)
        self.assertEqual(round(st.hM.Jmol, 0), 5512)
        self.assertEqual(round(st.sM.JmolK, 1), 133.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.87)
        self.assertEqual(round(st.cpM.JmolK, 2), 21.93)
        self.assertEqual(round(st.w, 1), 184.0)
        self.assertEqual(round(st.mu.muPas, 2), 8.24)
        self.assertEqual(round(st.k.WmK, 5), 0.00663)

        st = Ar(T=300, P=101325, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 3), 1.624)
        self.assertEqual(round(st.rhoM, 5), 0.04065)
        self.assertEqual(round(st.uM.Jmol, 0), 7223)
        self.assertEqual(round(st.hM.Jmol, 0), 9716)
        self.assertEqual(round(st.sM.JmolK, 1), 156.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.48)
        self.assertEqual(round(st.cpM.JmolK, 2), 20.83)
        self.assertEqual(round(st.w, 1), 322.7)
        self.assertEqual(round(st.mu.muPas, 1), 22.9)
        self.assertEqual(round(st.k.WmK, 4), 0.0179)

        st = Ar(T=90, P=2e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1377)
        self.assertEqual(round(st.rhoM, 2), 34.48)
        self.assertEqual(round(st.uM.Jmol, 0), -1100)
        self.assertEqual(round(st.hM.Jmol, 0), -1094)
        self.assertEqual(round(st.sM.JmolK, 2), 57.99)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.52)
        self.assertEqual(round(st.cpM.JmolK, 2), 44.89)
        self.assertEqual(round(st.w, 1), 805.0)
        self.assertEqual(round(st.mu.muPas, 0), 241)
        self.assertEqual(round(st.k.WmK, 3), 0.124)

        st = Ar(T=100, P=3e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 2), 15.52)
        self.assertEqual(round(st.rhoM, 4), 0.3884)
        self.assertEqual(round(st.uM.Jmol, 0), 4623)
        self.assertEqual(round(st.hM.Jmol, 0), 5395)
        self.assertEqual(round(st.sM.JmolK, 1), 123.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.77)
        self.assertEqual(round(st.cpM.JmolK, 2), 24.78)
        self.assertEqual(round(st.w, 1), 179.2)
        self.assertEqual(round(st.mu.muPas, 2), 8.34)
        self.assertEqual(round(st.k.WmK, 5), 0.00704)

        st = Ar(T=400, P=4e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 3), 4.805)
        self.assertEqual(round(st.rhoM, 4), 0.1203)
        self.assertEqual(round(st.uM.Jmol, 0), 8460)
        self.assertEqual(round(st.hM.Jmol, -1), 11790)
        self.assertEqual(round(st.sM.JmolK, 1), 151.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.49)
        self.assertEqual(round(st.cpM.JmolK, 2), 20.89)
        self.assertEqual(round(st.w, 1), 373.1)
        self.assertEqual(round(st.mu.muPas, 1), 29.0)
        self.assertEqual(round(st.k.WmK, 4), 0.0227)

        st = Ar(T=90, P=6e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1379)
        self.assertEqual(round(st.rhoM, 2), 34.51)
        self.assertEqual(round(st.uM.Jmol, 0), -1104)
        self.assertEqual(round(st.hM.Jmol, 0), -1087)
        self.assertEqual(round(st.sM.JmolK, 2), 57.93)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.54)
        self.assertEqual(round(st.cpM.JmolK, 2), 44.81)
        self.assertEqual(round(st.w, 1), 806.7)
        self.assertEqual(round(st.mu.muPas, 0), 242)
        self.assertEqual(round(st.k.WmK, 3), 0.124)

        st = Ar(T=130, P=8e5, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 2), 32.50)
        self.assertEqual(round(st.rhoM, 4), 0.8135)
        self.assertEqual(round(st.uM.Jmol, 0), 4926)
        self.assertEqual(round(st.hM.Jmol, 0), 5909)
        self.assertEqual(round(st.sM.JmolK, 1), 120.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.61)
        self.assertEqual(round(st.cpM.JmolK, 2), 25.41)
        self.assertEqual(round(st.w, 1), 203.7)
        self.assertEqual(round(st.mu.muPas, 1), 10.9)
        # self.assertEqual(round(st.k.WmK, 5), 0.00938)

        st = Ar(T=110, P=1e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1242)
        self.assertEqual(round(st.rhoM, 2), 31.10)
        self.assertEqual(round(st.uM.Jmol, 1), -183.7)
        self.assertEqual(round(st.hM.Jmol, 1), -151.6)
        self.assertEqual(round(st.sM.JmolK, 2), 67.18)
        self.assertEqual(round(st.cvM.JmolK, 2), 18.89)
        self.assertEqual(round(st.cpM.JmolK, 2), 48.89)
        self.assertEqual(round(st.w, 1), 679.0)
        self.assertEqual(round(st.mu.muPas, 0), 143)
        self.assertEqual(round(st.k.WmK, 4), 0.0961)

        st = Ar(T=280, P=1.5e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 2), 26.06)
        self.assertEqual(round(st.rhoM, 4), 0.6523)
        self.assertEqual(round(st.uM.Jmol, 0), 6878)
        self.assertEqual(round(st.hM.Jmol, 0), 9177)
        self.assertEqual(round(st.sM.JmolK, 1), 132.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.59)
        self.assertEqual(round(st.cpM.JmolK, 2), 21.66)
        self.assertEqual(round(st.w, 1), 312.7)
        self.assertEqual(round(st.mu.muPas, 1), 21.9)
        self.assertEqual(round(st.k.WmK, 4), 0.0175)

        st = Ar(T=128, P=2e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1089)
        self.assertEqual(round(st.rhoM, 2), 27.25)
        self.assertEqual(round(st.uM.Jmol, 1), 733.7)
        self.assertEqual(round(st.hM.Jmol, 1), 807.1)
        self.assertEqual(round(st.sM.JmolK, 2), 74.95)
        self.assertEqual(round(st.cvM.JmolK, 2), 17.70)
        self.assertEqual(round(st.cpM.JmolK, 2), 59.76)
        self.assertEqual(round(st.w, 1), 521.7)
        self.assertEqual(round(st.mu.muPas, 1), 91.9)
        self.assertEqual(round(st.k.WmK, 4), 0.0746)

        st = Ar(T=240, P=3e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 2), 62.93)
        self.assertEqual(round(st.rhoM, 3), 1.575)
        self.assertEqual(round(st.uM.Jmol, 0), 6224)
        self.assertEqual(round(st.hM.Jmol, 0), 8128)
        self.assertEqual(round(st.sM.JmolK, 1), 122.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.84)
        self.assertEqual(round(st.cpM.JmolK, 2), 23.49)
        self.assertEqual(round(st.w, 1), 288.7)
        self.assertEqual(round(st.mu.muPas, 1), 19.7)
        self.assertEqual(round(st.k.WmK, 4), 0.0164)

        st = Ar(T=150, P=4e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 209.7)
        self.assertEqual(round(st.rhoM, 3), 5.250)
        self.assertEqual(round(st.uM.Jmol, 0), 4317)
        self.assertEqual(round(st.hM.Jmol, 0), 5078)
        self.assertEqual(round(st.sM.JmolK, 1), 104.2)
        self.assertEqual(round(st.cvM.JmolK, 2), 17.85)
        self.assertEqual(round(st.cpM.JmolK, 2), 75.22)
        self.assertEqual(round(st.w, 1), 195.1)
        self.assertEqual(round(st.mu.muPas, 1), 15.6)
        self.assertEqual(round(st.k.WmK, 4), 0.0230)

        st = Ar(T=150.3, P=4.8e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 670.6)
        self.assertEqual(round(st.rhoM, 2), 16.79)
        self.assertEqual(round(st.uM.Jmol, 0), 2557)
        self.assertEqual(round(st.hM.Jmol, 0), 2843)
        self.assertEqual(round(st.sM.JmolK, 2), 88.49)
        self.assertEqual(round(st.cvM.JmolK, 2), 19.90)
        self.assertEqual(round(st.cpM.JmolK, 0), 1058)
        self.assertEqual(round(st.w, 1), 221.6)
        self.assertEqual(round(st.mu.muPas, 1), 35.8)
        self.assertEqual(round(st.k.WmK, 4), 0.0538)

        st = Ar(T=151.3, P=5e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 559.0)
        self.assertEqual(round(st.rhoM, 2), 13.99)
        self.assertEqual(round(st.uM.Jmol, 0), 2936)
        self.assertEqual(round(st.hM.Jmol, 0), 3293)
        self.assertEqual(round(st.sM.JmolK, 2), 91.38)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.18)
        self.assertEqual(round(st.cpM.JmolK, 0), 4808)
        self.assertEqual(round(st.w, 1), 195.4)
        self.assertEqual(round(st.mu.muPas, 1), 28.7)
        self.assertEqual(round(st.k.WmK, 4), 0.0759)

        st = Ar(T=90, P=6e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1395)
        self.assertEqual(round(st.rhoM, 2), 34.92)
        self.assertEqual(round(st.uM.Jmol, 0), -1165)
        self.assertEqual(round(st.hM.Jmol, 1), -992.7)
        self.assertEqual(round(st.sM.JmolK, 2), 57.25)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.75)
        self.assertEqual(round(st.cpM.JmolK, 2), 43.80)
        self.assertEqual(round(st.w, 1), 829.6)
        self.assertEqual(round(st.mu.muPas, 0), 255)
        self.assertEqual(round(st.k.WmK, 3), 0.127)

        st = Ar(T=186, P=7e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 252.0)
        self.assertEqual(round(st.rhoM, 3), 6.308)
        self.assertEqual(round(st.uM.Jmol, 0), 4723)
        self.assertEqual(round(st.hM.Jmol, 0), 5833)
        self.assertEqual(round(st.sM.JmolK, 1), 105.6)
        self.assertEqual(round(st.cvM.JmolK, 2), 15.32)
        self.assertEqual(round(st.cpM.JmolK, 2), 46.31)
        self.assertEqual(round(st.w, 1), 243.6)
        self.assertEqual(round(st.mu.muPas, 1), 19.6)
        self.assertEqual(round(st.k.WmK, 4), 0.0216)

        st = Ar(T=400, P=8e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 2), 95.74)
        self.assertEqual(round(st.rhoM, 3), 2.397)
        self.assertEqual(round(st.uM.Jmol, 0), 8152)
        self.assertEqual(round(st.hM.Jmol, -1), 11490)
        self.assertEqual(round(st.sM.JmolK, 1), 125.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.81)
        self.assertEqual(round(st.cpM.JmolK, 2), 22.75)
        self.assertEqual(round(st.w, 1), 386.9)
        self.assertEqual(round(st.mu.muPas, 1), 30.3)
        self.assertEqual(round(st.k.WmK, 4), 0.0250)

        st = Ar(T=188, P=9e6, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 351.3)
        self.assertEqual(round(st.rhoM, 3), 8.795)
        self.assertEqual(round(st.uM.Jmol, 0), 4365)
        self.assertEqual(round(st.hM.Jmol, 0), 5389)
        self.assertEqual(round(st.sM.JmolK, 1), 101.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 16.06)
        self.assertEqual(round(st.cpM.JmolK, 2), 58.59)
        self.assertEqual(round(st.w, 1), 252.1)
        self.assertEqual(round(st.mu.muPas, 1), 22.8)
        self.assertEqual(round(st.k.WmK, 4), 0.0255)

        st = Ar(T=145, P=1e7, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1016)
        self.assertEqual(round(st.rhoM, 2), 25.42)
        self.assertEqual(round(st.uM.Jmol, 0), 1302)
        self.assertEqual(round(st.hM.Jmol, 0), 1696)
        self.assertEqual(round(st.sM.JmolK, 2), 79.24)
        self.assertEqual(round(st.cvM.JmolK, 2), 17.43)
        self.assertEqual(round(st.cpM.JmolK, 2), 57.54)
        self.assertEqual(round(st.w, 1), 486.1)
        self.assertEqual(round(st.mu.muPas, 1), 76.6)
        self.assertEqual(round(st.k.WmK, 4), 0.0665)

        st = Ar(T=400, P=1.4e7, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 166.0)
        self.assertEqual(round(st.rhoM, 3), 4.155)
        self.assertEqual(round(st.uM.Jmol, 0), 7919)
        self.assertEqual(round(st.hM.Jmol, -1), 11290)
        self.assertEqual(round(st.sM.JmolK, 1), 120.3)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.03)
        self.assertEqual(round(st.cpM.JmolK, 2), 24.07)
        self.assertEqual(round(st.w, 1), 400.5)
        self.assertEqual(round(st.mu.muPas, 1), 31.7)
        self.assertEqual(round(st.k.WmK, 4), 0.0268)

        st = Ar(T=100, P=2e7, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 0), 1379)
        self.assertEqual(round(st.rhoM, 2), 34.52)
        self.assertEqual(round(st.uM.Jmol, 1), -894.8)
        self.assertEqual(round(st.hM.Jmol, 1), -315.4)
        self.assertEqual(round(st.sM.JmolK, 2), 60.14)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.78)
        self.assertEqual(round(st.cpM.JmolK, 2), 42.43)
        self.assertEqual(round(st.w, 1), 846.3)
        self.assertEqual(round(st.mu.muPas, 0), 226)
        self.assertEqual(round(st.k.WmK, 3), 0.121)

        st = Ar(T=310, P=4e7, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 577.0)
        self.assertEqual(round(st.rhoM, 2), 14.44)
        self.assertEqual(round(st.uM.Jmol, 0), 5354)
        self.assertEqual(round(st.hM.Jmol, 0), 8123)
        self.assertEqual(round(st.sM.JmolK, 1), 102.4)
        self.assertEqual(round(st.cvM.JmolK, 2), 14.12)
        self.assertEqual(round(st.cpM.JmolK, 2), 32.18)
        self.assertEqual(round(st.w, 1), 471.1)
        self.assertEqual(round(st.mu.muPas, 1), 41.9)
        self.assertEqual(round(st.k.WmK, 4), 0.0393)

        st = Ar(T=400, P=1e8, eq="younglove", visco=2, thermal=1)
        self.assertEqual(round(st.rho, 1), 787.4)
        self.assertEqual(round(st.rhoM, 2), 19.71)
        self.assertEqual(round(st.uM.Jmol, 0), 5999)
        self.assertEqual(round(st.hM.Jmol, -1), 11070)
        self.assertEqual(round(st.sM.JmolK, 1), 100.9)
        self.assertEqual(round(st.cvM.JmolK, 2), 14.82)
        self.assertEqual(round(st.cpM.JmolK, 2), 28.48)
        self.assertEqual(round(st.w, 1), 687.3)
        self.assertEqual(round(st.mu.muPas, 1), 62.1)
        self.assertEqual(round(st.k.WmK, 4), 0.0606)

    def test_Stewart(self):
        # Saturation pressures from Table 12, pag. 675
        st = Ar(T=84, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 0.0705)
        st = Ar(T=90, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 0.1336)
        st = Ar(T=95, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 0.2132)
        st = Ar(T=100, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 0.3240)
        st = Ar(T=105, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 0.47258)
        st = Ar(T=110, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 0.66574)
        st = Ar(T=115, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 0.91046)
        st = Ar(T=120, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 1.21391)
        st = Ar(T=125, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 1.5835)
        st = Ar(T=130, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 2.02700)
        st = Ar(T=135, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 2.55295)
        st = Ar(T=140, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 3.17100)
        st = Ar(T=145, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 3.89294)
        st = Ar(T=150, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 4.73599)

        # Selected point from Table 14, Pag 679, saturation states
        st = Ar(T=83.804, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 0.06896)
        self.assertEqual(round(st.Liquido.rhoM, 3), 35.475)
        self.assertEqual(round(st.Liquido.hM.kJkmol, 1), -4835.9)
        self.assertEqual(round(st.Liquido.sM.kJkmolK, 2), 53.29)
        self.assertEqual(round(st.Liquido.cvM.kJkmolK, 2), 21.34)
        self.assertEqual(round(st.Liquido.cpM.kJkmolK, 2), 42.61)
        self.assertEqual(round(st.Liquido.w, 0), 853)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.10154)
        self.assertEqual(round(st.Gas.hM.kJkmol, 1), 1701.4)
        self.assertEqual(round(st.Gas.sM.kJkmolK, 1), 131.3)
        self.assertEqual(round(st.Gas.w, 0), 209)

        st = Ar(T=90, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 5), 0.13361)
        self.assertEqual(round(st.Liquido.rhoM, 3), 34.538)
        self.assertEqual(round(st.Liquido.hM.kJkmol, 1), -4568.1)
        self.assertEqual(round(st.Liquido.sM.kJkmolK, 2), 56.35)
        self.assertEqual(round(st.Liquido.cvM.kJkmolK, 2), 20.59)
        self.assertEqual(round(st.Liquido.cpM.kJkmolK, 2), 43.49)
        self.assertEqual(round(st.Liquido.w, 0), 812)
        self.assertEqual(round(st.Gas.rhoM, 5), 0.18649)
        self.assertEqual(round(st.Gas.hM.kJkmol, 1), 1777.5)
        self.assertEqual(round(st.Gas.sM.kJkmolK, 2), 126.86)
        self.assertEqual(round(st.Gas.w, 0), 186)

        st = Ar(T=120, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 1.2139)
        self.assertEqual(round(st.Liquido.rhoM, 3), 29.123)
        self.assertEqual(round(st.Liquido.hM.kJkmol, 1), -3138.4)
        self.assertEqual(round(st.Liquido.sM.kJkmolK, 2), 69.68)
        self.assertEqual(round(st.Liquido.cvM.kJkmolK, 2), 18.16)
        self.assertEqual(round(st.Liquido.cpM.kJkmolK, 2), 53.27)
        self.assertEqual(round(st.Liquido.w, 0), 587)
        self.assertEqual(round(st.Gas.rhoM, 4), 1.5089)
        self.assertEqual(round(st.Gas.hM.kJkmol, 1), 1917.7)
        self.assertEqual(round(st.Gas.sM.kJkmolK, 2), 111.81)
        self.assertEqual(round(st.Gas.cvM.kJkmolK, 2), 16.75)
        self.assertEqual(round(st.Gas.cpM.kJkmolK, 2), 36.15)
        self.assertEqual(round(st.Gas.w, 0), 182)

        st = Ar(T=149, x=0.5, eq="stewart")
        self.assertEqual(round(st.P.MPa, 4), 4.5560)
        self.assertEqual(round(st.Liquido.rhoM, 3), 18.235)
        self.assertEqual(round(st.Liquido.hM.kJkmol, 1), -911.4)
        self.assertEqual(round(st.Liquido.sM.kJkmolK, 2), 84.99)
        self.assertEqual(round(st.Liquido.cvM.kJkmolK, 2), 22.45)
        self.assertEqual(round(st.Liquido.cpM.kJkmolK, 2), 366.38)
        self.assertEqual(round(st.Liquido.w, 0), 218)
        self.assertEqual(round(st.Gas.rhoM, 3), 8.564)
        self.assertEqual(round(st.Gas.hM.kJkmol, 1), 718.5)
        self.assertEqual(round(st.Gas.sM.kJkmolK, 2), 95.93)
        self.assertEqual(round(st.Gas.cvM.kJkmolK, 2), 24.79)
        self.assertEqual(round(st.Gas.cpM.kJkmolK, 1), 472.2)
        self.assertEqual(round(st.Gas.w, 0), 173)

        # Table 15, Pag 684, Single phase points
        st = Ar(T=84, P=8e4, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 35.447)
        self.assertEqual(round(st.uM.kJkmol, 1), -4829.6)
        self.assertEqual(round(st.hM.kJkmol, 1), -4827.4)
        self.assertEqual(round(st.sM.kJkmolK, 2), 53.39)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 21.31)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 42.64)
        self.assertEqual(round(st.w, 0), 852)

        st = Ar(T=1200, P=8e4, eq="stewart")
        self.assertEqual(round(st.rhoM, 5), 0.00802)
        self.assertEqual(round(st.uM.kJkmol, 0), 14965)
        self.assertEqual(round(st.hM.kJkmol, 0), 24944)
        self.assertEqual(round(st.sM.kJkmolK, 2), 185.64)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 12.47)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 20.79)
        self.assertEqual(round(st.w, 0), 645)

        st = Ar(T=84, P=1e5, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 35.448)
        self.assertEqual(round(st.uM.kJkmol, 1), -4829.8)
        self.assertEqual(round(st.hM.kJkmol, 1), -4827.0)
        self.assertEqual(round(st.sM.kJkmolK, 2), 53.39)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 21.31)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 42.63)
        self.assertEqual(round(st.w, 0), 852)

        st = Ar(T=150, P=1.5e5, eq="stewart")
        self.assertEqual(round(st.rhoM, 5), 0.12154)
        self.assertEqual(round(st.uM.kJkmol, 1), 1845.3)
        self.assertEqual(round(st.hM.kJkmol, 1), 3079.4)
        self.assertEqual(round(st.sM.kJkmolK, 2), 137.02)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 12.58)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 21.23)
        self.assertEqual(round(st.w, 0), 227)

        st = Ar(T=1200, P=2e5, eq="stewart")
        self.assertEqual(round(st.rhoM, 5), 0.02004)
        self.assertEqual(round(st.uM.kJkmol, 0), 14964)
        self.assertEqual(round(st.hM.kJkmol, 0), 24945)
        self.assertEqual(round(st.sM.kJkmolK, 2), 178.02)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 12.47)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 20.79)
        self.assertEqual(round(st.w, 0), 645)

        st = Ar(T=100, P=5e5, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 32.937)
        self.assertEqual(round(st.uM.kJkmol, 1), -4133.5)
        self.assertEqual(round(st.hM.kJkmol, 1), -4118.4)
        self.assertEqual(round(st.sM.kJkmolK, 2), 60.98)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 19.56)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 45.42)
        self.assertEqual(round(st.w, 0), 744)

        st = Ar(T=116, P=1e6, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 29.967)
        self.assertEqual(round(st.uM.kJkmol, 1), -3380.9)
        self.assertEqual(round(st.hM.kJkmol, 1), -3347.5)
        self.assertEqual(round(st.sM.kJkmolK, 2), 67.97)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 18.37)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 50.98)
        self.assertEqual(round(st.w, 0), 621)

        st = Ar(T=130, P=2e6, eq="stewart")
        self.assertEqual(round(st.rhoM, 4), 2.5433)
        self.assertEqual(round(st.uM.kJkmol, 1), 1030.3)
        self.assertEqual(round(st.hM.kJkmol, 1), 1816.7)
        self.assertEqual(round(st.sM.kJkmolK, 2), 107.82)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 17.89)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 47.16)
        self.assertEqual(round(st.w, 0), 183)

        st = Ar(T=150, P=5e6, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 18.975)
        self.assertEqual(round(st.uM.kJkmol, 1), -1229.6)
        self.assertEqual(round(st.hM.kJkmol, 2), -966.13)
        self.assertEqual(round(st.sM.kJkmolK, 2), 84.46)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 21.06)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 209.4)
        self.assertEqual(round(st.w, 0), 245)

        st = Ar(T=100, P=1e7, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 33.825)
        self.assertEqual(round(st.uM.kJkmol, 1), -4265.2)
        self.assertEqual(round(st.hM.kJkmol, 1), -3969.6)
        self.assertEqual(round(st.sM.kJkmolK, 2), 59.62)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 19.92)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 42.91)
        self.assertEqual(round(st.w, 0), 802)

        st = Ar(T=310, P=5e7, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 16.752)
        self.assertEqual(round(st.uM.kJkmol, 1), 1583.6)
        self.assertEqual(round(st.hM.kJkmol, 1), 4568.4)
        self.assertEqual(round(st.sM.kJkmolK, 2), 98.28)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 14.54)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 32.29)
        self.assertEqual(round(st.w, 0), 519)

        st = Ar(T=106.75, P=1e8, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 37.796)
        self.assertEqual(round(st.uM.kJkmol, 1), -4624.0)
        self.assertEqual(round(st.hM.kJkmol, 1), -1978.2)
        self.assertEqual(round(st.sM.kJkmolK, 2), 54.71)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 22.11)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 36.10)
        self.assertEqual(round(st.w, 0), 1064)

        st = Ar(T=850, P=5e8, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 27.468)
        self.assertEqual(round(st.uM.kJkmol, 1), 8752.3)
        self.assertEqual(round(st.hM.kJkmol, 0), 26955)
        self.assertEqual(round(st.sM.kJkmolK, 2), 103.80)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 15.18)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 24.23)
        self.assertEqual(round(st.w, 0), 1304)

        st = Ar(T=253.52, P=1e9, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 45.826)
        self.assertEqual(round(st.uM.kJkmol, 1), -1222.0)
        self.assertEqual(round(st.hM.kJkmol, 0), 20600)
        self.assertEqual(round(st.sM.kJkmolK, 2), 62.61)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 25.18)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 33.75)
        self.assertEqual(round(st.w, 0), 1937)

        st = Ar(T=1200, P=1e9, eq="stewart")
        self.assertEqual(round(st.rhoM, 3), 32.059)
        self.assertEqual(round(st.uM.kJkmol, 0), 14382)
        self.assertEqual(round(st.hM.kJkmol, 0), 45574)
        self.assertEqual(round(st.sM.kJkmolK, 2), 105.67)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 15.44)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 23.57)
        self.assertEqual(round(st.w, 0), 1715)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = Ar(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.5203)
        self.assertEqual(round(st.P.MPa, 3), 31.922)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.5630)

        st2 = Ar(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 25.97)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.18479)

    def test_LemmonTransport(self):
        # Table V, pag 28
        # Viscosity
        self.assertEqual(round(Ar(T=100, rhom=0).mu.muPas, 5), 8.18940)
        self.assertEqual(round(Ar(T=300, rhom=0).mu.muPas, 4), 22.7241)
        self.assertEqual(round(Ar(T=100, rhom=33).mu.muPas, 3), 184.232)
        self.assertEqual(round(Ar(T=200, rhom=10).mu.muPas, 4), 25.5662)
        self.assertEqual(round(Ar(T=300, rhom=5).mu.muPas, 4), 26.3706)
        self.assertEqual(round(Ar(T=150.69, rhom=13.4).mu.muPas, 4), 27.6101)

        # Thermal Conductivity
        self.assertEqual(round(Ar(rhom=0, T=100).k.mWmK, 5), 6.36587)
        self.assertEqual(round(Ar(rhom=0, T=300).k.mWmK, 4), 17.8042)
        self.assertEqual(round(Ar(rhom=33, T=100).k.mWmK, 3), 111.266)
        self.assertEqual(round(Ar(rhom=10, T=200).k.mWmK, 4), 26.1377)
        self.assertEqual(round(Ar(rhom=5, T=300).k.mWmK, 4), 23.2302)
        self.assertEqual(round(Ar(rhom=13.4, T=150.69).k.mWmK, 1), 856.8)

    def test_YoungloveHanley(self):
        st = Ar(T=90, P=1e5, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 7.2)
        self.assertEqual(round(st.k.mWmK, 1), 5.5)

        st = Ar(T=90, P=5e6, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 249.6)
        self.assertEqual(round(st.k.mWmK, 1), 126.7)

        st = Ar(T=200, P=1e6, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 16.1)
        self.assertEqual(round(st.k.mWmK, 1), 12.7)

        st = Ar(T=300, P=5e5, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 22.9)
        self.assertEqual(round(st.k.mWmK, 1), 18.0)

        st = Ar(T=400, P=2e6, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 29.1)
        self.assertEqual(round(st.k.mWmK, 1), 23.1)

        st = Ar(T=500, P=3e6, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 34.6)
        self.assertEqual(round(st.k.mWmK, 1), 27.5)

        st = Ar(T=90, P=6e6, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 251.8)
        self.assertEqual(round(st.k.mWmK, 1), 127.2)

        st = Ar(T=200, P=5e7, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 78.4)
        self.assertEqual(round(st.k.mWmK, 1), 69.9)

        st = Ar(T=150, P=3e7, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 103.3)
        self.assertEqual(round(st.k.mWmK, 1), 83.0)

        st = Ar(T=320, P=1e8, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 71.1)
        self.assertEqual(round(st.k.mWmK, 1), 66.4)

        st = Ar(T=370, P=2e8, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 102.2)
        self.assertEqual(round(st.k.mWmK, 1), 92.6)

        st = Ar(T=500, P=1e7, eq="younglove", visco=1, thermal=3)
        self.assertEqual(round(st.mu.muPas, 1), 35.6)
        self.assertEqual(round(st.k.mWmK, 1), 29.4)
