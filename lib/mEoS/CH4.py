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

from lib import unidades
from lib.meos import MEoS


class CH4(MEoS):
    """Multiparameter equation of state for methane"""
    name = "methane"
    CASNumber = "74-82-8"
    formula = "CH4"
    synonym = "R-50"
    _refPropName = "METHANE"
    _coolPropName = "Methane"
    rhoc = unidades.Density(162.66)
    Tc = unidades.Temperature(190.564)
    Pc = unidades.Pressure(4599.2, "kPa")
    M = 16.0428  # g/mol
    Tt = unidades.Temperature(90.694)
    Tb = unidades.Temperature(111.667)
    f_acent = 0.01142
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 2
    _Tr = unidades.Temperature(186.659809)
    _rhor = unidades.Density(163.413536)
    _w = 0.010528102

    Fi1 = {"ao_log": [1, 3.00160],
           "pow": [0, 1], "ao_pow": [9.91243972, -6.33270087],
           "ao_exp": [0.008449, 4.6942, 3.4865, 1.6572, 1.4115],
           "titao": [648/Tc, 1957/Tc, 3895/Tc, 5705/Tc, 15080/Tc]}

    Fi2 = {"R": 8.314510,
           "ao_log": [1, 3.00088],
           "pow": [0, 1], "ao_pow": [19.597508817, -83.959667892],
           "ao_exp": [], "titao": [],
           "ao_sinh": [0.76315, 8.74432], "sinh": [4.306474465, 5.577233895],
           "ao_cosh": [0.0046, -4.46921], "cosh": [0.936220902, 5.722644361]}

    Fi3 = {"ao_log": [1, 2.5998324],
           "pow": [0, -1./3, -2./3, -1],
           "ao_pow": [-10.413865, -3.3854083, 1.6900979, -0.3911541],
           "ao_exp": [4.7206715], "titao": [10.543907]}

    CP4 = {"ao": 0.15438149595e2,
           "an": [-0.18044750507e7, 0.77426666393e5, -0.13241658754e4,
                  -0.51479005257e-1, 0.10809172196e-3, -0.65501783437e-7],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-0.67490056171e1], "exp": [3000]}

    setzmann = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Setzmann and "
                    "Wagner (1991)",
        "__doi__": {
            "autor": "Setzmann, U., Wagner, W.",
            "title": "A New Equation of State and Tables of Thermodynamic "
                     "Properties for Methane Covering the Range from the "
                     "Melting Line to 625 K at Pressures up to 1000 MPa",
            "ref": "J. Phys. Chem. Ref. Data, 20(6) (1991) 1061-1155",
            "doi": "10.1063/1.555898"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072,

        "nr1": [0.43679010280e-1, 0.67092361990, -0.17655778590e01,
                0.85823302410, -0.12065130520e01, 0.51204672200,
                -0.40000107910e-3, -0.12478424230e-1, 0.31002697010e-1,
                0.17547485220e-2, -0.31719216050e-5, -0.22403468400e-5,
                0.29470561560e-6],
        "d1": [1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 8, 9, 10],
        "t1": [-0.5, 0.5, 1., 0.5, 1., 1.5, 4.5, 0., 1., 3., 1., 3., 3.],

        "nr2": [0.18304879090, 0.15118836790, -0.42893638770, 0.68940024460e-1,
                -0.14083139960e-1, -0.30630548300e-1, -0.29699067080e-1,
                -0.19320408310e-1, -0.11057399590, 0.99525489950e-1,
                0.85484378250e-2, -0.61505556620e-1, -0.42917924230e-1,
                -0.18132072900e-1, 0.34459047600e-1, -0.23859194500e-2,
                -0.11590949390e-1, 0.66416936020e-1, -0.23715495900e-1,
                -0.39616249050e-1, -0.13872920440e-1, 0.33894895990e-1,
                -0.29273787530e-2],
        "d2": [1, 1, 1, 2, 4, 5, 6, 1, 2, 3, 4, 4, 3, 5, 5, 8, 2, 3, 4, 4, 4,
               5, 6],
        "t2": [0., 1., 2., 0., 0., 2., 2., 5., 5., 5., 2., 4., 12., 8., 10.,
               10., 10., 14., 12., 18., 22., 18., 14.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4,
               4, 4],
        "gamma2": [1]*23,

        "nr3": [9.324799946e-5, -6.287171518, 12.71069467, -6.423953466],
        "d3": [2, 0, 0, 0],
        "t3": [2., 0., 1., 2.],
        "alfa3": [20, 40, 40, 40],
        "beta3": [200, 250, 250, 250],
        "gamma3": [1.07, 1.11, 1.11, 1.11],
        "epsilon3": [1]*4}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethane of Younglove and Ely "
                    "(1987)",
        "__doi__": {
            "autor": "Younglove, B.A., Ely, J.F.",
            "title": "Thermophysical Properties of Fluids. II. Methane, "
                     "Ethane, Propane, Isobutane, and Normal Butane",
            "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
            "doi": "10.1063/1.555785"},

        "Tmin": 90.68, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 36.2029,

        "R": 8.31434,
        "M": 16.043, "Tt": 90.68, "Tc": 190.53, "Pc": 4.59797, "rhoc": 10.15,

        "cp": CP4,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 10068.3, "so": 186.447},

        "b": [None, 0.9898937956e-4, 0.2199608275, -0.5322788000e1,
              0.2021657962e3, -0.2234398926e5, 0.106794028e-3, 0.1457922469e-2,
              -9.265816666, 0.2915364732e4, 0.2313546209e-5, 0.1387214274e-2,
              0.4780467451e-1, 0.1176103833e-3, -0.198209673e-2, -0.2512887756,
              0.9748899826e-4, -0.1202192137e-5, 0.4128353939e-3,
              -0.7215842918e-5, 0.5081738255e4, -0.9198903192e6, -27.32264677,
              0.7499024351e6, 0.01114060908, 0.1083955159e2, -0.4490960312e-3,
              -0.1380337847e2, -0.2371902232e-6, 0.3761652197e-3,
              -0.2375166954e-8, -0.1237640790e-6, 0.6766926453e-5]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {
            "autor": "Kunz, O., Wagner, W.",
            "title": "The GERG-2008 Wide-Range Equation of State for Natural "
                     "Gases and Other Mixtures: An Expansion of GERG-2004",
            "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
            "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": 90.6941, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072,

        "nr1": [0.57335704239162, -0.16760687523730e1, 0.23405291834916,
                -0.21947376343441, 0.16369201404128e-1, 0.15004406389280e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.98990489492918e-1, 0.58382770929055, -0.7478686756039,
                0.30033302857974, 0.20985543806568, -0.18590151133061e-1,
                -0.15782558339049, 0.12716735220791, -0.32019743894346e-1,
                -0.68049729364536e-1, 0.24291412853736e-1, 0.51440451639444e-2,
                -0.019084949733532, 0.55229677241291e-2, -0.44197392976085e-2,
                0.040061416708429, -0.33752085907575e-1, -0.25127658213357e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1]*18}

    friend = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Friend (1989)",
        "__doi__": {"autor": "Friend, D.G., Ely, J.F., Ingham, H.",
                    "title": "Thermophysical Properties of Methane",
                    "ref": "J. Phys. Chem. Ref. Data 18(2) (1989) 583-638",
                    "doi": "10.1063/1.555828"},

        "R": 8.31451,
        "cp": Fi3,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 10017.7, "so": 186.266},
        "Tt": 90.6854, "Tc": 190.551, "rhoc": 10.139, "M": 16.043,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 100000.0, "rhomax": 29.714,

        "nr1": [0.384436099659, -0.179692598800e1, 0.329444947369,
                0.226312728442e-1, 0.759236768798e-1, 0.693758447259e-1,
                0.241163263947e-1, 0.107009920854e-1, -0.380933275164e-1,
                0.471537561143e-3, 0.556607678805e-3, 0.548759346533e-6,
                -0.999632699967e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.128087979280, 0.380198873377e-1, 0.139226650551,
                -0.874996348859e-1, -0.334894165760e-2, -0.517576297122e-1,
                0.252835179116e-1, 0.518703205950e-3, -0.166770594525e-2,
                -0.607401927389e-3, -0.972915359991e-4, -0.298844010462e-4,
                -0.130940111124e-1, 0.198175833798e-1, 0.208465762327e-1,
                -0.358025052631e-1, -0.203486851741, 0.215964755088,
                -0.429340628249e-2],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*19}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for methane of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",
        "M": 16.043,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 40.072,

        "nr1": [0.89269676, -0.25438282e1, 0.64980978, 0.20793471e-1,
                0.70189104e-1, 0.23700378e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.16653334, -0.43855669e-1, -0.1572678, -0.35311675e-1,
                -0.29570024e-1, 0.14019842e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072,

        "nr1": [1.25595787, 8.48007435e-1, -3.00939285, 5.99544996e-2,
                2.57003062e-4, -2.85914246e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-6.83210861e-2, -3.47523515e-2, 1.04637008e-1, -1.09884198e-2,
                -0.125124331, -5.53450960e-3, -1.51182884e-2, -2.04800000e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = setzmann, younglove, GERG, friend, shortSpan, sun
    _PR = [-0.5997, -12.094]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [0.09313, -0.05714], "exp": [1.367, 1.517]}
    _dielectric = {
        "eq": 1,
        "a": [6.5443, 0.0133], "b": [8.4578, 3.7196], "c": [-352.97, -100.65],
        "Au": 0, "D": 2}

    _melting = {
        "eq": 1,
        "__doi__": setzmann["__doi__"],

        "Tmin": Tt, "Tmax": 625.0,
        "Tref": Tt, "Pref": 11696,
        "a0": 1,
        "a2": [2.47568e4, -7.36602e3], "exp2": [1.85, 2.1]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.036219, 1.409353, -0.4945199, -1.443048],
        "t": [1, 1.5, 2, 4.5]}
    _liquid_Density = {
        "eq": 2,
        "n": [1.9906389, -0.78756197, 0.036976723],
        "t": [0.354, 0.5, 2.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.880284, -2.8526531, -3.000648, -5.251169, -13.191859,
              -37.553961],
        "t": [0.354, 5/6, 3/2, 5/2, 25/6, 47/6]}

    visco0 = {"__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 190.564,
              "no": [2.60536, -18.5247, 23.4216],
              "to": [0, 0.25, 0.5],

              "a": [-3.12118e-5, 1.99422e-7, 0],
              "b": [5.98858e-5, -4.91143e-5, 0],
              "c": [3.49668e-5, -1.73176e-5, 0],
              "A": [-8.52992e-10, -3.58009e-10, 0],
              "B": [1.60099e-8, 8.50221e-10, 0],
              "C": [-3.55631e-7, 2.80326e-7, 0]}

    visco1 = {"__name__": "Vogel (2000)",
              "__doi__": {
                  "autor": "Vogel, E., Wilhelm, J., Küchenmeister, C., "
                           "Jaesche, M.",
                  "title": "High-precision viscosity measurements on methane",
                  "ref": "High Temperatures-High Pressures 32(1) (2000) 73-81",
                  "doi": "10.1068/htwu359"},

              "eq": 1, "omega": 1,

              "ek": 160.78, "sigma": 0.37333,
              "n_chapman": 0.021357,
              "collision": [0.215309028, -0.46256942, 0.051313823,
                            0.030320660, -0.0070047029],

              "Tref_virial": 159.7,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 190.564, "rhoref_res": 10.139*M,
              "nr": [-3.02256904347, 17.6965130175, 3.11150846518,
                     -21.5685107769, 0.672852409238, 10.2387524315,
                     -1.09330775541, -1.20030749419],
              "tr": [0, 1, 0, 1, 0, 1, 0, 1],
              "dr": [2, 2, 3, 3, 4, 4, 5, 5],
              "gr": [0, 0, 0, 0, 0, 0, 0, 0],
              "cr": [0, 0, 0, 0, 0, 0, 0, 0],

              "CPf": -21.1009923406,
              "CPg1": 5.10860501398,
              "CPgi": [], "CPti": []}

    visco2 = {"__name__": "Younglove (1987)",
              "__doi__": {
                  "autor": "Younglove, B.A., Ely, J.F.",
                  "title": "Thermophysical Properties of Fluids. II. Methane, "
                           "Ethane, Propane, Isobutane, and Normal Butane",
                  "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                  "doi": "10.1063/1.555785"},

              "eq": 2, "omega": 2,

              "ek": 168., "sigma": 0.368,

              "F": [0.16969859271, -0.13337234608e-1, 1.4, 168],
              "E": [-0.1620427429e2, 0.4270589027e3, 0.1402596278e2,
                    -0.3916837745e4, -0.347709909e-1, 0.2136542674e2,
                    0.1436802482e4],
              "rhoc": 10.15}

    visco3 = {"__name__": "Friend (1989)",
              "__doi__": {
                  "autor": "Friend, D.G., Ely, J.F., Ingham, H.",
                  "title": "Thermophysical Properties of Methane",
                  "ref": "J. Phys. Chem. Ref. Data 18(2) (1989) 583-638",
                  "doi": "10.1063/1.555828"},

              "eq": 1, "omega": 2,

              "Tref": 174., "muref": 10.0,
              "ek": 174., "sigma": 0.36652,
              "n_chapman": 10.50/M**0.5*0.36652**2,

              "Tref_res": 190.551, "rhoref_res": 10.139*M, "muref_res": 12.149,
              "nr_num": [0.41250137, -0.14390912, 0.10366993, 0.40287464,
                         -0.24903524, -0.12953131, 0.06575776, 0.02566628,
                         -0.03716526],
              "tr_num": [0, 1, 0, 1, 1.5, 0, 2, 0, 1],
              "dr_num": [1, 1, 2, 2, 2, 3, 3, 4, 4],
              "gr_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "cr_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "nr_den": [1.0, -0.38798341, 0.03533815],
              "tr_den": [0, 0, 1],
              "dr_den": [0, 1, 1],
              "gr_den": [0, 0, 0],
              "cr_den": [0, 0, 0]}

    _viscosity = visco0, visco1, visco2, visco3

    thermo0 = {"__name__": "Friend (1989)",
               "__doi__": {
                   "autor": "Friend, D.G., Ely, J.F., Ingham, H.",
                   "title": "Thermophysical Properties of Methane",
                   "ref": "J. Phys. Chem. Ref. Data 18(2) (1989) 583-638",
                   "doi": "10.1063/1.555828"},

               "eq": 1,

               "Toref": 174., "koref": 1e-3,
               "visco": visco3,
               "no_viscoCp": [1.45885, -0.4377162],
               "to_viscoCp": [0, -1],

               "Tref_res": 190.551, "rhoref_res": 10.139*16.043,
               "kref_res": 6.29638e-3,
               "nr": [2.4149207, 0.55166331, -0.52837734, 0.073809553,
                      0.24465507, -0.047613626],
               "tr": [0, 0, 0, 1, 0, 1],
               "dr": [1, 3, 4, 4, 5, 5],

               "nr_s": [1.5554612], "tr_s": [0], "dr_s": [2],

               "critical": "_ThCritical"}

    def _ThCritical(self, rho, T, fase):
        """Special critical enhancement for Friend et al. thermal conductivity
        correlation"""
        tau = self.Tc/T
        delta = rho/self.rhoc

        T_ = 1-1/tau                                                    # Eq 21
        d_ = 1-delta                                                    # Eq 22

        xt = 0.28631*delta*tau/(1 + 2*delta*fase.fird + delta**2*fase.firdd)

        F = exp(-2.646*abs(T_)**0.5 - 2.678*d_**2 + 0.637*d_)           # Eq 20

        # Eq 19
        tc = 91.855/fase.mu.muPas/tau**2*xt**0.4681*F * \
            (1 + delta*fase.fird - delta*tau*fase.firdt)**2
        return tc*1e-3

    thermo1 = {"__name__": "Younglove (1987)",
               "__doi__": {
                   "autor": "Younglove, B.A., Ely, J.F.",
                   "title": "Thermophysical Properties of Fluids. II. Methane,"
                            " Ethane, Propane, Isobutane, and Normal Butane",
                   "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                   "doi": "10.1063/1.555785"},

               "eq": 3,

               "ek": 168.,
               "G": [1.346953698, -0.3254677753],
               "E": [0.2325800819e-2, -0.2477927999, 0.3880593713e2,
                     -0.1579519146e-6, 0.3717991328e-2, -0.9616989434,
                     -0.3017352774e-1, 0.4298153386],

               "critical": 2,
               "Tc": 190.53, "rhoc": 10.15*16.043,
               "X": [37.42368, 3.16714, 0.78035, 0.60103],
               "Z": 6.512707e-10}

    _thermal = thermo0, thermo1


class Test(TestCase):
    """Testing"""
    def test_setzmann(self):
        """Selected point from Table 39, Pag 1117, saturation state"""
        st = CH4(T=90.694, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.011696)
        self.assertEqual(round(st.Liquido.rho, 2), 451.48)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -982.76)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -7.3868)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.1677)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.3678)
        self.assertEqual(round(st.Liquido.w, 1), 1538.6)
        self.assertEqual(round(st.Gas.rho, 5), 0.25074)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -438.50)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.3857)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.5735)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.1100)
        self.assertEqual(round(st.Gas.w, 2), 249.13)

        st = CH4(T=100, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.034376)
        self.assertEqual(round(st.Liquido.rho, 2), 438.89)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -951.21)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -7.0562)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.1136)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.4084)
        self.assertEqual(round(st.Liquido.w, 1), 1452.0)
        self.assertEqual(round(st.Gas.rho, 5), 0.67457)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -420.73)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -1.7514)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.5887)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.1458)
        self.assertEqual(round(st.Gas.w, 2), 260.09)

        st = CH4(T=110, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.088130)
        self.assertEqual(round(st.Liquido.rho, 2), 424.78)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -916.75)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -6.7290)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.0642)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.4692)
        self.assertEqual(round(st.Liquido.w, 1), 1354.7)
        self.assertEqual(round(st.Gas.rho, 4), 1.5982)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -402.92)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -2.0578)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.6108)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.2053)
        self.assertEqual(round(st.Gas.w, 2), 270.01)

        st = CH4(T=120, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.19143)
        self.assertEqual(round(st.Liquido.rho, 2), 409.90)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -881.54)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -6.4248)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.0196)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.5493)
        self.assertEqual(round(st.Liquido.w, 1), 1253.5)
        self.assertEqual(round(st.Gas.rho, 4), 3.2619)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -386.93)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -2.3030)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.6390)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.2930)
        self.assertEqual(round(st.Gas.w, 2), 277.76)

        st = CH4(T=140, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.64118)
        self.assertEqual(round(st.Liquido.rho, 2), 376.87)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -807.74)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -5.8653)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.9452)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 3.8129)
        self.assertEqual(round(st.Liquido.w, 1), 1037.7)
        self.assertEqual(round(st.Gas.rho, 3), 10.152)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -362.60)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -2.6857)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.7172)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.6108)
        self.assertEqual(round(st.Gas.w, 2), 285.93)

        st = CH4(T=160, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 1.5921)
        self.assertEqual(round(st.Liquido.rho, 2), 336.31)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -726.14)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -5.3391)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.9037)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 4.4354)
        self.assertEqual(round(st.Liquido.w, 2), 795.43)
        self.assertEqual(round(st.Gas.rho, 3), 25.382)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -353.87)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -3.0124)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.8473)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 3.4189)
        self.assertEqual(round(st.Gas.w, 2), 283.01)

        st = CH4(T=180, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 3.2852)
        self.assertEqual(round(st.Liquido.rho, 2), 276.23)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -625.00)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -4.7778)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.9669)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 7.2923)
        self.assertEqual(round(st.Liquido.w, 2), 497.01)
        self.assertEqual(round(st.Gas.rho, 3), 61.375)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -378.11)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -3.4062)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 2.1404)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 7.5740)
        self.assertEqual(round(st.Gas.w, 2), 266.04)

        st = CH4(T=190, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 4.5186)
        self.assertEqual(round(st.Liquido.rho, 2), 200.78)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -532.67)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -4.3082)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 2.6022)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 94.012)
        self.assertEqual(round(st.Liquido.w, 2), 250.31)
        self.assertEqual(round(st.Gas.rho, 2), 125.18)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -451.91)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), -3.8831)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 2.8546)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 140.81)
        self.assertEqual(round(st.Gas.w, 2), 238.55)

        # Selected points from Table 40, Pag 1119
        st = CH4(T=90.698, P=2.5e4)
        self.assertEqual(round(st.rho, 2), 451.48)
        self.assertEqual(round(st.u.kJkg, 2), -982.78)
        self.assertEqual(round(st.h.kJkg, 2), -982.73)
        self.assertEqual(round(st.s.kJkgK, 4), -7.3868)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.1677)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.3677)
        self.assertEqual(round(st.w, 1), 1538.7)

        st = CH4(T=115, P=1e5)
        self.assertEqual(round(st.rho, 4), 1.7341)
        self.assertEqual(round(st.u.kJkg, 2), -450.34)
        self.assertEqual(round(st.h.kJkg, 2), -392.68)
        self.assertEqual(round(st.s.kJkgK, 4), -2.0301)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6041)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.1957)
        self.assertEqual(round(st.w, 2), 276.16)

        st = CH4(T=110, P=101325)
        self.assertEqual(round(st.rho, 2), 424.79)
        self.assertEqual(round(st.u.kJkg, 2), -916.97)
        self.assertEqual(round(st.h.kJkg, 2), -916.74)
        self.assertEqual(round(st.s.kJkgK, 4), -6.7291)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.0642)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.4691)
        self.assertEqual(round(st.w, 1), 1354.8)

        st = CH4(T=620, P=5e5)
        self.assertEqual(round(st.rho, 4), 1.5546)
        self.assertEqual(round(st.u.kJkg, 2), 564.67)
        self.assertEqual(round(st.h.kJkg, 2), 886.30)
        self.assertEqual(round(st.s.kJkgK, 4), 1.1384)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.8272)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.3497)
        self.assertEqual(round(st.w, 2), 617.61)

        st = CH4(T=150, P=1e6)
        self.assertEqual(round(st.rho, 3), 15.536)
        self.assertEqual(round(st.u.kJkg, 2), -418.03)
        self.assertEqual(round(st.h.kJkg, 2), -353.67)
        self.assertEqual(round(st.s.kJkgK, 4), -2.8198)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.7556)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.8402)
        self.assertEqual(round(st.w, 2), 287.76)

        st = CH4(T=200, P=5e6)
        self.assertEqual(round(st.rho, 3), 87.764)
        self.assertEqual(round(st.u.kJkg, 2), -423.51)
        self.assertEqual(round(st.h.kJkg, 2), -366.54)
        self.assertEqual(round(st.s.kJkgK, 4), -3.4670)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.9965)
        self.assertEqual(round(st.cp.kJkgK, 4), 7.2726)
        self.assertEqual(round(st.w, 2), 291.29)

        st = CH4(T=100, P=1e7)
        self.assertEqual(round(st.rho, 2), 446.02)
        self.assertEqual(round(st.u.kJkg, 2), -957.84)
        self.assertEqual(round(st.h.kJkg, 2), -935.42)
        self.assertEqual(round(st.s.kJkgK, 4), -7.1235)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.1389)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.3437)
        self.assertEqual(round(st.w, 1), 1525.7)

        st = CH4(T=600, P=5e7)
        self.assertEqual(round(st.rho, 2), 134.47)
        self.assertEqual(round(st.u.kJkg, 2), 400.24)
        self.assertEqual(round(st.h.kJkg, 2), 772.07)
        self.assertEqual(round(st.s.kJkgK, 4), -1.5109)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.8062)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.5887)
        self.assertEqual(round(st.w, 2), 786.69)

        st = CH4(T=420, P=1e8)
        self.assertEqual(round(st.rho, 2), 277.63)
        self.assertEqual(round(st.u.kJkg, 2), -165.97)
        self.assertEqual(round(st.h.kJkg, 2), 194.23)
        self.assertEqual(round(st.s.kJkgK, 4), -3.1391)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.2570)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.1894)
        self.assertEqual(round(st.w, 1), 1096.5)

        st = CH4(T=620, P=1e9)
        self.assertEqual(round(st.rho, 2), 503.16)
        self.assertEqual(round(st.u.kJkg, 2), 350.14)
        self.assertEqual(round(st.h.kJkg, 2), 2337.56)
        self.assertEqual(round(st.s.kJkgK, 4), -3.1888)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.2951)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.8115)
        self.assertEqual(round(st.w, 1), 2915.6)

    def test_younglove(self):
        """Selected point from Appendix E, Pag 588, single phase region"""
        # The saturation point use the ancillary equation for calculate
        # pressure and density, so the values differ of values give by mBWR,
        # so not used for testing
        kw = {"eq": "younglove", "visco": 2, "thermal": 1}

        st = CH4(T=360, P=1e4, **kw)
        self.assertEqual(round(st.rho, 4), 0.05360)
        self.assertEqual(round(st.rhoM, 6), 0.003341)
        self.assertEqual(round(st.uM.kJkmol, 0), 9313)
        self.assertEqual(round(st.hM.kJkmol, -1), 12310)
        self.assertEqual(round(st.sM.kJkmolK, 1), 212.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 30.15)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 38.47)
        self.assertEqual(round(st.w, 1), 487.9)
        self.assertEqual(round(st.mu.muPas, 1), 13.1)
        self.assertEqual(round(st.k, 4), 0.0427)

        st = CH4(T=100, P=5e4, **kw)
        self.assertEqual(round(st.rho, 1), 439.1)
        self.assertEqual(round(st.rhoM, 2), 27.37)
        self.assertEqual(round(st.uM.kJkmol, 0), -5231)
        self.assertEqual(round(st.hM.kJkmol, 0), -5229)
        self.assertEqual(round(st.sM.kJkmolK, 2), 73.16)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 33.40)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 54.13)
        self.assertEqual(round(st.w, 0), 1450)
        self.assertEqual(round(st.mu.muPas, 0), 156)
        self.assertEqual(round(st.k, 3), 0.208)

        st = CH4(T=115, P=1e5, rho0=1, **kw)
        self.assertEqual(round(st.rho, 3), 1.735)
        self.assertEqual(round(st.rhoM, 4), 0.1081)
        self.assertEqual(round(st.uM.kJkmol, 0), 2792)
        self.assertEqual(round(st.hM.kJkmol, 0), 3717)
        self.assertEqual(round(st.sM.kJkmolK, 1), 153.7)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 25.65)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 35.14)
        self.assertEqual(round(st.w, 1), 276.2)
        self.assertEqual(round(st.mu.muPas, 2), 4.59)
        self.assertEqual(round(st.k, 4), 0.0123)

        st = CH4(T=300, P=101325, **kw)
        self.assertEqual(round(st.rho, 4), 0.6528)
        self.assertEqual(round(st.rhoM, 5), 0.04069)
        self.assertEqual(round(st.uM.kJkmol, 0), 7578)
        self.assertEqual(round(st.hM.kJkmol, -1), 10070)
        self.assertEqual(round(st.sM.kJkmolK, 1), 186.4)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 27.48)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 35.87)
        self.assertEqual(round(st.w, 1), 449.7)
        self.assertEqual(round(st.mu.muPas, 1), 11.2)
        self.assertEqual(round(st.k, 4), 0.0341)

        st = CH4(T=120, P=2e5, **kw)
        self.assertEqual(round(st.rho, 1), 410.1)
        self.assertEqual(round(st.rhoM, 2), 25.56)
        self.assertEqual(round(st.uM.kJkmol, 0), -4126)
        self.assertEqual(round(st.hM.kJkmol, 0), -4118)
        self.assertEqual(round(st.sM.kJkmolK, 2), 83.23)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 32.56)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 56.80)
        self.assertEqual(round(st.w, 0), 1248)
        self.assertEqual(round(st.mu.muPas, 1), 97.7)
        self.assertEqual(round(st.k, 3), 0.175)

        st = CH4(T=600, P=3e5, **kw)
        self.assertEqual(round(st.rho, 4), 0.9643)
        self.assertEqual(round(st.rhoM, 5), 0.06011)
        self.assertEqual(round(st.uM.kJkmol, -1), 18130)
        self.assertEqual(round(st.hM.kJkmol, -1), 23120)
        self.assertEqual(round(st.sM.kJkmolK, 1), 206.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 42.94)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 51.29)
        self.assertEqual(round(st.w, 1), 609.8)
        self.assertEqual(round(st.mu.muPas, 1), 19.4)
        self.assertEqual(round(st.k, 4), 0.0842)

        st = CH4(T=135, P=4e5, rho0=1, **kw)
        self.assertEqual(round(st.rho, 3), 6.255)
        self.assertEqual(round(st.rhoM, 4), 0.3899)
        self.assertEqual(round(st.uM.kJkmol, 0), 3167)
        self.assertEqual(round(st.hM.kJkmol, 0), 4193)
        self.assertEqual(round(st.sM.kJkmolK, 1), 146.6)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 26.53)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 38.29)
        self.assertEqual(round(st.w, 1), 289.2)
        self.assertEqual(round(st.mu.muPas, 2), 5.45)
        self.assertEqual(round(st.k, 4), 0.0153)

        st = CH4(T=350, P=5e5, **kw)
        self.assertEqual(round(st.rho, 3), 2.769)
        self.assertEqual(round(st.rhoM, 4), 0.1726)
        self.assertEqual(round(st.uM.kJkmol, 0), 8969)
        self.assertEqual(round(st.hM.kJkmol, -1), 11870)
        self.assertEqual(round(st.sM.kJkmolK, 1), 178.8)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 29.70)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 38.27)
        self.assertEqual(round(st.w, 1), 481.3)
        self.assertEqual(round(st.mu.muPas, 1), 12.8)
        self.assertEqual(round(st.k, 4), 0.0415)

        st = CH4(T=600, P=6e5, **kw)
        self.assertEqual(round(st.rho, 3), 1.928)
        self.assertEqual(round(st.rhoM, 4), 0.1202)
        self.assertEqual(round(st.uM.kJkmol, -1), 18120)
        self.assertEqual(round(st.hM.kJkmol, -1), 23110)
        self.assertEqual(round(st.sM.kJkmolK, 1), 201.2)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 42.95)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 51.35)
        self.assertEqual(round(st.w, 1), 610.3)
        self.assertEqual(round(st.mu.muPas, 1), 19.5)
        self.assertEqual(round(st.k, 4), 0.0843)

        st = CH4(T=100, P=8e5, **kw)
        self.assertEqual(round(st.rho, 1), 439.6)
        self.assertEqual(round(st.rhoM, 2), 27.40)
        self.assertEqual(round(st.uM.kJkmol, 0), -5239)
        self.assertEqual(round(st.hM.kJkmol, 0), -5210)
        self.assertEqual(round(st.sM.kJkmolK, 2), 73.07)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 33.44)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 54.04)
        self.assertEqual(round(st.w, 0), 1456)
        self.assertEqual(round(st.mu.muPas, 0), 158)
        self.assertEqual(round(st.k, 3), 0.209)

        st = CH4(T=150, P=1e6, rho0=1, **kw)
        self.assertEqual(round(st.rho, 2), 15.54)
        self.assertEqual(round(st.rhoM, 4), 0.9689)
        self.assertEqual(round(st.uM.kJkmol, 0), 3309)
        self.assertEqual(round(st.hM.kJkmol, 0), 4341)
        self.assertEqual(round(st.sM.kJkmolK, 0), 141)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 27.92)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 45.17)
        self.assertEqual(round(st.w, 1), 287.8)
        self.assertEqual(round(st.mu.muPas, 2), 6.24)
        self.assertEqual(round(st.k, 4), 0.0186)

        st = CH4(T=170, P=2e6, rho0=1, **kw)
        self.assertEqual(round(st.rho, 2), 30.57)
        self.assertEqual(round(st.rhoM, 3), 1.905)
        self.assertEqual(round(st.uM.kJkmol, 0), 3496)
        self.assertEqual(round(st.hM.kJkmol, 0), 4545)
        self.assertEqual(round(st.sM.kJkmolK, 1), 137.8)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 29.07)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 55.44)
        self.assertEqual(round(st.w, 1), 290.3)
        self.assertEqual(round(st.mu.muPas, 2), 7.41)
        self.assertEqual(round(st.k, 4), 0.0241)

        st = CH4(T=196, P=3e6, **kw)
        self.assertEqual(round(st.rho, 2), 38.88)
        self.assertEqual(round(st.rhoM, 3), 2.424)
        self.assertEqual(round(st.uM.kJkmol, 0), 4056)
        self.assertEqual(round(st.hM.kJkmol, 0), 5294)
        self.assertEqual(round(st.sM.kJkmolK, 1), 139.3)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 28.30)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 52.39)
        self.assertEqual(round(st.w, 1), 319.7)
        self.assertEqual(round(st.mu.muPas, 2), 8.68)
        self.assertEqual(round(st.k, 4), 0.0289)

        st = CH4(T=185.5, P=4e6, **kw)
        self.assertEqual(round(st.rho, 1), 250.9)
        self.assertEqual(round(st.rhoM, 2), 15.64)
        self.assertEqual(round(st.uM.kJkmol, 1), 339.4)
        self.assertEqual(round(st.hM.kJkmol, 1), 595.1)
        self.assertEqual(round(st.sM.kJkmolK, 1), 112.7)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 31.54)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 176.5)
        self.assertEqual(round(st.w, 1), 416.7)
        self.assertEqual(round(st.mu.muPas, 1), 25.7)
        self.assertEqual(round(st.k, 4), 0.1110)

        st = CH4(T=600, P=5e6, **kw)
        self.assertEqual(round(st.rho, 2), 15.93)
        self.assertEqual(round(st.rhoM, 4), 0.9930)
        self.assertEqual(round(st.uM.kJkmol, -1), 17930)
        self.assertEqual(round(st.hM.kJkmol, -1), 22970)
        self.assertEqual(round(st.sM.kJkmolK, 1), 183.2)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 43.09)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 52.11)
        self.assertEqual(round(st.w, 1), 619.3)
        self.assertEqual(round(st.mu.muPas, 1), 19.8)
        self.assertEqual(round(st.k, 4), 0.0861)

        st = CH4(T=100, P=6e6, **kw)
        self.assertEqual(round(st.rho, 1), 443.5)
        self.assertEqual(round(st.rhoM, 2), 27.64)
        self.assertEqual(round(st.uM.kJkmol, 0), -5296)
        self.assertEqual(round(st.hM.kJkmol, 0), -5079)
        self.assertEqual(round(st.sM.kJkmolK, 2), 72.50)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 33.64)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 53.51)
        self.assertEqual(round(st.w, 0), 1496)
        self.assertEqual(round(st.mu.muPas, 0), 168)
        self.assertEqual(round(st.k, 3), 0.213)

        st = CH4(T=410, P=1e7, **kw)
        self.assertEqual(round(st.rho, 2), 48.25)
        self.assertEqual(round(st.rhoM, 3), 3.008)
        self.assertEqual(round(st.uM.kJkmol, -1), 10150)
        self.assertEqual(round(st.hM.kJkmol, -1), 13480)
        self.assertEqual(round(st.sM.kJkmolK, 1), 158.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 33.45)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 45.50)
        self.assertEqual(round(st.w, 1), 529.0)
        self.assertEqual(round(st.mu.muPas, 1), 16.0)
        self.assertEqual(round(st.k, 4), 0.0572)

        st = CH4(T=120, P=2e7, **kw)
        self.assertEqual(round(st.rho, 1), 428.5)
        self.assertEqual(round(st.rhoM, 2), 26.71)
        self.assertEqual(round(st.uM.kJkmol, 0), -4405)
        self.assertEqual(round(st.hM.kJkmol, 0), -3657)
        self.assertEqual(round(st.sM.kJkmolK, 2), 80.77)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 33.26)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 53.59)
        self.assertEqual(round(st.w, 0), 1437)
        self.assertEqual(round(st.mu.muPas, 0), 123)
        self.assertEqual(round(st.k, 3), 0.196)

        st = CH4(T=540, P=3e7, **kw)
        self.assertEqual(round(st.rho, 2), 98.78)
        self.assertEqual(round(st.rhoM, 3), 6.157)
        self.assertEqual(round(st.uM.kJkmol, -1), 14310)
        self.assertEqual(round(st.hM.kJkmol, -1), 19180)
        self.assertEqual(round(st.sM.kJkmolK, 1), 161.2)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 41.16)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 53.72)
        self.assertEqual(round(st.w, 1), 675.5)
        self.assertEqual(round(st.mu.muPas, 1), 21.8)
        self.assertEqual(round(st.k, 4), 0.0888)

        st = CH4(T=200, P=5e7, **kw)
        self.assertEqual(round(st.rho, 1), 367.2)
        self.assertEqual(round(st.rhoM, 2), 22.89)
        self.assertEqual(round(st.uM.kJkmol, 1), -953.0)
        self.assertEqual(round(st.hM.kJkmol, 0), 1231)
        self.assertEqual(round(st.sM.kJkmolK, 1), 104.3)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 30.40)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 51.46)
        self.assertEqual(round(st.w, 0), 1208)
        self.assertEqual(round(st.mu.muPas, 1), 59.1)
        self.assertEqual(round(st.k, 3), 0.145)

        st = CH4(T=600, P=1e8, **kw)
        self.assertEqual(round(st.rho, 0), 214)
        self.assertEqual(round(st.rhoM, 2), 13.34)
        self.assertEqual(round(st.uM.kJkmol, -1), 15460)
        self.assertEqual(round(st.hM.kJkmol, -1), 22960)
        self.assertEqual(round(st.sM.kJkmolK, 1), 155.3)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 44.99)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 57.65)
        self.assertEqual(round(st.w, 0), 1023)
        self.assertEqual(round(st.mu.muPas, 1), 32.1)
        self.assertEqual(round(st.k, 3), 0.129)

        st = CH4(T=150, P=2e8, **kw)
        self.assertEqual(round(st.rho, 1), 484.6)
        self.assertEqual(round(st.rhoM, 2), 30.21)
        self.assertEqual(round(st.uM.kJkmol, 0), -4018)
        self.assertEqual(round(st.hM.kJkmol, 0), 2603)
        self.assertEqual(round(st.sM.kJkmolK, 2), 80.58)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 36.33)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 49.56)
        self.assertEqual(round(st.w, 0), 2110)
        self.assertEqual(round(st.mu.muPas, 0), 197)
        self.assertEqual(round(st.k, 3), 0.336)

    def test_friend(self):
        """Selected point from Table A1, Pag 630, ideal gas"""

        # The point really are not the gas ideal at zero preesure as it's
        # calculated here, it's the value at 1 bar pressure so the values of
        # Helmholtz free energy and entropy differ of ideal properties
        st = CH4(T=100, rho=0, eq="friend", visco=3)
        self.assertEqual(round(st.hM.kJmol, 3), 3.311)
        self.assertEqual(round(st.cpM.JmolK, 3), 33.277)
        self.assertEqual(round(st.mu.muPas, 2), 3.95)
        self.assertEqual(round(st.k.mWmK, 2), 9.83)

        st = CH4(T=400, rho=0, eq="friend", visco=3)
        self.assertEqual(round(st.hM.kJmol, 3), 13.888)
        self.assertEqual(round(st.cpM.JmolK, 3), 40.613)
        self.assertEqual(round(st.mu.muPas, 2), 14.23)
        self.assertEqual(round(st.k.mWmK, 2), 49.96)

        # Selected point from Table A2, Pag 631, saturation state
        # This table use ancillary equation for calculate presure and densities
        # of two phases so the result may differ
        st = CH4(T=92, x=0.5, eq="friend", visco=3)
        self.assertEqual(round(st.P.MPa, 3), 0.014)
        self.assertEqual(round(st.Liquido.rhoM, 2), 28.04)
        self.assertEqual(round(st.Gas.rhoM, 3), 0.018)
        self.assertEqual(round(st.Liquido.w, 1), 1532.8)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 194.92)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 209.7)

        st = CH4(T=120, x=0.5, eq="friend", visco=3)
        self.assertEqual(round(st.P.MPa, 3), 0.192)
        self.assertEqual(round(st.Liquido.rhoM, 2), 25.56)
        self.assertEqual(round(st.Gas.rhoM, 3), 0.204)
        self.assertEqual(round(st.Liquido.w, 1), 1247.5)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 98.09)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 172.2)

        st = CH4(T=160, x=0.5, eq="friend", visco=3)
        self.assertEqual(round(st.P.MPa, 3), 1.593)
        self.assertEqual(round(st.Liquido.rhoM, 2), 20.96)
        self.assertEqual(round(st.Gas.rhoM, 3), 1.584)
        self.assertEqual(round(st.Liquido.w, 1), 791.5)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 46.70)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 115.0)

        st = CH4(T=190, x=0.5, eq="friend", visco=3)
        self.assertEqual(round(st.P.MPa, 3), 4.521)
        self.assertEqual(round(st.Liquido.rhoM, 2), 12.55)
        self.assertEqual(round(st.Gas.rhoM, 3), 7.822)
        self.assertEqual(round(st.Liquido.w, 1), 266.1)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 19.84)
        # The paper use a scaled equation of state near of critical point
        # self.assertEqual(round(st.Liquido.k.mWmK, 1), 100.3)

        # Selected point from Table A3, Pag 632
        st = CH4(T=100, P=1e5, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 27.37)
        self.assertEqual(round(st.hM.kJmol, 3), -5.242)
        self.assertEqual(round(st.sM.JmolK, 2), 73.05)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.08)
        self.assertEqual(round(st.cpM.JmolK, 2), 54.64)
        self.assertEqual(round(st.w, 1), 1446.8)
        self.assertEqual(round(st.mu.muPas, 2), 157.49)
        self.assertEqual(round(st.k.mWmK, 1), 199.8)

        st = CH4(T=150, P=5e6, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 22.86)
        self.assertEqual(round(st.hM.kJmol, 3), -2.274)
        self.assertEqual(round(st.sM.JmolK, 2), 95.48)
        self.assertEqual(round(st.cvM.JmolK, 2), 30.84)
        self.assertEqual(round(st.cpM.JmolK, 2), 61.25)
        self.assertEqual(round(st.w, 1), 996.2)
        self.assertEqual(round(st.mu.muPas, 2), 60.01)
        self.assertEqual(round(st.k.mWmK, 1), 136.5)

        st = CH4(T=200, P=2e6, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 1.40)
        self.assertEqual(round(st.hM.kJmol, 3), 5.939)
        self.assertEqual(round(st.sM.JmolK, 2), 145.34)
        self.assertEqual(round(st.cvM.JmolK, 2), 26.77)
        self.assertEqual(round(st.cpM.JmolK, 2), 41.60)
        self.assertEqual(round(st.w, 1), 343.0)
        self.assertEqual(round(st.mu.muPas, 2), 8.28)
        self.assertEqual(round(st.k.mWmK, 1), 24.5)

        st = CH4(T=250, P=5e7, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 19.80)
        self.assertEqual(round(st.hM.kJmol, 3), 3.809)
        self.assertEqual(round(st.sM.JmolK, 2), 115.82)
        self.assertEqual(round(st.cvM.JmolK, 2), 30.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 51.95)
        self.assertEqual(round(st.w, 1), 999.3)
        self.assertEqual(round(st.mu.muPas, 2), 42.74)
        self.assertEqual(round(st.k.mWmK, 1), 114.8)

        st = CH4(T=300, P=1e5, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 0.04)
        self.assertEqual(round(st.hM.kJmol, 3), 10.068)
        self.assertEqual(round(st.sM.JmolK, 2), 186.56)
        self.assertEqual(round(st.cvM.JmolK, 2), 27.47)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.86)
        self.assertEqual(round(st.w, 1), 449.7)
        self.assertEqual(round(st.mu.muPas, 2), 11.20)
        self.assertEqual(round(st.k.mWmK, 1), 34.4)

        st = CH4(T=350, P=1e6, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 0.35)
        self.assertEqual(round(st.hM.kJmol, 3), 11.807)
        self.assertEqual(round(st.sM.JmolK, 2), 172.86)
        self.assertEqual(round(st.cvM.JmolK, 2), 29.75)
        self.assertEqual(round(st.cpM.JmolK, 2), 38.60)
        self.assertEqual(round(st.w, 1), 480.8)
        self.assertEqual(round(st.mu.muPas, 2), 12.89)
        self.assertEqual(round(st.k.mWmK, 1), 42.3)

        st = CH4(T=400, P=5e7, eq="friend", visco=3)
        self.assertEqual(round(st.rhoM, 2), 12.68)
        self.assertEqual(round(st.hM.kJmol, 3), 11.583)
        self.assertEqual(round(st.sM.JmolK, 2), 140.19)
        self.assertEqual(round(st.cvM.JmolK, 2), 34.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 51.82)
        self.assertEqual(round(st.w, 1), 778.9)
        self.assertEqual(round(st.mu.muPas, 2), 26.93)
        self.assertEqual(round(st.k.mWmK, 1), 88.7)

    def test_shortSpan(self):
        """Table III, Pag 46"""
        st = CH4(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.6278)
        self.assertEqual(round(st.P.MPa, 3), 108.108)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.9282)

        st2 = CH4(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 142.73)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.78166)
