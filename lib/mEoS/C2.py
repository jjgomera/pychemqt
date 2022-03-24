#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


class C2(MEoS):
    """Multiparameter equation of state for ethane"""
    name = "ethane"
    CASNumber = "74-84-0"
    formula = "CH3CH3"
    synonym = "R-170"
    _refPropName = "ETHANE"
    _coolPropName = "Ethane"
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
           "ao_sinh": [4.33939, 13.1974],
           "sinh": [559.314/Tc, 1031.38/Tc],
           "ao_cosh": [1.23722, -6.01989],
           "cosh": [223.284/Tc, 1071.29/Tc]}

    Fi3 = {"ao_log": [1, 3.8159476],
           "pow": [0, -1./3, -2./3, -1],
           "ao_pow": [-23.446765, 8.6021299, -3.3075735, -0.55956678],
           "ao_exp": [5.0722267], "titao": [5.5074874]}

    CP5 = {"ao": 9.9507922459,
           "an": [-6.9341406909e5, 3.1534834135e4, -6.103375287e2,
                  -2.8657877948e-2, 9.0922897821e-5, -5.2750109915e-8],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-1.4243593411e1], "exp": [3000]}

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

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethane of Younglove and Ely "
                    "(1987)",
        "__doi__": {"autor": "Younglove, B.A., Ely, J.F.",
                    "title": "Thermophysical Properties of Fluids. II. "
                             "Methane, Ethane, Propane, Isobutane, and Normal "
                             "Butane",
                    "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                    "doi": "10.1063/1.555785"},

        "R": 8.31434,
        "M": 30.07, "Tt": 90.348, "Tc": 305.34, "Pc": 4871.43, "rhoc": 6.875,

        "cp": CP5,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 11913.6, "so": 229.35},

        "Tmin": 90.348, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 21.68,

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
        "__doi__": {"autor": "Friend, D.G., Ingham, H., Ely, J.F.",
                    "title": "Thermophysical Properties of Ethane",
                    "ref": "J. Phys. Chem. Ref. Data 20(2) (1991) 275-347",
                    "doi": "10.1063/1.555881"},

        "R": 8.31451,
        "cp": Fi3,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 11874, "so": 229.12},
        "Tt": 90.352, "Tc": 305.33, "Pc": 4871.8, "rhoc": 6.87, "M": 30.07,

        "Tmin": 90.352, "Tmax": 625.0, "Pmax": 70000.0, "rhomax": 22.419,

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
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 675.0, "Pmax": 900000.0, "rhomax": 22.419,

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

    eq = buecker, younglove, GERG, friend, shortSpan, sun
    _PR = [-0.2955, -14.9384]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [-0.01484, 0.0619], "exp": [1.869, 1.285]}
    _dielectric = {
        "eq": 1,
        "a": [11.1552, 0.0112], "b": [36.759, 23.639], "c": [-808.03, -378.84],
        "Au": 0, "D": 1.75}

    _melting = {
        "eq": 1,
        "__doi__": buecker["__doi__"],

        "Tmin": Tt, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 1.14,
        "a0": 1,
        "a2": [1, 2.23626315e8, 1.05262374e8], "exp2": [0, 1, 2.55]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.48647577, 1.47010078, -1.66261122, 3.57898378, -4.79105705],
        "t": [1, 1.5, 2.5, 3.5, 4]}
    _liquid_Density = {
        "eq": 2,
        "n": [1.56138026, -0.381552776, 0.078537204, 0.0370315089],
        "t": [0.329, 4/6, 8/6, 19/6]}
    _vapor_Density = {
        "eq": 3,
        "n": [-1.89879145, -3.65459262, 0.850562745, 0.363965487, -1.50005943,
              -2.26690389],
        "t": [0.346, 5/6, 1, 2, 3, 5]}

    visco0 = {"__name__": "Herrmann (2018)",
              "__doi__": {
                  "autor": "Herrmann, S., Hellmann, R., Vogel, E.",
                  "title": "Update: Reference Correlation for the Viscosity "
                           "of Ethane [J. Phys. Chem. Ref. Data 44, 043101 "
                           "(2015)]",
                  "ref": "J. Phys. Chem. Ref. Data 44(4) (2015) 043101",
                  "doi": "10.1063/1.4930838"},

              "eq": 1, "omega": 0,

              "special0": "_mu0",

              "Tref_res": 305.322,
              "rhoref_res": 206.18,
              "nr": [2.0101502550505e1, -5.3310495236264, -9.1932109313455e-3,
                     -3.2692750337240e-1, 2.8850975314454e-1,
                     1.8714822381408e-4, 6.5739794167039e-6,
                     -1.3249185654669e1, -2.5958989867061e1, 2.3489171321456e1,
                     7.8770433769273],
              "tr": [0, 1, 7, 2, 0, 7, 1, 0, 1, 2, 0],
              "dr": [1, 1, 2, 4, 6, 6, 15, 1, 3, 3, 4],
              "gr": [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
              "cr": [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],

              "nr_gaus": [6.4266651197353e-1, 7.0693316317088e-1],
              "br_gaus": [90, 50],
              "er_gaus": [100, 250]}

    def _mu0(self, T):
        """Special zero-density correlation for the Herrmann viscosity
        correlation"""
        # Eq 6
        X = 0.79330 + 262.946/exp(T**(1/3)) + \
            (13.8366+1339.77/exp(T**(1/3)))*(1/T)**0.5 + \
            -322.242*T/exp(2*T**(1/3))

        # Eq 5
        muo = T**0.5/X

        return muo

    visco1 = {"__name__": "Vogel (2015)",
              "__doi__": {
                  "autor": "Vogel, E., Span, R., Herrmann, S.",
                  "title": "Reference Correlation for the Viscosity of Ethane",
                  "ref": "J. Phys. Chem. Ref. Data 44(4) (2015) 043101",
                  "doi": "10.1063/1.4930838"},

              "eq": 1, "omega": 0,

              "Toref": 305.322,
              "no": [9.6634694892149, -2.2985582151676e-1],
              "to": [1, 3],

              "Tref_res": 305.322,
              "rhoref_res": 206.18,
              "nr": [6.6687966976352, -4.6983342709702, 1.9688847427047e1,
                     -9.5399537393789, 6.3640646131666e-2, 7.9981217444542e-3,
                     7.0489675750657e-8, -2.2734655865556e1, 2.2124096051632e1,
                     -3.0986358885564e-1],
              "tr": [0, 1, 0, 1, 0, 1, 3, 0, 2, 5],
              "dr": [1, 1, 2, 2, 7, 8, 17, 3, 3, 3],
              "gr": [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
              "cr": [0, 0, 0, 0, 0, 0, 0, 1, 1, 1],

              "nr_gaus": [6.4034200732045e-1, 7.0437620805249e-1],
              "br_gaus": [90, 50],
              "er_gaus": [100, 250]}

    visco2 = {"__name__": "Friend (1991)",
              "__doi__": {
                  "autor": "Friend, D.G., Ingham, H., Ely, J.F.",
                  "title": "Thermophysical Properties of Ethane",
                  "ref": "J. Phys. Chem. Ref. Data 20(2) (1991) 275-347",
                  "doi": "10.1063/1.555881"},

              "eq": 1, "omega": 2,

              "ek": 245.0, "sigma": 0.43682,
              "Tref": 245.0,
              "n_chapman": 12.0085/M**0.5*0.43682**2,

              "muref_res": 15.977,
              "nr_num": [0.47177003, -0.23950311, 0.39808301, -0.27343335,
                         0.35192260, -0.21101308, -0.00478579, 0.07378129,
                         -0.030425255],
              "tr_num": [0, 1, 0, 1, 1.5, 0, 2, 0, 1],
              "dr_num": [1, 1, 2, 2, 2, 3, 3, 4, 4],
              "gr_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "cr_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "nr_den": [1., -0.30435286, 0.001215675],
              "tr_den": [0, 0, 1],
              "dr_den": [0, 1, 1],
              "gr_den": [0, 0, 0],
              "cr_den": [0, 0, 0]
              }

    visco3 = {"__name__": "Younglove (1987)",
              "__doi__": {
                  "autor": "Younglove, B.A., Ely, J.F.",
                  "title": "Thermophysical Properties of Fluids. II. Methane, "
                           "Ethane, Propane, Isobutane, and Normal Butane",
                  "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                  "doi": "10.1063/1.555785"},

              "eq": 2, "omega": 2,
              "ek": 240.0, "sigma": 0.440110,

              "F": [0.2102436247e1, -0.1065920192e1, 1.4, 305.33],
              "E": [-0.1903481042e2, 0.1799260494e4, 0.1561316986e2,
                    -0.1497221136e5, 0.1130374601, -0.2186440756e2,
                    0.8235954037e4],
              "rhoc": 6.875}

    visco4 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 305.322,
              "no": [15.9252, -49.7734, 43.4368],
              "to": [0, 0.25, 0.5],

              "a": [-7.50686e-6, -1.50327e-6, 0],
              "b": [6.72862e-5, -4.36451e-5, 0],
              "c": [3.88040e-5, -1.38524e-5, 0],
              "A": [7.68043e-10, -1.32048e-10, 0.0],
              "B": [9.15407e-9, 4.13028e-10, 0.0],
              "C": [-1.45842e-7, 2.39764e-7, 0.0]}

    _viscosity = visco0, visco1, visco2, visco3, visco4

    thermo0 = {"__name__": "Friend (1991)",
               "__doi__": {
                   "autor": "Friend, D.G., Ingham, H., Ely, J.F.",
                   "title": "Thermophysical Properties of Ethane",
                   "ref": "J. Phys. Chem. Ref. Data 20(2) (1991) 275-347",
                   "doi": "10.1063/1.555881"},

               "eq": 1,

               "Toref": 245.0, "koref": 1e-3,
               "no_viscoCp": [1.7104147, -0.6936482],
               "to_viscoCp": [0, -1],

               "Tref_res": 305.33, "rhoref_res": 6.87*30.07,
               "kref_res": 4.41786e-3,
               "nr": [0.96084322, 2.7500235, -0.026609289, -0.078146729,
                      0.21881339, 2.3849563, -0.75113971],
               "tr": [0, 0, 0, 0, 0, 1.5, 1],
               "dr": [1, 2, 3, 4, 5, 1, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.242, "R0": 1.01,
               "Xio": 0.19e-9, "gam0": 0.0563, "qd": 0.545e-9, "Tcref": 610.66}

    thermo1 = {"__name__": "Younglove (1987)",
               "__doi__": {
                   "autor": "Younglove, B.A., Ely, J.F.",
                   "title": "Thermophysical Properties of Fluids. II. Methane,"
                            " Ethane, Propane, Isobutane, and Normal Butane",
                   "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                   "doi": "10.1063/1.555785"},

               "eq": 3,

               "ek": 240,
               "G": [1.545691277, -0.5086287855],
               "E": [0.2863803648e-2, -0.459858003, 0.7772750057e2,
                     0.138460594e-5, 0.1874040714e-1, -0.3009947821e1,
                     -0.4225741011e-1, 0.1028764297e1],

               "critical": 2,
               "Tc": 305.34, "rhoc": 6.875*30.07,
               "X": [0.225388, 10.51088, 0.45, 1],
               "Z": 7.42399e-10}

    _thermal = thermo0, thermo1


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

    def test_younglove(self):
        # The saturation point use the ancillary equation for calculate
        # pressure and density, so the values differ of values give by mBWR,
        # so not used for testing
        kw = {"eq": "younglove", "visco": 3, "thermal": 1}

        # Selected point from Appendix F, Pag 642, single phase region
        st = C2(T=100, P=1e4, **kw)
        self.assertEqual(round(st.rho, 1), 641.1)
        self.assertEqual(round(st.rhoM, 2), 21.32)
        self.assertEqual(round(st.uM.kJkmol, -1), -14210)
        self.assertEqual(round(st.hM.kJkmol, -1), -14210)
        self.assertEqual(round(st.sM.kJkmolK, 2), 83.82)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 47.14)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 69.32)
        self.assertEqual(round(st.w, 0), 1943)
        self.assertEqual(round(st.mu.muPas, 0), 797)
        self.assertEqual(round(st.k, 3), 0.248)

        st = C2(T=175, P=5e4, rho0=1, **kw)
        self.assertEqual(round(st.rho, 3), 1.054)
        self.assertEqual(round(st.rhoM, 5), 0.03505)
        self.assertEqual(round(st.uM.kJkmol, 0), 4708)
        self.assertEqual(round(st.hM.kJkmol, 0), 6135)
        self.assertEqual(round(st.sM.kJkmolK, 1), 210.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 32.53)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 41.58)
        self.assertEqual(round(st.w, 1), 243.8)
        self.assertEqual(round(st.mu.muPas, 2), 5.58)
        self.assertEqual(round(st.k, 5), 0.00890)

        st = C2(T=600, P=1e5, **kw)
        self.assertEqual(round(st.rho, 4), 0.6031)
        self.assertEqual(round(st.rhoM, 5), 0.02006)
        self.assertEqual(round(st.uM.kJkmol, -1), 28430)
        self.assertEqual(round(st.hM.kJkmol, -1), 33420)
        self.assertEqual(round(st.sM.kJkmolK, 1), 277.7)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 80.90)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 89.26)
        self.assertEqual(round(st.w, 1), 427.6)
        self.assertEqual(round(st.mu.muPas, 1), 17.0)
        self.assertEqual(round(st.k, 4), 0.0695)

        st = C2(T=240, P=101325, rho0=1, **kw)
        self.assertEqual(round(st.rho, 3), 1.550)
        self.assertEqual(round(st.rhoM, 5), 0.05155)
        self.assertEqual(round(st.uM.kJkmol, 0), 6962)
        self.assertEqual(round(st.hM.kJkmol, 0), 8928)
        self.assertEqual(round(st.sM.kJkmolK, 1), 218.3)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 37.90)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 46.74)
        self.assertEqual(round(st.w, 1), 281.8)
        self.assertEqual(round(st.mu.muPas, 2), 7.56)
        self.assertEqual(round(st.k, 4), 0.0145)

        st = C2(T=200, P=2e5, rho0=1, **kw)
        self.assertEqual(round(st.rho, 3), 3.822)
        self.assertEqual(round(st.rhoM, 4), 0.1271)
        self.assertEqual(round(st.uM.kJkmol, 0), 5407)
        self.assertEqual(round(st.hM.kJkmol, 0), 6981)
        self.assertEqual(round(st.sM.kJkmolK, 1), 203.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 35.19)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 45.65)
        self.assertEqual(round(st.w, 1), 253.1)
        self.assertEqual(round(st.mu.muPas, 2), 6.32)
        self.assertEqual(round(st.k, 4), 0.0111)

        st = C2(T=200, P=3e5, **kw)
        self.assertEqual(round(st.rho, 1), 524.3)
        self.assertEqual(round(st.rhoM, 2), 17.44)
        self.assertEqual(round(st.uM.kJkmol, 0), -7107)
        self.assertEqual(round(st.hM.kJkmol, 0), -7090)
        self.assertEqual(round(st.sM.kJkmolK, 1), 132.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 42.85)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 75.14)
        self.assertEqual(round(st.w, 0), 1203)
        self.assertEqual(round(st.mu.muPas, 0), 136)
        self.assertEqual(round(st.k, 3), 0.152)

        st = C2(T=420, P=4e5, **kw)
        self.assertEqual(round(st.rho, 3), 3.478)
        self.assertEqual(round(st.rhoM, 4), 0.1157)
        self.assertEqual(round(st.uM.kJkmol, -1), 15630)
        self.assertEqual(round(st.hM.kJkmol, -1), 19090)
        self.assertEqual(round(st.sM.kJkmolK, 1), 238.0)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 59.86)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 68.61)
        self.assertEqual(round(st.w, 1), 361.2)
        self.assertEqual(round(st.mu.muPas, 1), 12.8)
        self.assertEqual(round(st.k, 4), 0.0383)

        st = C2(T=200, P=5e5, **kw)
        self.assertEqual(round(st.rho, 1), 524.5)
        self.assertEqual(round(st.rhoM, 2), 17.44)
        self.assertEqual(round(st.uM.kJkmol, 0), -7113)
        self.assertEqual(round(st.hM.kJkmol, 0), -7084)
        self.assertEqual(round(st.sM.kJkmolK, 1), 132.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 42.86)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 75.09)
        self.assertEqual(round(st.w, 0), 1205)
        self.assertEqual(round(st.mu.muPas, 0), 136)
        self.assertEqual(round(st.k, 3), 0.152)

        st = C2(T=230, P=6e5, rho0=1, **kw)
        self.assertEqual(round(st.rho, 2), 10.61)
        self.assertEqual(round(st.rhoM, 4), 0.3528)
        self.assertEqual(round(st.uM.kJkmol, 0), 6212)
        self.assertEqual(round(st.hM.kJkmol, 0), 7912)
        self.assertEqual(round(st.sM.kJkmolK, 1), 199.8)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 38.94)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 52.39)
        self.assertEqual(round(st.w, 1), 258.2)
        self.assertEqual(round(st.mu.muPas, 2), 7.25)
        self.assertEqual(round(st.k, 4), 0.0143)

        st = C2(T=360, P=8e5, **kw)
        self.assertEqual(round(st.rho, 3), 8.316)
        self.assertEqual(round(st.rhoM, 4), 0.2766)
        self.assertEqual(round(st.uM.kJkmol, -1), 12130)
        self.assertEqual(round(st.hM.kJkmol, -1), 15030)
        self.assertEqual(round(st.sM.kJkmolK, 1), 221.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 52.39)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 62.08)
        self.assertEqual(round(st.w, 1), 331.8)
        self.assertEqual(round(st.mu.muPas, 1), 11.3)
        self.assertEqual(round(st.k, 4), 0.0295)

        st = C2(T=245, P=1e6, rho0=1, **kw)
        self.assertEqual(round(st.rho, 2), 17.56)
        self.assertEqual(round(st.rhoM, 4), 0.5838)
        self.assertEqual(round(st.uM.kJkmol, 0), 6532)
        self.assertEqual(round(st.hM.kJkmol, 0), 8245)
        self.assertEqual(round(st.sM.kJkmolK, 1), 197.6)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 41.37)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 58.26)
        self.assertEqual(round(st.w, 1), 255.6)
        self.assertEqual(round(st.mu.muPas, 2), 7.78)
        self.assertEqual(round(st.k, 4), 0.0162)

        st = C2(T=260, P=2e6, **kw)
        self.assertEqual(round(st.rho, 1), 430.1)
        self.assertEqual(round(st.rhoM, 2), 14.30)
        self.assertEqual(round(st.uM.kJkmol, 0), -2248)
        self.assertEqual(round(st.hM.kJkmol, 0), -2108)
        self.assertEqual(round(st.sM.kJkmolK, 1), 154.1)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 46.38)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 95.30)
        self.assertEqual(round(st.w, 1), 713.4)
        self.assertEqual(round(st.mu.muPas, 1), 68.7)
        self.assertEqual(round(st.k, 3), 0.101)

        st = C2(T=150, P=3e6, **kw)
        self.assertEqual(round(st.rho, 1), 587.5)
        self.assertEqual(round(st.rhoM, 2), 19.54)
        self.assertEqual(round(st.uM.kJkmol, -1), -10760)
        self.assertEqual(round(st.hM.kJkmol, -1), -10610)
        self.assertEqual(round(st.sM.kJkmolK, 1), 111.8)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 43.54)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 70.18)
        self.assertEqual(round(st.w, 0), 1602)
        self.assertEqual(round(st.mu.muPas, 0), 258)
        self.assertEqual(round(st.k, 3), 0.204)

        st = C2(T=300, P=4e6, **kw)
        self.assertEqual(round(st.rho, 2), 84.99)
        self.assertEqual(round(st.rhoM, 3), 2.827)
        self.assertEqual(round(st.uM.kJkmol, 0), 6868)
        self.assertEqual(round(st.hM.kJkmol, 0), 8283)
        self.assertEqual(round(st.sM.kJkmolK, 1), 189.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 53.64)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 144.3)
        self.assertEqual(round(st.w, 1), 222.9)
        self.assertEqual(round(st.mu.muPas, 1), 11.8)
        self.assertEqual(round(st.k, 4), 0.0308)

        st = C2(T=150, P=5e6, **kw)
        self.assertEqual(round(st.rho, 1), 588.7)
        self.assertEqual(round(st.rhoM, 2), 19.58)
        self.assertEqual(round(st.uM.kJkmol, -1), -10790)
        self.assertEqual(round(st.hM.kJkmol, -1), -10530)
        self.assertEqual(round(st.sM.kJkmolK, 1), 111.6)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 43.66)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 70.02)
        self.assertEqual(round(st.w, 0), 1612)
        self.assertEqual(round(st.mu.muPas, 0), 262)
        self.assertEqual(round(st.k, 3), 0.205)

        st = C2(T=311, P=6e6, **kw)
        self.assertEqual(round(st.rho, 1), 280.8)
        self.assertEqual(round(st.rhoM, 3), 9.339)
        self.assertEqual(round(st.uM.kJkmol, 0), 3265)
        self.assertEqual(round(st.hM.kJkmol, 0), 3908)
        self.assertEqual(round(st.sM.kJkmolK, 1), 173.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 55.44)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 255.5)
        self.assertEqual(round(st.w, 1), 280.7)
        self.assertEqual(round(st.mu.muPas, 1), 30.4)
        self.assertEqual(round(st.k, 4), 0.0699)

        st = C2(T=580, P=7e6, **kw)
        self.assertEqual(round(st.rho, 2), 45.23)
        self.assertEqual(round(st.rhoM, 3), 1.504)
        self.assertEqual(round(st.uM.kJkmol, -1), 25890)
        self.assertEqual(round(st.hM.kJkmol, -1), 30550)
        self.assertEqual(round(st.sM.kJkmolK, 1), 237.8)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 79.51)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 91.35)
        self.assertEqual(round(st.w, 1), 416.2)
        self.assertEqual(round(st.mu.muPas, 1), 19.2)
        self.assertEqual(round(st.k, 4), 0.0696)

        st = C2(T=350, P=8e6, **kw)
        self.assertEqual(round(st.rho, 1), 149.3)
        self.assertEqual(round(st.rhoM, 3), 4.966)
        self.assertEqual(round(st.uM.kJkmol, 0), 8078)
        self.assertEqual(round(st.hM.kJkmol, 0), 9689)
        self.assertEqual(round(st.sM.kJkmolK, 1), 190.6)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 58.59)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 138.3)
        self.assertEqual(round(st.w, 1), 259.1)
        self.assertEqual(round(st.mu.muPas, 1), 18.1)
        self.assertEqual(round(st.k, 4), 0.0513)

        st = C2(T=150, P=1e7, **kw)
        self.assertEqual(round(st.rho, 1), 591.7)
        self.assertEqual(round(st.rhoM, 2), 19.68)
        self.assertEqual(round(st.uM.kJkmol, -1), -10860)
        self.assertEqual(round(st.hM.kJkmol, -1), -10350)
        self.assertEqual(round(st.sM.kJkmolK, 1), 111.2)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 43.93)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 69.65)
        self.assertEqual(round(st.w, 0), 1637)
        self.assertEqual(round(st.mu.muPas, 0), 273)
        self.assertEqual(round(st.k, 3), 0.208)

        st = C2(T=530, P=2e7, **kw)
        self.assertEqual(round(st.rho, 1), 149.4)
        self.assertEqual(round(st.rhoM, 3), 4.970)
        self.assertEqual(round(st.uM.kJkmol, -1), 19960)
        self.assertEqual(round(st.hM.kJkmol, -1), 23980)
        self.assertEqual(round(st.sM.kJkmolK, 1), 217.7)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 75.23)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 96.67)
        self.assertEqual(round(st.w, 1), 426.0)
        self.assertEqual(round(st.mu.muPas, 1), 25.3)
        self.assertEqual(round(st.k, 4), 0.0779)

        st = C2(T=200, P=4e7, **kw)
        self.assertEqual(round(st.rho, 1), 560.2)
        self.assertEqual(round(st.rhoM, 2), 18.63)
        self.assertEqual(round(st.uM.kJkmol, 0), -7940)
        self.assertEqual(round(st.hM.kJkmol, 0), -5793)
        self.assertEqual(round(st.sM.kJkmolK, 1), 128.4)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 44.75)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 69.62)
        self.assertEqual(round(st.w, 0), 1488)
        self.assertEqual(round(st.mu.muPas, 0), 187)
        self.assertEqual(round(st.k, 3), 0.184)

        st = C2(T=600, P=6e7, **kw)
        self.assertEqual(round(st.rho, 1), 281.9)
        self.assertEqual(round(st.rhoM, 3), 9.374)
        self.assertEqual(round(st.uM.kJkmol, -1), 22940)
        self.assertEqual(round(st.hM.kJkmol, -1), 29340)
        self.assertEqual(round(st.sM.kJkmolK, 1), 217.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 83.63)
        self.assertEqual(round(st.cpM.kJkmolK, 0), 103)
        self.assertEqual(round(st.w, 1), 702.9)
        self.assertEqual(round(st.mu.muPas, 1), 41.7)
        self.assertEqual(round(st.k, 3), 0.115)

    def test_friend(self):
        # Selected point from Table A1, Pag 336, ideal gas
        # The point really are not the gas ideal at zero preesure as it's
        # calculated here, it's the value at 1 bar pressure so the values of
        # Helmholtz free energy and entropy differ of ideal properties

        st = C2(T=100, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.hM.kJmol, 3), 3.384)
        self.assertEqual(round(st.cpM.JmolK, 3), 35.698)
        self.assertEqual(round(st.mu.muPas, 2), 3.32)
        self.assertEqual(round(st.k.mWmK, 2), 3.46)

        st = C2(T=500, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.hM.kJmol, 3), 25.059)
        self.assertEqual(round(st.cpM.JmolK, 3), 77.987)
        self.assertEqual(round(st.mu.muPas, 2), 14.76)
        self.assertEqual(round(st.k.mWmK, 2), 53.78)

        # Selected point from Table A2, Pag 337, saturation state
        # This table use ancillary equation for calculate presure and densities
        # of two phases so the result may differ

        # TODO: Add heat capacity along the saturated boundary to calculated
        # propeties, the table show this value
        st = C2(T=100, x=0.5, eq="friend", visco=2)
        self.assertEqual(round(st.P.MPa, 6), 0.11e-4)
        self.assertEqual(round(st.Liquido.rhoM, 2), 21.32)
        self.assertEqual(round(st.Gas.rhoM, 7), 0.133e-4)
        self.assertEqual(round(st.Liquido.w, 1), 1938.3)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 877.99)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 248.2)

        st = C2(T=200, x=0.5, eq="friend", visco=2)
        self.assertEqual(round(st.P.MPa, 3), 0.217)
        self.assertEqual(round(st.Liquido.rhoM, 2), 17.42)
        self.assertEqual(round(st.Gas.rhoM, 3), 0.139)
        self.assertEqual(round(st.Liquido.w, 1), 1194.8)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 138.23)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 152.6)

        st = C2(T=300, x=0.5, eq="friend", visco=2)
        self.assertEqual(round(st.P.MPa, 3), 4.357)
        self.assertEqual(round(st.Liquido.rhoM, 2), 10.08)
        self.assertEqual(round(st.Gas.rhoM, 3), 3.808)
        self.assertEqual(round(st.Liquido.w, 1), 276.8)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 34.92)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 71.4)

        # Selected point from Table A3, Pag 339, single phase region
        st = C2(T=100, P=1e5, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 21.33)
        self.assertEqual(round(st.hM.kJmol, 3), -14.221)
        self.assertEqual(round(st.sM.JmolK, 2), 83.60)
        self.assertEqual(round(st.cvM.JmolK, 2), 48.15)
        self.assertEqual(round(st.cpM.JmolK, 2), 70.11)
        self.assertEqual(round(st.w, 1), 1938.7)
        self.assertEqual(round(st.mu.muPas, 2), 878.69)
        self.assertEqual(round(st.k.mWmK, 1), 248.2)

        st = C2(T=130, P=1e6, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 20.24)
        self.assertEqual(round(st.hM.kJmol, 3), -12.071)
        self.assertEqual(round(st.sM.JmolK, 2), 102.03)
        self.assertEqual(round(st.cvM.JmolK, 2), 45.01)
        self.assertEqual(round(st.cpM.JmolK, 2), 70.10)
        self.assertEqual(round(st.w, 1), 1726.9)
        self.assertEqual(round(st.mu.muPas, 2), 392.40)
        self.assertEqual(round(st.k.mWmK, 1), 221.3)

        st = C2(T=160, P=5e7, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 20.06)
        self.assertEqual(round(st.hM.kJmol, 3), -8.163)
        self.assertEqual(round(st.sM.JmolK, 2), 112.25)
        self.assertEqual(round(st.cvM.JmolK, 2), 44.94)
        self.assertEqual(round(st.cpM.JmolK, 2), 67.54)
        self.assertEqual(round(st.w, 1), 1763.8)
        self.assertEqual(round(st.mu.muPas, 2), 323.48)
        self.assertEqual(round(st.k.mWmK, 1), 223.9)

        st = C2(T=200, P=1e7, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 17.79)
        self.assertEqual(round(st.hM.kJmol, 3), -6.804)
        self.assertEqual(round(st.sM.JmolK, 2), 131.51)
        self.assertEqual(round(st.cvM.JmolK, 2), 43.41)
        self.assertEqual(round(st.cpM.JmolK, 2), 73.00)
        self.assertEqual(round(st.w, 1), 1281.7)
        self.assertEqual(round(st.mu.muPas, 2), 151.38)
        self.assertEqual(round(st.k.mWmK, 1), 161.5)

        st = C2(T=250, P=5e5, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 0.26)
        self.assertEqual(round(st.hM.kJmol, 3), 9.048)
        self.assertEqual(round(st.sM.JmolK, 2), 205.92)
        self.assertEqual(round(st.cvM.JmolK, 2), 39.92)
        self.assertEqual(round(st.cpM.JmolK, 2), 51.02)
        self.assertEqual(round(st.w, 1), 276.3)
        self.assertEqual(round(st.mu.muPas, 2), 8.01)
        self.assertEqual(round(st.k.mWmK, 1), 16.0)

        st = C2(T=300, P=3e7, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 14.78)
        self.assertEqual(round(st.hM.kJmol, 3), 1.409)
        self.assertEqual(round(st.sM.JmolK, 2), 159.75)
        self.assertEqual(round(st.cvM.JmolK, 2), 50.85)
        self.assertEqual(round(st.cpM.JmolK, 2), 81.83)
        self.assertEqual(round(st.w, 1), 905.4)
        self.assertEqual(round(st.mu.muPas, 2), 77.56)
        self.assertEqual(round(st.k.mWmK, 1), 112.6)

        st = C2(T=350, P=1e5, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 0.03)
        self.assertEqual(round(st.hM.kJmol, 3), 14.720)
        self.assertEqual(round(st.sM.JmolK, 2), 238.06)
        self.assertEqual(round(st.cvM.JmolK, 2), 50.76)
        self.assertEqual(round(st.cpM.JmolK, 2), 59.24)
        self.assertEqual(round(st.w, 1), 334.6)
        self.assertEqual(round(st.mu.muPas, 2), 10.84)
        self.assertEqual(round(st.k.mWmK, 1), 28.1)

        st = C2(T=400, P=5e6, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 1.77)
        self.assertEqual(round(st.hM.kJmol, 3), 16.051)
        self.assertEqual(round(st.sM.JmolK, 2), 210.58)
        self.assertEqual(round(st.cvM.JmolK, 2), 59.05)
        self.assertEqual(round(st.cpM.JmolK, 2), 76.57)
        self.assertEqual(round(st.w, 1), 322.4)
        self.assertEqual(round(st.mu.muPas, 2), 13.91)
        self.assertEqual(round(st.k.mWmK, 1), 40.0)

        st = C2(T=450, P=1e5, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 0.03)
        self.assertEqual(round(st.hM.kJmol, 3), 21.285)
        self.assertEqual(round(st.sM.JmolK, 2), 254.49)
        self.assertEqual(round(st.cvM.JmolK, 2), 63.59)
        self.assertEqual(round(st.cpM.JmolK, 2), 71.99)
        self.assertEqual(round(st.w, 1), 374.6)
        self.assertEqual(round(st.mu.muPas, 2), 13.52)
        self.assertEqual(round(st.k.mWmK, 1), 44.6)

        st = C2(T=500, P=1e7, eq="friend", visco=2)
        self.assertEqual(round(st.rhoM, 2), 2.68)
        self.assertEqual(round(st.hM.kJmol, 3), 22.854)
        self.assertEqual(round(st.sM.JmolK, 2), 220.68)
        self.assertEqual(round(st.cvM.JmolK, 2), 71.12)
        self.assertEqual(round(st.cpM.JmolK, 2), 88.12)
        self.assertEqual(round(st.w, 1), 378.1)
        self.assertEqual(round(st.mu.muPas, 2), 17.95)
        self.assertEqual(round(st.k.mWmK, 1), 59.6)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = C2(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 3), 3.299)
        self.assertEqual(round(st.P.MPa, 3), 44.781)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.6276)

        st2 = C2(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 209.07)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.50714)

    def test_Herrmann(self):
        # Table 3, Pag 5
        self.assertEqual(round(C2(T=100, rho=650).mu.muPas, 3), 1050.774)
        self.assertEqual(round(C2(T=300, rho=0).mu.muPas, 6), 9.285785)
        self.assertEqual(round(C2(T=300, rho=1).mu.muPas, 6), 9.293015)
        self.assertEqual(round(C2(T=300, rho=100).mu.muPas, 5), 12.55151)
        self.assertEqual(round(C2(T=300, rho=500).mu.muPas, 4), 114.3985)
        self.assertEqual(round(C2(T=C2.Tc, rho=206.18).mu.muPas, 5), 22.63285)
        self.assertEqual(round(C2(T=310, rho=1).mu.muPas, 6), 9.592824)
        self.assertEqual(round(C2(T=310, rho=100).mu.muPas, 5), 12.88843)
        self.assertEqual(round(C2(T=310, rho=500).mu.muPas, 4), 114.8683)
        self.assertEqual(round(C2(T=500, rho=1).mu.muPas, 5), 14.85686)
        self.assertEqual(round(C2(T=500, rho=100).mu.muPas, 5), 18.82032)
        self.assertEqual(round(C2(T=500, rho=400).mu.muPas, 5), 66.15890)
        self.assertEqual(round(C2(T=675, rho=0).mu.muPas, 5), 18.97242)
        self.assertEqual(round(C2(T=675, rho=1).mu.muPas, 5), 18.99427)
        self.assertEqual(round(C2(T=675, rho=100).mu.muPas, 5), 23.37711)
        self.assertEqual(round(C2(T=675, rho=300).mu.muPas, 5), 45.90524)

    def test_Vogel(self):
        # Table 6, Pag 22
        self.assertEqual(round(
            C2(T=100, rho=650, visco=1).mu.muPas, 9), 1058.240219296)
        self.assertEqual(round(
            C2(T=300, rho=1, visco=1).mu.muPas, 9), 9.286370770)
        self.assertEqual(round(
            C2(T=300, rho=100, visco=1).mu.muPas, 9), 12.529249197)
        self.assertEqual(round(
            C2(T=300, rho=500, visco=1).mu.muPas, 9), 113.730164659)
        self.assertEqual(round(
            C2(T=305.322, rho=206.18, visco=1).mu.muPas, 9), 22.630721989)
        self.assertEqual(round(
            C2(T=310, rho=1, visco=1).mu.muPas, 9), 9.581087691)
        self.assertEqual(round(
            C2(T=310, rho=100, visco=1).mu.muPas, 9), 12.873112150)
        self.assertEqual(round(
            C2(T=310, rho=500, visco=1).mu.muPas, 9), 114.131039291)
        self.assertEqual(round(
            C2(T=500, rho=1, visco=1).mu.muPas, 9), 14.834335193)
        self.assertEqual(round(
            C2(T=500, rho=100, visco=1).mu.muPas, 9), 18.900768605)
        self.assertEqual(round(
            C2(T=500, rho=400, visco=1).mu.muPas, 9), 66.704548052)
        self.assertEqual(round(
            C2(T=675, rho=1, visco=1).mu.muPas, 9), 18.902546759)
        self.assertEqual(round(
            C2(T=675, rho=100, visco=1).mu.muPas, 9), 23.421278317)
        self.assertEqual(round(
            C2(T=675, rho=300, visco=1).mu.muPas, 9), 45.895671715)

    def test_friendThermo(self):
        # Selected point from Table A1, Pag 336, ideal gas
        st = C2(T=100, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 3.32)
        self.assertEqual(round(st.k.mWmK, 2), 3.46)

        st = C2(T=200, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 6.35)
        self.assertEqual(round(st.k.mWmK, 2), 10.49)

        st = C2(T=300, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 9.39)
        self.assertEqual(round(st.k.mWmK, 2), 21.13)

        st = C2(T=400, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 12.19)
        self.assertEqual(round(st.k.mWmK, 2), 35.95)

        st = C2(T=500, rho=0, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 14.76)
        self.assertEqual(round(st.k.mWmK, 2), 53.78)

        # Selected point from Table A2, Pag 337, saturation state
        # This table has tiny desviation in saturation calculation
        # st = C2(T=304, x=0.5, eq="friend", visco=2)
        # self.assertEqual(round(st.Liquido.mu.muPas, 2), 28.97)
        # self.assertEqual(round(st.Liquido.k.mWmK, 1), 79.0)

        # Selected point from Table A3, Pag 339, single phase region
        st = C2(T=100, P=1e5, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 878.69)
        self.assertEqual(round(st.k.mWmK, 1), 248.2)

        st = C2(T=170, P=6e7, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 301.71)
        self.assertEqual(round(st.k.mWmK, 1), 221.8)

        st = C2(T=260, P=5e6, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 74.61)
        self.assertEqual(round(st.k.mWmK, 1), 106.5)

        st = C2(T=330, P=5e5, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 10.37)
        self.assertEqual(round(st.k.mWmK, 1), 25.6)

        st = C2(T=380, P=1e6, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 11.88)
        self.assertEqual(round(st.k.mWmK, 1), 33.3)

        st = C2(T=420, P=4e7, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 46.86)
        self.assertEqual(round(st.k.mWmK, 1), 89.7)

        st = C2(T=480, P=1e5, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 14.28)
        self.assertEqual(round(st.k.mWmK, 1), 50.1)

        st = C2(T=500, P=6e7, eq="friend", visco=2)
        self.assertEqual(round(st.mu.muPas, 2), 48.34)
        self.assertEqual(round(st.k.mWmK, 1), 101.4)
