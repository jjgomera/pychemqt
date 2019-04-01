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

from lib import unidades
from lib.meos import MEoS


class iC4(MEoS):
    """Multiparameter equation of state for isobutane"""
    name = "isobutane"
    CASNumber = "75-28-5"
    formula = "CH(CH3)3"
    synonym = "R-600a"
    _refPropName = "ISOBUTAN"
    _coolPropName = "IsoButane"
    rhoc = unidades.Density(225.5)
    Tc = unidades.Temperature(407.81)
    Pc = unidades.Pressure(3629.0, "kPa")
    M = 58.1222  # g/mol
    Tt = unidades.Temperature(113.73)
    Tb = unidades.Temperature(261.401)
    f_acent = 0.184
    momentoDipolar = unidades.DipoleMoment(0.132, "Debye")
    id = 5
    _Tr = unidades.Temperature(390.355535)
    _rhor = unidades.Density(228.302484)
    _w = 0.178714317

    Fi1 = {"ao_log": [1, 3.05956619],
           "pow": [0, 1],
           "ao_pow": [11.60865546, -5.29450411],
           "ao_exp": [4.94641014, 4.09475197, 15.6632824, 9.73918122],
           "titao": [0.9512779015, 2.3878958853, 4.3469042691, 10.3688586351]}

    Fi2 = {"ao_log": [1, 3.06714],
           "pow": [0, 1],
           "ao_pow": [20.413726078, -94.467620036],
           "ao_exp": [], "titao": [],
           "ao_sinh": [8.97575, 25.1423], "sinh": [1.074673199, 4.671261865],
           "ao_cosh": [5.25156, 16.1388], "cosh": [0.485556021, 2.19158348]}

    Fi3 = {"ao_log": [1, 3.059347],
           "ao_pow": [-6.026745, -5.035251],
           "pow": [0, 1],
           "ao_exp": [4.940314, 4.090139, 9.739581, 15.68832],
           "titao": [0.9508183, 2.383449, 10.38655, 4.347095]}

    CP4 = {"ao": -1.7231723278e1,
           "an": [1.7027919006e7, -4.7269724737e5, 4.7301406581e3,
                  5.8491344291e-2, 8.9440351886e-6, -1.8274599197e-8],
           "pow": [-3, -2, -1, 1, 2, 3],
           "ao_exp": [-1.9283021962e1], "exp": [3000]}

    CP5 = {"ao": 4.06714,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_sinh": [8.97575, 25.1423], "sinh": [438.27, 1905.02],
           "ao_cosh": [5.25156, 16.1388], "cosh": [198.018, 893.765]}

    CP6 = {"ao": 0.397893/8.3143*58.124,
           "an": [0.412501e-2/8.3143*58.124, -0.196195e-6/8.3143*58.124,
                  0.380185e-8/8.3143*58.124, -0.523950e-11/8.3143*58.124],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": []}

    buecker = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Buecker and "
                    "Wagner (2006)",
        "__doi__": {"autor": "Bücker, D., Wagner, W.",
                    "title": "Reference Equations of State for the "
                             "Thermodynamic Properties of Fluid Phase "
                             "n-Butane and Isobutane",
                    "ref": "J. Phys. Chem. Ref. Data 35(2) (2006) 929-1019",
                    "doi": "10.1063/1.1901687"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 35000.0, "rhomax": 12.9,
        "Pmin": 0.0000219, "rhomin": 12.74,

        "nr1":  [0.20686820727966e1, -0.36400098615204e1, 0.51968754427244,
                 0.17745845870123, -0.12361807851599, 0.45145314010528e-1,
                 0.30476479965980e-1],
        "d1": [1, 1, 1, 2, 3, 4, 4],
        "t1": [0.50, 1.00, 1.50, 0.00, 0.50, 0.50, 0.75],
        "nr2": [0.75508387706302, -0.85885381015629, 0.36324009830684e-1,
                -0.01954879945055, -0.44452392904960e-2, 0.46410763666460e-2,
                -0.71444097992825e-1, -0.80765060030713e-1, 0.15560460945053,
                0.20318752160332e-2, -0.10624883571689, 0.39807690546305e-1,
                0.16371431292386e-1, 0.53212200682628e-3, -0.78681561156387e-2,
                -0.30981191888963e-2],
        "d2": [1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6],
        "t2": [2.00, 2.50, 2.50, 1.50, 1.00, 1.50, 4.00, 7.00, 3.00, 7.00,
               3.00, 1.00, 6.00, 0.00, 6.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*16,

        "nr3": [-0.42276036810382e-1, -0.53001044558079e-2],
        "d3": [1, 2],
        "t3": [2., 0.],
        "alfa3": [10, 10],
        "beta3": [150, 200],
        "gamma3": [1.16, 1.13],
        "epsilon3": [0.85, 1.]}

    younglove = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for isobutane of Younglove and "
                    "Ely (1987)",
        "__doi__": {"autor": "Younglove, B.A. and Ely, J.F.",
                    "title": "Thermophysical Properties of Fluids. II. "
                             "Methane, Ethane, Propane, Isobutane, and Normal "
                             "Butane",
                    "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                    "doi": "10.1063/1.555785"},

        "R": 8.31434,
        "M": 58.125, "Tt": 113.55, "Tc": 407.85, "Pc": 3640, "rhoc": 3.86,

        "cp": CP4,
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 17896.45, "so": 295.48},

        "Tmin": 113.55, "Tmax": 600.0, "Pmax": 35000.0, "rhomax": 12.89,
        "Pmin": 1.948e-5, "rhomin": 12.755,

        "b": [None, 0.1307325972e-1, 0.3927802742, -0.3185427394e2,
              0.7608825192e4, -0.1753919859e7, -0.2090019755e-2, 8.959557971,
              -0.6816710130e4, -0.1111271045e7, 0.3248737572e-3, -1.046526456,
              0.6536598969e3, 0.3726503734e-1, 0.8553649395e1, 0.2109987236e4,
              -0.1401267363e1, 0.5213089327e-1, -0.1925026382e2, 0.7640067895,
              0.3425854273e7, -0.3373475924e9, 0.1180683444e6, 0.1529683738e10,
              0.3323837416e4, 0.6423169487e5, 0.3891706042e2, -0.1494755736e7,
              -0.1720240173e-1, 0.2894195375e3, 0.2005086329e-2, -0.4448393005,
              0.8028488415e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 35000.0, "rhomax": 12.9,
        "Pmin": 7.36, "rhomin": 38.2,

        "nr1":  [0.10429331589100e1, -0.28184272548892e1, 0.86176232397850,
                 -0.10613619452487, 0.98615749302134e-1, 0.23948208682322e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30330004856950, -0.041598156135099, -0.29991937470058,
                -0.080369342764109, -0.029761373251151, 0.013059630303140],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    miyamoto = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Miyamoto "
                    "and Watanabe (2001)",
        "__doi__": {"autor": "Miyamoto, H. and Watanabe, K.",
                    "title": "A Thermodynamic Property Model for Fluid-Phase "
                             "Isobutane",
                    "ref": "Int. J. Thermophys., 23(2) (2002) 477-499",
                    "doi": "10.1023/A:1015161519954"},

        "R": 8.314472,
        "M": 58.1222, "Tc": 407.817, "rhoc": 224.36/58.1222,
        "cp": Fi3,
        "ref": "IIR",

        "Tmin": 113.56, "Tmax": 573.0, "Pmax": 35000.0, "rhomax": 12.9,
        "Pmin": 0.000021, "rhomin": 12.738,

        "nr1":  [2.892737e-1, -1.342570, -7.976713e-3, 2.025793e-1,
                 -4.241612e-2, 2.617971e-3, 5.068955e-5, -1.144596e-6],
        "d1": [1, 1, 2, 2, 3, 5, 8, 8],
        "t1": [-0.25, 1.5, -0.75, 0, 1.25, 1.5, 0.5, 2.5],

        "nr2": [-1.930153, 1.982609, 2.076533e-3, -4.958752e-3, 1.377372e-3,
                -1.582662e-1, -4.961892e-2, 9.451030e-4, -3.037276e-2,
                -1.382675e-2, 8.876254e-5],
        "d2": [3, 3, 8, 5, 6, 1, 5, 7, 2, 3, 15],
        "t2": [1.5, 1.75, -0.25, 3, 3, 4, 2, -1, 2, 19, 5],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*11}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isobutane of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP5,
        "ref": "OTO",
        "M": 58.123, "Tc": 407.817, "rhoc": 224.36/58.123,

        "Tmin": 113.55, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 12.89,
        "Pmin": 0.000020860, "rhomin": 12.784,

        "nr1":  [0.10429332e1, -0.28184273e1, 0.86176232, -0.10613619,
                 0.986157490e-1, 0.23948209e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30330005, -0.41598156e-1, -0.29991937, -0.80369343e-1,
                -0.29761373e-1, 0.1305963e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP6,
        "ref": "NBP",

        "Tmin": 120.0, "Tmax": 498.0, "Pmax": 35000.0, "rhomax": 12.89,
        "Pmin": 0.46491e-4, "rhomin": 12.649,

        "nr1":  [-0.958589873652, 0.818846326211, -0.115814967179,
                 0.345513148715, -0.168751721524e1, 0.936693300209,
                 -0.106644545724e1, 0.980958295776e-1, 0.495941129005,
                 -0.261313404262, 0.485109471188, -0.177275820736,
                 -0.209415485311e-1, 0.788178884079e-1, -0.102751671767,
                 0.178645875838e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.958589873652, -0.818846326211, 0.115814967179,
                0.537585249054, -0.71942446879, 0.245830118086],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.0071072]*6}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Sun and Ely "
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

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.18083775, 9.46903331e-1, -2.90618044, 8.51346220e-2,
                2.79868503e-4, -1.68266335e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-2.01202825e-1, -3.32570120e-2, 2.42967225e-1, -4.20931100e-3,
                -0.224528572, -1.41307663e-2, -5.93401702e-2, -2.27862942e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = buecker, younglove, GERG, miyamoto, shortSpan, polt, sun

    _surface = {"sigma": [-0.01639, 0.06121], "exp": [2.102, 1.304]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.388417],  "expt0": [-1.], "expd0": [1.],
                   "a1": [20.534, 0.02], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [126.25, 52.91, -7501.4, -2672.9],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.9, 2.9]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.000022891,
                "Tmin": Tt, "Tmax": 575.0,
                "a1": [-1953637129., 1953637130.], "exp1": [0, 6.12],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.85093103, 1.36543198, -1.32542691, -2.56190994],
        "t": [1, 1.5, 2.5, 4.5]}
    _liquid_Density = {
        "eq": 1,
        "n": [2.04025104, 0.850874089, -0.479052281, 0.348201252],
        "t": [0.355, 1, 4/3, 7/3]}
    _vapor_Density = {
        "eq": 3,
        "n": [-2.12933323, -2.93790085, -0.89441086, -3.46343707],
        "t": [0.355, 5/6, 19/6, 26/6]}

    visco0 = {"__name__": "Herrmann (2018)",
              "__doi__": {
                  "autor": "Herrmann, S., Vogel, E.",
                  "title": "New Formulation for the Viscosity of Isobutane",
                  "ref": "J. Phys. Chem. Ref. Data 47(4) (2018) 043103",
                  "doi": "10.1063/1.5057413"},

              "eq": 1, "omega": 0,

              "special0": "_mu0",

              "Tref_virial": 407.81,
              # Special term of virial coefficient, with δ term and μPa·s
              "muref_virial": 2.961503647e-1/225.5*M/1e3,
              "n_virial": [-1.9572881e1, 2.035770736e2, -8.714522319e2,
                           1.964871902e3, -2.486422559e3, 1.700542778e3,
                           -4.9777979e2, 6.560977477, -6.455286496e-2],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 407.81, "rhoref_res": 225.5,
              "nr": [2.871640981e1, -1.636506531e1, 7.67015053, -5.790243964e1,
                     3.891832708e1, 2.730660102e1, -1.658899499e1,
                     2.246568982e-1, -2.343800854e-1, 4.328758428e-2],
              "tr": [0, 3, 4, 0, 1, 0, 1, 4, 4, 2],
              "dr": [2, 2, 2, 2.5, 2.5, 3, 3, 3, 5, 8],

              "nr_gaus": [2.75e-1, 5.5e-1],
              "br_gaus": [5, 50],
              "er_gaus": [30, 100],

              "special": "_mur"}

    def _mu0(self):
        """Special term for zero-density viscosity for Herrmann correlation"""
        Tc = 407.81
        tau = Tc/self.T

        npf = 1.032954077e3
        n00 = 1.6520217e3
        n10 = 1.95822722e1
        n20 = 2.055051549e3

        # Eq 9
        muo = npf/tau**0.5/(
            n00*tau**0.5/Tc**0.5 + n10 + n20*exp(-(Tc/tau)**(1/3)))
        return muo

    def _mur(self, rho, T, fase):
        """Special term of residual viscosity for Herrmann correlation"""
        if rho:
            tau = 407.81/T
            delta = rho/225.5
            mur = tau**0.5/delta**(2/3)*1.693411655e-10*(delta**5.6*tau)**4
        else:
            mur = 0
        return mur

    visco1 = {"__name__": "Vogel (2000)",
              "__doi__": {
                  "autor": "Vogel, E., Küchenmeister, C., Bich, E.",
                  "title": "Viscosity Correlation for Isobutane over Wide "
                           "Ranges of the Fluid Region",
                  "ref": "Int. J. Thermophys 21(2) (2000) 343-356",
                  "doi": "10.1023/A:1006623310780"},

              "eq": 1, "omega": 1,
              "ek": 307.55, "sigma": 0.46445,
              "n_chapman": 0.021357,
              "collision": [0.53583008, -0.45629630, 0.049911282],

              "Tref_virial": 307.55,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 407.817, "rhoref_res": 3.86*M, "etaref_res": 1,
              "nr": [103.511763411, -312.670896234, 145.253750239,
                     -210.649894193, 386.269696509, -214.963015527,
                     112.58036092, -223.242033154, 119.114788598, -18.19097459,
                     36.0438957232, -21.3960184050],
              "tr": [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2],
              "dr": [2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5],

              "CPf": 1940.3760699,
              "CPg1": 2.33859774637,
              "CPgi": [1.00596672174],
              "CPti": [-0.5]}

    visco2 = {"__name__": "Younglove (1987)",
              "__doi__": {
                  "autor": "Younglove, B.A., Ely, J.F.",
                  "title": "Thermophysical Properties of Fluids. II. Methane, "
                           "Ethane, Propane, Isobutane, and Normal Butane",
                  "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                  "doi": "10.1063/1.555785"},

              "eq": 2, "omega": 2,

              "ek": 418.0, "sigma": 0.509217,

              "F": [1.687838652, 0.0, 1.40, 407.85],
              "E": [-0.2055498053e2, 0.1357076181e4, 0.1893774336e2,
                    -0.1822277344e5, -0.4599387773e-2, 0.6305247065e2,
                    0.1282253921e5],
              "rhoc": 3.86}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"__name__": "Perkins (2002)",
               "__doi__": {
                   "autor": "Perkins, R.A.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of Isobutane from 114 K to 600 K at "
                            "Pressures to 70 MPa",
                   "ref": "J. Chem. Eng. Data 47(5) (2002) 1272-1279",
                   "doi": "10.1021/je010121u"},

               "eq": 1,

               "Toref": 407.82, "koref": 1,
               "no": [-2.37901e-3, 1.06601e-2, 2.15811e-2],
               "to": [0, 1, 2],

               "Tref_res": 407.82, "rhoref_res": 3.86*M, "kref_res": 1,
               "nr": [-4.11789e-2, 4.76346e-2, 1.46805e-1, -1.28445e-1,
                      -1.19190e-1, 1.07565e-1, 4.10226e-2, -3.85968e-2,
                      -4.88704e-3, 5.20901e-3],
               "tr": [0, -1]*5,
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 0.657661e-9, "Tcref": 611.73}

    thermo1 = {"__name__": "Younglove (1987)",
               "__doi__": {
                   "autor": "Younglove, B.A., Ely, J.F.",
                   "title": "Thermophysical Properties of Fluids. II. Methane,"
                            " Ethane, Propane, Isobutane, and Normal Butane",
                   "ref": "J. Phys. Chem. Ref. Data 16(4) (1987) 577-798",
                   "doi": "10.1063/1.555785"},

               "eq": 3,

               "ek": 418.0,
               "G": [0.1449797353e1, -0.1685643887],
               "E": [0.4307008989e-2, -0.1509010974e1, 0.4693712392e3,
                     -0.3554280979e-3, 0.1841552874, -0.3892338766e2,
                     -0.9354624917e-1, 0.7114330590e1],

               "critical": 2,
               "Tc": 407.85, "rhoc": 3.86*58.125,
               "X": [0.0034718, 10.1207, 0.466392, 1.00344],
               "Z": 9.10218e-10}

    _thermal = thermo0, thermo1


class Test(TestCase):

    def test_buecker(self):
        # Selected point from Table 46, Pag 996, saturation state
        st = iC4(T=114, x=0.5)
        self.assertEqual(round(st.P.MPa, 8), 0.00000002)
        self.assertEqual(round(st.Liquido.rho, 5), 740.08373)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -713.59)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -3.195)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.175)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.690)
        self.assertEqual(round(st.Liquido.w, 2), 1997.42)
        self.assertEqual(round(st.Gas.rho, 7), 0.0000015)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -233.11)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 1.020)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.738)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.881)
        self.assertEqual(round(st.Gas.w, 2), 139.53)

        st = iC4(T=150, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000024)
        self.assertEqual(round(st.Liquido.rho, 3), 706.038)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -650.69)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -2.716)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.251)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.805)
        self.assertEqual(round(st.Liquido.w, 2), 1714.39)
        self.assertEqual(round(st.Gas.rho, 5), 0.00111)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -198.52)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 0.299)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.894)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.037)
        self.assertEqual(round(st.Gas.w, 2), 157.76)

        st = iC4(T=200, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.003814)
        self.assertEqual(round(st.Liquido.rho, 3), 657.706)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -556.44)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -2.175)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.363)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.968)
        self.assertEqual(round(st.Liquido.w, 2), 1389.12)
        self.assertEqual(round(st.Gas.rho, 5), 0.13378)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -142.03)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.103)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.095)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.241)
        self.assertEqual(round(st.Gas.w, 2), 179.38)

        st = iC4(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.37000)
        self.assertEqual(round(st.Liquido.rho, 2), 548.32)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -338.17)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -1.299)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.690)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 2.442)
        self.assertEqual(round(st.Liquido.w, 2), 810.25)
        self.assertEqual(round(st.Gas.rho, 4), 9.6096)
        self.assertEqual(round(st.Gas.h.kJkg, 2), -11.29)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.209)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.578)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.810)
        self.assertEqual(round(st.Gas.w, 2), 197.74)

        st = iC4(T=400, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 3.1856)
        self.assertEqual(round(st.Liquido.rho, 2), 341.03)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -36.27)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -0.459)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 2.250)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 6.349)
        self.assertEqual(round(st.Liquido.w, 2), 184.38)
        self.assertEqual(round(st.Gas.rho, 2), 118.39)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 81.59)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.164)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 2.354)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 7.555)
        self.assertEqual(round(st.Gas.w, 2), 128.90)

        st = iC4(T=407, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 3.5801)
        self.assertEqual(round(st.Liquido.rho, 2), 276.83)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 7.30)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), -0.354)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 2.465)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 42.319)
        self.assertEqual(round(st.Liquido.w, 2), 116.97)
        self.assertEqual(round(st.Gas.rho, 2), 173.46)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 59.64)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), -0.225)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 2.575)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 61.185)
        self.assertEqual(round(st.Gas.w, 2), 113.83)

        # Selected point from Table 47, Pag 1003
        st = iC4(T=115, P=1e5)
        self.assertEqual(round(st.rho, 2), 739.18)
        self.assertEqual(round(st.u.kJkg, 2), -711.91)
        self.assertEqual(round(st.h.kJkg, 2), -711.78)
        self.assertEqual(round(st.s.kJkgK, 4), -3.1801)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.1769)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.6927)
        self.assertEqual(round(st.w, 2), 1988.81)

        st = iC4(T=320, P=5e5)
        self.assertEqual(round(st.rho, 3), 12.322)
        self.assertEqual(round(st.u.kJkg, 3), -20.015)
        self.assertEqual(round(st.h.kJkg, 3), 20.563)
        self.assertEqual(round(st.s.kJkgK, 5), -0.14464)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.6742)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.9171)
        self.assertEqual(round(st.w, 2), 201.39)

        st = iC4(T=330, P=1e6)
        self.assertEqual(round(st.rho, 2), 508.12)
        self.assertEqual(round(st.u.kJkg, 2), -263.45)
        self.assertEqual(round(st.h.kJkg, 2), -261.48)
        self.assertEqual(round(st.s.kJkgK, 4), -1.0589)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.8167)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.6701)
        self.assertEqual(round(st.w, 2), 644.33)

        st = iC4(T=575, P=2e6)
        self.assertEqual(round(st.rho, 3), 25.769)
        self.assertEqual(round(st.u.kJkg, 2), 537.11)
        self.assertEqual(round(st.h.kJkg, 2), 614.72)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0085)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.7190)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.9159)
        self.assertEqual(round(st.w, 2), 280.39)

        st = iC4(T=190, P=5e6)
        self.assertEqual(round(st.rho, 2), 670.90)
        self.assertEqual(round(st.u.kJkg, 2), -577.98)
        self.assertEqual(round(st.h.kJkg, 2), -570.53)
        self.assertEqual(round(st.s.kJkgK, 4), -2.2856)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.3443)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.9257)
        self.assertEqual(round(st.w, 2), 1477.48)

        st = iC4(T=120, P=1e7)
        self.assertEqual(round(st.rho, 2), 738.15)
        self.assertEqual(round(st.u.kJkg, 2), -705.41)
        self.assertEqual(round(st.h.kJkg, 2), -691.86)
        self.assertEqual(round(st.s.kJkgK, 4), -3.1247)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.2021)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.7037)
        self.assertEqual(round(st.w, 2), 1969.38)

        st = iC4(T=230, P=3e7)
        self.assertEqual(round(st.rho, 2), 652.07)
        self.assertEqual(round(st.u.kJkg, 2), -510.55)
        self.assertEqual(round(st.h.kJkg, 2), -464.54)
        self.assertEqual(round(st.s.kJkgK, 4), -1.9602)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.4740)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.0224)
        self.assertEqual(round(st.w, 2), 1398.87)

    def test_younglove(self):
        # The saturation point use the ancillary equation for calculate
        # pressure and density, so the values differ of values give by mBWR,
        # so not used for testing
        kw = {"eq": "younglove", "visco": 2, "thermal": 1}

        # Selected point from Appendix G, Pag 688, single phase region
        st = iC4(T=120, P=1e4, **kw)
        self.assertEqual(round(st.rho, 1), 735.1)
        self.assertEqual(round(st.rhoM, 2), 12.65)
        self.assertEqual(round(st.uM.kJkmol, -1), -23160)
        self.assertEqual(round(st.hM.kJkmol, -1), -23160)
        self.assertEqual(round(st.sM.kJkmolK, 1), 113.9)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 74.26)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 100.1)
        self.assertEqual(round(st.w, 0), 1698)
        self.assertEqual(round(st.mu.muPas, -1), 6160)
        self.assertEqual(round(st.k, 3), 0.116)

        st = iC4(T=200, P=5e4, **kw)
        self.assertEqual(round(st.rho, 1), 657.6)
        self.assertEqual(round(st.rhoM, 2), 11.31)
        self.assertEqual(round(st.uM.kJkmol, -1), -14600)
        self.assertEqual(round(st.hM.kJkmol, -1), -14590)
        self.assertEqual(round(st.sM.kJkmolK, 1), 168.3)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 82.19)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 114.9)
        self.assertEqual(round(st.w, 0), 1319)
        self.assertEqual(round(st.mu.muPas, 0), 545)
        self.assertEqual(round(st.k, 3), 0.135)

        st = iC4(T=520, P=1e5, **kw)
        self.assertEqual(round(st.rho, 3), 1.350)
        self.assertEqual(round(st.rhoM, 5), 0.02323)
        self.assertEqual(round(st.uM.kJkmol, -1), 41570)
        self.assertEqual(round(st.hM.kJkmol, -1), 45880)
        self.assertEqual(round(st.sM.kJkmolK, 1), 364.1)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 145.4)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 153.9)
        self.assertEqual(round(st.w, 1), 279.4)
        self.assertEqual(round(st.mu.muPas, 1), 12.8)
        self.assertEqual(round(st.k, 4), 0.0452)

        st = iC4(T=300, P=101325, **kw)
        self.assertEqual(round(st.rho, 3), 2.423)
        self.assertEqual(round(st.rhoM, 5), 0.04169)
        self.assertEqual(round(st.uM.kJkmol, -1), 15470)
        self.assertEqual(round(st.hM.kJkmol, -1), 17900)
        self.assertEqual(round(st.sM.kJkmolK, 1), 295.5)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 89.80)
        self.assertEqual(round(st.cpM.kJkmolK, 2), 99.19)
        self.assertEqual(round(st.w, 1), 212.1)
        self.assertEqual(round(st.mu.muPas, 2), 7.56)
        self.assertEqual(round(st.k, 4), 0.0161)

        st = iC4(T=250, P=2e5, **kw)
        self.assertEqual(round(st.rho, 1), 606.4)
        self.assertEqual(round(st.rhoM, 2), 10.43)
        self.assertEqual(round(st.uM.kJkmol, 0), -8563)
        self.assertEqual(round(st.hM.kJkmol, 0), -8543)
        self.assertEqual(round(st.sM.kJkmolK, 1), 195.2)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 90.00)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 127.2)
        self.assertEqual(round(st.w, 0), 1070)
        self.assertEqual(round(st.mu.muPas, 0), 263)
        self.assertEqual(round(st.k, 3), 0.117)

        st = iC4(T=600, P=3e5, **kw)
        self.assertEqual(round(st.rho, 3), 3.523)
        self.assertEqual(round(st.rhoM, 5), 0.06062)
        self.assertEqual(round(st.uM.kJkmol, -1), 53820)
        self.assertEqual(round(st.hM.kJkmol, -1), 58770)
        self.assertEqual(round(st.sM.kJkmolK, 1), 378.1)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 162.1)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 170.8)
        self.assertEqual(round(st.w, 1), 298.3)
        self.assertEqual(round(st.mu.muPas, 1), 14.6)
        self.assertEqual(round(st.k, 4), 0.0578)

        st = iC4(T=330, P=5e5, **kw)
        self.assertEqual(round(st.rho, 2), 11.77)
        self.assertEqual(round(st.rhoM, 4), 0.2024)
        self.assertEqual(round(st.uM.kJkmol, -1), 17780)
        self.assertEqual(round(st.hM.kJkmol, -1), 20250)
        self.assertEqual(round(st.sM.kJkmolK, 1), 290.4)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 100.4)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 113.8)
        self.assertEqual(round(st.w, 1), 206.9)
        self.assertEqual(round(st.mu.muPas, 2), 8.61)
        self.assertEqual(round(st.k, 4), 0.0200)

        st = iC4(T=320, P=8e5, **kw)
        self.assertEqual(round(st.rho, 1), 522.3)
        self.assertEqual(round(st.rhoM, 3), 8.986)
        self.assertEqual(round(st.uM.kJkmol, 0), 1074)
        self.assertEqual(round(st.hM.kJkmol, 0), 1163)
        self.assertEqual(round(st.sM.kJkmolK, 1), 229.1)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 105.0)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 151.7)
        self.assertEqual(round(st.w, 1), 692.0)
        self.assertEqual(round(st.mu.muPas, 0), 120)
        self.assertEqual(round(st.k, 4), 0.0866)

        st = iC4(T=460, P=1e6, **kw)
        self.assertEqual(round(st.rho, 2), 16.26)
        self.assertEqual(round(st.rhoM, 4), 0.2798)
        self.assertEqual(round(st.uM.kJkmol, -1), 32750)
        self.assertEqual(round(st.hM.kJkmol, -1), 36330)
        self.assertEqual(round(st.sM.kJkmolK, 1), 325.9)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 132.6)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 143.7)
        self.assertEqual(round(st.w, 1), 249.2)
        self.assertEqual(round(st.mu.muPas, 1), 11.9)
        self.assertEqual(round(st.k, 4), 0.0371)

        st = iC4(T=372, P=2e6, **kw)
        self.assertEqual(round(st.rho, 1), 432.9)
        self.assertEqual(round(st.rhoM, 3), 7.448)
        self.assertEqual(round(st.uM.kJkmol, 0), 9576)
        self.assertEqual(round(st.hM.kJkmol, 0), 9844)
        self.assertEqual(round(st.sM.kJkmolK, 1), 253.7)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 119.3)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 193.0)
        self.assertEqual(round(st.w, 1), 388.4)
        self.assertEqual(round(st.mu.muPas, 1), 66.1)
        self.assertEqual(round(st.k, 4), 0.0674)

        st = iC4(T=400, P=3e6, **kw)
        self.assertEqual(round(st.rho, 2), 96.31)
        self.assertEqual(round(st.rhoM, 3), 1.657)
        self.assertEqual(round(st.uM.kJkmol, -1), 22000)
        self.assertEqual(round(st.hM.kJkmol, -1), 23810)
        self.assertEqual(round(st.sM.kJkmolK, 1), 289.0)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 130.2)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 250.2)
        self.assertEqual(round(st.w, 1), 145.1)
        self.assertEqual(round(st.mu.muPas, 1), 13.9)
        # The point near critical point differ in thermal conductivity
        # because the reference calculate compresibility (Kt) with a special
        # scaled equation
        # self.assertEqual(round(st.k, 4), 0.0404)

        st = iC4(T=438, P=4e6, **kw)
        self.assertEqual(round(st.rho, 1), 103.9)
        self.assertEqual(round(st.rhoM, 3), 1.788)
        self.assertEqual(round(st.uM.kJkmol, -1), 26780)
        self.assertEqual(round(st.hM.kJkmol, -1), 29010)
        self.assertEqual(round(st.sM.kJkmolK, 1), 300.0)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 135.9)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 199.2)
        self.assertEqual(round(st.w, 1), 173.6)
        self.assertEqual(round(st.mu.muPas, 1), 15.3)
        self.assertEqual(round(st.k, 4), 0.0413)

        st = iC4(T=200, P=5e6, **kw)
        self.assertEqual(round(st.rho, 1), 661.5)
        self.assertEqual(round(st.rhoM, 2), 11.38)
        self.assertEqual(round(st.uM.kJkmol, -1), -14720)
        self.assertEqual(round(st.hM.kJkmol, -1), -14290)
        self.assertEqual(round(st.sM.kJkmolK, 1), 167.6)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 82.21)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 114.4)
        self.assertEqual(round(st.w, 0), 1352)
        self.assertEqual(round(st.mu.muPas, 0), 577)
        self.assertEqual(round(st.k, 3), 0.137)

        st = iC4(T=440, P=6e6, **kw)
        self.assertEqual(round(st.rho, 1), 255.1)
        self.assertEqual(round(st.rhoM, 3), 4.389)
        self.assertEqual(round(st.uM.kJkmol, -1), 22760)
        self.assertEqual(round(st.hM.kJkmol, -1), 24130)
        self.assertEqual(round(st.sM.kJkmolK, 1), 287.2)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 139.4)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 298.3)
        self.assertEqual(round(st.w, 1), 172.8)
        self.assertEqual(round(st.mu.muPas, 1), 27.8)
        self.assertEqual(round(st.k, 4), 0.0589)

        st = iC4(T=200, P=7e6, **kw)
        self.assertEqual(round(st.rho, 1), 663.0)
        self.assertEqual(round(st.rhoM, 2), 11.41)
        self.assertEqual(round(st.uM.kJkmol, -1), -14770)
        self.assertEqual(round(st.hM.kJkmol, -1), -14160)
        self.assertEqual(round(st.sM.kJkmolK, 1), 167.4)
        self.assertEqual(round(st.cvM.kJkmolK, 2), 82.24)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 114.2)
        self.assertEqual(round(st.w, 0), 1364)
        self.assertEqual(round(st.mu.muPas, 0), 589)
        self.assertEqual(round(st.k, 3), 0.137)

        st = iC4(T=530, P=8e6, **kw)
        self.assertEqual(round(st.rho, 1), 148.0)
        self.assertEqual(round(st.rhoM, 3), 2.546)
        self.assertEqual(round(st.uM.kJkmol, -1), 38730)
        self.assertEqual(round(st.hM.kJkmol, -1), 41870)
        self.assertEqual(round(st.sM.kJkmolK, 1), 322.8)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 152.7)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 192.8)
        self.assertEqual(round(st.w, 1), 231.4)
        self.assertEqual(round(st.mu.muPas, 1), 20.2)
        self.assertEqual(round(st.k, 4), 0.0594)

        st = iC4(T=420, P=1e7, **kw)
        self.assertEqual(round(st.rho, 1), 405.0)
        self.assertEqual(round(st.rhoM, 3), 6.968)
        self.assertEqual(round(st.uM.kJkmol, -1), 16310)
        self.assertEqual(round(st.hM.kJkmol, -1), 17750)
        self.assertEqual(round(st.sM.kJkmolK, 1), 270.9)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 129.2)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 186.4)
        self.assertEqual(round(st.w, 1), 411.3)
        self.assertEqual(round(st.mu.muPas, 1), 57.5)
        self.assertEqual(round(st.k, 4), 0.0667)

        st = iC4(T=560, P=2e7, **kw)
        self.assertEqual(round(st.rho, 1), 303.1)
        self.assertEqual(round(st.rhoM, 3), 5.215)
        self.assertEqual(round(st.uM.kJkmol, -1), 39260)
        self.assertEqual(round(st.hM.kJkmol, -1), 43090)
        self.assertEqual(round(st.sM.kJkmolK, 1), 319.4)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 157.9)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 195.7)
        self.assertEqual(round(st.w, 1), 380.8)
        self.assertEqual(round(st.mu.muPas, 1), 37.4)
        self.assertEqual(round(st.k, 4), 0.0772)

        st = iC4(T=600, P=3.5e7, **kw)
        self.assertEqual(round(st.rho, 1), 358.6)
        self.assertEqual(round(st.rhoM, 3), 6.170)
        self.assertEqual(round(st.uM.kJkmol, -1), 44210)
        self.assertEqual(round(st.hM.kJkmol, -1), 49880)
        self.assertEqual(round(st.sM.kJkmolK, 1), 326.6)
        self.assertEqual(round(st.cvM.kJkmolK, 1), 165.6)
        self.assertEqual(round(st.cpM.kJkmolK, 1), 194.4)
        self.assertEqual(round(st.w, 1), 532.1)
        self.assertEqual(round(st.mu.muPas, 1), 48.0)
        self.assertEqual(round(st.k, 4), 0.0857)

    def test_shortSpan(self):
        # Table III, Pag 46
        st = iC4(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.2393)
        self.assertEqual(round(st.P.MPa, 3), 19.108)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.5575)

        st2 = iC4(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 210.32)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.37469)

    def test_herrmann(self):
        # Table 6, Pag 17

        # This point isn't real, it's in two phases region so need force
        # calculation
        st = iC4(T=115, rho=739)
        mu = st._Viscosity(739, 115, None)
        self.assertEqual(round(mu.muPas, 3), 8255.646)

        self.assertEqual(round(iC4(T=300, rho=0).mu.muPas, 6), 7.539691)
        self.assertEqual(round(iC4(T=300, rho=1).mu.muPas, 6), 7.528194)
        self.assertEqual(round(iC4(T=300, rho=9).mu.muPas, 6), 7.454569)
        self.assertEqual(round(iC4(T=300, rho=550).mu.muPas, 4), 148.1347)
        self.assertEqual(round(iC4(T=300, rho=625).mu.muPas, 4), 310.3424)
        self.assertEqual(round(iC4(T=400, rho=1).mu.muPas, 6), 9.889178)
        self.assertEqual(round(iC4(T=400, rho=118).mu.muPas, 5), 13.52154)
        self.assertEqual(round(iC4(T=400, rho=345).mu.muPas, 5), 41.45875)
        self.assertEqual(round(iC4(T=400, rho=560).mu.muPas, 4), 159.6867)
        self.assertEqual(round(iC4(T=407.81, rho=225.5).mu.muPas, 5), 23.45602)
        self.assertEqual(round(iC4(T=500, rho=1).mu.muPas, 5), 12.15232)
        self.assertEqual(round(iC4(T=500, rho=100).mu.muPas, 5), 15.45786)
        self.assertEqual(round(iC4(T=500, rho=500).mu.muPas, 4), 102.1956)
        self.assertEqual(round(iC4(T=650, rho=0).mu.muPas, 5), 15.39019)
        self.assertEqual(round(iC4(T=650, rho=1).mu.muPas, 5), 15.40515)
        self.assertEqual(round(iC4(T=650, rho=100).mu.muPas, 5), 19.01606)
        self.assertEqual(round(iC4(T=650, rho=420).mu.muPas, 5), 63.72888)
