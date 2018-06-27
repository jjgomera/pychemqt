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
from copy import copy

from scipy import exp
from iapws import _Viscosity, _ThCond, _Dielectric, _Tension

from lib import unidades
from lib.meos import MEoS


class H2O(MEoS):
    """Multiparameter equation of state for water (including IAPWS95)"""

    name = "water"
    CASNumber = "7732-18-5"
    formula = "H2O"
    synonym = "R-718"
    _refPropName = "WATER"
    Tc = unidades.Temperature(647.096)
    rhoc = unidades.Density(322.)
    Pc = unidades.Pressure(22064.0, "kPa")
    M = 18.015268  # g/mol
    Tt = unidades.Temperature(273.16)
    Tb = unidades.Temperature(373.1243)
    f_acent = 0.3443
    momentoDipolar = unidades.DipoleMoment(1.855, "Debye")
    id = 62

    Fi1 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.3204464837497, 6.6832105275932],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.28728967, 3.53734222, 7.74073708, 9.24437796,
                     27.5075105]}

    Fi2 = {"ao_log": [1, 3.00392],
           "pow": [0, 1],
           "ao_pow": [8.203520690, -11.996306443],
           "ao_exp": [], "titao": [],
           "ao_hyp": [0.01059, -0.98763, 3.06904, 0],
           "hyp": [0.415386589, 1.763895929, 3.874803739, 0]}

    Fi3 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.318441, 6.681816],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.287202151, 3.537101709, 7.740210774, 9.243749421,
                     27.5056402]}

    Fi4 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.3177095, 6.6815049],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.287202151, 3.537101709, 7.740210774, 9.243749421,
                     27.5056402]}

    iapws = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for water of Wagner and"
                    u"Pruß (2002).",
        "__doi__": {"autor": u"Wagner, W., Pruß, A.",
                    "title": "The IAPWS Formulation 1995 for the Thermodynamic"
                             "Properties of Ordinary Water Substance for"
                             "General and Scientific Use",
                    "ref": "J. Phys. Chem. Ref. Data 31, 387 (2002)",
                    "doi": "10.1063/1.1461829"},

        "R": 8.314371357587,
        "cp": Fi1,
        "ref": {"Tref": Tt, "Pref": 0.611655, "ho": 0.611872, "so": 0},

        "Tmin": Tt, "Tmax": 2000., "Pmax": 2000000.0, "rhomax": 73.96,
        "Pmin": 0.61248, "rhomin": 55.49696,

        "nr1": [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
                0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
                0.88089493102134e-2],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1],

        "nr2": [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
                -0.19232721156002, -0.25709043003438, 0.16074868486251,
                -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
                0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
                .36582165144204e-6, -.13251180074668e-11, -0.62639586912454e-9,
                -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
                -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
                -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
                0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
                0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
                -0.20393486513704e-1, -.16554050063734e-2, 0.19955571979541e-2,
                0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
                0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
                -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
                0.31777497330738, -0.11841182425981],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,
               6, 6],
        "d2": [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
               4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14,
               3, 6, 6, 6],
        "t2": [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10,
               10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,
               23, 10, 50, 44, 46, 50],
        "gamma2": [1]*44,

        "nr3": [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4],
        "d3": [3]*3,
        "t3": [0, 1, 4],
        "alfa3": [20]*3,
        "beta3": [150, 150, 250],
        "gamma3": [1.21, 1.21, 1.25],
        "epsilon3": [1.]*3,

        "nr4": [-0.14874640856724, 0.31806110878444],
        "a4": [3.5, 3.5],
        "b4": [0.85, 0.95],
        "B": [0.2, 0.2],
        "C": [28, 32],
        "D": [700, 800],
        "A": [0.32, .32],
        "beta4": [0.3, 0.3]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Kunz and Wagner"
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

        "Tmin": Tt, "Tmax": 1350.0, "Pmax": 1000000.0, "rhomax": 73.96,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.82728408749586, -0.18602220416584e1, -0.11199009613744e1,
                0.15635753976056, 0.87375844859025, -0.36674403715731,
                0.53987893432436e-1],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [0.5, 1.25, 1.875, 0.125, 1.5, 1, 0.75],

        "nr2": [0.10957690214499e1, 0.53213037828563e-1, 0.13050533930825e-1,
                -0.41079520434476, 0.14637443344120, -0.55726838623719e-1,
                -0.0112017741438, -0.66062758068099e-2, .46918522004538e-2],
        "c2": [1, 1, 1, 2, 2, 2, 3, 5, 5],
        "d2": [1, 5, 5, 1, 2, 4, 4, 1, 1],
        "t2": [1.5, 0.625, 2.625, 5, 4, 4.5, 3, 4, 6],
        "gamma2": [1]*9}

    saul = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Saul and"
                    "Wagner-58 coeff (1989).",
        "__doi__": {"autor": "Saul, A. and Wagner, W.",
                    "title": "A Fundamental Equation for Water Covering the"
                             "Range from the Melting Line to 1273 K at"
                             "Pressures up to 25000 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 18, 1537 (1989)",
                    "doi": "10.1063/1.555836"},

        "R": 8.31434,
        "cp": Fi3,
        "ref": {"Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0},

        "Tmin": Tt, "Tmax": 1273., "Pmax": 400000.0, "rhomax": 55.49,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.8216377478, -0.2543894379, -0.08830868648, -0.8903097248e-6,
                -0.1241333357e-5, 0.2895590286e-8, 0.1403610309e-10,
                0.8183943371e-12, -0.2397905287e-12],
        "d1": [1, 1, 2, 5, 8, 11, 11, 13, 13],
        "t1": [0, 2, 0, 9, 0, 0, 12, 7, 13],

        "nr2": [-0.7519743341, -0.4151278588, -0.103051374e1, -0.1648036888e1,
                -0.4686350251, 0.3560258142, -0.6364658294, 0.2227482363,
                -0.8954849939e-1, 0.1557686788e-2, 0.1347719088e-2,
                -0.1301353385e-2, 0.9987368673e-6, 0.2263629476e-3,
                0.289330495e-5, 0.1995437169, -0.2707767662e-1,
                0.1849068216e-1, -0.4402394357e-2, -0.8546876737e-1,
                0.1220538576, -0.2562237041, 0.2555034636, -0.6323203907e-1,
                0.3351397575e-4, -0.6152834985e-1, -0.3533048208e-3,
                0.3146309259e-1, -0.2261795983e-2, 0.18689702e-3,
                -0.1384614556e-2, 0.2713160073e-2, -0.4866118539e-2,
                0.3751789129e-2, -0.5692669373e-3, -0.5876414555, 0.5687838346,
                -0.1642158198, 0.5878635885, -0.2844301931, -0.2049198337,
                -0.4039233716e-2, 0.5459049594e-1, -0.8914260146e-2,
                0.4974411254e-2],
        "c2": [1]*15+[2]*20+[3]*10,
        "d2": [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 11, 1, 1, 1, 1, 2, 2,
               4, 5, 6, 6, 7, 7, 8, 10, 10, 11, 11, 11, 11, 11, 2, 2, 2, 3, 3,
               4, 4, 5, 5, 5],
        "t2": [0, 1, 3, 1, 5, 5, 2, 3, 5, 6, 4, 1, 8, 0, 1, 0, 9, 10, 11, 0,
               8, 5, 4, 2, 12, 3, 10, 3, 2, 8, 0, 1, 3, 4, 6, 13, 14, 15, 14,
               16, 13, 26, 15, 23, 25],
        "gamma2": [1]*45,

        "nr5": [-0.709318338e-2, 0.1718796342e-1, -0.1482653038e-1,
                0.4517292884e-2],
        "d5": [1, 2, 3, 4],
        "t5": [50, 40, 32, 26]
        }

    saul2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Saul and"
                    "Wagner-38 coeff (1989).",
        "__doi__": {"autor": "Saul, A. and Wagner, W.",
                    "title": "A Fundamental Equation for Water Covering the"
                             "Range from the Melting Line to 1273 K at"
                             "Pressures up to 25000 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 18, 1537 (1989)",
                    "doi":  "10.1063/1.555836"},

        "R": 8.31434,
        "cp": Fi4,
        "ref": {"Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0},

        "Tmin": Tt, "Tmax": 1273., "Pmax": 400000.0, "rhomax": 55.49,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1": [0.2330009013, -0.1402091128e1, 0.1172248041, -0.1850749499,
                0.1770110422, .5525151794e-1, -0.341325738e-3, 0.8557274367e-3,
                0.3716900685e-3, -0.1308871233e-3, 0.3216895199e-4,
                0.2785881034e-6],
        "d1": [1, 1, 2, 2, 2, 2, 3, 5, 5, 6, 7, 8],
        "t1": [0, 2, 0, 1, 2, 3, 5, 0, 1, 3, 2, 5],

        "nr2": [-0.352151113, 0.7881914536e-1, -0.151966661e-1, -0.1068458586,
                -0.2055046288, 0.9146198012, 0.3213343569e-3, -0.1133591391e1,
                -0.3107520749, 0.1217901527e1, -0.4481710831, 0.5494218772e-1,
                -0.8665222096e-4, 0.3844084088e-1, 0.9853044884e-2,
                -0.1767598472e-1, 0.1488549222e-2, -0.3070719069e-2,
                0.388080328e-2, -.2627505215e-2, .5258371388e-3, -0.1716396901,
                .7188823624e-1, .5881268357e-1, -.145593888e-1, -.12161394e-1],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 3, 3, 3, 3, 3],
        "d2": [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 11, 11,
               11, 2, 2, 3, 3, 5],
        "t2": [5, 7, 9, 5, 4, 6, 13, 5, 2, 3, 2, 0, 11, 1, 4, 0, 0, 3, 5, 6,
               7, 13, 14, 15, 24, 15],
        "gamma2": [1]*26}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Sun and Ely"
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering"
                             "application: Algorithm and  application to"
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314371357587,
        "cp": Fi1,
        "ref": {"name": "CUSTOM",
                "Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [3.46821920e-1, 5.03423025e-1, -3.51059570e-1, 5.07004866e-2,
                1.99939129e-4, -5.69888763e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.96198912e-1, -2.02509554, -1.09353609, 7.25785202e-2,
                2.16072642e-1, -1.01542630e-1, 7.46926106e-2, 2.18830463e-3],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = iapws, GERG, saul, saul2, sun
    _PR = 0.0043451

    _surface = {"sigma": [-0.1306, 0.2151], "exp": [2.471, 1.233]}
    _melting = {"Tmin": 251.165, "Tmax": 370.0}
    _sublimation = {"Tmin": 50.0, "Tmax": Tt}

    _vapor_Pressure = {
        "eq": 6,
        "ao": [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719,
               1.80122502],
        "exp": [2, 3, 6, 7, 8, 15]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352,
               -6.74694450e5],
        "exp": [1, 2, 5, 16, 43, 110]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581,
               -63.9201063],
        "exp": [1, 2, 4, 9, 18.5, 35.5]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (1997)",
              "__doi__": {
                  "autor": "Huber, M.L., Perkins, R.A., Lasecke, A., Friend, "
                           "D.G., Sengers, J.V., Assael, M.J., Metaxa, I.N.,"
                           "Vogel, E., Mareš, R., Miyagawa, K.",
                  "title": "New International Formulation for the Viscosity"
                           "of H2O",
                  "ref": "J. Phys. Chem. Ref. Data 38(2) (2009) 101-125",
                  "doi": "10.1063/1.3088050"},
              "__code__": (_Viscosity, )}

    def _visco0(self, rho, T, fase):
        ref = H2O()
        ref._ref("OTO")
        estado = ref._Helmholtz(rho, 1.5*self.Tc)
        drho = 1/estado["dpdrho"]*1e6

        # convert ∂ρ/∂P]τ to IAPWS units, [kg/m³·MPa]
        if fase:
            fase = copy(fase)
            fase.drhodP_T *= 1e6

        mu = _Viscosity(rho, T, fase=fase, drho=drho)
        return unidades.Viscosity(mu)

    visco1 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 647.096,
              "no": [151.138, -444.318, 398.262, -81.7008],
              "to": [0, 0.25, 0.5, 0.75],

              "a": [-1.17407e-5, -3.78855e-7, 3.56743e-8],
              "b": [1.62216e-6, -8.36595e-6, 9.10863e-8],
              "c": [1.92707e-5, -1.2868e-5, 0.0],
              "A": [-3.30145e-10, 0.0, 1.02931e-11],
              "B": [5.03139997945133e-10, 1.82304e-10, 0.0],
              "C": [8.01449e-10, 5.65614e-9, 1.10163e-10]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1997)",
               "__doi__": {
                   "autor": "Huber, M.L., Perkins, R.A., Friend, D.G., "
                            "Sengers, J.V., Assael, M.J., Metaxa, I.N., "
                            "Miyagawa, K., Hellmann, R., Vogel, E.",
                   "title": "New International Formulation for the Thermal"
                            "Conductivity of H20",
                   "ref": "J. Phys. Chem. Ref. Data 41(3) (2012) 033102",
                   "doi": "10.1063/1.4738955"},
               "__code__": (_ThCond, )}

    _thermal = thermo0,

    def _thermo0(self, rho, T, fase):
        ref = H2O()
        ref._ref("OTO")
        estado = ref._Helmholtz(rho, 1.5*self.Tc)
        drho = 1/estado["dpdrho"]*1e6

        # convert values to IAPWS units
        # ∂ρ/∂P]τ, [kg/m³·MPa]
        # cp, [kJ/kg]
        if fase:
            fase = copy(fase)
            fase.drhodP_T *= 1e6
            fase.cp /= 1e3

        return _ThCond(rho, T, fase, drho)

    def _Dielectric(self, rho, T):
        try:
            nu = _Dielectric(rho, T)
        except NotImplementedError:
            nu = None
        return unidades.Dimensionless(nu)

    def _Surface(self, T):
        """Equation for the surface tension"""
        try:
            s = _Tension(T)
        except NotImplementedError:
            s = None
        return unidades.Tension(s)

    @classmethod
    def _Melting_Pressure(cls, T):
        if 251.165 <= T <= 256.164:
            Tref = 251.165
            Pref = 208566.
            Tita = T/Tref
            P2 = Pref*(1-0.299948*(1-Tita**60.))
        elif 256.164 < T <= 273.31:
            Tref = 256.164
            Pref = 350100.
            Tita = T/Tref
            P2 = Pref*(1-1.18721*(1-Tita**8.))
        elif 273.31 < T <= 355:
            Tref = 273.31
            Pref = 632400.
            Tita = T/Tref
            P2 = Pref*(1-1.07476*(1-Tita**4.6))
        elif 355. < T:
            Tref = 355
            Pref = 2216000.
            Tita = T/Tref
            P2 = Pref*exp(1.73683*(1-1./Tita)-0.544606e-1*(1-Tita**5) +
                          0.806106e-7*(1-Tita**22))
        return unidades.Pressure(P2, "kPa")

    @classmethod
    def _Sublimation_Pressure(cls, T):
        Pref = 611.657
        Tita = T/cls.Tt
        a = [-0.212144006e2, 0.273203819e2, -0.61059813e1]
        expo = [0.333333333e-2, 1.20666667, 1.70333333]
        suma = 0
        for ai, expi in zip(a, expo):
            suma += ai*Tita**expi
        return unidades.Pressure(exp(suma/Tita)*Pref)


class Test(TestCase):

    def test_iapws(self):
        # Table 6.6, pag 436"""
        fluid = H2O()

        delta = 838.025/fluid.rhoc
        tau = fluid.Tc/500
        ideal = fluid._phi0(fluid._constants["cp"], tau, delta)
        self.assertEqual(round(ideal["fio"], 8), 2.04797733)
        self.assertEqual(round(ideal["fiod"], 9), 0.384236747)
        self.assertEqual(round(ideal["fiodd"], 9), -0.147637878)
        self.assertEqual(round(ideal["fiot"], 8), 9.04611106)
        self.assertEqual(round(ideal["fiott"], 8), -1.93249185)
        self.assertEqual(round(ideal["fiodt"], 8), 0.0)
        res = fluid._phir(tau, delta)
        self.assertEqual(round(res["fir"], 8), -3.42693206)
        self.assertEqual(round(res["fird"], 9), -0.364366650)
        self.assertEqual(round(res["firdd"], 9), 0.856063701)
        self.assertEqual(round(res["firt"], 8), -5.81403435)
        self.assertEqual(round(res["firtt"], 8), -2.23440737)
        self.assertEqual(round(res["firdt"], 8), -1.12176915)

        delta = 358/fluid.rhoc
        tau = fluid.Tc/647
        ideal = fluid._phi0(fluid._constants["cp"], tau, delta)
        self.assertEqual(round(ideal["fio"], 8), -1.56319605)
        self.assertEqual(round(ideal["fiod"], 9), 0.899441341)
        self.assertEqual(round(ideal["fiodd"], 9), -0.808994726)
        self.assertEqual(round(ideal["fiot"], 8), 9.80343918)
        self.assertEqual(round(ideal["fiott"], 8), -3.43316334)
        self.assertEqual(round(ideal["fiodt"], 8), 0.0)
        res = fluid._phir(tau, delta)
        self.assertEqual(round(res["fir"], 8), -1.21202657)
        self.assertEqual(round(res["fird"], 9), -0.714012024)
        self.assertEqual(round(res["firdd"], 9), 0.475730696)
        self.assertEqual(round(res["firt"], 8), -3.21722501)
        self.assertEqual(round(res["firtt"], 8), -9.96029507)
        self.assertEqual(round(res["firdt"], 8), -1.33214720)

        # Table 13.1, Pag 486, Saturation state
        st = H2O(T=H2O.Tt, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000612)
        self.assertEqual(round(st.Liquido.rho, 3), 999.793)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 0.001)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.0000)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 4.2174)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 4.2199)
        self.assertEqual(round(st.Liquido.w, 1), 1402.3)
        self.assertEqual(round(st.Gas.rho, 5), 0.00485)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2500.92)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 9.1555)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.4184)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.8844)
        self.assertEqual(round(st.Gas.w, 2), 409.00)

        st = H2O(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.003537)
        self.assertEqual(round(st.Liquido.rho, 3), 996.513)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 112.565)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.3931)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 4.1305)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 4.1809)
        self.assertEqual(round(st.Liquido.w, 1), 1501.4)
        self.assertEqual(round(st.Gas.rho, 5), 0.02559)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2549.85)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 8.5174)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.4422)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.9141)
        self.assertEqual(round(st.Gas.w, 2), 427.89)

        st = H2O(T=400, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.24577)
        self.assertEqual(round(st.Liquido.rho, 3), 937.486)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 532.953)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.6013)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 3.6324)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 4.2555)
        self.assertEqual(round(st.Liquido.w, 1), 1509.5)
        self.assertEqual(round(st.Gas.rho, 4), 1.3694)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2715.70)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 7.0581)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.6435)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 2.2183)
        self.assertEqual(round(st.Gas.w, 2), 484.67)

        st = H2O(T=500, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 2.6392)
        self.assertEqual(round(st.Liquido.rho, 3), 831.313)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 975.431)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 2.5810)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 3.2255)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 4.6635)
        self.assertEqual(round(st.Liquido.w, 1), 1239.6)
        self.assertEqual(round(st.Gas.rho, 3), 13.199)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2802.48)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 6.2351)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 2.2714)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 3.4631)
        self.assertEqual(round(st.Gas.w, 2), 504.55)

        st = H2O(T=600, x=0.5)
        self.assertEqual(round(st.P.MPa, 3), 12.345)
        self.assertEqual(round(st.Liquido.rho, 3), 649.411)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 1505.36)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 3.5190)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 3.0475)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 6.9532)
        self.assertEqual(round(st.Liquido.w, 2), 749.57)
        self.assertEqual(round(st.Gas.rho, 3), 72.842)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2677.81)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 5.4731)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 3.3271)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 9.1809)
        self.assertEqual(round(st.Gas.w, 2), 457.33)

        st = H2O(T=647, x=0.5)
        self.assertEqual(round(st.P.MPa, 3), 22.038)
        self.assertEqual(round(st.Liquido.rho, 2), 357.34)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 2029.44)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 4.3224)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 6.2344)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 1), 3905.2)
        self.assertEqual(round(st.Liquido.w, 2), 251.19)
        self.assertEqual(round(st.Gas.rho, 2), 286.51)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 2148.56)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 4.5065)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.2740)
        self.assertEqual(round(st.Gas.cp.kJkgK, 1), 5334.1)
        self.assertEqual(round(st.Gas.w, 2), 285.32)

        # Table 13.2, Pag 495, Single phase region
        st = H2O(T=275, P=5e4)
        self.assertEqual(round(st.rho, 3), 999.912)
        self.assertEqual(round(st.u.kJkg, 3), 7.760)
        self.assertEqual(round(st.h.kJkg, 3), 7.810)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0283)
        self.assertEqual(round(st.cv.kJkgK, 4), 4.2130)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.2137)
        self.assertEqual(round(st.w, 1), 1411.4)

        st = H2O(T=370, P=1e5)
        self.assertEqual(round(st.rho, 3), 960.591)
        self.assertEqual(round(st.u.kJkg, 3), 405.787)
        self.assertEqual(round(st.h.kJkg, 3), 405.891)
        self.assertEqual(round(st.s.kJkgK, 4), 1.2715)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.7845)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.2121)
        self.assertEqual(round(st.w, 1), 1545.8)

        st = H2O(T=430, P=5e5)
        self.assertEqual(round(st.rho, 4), 2.6297)
        self.assertEqual(round(st.u.kJkg, 2), 2569.90)
        self.assertEqual(round(st.h.kJkg, 2), 2760.04)
        self.assertEqual(round(st.s.kJkgK, 4), 6.8486)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.7170)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.3469)
        self.assertEqual(round(st.w, 2), 497.85)

        st = H2O(T=275, P=1e6)
        self.assertEqual(round(st.rho, 2), 1000.39)
        self.assertEqual(round(st.u.kJkg, 3), 7.769)
        self.assertEqual(round(st.h.kJkg, 3), 8.768)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0283)
        self.assertEqual(round(st.cv.kJkgK, 4), 4.2087)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.2093)
        self.assertEqual(round(st.w, 1), 1413.0)

        st = H2O(T=1273, P=2e6)
        self.assertEqual(round(st.rho, 4), 3.4085)
        self.assertEqual(round(st.u.kJkg, 2), 4049.90)
        self.assertEqual(round(st.h.kJkg, 2), 4636.67)
        self.assertEqual(round(st.s.kJkgK, 4), 8.5933)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.0205)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.4901)
        self.assertEqual(round(st.w, 2), 849.82)

        st = H2O(T=530, P=5e6)
        self.assertEqual(round(st.rho, 3), 789.208)
        self.assertEqual(round(st.u.kJkg, 2), 1112.95)
        self.assertEqual(round(st.h.kJkg, 2), 1119.28)
        self.assertEqual(round(st.s.kJkgK, 4), 2.8547)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.1377)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.9386)
        self.assertEqual(round(st.w, 1), 1122.2)

        st = H2O(T=375, P=1e7)
        self.assertEqual(round(st.rho, 3), 961.618)
        self.assertEqual(round(st.u.kJkg, 3), 423.977)
        self.assertEqual(round(st.h.kJkg, 3), 434.376)
        self.assertEqual(round(st.s.kJkgK, 4), 1.3203)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.7438)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.1955)
        self.assertEqual(round(st.w, 1), 1562.6)

        st = H2O(T=630, P=2e7)
        self.assertEqual(round(st.rho, 3), 567.644)
        self.assertEqual(round(st.u.kJkg, 2), 1671.50)
        self.assertEqual(round(st.h.kJkg, 2), 1706.73)
        self.assertEqual(round(st.s.kJkgK, 4), 3.8259)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.1104)
        self.assertEqual(round(st.cp.kJkgK, 4), 9.8619)
        self.assertEqual(round(st.w, 2), 587.65)

        st = H2O(T=900, P=3e7)
        self.assertEqual(round(st.rho, 3), 82.840)
        self.assertEqual(round(st.u.kJkg, 2), 3167.87)
        self.assertEqual(round(st.h.kJkg, 2), 3530.02)
        self.assertEqual(round(st.s.kJkgK, 4), 6.3313)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.0283)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.0408)
        self.assertEqual(round(st.w, 2), 687.91)

        st = H2O(T=1100, P=5e7)
        self.assertEqual(round(st.rho, 2), 106.26)
        self.assertEqual(round(st.u.kJkg, 2), 3534.70)
        self.assertEqual(round(st.h.kJkg, 2), 4005.26)
        self.assertEqual(round(st.s.kJkgK, 4), 6.5956)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.0843)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.9377)
        self.assertEqual(round(st.w, 2), 789.83)

        st = H2O(T=600, P=1e8)
        self.assertEqual(round(st.rho, 3), 791.493)
        self.assertEqual(round(st.u.kJkg, 2), 1322.23)
        self.assertEqual(round(st.h.kJkg, 2), 1448.58)
        self.assertEqual(round(st.s.kJkgK, 4), 3.2256)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.9295)
        self.assertEqual(round(st.cp.kJkgK, 4), 4.5019)
        self.assertEqual(round(st.w, 1), 1300.4)

        st = H2O(T=400, P=2e8)
        self.assertEqual(round(st.rho, 2), 1017.12)
        self.assertEqual(round(st.u.kJkg, 3), 480.471)
        self.assertEqual(round(st.h.kJkg, 3), 677.105)
        self.assertEqual(round(st.s.kJkgK, 4), 1.4520)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.4513)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.9363)
        self.assertEqual(round(st.w, 1), 1882.7)

        st = H2O(T=305, P=1e9)
        self.assertEqual(round(st.rho, 2), 1234.94)
        self.assertEqual(round(st.u.kJkg, 3), 93.201)
        self.assertEqual(round(st.h.kJkg, 3), 902.957)
        self.assertEqual(round(st.s.kJkgK, 4), 0.0895)
        self.assertEqual(round(st.cv.kJkgK, 4), 3.4296)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.7845)
        self.assertEqual(round(st.w, 1), 2723.2)

        st = H2O(T=1273, P=1e9)
        self.assertEqual(round(st.rho, 2), 809.28)
        self.assertEqual(round(st.u.kJkg, 2), 3097.36)
        self.assertEqual(round(st.h.kJkg, 2), 4333.03)
        self.assertEqual(round(st.s.kJkgK, 4), 5.2048)
        self.assertEqual(round(st.cv.kJkgK, 4), 2.6446)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.4245)
        self.assertEqual(round(st.w, 1), 2095.5)

    def test_Viscosity(self):
        # Table 6, pag 116
        self.assertEqual(round(_Viscosity(998, 298.15)*1e6, 6), 889.735100)
        self.assertEqual(round(_Viscosity(1200, 298.15)*1e6, 6), 1437.649467)
        self.assertEqual(round(_Viscosity(1000, 373.15)*1e6, 6), 307.883622)
        self.assertEqual(round(_Viscosity(1, 433.15)*1e6, 6), 14.538324)
        self.assertEqual(round(_Viscosity(1000, 433.15)*1e6, 6), 217.685358)
        self.assertEqual(round(_Viscosity(1, 873.15)*1e6, 6), 32.619287)
        self.assertEqual(round(_Viscosity(100, 873.15)*1e6, 6), 35.802262)
        self.assertEqual(round(_Viscosity(600, 873.15)*1e6, 6), 77.430195)
        self.assertEqual(round(_Viscosity(1, 1173.15)*1e6, 6), 44.217245)
        self.assertEqual(round(_Viscosity(100, 1173.15)*1e6, 6), 47.640433)
        self.assertEqual(round(_Viscosity(400, 1173.15)*1e6, 6), 64.154608)

        # Table 7, pag 116
        fluid = H2O(rho=122, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 25.520677)
        fluid = H2O(rho=222, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 31.337589)
        fluid = H2O(rho=272, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 36.228143)
        fluid = H2O(rho=322, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 42.961579)
        fluid = H2O(rho=372, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 45.688204)
        fluid = H2O(rho=422, T=647.35)
        self.assertEqual(round(fluid.mu.muPas, 6), 49.436256)

    def test_ThCond(self):
        # Table 7, pag 12
        self.assertEqual(round(_ThCond(0, 298.15)*1000, 7), 18.4341883)
        self.assertEqual(round(_ThCond(998, 298.15)*1000, 6), 607.712868)
        self.assertEqual(round(_ThCond(1200, 298.15)*1000, 6), 799.038144)
        self.assertEqual(round(_ThCond(0, 873.15)*1000, 7), 79.1034659)

        # Table 8, pag 13
        fluid = H2O(rho=1, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 7), 51.9298924)
        fluid = H2O(rho=122, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 130.922885)
        fluid = H2O(rho=222, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 367.787459)
        fluid = H2O(rho=272, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 757.959776)
        fluid = H2O(rho=322, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 5), 1443.75556)
        fluid = H2O(rho=372, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 650.319402)
        fluid = H2O(rho=422, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 448.883487)
        fluid = H2O(rho=750, T=647.35)
        self.assertEqual(round(fluid.k.mWmK, 6), 600.961346)
