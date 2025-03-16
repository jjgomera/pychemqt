#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


from copy import copy
from unittest import TestCase

from iapws._iapws import _D2O_Viscosity, _D2O_ThCond, _D2O_Tension
from iapws._iapws import _D2O_Melting_Pressure, _D2O_Sublimation_Pressure

from lib import unidades
from lib.meos import MEoS


class D2O(MEoS):
    """Multiparameter equation of state for heavy water"""
    name = "heavy water"
    CASNumber = "7789-20-0"
    formula = "D2O"
    synonym = "deuterium oxide"
    _refPropName = "D2O"
    _coolPropName = "HeavyWater"
    Tc = unidades.Temperature(643.847)
    rhoc = unidades.Density(356)
    Pc = unidades.Pressure(21661.8, "kPa")
    M = 20.027508  # g/mol
    Tt = unidades.Temperature(276.969)
    Tb = unidades.Temperature(374.549)
    f_acent = 0.364
    momentoDipolar = unidades.DipoleMoment(1.9, "Debye")

    Fi0 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-8.670994022646, 6.96033578458778],
           "ao_exp": [0.010633, 0.99787, 2.1483, 0.3549],
           "titao": [308/Tc, 1695/Tc, 3949/Tc, 10317/Tc]}

    Fi1 = {"ao_log": [0.5399322597e-2, 0],
           "pow": [0, 1, 2, 3, 4, 5],
           "ao_pow": [0.3087155964e2, -.3827264031e2, 0.4424799189,
                      -.1256336874e1, 0.2843343470, -.2401555088e-1],
           "tau*logtau": -.1288399716e2,
           "tau*logdelta": 0.4415884023e1}

    herrig = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heavy water of Herrig"
        " (2018).",
        "__doi__": {
            "autor": "Herrig, S., Thol, M., Harvey, A.H., Lemmon, E.W.",
            "title": "A Reference Equation of State for Heavy Water",
            "ref": "J. Phys. Chem. Ref. Data 47(4) (2018) 043102",
            "doi": "10.1063/1.5053993"},

        "R": 8.3144598,
        "rhoc": 17.77555,
        "cp": Fi0,
        "ref": {"Tref": 276.95, "Pref": 0.660096, "ho": 0.598, "so": 0},

        "Tmin": Tt, "Tmax": 800.0, "Pmax": 100000.0, "rhomax": 65.,

        "nr1": [0.122082060e-1, 0.296956870e1, -0.379004540e1, 0.941089600,
                -0.922466250, -0.139604190e-1],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0000, 0.6555, 0.9369, 0.5610, 0.7017, 1.0672],

        "nr2": [-0.125203570, -0.555391500e1, -0.493009740e1, -0.359470240e-1,
                -0.936172870e1, -0.691835150],
        "c2": [1, 2, 2, 1, 2, 2],
        "d2": [1, 1, 3, 2, 2, 1],
        "t2": [3.9515, 4.6000, 5.1590, 0.2000, 5.4644, 2.3660],
        "gamma2": [1]*6,

        "nr3": [-0.456110600e-1, -0.224513300e1, 0.860006070e1, -0.248410420e1,
                0.164476900e2, 0.270393360e1, 0.375637470e2, -0.177607760e1,
                0.220924640e1, 0.519652000e1, 0.421097400, -0.391921100],
        "t3": [3.4553, 1.4150, 1.5745, 3.4540, 3.8106, 4.8950, 1.4300, 1.5870,
               3.7900, 2.6200, 1.9000, 4.3200],
        "d3": [1, 3, 1, 3, 1, 1, 2, 2, 2, 1, 1, 1],
        "alfa3": [0.6014, 1.4723, 1.5305, 2.4297, 1.3086, 1.3528, 3.4456,
                  1.2645, 2.5547, 1.2148, 18.738, 18.677],
        "beta3": [0.4200, 2.4318, 1.2888, 8.2710, 0.3673, 0.9504, 7.8318,
                  3.3281, 7.1753, 0.9465, 1177.0, 1167.0],
        "epsilon3": [1.8663, 0.2895, 0.5803, 0.2236, 0.6815, 0.9495, 1.1158,
                     0.1607, 0.4144, 0.9683, 0.9488, 0.9487],
        "gamma3": [1.5414, 1.3794, 1.7385, 1.3045, 2.7242, 3.5321, 2.4552,
                   0.8319, 1.3500, 2.5617, 1.0491, 1.0486]}

    hill = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heavy water of Hill et "
                    "al. (1982).",
        "__doi__": {"autor": "Hill, P.G., MacMillan, R.D.C., Lee, V.",
                    "title": "A Fundamental Equation of State for Heavy Water",
                    "ref": "J. Phys. Chem. Ref. Data 11, 1 (1982)",
                    "doi": "10.1063/1.555661"},

        "R": 8.3143565, "rhoc": 17.875414, "Tt": 276.97,
        "cp": Fi1,
        "ref": {"Tref": 276.95, "Pref": 0.660096, "ho": 0.598, "so": 0},

        "Tmin": Tt, "Tmax": 800.0, "Pmax": 100000.0, "rhomax": 65.,

        "nr1": [-0.384820628204e3, 0.108213047259e4, -0.110768260635e4,
                0.164668954246e4, -0.137959852228e4, 0.598964185629e3,
                -0.100451752702e3, 0.419192736351e3, -0.107279987867e4,
                0.653852283544e3, -0.984305985655e3, 0.845444459339e3,
                -0.376799930490e3, 0.644512590492e2, -0.214911115714e3,
                0.531113962967e3, -0.135454224420e3, 0.202814416558e3,
                -0.178293865031e3, 0.818739394970e2, -0.143312594493e2,
                0.651202383207e2, -0.171227351208e3, 0.100859921516e2,
                -0.144684680657e2, 0.128871134847e2, -0.610605957134e1,
                0.109663804408e1, -0.115734899702e2, 0.374970075409e2,
                0.897967147669, -0.527005883203e1, 0.438084681795e-1,
                0.406772082680, -0.965258571044e-2, -0.119044600379e-1],
        "d1": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
               4, 4, 4, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8],
        "t1": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
               0, 1, 2, 3, 4, 5, 6, 0, 1, 0, 1, 0, 1, 0, 1],

        "nr2": [0.382589102341e3, -0.106406466204e4, 0.105544952919e4,
                -0.157579942855e4, 0.132703387531e4, -0.579348879870e3,
                0.974163902526e2, 0.286799294226e3, -0.127543020847e4,
                0.275802674911e4, -0.381284331492e4, 0.293755152012e4,
                -0.117858249946e4, 0.186261198012e3],
        "c2": [1]*14,
        "d2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2],
        "t2": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6],
        "gamma2": [1.5394]*14}

    # Unorthodox formulation, the IAPWS 1984
    # Kestin, J., Sengers, J.V.
    # New International Formulations for the Thermodynamic Properties of Light
    # and Heavy Water
    # J. Phys. Chem. Ref. Data 15(1) (1986) 305-320
    # doi: 10.1063/1.555772

    eq = herrig, hill

    _melting = {"Tmin": 254.415, "Tmax": 315}

    @classmethod
    def _Melting_Pressure(cls, T, melting=None):
        try:
            Pm = _D2O_Melting_Pressure(T)
        except NotImplementedError:
            Pm = None
        return unidades.Pressure(Pm, "MPa")

    @classmethod
    def _Sublimation_Pressure(cls, T):
        try:
            Ps = _D2O_Sublimation_Pressure(T)
        except NotImplementedError:
            Ps = None
        return unidades.Pressure(Ps, "MPa")

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.794440e1, 0.194340e1, -0.243530e1, -0.342000e1, 0.355000e2,
              -0.302000e3],
        "t": [1.0, 1.5, 2.44, 5.3, 14, 20]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.166200e1, 0.901130e1, -0.154210e2, 0.115760e2, -0.516940e1,
              -0.236240e3],
        "t": [0.29, 1, 1.3, 1.77, 2.5, 16]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.247140e1, -0.266744e2, 0.531080e2, -0.480150e2, -0.576230e2,
              -0.371720e3],
        "t": [0.33, 1.29, 1.68, 2.09, 6.1, 17]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (2021)",
              "__doi__": {
                  "autor": "Assael, M.J., Monogenidou, S.A., Huber, M.L., "
                           "Perkins, R.A., Sengers, J.V.",
                  "title": "New International Formulation for the Viscosity "
                           "of Heavy Water",
                  "ref": "J. Phys. Chem. Ref. Data 50(3) (2021) 033102",
                  "doi": "10.1063/5.0048711"},
              "__code__": (_D2O_Viscosity, )}

    def _visco0(self, rho, T, fase):
        ref = D2O()
        ref._ref(False)
        estado = ref._eq(rho, 1.5*self.Tc)
        delta = estado["delta"]
        fird = estado["fird"]
        firdd = estado["firdd"]
        dpdrho = self.R*estado["T"]*(1+2*delta*fird+delta**2*firdd)
        drho = 1/dpdrho*1e6

        # convert ∂ρ/∂P]τ to IAPWS units, [kg/m³·MPa]
        if rho and fase:
            fase = copy(fase)
            fase.drhodP_T *= 1e6

        mu = _D2O_Viscosity(rho, T, fase=fase, drho=drho)
        return unidades.Viscosity(mu)

    _viscosity = (visco0, )

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (2021)",
               "__doi__": {
                   "autor": "Huber, M.L., Perkins, R.A., Assael, M.J., "
                            "Monogenidou, S.A., Hellmann, R., Sengers, J.V.",
                   "title": "New International Formulation for the Thermal "
                            "Conductivity of Heavy Water",
                   "ref": "J. Phys. Chem. Ref. Data 51(1) (2022) 013102",
                   "doi": "10.1063/5.0084222"},
               "__code__": (_D2O_ThCond, )}

    def _thermo0(self, rho, T, fase):
        ref = D2O()
        ref._ref(False)
        estado = ref._eq(rho, 1.5*self.Tc)
        delta = estado["delta"]
        fird = estado["fird"]
        firdd = estado["firdd"]
        dpdrho = self.R*estado["T"]*(1+2*delta*fird+delta**2*firdd)
        drho = 1/dpdrho*1e6

        # convert ∂ρ/∂P]τ to IAPWS units, [kg/m³·MPa]
        # convert cp to [J/kmol·K]
        if rho and fase:
            fase = copy(fase)
            fase.drhodP_T *= 1e6
            fase.cp /= 1e3

        k = _D2O_ThCond(rho, T, fase=fase, drho=drho)
        return unidades.ThermalConductivity(k)

    _thermal = (thermo0, )

    def _Surface(self, T):
        """Equation for the surface tension"""
        try:
            s = _D2O_Tension(T)
        except NotImplementedError:
            s = None
        return unidades.Tension(s)


class Test(TestCase):
    """Testing"""
    def test_herrig(self):
        """Table 6, pag 12"""
        fluid = D2O()

        delta = 46.26*fluid.M/fluid.rhoc
        tau = fluid.Tc/500
        ideal = fluid._phi0(fluid._constants["cp"], tau, delta)
        self.assertEqual(round(ideal["fio"], 8), 1.96352717)
        self.assertEqual(round(ideal["fiod"], 9), 0.384253134)
        self.assertEqual(round(ideal["fiodd"], 9), -0.147650471)
        self.assertEqual(round(ideal["fiot"], 8), 9.39259413)
        self.assertEqual(round(ideal["fiott"], 8), -2.09517144)
        self.assertEqual(round(ideal["fiodt"], 8), 0)
        res = fluid._Helmholtz(tau, delta)
        self.assertEqual(round(res["fir"], 8), -3.42291092)
        self.assertEqual(round(res["fird"], 9), -0.367562780)
        self.assertEqual(round(res["firdd"], 9), 0.835183806)
        self.assertEqual(round(res["firt"], 8), -5.89707436)
        self.assertEqual(round(res["firtt"], 8), -2.45187285)
        self.assertEqual(round(res["firdt"], 8), -1.13178440)

        # Table 7, Pag 12, Single phase region
        st = D2O(T=300, rhom=55.126)
        self.assertEqual(round(st.P.MPa, 10), 0.0529123711)
        self.assertEqual(round(st.cvM.JmolK, 7), 83.3839128)
        self.assertEqual(round(st.w, 5), 1403.74625)
        self.assertEqual(round(st.sM.JmolK, 8), 6.73910582)

        st = D2O(T=300, rhom=60)
        self.assertEqual(round(st.P.MPa, 6), 238.222326)
        self.assertEqual(round(st.cvM.JmolK, 7), 73.8561038)
        self.assertEqual(round(st.w, 5), 1772.79674)
        self.assertEqual(round(st.sM.JmolK, 8), 5.40117148)

        st = D2O(T=300, rhom=65)
        self.assertEqual(round(st.P.MPa, 6), 626.176781)
        self.assertEqual(round(st.cvM.JmolK, 7), 69.9125978)
        self.assertEqual(round(st.w, 5), 2296.97942)
        self.assertEqual(round(st.sM.JmolK, 8), 2.71566150)

        st = D2O(T=500, rhom=0.05)
        self.assertEqual(round(st.P.MPa, 9), 0.206052588)
        self.assertEqual(round(st.cvM.JmolK, 7), 29.4298102)
        self.assertEqual(round(st.w, 6), 514.480413)
        self.assertEqual(round(st.sM.JmolK, 6), 140.879085)

        st = D2O(T=500, rhom=0.5)
        self.assertEqual(round(st.P.MPa, 8), 1.88967446)
        self.assertEqual(round(st.cvM.JmolK, 7), 36.6460545)
        self.assertEqual(round(st.w, 6), 489.633254)
        self.assertEqual(round(st.sM.JmolK, 6), 120.227024)

        st = D2O(T=500, rhom=46.26)
        self.assertEqual(round(st.P.MPa, 8), 8.35329492)
        self.assertEqual(round(st.cvM.JmolK, 7), 62.6885994)
        self.assertEqual(round(st.w, 5), 1178.88631)
        self.assertEqual(round(st.sM.JmolK, 7), 49.5587000)

        st = D2O(T=500, rhom=50)
        self.assertEqual(round(st.P.MPa, 6), 107.462884)
        self.assertEqual(round(st.cvM.JmolK, 7), 61.7372286)
        self.assertEqual(round(st.w, 5), 1483.74868)
        self.assertEqual(round(st.sM.JmolK, 7), 46.9453826)

        st = D2O(T=500, rhom=60)
        self.assertEqual(round(st.P.MPa, 6), 721.798322)
        self.assertEqual(round(st.cvM.JmolK, 7), 57.6860681)
        self.assertEqual(round(st.w, 5), 2413.93520)
        self.assertEqual(round(st.sM.JmolK, 7), 39.3599094)

        st = D2O(T=643.8, rhom=20)
        self.assertEqual(round(st.P.MPa, 7), 21.6503820)
        self.assertEqual(round(st.cvM.JmolK, 7), 99.2661842)
        self.assertEqual(round(st.w, 6), 256.043612)
        self.assertEqual(round(st.sM.JmolK, 7), 81.7656125)

        st = D2O(T=800, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 10), 0.0664864175)
        self.assertEqual(round(st.cvM.JmolK, 7), 34.0033604)
        self.assertEqual(round(st.w, 6), 642.794634)
        self.assertEqual(round(st.sM.JmolK, 6), 169.067586)

        st = D2O(T=800, rhom=0.25)
        self.assertEqual(round(st.P.MPa, 8), 1.64466177)
        self.assertEqual(round(st.cvM.JmolK, 7), 34.4327932)
        self.assertEqual(round(st.w, 6), 639.281410)
        self.assertEqual(round(st.sM.JmolK, 6), 142.125615)

        # Table 8, Pag 13, Saturation state
        st = D2O(T=280, x=0.5)
        self.assertEqual(round(st.P.MPa, 12), 0.000823054058)
        self.assertEqual(round(st.Liquido.rhoM, 7), 55.2072786)
        self.assertEqual(round(st.Gas.rhoM, 12), 0.000353747143)
        self.assertEqual(round(st.Liquido.hM.Jmol, 6), 257.444444)
        self.assertEqual(round(st.Gas.hM.Jmol, 4), 46610.6716)
        self.assertEqual(round(st.Liquido.sM.JmolK, 9), 0.924406091)
        self.assertEqual(round(st.Gas.sM.JmolK, 6), 166.471646)

        st = D2O(T=450, x=0.5)
        self.assertEqual(round(st.P.MPa, 9), 0.921212105)
        self.assertEqual(round(st.Liquido.rhoM, 7), 49.2937575)
        self.assertEqual(round(st.Gas.rhoM, 9), 0.264075691)
        self.assertEqual(round(st.Liquido.hM.Jmol, 4), 14512.7149)
        self.assertEqual(round(st.Gas.hM.Jmol, 4), 51501.9146)
        self.assertEqual(round(st.Liquido.sM.JmolK, 7), 40.6584121)
        self.assertEqual(round(st.Gas.sM.JmolK, 6), 122.856634)

        st = D2O(T=625, x=0.5)
        self.assertEqual(round(st.P.MPa, 7), 17.2118129)
        self.assertEqual(round(st.Liquido.rhoM, 7), 30.6770554)
        self.assertEqual(round(st.Gas.rhoM, 8), 6.94443339)
        self.assertEqual(round(st.Liquido.hM.Jmol, 4), 32453.3556)
        self.assertEqual(round(st.Gas.hM.Jmol, 4), 47246.0343)
        self.assertEqual(round(st.Liquido.sM.JmolK, 7), 73.1042291)
        self.assertEqual(round(st.Gas.sM.JmolK, 7), 96.7725149)

        # Sublimation-pressure equation
        # Inline point in section 3.5, pag 9
        P = D2O._Sublimation_Pressure(245).MPa
        self.assertEqual(round(P, 13), 3.27390934e-5)

        # Melting-pressure equation
        # Inline point in section 3.4, pag 7
        P = D2O._Melting_Pressure(270).MPa
        self.assertEqual(round(P, 7), 83.7888413)

    def test_D2O_Viscosity(self):
        """Table 4, pag 11"""
        self.assertEqual(round(_D2O_Viscosity(0, 298.15)*1e6, 6), 10.035938)
        self.assertEqual(round(_D2O_Viscosity(1105, 298.15)*1e6, 4), 1092.6424)
        self.assertEqual(round(_D2O_Viscosity(1130, 298.15)*1e6, 4), 1088.3626)
        self.assertEqual(round(_D2O_Viscosity(1064, 373.15)*1e6, 5), 326.63791)
        self.assertEqual(round(_D2O_Viscosity(1, 775)*1e6, 6), 29.639474)
        self.assertEqual(round(_D2O_Viscosity(100, 775)*1e6, 6), 31.930085)
        self.assertEqual(round(_D2O_Viscosity(400, 775)*1e6, 6), 53.324172)

        # Table 5, pag 12
        fluid = D2O(rho=145, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 26.640959)
        fluid = D2O(rho=245, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 32.119967)
        fluid = D2O(rho=295, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 36.828275)
        fluid = D2O(rho=345, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 43.225016)
        fluid = D2O(rho=395, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 47.193530)
        fluid = D2O(rho=445, T=644.101)
        self.assertEqual(round(fluid.mu*1e6, 6), 50.241640)

    def test_D2O_ThCond(self):
        """Table 6, pag 10"""
        self.assertEqual(round(_D2O_ThCond(0, 298.15)*1e3, 4), 17.7498)
        self.assertEqual(round(_D2O_ThCond(1104.5, 298.15)*1e3, 3), 599.557)
        self.assertEqual(round(_D2O_ThCond(1200, 298.15)*1e3, 3), 690.421)
        self.assertEqual(round(_D2O_ThCond(0, 825)*1e3, 4), 76.4492)

        # Table 7, pag 10
        fluid = D2O(rho=1, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 4), 52.4527)
        fluid = D2O(rho=106, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 103.342)
        fluid = D2O(rho=256, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 394.612)
        fluid = D2O(rho=306, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 801.382)
        fluid = D2O(rho=356, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 1278.423)
        fluid = D2O(rho=406, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 670.833)
        fluid = D2O(rho=456, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 423.603)
        fluid = D2O(rho=750, T=644.1)
        self.assertEqual(round(fluid.k*1e3, 3), 454.846)

    def test_D2O_Tension(self):
        """Selected values from table 1 in IAPWS R5-85"""
        st = D2O(T=283.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 74.06)

        st = D2O(T=373.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 58.93)

        st = D2O(T=643.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 0.05)
