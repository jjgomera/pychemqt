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

from lib import unidades
from lib.meos import MEoS
from lib.mEoS import C3


class R125(MEoS):
    """Multiparameter equation of state for R125"""
    name = "pentafluoroethane"
    CASNumber = "354-33-6"
    formula = "CHF2CF3"
    synonym = "R125"
    _refPropName = "R125"
    _coolPropName = "R125"
    rhoc = unidades.Density(573.5822706)
    Tc = unidades.Temperature(339.173)
    Pc = unidades.Pressure(3617.7, "kPa")
    M = 120.0214  # g/mol
    Tt = unidades.Temperature(172.52)
    Tb = unidades.Temperature(225.06)
    f_acent = 0.3052
    momentoDipolar = unidades.DipoleMoment(1.563, "Debye")
    id = 1231

    CP1 = {"ao": 0,
           "an": [3.063], "pow": [0.1],
           "ao_exp": [2.303, 5.086, 7.3], "exp": [314, 756, 1707]}

    Fi2 = {"ao_log": [1, 4.911212],
           "pow": [0, 1],
           "ao_pow": [-12.9469, 8.512891],
           "ao_exp": [6.856764, 4.885985, 3.292859],
           "titao": [1.9757425, 4.7965398, 5.4932421]}

    Fi3 = {"ao_log": [1, 11.990267052],
           "pow": [0, 1, 0.25],
           "ao_pow": [13.79478971, 9.231669075, -27.87317349],
           "ao_exp": [7.028445731, 4.586635360],
           "titao": [4.907126427, 2.080818176]}

    CP3 = {"ao": 4.3987,
           "an": [0.0242728, -4.099e-6], "pow": [1, 2],
           "ao_exp": [], "exp": []}

    CP5 = {"ao": 3.111514,
           "an": [10.982115/339.33, -1.843797/339.33**2, 0.019273/339.33**3],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": []}

    # Expression in tau term, dividing by Tc in all terms
    CP4 = {"ao": 2.9072,
           "an": [11.5586/339.33, -2.1135/339.33**2],
           "pow": [1, 2],
           "ao_exp": [], "exp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Lemmon (2005)",
        "__doi__": {"autor": "Lemmon, E.W. Jacobsen, R.T",
                    "title": "A New Functional Form and New Fitting Techniques"
                             " for Equations of State with Application to "
                             "Pentafluoroethane (HFC-125)",
                    "ref": "J. Phys. Chem. Ref. Data 34(1) (2005) 69-108",
                    "doi": "10.1063/1.1797813"},

        "R": 8.314472,
        "rhoc": 4.779,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.09,

        "nr1": [5.28076, -8.67658, 0.7501127, 0.7590023, 0.01451899],
        "d1": [1, 1, 1, 2, 4],
        "t1": [0.669, 1.05, 2.75, 0.956, 1.],

        "nr2": [4.777189, -3.330988, 3.775673, -2.290919, 0.8888268,
                -0.6234864, -0.04127263, -0.08455389, -0.1308752, 0.008344962],
        "d2": [1, 1, 2, 2, 3, 4, 5, 1, 5, 1],
        "t2": [2, 2.75, 2.38, 3.37, 3.47, 2.63, 3.45, 0.72, 4.23, 0.2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
        "gamma2": [1]*10,

        "nr3": [-1.532005, -0.05883649, 0.02296658],
        "t3": [4.5, 29, 24],
        "d3": [2, 3, 5],
        "alfa3": [1]*3,
        "beta3": [1]*3,
        "exp1": [2, 3, 3],
        "exp2": [1.7, 7, 6],
        "gamma3": [0]*3,
        "epsilon3": [0]*3}

    outcalt = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-125 of Outcalt and McLinden "
                    "(1995)",
        "__doi__": {"autor": "Outcalt, S.L., McLinden, M.O.",
                    "title": "Equations of State for the Thermodynamic "
                             "Properties of R32 (Difluoromethane) and R125 "
                             "(Pentafluoroethane)",
                    "ref": "Int. J. Thermophysics 16(1) (1995) 79-89.",
                    "doi": "10.1007/BF01438959"},

        "R": 8.314471,
        "cp": CP5,
        "ref": "IIR",
        "Tc": 339.33, "Pc": 3629., "rhoc": 4.75996,

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.10,

        "b": [None, -0.523369607050e-1, 0.378761878904e1, -0.807152818990e2,
              0.115654605248e5, -0.152175619161e7, 0.597541484451e-2,
              -0.145990589966e1, -0.992338995652e3, -0.399180535687e6,
              -0.722591037504e-3, 0.358108080969, -0.108627994573e3,
              0.229821626570e-1, 0.149537670449e1, 0.911199833952e3,
              -0.254479949722, 0.102433894096e-1, -0.645583164735e1,
              0.218649963191, 0.114748721552e7, -0.118389825386e9,
              0.306539775027e5, 0.542870289406e9, 0.903502635609e3,
              -0.153646507435e6, 0.314617903718e1, 0.429297546671e6,
              0.109652021582, -0.329350271819e2, -0.338796950505e-3,
              0.384533651902, -0.491511706857e2]}

    sunaga = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Sunaga (1998)",
        "__doi__": {"autor": "Sunaga, H., Tillner-Roth, R., Sato, H., "
                             "Watanabe, K.",
                    "title": "A Thermodynamic Equation of State for "
                             "Pentafluoroethane (R-125)",
                    "ref": "Int. J. Thermophys., 19(6) (1998) 1623-1635",
                    "doi": "10.1007/BF03344914"},

        "R": 8.314471,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.09,

        "nr1": [0.12439220, 0.27922179, -1.1822597, 0.23616512, -0.01157181],
        "d1": [1, 2, 2, 3, 2],
        "t1": [-0.5, 0, 1.5, 1.5, 3],

        "nr2": [0.1225177e1, -0.2147964e1, -0.298138, 0.3391211, -0.6322995e-3,
                0.1271747e-3, 0.5026962e-5, -0.1667058, -0.733275e-1,
                -0.637878e-1, 0.683311e-5, -0.1995426e-1, 0.1260026e-1],
        "d2": [1, 1, 1, 3, 8, 10, 12, 1, 2, 4, 15, 3, 4],
        "t2": [0.5, 1, 3, 2.75, 2, -1, 1.25, 4, 4, 3, 0.25, 23, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*13}

    piao = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Piao (1998)",
        "__doi__": {"autor": "Piao, C.-C., Noguchi, M.",
                    "title": "An International Standard Equation of State for "
                             "the Thermodynamic Properties of HFC-125 "
                             "(Pentafluoroethane)",
                    "ref": "J. Phys. Chem. Ref. Data, 27(4) (1998) 775-806",
                    "doi": "10.1063/1.556021"},

        "R": 8.314471,
        "M": 120.022, "Tc": 339.165, "Pc": 3617.5, "rhoc": 568/120.022,
        "cp": CP3,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 14.11,

        "nr1": [0.85393382372e-1, -0.133260499658, 0.257817782488,
                -0.735018179542, -0.787454743426, -0.190320468891e-1,
                .388329449013, -0.631901774641, 0.623842653447, 0.109925047828,
                -0.993099630896e-1, -0.104601585904e-1, -0.769998709731e-1,
                0.149829594347e-1, 0.166640927925e-1, -0.181492321758e-2],
        "d1": [0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6],
        "t1": [1, 2, 0, 1, 3, 4, 0, 1, 3, 1, 3, 1, 2, 1, 2, 1],

        "nr2": [-.085393382372, .133260499658, .410983574575, -.45298892633],
        "d2": [0, 0, 2, 2],
        "t2": [1, 2, 1, 2],
        "c2": [2, 2, 2, 2],
        "gamma2": [1]*4}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-125 of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": CP3,
        "ref": "IIR",
        "M": 120.022, "Tc": 339.33, "rhoc": 571.3/120.022,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 14.1,

        "nr1": [0.11290996e1, -0.28349269e1, 0.29968733, 0.87282204e-1,
                0.26347747e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.61056963, 0.90073581, -0.68788457e-2, -0.44211186,
                -0.35041493e-1, -0.1269863, -0.25185874e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    astina = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Astina (2004)",
        "__doi__": {"autor": "Astina, I.M., Sato, H.",
                    "title": "A Rational Fundamental Equation of State for "
                             "Pentafluoroethane with Theoretical and "
                             "Experimental Bases",
                    "ref": "Int. J. Thermophys., 25(1) (2004) 113-131",
                    "doi": "10.1023/B:IJOT.0000022330.46522.68"},

        "R": 8.314472,
        "cp": Fi3,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 70000.0, "rhomax": 14.1,

        "nr1": [1.51628822, -1.4959805, -1.2893965, 1.47295195, -2.22976436,
                1.02082011, -9.61695881e-3, 4.14142522e-2],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4],
        "t1": [0.5, 0.75, 2.25, 0.5, 0.875, 2, 3, 0.5],

        "nr2": [1.46217490e-1, -6.56486371e-2, -9.18319727e-2, -2.90343386e-2,
                -1.74343357e-2, -8.77406498e-4, -5.10648362e-3, 3.52425947e-3,
                4.98022850e-4],
        "d2": [3, 6, 4, 2, 4, 4, 3, 5, 7],
        "t2": [4, 2, 3.25, 9.5, 4.5, 10.5, 25, 5, 28],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [7.41057508e-1, 1.13555445, -3.12563760, 9.32031442e-2,
                2.76844975e-4, -5.64403707e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [9.63969526e-3, 4.30480259e-1, 7.65668079e-1, -1.13913859e-2,
                -0.441468178, -2.00943884e-2, -1.26041587e-1, -2.32331768e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    vasserman = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-125 of Vasserman and"
                    "Fominsky (2001)",
        "__doi__": {"autor": "Vasserman A.A., Fominsky D.V.",
                    "title": "Equations of State for the Ozone-Safe "
                             "Refrigerants R32 and R125",
                    "ref": "Int. J. Thermophysics 22(4) (2001) 1089-1098",
                    "doi": "10.1023/a_1010699806169"},

        "R": 0.069275*M,
        "Tc": 339.33, "rhoc": 571.29/M,
        "cp": CP4,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [2.825627e-1, -6.854910e-1, -9.975127e-1, 6.445149e-1,
                4.024926e-1, -7.060325e-1, 2.621465e-1, 1.118047e-1,
                -1.217283e-1, -3.381068e-2, 7.828134e-3, 2.074012e-3,
                -5.912590e-5, 8.017074e-5, -5.547434e-5],
        "d1": [1, 1, 1, 2, 2, 3, 4, 5, 5, 6, 7, 7, 7, 8, 10],
        "t1": [0, 1, 3, 1, 3, 1, 2, 0, 2, 0, 1, 2, 4, 3, 1],

        "nr2": [1.203782e-1, 1.687868e1, -6.643158e1, 9.473412e1, -5.988361e1,
                1.396914e1, 6.765262e-1, -5.880934e-2, -5.648372e-1,
                1.643516e-1, 6.900750e-1, -1.066543, 3.699827e-1, -4.472888e-2,
                6.506423e-3, -8.564364e-3],
        "d2": [1, 2, 2, 2, 2, 2, 3, 4, 4, 5, 6, 6, 6, 7, 10, 10],
        "t2": [4, 1, 2, 3, 4, 5, 2, 1, 3, 4, 1, 2, 4, 5, 1, 3],
        "c2": [2]*16,
        "gamma2": [1]*16}

    eq = lemmon, outcalt, sunaga, piao, shortSpan, astina, vasserman, sun
    _PR = -0.00247

    _surface = {"sigma": [0.05252], "exp": [1.237]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.5295, 1.9026, -2.2966, -3.448],
        "t": [1, 1.5, 2.3, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.6684, 0.88415, 0.44383],
        "t": [1/3, 0.6, 2.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.8403, -7.2738, -21.89, -58.825],
        "t": [0.38, 1.22, 3.3, 6.9]}

    visco0 = {"__name__": "Huber (2006)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A.",
                  "title": "Correlation for the Viscosity of Pentafluoroethane"
                           " (R125) from the Triple Point to 500 K at "
                           "Pressures up to 60 MPa",
                  "ref": "Ind. Eng. Chem. Res. 45(12) (2006) 4447-4453",
                  "doi": "10.1021/ie051367l"},

              "eq": 1, "omega": 5,

              "ek": 237.077, "sigma": 0.5235,

              "Tref_virial": 237.077,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 339.173, "rhoref_res": 4.779*M, "muref_res": 1000,
              "nr": [5.677448e-3, -5.096662e-3],
              "tr": [2, 1],
              "dr": [2, 3],

              "CPf": 141.2564,
              "CPg1": 3.033797,
              "CPgi": [0.2992464/3.033797],
              "CPti": [-0.5]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",
              "thermo": "thermo0",

              "ek": 249, "sigma": 0.519, "omega": 5,

              "psi": [1.00907, 1.93968e-2], "psi_d": [0, 1],
              "fint": [1.18189e-3, 6.63334e-7], "fint_t": [0, 1],
              "chi": [1.2159, -5.6531e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
              "gam0": 0.0496, "qd": 7.48281e-10, "Tcref": 1.5*Tc}

    _viscosity = visco0, trnECS

    thermo0 = {"__name__": "Perkins (2006)",
               "__doi__": {
                   "autor": "Perkins, R.A., Huber, M.L.",
                   "title": "Measurement and Correlation of the Thermal "
                            "Conductivity of Pentafluoroethane (R125) from "
                            "190 K to 512 K at Pressures to 70 MPa",
                   "ref": "J. Chem. Eng. Data 51(3) (2006) 898-904",
                   "doi": "10.1021/je050372t"},

               "eq": 1,

               "Toref": 339.17, "koref": 1.,
               "no": [-0.460820e-2, 0.168688e-1, 0.488345e-2],
               "to": [0, 1, 2],

               "Tref_res": 339.17, "rhoref_res": 4.779*M, "kref_res": 1.,
               "nr": [-7.29410e-3, 1.10497e-2, 4.16339e-2, -2.89236e-2,
                      -3.11487e-2, 2.78399e-2, 1.12682e-2, -1.21100e-2,
                      -1.38322e-3, 2.11196e-3],
               "tr": [0, -1]*5,
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 5.834646e-10, "Tcref": 508.755}

    _thermal = thermo0, trnECS


class Test(TestCase):

    def test_lemmon(self):
        # Enthalpy-Entropy reference state
        st = R125(T=273.15, x=0)
        self.assertEqual(round(st.h.kJkg, 3), 200)
        self.assertEqual(round(st.s.kJkgK, 5), 1)

        # Table 12, Pag 104
        st = R125(T=200, rhom=14)
        self.assertEqual(round(st.P.MPa, 6), 42.302520)
        self.assertEqual(round(st.cvM.JmolK, 6), 85.816305)
        self.assertEqual(round(st.cpM.JmolK, 5), 123.53641)
        self.assertEqual(round(st.w, 5), 968.67194)

        st = R125(T=300, rhom=10)
        self.assertEqual(round(st.P.MPa, 7), 2.9023498)
        self.assertEqual(round(st.cvM.JmolK, 6), 99.919660)
        self.assertEqual(round(st.cpM.JmolK, 5), 164.16914)
        self.assertEqual(round(st.w, 5), 345.91235)

        st = R125(T=300, rhom=0.7)
        self.assertEqual(round(st.P.MPa, 7), 1.3245058)
        self.assertEqual(round(st.cvM.JmolK, 6), 94.823171)
        self.assertEqual(round(st.cpM.JmolK, 5), 124.96009)
        self.assertEqual(round(st.w, 5), 120.56007)

        st = R125(T=400, rhom=5.0)
        self.assertEqual(round(st.P.MPa, 7), 9.0495658)
        self.assertEqual(round(st.cvM.JmolK, 5), 114.41819)
        self.assertEqual(round(st.cpM.JmolK, 5), 198.11792)
        self.assertEqual(round(st.w, 5), 151.53060)

        st = R125(T=339.2, rhom=4.8)
        self.assertEqual(round(st.P.MPa, 7), 3.6201215)
        self.assertEqual(round(st.cvM.JmolK, 5), 130.63650)
        self.assertEqual(round(st.cpM.JmolK, 2), 274863.02)
        self.assertEqual(round(st.w, 6), 78.735928)

        # Selected point from Table 9, Appendix B, pag 104
        st = R125(T=R125.Tt, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00291)
        self.assertEqual(round(st.Liquido.rho, 1), 1690.7)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 87.130)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.49022)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.6776)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.035)
        self.assertEqual(round(st.Liquido.w, 1), 932.6)
        self.assertEqual(round(st.Gas.rho, 5), 0.24462)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 277.39)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5931)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.4984)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.5689)
        self.assertEqual(round(st.Gas.w, 1), 116.4)

        st = R125(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.67052)
        self.assertEqual(round(st.Liquido.rho, 1), 1319.8)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Liquido.s.kJkgK, 3), 1.000)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7948)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.255)
        self.assertEqual(round(st.Liquido.w, 1), 448.0)
        self.assertEqual(round(st.Gas.rho, 3), 42.070)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 333.16)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.4875)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.7240)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.8797)
        self.assertEqual(round(st.Gas.w, 1), 125.8)

        st = R125(T=66+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 3.61607)
        self.assertEqual(round(st.Liquido.rho, 2), 612.97)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 314.86)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.3602)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.078)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 1), 818.2)
        self.assertEqual(round(st.Liquido.w, 1), 78.9)
        self.assertEqual(round(st.Gas.rho, 2), 534.59)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 321.41)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.3795)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.093)
        self.assertEqual(round(st.Gas.cp.kJkgK, 1), 923.0)
        self.assertEqual(round(st.Gas.w, 1), 79.2)

    def test_piao(self):
        # Selected properties of Table 12, pag 790, saturation state
        st = R125(T=R125.Tt, x=0.5, eq="piao")
        self.assertEqual(round(st.P.kPa, 1), 3.0)
        self.assertEqual(round(st.Liquido.rho, 1), 1692.4)
        self.assertEqual(round(st.Gas.rho, 5), 0.24816)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 87.50)
        self.assertEqual(round(st.Hvap.kJkg, 2), 189.29)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 276.78)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.49158)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5888)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.0589)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.58880)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.66230)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.51827)
        self.assertEqual(round(st.Liquido.w, 2), 950.76)
        self.assertEqual(round(st.Gas.w, 2), 116.14)

        st = R125(T=-50+273.15, x=0.5, eq="piao")
        self.assertEqual(round(st.P.kPa, 1), 92.2)
        self.assertEqual(round(st.Liquido.rho, 1), 1519.8)
        self.assertEqual(round(st.Gas.rho, 4), 6.2172)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 141.55)
        self.assertEqual(round(st.Hvap.kJkg, 2), 164.80)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 306.35)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 0.76579)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5043)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.1012)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.69560)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.73077)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.60994)
        self.assertEqual(round(st.Liquido.w, 2), 674.24)
        self.assertEqual(round(st.Gas.w, 2), 127.30)

        st = R125(T=273.15, x=0.5, eq="piao")
        self.assertEqual(round(st.P.kPa, 1), 670.7)
        self.assertEqual(round(st.Liquido.rho, 1), 1320.9)
        self.assertEqual(round(st.Gas.rho, 3), 42.098)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Hvap.kJkg, 2), 133.15)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 333.15)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.4875)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.2524)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.87552)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.79706)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.72077)
        self.assertEqual(round(st.Liquido.w, 2), 447.26)
        self.assertEqual(round(st.Gas.w, 2), 125.88)

        st = R125(T=50+273.15, x=0.5, eq="piao")
        self.assertEqual(round(st.P.kPa, 1), 2536.2)
        self.assertEqual(round(st.Liquido.rho, 1), 1000.8)
        self.assertEqual(round(st.Gas.rho, 2), 196.27)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 271.03)
        self.assertEqual(round(st.Hvap.kJkg, 2), 76.47)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 347.51)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.2324)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.4691)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.9262)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.6708)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.87922)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.86250)
        self.assertEqual(round(st.Liquido.w, 2), 197.51)
        self.assertEqual(round(st.Gas.w, 2), 100.44)

        st = R125(T=64+273.15, x=0.5, eq="piao")
        self.assertEqual(round(st.P.kPa, 1), 3461.7)
        self.assertEqual(round(st.Liquido.rho, 2), 779.19)
        self.assertEqual(round(st.Gas.rho, 2), 366.82)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 301.01)
        self.assertEqual(round(st.Hvap.kJkg, 2), 36.88)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 337.88)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.3199)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.4293)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 6.1404)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 6.8345)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.92205)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.92240)
        self.assertEqual(round(st.Liquido.w, 2), 109.97)
        self.assertEqual(round(st.Gas.w, 2), 87.06)

        # Selected point of Table 13, pag 792, sinfle phase region
        st = R125(T=-90+273.15, P=50e3, eq="piao")
        self.assertEqual(round(st.rho, 1), 1656.5)
        self.assertEqual(round(st.h.kJkg, 2), 98.72)
        self.assertEqual(round(st.s.kJkgK, 5), 0.55457)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0524)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.67489)
        self.assertEqual(round(st.w, 2), 882.70)

        st = R125(T=-50+273.15, P=100e3, eq="piao")
        self.assertEqual(round(st.rho, 1), 1519.8)
        self.assertEqual(round(st.h.kJkg, 2), 141.55)
        self.assertEqual(round(st.s.kJkgK, 5), 0.76578)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1012)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.73077)
        self.assertEqual(round(st.w, 2), 674.29)

        st = R125(T=200+273.15, P=101325, eq="piao")
        self.assertEqual(round(st.rho, 4), 3.1008)
        self.assertEqual(round(st.h.kJkg, 2), 521.84)
        self.assertEqual(round(st.s.kJkgK, 4), 2.1294)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0385)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.96815)
        self.assertEqual(round(st.w, 2), 186.93)

        st = R125(T=100+273.15, P=150e3, eq="piao")
        self.assertEqual(round(st.rho, 4), 5.8639)
        self.assertEqual(round(st.h.kJkg, 2), 424.74)
        self.assertEqual(round(st.s.kJkgK, 4), 1.8725)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.89891)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.82599)
        self.assertEqual(round(st.w, 2), 165.97)

        st = R125(T=-40+273.15, P=200e3, eq="piao")
        self.assertEqual(round(st.rho, 1), 1483.9)
        self.assertEqual(round(st.h.kJkg, 2), 152.70)
        self.assertEqual(round(st.s.kJkgK, 5), 0.81436)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1229)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.74423)
        self.assertEqual(round(st.w, 2), 628.38)

        st = R125(T=273.15, P=300e3, eq="piao")
        self.assertEqual(round(st.rho, 3), 16.956)
        self.assertEqual(round(st.h.kJkg, 2), 339.48)
        self.assertEqual(round(st.s.kJkgK, 4), 1.5606)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.78927)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.69305)
        self.assertEqual(round(st.w, 2), 136.99)

        st = R125(T=200+273.15, P=400e3, eq="piao")
        self.assertEqual(round(st.rho, 3), 12.353)
        self.assertEqual(round(st.h.kJkg, 2), 520.78)
        self.assertEqual(round(st.s.kJkgK, 4), 2.0327)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0438)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.97016)
        self.assertEqual(round(st.w, 2), 185.53)

        st = R125(T=-10+273.15, P=500e3, eq="piao")
        self.assertEqual(round(st.rho, 1), 1365.4)
        self.assertEqual(round(st.h.kJkg, 2), 187.68)
        self.assertEqual(round(st.s.kJkgK, 5), 0.95452)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2114)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.78368)
        self.assertEqual(round(st.w, 2), 492.83)

        st = R125(T=-90+273.15, P=800e3, eq="piao")
        self.assertEqual(round(st.rho, 1), 1658.0)
        self.assertEqual(round(st.h.kJkg, 2), 99.01)
        self.assertEqual(round(st.s.kJkgK, 5), 0.55366)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0518)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.67519)
        self.assertEqual(round(st.w, 2), 886.64)

        st = R125(T=10+273.15, P=1e6, eq="piao")
        self.assertEqual(round(st.rho, 1), 1273.6)
        self.assertEqual(round(st.h.kJkg, 2), 212.74)
        self.assertEqual(round(st.s.kJkgK, 4), 1.0449)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.3011)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.81077)
        self.assertEqual(round(st.w, 2), 402.37)

        st = R125(T=120+273.15, P=3e6, eq="piao")
        self.assertEqual(round(st.rho, 2), 135.20)
        self.assertEqual(round(st.h.kJkg, 2), 424.87)
        self.assertEqual(round(st.s.kJkgK, 4), 1.6783)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0810)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.89894)
        self.assertEqual(round(st.w, 2), 146.62)

        st = R125(T=273.15, P=5e6, eq="piao")
        self.assertEqual(round(st.rho, 1), 1351.1)
        self.assertEqual(round(st.h.kJkg, 2), 200.30)
        self.assertEqual(round(st.s.kJkgK, 5), 0.98923)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2088)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.79411)
        self.assertEqual(round(st.w, 2), 494.19)

        st = R125(T=200+273.15, P=8e6, eq="piao")
        self.assertEqual(round(st.rho, 2), 298.28)
        self.assertEqual(round(st.h.kJkg, 2), 491.88)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7783)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2466)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0182)
        self.assertEqual(round(st.w, 2), 171.00)

        st = R125(T=-90+273.15, P=10e6, eq="piao")
        self.assertEqual(round(st.rho, 1), 1675.2)
        self.assertEqual(round(st.h.kJkg, 2), 102.56)
        self.assertEqual(round(st.s.kJkgK, 5), 0.54290)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0463)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.67920)
        self.assertEqual(round(st.w, 2), 932.18)

        st = R125(T=60+273.15, P=15e6, eq="piao")
        self.assertEqual(round(st.rho, 1), 1177.3)
        self.assertEqual(round(st.h.kJkg, 2), 275.32)
        self.assertEqual(round(st.s.kJkgK, 4), 1.2112)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2853)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.86279)
        self.assertEqual(round(st.w, 2), 386.86)

        st = R125(T=200+273.15, P=20e6, eq="piao")
        self.assertEqual(round(st.rho, 2), 722.89)
        self.assertEqual(round(st.h.kJkg, 2), 460.06)
        self.assertEqual(round(st.s.kJkgK, 4), 1.6602)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.3664)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0312)
        self.assertEqual(round(st.w, 2), 239.78)

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R125(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.0745)
        self.assertEqual(round(st.P.MPa, 3), 14.620)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.3219)

        st2 = R125(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 151.30)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.35860)

    def test_Huber(self):
        # Data in text, pag 4452
        st = R125(T=300, P=1e7)
        self.assertEqual(round(st.mu.muPas, 2), 177.38)
        st = R125(T=400, P=101325)
        self.assertEqual(round(st.mu.muPas, 3), 17.063)
