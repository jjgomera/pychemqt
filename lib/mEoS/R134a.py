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


class R134a(MEoS):
    """Multiparameter equation of state for R134a"""
    name = "1,1,1,2-tetrafluoroethane"
    CASNumber = "811-97-2"
    formula = "CF3CH2F"
    synonym = "R134a"
    _refPropName = "R134A"
    rhoc = unidades.Density(511.9)
    Tc = unidades.Temperature(374.21)
    Pc = unidades.Pressure(4059.28, "kPa")
    M = 102.032  # g/mol
    Tt = unidades.Temperature(169.85)
    Tb = unidades.Temperature(247.076)
    f_acent = 0.32684
    momentoDipolar = unidades.DipoleMoment(2.058, "Debye")
    id = 1235

    Fi1 = {"R": 8.314471,
           "ao_log": [1, -1.629789],
           "pow": [0, 1, -0.5, -0.75],
           "ao_pow": [-1.019535, 9.047135, -9.723916, -3.92717],
           "ao_exp": [], "titao": []}

    Fi2 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.25],
           "ao_pow": [10.78497786, 8.612977410, -24.37548384],
           "ao_exp": [7.451784998, -4.239239505, 6.445739825],
           "titao": [-4.103830338, -2.561528683, -2.084607363]}

    CP2 = {"ao": 19.4006,
           "an": [0.258531, -1.29665e-4], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Tillner-Roth "
                    "(1994).",
        "__doi__": {"autor": "Tillner-Roth, R., Baehr, H.D.",
                    "title": "An International Standard Formulation for the "
                             "Thermodynamic Properties of 1,1,1,2-"
                             "Tetrafluoroethane (HFC-134a) for Temperatures "
                             "from 170 K to 455 K at Pressures up to 70 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 23(5) (1994) 657-729",
                    "doi": "10.1063/1.555958"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",
        "Tc": 374.18, "rhoc": 508/M,

        "Tmin": Tt, "Tmax": 465.0, "Pmax": 70000.0, "rhomax": 15.60,
        "Pmin": 0.3896, "rhomin": 15.5942,

        "nr1": [0.5586817e-1, 0.498223, 0.2458698e-1, 0.8570145e-3,
                0.4788584e-3, -0.1800808e1, 0.2671641, -0.4781652e-1],
        "d1": [2, 1, 3, 6, 6, 1, 1, 2],
        "t1": [-0.5, 0, 0, 0, 1.5, 1.5, 2, 2],

        "nr2": [0.1423987e-1, 0.3324062, -0.7485907e-2, 0.1017263e-3,
                -0.5184567, -0.8692288e-1, 0.2057144, -0.5000457e-2,
                0.4603262e-3, -0.3497836e-2, 0.6995038e-2, -0.1452184e-1,
                -0.1285458e-3],
        "d2": [5, 2, 2, 4, 1, 4, 1, 2, 4, 1, 5, 3, 10],
        "t2": [1, 3, 5, 1, 5, 5, 6, 10, 10, 10, 18, 22, 50],
        "c2": [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4],
        "gamma2": [1]*13}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-134a of Huber and McLinden "
                    "(1992)",
        "__doi__": {"autor": "Huber, M.L., McLinden, M.O.",
                    "title": "Thermodynamic Properties of R134a "
                             "(1,1,1,2-tetrafluoroethane)",
                    "ref": "International Refrigeration and Air Condictioning "
                           "Conference, Paper 184, pag. 453-462, 1992.",
                    "doi": ""},

        "R": 8.314471,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 15.60,
        "Pmin": 0.3922, "rhomin": 15.60,

        "b": [None, 0.965209362217e-1, -0.401824768889e1, 0.395239532858e2,
              0.134532868960e4, -0.139439741347e7, -0.309281355175e-2,
              0.292381512283e1, -0.165146613555e4, 0.150706003118e7,
              0.534973948313e-4, 0.543933317622, -0.211326049762e3,
              -0.268191203847e-1, -0.541067125950, -0.851731779398e3,
              0.205188253646, -0.733050188093e-2, 0.380655963862e1,
              -0.105832087589, -0.679243084424e6, -0.126998378601e9,
              -0.426234431829e5, 0.101973338234e10, -0.186699526782e3,
              -0.933426323419e5, -0.571735208963e1, -0.176762738787e6,
              -0.397282752308e-1, 0.143016844796e2, 0.803085294260e-4,
              -0.171959073552, 0.226238385661e1]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-134a of Span and "
                    "Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",
        "Tc": 374.18, "rhoc": 508/M,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 15.6,
        "Pmin": 0.38818, "rhomin": 15.588,

        "nr1": [0.106631890000e1, -0.244959700000e1, 0.446457180000e-1,
                0.756568840000e-1, 0.206520890000e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.42006912, 0.76739111, 0.17897427e-2, -0.36219746,
                -0.6780937e-1, -0.10616419, -0.18185791e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    astina = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Astina (2004)",
        "__doi__": {"autor": "Astina, I.M. Sato, H.",
                    "title": "A Fundamental Equation of State for "
                             "1,1,1,2-tetrafluoroethane with an Intermolecular"
                             " Potential Energy Background and Relialbe "
                             "Ideal-Gas Properties",
                    "ref": "Fluid Phase Equilib., 221 (2004) 103-111",
                    "doi": "10.1016/j.fluid.2004.03.004"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 460.0, "Pmax": 70000.0, "rhomax": 15.58,
        "Pmin": 0.327, "rhomin": 15.58,

        "nr1": [1.832124209, -2.940698861, 5.156071823e-1, 2.756965911e-1,
                1.225264939, -6.486749497e-1, -9.286738053e-1, 3.920381291e-1,
                1.056692108e-1],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4],
        "t1": [0.5, 1.125, 3.25, 0.5, 1.875, 2.75, 1.625, 2.125, 1.125],

        "nr2": [-7.586523371e-1, -1.198140136, -2.878260390e-1,
                -9.723032379e-2, 5.307113358e-2, -4.681610582e-2,
                -9.604697902e-3, 6.668035048e-3, 2.361266290e-3],
        "d2": [1, 2, 3, 2, 3, 4, 4, 5, 6],
        "t2": [3.75, 1.5, 0.75, 9, 8.5, 5.5, 32, 23, 31],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.08605179, 1.03772416, -2.92069735, 9.15573346e-2,
                2.40541430e-4, -2.00239570e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.61424796e-2, -2.15499979e-1, 3.11819936e-1, 1.12867938e-3,
                -0.283454532, -4.21157950e-2, -8.08314045e-2, -1.59762784e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Huber and Eli"
                    "(1992).",
        "__doi__": {"autor": "Huber, M.L., Ely, J.F.",
                    "title": "An Equation of State Formulation of the "
                             "Thermodynamic Properties of R134a ("
                             "1,1,1,2-Tetrafluoroethane)",
                    "ref": "Int. J. Refrig. 15(6) (1992) 393-400",
                    "doi": "10.1016/0140-7007(92)90024-o"},

        "R": 8.314471,
        "cp": CP2,
        "ref": "IIR",
        "Tc": 374.179, "rhoc": 5.0308, "Pc": 4058.59,

        "Tmin": Tt, "Tmax": 465.0, "Pmax": 70000.0, "rhomax": 15.60,
        "Pmin": 0.3896, "rhomin": 15.5942,

        "nr1": [6.81716385385e-1, -2.35124614105, 6.70216482859e-1,
                -3.07204611902e-2, 3.74529023556e-1, -1.57205367415e-1,
                6.52988383109e-2, -5.10116156742e-2, -5.69183659026e-2,
                6.45310700471e-4, 1.02593424592e-3, 6.77375367275e-7,
                -1.92870222869e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-4.95254825047e-1, 1.28070070661e-1, 2.76305386558e-1,
                -1.53983381830e-1, -2.11838190838e-1, -2.39896004684e-2,
                -3.27379569918e-3, -9.27516738026e-4, -2.09645193939e-2,
                2.14330093737e-3, -5.41732277806e-4, 3.47165872395e-3,
                4.91210193371e-2, -3.69286578727e-2, -6.94084047023e-2,
                4.73399474790e-2, 6.55276251860e-1, -6.87628059906e-1,
                4.30311999742e-2],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2]*11+[4]*8,
        "gamma2": [1]*19}

    eq = tillner, MBWR, shortSpan, astina, sun, huber
    _PR = 0.001032

    _surface = {"sigma": [0.05801], "exp": [1.241]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77513e1, 0.29263e1, -0.26622e1, -0.39711e1],
        "exp": [1.0, 1.5, 1.9, 4.25]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.12449e2, -0.41023e2, 0.73641e2, -0.64635e2, 0.22551e2],
        "exp": [0.5, 0.7, 0.9, 1.1, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.29174e1, -0.72542e1, -0.23306e2, 0.59840e1, -0.71821e2],
        "exp": [0.383, 1.21, 3.3, 5.6, 7.0]}

    visco0 = {"__name__": "Huber (2003)",
              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": 1, "omega": 1,

              "M": 102.031, "ek": 299.363, "sigma": 0.46893,
              "n_chapman": 0.021357,
              "collision": [0.355404, -0.464337, 0.257353e-1],

              "Tref_virial": 299.363,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 374.21, "rhoref_res": 511.9, "muref_res": 1000,
              "nr": [-0.206900719e-1, 0.356029549e-3, 0.211101816e-2,
                     0.139601415e-1, -0.456435020e-2, -0.351593275e-2],
              "tr": [0, 6, 2, 0.5, -2, 0],
              "dr": [1, 2, 2, 2, 2, 3],

              "special": "_closePacked"}

    def _closePacked(self, rho, T, fase):
        """Special closed packed term"""
        c10 = 3.163695636
        c9 = 0.100035295
        c8 = -0.890173375e-1
        c7 = 0.214763320

        tau = T/374.21
        delta = rho/511.9
        delta0 = c10/(1+c8*tau+c9*tau**2)
        muCP = c7/(delta0-delta) - c7/delta0
        return muCP * 1000

    visco1 = {"__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 374.21,
              "no": [31.2515, -89.6122, 73.0823],
              "to": [0, 0.25, 0.5],

              "a": [1.07271e-4, -4.41655e-5, 0.0],
              "b": [1.66457e-4, -4.80293e-5, 0.0],
              "c": [8.08333e-5, -4.90360e-5, 0.0],
              "A": [-2.12476e-8, 2.81647e-9, 0.0],
              "B": [1.35594e-8, 0.0, 3.17550e-10],
              "C": [0.0, 4.81769e-7, -1.17149e-7]}

    _viscosity = visco0, visco1,

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2000)",
               "__doi__": {"autor": "Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L.",
                           "title": "Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)",
                           "ref": "NIST Interagency/Internal Report (NISTIR) - 6605",
                           "doi": ""},

               "Tref": 1., "kref": 1.,
               "no": [-1.05248e-2, 8.00982e-5],
               "co": [0, 1],

               "Trefb": 339.173, "rhorefb": 5.049886, "krefb": 2.055e-3,
               "nb": [1.836526, 5.126143, -1.436883, 6.261441e-1],
               "tb": [0]*4,
               "db": [1, 2, 3, 4],
               "cb": [0]*4,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.285356e-10, "Tcref": 561.411}

    _thermal = thermo0,


class Test(TestCase):

    def test_tillner(self):
        # The tables have fairly values that differ in the last decimal place,
        # always a value below the one obtained, perhaps a problem of rounding

        # Selected point from Table 8, Pag 696, saturation state
        st = R134a(T=169.85, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00039)
        self.assertEqual(round(st.Liquido.rho, 1), 1591.1)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 71.455)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.4126)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.7922)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.1838)
        self.assertEqual(round(st.Liquido.w, 1), 1120.0)
        self.assertEqual(round(st.Gas.rho, 5), 0.02817)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 334.94)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.9639)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.5030)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.5853)
        self.assertEqual(round(st.Gas.w, 2), 126.79)

        st = R134a(T=200, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.00631)
        self.assertEqual(round(st.Liquido.rho, 1), 1510.5)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 107.39)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.6073)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8016)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.2058)
        self.assertEqual(round(st.Liquido.w, 2), 967.61)
        self.assertEqual(round(st.Gas.rho, 5), 0.38977)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 353.06)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.8356)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.5732)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.6586)
        self.assertEqual(round(st.Gas.w, 2), 135.98)

        st = R134a(T=250, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.11561)
        self.assertEqual(round(st.Liquido.rho, 1), 1367.9)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 169.57)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.8841)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8515)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.2865)
        self.assertEqual(round(st.Liquido.w, 2), 728.39)
        self.assertEqual(round(st.Gas.rho, 4), 5.9546)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 384.60)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.7443)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.6962)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.8044)
        self.assertEqual(round(st.Gas.w, 2), 145.98)

        st = R134a(P=1e5, x=0.5)
        self.assertEqual(round(st.T, 2), 246.79)
        self.assertEqual(round(st.Liquido.rho, 1), 1377.5)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 165.44)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.8676)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.8478)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.2800)
        self.assertEqual(round(st.Liquido.w, 2), 743.31)
        self.assertEqual(round(st.Gas.rho, 4), 5.1932)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 382.60)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.7475)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.6876)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 0.7932)
        self.assertEqual(round(st.Gas.w, 2), 145.63)

        st = R134a(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 5), 0.70282)
        self.assertEqual(round(st.Liquido.rho, 1), 1199.7)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 237.19)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.1287)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 0.9144)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.4324)
        self.assertEqual(round(st.Liquido.w, 2), 497.89)
        self.assertEqual(round(st.Gas.rho, 3), 34.193)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 413.27)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.7156)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 0.8426)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.0438)
        self.assertEqual(round(st.Gas.w, 2), 143.88)

        st = R134a(T=350, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 2.4611)
        self.assertEqual(round(st.Liquido.rho, 2), 951.32)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 316.50)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.3674)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.0037)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.9614)
        self.assertEqual(round(st.Liquido.w, 2), 254.06)
        self.assertEqual(round(st.Gas.rho, 2), 140.99)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 429.03)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.6889)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.0301)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.8494)
        self.assertEqual(round(st.Gas.w, 2), 120.33)

        st = R134a(T=374, x=0.5)
        self.assertEqual(round(st.P.MPa, 4), 4.0416)
        self.assertEqual(round(st.Liquido.rho, 2), 587.91)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 380.86)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.5387)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.2121)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 101.66)
        self.assertEqual(round(st.Liquido.w, 3), 92.401)
        self.assertEqual(round(st.Gas.rho, 2), 434.06)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 399.51)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.5886)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.2409)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 137.24)
        self.assertEqual(round(st.Gas.w, 3), 91.390)

        # Selected point from Table 9, Pag 700, single phase region
        st = R134a(T=205, P=1e4)
        self.assertEqual(round(st.rho, 1), 1496.8)
        self.assertEqual(round(st.h.kJkg, 3), 113.44)
        self.assertEqual(round(st.s.kJkgK, 4), 0.6372)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8055)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2120)
        self.assertEqual(round(st.w, 2), 942.87)

        st = R134a(T=325, P=2e4)
        self.assertEqual(round(st.rho, 5), 0.75738)
        self.assertEqual(round(st.h.kJkg, 2), 449.09)
        self.assertEqual(round(st.s.kJkgK, 4), 2.1106)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8000)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.8828)
        self.assertEqual(round(st.w, 2), 170.45)

        st = R134a(T=460, P=6e4)
        self.assertEqual(round(st.rho, 4), 1.6045)
        self.assertEqual(round(st.h.kJkg, 2), 582.87)
        self.assertEqual(round(st.s.kJkgK, 4), 2.3634)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0172)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0998)
        self.assertEqual(round(st.w, 2), 200.83)

        st = R134a(T=170, P=1e5)
        self.assertEqual(round(st.rho, 1), 1590.8)
        self.assertEqual(round(st.h.kJkg, 3), 71.678)
        self.assertEqual(round(st.s.kJkgK, 4), 0.4136)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7922)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1837)
        self.assertEqual(round(st.w, 1), 1119.6)

        st = R134a(T=265, P=2e5)
        self.assertEqual(round(st.rho, 4), 9.9171)
        self.assertEqual(round(st.h.kJkg, 2), 394.26)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7396)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.7324)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.8523)
        self.assertEqual(round(st.w, 2), 147.61)

        st = R134a(T=240, P=5e5)
        self.assertEqual(round(st.rho, 1), 1398.8)
        self.assertEqual(round(st.h.kJkg, 2), 156.94)
        self.assertEqual(round(st.s.kJkgK, 4), 0.8314)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8403)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2659)
        self.assertEqual(round(st.w, 2), 777.64)

        st = R134a(T=310, P=1e6)
        self.assertEqual(round(st.rho, 1), 1160.4)
        self.assertEqual(round(st.h.kJkg, 2), 251.72)
        self.assertEqual(round(st.s.kJkgK, 4), 1.1755)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9288)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.4795)
        self.assertEqual(round(st.w, 2), 452.15)

        st = R134a(T=345, P=2e6)
        self.assertEqual(round(st.rho, 2), 102.38)
        self.assertEqual(round(st.h.kJkg, 2), 434.73)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7164)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.9701)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.4224)
        self.assertEqual(round(st.w, 2), 131.78)

        st = R134a(T=460, P=5e6)
        self.assertEqual(round(st.rho, 2), 169.51)
        self.assertEqual(round(st.h.kJkg, 2), 549.57)
        self.assertEqual(round(st.s.kJkgK, 4), 1.9474)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.0760)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.3429)
        self.assertEqual(round(st.w, 2), 169.93)

        st = R134a(T=270, P=1e7)
        self.assertEqual(round(st.rho, 1), 1337.4)
        self.assertEqual(round(st.h.kJkg, 2), 198.54)
        self.assertEqual(round(st.s.kJkgK, 4), 0.9675)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8750)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2947)
        self.assertEqual(round(st.w, 2), 707.43)

        st = R134a(T=460, P=2e7)
        self.assertEqual(round(st.rho, 2), 767.68)
        self.assertEqual(round(st.h.kJkg, 2), 475.99)
        self.assertEqual(round(st.s.kJkgK, 4), 1.7125)
        self.assertEqual(round(st.cv.kJkgK, 4), 1.1073)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.6079)
        self.assertEqual(round(st.w, 2), 281.20)

        st = R134a(T=185, P=5e7)
        self.assertEqual(round(st.rho, 1), 1608.6)
        self.assertEqual(round(st.h.kJkg, 2), 112.12)
        self.assertEqual(round(st.s.kJkgK, 4), 0.4657)
        self.assertEqual(round(st.cv.kJkgK, 4), 0.8100)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1623)
        self.assertEqual(round(st.w, 1), 1208.3)

    def test_MBWR(self):
        # TODO: Add testing when fix MBWR equations
        pass

        # Selected points from Table 9, Pag 459, saturation state
        # st = R134a(T=-100+273.15, x=0.5)
        # self.assertEqual(round(st.P.MPa, 5), 0.00056)
        # self.assertEqual(round(st.Liquido.rho, 1), 1581.9)
        # self.assertEqual(round(st.Gas.rho, 3), 0.040)
        # self.assertEqual(round(st.Liquido.h.kJkg, 2), 75.71)
        # self.assertEqual(round(st.Gas.h.kJkg, 2), 337.00)
        # self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.4366)
        # self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.9456)
        # self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.168)
        # self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.592)
        # self.assertEqual(round(st.Liquido.w, 0), 1111)
        # self.assertEqual(round(st.Gas.w, 0), 128)

        # st = R134a(T=+273.15, x=0.5)
        # self.assertEqual(round(st.P.MPa, 5), )
        # self.assertEqual(round(st.Liquido.rho, 1), )
        # self.assertEqual(round(st.Gas.rho, 3), )
        # self.assertEqual(round(st.Liquido.h.kJkg, 2), )
        # self.assertEqual(round(st.Gas.h.kJkg, 2), )
        # self.assertEqual(round(st.Liquido.s.kJkgK, 4), )
        # self.assertEqual(round(st.Gas.s.kJkgK, 4), )
        # self.assertEqual(round(st.Liquido.cp.kJkgK, 3), )
        # self.assertEqual(round(st.Gas.cp.kJkgK, 3), )
        # self.assertEqual(round(st.Liquido.w, 0), )
        # self.assertEqual(round(st.Gas.w, 0), )

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R134a(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.1577)
        self.assertEqual(round(st.P.MPa, 3), 14.656)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.6129)

        st2 = R134a(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 181.97)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.41386)
