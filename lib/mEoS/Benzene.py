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


from lib.meos import MEoS
from lib import unidades


class Benzene(MEoS):
    """Multiparameter equation of state for benzene"""
    name = "benzene"
    CASNumber = "71-43-2"
    formula = "C6H6"
    synonym = ""
    rhoc = unidades.Density(304.79239968)
    Tc = unidades.Temperature(562.02)
    Pc = unidades.Pressure(4906.3, "kPa")
    M = 78.11184  # g/mol
    Tt = unidades.Temperature(278.674)
    Tb = unidades.Temperature(353.22)
    f_acent = 0.211
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 40

    Fi1 = {"ao_log": [1, 2.94645],
           "pow": [0, 1],
           "ao_pow": [-0.6740687105, 2.5560188958],
           "ao_exp": [7.36374, 18.6490, 4.01834],
           "titao": [4116/Tc, 1511/Tc, 630/Tc]}

    CP2 = {"ao": -0.478176/8.3143*78.108,
           "an": [0.618649e-2/8.3143*78.108, -0.380363e-5/8.3143*78.108,
                  0.699648e-9/8.3143*78.108, 0.42661e-13/8.3143*78.108],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Thol et al. (2013).",
        "__doi__": {"autor": "Thol M., Lemmon E.W., Span R.",
                    "title": "Equation of state for benzene for temperatures from the melting line up to 750 K and pressures up to 500 MPa (final)",
                    "ref": "to be published, 2013",
                    "doi": ""},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750., "Pmax": 500000.0, "rhomax": 11.45,
        "Pmin": 4.78, "rhomin": 11.45,

        "nr1": [0.03512459, 2.2338, -3.10542612, -0.577233, 0.25101],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.29, 0.696, 1.212, 0.595],

        "nr2": [-0.705518, -0.139648, 0.83494, -0.331456, -0.0279953],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.51, 3.96, 1.24, 1.83, 0.82],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1.]*5,

        "nr3": [0.7099766, -0.3732185, -0.0629985, -0.803041],
        "d3": [1, 1, 3, 3],
        "t3": [0.57, 2.04, 3.2, 0.78],
        "alfa3": [1.032, 1.423, 1.071, 14.2],
        "beta3": [1.864, 1.766, 1.825, 297.9],
        "gamma3": [1.118, 0.639, 0.654, 1.164],
        "epsilon3": [0.729, 0.907, 0.765, 0.87],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Thol et al. (2012).",
        "__doi__": {"autor": "Thol M., Lemmon E.W., Span R.",
                    "title": "Equation of state for benzene for temperatures from the melting line up to 750 K and pressures up to 500 MPa",
                    "ref": "High Temperatures-High Pressures 01/2012; 41:81.",
                    "doi": ""},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750., "Pmax": 500000.0, "rhomax": 11.45,
        "Pmin": 4.78, "rhomin": 11.45,

        "nr1": [0.3513062e-1, 0.2229707e1, -0.3100459e1, -0.5763224, 0.2504179],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.744, 1.174, 0.68],

        "nr2": [-0.7049091, -0.1393433, 0.8319673, -0.3310741, -0.2793578e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.5, 3.67, 1.26, 2.6, 0.95],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1.]*5,

        "nr3": [0.7087408, -0.3723906, -0.6267414e-1, -0.8629500],
        "d3": [1, 1, 3, 3],
        "t3": [1, 2.47, 3.35, 0.75],
        "alfa3": [1.032, 1.423, 1.071, 14.35],
        "beta3": [1.867, 1.766, 1.824, 297.5],
        "gamma3": [1.1180, 0.6392, 0.6536, 1.1640],
        "epsilon3": [0.7289, 0.9074, 0.7655, 0.8711],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Polt et al. (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 278.7, "Tmax": 635.0, "Pmax": 78000.0, "rhomax": 11.45,
        "Pmin": 6.0329, "rhomin": 11.385,

        "nr1": [-0.918572178424, 0.155357491575e1, -0.356149241161,
                0.817273664265, -0.331303917534e1, 0.335336626528e1,
                -0.256976312022e1, 0.427304812515, 0.406483484297,
                -0.329744378187, 0.208907540720, 0.777471199254e-1,
                -0.202621443063, -0.148580350700e-1, 0.503167715817e-1,
                0.293012717053e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.918572178424, -0.155357491575e1, 0.356149241161,
                -0.447029533153e-1, 0.957712367542, -0.114688433057e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.95481]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.76284970, 1.02610647, -3.74263321, 9.57682041e-2,
                2.59179321e-4, -1.03082188e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [1.07359246e-1, -1.12562310e-1, 3.18737987e-1, -3.07549016e-2,
                -3.25082386e-1, 2.28099159e-2, -7.07431076e-2, -1.96809158e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2, helmholtz3, helmholtz4

    _surface = {"sigma": [0.07298, -0.0007802, -0.0001756],
                "exp": [1.232, 0.8635, 0.3065]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71661e1, 0.21551e1, -0.20297e1, -0.40668e1, 0.38092],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18160e2, -0.56879e2, 0.87478e2, -0.64365e2, 0.18500e2],
        "exp": [0.534, 0.686, 0.84, 1.0, 1.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31147e1, -0.46689e1, -0.16161e2, -0.14650e3, 0.51887e3, -0.82772e3],
        "exp": [0.419, 1.12, 2.8, 7.3, 10., 12.]}

    visco0 = {"eq": 5, "omega": 3,
              "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E",
                          "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties",
                          "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                          "doi": "10.1021/ie00076a024"},
              "__name__": "Chung (1988)",
              "w": 0.5693, "mur": 0.3209, "k": 0.0642}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Assael (2012)",
               "__doi__": {"autor": "Assael, M.J., Mihailidou, E., Huber, M.L. and Perkins, R.A.",
                           "title": "Reference Correlation of the Thermal Conductivity of Benzene from the Triple Point to 725 K and up to 500 MPa",
                           "ref": "J. Phys. Chem. Ref. Data 41, 043102 (2012)",
                           "doi": "10.1063/1.4755781"},

               "Tref": 1., "kref": 1e-3,
               "no": [56991.07, -521.44, 1.5449],
               "co": [0, 1, 2],
               "noden": [562.02, 9.714, 0.0026102],
               "coden": [0, 1, 2],

               "Trefb": Tc, "rhorefb": 2.3153, "krefb": 1e-3,
               "nb": [.282489e-1, -.773415e-1, .714001e-1, -.236798e-1,
                      .300875e-2, -.119268e-1, .833389e-1, -.898176e-1,
                      .363025e-1, -.490052e-2],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.216-9, "gam0": 0.0569, "qd": 0.62e-9, "Tcref": 843}

    _thermal = thermo0,
