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


from unittest import TestCase

from lib import unidades
from lib.meos import MEoS
from lib.mEoS import R134a


class R1224ydZ(MEoS):
    """Multiparameter equation of state for R1224yd(Z)"""
    name = "cis-1-chloro-2,3,3,3-tetrafluoropropene"
    CASNumber = "111512-60-8"
    formula = "CHCl=CFCF3"
    synonym = "R-1224yd(Z)"
    rhoc = unidades.Density(527.127785)
    Tc = unidades.Temperature(428.69)
    Pc = unidades.Pressure(3337, "kPa")
    M = 148.4867  # g/mol
    Tt = unidades.Temperature(263)
    Tb = unidades.Temperature(287.767)
    f_acent = 0.322
    momentoDipolar = unidades.DipoleMoment(1.47, "Debye")

    CP1 = {"ao": 4,
           "ao_exp": [15.1, 0.703],
           "exp": [727, 1500]}

    akasaka = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1243zf of Akasaka "
                    "(2016)",
        "__doi__": {"autor": "Akasaka, R., Fukushima, M., Lemmon, E.W.",
                    "title": "A Helmholtz Energy Equation of State for "
                             "Cis-1-chloro-2,3,3,3-tetrafluoropropene "
                             "(R-1224yd(Z))",
                    "ref": "21st European Conference on Thermophysical "
                           "Properties, Graz, Austria, September 3-8, 2017",
                    "doi": ""},
        # data from the .fld file courteously provided by Mr. Akasaka

        "R": 8.3144598,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": 263, "Tmax": 473.15, "Pmax": 40000,

        "nr1": [0.04437748, 1.234842, -1.23132, -0.8495978, 0.161964],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.187, 1, 1, 0.36],

        "nr2": [-3.257722, -1.884409, 0.6828509, -1.470755, -0.01296567],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.564, 1.966, 0.86, 1.835, 1.134],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.217795, 1.619814, -0.7848116, -1.22362, -0.4546325],
        "d3": [1, 1, 3, 2, 3],
        "t3": [1.51, 0.922, 0.345, 1.225, 1.2],
        "alfa3": [0.713, 1.614, 1.496, 0.7278, 1.599],
        "beta3": [0.878, 1.44, 2.41, 0.86, 0.986],
        "gamma3": [1.075, 1.124, 1.084, 1.093, 1.07],
        "epsilon3": [0.763, 0.884, 0.265, 0.696, 0.681]}

    eq = (akasaka, )

    _surface = {
        "__doi__": {"autor": "Huber, M.L.",
                    "title": "Models for Viscosity, Thermal Conductivity, and "
                             "Surface Tension of Selected Pure Fluids as "
                             "Implemented in REFPROP v10.0",
                    "ref": "NISTIR 8209",
                    "doi": "10.6028/NIST.IR.8209"},
        "sigma": [0.06195], "exp": [1.277]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.5822, 1.6998, -2.6426, -3.5124],
        "t": [1.0, 1.5, 2.52, 5]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.939, 1.7931, 0.043779, -0.097119],
        "t": [0.22, 0.62, 1, 1.8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.7785, -7.1173, -19.416, -43.38],
        "t": [0.284, 1, 3, 6]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": R134a,

              "ek": 340.42, "sigma": 0.53, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 0.92,

              "psi": [0.712387, 0.186976, -0.0316058], "psi_d": [0, 1, 2],
              "fint": [0.00125], "fint_t": [0],
              "chi": [1.04], "chi_d": [0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
              "Xio": 0.214e-9, "gam0": 0.058, "qd": 0.646e-9, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


class Test(TestCase):
    """Test class"""

    def test_Huber(self):
        """Table 7, pag 266"""
        st = R1224ydZ(T=385.8, rhom=7.275)
        self.assertEqual(round(st.mu.muPas, 4), 113.9877)
        self.assertEqual(round(st.k.mWmK, 4), 53.9324)

    def test_Surface(self):
        """Table 10, pag 271"""
        self.assertEqual(round(R1224ydZ(T=385.8, x=0.5).sigma, 7), 0.0032758)


if __name__ == "__main__":
    # Trying to generate Fig.7 in paper for pseudo-testing of code
    import numpy as np
    from matplotlib import pyplot

    Ti = np.logspace(np.log10(280), np.log10(R1224ydZ.Tc-1e-5), 250)
    M = R1224ydZ.M

    Tl = []
    rhol = []
    Tg = []
    rhog = []
    for t in Ti:
        stl = R1224ydZ(T=t, x=0)
        if stl.status:
            Tl.append(stl.T)
            rhol.append(stl.rhoM)
        stg = R1224ydZ(T=t, x=1)
        if stg.status:
            Tg.append(stg.T)
            rhog.append(stg.rhoM)
    pyplot.plot(rhol, Tl, ls="-", lw=0.8, color="black")
    pyplot.plot(rhog, Tg, ls="-", lw=0.8, color="black")

    # Data from
    # Fukushima, M., Hayamizu, H., Hashimoto, M.
    # Thermodynamic Properties of Low-GWP Refrigerant for Centrifugal Chiller
    # International Refrigeration and Air Conditioning Conference. Paper 1631.
    # http://docs.lib.purdue.edu/iracc/1631

    # Table 3, saturated liquid
    Ts = np.r_[286.26, 290.26, 295.44, 300.55, 308.11, 314.18, 317.59, 325.87,
               428.50, 422.56]
    rhols = np.r_[1395, 1384, 1371, 1354, 1336, 1320, 1304, 1281, 590, 789]/M
    pyplot.plot(rhols, Ts, ls="", marker="o", mec="black", mfc="black")

    # Table 1, critical region
    Tc = np.r_[427.71, 428.25, 428.66, 429.19, 429.17, 429.02, 427.66]
    rhoc = np.r_[371.20, 414.70, 417.01, 469.29, 530.18, 600.00, 674.13]/M
    pyplot.plot(rhoc, Tc, ls="", marker="s", mec="black", mfc="#00000000")

    # Data from
    # Sakoda, N., Higashi, Y.
    # Measurements of PvT Properties, Vapor Pressures, Saturated Densities,
    # and Critical Parameters for cis-1-Chloro-2,3,3,3-tetrafluoropropene
    # (R1224yd(Z))
    # J. Chem. Eng. Data 64(9) (2019) 3983-3987
    # doi: 10.1021/acs.jced.9b00374

    # Table 3, Visual observation of the V-L meniscus disappearance
    Ts = np.r_[424.261, 427.532, 428.216, 428.671, 428.688, 428.690, 427.089,
               425.294, 415.835, 408.745]
    rhos = np.r_[297.0, 387.5, 415.1, 490.3, 525.2, 541.7, 685.5, 734.3,
                 867.9, 930.0]/M
    pyplot.plot(rhos, Ts, ls="", marker="o", mec="blue", mfc="blue")

    # Table 4, Intersections with the vapor pressure curves
    Ts = np.r_[325.814, 355.441, 361.676, 387.362, 401.187, 374.157, 351.037]
    rhos = np.r_[21.8, 47.9, 55.4, 103.5, 970.7, 1107.5, 1197.8]/M
    pyplot.plot(rhos, Ts, ls="", marker="^", mec="blue", mfc="#00000000")

    pyplot.xlabel("ρ (mol/l)")
    pyplot.ylabel("T (K)")

    pyplot.show()
