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


from lib import unidades
from lib.meos import MEoS
from lib.mEoS import C3


class R13(MEoS):
    """Multiparameter equation of state for R13"""
    name = "chlorotrifluoromethane"
    CASNumber = "75-72-9"
    formula = "CClF3"
    synonym = "R13"
    _refPropName = "R13"
    _coolPropName = "R13"
    rhoc = unidades.Density(582.88)
    Tc = unidades.Temperature(302.0)
    Pc = unidades.Pressure(3879.0, "kPa")
    M = 104.459  # g/mol
    Tt = unidades.Temperature(92.0)
    Tb = unidades.Temperature(191.67)
    f_acent = 0.1723
    momentoDipolar = unidades.DipoleMoment(0.51, "Debye")
    id = 215

    # Cp⁰ coefficient for Tr term, so divide by Tc
    CP1 = {"ao": 1.86012334,
           "an": [8.07314520/Tc, -1.87713639/Tc**2, 3.17242858e-2/Tc**3],
           "pow": [1, 2, 3]}

    CP2 = {"ao": 2.4766458,
           "an": [0.018074269, 2.1945535e-5, -8.5810657e-8, 6.3199171e-11],
           "pow": [1, 2, 3, 4]}

    magee = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-13 of Magee et al. (2000)",
        "__doi__": {"autor": "Magee, J.W., Outcalt, S.L., Ely, J.F.",
                    "title": "Molar Heat Capacity Cv, Vapor Pressure, and "
                             "(p, ρ, T) Measurements from 92 to 350 K at "
                             "Pressures to 35 MPa and a New Equation of State "
                             "for Chlorotrifluoromethane (R13)",
                    "ref": "Int. J. Thermophys., 21(5):1097-1121, 2000.",
                    "doi": "10.1023/A:1026446004383"},

        "R": 8.314471,
        "Tc": 302, "Pc": 3879., "rhoc": 5.58, "M": 104.459,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 403.0, "Pmax": 35000.0, "rhomax": 17.85,

        "b": [None, 0.427710490378e-2, 0.106603397093e1, -0.383065097813e2,
              0.661580211522e4, -0.800160780370e6, -0.406405755462e-2,
              0.561380767634e1, -0.247694806929e4, -0.639834580892e5,
              0.198818486764e-3, -0.206916891385, 0.749317872337e2,
              -0.431471653965e-2, 0.181741326553e1, -0.206066849491e2,
              -0.136681208829, 0.260496240940e-2, 0.287244312242,
              -0.105459756169e-1, 0.582404815872e6, -0.455721947029e8,
              0.114174177352e5, 0.265590236008e6, 0.135249873550e3,
              0.128289104267e4, 0.800900540368, -0.703307137789e4,
              0.235567665577e-2, 0.131830636112e1, -0.115187941781e-4,
              0.564530387616e-2, 0.336242130107]}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-13 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},
        "R": 8.31451,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 50000.0, "rhomax": 17.699806,

        "nr1": [-0.628346559920, 0.792797111341, -0.134038992692,
                0.761143010172, -0.194465098795e1, 0.940938700406,
                -0.108107050239e1, 0.117501564976, 0.228305167217,
                -0.403338888789, 0.375585713420, -0.617543677315e-1,
                0.170326226881, 0.536612457231e-1, -0.151603010301,
                0.252033265074e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.628346559920, -0.792797111341, 0.134038992692,
                -0.399863840975e-1, 0.436410910529, -0.448724904991],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.98230055]*6}

    eq = magee, platzer
    _PR = [-0.2395, -15.5913]

    _surface = {"sigma": [0.05045], "exp": [1.269]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.69311e1, 0.18281e1, -0.21901e1, -0.38177e1, 0.20803e1],
        "t": [1.0, 1.5, 2.5, 6.0, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.95469e1, -0.24017e2, 0.33365e2, -0.26837e2, 0.10638e2],
        "t": [0.51, 0.72, 0.94, 1.2, 1.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.31949e1, -0.73425e1, -0.21966e2, -0.51459e2, -0.85359e2],
        "t": [0.414, 1.41, 3.7, 7.7, 15.0]}

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

              "ek": 204.0, "sigma": 0.4971, "omega": 5,

              "psi": [0.976177, 1.48047e-2], "psi_d": [0, 1],
              "fint": [1.07447e-3, 6.42373e-7], "fint_t": [0, 1],
              "chi": [1.1394, -3.65562e-2], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
              "gam0": 0.0496, "qd": 3.49636e-10, "Tcref": 453.0}

    _viscosity = (trnECS, )
    _thermal = (trnECS, )


if __name__ == "__main__":
    st = R13(P=256, T=92)
    print(st.status, st.msg)
    print(st.P)
    print(st.v, st.u)
