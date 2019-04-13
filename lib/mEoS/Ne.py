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


from scipy import log10
from scipy.constants import pi, Avogadro

from lib import unidades
from lib.meos import MEoS


class Ne(MEoS):
    """Multiparameter equation of state for neon"""
    name = "neon"
    CASNumber = "7440-01-9"
    formula = "Ne"
    synonym = "R-720"
    _refPropName = "NEON"
    _coolPropName = "Neon"
    rhoc = unidades.Density(481.914888)
    Tc = unidades.Temperature(44.4918)
    Pc = unidades.Pressure(2678.6, "kPa")
    M = 20.179  # g/mol
    Tt = unidades.Temperature(24.556)
    Tb = unidades.Temperature(27.104)
    f_acent = -0.0387
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 107

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [], "exp": []}

    katti = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for neon of Katti (1986)",
        "__doi__": {"autor": "Katti, R.S., Jacobsen, R.T, Stewart, R.B., "
                             "Jahangiri, M.",
                    "title": "Thermodynamic Properties of Neon for "
                             "Temperatures from the Triple Point to 700 K at "
                             "Pressures up to 700 MPa",
                    "ref": "Adv. Cryo. Eng. 31 (1986) 1189-1197",
                    "doi": "10.1007/978-1-4613-2213-9_132"},

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 6179, "so": 146.214},

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 700000.0, "rhomax": 90.56,

        "nr1": [0.3532653449e1, -0.4513954384e1, -0.1524027959, 0.2188568609e1,
                -7.44299997, 0.7755627402e1, -0.3122553128e1, 0.1014206899e1,
                -0.5289214086e-1, 0.1566849239, -0.222852705, -0.1410150942e-1,
                0.7036229719e-1, -0.5882048367e-1, 0.1571172741e-1,
                0.1292202769e-2, 0.7902035603e-3, -0.3794403616e-3],
        "d1": [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6],
        "t1": [0.5, 0.75, 3.5, 0.5, 0.75, 1, 1.5, 2.5, 0.25, 0.5, 2.5, 1, 3, 4,
               5, 1, 5, 6],

        "nr2": [0.04652799333,  0.04524001818, -0.2383421991, 0.629359013e-2,
                -0.1272313644e-2, -0.175235256e-6, 0.7188419232e-2,
                -0.5403006914e-1, 0.7578222187e-1, -0.3808588254e-1,
                0.6034022431e-2],
        "d2": [1, 2, 2, 2, 2, 2, 4, 8, 8, 8, 8],
        "t2": [4, 1, 5, 8, 12, 32, 10, 6, 7, 8, 9],
        "c2": [3, 2, 2, 4, 6, 6, 2, 2, 2, 2, 2],
        "gamma2": [1]*11,

        "nr3": [],
        "nr4": []}

    eq = katti,

    _surface = {"sigma": [0.012254, 0.02728, -0.025715],
                "exp": [1.4136, 1.4517, 1.6567]}
    _dielectric = {
        "eq": 1,
        "a": [0.9969, 0], "b": [-0.109, 0.0708], "c": [-2.88, -1.0],
        "Au": 0, "D": 2}

    _melting = {
        "eq": 2,
        "__doi__": {
            "autor": "Santamaría-Pérez, D., Mukherjee, G.D., Schwager, B., "
                     "Boehler, R.",
            "title": "High-pressure melting curve of helium and neon: "
                     "Deviations from corresponding states theory",
            "ref": "Physical Review B 81 (2010) 214101",
            "doi": "10.1103/PhysRevB.81.214101"},

        "Tmin": 24.4, "Tmax": 700.0,
        "Tref": 24.4, "Pref": 101325,

        "a2": [0.17e9], "exp2": [1/0.77]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.55805e1, 0.68795e-1, 0.54840e1, -0.83760e1, 0.34276e1],
        "t": [1, 1.5, 2.3, 2.8, 3.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.0601, 120.76, -385.53, 816.55, -899.07, 354.66],
        "t": [0.33, 1.4, 1.7, 2.2, 2.6, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.23338e1, -0.36834e1, -0.85368e2, 0.22769e3, -0.17290e3],
        "t": [0.444, 0.95, 3.5, 4.1, 4.5]}

    visco0 = {"__name__": "Rabinovich (1988)",
              "__doi__": {
                  "autor": "Rabinovich, V.A., Vasserman, A.A., Nedostup, V.I.,"
                           " Veksler, L.S.",
                  "title": "Thermophysical Properties of Neon, Argon, "
                           "Krypton, and Xenon",
                  "ref": "Hemisphere Publishing Corp., 1988.",
                  "doi": ""},

              "eq": 0,
              "method": "_visco0"}

    _viscosity = visco0,

    def _visco0(self, rho, T, fase=None):
        a = [17.67484, -2.78751, 311498.7, -48826500, 3938774000, -1.654629e11,
             2.86561e12]
        Tr = T/0.29944
        y = 0.68321*(a[0] + a[1]*log10(Tr) + a[2]/Tr**2 + a[3]/Tr**3 +
                     a[4]/Tr**4 + a[5]/Tr**5 + a[6]/Tr**6)
        nt = 266.93*(T*self.M)**0.5/y
        om = rho/1673.0
        c = [1.03010, -0.99175, 2.47127, -3.11864, 1.57066]
        b = [0.48148, -1.18732, 2.80277, -5.41058, 7.04779, -3.76608]
        sum1 = sum([ci*om**i for i, ci in enumerate(c)])
        sum2 = sum([bi*om**i for i, bi in enumerate(b)])
        sigma = 3.05e-10*(sum1-sum2*log10(T/122.1))
        br = 2.0/3.0*pi*Avogadro*sigma**3
        brho = rho/self.M*1000*br
        d = [1, 0.27676, 0.014355, 2.6480, -1.9643, 0.89161]
        nd = sum([di*brho**i for i, di in enumerate(d)])
        return unidades.Viscosity(nd*nt/100, "muPas")
