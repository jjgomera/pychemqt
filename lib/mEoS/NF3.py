#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


class NF3(MEoS):
    """Multiparameter equation of state for nitrogen trifluoride"""
    name = "nitrogen trifluoride"
    CASNumber = "7783-54-2"
    formula = "NF3"
    synonym = ""
    rhoc = unidades.Density(562.47)
    Tc = unidades.Temperature(234.0)
    Pc = unidades.Pressure(4460.7, "kPa")
    M = 71.019  # g/mol
    Tt = unidades.Temperature(66.36)
    Tb = unidades.Temperature(144.138)
    f_acent = 0.126
    momentoDipolar = unidades.DipoleMoment(0.235, "Debye")
    # id = 951
    id = 60

    CP1 = {"ao": -7.140693612211,
           "an": [0.7427518245951e6, -0.4389825372134e5, 0.1012629224351e4,
                  0.5481339146452e-1, -0.7677196006769e-4, 0.4203630864340e-7],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-0.6328752997967],
           "exp": [3000],
           "ao_hyp": [], "hyp": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for nitrogen trifluoride of Younglove (1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.",
                    "doi": ""},

        "R": 8.314471,
        "cp": CP1,

        "Tmin": 85.0, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 25.3,
        "Pmin": 0.000186, "rhomin": 26.32,

        "b": [None, 0.1774353868e-1, -0.5409379418, 0.3976634466e1,
              -0.5209476694e3, -0.3286322888e5, -0.5990517411e-3, 0.9217525601,
              -0.4848977075e3, -0.4235892691e7, -0.9824248063e-5, 0.5432235989e-1,
              -0.1462388500e2, -0.3366180440e-2, 0.2801374599, 0.8435288597e1,
              -0.1324421452e-1, 0.1875604377e-3, 0.2959643991, -0.700997687e-2,
              0.4365820912e7, -0.1111397536e8, 0.2411866612e5, 0.3179136276e7,
              0.6166849090e2, 0.4260854720e2, 0.1090598789, -0.3340951059e2,
              0.8597429644e-4, 0.1240544214e-2, 0.1286224248e-6,
              -0.8941104276e-6, 0.3353054595e-4]}

    eq = MBWR,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.66672e1, 0.33864e1, -0.28222e1, -0.50602e1, 0.32481e1],
        "exp": [1.0, 1.5, 1.7, 5.5, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.22080e1, 0.35709e2, -0.92868e2, 0.66666e2, -0.93589e1],
        "exp": [0.35, 2.4, 2.7, 3.0, 4.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.3061e1, -0.80541e1, -0.19619e2, -0.13432e2, -0.3276e2, -0.67907e2],
        "exp": [0.421, 1.48, 3.9, 7.0, 8.0, 15.0]}
