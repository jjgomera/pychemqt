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


class oH2(MEoS):
    """Multiparameter equation of state for hydrogen (orto)"""
    name = "ortohydrogen"
    CASNumber = "1333-74-0o"
    formula = "H2"
    synonym = "R-702o"
    rhoc = unidades.Density(31.1361933)
    Tc = unidades.Temperature(33.22)
    Pc = unidades.Pressure(1310.65, "kPa")
    M = 2.01594  # g/mol
    Tt = unidades.Temperature(14.008)
    Tb = unidades.Temperature(20.4)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4675442336, 1.8845068862],
           "ao_exp": [2.54151, -2.3661, 1.00365, 1.22447],
           "titao": [856/Tc, 1444/Tc, 2194/Tc, 6968/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ortohydrogen of Leachman et al. (2007)",
        "__doi__": {"autor": "Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.",
                    "title": "Fundamental equations of state for parahydrogen, normal hydrogen, and orthohydrogen",
                    "ref": "J. Phys. Chem. Ref. Data, 38 (2009), 721 – 748",
                    "doi": "10.1063/1.3160306"},
        "__test__": """
            >>> st=oH2(T=14.008, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            14.008 7.5600 77.010 0.13273 −53.820 400.77 −3.0625 29.390 5.1746 6.2707 7.1448 10.557 1264.7 307.38
            """, # Table 15, Pag 746

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 38.2,
        "Pmin": 7.461, "rhomin": 38.2,

        "nr1": [-6.83148, 0.01, 2.11505, 4.38353, 0.211292, -1.00939, 0.142086],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.7333, 1, 1.1372, 0.5136, 0.5638, 1.6248, 1.829],

        "nr2": [-0.87696, 0.804927],
        "d2": [1, 3],
        "t2": [2.404, 2.105],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.710775, 0.0639688, 0.0710858, -0.087654, 0.647088],
        "d3": [2, 1, 3, 1, 1],
        "t3": [4.1, 7.658, 1.259, 7.589, 3.946],
        "alfa3": [1.169, 0.894, 0.04, 2.072, 1.306],
        "beta3": [0.4555, 0.4046, 0.0869, 0.4415, 0.5743],
        "gamma3": [1.5444, 0.6627, 0.763, 0.6587, 1.4327],
        "epsilon3": [0.6366, 0.3876, 0.9437, 0.3976, 0.9626],
        "nr4": []}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.488684e1, 0.105310e1, 0.856947, -0.185355],
        "exp": [1.0, 1.5, 2.7, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.43911e1, -0.75872e1, 0.10402e2, -0.72651e1, 0.18302e1],
        "exp": [0.53, 0.93, 1.35, 1.8, 2.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31463e1, -0.16183e2, 0.31803e2, -0.21961e3, 0.43123e3, -0.25591e3],
        "exp": [0.491, 2.1, 2.9, 4.4, 5.0, 5.5]}
