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


class Cyclopentane(MEoS):
    """Multiparameter equation of state for cyclopentane"""
    name = "cyclopropane"
    CASNumber = "287-92-3"
    formula = "C5H10"
    synonym = ""
    rhoc = unidades.Density(274.920968)
    Tc = unidades.Temperature(511.72)
    Pc = unidades.Pressure(4571.2, "kPa")
    M = 70.1329  # g/mol
    Tt = unidades.Temperature(179.7)
    Tb = unidades.Temperature(322.405)
    f_acent = 0.201
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 36

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-0.3946233253, 2.4918910143],
           "ao_exp": [1.34, 13.4, 17.4, 6.65],
           "titao": [230/Tc, 1180/Tc, 2200/Tc, 5200/Tc],
           "ao_hyp": [], "hyp": []}

    CP1 = {"ao": 3.263,
           "an": [], "pow": [],
           "ao_exp": [2.151, 19.55, 14.45, 3.594],
           "exp": [179.0, 1336.0, 2911.0, 6420.0],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclopentane of Gedanitz et al. (2015).",
        "__doi__": {"autor": "Gedanitz, H., Davila, M.J., Lemmon, E.W.",
                    "title": "Speed of Sound Measurements and a Fundamental Equation of State for Cyclopentane",
                    "ref": "J. Chem. Eng. Data, 2015, 60 (5), pp 1331–1337",
                    "doi": "10.1021/je5010164"},
        "__test__": """
            >>> st=Cyclopentane(T=330, rhom=0.01)
            >>> print "%0.0f %0.7g %0.2f %0.7g %0.7g %0.7g" % (\
                st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            330 0.02720379 0.01 86.25843 94.88857 205.6768
            >>> st=Cyclopentane(T=330, rhom=11)
            >>> print "%0.0f %0.7g %0.1f %0.7g %0.7g %0.7g" % (\
                st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            330 75.55974 11.0 100.6003 130.5568 1471.842
            >>> st=Cyclopentane(T=530, rhom=1)
            >>> print "%0.0f %0.7g %0.1f %0.7g %0.7g %0.7g" % (\
                st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            530 3.240235 1.0 156.4924 187.7243 195.4293
            >>> st=Cyclopentane(T=512, rhom=4)
            >>> print "%0.0f %0.7g %0.1f %0.7g %0.7g %0.7g" % (\
                st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            512 4.601539 4.0 161.5786 22857.91 113.0171
            >>> st=Cyclopentane(T=520, rhom=6)
            >>> print "%0.0f %0.7g %0.1f %0.7g %0.7g %0.7g" % (\
                st.T, st.P.MPa, st.rhoM, st.cvM.JmolK, st.cpM.JmolK, st.w)
            520 6.522373 6.0 159.2304 276.7530 234.2660
            """, # Table 5, Pag F

        "R": 8.314472, "rhoc": 3.92,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 250000.0, "rhomax": 12.11,
        "Pmin": 0.008854, "rhomin": 12.1,

        "nr1": [0.0630928, 1.50365, -2.37099, -0.484886, 0.191843],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.29, 0.85, 1.185, 0.45],

        "nr2": [-0.835582, -0.435929, 0.545607, -0.209741, -0.0387635],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.28, 1.8, 1.5, 2.9, 0.93],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.677674, -0.137043, -0.0852862, -0.128085, -0.00389381],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.05, 4.0, 2.33, 1.5, 1.0],
        "alfa3": [0.86, 0.85, 0.86, 1.53, 5.13],
        "beta3": [0.63, 2.8, 0.5, 0.95, 0.23],
        "gamma3": [1.22, 0.32, 0.22, 1.94, 1.21],
        "epsilon3": [0.684, 0.7, 0.77, 0.625, 0.42]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclopentane of Lemmon et al. (2008).",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "",
                    "ref": "unpublished equation",
                    "doi": ""},
        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 12.2,
        "Pmin": 0.0089, "rhomin": 12.1,

        "nr1": [0.4909331e-1, 0.1244679e1, -0.1990222e1, -0.5245596, 0.1764215],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.23, 0.94, 1.08, 0.53],

        "nr2": [-0.1066798e1, -0.5028152, 0.8484762, -0.4547443, -0.2767817e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.67, 1.8, 1.3, 2.5, 1.0],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.9455318, -0.3014822, -0.1675668, -0.637707],
        "d3": [1, 1, 3, 3],
        "t3": [0.87, 1.4, 2.4, 1.3],
        "alfa3": [1.023, 1.383, 0.996, 7.038],
        "beta3": [1.7, 1.55, 1.07, 87.17],
        "gamma3": [1.1, 0.64, 0.5, 1.26],
        "epsilon3": [0.713, 0.917, 0.688, 0.748]}

    eq = helmholtz1, helmholtz2

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.1905, 1.8637, -1.6442, -2.72],
        "exp": [1.0, 1.5, 5.5, 2.9]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.0741, 81.968, 173.88, -68.519, -184.74],
        "exp": [0.1, 0.9, 1.25, 1.4, 1.05]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.0559, -6.4211, -46.926, 28.082, -70.838],
        "exp": [0.1, 0.65, 3.2, 3.55, 7.5]}
