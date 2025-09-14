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


from math import exp

from lib import unidades
from lib.newComponent._base import GroupContribution


class Joback(GroupContribution):
    """
    Group contribution for definition of unknown component using the Joback
    procedure (1987). This method is fairly complete, can calculate critical
    properties, boiling and melting temperature, enthalpy and gibbs free energy
    of formation, vaporization and melting heat, ideal gas heat capacity
    dependence with temperature and viscosity.

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound, using general properties for calculation of other
    mandatory properties don't defined by the method.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float, optional
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    M, Tb and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    p-dichlorobenzene example in [2]_, Table V

    >>> cmp = Joback(group=[16, 13, 14], contribution=[2, 4, 2])
    >>> "%0.1f %0.0f %0.0f %0.1f" % (cmp.Tb, cmp.Tf, cmp.Tc, cmp.Pc.bar)
    '443.4 256 675 41.5'
    >>> "%0.0f" % (cmp.Vc.ccg*cmp.M)
    '362'
    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '26.41 78.56'
    >>> "%0.0f %0.0f" % (cmp._Cp0(298).JgK*cmp.M, cmp._Cp0(400).JgK*cmp.M)
    '112 139'
    >>> "%0.0f %0.0f" % (cmp._Cp0(800).JgK*cmp.M, cmp._Cp0(1000).JgK*cmp.M)
    '206 224'
    >>> "%0.2f %0.1f" % (cmp.Hv.kJg*cmp.M, cmp.Hm.kJg*cmp.M)
    '40.66 13.3'
    >>> "%0.2e %0.2e" % (cmp._Visco(333.8), cmp._Visco(374.4))
    '7.26e-04 4.92e-04'
    >>> "%0.2e %0.1e" % (cmp._Visco(403.1), cmp._Visco(423.3))
    '3.91e-04 3.4e-04'

    Example 2-1 in [1]_, 2-ethylphenol critical properties

    >>> cmp = Joback(group=[0, 1, 13, 14, 20], contribution=[1, 1, 4, 2, 1])
    >>> "%0.2f %0.1f %0.2f" % (cmp.Tb, cmp.Tc, cmp.Pc.bar)
    '489.94 716.0 44.09'
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '341.5'
    >>> cmp.formula
    'C8H10O'

    Example 3-1 in [1]_, 2,4 dimethylphenol ΔH and ΔG

    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '-149.23 -25.73'
    >>> "%0.1f" % (cmp._Cp0(700).JgK*cmp.M)
    '281.2'

    Example 2-10 in [1]_, 2,4 dimethylphenol Tb and Tf

    >>> cmp = Joback(group=[0, 13, 14, 20], contribution=[2, 3, 3, 1])
    >>> "%0.2f %0.2f" % (cmp.Tf, cmp.Tb)
    '330.58 494.92'

    Example in http://en.wikipedia.org/wiki/Joback_method, acetone

    >>> cmp = Joback(group=[0, 23], contribution=[2, 1])
    >>> "%0.3f %0.3f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Tb, cmp.Tf)
    '500.559 48.025 322.11 173.5'
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '209.5'
    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '-217.83 -154.54'
    >>> "%0.4f" % (cmp._Cp0(300).JgK*cmp.M)
    '75.3264'
    >>> "%0.2f %0.2f" % (cmp.Hm.kJg*cmp.M, cmp.Hv.kJg*cmp.M)
    '5.12 29.02'
    >>> "%0.7f" % (cmp._Visco(300))
    '0.0002942'

    Example in [3]_ pag. 2-470, o-xylene

    >>> cmp = Joback(group=[13, 14, 0], contribution=[4, 2, 2], Tb=417.58)
    >>> "%0.2f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '630.37 35.86 375.5'

    Example in [3]_ pag. 2-470, sec-butanol

    >>> cmp = Joback(group=[0, 1, 2, 19], contribution=[2, 1, 1, 1], Tb=372.7)
    >>> "%0.1f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '534.1 44.33 272.5'
    """
    __title__ = "Joback-Reid (1987)"
    __doi__ = {
      1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
      2:
        {"autor": "Joback, K.G., Reid, R.C.",
         "title": "Estimation of Pure-Component Properties from "
                  "Group-Contributions.",
         "ref": "Chemical Engineering Communications, 57 (1987) 233-243",
         "doi": "10.1080/00986448708960487"},
      3:
        {"autor": "",
         "title": "Perry's Chemical Engineers' Handbook 8th Edition",
         "ref": "McGraw Hill (2008)",
         "doi": ""}}

    __coeff__ = {
        # Table III
        "tc": [0.0141, 0.0189, 0.0164, 0.0067, 0.0113, 0.0129, 0.0117, 0.0026,
               0.0027, 0.002, 0.01, 0.0122, 0.0042, 0.0082, 0.0143, 0.0111,
               0.0105, 0.0133, 0.0068, 0.0741, 0.024, 0.0168, 0.0098, 0.038,
               0.0284, 0.0379, 0.0791, 0.0481, 0.0143, 0.0243, 0.0295, 0.0130,
               .0169, .0255, .0085, .0, .0496, .0437, 0.0031, 0.0119, 0.0019],
        "Pc": [-0.0012, 0.0, 0.002, 0.0043, -0.0028, -0.0006, 0.0011, 0.0028,
               -0.0008, 0.0016, 0.0025, 0.0004, 0.0061, 0.0011, 0.0008, -.0057,
               -0.0049, 0.0057, -0.0034, 0.0112, 0.0184, 0.0015, 0.0048, .0031,
               0.0028, 0.0030, 0.0077, 0.0005, 0.0101, 0.0109, 0.0077, 0.0114,
               0.0074, -0.0099, .0076, .0, -.0101, .0064, .0084, .0049, .0051],
        "vc": [65, 56, 41, 27, 56, 46, 38, 36, 46, 37, 48, 38, 27, 41, 32, 27,
               58, 71, 97, 28, -25, 18, 13, 62, 55, 82, 89, 82, 36, 38, 35, 29,
               9, 0, 34, 0, 91, 91, 63, 54, 38],
        "tb": [23.58, 22.88, 21.74, 18.25, 18.18, 24.96, 24.14, 26.15, 9.20,
               27.38, 27.15, 21.78, 21.32, 26.73, 31.01, -0.03, 38.13, 66.86,
               93.84, 92.88, 76.34, 22.42, 31.22, 76.75, 94.97, 72.24, 169.09,
               81.10, -10.5, 73.23, 50.17, 52.82, 11.74, 74.6, 57.55, 0.0,
               125.66, 152.54, 63.56, 68.78, 52.10],
        "tf": [-5.10, 11.27, 12.64, 46.43, -4.32, 8.73, 11.14, 17.78, -11.18,
               64.32, 7.75, 19.88, 60.15, 8.13, 37.02, -15.78, 13.55, 43.43,
               41.69, 44.45, 82.83, 22.23, 23.05, 61.20, 75.97, 36.9, 155.5,
               53.6, 2.08, 66.89, 52.66, 101.51, 48.84, 0, 68.4, 0.0, 59.89,
               127.24, 20.09, 34.4, 79.93],
        "hf": [-76.45, -20.64, 29.89, 82.23, -9.63, 37.97, 83.99, 142.14,
               79.30, 115.51, -26.8, 8.67, 79.72, 2.09, 46.43, -251.92, -71.55,
               -29.48, 21.06, -208.04, -221.65, -132.22, -138.16, -133.22,
               -164.50, -162.03, -426.72, -337.92, -247.61, -22.02, 53.47,
               31.65, 123.34, 23.61, 55.52, 93.7, 88.43, -66.57, -17.33, 41.87,
               39.1],
        "gf": [-43.96, 8.42, 58.36, 116.02, 3.77, 48.53, 92.36, 136.7, 77.71,
               109.82, -3.68, 40.99, 87.88, 11.30, 54.05, -247.19, -64.31,
               -38.06, 5.74, -189.2, -197.37, -105.0, -98.22, -120.50, -126.27,
               -143.48, -387.87, -301.95, -250.83, 14.07, 89.39, 75.61, 163.16,
               0.0, 79.93, 119.66, 89.22, -16.83, -22.99, 33.12, 27.73],
        "hv": [567, 532, 404, 152, 412, 527, 511, 636, 276, 789, 573, 464, 154,
               608, 731, -160, 1083, 1573, 2275, 4021, 2987, 576, 1119, 2144,
               1588, 2173, 4669, 2302, 1412, 2578, 1538, 1656, 453, 797, 1560,
               2908, 3071, 4000, 1645, 1629, 1430],
        "hm": [217, 619, 179, -349, -113, 643, 732, 1128, 555, 992, 117, 775,
               -328, 263, 572, 334, 601, 861, 651, 575, 1073, 284, 1405, 1001,
               0, 764, 2641, 1663, 866, 840, 1197, 1790, 1124, 0, 872, 0, 577,
               2313, 564, 987, 372],
        "cpa": [19.5, -0.909, -23.0, -66.2, -23.6, -8.0, -28.1, 27.4, 24.5,
                7.87, -6.03, 8.67, -90.9, -2.14, -8.25, 26.5, 33.3, 28.6, 32.1,
                25.7, -2.81, 25.5, 12.2, 6.45, 30.4, 30.9, 24.1, 24.5, 6.82,
                26.9, -1.21, 11.8, -31.1, 0.0, 8.83, 5.69, 36.5, 25.9, 35.3,
                19.6, 16.7],
        "cpb": [-8.08e-3, 9.5e-2, 2.04e-1, 4.27e-1, -3.81e-2, 1.05e-1, 2.08e-1,
                -5.57e-2, -2.71e-2, 2.01e-2, 8.54e-2, 1.62e-1, 5.57e-1,
                5.74e-2, 1.01e-1, -9.13e-2, -9.63e-2, -6.49e-2, -6.41e-2,
                -6.91e-2, 1.11e-1, -6.32e-2, -1.26e-2, 6.7e-2, -8.29e-2,
                -3.36e-2, 4.27e-2, 4.02e-2, 1.96e-2, -4.12e-2, 7.62e-2,
                -2.3e-2, 2.27e-1, 0.0, -3.84e-3, -4.12e-3, -7.33e-2, -3.74e-3,
                -7.58e-2, -5.61e-3, 4.81e-3],
        "cpc": [1.53e-4, -5.44e-5, -2.65e-4, -6.41e-4, 1.72e-4, -9.63e-5,
                -3.06e-4, 1.01e-4, 1.11e-4, -8.33e-6, -8.0e-6, -1.6e-4,
                -9.0e-4, -1.64e-6, -1.42e-4, 1.91e-4, 1.87e-4, 1.36e-4,
                1.26e-4, 1.77e-4, -1.16e-4, 1.11e-4, 6.03e-5, -3.57e-5,
                2.36e-4, 1.6e-4, 8.04e-5, 4.02e-5, 1.27e-5, 1.64e-4, -4.86e-5,
                1.07e-4, -3.2e-4, 0.0, 4.35e-5, 1.28e-4, 1.84e-4, 1.29e-4,
                1.85e-4, 4.02e-5, 2.77e-5],
        "cpd": [-9.67e-8, 1.19e-8, 1.2e-7, 3.01e-7, -1.03e-7, 3.56e-8, 1.46e-7,
                -5.02e-8, -6.78e-8, 1.39e-9, -1.8e-8, 6.24e-8, 4.69e-7,
                -1.59e-8, 6.78e-8, -1.03e-7, -9.96e-8, -7.45e-8, -6.87e-8,
                -9.88e-8, 4.94e-8, -5.48e-8, -3.86e-8, 2.86e-9, -1.31e-7,
                -9.88e-8, -6.87e-8, -4.52e-8, -1.78e-8, -9.76e-8, 1.05e-8,
                -6.28e-8, 1.46e-7, 0.0, -2.6e-8, -8.88e-8, -1.03e-7, -8.88e-8,
                -1.03e-7, -2.76e-8, -2.11e-8],
        "mua": [548.29, 94.16, -322.15, -573.56, 495.01, 82.28, 0, 0, 0, 0,
                307.53, -394.29, 0, 259.65, -245.74, 0, 625.45, 738.91, 809.55,
                2173.72, 3018.17, 122.09, 440.24, 340.35, 0, 740.92, 1317.23,
                483.88, 675.24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "mub": [-1.719, -0.199, 1.187, 2.307, -1.539, -0.242, 0, 0, 0, 0,
                -0.798, 1.251, 0, -0.702, 0.912, 0, -1.814, -2.038, -2.224,
                -5.057, -7.314, -0.386, -0.953, -0.35, 0, -1.713, -2.578,
                -0.966, -1.34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],

        # Name and element composition
        "txt": [("CH3", ),
                ("CH2", ),
                ("CH", ),
                ("C", ),
                ("=CH2", ),
                ("=CH", ),
                ("=C", ),
                ("=C=", ),
                ("≡CH", ),
                ("≡C", ),
                ("CH2 (cyclic)", ),
                ("CH (cyclic)", ),
                ("C (cyclic)", ),
                ("-CH (Aromatic)", ),
                ("=C (Aromatic)", ),
                ("F", ),
                ("Cl", ),
                ("Br", ),
                ("I", ),
                ("-OH", ),
                ("-OH (Aromatic)", ),
                ("-O-", ),
                ("-O- (cyclic)", ),
                ("C=O", ),
                ("C=O (cyclic)", ),
                ("CH=O", ),
                ("COOH", ),
                ("COO", ),
                ("=O", ),
                ("NH2", ),
                ("NH", ),
                ("NH (cyclic)", ),
                ("N", ),
                ("=N-", ),
                ("=N- (cyclic)", ),
                ("=NH", ),
                ("CN", ),
                ("NO2", ),
                ("SH", ),
                ("S", ),
                ("S (cyclic)", )]}

    FirstOrder = 41

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        if self.kwargs["Tb"]:
            Tb = self.kwargs["Tb"]
        else:
            # Eq 2
            Tb = 198.2+sum([c*self.__coeff__["tb"][i] for i, c in zip(
                self.kwargs["group"], self.kwargs["contribution"])])
        self.Tb = unidades.Temperature(Tb)

        # Equations of Table II
        self.Na = self._atoms()
        tcsuma, pcsuma, vcsuma = 0, 0, 0
        Tf = 122.5
        Hf = 68.29
        Gf = 53.88
        Hv = 15.3
        Hm = -0.88
        cpa = -37.93
        cpb = 0.21
        cpc = -3.91e-4
        cpd = 2.06e-7
        mua, mub = 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            Tf += c*self.__coeff__["tf"][i]                             # Eq 3
            tcsuma += c*self.__coeff__["tc"][i]                         # Eq 4
            pcsuma += c*self.__coeff__["Pc"][i]                         # Eq 5
            vcsuma += c*self.__coeff__["vc"][i]                         # Eq 6
            Hf += c*self.__coeff__["hf"][i]                             # Eq 7
            Gf += c*self.__coeff__["gf"][i]                             # Eq 8
            Hv += c*self.__coeff__["hv"][i]*0.004184                    # Eq 10
            Hm += c*self.__coeff__["hm"][i]*0.004184                    # Eq 11
            cpa += c*self.__coeff__["cpa"][i]
            cpb += c*self.__coeff__["cpb"][i]
            cpc += c*self.__coeff__["cpc"][i]
            cpd += c*self.__coeff__["cpd"][i]
            mua += c*self.__coeff__["mua"][i]
            mub += c*self.__coeff__["mub"][i]
        self.Tf = unidades.Temperature(Tf)
        self.Tc = unidades.Temperature(self.Tb/(0.584+0.965*tcsuma-tcsuma**2))
        self.Pc = unidades.Pressure((0.113+0.0032*self.Na-pcsuma)**-2, "bar")
        self.Vc = unidades.SpecificVolume((vcsuma+17.5)/1000/self.M)
        self.Hf = unidades.Enthalpy(Hf/self.M, "kJg")
        self.Gf = unidades.Enthalpy(Gf/self.M, "kJg")
        self.Hv = unidades.Enthalpy(Hv/self.M, "kJg")
        self.Hm = unidades.Enthalpy(Hm/self.M, "kJg")
        self.cp = [cpa, cpb, cpc, cpd, 0, 0]

        # Adjust the viscosity correlation with the parametric viscosity
        self.mul = [mua-597.82, mub-11.202]

        GroupContribution.calculo(self)

    def _Visco(self, T):
        """Viscosity calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        mu : float
            Viscosity, [Pas]
        """
        mu = self.M*exp(self.mul[0]/T + self.mul[1])
        return unidades.Viscosity(mu)

    def _Cp0(self, T):
        """Ideal gas specific heat calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        cp0 : float
            Ideal gas specific heat, [J/kgK]
        """
        cp0 = 0
        for i, c in enumerate(self.cp):
            cp0 += c*T**i

        return unidades.SpecificHeat(cp0/self.M, "JgK")
