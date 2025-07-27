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


from tools.qt import translate

from lib import unidades
from lib.newComponent._base import GroupContribution


class Lydersen(GroupContribution):
    """
    Group contribution for definition of unknown component using the Lydersen
    procedure (1955). This method is able to calculate the critical properties.

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
    Tb : float
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2 in [1]_ pag 2-343, 2-butanol critical properties

    >>> cmp = Lydersen(Tb=372.7, group=[0, 1, 2, 18],
    ... contribution=[2, 1, 1, 1])
    >>> "%0.1f %0.3f" % (cmp.Tc, cmp.Pc.MPa)
    '534.5 4.506'
    >>> cmp.formula
    'C4H10O'

    Example in http://en.wikipedia.org/wiki/Lydersen_method, acetone

    >>> cmp = Lydersen(Tb=329.25, group=[0, 22], contribution=[2, 1])
    >>> "%0.0f" % (cmp.Vc.ccg*cmp.M)
    '210'
    """
    __title__ = "Lydersen (1955)"

    __doi__ = {
      1:
        {"autor": "Maloney, J.O.",
         "title": "Perry's Chemical Engineers' Handbook 8th Edition",
         "ref": "McGraw Hill (2008)",
         "doi": ""},
      2:
        {"autor": "Lydersen, A. L.",
         "title": "Estimation of Critical Properties of Organic Compounds",
         "ref": "Coll. Eng. Univ. Wisconsin, Engineering Experimental Station "
                "Rept. 3, Madison, WI (1955)",
         "doi": ""}}

    __coeff__ = {
        # Table III
        "tc": [0.020, 0.020, 0.012, 0.00, 0.018, 0.018, 0.00, 0.00, 0.005,
               0.005, 0.013, 0.012, -0.007, 0.011, 0.011, 0.011, 0.066, 0.066,
               0.082, 0.031, 0.021, 0.014, 0.040, 0.033, 0.048, 0.085, 0.047,
               0.02, 0.018, 0.017, 0.010, 0.012, 0.031, 0.031, 0.024, 0.014,
               0.007, 0.060, 0.055, 0.015, 0.015, 0.008, 0.003, 0.026, 0.040,
               0.027, 0.025, 0.027],
        "Pc": [0.227, 0.227, 0.210, 0.210, 0.198, 0.198, 0.198, 0.198, 0.153,
               0.153, 0.184, 0.192, 0.154, 0.154, 0.154, 0.154, 0.924, 0.924,
               0.06, -0.02, 0.16, 0.12, 0.29, 0.2, 0.33, 0.4, 0.47, 0.12,
               0.224, 0.320, 0.50, 0.83, 0.095, 0.135, 0.09, 0.17, 0.13, 0.36,
               0.42, 0.27, 0.27, 0.24, 0.24, 0.468, 0.513, 0, 0.730, 0.668],
        "vc": [0.055, 0.055, 0.051, 0.041, 0.045, 0.045, 0.036, 0.036, 0.036,
               0.036, 0.0445, 0.046, 0.031, 0.036, 0.036, 0.037, 0, 0, 0.018,
               0.003, 0.020, 0.008, 0.060, 0.050, 0.073, 0.080, 0.080, 0.011,
               0.018, 0.049, 0.070, 0.095, 0.028, 0.037, 0.027, 0.042, 0.032,
               0.080, 0.078, 0.055, 0.055, 0.045, 0.047, 0, 0, 0, 0, 0],

        # Name and element composition
        "txt": [("CH3-", ),                     # 0
                ("-CH2-", ),
                ("-CH<", ),
                (">C<", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=C<", ),
                ("=C=", ),
                ("≡CH", ),
                ("≡C-", ),
                ("-CH2- (cyclic)", ),           # 10
                ("-CH< (cyclic)", ),
                (">C< (cyclic)", ),
                ("=C< (cyclic)", ),
                ("=C= (cyclic)", ),
                ("=CH- (cyclic)", ),
                ("-CH= (Aromatic)", ),
                ("=C< (Aromatic)", ),
                ("-OH", ),
                ("-OH (Aromatic)", ),
                ("-O-", ),                      # 20
                ("-O- (cyclic)", ),
                (">C=O", ),
                (">C=O (cyclic)", ),
                ("-CH=O", ),
                ("-COOH", ),
                ("-COO-", ),
                ("=O (other)", ),
                ("F", ),
                ("Cl", ),
                ("Br", ),                       # 30
                ("I", ),
                ("-NH2", ),
                (">NH", ),
                (">NH (cyclic)", ),
                (">N-", ),
                (">N- (cyclic)", ),
                ("-CN", ),
                ("NO2", ),
                ("-SH", ),
                ("-S-", ),                      # 40
                ("-S- (cyclic)", ),
                ("=S", ),
                (">Si<", ),
                ("-SiH<", ),
                ("-SiH3", ),
                (">SiO-", ),
                (">SiO- (cyclic)", ),

                ]}

    FirstOrder = 47

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = translate("newComponent", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self.__coeff__["tc"][i]
            pc += c*self.__coeff__["Pc"][i]
            vc += c*self.__coeff__["vc"][i]

        self.Tc = unidades.Temperature(self.Tb/(0.567+tc-tc**2))
        self.Pc = unidades.Pressure(self.M/(0.34+pc)**2, "atm")
        self.Vc = unidades.SpecificVolume((0.04+vc)/self.M)

        GroupContribution.calculo(self)
