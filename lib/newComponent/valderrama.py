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
from lib.newComponent._base import GroupContribution


class Valderrama(GroupContribution):
    """
    Group contribution for definition of unknown component using the Valderrama
    procedure (2006). This method is able to calculate the critical
    properties.

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
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example from in [1]_ Table 8, phenantrene

    >>> cmp = Valderrama(Tb=372.7, group=[30, 32], contribution=[10, 4])
    >>> "%0.0f %0.1f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '871 31.6 550.55'

    2-methyl-1-pentanol

    >>> cmp = Valderrama(Tb=401.85,
    ... group=[0, 1, 2, 10], contribution=[2, 3, 1, 1])
    >>> "%0.0f %0.1f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '602 33.5 381.53'
    """
    __title__ = "Valderrama (2006)"

    __doi__ = (
        {"autor": "Valderrama, J.O., Álvarez, V.H.",
         "title": "A New Group Contribution Method Based on Equation of State "
                  "Parameters to Evaluate the Critical Properties of Simple "
                  "and Complex Molecules",
         "ref": "Can. J. Chem. Eng. 84(4) (2006) 431-446",
         "doi": "10.1002/cjce.5450840404"}, )

    __coeff__ = {
        # Table III
        "tc": [8.26, 20.07, 27.11, 35.81, 3.27, 18.81, 26.64, 15.40, 6.46,
               6.46, 11.21, 9.93, 34.18, 21.98, 50.52, 40.07, 28.83, -7.64,
               13.70, 15.00, 24.66, 26.42, 49.30, 48.83, -7.25, 8.28, 16.07,
               41.78, 15.09, 21.92, 11.23, 31.39, 27.40, 2.13, 1.17, 27.48,
               4.28, 33.80, 9.23],
        "Pc": [0.8446, 1.2910, 1.4806, 1.6350, 0.4418, 1.1606, 1.4472, 0.6119,
               0.2734, 0.2734, 0.0881, 0.3628, 1.5195, 0.7953, 1.9871, 1.9804,
               1.3779, -0.9769, 0.4922, 0.2946, 0.9936, 1.8782, 2.5229, 2.0970,
               0.0210, 0.5375, 0.5006, 1.7649, 0.9240, 1.0535, 0.6178, 1.3134,
               1.2796, -0.0201, -0.4922, 1.0798, -0.4429, 1.2695, -0.0791],
        "vc": [45.42, 37.60, 30.00, 13.28, 41.36, 33.48, 21.93, 21.56, 28.60,
               28.60, 20.75, 10.81, 47.41, 52.80, 60.10, 56.69, 65.87, 26.61,
               32.98, 51.67, 17.34, 32.62, 58.84, 83.38, 21.46, 41.72, 51.66,
               68.45, 34.69, 20.42, 28.49, 11.00, 20.63, 10.59, -11.19, 40.67,
               19.04, 16.14, 29.24],

        # Name and element composition
        "txt": [("CH3-", ),                     # 0
                ("-CH2-", ),
                (">CH-", ),
                (">C<", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=C<", ),
                ("=C=", ),
                ("=(-)CH", ),
                ("=(-)C-", ),
                ("-OH", ),                      # 10
                ("-O-", ),
                (">C=O", ),
                ("-CHO", ),
                ("-COOH", ),
                ("-COO-", ),
                ("HCOO-", ),
                ("=O (Any_other)", ),
                ("-NH2", ),
                ("-NH-", ),
                (">N-", ),                      # 20
                ("=N-", ),
                ("-CN", ),
                ("-NO2", ),
                ("-F", ),
                ("-Cl", ),
                ("-Br", ),
                ("-I", ),
                ("-CH2- (ring)", ),
                (">CH- (ring)", ),
                ("=CH- (ring)", ),              # 30
                (">C< (ring)", ),
                ("=C< (ring)", ),
                ("-O- (ring)", ),
                ("-OH(Phenols) (ring)", ),
                (">C=O (ring)", ),
                ("-NH- (ring)", ),
                (">N- (ring)", ),
                ("=N- (ring)", )
                ]}

    FirstOrder = 38

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
            self.Tb = unidades.Temperature(self.kwargs["Tb"])

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self.__coeff__["tc"][i]
            pc += c*self.__coeff__["Pc"][i]
            vc += c*self.__coeff__["vc"][i]

        g1 = 38.91+tc**0.88
        g2 = 5.84+pc**1.27
        self.Pc = unidades.Pressure((g1/g2)**2, "bar")
        self.Tc = unidades.Temperature(self.Pc.bar*g2)
        self.Vc = unidades.SpecificVolume((26.86+vc**1.06)/self.M, "ccg")

        GroupContribution.calculo(self)
