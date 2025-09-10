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

from tools.qt import translate

from lib import unidades
from lib.newComponent._base import GroupContribution


class Wilson(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Wilson-Jasperson procedure (1994). This method is able to calculate only
    critical temperature and pressure.

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound, using general properties for calculation of other
    mandatory properties don't defined by the method.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    ring : int
        Ring in the atom, [-]
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2-4 in [1]_, 2-ethylphenol critical properties

    >>> c1 = Wilson(Tb=477.67, ring=1)
    >>> c2 = Wilson(Tb=477.67, ring=1)
    >>> c1(group=[3, 0, 5], contribution=[8, 10, 1])
    >>> c2(group=[3, 0, 5, 42], contribution=[8, 10, 1, 1])
    >>> "%0.1f %0.2f %0.1f" % (c1.Tc, c1.Pc.bar, c2.Tc)
    '702.9 37.94 693.6'
    >>> c1.formula
    'C8H10O'
    """
    __title__ = "Wilson-Jasperson (1996)"
    __doi__ = {
      1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
      2:
        {"autor": "Wilson, G.M. Jasperson, L.V.",
         "title": "Critical constants Tc and Pc, estimation based on zero, "
                  "first and second order methods",
         "ref": "Paper given at AIChE Spring National Meeting, New Orleans, "
                "LA, USA, February 25-29, 1996.",
         "doi": ""}}

    kwargs = GroupContribution.kwargs.copy()
    kwargs["ring"] = 0

    __coeff__ = {
        "tc": [0.002793, 0.320000, 0.019000, 0.008532, 0.019181, 0.020341,
               0.008810, 0.036400, 0.088000, 0.020000, 0.012000, 0.007271,
               0.011151, 0.016800, 0.014000, 0.018600, 0.059000, 0.031000,
               0.007000, 0.010300, 0.012447, 0.013300, -0.027000, 0.175000,
               0.017600, 0.007000, 0.020000, 0.010000, 0.000000, 0.005900,
               0.017000, -0.027500, 0.219000, 0.013000, 0.011000, 0.014000,
               -0.050000, 0.000000, 0.000000, 0.007000, 0.015000, 0.0350,
               0.0100, -0.0075, -0.0040, 0.0000, -0.0550, 0.0170, -0.0150,
               0.0170, -0.0200, 0.0020, 0.0000, -0.0250],
        "Pc": [0.12660, 0.43400, 0.91000, 0.72983, 0.44805, 0.43360, 0.32868,
               0.12600, 6.05000, 1.34000, 1.22000, 1.04713, 0.97711, 0.79600,
               1.19000, 0.0, 0.0, 1.42000, 2.68000, 1.20000, 0.97151, 1.11000,
               0.0, 1.11000, 2.71000, 1.69000, 1.95000, 0.0, 0.43000, 1.315930,
               1.66000, 6.33000, 1.07000, 0.0, 1.08000, 0.0, 0.0, -0.08000,
               0.69000, 2.05000, 2.04000, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00,
               0.50, 0.00, 1.50, 1.00, 0.00, 0.00, -0.50],

        "txt": [("H",),
                ("He",),
                ("B",),
                ("C",),
                ("N",),
                ("O",),
                ("F",),
                ("Ne",),
                ("Al",),
                ("Si",),
                ("P",),
                ("S",),
                ("Cl",),
                ("Ar",),
                ("Ti",),
                ("V",),
                ("Ga",),
                ("Ge",),
                ("As",),
                ("Se",),
                ("Br",),
                ("Kr",),
                ("Rb",),
                ("Zr",),
                ("Nb",),
                ("Mo",),
                ("Sn",),
                ("Sb",),
                ("Te",),
                ("I",),
                ("Xe",),
                ("Cs",),
                ("Hf",),
                ("Ta",),
                ("W",),
                ("Re",),
                ("Os",),
                ("Hg",),
                ("Bi",),
                ("Rn",),
                ("U",),

                # 2nd Order term
                ("-OH, C4 or less",),
                ("-OH, C5 or more",),
                ("-O-",),
                ("-NH2, >NH, >N-",),
                ("-CHO",),
                (">CO",),
                ("-COOH",),
                ("-COO-",),
                ("-CN",),
                ("-NO2",),
                ("Organic Halides (once / molecule)",),
                ("-SH, -S-, -SS-",),
                ("Siloxane bond",)]}

    FirstOrder = 41
    SecondOrder = 54

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = translate("newComponent", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        tc = Pc = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self.__coeff__["tc"][i]
            Pc += c*self.__coeff__["Pc"][i]

        Nr = self.kwargs["ring"]
        self.Tc = unidades.Temperature(self.Tb/(0.048271-0.019846*Nr+tc)**0.2)
        self.Pc = unidades.Pressure(0.0186233*self.Tc/(-0.96601+exp(
            -0.00922295-0.0290403*Nr+0.041*Pc)), "bar")

        GroupContribution.calculo(self)
