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


class Li(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Li-Xia-Xiang procedure (2016)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
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
    Example A.1 in [21]_, critical temperature of hexanal

    >>> cmp = Li(Tb=401.45, group=[0, 1, 2, 10, 14, 17],
    ... contribution=[6, 12, 1, 5, 1, 12])
    >>> "%0.2f" % cmp.Tc
    '588.90'
    >>> cmp.formula
    'C6H12O'

    Example A.2 in [21]_, critical temperature of chlorotrimethylsilane

    >>> cmp = Li(Tb=330.75, group=[0, 1, 6, 9, 24, 25, 17],
    ... contribution=[3, 9, 1, 1, 3, 1, 9])
    >>> "%0.2f" % cmp.Tc
    '497.69'

    Example A.3 in [21]_, critical pressure of n-propyl formate

    >>> cmp = Li(Tb=353.97, group=[0, 1, 2, 10, 14, 17, 13],
    ... contribution=[4, 8, 2, 2, 1, 8, 2])
    >>> "%0.3f" % cmp.Pc.MPa
    '4.033'

    Example A.4 in [21]_, critical volume of acetic acid
    The paper has a bug, the calculated value is erroneous

    >>> cmp = Li(Tb=351.44, group=[0, 1, 2, 10, 14, 17, 16],
    ... contribution=[2, 4, 2, 1, 1, 3, 1])
    >>> "%0.2f" % (cmp.Vc.ccg*cmp.M)
    '180.65'
    """
    __title__ = "Li-Xia-Xiang (2016)"

    __doi__ = (
        {"autor": "Li, J., Xia, L., Xiang, S.",
         "title": "A New Method Based on Elements and Chemical Bonds for "
                  "Organic Compounds Critical Properties Estimation",
         "ref": "Fluid Phase Equil. 417 (2016) 1-6",
         "doi": "10.1016/j.fluid.2016.01.008"})

    _coeff = {
        "tc": [-0.0003, -0.0016, -0.0472, -0.0533, 0.0179, -0.0478, -0.0380,
               0.0071, 0.0091, 0.0224, -0.0168, -0.0080, -0.0079, 0.0059,
               0.0200, 0, -0.0114, -0.0001, 0.0285, 0.0358, 0.0071, 0.0091,
               -0.0153, 0.0360, -0.0316, 0.0209, 0.0268, 0.0029, -0.0192,
               0.0152, -0.0024, 0, 0.0180, 0.0171, 0, 0.0528, -0.0061, -0.0055,
               0.0238, 0.0022, 0.0209, 0.0285, -0.0398],
        "Pc": [0.3304, 0.0118, 1.0330, 0.4774, 0.5733, 0.9346, 1.1961, 0.9630,
               1.4724, 0.0654, 0.3317, 0.2350, 0.0604, -0.1030, -0.6309,
               0.0000, -0.6777, 0.0179, -0.2329, -0.1194, 0.9630, 1.4724,
               0.2999, 0.2856, 0.6898, 0.1418, 0.1493, 0.0000, 0.0000, 0.0973,
               0.2519, 0.0000, -0.0060, -0.1233, 0.0000, -0.6638, 0.1793,
               0.0368, -0.2843, 0.1783, -0.0381, -0.1360, 0.4207],
        "vc": [60.0289, 5.6116, 60.5913, -58.5859, 58.2972, 51.3276, 90.5686,
               88.0064, 128.9134, 30.6559, 5.3217, 3.3178, -5.4347, -5.3629,
               -15.8297, 0.0000, -12.3706, -1.3769, -9.4243, 8.1548, 88.0064,
               128.9134, 29.1266, 29.0457, 27.0244, 46.3533, 49.2458, 0, 0,
               19.6420, 48.3117, 0, 116.9220, 39.4229, 0, 23.5185, -2.1851,
               -15.5277, -14.5939, 9.7344, 28.5488, 91.7583, -50.7082],

        "txt": [("C",),                          # 0
                ("H",),
                ("O",),
                ("N",),
                ("S",),
                ("F",),
                ("Cl",),
                ("Br",),
                ("I",),
                ("Si",),

                # 2nd Order term
                ("C-C",),                        # 10
                ("C=C",),
                ("C≡C",),
                ("C-O",),
                ("C=O",),
                ("O-O",),
                ("O-H",),
                ("C-H",),
                ("C-F",),
                ("C-Cl",),
                ("C-Br",),                       # 20
                ("C-I",),
                ("C-S",),
                ("C=S",),
                ("C-Si",),
                ("Si-Cl",),
                ("Si-O",),
                ("Si-H",),
                ("S-S",),
                ("S-H",),
                ("C-N",),                        # 30
                ("C=N",),
                ("C≡N",),
                ("N-H",),
                ("N-N",),
                ("N=O",),
                ("C-C [r]",),
                ("C=C [r]",),
                ("C-O [r]",),
                ("C-S [r]",),
                ("C-N [r]",),                    # 40
                ("C=N [r]",),
                ("benzene",)]}

    FirstOrder = 10
    SecondOrder = 43

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

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]

        self.Tc = unidades.Temperature(self.Tb*(1.5530+tc)+18.9999)
        self.Pc = unidades.Pressure(self.M/(1.2220+pc)**2, "MPa")
        self.Vc = unidades.SpecificVolume(
                (19.6531-1.2603*self.M+vc)/self.M, "ccg")

        GroupContribution.calculo(self)
