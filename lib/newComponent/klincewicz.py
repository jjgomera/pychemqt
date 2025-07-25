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


class Klincewicz(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Klincewicz-Reid procedure (1984). This method is able to calculate the
    critical properties.

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound, using general properties for calculation of other
    mandatory properties don't defined by the method.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    nogroup : boolean
        Use the simple correlation without group contribution
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]
    atoms : int, optional
        Atoms count, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation. This
    method has a alternate not group contribution procedure known only Tb and M
    Na is only necessary if the nogroup opcion is enabled

    Examples
    --------
    Example in https://en.wikipedia.org/wiki/Klincewicz_method, acetone

    >>> cmp = Klincewicz(Tb=329.25, group=[0, 18], contribution=[2, 1])
    >>> "%0.2f %0.2f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '510.48 45.69 213.52'

    Same example without group contribution

    >>> cmp = Klincewicz(Tb=329.25, M=58.08, atoms=10, nogroup=True)
    >>> "%0.2f %0.4f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '505.15 52.9098 205.2'
    """
    __title__ = "Klincewicz (1984)"

    __doi__ = (
        {"autor": "Klincewicz, K.M., Reid, R.C.",
         "title": "Estimation of Critical Properties with Group Contribution "
                  "Methods",
         "ref": "AIChE J. 30(1) (1984) 137-142",
         "doi": "10.1002/aic.690300119"}, )

    kwargs = GroupContribution.kwargs.copy()
    kwargs["nogroup"] = False
    kwargs["atoms"] = 0

    __coeff__ = {
        # Table 6
        "tc": [-2.433, 0.353, 4.253, 6.266, -0.335, 16.416, 12.435, -0.991,
               3.786, 3.373, 7.169, 7.169, 5.623, -4.561, 7.341, -28.930,
               5.389, 7.127, 4.332, 4.332, -25.085, 8.890, -4.153, 2.005,
               2.773, 12.253, 8.239, -10.381, 28.529, 23.905, 31.537, 5.191,
               18.353, 53.456, 94.186, -1.770, -1.770, -1.770, -1.770, 11.709],
        "Pc": [0.026, -0.015, -0.046, -0.083, -0.027, -0.136, -0.111, -0.015,
               -0.050, -0.066, -0.067, -0.067, -0.089, -0.056, -0.112, -0.190,
               -0.143, -0.116, -0.196, -0.196, -0.251, -0.277, -0.127, -0.180,
               -0.172, -0.163, -0.104, -0.064, -0.303, -0.311, -0.208, -0.067,
               -0.244, -0.692, -1.051, 0.032, 0.032, 0.032, 0.032, -0.325],
        "vc": [16.2, 16.1, 8.2, 12.1, 7.4, 8.95, -6.6, 13.9, 9.8, 5.1, 2.7,
               2.7, 0.2, 7.5, 3.0, -24.0, -26.1, -36.6, -6.7, -6.7, -37.0,
               -28.2, -0.1, 53.7, -8.0, -0.7, -18.4, 12.0, -27.7, -27.3, -61.9,
               -34.1, -47.4, -148.1, -270.6, 0.8, 0.8, 0.8, 0.8, -39.2],

        "txt": [("-CH3", ),                     # 0
                ("-CH2-", ),
                ("-CH2- (ring)", ),
                (">CH-", ),
                (">CH- (ring)", ),
                (">C<", ),
                (">C< (ring)", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=CH- (ring)", ),
                (">c=", ),                      # 10
                ("=c=", ),
                (">C= (ring)", ),
                ("=CH", ),
                ("s-", ),
                ("-OH", ),
                ("-0-", ),
                ("-0- (ring)", ),
                (">CO", ),
                ("-CHO", ),
                ("-COOH", ),                    # 20
                ("-CO-0-", ),
                ("-NH2", ),
                (">NH", ),
                (">NH (ring)", ),
                (">N-", ),
                ("=N- (ring)", ),
                ("-CN", ),
                ("-SH", ),
                ("-S-", ),
                ("-S- (ring)", ),               # 30
                ("-F", ),
                ("-C1", ),
                ("-Br", ),
                ("-I", ),
                ("-CF2", ),
                ("-CCl2", ),
                ("-CBr2", ),
                ("-CI2", ),
                ("-NO2", )]}

    FirstOrder = 40

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = translate("newComponent", "undefined boiling point")
            self.status = 0
        elif self.kwargs["nogroup"]:
            self.group = []
            if not self.kwargs["M"]:
                self.msg = translate(
                    "newComponent", "undefined molecular weight")
                self.status = 0
            elif not self.kwargs["atoms"]:
                self.msg = translate(
                    "newComponent", "undefined atoms number of molecule")
                self.status = 0
            else:
                self.status = 1
                self.msg = ""
                return True
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        if self.kwargs["nogroup"]:
            # No group calculation
            self.Na = self.kwargs["atoms"]

            # Eq 13-15
            self.Tc = unidades.Temperature(50.2-0.16*self.M+1.41*self.Tb)
            self.Pc = unidades.Pressure(
                    self.M/(0.335+0.009*self.M+0.019*self.Na)**2, "bar")
            self.Vc = unidades.SpecificVolume(
                    (20.1+0.88*self.M+13.4*self.Na)/self.M, "ccg")

        else:
            Pc = tc = vc = 0
            for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
                tc += c*self.__coeff__["tc"][i]
                Pc += c*self.__coeff__["Pc"][i]
                vc += c*self.__coeff__["vc"][i]

            # Eq 10-12
            self.Tc = unidades.Temperature(45.4-0.77*self.M+1.55*self.Tb+tc)
            self.Pc = unidades.Pressure(
                    self.M/(0.348+0.0159*self.M+Pc)**2, "bar")
            self.Vc = unidades.SpecificVolume(
                    (25.2+2.8*self.M+vc)/self.M, "ccg")

        GroupContribution.calculo(self)
