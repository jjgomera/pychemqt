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


class Ambrose(GroupContribution):
    """
    Group contribution for definition of unknown component using the Ambrose
    procedure as use in API Technical Databook, procedure 4A1.1 with aditional
    term from Perry's Handbook. This method is able to calculate the critical
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
    platt : int
        ΔPlatt number, [-]
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation
    The Platt number is the number of pairs of carbon atoms which are
    separated by three carbon-carbon bonds and is an indicator of the degree
    of branching in the molecule. The Platt number of an n-alkane is equal to
    the number of carbons minus three.

    Examples
    --------
    Example 1 in [2]_, 2,2,3-Trimethylpentane

    >>> Tb = unidades.Temperature(229.72, "F")
    >>> cmp = Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1],
    ... Tb=Tb, platt=3)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '555.83 400.74 0.0639'
    >>> cmp.formula
    'C8H18'

    Example 2 in [2]_, 2-Methyl-1-butene

    >>> Tb = unidades.Temperature(88.09, "F")
    >>> cmp = Ambrose(group=[0, 1, 4, 6], contribution=[2, 1, 1, 1],
    ... Tb=Tb, platt=0)
    >>> "%0.2f %0.1f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '385.95 520.3 0.0657'

    Example 3 in [2]_, cis-Decalin

    >>> Tb = unidades.Temperature(384.47, "F")
    >>> cmp = Ambrose(group=[10, 12], contribution=[8, 2], Tb=Tb, platt=0)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '801.95 430.06 0.0562'

    Example 4 in [2]_, tert-Butyl benzene

    >>> Tb = unidades.Temperature(336.41, "F")
    >>> cmp = Ambrose(group=[0, 3, 17], contribution=[3, 1, 1],
    ... Tb=Tb, platt=-1)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '705.82 415.97 0.0555'

    Example 5 in [2]_, Anthracene

    >>> Tb = unidades.Temperature(646.16, "F")
    >>> cmp = Ambrose(group=[18, 28], contribution=[1, 2], M=178.23, Tb=Tb,
    ... platt=0)
    >>> "%0.1f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '1165.3 504.64 0.0502'

    Example from [3]_, 2,2,4-trimethylpentane

    >>> cmp = Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1],
    ... Tb=372.39, platt=0)
    >>> "%0.1f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '543.0 25.63 455.8'
    """
    __title__ = "Ambrose (1980)"
    __doi__ = (
        {"autor": "Ambrose, D.",
         "title": "Correlation and Estimation of Vapor-Liquid Critical "
                  "Properties: I. Critical Temperatures of Organic Compounds",
         "ref": "National Physical Laboratory, Teddington, NPL Rep. Chern.  "
                "92, 1978, corrected 1980.",
         "doi": ""},
        {"autor": "Ambrose, D.",
         "title": "Correlation and Estimation of Vapor-Liquid Critical "
                  "Properties: II. Critical Pressures and Volumes of Organic "
                  "Compounds",
         "ref": "National Physical Laboratory, Teddington, NPL Rep. 98, 1979",
         "doi": ""},
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
        {"autor": "Maloney, J.O.",
         "title": "Perry's Chemical Engineers' Handbook 8th Edition",
         "ref": "McGraw Hill (2008)",
         "doi": ""})

    kwargs = GroupContribution.kwargs.copy()
    kwargs["platt"] = 0

    __coeff__ = {
        "Pc": [0.2260, 0.2260, 0.22, 0.1960, 0.1935, 0.1935, 0.1875, 0.1610,
               0.1410, 0.1410, 0.1820, 0.1820, 0.1820, 0.1820, 0.1495, 0.1495,
               0.1170, 0.9240, 0.8940, 0.9440, 0.9440, 0.8640, 0.9140, 0.8340,
               0.8840, 0.8840, 0.8040, 0.7240, 0.5150, None, 0.160, 0.282,
               0.220, 0.450, 0.900, 0.470, 0.420, 0.095, 0.135, 0.170, 0.360,
               0.270, 0.270, 0.461, 0.507, 0, 0.725, 0.663, 0.223, 0.318,
               0.500, -0.025, 0.515, 0.183, 0.318, 0.600, 0.850, -0.065,
               -0.170, 0.924, 0.850, 0, 0.020, -0.050, 0],
        "tc": [0.138, 0.138, 0.095, 0.018, 0.113, 0.113, 0.070, 0.088, 0.038,
               0.038, 0.09, 0.09, 0.03, 0.09, 0.075, 0.075, 0.06, 0.458, 0.448,
               0.488, 0.488, 0.438, 0.478, 0.428, .468, .468, .418, .368, .22,
               None, 0.138, 0.220, 0.220, 0.578, 1.156, 0.330, 0.370, 0.208,
               0.208, 0.088, 0.423, 0.105, 0.090, 0.138, 0.371, 0.195, 0.159,
               0.131, 0.180, 0.110, 0.110, 0.198, 0.220, 0.080, 0.080, 0.080,
               0.010, -0.050, -0.200, 0.448, 0.448, 0.030, -0.040, -0.080],
        "vc": [55.1, 55.1, 47.1, 38.1, 45.1, 45.1, 37.1, 35.1, 35.1, 35.1,
               44.5, 44.5, 44.5, 44.5, 37, 37, 29.5, 222, 222, 222, 222, 222,
               222, 222, 222, 222, 222, 222, 148, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0],

        "txt": [("CH3", {"C": 1, "H": 3}),                                # 0
                ("CH2", {"C": 1, "H": 2}),
                ("CH", {"C": 1, "H": 1}),
                ("C", {"C": 1}),
                ("=CH2", {"C": 1, "H": 2}),
                ("=CH-", {"C": 1, "H": 1}),
                ("=C", {"C": 1}),
                ("=C=", {"C": 1}),
                ("≡CH", {"C": 1, "H": 1}),
                ("≡C-", {"C": 1}),
                ("-CH2- (Cyclic)", {"C": 1, "H": 2}),                     # 10
                ("-CH< (Ciclic)", {"C": 1, "H": 1}),
                ("-CH< (in fused ring)", {"C": 1, "H": 1}),
                (">C< (Ciclic)", {"C": 1}),
                ("=CH- (Cyclic)", {"C": 1, "H": 1}),
                ("=C< (Cyclic)", {"C": 1}),
                ("=C= (Cyclic)", {"C": 1}),
                ("Phenyl- ", {"C": 6, "H": 5}),
                ("o-Phenyl- ", {"C": 6, "H": 4}),
                ("m-Phenyl- ", {"C": 6, "H": 4}),
                ("p-Phenyl- ", {"C": 6, "H": 4}),                         # 20
                ("1,2,3-Phenyl- ", {"C": 6, "H": 3}),
                ("1,2,4-Phenyl- ", {"C": 6, "H": 3}),
                ("1,2,3,4-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,3,5-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,4,5-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,3,4,5-Phenyl- ", {"C": 6, "H": 1}),
                ("1,2,4,5,6-Phenyl- ", {"C": 6}),
                ("=CH-CH= (in fused Aromatic ring)", {"C": 2, "H": 2}),
                ("-OH", {"H": 1, "O": 1}),
                ("-O-", {"O": 1}),                                        # 30
                ("-CO-", {"C": 1, "O": 1}),
                ("-CHO", {"C": 1, "H": 1, "O": 1}),
                ("-COOH", {"C": 1, "H": 1, "O": 2}),
                ("-CO-O-OC-", {"C": 2, "O": 3}),
                ("-CO-O-", {"C": 1, "O": 2}),
                ("-NO2", {"N": 1, "O": 2}),
                ("-NH2", {"N": 1, "H": 2}),
                ("-NH-", {"N": 1, "H": 1}),
                ("-N<", {"N": 1}),
                ("-CN", {"N": 1, "C": 1}),                                # 40
                ("-S-", {"S": 1}),
                ("-SH", {"S": 1, "H": 1}),
                (">Si<", {"Si": 1}),
                (">SiH-", {"Si": 1, "H": 1}),
                ("-SiH3", {"Si": 1, "H": 3}),
                (">SiO-", {"Si": 1, "O": 1}),
                (">SiO- (cyclic)", {"Si": 1, "O": 1}),
                ("-F", {"F": 1}),
                ("-Cl", {"Cl": 1}),
                ("-Br", {"Br": 1}),                                       # 50
                ("-OH (Aromatic)", {"O": 1, "H": 1}),
                ("C4H4 (fused ring)", {"C": 4, "H": 4}),
                ("-F (Aromatic)", {"F": 1}),
                ("-Cl (Aromatic)", {"Cl": 1}),
                ("-Br (Aromatic)", {"Br": 1}),
                ("-I (Aromatic)", {"I": 1}),

                # 2nd Order terms
                ("Double Bond", ),
                ("Triple Bond", ),
                ("Benzene", ),
                ("Pyridine", ),                                           # 60
                ("The single or first substituent on an aromatic ring", ),
                ("The second or subsequent substituent on an aromatic ring", ),
                ("Each pair of rign substituents in ortho positions", ),
                ("If one of the ortho pair is -OH", )]}

    FirstOrder = 57
    SecondOrder = 65

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = translate("newComponent", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def _group(self):
        """From group contribution desglose the chemical composition"""
        group = []
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            # Only the first order term count for this
            if i < self.FirstOrder:
                grp = self.__coeff__["txt"][i][1]
                for x in range(c):
                    group.append(grp)
        return group

    def calculo(self):
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        Pc = tc = vc = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i == 29:
                # Use special case for -OH group
                n = 0
                for g in self.group:
                    n += g.get("C", 0)

                tc += 0.87-0.11*n+0.003*n**2
                Pc += 0.1-0.013*n
                vc += 0
            else:
                tc += c*self.__coeff__["tc"][i]
                Pc += c*self.__coeff__["Pc"][i]
                vc += c*self.__coeff__["vc"][i]

        Pt = self.kwargs["platt"]
        self.Tc = unidades.Temperature(self.Tb*(1+1/(1.242+tc-0.023*Pt)))
        self.Pc = unidades.Pressure(14.5*self.M/(0.339+Pc-0.026*Pt)**2, "psi")
        self.Vc = unidades.SpecificVolume(0.01602*(40+vc)/self.M, "ft3lb")

        GroupContribution.calculo(self)
