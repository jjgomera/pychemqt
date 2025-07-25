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
from lib.compuestos import atomic_decomposition
from lib.newComponent._base import GroupContribution


class Wen(GroupContribution):
    """
    Group contribution for definition of unknown component using the Wen-Qiang
    procedure (2001). This method is able to calculate the critical properties.

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
    Example 1 in [1]_, Tc of n-Butylaniline
    The last group containing a carbon-adjacent atom pair has a typo, must be
    =C<[r]/>N-

    >>> cmp = Wen(Tb=513.9, group=[1, 11, 18, 87, 96, 99, 135],
    ... contribution=[1, 5, 1, 10, 2, 1, 1])
    >>> "%0.3f" % cmp.Tc
    '721.275'
    >>> cmp.formula
    'C10H15N'

    Example 2 in [1]_, Pc of Benzoic acid

    >>> cmp = Wen(group=[87, 93, 96, 114], contribution=[10, 1, 1, 1])
    >>> "%0.3f" % cmp.Pc.MPa
    '4.547'

    Example 3 in [1]_, Vc of chloropentafluorobenzene

    >>> cmp = Wen(group=[96, 101, 102], contribution=[12, 5, 1])
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '374.6'
    >>> cmp.formula
    'C6F5Cl'
    """
    __title__ = "Wen-Qiang (2001)"

    __doi__ = (
        {"autor": "Wen, X., Quiang, Y.",
         "title": "A New Group Contribution Method for Estimating Critical "
                  "Properties of Orgnic Compounds",
         "ref": "Ind. Eng. Chem. Res. 40(26) (2001) 6245-6250.",
         "doi": "10.1021/ie010374g"},)

    __coeff__ = {
        # Table III
        "tc": [-2.885, 2.424, 0.048, 22.766, -3.404, 2.495, 2.275, 2.602,
               -1.601, 0.000, 35.848, 2.124, -0.708, 22.576, -3.085, 1.578,
               -0.030, 2.256, -2.322, 2.549, 13.769, 20.882, 25.177, 1.934,
               -3.377, 0.000, 0.000, -1.765, -3.163, -5.588, 0.830, 11.483,
               2.183, -7.415, 0.007, -0.651, 0.178, -4.384, -0.502, 6.664,
               9.639, 2.262, -1.648, 0.000, 5.838, 2.218, 10.659, -2.228,
               0.017, 3.541, -0.748, -1.248, 9.056, 4.564, 2.737, 0.007, 1.136,
               -0.008, 13.166, -0.009, -2.298, 9.242, 1.143, -11.918, 0.007,
               -0.953, 0.000, 0.000, 5.993, 3.162, 0.000, 2.707, 3.180, 6.070,
               13.625, 7.842, 5.897, 0.004, -781.237, 2.367, 0.000, -7.274,
               3.798, 5.571, 2.446, -1.223, 0.754, 2.852, 2.013, 7.937, 14.661,
               10.141, -0.603, 2.172, 0.000, 0.009, 4.660, -2.465, 3.701,
               1.700, 6.344, -5.547, 5.600, 12.840, 28.472, 10.144, 18.220,
               15.436, 20.655, 11.326, 13.047, 34.349, 35.591, 34.476, 35.009,
               7.041, 9.909, 9.193, 1.933, 16.126, 23.640, 9.697, 15.542,
               3.835, 11.598, 4.679, 22.777, 18.847, 8.357, 0.003, 11.032,
               14.462, 16.658, 1.251, 15.982, 17.609, 16.935, 16.840, 0.001,
               0.004],
        "Pc": [9.361, 4.035, 2.682, 7.267, 9.548, 3.297, 1.286, 2.797, 5.280,
               0.000, 2.517, 1.454, 0.314, 5.065, 9.483, 2.239, 0.038, 1.497,
               3.545, 1.851, 1.707, 1.517, 3.345, 0.460, 0.358, 0.000, 6.274,
               1.087, 0.187, 3.501, 2.770, 1.318, -0.282, 0.141, 6.690, -0.032,
               -0.013, 0.214, 2.885, 2.880, 2.604, 3.648, 4.393, 0.000, 2.686,
               0.731, 5.158, 3.462, 0.000, -0.290, -1.534, 3.842, 2.532,
               1.624, 0.305, 0.000, 2.352, 0.000, -0.864, -0.854, 3.177,
               2.036, 1.668, 0.000, -0.001, -3.595, 0.000, 0.000, -1.661,
               -0.802, 5.880, 1.006, 2.175, 0.535, 0.000, -0.152, -4.382,
               0.000, 0.000, 0.464, 0.000, 4.682, 1.848, -5.905, 0.883, 1.355,
               -0.296, 0.680, 0.853, 0.275, -0.001, 0.189, 1.671, 2.011, 4.716,
               -0.004, -0.064, -0.225, -1.345, 5.744, -0.319, 3.029, 3.476,
               4.375, 3.073, 1.662, 2.422, 1.621, 0.037, 1.516, -0.543, -2.792,
               1.524, -0.406, 0.855, -0.026, 2.966, 1.840, 2.188, 1.785,
               -0.332, -0.622, 5.780, 0.078, 0.001, 1.800, 0.000, -1.320,
               4.158, 0.005, -1.452, -2.250, -0.001, 0.000, -3.152, -4.585,
               -7.394, 0.000, 0.000, -0.002],
        "vc": [125.58, 86.72, 70.64, 105.31, 45.38, 91.08, 62.82, 113.65,
               49.86, 0.00, 200.24, 28.22, 14.98, 44.98, 0.07, 46.54, 13.89,
               57.95, -20.09, 81.10, 103.74, 127.10, 0.00, 5.53, 2.54, 0.00,
               0.00, 16.98, 1.47, -46.51, 55.26, 84.21, -4.52, -16.38, 0.00,
               0.00, -0.08, -3.45, 44.84, 75.12, 98.04, 78.05, 80.34, 0.00,
               38.30, 23.00, 136.00, 28.40, 0.00, 0.00, 21.77, 68.06, 98.13,
               20.29, 5.01, 0.00, 3.48, 0.00, 62.65, -99.32, 48.52, 77.32,
               69.56, 0.01, -95.79, 18.23, 0.00, 0.00, 14.19, 77.52, 0.00,
               27.84, 46.55, 38.22, 0.00, 50.59, -0.02, 0.00, 0.00, 34.33,
               0.00, 0.00, 36.94, 0.04, 36.83, 4.71, 0.00, 23.60, 51.71, 52.22,
               0.00, 45.09, -2.49, 112.28, 0.00, 0.00, 9.87, -49.16, 5.85,
               -64.77, 5.54, 40.93, 78.55, 95.28, 122.28, 0.00, 113.89, 0.00,
               0.00, 46.35, -81.99, 0.00, 127.31, 112.99, 0.00, 74.01, 60.54,
               39.13, -68.36, 0.00, -28.91, 17.46, 165.84, 94.53, -0.01, 0.00,
               0.00, 58.88, 46.41, 0.01, 101.11, 59.67, 0.00, 0.00, 110.06,
               138.33, 131.47, 0.00, 0.00, -0.01, ],
        "tc_": [-8.8072, -1.1863, 2.0695, 0.8880, 0.4312, -2.9673, -5.6886,
                -1.8098, -0.5794, -0.0174, 4.5913, 2.2116, 5.2478, 0.5832,
                3.3562, -0.2371, -0.0296, 2.8167, 2.9084, 1.0111, -0.8130,
                -1.3763, -0.0022, 2.7600, 7.5553, 0.0000, 0.0000, 0.3694,
                0.9548, 3.9991, 1.6760, 0.5273, 2.7770, 5.3585, -0.0076,
                0.8857, -0.0093, 0.7629, 0.7647, 1.6201, 0.4661, 1.7766,
                4.6760, 0.0000, -1.4353, 1.5808, -5.4993, 13.0182, -0.1973,
                -5.5762, -2.5308, -1.7201, -0.9027, 0.4717, 0.9585, -0.0014,
                6.5524, 0.0038, -4.2721, 0.0060, 1.3320, 0.6099, -1.1118,
                -6.0111, -0.0038, 4.4394, 0.0000, 0.0000, -6.3639, -0.8632,
                0.0000, 0.5719, 1.2736, -0.1301, -0.9306, 0.1783, 2.0133,
                0.0204, 0.0000, 0.0094, 0.0000, 0.7476, 2.4465, 4.1358,
                -1.5784, 4.0418, -1.1378, 0.3159, -0.1566, -1.1410, 0.4046,
                -0.0075, 5.6200, 7.6978, 0.0000, -0.0180, 1.2080, -4.5849,
                0.8675, 7.0762, 2.773, 0.7764, 1.3463, 0.2835, 0.0770, 0.6733,
                -1.1847, 5.4387, -4.3184, 1.0536, -3.2089, 6.8866, 12.5998,
                14.3778, 10.1056, 2.6863, 5.4168, 11.107, 8.7813, 15.4803,
                -0.4178, 2.5145, 9.8751, 9.4595, 15.4000, 1.3670, 7.2107,
                12.1779, 11.7647, 0.0013, 6.3813, -4.7655, 2.7693, -0.5727,
                1.4762, 4.2350, 7.5331, 5.1503, 0.0019, 0.0590],

        # Name and element composition
        "txt": [("CH3- & H",),                          # 0
                ("CH3- & >C<",),
                ("CH3- & =C<",),
                ("CH3- & ≡C-",),
                ("CH3- & >C< [r]",),
                ("CH3- & =C< [r]",),
                ("CH3- & -O-",),
                ("CH3- & -S-",),
                ("CH3- & >N-",),
                ("CH3- & =N-",),
                ("CH3- & -NO2",),                        # 10
                ("-CH2- & >C<",),
                ("-CH2- & =C<",),
                ("-CH2- & ≡C-",),
                ("-CH2- & >C< [r]",),
                ("-CH2- & =C< [r]",),
                ("-CH2- & -O-",),
                ("-CH2- & -S-",),
                ("-CH2- & >N-",),
                ("-CH2- & F-",),
                ("-CH2- & Cl-",),                        # 20
                ("-CH2- & Br-",),
                ("-CH2- & I-",),
                (">CH- & >C<",),
                (">CH- & =C<",),
                (">CH- & ≡C-",),
                (">CH- & >C< [r]",),
                (">CH- & =C< [r]",),
                (">CH- & -O-",),
                (">CH- & >N-",),
                (">CH- & F-",),                          # 30
                (">CH- & Cl-",),
                (">C< & >C<",),
                (">C< & =C<",),
                (">C< & ≡C-",),
                (">C< & >C< [r]",),
                (">C< & =C< [r]",),
                (">C< & -O-",),
                (">C< & F-",),
                (">C< & Cl-",),
                (">C< & Br-",),                          # 40
                ("=CH2 & =C<",),
                ("=CH2 & =C=",),
                ("=CH2 & =C< [r]",),
                ("=CH- & >C<",),
                ("=CH- & =C<",),
                ("=CH- & ≡C-",),
                ("=CH- & =C=",),
                ("=CH- & >C< [r]",),
                ("=CH- & =C< [r]",),
                ("=CH- & -O-",),                         # 50
                ("=CH- & F-",),
                ("=CH- & Cl-",),
                ("=C< & >C<",),
                ("=C< & =C<",),
                ("=C< & ≡C-",),
                ("=C< & =C=",),
                ("=C< & =C< [r]",),
                ("=C< & =O",),
                ("=C< & >N-",),
                ("=C< & F-",),                           # 60
                ("=C< & Cl-",),
                ("≡CH & ≡C-",),
                ("≡C- & >C<",),
                ("≡C- & =C<",),
                ("≡C- & ≡C-",),
                ("≡C- & >C< [r]",),
                ("≡C- & =C< [r]",),
                ("=C= [r] & =C<",),
                ("=C= [r] & =O",),
                ("=C= [r] & =N-",),                      # 70
                ("-CH2- [r] & >C< [r]",),
                ("-CH2- [r] & =C< [r]",),
                ("-CH2- [r] & -O- [r]",),
                ("-CH2- [r] & -S- [r]",),
                ("-CH2- [r] & >N- [r]",),
                (">CH- [r] & >C<",),
                (">CH- [r] & =C<",),
                (">CH- [r] & ≡C-",),
                (">CH- [r] & >C< [r]",),
                (">CH- [r] & =C< [r]",),                  # 80
                (">CH- [r] & -O-",),
                (">CH- [r] & -O- [r]",),
                (">C< [r] & >C<",),
                (">C< [r] & >C< [r]",),
                (">C< [r] & F-",),
                ("=CH- [r] & >C< [r]",),
                ("=CH- [r] & =C< [r]",),
                ("=CH- [r] & -O- [r]",),
                ("=CH- [r] & -S- [r]",),
                ("=CH- [r] & >N- [r]",),                  # 90
                ("=CH- [r] & =N- [r]",),
                ("=C< [r] & >C<",),
                ("=C< [r] & =C<",),
                ("=C< [r] & ≡C-",),
                ("=C< [r] & >C< [r]",),
                ("=C< [r] & =C< [r]",),
                ("=C< [r] & -O-",),
                ("=C< [r] & -O- [r]",),
                ("=C< [r] & >N-",),
                ("=C< [r] & >N- [r]",),                   # 100
                ("=C< [r] & F-",),
                ("=C< [r] & Cl-",),
                ("=C< [r] & Br-",),
                ("=C< [r] & I-",),
                ("-CHO & -H",),
                ("-CHO & >C<",),
                ("-CHO & =C<",),
                ("-CHO & =C< [r]",),
                ("-CO & >C<",),
                ("-CO & =C< [r]",),                       # 110
                ("-COOH & -H",),
                ("-COOH & >C<",),
                ("-COOH & =C<",),
                ("-COOH & =C< [r]",),
                ("-COO- & -H",),
                ("-COO- & >C<",),
                ("-COO- & =C<",),
                ("-COO- & =C< [r]",),
                (">C2O3 & =C<",),
                (">C2O3 & =C< [r]",),                     # 120
                (">CO [r] & >C< [r]",),
                ("-CN & -H",),
                ("-CN & >C<",),
                ("-CN & =C<",),
                ("-CN & ≡C-",),
                ("-CN & =C< [r]",),

                # Group Structural Unit Parameter
                ("-OH", ),
                ("-O-", ),
                ("-O- (cyclic)", ),
                ("-OH (Aromatic)", ),                     # 130
                ("-SH", ),
                ("-S-", ),
                ("-S- (cyclic)", ),
                ("-NH2", ),
                (">NH", ),
                (">N-", ),
                ("=N-", ),
                (">NH (cyclic)", ),
                (">N- (cyclic)", )]}

    FirstOrder = 127
    SecondOrder = 140

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        tc, pc, vc, tc_ = 0, 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self.__coeff__["tc"][i]
            tc_ += c*self.__coeff__["tc_"][i]
            pc += c*self.__coeff__["Pc"][i]
            vc += c*self.__coeff__["vc"][i]

        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
            Tc = (((127.754+tc_)/1e2)**-2+1)*self.Tb
            self.Tc = unidades.Temperature(Tc)
        else:
            self.Tc = unidades.Temperature(((4.72+tc)*1e6)**(1/2.747))
        self.Pc = unidades.Pressure(((37.293+pc)/1e2)**-2, "MPa")
        self.Vc = unidades.SpecificVolume((-27.04+vc)/self.M, "ccg")

        GroupContribution.calculo(self)

    def _group(self):
        """Specific procedure to calculate the molecular weight of compound
        from group contribution"""
        # TODO: Add group check coherent definition
        group = []
        rest = {}
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if "&" in self.__coeff__["txt"][i][0]:
                grp = self.__coeff__["txt"][i][0].split(" & ")[0]
                second = self.__coeff__["txt"][i][0].split(" & ")[1]

                # Discard second term with carbons and add the heteroatoms term
                grp2 = atomic_decomposition(second)
                if "H" in grp2 or "F" in grp2 or "Cl" in grp2 or \
                        "Br" in grp2 or "I" in grp2:
                    for x in range(c):
                        group.append(grp2)
            else:
                grp = self.__coeff__["txt"][i][0]

            for x in range(c):
                cmp = grp
                molecule = grp
                # Clean additional comment of group, ring flag and other,
                # separated of main group by spaces
                if " " in cmp:
                    cmp = cmp.split(" ")[0]

                # Calculate the remain appeareances of group contributions
                # X- : X is one connected so don't appearance reamaining
                # -X- : X has 2 connection so must be stay in another group
                # -X<: X has 3 connection so must be stay in other 2 group
                # >X<: X has 4 connection so must be stay in other 3 group
                restLink = -1
                if cmp[0] in ("-=≡"):
                    restLink += 1
                    molecule = molecule[1:]
                elif cmp[0] == ">":
                    restLink += 2
                    molecule = molecule[1:]
                if cmp[-1] in ("-=≡"):
                    restLink += 1
                    molecule = molecule[:-1]
                elif cmp[-1] == "<":
                    restLink += 2
                    molecule = molecule[:-1]

                if grp in rest and rest[grp]:
                    rest[grp] -= 1
                else:
                    rest[grp] = restLink
                    group.append(atomic_decomposition(molecule))
        return group
