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
from lib.compuestos import atomic_decomposition
from lib.newComponent._base import GroupContribution


class Marrero(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Marrero-Pardillo procedure (1999). This method is able to calculate the
    critical properties, boiling temperature and viscosity.

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound, using general properties for calculation of other
    mandatory properties don't defined by the method.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    Tb : float, optional
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2-6 in [1]_, 2-ethylphenol critical properties

    >>> cmp = Marrero(Tb=477.67, group=[1, 36, 129, 130, 132, 140, 148],
    ... contribution=[1, 1, 1, 2, 2, 1, 1])
    >>> "%0.1f %0.1f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '699.8 42.2 378.8'
    >>> cmp.formula
    'C8H10O'

    Example 2-7 in [1]_, butanols

    >>> b1 = Marrero(group=[1, 28, 41], contribution=[1, 2, 1])
    >>> b2m1 = Marrero(group=[2, 29, 41], contribution=[2, 1, 1])
    >>> b2m2 = Marrero(group=[3, 79], contribution=[3, 1])
    >>> b2 = Marrero(group=[1, 2, 29, 62], contribution=[1, 1, 1, 1])

    >>> "%0.2f %0.2f %0.2f %0.2f" % (b1.Tc, b2m1.Tc, b2m2.Tc, b2.Tc)
    '565.67 555.80 504.45 534.96'
    >>> "%0.2f %0.2f %0.2f %0.2f" % (
    ... b1.Pc.bar, b2m1.Pc.bar, b2m2.Pc.bar, b2.Pc.bar)
    '44.86 45.04 41.30 43.40'
    >>> "%0.1f %0.1f %0.1f %0.1f" % (b1.Vc.ccg*b1.M, b2m1.Vc.ccg*b2m1.M,
    ...  b2m2.Vc.ccg*b2m2.M, b2.Vc.ccg*b2.M)
    '272.1 267.3 275.1 277.9'

    Example 2-13 in [1]_, 2,4-dimethylphenol

    >>> cmp = Marrero(group=[9, 130, 132, 133, 140, 148],
    ... contribution=[2, 3, 1, 1, 1, 1], M=122.167)
    >>> "%0.2f" % cmp.Tb
    '489.39'

    Example 2-14 in [1]_, cycloalkanes

    >>> c7 = Marrero(group=[111], contribution=[7])
    >>> mc6 = Marrero(group=[7, 111, 112], contribution=[1, 4, 2])
    >>> ec5 = Marrero(group=[1, 34, 111, 112], contribution=[1, 1, 3, 2])
    >>> c5 = Marrero(group=[7, 111, 112], contribution=[2, 1, 4])
    >>> t5 = Marrero(group=[7, 111, 112], contribution=[2, 1, 4])
    >>> "%0.2f %0.2f %0.2f %0.2f %0.2f" % (c7.Tb, mc6.Tb, ec5.Tb, c5.Tb, t5.Tb)
    '377.53 375.85 384.47 374.17 374.17'

    Table 8a in [2]_ for 1,3,5-Trichlorotrifluorobenzene

    >>> cmp = Marrero(group=[138, 140, 144, 145], contribution=[3, 3, 3, 3])
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '471.9 685.0 32.4 452.3'

    Table 8b in [2]_ for ethyl acrylate

    >>> cmp = Marrero(group=[84, 1, 46, 100], contribution=[1, 1, 1, 1])
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '373.9 553.6 36.6 325.1'

    Table 8c in [2]_ for isoquinoline

    >>> cmp = Marrero(group=[129, 138, 131, 136, 132, 133],
    ... contribution=[3, 1, 1, 1, 1, 4])
    >>> "%0.1f %0.1f %0.1f" % (cmp.Tb, cmp.Tc, cmp.Vc.ccg*cmp.M)
    '519.5 787.3 405.1'

    Table 8d in [2]_ for m-Terphenyl

    >>> cmp = Marrero(group=[129, 130, 141, 132, 133],
    ... contribution=[5, 4, 2, 5, 4], Tb=638)
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '638.0 907.3 33.9 764.3'

    Table 3 in [3]_ for o-phthalate

    >>> cmp = Marrero(group=[129, 132, 133, 138, 153, 46, 28, 1],
    ... contribution=[2, 1, 2, 1, 2, 2, 4, 2])
    >>> "%0.2f" % cmp.mu.muPas
    '18.58'
    """
    __title__ = "Marrero-Pardillo (1999)"
    __doi__ = (
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
        {"autor": "Marrero-Morejón, J., Pardillo-Fontdevila, F.",
         "title": "Estimation of Pure Compound Properties Using "
                  "Group-Interaction C9ontributions",
         "ref": "AIChE J., 45(3) (1999) 615-621",
         "doi": "10.1002/aic.690450318"},
        {"autor": "Marrero-Morejon, J., Pardillo-Fontdevila, E.",
         "title": "Estimation of Liquid Viscosity at Ambient Temperature of "
                  "Pure Organic Compounds by Using Group-Interaction "
                  "Contributions",
         "ref": "Chemical Engineering Journal 79 (2000) 69-72",
         "doi": "10.1016/s1385-8947(99)00173-4"})

    __coeff__ = {
        # Table 5 from [2]_
        "tc": [-0.0213, -0.0227, -0.0223, -0.0189, 0.8526, 0.1792, 0.3818,
               -0.0214, 0.1117, 0.0987, -0.0370, -0.9141, -0.9166, -0.9146,
               -0.0876, -0.0205, -0.0362, -0.0606, -0.0890, 0.0267, -0.0974,
               -0.0397, -0.0313, -0.0199, -0.0766, -0.0591, -0.9192, -0.0181,
               -0.0206, -0.0134, -0.0098, 0.8636, 0.1874, 0.4160, -0.0149,
               0.1193, 0.1012, -0.0255, -0.0162, -0.0205, -0.0210, -0.0786,
               -0.0205, -0.0256, -0.0267, -0.0932, 0.0276, -0.0993, -0.0301,
               -0.0248, -0.0161, -0.0654, -0.0137, -0.0192, -0.0039, 0.0025,
               0.8547, 0.1969, 0.0025, 0.1187, -0.0200, -0.0142, -0.0757,
               -0.0162, -0.0194, -0.0406, -0.0918, -0.1054, -0.0286, -0.0158,
               0.0084, 0.8767, 0.2061, 0.0207, 0.0049, 0.1249, -0.0176,
               -0.0133, -0.0084, -0.0780, -0.0156, -0.0114, -0.1008, -0.9129,
               -0.8933, -0.4158, -0.0123, -1.7660, -1.2909, -0.8945, 1.7377,
               1.0731, 1.2865, 0.9929, 0.8623, 0.8613, 0.8565, 0.8246, 0.7862,
               0.8818, 0.7780, 0.8122, -0.8155, -0.4009, 0.3043, 0.1868,
               0.1886, -0.0159, -0.0288, -0.4222, -0.7958, -0.0098, -0.0093,
               -0.1386, 0.0976, 0.1089, -0.0092, -0.0148, -0.0139, -0.0071,
               -0.0055, -0.1341, 0.0, 0.0, -0.0218, -0.0737, 0.0329, 0.0,
               -0.0314, -0.2246, -0.3586, 0.3913, 0.2089, 0.2190, 0.1000,
               0.0947, -0.4067, 0.1027, -0.4848, 0.2541, 0.2318, 0.2424,
               0.1104, -0.3972, 0.1069, 0.1028, 0.1060, 0.1075, 0.0931, 0.0997,
               0.1112, 0.0919, 0.0313, 0.0241, 0.0830, 0.0978, 0.0938, 0.0768,
               -0.0191, -0.1926, -0.5728, -0.3553, -0.0422, -0.0690, -0.0781,
               -0.0301, -0.0124],
        "Pc": [-0.0618, -0.0430, -0.0376, -0.0354, 0.0654, 0.0851, -0.2320,
               -0.0396, -0.0597, -0.0746, -0.0345, -0.0231, -0.0239, -0.0241,
               -0.0180, -0.0321, -0.0363, -0.0466, -0.0499, 0.1462, -0.2290,
               -0.0288, -0.0317, -0.0348, -0.0507, -0.0385, -0.0244, -0.0305,
               -0.0272, -0.0219, -0.0162, 0.0818, 0.1010, -0.2199, -0.0265,
               -0.0423, -0.0626, -0.0161, -0.0150, -0.0140, -0.0214, -0.0119,
               -0.0184, -0.0204, -0.0210, -0.0253, 0.1561, -0.2150, -0.0214,
               -0.0203, -0.0170, -0.0329, -0.0163, -0.0173, -0.0137, -0.0085,
               0.0816, 0.1080, -0.0168, -0.0556, -0.0147, -0.0131, -0.0093,
               -0.0155, -0.0112, -0.0280, -0.2098, -0.0358, -0.0212, -0.0162,
               0.0002, 0.0953, 0.1109, 0.0213, -0.0111, -0.0510, -0.0161,
               -0.0129, -0.0121, -0.0094, -0.0103, -0.0085, -0.0455, -0.0476,
               -0.1378, -0.2709, -0.0239, -0.2291, -0.3613, -0.1202, 0.1944,
               0.2146, -0.1087, 0.0533, 0.0929, 0.0919, 0.0947, 0.0801, 0.0806,
               0.2743, -0.1007, 0.0771, -0.4920, -0.2502, 0.0705, 0.1064,
               0.1102, -0.0010, -0.0226, 0.1860, 0.3933, -0.0221, -0.0181,
               0.0081, -0.1034, -0.0527, -0.0119, -0.0177, -0.0127, 0.0,
               -0.0088, 0.0162, 0.0, 0.0, -0.0091, -0.0220, -0.0071, 0.0,
               -0.0119, 0.1542, 0.1490, 0.1356, -0.1822, -0.1324, -0.0900, 0.0,
               -0.1491, -0.0916, 0.1432, 0.0, -0.0809, -0.0792, -0.0374,
               -0.0971, -0.0504, -0.0512, -0.0548, -0.0514, -0.0388, -0.0523,
               -0.0528, -0.0597, -0.0684, -0.2573, -0.0579, -0.0471, -0.0462,
               -0.0625, -0.0125, -0.0878, 0.0, -0.0176, -0.0123, 0.0, -0.1878,
               0.0, 0.0],
        "vc": [123.2, 88.6, 78.4, 69.8, 81.5, 57.7, 65.8, 58.3, 49.0, 71.7,
               88.1, 113.8, 0.0, 0.0, 92.9, 66.0, 88.9, 128.9, 145.9, 93.3,
               108.2, 0.0, 0.0, 76.3, 147.9, 148.1, 119.7, 87.9, 56.6, 40.2,
               32.0, 50.7, 24.0, 33.9, 31.9, 0.0, 52.1, 49.3, 80.8, 101.3, 0.0,
               45.2, 34.5, 62.3, 106.1, 114.0, 69.9, 79.1, 63.3, 49.4, 32.7,
               113.5, 93.3, 57.9, 18.3, 8.6, 48.9, 4.3, 0.0, 0.0, 37.7, 68.6,
               45.6, 23.7, 39.3, 92.2, 72.3, 110.2, 39.2, 0.0, 22.7, 23.4, 8.8,
               0.0, 0.0, 0.0, 30.0, 63.7, 85.7, 40.6, 40.8, 62.1, 89.0, 105.3,
               77.4, 99.2, 68.4, 47.8, 73.6, 43.6, 42.1, 16.6, 0.0, 0.0, 41.4,
               68.7, 36.4, 0.0, 107.4, 55.2, 64.1, 107.4, 93.7, 58.1, 0.0,
               14.6, 43.3, 51.4, 87.6, 73.1, 64.3, 47.2, 47.5, 49.9, 42.5, 0.0,
               29.2, 50.7, 38.8, 0.0, 33.9, 0.0, 0.0, 0.0, 19.2, 0.0, 36.2,
               0.0, 18.4, 36.5, 34.4, 8.3, 39.3, 29.8, 40.3, 0.0, 65.9, 40.8,
               37.8, 0.0, 20.6, 51.7, -0.3, 35.6, 23.7, 60.3, 83.2, 110.2, 8.5,
               0.0, 46.3, 0.0, 100.2, 55.2, 33.2, 0.0, 0.0, 0.0, 84.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 51.2, 0.0, 0.0],
        "tb": [113.12, 194.25, 194.27, 186.41, 137.18, 182.20, 194.40, 176.16,
               180.60, 145.56, 160.83, 453.70, 758.44, 1181.44, 736.93, 228.01,
               445.61, 636.49, 1228.84, 456.92, 510.65, 443.76, 293.86, 207.75,
               891.15, 1148.58, 588.31, 409.85, 244.88, 244.14, 273.26, 201.80,
               242.47, 207.49, 238.81, 260.00, 167.85, 166.59, 517.62, 875.85,
               1262.80, 673.24, 243.37, 451.27, 648.70, 1280.39, 475.65,
               541.29, 452.30, 314.71, 240.08, 869.18, 612.31, 451.03, 291.41,
               344.06, 179.96, 249.10, 295.33, 132.66, 68.80, 438.47, 585.99,
               215.94, 434.45, 630.07, 497.58, 1270.16, 388.44, 260.32, 411.56,
               286.30, 286.42, 456.90, 340.00, 188.99, -16.64, 360.79, 610.26,
               540.38, 267.26, 373.71, 1336.54, 51.13, 205.73, 245.27, 183.55,
               334.64, 354.41, 316.46, 174.18, 228.38, 174.39, 184.20, 5.57,
               370.60, 204.81, 658.53, 1245.86, 423.86, 525.35, 761.36, 399.58,
               321.02, 250.88, -37.99, 367.05, 160.42, 120.85, 222.40, 333.26,
               201.89, 209.40, 182.74, 218.07, 106.21, 225.52, 451.74, 283.55,
               424.13, 210.66, 220.24, 254.50, 184.36, 169.17, 597.82, 348.23,
               111.51, -41.35, 112.00, 291.15, 221.55, 285.07, 237.22, 171.59,
               420.54, 321.44, 348.00, 477.77, 334.09, 180.07, 123.05, 134.23,
               174.31, -48.79, 347.33, 716.23, 1294.98, 456.25, 199.70, 437.51,
               700.06, 1232.55, 437.78, 517.75, 411.29, 422.51, 682.19, 532.24,
               1012.51, 382.25, 385.36, 387.17, 1022.45, 298.12, 673.59,
               597.59],

        # Table 1 and 2 from [3]_
        "mu": [0.4712, 0.2588, 0.2472, 0.1833, 0.0834, 0.1670, 0, -0.0926, 0,
               0.0314, 0, 0.6294, 0.7161, 0.6086, 2.6902, 0.4020, 0.5475,
               1.3881, 2.6110, 0.6410, 0.5904, 0, 1.3476, 0.1695, 1.9082,
               1.9476, 0, 0.5587, 0.0975, -0.0543, -0.0417, -0.0775, -0.0488,
               0, -0.4464, 0, -0.4327, 0, 0.7214, 0.5799, 0.3971, 2.9295,
               0.1089, 0.2605, 1.0415, 1.9798, 0.3632, 0.2089, 1.5471, 0.2857,
               -0.0047, 1.9082, 0.7255, 0.2785, -0.1470, 0.1181, 0, 0, 0,
               -0.0701, 0, 0.3633, 2.9212, -0.1711, 0.2242, 0, -0.0606, 1.6883,
               0.6492, 0, 0, 0, 0, 0, 0, 0, -0.1491, 0.2947, 0, 2.8610, 0, 0,
               0, 0, 0.2894, 0.3306, 0, 0.3096, 0.1060, 0, 0, -0.0465, 0,
               -0.2144, 0, 0.3593, -0.0982, 1.1338, 0, 0.1966, 1.1345, 1.2214,
               0.2482, 0, 0, 0, 0.1831, 0, 0, 0, 0, 0.2980, 0.2136, 0.1692,
               0.2675, -0.0704, 0.3781, 0.6534, 0.6341, 0.5436, -0.3205, 0, 0,
               0, 0.8726, 3.1712, 0, 0, 0, 0.0952, 0.2486, 0.5734, 0.3422,
               0.2130, 0.3708, 1.0147, 0.5751, 0.5128, 0.4700, 0.7030, 0.2292,
               0, 0, 0.1880, -0.1932, 0.0404, -0.0847, 0.0139, 4.3481, 0.0963,
               0.2243, 0.7089, 0, 0.1514, 1.7430, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0],

        # Custom dict for molecular weight and empiric formula calculation
        "txt": [("CH3- & CH3-",),          # 0
                ("CH3- & -CH2-",),
                ("CH3- & >CH-",),
                ("CH3- & >C<",),
                ("CH3- & =CH-",),
                ("CH3- & =C<",),
                ("CH3- & ≡C-",),
                ("CH3- & >CH- [r]",),
                ("CH3- & >C< [r]",),
                ("CH3- & =C< [r]",),
                ("CH3- & F-",),            # 10
                ("CH3- & Cl-",),
                ("CH3- & Br-",),
                ("CH3- & I-",),
                ("CH3- & -OH",),
                ("CH3- & -O-",),
                ("CH3- & >CO",),
                ("CH3- & -CHO",),
                ("CH3- & -COOH",),
                ("CH3- & -COO[-]",),
                ("CH3- & [-]COO-",),       # 20
                ("CH3- & -NH2",),
                ("CH3- & -NH-",),
                ("CH3- & >N-",),
                ("CH3- & -CN",),
                ("CH3- & -NO2",),
                ("CH3- & -SH",),
                ("CH3- & -S-",),
                ("-CH2- & -CH2-",),
                ("-CH2- & >CH-",),
                ("-CH2- & >C<",),          # 30
                ("-CH2- & =CH-",),
                ("-CH2- & =C<",),
                ("-CH2- & ≡C-",),
                ("-CH2- & >CH- [r]",),
                ("-CH2- & >C< [r]",),
                ("-CH2- & =C< [r]",),
                ("-CH2- & F-",),
                ("-CH2- & Cl-",),
                ("-CH2- & Br-",),
                ("-CH2- & I-",),           # 40
                ("-CH2- & -OH",),
                ("-CH2- & -O-",),
                ("-CH2- & >CO",),
                ("-CH2- & -CHO",),
                ("-CH2- & -COOH",),
                ("-CH2- & -COO[-]",),
                ("-CH2- & [-]COO-",),
                ("-CH2- & -NH2",),
                ("-CH2- & -NH-",),
                ("-CH2- & >N-",),          # 50
                ("-CH2- & -CN",),
                ("-CH2- & -SH",),
                ("-CH2- & -S-",),
                (">CH- & >CH-",),
                (">CH- & >C<",),
                (">CH- & =CH-",),
                (">CH- & =C<",),
                (">CH- & >CH- [r]",),
                (">CH- & =C< [r]",),
                (">CH- & F-",),            # 60
                (">CH- & Cl-",),
                (">CH- & -OH",),
                (">CH- & -O-",),
                (">CH- & >CO",),
                (">CH- & -CHO",),
                (">CH- & [-]COO-",),
                (">CH- & -COOH",),
                (">CH- & -NH2",),
                (">CH- & -NH-",),
                (">C< & >C<",),            # 70
                (">C< & =CH-",),
                (">C< & =C<",),
                (">C< & >C< [r]",),
                (">C< & >CH- [r]",),
                (">C< & =C< [r]",),
                (">C< & F-",),
                (">C< & Cl-",),
                (">C< & Br-",),
                (">C< & -OH",),
                (">C< & -O-",),            # 80
                (">C< & >CO",),
                (">C< & -COOH",),
                ("[=]CH2 & [=]CH2",),
                ("[=]CH2 & [=]CH-",),
                ("[=]CH2 & [=]C<",),
                ("[=]CH2 & =C[=]",),
                ("[=]CH- & [=]CH-",),
                ("[=]CH- & [=]C<",),
                ("[=]CH- & =C[=]",),
                ("=CH- & =CH-",),          # 90
                ("=CH- & =C<",),
                ("=CH- & ≡C-",),
                ("=CH- & =C< [r]",),
                ("=CH- & F-",),
                ("=CH- & Cl-",),
                ("=CH- & -O-",),
                ("=CH- & -CHO",),
                ("=CH- & -COOH",),
                ("=CH- & -COO[-]",),
                ("=CH- & [-]COO-",),       # 100
                ("=CH- & -CN",),
                ("[=]C< & [=]C<",),
                ("[=]C< & =C[=]",),
                ("=C< & =C< [r]",),
                ("=C< & F-",),
                ("=C< & Cl-",),
                ("=C[=] & O[=]",),
                ("CH[≡] & CH[≡]",),
                ("CH[≡] & -C[≡]",),
                ("-C[≡] & -C[≡]",),             # 110
                ("-CH2- [r] & -CH2- [r]",),
                ("-CH2- [r] & >CH- [r]",),
                ("-CH2- [r] & >C< [r]",),
                ("-CH2- [r] & =CH- [r]",),
                ("-CH2- [r] & =C< [r]",),
                ("-CH2- [r] & -O- [r]",),
                ("-CH2- [r] & >CO [r]",),
                ("-CH2- [r] & -NH- [r]",),
                ("-CH2- [r] & -S- [r]",),
                (">CH- [r] & >CH- [r]",),       # 120
                (">CH- [r] & >C< [r]",),
                (">CH- [r] & >CH- [r]",),
                (">CH- [r] & [=]C< [r]",),
                (">CH- [r] & -O- [r]",),
                (">CH- [r] & -OH",),
                (">C< [r] & >C< [r]",),
                (">C< [r] & =C< [r]",),
                (">C< [r] & F-",),
                ("[=]CH- [r] & [=]CH- [r]",),
                ("[=]CH- [r] & [=]C< [r]",),    # 130
                ("[=]CH- [r] & [=]N- [r]",),
                ("=CH- [r] & =CH- [r]",),
                ("=CH- [r] & =C< [r]",),
                ("=CH- [r] & -O- [r]",),
                ("=CH- [r] & -NH- [r]",),
                ("=CH- [r] & =N- [r]",),
                ("=CH- [r] & -S- [r]",),
                ("[=]C< [r] & [=]C< [r]",),
                ("[=]C< [r] & -N[=] [r]",),
                ("=C< [r] & =C< [r]",),         # 140
                ("=C< [r] & =C< [r]",),
                ("=C< [r] & -O- [r]",),
                ("=C< [r] & =N- [r]",),
                ("=C< [r] & F-",),
                ("=C< [r] & Cl-",),
                ("=C< [r] & Br-",),
                ("=C< [r] & I-",),
                ("=C< [r] & -OH",),
                ("=C< [r] & -O-",),
                ("=C< [r] & >CO",),             # 150
                ("=C< [r] & -CHO",),
                ("=C< [r] & -COOH",),
                ("=C< [r] & [-]COO-",),
                ("=C< [r] & -NH2",),
                ("=C< [r] & -NH-",),
                ("=C< [r] & >N-",),
                ("=C< [r] & -CN",),
                ("Cl- & >CO",),
                ("[-]COO- & [-]COO-",),
                ("-O- [r] & =N- [r]",),         # 160
                (">CO & -O-",),
                ("-H & -CHO",),
                ("-H & -COOH",),
                ("-H & [-]COO-",),
                ("-NH- & -NH2",),
                ("-S- & -S-",)]}

    FirstOrder = 29

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        group, rest = self._decomposition()
        self.group = group

        resto = 0
        for grp, r in rest.items():
            resto += r

        if resto:
            self.msg = translate(
                "newComponent",
                "Bad definition, check input group and contribution")
            self.status = 0
        else:
            self.msg = ""
            self.status = 1
            return True

    def _decomposition(self):
        """Specific procedure to calculate the molecular weight of compound
        from group contribution"""
        group = []
        rest = {}
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            A, B = self.__coeff__["txt"][i][0].split(" & ")
            for x in range(c):
                for grp in (A, B):
                    # Remove the optional corchetes in multiple link definition
                    # TODO: Add support for structure check for group with
                    # several different connection, in that case the link used
                    # are the corchetes signaled
                    if "[-]" in grp:
                        grp = grp.replace("[-]", "-")
                    if "[=]" in grp:
                        grp = grp.replace("[=]", "=")
                    if "[≡]" in grp:
                        grp = grp.replace("[≡]", "≡")

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
        return group, rest

    def calculo(self):
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        self.Na = self._atoms()

        tc = Pc = vc = tb = mu = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb += c*self.__coeff__["tb"][i]
            tc += c*self.__coeff__["tc"][i]
            Pc += c*self.__coeff__["Pc"][i]
            vc += c*self.__coeff__["vc"][i]
            mu += c*self.__coeff__["mu"][i]

        # Table 1 equations in [6]_
        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb = unidades.Temperature(self.M**-0.404*tb+156)
        self.Tc = unidades.Temperature(self.Tb/(0.5851-0.9286*tc-tc**2))
        self.Pc = unidades.Pressure((0.1285-0.0059*self.Na-Pc)**-2, "bar")
        self.Vc = unidades.SpecificVolume((25.1+vc)/self.M, "ccg")

        # Eq 1 in [7]_
        self.mu = unidades.Viscosity(self.M**1.279*exp(mu-7.6425), "muPas")

        GroupContribution.calculo(self)
