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


from itertools import permutations
from math import exp

from lib import unidades
from lib.newComponent._base import GroupContribution


class Nannoolal(GroupContribution):
    """
    Group contribution for definition of unknown component using the Nannoolal
    procedure (2007). This method can calculate critical properties, boiling
    temperature, and vapor presssure.

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
    Table 17a in [1]_, 3,3,4,4-tetramethylhexane

    >>> cmp = Nannoolal(group=[0, 3, 5, 131], contribution=[6, 2, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '429.5'
    >>> cmp.formula
    'C10H22'

    Table 17b in [1]_, di-isopropanolamine, find too in [10]_

    >>> cmp = Nannoolal(group=[0, 6, 33, 41], contribution=[2, 4, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '509.3'

    Table 17c in [1]_, perfluoro-2-propanone

    >>> cmp = Nannoolal(group=[6, 20, 50, 118, 119, 121],
    ... contribution=[2, 6, 1, 1, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '246.3'

    Table 17d in [1]_, methyl m-toluate

    >>> cmp = Nannoolal(group=[1, 2, 14, 15, 44, 126, 132],
    ... contribution=[1, 1, 4, 2, 1, 1, 1])
    >>> "%0.1f" % cmp.Tb
    '490.7'

    Table 42a in [2]_, 2,2,3,3-tetramethylbutane

    >>> cmp = Nannoolal(Tb=379.6, group=[0, 5, 131], contribution=[6, 2, 1])
    >>> "%0.1f" % cmp.Tc
    '566.8'

    Table 42b in [2]_, diethylene glycol monomethyl ether

    >>> cmp = Nannoolal(M=120.15, group=[1, 6, 34, 37],
    ... contribution=[1, 4, 1, 2])
    >>> "%0.1f" % cmp.Pc.kPa
    '3689.9'

    Table 42c in [2]_, trichloro silane

    >>> cmp = Nannoolal(group=[26, 111, 122, 133], contribution=[3, 1, 1, 1])
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '262.8'

    Table 42d in [2]_, perfluoro-2-propanone

    >>> cmp = Nannoolal(Tb=245.9, group=[6, 20, 50, 118, 119, 121],
    ... contribution=[2, 6, 1, 1, 2, 1])
    >>> "%0.1f" % cmp.Tc
    '358.0'

    Table 15c in [3]_

    >>> "%0.7f" % cmp.db
    '0.1179392'
    >>> "%0.2f" % cmp._Pv(210.16).kPa
    '14.63'

    Table 15a in [3]_, α-pinene

    >>> cmp = Nannoolal(Tb=429, group=[0, 8, 9, 10, 61, 123, 130],
    ... contribution=[3, 2, 2, 1, 1, 1, 2])
    >>> "%0.7f" % cmp.db
    '0.1078596'
    >>> "%0.2f" % cmp._Pv(388.15).kPa
    '31.03'

    Table 15b in [3]_, 1,2-ethanediol

    >>> cmp = Nannoolal(Tb=470.5, group=[6, 35], contribution=[2, 2])
    >>> "%0.7f" % cmp.db
    '1.1310491'
    >>> "%0.2f" % cmp._Pv(410.65).kPa
    '13.05'

    Table 15d in [3]_, acrylic acid

    >>> cmp = Nannoolal(Tb=413.6, group=[43, 60, 132], contribution=[1, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.9163297'
    >>> "%0.2f" % cmp._Pv(344.15).kPa
    '6.52'

    Table 15e in [3]_, glycol monoacetate

    >>> cmp = Nannoolal(Tb=458.65, group=[0, 6, 35, 44],
    ... contribution=[1, 2, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.5111422'
    >>> "%0.2f" % cmp._Pv(352.65).kPa
    '2.24'

    Table 15f in [3]_, dipropyl succinate

    >>> cmp = Nannoolal(group=[0, 3, 6, 44], contribution=[2, 4, 2, 2])
    >>> "%0.7f" % cmp.db
    '0.9878298'

    Table 15g in [3]_, diethanolamine
    The example in paper has a error, the group interaction contribution has
    the sign changed, fixing this the value get is near to experimental value

    >>> cmp = Nannoolal(Tb=541.15, group=[6, 35, 41], contribution=[4, 2, 1])
    >>> "%0.5f" % cmp._Pv(401.13).kPa
    '0.40405'

    Table 15h in [3]_, R122

    >>> cmp = Nannoolal(Tb=344.25, group=[6, 20, 25, 26, 119, 122],
    ... contribution=[2, 2, 2, 1, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.0447374'
    >>> "%0.4f" % cmp._Pv(297.46).kPa
    '17.5093'

    Table 28a in [4]_, N,N-diethylamine

    >>> cmp = Nannoolal(Tb=329, group=[0, 6, 41], contribution=[2, 2, 1])
    >>> "%0.7f %0.2f" % (cmp.dBv, cmp.Tv)
    '4.7713112 210.38'
    >>> "%0.4f" % cmp._Visco(308.15).mPas
    '0.2674'

    Table 28b in [4]_, ethylene glycol monopropyl

    >>> cmp = Nannoolal(Tb=424.5, group=[0, 3, 6, 34, 37],
    ... contribution=[1, 1, 3, 1, 1])
    >>> "%0.7f %0.3f" % (cmp.dBv, cmp.Tv)
    '6.0699241 301.197'
    >>> "%0.6f" % cmp._Visco(318.15).mPas
    '0.921814'

    Table 28c in [4]_, monoethanolamine

    >>> cmp = Nannoolal(Tb=443.45, group=[6, 35, 39], contribution=[2, 1, 1])
    >>> "%0.7f %0.3f" % (cmp.dBv, cmp.Tv)
    '12.3872997 382.275'
    >>> "%0.4f" % cmp._Visco(363.15).mPas
    '2.6135'
    """
    __title__ = "Nannoolal (2007)"

    __doi__ = {
      1:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D., Cordes, W.",
         "title": "Estimation of Pure Component Properties 1. Estimation of "
                  "the Normal Boiling Point of Non-electrolyte Organic "
                  "Compounds Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 226 (2004) 45-63",
         "doi": "10.1016/j.fluid.2004.09.001"},
      2:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 2. Estimation of "
                  "Critical Property Data by Group Contribution",
         "ref": "Fluid Phase Equilib., 252 (2007) 1-27",
         "doi": "10.1016/j.fluid.2006.11.014"},
      3:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 3. Estimation of "
                  "the Vapor Pressure of Non-Electrolyte Organic Compounds "
                  "Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 269 (2008) 117-133",
         "doi": "10.1016/j.fluid.2008.04.020"},
      4:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 4. Estimation of "
                  "the Saturated Liquid Viscosity of Non-Electrolyte Organic "
                  "Compounds Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 281 (2009) 97-119",
         "doi": "10.1016/j.fluid.2009.02.016"},
      5:
        {"autor": "",
         "title": "Perry's Chemical Engineers' Handbook 9th Edition",
         "ref": "McGraw-Hill (2019)",
         "doi": ""}
        }

    __coeff__ = {
        # Be careful, there are several changes in group between Tb paper and
        # the other paper and several index ajust
        #   - 78 : Two different group with equal parameter but different M,
        #     the Si related copy to 110
        #   - 71, 93: Change definition, use the critical paper definition
        #   - Mising 98, use the 214 defined in critical (be careful to delete
        #     in table
        #   - Mising 112, use the 215 defined in critical
        #   - Mising 114, use the 216 defined in critical
        #   - 217 set at end of second order as 135

        # Table 3 & 4 in [1]_
        "tb": [177.3066, 251.8338, 157.9527, 239.4531, 240.6785, 249.5809,
               266.8769, 201.0115, 239.4957, 222.1163, 209.9749, 250.9584,
               291.2291, 244.3581, 235.3462, 315.4128, 348.2779, 367.9649,
               106.5492, 49.2701, 53.1871, 78.7578, 103.5672, -19.5575,
               330.9117, 287.1863, 267.4170, 205.7363, 292.5816, 419.4959,
               377.6775, 556.3944, 349.9409, 390.2446, 443.8712, 488.0819,
               361.4775, 146.4836, 820.7118, 321.1759, 441.4388, 223.0992,
               126.2952, 1080.3139, 636.2020, 642.0427, 1142.6119, 1052.6072,
               1364.5333, 1487.4109, 618.9782, 553.8090, 434.0811, 461.5784,
               864.5074, 304.3321, 719.2462, 475.7958, 586.1413, 500.2434,
               412.6276, 475.9623, 512.2893, 422.2307, 37.1936, 453.3397,
               306.7139, 866.5843, 821.4141, 282.0181, 207.9312, 920.3617,
               1153.1344, 494.2668, 1041.0851, 1251.2675, 778.9151, 540.0895,
               879.7062, 660.4645, 1018.4865, 1559.9840, 510.4223, 1149.9670,
               1209.2972, 347.7717, 664.0903, 957.6388, 928.9954, 560.1024,
               229.2288, 606.1797, 215.3416, 273.1755, 1218.1878, 2082.3288,
               201.3224, 0, 886.7613, 1045.0343, -109.6269, 111.0590,
               1573.3769, 1483.1289, 1506.8136, 484.6371, 1379.4485, 659.7336,
               492.0707, 540.0895, 971.0365, 0, 428.8911, 0, 612.9506,
               562.1791, 761.6006,
               -82.2328, -247.8893, -20.3996, 15.4720,
               -172.4201, -99.8035, -62.3740, -40.0058, -27.2705, -3.5075,
               16.1061, 25.8348, 35.8330, 51.9098, 111.8372, 40.4205, 0],

        # Table 7 & 8 in [2]_
        "tc": [41.8682, 33.1371, -1.0710, 40.0977, 30.2069, -3.8778, 52.8003,
               9.4422, 21.2898, 26.3513, -17.0459, 51.7974, 18.9549, -29.1568,
               16.1154, 68.2045, 68.1923, 29.8039, 15.6068, 11.0757, 18.1302,
               19.1772, 20.8519, -24.0220, -1.3329, 2.6113, 15.5010, -16.1905,
               60.1907, 5.2621, -21.5199, -8.6881, 84.8567, 79.3047, 49.5968,
               130.1320, 14.0159, 12.5082, 41.3490, 18.3404, -50.6419, 17.1780,
               -0.5820, 199.9042, 75.7089, 58.0782, 109.1930, 102.1024, 0, 0,
               56.1572, 44.2000, -7.1070, 0.5887, 0, -7.7181, 117.1330,
               45.1531, 0, 67.9821, 45.4406, 56.4059, -19.9737, 36.0883,
               10.4146, 18.9903, 10.9495, 82.6239, 0, 25.4209, 72.5587, 0, 0,
               0, 0, 164.3355, 0, 157.3401, 97.2830, 153.7225, 0, 90.9726,
               62.3642, 0, 0, 0, 53.6350, 24.7302, 0, 38.4681, 0, 63.6504,
               34.2058, 0, 0, 0, 27.3441, 48.1680, 0, 0, 0, 1.3231, 764.9595,
               0, 0, 0, 0, 0, 0, 157.3401, 0, 0.2842, 0, -0.6536, 36.0361, 0,
               0, 32.1829, 11.4437, -1.3023, -34.3037, -1.3798, -2.7180,
               11.3251, -4.7516, 1.2823, 6.7099, 0, -33.8201, -18.4815,
               -23.6024, -24.5802, -35.6113, 62.0286],

        # Table 10 & 11 in [2]_
        "Pc": [8.1620, 5.5262, 4.1660, 5.2623, 2.3009, -2.9925, 3.4310, 2.3665,
               3.4027, 3.6162, -5.1299, 4.1421, 0.8765, -0.1320, 2.1064,
               4.1826, 3.5500, 1.0997, 0.7328, 4.3757, 3.4933, 2.6558, 1.6547,
               0.5236, -2.2611, -1.4992, 0.4883, -0.9280, 11.8687, -4.3170,
               -2.2409, -4.7841, -7.4244, -4.4735, -1.8153, -6.8991, -12.1664,
               2.0592, 0.1759, -4.4164, -9.0065, -0.4086, 2.3625, 3.9873,
               4.3592, 1.0266, 0.4329, 0.5172, 0, 0, 0.1190, -2.3615, -9.4154,
               -8.2595, 0, -4.9259, 5.1666, 7.1581, 0, -6.2791, 9.6413, 3.4731,
               -2.2718, 2.4489, -0.5403, 8.3052, -4.7101, -5.0929, 0, 5.7270,
               2.7602, 0, 0, 0, 0, 4.0458, 0, 12.6786, 0.2822, 0, 0, -23.9221,
               0.7043, 0, 0, 0, 12.6128, -10.2451, 0, -4.0133, 0, -5.0403,
               3.2023, 0, 0, 0, -4.3834, 1.5574, 0, 0, 0, 3.3971, 58.9190, 0,
               0, 0, 0, 0, 0, 12.6786, 0, 3.8751, 0, 4.4882, -5.1116, 0, 0,
               7.3149, 4.1439, 0.4387, -4.2678, 4.8944, 2.8103, -0.3035,
               0.0930, 0.7061, -0.7246, 0, -8.8457, -2.2542, -3.2460, -5.3113,
               1.0934, 8.6126],

        # Table 13 & 14 in [2]_
        "vc": [28.7855, 28.8811, 26.7237, 32.0493, 32.1108, 28.0534, 33.7577,
               28.8792, 24.8517, 30.9323, 5.9550, 29.5901, 20.2325, 10.5669,
               19.4020, 25.0434, 5.6704, 16.4118, -5.0331, 1.5646, 3.3646,
               1.0897, 1.1084, 19.3190, 22.0457, 23.9279, 26.2582, 36.7624,
               34.4110, 36.0223, 30.7004, 48.2989, 10.6790, 5.6645, 2.0869,
               3.7778, 25.6584, 11.6284, 46.7680, 13.2571, 73.7444, 20.5722,
               6.0178, 40.3909, 42.6733, 36.1286, 0, 64.3506, 0, 0, 30.9229,
               25.5034, 34.7699, 38.0185, 0, 20.3127, 43.7983, 0, 0, 51.0710,
               48.1957, 34.1240, 40.9263, 29.8612, 4.7476, -25.3680, 23.6094,
               34.8472, 0, 75.7193, 69.5645, 0, 0, 0, 0, 0, 0, 0, 52.8789,
               27.1026, 0, 68.0701, 0, 0, 0, 0, 0, 64.4616, 0, 20.0440, 0,
               28.7127, 55.3822, 0, 0, 0, 29.3068, 16.3122, 0, 0, 0, 1.3597, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 37.0423, 0, 55.7432, 16.2688, 0, 0,
               -3.8033, 27.5326, 1.5807, -2.6235, -5.3091, -6.1909, 3.2219,
               -6.3900, -3.5964, 1.5196, 0, -4.6483, -5.0563, -6.3267, 4.9392,
               2.8889, 19.4348],

        # Table 4 & 5 in [3]_
        "Pv": [13.3063, 91.8000, 50.1939, 54.6564, 45.7437, -31.7531, 37.8485,
               96.1386, 22.2573, 32.8162, 4.8500, 23.6411, 49.8237, -3.6950,
               32.7177, 69.8796, 41.5534, 43.7191, 79.5429, 51.2880, 42.0887,
               56.9998, 142.1060, 45.9652, 93.6679, 67.8082, 55.9304, 46.0435,
               84.9162, 104.9291, -40.1837, 134.3501, 719.3666, 758.4218,
               700.7226, 756.0824, 441.8437, 108.4964, 286.9731, 251.9212,
               361.7760, 193.7667, -102.7252, 1074.1000, 355.7381, 350.5184,
               292.8046, 269.2471, 736.9540, 1216.0700, 255.8480, 252.9059,
               123.2143, 127.3380, 222.2789, 20.1604, 226.1873, 86.4601,
               224.2238, 134.9382, 34.2541, 97.4210, 206.6665, 128.0247,
               48.8839, 375.0486, 126.3340, 375.8217, 238.2066, 2.8992, 9.3624,
               603.5347, 662.0582, 510.9666, 1317.4360, 681.3525, 564.1116,
               391.3697, 318.2350, 435.8446, 218.6012, 381.5945, 80.2735,
               231.3919, 186.9204, 48.5026, 168.3007, 108.5277, 213.7165,
               183.1130, 1178.1950, 158.3258, -47.3420, 186.7950, 392.2006,
               595.1778, 268.7081, 86.9450, 612.9546, 258.9924, -316.4392,
               64.0566, 660.2247, 554.7941, 420.7591, 0, -237.2124, 0, 0,
               391.3697, 319.4879, -66.5670, 118.8412, -81.1543, 305.1341,
               191.5058, 423.5251, 34.3545, 2.5030, -83.3326, -64.4854,
               -125.9208, -47.2962, 33.9765, -7.0982, -45.0531, -3.2036, 0,
               -20.6706, -36.3170, -1.1994, 123.7433, -15.9694, 36.7574],

        # Table 5 & 6 in [4]_
        "dBv": [13.9133, 11.7002, -11.0660, 2.1727, 4.5878, 37.0296, 21.3473,
                5.9452, 10.8799, -7.2202, 142.1976, 61.0811, 28.7351, -12.3456,
                -2.7840, 45.9403, 80.5124, 37.4124, 5.1640, 0, 2.8323, 0.7129,
                -36.3189, -61.9434, 4.7579, 5.8228, 4.6555, -67.3989, -9.6209,
                0.5164, -18.1984, -17.3110, 336.8834, 365.8067, 249.0118,
                218.8000, 160.8315, -35.3055, 85.3693, 58.9131, 44.0698,
                13.6479, -58.7354, 54.7891, 17.6757, 0.4267, 47.6109, 12.6717,
                129.8293, 202.2864, 24.2524, 18.4961, 30.5022, -0.0276,
                -13.4614, 18.5507, 23.1459, 9.5809, 152.2693, 18.5983, 21.4560,
                19.7836, -165.0071, 13.3585, 42.7958, 151.9493, 52.5900,
                -34.3948, -6.5626, -25.5950, -28.3943, -30.6156, 45.9972,
                -7.3298, 369.8367, 16.3525, -2.6553, -61.2368, -7.5067, 4.1408,
                -46.5613, 122.6902, 0, 0, 0, -70.9713, 0, 29.0985, 125.0861,
                -14.2823, 0, -8.5352, 9.9037, 0, 0, 102.0816, 74.0520, 37.5669,
                0, 54.4769, 0, 5.7765, 95.6531, 56.9133, 64.7133, 0, 22.9969,
                23.2473, 0, -61.2368, 0, 64.6600, 0, 68.4952, 0, 0, 0, 0.3041,
                0, -6.1420, -26.4635, -14.9636, -25.9017, -57.3789, -21.2204,
                -20.1917, -34.5860, 0, -110.7391, 2.4859, -59.3670, 0, 13.1413,
                -76.1631],

        # Table 8 & 9 in [4]_
        "Tv": [89.0803, 216.0226, 80.9698, 60.3316, 24.2637, 244.4643,
               103.4109, -16.5212, 174.1316, -37.7584, 252.0190, 251.9299,
               330.7100, 294.3323, 113.9028, -26.6195, 133.5499, 128.4739,
               208.3258, 0, 35.2688, 207.3562, -15.8544, 112.1172, 329.0113,
               313.1106, 194.6060, -8.6247, 182.7067, 456.3713, 391.6060,
               499.2149, 1199.4010, 1198.1040, 1078.0840, 1284.7450, 1134.1640,
               -34.9892, 612.7222, 458.7425, 705.1250, 159.5146, -284.4707,
               1446.0240, 325.5736, 454.1671, 374.6477, 289.9690, 1150.8290,
               1619.1650, 304.5982, 394.7932, 294.7319, 206.6432, 292.3613,
               302.2321, 346.9998, -23.9801, 238.3242, 137.5408, 74.4489,
               304.9257, -32.4179, 57.8131, 279.2114, 662.0051, 277.5038,
               369.4221, 488.1136, -10.6146, -181.7627, 351.0623, -15.2801,
               174.3672, 1098.1570, 549.1481, 394.5776, 10.3752, 365.8081,
               164.8904, 197.1806, 1297.7560, 0, 0, 0, 35.4672, 0, 495.5141,
               551.9254, 490.7224, 0, 669.0158, 256.5078, 0, 0, 1787.0390,
               220.0803, 192.1303, 0, 131.2253, 0, 53.2507, 288.4140, 542.6641,
               714.0494, 0, 797.2271, 253.5303, 0, 10.3752, 0, 377.7146, 0,
               806.8125, 0, 0, 0, -180.3686, 0, 241.8968, 138.6555, -71.1647,
               -115.0418, -96.7544, -153.8442, -22.1041, 24.7835, 0, 224.2439,
               24.2539, 137.8708, 0, -54.1782, -726.4291],

        # Name and group composition, Table 1 in [1]_
        "txt": [                                           # 0
                ("CH3-", ),
                ("CH3- (N, O, F, Cl)", ),
                ("CH3- (aromatic)", ),
                ("-CH2- ", ),
                (">CH- ", ),
                (">C< ", ),
                (">C< (N, O, F, Cl)", ),
                (">C< (aromatic)", ),
                ("-CH2- (ring)", ),
                (">CH- (ring)", ),                         # 10
                (">C< (ring) ", ),
                (">C< (ring) (outer N, O, Cl, F)", ),
                (">C< (ring N, O)", ),
                (">C< (aromatic)", ),
                ("=CH- (aromatic) ", ),
                ("=C< (aromatic)", ),
                ("=C< (aromatic) (N, O, Cl, F)", ),
                ("=C< (full aromatic)", ),
                ("F– (C, Si)", ),
                ("-CF=C< ", ),                             # 20
                ("F– (CF)", ),
                ("F– (CF, CCl)", ),
                ("F– (CF2, CCl2)", ),
                ("F- (aromatic) ", ),
                ("Cl– (C, Si)", ),
                ("Cl– (CF, CCl)", ),
                ("Cl– (CF2, CCl2)", ),
                ("Cl– (aromatic)", ),
                ("-CCl=C< ", ),
                ("Br-", ),                                 # 30
                ("Br- (aromatic)", ),
                ("I-", ),
                ("-OH (>C< no H)", ),
                ("-OH (-CH2-)", ),
                ("-OH (>C5 chain)", ),
                ("-OH (<C5 chain)", ),
                ("-OH (aromatic)", ),
                ("-O-", ),
                (">OC2<", ),
                ("NH2-", ),                                # 40
                ("NH2- (aromatic)", ),
                ("-NH-", ),
                (">N-", ),
                ("-COOH", ),
                ("-COO-", ),
                ("HCOO-", ),
                ("-COO– (ring)", ),
                ("-CON<", ),
                ("-CONH-", ),
                ("-CONH2", ),                              # 50
                (">C=O", ),
                ("-CHO", ),
                ("-SH", ),
                ("-S-", ),
                ("-S-S-", ),
                ("-S- (aromatic) ", ),
                ("-C≡N", ),
                (">C=C< ", ),
                (">C=C< (aromatic)", ),
                (">C=C< (N, O, Cl, F)", ),                 # 60
                ("H2C=C<", ),
                (">C=C< (ring) ", ),
                ("-C≡C-", ),
                ("HC≡C-", ),
                ("-O- (aromatic)", ),
                ("=N- (ring C5)", ),
                ("=N- (ring C6)", ),
                ("NO2-", ),
                ("NO2- (aromatic)", ),
                (">Si<", ),                                # 70
                (">Si< (3CH)", ),  # The definiton in Tb is >Si< (O)
                ("NO3-", ),
                (">PO4-", ),
                ("O=NO-", ),
                ("ONC-", ),
                ("-C=O-O-C=O-", ),
                ("COCl-", ),
                (">BO3-", ),
                ("COOO<", ),
                ("OCN-", ),                                # 80
                ("SCN-", ),
                ("-SO2-", ),
                (">Sn<", ),
                ("AsCl2-", ),
                ("-GeCl3 ", ),
                (">Ge<", ),
                (">C=C=C<", ),
                (">C=C-C=C<", ),
                (">C=C-C=C<", ),
                ("-CHO (aromatic)", ),                     # 90
                ("=N-", ),
                (">C=O (aromatic)", ),
                (">Si< (2CH)", ),  # The definiton in Tb is >Si< (F, Cl)
                ("-O-O-", ),
                (">C≡C-C≡C<", ),
                ("-C=O-O-C=O- (ring)", ),
                ("-NH- (aromatic)", ),
                ("=C< (aromatic)", ),  # Missing 98, set 214
                ("-OCON<", ),
                (">N-C=O-N<", ),                          # 100
                (">N<", ),
                ("F– (CCl)", ),
                ("-OCOO-", ),
                (">SO4", ),
                ("-SO2N<", ),
                ("=CNC=NC=", ),
                (">S=O", ),
                ("(S) -C≡N", ),
                (">N-C=O", ),
                # Missing 110, use a dupplicate of 78    # 110
                # with different molecular weigh
                (">Si< (F, Cl)", ),
                ("(N) -C≡N", ),
                (">Si< (1CH)", ),  # Missing 112, set 215
                (">P<", ),
                (">Si< (no CH)", ),  # Missing 114, set 216
                ("-ON=", ),
                (">Se<", ),
                (">Al<", ),
                # 118 in Tb equal to 134 in other, delete reference
                # So down to this is 1 index below

                # Second order contributions
                ("CO-CXy (X: F,Cl) (y>2)", ),
                ("CXy-CO-CXy (X: F,Cl) (y>2)", ),         # 120
                ("-CF3, -CCl3", ),
                (">CF2, >CCl2", ),
                ("No hydrogen", ),
                ("One hydrogen", ),
                ("3/4 Ring", ),
                ("Five-ring", ),
                ("Ortho pair(s)", ),
                ("Meta pair(s)", ),
                ("Para pair(s)", ),
                ("(C)(C=)>C<", ),                         # 130
                ("C2C-CC2", ),
                ("C3C-CC2", ),
                ("C3C-CC3", ),
                ("C=C-C=O", ),  # Same as 118 in Tb
                (">Si< (F,Cl,Br,I)", )]}  # 217 set to 135

    # Group interactions definition, Table 9 in [1]_
    # See Table 5 & 6 in [2]_ to check changes
    GI = {
        33: "A", 34: "A", 35: "A",
        36: "B",
        43: "C",
        37: "D",
        38: "E",
        44: "F", 45: "F", 46: "F",
        50: "G", 91: "G",
        51: "H", 89: "H",
        64: "I",
        53: "J",
        55: "K",
        52: "L",
        39: "M", 40: "M",
        41: "N", 96: "N",
        79: "O",
        56: "P",
        68: "Q",
        65: "R",
        66: "S"}

    # Table 5 in [1]_
    GI_Tb = {
        "AA": 291.7985, "AM": 314.6126, "AN": 286.9698, "AL": 38.6974,
        "AC": 146.7286, "AD": 135.3991, "AE": 226.4980, "AF": 211.6814,
        "AG": 46.3754, "AJ": -74.0193, "AP": 306.3979, "AI": 435.0923,
        "AS": 1334.6747, "BB": 288.6155, "BM": 797.4327, "BC": -1477.9671,
        "BD": 130.3742, "BF": -1184.9784, "BG": 0, "BH": 43.9722,
        "BQ": -1048.1236, "BS": -614.3624, "MM": 174.0258, "MN": 510.3473,
        "MD": 124.3549, "MF": 182.6291, "MJ": -562.3061, "MQ": 663.8009,
        "MI": 395.4093, "MS": 27.2735, "NN": 239.8076, "ND": 101.8475,
        "NF": 317.0200, "NG": -215.3532, "NS": 758.9855, "LL": 217.6360,
        "LA": 501.2778, "CC": 117.2044, "CD": 612.8821, "CF": -183.2986,
        "CG": -55.9871, "OO": -356.5017, "OQ": -263.0807, "DD": 91.4997,
        "DE": 178.7845, "DF": 322.5671, "DG": 15.6980, "DH": 17.0400,
        "DJ": 394.5505, "DQ": 963.6518, "DP": 293.5974, "DI": 329.0050,
        "EE": 1006.3880, "EH": 163.5475, "FF": 431.0990, "FG": 22.5208,
        "FQ": -205.6165, "FP": 517.0677, "FI": 707.9404, "GG": -303.9653,
        "GH": -391.3690, "GQ": -3628.9026, "GK": 381.0107, "GP": -574.2230,
        "GI": 176.5481, "GS": 124.1943, "HH": 582.1763, "HQ": 140.9644,
        "HK": 397.5750, "HI": 674.6858, "JJ": -11.9406, "QQ": 65.1432,
        "KP": -101.2319, "KR": -348.7400, "PS": -370.9729, "IR": -888.6123,
        "RR": 0, "SS": -271.9449}

    # Table 9 in [2]_
    GI_Tc = {
        "AA": -434.8568, "AM": 120.9166, "AN": -30.4354, "AD": -146.7881,
        "BB": 144.4697, "MM": -60.9217, "DM": -738.0515, "NN": -49.7641,
        "OO": -1866.0970, "DD": 162.6878, "DE": 707.4116, "DF": 128.2740,
        "DJ": -654.1363, "DP": 741.8565, "FF": 366.2663, "GG": 1605.5640,
        "JJ": -861.1528, "KR": 131.7924, "IR": 24.0243, "SS": -32.3208}

    # Table 12 in [2]_
    GI_Pc = {
        "AA": -5.6023, "AM": 69.8200, "AN": 6.1331, "AD": 7.3373,
        "BB": 57.8350, "MM": -0.6754, "DM": -125.5983, "NN": 22.1871,
        "DD": 2.6751, "DE": 88.8752, "DF": -1.0295, "DJ": 25.8246,
        "FF": 0.5195, "GG": -78.2743, "JJ": 43.9001, "KR": -19.7033,
        "IR": -35.1998, "SS": 12.5371}

    # Table 15 in [2]_
    GI_Vc = {
        "AM": -8.0423, "AD": 19.7707, "BB": 97.5425, "MM": -57.1233,
        "OO": 44.1062, "DD": -23.6366, "DE": -329.5074, "DF": -55.5112,
        "DJ": -37.2468, "FF": -74.8680, "GG": -413.3976, "JJ": -403.1196,
        "KR": 164.2930, "IR": 217.9243, "SS": -26.4556}

    # Table 6 in [3]_
    GI_Pv = {
        "AA": -561.5153, "AM": 1067.6660, "AN": 42.4825, "AD": -799.5332,
        "AE": -618.2760, "AF": -1797.6930, "AG": -1181.5990, "AP": 1431.2430,
        "AI": 1132.0400, "BB": -97.6205, "BD": -751.6676, "BF": 548.4352,
        "MM": 1085.8320, "MN": -206.7811, "MD": -198.2791, "MF": -1676.4770,
        "MQ": 1659.0340, "NN": -307.1018, "ND": 65.4421, "LL": 240.3037,
        "CC": -2601.7090, "CG": -787.8563, "OO": -3929.1300, "DD": 144.6074,
        "DE": 1118.9580, "DF": -225.7802, "DG": 1981.2980, "DH": 362.7540,
        "DJ": -1425.0170, "DP": 743.3353, "EE": -3748.8180, "FF": 920.3138,
        "FG": 1594.1640, "FP": 108.1305, "FI": 1590.3210, "GG": -1270.0830,
        "HH": 946.7309, "HI": 705.3049, "JJ": 838.3372, "QQ": -1501.3550,
        "KR": 675.0414, "PS": 994.4996, "IR": 135.5896, "SS": -29.6785}

    # Table 7 in [4]_
    GI_dBv = {
        "AA": -112.4939, "AM": 1031.5920, "AN": 853.2318, "AD": -423.9834,
        "AP": -683.0189, "AI": -557.5079, "BB": -1186.0500, "BD": -333.5638,
        "BQ": -878.0615, "MM": 135.3183, "MD": 219.9701, "ND": -134.4625,
        "DD": 132.0275, "DF": 44.8702, "DG": -219.5265, "DH": 546.5846,
        "DQ": -59.3635, "FF": 964.0840, "FG": 126.0380, "FP": 539.2401,
        "GG": 3705.4400, "HI": 50.1063, "QQ": 896.3606, "PS": -196.6361}

    # Table 10 in [4]_
    GI_Tv = {
        "AA": -1313.5690, "AM": -41.9608, "AN": -1868.6060, "AD": -643.4378,
        "AP": -345.7844, "AI": 50.2582, "BB": -1146.1070, "BD": -229.2406,
        "BQ": 515.1511, "MM": 86.7249, "MD": -57.1437, "ND": 54.2025,
        "DD": 156.7495, "DF": 273.6616, "DG": -339.6071, "DH": 1050.3190,
        "DQ": 355.0508, "FF": 167.7204, "FG": 244.0583, "FP": 334.4856,
        "GG": 1985.8270, "HI": 161.7447, "QQ": 1839.2630, "PS": 718.1262}

    FirstOrder = 115
    SecondOrder = 134

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        nh = self._atomX("H")
        na = self._atoms()
        n = na-nh  # Total atoms of compound except hydrogen

        tb, tc, pc, vc, pv, dbv, tv = 0, 0, 0, 0, 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb += c*self.__coeff__["tb"][i]
            tc += c*self.__coeff__["tc"][i]*1e-3
            pc += c*self.__coeff__["Pc"][i]*1e-4
            vc += c*self.__coeff__["vc"][i]
            pv += c*self.__coeff__["Pv"][i]*1e-3
            dbv += c*self.__coeff__["dBv"][i]*1e-3
            tv += c*self.__coeff__["Tv"][i]

        # Group interaction calculation
        GI = []
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i in self.GI:
                for x in range(c):
                    GI.append(self.GI[i])
        m = len(GI)

        for pair in permutations(GI, 2):
            key = "".join(sorted(pair))
            tb += self.GI_Tb.get(key, 0)/n/(m-1)
            tc += self.GI_Tc.get(key, 0)*1e-3/n/(m-1)
            pc += self.GI_Pc.get(key, 0)*1e-4/n/(m-1)
            vc += self.GI_Vc.get(key, 0)/n/(m-1)
            pv += self.GI_Pv.get(key, 0)*1e-3/n/(m-1)
            dbv += self.GI_dBv.get(key, 0)*1e-3/n/(m-1)
            tv += 2*self.GI_Tv.get(key, 0)/n/(m-1)

        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb = unidades.Temperature(tb/(n**0.6583+1.6868)+84.3395)
        self.Tc = unidades.Temperature(self.Tb*(0.699+1/(0.9889+tc**0.8607)))
        self.Pc = unidades.Pressure(self.M**-0.14041/(0.00939+pc)**2, "kPa")
        self.Vc = unidades.SpecificVolume((vc/n**-.2266+86.1539)/self.M, "ccg")

        # Eq 7 in [3]_
        self.db = pv - 0.176055

        # Eq 7 in [4]_
        self.dBv = dbv/(n**-2.5635+0.0685)+3.777

        # Eq 8 in [4]_
        Tv = 21.8444*self.Tb**0.5+tv**0.9315/(n**0.6577+4.9259)-231.1361
        self.Tv = unidades.Temperature(Tv)

        GroupContribution.calculo(self)

    def _Pv(self, T):
        """Vapor pressure calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        Pv : float
            Vapor Pressure, [Pa]
        """
        # Eq 6 in [3]_
        Trb = T/self.Tb
        Pv = 10**((4.1012+self.db)*((Trb-1)/(Trb-0.125)))
        return unidades.Pressure(Pv, "atm")

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
        mu = 1.3*exp(-self.dBv*(T-self.Tv)/(T-self.Tv/16))
        return unidades.Viscosity(mu, "cP")
