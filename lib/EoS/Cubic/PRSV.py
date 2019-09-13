#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from scipy.constants import R

from lib.EoS.cubic import Cubic


# K1 values from Table 1 from [1]_
# K2 and K3 values from Table 2 from reference for PRSV2
dat = {
    46: (0.01996, 0.3162, 0.535),
    47: (0.01512, -0.0090, 0.490),
    49: (0.04285, 0.0000, 0.000),
    63: (0.00100, -0.1265, 0.510),
    62: (-0.06635, 0.0199, 0.443),
    104: (0.01989, -0.0036, 0.310),

    2: (-0.00159, 0.1521, 0.517),
    3: (0.02669, 0.1358, 0.424),
    22: (0.04191, 0.2610, 0.424),
    4: (0.03136, 0.2757, 0.447),
    6: (0.03443, 0.6767, 0.461),
    8: (0.03946, 0.3940, 0.457),
    9: (0.04303, 0.8697, 0.615),
    10: (0.05104, 0.8634, 0.460),
    11: (0.04648, 0.9331, 0.496),
    12: (0.04464, 0.6214, 0.509),
    13: (0.04104, 0.6621, 0.519),
    14: (0.04510, 0.8549, 0.527),
    15: (0.02919, 1.3288, 0.568),
    16: (0.05426, 0.8744, 0.505),
    17: (0.04157, 0.9387, 0.528),
    18: (0.02686, 0.9408, 0.528),
    19: (0.01892, 1.0908, 0.559),
    20: (0.02665, 0.0334, 0.767),
    21: (0.04048, 2.9805, 0.571),
    90: (0.08291, 4.1441, 0.577),
    38: (0.07023, 0.6146, 0.530),
    # 888 : (0.01805, 3.0438, 0.606),

    40: (0.07019, 0.7939, 0.523),
    41: (0.03849, 0.5261, 0.510),
    45: (0.03994, 0.5342, 0.519),
    44: (0.01277, 0.5963, 0.524),
    182: (0.01173, 0.9246, 0.548),
    70: (0.02715, 0.7310, 0.530),
    75: (-0.01384, 0.4777, 0.538),
    185: (0.03297, 0.6634, 0.510),
    191: (-0.01842, -0.8140, 0.577),
    192: (-0.01639, -0.0750, 0.597),
    194: (0.11487, 0.1077, 0.407),
    406: (0.05955, 0.3703, 0.587),
    # : (-0.01393, -5.8256, 0.587),  # 9,10-Dihydrophenanthrene

    140: (-0.00888, 0.2871, 0.537),
    448: (0.00554, 0.6847, 0.613),
    304: (0.01681, 0.7787, 0.620),
    165: (0.03558, 0.6180, 0.610),
    305: (0.04113, 0.5138, 0.615),
    771: (0.00984, 0.8448, 0.550),
    # 1356: (0.02321, 0.8891, 0.619),
    # 1354: (0.04005, 0.7451, 0.535),
    # 809: (0.02731, 0.9236, 0.561),
    610: (0.02002, 0.9649, 0.579),

    117: (-0.16816, -1.3400, 0.588),
    134: (-0.03374, -2.6846, 0.592),
    146: (0.21419, -3.6816, 0.640),
    145: (0.23264, -3.5578, 0.652),
    160: (0.33431, -1.1748, 0.642),
    159: (0.39045, 0.0026, 0.676),
    314: (0.37200, -1.2792, 0.642),
    316: (0.43099, -0.0480, 0.658),
    313: (0.36781, 0.2918, 0.621),
    335: (-0.00237, -3.3938, 0.551),
    360: (0.82940, 2.7372, 0.543),
    396: (0.80898, 1.4978, 0.470),

    133: (0.05717, -0.1211, 0.481),
    486: (0.16948, 0.0515, 0.768),
    # 1287: 0.02300, 0.9179, 0.558),
    # 1286: 0.04123, 0.3833, 0.562),
    # 1317: 0.01622, 0.6140, 0.548),
    456: (0.05129, -0.2022, 0.585),
    318: (-0.01668, 1.1679, 0.553),
    778: (-0.03162, 1.4094, 0.577),
    337: (0.03751, 0.8810, 0.590),
    344: (0.01610, 1.0478, 0.616),

    226: (-0.10299, 0.5905, 0.463),
    125: (-0.13991, -0.3777, 0.522),
    130: (-0.19724, 0.8136, 0.541),
    144: (0.18999, -0.1857, 0.470),
    688: (-0.42503, 0.0020, 0.407),
    270: (0.14326, 1.1002, 0.614),
    271: (0.06001, 2.2851, 0.603),
    # 1285: (-0.09508, 0.0033, 0.931),
    # : (0.13440, 0.9308, 0.635),     # 2-Methyl-2-Propylamine
    281: (0.03961, 0.5029, 0.573),
    295: (0.06946, 0.4241, 0.575),
    164: (-0.03471, 4.5512, 0.539),
    497: (0.11367, -0.7871, 0.537),
    319: (0.02752, 0.8172, 0.565),
    796: (-0.00901, 0.5425, 0.582),
    346: (0.24705, 5.6910, 0.597),
    # : (0.06043, -0.1210, 0.592),   # Thianaphthene

    # Aditional parameter from Table 1 in [2]_, only for PRSV equation
    # Inorganic
    51: (0.03962, 0, 0),
    48: (0.04279, 0, 0),
    50: (0.03160, 0, 0),
    105: (0.03641, 0, 0),
    110: (0.06325, 0, 0),

    # Hydrocarbons
    65: (0.09919, 0, 0),
    5: (-0.00238, 0, 0),
    24: (0.02222, 0, 0),
    68: (-0.00516, 0, 0),
    25: (0.04457, 0, 0),
    26: (0.00435, 0, 0),
    28: (0.02327, 0, 0),
    7: (0.02840, 0, 0),
    29: (0.03521, 0, 0),
    82: (0.03363, 0, 0),
    36: (0.04602, 0, 0),
    39: (0.04859, 0, 0),
    42: (0.02291, 0, 0),
    43: (0.01610, 0, 0),
    71: (0.05509, 0, 0),
    178: (0.16109, 0, 0),

    # Acids
    143: (-0.03093, 0, 0),
    154: (-0.14012, 0, 0),
    306: (-0.08171, 0, 0),

    # Aldehydes
    261: (-0.03642, 0, 0),

    # Esters
    141: (0.11372, 0, 0),
    142: (0.05791, 0, 0),
    155: (0.06464, 0, 0),
    166: (0.06531, 0, 0),
    330: (0.07587, 0, 0),
    309: (0.06362, 0, 0),

    # Ethers
    162: (0.05004, 0, 0),

    # Ketones
    327: (-0.29473, 0, 0),

    # Halogen Compounds
    100: (0.05212, 0, 0),
    218: (0.02136, 0, 0),
    636: (0.04014, 0, 0),
    215: (0.05054, 0, 0),
    216: (0.04722, 0, 0),
    217: (0.03708, 0, 0),
    112: (0.02899, 0, 0),
    642: (0.07915, 0, 0),
    220: (0.04513, 0, 0),
    643: (0.00535, 0, 0),
    645: (-0.02874, 0, 0),
    115: (0.01160, 0, 0),
    225: (-0.00374, 0, 0),
    229: (0.03220, 0, 0),
    231: (0.06095, 0, 0),
    232: (0.05596, 0, 0),
    126: (-0.00984, 0, 0),
    127: (-0.04659, 0, 0),
    132: (0.03982, 0, 0),
    172: (0.03123, 0, 0),

    # Nitrogen Compounds
    118: (0.10475, 0, 0),
    138: (0.03536, 0, 0),
    249: (0.20709, 0, 0),
    147: (0.06414, 0, 0),
    294: (0.11210, 0, 0),
    339: (0.02951, 0, 0),
    175: (0.00316, 0, 0),
    353: (0.03690, 0, 0),
    354: (0.08293, 0, 0),

    # Oxygen Compounds
    129: (0.01757, 0, 0),
    444: (0.03693, 0, 0),
    174: (0.02100, 0, 0),
    177: (-0.00849, 0, 0),
    347: (-0.04098, 0, 0),

    # Sulful Compounds
    149: (0.07654, 0, 0)}


class PRSV(Cubic):
    r"""Peng-Robinson cubic equation of state with a modified dependence of
    temperature by Stryjek-Vera [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + k\left(1-Tr^{0.5}\right)\\
        k = k_0+k_1\left(1+\sqrt{T_r}\right)\left(0.7-T_r\right)\\
        k_0 = 0.378893+1.4897153\omega-0.17131848\omega^2+0.0196554\omega^3\\
        \end{array}

    :math:`k_1` is a parameter characteristic to each compound given in [1]_
    and [2]_

    Examples
    --------
    Example 4.3 from [3]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PRSV(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2040 86.8'
    >>> eq = PRSV(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '84.2'
    """

    __title__ = "PRSV (1986)"
    __status__ = "PRSV"
    __doi__ = (
     {
        "autor": "Stryjek, R., Vera, J.H.",
        "title": "PRSV: An Improved Peng—Robinson Equation of State for Pure "
                 "Compounds and Mixtures",
        "ref": "Can. J. Chem. Eng. 64 (1986) 323-333",
        "doi": "10.1002/cjce.5450640224"},
     {
        "autor": "Proust, P., Vera, J.H.",
        "title": "PRSV: The Stryjek-Vera Modification of the Peng-Robinson"
                 "Equation of Stte. Parameters for Other Pure Compounds of "
                 "Industrial Interest",
        "ref": "Can. J. Chem. Eng. 67 (1989) 170-173",
        "doi": "10.1002/cjce.5450670125"},
     {
        "autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
        "title": "The Properties of Gases and Liquids 5th Edition",
        "ref": "McGraw-Hill, New York, 2001",
        "doi": ""})

    def __init__(self, T, P, mezcla):

        self.T = T
        self.P = P
        self.mezcla = mezcla

        ai = []
        bi = []
        for cmp in mezcla.componente:
            Tr = T/cmp.Tc

            k = self._k(cmp, Tr)

            alfa = (1 + k*(1-Tr**0.5))**2
            ao = 0.457235*R**2*cmp.Tc**2/cmp.Pc
            b = 0.077796*R*cmp.Tc/cmp.Pc

            ai.append(ao*alfa)
            bi.append(b)

        am, bm = self._mixture(None, mezcla.ids, [ai, bi])

        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = 2*bm
        self.epsilon = -bm**2

        super(PRSV, self).__init__(T, P, mezcla)

    def _k(self, cmp, Tr):
        # Eq 11
        ko = 0.378893 + 1.4897153*cmp.f_acent - \
            0.17131848*cmp.f_acent**2 + 0.0196554*cmp.f_acent**3

        if cmp.id in dat and Tr < 1:
            k1 = dat[cmp.id][0]
        else:
            k1 = 0
        k = ko + k1*(1+Tr**0.5)*(0.7-Tr)
        return k


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRSV(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRSV(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
