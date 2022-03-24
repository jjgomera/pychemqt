#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy import roots
from scipy.constants import R

from lib.EoS.cubic import Cubic


# Table 1
dat = {
    98: (0.3069, 0.5542, 0.1182, 24.58),
    105: (0.2796, 0.4462, 0.7858, 48.74),
    208: (0.3118, 0.8276, 0.3548, 17.99),
    951: (0.3014, 0.9671, 0.4684, 34.2),
    212: (0.3311, 0.4566, -0.0954, 12.69),
    46: (0.3166, 0.8979, 0.1504, 18.61),
    110: (0.2916, 0.8304, -0.1884, 31.06),
    107: (0.3121, 0.4599, -0.006, 14.54),
    47: (0.3103, 0.7317, 0.2532, 20.19),
    51: (0.2872, 1.12, 0.0178, 33.46),
    111: (0.2781, 1.313, 0.7164, 48.17),
    104: (0.2692, 0.542, 0.5979, 32.16),
    1: (0.3237, -0.0152, -0.1086, 19.64),
    62: (0.2496, 1.1172, 0.7625, 14.13),
    50: (0.3003, 0.7915, 0.2366, 30.89),
    63: (0.2705, 1.1351, 0.8093, 14.08),
    # 962: (0.2824, 0.3701, -0.0181, 51),
    215: (0.2925, 0.9334, 0.8316, 60.89),
    216: (0.2914, 0.9599, 0.6748, 70.59),
    101: (0.2952, 1.0254, 0.4435, 63.41),
    217: (0.2795, 0.8145, 1.0604, 86.66),
    100: (0.2802, 0.8039, 0.5362, 105.38),
    218: (0.2873, 0.8498, 0.8453, 51.55),
    48: (0.3154, 0.8881, 0.1842, 22.48),
    49: (0.2897, 1.0099, 0.2028, 31.19),
    220: (0.2813, 0.9525, 1.1588, 54.98),
    642: (0.2858, 0.9759, 1.0158, 65.51),
    115: (0.288, 0.9167, 0.0415, 35.85),
    2: (0.3111, 0.7197, 0.1406, 25.2),
    117: (0.2317, 1.2548, 1.4833, 44.19),
    227: (0.2818, 0.7023, 0.9248, 59.58),
    65: (0.287, 0.91, 0.4829, 39.1),
    125: (0.2129, 0.764, 0.77, 30.71),
    22: (0.3057, 0.8972, 0.5696, 34.28),
    130: (0.2448, 1.2171, 0.6035, 57.09),
    131: (0.2636, 0.8285, 1.5152, 65.4),
    132: (0.2781, 0.8472, 0.4738, 66.7),
    3: (0.3093, 0.9918, 0.4546, 34.86),
    133: (0.2759, 0.7871, 1.1618, 65.48),
    134: (0.2646, 1.8432, 1.8157, 59.67),
    137: (0.2843, 0.9228, 0.3603, 68.8),
    136: (0.2738, 0.7885, 0.2713, 67.68),
    66: (0.2875, 0.97, -0.0147, 50.08),
    258: (0.2935, 0.8559, 0.0702, 45.87),
    23: (0.2895, 0.7301, -0.2315, 66.48),
    140: (0.2528, 0.8495, 0.8263, 85.96),
    141: (0.2655, 0.9159, 1.5263, 86.97),
    142: (0.2623, 0.9893, 1.3609, 87.46),
    4: (0.3045, 1.1084, 0.6398, 48.77),
    146: (0.2559, 1.6081, 4.1494, 86.84),
    24: (0.2939, 0.9444, 0.4361, 82.55),
    155: (0.2624, 1.1093, 1.5612, 105.67),
    156: (0.2640, 1.0840, 1.3797, 104.77),
    157: (0.2687, 1.0333, 1.1352, 106.76),
    6: (0.2985, 1.1527, 0.288, 69.81),
    5: (0.3032, 1.1856, 0.1648, 63.82),
    162: (0.2707, 0.9429, 1.2132, 106.34),
    290: (0.281, 1.1586, 1.2374, 99.25),
    294: (0.2778, 1.1409, 0.8178, 99.1),
    166: (0.2585, 1.133, 1.5179, 129.6),
    309: (0.2622, 1.1438, 1.424, 131.81),
    310: (0.2628, 1.0929, 1.3114, 133.76),
    311: (0.2651, 1.0882, 1.3618, 131.14),
    8: (0.2931, 1.2328, 0.2250, 88.87),
    7: (0.2798, 0.8820, 0.4085, 118.06),
    318: (0.273, 0.974, 1.1706, 138.56),
    171: (0.2632, 0.7355, 1.269, 132.69),
    172: (0.2688, 0.8204, 1.1676, 123.49),
    322: (0.2721, 0.8713, 0.9878, 105.38),
    173: (0.2669, 0.8009, 1.1559, 137.62),
    40: (0.2923, 1.1336, 0.8214, 70.01),
    38: (0.2861, 0.9236, 0.8187, 114.62),
    10: (0.2904, 1.3817, 0.8945, 98.82),
    11: (0.2768, 1.2947, 0.8374, 136.39),
    12: (0.2739, 1.353, 1.1538, 164.39),
    13: (0.2547, 1.2763, 1.7407, 184.27),
    185: (0.2624, 0.9845, 1.2983, 136.67),
    14: (0.2528, 1.3683, 1.9376, 203.6),
    15: (0.2488, 1.4286, 2.1253, 222.93),
    16: (0.2461, 1.5018, 2.3166, 242.27),
    17: (0.2425, 1.5602, 2.5753, 261.6),
    18: (0.2384, 1.6105, 2.7125, 280.93),
    19: (0.2257, 1.5396, 2.9634, 300.27),
    20: (0.2193, 1.5503, 3.1788, 319.6),
    21: (0.2305, 1.7846, 3.1922, 338.93),
    90: (0.224, 1.7718, 3.3393, 358.27),
    91: (0.2063, 1.6075, 3.651, 377.6),
    92: (0.2069, 1.7712, 4.4157, 396.93)}


class TBS(Cubic):
    r"""Salim modification of Trebble-Bishnoi cubic equation of state

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V^2+(b+c)V+bc+d^2}\\
        a = a_c \alpha\\
        \alpha = 1 + m\left(1-\sqrt{T_r}\right) +
        p\left(\sqrt{0.7}-\sqrt{T_r}\right)\left(1-\sqrt{T_r}\right)\\
        c = \frac{RT_c}{P_c}\left(1-3\zeta_c\right)\\
        \end{array}

    d, ζc, m and p are compound specific parameter available from [1]_ for
    several compounds, for other compound there are generalization correlation

    .. math::
        \begin{array}[t]{l}
        d = v_c/3\\
        \zeta_c = 1.063Z_c\\
        m = 0.662 + 3.12\omega - 0.854\omega^2 + 9.3\left(Z_c-0.3\right)\\
        p = 0.475 + 2.0\omega \; for M <= 128\\
        p = 0.613 + 0.62\omega + 4.06\omega^2 \; for M > 128\\
        \end{array}

    This EoS remove the temperature dependence of co-volume parameter, b, to
    avoid negative values in isochoric heat capacity.
    """

    __title__ = "Trebble-Bishnoi-Salim (1991)"
    __status__ = "TBS"
    __doi__ = {
        "autor": "Salim, P.H., Trebble, M.A.",
        "title": "A Modified Trebble-Bishnoi Equation of State: Thermodynamic "
                 "Consistency Revisited",
        "ref": "Fluid Phase Equilibria 65 (1991) 59-71",
        "doi": "10.1016/0378-3812(91)87017-4"},

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""

        ai = []
        bi = []
        ci = []
        di = []
        for cmp in self.componente:
            a, b, c, d = self.__lib(cmp, T)

            ai.append(a)
            bi.append(b)
            ci.append(c)
            di.append(d)

        self.ai = ai
        self.bi = bi
        self.ci = ci
        self.di = di

    def _GEOS(self, xi):
        coef = [self.ai, self.bi, self.ci, self.di]
        am, bm, cm, dm = self._mixture(None, xi, coef)

        delta = bm+cm
        epsilon = -bm*cm - dm**2

        return am, bm, delta, epsilon

    def __lib(self, cmp, T):
        Tr = T/cmp.Tc

        if cmp.id in dat:
            Xc, m, p, d = dat[cmp.id]
        else:
            Zc = cmp.Pc.kPa*cmp.Vc*cmp.M/R/cmp.Tc

            d = cmp.Vc.ccg*cmp.M/3                                     # Eq 11
            Xc = 1.063*Zc                                              # Eq 12

            # Eq 13
            m = 0.662 + 3.12*cmp.f_acent - 0.854*cmp.f_acent**2 + 9.3*(Zc-0.3)

            if cmp.M <= 128:
                p = 0.475 + 2*cmp.f_acent                              # Eq 14a
            else:
                p = 0.613 + 0.62*cmp.f_acent + 4.06*cmp.f_acent**2     # Eq 14b

        alfa = 1 + m*(1-Tr**0.5) + p*(0.7**0.5-Tr**0.5)*(1-Tr**0.5)

        # From here common functionality of original TB EoS
        # Convert d from cc/mol to m3/mol
        d /= 1e6

        Cc = 1-3*Xc                                                    # Eq 6
        Dc = d*cmp.Pc/R/cmp.Tc                                         # Eq 12

        # Bc calculated ad the smallest positive root of Eq 8
        B = roots([1, 2-3*Xc, 3*Xc**2, -Dc**2-Xc**3])
        Bpositivos = []
        for Bi in B:
            if Bi > 0:
                Bpositivos.append(Bi)
        Bc = min(Bpositivos).real

        Ac = 3*Xc**2 + 2*Bc*Cc + Bc + Cc + Bc**2 + Dc**2                # Eq 7
        ac = Ac*R**2*cmp.Tc**2/cmp.Pc
        a = ac*alfa                                                     # Eq 23

        b = Bc*R*cmp.Tc/cmp.Pc
        c = R*cmp.Tc/cmp.Pc*(1-3*Xc)                                    # Eq 10

        return a, b, c, d


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = TBS(300, 9.9742e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
    eq = TBS(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vg.ccmol))
