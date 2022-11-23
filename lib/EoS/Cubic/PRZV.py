#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


from math import log

from lib.EoS.Cubic.PR import PR


# Compound specific parameters given in Table 3
dat = {
    1: (-0.23500, 0.00118),
    46: (-0.26850, 0.00009),
    47: (0.63500, -0.00879),
    107: (-0.06690, -0.00003),
    105: (0.62930, -0.00940),
    98: (-0.65220, 0.00012),
    99: (-0.57610, -0.00024),
    213: (-1.01540, 0.00016),
    63: (0.20980, -0.00027),
    62: (1.07890, -0.00807),
    50: (-0.70300, -0.00004),
    104: (-0.05150, 0.00019),
    106: (-0.42120, 0.00003),
    48: (-0.25670, 0.00155),
    51: (-0.10090, -0.00090),
    2: (-0.56500, 0.00005),
    3: (0.48530, -0.01088),
    4: (0.37320, -0.00114),
    6: (0.74250, -0.01209),
    8: (0.68330, -0.00860),
    10: (0.20780, -0.00924),
    11: (1.11240, -0.01066),
    12: (1.21350, -0.01169),
    13: (1.09110, -0.00938),
    14: (1.37690, -0.01064),
    15: (1.06960, -0.00729),
    16: (1.15280, -0.00670),
    17: (1.40340, -0.00500),
    18: (1.60370, -0.00536),
    19: (1.85870, -0.00313),
    20: (1.53570, -0.00146),
    21: (1.82780, -0.00116),
    90: (2.12220, -0.00113),
    38: (-0.18820, -0.00007),
    9: (-0.88220, 0.00011),
    888: (0.65650, -0.00565),
    22: (0.49930, -0.01076),
    23: (0.30850, -0.00268),
    40: (-0.83570, -0.00017),
    41: (0.78220, -0.01074),
    45: (0.50850, -0.00269),
    70: (0.42120, -0.00075),
    42: (0.82100, -0.01235),
    43: (0.99940, -0.01149),
    44: (1.11130, -0.00692),
    75: (1.22690, -0.02824),
    182: (0.83860, -0.01082),
    178: (1.27390, -0.01432),
    194: (0.49530, -0.01263),
    406: (0.81050, -0.00801),
    185: (0.69880, -0.01145),
    191: (1.62110, -0.01164),
    192: (-0.07110, -0.00887),
    200: (0.50500, -0.00318),
    133: (1.20930, -0.01463),
    486: (-0.00790, -0.01161),
    162: (0.36140, -0.01169),
    # 1287: (1.15620, -0.00930),
    # 1286: (0.80540, -0.03439),
    # 1317: (0.30500, -0.00667),
    456: (0.71820, -0.03442),
    318: (-0.87620, -0.00055),
    778: (0.17190, -0.00068),
    337: (2.43180, -0.01728),
    140: (1.17340, -0.01182),
    153: (1.41690, -0.00763),
    165: (1.36050, -0.01400),
    304: (2.16220, -0.01295),
    771: (1.19260, -0.03419),
    # 1356: (1.02970, -0.01223),
    # 1400: (2.74400, -0.03240),
    610: (-0.08710, -0.04217),
    # 1354: (0.52930, -0.03585),
    117: (2.19680, -0.00505),
    134: (2.35920, -0.00017),
    146: (2.44650, -0.00003),
    160: (2.48680, -0.00255),
    313: (1.01740, -0.00419),
    335: (3.32420, -0.02722),
    360: (-0.68880, -0.08422),
    396: (0.80060, -0.02099),
    417: (2.76770, -0.07755),
    422: (4.90170, -0.04306),
    145: (2.73290, -0.00506),
    159: (2.03360, -0.00159),
    450: (2.02440, -0.00130),
    161: (-2.16510, -0.00025),
    130: (0.87350, -0.00213),
    143: (3.67200, -0.02336),
    270: (1.78050, -0.01903),
    271: (1.62350, -0.01575),
    292: (2.22210, -0.02147),
    718: (2.12430, -0.01938),
    717: (0.95570, -0.03363),
    293: (0.16670, -0.04010),
    294: (1.98150, -0.01714),
    112: (0.83430, -0.00744),
    115: (0.31930, -0.00937),
    126: (0.45640, -0.01008),
    127: (1.44510, -0.01047),
    226: (0.61700, -0.00754),
    796: (2.35030, -0.00489),
    797: (2.02630, -0.01577),
    798: (2.86580, -0.00661),
    344: (0.14750, -0.03986),
    281: (0.60340, -0.01238),
    155: (0.72750, -0.03795),
    164: (1.18820, -0.00466),
    346: (1.59520, -0.01114),
    125: (1.30280, -0.00874),
    295: (0.80770, -0.01428),
    739: (0.71150, -0.03169),
    # 1285: (2.26580, -0.01101),
    144: (-5.70320, -0.06072),
    # 1410: (0.46880, -0.03572),
    319: (-3.92610, -0.00081)}


class PRZV(PR):
    r"""Peng-Robinson cubic equation of state with a modified dependence of
    temperature by Zaboloy-Vera

    .. math::
        \begin{array}[t]{l}
        \alpha = f_1\left(T_{r_{bp}},T_{r_{tp}}\right)\left(
        f_2\left(T_r, T_{r_{bp}},T_{r_{tp}}, \alpha_{tp}, \alpha_{bp},
        \alpha_{\omega}\right) +
        k_1 f_3\left(2, T_r, T_{r_{bp}}, T_{r_{tp}}\right) +
        k_2 f_3\left(11, T_r, T_{r_{bp}}, T_{r_{tp}}\right)\right)\\
        \end{array}

    where k₁ and k₂ are adjustable parameters compound specific
    """

    __title__ = "PR ZV (1996)"
    __status__ = "PRZV"
    __doi__ = {
        "autor": "Zaboloy, M.S., Vera, J.H.",
        "title": "Cubic Equation of State for Pure Compound Vapor Pressure "
                 "from the Triple Point to the Critical Point",
        "ref": "Ind. Eng. Chem. Res. 35(3) (1996) 829-836",
        "doi": "10.1021/ie950306s"},

    def _alfa(self, cmp, T):
        # m,n parameters compound specific
        if cmp.id in dat:
            k1, k2 = dat[cmp.id]

            Tr = T/cmp.Tc
            Trbp = cmp.Tb/cmp.Tc

            # TODO: The correlation use the triple point temperature dont
            # available in database, using normal melting point
            Trtp = cmp.Tf/cmp.Tc

            alfatp = PR._alfa(self, cmp, cmp.Tf)[1]
            alfabp = PR._alfa(self, cmp, cmp.Tb)[1]
            alfaw = PR._alfa(self, cmp, 0.7*cmp.Tc)[1]

            # Eq B-4
            def f4(Tr):
                return Tr*log(Tr)+7/3*(Tr-1)*log(0.7)

            # Eq B-5
            def f5(q, Tr):
                return 10/3*(0.7**(q+1)-1)*(Tr-1)+Tr**(q+1)-1

            # Eq B-6
            def f6(alfaw, Tr):
                return 10/3*(1-alfaw)*(Tr-1)+1

            # Eq B-7
            def f7(q, Tr, Trbp):
                return f5(q, Tr)*f4(Trbp) - f5(q, Trbp)*f4(Tr)

            # Eq B-8
            def f8(Tr, Trbp, alfabp, alfaw):
                return f6(Tr, alfaw)*f4(Trbp) + f4(Tr)*(alfabp-f6(Trbp, alfaw))

            f1 = 1/f4(Trbp)/f7(1, Trtp, Trbp)
            f2 = f7(1, Tr, Trbp) * \
                (alfatp*f4(Trbp) - f8(Trtp, Trbp, alfabp, alfaw)) + \
                f8(Tr, Trbp, alfabp, alfaw)*f7(1, Trtp, Trbp)

            # Eq B-3
            def f3(q, Tr, Trbp, Trtp):
                return f7(q, Tr, Trbp)*f7(1, Trtp, Trbp) - \
                    f7(q, Trtp, Trbp)*f7(1, Tr, Trbp)

            alfa = f1*(f2+k1*f3(2, Tr, Trbp, Trtp)+k2*f3(11, Tr, Trbp, Trtp))

        else:
            m, alfa = PR._alfa(self, cmp, T)
        return 0, alfa


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRZV(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRZV(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
