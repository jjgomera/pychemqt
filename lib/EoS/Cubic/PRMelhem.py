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


from math import exp

from lib.EoS.Cubic.PR import PR


dat = {
    126 : (0.8313, -0.5701),
    127 : (0.7311, 0.3594),
    28 : (0.6644, 0.2885),
    24 : (0.6396, 0.4639),
    55 : (0.7303, 0.3486),
    159 : (1.0223, 1.7081),
    52 : (0.7825, 0.1751),
    32 : (0.7366, 0.0216),
    128 : (0.8317, -0.0149),
    130 : (1.0791, -0.5313),
    447 : (2.1183, -6.6355),
    140 : (0.8283, 0.1495),
    125 : (0.9813, 0.9529),
    65 : (0.6640, 0.2062),
    139 : (0.9016, -0.2409),
    63 : (0.7393, 0.2837),
    175 : (0.9227, 0.1891),
    98 : (0.3801, 0.2109),
    342 : (0.7979, 0.4790),
    40 : (0.6670, 0.4720),
    345 : (0.7357, 1.8414),
    99 : (0.5082, 0.5011),
    171 : (0.7497, 0.2442),
    49 : (0.6877, 0.3813),
    102 : (0.5924, -0.1033),
    48 : (0.4489, 0.3207),
    100 : (0.6433, 0.4615),
    105 : (0.5014, 0.3624),
    172 : (0.7511, 0.1902),
    112 : (0.6705, 0.4636),
    38 : (0.6651, 0.4814),
    36 : (0.6366, 0.5328),
    222 : (0.6955, 0.2458),
    163 : (2.4837, -6.5968),
    294 : (0.7712, 0.6816),
    133 : (0.6572, 0.3478),
    3 : (0.5336, 0.2431),
    134 : (1.2417, 0.1381),
    22 : (0.4873, 0.4570),
    135 : (2.2459, -4.9603),
    129 : (0.6499, 0.4672),
    155 : (0.8758, 0.4057),
    45 : (0.8084, 0.3129),
    162 : (0.7740, 0.3761),
    137 : (0.6495, 0.3559),
    197 : (0.4516, 0.3038),
    212 : (-0.1603, -0.9866),
    1 : (0.0829, -0.4780),
    209 : (0.6088, -0.6640),
    104 : (0.6018, -0.0507 ),
    113 : (1.1635, -2.0723),
    210 : (1.1982, -2.8266),
    106 : (0.5193, -0.1840),
    50 : (0.5050, 0.4719),
    213 : (0.5024, 0.5590),
    5 : (0.6384, 0.3959),
    24 : (0.6576, 0.3844),
    7 : (0.7052, 0.3361),
    145 : (1.1979, 0.8456),
    337 : (0.8718, 0.1339),
    159 : (1.0553, 1.5978),
    2 : (0.4045, 0.1799),
    117 : (1.2107, -0.5292),
    142 : (0.8555, 0.1447),
    118 : (0.7699, 0.5226),
    115 : (0.6068, 0.3011),
    153 : (0.8455, 0.2145),
    225 : (0.6730, 0.2581),
    131 : (0.7202, 0.5485),
    116 : (0.5815, 0.3936),
    43 : (0.8477, 0.2164),
    108 : (1.3541, -1.5199),
    46 : (0.4274, 0.3359),
    109 : (1.8386, -3.6754),
    226 : (0.9585, 0.6643),
    110 : (0.5981, 0.4349),
    6 : (0.6543, 0.4308),
    160 : (1.0705, 1.3402),
    14 : (1.0763, 0.0164),
    16 : (1.1778, 0.0531),
    11 : (0.8620, 0.3973),
    10 : (0.7939, 0.4116),
    13 : (1.0203, 0.0654),
    12 : (0.9273, 0.3748),
    8 : (0.7362, 0.3518),
    17 : (1.2243, 0.1019),
    15 : (1.1394, -0.0157),
    47 : (0.4137, 0.2784),
    214 : (0.6838, 0.2261),
    42 : (0.8283, 0.2267),
    174 : (1.0205, -0.0250),
    4 : (0.5909, 0.3997),
    146 : (1.1505, 0.8075),
    261 : (0.7630, 0.1883),
    23 : (0.5743, 0.4423),
    294 : (0.7124, 0.4424),
    44 : (0.8478, 0.1643),
    51 : (0.7238, 0.4537),
    111 : (0.4072, 6.1502),
    161 : (1.0636, 1.7938),
    41 : (0.7563, 0.3125),
    62 : (0.8893, 0.0151)}


class PRMelhem(PR):
    """Peng-Robinson cubic equation of state with a modified dependence of 
    temperature by Melhem

    .. math::
        \alpha = \exp\left(m\left(1-T_r\right)
        \left(1-\sqrt{T_r}\right)^2\right)
    
    where n and m are adjustable parameters compound specific
    """

    __title__="PR Melhem (1989)"
    __status__="PRMelhem"
    __doi__ = {
        "autor": "Melhem, G.A., Saini, R., Goodwin, B.M.",
        "title": "A Modified Peng-Robinson Equation of State",
        "ref": "Fluid Phase Equilibria 47 (1989) 189-237",
        "doi": "10.1016_0378-3812(89)80176-1"},

    def _alfa(self, cmp, T):
        # m,n parameters per compound, add a dict with values
        if cmp.id in dat:
            m, n = dat[cmp.id]
            alfa = exp(m*(1-T/cmp.Tc) + n*(1-(T/cmp.Tc)**0.5)**2)       # Eq 1
        else:
            alfa = PR._alfa(self, cmp, T)
        return 0, alfa


if __name__ == "__main__":
    
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRMelhem(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRMelhem(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))

