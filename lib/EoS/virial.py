#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.



Virial equation of state implementation

:class:`lib.EoS.virial.Virial`: The main class with all integrated
functionality

Second virial coefficient correlations:
    * :func:`lib.EoS.virial.B_Tsonopoulos`
    * :func:`lib.EoS.virial.B_IglesiasSilva`
    * :func:`lib.EoS.virial.B_Meng`
    * :func:`lib.EoS.virial.B_Orbey`
    * :func:`lib.EoS.virial.B_Tarakad`

Third virial coefficient correlations:
    * :func:`lib.EoS.virial.C_Orbey_Vera`
    * :func:`lib.EoS.virial.C_Liu_Xiang`
    * :func:`lib.EoS.virial.C_Meng`

.. include:: virial.rst

'''


from numpy import exp
from numpy.lib.scimath import log
from scipy import zeros
from scipy.constants import R
from scipy.misc import derivative

from lib.eos import EoS
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Dymond, J.H., Marsh, K.N., Wilhoit, R.C., Wong, K.C., "
                  "Frenkel, M.",
         "title": "Virial Coefficients of Pure Gases (Landolt-Börnstein - "
                  "Group IV Physical Chemistry 21A)",
         "ref": "Springer-Verlag",
         "doi": "10.1007/10693952_16"},
    2:
        {"autor": "Dymond, J.H., Marsh, K.N., Wilhoit, R.C., Wong, K.C., "
                  "Frenkel, M.",
         "title": "Virial Coefficients of Mixtures "
                  "(Landolt-Börnstein - Group IV Physical Chemistry 21B)",
         "ref": "Springer-Verlag",
         "doi": "10.1007/10754889"},
    3:
        {"autor": "Tsonopoulos, C.",
         "title": "An Empirical Correlation of Second Virial Coefficients",
         "ref": "AIChE Journal 20(2) (1974) 263-272",
         "doi": "10.1002/aic.690200209"},
    4:
        {"autor": "Orbey, H., Vera, J.H.",
         "title": "Correlation for the Third Virial Coefficient Using Tc, Pc "
                  "and ω as Parameters",
         "ref": "AIChE Journal 29(1) (1983) 107-113",
         "doi": "10.1002/aic.690290115"},
    5:
        {"autor": "Liu, D.X., Xiang, H.W.",
         "title": "Corresponding-States Correlation and Prediction of Third "
                  "Virial Coefficients for a Wide Range of Substances",
         "ref": "Int. J. Thermophysics 24(6) (2003) 1667-1680",
         "doi": "10.1023/b_ijot.0000004098.98614.38"},
    6:
        {"autor": "Iglesias-Silva, G.A., Hall K.R.",
         "title": "An Equation for Prediction and/or Correlation of Second "
                  "Virial Coefficients",
         "ref": "Ind. Eng. Chem. Res. 40(8) (2001) 1968-1974",
         "doi": "10.1021/ie0006817"},
    7:
        {"autor": "Meng, L., Duan, Y.Y. Li, L.",
         "title": "Correlations for Second and Third Virial Coefficients of "
                  "Pure Fluids",
         "ref": "Fluid Phase Equilibria 226 (2004) 109-120",
         "doi": "10.1016/j.fluid.2004.09.023"},
    8:
        {"autor": "Orbey, H.",
         "title": "A Four Parameter Pitzer-Curl Type Correlation of Second "
                  "Virial Coefficients",
         "ref": "Chemical Engineering Comm. 65(1) (1988) 1-19",
         "doi": "10.1080/00986448808940239"},
    9:
        {"autor": "Tarakad, R.R., Danner, R.P.",
         "title": "An Improved Corresponding States Method for Polar Fluids: "
                  "Correlation of Second Virial Coefficients",
         "ref": "AIChE J. 23(5) (1977) 685-695",
         "doi": "10.1002/aic.690230510"}}


# Parameters for second virial coefficient for pure compounds, ref [1]_
B_Database = {
    98: [3.4162e1, -1.2087e4, -7.6702e5, -1.96e7],
    630: [9.1039e1, -5.9081e4, 1.0478e7, 3.0463e9],
    113: [-6.0611e3, 4.6389e6, 9.8352e8],
    48: [4.8826e1, -1.5614e4, -2.757e5, -4.7684e7],
    49: [5.74e1, -3.8829e4, 4.2899e5, -1.4661e9],
    102: [2.8415e3, -2.8985e6, 9.583e8, -1.2515e11],
    # Chlorine triluoride, pag 36
    # 215: [-1.0754e3, 9.3227e5, -2.7872e8],
    104: [8.1566e1, -9.7274e4, 2.7321e7, -5.332e9],
    105: [1.3171e1, 5.1055e3, -2.9404e7, -5.0383e8],
    208: [3.3609e1, -1.0625e4, -6.078e5, -2.22759e7],
    # Inorganic fluorides pag 39-46
    62: [1.5883e2, -3.0107e5, 1.8189e8, -5.6932e10],
    50: [3.6507e1, -7.9791e4, 2.5927e7, -6.7852e9],
    63: [1.5198e2, 2.0281e5, 8.1489e7, -1.7248e10],
    # Phosphine pag 61
    # Krypton pag 69
    108: [9.7906e1, -6.2888e4, 1.1825e7, -1.0802e9],
    46: [4.0286e1, -9.3378e3, -1.4164e6, 6.1253e7, -2.7198e9],
    110: [5.6654e1, 5.4768e4, 9.4786e6, -2.992e9],
    107: [1.5894e1, -9.9406e2, -1.2641e5, 2.2721e6],
    47: [4.2859e1, -1.7696e4, 5.2007e5, -1.6393e8, 5.0855e9],
    51: [2.4794e1, -2.3892e4, -4.9786e6, -8.2606e9],
    # Xenon pag 84
    636: [-3.0298e2, 1.787e5, -5.2361e7],
    215: [8.526e1, -3.6077e4, -1.7214e7, 1.0634e8],
    216: [2.0506e2, -2.1105e5, 3.941e7, -1.0597e10],
    217: [2.7564e2, -2.5324e5, -7.489e6],
    100: [-1.3094e3, 9.9093e5, -3.1965e8],
    # trifluoroiodomethane pag 92
    218: [8.2106e1, -3.0923e4, -5.9874e6, 9.5915e6],
    220: [9.8669e1, -9.879e4, 1.5682e7, -7.9366e9],
    642: [4.1796e2, -3.2472e5, 8.1294e6, 2.9564e8],
    112: [-3.2242e4, 3.3904e7, -1.1904e10, 1.3609e12],
    643: [4.3496e1, -2.6827e4, -3.6456e6, -2.6677e9],
    222: [6.6308e2, -9.0425e5, 3.6968e8, -7.0487e10],
    645: [1.5688e2, -1.0739e5, 1.2759e7, -6.5037e9],
    224: [1.6048e3, -1.751e6, 6.1379e8, -8.4605e10],
    115: [2.332e2, -1.9046e5, 1.9734e7, -5.576e9],
    225: [4.6669e2, -4.7431e5, 1.4966e8, 2.0299e10],
    116: [1.1563e3, -5.7779e5, 3.1037e5],
    226: [-2.5351e4, 1.9039e7, -3.8114e9],
    2: [4.4344e1, -1.6608e4, -3.543e6, 2.9832e8, 2.3448e10],
    117: [2.1043e3, -3.2372e6, 1.6677e9, -3.1011e11],
    118: [1.021e2, -1.6761e5, 5.4579e7, -1.6593e10],
    229: [-1.9198e2, 2.1731e5, -9.2469e7],
    230: [7.9255e2, -4.6603e5, -1.4163e6],
    231: [2.3066e2, -2.3e5, 2.1747e7, -1.3611e10],
    232: [5.1747e2, -4.6798e5, 4.726e6],
    235: [1.3583e2, -9.2618e4],
    236: [2.1727e3, -1.6683e6, 4.0046e8, -3.5591e10],
    654: [2.8867e2, -2.1746e5],
    # 1-cloro-1,2,2,2-tetrafluoroetano pag 118
    # 1,1-dicloro-2,2,2-tetrafluoroetano pag 119
    119: [1.3094e3, -8.7959e5],
    # Pentafluroetano pag 120
    # Difluromethoxy-triflurometano pag 122
    65: [-7.8931e2, 5.2035e5, -1.0229e8, 1.2008e9],
    239: [1.4779e2, -9.305e4, -2.3775e6],
    # 1,1,1,2-tetrafluroetano pag 125
    122: [2.174e2, -1.945e5],
    241: [-1.3226e2, 5.4763e4, -1.875e7, -1.252e10],
    # 1,1-dicloro-1-fluoroetano pag 129
    479: [6.5521e2, -6.0452e5],
    243: [7.0822e2, -6.8308e5, 2.0059e8, -2.8667e10],
    # 2,2,2-Trifluroetanol pag 132
    125: [3.6936e3, -6.0221e6, 3.2337e9, -6.6701e11],
    22: [8.7776e1, -6.0061e4, 1.2857e6, -1.0861e9],
    126: [-5.2377e2, 4.5869e5, -1.9557e8],
    127: [9.0835e2, -1.2761e6, 5.2426e8, -1.0654e11],
    245: [1.3654e3, -1.471e6, 5.0537e8, -6.9611e10],
    128: [2.8626e3, -3.53e6, 1.4489e9, -2.2614e11],
    129: [1.9619e2, -8.6867e4, -4.5348e7],
    130: [-5.2743e4, 6.2086e7, -1.8447e10],
    131: [-3.0469e3, 2.3507e6, -5.2465e8],
    246: [1.08e3, -5.5106e5],
    132: [-5.8568e2, 7.9066e5, -3.5796e8, 3.2306e10],
    3: [1.0773e2, -8.2548e4, 5.2387e6, -1.9764e9],
    133: [-7.9769e2, 5.3083e5, -1.2938e8],
    134: [9.6838e3, -1.3575e7, 6.3248e9, -1.0114e12],
    136: [1.187e3, -6.263e5],
    138: [2.9097e3, -2.9479e6, 9.6424e8, -1.2336e11],
    249: [2.5017e2, -1.0569e5, -4.9303e7],
    237: [1.2333e2, -2.4851e4, -3.7728e7],
    # 2,2,4,4,5,5-Hexafluro-1,3-dioxonale pag 153
    671: [9.1463e2, -4.4617e5],
    # Heptafluropropanes pag 153
    66: [-4.9112e2, 5.4036e5, -2.1854e8, 1.9163e10],
    57: [2.6664e3, -2.2784e6, 6.0977e8, -6.0081e10],
    258: [-1.2251e3, 1.337e6, -5.1254e8, 5.6467e10],
    23: [1.0101e2, -7.5735e4, -7.9502e6, -2.7987e9],
    # 2,2-dicloropropane pag 159
    484: [-3.1654e3, 2.2402e6, -4.78e8],
    140: [5.1325e3, -6.4311e6, 2.6952e9, -4.2124e11],
    261: [-2.3257e4, 1.5352e7, -2.6699e9],
    141: [-5.6624e3, 4.3602e6, -9.3449e8],
    142: [-2.7935e2, 5.7525e5, -2.8774e8],
    682: [-9.2253e2, 7.2399e5, -2.5425e8],
    263: [2.6714e2, -3.2743e5, 9.285e7, -3.4217e10],
    264: [6.6635e2, -4.3329e5, -1.4434e7],
    4: [1.0971e2, -8.4673e4, -8.1215e6, -3.4382e9],
    486: [1.0121e3, -5.6471e5, 2.5578e7],
    146: [5.6134e3, -8.3867e6, 4.1544e9, -7.2489e11],
    145: [1.0296e4, -1.414e7, 6.4638e9, -1.0248e12],
    # 2-propanethiol pag 169
    147: [-1.4209e3, 1.0335e6, -2.4971e8],
    # Trimetilborato pag 170
    692: [4.2109e2, -4.3943e5, 1.2859e8, -3.1178e10],
    693: [2.0874e3, -1.223e6, 9.7474e7],
    # 1,1,1,2,2,3,3,4-octaflurobutano pag 172
    # 2-clorotiophene pag 173
    273: [6.2291e3, -6.1819e6, 1.9793e9, -2.2748e11],
    149: [-3.0669e3, 2.2018e6, -5.088e8],
    28: [1.7027e4, -1.6558e7, 5.3391e9, -5.8813e11],
    151: [7.8919e2, 4.9648e5],
    # 2,5-dihydrofuran pag 175
    278: [-6.8713e4, 5.4024e7, -1.0864e10],
    492: [3.8855e3, -2.156e6],
    24: [1.8449e3, -2.4771e6, 1.1518e9, -2.4483e11, 1.6568e13],
    26: [-1.2888e3, 8.6278e5, -2.07e8],
    25: [-2.461e3, 1.6351e6, -3.3194e8],
    27: [-1.5536, 2.6336e6, -1.6613e9, 4.2568e11, -4.187e13],
    448: [-2.4919e3, 1.7921e6, -4.5499e8],
    153: [-7.3871e3, 5.8855e6, -1.3266e9],
    282: [-6.9621e2, 6.8028e5, -3.0414e8],
    155: [-6.583e3, 5.1996e6, -1.167e9],
    156: [-3.5353e3, 2.975e6, -7.697e8],
    157: [-4.6773e3, 3.8319e6, -9.1872e8],
    712: [1.2098e3, -8.7282e5],
    284: [8.0474e2, -1.114e6, 4.4699e8, -9.9235e10],
    6: [2.272e2, -2.2797e5, 2.9855e7, -1.3706e10],
    5: [1.1625e2, -1.0293e5, -1.2475e7, -7.049e9],
    160: [1.2021e2, 1.589e6, 1.3878e9, -3.6783e11],
    161: [2.8091e3, -1.4181e6],
    162: [-2.3331e3, 1.9751e6, -5.4039e8, 1.5345e10],
    159: [2.6055e3, -1.431e6],
    450: [2.9765e3, -1.5529e6],
    294: [7.4083e3, -7.3873e6, 2.3796e9, -2.8822e11],
    # Tetrametilsilane pag 193
    # Dodecafluoropentano pag 195
    295: [3.4839e3, -2.4252e6, 2.68e8],
    # 2-metilfuran pag 196
    # metiltiofenos pag 197
    69: [-2.6712e3, 1.9665e6, -4.4974e8],
    297: [2.1158e3, -1.1421e6],
    724: [1.2261e3, -8.0976e5],
    29: [1.95e4, -1.1193e7, 3.9052e9, -4.7699e11],
    165: [-9.9989e3, 8.5387e6, -1.999e9],
    # Cyclopentanethiol pag 202
    # 2-cloro-2-metylbutane pag 203
    9: [4.9697e2, -4.7244e5, 1.031e8, -2.3475e10],
    7: [-2.0233e2, 1.8489e5, -1.2348e8, -2.3136e9],
    8: [2.9122e2, -3.3475e5, 6.5631e7, -2.8478e10],
    # butyl methyl ether pag 208
    743: [2.7106e3, -1.6467e6],
    # 3,3-dimetil-2-tiabutano pag 209
    # 2-methyl-1-butanethiol pag 210
    319: [1.0224e3, -1.5575e6, 7.5148e8, -1.8001e11],
    320: [-5.0538e1, 2.2888e5, -2.1502e8],
    # 2,3-bis(Trifluromethyl)-perfluorobutane pag 212
    321: [-3.7991e2, 6.1865e5, -3.4962e8],
    # Undecafluoro-2-(trifluoromethyl)-pentane pag 213
    172: [-4.7587e2, 6.75e5, -3.9453e8],
    322: [3.782e2, -4.2848e5, 7.0203e7, -3.4266e10],
    40: [4.7946e2, -6.8047e5, 2.3851e8, -6.2693e10],
    174: [-1.8098e3, 2.4621e6, -1.1103e9, 5.3333e10],
    622: [-2.0978e4, 2.4607e7, -9.5935e9, 1.166e12],
    623: [2.7217e3, -1.6134e6, -2.9073e7],
    323: [7.3393e2, -4.9526e4, -3.3521e8],
    38: [7.3023e1, -1.2813e5, -1.3635e7, -2.8581e10],
    35: [-4.1329e3, 4.2605e6, -1.5173e9, 1.3838e11],
    37: [1.3965e3, -8.6726e5],
    54: [5.1719e2, -3.8352e5, -4.0872e7],
    55: [7.7373e2, -5.3442e5, -4.6235e7],
    52: [-7.0489e2, 7.7395e5, -3.2471e8],
    53: [7.5391e2, -3.9836e5, -9.7628e7],
    10: [1.4421e3, -1.9714e6, 8.1722e8, -1.589e11],
    339: [1.2207e3, -6.9533e5, -7.1981e7],
    # Hexamethyldisiloxane pag 230
    340: [1.1512e4, -1.5192e7, 6.6591e9, -1.0255e12],
    524: [2.7009e3, -1.5608e6],
    343: [-1.3687e5, 1.5655e8, 4.5344e10],
    # fluoro-metil benzenos pag 233
    41: [1.9007e3, -2.7592e6, 1.2424e9, -2.4902e11],
    344: [-2.2836e3, 2.6153e6, -9.5452e8],
    348: [3.1208e3, -2.0276e6],
    # 2,4-dimethylpyridine pag 238
    349: [2.724e3, -1.8369e6],
    # 2,6-dimethylpyridine pag 238
    56: [-4.3752e2, 6.8995e5, -3.9007e8],
    11: [8.7046e2, -1.3176e6, 5.6372e8, -1.4992e11],
    42: [-1.0007e4, 8.9315e6, -2.2448e9],
    43: [-6.3438e3, 6.3946e6, -1.8133e9],
    44: [3.5344e3, -5.2399e6, 2.5224e9, -5.013e11],
    45: [-3.4316e3, 3.2595e6, -9.9758e8],
    83: [4.5025e3, -2.9103e6, 1.8712e8],
    582: [1.5963e3, -1.0012e6, -7.3728e7],
    600: [1.5335e3, -8.8627e5, -1.1686e8],
    # Octamethyltrisilosane pag 247
    12: [2.5227e3, -3.8225e6, 1.8393e9, -3.9652e11]}


# Parameters for second cross virial coefficient, ref [2]_
Bij_Database = {
    "98-48": [3.1874e1, 8.687e3, -1.4896e6, 2.5386e7],
    "98-49": [-7.2525e2, 9.6341e5, -4.5827e8, 9.1873e10, -6.7735e12],
    "98-104": [1.6897e1, -8.7669e2, -4.8549e6, 1.5096e8],
    # Ar-F6S, pag 23
    "98-1": [3.3749e1, 8.2782e3, -7.4655e4, 7.1912e7, -5.7163e9],
    "98-62": [1.7025e3, -1.8921e6, 6.9275e8, -8.4515e10],
    "98-63": [-1.1969e1, 2.1313e4, -8.7141e6],
    "98-212": [1.4965e1, 6.8188e3, -2.6724e6, 3.2953e8, -1.3616e10],
    # Ar-Kr, pag 29
    "98-46": [3.9467e1, -1.4212e4, -1.8873e5, -4.7010e7],
    "98-107": [2.9286, 1.2477e4, -4.1259e6, 3.7687e8, -1.2338e10],
    "98-47": [4.8775e1, -2.3064e4, 1.7036e6, -1.6031e8],
    "98-115": [-7.6625e2, 9.4837e5, -3.5161e8, 3.9333e10],
    "98-2": [-1.2034e1, 3.569e4, -1.7184e7, 1.9559e9, -8.5469e10],
    "98-3": [9.2713e1, -5.6437e4, 7.0318e6, -9.3813e8],
    "98-4": [1.3285e2, -9.7467e4, 1.5052e7, -1.8464e8],
    "98-6": [-2.8383e2, 3.3643e5, -1.3313e8, 1.4337e10],
    # Ar-C4H12Si, pag 49
    "98-9": [9.6764e1, -5.1023e4, -3.848e6],
    "48-49": [-2.4271, 2.7027e4, -1.7561e7, 1.6094e9],
    "48-1": [4.4409e1, -9.8417e3, 1.7911e5, -7.7626e6],
    "48-212": [3.4426e1, -5.6592e3, 7.7493e5, -4.7919e7],
    "48-107": [3.3843e1, 6.4861e3, 2.5955e5, -2.4332e7],
    "48-2": [2.2044e3, -1.9304e6, 5.6339e8, -5.5377e10],
    "48-40": [2.4487e2, -2.4165e5, 7.124e7, -9.9105e9],
    "48-12": [2.4481e2, -1.6544e5, 1.32e7],
    "49-1": [2.2598e2, -1.7768e5, 4.7952e7, -4.5198e9],
    "49-62": [-1.0744e2, 1.1123e5, -4.0394e7],
    "49-212": [4.8515e1, -7.9326e3],
    "49-46": [1.8683e1, 1.9172e4, -2.0167e7, 3.425e9, -2.1815e11],
    "49-2": [1.6073e2, -1.3829e5, 3.5911e7, -4.3948e9],
    "49-3": [-2.8002e3, 2.458e6, -7.2439e8, 6.8408e10],
    "49-40": [-1.7625e1, -8.5149e4],
    "49-38": [-1.9488e2, 3.9274e5, -2.4825e8, 3.8182e10],
    # F6S-He, pag 114
    # F6S-Kr, pag 115
    # F6S-CF4, pag 118
    # F6S-CH4, pag 121
    # F6S-C4H12Si, pag 122
    "1-212": [1.8618e1, -7.5727e2, 6.8549e2, -3.1665e3],
    # H2-Kr, pag 125
    "1-46": [2.2959e1, -2.6461e3, 4.9944e4, -4.5678e5],
    "1-107": [3.5658e1, -7.4422e3, 1.4167e5, -7.5396e6],
    # H2-Xe, pag 130
    "1-100": [5.8322e3, -5.0461e6, 1.4589e9, -1.4108e11],
    "1-218": [6.6438e1, -1.9787e4, 1.4512e6, -1.1631e8],
    "1-2": [-1.0011e2, 7.6037e4, -1.2943e7],
    "1-40": [7.1836e1, -2.5863e4],
    "1-12": [1.0789e2, -3.2428e4],
    "62-46": [2.1501e2, -2.1529e5, 7.2238e7, -9.0135e9],
    "62-47": [3.6866e1, -1.014e4, -3.3711e6, 2.1458e8],
    "62-2": [6.7399e1, -3.5456e4],
    "62-117": [2.7643e3, -3.7597e6, 1.7211e9, -2.8931e11],
    "62-3": [-2.9694e1, 2.9111e4, -1.5954e7],
    "62-134": [-1.6626e3, 1.6933e6, -4.7847e8],
    "62-4": [1.9447e2, -1.0179e5],
    "62-322": [-8.1996e2, 7.8141e5, -2.1661e8],
    "62-10": [1.6413e2, -1.0807e5],
    "62-41": [5.6094e3, -8.2825e6, 4.0104e9, -6.5703e11],
    "62-11": [-5.7138e2, 4.767e5, -1.2224e8],
    "62-12": [1.7899e2, -1.2172e5],
    "50-2": [5.4062e1, -2.1818e4, -2.1814e6],
    "63-212": [3.6769e1, -8.119e3],
    "63-46": [-1.4663e1, 2.5143e4, +9.6702e6],
    "63-2": [-7.5406e1, 6.6052e4, -1.8725e7],
    # He-Kr, pag 177
    "212-108": [2.5168e1, -4.8809e2, -2.3549e5],
    "212-46": [1.8417e1, 3.045e3, -7.3264e5, 3.6364e7, -6.3776e8],
    "212-110": [4.4476e1, -9.8008e3, 7.9041e5],
    "212-107": [1.3558e1, 2.5654e2, -2.8502e4, 2.206e5],
    "212-47": [2.6501e1, -2.5239e3, -1.9725e4, -2.9539e4],
    # He-Xe, pag 187
    "212-218": [9.6369e1, -3.6626e4, 9.6817e6, -9.1723e8],
    "212-2": [2.9872e1, -9.6196e2, -2.3581e5, 1.7281e7],
    "212-4": [6.2759e1, -1.4565e4, 2.1194e6, 1.9238e8],
    "212-6": [-1.2344e1, 5.8166e4, -1.3241e7],
    # Kr-Ne, pag 202
    # Kr-Xe, pag 204
    # Kr-CH4, pag 206
    "46-107": [3.2804e1, -5.4303e3, -1.2152e5, 8.2859e6],
    "46-47": [3.8110e1, -1.2265e4, 6.6992e5, -3.2660e7, 1.1664e9],
    "46-222": [8.1349e2, -9.0416e5, 3.2388e8, -4.0768e10],
    "46-2": [2.3092e1, -1.6786e3, -3.3822e6, 7.1634e7],
    "46-117": [-3.4678e2, 2.6522e5, -6.9489e7, 4.2123e9],
    "46-4": [4.3087e2, -4.013e5, 1.2091e8, -1.3638e10],
    "46-40": [1.307e2, -7.2389e4],
    "46-38": [-2.5418e2, 2.1825e5, -5.4381e7],
    "46-41": [-1.7879e2, 1.2731e5, -3.3317e7],
    "46-11": [1.6425e3, -1.5234e6, 4.6319e8, -4.9785e10],
    "107-47": [3.2788e1, -8.0888e3, 2.9009e5, -1.6447e7],
    # Ne-Xe, pag 243
    "107-2": [3.5364e1, -5.7379e3, -1.4498e5, 6.8207e6],
    # Xe-C2H6, pag 255
    "218-2": [8.6407e1, -5.9805e4, 1.8721e7, -5.6828e9, 5.4603e11],
    "220-231": [-2.1738e1, 4.0729e4, -5.2896e7],
    "112-140": [-1.6239e4, 1.2491e7, -2.5216e9],
    "112-162": [-4.1677e3, 3.8019e6, -9.5242e8],
    "112-40": [1.5715e3, -9.3146e5],
    "222-40": [-5.7315e2, 5.2812e5, -2.1972e8],
    "2-117": [-1.286e3, 8.0897e5, -1.372e8],
    "2-22": [7.3929e1, -4.4007e4],
    "2-3": [4.9899e1, -2.227e4, -5.8719e6],
    "2-4": [1.4877e2, -8.3334e4],
    "2-6": [1.2142e2, -8.589e4],
    # CH4-C4H12Si, pag 296
    "2-8": [-5.3644e2, 3.9304e5, -8.9133e7],
    "2-9": [1.3314e2, -8.4264e4, -3.3841e6],
    "2-10": [3.0763e2, -1.6881e5],
    "2-11": [-1.6715e2, 1.4053e5, -5.3187e7],
    "117-40": [6.4967e1, 9.8356e3, -6.7332e7],
    "117-38": [9.6586e1, -4.5832e4, -3.4819e7],
    "117-10": [-2.7957e2, 3.3896e5, -1.2459e8],
    "3-4": [2.7832e2, -1.6449e5],
    "3-6": [2.3606e2, -1.7765e5],
    "3-8": [-9.987e2, 6.7418e5, -1.5158e8],
    "134-40": [4.5376e2, -3.5275e5],
    "134-38": [7.8868e1, -2.2292e4, -5.9299e7],
    "140-40": [1.4272e3, -7.9349e5],
    "4-6": [1.6045e2, -5.6412e4, -4.5336e7],
    "4-10": [-4.1284e2, 3.6086e5, -1.4408e8],
    "693-6": [9.3979e2, -4.8933e5],
    "40-38": [1.8243e3, -9.9756e5]}


def _B_Database(T, args):
    """Calculate second virial coefficient, its 1st and 2nd temperature
    derivarives from coefficient in database

    Parameters
    ----------
    T : float
        Temperature, [K]
    args : list
        Array with polinomial coefficient

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]
    """
    B = 0
    Bt = 0
    Btt = 0
    for i, a in enumerate(args):
        B += a/T**i
        Bt -= a*i*T**(i-1)
        Btt += a*i*(i-1)*T**(i-2)
    return B/1000, Bt/1000, Btt/1000


@refDoc(__doi__, [3])
def B_Tsonopoulos(T, Tc, Pc, w, mu=None):
    r"""Calculate the 2nd virial coefficient using the Tsonopoulos correlation

    .. math::
        \frac{BP_c}{RT_c} = f^{(0)}+\omega f^{(1)}+f^{(2)}+f^{(3)}

    .. math::
        f^{(0)} = 0.1445 - \frac{0.33}{T_r} - \frac{0.1385}{T_r^2} -
        \frac{0.0121}{T_r^3} - \frac{0.000607}{T_r^8}

    .. math::
        f^{(1)} = 0.0637 - \frac{0.331}{T_r^2} - \frac{0.423}{T_r^3} -
        \frac{0.008}{T_r^8}

    .. math::
        f^{(2)} = \frac{a}{T_r^6}

    .. math::
        f^{(3)} = \frac{b}{T_r^8}

    .. math::
        a = -2.14e-4\mu_r-4.308e-21*\mu_r^8

    .. math::
        b = 0.00908+0.0006957*\mu_r

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    w : float
        Acentric factor [-]
    mu : float, optional
        dipole moment [debye]

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]

    Notes
    -----
    With the B_database this correlation is only for completeness, the a, b
    correlations are general and possibly not applicable to the compounds
    not availables in the B_database

    Examples
    --------
    Selected date from Table 2, pag 74 for neon

    >>> from lib.mEoS import Ne
    >>> D = Ne.momentoDipolar
    >>> "%0.4f" % B_Tsonopoulos(262, Ne.Tc, Ne.Pc, Ne.f_acent, D)[0]
    '0.0113'
    """
    Tr = T/Tc
    if mu:
        mur = mu**2*Pc/1.01325/Tc**2
        a = -2.14e-4*mur-7.831e-21*mur**8
        b = 0.00908+0.0006957*mur
    else:
        a, b = 0, 0

    def f(T):
        f0 = 0.1445-0.33/Tr-0.1385/Tr**2-0.0121/Tr**3-0.000607/Tr**8
        f1 = 0.0637+0.331/Tr**2-0.423/Tr**3-0.008/Tr**8
        f2 = 1/Tr**6
        f3 = -1/Tr**8
        f = f0 + w*f1 + a*f2 + b*f3
        return f*R*Tc/Pc*1e3

    B = f(T)
    B1 = derivative(f, T, n=1)
    B2 = derivative(f, T, n=2)
    return B, B1, B2


@refDoc(__doi__, [6])
def B_IglesiasSilva(T, Tc, Pc, Vc, w, D):
    r"""Calculate the 2nd virial coefficient using the Iglesias-Silva Hall
    correlation

    .. math::
        \frac{B}{b_o} = \left(\frac{T_B}{T}\right)^{0.2}
        \left[1-\left(\frac{T_B}{T}\right)^{0.8}\right] \left[\frac{B_c}
        {b_o\left(\left(T_B/T_C\right)^{0.2}-\left(T_B/T_C\right)\right)}
        \right]^{\left(T_c/T\right)^n}

    .. math::
        \frac{B_C}{V_C} = -1.1747 - 0.3668\omega - 0.00061\mu_R

    .. math::
        n = 1.4187 + 1.2058\omega

    .. math::
        \frac{b_o}{V_C} = 0.1368 - 0.4791\omega + 13.81\left(T_B/T_C\right)^2
        \exp\left(-1.95T_B/T_C\right)

    .. math::
        \frac{T_B}{T_C} = 2.0525 + 0.6428\exp\left(-3.6167\omega\right)

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure, [Pa]
    Vc : float
        Critical specific volume, [m³/mol]
    w : float
        Acentric factor [-]
    D : float
        dipole moment [debye]

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]

    Examples
    --------
    Selected date from Table 2, pag 74 for neon

    >>> from lib.mEoS import Ne
    >>> Vc = Ne.M/Ne.rhoc
    >>> D = Ne.momentoDipolar.Debye
    >>> "%0.4f" % B_IglesiasSilva(262, Ne.Tc, Ne.Pc, Vc, Ne.f_acent, D)[0]
    '0.0102'
    """
    # Reduced dipole moment
    muR = D**2*Pc/101325/Tc**2

    TB = Tc*(2.0525 + 0.6428*exp(-3.6167*w))                           # Eq 19
    n = 1.4187 + 1.2058*w                                              # Eq 14
    bo = Vc*(0.1368-0.4791*w+13.81*(TB/Tc)**2*exp(-1.95*TB/Tc))        # Eq 15
    Bc = Vc*(-1.1747-0.3668*w-0.00061*muR)

    def f(T):
        # Eq 12
        return bo*(TB/T)**.2 * (1-(TB/T)**.8) * \
            (Bc/bo/((TB/Tc)**.2-TB/Tc))**((Tc/T)**n)

    B = f(T)
    B1 = derivative(f, T, n=1)
    B2 = derivative(f, T, n=2)

    return B, B1, B2


@refDoc(__doi__, [7])
def B_Meng(T, Tc, Pc, w, D):
    r"""Calculate the 2nd virial coefficient using the Meng-Duan-Li correlation

    .. math::
        \frac{BP_c}{RT_c} = f^{(0)}+\omega f^{(1)} + f^{(2)}

    .. math::
        f^{(0)} = 0.13356 - \frac{0.30252}{T_r} - \frac{0.15668}{T_r^2} -
        \frac{0.00724}{T_r^3} - \frac{0.00022}{T_r^8}

    .. math::
        f^{(1)} = 0.17404 - \frac{0.15581}{T_r^2} + \frac{0.38183}{T_r^3} -
        \frac{0.44044}{T_r^3} - \frac{0.00541}{T_r^8}

    .. math::
        f^{(2)} = \frac{a}{T_r^6}

    .. math::
        a = -3.0309e-6\mu_r^2+9.503e-11\mu_r^4-1.2469e-15\mu_r^6

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    w : float
        Acentric factor [-]
    D : float
        dipole moment [debye]

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]

    Examples
    --------
    Selected date from Table 2, pag 74 for neon

    >>> from lib.mEoS import Ne
    >>> D = Ne.momentoDipolar.Debye
    >>> "%0.4f" % (B_Meng(262, Ne.Tc, Ne.Pc, Ne.f_acent, D)[0])
    '0.0099'
    """

    mur = D**2*Pc/1.01325/Tc**2
    a = -3.0309e-6*mur**2 + 9.503e-11*mur**4 - 1.2469e-15*mur**6

    def f(T):
        Tr = T/Tc
        f0 = .13356 - .30252/Tr - .15668/Tr**2 - .00724/Tr**3 - .00022/Tr**8
        f1 = .17404 - .15581/Tr + .38183/Tr**2 - .44044/Tr**3 - .00541/Tr**8
        f2 = 1/Tr**6
        f = f0 + w*f1 + a*f2
        return f*R*Tc/Pc*1e3

    # There are a ampliation to associating fluid, only applicable to alcohols,
    # amines and water, and with the necessity of chemical type in database
    # Meng, L., Duan, Y.-Y.
    # An Extended Correlation for Second Virial Coefficients of Associated and
    # Quantum Fluids
    # Fluid Phase Equilibria 258 (2007) 29-33
    # doi: 10.1016/j.fluid.2007.05.010

    B = f(T)
    B1 = derivative(f, T, n=1)
    B2 = derivative(f, T, n=2)

    return B, B1, B2


@refDoc(__doi__, [8])
def B_Orbey(T, Tc, Pc, w, id=None):
    r"""Calculate the 2nd virial coefficient using the Orbey correlation

    .. math::
        \frac{BP_c}{RT_c} = f^{(0)} + \omega f^{(1)} + f^{(2)}

    .. math::
        f^{(0)} = 0.1479 - \frac{0.3321}{T_r} - \frac{0.02907}{T_r^2} -
        \frac{0.06849}{T_r^3}

    .. math::
        f^{(1)} = 0.2473 - \frac{0.6092}{T_r} + \frac{1.0749}{T_r^2} -
        \frac{0.7569}{T_r^3}

    .. math::
        f^{(2)} = \frac{a}{T_r^7}

    where a is compound dependend parameters

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    w : float
        Acentric factor [-]
    id : integer
        Index of compound in database

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]

    Examples
    --------
    Using same value of Tsonopoulos correlation as reference

    >>> from lib.mEoS import Ne
    >>> "%0.4f" % B_Orbey(262, Ne.Tc, Ne.Pc, Ne.f_acent)[0]
    '0.0104'
    """
    Tr = T/Tc

    par = {
        100: -0.0021,
        220: -0.0021,
        224: -0.0060,
        115: -0.0065,
        225: -0.0550,
        226: -0.0320,
        117: -0.0350,
        227: -0.0038,
        118: -0.0230,
        102: -0.0010,
        235: -0.0030,
        236: -0.0069,
        119: -0.0050,
        65: -0.0300,
        125: -0.0730,
        126: -0.0031,
        128: -0.0340,
        129: -0.0092,
        131: -0.0140,
        132: -0.0085,
        133: -0.0120,
        134: -0.0040,
        137: -0.0015,
        249: -0.0115,
        138: -0.0130,
        66: -0.0010,
        140: -0.0230,
        141: -0.0045,
        142: -0.0145,
        264: -0.0030,

        486: -0.0054,
        304: -0.0150,
        305: -0.0083,
        165: -0.0140,
        166: -0.0080,
        309: -0.0066,
        7: -0.0020,
        318: -0.0012,
        319: -0.0077,
        322: -0.0013,
        174: -0.0100,
        323: -0.0020,
        # Methyl isobutyl ketone: -0.0210,
        41: -0.0025,
        45: -0.0023,
        42: -0.0050,
        43: -0.0055,
        145: -0.0060,
        147: -0.0030,
        270: -0.0116,
        67: -0.0046,
        27: -0.0020,
        25: -0.0020,
        26: -0.0020,
        153: -0.0130,
        157: -0.0045,
        155: -0.0065,
        156: -0.0082,
        162: -0.0005,
        290: -0.0015,
        294: -0.0040,
        29: -0.0080,
        34: -0.0018,
        44: -0.0066,

        71: -0.0066,
        72: -0.0022,
        74: -0.0022,
        75: -0.0022,
        76: -0.0021,
        77: -0.0013,
        208: -0.0090,
        62: -0.0090,
        63: -0.0220,
        51: -0.0060}

    if id and id in par:
        a = par[id]
    else:
        a = 0

    def f(T):
        f0 = 0.1479 - 0.3821/Tr - 0.02907/Tr**2 - 0.06849/Tr**3
        f1 = 0.2473 - 0.6092/Tr + 1.0749/Tr**2 - 0.7569/Tr**3
        f2 = a/Tr**6
        f = f0 + w*f1 + f2
        return f*R*Tc/Pc*1e3

    B = f(T)
    B1 = derivative(f, T, n=1)
    B2 = derivative(f, T, n=2)
    return B, B1, B2


@refDoc(__doi__, [9])
def B_Tarakad(T, Tc, Pc, id):
    r"""Calculate the 2nd virial coefficient using the Tarakad-Danner
    correlation

    .. math::
        B^o = \frac{BP_c}{RT_c} = B_{simple}^o + B_{size-shape}^o + B_{polar}^o

    .. math::
        B_{simple}^o = 0.1445 - \frac{0.33}{T_r} - \frac{0.1385}{T_r^2} -
        \frac{0.0121}{T_r^3} - \frac{0.000607}{T_r^8}

    .. math::
        B_{size-shape}^o = \left(-0.00787 + \frac{0.0812}{T_R^2} -
        \frac{0.0646}{T_R^3}\right)R - \left(\frac{0.00347}{T_R^2} -
        \frac{0.000149}{T_R^7}\right)R^2

    .. math::
        B_{polar}^o = -\frac{0.028}{T_r^7}\Psi

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    id : integer
        Index of compound in database

    Returns
    -------
    B : float
        Second virial coefficient [dm³/mol]
    Bt : float
        T(∂B/∂T) [dm³/molK]
    Btt : float
        T²(∂²B/∂T²) [dm³/molK²]

    Examples
    --------
    Using same value of Tsonopoulos correlation as reference

    >>> from lib.mEoS import Ne
    >>> "%0.4f" % B_Tarakad(262, Ne.Tc, Ne.Pc, Ne.id)[0]
    '0.0117'
    """
    Tr = T/Tc

    par = {
        65: (1.094, 0.439),
        66: (1.886, 0.133),
        67: (2.713, 0),
        40: (3.004, 0),
        41: (3.443, 0.105),
        45: (3.821, 0.101),
        42: (3.789, 0.229),
        43: (3.897, 0.329),
        44: (3.796, 0.308),
        71: (4.187, 0.199),
        72: (4.130, 0.049),
        73: (4.285, 0),
        74: (4.166, 0.066),
        75: (4.100, 0.179),
        76: (4.168, 0.152),
        77: (4.341, 0.127),
        225: (1.419, 0.867),
        643: (2.319, 0.250),
        218: (2.635, 0),
        247: (2.176, 0.359),
        245: (2.498, 0.425),
        243: (2.672, 0),
        235: (3.380, 0),
        691: (4.443, 0),
        321: (5.521, 0),
        326: (6.022, 0),
        320: (4.683, 0),
        236: (3.150, 0),
        322: (3.345, 0.04),
        319: (4.658, 0.323),
        642: (2.880, 0),
        220: (2.569, 0.119),
        216: (3.026, 0),
        215: (2.810, 0),
        217: (3.304, 0),
        115: (1.450, 0.436),
        222: (2.342, 0.050),
        112: (3.178, 0),
        100: (3.458, 0),
        132: (2.281, 0.396),
        126: (2.995, 0.151),
        127: (2.851, 0.071),
        479: (3.357, 0.229),
        122: (2.122, 0),
        658: (3.013, 0.195),
        659: (3.013, 0.041),
        119: (3.759, 0.061),
        224: (1.180, 0.348),
        246: (1.988, 0),
        116: (1.046, 0.118),
        133: (2.127, 0.460),
        486: (2.641, 0.255),
        162: (3.139, 0.157),
        318: (3.547, 0.129),
        337: (3.681, 0.154),
        280: (2.990, 0.200),
        129: (1.901, 0.422),
        117: (1.536, 2.146),
        134: (2.250, 2.226),
        146: (2.736, 1.431),
        62: (0.615, 1.220),
        145: (2.726, 1.696),
        159: (3.182, 1.530),
        161: (3.019, 1.749),
        160: (3.225, 0.611),
        450: (3.151, 1.349),
        174: (3.550, 0.634),
        140: (2.740, 0.966),
        153: (3.139, 0.557),
        304: (3.627, 0.797),
        305: (3.415, 0.421),
        165: (3.482, 0.558),
        329: (3.572, 0.925),
        131: (2.360, 0.683),
        141: (2.870, 0.335),
        157: (3.419, 0.260),
        142: (2.862, 0.630),
        155: (3.348, 0.454),
        166: (3.778, 0.425),
        156: (3.304, 0.454),
        309: (3.712, 0.386),
        136: (2.372, 0.077),
        269: (2.809, 0.559),
        290: (3.207, 0.139),
        # Methyl propyl sulfide: 3.317, 0.11,
        # Methyl isopropyl sulfide: 3.236, 0,
        227: (1.611, 0.197),
        137: (2.341, 0.094),
        # 2-methyl-2-propanethiol: 3.110, 0,
        # 1-pentanethiol: 3.791, 0.097,
        113: (0.650, 2.065),
        125: (1.821, 3.008),
        118: (1.662, 0.903),
        249: (2.264, 0.552),
        147: (2.736, 0.153),
        138: (2.309, 0.660),
        294: (3.161, 0.197),
        270: (2.789, 0.452),
        226: (2.306, 1.398),
        108: (0.530, 0.454),
        110: (1.191, 0.231),
        48: (0.558, 0),
        49: (0.992, 0.152),
        51: (1.674, 0.512),
        50: (0.638, 0.301),
        102: (1.424, 0.074),
        63: (0.853, 0.917),
        630: (2.356, 0),
        629: (3.279, 0),
        237: (1.423, 0)}

    if id and id in par:
        R_, fi = par[id]
    else:
        R_, fi = 0, 0

    def f(T):
        Bsimple = .1445 - 0.33/Tr - .1385/Tr**2 - .0121/Tr**3 - .000607/Tr**8
        Bshape = (-0.00787 + 0.0812/Tr**2 - 0.0646/Tr**3)*R_ - \
            (0.00347/Tr**2-0.000149/Tr**7)*R_**2
        Bpolar = -0.028/Tr**7*fi

        B = Bsimple + Bshape + Bpolar
        return B*R*Tc/Pc*1e3

    B = f(T)
    B1 = derivative(f, T, n=1)
    B2 = derivative(f, T, n=2)
    return B, B1, B2


@refDoc(__doi__, [4])
def C_OrbeyVera(T, Tc, Pc, w):
    r"""Calculate the third virial coefficient using the Orbey-Vera correlation

    .. math::
        \frac{CP_c^2}{R^2T_c^2} = f^{(0)} + f^{(1)}\omega

    .. math::
        f^{(0)} = 0.01407+\frac{0.02432}{T_r^{2.8}}-\frac{0.00313}{T_r^{10.5}}

    .. math::
        f^{(1)} = -0.02676+\frac{0.0177}{T_r^{2.8}}+\frac{0.04}{T_r^3}
        -\frac{0.003}{T_r^6}-\frac{0.00228}{T_r^{10.5}}

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    w : float
        Acentric factor [-]

    Returns
    -------
    C : float
        Third virial coefficient [dm⁶/mol²]
    C1 : float
        T(∂C/∂T) [dm⁶/mol²]
    C2 : float
        T²(∂²C/∂T²) [dm⁶/mol²]

    Examples
    --------
    Selected points from Table 2 of paper

    >>> from lib.mEoS.Benzene import Benzene as Bz
    >>> "%.1f" % (C_OrbeyVera(0.877*Bz.Tc, Bz.Tc, Bz.Pc, Bz.f_acent)[0]*1e9)
    '41.7'
    >>> "%.1f" % (C_OrbeyVera(1.019*Bz.Tc, Bz.Tc, Bz.Pc, Bz.f_acent)[0]*1e9)
    '36.0'
    """
    def f(T):
        Tr = T/Tc
        g0 = 0.01407+0.02432/Tr**2.8-0.00313/Tr**10.5
        g1 = -0.02676+0.0177/Tr**2.8+0.04/Tr**3-0.003/Tr**6-0.00228/Tr**10.5
        g = g0+w*g1
        return g*R**2*Tc**2/Pc**2

    C = f(T)
    C1 = derivative(f, T, n=1)
    C2 = derivative(f, T, n=2)

    return C, C1, C2


@refDoc(__doi__, [5])
def C_LiuXiang(T, Tc, Pc, w, Zc):
    r"""Calculate the third virial coefficient using the Liu-Xiang correlation

    .. math::
        \frac{CP_c^2}{R^2T_c^2} = f^{(0)} + \omega f^{(1)} + \theta f^{(2)}

    .. math::
        \theta = \left(Z_c-0.29\right)^2

    .. math::
        f^{(0)} = a_{00} + \frac{a_{10}}{T_r^3} + \frac{a_20}{T_r^6} +
        \frac{a_{30}}{T_r^{11}}

    .. math::
        f^{(1)} = a_{01} + \frac{a_{11}}{T_r^3} + \frac{a_21}{T_r^6} +
        \frac{a_{31}}{T_r^{11}}

    .. math::
        f^{(2)} = a_{02} + \frac{a_{12}}{T_r^3} + \frac{a_22}{T_r^6} +
        \frac{a_{32}}{T_r^{11}}

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    w : float
        Acentric factor [-]
    Zc : float
        Critical compresibility factor [-]

    Returns
    -------
    C : float
        Third virial coefficient [dm⁶/mol²]
    C1 : float
        T(∂C/∂T) [dm⁶/mol²]
    C2 : float
        T²(∂²C/∂T²) [dm⁶/mol²]
    """
    X = (Zc-0.29)**2

    def f(T):
        Tr = T/Tc
        g0 = 0.1623538 + 0.3087440/Tr**3 - 0.01790184/Tr**6 - 0.02789157/Tr**11
        g1 = -0.5390344 + 1.783526/Tr**3 - 1.055391/Tr**6 + 0.09955867/Tr**11
        g2 = 34.22804 - 74.76559/Tr**3 + 279.9220/Tr**6 - 62.85431/Tr**11
        g = g0+w*g1+X*g2
        return g*R**2*Tc**2/Pc**2*1e-1

    C = f(T)
    C1 = derivative(f, T, n=1)
    C2 = derivative(f, T, n=2)

    return C, C1, C2


@refDoc(__doi__, [7])
def C_Meng(T, Tc, Pc, D, B):
    r"""Calculate the 3rd virial coefficient using the Meng-Duan-Li correlation

    .. math::
        \frac{CP_c^2}{R^2T_c^2} = 5.476e-3 + \left(\frac{BP_c}{RT_c}-0.0936
        \right)^2 \left(f^{(0)} + \mu_r^4 f^{(1)}x10^{-10}\right)

    .. math::
        f^{(0)} = 1094.051 - \frac{3334.145}{T_r^{0.1}} +
        \frac{3389.848}{T_r^{0.2}} - \frac{1149.58}{T_r^{0.3}}

    .. math::
        f^{(1)} = 2.0243 - \frac{0.85902}{T_r}

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure
    D : float
        dipole moment [debye]
    B : list
        Second virial coefficient tuple with B, ∂B/∂T, ∂²B/∂T²

    Returns
    -------
    C : float
        Third virial coefficient [dm⁶/mol²]
    C1 : float
        T(∂C/∂T) [dm⁶/mol²]
    C2 : float
        T²(∂²C/∂T²) [dm⁶/mol²]
    """
    mur = D**2*Pc/1.01325/Tc**2

    def f(T, *args):
        B = args[-1]*1e-3
        Br = B*Pc/R/Tc
        Tr = T/Tc
        f0 = 1094.051 - 3334.145/Tr**0.1 + 3389.848/Tr**0.2 - 1149.58/Tr**0.3
        f1 = (2.0243-0.85902/Tr)*1e-10
        Cr = (5.476e-3 + (Br-0.0936)**2*(f0+mur**4*f1))
        return Cr/Pc**2*R**2*Tc**2

    C = f(T, B[0])
    C1 = derivative(f, T, n=1, args=(T, B[1]))
    C2 = derivative(f, T, n=2, args=(T, B[2]))

    return C, C1, C2


class Virial(EoS):
    r"""Class to model the virial equation of state

    .. math::
        Z = \frac{PV}{RT} = 1 + \frac{B}{V} + \frac{C}{V^2} + \cdots

    The implementation use the form truncated at third term using the virial
    coefficient from database or try to predicted from reference correlations.
    This equation is only appropiate for single-phase gas systems.
    """
    __title__ = "Virial"
    __status__ = "Virial"

    METHODS_B = ["Tsonopoulos (1974)", "Iglesias-Silva (2001)", "Meng (2004)"]
    METHODS_C = ["Orbey-Vera (1983)", "Liu-Xiang (2003)", "Meng (2004)"]

    def __init__(self, *args, **kwargs):
        EoS.__init__(self, *args, **kwargs)
        self._physics(*args)
        self.x, self.xi, self.yi, self.Ki = self._Flash()

    def _Bi(self, T):
        """Second virial coefficient muxture contributions"""
        B = []
        Bt = []
        Btt = []
        for cmp in self.componente:
            if cmp.id in B_Database:
                if cmp.id == 1:  # Hydrogen special case
                    if self.T < 60:
                        coef = [2.0375e1, -2-2113e3, -2.0892e4, -6.5299e4]
                    else:
                        coef = [1.7472e1, 1.2926e2, -2.6988e5, 8.0282e6]
                elif cmp.id == 212:   # Helium special case
                    if self.T < 35.1:
                        coef = [1.5943e1, -3.4601e2, -5.9545e2, 1.9929e3,
                                2.2269e3]
                    else:
                        coef = [9.2479, 1.0876e3, -1.088e5, 2.3869e6]
                else:
                    coef = B_Database[cmp.id]
                Bi, Bit, Bitt = _B_Database(T, coef)
            else:
                # Use general correlation
                Bi, Bit, Bitt = B_Tsonopoulos(
                    T, cmp.Tc, cmp.Pc, cmp.f_acent, cmp.dipole)
            B.append(Bi)
            Bt.append(Bit)
            Btt.append(Bitt)
        return B, Bt, Btt

    def _Ci(self, T):
        """Third virial coefficient mixture contributions"""
        C = []
        Ct = []
        Ctt = []
        for comp in self.componente:
            if self.kwargs.get("C", 0):
                Ci, Cit, Citt = C_OrbeyVera(T, comp.Tc, comp.Pc, comp.f_acent)
            else:
                Ci, Cit, Citt = C_LiuXiang(
                    T, comp.Tc, comp.Pc, comp.f_acent, comp.Zc)
            C.append(Ci)
            Ct.append(Cit)
            Ctt.append(Citt)
        return C, Ct, Ctt

    def B(self, T):
        """Second virial coefficient calculation"""
        Bi, Bit, Bitt = self._Bi(T)
        B, Bt, Btt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                if i == j:
                    Bij = Bi[i]
                    Bijt = Bit[i]
                    Bijtt = Bitt[i]
                elif "%i-%i" % (i, j) in Bij_Database:
                    id = "%i-%i" % (i, j)
                    coef = Bij_Database[id]
                    Bij, Bijt, Bijtt = _B_Database(T, coef)
                elif "%i-%i" % (j, i) in Bij_Database:
                    id = "%i-%i" % (j, i)
                    coef = Bij_Database[id]
                    Bij, Bijt, Bijtt = _B_Database(T, coef)
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/(ci.Vc**(1/3)+cj.Vc**(1/3))**3
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    Bij, Bijt, Bijtt = B_Tsonopoulos(self.T, Tcij, Pcij, wij)
                B += xi*xj*Bij
                Bt += xi*xj*Bijt
                Btt += xi*xj*Bijtt
        return B, Bt, Btt

    def C(self, T):
        """Third virial coefficient calculation"""
        Ci, Cit, Citt = self._Ci(T)
        Cij = zeros((len(Ci), len(Ci)))
        Cijt = zeros((len(Ci), len(Ci)))
        Cijtt = zeros((len(Ci), len(Ci)))
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                if i == j:
                    Cij[i, j] = Ci[i]
                    Cijt[i, j] = Cit[i]
                    Cijtt[i, j] = Citt[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Vcij = (ci.Vc**(1./3)+cj.Vc**(1./3))**3
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/Vcij
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    if self.kwargs.get("C", 0):
                        Cij[i, j] = C_OrbeyVera(self.T, Tcij, Pcij, wij)
                    else:
                        Zcij = Pcij*Vcij/R/Tcij
                        Cij[i, j] = C_LiuXiang(
                            self.T, Tcij, Pcij, wij, Zcij)[0]

        C, Ct, Ctt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                for k, xk in enumerate(self.fraccion):
                    C += xi*xj*xk*(Cij[i, j]*Cij[j, k]*Cij[j, k])**(1./3)
        return C, Ct, Ctt

    def _physics(self, T, P, mezcla):
        """Properties of Gases calculation. Explanation in [1]_ section 1.4"""
        B, B1, B2 = self.B(T)
        C, C1, C2 = self.C(T)

        self.Z = 1+B*(P/R/T)+(C-B**2)*(P/R/T)**2
        V = self.Z*R*T/P
        self.U_exc = -R*T*(B1/V+C1/2/V**2)
        self.H_exc = R*T*((B-B1)/V+(2*C-C1)/2/V**2)
        self.Cv_exc = -R*((2*B1+B2)/V+(2*C1+C2)/2/V**2)
        self.Cp_exc = -R*(B2/V-((B-B1)**2-(C-C1)-C2/2)/V**2)
        self.S_exc = -R*(log(P)+B1/V+(B**2-C+C1)/2/V**2)
        self.A_exc = R*T*(log(P)+(B**2-C/2/V**2))
        self.G_exc = R*T*(log(P)+B/V+(B**2+C)/2/V**2)

        self.fug = P*exp(B/V+(C+B**2)/2/V**2)

    def _fug(self, Z, xi):
        """Calculate partial fugacities coefficieint of components, Eq 6-7.9"""
        pass


_all = [Virial]


if __name__ == "__main__":
    from lib.mezcla import Mezcla

    mezcla = Mezcla(2, ids=[10, 38, 22, 61],
                    caudalUnitarioMolar=[0.3, 0.5, 0.05, 0.15])
    eq = Virial(340, 1., mezcla)
    print(dir(eq))
    print(eq.x)
