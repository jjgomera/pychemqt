#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Virial equation of state implementation
###############################################################################

from scipy import zeros, log, exp
from scipy.constants import R

from lib.eos import EoS


__doi__ = {
    1:
        {"autor": "Dymond, J.H., Marsh, K.N., Wilhoit, R.C., Wong, K.C., "
                  "Frenkel, M.,",
         "title": "Virial Coefficients of Pure Gases (Landolt-Börnstein - "
                  "Group IV Physical Chemistry 21A)",
         "ref": "Springer-Verlag",
         "doi": "10.1007/10693952_16"},
    2:
        {"autor": "Tsonopoulos, C.",
         "title": "An empirical correlation of second virial coefficients",
         "ref": "AICHE Journal 20, pp 263 (1974)",
         "doi": "10.1002/aic.690200209"},
    3:
        {"autor": "Orbey, H., Vera, J.H.",
         "title": "Correlation for the third virial coefficient using Tc, Pc "
                  "and ω as parameters",
         "ref": "AIChE Journal 29, 107 (1983)",
         "doi": "10.1002/aic.690290115"},
    4:
        {"autor": "Liu, D.X., Xiang, H.W.",
         "title": "Corresponding-States Correlation and Prediction of Third "
                  "Virial Coefficients for a Wide Range of Substances",
         "ref": "Int. J. Thermophysics 24(6), 1667-1680 (2003)",
         "doi": "10.1023/b:ijot.0000004098.98614.38"},
}


B_Database = {
    98: [3.4162e1, -1.2087e4, -7.6702e5, -1.96e7],
    630: [9.1039e1, -5.9081e4, 1.0478e7, 3.0463e9],
    113: [-6.0611e3, 4.6389e6, 9.8352e8],
    48: [4.8826e1, -1.5614e4, -2.757e5, -4.7684e7],
    49: [5.74e1, -3.8829e4, 4.2899e5, -1.4661e9],
    102: [2.8415e3, -2.8985e6, 9.583e8, -1.2515e11],
    215: [-1.0754e3, 9.3227e5, -2.7872e8],
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
    225: [1.6048e3, -1.751e6, 6.1379e8, -8.4605e10],
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
    12: [2.5227e3, -3.8225e6, 1.8393e9, -3.9652e11],
    # Octamethyltrisilosane pag 247
}


def B_Tsonopoulos(T, Tc, Pc, w, mu=None):
    """Calculate the 2nd virial coefficient using the Tsonopoulos correlation

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
        Second virial coefficient [m³/mol]
    B1 : float
        T(dB/dT) [m³/mol]
    B2 : float
        T²(d²B/dT²) [m³/mol]

    Notes
    -----
    With the B_database this correlation is only for completeness, the a, b
    correlations are general and possibly not applicable to the compounds
    not availables in the B_database

    References
    ----------
    .. [2] Tsonopoulos, C. An empirical correlation of second virial
        coefficients. AICHE Journal 20, pp 263 (1974)
    """
    Tr = T/Tc
    if mu:
        mur = mu**2*Pc/1.01325/Tc**2
        a = -2.14e-4*mur-7.831e-21*mur**8
        b = 0.00908+0.0006957*mur
    else:
        a, b = 0, 0

    f0 = 0.1445-0.33/Tr-0.1385/Tr**2-0.0121/Tr**3-0.000607/Tr**8
    f1 = 0.0637+0.331/Tr**2-0.423/Tr**3-0.008/Tr**8
    f2 = 1/Tr**6
    f3 = -1/Tr**8
    f = f0 + w*f1 + a*f2 + b*f3
    B = f*R*Tc/Pc
    f0t = 0.33*Tc/T**2 + 2*0.1385*Tc**2/T**3 + 3*0.0121*Tc**3/Tr**4 + \
        8*0.000607*Tc**8/T**9
    f1t = -2*0.331*Tc**2/T**3+3*0.423*Tc**3/T**4+8*0.008+Tc**8/T**9
    f2t = -6*Tc**6/T**7
    f3t = 8*Tc**8/T**9
    ft = f0t + w*f1t + a*f2t + b*f3t
    B1 = ft*R*Tc/Pc
    f0tt = -2*0.33*Tc/T**3 - 2*3*0.1385*Tc**2/T**4 - 3*4*0.0121*Tc**3/Tr**5 - \
        8*9*0.000607*Tc**8/T**10
    f1tt = -2*0.331*Tc**2/T**3+3*0.423*Tc**3/T**4+8*0.008+Tc**8/T**9
    f2tt = 6*7*Tc**6/T**8
    f3tt = -8*9*Tc**8/T**10
    ftt = f0tt + w*f1tt + a*f2tt + b*f3tt
    B2 = ftt*R*Tc/Pc

    return B, B1, B2


def C_Orbey_Vera(T, Tc, Pc, w):
    """Calculate the third virial coefficient using the Orbey-Vera correlation

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
        Third virial coefficient [m⁶/mol²]
    C1 : float
        T(dC/dT) [m⁶/mol²]
    C2 : float
        T²(d²C/dT²) [m⁶/mol²]

    Examples
    --------
    Selected points from Table 2 of paper
    >>> from lib.mEoS.Benzene import Benzene as Bz
    >>> "%.1f" % (C_Orbey_Vera(0.877*Bz.Tc, Bz.Tc, Bz.Pc, Bz.f_acent)[0]*1e9)
    '41.5'
    >>> "%.1f" % (C_Orbey_Vera(1.019*Bz.Tc, Bz.Tc, Bz.Pc, Bz.f_acent)[0]*1e9)
    '35.8'

    References
    ----------
    [3] .. Orbey, H., Vera, J.H.: Correlation for the third virial coefficient
        using Tc, Pc and ω as parameters, AIChE Journal 29, 107 (1983)
    """
    Tr = T/Tc
    g0 = 0.01407+0.02432/Tr**2.8-0.00313/Tr**10.5
    g1 = -0.02676+0.0177/Tr**2.8+0.04/Tr**3-0.003/Tr**6-0.00228/Tr**10.5
    g = g0+w*g1
    C = g*R**2*Tc**2/Pc**2
    g0t = -2.8*0.02432*Tc**2.8/T**3.8 + 10.5*0.00313*Tc**10.5/T**11.5
    g1t = -2.8*0.0177*Tc**2.8/T**3.8 - 3*0.04*Tc**3/Tr**4 - \
        0.003*Tc**6/T**7 + 10.5*0.00228*Tc**10.5/T**11.5
    gt = g0t+w*g1t
    C1 = gt*R**2*Tc**2/Pc**2*T
    g0tt = 2.8*3.8*0.02432*Tc**2.8/T**4.8 - 10.5*11.5*0.00313*Tc**10.5/T**12.5
    g1tt = 2.8*3.8*0.0177*Tc**2.8/T**4.8 + 3*4*0.04*Tc**3/Tr**5 - \
        0.003*Tc**6/T**8 - 10.5*11.5*0.00228*Tc**10.5/T**12.5
    gtt = g0tt+w*g1tt
    C2 = gtt*R**2*Tc**2/Pc**2*T**2

    return C, C1, C2


def C_Liu_Xiang(T, Tc, Pc, w, Zc):
    """Calculate the third virial coefficient using the Liu-Xiang correlation

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
        Third virial coefficient [m⁶/mol²]
    C1 : float
        T(dC/dT) [m⁶/mol²]
    C2 : float
        T²(d²C/dT²) [m⁶/mol²]

    References
    ----------
    [4] .. Liu, D.X., Xiang, H.W.: Corresponding-States Correlation and
        Prediction of Third Virial Coefficients for a Wide Range of Substances.
        International Journal of Thermophysics, November 2003, Volume 24,
        Issue 6, pp 1667-1680
    """
    Tr = T/Tc
    X = (Zc-0.29)**2

    g0 = 0.1623538 + 0.3087440/Tr**3 - 0.01790184/Tr**6 - 0.02789157/Tr**11
    g1 = -0.5390344 + 1.783526/Tr**3 - 1.055391/Tr**6 + 0.09955867/Tr**11
    g2 = 34.22804 - 74.76559/Tr**3 + 279.9220/Tr**6 - 62.85431/Tr**11
    g = g0+w*g1+X*g2
    C = g*R**2*Tc**2/Pc**2

    g0t = -3*0.3087440*Tc**3/T**4 + 6*0.01790184*Tc**6/T**7 + \
        11*0.02789157*Tc**11/T**12
    g1t = -3*1.783526*Tc**3/T**4 + 6*1.055391*Tc**6/T**7 - \
        11*0.09955867*Tc**11/T**12
    g2t = 3*74.76559*Tc**3/Tr**4 - 6*279.9220*Tc**6/T**7 + \
        11*62.85431*Tc**11/T**12
    gt = g0t+w*g1t+X*g2t
    C1 = gt*R**2*Tc**2/Pc**2*T

    g0tt = 3*4*0.3087440*Tc**3/T**5 - 6*7*0.01790184*Tc**6/T**8 - \
        11*12*0.02789157*Tc**11/T**13
    g1tt = 3*4*1.783526*Tc**3/T**5 - 6*7*1.055391*Tc**6/T**8 + \
        11*12*0.09955867*Tc**11/T**13
    g2tt = -3*4*74.76559*Tc**3/Tr**5 + 6*7*279.9220*Tc**6/T**8 - \
        11*12*62.85431*Tc**11/T**13
    gtt = g0tt+w*g1tt+X*g2tt
    C2 = gtt*R**2*Tc**2/Pc**2*T**2

    return C, C1, C2


class Virial(EoS):
    """Class to model virial equation of state"""

    def __init__(self, *args, **kwargs):
        EoS.__init__(self, *args, **kwargs)
        self._physics(*args)

    def _Bi(self):
        """Second virial coefficient muxture contributions"""
        B = []
        Bt = []
        Btt = []
        for comp in self.componente:
            if comp.indice in B_Database:
                Bi = 0
                Bit = 0
                Bitt = 0
                if comp.indice == 1:  # Hydrogen special case
                    if self.T < 60:
                        coef = [2.0375e1, -2-2113e3, -2.0892e4, -6.5299e4]
                    else:
                        coef = [1.7472e1, 1.2926e2, -2.6988e5, 8.0282e6]
                elif comp.indice == 212:   # Helium special case
                    if self.T < 35.1:
                        coef = [1.5943e1, -3.4601e2, -5.9545e2, 1.9929e3,
                                2.2269e3]
                    else:
                        coef = [9.2479, 1.0876e3, -1.088e5, 2.3869e6]
                else:
                    coef = B_Database[comp.indice]
                for i, a in enumerate(coef):
                    Bi += a/self.T**i
                    Bit += -a*i*self.T**(i-1)
                    Bitt += a*i*(i-1)*self.T**(i-2)
            else:
                # Use general correlation
                Bi, Bit, Bitt = B_Tsonopoulos(
                        self.T, comp.Tc, comp.Pc, comp.f_acent,
                        comp.dipole)
        B.append(Bi)
        Bt.append(Bit)
        Btt.append(Bitt)
        return B, Bt, Btt

    def _Ci(self):
        """Third virial coefficient muxture contributions"""
        C = []
        Ct = []
        Ctt = []
        for comp in self.componente:
            if self.kwargs.get("C", 0):
                Ci, Cit, Citt = C_Orbey_Vera(
                    self.T, comp.Tc, comp.Pc, comp.f_acent)
            else:
                Ci, Cit, Citt = C_Liu_Xiang(
                    self.T, comp.Tc, comp.Pc, comp.f_acent, comp.Zc)
            C.append(Ci)
            Ct.append(Cit)
            Ctt.append(Citt)
        return C, Ct, Ctt

    def B(self):
        """Second virial coefficient calculation"""
        Bi, Bit, Bitt = self._Bi()
        B, Bt, Btt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                if i == j:
                    Bij = Bi[i]
                    Bijt = Bit[i]
                    Bijtt = Bitt[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/(ci.Vc**(1/3)+cj.Vc**(1/3))**3
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    Bij, Bijt, Bijtt = self._B_Tsonopoulos(Tcij, Pcij, wij)
                B += xi*xj*Bij
                Bt += xi*xj*Bijt
                Btt += xi*xj*Bijtt
        return B, Bt, Btt

    def C(self):
        """Third virial coefficient calculation"""
        Ci, Cit, Citt = self._Ci()
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
                        Cij[i, j] = self._C_Orbey_Vera(Tcij, Pcij, wij)
                    else:
                        Zcij = Pcij*Vcij/R/Tcij
                        Cij[i, j] = self._C_Liu_Xiang(Tcij, Pcij, wij, Zcij)

        C, Ct, Ctt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                for k, xk in enumerate(self.fraccion):
                    C += xi*xj*xk*(Cij[i, j]*Cij[j, k]*Cij[j, k])**(1./3)
        return C, Ct, Ctt

    def _physics(self, T, P, mezcla):
        """Properties of Gases calculation. Explanation in [1]_ section 1.4"""
        B, B1, B2 = self.B()
        C, C1, C2 = self.C()

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
