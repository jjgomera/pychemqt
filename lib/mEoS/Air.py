#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


from lib.meos import MEoSBlend
from lib import unidades


class Air(MEoSBlend):
    """Multiparameter equation of state for Air as pseudocomponent"""
    name = "air"
    CASNumber = "1"
    formula = "N2+Ar+O2"
    synonym = "R-729"
    rhoc = unidades.Density(342.60456)
    Tc = unidades.Temperature(132.6306)
    Pc = unidades.Pressure(3786.0, "kPa")
    M = 28.96546  # g/mol
    Tt = unidades.Temperature(59.75)
    Tb = unidades.Temperature(78.903)
    f_acent = 0.0335
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 475

    Fi1 = {"ao_log": [1, 2.490888032],
           "pow": [-3, -2, -1, 0, 1, 1.5],
           "ao_pow": [0.6057194e-7, -0.210274769e-4, -0.158860716e-3,
                      -13.841928076, 17.275266575, -0.19536342e-3],
           "ao_exp": [0.791309509, 0.212236768],
           "titao": [25.36365, 16.90741],
           "ao_exp2": [-0.197938904],
           "titao2": [87.31279],
           "sum2": [2./3]
           }

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Lemmon et al. (2000)",
        "__doi__": {"autor": "Lemmon, E.W., Jacobsen, R.T, Penoncello, S.G., and Friend, D.G.",
                    "title": "Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon, and Oxygen From 60 to 2000 K at Pressures to 2000 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29, 331 (2000)",
                    "doi":  "10.1063/1.1285884"},
        "__test__":
            # Table A1, Pag 363
            """
            >>> print "%0.6f %0.5f" % (Air._bubbleP(59.75).MPa, Air._dewP(59.75).MPa)
            0.005265 0.00243
            >>> print "%0.5f %0.5f" % (Air._bubbleP(70).MPa, Air._dewP(70).MPa)
            0.03191 0.01943
            >>> print "%0.5f %0.5f" % (Air._bubbleP(80).MPa, Air._dewP(80).MPa)
            0.11462 0.08232
            >>> print "%0.5f %0.5f" % (Air._bubbleP(100).MPa, Air._dewP(100).MPa)
            0.66313 0.56742
            >>> print "%0.5f %0.5f" % (Air._bubbleP(120).MPa, Air._dewP(120).MPa)
            2.15573 2.00674
            >>> print "%0.5f %0.5f" % (Air._bubbleP(130).MPa, Air._dewP(130).MPa)
            3.42947 3.30835
            """
            # Table A2, Pag 366
            """
            >>> st=Air(T=100, P=101325)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            100 0.12449 2028.2 2784.1 166.61 21.09 30.13 198.2
            >>> st=Air(T=2000, P=101325)
            >>> print "%0.0f %0.4g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            2000 0.006092 48610 65242 259.62 27.90 36.21 863.5
            >>> st=Air(T=500, P=2e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            500 0.048077 10418 14578 208.2 21.51 29.84 446.6
            >>> st=Air(T=300, P=5e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 0.20075 6179.8 8670.5 185.5 20.82 29.33 347.8
            >>> st=Air(T=130, P=1e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            130 1.0295 2461.1 3432.5 153.79 22.058 34.69 216.8
            >>> st=Air(T=70, P=5e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            70 31.895 -4198 -4041.2 78.907 32.17 54.57 974.6
            >>> st=Air(T=2000, P=1e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            2000 0.59094 48600 65522 221.44 27.93 36.25 878.6
            >>> st=Air(T=130, P=5e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            130 27.946 -1831.3 -42.096 104.84 27.68 48.19 878.8
            >>> st=Air(T=100, P=1e8)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            100 33.161 -3403.9 -388.34 87.644 31.98 48.22 1192.4
            >>> st=Air(T=1000, P=1e9)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.2f %0.2f %0.1f" % (\
                st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            1000 30.791 21944 54421 156.83 29.07 36.77 1966.3
            """,

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 8649.34, "so": 194.},

        "M": 28.9586, "Tc": 132.6312, "rhoc": 10.4477,

        "Tmin": Tt, "Tmax": 2000., "Pmax": 2000000.0, "rhomax": 53.73,
        "Pmin": 5.2646, "rhomin": 33.067,

        "Tj": 132.6312, "Pj": 3.78502,
        "dew": {"i": [1, 2, 5, 8],
                "n": [-0.1567266, -5.539635, 0.7567212, -3.514322]},
        "bubble": {"i": [1, 2, 3, 4, 5, 6],
                   "n": [0.2260724, -7.080499, 5.700283, -12.44017, 17.81926,
                         -10.81364]},

        "nr1": [0.118160747229, 0.713116392079, -0.161824192067e1,
                0.714140178971e-1, -0.865421396646e-1, 0.134211176704,
                0.112626704218e-1, -0.420533228842e-1, 0.349008431982e-1,
                0.164957183186e-3],
        "d1": [1, 1, 1, 2, 3, 3, 4, 4, 4, 6],
        "t1": [0, 0.33, 1.01, 0, 0, 0.15, 0, 0.2, 0.35, 1.35],

        "nr2": [-0.101365037912, -0.173813690970, -0.472103183731e-1,
                -0.122523554253e-1, -0.146629609713, -0.316055879821e-1,
                0.233594806142e-3, 0.148287891978e-1, -0.938782884667e-2],
        "d2": [1, 3, 5, 6, 1, 3, 11, 1, 3],
        "t2": [1.6, 0.8, 0.95, 1.25, 3.6, 6, 3.25, 3.5, 15],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3, 3],
        "gamma2": [1]*9}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Jacobsen et al. (1992)",
        "__doi__": {"autor": "Jacobsen, R.T, Penoncello, S.G., Beyerlein, S.W., Clarke, W.P., and Lemmon, E.W.",
                    "title": "A Thermodynamic Property Formulation for Air",
                    "ref": "Fluid Phase Equilibria, 79:113-124, 1992.",
                    "doi":  "10.1016/0378-3812(92)85124-Q"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 8649.34, "so": 194.},
        "M": 28.9586, "Tc": 132.6312, "rhoc": 10.4477,

        "Tmin": Tt, "Tmax": 870.0, "Pmax": 70000.0, "rhomax": 34.628,
        "Pmin": 6.2545, "rhomin": 33.073,

        "Tj": 132.61738, "Pj": 3.78502,
        "dew": {"i": [1, 2, 10, 11, 13, 14],
                "n": [-0.1537763029, -5.544542064, 312.7182733, -895.9553274,
                      1834.176566, -1321.892808]},
        "bubble": {"i": [1, 2, 4, 5, 6, 7, 12],
                   "n": [0.2095592444, -6.654905539, 22.13718815, -84.14553609,
                         135.9753732, -83.66895082, 17.97856602]},

        "nr1": [0.206604930965, 0.367099749382, -0.943192015369,
                0.382519513142e-2, -0.865385542309e-1, 0.323019987452,
                0.608695449299e-2, 0.128352106296e-3, -0.400058181940e-5],
        "d1": [1, 1, 1, 1, 2, 2, 4, 6, 7],
        "t1": [0, 0.25, 1, 3.5, 0, 0.25, 0.5, 2, 3],

        "nr2": [-0.544697915817, -0.526471065792, -0.608529300347,
                -0.124174277875, -0.595578533411e-2, -0.157523548353,
                -0.346463251040e-2, 0.837023084176e-2, -0.316701981142e-1,
                -0.721856676857e-2, 0.276838040645e-3, 0.160877459321e-4,
                0.409235806738e-1, 0.652776125216e-3, -0.952903961290e-2,
                -0.100337820004e-1, 0.701111041628e-2, -0.472754336912e-2,
                0.399257638569e-2, 0.968453675994e-2, -0.106826283630e-1,
                -0.489679885832e-2],
        "d2": [1, 2, 3, 5, 6, 1, 1, 2, 2, 3, 11, 11, 1, 1, 2, 3, 7, 8, 2, 4, 5, 2],
        "t2": [1.5, 1, 1, 1, 2, 3, 8, 0.5, 5.5, 9, 3, 6, 3, 9, 2, 13, 11, 11, 8,
               22, 23, 11],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5],
        "gamma2": [1]*22}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.03046], "exp": [1.28]}
    _melting = {"eq": 1, "Tref": Tb, "Pref": 5.265,
                "Tmin": 59.75, "Tmax": 2000.0,
                "a1": [1, 0.354935e5, -0.354935e5],
                "exp1": [0, 0.178963e1, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.20466e1, -0.4752e1, -0.13259e2, -0.47652e2],
        "exp": [0.41, 1, 2.8, 6.5]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.1567266, -0.5539635e1, 0.7567212, -0.3514322e1],
        "exp": [0.5, 1, 2.5, 4]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Lemmon (2004)",
              "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                          "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air",
                          "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                          "doi": "10.1023/B:IJOT.0000022327.04529.f3"},
              "__test__": """
                  >>> st=Air(T=100, rhom=0)
                  >>> print "%0.5f" % st.mu.muPas
                  7.09559
                  >>> st=Air(T=300, rhom=0)
                  >>> print "%0.4f" % st.mu.muPas
                  18.5230
                  >>> st=Air(T=100, rhom=28)
                  >>> print "%0.3f" % st.mu.muPas
                  107.923
                  >>> st=Air(T=200, rhom=10)
                  >>> print "%0.4f" % st.mu.muPas
                  21.1392
                  >>> st=Air(T=300, rhom=5)
                  >>> print "%0.4f" % st.mu.muPas
                  21.3241
                  >>> st=Air(T=132.64, rhom=10.4)
                  >>> print "%0.4f" % st.mu.muPas
                  17.7623
                  """, # Table V, Pag 28

            "ek": 103.3, "sigma": 0.36,
            "Tref": 1, "rhoref": 1.*M,

            "Tref_res": 132.6312, "rhoref_res": 10.4477*M,
            "n_poly": [10.72, 1.122, 0.002019, -8.876, -0.02916],
            "t_poly": [.2, .05, 2.4, .6, 3.6],
            "d_poly": [1, 4, 9, 1, 8],
            "g_poly": [0, 0, 0, 1, 1],
            "c_poly": [0, 0, 0, 1, 1]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Lemmon (2004)",
               "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                            "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air",
                            "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                            "doi": "10.1023/B:IJOT.0000022327.04529.f3"},
               "__test__": """
                    >>> st=Air(T=100, rhom=0)
                    >>> print "%0.5f" % st.k.mWmK
                    9.35902
                    >>> st=Air(T=300, rhom=0)
                    >>> print "%0.4f" % st.k.mWmK
                    26.3529
                    >>> st=Air(T=100, rhom=28)
                    >>> print "%0.3f" % st.k.mWmK
                    119.221
                    >>> st=Air(T=200, rhom=10)
                    >>> print "%0.4f" % st.k.mWmK
                    35.3185
                    >>> st=Air(T=300, rhom=5)
                    >>> print "%0.4f" % st.k.mWmK
                    32.6062
                    >>> st=Air(T=132.64, rhom=10.4)
                    >>> print "%0.4f" % st.k.mWmK
                    75.6231
                    """, # Table V, Pag 28

               "Tref": 132.6312, "kref": 1e-3,
               "no": [1.308, 1.405, -1.036],
               "co": [-97, -1.1, -0.3],

               "Trefb": 132.6312, "rhorefb": 10.4477, "krefb": 1e-3,
               "nb": [8.743, 14.76, -16.62, 3.793, -6.142, -0.3778],
               "tb": [0.1, 0, 0.5, 2.7, 0.3, 1.3],
               "db": [1, 2, 3, 7, 7, 11],
               "cb": [0, 0, 2, 2, 2, 2],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.11e-9, "gam0": 0.55e-1, "qd": 0.31e-9, "Tcref": 265.262}

    _thermal = thermo0,
