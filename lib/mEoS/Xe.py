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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from math import log, exp
from unittest import TestCase

from lib import unidades
from lib.meos import MEoS


class Xe(MEoS):
    """Multiparameter equation of state for xenon"""
    name = "xenon"
    CASNumber = "7440-63-3"
    formula = "Xe"
    synonym = ""
    _refPropName = "XENON"
    _coolPropName = "Xenon"
    rhoc = unidades.Density(1102.8612)
    Tc = unidades.Temperature(289.733)
    Pc = unidades.Pressure(5842.0, "kPa")
    M = 131.293  # g/mol
    Tt = unidades.Temperature(161.405)
    Tb = unidades.Temperature(165.05)
    f_acent = 0.00363
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 994

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-3.8227178129, 3.8416395351]}

    CP1 = {"ao": 2.5}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for xenon of Lemmon "
                    "and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi": "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 700000.0, "rhomax": 28.78,

        "nr1": [0.83115, -2.3553, 0.53904, 0.014382, 0.066309, 0.00019649],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.14996, -0.035319, -0.15929, -0.027521, -0.023305, 0.0086941],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = (lemmon, )
    _PR = [-0.5639, -12.9972]

    # Other equation with unorthodox formulations not implemented
    # S̆ifner, O., Klomfar, J.
    # Thermodynamic Properties of Xenon from the Triple Point to 800 K with
    # Pressures up to 350 MPa.
    # J. Phys. Chem. Ref. Data 23(1) (1994) 63-152
    # http:/dx.doi.org/10.1063/1.555956

    _surface = {"sigma": [-0.11538, 0.16598], "exp": [1.0512, 1.098]}
    _dielectric = {
        "eq": 1,
        "a": [10.122, 0], "b": [31.97, 46.97], "c": [-948.4, 0],
        "Au": 0, "D": 1.7}

    _melting = {
        "eq": 2,
        "__doi__": {"autor": "Michels, A., Prins, C.",
                    "title": "The Melting Lines of Argon, Krypton and Xenon "
                             "up to 1500 Atm; Representation of the Results "
                             "by a Law of Corresponding States",
                    "ref": "Physica 28 (1962) 101-116",
                    "doi": "10.1016/0031-8914(62)90096-4"},

        "Tmin": Tt, "Tmax": 1300.0,
        "Tref": 1, "Pref": -2576*101325,
        "a1": [0.7983277027965369*101325], "exp1": [1.589165]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.60231e1, 0.14989e1, -0.74906, -0.12194e1, -0.44905],
        "t": [1., 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.13570e2, -0.47545e2, 0.63876e2, -0.39983e2, 0.12701e2],
        "t": [0.56, 0.8, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-3.0026, -6.056, -0.60339e2, 0.48838e3, -0.81974e3, 0.47287e3],
        "t": [0.435, 1.4, 4.4, 6.2, 7.0, 8.6]}

    visco0 = {"__name__": "Velliadou (2021)",
              "__doi__": {
                  "autor": "Velliadou, D., Tasidou, K.A., Antoniadis, K.D., "
                           "Assael, M.J., Perkins, R.A., Huber, M.L.,",
                  "title": "Reference Correlation for the Viscosity of Xenon "
                           "from the Triple Point to 750 K and up to 86 MPa",
                  "ref": "Int. J. Thermophysics 42 (2021) 74",
                  "doi": "10.1007/s10765-021-02818-9"},

              "eq": 1, "omega": 0,

              "special0": "_mu0",

              "ek": 250, "sigma": 0.396,

              "Tref_virial": 250,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "special": "_mur",

              "critical": 0, "Xic": 0.06e-9,
              "xu": 0.068, "gnu": 0.63, "gamma": 1.239, "Xio": 0.184e-9,
              "gam0": 0.058, "qc": 3.6e-9, "qd": 1.15e-9, "Tcref": 1.5*Tc}

    def _mu0(self, T):
        """Special term for zero-density viscosity for Velliadou correlation"""
        ai = [9.652514e-1, -5.237199e-2, -6.758414e-2, 2.855787e-2,
              1.002789e-2, -9.639621e-3, 1.329770e-3, 1.114305e-3,
              -5.992234e-4, 1.224218e-4, -9.584978e-6]

        # Eq 2
        suma = 0
        for i, a in enumerate(ai, 1):
            suma += a*log(T/298.15)**i
        return 23.0183*exp(suma)

    def _mur(self, rho, T, fase):
        """Special term of residual viscosity for Velliadou correlation"""
        Tr = T/self.Tc
        rhor = rho/self.rhoc

        # Eq 6
        mur = rhor**(2/3)*Tr**0.5 * (
            Tr + 1.396328251*Tr*rhor**4 + 5.418871011e-4*rhor**12/Tr
            + (4.478809952+24.91698858*rhor)/Tr**2)
        return mur

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Velliadou (2021)",
               "__doi__": {
                   "autor": "Velliadou, D., Assael, M.J., Antoniadis, K.D., "
                            "Huber, M.L.",
                   "title": "Reference Correlations for the Thermal "
                            "Conductivity of Xenon from the Triple Point to "
                            "606 K and Pressures up to 400 MPa",
                   "ref": "Int. J. Thermophysics 42 (2021) 51",
                   "doi": "10.1007/s10765-021-02803-2"},

               "eq": 1,

               "special": "_thermo0",

               "Tref_res": Tc, "rhoref_res": rhoc, "kref_res": 1,
               "nr": [0.694552e-2, -0.732747e-4, 0.876111e-2, -0.268366e-2,
                      -0.119900e-1, 0.563598e-2, 0.684476e-2, -0.314076e-2,
                      -0.102229e-2, 0.605394e-3],
               "tr": [0, -1, 0, -1, 0, -1, 0, -1, 0, -1],
               "dr": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.182e-9, "gam0": 0.058, "qd": 0.479e-9, "Tcref": 434.6}

    _thermal = (thermo0, )

    def _thermo0(self, rho, T, fase):
        """Custom dilute-gas limit thermal conductivity for Velliadou method"""
        Tr = T/298.15
        bi = [9.65520e-1, -5.12353e-2, -6.70913e-2, 2.88938e-2, 9.25546e-3,
              -9.72175e-3, 1.69364e-3, 9.96803e-4, -6.10466e-4, 1.33327e-4,
              -1.09858e-5]

        # Eq 2
        ko = 5.4666*exp(sum(b*log(Tr)**(i+1) for i, b in enumerate(bi)))

        return ko*1e-3


class Test(TestCase):
    """Testing"""
    def test_shortLemmon(self):
        """Table 10, Pag 842"""
        st = Xe(T=291, rhom=8)
        self.assertEqual(round(st.P.kPa, 3), 5986.014)
        self.assertEqual(round(st.hM.kJkmol, 3), 9193.668)
        self.assertEqual(round(st.sM.kJkmolK, 3), 36.895)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 28.692)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3063.309)
        self.assertEqual(round(st.w, 3), 125.648)

    def test_Michels(self):
        """Table I, pag 105"""
        self.assertEqual(round(Xe._Melting_Pressure(161.554).atm, 2), 3.66)
        self.assertEqual(round(Xe._Melting_Pressure(167.154).atm, 2), 147.21)
        self.assertEqual(round(Xe._Melting_Pressure(191.144).atm, 2), 794.08)
        self.assertEqual(round(Xe._Melting_Pressure(215.264).atm, 2), 1494.59)

    def test_VelliadouVisco(self):
        """Point data given in Section 4"""
        self.assertEqual(round(Xe(T=300, rho=0).mu.muPas, 4), 23.1561)
        self.assertEqual(round(Xe(T=300, rho=6).mu.muPas, 4), 23.3186)
        self.assertEqual(round(Xe(T=300, rho=2500).mu.muPas, 3), 206.449)
        self.assertEqual(round(Xe(T=292.711322, rho=0).mu.muPas, 4), 22.6125)

        # The critical enhancement fail
        # self.assertEqual(round(
        #     Xe(T=292.711322, rho=1102.9).mu.muPas, 5), 52.82074)

    def test_VelliadouThermo(self):
        """Point data given in Section 4"""
        self.assertEqual(round(Xe(T=300, rho=0).k.mWmK, 4), 5.4993)
        self.assertEqual(round(Xe(T=300, rho=1200).k.mWmK, 4), 22.7675)


if __name__ == "__main__":
    # Generate Fig. 6 of viscosity paper
    # To get work, you must enable viscosity critical enhancement parameter for
    # Xenon, change 0 by 1 in critical viscosity parameters
    from numpy import logspace, r_
    import matplotlib.pyplot as plt

    t = logspace(-4.8, -1.2, 50)
    Ti = t*Xe.Tc+Xe.Tc

    mu = []
    mu1 = []
    mu2 = []
    for ti in Ti:
        state = Xe(T=ti, rho=Xe.rhoc)
        mu.append(state.mu.muPas)
        state = Xe(T=ti, rho=Xe.rhoc, viscocriticallineal=True)
        mu1.append(state.mu.muPas)
        state = Xe(T=ti, rho=Xe.rhoc, viscocritical=False)
        mu2.append(state.mu.muPas)

    plt.plot(t, mu, color="k", ls="-")
    plt.plot(t, mu1, color="k", ls="--")
    plt.plot(t, mu2, color="k", ls=":")

    # Experimental data from
    # Berg, R.F., Moldover, M.R., Zimmerli, G.A.
    # Frequency-dependent viscosity of xenon near the critical point
    # Physycal Review E. 60(4) (1999) 4079-4098
    # doi: 10.1103/physreve.60.4079
    t3 = r_[1029.0, 344.3, 104.1, 100.0, 96.35, 65.51, 59.60, 43.36, 35.89,
            35.48, 34.76, 29.69, 21.18, 20.76, 17.07, 14.93, 11.73, 11.03,
            8.645, 8.620, 6.808, 5.886, 4.978, 4.381, 3.142, 2.960, 2.211]*1e-5
    mu3 = [52.241, 53.009, 54.509, 54.585, 54.576, 55.213, 55.275, 55.923,
           56.305, 56.339, 56.349, 56.709, 57.359, 57.389, 57.824, 58.080,
           58.658, 58.784, 59.291, 59.351, 59.876, 60.246, 60.610, 60.937,
           61.765, 61.930, 62.713]
    plt.plot(t3, mu3, ls="None", marker="^", mec="k", mfc="k")

    t4 = r_[5527, 3976, 2420, 1730, 1729, 1729, 1729, 1729, 1039, 1034, 354.2,
            338.0, 220.9, 140.9, 105.6, 105.4]*1e-5
    mu4 = [52.623, 52.316, 52.073, 52.047, 52.048, 52.048, 52.047, 52.053,
           52.178, 52.162, 52.928, 52.974, 53.446, 54.011, 54.375, 54.360]
    plt.plot(t4, mu4, ls="None", marker="o", mec="k", mfc="None")

    plt.xscale("log")
    plt.show()
