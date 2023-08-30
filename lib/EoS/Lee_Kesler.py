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


###############################################################################
# Lee-Kesler equation of state implementation
###############################################################################


from numpy import exp, r_
from numpy.lib.scimath import log
from scipy.constants import R
from scipy.optimize import fsolve

from lib.compuestos import RhoL_Costald

from lib.bip import Kij
from lib.eos import EoS
from lib.physics import R_atml, factor_acentrico_octano


# Table 1
coef = {
    "b1": (0.1181193, 0.2026579),
    "b2": (0.265728, 0.331511),
    "b3": (0.154790, 0.027655),
    "b4": (0.030323, 0.203488),
    "c1": (0.0236744, 0.0313385),
    "c2": (0.0186984, 0.0503618),
    "c3": (0.0, 0.016901),
    "c4": (0.042724, 0.041577),
    "d1": (0.155488e-4, 0.48736e-4),
    "d2": (0.623689e-4, 0.0740336e-4),
    "beta": (0.65392, 1.226),
    "gamma": (0.060167, 0.03754)}


# Table I from Plöcker
bip = {
    "2-3": 1.052,
    "2-22": 1.014,
    "2-4": 1.113,
    "2-23": 1.089,
    "2-6": 1.171,
    "2-5": 1.155,
    "2-8": 1.240,
    "2-7": 1.228,
    "2-10": 1.304,
    "2-38": 1.269,
    "2-40": 1.234,
    "2-11": 1.367,
    "2-12": 1.423,
    "2-13": 1.484,
    "2-14": 1.533,
    "3-22": 0.991,
    "3-4": 1.010,
    "3-23": 1.002,
    "3-6": 1.029,
    "3-5": 1.036,
    "3-8": 1.064,
    "3-7": 1.070,
    "3-10": 1.106,
    "3-38": 1.081,
    "3-40": 1.066,
    "3-11": 1.143,
    "3-12": 1.165,
    "3-13": 1.214,
    "3-14": 1.237,
    "22-6": 0.998,
    "22-40": 1.094,
    "22-11": 1.163,
    "65-22": 0.948,
    "4-23": 0.992,
    "4-6": 1.003,
    "4-5": 1.003,
    "4-8": 1.006,
    "4-7": 1.009,
    "4-10": 1.047,
    "4-38": 1.037,
    "4-40": 1.011,
    "4-11": 1.067,
    "4-12": 1.090,
    "4-13": 1.115,
    "4-14": 1.139,
    "23-6": 1.010,
    "23-5": 1.009,
    "23-27": 1.006,
    "6-5": 1.001,
    "6-8": 0.994,
    "6-7": 0.998,
    "6-10": 1.018,
    "6-38": 1.008,
    "6-40": 0.999,
    "6-11": 1.027,
    "6-12": 1.046,
    "6-13": 1.064,
    "6-14": 1.078,
    "8-7": 0.987,
    "8-10": 0.996,
    "8-38": 0.996,
    "8-40": 0.977,
    "8-11": 1.004,
    "8-12": 1.020,
    "8-13": 1.033,
    "8-14": 1.045,
    "10-38": 0.998,
    "10-40": 0.978,
    "10-11": 1.008,
    "10-12": 1.005,
    "10-13": 1.015,
    "10-14": 1.025,
    "40-38": 0.979,
    "40-11": 0.985,
    "40-12": 0.987,
    "40-82": 0.982,
    "40-13": 1.034,
    "40-14": 1.047,
    "38-11": 0.999,
    "38-12": 1.010,
    "38-13": 1.021,
    "38-14": 1.032,
    "11-12": 0.993,
    "11-82": 1.002,
    "11-13": 1.002,
    "11-14": 1.010,
    "12-13": 0.993,
    "12-14": 0.999,
    "13-14": 0.991,
    "46-2": 0.977,
    "46-22": 1.032,
    "46-3": 1.082,
    "46-4": 1.177,
    "46-23": 1.151,
    "46-6": 1.276,
    "46-8": 1.372,
    "46-10": 1.442,
    "46-47": 0.997,
    "46-48": 0.987,
    "46-98": 0.988,
    "46-50": 0.983,
    "46-49": 1.110,
    "46-110": 1.073,
    "46-63": 1.033,
    "49-2": 0.975,
    "49-3": 0.938,
    "49-4": 0.925,
    "49-6": 0.955,
    "49-5": 0.946,
    "49-8": 1.002,
    "49-10": 1.018,
    "49-38": 1.054,
    "49-40": 1.018,
    "49-11": 1.058,
    "49-12": 1.090,
    "49-13": 1.126,
    "49-14": 1.160,
    "49-50": 0.922,
    "49-216": 0.969,
    "49-117": 1.069,
    "1-2": 1.216,
    "1-3": 1.604,
    "1-22": 1.498,
    "1-4": 1.826,
    "1-6": 2.093,
    "1-8": 2.335,
    "1-10": 2.456,
    "1-11": 2.634,
    "1-46": 1.080,
    "1-48": 1.085,
    "1-49": 1.624,
    "98-47": 0.985,
    "98-2": 1.010,
    "98-63": 0.984,
    "47-110": 1.057,
    "48-2": 0.974,
    "971-47": 0.989,
    "50-5": 0.947,
    "110-2": 1.017,
    "62-49": 0.920,
    "62-63": 1.152,
    "62-117": 0.979}


class Lee_Kesler(EoS):
    r"""
    Corresponding state equation of state of Lee-Kesler

    .. math::
        \begin{array}[t]{l}
        Z = Z^{(0)} + \omega Z^{(1)}\\
        Z = \frac{P_rV_r}{T_r} = 1 + \frac{B}{V_r} + \frac{C}{V_r^2} +
        \frac{D}{V_r^5} + \frac{c_4}{T_r^3V_r^2}\left(\beta+\frac{\gamma}
        {V_r^2}\right)\exp{\left(-\frac{\gamma}{V_r^2}\right)}\\
        B = b_1 - \frac{b_2}{T_r} - \frac{b_3}{T_r^2}-\frac{b_4}{T_r^3}\\
        C = c_1 - \frac{c_2}{T_r} - \frac{c_3}{T_r^3}\\
        D = d_1 - \frac{d_2}{T_r}\\
        \end{array}

    Using the mixing rules defined by Plöcker [2]_

    .. math::
        \begin{array}[t]{l}
        T_{CM} = \frac{1}{v_{CM}^\nu \sum_j \sum_k Z_j Z_k v_{Cjk}^\nu T_{Cjk}\\
        v_{CM} = \sum_j \sum_k Z_j Z_k v_{Cjk}\\
        \omega_M = \sum_j Z_j \omega_j\\
        P_{CM} = \left(0.2905-0.085\omega_m\right) R \frac{T_{CM}}{v_{CM}\\
        \end{array}

    with the critical cross parameters

    .. math::
        \begin{array}[t]{l}
        T_{Cjk} = \left(T_{Cj} T_{Ck}\right)^{1/2) k_{jk}\\
        v_{Cjk} = \frac{1}{8} \left(v_{Cj}^{1/3} + v_{Ck}^{1/3}\right)^3\\
        \end{array}

    Examples
    --------
    Example 1.17 from [3]_, Propane compressibility

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(1, ids=[8], caudalUnitarioMasico=[1.])
    >>> eq = Lee_Kesler(0.8*mix.Tc, 0.4*mix.Pc, mix)
    >>> '%0.4f' % (eq.Z[0])
    '0.0592'
    >>> eq = Lee_Kesler(0.9*mix.Tc, 0.4*mix.Pc, mix)
    >>> '%0.4f' % (eq.Z[1])
    '0.7520'
    >>> eq = Lee_Kesler(1*mix.Tc, 0.4*mix.Pc, mix)
    >>> '%0.4f' % (eq.Z[1])
    '0.8437'
    >>> eq = Lee_Kesler(1.1*mix.Tc, 0.4*mix.Pc, mix)
    >>> '%0.4f' % (eq.Z[1])
    '0.8939'
    >>> eq = Lee_Kesler(1.2*mix.Tc, 0.4*mix.Pc, mix)
    >>> '%0.4f' % (eq.Z[1])
    '0.9253'
    """

    __title__ = "Lee Kesler"
    __status__ = "LK"

    __doi__ = (
        {"autor": "Lee, B.I., Kesler, M.G.",
         "title": "A Generalized Thermodynamic Correlation Based on "
                  "Three-Parameter Corresponding States",
         "ref": "AIChE Journal 21(3) (1975) 510-527",
         "doi": "10.1002/aic.690210313"},
        {"autor": "Plöcker, U., Knapp, H., Prausnitz, J.",
         "title": "Calculation of High-Pressure Vapor-Liquid Equilibria from "
                  "a Corresponding-States Correlation with Emphasis on "
                  "Asymmetric Mixtures",
         "ref": "Ind. Eng. Chem. Process Des. Dev. 17(3) (1978) 324-332",
         "doi": "10.1021/i260067a020"},
        {"autor": "Walas, S.M.",
         "title": "Phase Equilibria in Chemical Engineering",
         "ref": "Butterworth, 1985",
         "doi": ""},
        {"autor": "Joffe, J.",
         "title": "Vapor-Liquid Equilibria by the Pseudocritical Method",
         "ref": "Ind. Eng. Chem. Fundam. 15(4) (1976) 298-303",
         "doi": "10.1021/i160060a013"},

        )

    def __init__(self, T, P, mezcla):
        EoS.__init__(self, T, P, mezcla)
        self.kij = Kij(mezcla.ids, "LK", bip)

        self.Zci = [0.2905-0.085*cmp.f_acent for cmp in self.componente]
        self.Vci = [zc*R*cmp.Tc/cmp.Pc.kPa
                    for zc, cmp in zip(self.Zci, self.componente)]
        self.Tci = [cmp.Tc for cmp in self.componente]
        self.wi = mezcla._arraylize("f_acent")

        # Cross critical coefficient
        # Eq 12 from Plöcker
        self.Tcij = []
        for tci, kiji in zip(self.Tci, self.kij):
            tcij = []
            for tcj, kij in zip(self.Tci, kiji):
                tcij.append((tci*tcj)**0.5*kij)
            self.Tcij.append(tcij)

        # Eq 13 from Plöcker
        self.Vcij = []
        for vci in self.Vci:
            vcij = []
            for vcj in self.Vci:
                vcij.append((vci**(1/3)+vcj**(1/3))**3/8)
            self.Vcij.append(vcij)

        self.Z = self._Z(self.zi, T, P)
        self.V = self.Z*R_atml*self.T/self.P.atm  # mol/l
        self._fug(self.zi, self.zi, T, P)
        # self.H_exc = r_[Hv, Hl]
        self.x, self.xi, self.yi, self.Ki = self._Flash()

    def _mix(self, zi):
        """Mixing rules, Eq 20-25"""
        # Using the binary interaction parameters defined in [2]_
        sumV = 0
        sumT = 0
        for xj, Vcj, Tcj, kiji in zip(zi, self.Vci, self.Tci, self.kij):
            for xk, Vck, Tck, kij in zip(zi, self.Vci, self.Tci, kiji):
                sumV += xj*xk*(Vcj**(1/3)+Vck**(1/3))**3
                sumT += xj*xk*(Vcj**(1/3)+Vck**(1/3))**3*(Tcj*Tck)**0.5*kij
        Vc = sumV/8
        Tc = sumT/8/Vc
        Pc = (0.2905-0.085*self.mezcla.f_acent)*R_atml*Tc/Vc*101325

        return Tc, Pc, Vc

    def _lib(self, zi=None, T=None, P=None, rho0=None):

        Tc, Pc, Vc = self._mix(zi)

        Tr = T/Tc
        Pr = P/Pc

        # Table 1
        b1 = 0.1181193, 0.2026579
        b2 = 0.265728, 0.331511
        b3 = 0.154790, 0.027655
        b4 = 0.030323, 0.203488
        c1 = 0.0236744, 0.0313385
        c2 = 0.0186984, 0.0503618
        c3 = 0.0, 0.016901
        c4 = 0.042724, 0.041577
        d1 = 0.155488e-4, 0.48736e-4
        d2 = 0.623689e-4, 0.0740336e-4
        beta = 0.65392, 1.226
        gamma = 0.060167, 0.03754


        Bo = b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
        Co = c1[0]-c2[0]/Tr+c3[0]/Tr**3
        Do = d1[0]+d2[0]/Tr

        def Vr(V):
            Vr = 1 + Bo/V + Co/V**2 + Do/V**5 + c4[0]/Tr**3/V**2 * \
                (beta[0]+gamma[0]/V**2) * exp(-gamma[0]/V**2)-Pr*V/Tr
            return Vr

        Bh = b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
        Ch = c1[1]-c2[1]/Tr+c3[1]/Tr**3
        Dh = d1[1]+d2[1]/Tr

        def Vrh(V):
            Vrh = 1 + Bh/V + Ch/V**2 + Dh/V**5 + c4[1]/Tr**3/V**2 * \
                (beta[1]+gamma[1]/V**2) * exp(-gamma[1]/V**2)-Pr*V/Tr
            return Vrh

        # Used initial values for iteration
        Vlo = RhoL_Costald(T, Tc, self.mezcla.f_acent, Vc)
        Vgo = R_atml*T/P

        vr0v = fsolve(Vr, Vgo)[0]
        vrhv = fsolve(Vrh, Vgo)[0]
        vr0l = fsolve(Vr, Vlo)[0]
        vrhl = fsolve(Vrh, Vlo)[0]

        return Tr, Pr, vr0v, vrhv, vr0l, vrhl

    def _Z(self, zi=None, T=None, P=None, rho0=None):
        Tr, Pr, vr0v, vrhv, vr0l, vrhl = self._lib(zi, T, P, rho0)

        z0l = Pr*vr0l/Tr
        zhl = Pr*vrhl/Tr
        z0v = Pr*vr0v/Tr
        zhv = Pr*vrhv/Tr
        return r_[z0v+self.mezcla.f_acent/factor_acentrico_octano*(zhv-z0v),
                  z0l+self.mezcla.f_acent/factor_acentrico_octano*(zhl-z0l)]

    def _fug(self, xi, yi, T, P):
        """Fugacities of component in mixture calculation

        Parameters
        ----------
        xi : list
            Molar fraction of component in liquid phase, [-]
        yi : list
            Molar fraction of component in vapor phase, [-]

        Returns
        -------
        tital : list
            List with liquid phase component fugacities
        titav : list
            List with vapour phase component fugacities
        """
        Trl, Prl, vr0v_, vrhv_, vr0l, vrhl = self._lib(xi, T, P)
        Trv, Prv, vr0v, vrhv, vr0l_, vrhl_ = self._lib(yi, T, P)
        Tcl, Pcl, Vcl = self._mix(xi)
        Tcv, Pcv, Vcv = self._mix(yi)

        z0l = Prl*vr0l/Trl
        zhl = Prl*vrhl/Trl
        z0v = Prv*vr0v/Trv
        zhv = Prv*vrhv/Trv
        zl = z0l+self.mezcla.f_acent/factor_acentrico_octano*(zhl-z0l)
        zv = z0v+self.mezcla.f_acent/factor_acentrico_octano*(zhv-z0v)

        f0v = self._fugacity(0, Trv, Prv, vr0v)
        fhv = self._fugacity(1, Trv, Prv, vrhv)
        f0l = self._fugacity(0, Trl, Prl, vr0l)
        fhl = self._fugacity(1, Trl, Prl, vrhl)
        fv = f0v+self.mezcla.f_acent/factor_acentrico_octano*(fhv-f0v)
        fl = f0l+self.mezcla.f_acent/factor_acentrico_octano*(fhl-f0l)

        # Eq A8
        dfvdw = 1/factor_acentrico_octano*(fhv-f0v)
        dfldw = 1/factor_acentrico_octano*(fhl-f0l)

        H0v = self._Hexc_lib(0, Trv, Prv, vr0v)
        Hhv = self._Hexc_lib(1, Trv, Prv, vrhv)
        H0l = self._Hexc_lib(0, Trl, Prl, vr0l)
        Hhl = self._Hexc_lib(1, Trl, Prl, vrhl)
        hv = H0v+self.mezcla.f_acent/factor_acentrico_octano*(Hhv-H0v)
        hl = H0l+self.mezcla.f_acent/factor_acentrico_octano*(Hhl-H0l)

        # Eq A13
        dwZ = []
        dzZ = []
        for wi in self.wi:
            dwZi = []
            dzZi = []
            for wj in self.wi:
                dwZi.append(wj-wi)
                dzZi.append(-0.085*(wj-wi))
            dwZ.append(dwZi)
            dzZ.append(dzZi)

        # Eq A12 for both phases
        dVcx = []
        dVcy = []
        for vci in self.Vci:
            dVcxi = 0
            dVcyi = 0
            for xj, yj, vcj in zip(xi, yi, self.Vci):
                dVcxi += xj*(vcj-vci)
                dVcyi += yj*(vcj-vci)
            dVcx.append(2*dVcxi)
            dVcy.append(2*dVcyi)

        # Eq A9 for both phases
        dTcZl = []
        dTcZv = []
        for x, y, vci, tci, vcij, tcij, dvcx, dvcy in \
                zip(xi, yi, self.Vci, self.Tci, self.Vcij, self.Tcij, dVcx, dVcy):
            suml = 0
            sumv = 0
            for vcj, tcj in zip(vcij, tcij):
                suml += x*(vcj**0.25*tcj-vci**0.25*tci)
                sumv += y*(vcj**0.25*tcj-vci**0.25*tci)
            dTcZl.append((2*suml-0.25*Vcl**-0.75*dvcx*Tcl)/Vcl**0.25)
            dTcZv.append((2*sumv-0.25*Vcv**-0.75*dvcy*Tcv)/Vcv**0.25)

        # Eq A9 for both phases
        dPcZl = []
        dPcZv = []
        for x, y, dzi, dTli, dTvi, dvcx, dvcy in \
                zip(xi, yi, dzZ, dTcZl, dTcZv, dVcx, dVcy):
            dpli = []
            dpvi = []
            for dz, dTl, dTv in zip(dzi, dTli, dTvi):
                dPli.append(Pcl*(dz/zl+dTl/Tcl-dvcx/Vcl))
                dPvi.append(Pcv*(dz/zv+dTv/Tcv-dvcy/Vcv))
            dPcZl.append(dPli)
            dPcZv.append(dPvi)
        print("dPc", dPcZl)

        sumTcx = []
        sumTcy = []
        for i in range(len(xi)):
            sumx = 0
            sumy = 0
            for j, dTcl in enumerate(dTcZl):
                if i != j:
                    sumx += xi[j]*dTcl
                    sumy += yi[j]*dTcZv[j]
            sumTcx.append(sumx)
            sumTcy.append(sumy)

        sumPcx = []
        sumPcy = []
        for i in range(len(xi)):
            sumx = 0
            sumy = 0
            for j, dPcl in enumerate(dPcZl):
                if i != j:
                    sumx += xi[j]*dPcl
                    sumy += yi[j]*dPcZv[j]
            sumPcx.append(sumx)
            sumPcy.append(sumy)

        # print("sumPcx", sumPcx)

        sumw = []
        for i in range(len(xi)):
            sum = 0
            for j, dw in enumerate(dwZ[i]):
                if i != j:
                    sum += xi[j]*dw
            sumw.append(sum)

        # Eq A7
        fil = []
        for sTc, sPc, sw in zip(sumTcx, sumPcx, sumw):
            fil.append(fl-hl/R/T/Tcl*sTc + (zl-1)/Pcl*sPc - dfldw*sw)
        fiv = []
        for sTc, sPc, sw in zip(sumTcy, sumPcy, sumw):
            fiv.append(fv-hv/R/T/Tcv*sTc + (zv-1)/Pcv*sPc - dfvdw*sw)

        print("fil", fil)
        print("fiv", fiv)
        return fil, fiv


    def _fugacity(self, ref, Tr, Pr, vr):
        """Calculate fugacity coefficient of mixture

        Parameters
        ----------
        ref : integer
            Index of parameter, 0 for simple fluid and 1 for reference fluid
        Tr : float
            Reduced temperature
        Pr : float
            Reduced pressure
        vr : float
            Reduced volume
        """
        b1 = coef["b1"][ref]
        b2 = coef["b2"][ref]
        b3 = coef["b3"][ref]
        b4 = coef["b4"][ref]
        c1 = coef["c1"][ref]
        c2 = coef["c2"][ref]
        c3 = coef["c3"][ref]
        c4 = coef["c4"][ref]
        d1 = coef["d1"][ref]
        d2 = coef["d2"][ref]
        gamma = coef["gamma"][ref]
        beta = coef["beta"][ref]

        z = Pr*vr/Tr

        E = c4/(2*Tr**3*gamma)*(
            beta + 1 - (beta+1+gamma/vr**2)*exp(-gamma/vr**2))
        f = z - 1 - log(z) + (b1 - b2/Tr + b3/Tr**2 + b4/Tr**3)/vr \
            + (c1 - c2/Tr + c3/Tr**3)/2/vr**2 + (d1 + d2)/5/vr**5 + E
        return f

    def _Hexc(self, zi, T, P):

        Tr, Pr, vr0v, vrhv, vr0l, vrhl = self._lib(zi, T, P)
        z0l = Pr*vr0l/Tr
        zhl = Pr*vrhl/Tr
        z0v = Pr*vr0v/Tr
        zhv = Pr*vrhv/Tr

        H0v = self._Hexc_lib(0, Tr, Pr, vr0v)
        Hhv = self._Hexc_lib(1, Tr, Pr, vrhv)
        H0l = self._Hexc_lib(0, Tr, Pr, vr0l)
        Hhl = self._Hexc_lib(1, Tr, Pr, vrhl)

        Hv = H0v+mezcla.f_acent/factor_acentrico_octano*(Hhv-H0v)
        Hl = H0l+mezcla.f_acent/factor_acentrico_octano*(Hhl-H0l)
        return Hv, Hl

    def _Hexc_lib(self, ref, Tr, Pr, vr):
        """Calculate enthalpy excess of mixture

        Parameters
        ----------
        ref : integer
            Index of parameter, 0 for simple fluid and 1 for reference fluid
        Tr : float
            Reduced temperature
        Pr : float
            Reduced pressure
        vr : float
            Reduced volume
        """
        b1 = coef["b1"][ref]
        b2 = coef["b2"][ref]
        b3 = coef["b3"][ref]
        b4 = coef["b4"][ref]
        c1 = coef["c1"][ref]
        c2 = coef["c2"][ref]
        c3 = coef["c3"][ref]
        c4 = coef["c4"][ref]
        d1 = coef["d1"][ref]
        d2 = coef["d2"][ref]
        gamma = coef["gamma"][ref]
        beta = coef["beta"][ref]

        z = Pr*vr/Tr

        E = c4/(2*Tr**3*gamma)*(
            beta + 1 - (beta+1+gamma/vr**2)*exp(-gamma/vr**2))
        H0 = -Tr*(z - 1 - (b2 + 2*b3/Tr + 3*b4/Tr**2)/Tr/vr
                  - c2/Tr/2/vr**2 + d2/5/Tr/vr**5 + 3*E)
        return H0


    # def Cp_Lee_Kesler(self, T, P, fase=None):
        # """Método alternativo para el cálculo de la capacidad calorífica
        # Procedure API 7D3.6 Pag.711"""
        # Tr=self.tr(T)
        # if fase==None:
            # fase=self.Fase(T, P)
        # Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(T, P, fase)

        # B=0.1181193-0.265728/Tr-0.154790/Tr**2-0.030323/Tr**3
        # C=0.0236744-0.0186984/Tr
        # D=0.155488e-4+0.623689e-4/Tr
        # dpdt_0=1/vr0*(1+(0.1181193+0.154790/Tr**2+2*0.030323/Tr**3)/vr0+0.0236744/vr0**2+0.155488e-4/vr0**5-2*0.042724/Tr**3/vr0**2*((0.65392+0.060167/vr0**2)*exp(-0.060167/vr0**2)))
        # dpdv_0=-Tr/vr0**2*(1+2*B/vr0+3*C/vr0**2+6*D/vr0**5+0.042724/Tr**3/vr0**2*(3*0.65392+(5-2*(0.65392+0.060167/vr0**2))*0.060167/vr0**2)*exp(-0.060167/vr0**2))
        # Cp0=1+Tr*dpdt_0**2/dpdv_0+Cv0

        # B=0.2026579-0.331511/Tr-0.027655/Tr**2-0.203488/Tr**3
        # C=0.0313385-0.0503618/Tr+0.016901/Tr**3
        # D=0.48736e-4+0.0740336e-4/Tr
        # dpdt_h=1/vrh*(1+(0.2026579+0.027655/Tr**2+2*0.203488/Tr**3)/vrh+(0.0313385-2*0.016901/Tr**3)/vrh**2+0.48736e-4/vrh**5-2*0.041577/Tr**3/vrh**2*((1.226+0.03754/vrh**2)*exp(-0.03754/vrh**2)))
        # dpdv_h=-Tr/vrh**2*(1+2*B/vrh+3*C/vrh**2+6*D/vrh**5+0.041577/Tr**3/vrh**2*(3*1.226+(5-2*(1.226+0.03754/vrh**2))*0.03754/vrh**2)*exp(-0.03754/vrh**2))
        # Cph=1+Tr*dpdt_h**2/dpdv_h+Cvh

        # Cp_adimensional=Cp0+self.f_acent/factor_acentrico_octano*(Cph-Cp0)
        # return unidades.SpecificHeat(self._Cpo(T).JgK-R/self.M*Cp_adimensional, "JgK")

    # def Cv_Lee_Kesler(self, T, P, fase=None):
        # """Método de cálculo de la capacidad calorífica a volumen constante
        # Procedure API 7E1.6 Pag.726"""
        # #FIXME: No sale, un factor de 100 tengo que añadir no sé de donde
        # Pr=P/self.Pc
        # Tr=T/self.Tc
        # if fase==None:
            # fase=self.Fase(T, P)
        # Cpo=self._Cpo(T)
        # Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr, fase)
        # Cv_adimensional=Cv0+self.f_acent/factor_acentrico_octano*(Cvh-Cv0)
        # return unidades.SpecificHeat(100*(Cpo.JgK-R/self.M*(1+Cv_adimensional)), "JgK")


    # def Cp_Cv_Lee_Kesler(self, T, P):
        # """Método de cálculo de la capacidad calorífica a volumen constante
        # Procedure API 7E1.6 Pag.726"""
        # Cv=self.Cv_Lee_Kesler(T, P.atm)
        # Cp=self.Cp_Lee_Kesler(T, P.atm)
# #        print Cp.BtulbF, Cv
        # return Cp/Cv


    # def Fugacidad_Lee_Kesler(self, T, P):
        # """Método de cálculo de la fugacidad
        # Procedure API 7G1.8 Pag.752"""
        # Tr=T/self.Tc
        # Pr=P/self.Pc
        # f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P))
        # return unidades.Pressure(P*exp(f), "atm")

    # def Entropia_Lee_Kesler(self, T, P):
        # """Método de cálculo de la entropia
        # Procedure API 7F1.7 Pag.739"""
        # Tr=T/self.Tc
        # Pr=P/self.Pc
        # S0=self._so(T)
        # H_adimensional=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        # f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        # S=H_adimensional+f+log(P/101325)

        # return unidades.SpecificHeat(S0.JgK-R*S/self.M, "JgK")

    # def Hv_Lee_Kesler(self, T):
        # """Método alternativo para el cálculo del calor de vaporización haciendo uso de las propiedades críticas
        # Procedure API 7C1.16 Pag.680
        # Valor en J/mol"""
        # Pv=self.Pv_DIPPR(T)
        # Tr=T/self.Tc
        # Pr=Pv/self.Pc
        # H_adimensional_vapor=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 1)
        # H_adimensional_liquido=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 0)
        # return unidades.Enthalpy(R*self.Tc/self.M*(H_adimensional_vapor-H_adimensional_liquido), "Jg")


    # def RhoG_Lee_Kesler(self, T, P):
        # a, b=eos.SRK_lib(self, T)
        # Z_srk=eos.Z_Cubic_EoS(T, P, b, a, b, 0, b)
        # Vvo=Z_srk[0]*R_atml*T/P

        # vr0v, vrhv, vr0l, vrhl=eos.Lee_Kesler_lib(T/self.Tc, P/self.Pc.atm, fase=1, Vvo=Vvo)
        # z0v=P/self.Pc.atm*vr0v/T*self.Tc
        # zhv=P/self.Pc.atm*vrhv/T*self.Tc
        # z=z0v+self.f_acent/factor_acentrico_octano*(zhv-z0v)
        # return P/z/R_atml/T




b1=0.1181193, 0.2026579
b2=0.265728, 0.331511
b3=0.154790, 0.027655
b4=0.030323, 0.203488
c1=0.0236744, 0.0313385
c2=0.0186984, 0.0503618
c3=0.0, 0.016901
c4=0.042724, 0.041577
d1=0.155488e-4, 0.48736e-4
d2=0.623689e-4, 0.0740336e-4
beta=0.65392, 1.226
gamma=0.060167, 0.03754

def Lee_Kesler_lib(Tr, Pr, fase=2, Vvo=0.0001, Vlo=5):
    """Librería para el cálculo de la EoS de Lee-Kesler
    Procedure API 6B1.8 pag 518
    Perry pag 2-358
    fase: fase para la que se realiza el cálculo
        0   -   Liquido
        1   -   Vapor
        2   -   Ambas"""

    Bo=b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
    Co=c1[0]-c2[0]/Tr+c3[0]/Tr**3
    Do=d1[0]+d2[0]/Tr
    def Vr(V):
        return 1+Bo/V+Co/V**2+Do/V**5+c4[0]/Tr**3/V**2*(beta[0]+gamma[0]/V**2)*exp(-gamma[0]/V**2)-Pr*V/Tr

    Bh=b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
    Ch=c1[1]-c2[1]/Tr+c3[1]/Tr**3
    Dh=d1[1]+d2[1]/Tr
    def Vrh(V):
        return 1+Bh/V+Ch/V**2+Dh/V**5+c4[1]/Tr**3/V**2*(beta[1]+gamma[1]/V**2)*exp(-gamma[1]/V**2)-Pr*V/Tr

    vr0v=vr0l=vrhv=vrhl=None
    if fase!=0:
        vr0v=fsolve(Vr, Vvo)
        vrhv=fsolve(Vrh, Vvo)
    elif fase!=1:
        vr0l=fsolve(Vr, Vlo)
        vrhl=fsolve(Vrh, Vlo)

    return vr0v, vrhv, vr0l, vrhl


def Z_Lee_Kesler(T, P, mezcla):
    """Factor de compresibilidad según la ecuación de estado de Lee-Kesler"""
    a, b, ai, bi, tdadt=SRK_lib(mezcla)
    Z_srk=Z_Cubic_EoS(b, a, b, 0, b)
#    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z_srk, b, a, ai, bi, 1, 0)

    Vvo=Z_srk[0]*R_atml*self.T/self.P.atm
    Vlo=Z_srk[1]*R_atml*self.T/self.P.atm

    Tmc, Pmc, f_acent, Vmc=Mix_Lee_Kesler(mezcla.fraccion, mezcla.componentes, mezcla.f_acent)
    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(T/Tmc, P/Pmc, Vvo, Vlo)
    z0l=self.P.atm/Pmc*vr0l/self.T*Tmc
    zhl=self.P.atm/Pmc*vrhl/self.T*Tmc
    z0v=self.P.atm/Pmc*vr0v/self.T*Tmc
    zhv=self.P.atm/Pmc*vrhv/self.T*Tmc

    return z0v+f_acent/factor_acentrico_octano*(zhv-z0v), z0l+f_acent/factor_acentrico_octano*(zhl-z0l)


def Lee_Kesler_lib_Cp(Tr, Pr, fase=1):
    """Librería para el cálculo de capacidades calorificas, usada a continuación en diferentes funciones
    Procedure API 7E1.6 Pag.726"""
    #FIXME: No concuerdan mucho los valores de cp y cv con los valores por DIPPR
    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(Tr, Pr, fase)
    if fase:
        E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        Cv0=-2*(b3[0]+3*b4[0]/Tr)/Tr**2/vr0v+3*c3[0]/Tr**3/vr0v**2+6*E
        E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Cvh=-2*(b3[1]+3*b4[1]/Tr)/Tr**2/vrhv+3*c3[1]/Tr**3/vrhv**2+6*E
        vr0=vr0v
        vrh=vrhv
    else:
        E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        Cv0=-2*(b3[0]+3*b4[0]/Tr)/Tr**2/vr0l+3*c3[0]/Tr**3/vr0l**2+6*E
        E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Cvh=-2*(b3[1]+3*b4[1]/Tr)/Tr**2/vrhl+3*c3[1]/Tr**3/vrhl**2+6*E
        vr0=vr0l
        vrh=vrhl

    return Cv0, Cvh, vr0, vrh


def Lee_Kesler_Entalpia_lib(Tr, Pr, w, fase=1):
    """Librería para el cálculo del factor adimensional de influencia de la presión sobre la temperatura.
    Usado en diversos métodos a continuación
    eq 7B3.7-1 pag 643"""

    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(Tr, Pr, fase)
    if fase:
        z0=Pr*vr0v/Tr
        zh=Pr*vrhv/Tr
        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        H0=-Tr*(z0-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0v-c2[0]/Tr/2/vr0v**2+d2[0]/5/Tr/vr0v**5+3*E)

        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Hh=-Tr*(zh-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhv-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhv**2+d2[1]/5/Tr/vrhv**5+3*E)
        H=H0+w/factor_acentrico_octano*(Hh-H0)
    else:
        z0=Pr*vr0l/Tr
        zh=Pr*vrhl/Tr
        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        H0=-Tr*(z0-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0l-c2[0]/Tr/2/vr0l**2+d2[0]/5/Tr/vr0l**5+3*E)

        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Hh=-Tr*(zh-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhl-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhl**2+d2[1]/5/Tr/vrhl**5+3*E)
        H=H0+w/factor_acentrico_octano*(Hh-H0)

    return H

def Lee_Kesler_Fugacidad_lib(Tr, Pr, w, fase=1):
    """Librería para el cálculo de la fugacidad, entropia...
    Procedure API 7G1.8 Pag.752"""
    vr0, vrh=self.Lee_Kesler_lib(Tr, Pr, fase)
    z0=Pr*vr0/Tr
    B=b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
    C=c1[0]-c2[0]/Tr+c3[0]/Tr**3
    D=d1[0]+d2[0]/Tr
    E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0**2)*exp(-gamma[0]/vr0**2))
    f0=z0-1-log(z0)+B/vr0+C/2/vr0**2+D/5/vr0**5+E

    zh=Pr*vrh/Tr
    B=b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
    C=c1[1]-c2[1]/Tr+c3[1]/Tr**3
    D=d1[1]+d2[1]/Tr
    E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrh**2)*exp(-gamma[1]/vrh**2))
    fh=zh-1-log(zh)+B/vrh+C/2/vrh**2+D/5/vrh**5+E

    return f0+w/factor_acentrico_octano*(fh-f0)

def Entalpia_Lee_Kesler(self):
    """Librería para el cálculo del factor adimensional de influencia de la presión sobre la temperatura.
    Usado en diversos métodos a continuación
    eq 7B3.7-1 pag 643"""
    Tmc, Pmc, f_acent, Vmc=self.mezcla.Mix_Lee_Kesler()

    a, b, ai, bi, tdadt=self.SRK_lib()
    Z_srk=self.Z_Cubic_EoS(b, a, b, 0, b)
    Vvo=Z_srk[0]*R_atml*self.T/self.P.atm
    Vlo=Z_srk[1]*R_atml*self.T/self.P.atm

    vr0, vrh, vr0l, vrhl=self.Lee_Kesler_lib(Tmc, Pmc, Vmc, Vvo, Vlo)
    Tr=self.T/Tmc
    Pr=self.P.atm/Pmc
    z0=Pr*vr0/self.T*Tmc
    E=0.041724/(2*Tr**3*0.060167)*(0.65392+1-(0.65392+1+0.060167/vr0**2)*exp(-0.060167/vr0**2))
    H0=-Tr*(z0-1-(0.2658728+2*0.154790/Tr+3*0.030323/Tr**2)/Tr/vr0-0.0186984/Tr/2/vr0**2+0.623689e-4/5/Tr/vr0**5+3*E)

    zh=Pr*vrh/Tr
    E=0.041577/(2*Tr**3*0.03754)*(1.226+1-(1.226+1+0.03754/vrh**2)*exp(-0.03754/vrh**2))
    Hh=-Tr*(zh-1-(0.331511+2*0.027655/Tr+3*0.203488/Tr**2)/Tr/vrh-(0.0503618-3*0.016901/Tr**2)/Tr/2/vrh**2+0.0740336e-4/5/Tr/vrh**5+3*E)
    return H0+f_acent/factor_acentrico_octano*(Hh-H0)


_all = [Lee_Kesler]

if __name__ == "__main__":
    from lib.corriente import Mezcla
    from lib import unidades
    # mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    # for T in [125, 135, 145, 165, 185, 205]:
        # eq = Lee_Kesler(T, 1, mezcla)
# #         print(eq.H_exc)
        # print(eq.Z)

    T = unidades.Temperature(120, "F")
    P = unidades.Pressure(485, "psi")
    mix = Mezcla(5, caudalMolar=5597, ids=[1, 2, 40, 41],
                 fraccionMolar=[0.3177, 0.5894, 0.0715, 0.0214])
    eq = Lee_Kesler(T, P, mix)
    print(eq.x, eq.xi, eq.yi)

    # mix = Mezcla(5, caudalMolar=5597, ids=[22],
                 # fraccionMolar=[1])
    # eq = Lee_Kesler(40+273, 9e6, mix)
    # # print(eq.x, eq.xi, eq.yi)
    # print(eq.Z)
