#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from math import exp, log

from scipy.optimize import fsolve
from scipy import constants as k

from lib import unidades
from lib.eos import EoS
from lib.bip import Kij


class BWRSoave(EoS):
    r"""Soave modification of Benedict-Webb-Rubin equation of state as give in
    [1]_

    .. math::
        Z = 1 + B\rho + D\rho^4 + E\rho^2\left(1+F\rho^2\right)
        \exp(-F\rho^2)

    .. math::
        \begin{array}[t]{l}
        \beta = B \left(\frac{P_c}{RT_c}\right) = \beta_c +
        0.422\left(1-\frac{1}{T_r^{1.6}}\right) +
        0.234 \omega \left(1-\frac{1}{T_r^3}\right) \\
        \delta = D \left(\frac{P_c}{RT_c}\right)^4 = \delta_c \left[
        1 + d_1\left(\frac{1}{T_r}-1\right) +
        d_2\left(\frac{1}{T_r}-1\right)^2\right] \\
        \epsilon = E \left(\frac{P_c}{RT_c}\right)^2 = \epsilon_c +
        e_1\left(\frac{1}{T_r}-1\right) +
        e_2\left(\frac{1}{T_r}-1\right)^2 +
        e_3\left(\frac{1}{T_r}-1\right)^3 \\
        \phi = F \left(\frac{P_c}{RT_c}\right)^2 = f Z_c^2
        \end{array}

    where:

    .. math::
        \begin{array}[t]{l}
        d_1 = 0.4912 + 0.6478\omega\\
        d_2 = 0.3 + 0.3619\omega\\
        e_1 = 0.0841 + 0.1318\omega + 0.0018\omega^2\\
        e_2 = 0.075 + 0.2408\omega + 0.014\omega^2\\
        e_3 = -0.0065 + 0.1798\omega + 0.0078\omega^2\\
        f = 0.77\\
        \beta_c = bZ_c\\
        \delta_c = dZ_c^4\\
        \epsilon_c= eZ_c^2\\
        e = \frac{2-5Z_c}{\left(1+f+3f^2-2f^3\right)\exp(-f)}\\
        d = \frac{1}{3} \left[1-2Z_c-e\left(1+f-2f^2\right)\exp(-f)\right]\\
        b = Z_c - 1 - d - e\left(1+f\right)\exp(-f)
        \end{array}

    This equation use critical parameters as parameters including Zc and
    acentric factor.

    The mixture are treated as fictitious pure compounds with the critical
    parameters calculated by this mixing rules:

    .. math::
        \begin{array}[t]{l}
        T_{cm} = \frac{S_1}{\left(\sqrt{S_2}+\sqrt{S_3}\right)^2} \\
        P_{cm} = \frac{T_{cm}}{S_3} \\
        m_m = \sqrt{\frac{S_2}{S_3}} \\
        Z_{cm} = \frac{\sum_ix_iZ_{ci}T{_ci}/P_{ci}}{\sum_ix_i T_{ci}/P_{ci}}\\
        S_1 = \sum_i \sum_j x_ix_j\left(1-k_{ij}\right)\frac{T_{ci}T_{cj}}
        {\sqrt{P_{ci}P_{cj}}} \left(1+m_i\right)\left(1+m_j\right) \\
        S_2 = \sum_i \sum_j x_ix_j\left(1-k_{ij}\right)\frac{T_{ci}T_{cj}}
        {\sqrt{P_{ci}P_{cj}}} m_i m_j \\
        S_3 = \sum_i x_i \frac{T_{ci}}{P_{ci}} \\
        \end{array}

    """
    __title__ = "Benedict-Webb-Rubin-Soave  (1999)"
    __status__ = "BWRSoave"

    __doi__ = (
      {
        "autor": "Soave, G.S.",
        "title": "An Effective Modification of the Benedict-Webb-Rubin "
                 "Equation of State",
        "ref": "Fluid Phase Equilibria 164(2) (1999) 157-172",
        "doi": "10.1016/s0378-3812(99)00252-6"},
      {
        "autor": "Soave, G.S.",
        "title": "A Noncubic Equation of State for the Tretament of "
                 "Hydrocarbon Fluids at Rerservoir Conditions",
        "ref": "Ind. Eng. Chem. Res. 34(11) (1995) 3981-3994",
        "doi": "10.1021/ie00038a039"},
      {
        "autor": "",
        "title": "",
        "ref": "",
        "doi": ""},
    )

    def __init__(self, T, P, mezcla):
        EoS.__init__(self, T, P, mezcla)

        self.kij = Kij(mezcla.ids, "BWRSoave")

        # print(self.mezcla.Tc, T)
        # if self.mezcla.Tc < T:
            # self.x = 1
            # self.xi = self.zi
            # self.yi = self.zi
            # self.Zg = self._Z()
            # self.Zl = None

        # else:
        self.x, self.Zl, self.Zg, self.xi, self.yi, self.Ki = self._Flash()
            # # print("q = ", self.x)
            # # print("x = ", self.xi)
            # # print("y = ", self.yi)
            # # print("K = ", self.Ki)

        if self.Zl:
            self.Vl = unidades.MolarVolume(self.Zl*k.R*T/P, "m3mol")   # l/mol
            self.rhoL = unidades.MolarDensity(1/self.Vl)
        else:
            self.Vl = None
            self.rhoL = None

        if self.Zg:
            self.Vg = unidades.MolarVolume(self.Zg*k.R*T/P, "m3mol")  # l/mol
            self.rhoG = unidades.MolarDensity(1/self.Vg)

        else:
            self.Vg = None
            self.rhoG = None

    def _lib(self, T, **kw):
        """Library for compound specified parameters"""

        Tc = kw["Tc"]
        w = kw["w"]
        Zc = kw["Zc"]

        Tr = T/Tc
        tau = 1/Tr-1

        f = 0.77                                                        # Eq 12
        e = (2 - 5*Zc)/((1 + f + 3*f**2 - 2*f**3)*exp(-f))              # Eq 14
        d = (1 - 2*Zc - e*(1 + f - 2*f**2)*exp(-f))/3
        b = Zc - 1 - d - e*(1+f)*exp(-f)
        beta_c = b*Zc                                                   # Eq 13
        delta_c = d*Zc**4
        epsilon_c = e*Zc**2

        d1 = 0.4912 + 0.6478*w                                          # Eq 10
        d2 = 0.3 + 0.3619*w
        e1 = 0.0841 + 0.1318*w + 0.0018*w**2                            # Eq 11
        e2 = 0.0750 + 0.2408*w - 0.0140*w**2
        e3 = -0.0065 + 0.1798*w - 0.0078*w**2

        # Eq 6
        beta = beta_c + 0.422*(1-1/Tr**1.6) + 0.234*w*(1-1/Tr**3)

        # Eq 7
        delta = delta_c*(1 + d1*tau + d2*tau**2)

        # Eq 8
        epsilon = epsilon_c + e1*tau + e2*tau**2 + e3*tau**3

        # Eq 9
        phi = f*Zc**2

        kw = {}
        kw["beta"] = beta
        kw["delta"] = delta
        kw["epsilon"] = epsilon
        kw["phi"] = phi
        return kw

    def _mixture(self, zi):
        """Mixing rules"""
        # Linear correlation for acentric factor
        m = [1.25*cmp.f_acent for cmp in self.componente]

        S1 = S2 = S3 = 0
        for ci, kiji, mi, xi in zip(self.componente, self.kij, m, zi):
            for cj, ki, mj, xj in zip(self.componente, kiji, m, zi):
                S1 += xi*xj*(1-ki)*ci.Tc*cj.Tc/(ci.Pc*cj.Pc)**0.5*(1+mi)*(1+mj)
                S2 += xi*xj*(1-ki)*(ci.Tc*cj.Tc/ci.Pc/cj.Pc)**0.5*mi*mj
            S3 += xi*ci.Tc/ci.Pc

        Tcm = S1/(S2**0.5+S3**0.5)**2
        Pcm = Tcm/S3
        mm = (S2/S3)**0.5

        # Eq 23
        num = 0
        den = 0
        for xi, cmp in zip(zi, self.componente):
            num += xi*cmp.Zc*cmp.Tc/cmp.Pc
            den += xi*cmp.Tc/cmp.Pc
        Zcm = num/den

        kw = {}
        kw["Tc"] = Tcm
        kw["Pc"] = Pcm
        kw["w"] = mm/1.25
        kw["Zc"] = Zcm
        kw["S1"] = S1
        kw["S2"] = S2
        kw["S3"] = S3
        return kw

    def _Z(self, zi=None, T=None, P=None, crit=None, phase=0):

        if not zi:
            zi = self.zi
        if not T:
            T = self.T
        if not P:
            P = self.P

        if not crit:
            crit = self._mixture(zi)
        kw = self._lib(T, **crit)

        b = kw["beta"]
        d = kw["delta"]
        e = kw["epsilon"]
        ph = kw["phi"]
        Tr = T/crit["Tc"]
        Pr = P/crit["Pc"]

        def func(ps):
            rlh = ps*(1+b*ps+d*ps**4+e*ps**2*(1+ph*ps**2)*exp(-ph*ps**2))
            return rlh-Pr/Tr

        if Tr > 1:
            y0 = Pr/Tr/(1+b*Pr/Tr)
        elif phase == 0:
            y0 = Pr/Tr/(1+b*Pr/Tr)
        else:
            y0 = 1/crit["Zc"]**(1+(1-Tr)**(2/7))

        rinput = fsolve(func, y0, full_output=True)
        y = rinput[0]
        Z = Pr/y/Tr
        # print(Z, rinput)

        return Z

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
        critL = self._mixture(xi)
        Zl = self._Z(xi, T, P, crit=critL, phase=1)[0]
        kwL = self._lib(T, **critL)
        kwL.update(critL)
        tital = self._fugacity(Zl, xi, **kwL)

        critG = self._mixture(yi)
        Zv = self._Z(yi, T, P, crit=critG, phase=0)[0]
        kwG = self._lib(T, **critG)
        kwG.update(critG)
        titav = self._fugacity(Zv, yi, **kwG)
        return tital, titav

    def _fugacity(self, Z, x, **kw):
        # Fugacity coefficient in a mixture, Appendix A
        print(Z, x)
        f = 0.77
        Tc = kw["Tc"]
        Pc = kw["Pc"]
        Zc = kw["Zc"]
        w = kw["w"]
        beta = kw["beta"]
        delta = kw["delta"]
        epsilon = kw["epsilon"]
        phi = kw["phi"]
        S1 = kw["S1"]
        S2 = kw["S2"]
        S3 = kw["S3"]

        Tr = self.T/Tc
        Pr = self.P/Pc
        psi = 1/Z*Pr/Tr

        f = 0.77
        e = (2 - 5*Zc)/((1 + f + 3*f**2 - 2*f**3)*exp(-f))
        d = (1 - 2*Zc - e*(1+f-2*f**2)*exp(-f))/3
        b = Zc - 1 - d - e*(1+f)*exp(-f)
        # betac = b*Zc
        deltac = d*Zc**4
        # epsilonc = e*Zc**2

        d1 = 0.4912 + 0.6478*w
        d2 = 0.3 + 0.3619*w
        e1 = 0.0841 + 0.1318*w + 0.0018*w**2
        e2 = 0.0750 + 0.2408*w - 0.0140*w**2
        e3 = -0.0065 + 0.1798*w - 0.0078*w**2

        dZcdn = []
        for cmp in self.componente:
            dZcdn.append(cmp.Tc/cmp.Pc/Tc*Pc*(cmp.Zc-Zc))

        # Appendix B
        S1i = []
        S2i = []
        for cmp in self.componente:
            S1i.append(cmp.Tc**2/cmp.Pc*(1+1.4*cmp.f_acent)**2)
            S2i.append(cmp.Tc/cmp.Pc*(1.4*cmp.f_acent)**2)

        S1ij = []
        S2ij = []
        for s1i, s2i, kiji in zip(S1i, S2i, self.kij):
            S1iji = []
            S2iji = []
            for s1j, s2j, kij in zip(S1i, S2i, kiji):
                S1iji.append((1-kij)*(s1i*s1j)**0.5)
                S2iji.append((1-kij)*(s2i*s2j)**0.5)
            S1ij.append(S1iji)
            S2ij.append(S2iji)

        from pprint import pprint
        pprint(S1ij)
        pprint(S1i)
        sum1 = [xi*s1ij[i] for i, (xi, s1ij) in enumerate(zip(x, S1ij))]
        sum2 = [xi*s2ij[i] for i, (xi, s2ij) in enumerate(zip(x, S2ij))]
        print(sum1)
        S1m = sum([x * s1 for x, s1 in zip(self.zi, sum1)])
        S2m = sum([x * s2 for x, s2 in zip(self.zi, sum2)])

        # Appendix C, section B1, from [2]_
        ds1dn = [2*s1i - 2*S1m for s1i in sum1]
        ds2dn = [2*s2i - 2*S2m for s2i in sum2]
        ds3dn = [cmp.Tc/cmp.Pc-Tc/Pc for cmp in self.componente]

        # Appendix C, section B3, from [2]_
        dwdn = []
        for ds2, ds3 in zip(ds2dn, ds3dn):
            dwdn.append(w/2*(ds2/S2-ds3/S3))

        # Appendix C, section B2, from [2]_
        dTcdn = []
        for ds1, ds2, ds3 in zip(ds1dn, ds2dn, ds3dn):
            dTcdn.append(ds1/S1-(ds2/S2**0.5+ds3/S3**0.5)/(S2**0.5+S3**0.5))

        dedn = [-5*dZ/(1+f+3*f**2-2*f**3)/exp(-f) for dZ in dZcdn]
        dddn = [(1-Z-e*(1+f-2*f**2)*exp(-f))/3 for Z, e in zip(dZcdn, dedn)]
        dbdn = [Z-d-e*(1+4)*exp(-f) for Z, d, e in zip(dZcdn, dddn, dedn)]
        dbetacdn = [Zc*db + b*dZ for db, dZ in zip(dbdn, dZcdn)]
        ddeltacdn = [Zc**4*dd + 4*d*Zc**3*dZ for dd, dZ in zip(dddn, dZcdn)]
        depscdn = [Zc**2*de + 2*e*Zc*dZ for de, dZ in zip(dedn, dZcdn)]

        dbbcdTc = -1.6*0.422/Tr**1.6 - w*3*0.234/Tr**3

        dbbcdw = 0.422/Tr**0.234

        dbetadn = []
        for dbc, dtc, dwn in zip(dbetacdn, dTcdn, dwdn):
            dbetadn.append(dbc + dbbcdTc*dtc + dbbcdw*dwn)

        dddcdTc = (d1+d2*(1/Tr-1))/Tr
        dddcdw = 0.6478*(1/Tr-1) + 0.3619*(1/Tr-1)**2

        ddeltadn = []
        for dd, dt, dw in zip(ddeltacdn, dTcdn, dwdn):
            ddeltadn.append(delta/deltac*dd + deltac*(dddcdTc*dt + dddcdw*dw))

        deTcterm = (e1 + 2*e2*(1/Tr-1) + 3*e3*(1/Tr-1)**2)/Tr
        depsdw = (0.1318+2*0.0018*w)*(1/Tr-1) + \
            (0.2408-2*0.014*w)*(1/Tr-1)**2 + (0.1798-2*0.0078*w)*(1/Tr-1)**3

        depsilondn = []
        for dec, dtc, dwn in zip(depscdn, dTcdn, dwdn):
            depsilondn.append(dec + deTcterm*dtc + depsdw*dwn)

        dphidn = [2*f*Zc*dZn for dZn in dZcdn]

        dBn2 = []
        for dbet, cmp in zip(dbetadn, self.componente):
            dBn2.append(1 + dbet/beta + cmp.Tc/cmp.Pc/Tc*Pc)

        dDn5 = []
        for ddel, cmp in zip(ddeltadn, self.componente):
            dDn5.append(1 + ddel/delta + 4*cmp.Tc/cmp.Pc/Tc*Pc)

        dEn3 = []
        for dep, cmp in zip(depsilondn, self.componente):
            dEn3.append(dep/epsilon + 2*cmp.Tc/cmp.Pc/Tc*Pc)

        dFn2 = []
        for dp, cmp in zip(dphidn, self.componente):
            dFn2.append(dp/phi + 2*cmp.Tc/cmp.Pc/Tc*Pc)

        tita = []
        for db, dd, de, df in zip(dBn2, dDn5, dEn3, dFn2):
            fug = -log(Z) + db*beta*psi + dd*delta*psi**4/4 + \
                de*epsilon/phi*(1-(1+phi*psi**2/2)*exp(-phi*psi**2)) + \
                df*epsilon/phi*(
                    (1+phi*psi**2+phi**2*psi**4/2) * exp(-phi*psi**2)-1)
            tita.append(fug)

        return tita


_all = [BWRSoave]


if __name__ == "__main__":
    from lib.mezcla import Mezcla

    from lib.EoS.BWRS import BWRS
    # from lib.mEoS import O2
    # mezcla = Mezcla(3, ids=[47], caudalMasico=1, fraccionMolar=[1])
    # eq = BWRSoave(160, 1e7, mezcla)
    # # print(eq.Zg)
    # eq = BWRS(160, 1e7, mezcla)
    # print("BWRS", eq.Zg, eq.rhoG)
    # eq = O2(T=160, P=1e7)
    # print("MEoS", eq.Z, eq.rhoM)

    # from matplotlib.pyplot import *
    # from numpy import logspace
    # P = logspace(1, 3.5, 200)
    # Z1 = []
    # Z2 = []
    # for p in P:
        # eq = BWRSoave(160, p*1e5, mezcla)
        # if eq.x == 1:
            # Z1.append(eq.Zg)
        # else:
            # Z1.append(eq.Zl)
        # eq = BWRS(160, p*1e5, mezcla)
        # # eq = O2(T=160, P=p*1e5)
        # # Z2.appendg(eq.Z)
        # if eq.x == 1:
            # Z2.append(eq.Zg)
        # else:
            # Z2.append(eq.Zl)

    # plot(P, Z1)
    # plot(P, Z2)
    # yscale("log")
    # xscale("log")
    # xlim(10, 5000)
    # ylim(0.1, 10)
    # show()

    mezcla = Mezcla(3, ids=[2, 3, 14], caudalMasico=1,
                    fraccionMolar=[0.698, 0.297, 0.005])
    eq = BWRSoave(230, 4e6, mezcla)
    print(eq.rhoL.str)
