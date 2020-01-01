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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Library with the implementantion of a generic cubic equation of state with the
form

.. math::
    P = \frac{RT}{V-b}-\frac{\alpha(T)}{V^2+\delta V+\epsilon}

Expressing as a cubic polynomy in compressibility factor easy to solve

.. math::
    Z^3 + \left(\delta'-B'-1\right)Z^2 +
    \left(a'+\epsilon'-\delta'\left(b'+1\right)\right)Z -
    \left(\epsilon'\left(b'+1\right)+a'b'\right) = 0

using the adimensional parameters

.. math::
    \begin{array}[t]{l}
    a' = \frac{aP}{RT}\\
    b' = \frac{bP}{RT}\\
    \delta' = \frac{\delta P}{RT}\\
    \epsilon' = \frac{\epsilon P}{RT}\\
    \end{array}

Each cubic EoS implemented here would a specific form of this general
expression changing the values of δ, ε and the expresion of α(T)

Each equation is specially suitable for different compounds, for example, the
Schmidt-Wenzel (SW) equation (1980) and the Adachi-Lu-Sugie (ALS) equation
(1983) are good for methane to n-decane. The Yu-Lu (YL) equation (1987) was
designed for asymmetric nonpolar mixtures, but not for polar substances. The
Iwai-Margerum-Lu (IML) equation ( 1987) was developed for polar substances, but
not suitable for nonpolar substances with large molecular weight.
"""


from scipy import roots, r_, log, exp, sqrt
from scipy.constants import R

from PyQt5.QtWidgets import QApplication

from lib import unidades
# from lib.corriente import Mezcla
from lib.eos import EoS
from lib.physics import R_atml
from lib.bip import Kij, Mixing_Rule
from lib.utilities import refDoc


# TODO: Add parameters, file
# Melhem, Almeida - A data Bank of Parameters for the Attractive-Aznar Telles
# self.Almeida = [0, 0]

# TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/alfas.pdf
# self.Mathias = 0
# self.Adachi = [0, 0]
# self.Andoulakis = [0, 0, 0]

__doi__ = {
    1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
    2:
        {"autor": "Ahmed, T.",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},
    3:
        {"autor": "Bell, I.H., Jäger, A.",
         "title": "Helmholtz Energy Transformations of Common Cubic Equations "
                  "of State for Use with Pure Fluids and Mixtures",
         "ref": "J. Res. of NIST 121 (2016) 236-263",
         "doi": "10.6028/jres.121.011"},

    4:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
        }


alfa = (QApplication.translate("pychemqt", "Original"),
        "Boston-Mathias",
        "Twu",
        "Doridon")


@refDoc(__doi__, [3])
def CubicHelmholtz(tau, delta, **kw):
    r"""Residual contribution to the free Helmholtz energy from a generic cubic
    equation of state with the form:

    .. math::
        P = \frac{RT}{V-b}-\frac{\alpha(T)}{\left(v+\Delta_1b\right)
        \left(v+\Delta_2b\right)}

    From this formulation it's possible calculate the Helmholtz free energy
    with the equation:

    .. math::
        \alpha^r = \phi^{(-)}-\frac{\tau\alpha}{RT_c}\phi^{(+)}

    Parameters
    ----------
    tau : float
        Inverse reduced temperature, Tc/T [-]
    delta : float
        Reduced density, rho/rhoc [-]
    kw : list
        Aditional parameters specific of cubic equation of state
        The parameters include: rhoc, Tc, b, alfa, Delta1, Delta2

    Returns
    -------
    prop : dictionary with residual adimensional helmholtz energy and deriv
        fir  [-]
        firt: [∂fir/∂τ]δ,x  [-]
        fird: [∂fir/∂δ]τ,x  [-]
        firtt: [∂²fir/∂τ²]δ,x  [-]
        firdt: [∂²fir/∂τ∂δ]x  [-]
        firdd: [∂²fir/∂δ²]τ,x  [-]
    """
    b = kw["b"]
    a = kw["a"]
    dat = kw["dat"]
    datt = kw["datt"]
    dattt = kw["dattt"]
    Delta1 = kw["Delta1"]
    Delta2 = kw["Delta2"]

    rhoc = kw.get("rhoc", 1)
    Tc = kw.get("Tc", 1)
    phi1 = -log(1-b*delta*rhoc)
    if Delta1 == Delta2:
        # Special case using the l'Hôpital's rule
        phi2 = rhoc*delta
    else:
        phi2 = log((Delta1*b*rhoc*delta+1)/(Delta2*b*rhoc*delta+1)) / \
            b/(Delta1-Delta2)

    phi1d = b*rhoc/(1-b*delta*rhoc)
    phi1dd = b**2*rhoc**2/(1-b*delta*rhoc)**2
    phi1ddd = 2*b**3*rhoc**3/(1-b*delta*rhoc)**3

    PI12 = (1+Delta1*b*rhoc*delta) * (1+Delta2*b*rhoc*delta)
    PI12d = b*rhoc * (2*Delta1*Delta2*b*delta*rhoc + Delta1 + Delta2)
    PI12dd = 2*Delta1*Delta2*b**2*rhoc**2

    phi2d = rhoc/PI12
    phi2dd = -rhoc*PI12d/PI12**2
    phi2ddd = rhoc*(-PI12*PI12dd+2*PI12d**2)/PI12**3

    fir = phi1 - tau*a/R/Tc*phi2
    fird = phi1d - tau*a/R/Tc*phi2d
    firdd = phi1dd - tau*a/R/Tc*phi2dd
    firddd = phi1ddd - tau*a/R/Tc*phi2ddd

    # Eq 32
    dtat = tau*dat + a
    dtatt = tau*datt + 2*dat
    dtattt = tau*dattt + 3*datt

    firt = -dtat/R/Tc * phi2
    firtt = -dtatt/R/Tc * phi2
    firttt = -dtattt/R/Tc * phi2
    firdt = -dtat/R/Tc * phi2d
    firddt = -dtat/R/Tc * phi2dd
    firdtt = -dtatt/R/Tc * phi2d

    prop = {}
    prop["fir"] = fir
    prop["fird"] = fird
    prop["firt"] = firt
    prop["firdd"] = firdd
    prop["firdt"] = firdt
    prop["firtt"] = firtt
    prop["firddd"] = firddd
    prop["firddt"] = firddt
    prop["firdtt"] = firdtt
    prop["firttt"] = firttt

    prop["B"] = 0
    prop["C"] = 0
    prop["D"] = 0

    return prop


@refDoc(__doi__, [1, 2])
class Cubic(EoS):
    r"""Class to implement the common functionality of cubic equation of state

    This class implement a general cubic equation of state in the form:

    .. math::
        P = \frac{RT}{V-b}-\frac{\alpha(T)}{V^2+\delta V+\epsilon}

    .. math::
        P = \frac{RT}{V-b}-\frac{\alpha(T)}{\left(V+\delta_1b\right)
        \left(V+\delta_2b\right)}

    .. math::
        \delta_1 = -\frac{\sqrt{\delta^2-4\epsilon}-\delta}{2b}

    .. math::
        \delta_2 = -\frac{\sqrt{\delta^2-4\epsilon}+\delta}{2b}
    """

    def __init__(self, T, P, mezcla):
        P_atm = P/101325
        self.T = unidades.Temperature(T)
        self.P = unidades.Pressure(P)
        self.mezcla = mezcla
        self.componente = mezcla.componente
        self.zi = mezcla.fraccion

        self.B = self.b*P/R/T
        self.Tita = self.tita*P/(R*T)**2

        delta = self.delta*P/R/T
        epsilon = self.epsilon*(P/R/T)**2

        # δ1, δ2 calculated from polynomial factorization
        self.delta1 = ((self.delta**2-4*self.epsilon)**0.5-self.delta)/2/self.b
        self.delta2 = ((self.delta**2-4*self.epsilon)**0.5+self.delta)/2/self.b

        # Eq 4-6.3 in [1]_
        coeff = [1, delta-self.B-1, self.Tita+epsilon-delta*(self.B+1),
                 -epsilon*(self.B+1)-self.Tita*self.B]
        Z = roots(coeff)
        # print("Z", Z)
        # TODO: use the anallycal solution, Span, pag 50

        # Set the minimum and maximum root values as liquid and gas Z values
        self.Z = r_[Z[0].real, Z[2].real]
        self.Zl = min(Z).real
        self.Zg = max(Z).real

        self.V = self.Z*R_atml*T/P_atm  # l/mol
        self.rho = 1/self.V
        self.Vl = unidades.MolarVolume(self.Zl*R*T/P, "m3mol")   # l/mol
        self.Vg = unidades.MolarVolume(self.Zg*R*T/P, "m3mol")  # l/mol

        # tau = mezcla.Tc/T
        # delta = self.V[0]*mezcla.Vc
        # kw = {}
        # print(CubicHelmholtz(tau, delta, **kw))


        # print(coeff, Z)
        # self.x, self.xi, self.yi, self.Ki = self._Flash()

        # dep_v = self._departure(self.tita, self.b, self.delta, self.epsilon, self.dTitadT, self.V[0], T)
        # dep_l = self._departure(self.tita, self.b, self.delta, self.epsilon, self.dTitadT, self.V[1], T)

    def _mixture(self, eq, ids, par):
        self.kij = Kij(ids, eq)
        mixpar = Mixing_Rule(self.mezcla.fraccion, par, self.kij)
        return mixpar

    # def _PHIO(self, cp, Tc):
        # """Convert cp dict in phi0 dict when the cp expression isn't in
        # Helmholtz free energy terms"""
        # co = cp["ao"]-1
        # ti = []
        # ci = []
        # for n, t in zip(cp["an"], cp["pow"]):
            # ti.append(-t)
            # ci.append(-n/(t*(t+1))*Tc**t)

        # # The integration constant are difficult to precalculate as depend of
        # # resitual Helmholtz free energy. It's easier use a offset system
        # # saved the values in database and retrieve for each reference state
        # cI = 0
        # cII = 0

        # Fi0 = {"ao_log": [1,  co],
               # "pow": [0, 1] + ti,
               # "ao_pow": [cII, cI] + ci}

        # return Fi0

    # def _phi0i(self, cp, tau, delta):
        # r"""Ideal gas Helmholtz free energy and derivatives
        # The ideal gas specific heat can have different contributions

        # .. math::
            # \frac{C_p^o}{R} = c_o + \sum_i c_iT_r^i

        # The dict with the definition of ideal gas specific heat must define
        # the parameters:

            # * ao: Independent of temperature coefficient
            # * an: Polynomial term coefficient
            # * pow: Polynomial term temperature exponent

        # Parameters
        # ----------
        # cp : dict
            # Ideal gas properties parameters, can be in Cp term of directly in
            # helmholtz free energy
        # tau : float
            # Inverse reduced temperature, Tc/T [-]
        # delta : float
            # Reduced density, rho/rhoc [-]

        # Returns
        # -------
        # prop : dictionary with ideal adimensional helmholtz energy and deriv
            # fio  [-]
            # fiot: [∂fio/∂τ]δ  [-]
            # fiod: [∂fio/∂δ]τ  [-]
            # fiott: [∂²fio/∂τ²]δ  [-]
            # fiodt: [∂²fio/∂τ∂δ]  [-]
            # fiodd: [∂²fio/∂δ²]τ  [-]
        # """

        # Fi0 = self._PHIO(cp)

        # fio = Fi0["ao_log"][1]*log(tau)
        # fiot = Fi0["ao_log"][1]/tau
        # fiott = -Fi0["ao_log"][1]/tau**2

        # if delta:
            # fiod = 1/delta
            # fiodd = -1/delta**2
        # else:
            # fiod, fiodd = 0, 0
        # fiodt = 0

        # for n, t in zip(Fi0["ao_pow"], Fi0["pow"]):
            # fio += n*tau**t
            # if t != 0:
                # fiot += t*n*tau**(t-1)
            # if t not in [0, 1]:
                # fiott += n*t*(t-1)*tau**(t-2)

        # prop = {}
        # prop["fio"] = fio
        # prop["fiot"] = fiot
        # prop["fiott"] = fiott
        # prop["fiod"] = fiod
        # prop["fiodd"] = fiodd
        # prop["fiodt"] = fiodt
        # return prop

    # def _phi0(self, tau, delta):
        # fio = fiod = fiodd = fiot = fiott = fiodt = 0
        # for cmp, x in zip(self.mezcla.componente, self.mezcla.fraccion):
            # phii = self._phi0i(cp, tau, delta)
            # fio += x*phii["fio"] + x*log(x)
            # fiod += x*phii["fiod"] + x*log(x)
            # fiodd += x*phii["fiodd"] + x*log(x)
            # fiot += x*phii["fiot"] + x*log(x)
            # fiott += x*phii["fiott"] + x*log(x)
            # fiodt += x*phii["fiodt"] + x*log(x)

        # prop = {}
        # prop["fio"] = fio
        # prop["fiot"] = fiot
        # prop["fiott"] = fiott
        # prop["fiod"] = fiod
        # prop["fiodd"] = fiodd
        # prop["fiodt"] = fiodt
        # return prop

    # def _ideal(self, tau, delta):
        # fio = self._phi0(tau, delta)
        # p = {}
        # p["U"] = fio["fiot"]
        # p["H"] = 1 + fio["fiot"]
        # p["S"] = fio["fiot"]-fio["fio"]
        # p["cv"] = -fio["fiott"]

        # return p

    def _excess(self, tau, delta):
        fir = self.fir["fir"]
        fird = self.fir["fird"]
        firt = self.fir["firt"]
        firtt = self.fir["firtt"]
        p = {}
        p["Z"] = 1 + delta*fird
        p["H"] = tau*firt + delta*fird
        p["S"] = tau*firt - fir
        p["cv"] = -tau**2*firtt

        return p

    def _departure(self, a, b, d, e, TdadT, V, T):
        """Calculate departure function, Table 6-3 from [1]"""
        Z = 1 + b/(V-b) - a*V/R_atml/T/(V**2+d*V+e)

        # Numerador and denominator used in several expression
        K = (d**2-4*e)**0.5
        num = 2*V + d - K
        den = 2*V + d + K
        kw = {}
        kw["Z"] = Z
        if K:
            kw["H"] = 1 - (a+TdadT)/R_atml/T/K*log(num/den) - Z
            kw["S"] = TdadT/R_atml/K*log(num/den) - log(Z*(1-b/V))
            kw["A"] = -a/R_atml/T/K*log(num/den) + log(Z*(1-b/V))
            kw["f"] = a/R_atml/T/K*log(num/den) - log(Z*(1-b/V)) - (1-Z)
        else:
            kw["H"] = 1 - Z
            kw["S"] = -log(Z*(1-b/V))
            kw["A"] = log(Z*(1-b/V))
            kw["f"] = -log(Z*(1-b/V)) - (1-Z)
        return kw

    def _fug(self, Z, xi):
        """Calculate partial fugacities coefficieint of components

        References
        ----------
        mollerup, Chap 2, pag 64 and so
        """
        V = Z*R_atml*self.T/self.P

        g = log(V-self.b) - log(V)                                     # Eq 61
        f = 1/R_atml/self.b/(self.delta1-self.delta2) * \
            log((V+self.delta1*self.b)/(V+self.delta2*self.b))         # Eq 62

        gB = -1/(V-self.b)                                             # Eq 80
        
        An = -g                                                        # Eq 75
        AB = -n*gB-D/self.T*fB                                         # Eq 78
        AD = -f/self.T                                                 # Eq 79

        # Ch.3, Eq 66
        dAni = An+AB*Bi+AD*Di

        # Ch.2, Eq 13
        fi = dAni - log(Z)
        return

        

    def _fug(self, Z, xi):
        Ai=[]
        for i in range(len(self.componente)):
            suma=0
            for j in range(len(self.componente)):
                suma+=self.zi[j]*self.ai[j]**0.5*(1-self.kij[i][j])
            Ai.append(1/self.tita*2*self.ai[i]**0.5*suma)
        tita=[]
        for i in range(len(self.componente)):
            tita.append(exp(self.bi[i]/self.b*(Z-1)-log(Z-self.B)-self.Tita/self.B/sqrt(self.u**2-4*self.w)*(Ai[i]-self.bi[i]/self.b)*log((Z+self.B/2*(self.u+sqrt(self.u**2-4*self.w)))/(Z+self.B/2*(self.u-sqrt(self.u**2-4*self.w))))).real)
        return tita

# if __name__ == "__main__":
    # # from lib.mezcla import Mezcla
    # # mix = Mezcla(1, ids=[62], caudalUnitarioMasico=[1.])
    # # # mix = Mezcla(tipo=5, caudalMolar=1, ids=[2, 47, 98], fraccionMolar=[0.5, 0.3, 0.2])
    # # eq = SRK(300, 101325, mix)
    # # print(1/eq.V)
    # # # for T in [125, 135, 145, 165, 185, 205]:
        # # # eq = SRK(T, 1, mezcla)
        # # # print(eq.H_exc)
    # # from lib.mEoS import H2O
    # # st = H2O(T=300, P=101325)
    # # print(st.Z, st.rho)
