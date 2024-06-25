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


from math import log, exp

from scipy.constants import R
from tools.qt import translate

from lib import unidades
from lib.eos import EoS
from lib.physics import R_atml, cubicCardano
from lib.bip import Kij, Mixing_Rule
from lib.utilities import refDoc


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


alfa = (translate("EoS", "Original"),
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
    R = kw["R"]

    # This parameters are necessary only for multicomponent mixtures to
    # calculate fugacity coefficient
    bi = kw.get("bi", None)
    daxi = kw.get("daxi", None)

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

    # Virial coefficient
    phi1d0 = b*rhoc
    phi2d0 = rhoc
    fird0 = phi1d0 - tau*a/R/Tc*phi2d0
    phi1dd0 = b**2*rhoc**2
    PI12d0 = b*rhoc * (Delta1 + Delta2)
    phi2dd0 = -rhoc*PI12d0
    firdd0 = phi1dd0 - tau*a/R/Tc*phi2dd0
    prop["B"] = fird0
    prop["C"] = firdd0
    prop["D"] = 0

    if bi:
        # Composition derivatives for fugacity coefficient calculation
        c = 1/b
        dbxi = bi                                                      # Eq 132
        A = log((delta*rhoc*b*Delta1+1)/(delta*rhoc*b*Delta2+1))       # Eq 103

        dAxi = [delta*rhoc*db*(Delta1-Delta2)/PI12 for db in dbxi]     # Eq 104
        dcxi = [-db/b**2 for db in dbxi]                               # Eq 107

        phi1xi = [delta*rhoc*db/(1-delta*rhoc*b) for db in dbxi]       # Eq 80

        # Eq 111
        phi2xi = [(A*dc + c*dA)/(Delta1-Delta2) for dc, dA in zip(dcxi, dAxi)]

        dtaxi = [tau*da for da in daxi]

        # Eq 77
        phirxi = []
        for dt, p1x, p2x in zip(dtaxi, phi1xi, phi2xi):
            phirxi.append(p1x - 1/R/Tc*(dt*phi2 + tau*a*p2x))

        prop["firxi"] = phirxi

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

    def __init__(self, T, P, mezcla, **kwargs):
        EoS.__init__(self, T, P, mezcla, **kwargs)

        if "R" in kwargs:
            self.R = kwargs["R"]
        else:
            self.R = R

        self._cubicDefinition(T)

#         if self.mezcla.Tc < T:
#             self.x = 1
#             self.xi = self.zi
#             self.yi = self.zi
#             self.Zg = self._Z(self.zi, T, P)[-1]
#             self.Zl = None

#         else:
        self.x, self.Zl, self.Zg, self.xi, self.yi, self.Ki = self._Flash()
            # print("q = ", self.x)
            # print("x = ", self.xi)
            # print("y = ", self.yi)
            # print("K = ", self.Ki)

        if self.Zl:
            self.Vl = unidades.MolarVolume(self.Zl*self.R*T/P, "m3mol")
            rhoL = self.P/self.Zl/self.R/self.T
            self.rhoL = unidades.MolarDensity(rhoL, "molm3")
        else:
            self.Vl = None
            self.rhoL = None

        if self.Zg:
            self.Vg = unidades.MolarVolume(self.Zg*self.R*T/P, "m3mol")
            rhoG = self.P/self.Zg/self.R/self.T
            self.rhoG = unidades.MolarDensity(rhoG, "molm3")
        else:
            self.Vg = None
            self.rhoG = None

        self._volumeCorrection()

        # tau = mezcla.Tc/T
        # delta = self.V[-1]*mezcla.Vc
        # kw = {}
        # print(CubicHelmholtz(tau, delta, **kw))

        # dep_v = self._departure(self.tita, self.b, self.delta, self.epsilon, self.dTitadT, self.V[-1], T)
        # dep_l = self._departure(self.tita, self.b, self.delta, self.epsilon, self.dTitadT, self.V[0], T)

        # rho = self.rhoG.molm3
        # from pprint import pprint
        # pprint(self._phir(self.T, rho, self.yi))

    def _volumeCorrection(self):
        """Apply volume correction to the rhoL property"""
        pass

    def _cubicDefinition(self, T):
        """Definition of individual component parameters of generalized cubic
        equation of state, its calculation don't depend composition"""
        # TODO: Split fixed paremeters calculation from temperature dependences
        # to speed up
        pass

    def _GEOS(self, xi):
        """Definition of parameters of generalized cubic equation of state,
        each child class must define in this procedure the values of mixture
        a, b, delta, epsilon. The returned values are not dimensionless.

        Parameters
        ----------
        xi : list
            Molar fraction of component in mixture, [-]

        Returns
        -------
        parameters : list
            Mixture parameters of equation, a, b, c, d
        """
        pass

    def _Z(self, xi, T, P):
        """Calculate root of cubic polynomial in terms of GCEoS as give in
        [1]_.

        Parameters
        ----------
        xi : list
            Molar fraction of component in mixture, [-]
        T : float
            Temperature, [K]
        P : float
            Pressure, [Pa]

        Returns
        -------
        Z : list
            List with real root of equation
        """

        self._cubicDefinition(T)
        tita, b, delta, epsilon = self._GEOS(xi)
        B = b*P/self.R/T
        A = tita*P/(self.R*T)**2

        D = delta*P/self.R/T
        E = epsilon*(P/self.R/T)**2

        # Eq 4-6.3 in [1]_
        # η by default set to b to reduce terms, if any equations need that
        # term redefine this procedure
        coeff = (1, D-B-1, A+E-D*(B+1), -E*(B+1)-A*B)
        Z = cubicCardano(*coeff)

        # Sort Z values, if typeerror is raise return is because there is
        # complex root, so return only the real root
        try:
            Z = sorted(map(float, Z))
        except TypeError:
            Z = Z[0:1]

        return Z

    def _fug(self, xi, yi, T, P):
        """Fugacities of component in mixture calculation

        Parameters
        ----------
        xi : list
            Molar fraction of component in liquid phase, [-]
        yi : list
            Molar fraction of component in vapor phase, [-]
        T : float
            Temperature, [K]
        P : float
            Pressure, [Pa]

        Returns
        -------
        tital : list
            List with liquid phase component fugacities
        titav : list
            List with vapour phase component fugacities
        """
        self._cubicDefinition(T)
        Bi = [bi*P/self.R/T for bi in self.bi]
        Ai = [ai*P/(self.R*T)**2 for ai in self.ai]

        al, bl, deltal, epsilonl = self._GEOS(xi)
        Bl = bl*P/self.R/T
        Al = al*P/(self.R*T)**2
        Zl = self._Z(xi, T, P)[0]
        tital = self._fugacity(Zl, xi, Al, Bl, Ai, Bi)

        Zv = self._Z(yi, T, P)[-1]
        av, bv, deltav, epsilonv = self._GEOS(yi)
        Bv = bv*P/self.R/T
        Av = av*P/(self.R*T)**2
        titav = self._fugacity(Zv, yi, Av, Bv, Ai, Bi)
        return tital, titav

    def _fugacity(self, Z, zi, A, B, Ai, Bi):
        """Fugacity for individual components in a mixture using the GEoS in
        the Schmidt-Wenzel formulation, so the subclass must define the
        parameters u and w in the EoS

        Any other subclass with different formulation must overwrite this
        method
        """
        # Precalculation of inner sum in equation
        aij = []
        for ai, kiji in zip(Ai, self.kij):
            suma = 0
            for xj, aj, kij in zip(zi, Ai, kiji):
                suma += xj*(1-kij)*(ai*aj)**0.5
            aij.append(suma)

        tita = []
        for bi, aai in zip(Bi, aij):
            rhs = bi/B*(Z-1) - log(Z-B) + A/B/(self.u-self.w)*(
                    bi/B-2/A*aai) * log((Z+self.u*B)/(Z+self.w*B))
            tita.append(exp(rhs))

        return tita

    def _mixture(self, eq, xi, par):
        """Apply mixing rules to individual parameters to get the mixture
        parameters for EoS

        Although it possible use any of available mixing rules, for now other
        properties calculation as fugacity helmholtz free energy are defined
        using the vdW mixing rules.

        Parameters
        ----------
        eq : str
            codename of equation, PR, SRK...
        xi : list
            Molar fraction of component, [-]
        par : list
            list with individual parameters of equation, [-]

        Returns
        -------
        mixpar : list
            List with mixture parameters, [-]
        """
        self.kij = Kij(self.mezcla.ids, eq)
        mixpar = Mixing_Rule(xi, par, self.kij)
        return mixpar

    def _Tr(self):
        """Definition of reducing parameters"""
        if len(self.mezcla.componente) > 1:
            # Mixture as one-fluid
            Tr = 1
            rhor = 1
        else:
            # Pure fluid
            Tr = self.mezcla.Tc
            rhor = 1/self.mezcla.Vc/1000  # m3/mol

        return Tr, rhor

    def _phir(self, T, rho, xi):

        Tr, rhor = self._Tr()

        tau = Tr/T
        delta = rho/rhor
        a, b, d, e = self._GEOS(xi)

        kw = self._da(tau, xi)

        Tr, rhor = self._Tr()
        kw["rhoc"] = rhor
        kw["Tc"] = Tr
        kw["Delta1"] = self.u
        kw["Delta2"] = self.w
        kw["bi"] = self.bi
        kw["b"] = b
        kw["a"] = a
        kw["R"] = self.R

        fir = CubicHelmholtz(tau, delta, **kw)
        # print(self._excess(tau, delta, fir))
        # print("fir: ", fir["fir"])
        # print("fird: ", fir["fird"]*delta)
        # print("firt: ", fir["firt"]*tau)
        # print("firdd: ", fir["firdd"]*delta**2)
        # print("firdt: ", fir["firdt"]*delta*tau)
        # print("firtt: ", fir["firtt"]*tau**2)
        # print("firddd: ", fir["firddd"]*delta**3)
        # print("firddt: ", fir["firddt"]*delta**2*tau)
        # print("firdtt: ", fir["firdtt"]*delta*tau**2)
        # print("firttt: ", fir["firttt"]*tau**3)
        # T = Tr/tau
        # rho = rhor*delta
        # print("P", (1+delta*fir["fird"])*R*T*rho)
        # print(delta, fir["fird"], R, T, rho)

        return fir

    def _excess(self, tau, delta, phir):
        fir = phir["fir"]
        fird = phir["fird"]
        firt = phir["firt"]
        firtt = phir["firtt"]
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

    # def _fug2(self, Z, xi):
        # """Calculate partial fugacities coefficieint of components

        # References
        # ----------
        # mollerup, Chap 2, pag 64 and so
        # """
        # V = Z*R_atml*self.T/self.P

        # g = log(V-self.b) - log(V)                                     # Eq 61
        # f = 1/R_atml/self.b/(self.delta1-self.delta2) * \
            # log((V+self.delta1*self.b)/(V+self.delta2*self.b))         # Eq 62

        # gB = -1/(V-self.b)                                             # Eq 80

        # An = -g                                                        # Eq 75
        # AB = -n*gB-D/self.T*fB                                         # Eq 78
        # AD = -f/self.T                                                 # Eq 79

        # # Ch.3, Eq 66
        # dAni = An+AB*Bi+AD*Di

        # # Ch.2, Eq 13
        # fi = dAni - log(Z)
        # return fi
