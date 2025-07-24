#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


from math import log
import time

from scipy.constants import R
from tools.qt import translate

from lib import unidades
from lib.compuestos import atomic_decomposition, facent_LeeKesler, RhoL_Rackett
from lib.elemental import databank
from lib.physics import R_atml, R_cal
from lib.utilities import refDoc


# Get molecular weight of atomic elements
MW = {}
databank.execute("SELECT symbol, atomic_mass FROM ELEMENTS")
for symbol, m in databank:
    MW[symbol] = m


__doi__ = {
    1:
        {"autor": "Hurst, J.E., Harrison, B.K.",
         "title": "Estimation of Liquid and Solid Heat Capacities Using a "
                  "Modified Kopp's Rule",
         "ref": "Chem. Eng. Comm. 112 (1992) 21-30",
         "doi": "10.1080/00986449208935989"}}


class newComponente():
    """Base class with general new component definition,
    interaction with database"""

    def calculo(self):
        """Calculate procedure with common functionality and define the
        properties don't defined by the method
        The child class must implement the specific calculate procedure and
        call this method it is necessary to finish definition"""

        if self.kwargs["name"]:
            self.name = str(self.kwargs["name"])
        else:
            self.name = self.__class__.__name__ + "_" + \
                time.strftime("%d/%m/%Y-%H:%M:%S")

    def export2Component(self):
        """Return the new component data as a list to add to database"""
        ele = []
        ele.append(self.formula)
        ele.append(self.name)
        ele.append(self.M)
        ele.append(self.Tc)
        ele.append(self.Pc.atm)
        ele.append(self.Vc)
        ele.append(self.API)
        ele.append(self.cp)

        # Parametrics
        ele.append([])
        ele.append([])
        if self.mul:
            ele.append(self.mul)
        else:
            ele.append([])
        ele.append([])

        # DIPPR
        ele.append([])
        ele.append([])
        ele.append([])
        ele.append([])
        ele.append([1, self.cps, 0, 0, 0, 0, 0, self.Tc])
        ele.append([1, self.cpl, 0, 0, 0, 0, 0, self.Tc])
        ele.append([])
        ele.append([])
        ele.append([])
        ele.append([])
        ele.append([])
        ele.append([])

        # Others
        ele.append(0)
        ele.append(self.Vliq)
        ele.append(self.rackett)
        ele.append(self.SG)
        ele.append(self.f_acent)
        ele.append(self.Parametro_solubilidad)
        ele.append(0)

        ele.append(0)
        ele.append(self.Tb)
        ele.append(self.Tf)
        ele.append("")
        ele.append("")

        ele.append("[]")

        ele.append(0)
        ele.append(0)
        ele.append(0)
        ele.append(0)
        ele.append(0)
        ele.append(self.Hf)
        ele.append(self.Gf)
        ele.append(0)
        ele.append(0)
        ele.append(0)

        ele.append("")

        ele.append(0)
        ele.append(0)
        ele.append(0)
        ele.append(0)

        ele.append("")

        # Add antoine and wagner parameters for vapor pressure
        ele.append([])
        ele.append([])
        return ele


class GroupContribution(newComponente):
    """Base class for group contribution methods"""

    kwargs = {"group": [],
              "contribution": [],
              "M": 0.0,
              "Tb": 0.0,
              "SG": 0.0,
              "name": ""}

    status = 0
    _bool = False
    msg = ""

    FirstOrder = 0
    SecondOrder = 0
    ThirdOrder = 0
    cp = []
    Tb = 0
    Tf = 0
    Hf = 0
    Gf = 0
    mul = None
    Hm = 0

    def __init__(self, **kwargs):
        """Constructor with kwargs of derived class"""
        self.kwargs = self.kwargs.copy()
        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        """Callable instance support"""
        self.kwargs.update(kwargs)
        self._bool = True
        if self.isCalculable():
            self.calculo()

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        self.group = self._group()
        if not self.kwargs["group"] or not self.kwargs["contribution"]:
            self.msg = translate("newComponent", "undefined group")
            self.status = 0
        else:
            self.status = 1
            self.msg = ""
            return True

    def calculo(self):
        """Calculate procedure with common functionality and define the
        properties don't defined by the method
        The child class must implement the specific calculate procedure and
        call this method it is necessary to finish definition"""

        if "Tb" not in self.__dict__:
            self.Tb = self._Tb()
        if "f_acent" not in self.__dict__:
            self.f_acent = facent_LeeKesler(self.Tb, self.Tc, self.Pc)
        if "Hv" not in self.__dict__:
            self.Hv = self._Calor_vaporizacion()

        self.rackett = self._Rackett()
        if "Vliq" not in self.__dict__:
            self.Vliq = self._VLiq()
        self.Parametro_solubilidad = self._SolubilityParameter()

        if self.kwargs["SG"]:
            self.SG = self.kwargs["SG"]
        else:
            self.SG = self._SG()
        self.Kw = self.Tb.R**(1./3)/self.SG
        if "Vc" not in self.__dict__:
            self.Vc = self._Vc()
        if "cp" not in self.__dict__:
            self.cp = self._cp()

        # Liquid and solid specific heat from Hurst correlation
        cps, cpl = cpLS_Hurst(self.group)
        self.cps = unidades.SpecificHeat(cps/self.M)
        self.cpl = unidades.SpecificHeat(cpl/self.M)

        self.API = 141.5/self.SG-131.5
        self.txt, self.formula = self.EmpiricFormula()
        newComponente.calculo(self)

    def __bool__(self):
        return self._bool

    def clear(self):
        self.kwargs = self.__class__.kwargs
        self.__dict__.clear()
        self._bool = False

    def _atoms(self):
        """Procedure to calculate the atom number of atoms in molecule"""
        a = 0
        for grp in self.group:
            for x in grp.values():
                a += x
        return a

    def _atomX(self, x):
        """Calculate the atom count of element X of the molecule"""
        a = 0
        for grp in self.group:
            a += grp.get(x, 0)
        return a

    def _M(self):
        """Procedure to calculate the molecular weight"""
        M = 0
        for grp in self.group:
            for ele, c in grp.items():
                M += c*MW[ele]
        return M

    def _group(self):
        """From group contribution desglose the chemical composition"""
        group = []
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            # Only the first order term count for this
            if i < self.FirstOrder:
                # Clean additional comment of group, ring flag and other,
                # separated of main group by spaces
                cmp = self._coeff["txt"][i][0]
                if " " in cmp:
                    cmp = cmp.split(" ")[0]

                grp = atomic_decomposition(cmp)
                for x in range(c):
                    group.append(grp)
        return group

    def _Tb(self):
        """Calculate boiling temperature"""
        # TODO: Add more general correlation or specific by compound type
        # For now suppose inorganic compound
        Tb = unidades.Temperature(1.64*self.Tc)
        return Tb

    def _SG(self):
        # FIXME: Don't work
        # volumen = self.Vliq*(5.7+3*288.71/self.Tc)
        # return 1/volumen*18
        return 1

    def _cp(self):
        """Default method to calculate the temperature dependence of ideal
        specific heat for method don't define this property"""
        cpa = unidades.Enthalpy((0.036863384*self.Kw-0.4673722)*self.M,
                                "Btulb").kcalkg/self.M*1.8
        cpb = unidades.Enthalpy((3.1865e-5*self.Kw+0.001045186)*self.M,
                                "Btulb").kcalkg/self.M*1.8**2
        cpc = unidades.Enthalpy(-4.9572e-7*self.M,
                                "Btulb").kcalkg/self.M*1.8**3
        cp = [cpa, cpb, cpc, 0, 0, 0]
        return cp

    def _Vc(self):
        """Método de cálculo del volumen crítico"""
        if self.Tc.R < 536.67:
            D = 8.75+1.987*log(self.Tb.R)+self.Tb.R/1.8
        # elif self.Tc.R > 593:
            # # FIXME: Don't work exactly
            # D = (0.398907*self.SG*(self.f_acent-592.4439)/self.M)**0.5
        else:
            f = ((self.Tc.R-536.67)/(self.Tc.R-self.Tb.R))**0.38
            D = 8.75+1.987*log(self.Tb.R)+self.Tb.R/1.8*f

        Zc = 1/(3.43+6.7e-9*D**2)
        return Zc*self.Tc.R*10.73/self.Pc.psi

    def _Calor_vaporizacion(self):
        """Método de cálculo del calor de vaporización,
        ref. chemcad pag 60"""
        tbr = self.Tb/self.Tc
        return unidades.Enthalpy(1.093*R*1000*self.Tc*(tbr*(log(self.Pc.atm)-1)/(0.930-tbr))/self.M, "calg")

    def _Rackett(self):
        """ref 64"""
        return 0.29056-0.08775*self.f_acent

    def _VLiq(self):
        Tr = 298.15/self.Tc
        V = R_atml*1000*self.Tc/self.Pc.atm*self.rackett**(1+(1-Tr)**(2/7))
        return V/(5.7+1611/self.Tc)  # cm3/mol

    def _SolubilityParameter(self):
        r"""Calculation procedure for solubility parameter when the compound
        has no defined value for that using the definition equation:

        .. math::
            \delta=\sqrt{\Delta H_{v}-\frac{RT}{M}}\rho

        Referenced in API databook Table 8B1.7 comment pag 812
        """
        if self.Tb < 298.15:
            T = self.Tb
        else:
            T = 298.15

        rhoL = RhoL_Rackett(T, self.Tc, self.Pc, self.rackett, self.M)
        DHv = self.Hv.calg
        delta = DHv-R_cal*T/self.M*rhoL
        if delta < 0:
            delta = 0

        return unidades.SolubilityParameter(delta**0.5, "calcc")

    @staticmethod
    def _atomicComposition(group):
        """Calculate of atomic composition"""
        total = {}
        for g in group:
            for key, value in g.items():
                if key not in total:
                    total[key] = 0
                total[key] += value
        return total

    def EmpiricFormula(self):
        "Calculate the empiric formulae of compound from group contribution"""

        total = self._atomicComposition(self.group)
        string = ""
        formula = ""
        for element in ["C", "H", "N", "O", "S", "F",  "Cl", "Br", "I"]:
            if element not in total:
                continue
            if total[element] > 1:
                string += "%s<sub>%i</sub>" % (element, total[element])
                formula += "%s%i" % (element, total[element])
            elif total[element] == 1:
                string += "%s" % element
                formula += "%s" % element
            del total[element]

        # Add other element at end of formule
        for key, value in total.items():
            if value > 1:
                string += "%s<sub>%i</sub>" % (key, value)
                formula += "%s%i" % (key, value)
            elif value == 1:
                string += "%s" % key
                formula += "%s" % key

        return string, formula


@refDoc(__doi__, [1])
def cpLS_Hurst(group):
    """Calculate liquid and solid heat capacities using the Hurst-Harrison
    method

    Parameters
    ----------
    group : dict
        Atomic composition of component

    Returns
    -------
    cpL : float
        Liquid heat capacity, []
    cpS : float
        Soid heat capacity, []

    Examples
    --------
    Example in [20]_, GdF3

    >>> "%0.1f" % cpLS_Hurst(group=[{"Gd": 1, "F": 3}])[0]
    '105.1'
    """
    Solid = {"H": 7.56, "Li": 23.25, "Be": 12.47, "B": 10.10, "C": 10.89,
             "N": 18.74, "O": 13.42, "F": 26.16, "Na": 26.19, "Mg": 22.69,
             "Al": 18.07, "Si": 7.00, "S": 12.36, "CI": 24.69, "K": 28.78,
             "Ca": 28.25, "Ti": 27.24, "V": 29.36, "Mn": 28.06, "Fe": 29.08,
             "Co": 25.71, "Ni": 25.46, "Cu": 26.92, "Br": 25.36, "Sr": 28.41,
             "Zr": 26.82, "Mo": 29.44, "I": 25.29, "Ba": 32.37, "W": 30.87,
             "Hg": 27.87, "Pb": 31.60, "Misc": 26.63}
    Liquid = {"H": 9.20, "C": 13.08, "N": 30.19, "O": 16.00, "F": 19.47,
              "S": 32.05, "Cl": 31.91, "Br": 37.23, "Misc": 26.19}

    cmp = GroupContribution._atomicComposition(group)

    cpl, cps = 0, 0
    for key, c in cmp.items():
        cps += c*Solid.get(key, Solid["Misc"])
        cpl += c*Liquid.get(key, Liquid["Misc"])

    return cps, cpl
