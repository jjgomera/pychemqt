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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Module to implement new component by group contribution methods
    * :func:`newComponente`: Base class for new component definition
    * :func:`GroupContribution`: Common group contribution functionality
    * :func:`Joback`: Group contribution method of Joback
    * :func:`Constantinou`: Group contribution method of Constaninou and Gani
    * :func:`Wilson`: Group contribution method of Wilson and Jasperson
    * :func:`Marrero`: Group contribution method of Marrero and Pardillo
    * :func:`Elliott`: Group contribution method of Elliott
    * :func:`Ambrose`: Group contribution method of Ambrose
    * :func:`Klincewicz`: Group contribution method of Klincewicz
    * :func:`Lydersen`: Group contribution method of Lydersen
    * :func:`Valderrama`: Group contribution method of Valderrama
    * :func:`Nannoolal`: Group contribution method of Nannoolal
    * :func:`Wen`: Group contribution method of Wen-Qiang
    * :func:`Li`: Group contribution method of Li-Xia-Xiang

Other methods:
    * :func:`cpLS_Hurst`: Liquid-Solid heat capacities using the Hurst-Harrison
    method
'''


from itertools import permutations
from math import exp, log
import time

from scipy.constants import R
from PyQt5.QtWidgets import QApplication

from lib import unidades
from lib.compuestos import atomic_decomposition, facent_LeeKesler, RhoL_Rackett
from lib.physics import R_atml, R_cal
from lib.elemental import databank
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
    2:
        {"autor": "Joback, K.G., Reid, R.C.",
         "title": "Estimation of Pure-Component Properties from "
                  "Group-Contributions.",
         "ref": "Chemical Engineering Communications, 57 (1987) 233-243",
         "doi": "10.1080/00986448708960487"},
    3:
        {"autor": "Constantinou, L., Gani, R.",
         "title": "New Group Controbution Method for Estimating Properties of "
                  "Pure Compounds",
         "ref": "AIChE J. 40(10) (1994) 1697-1710",
         "doi": "10.1002/aic.690401011"},
    4:
        {"autor": "Constantinou, L., Gani, R., O’Connell, J.P.",
         "title": "Estimation of the Acentric Factor and the Liquid Molar "
                  "Volume at 298K Using a New Group Contribution Method",
         "ref": "Fluid Phase Equilibria 103 (1995) 11-22",
         "doi": "10.1016/0378-3812(94)02593-p"},
    5:
        {"autor": "Wilson, G.M. Jasperson, L.V.",
         "title": "Critical constants Tc and Pc, estimation based on zero, "
                  "first and second order methods",
         "ref": "Paper given at AIChE Spring National Meeting, New Orleans, "
                "LA, USA, February 25-29, 1996.",
         "doi": ""},
    6:
        {"autor": "Marrero-Morejón, J., Pardillo-Fontdevila, F.",
         "title": "Estimation of Pure Compound Properties Using "
                  "Group-Interaction C9ontributions",
         "ref": "AIChE J., 45(3) (1999) 615-621",
         "doi": "10.1002/aic.690450318"},
    7:
        {"autor": "Marrero-Morejon, J., Pardillo-Fontdevila, E.",
         "title": "Estimation of Liquid Viscosity at Ambient Temperature of "
                  "Pure Organic Compounds by Using Group-Interaction "
                  "Contributions",
         "ref": "Chemical Engineering Journal 79 (2000) 69-72",
         "doi": "10.1016/s1385-8947(99)00173-4"},
    8:
        {"autor": "Ambrose, D.",
         "title": "Correlation and Estimation of Vapor-Liquid Critical "
                  "Properties: I. Critical Temperatures of Organic Compounds",
         "ref": "National Physical Laboratory, Teddington, NPL Rep. Chern.  "
                "92, 1978, corrected 1980.",
         "doi": ""},
    9:
        {"autor": "Ambrose, D.",
         "title": "Correlation and Estimation of Vapor-Liquid Critical "
                  "Properties: II. Critical Pressures and Volumes of Organic "
                  "Compounds",
         "ref": "National Physical Laboratory, Teddington, NPL Rep. 98, 1979",
         "doi": ""},
    10:
        {"autor": "API",
         "title": "Technical Data book: Petroleum Refining 6th Edition",
         "ref": "",
         "doi": ""},
    11:
        {"autor": "Maloney, J.O.",
         "title": "Perry's Chemical Engineers' Handbook 8th Edition",
         "ref": "McGraw Hill (2008)",
         "doi": ""},
    12:
        {"autor": "Klincewicz, K.M., Reid, R.C.",
         "title": "Estimation of Critical Properties with Group Contribution "
                  "Methods",
         "ref": "AIChE J. 30(1) (1984) 137-142",
         "doi": "10.1002/aic.690300119"},
    13:
        {"autor": "Lydersen, A. L.",
         "title": "Estimation of Critical Properties of Organic Compounds",
         "ref": "Coll. Eng. Univ. Wisconsin, Engineering Experimental Station "
                "Rept. 3, Madison, WI (1955)",
         "doi": ""},
    14:
        {"autor": "Valderrama, J.O., Álvarez, V.H.",
         "title": "A New Group Contribution Method Based on Equation of State "
                  "Parameters to Evaluate the Critical Properties of Simple "
                  "and Complex Molecules",
         "ref": "Can. J. Chem. Eng. 84(4) (2006) 431-446",
         "doi": "10.1002/cjce.5450840404"},
    15:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D., Cordes, W.",
         "title": "Estimation of Pure Component Properties 1. Estimation of "
                  "the Normal Boiling Point of Non-electrolyte Organic "
                  "Compounds Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 226 (2004) 45-63",
         "doi": "10.1016/j.fluid.2004.09.001"},
    16:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 2. Estimation of "
                  "Critical Property Data by Group Contribution",
         "ref": "Fluid Phase Equilib., 252 (2007) 1-27",
         "doi": "10.1016/j.fluid.2006.11.014"},
    17:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 3. Estimation of "
                  "the Vapor Pressure of Non-Electrolyte Organic Compounds "
                  "Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 269 (2008) 117-133",
         "doi": "10.1016/j.fluid.2008.04.020"},
    18:
        {"autor": "Nannoolal, Y., Rarey, J., Ramjugernath, D.",
         "title": "Estimation of Pure Component Properties 4. Estimation of "
                  "the Saturted Liquid Viscosity of Non-Electrolyte Organic "
                  "Compounds Via Group Contributions and Group Interactions",
         "ref": "Fluid Phase Equilib., 281 (2009) 97-119",
         "doi": "10.1016/j.fluid.2009.02.016"},
    19:
        {"autor": "Wen, X., Quiang, Y.",
         "title": "A New Group Contribution Method for Estimating Critical "
                  "Properties of Orgnic Compounds",
         "ref": "Ind. Eng. Chem. Res. 40(26) (2001) 6245-6250.",
         "doi": "10.1021/ie010374g"},
    20:
        {"autor": "Hurst, J.E., Harrison, B.K.",
         "title": "Estimation of Liquid and Solid Heat Capacities Using a "
                  "Modified Kopp's Rule",
         "ref": "Chem. Eng. Comm. 112 (1992) 21-30",
         "doi": "10.1080/00986449208935989"},
    21:
        {"autor": "Li, J., Xia, L., Xiang, S.",
         "title": "A New Method Based on Elements and Chemical Bonds for "
                  "Organic Compounds Critical Properties Estimation",
         "ref": "Fluid Phase Equil. 417 (2016) 1-6",
         "doi": "10.1016/j.fluid.2016.01.008"},
    22:
        {"autor": "Marrero, J.; Gani, R.",
         "title": "Group-contribution based estimation of pure component "
                  "properties",
         "ref": "Fluid Phase Equilib. 183-184 (2001), 183-208.",
         "doi": "10.1016_s0378-3812(01)00431-9"}}


# Get molecular weight of atomic elements
MW = {}
databank.execute("SELECT symbol, atomic_mass FROM ELEMENTS")
for symbol, m in databank:
    MW[symbol] = m


class newComponente(object):
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
        ele.append(self.watson)

        ele.append([])

        ele.append(0)
        ele.append(self.Tb)
        ele.append(self.Tf)
        ele.append("")
        ele.append("")

        ele.append([])

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
        ele.append(0)
        ele.append(0)
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
            self.msg = QApplication.translate("pychemqt", "undefined group")
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
        """Calculation procedure for solubility parameter when the compound has
        no defined value for that using the definition equation:

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


@refDoc(__doi__, [1, 2, 11])
class Joback(GroupContribution):
    """
    Group contribution for definition of unknown component using the Joback
    procedure (1987)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float, optional
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    M, Tb and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    p-dichlorobenzene example in 2_, Table V

    >>> cmp = Joback(group=[16, 13, 14], contribution=[2, 4, 2])
    >>> "%0.1f %0.0f %0.0f %0.1f" % (cmp.Tb, cmp.Tf, cmp.Tc, cmp.Pc.bar)
    '443.4 256 675 41.5'
    >>> "%0.0f" % (cmp.Vc.ccg*cmp.M)
    '362'
    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '26.41 78.56'
    >>> "%0.0f %0.0f" % (cmp._Cp0(298).JgK*cmp.M, cmp._Cp0(400).JgK*cmp.M)
    '112 139'
    >>> "%0.0f %0.0f" % (cmp._Cp0(800).JgK*cmp.M, cmp._Cp0(1000).JgK*cmp.M)
    '206 224'
    >>> "%0.2f %0.1f" % (cmp.Hv.kJg*cmp.M, cmp.Hm.kJg*cmp.M)
    '40.66 13.3'
    >>> "%0.2e %0.2e" % (cmp._Visco(333.8), cmp._Visco(374.4))
    '7.26e-04 4.92e-04'
    >>> "%0.2e %0.1e" % (cmp._Visco(403.1), cmp._Visco(423.3))
    '3.91e-04 3.4e-04'

    Example 2-1 in [1]_, 2-ethylphenol critical properties

    >>> cmp = Joback(group=[0, 1, 13, 14, 20], contribution=[1, 1, 4, 2, 1])
    >>> "%0.2f %0.1f %0.2f" % (cmp.Tb, cmp.Tc, cmp.Pc.bar)
    '489.94 716.0 44.09'
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '341.5'
    >>> cmp.formula
    'C8H10O'

    Example 3-1 in [1]_, 2,4 dimethylphenol ΔH and ΔG

    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '-149.23 -25.73'
    >>> "%0.1f" % (cmp._Cp0(700).JgK*cmp.M)
    '281.2'

    Example 2-10 in [1]_, 2,4 dimethylphenol Tb and Tf

    >>> cmp = Joback(group=[0, 13, 14, 20], contribution=[2, 3, 3, 1])
    >>> "%0.2f %0.2f" % (cmp.Tf, cmp.Tb)
    '330.58 494.92'

    Example in http://en.wikipedia.org/wiki/Joback_method, acetone

    >>> cmp = Joback(group=[0, 23], contribution=[2, 1])
    >>> "%0.3f %0.3f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Tb, cmp.Tf)
    '500.559 48.025 322.11 173.5'
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '209.5'
    >>> "%0.2f %0.2f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '-217.83 -154.54'
    >>> "%0.4f" % (cmp._Cp0(300).JgK*cmp.M)
    '75.3264'
    >>> "%0.2f %0.2f" % (cmp.Hm.kJg*cmp.M, cmp.Hv.kJg*cmp.M)
    '5.12 29.02'
    >>> "%0.7f" % (cmp._Visco(300))
    '0.0002942'

    Example in [11]_ pag. 2-470, o-xylene

    >>> cmp = Joback(group=[13, 14, 0], contribution=[4, 2, 2], Tb=417.58)
    >>> "%0.2f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '630.37 35.86 375.5'

    Example in [11]_ pag. 2-470, sec-butanol

    >>> cmp = Joback(group=[0, 1, 2, 19], contribution=[2, 1, 1, 1], Tb=372.7)
    >>> "%0.1f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '534.1 44.33 272.5'
    """
    __title__ = "Joback-Reid (1987)"
    _coeff = {
        # Table III
        "tc": [0.0141, 0.0189, 0.0164, 0.0067, 0.0113, 0.0129, 0.0117, 0.0026,
               0.0027, 0.002, 0.01, 0.0122, 0.0042, 0.0082, 0.0143, 0.0111,
               0.0105, 0.0133, 0.0068, 0.0741, 0.024, 0.0168, 0.0098, 0.038,
               0.0284, 0.0379, 0.0791, 0.0481, 0.0143, 0.0243, 0.0295, 0.0130,
               .0169, .0255, .0085, .0, .0496, .0437, 0.0031, 0.0119, 0.0019],
        "Pc": [-0.0012, 0.0, 0.002, 0.0043, -0.0028, -0.0006, 0.0011, 0.0028,
               -0.0008, 0.0016, 0.0025, 0.0004, 0.0061, 0.0011, 0.0008, -.0057,
               -0.0049, 0.0057, -0.0034, 0.0112, 0.0184, 0.0015, 0.0048, .0031,
               0.0028, 0.0030, 0.0077, 0.0005, 0.0101, 0.0109, 0.0077, 0.0114,
               0.0074, -0.0099, .0076, .0, -.0101, .0064, .0084, .0049, .0051],
        "vc": [65, 56, 41, 27, 56, 46, 38, 36, 46, 37, 48, 38, 27, 41, 32, 27,
               58, 71, 97, 28, -25, 18, 13, 62, 55, 82, 89, 82, 36, 38, 35, 29,
               9, 0, 34, 0, 91, 91, 63, 54, 38],
        "tb": [23.58, 22.88, 21.74, 18.25, 18.18, 24.96, 24.14, 26.15, 9.20,
               27.38, 27.15, 21.78, 21.32, 26.73, 31.01, -0.03, 38.13, 66.86,
               93.84, 92.88, 76.34, 22.42, 31.22, 76.75, 94.97, 72.24, 169.09,
               81.10, -10.5, 73.23, 50.17, 52.82, 11.74, 74.6, 57.55, 0.0,
               125.66, 152.54, 63.56, 68.78, 52.10],
        "tf": [-5.10, 11.27, 12.64, 46.43, -4.32, 8.73, 11.14, 17.78, -11.18,
               64.32, 7.75, 19.88, 60.15, 8.13, 37.02, -15.78, 13.55, 43.43,
               41.69, 44.45, 82.83, 22.23, 23.05, 61.20, 75.97, 36.9, 155.5,
               53.6, 2.08, 66.89, 52.66, 101.51, 48.84, 0, 68.4, 0.0, 59.89,
               127.24, 20.09, 34.4, 79.93],
        "hf": [-76.45, -20.64, 29.89, 82.23, -9.63, 37.97, 83.99, 142.14,
               79.30, 115.51, -26.8, 8.67, 79.72, 2.09, 46.43, -251.92, -71.55,
               -29.48, 21.06, -208.04, -221.65, -132.22, -138.16, -133.22,
               -164.50, -162.03, -426.72, -337.92, -247.61, -22.02, 53.47,
               31.65, 123.34, 23.61, 55.52, 93.7, 88.43, -66.57, -17.33, 41.87,
               39.1],
        "gf": [-43.96, 8.42, 58.36, 116.02, 3.77, 48.53, 92.36, 136.7, 77.71,
               109.82, -3.68, 40.99, 87.88, 11.30, 54.05, -247.19, -64.31,
               -38.06, 5.74, -189.2, -197.37, -105.0, -98.22, -120.50, -126.27,
               -143.48, -387.87, -301.95, -250.83, 14.07, 89.39, 75.61, 163.16,
               0.0, 79.93, 119.66, 89.22, -16.83, -22.99, 33.12, 27.73],
        "hv": [567, 532, 404, 152, 412, 527, 511, 636, 276, 789, 573, 464, 154,
               608, 731, -160, 1083, 1573, 2275, 4021, 2987, 576, 1119, 2144,
               1588, 2173, 4669, 2302, 1412, 2578, 1538, 1656, 453, 797, 1560,
               2908, 3071, 4000, 1645, 1629, 1430],
        "hm": [217, 619, 179, -349, -113, 643, 732, 1128, 555, 992, 117, 775,
               -328, 263, 572, 334, 601, 861, 651, 575, 1073, 284, 1405, 1001,
               0, 764, 2641, 1663, 866, 840, 1197, 1790, 1124, 0, 872, 0, 577,
               2313, 564, 987, 372],
        "cpa": [19.5, -0.909, -23.0, -66.2, -23.6, -8.0, -28.1, 27.4, 24.5,
                7.87, -6.03, 8.67, -90.9, -2.14, -8.25, 26.5, 33.3, 28.6, 32.1,
                25.7, -2.81, 25.5, 12.2, 6.45, 30.4, 30.9, 24.1, 24.5, 6.82,
                26.9, -1.21, 11.8, -31.1, 0.0, 8.83, 5.69, 36.5, 25.9, 35.3,
                19.6, 16.7],
        "cpb": [-8.08e-3, 9.5e-2, 2.04e-1, 4.27e-1, -3.81e-2, 1.05e-1, 2.08e-1,
                -5.57e-2, -2.71e-2, 2.01e-2, 8.54e-2, 1.62e-1, 5.57e-1,
                5.74e-2, 1.01e-1, -9.13e-2, -9.63e-2, -6.49e-2, -6.41e-2,
                -6.91e-2, 1.11e-1, -6.32e-2, -1.26e-2, 6.7e-2, -8.29e-2,
                -3.36e-2, 4.27e-2, 4.02e-2, 1.96e-2, -4.12e-2, 7.62e-2,
                -2.3e-2, 2.27e-1, 0.0, -3.84e-3, -4.12e-3, -7.33e-2, -3.74e-3,
                -7.58e-2, -5.61e-3, 4.81e-3],
        "cpc": [1.53e-4, -5.44e-5, -2.65e-4, -6.41e-4, 1.72e-4, -9.63e-5,
                -3.06e-4, 1.01e-4, 1.11e-4, -8.33e-6, -8.0e-6, -1.6e-4,
                -9.0e-4, -1.64e-6, -1.42e-4, 1.91e-4, 1.87e-4, 1.36e-4,
                1.26e-4, 1.77e-4, -1.16e-4, 1.11e-4, 6.03e-5, -3.57e-5,
                2.36e-4, 1.6e-4, 8.04e-5, 4.02e-5, 1.27e-5, 1.64e-4, -4.86e-5,
                1.07e-4, -3.2e-4, 0.0, 4.35e-5, 1.28e-4, 1.84e-4, 1.29e-4,
                1.85e-4, 4.02e-5, 2.77e-5],
        "cpd": [-9.67e-8, 1.19e-8, 1.2e-7, 3.01e-7, -1.03e-7, 3.56e-8, 1.46e-7,
                -5.02e-8, -6.78e-8, 1.39e-9, -1.8e-8, 6.24e-8, 4.69e-7,
                -1.59e-8, 6.78e-8, -1.03e-7, -9.96e-8, -7.45e-8, -6.87e-8,
                -9.88e-8, 4.94e-8, -5.48e-8, -3.86e-8, 2.86e-9, -1.31e-7,
                -9.88e-8, -6.87e-8, -4.52e-8, -1.78e-8, -9.76e-8, 1.05e-8,
                -6.28e-8, 1.46e-7, 0.0, -2.6e-8, -8.88e-8, -1.03e-7, -8.88e-8,
                -1.03e-7, -2.76e-8, -2.11e-8],
        "mua": [548.29, 94.16, -322.15, -573.56, 495.01, 82.28, 0, 0, 0, 0,
                307.53, -394.29, 0, 259.65, -245.74, 0, 625.45, 738.91, 809.55,
                2173.72, 3018.17, 122.09, 440.24, 340.35, 0, 740.92, 1317.23,
                483.88, 675.24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "mub": [-1.719, -0.199, 1.187, 2.307, -1.539, -0.242, 0, 0, 0, 0,
                -0.798, 1.251, 0, -0.702, 0.912, 0, -1.814, -2.038, -2.224,
                -5.057, -7.314, -0.386, -0.953, -0.35, 0, -1.713, -2.578,
                -0.966, -1.34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],

        # Name and element composition
        "txt": [("CH3", ),
                ("CH2", ),
                ("CH", ),
                ("C", ),
                ("=CH2", ),
                ("=CH", ),
                ("=C", ),
                ("=C=", ),
                ("≡CH", ),
                ("≡C", ),
                ("CH2 (cyclic)", ),
                ("CH (cyclic)", ),
                ("C (cyclic)", ),
                ("-CH (Aromatic)", ),
                ("=C (Aromatic)", ),
                ("F", ),
                ("Cl", ),
                ("Br", ),
                ("I", ),
                ("-OH", ),
                ("-OH (Aromatic)", ),
                ("-O-", ),
                ("-O- (cyclic)", ),
                ("C=O", ),
                ("C=O (cyclic)", ),
                ("CH=O", ),
                ("COOH", ),
                ("COO", ),
                ("=O", ),
                ("NH2", ),
                ("NH", ),
                ("NH (cyclic)", ),
                ("N", ),
                ("=N-", ),
                ("=N- (cyclic)", ),
                ("=NH", ),
                ("CN", ),
                ("NO2", ),
                ("SH", ),
                ("S", ),
                ("S (cyclic)", )]}

    FirstOrder = 41

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        if self.kwargs["Tb"]:
            Tb = self.kwargs["Tb"]
        else:
            # Eq 2
            Tb = 198.2+sum([c*self._coeff["tb"][i] for i, c in zip(
                self.kwargs["group"], self.kwargs["contribution"])])
        self.Tb = unidades.Temperature(Tb)

        # Equations of Table II
        self.Na = self._atoms()
        tcsuma, pcsuma, vcsuma = 0, 0, 0
        Tf = 122.5
        Hf = 68.29
        Gf = 53.88
        Hv = 15.3
        Hm = -0.88
        cpa = -37.93
        cpb = 0.21
        cpc = -3.91e-4
        cpd = 2.06e-7
        mua, mub = 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            Tf += c*self._coeff["tf"][i]                                # Eq 3
            tcsuma += c*self._coeff["tc"][i]                            # Eq 4
            pcsuma += c*self._coeff["Pc"][i]                            # Eq 5
            vcsuma += c*self._coeff["vc"][i]                            # Eq 6
            Hf += c*self._coeff["hf"][i]                                # Eq 7
            Gf += c*self._coeff["gf"][i]                                # Eq 8
            Hv += c*self._coeff["hv"][i]*0.004184                       # Eq 10
            Hm += c*self._coeff["hm"][i]*0.004184                       # Eq 11
            cpa += c*self._coeff["cpa"][i]
            cpb += c*self._coeff["cpb"][i]
            cpc += c*self._coeff["cpc"][i]
            cpd += c*self._coeff["cpd"][i]
            mua += c*self._coeff["mua"][i]
            mub += c*self._coeff["mub"][i]
        self.Tf = unidades.Temperature(Tf)
        self.Tc = unidades.Temperature(self.Tb/(0.584+0.965*tcsuma-tcsuma**2))
        self.Pc = unidades.Pressure((0.113+0.0032*self.Na-pcsuma)**-2, "bar")
        self.Vc = unidades.SpecificVolume((vcsuma+17.5)/1000/self.M)
        self.Hf = unidades.Enthalpy(Hf/self.M, "kJg")
        self.Gf = unidades.Enthalpy(Gf/self.M, "kJg")
        self.Hv = unidades.Enthalpy(Hv/self.M, "kJg")
        self.Hm = unidades.Enthalpy(Hm/self.M, "kJg")
        self.cp = [cpa, cpb, cpc, cpd, 0, 0]

        # Adjust the viscosity correlation with the parametric viscosity
        self.mul = [mua-597.82, mub-11.202]

        GroupContribution.calculo(self)

    def _Visco(self, T):
        """Viscosity calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        mu : float
            Viscosity, [Pas]
        """
        mu = self.M*exp(self.mul[0]/T + self.mul[1])
        return unidades.Viscosity(mu)

    def _Cp0(self, T):
        """Ideal gas specific heat calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        cp0 : float
            Ideal gas specific heat, [J/kgK]
        """
        cp0 = 0
        for i, c in enumerate(self.cp):
            cp0 += c*T**i

        return unidades.SpecificHeat(cp0/self.M, "JgK")


@refDoc(__doi__, [1, 3, 4, 11])
class Constantinou(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Constantinou-Gani procedure (1994)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example from 3_ to distinguish between isomers, dimethylhexanes
    1st order only

    >>> cmp = Constantinou(group=[0, 1, 2], contribution=[4, 2, 2])
    >>> "%0.2f %0.2f" % (cmp.Tc, cmp.Tb)
    '557.91 385.92'

    Adding 2nd order contributions

    >>> c23 = Constantinou(group=[0, 1, 2, 80], contribution=[4, 2, 2, 1])
    >>> c24 = Constantinou(group=[0, 1, 2, 78], contribution=[4, 2, 2, 1])
    >>> c25 = Constantinou(group=[0, 1, 2, 78], contribution=[4, 2, 2, 2])
    >>> "%0.2f %0.2f %0.2f" % (c23.Tc, c24.Tc, c25.Tc)
    '566.60 553.41 548.80'
    >>> "%0.2f %0.2f %0.2f" % (c23.Tb, c24.Tb, c25.Tb)
    '391.41 382.32 378.64'

    Examples from Table 15a and 15b and 16 in 3_

    >>> cmp1 = Constantinou(
    ... group=[0, 1, 2, 3, 15], contribution=[3, 1, 1, 1, 2])
    >>> cmp2 = Constantinou(group=[0, 1, 2, 38], contribution=[2, 1, 1, 1])
    >>> cmp3 = Constantinou(
    ... group=[0, 1, 2, 3, 15, 105, 106], contribution=[3, 1, 1, 1, 2, 1, 1])
    >>> cmp4 = Constantinou(
    ... group=[0, 1, 2, 38, 78], contribution=[2, 1, 1, 1, 1])
    >>> "%0.2f %0.2f %0.2f %0.2f" % (cmp1.Tb, cmp2.Tb, cmp3.Tb, cmp4.Tb)
    '488.39 452.14 465.18 449.54'

    >>> cmp1 = Constantinou(group=[6, 4, 0], contribution=[1, 1, 1])
    >>> cmp2 = Constantinou(group=[0, 13, 12, 10], contribution=[1, 1, 1, 4])
    >>> cmp3 = Constantinou(group=[6, 4, 0, 88], contribution=[1, 1, 1, 1])
    >>> "%0.2f %0.3f %0.3f" % (
    ... cmp1.Gf.kJg*cmp1.M, cmp2.Gf.kJg*cmp2.M, cmp3.Gf.kJg*cmp3.M)
    '150.47 131.007 144.965'

    Examples from Table A1, A2 and A3 in [4]_

    >>> cmp1 = Constantinou(group=[0, 2], contribution=[5, 3])
    >>> cmp2 = Constantinou(group=[0, 17, 1, 2], contribution=[2, 1, 1, 1])
    >>> cmp3 = Constantinou(group=[0, 2, 80], contribution=[5, 3, 2])
    >>> cmp4 = Constantinou(
    ... group=[0, 17, 1, 2, 78, 95], contribution=[2, 1, 1, 1, 1, 1])
    >>> "%0.5f %0.5f %0.5f %0.5f" % (
    ... cmp1.Vliq/cmp1.M, cmp2.Vliq/cmp2.M, cmp3.Vliq/cmp3.M, cmp4.Vliq/cmp4.M)
    '0.16414 0.12446 0.16008 0.12549'
    >>> "%0.4f %0.4f %0.4f %0.4f" % (
    ... cmp1.f_acent, cmp2.f_acent, cmp3.f_acent, cmp4.f_acent)
    '0.3195 0.4430 0.3167 0.3863'

    Example 2-2 in [1]_, 2-ethylphenol critical properties

    >>> cmp = Constantinou(group=[0, 10, 13, 16], contribution=[1, 4, 1, 1])
    >>> "%0.1f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '718.6 42.97 371.9'
    >>> cmp.formula
    'C8H10O'

    Example 2-3 in [1]_, butanols critical properties

    >>> b1 = Constantinou(group=[0, 1, 15], contribution=[1, 3, 1])

    2-met-1-propanol CHOH Second order contribution is zero

    >>> b2m1 = Constantinou(
    ... group=[0, 1, 2, 15, 78], contribution=[2, 1, 1, 1, 1])
    >>> b2m2 = Constantinou(
    ... group=[0, 3, 15, 79, 106], contribution=[3, 1, 1, 1, 1])
    >>> b2 = Constantinou(
    ... group=[0, 1, 2, 15, 105], contribution=[2, 1, 1, 1, 1])
    >>> "%0.2f %0.2f %0.2f %0.2f" % (b1.Tc, b2m1.Tc, b2m2.Tc, b2.Tc)
    '558.91 543.31 497.46 521.57'
    >>> "%0.2f %0.2f %0.2f %0.2f" % (
    ... b1.Pc.bar, b2m1.Pc.bar, b2m2.Pc.bar, b2.Pc.bar)
    '41.97 41.66 42.32 44.28'
    >>> "%0.1f %0.1f %0.1f %0.1f" % (b1.Vc.ccg*b1.M, b2m1.Vc.ccg*b2m1.M,
    ...  b2m2.Vc.ccg*b2m2.M, b2.Vc.ccg*b2.M)
    '276.9 276.0 280.2 264.2'

    Example 2-9 in [1]_, 2,3,3-trimethylpentane

    >>> c1 = Constantinou(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1])
    >>> c2 = Constantinou(group=[0, 1, 2, 3, 81], contribution=[5, 1, 1, 1, 1])
    >>> "%0.3f %0.3f" % (c1.f_acent, c2.f_acent)
    '0.301 0.292'

    Example 2-11 in [1]_, 2,4-dimethylphenol

    >>> cmp = Constantinou(group=[10, 12, 16], contribution=[3, 2, 1])
    >>> "%0.2f %0.2f" % (cmp.Tf, cmp.Tb)
    '315.96 492.33'

    Example 2-12 in [1]_, cycloalkanes

    >>> c7 = Constantinou(group=[1, 87], contribution=[7, 1])
    >>> mc6 = Constantinou(
    ... group=[0, 1, 2, 86, 92], contribution=[1, 5, 1, 1, 1])
    >>> ec5 = Constantinou(
    ... group=[0, 1, 2, 85, 92], contribution=[1, 5, 1, 1, 1])
    >>> c5 = Constantinou(
    ... group=[0, 1, 2, 85], contribution=[2, 3, 2, 1])
    >>> t5 = Constantinou(
    ... group=[0, 1, 2, 85], contribution=[2, 3, 2, 1])
    >>> "%0.2f %0.2f %0.2f %0.2f %0.2f" % (c7.Tf, mc6.Tf, ec5.Tf, c5.Tf, t5.Tf)
    '266.15 146.46 122.14 166.79 166.79'
    >>> "%0.2f %0.2f %0.2f %0.2f %0.2f" % (c7.Tb, mc6.Tb, ec5.Tb, c5.Tb, t5.Tb)
    '391.93 377.81 377.69 364.27 364.27'

    Example 3-2 in [1]_, 2-ethylphenol ΔH and ΔG

    >>> cmp = Constantinou(group=[0, 10, 13, 16], contribution=[1, 4, 1, 1])
    >>> "%0.3f %0.3f" % (cmp.Hf.kJg*cmp.M, cmp.Gf.kJg*cmp.M)
    '-145.561 -23.595'
    >>> "%0.2f" % (cmp._Cp0(700).JgK*cmp.M)
    '286.35'

    Example 3-3 in [1]_, butanols formations properties
    Component definition above in example 2-3

    >>> "%0.2f %0.2f %0.2f %0.2f" % (
    ... b1.Hf.kJg*b1.M, b2m1.Hf.kJg*b2m1.M, b2m2.Hf.kJg*b2m2.M, b2.Hf.kJg*b2.M)
    '-278.82 -287.87 -316.77 -290.90'
    >>> "%0.2f %0.2f %0.2f %0.2f" % (
    ... b1.Gf.kJg*b1.M, b2m1.Gf.kJg*b2m1.M, b2m2.Gf.kJg*b2m2.M, b2.Gf.kJg*b2.M)
    '-156.75 -161.10 -180.70 -168.17'
    >>> "%0.1f %0.1f %0.1f %0.1f" % (b1._Cp0(298).JgK*b1.M,
    ...  b2m1._Cp0(298).JgK*b2m1.M, b2m2._Cp0(298).JgK*b2m2.M,
    ... b2._Cp0(298).JgK*b2.M)
    '110.5 109.8 111.9 111.7'

    Example in [11]_ pag. 2-471, 2,6-dimethylpyridine

    >>> cmp = Constantinou(group=[0, 36, 86], contribution=[2, 1, 1])
    >>> "%0.0f" % cmp.Tf
    '278'
    """
    __title__ = "Constantinou-Gani (1995)"
    _coeff = {
        # The Second order term are append to the first order in each table

        # Table 1 in [3]_
        "tc": [1.6781, 3.4920, 4.0330, 4.8823, 5.0146, 7.3691, 6.5081, 8.9582,
               11.3764, 9.9318, 3.7337, 14.6409, 8.2130, 10.3239, 10.4664,
               9.7292, 25.9145, 13.2896, 14.6273, 10.1986, 12.5965, 13.8116,
               11.6057, 6.4737, 6.0723, 5.0663, 9.5059, 12.1726, 10.2075,
               9.8544, 10.4677, 7.2121, 7.6924, 5.5172, 28.7570, 29.1528,
               27.9464, 20.3781, 23.7593, 11.0752, 10.8632, 11.3959, 16.3945,
               0, 18.5875, 14.1565, 24.7369, 23.2050, 34.5870, 13.8058,
               17.3947, 10.5371, 7.5433, 11.4501, 5.4334, 2.8977, 0, 2.4778,
               1.7399, 3.5192, 12.1084, 9.8408, 0, 4.8923, 1.5974, 65.1053, 0,
               0, 36.1403, 0, 0, 17.9668, 0, 14.3969, 17.7916, 0, 0, 0,
               -0.5334, -0.5143, 1.0699, 1.9886, 5.8254, -2.3305, -1.2978,
               -0.6785, 0.8479, 3.6714, 0.4402, 0.0167, -0.5231, -0.3850,
               2.1160, 2.0427, -1.5826, 0.2996, 0.5018, 2.9571, 1.1696,
               -1.7493, 6.1279, -1.3406, 2.5413, -2.7617, -3.4235, -2.8035,
               -3.5442, 5.4941, 0.3233, 5.4864, 2.0699, 2.1345, 1.0159,
               -5.3307, 4.4847, -0.4996, -1.9334, 0, -2.2974, 2.8907, 0],
        "Pc": [0.019904, 0.010558, 0.001315, -0.010404, 0.025014, 0.017865,
               0.022319, 0.012590, 0.002044, 0.031270, 0.007542, 0.002136,
               0.019360, 0.012200, 0.002769, 0.005148, -0.007444, 0.025073,
               0.017841, 0.014091, 0.029020, 0.021836, 0.013797, 0.020440,
               0.015135, 0.009857, 0.009011, 0.012558, 0.010694, 0.012589,
               0.010390, -0.000462, 0.015874, 0.004917, 0.001120, 0.029565,
               0.025653, 0.036133, 0.011507, 0.019789, 0.011360, 0.003086,
               0.026808, 0, 0.034935, 0.013135, 0.020974, 0.012241, 0.015050,
               0.013572, 0.002753, -0.001771, 0.014827, 0.004115, 0.016004,
               0.013027, 0, 0.044232, 0.012884, 0.004673, 0.011294, 0.035446,
               0, 0.039004, 0.014434, 0.004266, 0, 0, 0.040149, 0, 0, 0.025435,
               0, 0.016048, 0.011105, 0, 0, 0, 0.000488, 0.001410, -0.001849,
               -0.005198, -0.013230, 0.003714, 0.001171, 0.000424, 0.002257,
               -0.009799, 0.004186, -0.000183, 0.003538, 0.005675, -0.002546,
               0.005175, 0.003659, 0.001474, -0.002303, 0.003818, -0.002481,
               0.004920, 0.000344, 0.000659, 0.001067, -0.004877, -0.000541,
               -0.004393, 0.000178, 0.005052, 0.006917, 0.001408, 0.002148,
               -0.005947, -0.000878, -0.002249, 0, 0.000319, -0.004305,
               0.009027, 0.008247, 0],
        "vc": [0.07504, 0.05576, 0.03153, -0.00034, 0.11648, 0.09541, 0.09183,
               0.07327, 0.07618, 0.14831, 0.04215, 0.03985, 0.10364, 0.10099,
               0.07120, 0.03897, 0.03162, 0.13396, 0.11195, 0.08635, 0.15890,
               0.13649, 0.10565, 0.08746, 0.07286, 0.05865, 0.06858, 0.13128,
               0.07527, 0.12152, 0.09956, 0.09165, 0.12598, 0.06705, 0.06358,
               0.24831, 0.17027, 0.15831, 0.10188, 0.11564, 0.10350, 0.07922,
               0.16951, 0, 0.21031, 0.10158, 0.16531, 0.14227, 0.14258,
               0.10252, 0.10814, 0.08281, 0.09331, 0.07627, 0.05687, 0.05672,
               0, 0.11480, 0.09519, 0, 0.08588, 0.18212, 0, 0.14753, 0.03783,
               0.14431, 0, 0, 0.25031, 0, 0, 0.16754, 0, 0.13021, 0.11650, 0,
               0, 0, 0.00400, 0.00572, -0.00398, -0.01081, -0.02300, -0.00014,
               -0.00851, -0.00866, 0.01636, -0.02700, -0.00781, -0.00098,
               0.00281, 0.00826, -0.01755, 0.00227, -0.00664, -0.00510,
               -0.00122, -0.01966, 0.00664, 0.00559, -0.00415, -0.00293,
               -0.00591, -0.00144, 0.02605, -0.00777, 0.01511, 0.00397,
               -0.02297, 0.00433, 0.00580, -0.01380, 0.00297, -0.00045, 0,
               -0.00596, 0.00507, 0, -0.00832, -0.00341, 0],
        "tb": [0.8894, 0.9225, 0.6033, 0.2878, 1.7827, 1.8433, 1.7117, 1.7957,
               1.8881, 3.1243, 0.9297, 1.6254, 1.9669, 1.9478, 1.7444, 3.2152,
               4.4014, 3.5668, 3.8967, 2.8526, 3.6360, 3.3953, 3.1459, 2.2536,
               1.6249, 1.1557, 2.5892, 3.1656, 2.5983, 3.1376, 2.6127, 1.5780,
               2.1647, 1.2171, 5.4736, 6.2800, 5.9234, 5.0525, 5.8337, 2.9637,
               2.6948, 2.2073, 3.9300, 3.5600, 4.5797, 2.6293, 5.7619, 5.0767,
               6.0837, 3.2914, 3.6650, 2.6495, 2.3678, 2.5645, 1.7824, 0.9442,
               7.2644, 1.2880, 0.6115, 1.1739, 2.6446, 2.8881, 2.3086, 1.9163,
               1.0081, 10.3428, 0, 0, 7.6904, 0, 6.7822, 5.5566, 5.4248,
               3.6796, 3.6763, 2.6812, 5.7093, 5.8260, -0.1157, -0.0489,
               0.1798, 0.3189, 0.7273, 0.4745, 0.3563, 0.1919, 0.1957, 0.3489,
               0.1589, 0.0668, -0.1406, -0.0900, 0.0511, 0.6884, -0.1074,
               0.0224, 0.0920, 0.5580, 0.0735, -0.1552, 0.7801, -0.2383,
               0.4456, -0.1977, 0.0835, -0.5385, -0.6331, 1.4108, -0.0690,
               1.0682, 0.4247, 0.2499, 0.1134, -0.2596, 0.4408, -0.1168,
               -0.3201, -0.4453, -0.6776, -0.3678, 0],
        "tf": [0.4640, 0.9246, 0.3557, 1.6479, 1.6472, 1.6322, 1.7899, 2.0018,
               5.1175, 3.3439, 1.4669, 0.2098, 1.8635, 0.4177, -1.7567, 3.5979,
               13.7349, 4.8776, 5.6622, 4.2927, 4.0823, 3.5572, 4.2250, 2.9248,
               2.0695, 4.0352, 4.5047, 6.7684, 4.1187, 4.5341, 6.0609, 3.4100,
               4.0580, 0.9544, 10.1031, 0, 12.6275, 4.1859, 11.5630, 3.3376,
               2.9933, 9.8409, 5.1638, 0, 10.2337, 2.7336, 5.5424, 4.9738,
               8.4724, 3.0044, 4.6089, 3.7442, 3.9106, 9.5793, 1.5598, 2.5015,
               0, 3.2411, 0, 0, 3.4448, 7.4756, 0, 2.7523, 1.9623, 31.2786, 0,
               0, 11.3770, 0, 0, 0, 0, 5.0506, 3.1468, 0, 0, 0, 0.0381,
               -0.2355, 0.4401, -0.4923, 6.0650, 1.3772, 0, 0.6824, 1.5656,
               6.9709, 1.9913, 0.2476, -0.5870, -0.2361, -2.8298, 1.4880,
               2.0547, -0.2951, -0.2986, 0.7143, -0.6697, -3.1034, 28.4324,
               0.4838, 0.0127, -2.3598, -2.0198, -0.5480, 0.3189, 0.9124,
               9.5209, 2.7826, 2.5114, 1.0729, 0.2476, 0.1175, -0.2914,
               -0.0514, -1.6425, 0, 2.5832, -1.5511, 0],

        # Table 2 in [3]_
        "hf": [-45.947, -20.763, -3.766, 17.119, 53.712, 69.939, 64.145,
               82.528, 104.293, 197.322, 11.189, 27.016, -19.243, 9.404,
               27.671, -181.422, -164.609, -182.329, -164.410, -129.2,
               -389.737, -359.258, -332.822, -163.569, -151.143, -129.488,
               -140.313, -15.505, 3.320, 5.432, 23.101, 26.718, 54.929, 69.885,
               20.079, 134.062, 139.758, 88.298, -396.242, -73.568, -63.795,
               -57.795, -82.921, 0, -107.188, -16.752, -66.138, -59.142,
               -7.365, -8.253, 57.546, 1.834, 220.803, 227.368, -36.097,
               -161.740, 0, -679.195, 0, 0, -313.545, -258.960, 0, -446.835,
               -223.398, -203.188, -67.778, -182.096, -189.888, -46.562, 0,
               -344.125, 0, -2.084, 18.022, 0, 0, 0, -0.860, -1.338, 6.771,
               7.205, 14.271, 104.800, 99.455, 13.782, -9.660, 15.465, -8.392,
               0.474, 1.472, 4.504, 1.252, -2.792, -2.092, 0.975, 4.753,
               14.145, -3.173, 1.279, 12.245, -7.807, 37.462, -16.097, -9.874,
               -3.887, -24.125, 0.366, -16.333, -2.992, 2.855, 0.351, -8.644,
               1.532, -0.329, 0, 11.989, 0, 12.285, 11.207, 11.740],
        "gf": [-8.030, 8.231, 19.848, 37.977, 84.926, 92.900, 88.402, 93.745,
               116.613, 221.308, 22.533, 30.485, 22.505, 41.228, 52.948,
               -158.589, -132.097, -131.366, -132.386, -107.858, -318.616,
               -291.188, -288.902, -105.767, -101.563, -92.099, -90.883,
               58.085, 63.051, 82.471, 95.888, 85.001, 128.602, 132.756,
               68.861, 199.958, 199.288, 121.544, -349.439, -33.373, -31.502,
               -25.261, -35.814, 0, -53.332, -0.596, 17.963, 18.088, 60.161,
               16.731, 46.945, -1.721, 217.003, 216.328, -28.148, -144.549, 0,
               -626.580, 0, 0, -281.495, -209.337, 0, -392.975, -212.718,
               -136.742, 0, 0, -65.642, 0, 0, -241.373, 0, 30.222, 38.346, 0,
               0, 0, 0.297, -0.399, 6.342, 7.466, 16.224, 94.564, 92.573,
               5.733, -8.180, 20.597, -5.505, 0.950, 0.699, 1.013, 1.041,
               -1.062, -1.359, 0.075, 0, 23.539, -2.602, 2.149, 10.715, -6.208,
               29.181, -11.809, -7.415, -6.770, -20.770, 3.805, -5.487, -1.600,
               1.858, 8.846, -13.167, -0.654, -2.091, 0, 12.373, 0, 14.161,
               12.530, 0],
        "hv": [4.116, 4.650, 2.771, 1.284, 6.714, 7.370, 6.797, 8.178, 9.342,
               12.318, 4.098, 12.552, 9.776, 10.185, 8.834, 24.529, 40.246,
               18.999, 20.041, 12.909, 22.709, 17.759, 0, 10.919, 7.478, 5.708,
               11.227, 14.599, 11.876, 14.452, 14.481, 0, 6.947, 6.918, 28.453,
               31.523, 31.005, 23.340, 43.046, 13.780, 11.985, 9.818, 19.208,
               17.574, 0, 11.883, 30.644, 26.277, 0, 14.931, 14.364, 11.423,
               7.751, 11.549, 0, 4.877, 0, 8.901, 1.860, 8.901, 0, 13.322, 0,
               8.301, 0, 0, 0, 51.787, 0, 0, 0, 0, 0, 16.921, 17.117, 13.265,
               27.966, 0, 0.292, -0.720, 0.868, 1.027, 2.426, 0, 0, -0.568,
               -0.905, -0.847, 2.057, -0.073, -0.369, 0.345, -0.114, 0, 0.207,
               -0.668, 0.071, 0.744, -3.410, 0, 8.502, -3.345, 0, 1.517, 0,
               -1.398, 0.320, -3.661, 4.626, 0, 0, 2.311, 0, 0, 0.972, 0, 0, 0,
               -7.488, -4.864, 0],

        # Table 1 in [4]_
        "w": [0.29602, 0.14691, -0.07063, -0.35125, 0.40842, 0.25224, 0.22309,
              0.23492, -0.21017, 0.73865, 0.15188, 0.02725, 0.33409, 0.14598,
              -0.08807, 1.52370, 0.73657, 1.01522, 0.63264, 0.96265, 1.13257,
              0.75574, 0.76454, 0.52646, 0.44184, 0.21808, 0.50922, 0.79963, 0,
              0.95344, 0.55018, 0.38623, 0.38447, 0.07508, 0.79337, 0, 0, 0,
              1.67037, 0.57021, 0, 0, 0.71592, 0, 0.61662, 0, 0, 0, 0, 0,
              0.23323, 0.27778, 0.61802, 0, 0, 0.26254, 0, 0.50023, 0, 0, 0,
              0.50260, 0, 0.54685, 0.43796, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.42753,
              0, 0, 0, 0.01740, 0.01922, -0.00475, -0.02883, -0.08632, 0.17563,
              0.22216, 0.16284, -0.03065, -0.02094, 0.01648, 0.00619, -0.01150,
              0.02778, -0.11024, -0.11240, 0, -0.20789, -0.16571, 0, 0,
              0.08774, 0, -0.26623, 0, 0.91939, 0, 0.03654, 0.21106, 0, 0, 0,
              0, -0.13106, 0, 0, -0.01509, 0, 0, 0, -0.03078, 0.00001, 0],
        "vliq": [0.02614, 0.01641, 0.00711, -0.00380, 0.03727, 0.02692,
                 0.02697, 0.01610, 0.00296, 0.04340, 0.01317, 0.00440, 0.02888,
                 0.01916, 0.00993, 0.00551, 0.01133, 0.03655, 0.02816, 0.02002,
                 0.04500, 0.03567, 0.02667, 0.03274, 0.02311, 0.01799, 0.02059,
                 0.02646, 0.01952, 0.02674, 0.02318, 0.01813, 0.01913, 0.01683,
                 0.01365, 0.06082, 0.05238, 0.03313, 0.02232, 0.03371, 0.02663,
                 0.02020, 0.04682, 0, 0.06202, 0.02414, 0.03375, 0.02620,
                 0.02505, 0.03446, 0.02791, 0.02143, 0, 0.01451, 0.01533,
                 0.01727, 0, 0, 0, 0, 0.01917, 0.05384, 0, 0.05383, 0, 0, 0, 0,
                 0.05477, 0, 0, 0.04104, 0, 0.03484, 0.02732, 0, 0, 0, 0.00133,
                 0.00179, -0.00203, -0.00243, -0.00744, 0, 0, 0.00213, 0.00063,
                 -0.00519, -0.00188, 0.00009, 0.00012, 0.00142, -0.00107, 0,
                 -0.00009, -0.00030, -0.00108, -0.00111, -0.00036, -0.00050,
                 0.00777, 0.00083, 0.00036, 0.00198, 0.00001, -0.00092,
                 0.00175, 0.00235, -0.00250, 0.00046, 0, -0.00179, -0.00206,
                 0.01203, -0.00023, 0, -0.0058, 0, 0.00178, 0.00171, 0],

        # Table C-2 in [1]_
        "cpa": [35.1152, 22.6346, 8.9272, 0.3456, 49.2506, 35.2248, 37.6299,
                21.3528, 10.2797, 66.0574, 16.3794, 10.4283, 42.8569, 32.8206,
                19.9504, 27.2107, 39.7712, 59.3032, 0, 40.7501, 66.8423, 0,
                51.5048, 50.5604, 39.5784, 25.6750, 0, 57.6861, 44.1122,
                53.7012, 44.6388, 0, 41.4064, 30.1561, 47.1311, 84.7602, 0,
                58.2837, 46.5577, 48.4648, 36.5885, 29.1848, 60.8262, 56.1685,
                78.6054, 33.6450, 63.7851, 51.1442, 0, 58.2445, 29.1815,
                28.0260, 45.9768, 26.7371, 25.8094, 30.1696, 0, 63.2024,
                44.3567, 0, 0, 0, 0, 0, 22.2082, 0, 0, 0, 0, 0, 0, 0, 0,
                57.7670, 45.0314, 40.5275, 80.3010, 0, 0.5830, 0.3226, 0.9668,
                -0.3082, -0.1201, 8.5546, 3.1721, -5.9060, -3.9682, -3.2746,
                2.6142, -1.3913, 0.2630, 6.5145, 4.1707, 0, 0, 3.7978, 0, 0, 0,
                0, -15.7667, 0, 0, -6.4072, 0, 2.4484, -1.5252, 0, 0, 0, 0, 0,
                0, 0, -2.7407, 0, -1.6978, 0, -2.2923, -0.3162, 0],
        "cpb": [39.5923, 45.0933, 59.9786, 74.0368, 59.3840, 62.1924, 62.1285,
                66.3947, 65.5372, 69.3936, 32.7433, 25.3634, 65.6464, 70.4153,
                81.8764, 2.7609, 35.5676, 67.8149, 0, 19.6990, 102.4553, 0,
                44.4133, 38.9681, 41.8177, 24.7281, 0, 64.0768, 77.2155,
                71.7948, 68.5041, 0, 85.0996, 81.6814, 51.3326, 177.2513, 0,
                49.6388, 48.2322, 37.2370, 47.6004, 52.3817, 41.9908, 46.9337,
                32.1318, 23.2759, 83.4744, 94.2934, 0, 46.9958, -9.7846,
                -7.1651, 20.6417, 21.7676, -5.2241, 26.9738, 0, 51.9366,
                44.5875, 0, 0, 0, 0, 0, -2.8385, 0, 0, 0, 0, 0, 0, 0, 0,
                44.1238, 55.1432, 55.0141, 132.7786, 0, -1.2002, 2.1309,
                -2.0762, 1.8969, 4.2846, -22.9771, -10.0834, -1.8710, 17.7889,
                32.1670, 4.4511, -1.5496, -2.3428, -17.5541, -3.1964, 0, 0,
                -7.3251, 0, 0, 0, 0, -0.1174, 0, 0, 15.2583, 0, -0.0765,
                -7.6380, 0, 0, 0, 0, 0, 0, 0, 11.1033, 0, 1.0477, 0, 3.1142,
                2.3711, 0],
        "cpc": [-9.9232, -15.7033, -29.5143, -45.7878, -21.7908, -24.8156,
                -26.0637, -29.3703, -30.6057, -25.1081, -13.1692, -12.7283,
                -21.0670, -28.9361, -40.2864, 1.3060, -15.5875, -20.9948, 0,
                -5.4360, -43.3306, 0, -19.6155, -4.7799, -11.0837, 4.2419, 0,
                -21.0480, -33.5086, -22.9685, -26.7106, 0, -35.6318, -36.1441,
                -25.0276, -72.3213, 0, -15.6291, -20.4868, -13.0635, -22.8148,
                -30.8526, -20.4091, -31.3325, -19.4033, -12.2406, -35.1171,
                -45.2029, 0, -10.5106, 3.4554, 2.4332, -8.3297, -6.4481,
                1.4542, -13.3722, 0, -28.6308, -23.2820, 0, 0, 0, 0, 0, 1.2679,
                0, 0, 0, 0, 0, 0, 0, 0, -9.5565, -18.7776, -31.7190, -58.3241,
                0, -0.0584, -1.5728, 0.3148, -1.6454, -2.0262, 10.7278, 4.9674,
                4.2945, -3.3639, -17.8246, -5.9808, 2.5899, 0.8975, 10.6977,
                -1.1997, 0, 0, 2.5312, 0, 0, 0, 0, 6.1191, 0, 0, -8.3149, 0,
                0.1460, 8.1795, 0, 0, 0, 0, 0, 0, 0, -11.0878, 0, 0.2002, 0,
                -1.4995, -1.4825, -0.0584],

        # Custom dict for molecular weight and empiric formula calculation
        "txt": [("CH3", ),                  # 0
                ("CH2", ),
                ("CH", ),
                ("C", ),
                ("CH2=CH", ),
                ("CH=CH", ),
                ("CH2=C", ),
                ("CH=C", ),
                ("C=C", ),
                ("CH2=C=CH", ),
                ("-CH (Aromatic)", ),       # 10
                ("=C (Aromatic)", ),
                ("-CCH3 (Aromatic)", ),
                ("-CCH2 (Aromatic)", ),
                ("-CCH (Aromatic)", ),
                ("-OH", ),
                ("-COH (Aromatic)", ),
                ("CH3CO", ),
                ("CH2CO", ),
                ("CHO", ),
                ("CH3COO", ),               # 20
                ("CH2COO", ),
                ("HCOO", ),
                ("CH3O", ),
                ("CH2O", ),
                ("CH-O", ),
                ("FCH2O", ),
                ("CH2NH2", ),
                ("CHNH2", ),
                ("CH3NH", ),
                ("CH2NH", ),                # 30
                ("CHNH", ),
                ("CH3N", ),
                ("CH2N", ),
                ("=CNH2 (Aromatic)", ),
                ("C5H4N", ),
                ("C5H3N", ),
                ("CH2CN", ),
                ("COOH", ),
                ("CH2Cl", ),
                ("CHCl", ),                 # 40
                ("CCl", ),
                ("CHCl2", ),
                ("CCl3", ),
                ("CCl2", ),
                ("=CCl (Aromatic)", ),
                ("CH2NO2", ),
                ("CHNO2", ),
                ("=CNO2 (Aromatic)", ),
                ("CH2SH", ),
                ("I", ),                    # 50
                ("Br", ),
                ("CH≡C", ),
                ("C≡C", ),
                ("Cl-C=C", ),
                ("=CF (Aromatic)", ),
                ("HCONCH2CH2", ),
                ("CF3", ),
                ("CF2", ),
                ("CF", ),
                ("COO", ),                  # 60
                ("CCl2F", ),
                ("HCClF", ),
                ("CClF2", ),
                ("F (others)", ),
                ("CONH2", ),
                ("CONHCH3", ),
                ("CONHCH2", ),
                ("CONCH3CH3", ),
                ("CONCH2CH2", ),
                ("CONCH2CH2", ),            # 70
                ("C2H5O2", ),
                ("C2H4O2", ),
                ("CH3S", ),
                ("CH2S", ),
                ("CHS", ),
                ("C4H3S", ),
                ("C4H2S", ),

                # Second order
                ("CH(CH3)2", ),
                ("C(CH3)3", ),
                ("CHCH3CHCH3", ),           # 80
                ("CH(CH3)C(CH3)2", ),
                ("C(CH3)2C(CH3)2", ),
                ("3 membered ring", ),
                ("4 membered ring", ),
                ("5 membered ring", ),
                ("6 membered ring", ),
                ("7 membered ring", ),
                ("CHn=CHm-CHp=CHk (m, p (0,1); k, n (0,2)", ),
                ("CH3-CHm=CHn (m (0,1); n (0,2))", ),
                ("CH2-CHm=CHn (m (0,1); n (0,2))", ),           # 90
                ("CH-CHm=CHn (m (0,1); n (0,2))", ),
                ("Alicyclic side-chain CcyclicCm", ),
                ("CH3CH3", ),
                ("CHCHO or CCHO", ),
                ("CH3COCH2", ),
                ("CH3COCH or CH3COC", ),
                ("Ccyclic=O", ),
                ("ACCHO", ),
                ("CHCOOH or CCOOH", ),
                ("ACCOOH", ),                                   # 100
                ("CH3COOCH or CH3COOC", ),
                ("COCH2COO or COCHCOO or COCCOO", ),
                ("CO-O-CO", ),
                ("ACCOO", ),
                ("CHOH", ),
                ("COH", ),
                ("CHm(OH)CHn(OH) (0<m,n<2)", ),
                ("CHm cyclic-OH (0<m<1)", ),
                ("CHn(OH)CHm(NHp) (0<m<1); (0<n,p<2)", ),
                ("CHm(NH2)CHn(NH2) (0<m,n<2)", ),               # 110
                ("CHm cyclic-NHp-CHn cyclic (0<n,m,p<1)", ),
                ("CHn-O-CHm=CHp (0<m<1); (0<n,p<2)", ),
                ("AC-O-CHm (0<m<3)", ),
                ("CHm cyclic-S-CHn cyclic (0<n,m<1)", ),
                ("CHn=CHm-F (0<m<1); (0<n<2)", ),
                ("CHn=CHm-Br (0<m<1); (0<n<2)", ),
                ("CHn=CHm-I (0<m<1); (0<n<2)", ),
                ("ACBr", ),
                ("ACI", ),
                ("CHm(NH2)-COOH (0<m<2)", )]                    # 120
            }

    FirstOrder = 78
    SecondOrder = 121

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        tc = Pc = tf = tb = vc = w = gf = hf = hv = vliq = cpa = cpb = cpc = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            Pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]
            tf += c*self._coeff["tf"][i]
            tb += c*self._coeff["tb"][i]
            hf += c*self._coeff["hf"][i]
            gf += c*self._coeff["gf"][i]
            hv += c*self._coeff["hv"][i]
            w += c*self._coeff["w"][i]
            vliq += c*self._coeff["vliq"][i]
            cpa += c*self._coeff["cpa"][i]
            cpb += c*self._coeff["cpb"][i]
            cpc += c*self._coeff["cpc"][i]

        # Table 5 with functions
        self.Tc = unidades.Temperature(181.128*log(tc))
        self.Pc = unidades.Pressure((Pc+0.10022)**-2+1.3705, "bar")
        self.Vc = unidades.SpecificVolume((vc-0.00435)/self.M)
        self.Tf = unidades.Temperature(102.425*log(tf))
        self.Tb = unidades.Temperature(204.359*log(tb))
        self.Hf = unidades.Enthalpy((hf+10.835)/self.M, "kJg")
        self.Gf = unidades.Enthalpy((gf-14.828)/self.M, "kJg")
        self.Hv = unidades.Enthalpy((hv+6.829)/self.M, "kJg")

        # Eq 2 in [4]_
        self.f_acent = 0.4085*log(w+1.1507)**(1/0.5050)
        # Eq 3 in [4]_
        self.Vliq = unidades.SpecificVolume((vliq+0.01211)*self.M)

        # Ideal gas specific heat expression from [1]_
        self.cp = [cpa-19.7779, cpb+22.5981, cpc-10.7983, 0, 0, 0]

        GroupContribution.calculo(self)

    def _Cp0(self, T):
        """Ideal gas specific heat calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        cp0 : float
            Ideal gas specific heat, [J/kgK]
        """
        Tita = (T-298)/700
        cp0 = 0
        for i, c in enumerate(self.cp):
            cp0 += c*Tita**i

        return unidades.SpecificHeat(cp0/self.M, "JgK")


@refDoc(__doi__, [1, 5])
class Wilson(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Wilson-Jasperson procedure (1994)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    ring : int
        Ring in the atom, [-]
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2-4 in [1]_, 2-ethylphenol critical properties

    >>> c1 = Wilson(Tb=477.67, ring=1)
    >>> c2 = Wilson(Tb=477.67, ring=1)
    >>> c1(group=[3, 0, 5], contribution=[8, 10, 1])
    >>> c2(group=[3, 0, 5, 42], contribution=[8, 10, 1, 1])
    >>> "%0.1f %0.2f %0.1f" % (c1.Tc, c1.Pc.bar, c2.Tc)
    '702.9 37.94 693.6'
    >>> c1.formula
    'C8H10O'
    """
    __title__ = "Wilson-Jasperson (1996)"
    kwargs = GroupContribution.kwargs.copy()
    kwargs["ring"] = 0

    _coeff = {
        "tc": [0.002793, 0.320000, 0.019000, 0.008532, 0.019181, 0.020341,
               0.008810, 0.036400, 0.088000, 0.020000, 0.012000, 0.007271,
               0.011151, 0.016800, 0.014000, 0.018600, 0.059000, 0.031000,
               0.007000, 0.010300, 0.012447, 0.013300, -0.027000, 0.175000,
               0.017600, 0.007000, 0.020000, 0.010000, 0.000000, 0.005900,
               0.017000, -0.027500, 0.219000, 0.013000, 0.011000, 0.014000,
               -0.050000, 0.000000, 0.000000, 0.007000, 0.015000, 0.0350,
               0.0100, -0.0075, -0.0040, 0.0000, -0.0550, 0.0170, -0.0150,
               0.0170, -0.0200, 0.0020, 0.0000, -0.0250],
        "Pc": [0.12660, 0.43400, 0.91000, 0.72983, 0.44805, 0.43360, 0.32868,
               0.12600, 6.05000, 1.34000, 1.22000, 1.04713, 0.97711, 0.79600,
               1.19000, 0.0, 0.0, 1.42000, 2.68000, 1.20000, 0.97151, 1.11000,
               0.0, 1.11000, 2.71000, 1.69000, 1.95000, 0.0, 0.43000, 1.315930,
               1.66000, 6.33000, 1.07000, 0.0, 1.08000, 0.0, 0.0, -0.08000,
               0.69000, 2.05000, 2.04000, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00,
               0.50, 0.00, 1.50, 1.00, 0.00, 0.00, -0.50],

        "txt": [("H",),
                ("He",),
                ("B",),
                ("C",),
                ("N",),
                ("O",),
                ("F",),
                ("Ne",),
                ("Al",),
                ("Si",),
                ("P",),
                ("S",),
                ("Cl",),
                ("Ar",),
                ("Ti",),
                ("V",),
                ("Ga",),
                ("Ge",),
                ("As",),
                ("Se",),
                ("Br",),
                ("Kr",),
                ("Rb",),
                ("Zr",),
                ("Nb",),
                ("Mo",),
                ("Sn",),
                ("Sb",),
                ("Te",),
                ("I",),
                ("Xe",),
                ("Cs",),
                ("Hf",),
                ("Ta",),
                ("W",),
                ("Re",),
                ("Os",),
                ("Hg",),
                ("Bi",),
                ("Rn",),
                ("U",),

                # 2nd Order term
                ("-OH, C4 or less",),
                ("-OH, C5 or more",),
                ("-O-",),
                ("-NH2, >NH, >N-",),
                ("-CHO",),
                (">CO",),
                ("-COOH",),
                ("-COO-",),
                ("-CN",),
                ("-NO2",),
                ("Organic Halides (once / molecule)",),
                ("-SH, -S-, -SS-",),
                ("Siloxane bond",)]}

    FirstOrder = 41
    SecondOrder = 54

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        tc = Pc = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            Pc += c*self._coeff["Pc"][i]

        Nr = self.kwargs["ring"]
        self.Tc = unidades.Temperature(self.Tb/(0.048271-0.019846*Nr+tc)**0.2)
        self.Pc = unidades.Pressure(0.0186233*self.Tc/(-0.96601+exp(
            -0.00922295-0.0290403*Nr+0.041*Pc)), "bar")

        GroupContribution.calculo(self)


@refDoc(__doi__, [1, 6, 7])
class Marrero(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Marrero-Pardillo procedure (1999)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    Tb : float, optional
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2-6 in [1]_, 2-ethylphenol critical properties

    >>> cmp = Marrero(Tb=477.67, group=[1, 36, 129, 130, 132, 140, 148],
    ... contribution=[1, 1, 1, 2, 2, 1, 1])
    >>> "%0.1f %0.1f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '699.8 42.2 378.8'
    >>> cmp.formula
    'C8H10O'

    Example 2-7 in [1]_, butanols

    >>> b1 = Marrero(group=[1, 28, 41], contribution=[1, 2, 1])
    >>> b2m1 = Marrero(group=[2, 29, 41], contribution=[2, 1, 1])
    >>> b2m2 = Marrero(group=[3, 79], contribution=[3, 1])
    >>> b2 = Marrero(group=[1, 2, 29, 62], contribution=[1, 1, 1, 1])

    >>> "%0.2f %0.2f %0.2f %0.2f" % (b1.Tc, b2m1.Tc, b2m2.Tc, b2.Tc)
    '565.67 555.80 504.45 534.96'
    >>> "%0.2f %0.2f %0.2f %0.2f" % (
    ... b1.Pc.bar, b2m1.Pc.bar, b2m2.Pc.bar, b2.Pc.bar)
    '44.86 45.04 41.30 43.40'
    >>> "%0.1f %0.1f %0.1f %0.1f" % (b1.Vc.ccg*b1.M, b2m1.Vc.ccg*b2m1.M,
    ...  b2m2.Vc.ccg*b2m2.M, b2.Vc.ccg*b2.M)
    '272.1 267.3 275.1 277.9'

    Example 2-13 in [1]_, 2,4-dimethylphenol

    >>> cmp = Marrero(group=[9, 130, 132, 133, 140, 148],
    ... contribution=[2, 3, 1, 1, 1, 1], M=122.167)
    >>> "%0.2f" % cmp.Tb
    '489.39'

    Example 2-14 in [1]_, cycloalkanes

    >>> c7 = Marrero(group=[111], contribution=[7])
    >>> mc6 = Marrero(group=[7, 111, 112], contribution=[1, 4, 2])
    >>> ec5 = Marrero(group=[1, 34, 111, 112], contribution=[1, 1, 3, 2])
    >>> c5 = Marrero(group=[7, 111, 112], contribution=[2, 1, 4])
    >>> t5 = Marrero(group=[7, 111, 112], contribution=[2, 1, 4])
    >>> "%0.2f %0.2f %0.2f %0.2f %0.2f" % (c7.Tb, mc6.Tb, ec5.Tb, c5.Tb, t5.Tb)
    '377.53 375.85 384.47 374.17 374.17'

    Table 8a in [6]_ for 1,3,5-Trichlorotrifluorobenzene

    >>> cmp = Marrero(group=[138, 140, 144, 145], contribution=[3, 3, 3, 3])
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '471.9 685.0 32.4 452.3'

    Table 8b in [6]_ for ethyl acrylate

    >>> cmp = Marrero(group=[84, 1, 46, 100], contribution=[1, 1, 1, 1])
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '373.9 553.6 36.6 325.1'

    Table 8c in [6]_ for isoquinoline

    >>> cmp = Marrero(group=[129, 138, 131, 136, 132, 133],
    ... contribution=[3, 1, 1, 1, 1, 4])
    >>> "%0.1f %0.1f %0.1f" % (cmp.Tb, cmp.Tc, cmp.Vc.ccg*cmp.M)
    '519.5 787.3 405.1'

    Table 8d in [6]_ for m-Terphenyl

    >>> cmp = Marrero(group=[129, 130, 141, 132, 133],
    ... contribution=[5, 4, 2, 5, 4], Tb=638)
    >>> "%0.1f %0.1f %0.1f %0.1f" % (
    ... cmp.Tb, cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '638.0 907.3 33.9 764.3'

    Table 3 in [7]_ for o-phthalate

    >>> cmp = Marrero(group=[129, 132, 133, 138, 153, 46, 28, 1],
    ... contribution=[2, 1, 2, 1, 2, 2, 4, 2])
    >>> "%0.2f" % cmp.mu.muPas
    '18.58'
    """
    __title__ = "Marrero-Pardillo (1999)"

    _coeff = {
        # Table 5 from [6]_
        "tc": [-0.0213, -0.0227, -0.0223, -0.0189, 0.8526, 0.1792, 0.3818,
               -0.0214, 0.1117, 0.0987, -0.0370, -0.9141, -0.9166, -0.9146,
               -0.0876, -0.0205, -0.0362, -0.0606, -0.0890, 0.0267, -0.0974,
               -0.0397, -0.0313, -0.0199, -0.0766, -0.0591, -0.9192, -0.0181,
               -0.0206, -0.0134, -0.0098, 0.8636, 0.1874, 0.4160, -0.0149,
               0.1193, 0.1012, -0.0255, -0.0162, -0.0205, -0.0210, -0.0786,
               -0.0205, -0.0256, -0.0267, -0.0932, 0.0276, -0.0993, -0.0301,
               -0.0248, -0.0161, -0.0654, -0.0137, -0.0192, -0.0039, 0.0025,
               0.8547, 0.1969, 0.0025, 0.1187, -0.0200, -0.0142, -0.0757,
               -0.0162, -0.0194, -0.0406, -0.0918, -0.1054, -0.0286, -0.0158,
               0.0084, 0.8767, 0.2061, 0.0207, 0.0049, 0.1249, -0.0176,
               -0.0133, -0.0084, -0.0780, -0.0156, -0.0114, -0.1008, -0.9129,
               -0.8933, -0.4158, -0.0123, -1.7660, -1.2909, -0.8945, 1.7377,
               1.0731, 1.2865, 0.9929, 0.8623, 0.8613, 0.8565, 0.8246, 0.7862,
               0.8818, 0.7780, 0.8122, -0.8155, -0.4009, 0.3043, 0.1868,
               0.1886, -0.0159, -0.0288, -0.4222, -0.7958, -0.0098, -0.0093,
               -0.1386, 0.0976, 0.1089, -0.0092, -0.0148, -0.0139, -0.0071,
               -0.0055, -0.1341, 0.0, 0.0, -0.0218, -0.0737, 0.0329, 0.0,
               -0.0314, -0.2246, -0.3586, 0.3913, 0.2089, 0.2190, 0.1000,
               0.0947, -0.4067, 0.1027, -0.4848, 0.2541, 0.2318, 0.2424,
               0.1104, -0.3972, 0.1069, 0.1028, 0.1060, 0.1075, 0.0931, 0.0997,
               0.1112, 0.0919, 0.0313, 0.0241, 0.0830, 0.0978, 0.0938, 0.0768,
               -0.0191, -0.1926, -0.5728, -0.3553, -0.0422, -0.0690, -0.0781,
               -0.0301, -0.0124],
        "Pc": [-0.0618, -0.0430, -0.0376, -0.0354, 0.0654, 0.0851, -0.2320,
               -0.0396, -0.0597, -0.0746, -0.0345, -0.0231, -0.0239, -0.0241,
               -0.0180, -0.0321, -0.0363, -0.0466, -0.0499, 0.1462, -0.2290,
               -0.0288, -0.0317, -0.0348, -0.0507, -0.0385, -0.0244, -0.0305,
               -0.0272, -0.0219, -0.0162, 0.0818, 0.1010, -0.2199, -0.0265,
               -0.0423, -0.0626, -0.0161, -0.0150, -0.0140, -0.0214, -0.0119,
               -0.0184, -0.0204, -0.0210, -0.0253, 0.1561, -0.2150, -0.0214,
               -0.0203, -0.0170, -0.0329, -0.0163, -0.0173, -0.0137, -0.0085,
               0.0816, 0.1080, -0.0168, -0.0556, -0.0147, -0.0131, -0.0093,
               -0.0155, -0.0112, -0.0280, -0.2098, -0.0358, -0.0212, -0.0162,
               0.0002, 0.0953, 0.1109, 0.0213, -0.0111, -0.0510, -0.0161,
               -0.0129, -0.0121, -0.0094, -0.0103, -0.0085, -0.0455, -0.0476,
               -0.1378, -0.2709, -0.0239, -0.2291, -0.3613, -0.1202, 0.1944,
               0.2146, -0.1087, 0.0533, 0.0929, 0.0919, 0.0947, 0.0801, 0.0806,
               0.2743, -0.1007, 0.0771, -0.4920, -0.2502, 0.0705, 0.1064,
               0.1102, -0.0010, -0.0226, 0.1860, 0.3933, -0.0221, -0.0181,
               0.0081, -0.1034, -0.0527, -0.0119, -0.0177, -0.0127, 0.0,
               -0.0088, 0.0162, 0.0, 0.0, -0.0091, -0.0220, -0.0071, 0.0,
               -0.0119, 0.1542, 0.1490, 0.1356, -0.1822, -0.1324, -0.0900, 0.0,
               -0.1491, -0.0916, 0.1432, 0.0, -0.0809, -0.0792, -0.0374,
               -0.0971, -0.0504, -0.0512, -0.0548, -0.0514, -0.0388, -0.0523,
               -0.0528, -0.0597, -0.0684, -0.2573, -0.0579, -0.0471, -0.0462,
               -0.0625, -0.0125, -0.0878, 0.0, -0.0176, -0.0123, 0.0, -0.1878,
               0.0, 0.0],
        "vc": [123.2, 88.6, 78.4, 69.8, 81.5, 57.7, 65.8, 58.3, 49.0, 71.7,
               88.1, 113.8, 0.0, 0.0, 92.9, 66.0, 88.9, 128.9, 145.9, 93.3,
               108.2, 0.0, 0.0, 76.3, 147.9, 148.1, 119.7, 87.9, 56.6, 40.2,
               32.0, 50.7, 24.0, 33.9, 31.9, 0.0, 52.1, 49.3, 80.8, 101.3, 0.0,
               45.2, 34.5, 62.3, 106.1, 114.0, 69.9, 79.1, 63.3, 49.4, 32.7,
               113.5, 93.3, 57.9, 18.3, 8.6, 48.9, 4.3, 0.0, 0.0, 37.7, 68.6,
               45.6, 23.7, 39.3, 92.2, 72.3, 110.2, 39.2, 0.0, 22.7, 23.4, 8.8,
               0.0, 0.0, 0.0, 30.0, 63.7, 85.7, 40.6, 40.8, 62.1, 89.0, 105.3,
               77.4, 99.2, 68.4, 47.8, 73.6, 43.6, 42.1, 16.6, 0.0, 0.0, 41.4,
               68.7, 36.4, 0.0, 107.4, 55.2, 64.1, 107.4, 93.7, 58.1, 0.0,
               14.6, 43.3, 51.4, 87.6, 73.1, 64.3, 47.2, 47.5, 49.9, 42.5, 0.0,
               29.2, 50.7, 38.8, 0.0, 33.9, 0.0, 0.0, 0.0, 19.2, 0.0, 36.2,
               0.0, 18.4, 36.5, 34.4, 8.3, 39.3, 29.8, 40.3, 0.0, 65.9, 40.8,
               37.8, 0.0, 20.6, 51.7, -0.3, 35.6, 23.7, 60.3, 83.2, 110.2, 8.5,
               0.0, 46.3, 0.0, 100.2, 55.2, 33.2, 0.0, 0.0, 0.0, 84.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 51.2, 0.0, 0.0],
        "tb": [113.12, 194.25, 194.27, 186.41, 137.18, 182.20, 194.40, 176.16,
               180.60, 145.56, 160.83, 453.70, 758.44, 1181.44, 736.93, 228.01,
               445.61, 636.49, 1228.84, 456.92, 510.65, 443.76, 293.86, 207.75,
               891.15, 1148.58, 588.31, 409.85, 244.88, 244.14, 273.26, 201.80,
               242.47, 207.49, 238.81, 260.00, 167.85, 166.59, 517.62, 875.85,
               1262.80, 673.24, 243.37, 451.27, 648.70, 1280.39, 475.65,
               541.29, 452.30, 314.71, 240.08, 869.18, 612.31, 451.03, 291.41,
               344.06, 179.96, 249.10, 295.33, 132.66, 68.80, 438.47, 585.99,
               215.94, 434.45, 630.07, 497.58, 1270.16, 388.44, 260.32, 411.56,
               286.30, 286.42, 456.90, 340.00, 188.99, -16.64, 360.79, 610.26,
               540.38, 267.26, 373.71, 1336.54, 51.13, 205.73, 245.27, 183.55,
               334.64, 354.41, 316.46, 174.18, 228.38, 174.39, 184.20, 5.57,
               370.60, 204.81, 658.53, 1245.86, 423.86, 525.35, 761.36, 399.58,
               321.02, 250.88, -37.99, 367.05, 160.42, 120.85, 222.40, 333.26,
               201.89, 209.40, 182.74, 218.07, 106.21, 225.52, 451.74, 283.55,
               424.13, 210.66, 220.24, 254.50, 184.36, 169.17, 597.82, 348.23,
               111.51, -41.35, 112.00, 291.15, 221.55, 285.07, 237.22, 171.59,
               420.54, 321.44, 348.00, 477.77, 334.09, 180.07, 123.05, 134.23,
               174.31, -48.79, 347.33, 716.23, 1294.98, 456.25, 199.70, 437.51,
               700.06, 1232.55, 437.78, 517.75, 411.29, 422.51, 682.19, 532.24,
               1012.51, 382.25, 385.36, 387.17, 1022.45, 298.12, 673.59,
               597.59],

        # Table 1 and 2 from [7]_
        "mu": [0.4712, 0.2588, 0.2472, 0.1833, 0.0834, 0.1670, 0, -0.0926, 0,
               0.0314, 0, 0.6294, 0.7161, 0.6086, 2.6902, 0.4020, 0.5475,
               1.3881, 2.6110, 0.6410, 0.5904, 0, 1.3476, 0.1695, 1.9082,
               1.9476, 0, 0.5587, 0.0975, -0.0543, -0.0417, -0.0775, -0.0488,
               0, -0.4464, 0, -0.4327, 0, 0.7214, 0.5799, 0.3971, 2.9295,
               0.1089, 0.2605, 1.0415, 1.9798, 0.3632, 0.2089, 1.5471, 0.2857,
               -0.0047, 1.9082, 0.7255, 0.2785, -0.1470, 0.1181, 0, 0, 0,
               -0.0701, 0, 0.3633, 2.9212, -0.1711, 0.2242, 0, -0.0606, 1.6883,
               0.6492, 0, 0, 0, 0, 0, 0, 0, -0.1491, 0.2947, 0, 2.8610, 0, 0,
               0, 0, 0.2894, 0.3306, 0, 0.3096, 0.1060, 0, 0, -0.0465, 0,
               -0.2144, 0, 0.3593, -0.0982, 1.1338, 0, 0.1966, 1.1345, 1.2214,
               0.2482, 0, 0, 0, 0.1831, 0, 0, 0, 0, 0.2980, 0.2136, 0.1692,
               0.2675, -0.0704, 0.3781, 0.6534, 0.6341, 0.5436, -0.3205, 0, 0,
               0, 0.8726, 3.1712, 0, 0, 0, 0.0952, 0.2486, 0.5734, 0.3422,
               0.2130, 0.3708, 1.0147, 0.5751, 0.5128, 0.4700, 0.7030, 0.2292,
               0, 0, 0.1880, -0.1932, 0.0404, -0.0847, 0.0139, 4.3481, 0.0963,
               0.2243, 0.7089, 0, 0.1514, 1.7430, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0],

        # Custom dict for molecular weight and empiric formula calculation
        "txt": [("CH3- & CH3-",),          # 0
                ("CH3- & -CH2-",),
                ("CH3- & >CH-",),
                ("CH3- & >C<",),
                ("CH3- & =CH-",),
                ("CH3- & =C<",),
                ("CH3- & ≡C-",),
                ("CH3- & >CH- [r]",),
                ("CH3- & >C< [r]",),
                ("CH3- & =C< [r]",),
                ("CH3- & F-",),            # 10
                ("CH3- & Cl-",),
                ("CH3- & Br-",),
                ("CH3- & I-",),
                ("CH3- & -OH",),
                ("CH3- & -O-",),
                ("CH3- & >CO",),
                ("CH3- & -CHO",),
                ("CH3- & -COOH",),
                ("CH3- & -COO[-]",),
                ("CH3- & [-]COO-",),       # 20
                ("CH3- & -NH2",),
                ("CH3- & -NH-",),
                ("CH3- & >N-",),
                ("CH3- & -CN",),
                ("CH3- & -NO2",),
                ("CH3- & -SH",),
                ("CH3- & -S-",),
                ("-CH2- & -CH2-",),
                ("-CH2- & >CH-",),
                ("-CH2- & >C<",),          # 30
                ("-CH2- & =CH-",),
                ("-CH2- & =C<",),
                ("-CH2- & ≡C-",),
                ("-CH2- & >CH- [r]",),
                ("-CH2- & >C< [r]",),
                ("-CH2- & =C< [r]",),
                ("-CH2- & F-",),
                ("-CH2- & Cl-",),
                ("-CH2- & Br-",),
                ("-CH2- & I-",),           # 40
                ("-CH2- & -OH",),
                ("-CH2- & -O-",),
                ("-CH2- & >CO",),
                ("-CH2- & -CHO",),
                ("-CH2- & -COOH",),
                ("-CH2- & -COO[-]",),
                ("-CH2- & [-]COO-",),
                ("-CH2- & -NH2",),
                ("-CH2- & -NH-",),
                ("-CH2- & >N-",),          # 50
                ("-CH2- & -CN",),
                ("-CH2- & -SH",),
                ("-CH2- & -S-",),
                (">CH- & >CH-",),
                (">CH- & >C<",),
                (">CH- & =CH-",),
                (">CH- & =C<",),
                (">CH- & >CH- [r]",),
                (">CH- & =C< [r]",),
                (">CH- & F-",),            # 60
                (">CH- & Cl-",),
                (">CH- & -OH",),
                (">CH- & -O-",),
                (">CH- & >CO",),
                (">CH- & -CHO",),
                (">CH- & [-]COO-",),
                (">CH- & -COOH",),
                (">CH- & -NH2",),
                (">CH- & -NH-",),
                (">C< & >C<",),            # 70
                (">C< & =CH-",),
                (">C< & =C<",),
                (">C< & >C< [r]",),
                (">C< & >CH- [r]",),
                (">C< & =C< [r]",),
                (">C< & F-",),
                (">C< & Cl-",),
                (">C< & Br-",),
                (">C< & -OH",),
                (">C< & -O-",),            # 80
                (">C< & >CO",),
                (">C< & -COOH",),
                ("[=]CH2 & [=]CH2",),
                ("[=]CH2 & [=]CH-",),
                ("[=]CH2 & [=]C<",),
                ("[=]CH2 & =C[=]",),
                ("[=]CH- & [=]CH-",),
                ("[=]CH- & [=]C<",),
                ("[=]CH- & =C[=]",),
                ("=CH- & =CH-",),          # 90
                ("=CH- & =C<",),
                ("=CH- & ≡C-",),
                ("=CH- & =C< [r]",),
                ("=CH- & F-",),
                ("=CH- & Cl-",),
                ("=CH- & -O-",),
                ("=CH- & -CHO",),
                ("=CH- & -COOH",),
                ("=CH- & -COO[-]",),
                ("=CH- & [-]COO-",),       # 100
                ("=CH- & -CN",),
                ("[=]C< & [=]C<",),
                ("[=]C< & =C[=]",),
                ("=C< & =C< [r]",),
                ("=C< & F-",),
                ("=C< & Cl-",),
                ("=C[=] & O[=]",),
                ("CH[≡] & CH[≡]",),
                ("CH[≡] & -C[≡]",),
                ("-C[≡] & -C[≡]",),             # 110
                ("-CH2- [r] & -CH2- [r]",),
                ("-CH2- [r] & >CH- [r]",),
                ("-CH2- [r] & >C< [r]",),
                ("-CH2- [r] & =CH- [r]",),
                ("-CH2- [r] & =C< [r]",),
                ("-CH2- [r] & -O- [r]",),
                ("-CH2- [r] & >CO [r]",),
                ("-CH2- [r] & -NH- [r]",),
                ("-CH2- [r] & -S- [r]",),
                (">CH- [r] & >CH- [r]",),       # 120
                (">CH- [r] & >C< [r]",),
                (">CH- [r] & >CH- [r]",),
                (">CH- [r] & [=]C< [r]",),
                (">CH- [r] & -O- [r]",),
                (">CH- [r] & -OH",),
                (">C< [r] & >C< [r]",),
                (">C< [r] & =C< [r]",),
                (">C< [r] & F-",),
                ("[=]CH- [r] & [=]CH- [r]",),
                ("[=]CH- [r] & [=]C< [r]",),    # 130
                ("[=]CH- [r] & [=]N- [r]",),
                ("=CH- [r] & =CH- [r]",),
                ("=CH- [r] & =C< [r]",),
                ("=CH- [r] & -O- [r]",),
                ("=CH- [r] & -NH- [r]",),
                ("=CH- [r] & =N- [r]",),
                ("=CH- [r] & -S- [r]",),
                ("[=]C< [r] & [=]C< [r]",),
                ("[=]C< [r] & -N[=] [r]",),
                ("=C< [r] & =C< [r]",),         # 140
                ("=C< [r] & =C< [r]",),
                ("=C< [r] & -O- [r]",),
                ("=C< [r] & =N- [r]",),
                ("=C< [r] & F-",),
                ("=C< [r] & Cl-",),
                ("=C< [r] & Br-",),
                ("=C< [r] & I-",),
                ("=C< [r] & -OH",),
                ("=C< [r] & -O-",),
                ("=C< [r] & >CO",),             # 150
                ("=C< [r] & -CHO",),
                ("=C< [r] & -COOH",),
                ("=C< [r] & [-]COO-",),
                ("=C< [r] & -NH2",),
                ("=C< [r] & -NH-",),
                ("=C< [r] & >N-",),
                ("=C< [r] & -CN",),
                ("Cl- & >CO",),
                ("[-]COO- & [-]COO-",),
                ("-O- [r] & =N- [r]",),         # 160
                (">CO & -O-",),
                ("-H & -CHO",),
                ("-H & -COOH",),
                ("-H & [-]COO-",),
                ("-NH- & -NH2",),
                ("-S- & -S-",)]}

    FirstOrder = 29

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        group, rest = self._decomposition()
        self.group = group

        resto = 0
        for grp, r in rest.items():
            resto += r

        if resto:
            self.msg = QApplication.translate(
                    "pychemqt",
                    "Bad definition, check input group and contribution")
            self.status = 0
        else:
            self.msg = ""
            self.status = 1
            return True

    def _decomposition(self):
        """Specific procedure to calculate the molecular weight of compound
        from group contribution"""
        group = []
        rest = {}
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            A, B = self._coeff["txt"][i][0].split(" & ")
            for x in range(c):
                for grp in (A, B):
                    # Remove the optional corchetes in multiple link definition
                    # TODO: Add support for structure check for group with
                    # several different connection, in that case the link used
                    # are the corchetes signaled
                    if "[-]" in grp:
                        grp = grp.replace("[-]", "-")
                    if "[=]" in grp:
                        grp = grp.replace("[=]", "=")
                    if "[≡]" in grp:
                        grp = grp.replace("[≡]", "≡")

                    cmp = grp
                    molecule = grp
                    # Clean additional comment of group, ring flag and other,
                    # separated of main group by spaces
                    if " " in cmp:
                        cmp = cmp.split(" ")[0]

                    # Calculate the remain appeareances of group contributions
                    # X- : X is one connected so don't appearance reamaining
                    # -X- : X has 2 connection so must be stay in another group
                    # -X<: X has 3 connection so must be stay in other 2 group
                    # >X<: X has 4 connection so must be stay in other 3 group
                    restLink = -1
                    if cmp[0] in ("-=≡"):
                        restLink += 1
                        molecule = molecule[1:]
                    elif cmp[0] == ">":
                        restLink += 2
                        molecule = molecule[1:]
                    if cmp[-1] in ("-=≡"):
                        restLink += 1
                        molecule = molecule[:-1]
                    elif cmp[-1] == "<":
                        restLink += 2
                        molecule = molecule[:-1]

                    if grp in rest and rest[grp]:
                        rest[grp] -= 1
                    else:
                        rest[grp] = restLink
                        group.append(atomic_decomposition(molecule))
        return group, rest

    def calculo(self):
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        self.Na = self._atoms()

        tc = Pc = vc = tb = mu = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb += c*self._coeff["tb"][i]
            tc += c*self._coeff["tc"][i]
            Pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]
            mu += c*self._coeff["mu"][i]

        # Table 1 equations in [6]_
        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb = unidades.Temperature(self.M**-0.404*tb+156)
        self.Tc = unidades.Temperature(self.Tb/(0.5851-0.9286*tc-tc**2))
        self.Pc = unidades.Pressure((0.1285-0.0059*self.Na-Pc)**-2, "bar")
        self.Vc = unidades.SpecificVolume((25.1+vc)/self.M, "ccg")

        # Eq 1 in [7]_
        self.mu = unidades.Viscosity(self.M**1.279*exp(mu-7.6425), "muPas")

        GroupContribution.calculo(self)


class Elliott(GroupContribution):
    """Zuppo and Elliott, Ind. Eng. Chem. Res. Submitted (1999).
    ref, chemcad propiedades fisicas pag 62
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    M: peso molecular
    Tb: Temperatura de ebullición, opcional
    SG: gravedad específica, opcional

    # >>> elliot=Elliott(group=[0, 5], contribution=[4, 1], M=72)
    # >>> print elliot.Tb, elliot.Tc
    # 333.268576405 829.20395796
    """
    __title__ = "UNIFAC (1999)"

    _coeff = {
        "tc": [0.135, 0.131, 0.077, 0.073, 0.070, -0.015, 0.070, 0.169, 0.169,
               0.169, 0.169, 0.169, 0.338, 0.069, 0.099, 0.221, 0.207, 0.136,
               0.554, 0.0, 0.0, 0.278, 0.387, 0.383, 0.299, 0.457, 0.453,
               0.305, 0.234, 0.230, 0.175, 0.140, 0.0, 0.301, 0.247, 0.306,
               0.301, 0.247, 0.148, 0.144, 0.270, 0.0, 0.433, 0.433, 0.0,
               0.512, 0.615, 0.0, 0.236, 0.178, 0.090, 0.0, 0.283, 0.196, 0.0,
               0.326, 0.0, 0.165, 0.0, 0.440, 0.440, 0.440, 0.0, 0.0, 0.203,
               0.0, 0.0, 0.056, 0.056, 0.125, 0.125, 0.0, 0.0, 0.082, 0.147,
               0.0, 0.0, 0.340, 0.222, 0.103, 0.327, 0.209, 0.205, 0.151,
               0.144, 0.245, 0.245, 0.215, 0.148, 0.0, 0.314, 0.0, 0.209,
               0.327, 0.0, 0.0, 0.0, 0.0, 0.422, 0.557, 0.553, 0.670, 0.666,
               0.662, 0.839, 0.609, 0.207, 0.203, 0.149, 0.0, 0.0, 0.379,
               0.372, 0.0],
        "Pc": [0.232, 0.224, 0.177, 0.186, 0.195, 0.143, 0.204, 0.360, 0.360,
               0.360, 0.360, 0.360, 0.720, 0.153, 0.173, 0.375, 0.370, 0.356,
               0.075, 0.0, 0.0, 0.126, 0.513, 0.504, 0.324, 0.712, 0.704,
               0.455, 0.367, 0.358, 0.311, 0.249, 0.0, 0.316, 0.269, 0.324,
               0.316, 0.269, 0.313, 0.304, 0.211, 0.0, 0.869, 0.869, 0.0,
               0.564, 0.511, 0.0, 0.542, 0.504, 0.461, 0.0, 0.822, 0.779, 0.0,
               1.161, 0.0, 0.460, 0.0, 0.617, 0.617, 0.617, 0.0, 0.0, 0.476,
               0.0, 0.0, 0.816, 0.522, 0.274, 0.274, 0.0, 0.0, 0.318, 0.340,
               0.0, 0.0, 0.886, 0.638, 0.391, 0.485, 0.398, 0.298, 0.251,
               0.269, 0.675, 0.675, 0.645, 0.200, 0.0, 1.027, 0.0, 0.709,
               0.956, 0.0, 0.0, 0.0, 0.0, 0.372, 0.605, 0.596, 0.946, 0.937,
               0.929, 0.658, 0.761, 0.485, 0.476, 0.429, 0.0, 0.0, 0.960,
               0.978, 0.0],
        "vc": [40, 41, 25, 30, 37, 5, 55, 32, 32, 32, 32, 32, 64, 16, 87, 68,
               95, 107, -25, 0.0, 0.0, -20, 77, 78, -8, 102, 103, -6, 41, 42,
               27, -57, 0.0, 78, 62, 77, 78, 62, 111, 112, 24, 0.0, 107, 107,
               0.0, 27, -31, 0.0, 79, 68, 43, 0.0, 107, 82, 0.0, 124, 0.0, 47,
               0.0, 34, 34, 34, 0.0, 0.0, 65, 0.0, 0.0, -7, 6, -12, -12, 0.0,
               0.0, 23, 27, 0.0, 0.0, 188, 127, 66, 47, -6, 41, 25, 37, 108,
               108, 108, -15, 0.0, 143, 0.0, 104, 165, 0.0, 0.0, 0.0, 0.0, 73,
               114, 115, 101, 102, 103, 55, 109, 64, 65, 49, 0.0, 0.0, 125,
               137, 0.0],
        "tb": [123, 121, 138,  97, 107,  74,  20, 257, 257, 257, 257, 257, 514,
               124, 247, 282, 303, 191, 474, 0.0, 0.0, 525, 514, 512, 396, 451,
               573, 426, 288, 286, 262, 323, 0.0, 437, 412, 444, 442, 418, 293,
               291, 655, 0.0, 942, 942, 0.0, 794, 858, 0.0, 360, 336, 313, 0.0,
               575, 552, 0.0, 598, 0.0, 358, 0.0, 692, 668, 818, 0.0, 0.0, 515,
               0.0, 0.0, 525, 353, 288, 288, 0.0, 0.0, 190, 135, 0.0, 0.0, 141,
               108, 91, 338, 164, 164, 164, 164, 44, 44, 61, 225, 0.0, 569,
               0.0, 477, 348, 0.0, 0.0, 0.0, 17, 707, 835, 833, 862, 860, 858,
               830, 495, 473, 471, 447, 0.0, 0.0, 0, 0, 0.0],
        "hf": [-45.947, -20.763, -20.763, -3.766, -3.766, 17.119, 17.119,
               53.712, 69.939, 64.145, 82.528, 104.293, 197.322, 11.189,
               27.016, -19.243, 9.404, 27.671, -181.422, 0.0, 0.0, -164.609,
               -182.329, -164.41, -129.158, -389.737, -359.258, -332.822,
               -163.569, -151.143, -129.488, -140.313, 0.0, -15.505, 3.32,
               5.432, 23.101, 26.718, 54.929, 69.885, 20.079, 0.0, 134.062,
               139.758, 0.0, 88.298, -396.242, 0.0, -73.568, -63.795, -57.795,
               0.0, -82.921, 0.0, 0.0, -107.188, 0.0, -16.752, 0.0, -66.138,
               -59.142, -7.365, 0.0, 0.0, -8.253, 0.0, 0.0, 57.546, 1.834,
               220.803, 227.368, 0.0, 0.0, -36.097, -161.74, 0.0, 0.0,
               -679.195, 0.0, 0.0, -313.545, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, -258.96, 0.0, 0.0, -446.835, 0.0, 0.0, 0.0, -223.398,
               -203.188, -67.778, -182.005, -189.888, -46.562, 0.0, -344.125,
               0.0, -2.084, 18.022, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "gf": [-8.03, 8.231, 8.231, 19.848, 19.848, 37.977, 37.977, 84.926,
               92.9, 88.402, 93.745, 116.613, 221.308, 22.533, 30.485, 22.505,
               41.228, 52.948, -158.589, 0.0, 0.0, -132.097, -131.366,
               -132.386, -107.858, -318.616, -291.188, -288.902, -105.767,
               -101.563, -92.099, -90.883, 0.0, 58.085, 63.051, 82.471, 95.888,
               85.001, 128.602, 132.756, 68.861, 0.0, 199.958, 199.288, 0.0,
               121.544, -349.439, 0.0, -33.373, -31.502, -25.261, 0.0, -35.814,
               0.0, 0.0, -53.332, 0.0, -0.50, 0.0, 17.963, 18.088, 60.161, 0.0,
               0.0, 16.731, 0.0, 0.0, 46.945, -1.721, 217.003, 216.328, 0.0,
               0.0, -28.148, -144.549, 0.0, 0.0, -626.58, 0.0, 0.0, -281.495,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -209.337, 0.0, 0.0,
               -392.975, 0.0, 0.0, 0.0, 212.718, 136.742, 0.0, 0.0, -65.642,
               0.0, 0.0, 241.373, 0.0, 30.222, 38.346, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0],
        "hv": [4.116, 4.65, 4.65, 2.771, 2.771, 1.284, 1.284, 6.714, 7.37,
               6.797, 8.178, 9.342, 12.318, 4.098, 12.552, 9.776, 10.185,
               8.834, 24.529, 0.0, 0.0, 40.246, 18.999, 20.041, 12.909, 22.709,
               17.759, 0.0, 10.919, 7.478, 5.708, 11.227, 0.0, 14.599, 11.876,
               14.452, 14.481, 0.0, 6.947, 6.918, 28.453, 0.0, 31.523, 31.005,
               0.0, 23.34, 43.046, 0.0, 13.78, 11.985, 9.818, 0.0, 19.208,
               17.574, 0.0, 0.0, 0.0, 11.883, 0.0, 30.644, 26.277, 0.0, 0.0,
               0.0, 14.931, 0.0, 0.0, 14.364, 11.423, 7.751, 11.549, 0.0, 0.0,
               4.877, 0.0, 0.0, 8.901, 1.86, 8.901, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 13.322, 0.0, 0.0, 8.301, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 51.787, 0.0, 0.0, 0.0, 0.0, 0.0, 16.921, 17.117,
               13.265, 0.0, 0.0, 27.966, 0.0, 0.0],
        "txt": [("CH3-",),
                ("CH2<",),
                ("RCH2<",),
                ("CH",),
                (">RCH-",),
                (">C<",),
                (">RC<",),
                ("CH2=CH",),
                ("CH=CH",),
                ("CH2=C",),
                ("CH=C",),
                ("C=C",),
                ("CH2=C=CH",),
                ("ACH",),
                ("AC-",),
                ("ACCH3",),
                ("ACCH2",),
                ("ACCH",),
                ("OH",),
                ("CH3OH",),
                ("H2O",),
                ("ACOH",),
                ("CH3CO",),
                ("CH2CO",),
                ("CHO",),
                ("CH3COO",),
                ("CH2COO",),
                ("HCOO",),
                ("CH3O",),
                ("CH2O",),
                ("CH-O",),
                ("FCH2O",),
                ("CH3NH2",),
                ("CH2NH2",),
                ("CHNH2",),
                ("CH3NH",),
                ("CH2NH",),
                ("CHNH",),
                ("CH3-RN",),
                ("CH2-RN",),
                ("ACNH2",),
                ("C5H5N",),
                ("C5H4N",),
                ("C5H3N",),
                ("CH3CN",),
                ("CH2CN",),
                ("COOH",),
                ("HCOOH",),
                ("CH2CL",),
                ("CHCL",),
                ("CCL",),
                ("CH2CL2",),
                ("CHCL2",),
                ("CCL2",),
                ("CHCL3",),
                ("CCL3",),
                ("CCL4",),
                ("ACCL",),
                ("CH3NO2",),
                ("CH2NO2",),
                ("CHNO2",),
                ("ACNO2",),
                ("CS2",),
                ("CH3SH",),
                ("CH2SH",),
                ("FURFURAL",),
                ("<CH2OH>2",),
                ("I",),
                ("Br",),
                ("CH===C",),
                ("C===C",),
                ("ME2SO",),
                ("ACRY",),
                ("CL<C=C>",),
                ("ACF",),
                ("DMF-1",),
                ("DMF-2",),
                ("CF3",),
                ("CF2",),
                ("CF",),
                ("COO",),
                ("SiH3",),
                ("SiH2",),
                ("SiH",),
                ("Si",),
                ("SiH2O",),
                ("SiHO",),
                ("SiO",),
                ("TERT-N",),
                ("CCL3F",),
                ("CCL2F",),
                ("HCCL2F",),
                ("HCCLF",),
                ("CCLF2",),
                ("HCCLF2",),
                ("CCLF3",),
                ("CCL2F2",),
                ("F (exceptions)",),
                ("CONH2",),
                ("CONHCH3",),
                ("CONHCH2",),
                ("CON<CH3>2",),
                ("CONCH3CH2",),
                ("CON<CH2>2",),
                ("C2H5O2",),
                ("C2H4O2",),
                ("CH3S",),
                ("CH2S",),
                ("CHS",),
                ("MORPH",),
                ("C4H4S",),
                ("C4H3S",),
                ("C4H2S",),
                ("NMP",)]}

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["M"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined molecular weight")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.M = self.kwargs["M"]
        tc = Pc = vc = tb = hv = gf = hf = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb += c*self._coeff["tb"][i]
            tc += c*self._coeff["tc"][i]
            Pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]
            hv += c*self._coeff["hv"][i]
            gf += c*self._coeff["gf"][i]
            hf += c*self._coeff["hf"][i]

        if self.kwargs["Tb"]:
            Tb = self.kwargs["Tb"]
        else:
            Tb = 1000/(0.5+35.7/tb**0.5+1000/(142+tb))
        self.Tb = unidades.Temperature(Tb)
        self.Tc = unidades.Temperature(self.Tb*(1+(1.28*tc)**-1))
        self.Pc = unidades.Pressure(self.M/(0.346+Pc)**2, "bar")
        self.Vc = unidades.SpecificVolume((172+vc)/self.M, "ccg")
        self.Hv = unidades.Enthalpy((hv+6.829)/self.M, "kJg")
        self.Hf = unidades.Enthalpy((hf+10.835)/self.M, "kJg")
        self.Gf = unidades.Enthalpy((gf-14.828)/self.M, "kJg")

        GroupContribution.calculo(self)


@refDoc(__doi__, [8, 9, 10, 11])
class Ambrose(GroupContribution):
    """
    Group contribution for definition of unknown component using the Ambrose
    procedure as use in API Technical Databook, procedure 4A1.1 with aditional
    term from Perry's Handbook

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    platt : int
        ΔPlatt number, [-]
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation
    The Platt number is the number of pairs of carbon atoms which are
    separated by three carbon-carbon bonds and is an indicator of the degree
    of branching in the molecule. The Platt number of an n-alkane is equal to
    the number of carbons minus three.

    Examples
    --------
    Example 1 in [10]_, 2,2,3-Trimethylpentane

    >>> Tb = unidades.Temperature(229.72, "F")
    >>> cmp = Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1],
    ... Tb=Tb, platt=3)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '555.83 400.74 0.0639'
    >>> cmp.formula
    'C8H18'

    Example 2 in [10]_, 2-Methyl-1-butene

    >>> Tb = unidades.Temperature(88.09, "F")
    >>> cmp = Ambrose(group=[0, 1, 4, 6], contribution=[2, 1, 1, 1],
    ... Tb=Tb, platt=0)
    >>> "%0.2f %0.1f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '385.95 520.3 0.0657'

    Example 3 in [10]_, cis-Decalin

    >>> Tb = unidades.Temperature(384.47, "F")
    >>> cmp = Ambrose(group=[10, 12], contribution=[8, 2], Tb=Tb, platt=0)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '801.95 430.06 0.0562'

    Example 4 in [10]_, tert-Butyl benzene

    >>> Tb = unidades.Temperature(336.41, "F")
    >>> cmp = Ambrose(group=[0, 3, 17], contribution=[3, 1, 1],
    ... Tb=Tb, platt=-1)
    >>> "%0.2f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '705.82 415.97 0.0555'

    Example 5 in [10]_, Anthracene

    >>> Tb = unidades.Temperature(646.16, "F")
    >>> cmp = Ambrose(group=[18, 28], contribution=[1, 2], M=178.23, Tb=Tb,
    ... platt=0)
    >>> "%0.1f %0.2f %0.4f" % (cmp.Tc.F, cmp.Pc.psi, cmp.Vc.ft3lb)
    '1165.3 504.64 0.0502'

    Example from [11]_, 2,2,4-trimethylpentane

    >>> cmp = Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1],
    ... Tb=372.39, platt=0)
    >>> "%0.1f %0.2f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '543.0 25.63 455.8'
    """
    __title__ = "Ambrose (1980)"
    kwargs = GroupContribution.kwargs.copy()
    kwargs["platt"] = 0

    _coeff = {
        "Pc": [0.2260, 0.2260, 0.22, 0.1960, 0.1935, 0.1935, 0.1875, 0.1610,
               0.1410, 0.1410, 0.1820, 0.1820, 0.1820, 0.1820, 0.1495, 0.1495,
               0.1170, 0.9240, 0.8940, 0.9440, 0.9440, 0.8640, 0.9140, 0.8340,
               0.8840, 0.8840, 0.8040, 0.7240, 0.5150, None, 0.160, 0.282,
               0.220, 0.450, 0.900, 0.470, 0.420, 0.095, 0.135, 0.170, 0.360,
               0.270, 0.270, 0.461, 0.507, 0, 0.725, 0.663, 0.223, 0.318,
               0.500, -0.025, 0.515, 0.183, 0.318, 0.600, 0.850, -0.065,
               -0.170, 0.924, 0.850, 0, 0.020, -0.050, 0],
        "tc": [0.138, 0.138, 0.095, 0.018, 0.113, 0.113, 0.070, 0.088, 0.038,
               0.038, 0.09, 0.09, 0.03, 0.09, 0.075, 0.075, 0.06, 0.458, 0.448,
               0.488, 0.488, 0.438, 0.478, 0.428, .468, .468, .418, .368, .22,
               None, 0.138, 0.220, 0.220, 0.578, 1.156, 0.330, 0.370, 0.208,
               0.208, 0.088, 0.423, 0.105, 0.090, 0.138, 0.371, 0.195, 0.159,
               0.131, 0.180, 0.110, 0.110, 0.198, 0.220, 0.080, 0.080, 0.080,
               0.010, -0.050, -0.200, 0.448, 0.448, 0.030, -0.040, -0.080],
        "vc": [55.1, 55.1, 47.1, 38.1, 45.1, 45.1, 37.1, 35.1, 35.1, 35.1,
               44.5, 44.5, 44.5, 44.5, 37, 37, 29.5, 222, 222, 222, 222, 222,
               222, 222, 222, 222, 222, 222, 148, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0],

        "txt": [("CH3", {"C": 1, "H": 3}),                                # 0
                ("CH2", {"C": 1, "H": 2}),
                ("CH", {"C": 1, "H": 1}),
                ("C", {"C": 1}),
                ("=CH2", {"C": 1, "H": 2}),
                ("=CH-", {"C": 1, "H": 1}),
                ("=C", {"C": 1}),
                ("=C=", {"C": 1}),
                ("≡CH", {"C": 1, "H": 1}),
                ("≡C-", {"C": 1}),
                ("-CH2- (Cyclic)", {"C": 1, "H": 2}),                     # 10
                ("-CH< (Ciclic)", {"C": 1, "H": 1}),
                ("-CH< (in fused ring)", {"C": 1, "H": 1}),
                (">C< (Ciclic)", {"C": 1}),
                ("=CH- (Cyclic)", {"C": 1, "H": 1}),
                ("=C< (Cyclic)", {"C": 1}),
                ("=C= (Cyclic)", {"C": 1}),
                ("Phenyl- ", {"C": 6, "H": 5}),
                ("o-Phenyl- ", {"C": 6, "H": 4}),
                ("m-Phenyl- ", {"C": 6, "H": 4}),
                ("p-Phenyl- ", {"C": 6, "H": 4}),                         # 20
                ("1,2,3-Phenyl- ", {"C": 6, "H": 3}),
                ("1,2,4-Phenyl- ", {"C": 6, "H": 3}),
                ("1,2,3,4-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,3,5-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,4,5-Phenyl- ", {"C": 6, "H": 2}),
                ("1,2,3,4,5-Phenyl- ", {"C": 6, "H": 1}),
                ("1,2,4,5,6-Phenyl- ", {"C": 6}),
                ("=CH-CH= (in fused Aromatic ring)", {"C": 2, "H": 2}),
                ("-OH", {"H": 1, "O": 1}),
                ("-O-", {"O": 1}),                                        # 30
                ("-CO-", {"C": 1, "O": 1}),
                ("-CHO", {"C": 1, "H": 1, "O": 1}),
                ("-COOH", {"C": 1, "H": 1, "O": 2}),
                ("-CO-O-OC-", {"C": 2, "O": 3}),
                ("-CO-O-", {"C": 1, "O": 2}),
                ("-NO2", {"N": 1, "O": 2}),
                ("-NH2", {"N": 1, "H": 2}),
                ("-NH-", {"N": 1, "H": 1}),
                ("-N<", {"N": 1}),
                ("-CN", {"N": 1, "C": 1}),                                # 40
                ("-S-", {"S": 1}),
                ("-SH", {"S": 1, "H": 1}),
                (">Si<", {"Si": 1}),
                (">SiH-", {"Si": 1, "H": 1}),
                ("-SiH3", {"Si": 1, "H": 3}),
                (">SiO-", {"Si": 1, "O": 1}),
                (">SiO- (cyclic)", {"Si": 1, "O": 1}),
                ("-F", {"F": 1}),
                ("-Cl", {"Cl": 1}),
                ("-Br", {"Br": 1}),                                       # 50
                ("-OH (Aromatic)", {"O": 1, "H": 1}),
                ("C4H4 (fused ring)", {"C": 4, "H": 4}),
                ("-F (Aromatic)", {"F": 1}),
                ("-Cl (Aromatic)", {"Cl": 1}),
                ("-Br (Aromatic)", {"Br": 1}),
                ("-I (Aromatic)", {"I": 1}),

                # 2nd Order terms
                ("Double Bond", ),
                ("Triple Bond", ),
                ("Benzene", ),
                ("Pyridine", ),                                           # 60
                ("The single or first substituent on an aromatic ring", ),
                ("The second or subsequent substituent on an aromatic ring", ),
                ("Each pair of rign substituents in ortho positions", ),
                ("If one of the ortho pair is -OH", )]}

    FirstOrder = 57
    SecondOrder = 65

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def _group(self):
        """From group contribution desglose the chemical composition"""
        group = []
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            # Only the first order term count for this
            if i < self.FirstOrder:
                grp = self._coeff["txt"][i][1]
                for x in range(c):
                    group.append(grp)
        return group

    def calculo(self):
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        Pc = tc = vc = 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i == 29:
                # Use special case for -OH group
                n = 0
                for g in self.group:
                    n += g.get("C", 0)

                tc += 0.87-0.11*n+0.003*n**2
                Pc += 0.1-0.013*n
                vc += 0
            else:
                tc += c*self._coeff["tc"][i]
                Pc += c*self._coeff["Pc"][i]
                vc += c*self._coeff["vc"][i]

        Pt = self.kwargs["platt"]
        self.Tc = unidades.Temperature(self.Tb*(1+1/(1.242+tc-0.023*Pt)))
        self.Pc = unidades.Pressure(14.5*self.M/(0.339+Pc-0.026*Pt)**2, "psi")
        self.Vc = unidades.SpecificVolume(0.01602*(40+vc)/self.M, "ft3lb")

        GroupContribution.calculo(self)


@refDoc(__doi__, [12])
class Klincewicz(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Klincewicz-Reid procedure (1984)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    nogroup : boolean
        Use the simple correlation without group contribution
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]
    atoms : int, optional
        Atoms count, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation. This
    method has a alternate not group contribution procedure known only Tb and M
    Na is only necessary if the nogroup opcion is enabled

    Examples
    --------
    Example in https://en.wikipedia.org/wiki/Klincewicz_method, acetone

    >>> cmp = Klincewicz(Tb=329.25, group=[0, 18], contribution=[2, 1])
    >>> "%0.2f %0.2f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '510.48 45.69 213.52'

    Same example without group contribution

    >>> cmp = Klincewicz(Tb=329.25, M=58.08, atoms=10, nogroup=True)
    >>> "%0.2f %0.4f %0.1f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '505.15 52.9098 205.2'
    """
    __title__ = "Klincewicz (1984)"
    kwargs = GroupContribution.kwargs.copy()
    kwargs["nogroup"] = False
    kwargs["atoms"] = 0

    _coeff = {
        # Table 6
        "tc": [-2.433, 0.353, 4.253, 6.266, -0.335, 16.416, 12.435, -0.991,
               3.786, 3.373, 7.169, 7.169, 5.623, -4.561, 7.341, -28.930,
               5.389, 7.127, 4.332, 4.332, -25.085, 8.890, -4.153, 2.005,
               2.773, 12.253, 8.239, -10.381, 28.529, 23.905, 31.537, 5.191,
               18.353, 53.456, 94.186, -1.770, -1.770, -1.770, -1.770, 11.709],
        "Pc": [0.026, -0.015, -0.046, -0.083, -0.027, -0.136, -0.111, -0.015,
               -0.050, -0.066, -0.067, -0.067, -0.089, -0.056, -0.112, -0.190,
               -0.143, -0.116, -0.196, -0.196, -0.251, -0.277, -0.127, -0.180,
               -0.172, -0.163, -0.104, -0.064, -0.303, -0.311, -0.208, -0.067,
               -0.244, -0.692, -1.051, 0.032, 0.032, 0.032, 0.032, -0.325],
        "vc": [16.2, 16.1, 8.2, 12.1, 7.4, 8.95, -6.6, 13.9, 9.8, 5.1, 2.7,
               2.7, 0.2, 7.5, 3.0, -24.0, -26.1, -36.6, -6.7, -6.7, -37.0,
               -28.2, -0.1, 53.7, -8.0, -0.7, -18.4, 12.0, -27.7, -27.3, -61.9,
               -34.1, -47.4, -148.1, -270.6, 0.8, 0.8, 0.8, 0.8, -39.2],

        "txt": [("-CH3", ),                     # 0
                ("-CH2-", ),
                ("-CH2- (ring)", ),
                (">CH-", ),
                (">CH- (ring)", ),
                (">C<", ),
                (">C< (ring)", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=CH- (ring)", ),
                (">c=", ),                      # 10
                ("=c=", ),
                (">C= (ring)", ),
                ("=CH", ),
                ("s-", ),
                ("-OH", ),
                ("-0-", ),
                ("-0- (ring)", ),
                (">CO", ),
                ("-CHO", ),
                ("-COOH", ),                    # 20
                ("-CO-0-", ),
                ("-NH2", ),
                (">NH", ),
                (">NH (ring)", ),
                (">N-", ),
                ("=N- (ring)", ),
                ("-CN", ),
                ("-SH", ),
                ("-S-", ),
                ("-S- (ring)", ),               # 30
                ("-F", ),
                ("-C1", ),
                ("-Br", ),
                ("-I", ),
                ("-CF2", ),
                ("-CCl2", ),
                ("-CBr2", ),
                ("-CI2", ),
                ("-NO2", )]}

    FirstOrder = 40

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined boiling point")
            self.status = 0
        elif self.kwargs["nogroup"]:
            self.group = []
            if not self.kwargs["M"]:
                self.msg = QApplication.translate(
                        "pychemqt", "undefined molecular weight")
                self.status = 0
            elif not self.kwargs["atoms"]:
                self.msg = QApplication.translate(
                        "pychemqt", "undefined atoms number of molecule")
                self.status = 0
            else:
                self.status = 1
                self.msg = ""
                return True
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        if self.kwargs["nogroup"]:
            # No group calculation
            self.Na = self.kwargs["atoms"]

            # Eq 13-15
            self.Tc = unidades.Temperature(50.2-0.16*self.M+1.41*self.Tb)
            self.Pc = unidades.Pressure(
                    self.M/(0.335+0.009*self.M+0.019*self.Na)**2, "bar")
            self.Vc = unidades.SpecificVolume(
                    (20.1+0.88*self.M+13.4*self.Na)/self.M, "ccg")

        else:
            Pc = tc = vc = 0
            for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
                tc += c*self._coeff["tc"][i]
                Pc += c*self._coeff["Pc"][i]
                vc += c*self._coeff["vc"][i]

            # Eq 10-12
            self.Tc = unidades.Temperature(45.4-0.77*self.M+1.55*self.Tb+tc)
            self.Pc = unidades.Pressure(
                    self.M/(0.348+0.0159*self.M+Pc)**2, "bar")
            self.Vc = unidades.SpecificVolume(
                    (25.2+2.8*self.M+vc)/self.M, "ccg")

        GroupContribution.calculo(self)


@refDoc(__doi__, [11, 13])
class Lydersen(GroupContribution):
    """
    Group contribution for definition of unknown component using the Lydersen
    procedure (1955)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 2 in [11]_ pag 2-343, 2-butanol critical properties

    >>> cmp = Lydersen(Tb=372.7, group=[0, 1, 2, 18],
    ... contribution=[2, 1, 1, 1])
    >>> "%0.1f %0.3f" % (cmp.Tc, cmp.Pc.MPa)
    '534.5 4.506'
    >>> cmp.formula
    'C4H10O'

    Example in http://en.wikipedia.org/wiki/Lydersen_method, acetone

    >>> cmp = Lydersen(Tb=329.25, group=[0, 22], contribution=[2, 1])
    >>> "%0.0f" % (cmp.Vc.ccg*cmp.M)
    '210'
    """
    __title__ = "Lydersen (1955)"
    _coeff = {
        # Table III
        "tc": [0.020, 0.020, 0.012, 0.00, 0.018, 0.018, 0.00, 0.00, 0.005,
               0.005, 0.013, 0.012, -0.007, 0.011, 0.011, 0.011, 0.066, 0.066,
               0.082, 0.031, 0.021, 0.014, 0.040, 0.033, 0.048, 0.085, 0.047,
               0.02, 0.018, 0.017, 0.010, 0.012, 0.031, 0.031, 0.024, 0.014,
               0.007, 0.060, 0.055, 0.015, 0.015, 0.008, 0.003, 0.026, 0.040,
               0.027, 0.025, 0.027],
        "Pc": [0.227, 0.227, 0.210, 0.210, 0.198, 0.198, 0.198, 0.198, 0.153,
               0.153, 0.184, 0.192, 0.154, 0.154, 0.154, 0.154, 0.924, 0.924,
               0.06, -0.02, 0.16, 0.12, 0.29, 0.2, 0.33, 0.4, 0.47, 0.12,
               0.224, 0.320, 0.50, 0.83, 0.095, 0.135, 0.09, 0.17, 0.13, 0.36,
               0.42, 0.27, 0.27, 0.24, 0.24, 0.468, 0.513, 0, 0.730, 0.668],
        "vc": [0.055, 0.055, 0.051, 0.041, 0.045, 0.045, 0.036, 0.036, 0.036,
               0.036, 0.0445, 0.046, 0.031, 0.036, 0.036, 0.037, 0, 0, 0.018,
               0.003, 0.020, 0.008, 0.060, 0.050, 0.073, 0.080, 0.080, 0.011,
               0.018, 0.049, 0.070, 0.095, 0.028, 0.037, 0.027, 0.042, 0.032,
               0.080, 0.078, 0.055, 0.055, 0.045, 0.047, 0, 0, 0, 0, 0],

        # Name and element composition
        "txt": [("CH3-", ),                     # 0
                ("-CH2-", ),
                ("-CH<", ),
                (">C<", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=C<", ),
                ("=C=", ),
                ("≡CH", ),
                ("≡C-", ),
                ("-CH2- (cyclic)", ),           # 10
                ("-CH< (cyclic)", ),
                (">C< (cyclic)", ),
                ("=C< (cyclic)", ),
                ("=C= (cyclic)", ),
                ("=CH- (cyclic)", ),
                ("-CH= (Aromatic)", ),
                ("=C< (Aromatic)", ),
                ("-OH", ),
                ("-OH (Aromatic)", ),
                ("-O-", ),                      # 20
                ("-O- (cyclic)", ),
                (">C=O", ),
                (">C=O (cyclic)", ),
                ("-CH=O", ),
                ("-COOH", ),
                ("-COO-", ),
                ("=O (other)", ),
                ("F", ),
                ("Cl", ),
                ("Br", ),                       # 30
                ("I", ),
                ("-NH2", ),
                (">NH", ),
                (">NH (cyclic)", ),
                (">N-", ),
                (">N- (cyclic)", ),
                ("-CN", ),
                ("NO2", ),
                ("-SH", ),
                ("-S-", ),                      # 40
                ("-S- (cyclic)", ),
                ("=S", ),
                (">Si<", ),
                ("-SiH<", ),
                ("-SiH3", ),
                (">SiO-", ),
                (">SiO- (cyclic)", ),

                ]}

    FirstOrder = 47

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]

        self.Tc = unidades.Temperature(self.Tb/(0.567+tc-tc**2))
        self.Pc = unidades.Pressure(self.M/(0.34+pc)**2, "atm")
        self.Vc = unidades.SpecificVolume((0.04+vc)/self.M)

        GroupContribution.calculo(self)


@refDoc(__doi__, [14])
class Valderrama(GroupContribution):
    """
    Group contribution for definition of unknown component using the Valderrama
    procedure (2006)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float, optional
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example from in [14]_ Table 8, phenantrene

    >>> cmp = Valderrama(Tb=372.7, group=[30, 32], contribution=[10, 4])
    >>> "%0.0f %0.1f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '871 31.6 550.55'

    2-methyl-1-pentanol

    >>> cmp = Valderrama(Tb=401.85,
    ... group=[0, 1, 2, 10], contribution=[2, 3, 1, 1])
    >>> "%0.0f %0.1f %0.2f" % (cmp.Tc, cmp.Pc.bar, cmp.Vc.ccg*cmp.M)
    '602 33.5 381.53'
    """
    __title__ = "Valderrama (2006)"
    _coeff = {
        # Table III
        "tc": [8.26, 20.07, 27.11, 35.81, 3.27, 18.81, 26.64, 15.40, 6.46,
               6.46, 11.21, 9.93, 34.18, 21.98, 50.52, 40.07, 28.83, -7.64,
               13.70, 15.00, 24.66, 26.42, 49.30, 48.83, -7.25, 8.28, 16.07,
               41.78, 15.09, 21.92, 11.23, 31.39, 27.40, 2.13, 1.17, 27.48,
               4.28, 33.80, 9.23],
        "Pc": [0.8446, 1.2910, 1.4806, 1.6350, 0.4418, 1.1606, 1.4472, 0.6119,
               0.2734, 0.2734, 0.0881, 0.3628, 1.5195, 0.7953, 1.9871, 1.9804,
               1.3779, -0.9769, 0.4922, 0.2946, 0.9936, 1.8782, 2.5229, 2.0970,
               0.0210, 0.5375, 0.5006, 1.7649, 0.9240, 1.0535, 0.6178, 1.3134,
               1.2796, -0.0201, -0.4922, 1.0798, -0.4429, 1.2695, -0.0791],
        "vc": [45.42, 37.60, 30.00, 13.28, 41.36, 33.48, 21.93, 21.56, 28.60,
               28.60, 20.75, 10.81, 47.41, 52.80, 60.10, 56.69, 65.87, 26.61,
               32.98, 51.67, 17.34, 32.62, 58.84, 83.38, 21.46, 41.72, 51.66,
               68.45, 34.69, 20.42, 28.49, 11.00, 20.63, 10.59, -11.19, 40.67,
               19.04, 16.14, 29.24],

        # Name and element composition
        "txt": [("CH3-", ),                     # 0
                ("-CH2-", ),
                (">CH-", ),
                (">C<", ),
                ("=CH2", ),
                ("=CH-", ),
                ("=C<", ),
                ("=C=", ),
                ("=(-)CH", ),
                ("=(-)C-", ),
                ("-OH", ),                      # 10
                ("-O-", ),
                (">C=O", ),
                ("-CHO", ),
                ("-COOH", ),
                ("-COO-", ),
                ("HCOO-", ),
                ("=O (Any_other)", ),
                ("-NH2", ),
                ("-NH-", ),
                (">N-", ),                      # 20
                ("=N-", ),
                ("-CN", ),
                ("-NO2", ),
                ("-F", ),
                ("-Cl", ),
                ("-Br", ),
                ("-I", ),
                ("-CH2- (ring)", ),
                (">CH- (ring)", ),
                ("=CH- (ring)", ),              # 30
                (">C< (ring)", ),
                ("=C< (ring)", ),
                ("-O- (ring)", ),
                ("-OH(Phenols) (ring)", ),
                (">C=O (ring)", ),
                ("-NH- (ring)", ),
                (">N- (ring)", ),
                ("=N- (ring)", )
                ]}

    FirstOrder = 38

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]

        g1 = 38.91+tc**0.88
        g2 = 5.84+pc**1.27
        self.Pc = unidades.Pressure((g1/g2)**2, "bar")
        self.Tc = unidades.Temperature(self.Pc.bar*g2)
        self.Vc = unidades.SpecificVolume((26.86+vc**1.06)/self.M, "ccg")

        GroupContribution.calculo(self)


@refDoc(__doi__, [15, 16, 17, 18])
class Nannoolal(GroupContribution):
    """
    Group contribution for definition of unknown component using the Nannoolal
    procedure (2007)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float, optional
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Table 17a in [15]_, 3,3,4,4-tetramethylhexane

    >>> cmp = Nannoolal(group=[0, 3, 5, 131], contribution=[6, 2, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '429.5'
    >>> cmp.formula
    'C10H22'

    Table 17b in [15]_, di-isopropanolamine, find too in [10]_

    >>> cmp = Nannoolal(group=[0, 6, 33, 41], contribution=[2, 4, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '509.3'

    Table 17c in [15]_, perfluoro-2-propanone

    >>> cmp = Nannoolal(group=[6, 20, 50, 118, 119, 121],
    ... contribution=[2, 6, 1, 1, 2, 1])
    >>> "%0.1f" % cmp.Tb
    '246.3'

    Table 17d in [15]_, methyl m-toluate

    >>> cmp = Nannoolal(group=[1, 2, 14, 15, 44, 126, 132],
    ... contribution=[1, 1, 4, 2, 1, 1, 1])
    >>> "%0.1f" % cmp.Tb
    '490.7'

    Table 42a in [16]_, 2,2,3,3-tetramethylbutane

    >>> cmp = Nannoolal(Tb=379.6, group=[0, 5, 131], contribution=[6, 2, 1])
    >>> "%0.1f" % cmp.Tc
    '566.8'

    Table 42b in [16]_, diethylene glycol monomethyl ether

    >>> cmp = Nannoolal(M=120.15, group=[1, 6, 34, 37],
    ... contribution=[1, 4, 1, 2])
    >>> "%0.1f" % cmp.Pc.kPa
    '3689.9'

    Table 42c in [16]_, trichloro silane

    >>> cmp = Nannoolal(group=[26, 111, 122, 133], contribution=[3, 1, 1, 1])
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '262.8'

    Table 42d in [16]_, perfluoro-2-propanone

    >>> cmp = Nannoolal(Tb=245.9, group=[6, 20, 50, 118, 119, 121],
    ... contribution=[2, 6, 1, 1, 2, 1])
    >>> "%0.1f" % cmp.Tc
    '358.0'

    Table 15c in [17]_

    >>> "%0.7f" % cmp.db
    '0.1179392'
    >>> "%0.2f" % cmp._Pv(210.16).kPa
    '14.63'

    Table 15a in [17]_, α-pinene

    >>> cmp = Nannoolal(Tb=429, group=[0, 8, 9, 10, 61, 123, 130],
    ... contribution=[3, 2, 2, 1, 1, 1, 2])
    >>> "%0.7f" % cmp.db
    '0.1078596'
    >>> "%0.2f" % cmp._Pv(388.15).kPa
    '31.03'

    Table 15b in [17]_, 1,2-ethanediol

    >>> cmp = Nannoolal(Tb=470.5, group=[6, 35], contribution=[2, 2])
    >>> "%0.7f" % cmp.db
    '1.1310491'
    >>> "%0.2f" % cmp._Pv(410.65).kPa
    '13.05'

    Table 15d in [17]_, acrylic acid

    >>> cmp = Nannoolal(Tb=413.6, group=[43, 60, 132], contribution=[1, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.9163297'
    >>> "%0.2f" % cmp._Pv(344.15).kPa
    '6.52'

    Table 15e in [17]_, glycol monoacetate

    >>> cmp = Nannoolal(Tb=458.65, group=[0, 6, 35, 44],
    ... contribution=[1, 2, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.5111422'
    >>> "%0.2f" % cmp._Pv(352.65).kPa
    '2.24'

    Table 15f in [17]_, dipropyl succinate

    >>> cmp = Nannoolal(group=[0, 3, 6, 44], contribution=[2, 4, 2, 2])
    >>> "%0.7f" % cmp.db
    '0.9878298'

    Table 15g in [17]_, diethanolamine
    The example in paper has a error, the group interaction contribution has
    the sign changed, fixing this the value get is near to experimental value

    >>> cmp = Nannoolal(Tb=541.15, group=[6, 35, 41], contribution=[4, 2, 1])
    >>> "%0.5f" % cmp._Pv(401.13).kPa
    '0.40405'

    Table 15h in [17]_, R122

    >>> cmp = Nannoolal(Tb=344.25, group=[6, 20, 25, 26, 119, 122],
    ... contribution=[2, 2, 2, 1, 1, 1])
    >>> "%0.7f" % cmp.db
    '0.0447374'
    >>> "%0.4f" % cmp._Pv(297.46).kPa
    '17.5093'

    Table 28a in [18]_, N,N-diethylamine

    >>> cmp = Nannoolal(Tb=329, group=[0, 6, 41], contribution=[2, 2, 1])
    >>> "%0.7f %0.2f" % (cmp.dBv, cmp.Tv)
    '4.7713112 210.38'
    >>> "%0.4f" % cmp._Visco(308.15).mPas
    '0.2674'

    Table 28b in [18]_, ethylene glycol monopropyl

    >>> cmp = Nannoolal(Tb=424.5, group=[0, 3, 6, 34, 37],
    ... contribution=[1, 1, 3, 1, 1])
    >>> "%0.7f %0.3f" % (cmp.dBv, cmp.Tv)
    '6.0699241 301.197'
    >>> "%0.6f" % cmp._Visco(318.15).mPas
    '0.921814'

    Table 28c in [18]_, monoethanolamine

    >>> cmp = Nannoolal(Tb=443.45, group=[6, 35, 39], contribution=[2, 1, 1])
    >>> "%0.7f %0.3f" % (cmp.dBv, cmp.Tv)
    '12.3872997 382.275'
    >>> "%0.4f" % cmp._Visco(363.15).mPas
    '2.6135'
    """
    __title__ = "Nannoolal (2007)"
    _coeff = {
        # Be careful, there are several changes in group between Tb paper and
        # the other paper and several index ajust
        #   - 78 : Two different group with equal parameter but different M,
        #     the Si related copy to 110
        #   - 71, 93: Change definition, use the critical paper definition
        #   - Mising 98, use the 214 defined in critical (be careful to delete
        #     in table
        #   - Mising 112, use the 215 defined in critical
        #   - Mising 114, use the 216 defined in critical
        #   - 217 set at end of second order as 135

        # Table 3 & 4 in [15]_
        "tb": [177.3066, 251.8338, 157.9527, 239.4531, 240.6785, 249.5809,
               266.8769, 201.0115, 239.4957, 222.1163, 209.9749, 250.9584,
               291.2291, 244.3581, 235.3462, 315.4128, 348.2779, 367.9649,
               106.5492, 49.2701, 53.1871, 78.7578, 103.5672, -19.5575,
               330.9117, 287.1863, 267.4170, 205.7363, 292.5816, 419.4959,
               377.6775, 556.3944, 349.9409, 390.2446, 443.8712, 488.0819,
               361.4775, 146.4836, 820.7118, 321.1759, 441.4388, 223.0992,
               126.2952, 1080.3139, 636.2020, 642.0427, 1142.6119, 1052.6072,
               1364.5333, 1487.4109, 618.9782, 553.8090, 434.0811, 461.5784,
               864.5074, 304.3321, 719.2462, 475.7958, 586.1413, 500.2434,
               412.6276, 475.9623, 512.2893, 422.2307, 37.1936, 453.3397,
               306.7139, 866.5843, 821.4141, 282.0181, 207.9312, 920.3617,
               1153.1344, 494.2668, 1041.0851, 1251.2675, 778.9151, 540.0895,
               879.7062, 660.4645, 1018.4865, 1559.9840, 510.4223, 1149.9670,
               1209.2972, 347.7717, 664.0903, 957.6388, 928.9954, 560.1024,
               229.2288, 606.1797, 215.3416, 273.1755, 1218.1878, 2082.3288,
               201.3224, 0, 886.7613, 1045.0343, -109.6269, 111.0590,
               1573.3769, 1483.1289, 1506.8136, 484.6371, 1379.4485, 659.7336,
               492.0707, 540.0895, 971.0365, 0, 428.8911, 0, 612.9506,
               562.1791, 761.6006,
               -82.2328, -247.8893, -20.3996, 15.4720,
               -172.4201, -99.8035, -62.3740, -40.0058, -27.2705, -3.5075,
               16.1061, 25.8348, 35.8330, 51.9098, 111.8372, 40.4205, 0],

        # Table 7 & 8 in [16]_
        "tc": [41.8682, 33.1371, -1.0710, 40.0977, 30.2069, -3.8778, 52.8003,
               9.4422, 21.2898, 26.3513, -17.0459, 51.7974, 18.9549, -29.1568,
               16.1154, 68.2045, 68.1923, 29.8039, 15.6068, 11.0757, 18.1302,
               19.1772, 20.8519, -24.0220, -1.3329, 2.6113, 15.5010, -16.1905,
               60.1907, 5.2621, -21.5199, -8.6881, 84.8567, 79.3047, 49.5968,
               130.1320, 14.0159, 12.5082, 41.3490, 18.3404, -50.6419, 17.1780,
               -0.5820, 199.9042, 75.7089, 58.0782, 109.1930, 102.1024, 0, 0,
               56.1572, 44.2000, -7.1070, 0.5887, 0, -7.7181, 117.1330,
               45.1531, 0, 67.9821, 45.4406, 56.4059, -19.9737, 36.0883,
               10.4146, 18.9903, 10.9495, 82.6239, 0, 25.4209, 72.5587, 0, 0,
               0, 0, 164.3355, 0, 157.3401, 97.2830, 153.7225, 0, 90.9726,
               62.3642, 0, 0, 0, 53.6350, 24.7302, 0, 38.4681, 0, 63.6504,
               34.2058, 0, 0, 0, 27.3441, 48.1680, 0, 0, 0, 1.3231, 764.9595,
               0, 0, 0, 0, 0, 0, 157.3401, 0, 0.2842, 0, -0.6536, 36.0361, 0,
               0, 32.1829, 11.4437, -1.3023, -34.3037, -1.3798, -2.7180,
               11.3251, -4.7516, 1.2823, 6.7099, 0, -33.8201, -18.4815,
               -23.6024, -24.5802, -35.6113, 62.0286],

        # Table 10 & 11 in [16]_
        "Pc": [8.1620, 5.5262, 4.1660, 5.2623, 2.3009, -2.9925, 3.4310, 2.3665,
               3.4027, 3.6162, -5.1299, 4.1421, 0.8765, -0.1320, 2.1064,
               4.1826, 3.5500, 1.0997, 0.7328, 4.3757, 3.4933, 2.6558, 1.6547,
               0.5236, -2.2611, -1.4992, 0.4883, -0.9280, 11.8687, -4.3170,
               -2.2409, -4.7841, -7.4244, -4.4735, -1.8153, -6.8991, -12.1664,
               2.0592, 0.1759, -4.4164, -9.0065, -0.4086, 2.3625, 3.9873,
               4.3592, 1.0266, 0.4329, 0.5172, 0, 0, 0.1190, -2.3615, -9.4154,
               -8.2595, 0, -4.9259, 5.1666, 7.1581, 0, -6.2791, 9.6413, 3.4731,
               -2.2718, 2.4489, -0.5403, 8.3052, -4.7101, -5.0929, 0, 5.7270,
               2.7602, 0, 0, 0, 0, 4.0458, 0, 12.6786, 0.2822, 0, 0, -23.9221,
               0.7043, 0, 0, 0, 12.6128, -10.2451, 0, -4.0133, 0, -5.0403,
               3.2023, 0, 0, 0, -4.3834, 1.5574, 0, 0, 0, 3.3971, 58.9190, 0,
               0, 0, 0, 0, 0, 12.6786, 0, 3.8751, 0, 4.4882, -5.1116, 0, 0,
               7.3149, 4.1439, 0.4387, -4.2678, 4.8944, 2.8103, -0.3035,
               0.0930, 0.7061, -0.7246, 0, -8.8457, -2.2542, -3.2460, -5.3113,
               1.0934, 8.6126],

        # Table 13 & 14 in [16]_
        "vc": [28.7855, 28.8811, 26.7237, 32.0493, 32.1108, 28.0534, 33.7577,
               28.8792, 24.8517, 30.9323, 5.9550, 29.5901, 20.2325, 10.5669,
               19.4020, 25.0434, 5.6704, 16.4118, -5.0331, 1.5646, 3.3646,
               1.0897, 1.1084, 19.3190, 22.0457, 23.9279, 26.2582, 36.7624,
               34.4110, 36.0223, 30.7004, 48.2989, 10.6790, 5.6645, 2.0869,
               3.7778, 25.6584, 11.6284, 46.7680, 13.2571, 73.7444, 20.5722,
               6.0178, 40.3909, 42.6733, 36.1286, 0, 64.3506, 0, 0, 30.9229,
               25.5034, 34.7699, 38.0185, 0, 20.3127, 43.7983, 0, 0, 51.0710,
               48.1957, 34.1240, 40.9263, 29.8612, 4.7476, -25.3680, 23.6094,
               34.8472, 0, 75.7193, 69.5645, 0, 0, 0, 0, 0, 0, 0, 52.8789,
               27.1026, 0, 68.0701, 0, 0, 0, 0, 0, 64.4616, 0, 20.0440, 0,
               28.7127, 55.3822, 0, 0, 0, 29.3068, 16.3122, 0, 0, 0, 1.3597, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 37.0423, 0, 55.7432, 16.2688, 0, 0,
               -3.8033, 27.5326, 1.5807, -2.6235, -5.3091, -6.1909, 3.2219,
               -6.3900, -3.5964, 1.5196, 0, -4.6483, -5.0563, -6.3267, 4.9392,
               2.8889, 19.4348],

        # Table 4 & 5 in [17]_
        "Pv": [13.3063, 91.8000, 50.1939, 54.6564, 45.7437, -31.7531, 37.8485,
               96.1386, 22.2573, 32.8162, 4.8500, 23.6411, 49.8237, -3.6950,
               32.7177, 69.8796, 41.5534, 43.7191, 79.5429, 51.2880, 42.0887,
               56.9998, 142.1060, 45.9652, 93.6679, 67.8082, 55.9304, 46.0435,
               84.9162, 104.9291, -40.1837, 134.3501, 719.3666, 758.4218,
               700.7226, 756.0824, 441.8437, 108.4964, 286.9731, 251.9212,
               361.7760, 193.7667, -102.7252, 1074.1000, 355.7381, 350.5184,
               292.8046, 269.2471, 736.9540, 1216.0700, 255.8480, 252.9059,
               123.2143, 127.3380, 222.2789, 20.1604, 226.1873, 86.4601,
               224.2238, 134.9382, 34.2541, 97.4210, 206.6665, 128.0247,
               48.8839, 375.0486, 126.3340, 375.8217, 238.2066, 2.8992, 9.3624,
               603.5347, 662.0582, 510.9666, 1317.4360, 681.3525, 564.1116,
               391.3697, 318.2350, 435.8446, 218.6012, 381.5945, 80.2735,
               231.3919, 186.9204, 48.5026, 168.3007, 108.5277, 213.7165,
               183.1130, 1178.1950, 158.3258, -47.3420, 186.7950, 392.2006,
               595.1778, 268.7081, 86.9450, 612.9546, 258.9924, -316.4392,
               64.0566, 660.2247, 554.7941, 420.7591, 0, -237.2124, 0, 0,
               391.3697, 319.4879, -66.5670, 118.8412, -81.1543, 305.1341,
               191.5058, 423.5251, 34.3545, 2.5030, -83.3326, -64.4854,
               -125.9208, -47.2962, 33.9765, -7.0982, -45.0531, -3.2036, 0,
               -20.6706, -36.3170, -1.1994, 123.7433, -15.9694, 36.7574],

        # Table 5 & 6 in [18]_
        "dBv": [13.9133, 11.7002, -11.0660, 2.1727, 4.5878, 37.0296, 21.3473,
                5.9452, 10.8799, -7.2202, 142.1976, 61.0811, 28.7351, -12.3456,
                -2.7840, 45.9403, 80.5124, 37.4124, 5.1640, 0, 2.8323, 0.7129,
                -36.3189, -61.9434, 4.7579, 5.8228, 4.6555, -67.3989, -9.6209,
                0.5164, -18.1984, -17.3110, 336.8834, 365.8067, 249.0118,
                218.8000, 160.8315, -35.3055, 85.3693, 58.9131, 44.0698,
                13.6479, -58.7354, 54.7891, 17.6757, 0.4267, 47.6109, 12.6717,
                129.8293, 202.2864, 24.2524, 18.4961, 30.5022, -0.0276,
                -13.4614, 18.5507, 23.1459, 9.5809, 152.2693, 18.5983, 21.4560,
                19.7836, -165.0071, 13.3585, 42.7958, 151.9493, 52.5900,
                -34.3948, -6.5626, -25.5950, -28.3943, -30.6156, 45.9972,
                -7.3298, 369.8367, 16.3525, -2.6553, -61.2368, -7.5067, 4.1408,
                -46.5613, 122.6902, 0, 0, 0, -70.9713, 0, 29.0985, 125.0861,
                -14.2823, 0, -8.5352, 9.9037, 0, 0, 102.0816, 74.0520, 37.5669,
                0, 54.4769, 0, 5.7765, 95.6531, 56.9133, 64.7133, 0, 22.9969,
                23.2473, 0, -61.2368, 0, 64.6600, 0, 68.4952, 0, 0, 0, 0.3041,
                0, -6.1420, -26.4635, -14.9636, -25.9017, -57.3789, -21.2204,
                -20.1917, -34.5860, 0, -110.7391, 2.4859, -59.3670, 0, 13.1413,
                -76.1631],

        # Table 8 & 9 in [18]_
        "Tv": [89.0803, 216.0226, 80.9698, 60.3316, 24.2637, 244.4643,
               103.4109, -16.5212, 174.1316, -37.7584, 252.0190, 251.9299,
               330.7100, 294.3323, 113.9028, -26.6195, 133.5499, 128.4739,
               208.3258, 0, 35.2688, 207.3562, -15.8544, 112.1172, 329.0113,
               313.1106, 194.6060, -8.6247, 182.7067, 456.3713, 391.6060,
               499.2149, 1199.4010, 1198.1040, 1078.0840, 1284.7450, 1134.1640,
               -34.9892, 612.7222, 458.7425, 705.1250, 159.5146, -284.4707,
               1446.0240, 325.5736, 454.1671, 374.6477, 289.9690, 1150.8290,
               1619.1650, 304.5982, 394.7932, 294.7319, 206.6432, 292.3613,
               302.2321, 346.9998, -23.9801, 238.3242, 137.5408, 74.4489,
               304.9257, -32.4179, 57.8131, 279.2114, 662.0051, 277.5038,
               369.4221, 488.1136, -10.6146, -181.7627, 351.0623, -15.2801,
               174.3672, 1098.1570, 549.1481, 394.5776, 10.3752, 365.8081,
               164.8904, 197.1806, 1297.7560, 0, 0, 0, 35.4672, 0, 495.5141,
               551.9254, 490.7224, 0, 669.0158, 256.5078, 0, 0, 1787.0390,
               220.0803, 192.1303, 0, 131.2253, 0, 53.2507, 288.4140, 542.6641,
               714.0494, 0, 797.2271, 253.5303, 0, 10.3752, 0, 377.7146, 0,
               806.8125, 0, 0, 0, -180.3686, 0, 241.8968, 138.6555, -71.1647,
               -115.0418, -96.7544, -153.8442, -22.1041, 24.7835, 0, 224.2439,
               24.2539, 137.8708, 0, -54.1782, -726.4291],

        # Name and group composition, Table 1 in [15]_
        "txt": [                                           # 0
                ("CH3-", ),
                ("CH3- (N, O, F, Cl)", ),
                ("CH3- (aromatic)", ),
                ("-CH2- ", ),
                (">CH- ", ),
                (">C< ", ),
                (">C< (N, O, F, Cl)", ),
                (">C< (aromatic)", ),
                ("-CH2- (ring)", ),
                (">CH- (ring)", ),                         # 10
                (">C< (ring) ", ),
                (">C< (ring) (outer N, O, Cl, F)", ),
                (">C< (ring N, O)", ),
                (">C< (aromatic)", ),
                ("=CH- (aromatic) ", ),
                ("=C< (aromatic)", ),
                ("=C< (aromatic) (N, O, Cl, F)", ),
                ("=C< (full aromatic)", ),
                ("F– (C, Si)", ),
                ("-CF=C< ", ),                             # 20
                ("F– (CF)", ),
                ("F– (CF, CCl)", ),
                ("F– (CF2, CCl2)", ),
                ("F- (aromatic) ", ),
                ("Cl– (C, Si)", ),
                ("Cl– (CF, CCl)", ),
                ("Cl– (CF2, CCl2)", ),
                ("Cl– (aromatic)", ),
                ("-CCl=C< ", ),
                ("Br-", ),                                 # 30
                ("Br- (aromatic)", ),
                ("I-", ),
                ("-OH (>C< no H)", ),
                ("-OH (-CH2-)", ),
                ("-OH (>C5 chain)", ),
                ("-OH (<C5 chain)", ),
                ("-OH (aromatic)", ),
                ("-O-", ),
                (">OC2<", ),
                ("NH2-", ),                                # 40
                ("NH2- (aromatic)", ),
                ("-NH-", ),
                (">N-", ),
                ("-COOH", ),
                ("-COO-", ),
                ("HCOO-", ),
                ("-COO– (ring)", ),
                ("-CON<", ),
                ("-CONH-", ),
                ("-CONH2", ),                              # 50
                (">C=O", ),
                ("-CHO", ),
                ("-SH", ),
                ("-S-", ),
                ("-S-S-", ),
                ("-S- (aromatic) ", ),
                ("-C≡N", ),
                (">C=C< ", ),
                (">C=C< (aromatic)", ),
                (">C=C< (N, O, Cl, F)", ),                 # 60
                ("H2C=C<", ),
                (">C=C< (ring) ", ),
                ("-C≡C-", ),
                ("HC≡C-", ),
                ("-O- (aromatic)", ),
                ("=N- (ring C5)", ),
                ("=N- (ring C6)", ),
                ("NO2-", ),
                ("NO2- (aromatic)", ),
                (">Si<", ),                                # 70
                (">Si< (3CH)", ),  # The definiton in Tb is >Si< (O)
                ("NO3-", ),
                (">PO4-", ),
                ("O=NO-", ),
                ("ONC-", ),
                ("-C=O-O-C=O-", ),
                ("COCl-", ),
                (">BO3-", ),
                ("COOO<", ),
                ("OCN-", ),                                # 80
                ("SCN-", ),
                ("-SO2-", ),
                (">Sn<", ),
                ("AsCl2-", ),
                ("-GeCl3 ", ),
                (">Ge<", ),
                (">C=C=C<", ),
                (">C=C-C=C<", ),
                (">C=C-C=C<", ),
                ("-CHO (aromatic)", ),                     # 90
                ("=N-", ),
                (">C=O (aromatic)", ),
                (">Si< (2CH)", ),  # The definiton in Tb is >Si< (F, Cl)
                ("-O-O-", ),
                (">C≡C-C≡C<", ),
                ("-C=O-O-C=O- (ring)", ),
                ("-NH- (aromatic)", ),
                ("=C< (aromatic)", ),  # Missing 98, set 214
                ("-OCON<", ),
                (">N-C=O-N<", ),                          # 100
                (">N<", ),
                ("F– (CCl)", ),
                ("-OCOO-", ),
                (">SO4", ),
                ("-SO2N<", ),
                ("=CNC=NC=", ),
                (">S=O", ),
                ("(S) -C≡N", ),
                (">N-C=O", ),
                # Missing 110, use a dupplicate of 78    # 110
                # with different molecular weigh
                (">Si< (F, Cl)", ),
                ("(N) -C≡N", ),
                (">Si< (1CH)", ),  # Missing 112, set 215
                (">P<", ),
                (">Si< (no CH)", ),  # Missing 114, set 216
                ("-ON=", ),
                (">Se<", ),
                (">Al<", ),
                # 118 in Tb equal to 134 in other, delete reference
                # So down to this is 1 index below

                # Second order contributions
                ("CO-CXy (X: F,Cl) (y>2)", ),
                ("CXy-CO-CXy (X: F,Cl) (y>2)", ),         # 120
                ("-CF3, -CCl3", ),
                (">CF2, >CCl2", ),
                ("No hydrogen", ),
                ("One hydrogen", ),
                ("3/4 Ring", ),
                ("Five-ring", ),
                ("Ortho pair(s)", ),
                ("Meta pair(s)", ),
                ("Para pair(s)", ),
                ("(C)(C=)>C<", ),                         # 130
                ("C2C-CC2", ),
                ("C3C-CC2", ),
                ("C3C-CC3", ),
                ("C=C-C=O", ),  # Same as 118 in Tb
                (">Si< (F,Cl,Br,I)", )]}  # 217 set to 135

    # Group interactions definition, Table 9 in [15]_
    # See Table 5 & 6 in [16]_ to check changes
    GI = {
        33: "A", 34: "A", 35: "A",
        36: "B",
        43: "C",
        37: "D",
        38: "E",
        44: "F", 45: "F", 46: "F",
        50: "G", 91: "G",
        51: "H", 89: "H",
        64: "I",
        53: "J",
        55: "K",
        52: "L",
        39: "M", 40: "M",
        41: "N", 96: "N",
        79: "O",
        56: "P",
        68: "Q",
        65: "R",
        66: "S"}

    # Table 5 in [15]_
    GI_Tb = {
        "AA": 291.7985, "AM": 314.6126, "AN": 286.9698, "AL": 38.6974,
        "AC": 146.7286, "AD": 135.3991, "AE": 226.4980, "AF": 211.6814,
        "AG": 46.3754, "AJ": -74.0193, "AP": 306.3979, "AI": 435.0923,
        "AS": 1334.6747, "BB": 288.6155, "BM": 797.4327, "BC": -1477.9671,
        "BD": 130.3742, "BF": -1184.9784, "BG": 0, "BH": 43.9722,
        "BQ": -1048.1236, "BS": -614.3624, "MM": 174.0258, "MN": 510.3473,
        "MD": 124.3549, "MF": 182.6291, "MJ": -562.3061, "MQ": 663.8009,
        "MI": 395.4093, "MS": 27.2735, "NN": 239.8076, "ND": 101.8475,
        "NF": 317.0200, "NG": -215.3532, "NS": 758.9855, "LL": 217.6360,
        "LA": 501.2778, "CC": 117.2044, "CD": 612.8821, "CF": -183.2986,
        "CG": -55.9871, "OO": -356.5017, "OQ": -263.0807, "DD": 91.4997,
        "DE": 178.7845, "DF": 322.5671, "DG": 15.6980, "DH": 17.0400,
        "DJ": 394.5505, "DQ": 963.6518, "DP": 293.5974, "DI": 329.0050,
        "EE": 1006.3880, "EH": 163.5475, "FF": 431.0990, "FG": 22.5208,
        "FQ": -205.6165, "FP": 517.0677, "FI": 707.9404, "GG": -303.9653,
        "GH": -391.3690, "GQ": -3628.9026, "GK": 381.0107, "GP": -574.2230,
        "GI": 176.5481, "GS": 124.1943, "HH": 582.1763, "HQ": 140.9644,
        "HK": 397.5750, "HI": 674.6858, "JJ": -11.9406, "QQ": 65.1432,
        "KP": -101.2319, "KR": -348.7400, "PS": -370.9729, "IR": -888.6123,
        "RR": 0, "SS": -271.9449}

    # Table 9 in [16]_
    GI_Tc = {
        "AA": -434.8568, "AM": 120.9166, "AN": -30.4354, "AD": -146.7881,
        "BB": 144.4697, "MM": -60.9217, "DM": -738.0515, "NN": -49.7641,
        "OO": -1866.0970, "DD": 162.6878, "DE": 707.4116, "DF": 128.2740,
        "DJ": -654.1363, "DP": 741.8565, "FF": 366.2663, "GG": 1605.5640,
        "JJ": -861.1528, "KR": 131.7924, "IR": 24.0243, "SS": -32.3208}

    # Table 12 in [16]_
    GI_Pc = {
        "AA": -5.6023, "AM": 69.8200, "AN": 6.1331, "AD": 7.3373,
        "BB": 57.8350, "MM": -0.6754, "DM": -125.5983, "NN": 22.1871,
        "DD": 2.6751, "DE": 88.8752, "DF": -1.0295, "DJ": 25.8246,
        "FF": 0.5195, "GG": -78.2743, "JJ": 43.9001, "KR": -19.7033,
        "IR": -35.1998, "SS": 12.5371}

    # Table 15 in [16]_
    GI_Vc = {
        "AM": -8.0423, "AD": 19.7707, "BB": 97.5425, "MM": -57.1233,
        "OO": 44.1062, "DD": -23.6366, "DE": -329.5074, "DF": -55.5112,
        "DJ": -37.2468, "FF": -74.8680, "GG": -413.3976, "JJ": -403.1196,
        "KR": 164.2930, "IR": 217.9243, "SS": -26.4556}

    # Table 6 in [17]_
    GI_Pv = {
        "AA": -561.5153, "AM": 1067.6660, "AN": 42.4825, "AD": -799.5332,
        "AE": -618.2760, "AF": -1797.6930, "AG": -1181.5990, "AP": 1431.2430,
        "AI": 1132.0400, "BB": -97.6205, "BD": -751.6676, "BF": 548.4352,
        "MM": 1085.8320, "MN": -206.7811, "MD": -198.2791, "MF": -1676.4770,
        "MQ": 1659.0340, "NN": -307.1018, "ND": 65.4421, "LL": 240.3037,
        "CC": -2601.7090, "CG": -787.8563, "OO": -3929.1300, "DD": 144.6074,
        "DE": 1118.9580, "DF": -225.7802, "DG": 1981.2980, "DH": 362.7540,
        "DJ": -1425.0170, "DP": 743.3353, "EE": -3748.8180, "FF": 920.3138,
        "FG": 1594.1640, "FP": 108.1305, "FI": 1590.3210, "GG": -1270.0830,
        "HH": 946.7309, "HI": 705.3049, "JJ": 838.3372, "QQ": -1501.3550,
        "KR": 675.0414, "PS": 994.4996, "IR": 135.5896, "SS": -29.6785}

    # Table 7 in [18]_
    GI_dBv = {
        "AA": -112.4939, "AM": 1031.5920, "AN": 853.2318, "AD": -423.9834,
        "AP": -683.0189, "AI": -557.5079, "BB": -1186.0500, "BD": -333.5638,
        "BQ": -878.0615, "MM": 135.3183, "MD": 219.9701, "ND": -134.4625,
        "DD": 132.0275, "DF": 44.8702, "DG": -219.5265, "DH": 546.5846,
        "DQ": -59.3635, "FF": 964.0840, "FG": 126.0380, "FP": 539.2401,
        "GG": 3705.4400, "HI": 50.1063, "QQ": 896.3606, "PS": -196.6361}

    # Table 10 in [18]_
    GI_Tv = {
        "AA": -1313.5690, "AM": -41.9608, "AN": -1868.6060, "AD": -643.4378,
        "AP": -345.7844, "AI": 50.2582, "BB": -1146.1070, "BD": -229.2406,
        "BQ": 515.1511, "MM": 86.7249, "MD": -57.1437, "ND": 54.2025,
        "DD": 156.7495, "DF": 273.6616, "DG": -339.6071, "DH": 1050.3190,
        "DQ": 355.0508, "FF": 167.7204, "FG": 244.0583, "FP": 334.4856,
        "GG": 1985.8270, "HI": 161.7447, "QQ": 1839.2630, "PS": 718.1262}

    FirstOrder = 115
    SecondOrder = 134

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        nh = self._atomX("H")
        na = self._atoms()
        n = na-nh  # Total atoms of compound except hydrogen

        tb, tc, pc, vc, pv, dbv, tv = 0, 0, 0, 0, 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb += c*self._coeff["tb"][i]
            tc += c*self._coeff["tc"][i]*1e-3
            pc += c*self._coeff["Pc"][i]*1e-4
            vc += c*self._coeff["vc"][i]
            pv += c*self._coeff["Pv"][i]*1e-3
            dbv += c*self._coeff["dBv"][i]*1e-3
            tv += c*self._coeff["Tv"][i]

        # Group interaction calculation
        GI = []
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i in self.GI:
                for x in range(c):
                    GI.append(self.GI[i])
        m = len(GI)

        for pair in permutations(GI, 2):
            key = "".join(sorted(pair))
            tb += self.GI_Tb.get(key, 0)/n/(m-1)
            tc += self.GI_Tc.get(key, 0)*1e-3/n/(m-1)
            pc += self.GI_Pc.get(key, 0)*1e-4/n/(m-1)
            vc += self.GI_Vc.get(key, 0)/n/(m-1)
            pv += self.GI_Pv.get(key, 0)*1e-3/n/(m-1)
            dbv += self.GI_dBv.get(key, 0)*1e-3/n/(m-1)
            tv += 2*self.GI_Tv.get(key, 0)/n/(m-1)

        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb = unidades.Temperature(tb/(n**0.6583+1.6868)+84.3395)
        self.Tc = unidades.Temperature(self.Tb*(0.699+1/(0.9889+tc**0.8607)))
        self.Pc = unidades.Pressure(self.M**-0.14041/(0.00939+pc)**2, "kPa")
        self.Vc = unidades.SpecificVolume((vc/n**-.2266+86.1539)/self.M, "ccg")

        # Eq 7 in [17]_
        self.db = pv - 0.176055

        # Eq 7 in [18]_
        self.dBv = dbv/(n**-2.5635+0.0685)+3.777

        # Eq 8 in [18]_
        Tv = 21.8444*self.Tb**0.5+tv**0.9315/(n**0.6577+4.9259)-231.1361
        self.Tv = unidades.Temperature(Tv)

        GroupContribution.calculo(self)

    def _Pv(self, T):
        """Vapor pressure calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        Pv : float
            Vapor Pressure, [Pa]
        """
        # Eq 6 in [17]_
        Trb = T/self.Tb
        Pv = 10**((4.1012+self.db)*((Trb-1)/(Trb-0.125)))
        return unidades.Pressure(Pv, "atm")

    def _Visco(self, T):
        """Viscosity calculation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        mu : float
            Viscosity, [Pas]
        """
        mu = 1.3*exp(-self.dBv*(T-self.Tv)/(T-self.Tv/16))
        return unidades.Viscosity(mu, "cP")


@refDoc(__doi__, [19])
class Wen(GroupContribution):
    """
    Group contribution for definition of unknown component using the Wen-Qiang
    procedure (2001)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    M: float, optional
        Molecular weight, [-]
    Tb : float, optional
        Normal boiling temperature, [K]
    SG: float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 1 in [19]_, Tc of n-Butylaniline
    The last group containing a carbon-adjacent atom pair has a typo, must be
    =C<[r]/>N-

    >>> cmp = Wen(Tb=513.9, group=[1, 11, 18, 87, 96, 99, 135],
    ... contribution=[1, 5, 1, 10, 2, 1, 1])
    >>> "%0.3f" % cmp.Tc
    '721.275'
    >>> cmp.formula
    'C10H15N'

    Example 2 in [19]_, Pc of Benzoic acid

    >>> cmp = Wen(group=[87, 93, 96, 114], contribution=[10, 1, 1, 1])
    >>> "%0.3f" % cmp.Pc.MPa
    '4.547'

    Example 3 in [19]_, Vc of chloropentafluorobenzene

    >>> cmp = Wen(group=[96, 101, 102], contribution=[12, 5, 1])
    >>> "%0.1f" % (cmp.Vc.ccg*cmp.M)
    '374.6'
    >>> cmp.formula
    'C6F5Cl'
    """
    __title__ = "Wen-Qiang (2001)"
    _coeff = {
        # Table III
        "tc": [-2.885, 2.424, 0.048, 22.766, -3.404, 2.495, 2.275, 2.602,
               -1.601, 0.000, 35.848, 2.124, -0.708, 22.576, -3.085, 1.578,
               -0.030, 2.256, -2.322, 2.549, 13.769, 20.882, 25.177, 1.934,
               -3.377, 0.000, 0.000, -1.765, -3.163, -5.588, 0.830, 11.483,
               2.183, -7.415, 0.007, -0.651, 0.178, -4.384, -0.502, 6.664,
               9.639, 2.262, -1.648, 0.000, 5.838, 2.218, 10.659, -2.228,
               0.017, 3.541, -0.748, -1.248, 9.056, 4.564, 2.737, 0.007, 1.136,
               -0.008, 13.166, -0.009, -2.298, 9.242, 1.143, -11.918, 0.007,
               -0.953, 0.000, 0.000, 5.993, 3.162, 0.000, 2.707, 3.180, 6.070,
               13.625, 7.842, 5.897, 0.004, -781.237, 2.367, 0.000, -7.274,
               3.798, 5.571, 2.446, -1.223, 0.754, 2.852, 2.013, 7.937, 14.661,
               10.141, -0.603, 2.172, 0.000, 0.009, 4.660, -2.465, 3.701,
               1.700, 6.344, -5.547, 5.600, 12.840, 28.472, 10.144, 18.220,
               15.436, 20.655, 11.326, 13.047, 34.349, 35.591, 34.476, 35.009,
               7.041, 9.909, 9.193, 1.933, 16.126, 23.640, 9.697, 15.542,
               3.835, 11.598, 4.679, 22.777, 18.847, 8.357, 0.003, 11.032,
               14.462, 16.658, 1.251, 15.982, 17.609, 16.935, 16.840, 0.001,
               0.004],
        "Pc": [9.361, 4.035, 2.682, 7.267, 9.548, 3.297, 1.286, 2.797, 5.280,
               0.000, 2.517, 1.454, 0.314, 5.065, 9.483, 2.239, 0.038, 1.497,
               3.545, 1.851, 1.707, 1.517, 3.345, 0.460, 0.358, 0.000, 6.274,
               1.087, 0.187, 3.501, 2.770, 1.318, -0.282, 0.141, 6.690, -0.032,
               -0.013, 0.214, 2.885, 2.880, 2.604, 3.648, 4.393, 0.000, 2.686,
               0.731, 5.158, 3.462, 0.000, -0.290, -1.534, 3.842, 2.532,
               1.624, 0.305, 0.000, 2.352, 0.000, -0.864, -0.854, 3.177,
               2.036, 1.668, 0.000, -0.001, -3.595, 0.000, 0.000, -1.661,
               -0.802, 5.880, 1.006, 2.175, 0.535, 0.000, -0.152, -4.382,
               0.000, 0.000, 0.464, 0.000, 4.682, 1.848, -5.905, 0.883, 1.355,
               -0.296, 0.680, 0.853, 0.275, -0.001, 0.189, 1.671, 2.011, 4.716,
               -0.004, -0.064, -0.225, -1.345, 5.744, -0.319, 3.029, 3.476,
               4.375, 3.073, 1.662, 2.422, 1.621, 0.037, 1.516, -0.543, -2.792,
               1.524, -0.406, 0.855, -0.026, 2.966, 1.840, 2.188, 1.785,
               -0.332, -0.622, 5.780, 0.078, 0.001, 1.800, 0.000, -1.320,
               4.158, 0.005, -1.452, -2.250, -0.001, 0.000, -3.152, -4.585,
               -7.394, 0.000, 0.000, -0.002],
        "vc": [125.58, 86.72, 70.64, 105.31, 45.38, 91.08, 62.82, 113.65,
               49.86, 0.00, 200.24, 28.22, 14.98, 44.98, 0.07, 46.54, 13.89,
               57.95, -20.09, 81.10, 103.74, 127.10, 0.00, 5.53, 2.54, 0.00,
               0.00, 16.98, 1.47, -46.51, 55.26, 84.21, -4.52, -16.38, 0.00,
               0.00, -0.08, -3.45, 44.84, 75.12, 98.04, 78.05, 80.34, 0.00,
               38.30, 23.00, 136.00, 28.40, 0.00, 0.00, 21.77, 68.06, 98.13,
               20.29, 5.01, 0.00, 3.48, 0.00, 62.65, -99.32, 48.52, 77.32,
               69.56, 0.01, -95.79, 18.23, 0.00, 0.00, 14.19, 77.52, 0.00,
               27.84, 46.55, 38.22, 0.00, 50.59, -0.02, 0.00, 0.00, 34.33,
               0.00, 0.00, 36.94, 0.04, 36.83, 4.71, 0.00, 23.60, 51.71, 52.22,
               0.00, 45.09, -2.49, 112.28, 0.00, 0.00, 9.87, -49.16, 5.85,
               -64.77, 5.54, 40.93, 78.55, 95.28, 122.28, 0.00, 113.89, 0.00,
               0.00, 46.35, -81.99, 0.00, 127.31, 112.99, 0.00, 74.01, 60.54,
               39.13, -68.36, 0.00, -28.91, 17.46, 165.84, 94.53, -0.01, 0.00,
               0.00, 58.88, 46.41, 0.01, 101.11, 59.67, 0.00, 0.00, 110.06,
               138.33, 131.47, 0.00, 0.00, -0.01, ],
        "tc_": [-8.8072, -1.1863, 2.0695, 0.8880, 0.4312, -2.9673, -5.6886,
                -1.8098, -0.5794, -0.0174, 4.5913, 2.2116, 5.2478, 0.5832,
                3.3562, -0.2371, -0.0296, 2.8167, 2.9084, 1.0111, -0.8130,
                -1.3763, -0.0022, 2.7600, 7.5553, 0.0000, 0.0000, 0.3694,
                0.9548, 3.9991, 1.6760, 0.5273, 2.7770, 5.3585, -0.0076,
                0.8857, -0.0093, 0.7629, 0.7647, 1.6201, 0.4661, 1.7766,
                4.6760, 0.0000, -1.4353, 1.5808, -5.4993, 13.0182, -0.1973,
                -5.5762, -2.5308, -1.7201, -0.9027, 0.4717, 0.9585, -0.0014,
                6.5524, 0.0038, -4.2721, 0.0060, 1.3320, 0.6099, -1.1118,
                -6.0111, -0.0038, 4.4394, 0.0000, 0.0000, -6.3639, -0.8632,
                0.0000, 0.5719, 1.2736, -0.1301, -0.9306, 0.1783, 2.0133,
                0.0204, 0.0000, 0.0094, 0.0000, 0.7476, 2.4465, 4.1358,
                -1.5784, 4.0418, -1.1378, 0.3159, -0.1566, -1.1410, 0.4046,
                -0.0075, 5.6200, 7.6978, 0.0000, -0.0180, 1.2080, -4.5849,
                0.8675, 7.0762, 2.773, 0.7764, 1.3463, 0.2835, 0.0770, 0.6733,
                -1.1847, 5.4387, -4.3184, 1.0536, -3.2089, 6.8866, 12.5998,
                14.3778, 10.1056, 2.6863, 5.4168, 11.107, 8.7813, 15.4803,
                -0.4178, 2.5145, 9.8751, 9.4595, 15.4000, 1.3670, 7.2107,
                12.1779, 11.7647, 0.0013, 6.3813, -4.7655, 2.7693, -0.5727,
                1.4762, 4.2350, 7.5331, 5.1503, 0.0019, 0.0590],

        # Name and element composition
        "txt": [("CH3- & H",),                          # 0
                ("CH3- & >C<",),
                ("CH3- & =C<",),
                ("CH3- & ≡C-",),
                ("CH3- & >C< [r]",),
                ("CH3- & =C< [r]",),
                ("CH3- & -O-",),
                ("CH3- & -S-",),
                ("CH3- & >N-",),
                ("CH3- & =N-",),
                ("CH3- & -NO2",),                        # 10
                ("-CH2- & >C<",),
                ("-CH2- & =C<",),
                ("-CH2- & ≡C-",),
                ("-CH2- & >C< [r]",),
                ("-CH2- & =C< [r]",),
                ("-CH2- & -O-",),
                ("-CH2- & -S-",),
                ("-CH2- & >N-",),
                ("-CH2- & F-",),
                ("-CH2- & Cl-",),                        # 20
                ("-CH2- & Br-",),
                ("-CH2- & I-",),
                (">CH- & >C<",),
                (">CH- & =C<",),
                (">CH- & ≡C-",),
                (">CH- & >C< [r]",),
                (">CH- & =C< [r]",),
                (">CH- & -O-",),
                (">CH- & >N-",),
                (">CH- & F-",),                          # 30
                (">CH- & Cl-",),
                (">C< & >C<",),
                (">C< & =C<",),
                (">C< & ≡C-",),
                (">C< & >C< [r]",),
                (">C< & =C< [r]",),
                (">C< & -O-",),
                (">C< & F-",),
                (">C< & Cl-",),
                (">C< & Br-",),                          # 40
                ("=CH2 & =C<",),
                ("=CH2 & =C=",),
                ("=CH2 & =C< [r]",),
                ("=CH- & >C<",),
                ("=CH- & =C<",),
                ("=CH- & ≡C-",),
                ("=CH- & =C=",),
                ("=CH- & >C< [r]",),
                ("=CH- & =C< [r]",),
                ("=CH- & -O-",),                         # 50
                ("=CH- & F-",),
                ("=CH- & Cl-",),
                ("=C< & >C<",),
                ("=C< & =C<",),
                ("=C< & ≡C-",),
                ("=C< & =C=",),
                ("=C< & =C< [r]",),
                ("=C< & =O",),
                ("=C< & >N-",),
                ("=C< & F-",),                           # 60
                ("=C< & Cl-",),
                ("≡CH & ≡C-",),
                ("≡C- & >C<",),
                ("≡C- & =C<",),
                ("≡C- & ≡C-",),
                ("≡C- & >C< [r]",),
                ("≡C- & =C< [r]",),
                ("=C= [r] & =C<",),
                ("=C= [r] & =O",),
                ("=C= [r] & =N-",),                      # 70
                ("-CH2- [r] & >C< [r]",),
                ("-CH2- [r] & =C< [r]",),
                ("-CH2- [r] & -O- [r]",),
                ("-CH2- [r] & -S- [r]",),
                ("-CH2- [r] & >N- [r]",),
                (">CH- [r] & >C<",),
                (">CH- [r] & =C<",),
                (">CH- [r] & ≡C-",),
                (">CH- [r] & >C< [r]",),
                (">CH- [r] & =C< [r]",),                  # 80
                (">CH- [r] & -O-",),
                (">CH- [r] & -O- [r]",),
                (">C< [r] & >C<",),
                (">C< [r] & >C< [r]",),
                (">C< [r] & F-",),
                ("=CH- [r] & >C< [r]",),
                ("=CH- [r] & =C< [r]",),
                ("=CH- [r] & -O- [r]",),
                ("=CH- [r] & -S- [r]",),
                ("=CH- [r] & >N- [r]",),                  # 90
                ("=CH- [r] & =N- [r]",),
                ("=C< [r] & >C<",),
                ("=C< [r] & =C<",),
                ("=C< [r] & ≡C-",),
                ("=C< [r] & >C< [r]",),
                ("=C< [r] & =C< [r]",),
                ("=C< [r] & -O-",),
                ("=C< [r] & -O- [r]",),
                ("=C< [r] & >N-",),
                ("=C< [r] & >N- [r]",),                   # 100
                ("=C< [r] & F-",),
                ("=C< [r] & Cl-",),
                ("=C< [r] & Br-",),
                ("=C< [r] & I-",),
                ("-CHO & -H",),
                ("-CHO & >C<",),
                ("-CHO & =C<",),
                ("-CHO & =C< [r]",),
                ("-CO & >C<",),
                ("-CO & =C< [r]",),                       # 110
                ("-COOH & -H",),
                ("-COOH & >C<",),
                ("-COOH & =C<",),
                ("-COOH & =C< [r]",),
                ("-COO- & -H",),
                ("-COO- & >C<",),
                ("-COO- & =C<",),
                ("-COO- & =C< [r]",),
                (">C2O3 & =C<",),
                (">C2O3 & =C< [r]",),                     # 120
                (">CO [r] & >C< [r]",),
                ("-CN & -H",),
                ("-CN & >C<",),
                ("-CN & =C<",),
                ("-CN & ≡C-",),
                ("-CN & =C< [r]",),

                # Group Structural Unit Parameter
                ("-OH", ),
                ("-O-", ),
                ("-O- (cyclic)", ),
                ("-OH (Aromatic)", ),                     # 130
                ("-SH", ),
                ("-S-", ),
                ("-S- (cyclic)", ),
                ("-NH2", ),
                (">NH", ),
                (">N-", ),
                ("=N-", ),
                (">NH (cyclic)", ),
                (">N- (cyclic)", )]}

    FirstOrder = 127
    SecondOrder = 140

    def calculo(self):
        """Calculation procedure"""
        # Use the input properties
        # SG is defined in base class
        if self.kwargs["M"]:
            M = self.kwargs["M"]
        else:
            M = self._M()
        self.M = unidades.Dimensionless(M)

        tc, pc, vc, tc_ = 0, 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            tc_ += c*self._coeff["tc_"][i]
            pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]

        if self.kwargs["Tb"]:
            self.Tb = unidades.Temperature(self.kwargs["Tb"])
            Tc = (((127.754+tc_)/1e2)**-2+1)*self.Tb
            self.Tc = unidades.Temperature(Tc)
        else:
            self.Tc = unidades.Temperature(((4.72+tc)*1e6)**(1/2.747))
        self.Pc = unidades.Pressure(((37.293+pc)/1e2)**-2, "MPa")
        self.Vc = unidades.SpecificVolume((-27.04+vc)/self.M, "ccg")

        GroupContribution.calculo(self)

    def _group(self):
        """Specific procedure to calculate the molecular weight of compound
        from group contribution"""
        # TODO: Add group check coherent definition
        group = []
        rest = {}
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if "&" in self._coeff["txt"][i][0]:
                grp = self._coeff["txt"][i][0].split(" & ")[0]
                second = self._coeff["txt"][i][0].split(" & ")[1]

                # Discard second term with carbons and add the heteroatoms term
                grp2 = atomic_decomposition(second)
                if "H" in grp2 or "F" in grp2 or "Cl" in grp2 or \
                        "Br" in grp2 or "I" in grp2:
                    for x in range(c):
                        group.append(grp2)
            else:
                grp = self._coeff["txt"][i][0]

            for x in range(c):
                cmp = grp
                molecule = grp
                # Clean additional comment of group, ring flag and other,
                # separated of main group by spaces
                if " " in cmp:
                    cmp = cmp.split(" ")[0]

                # Calculate the remain appeareances of group contributions
                # X- : X is one connected so don't appearance reamaining
                # -X- : X has 2 connection so must be stay in another group
                # -X<: X has 3 connection so must be stay in other 2 group
                # >X<: X has 4 connection so must be stay in other 3 group
                restLink = -1
                if cmp[0] in ("-=≡"):
                    restLink += 1
                    molecule = molecule[1:]
                elif cmp[0] == ">":
                    restLink += 2
                    molecule = molecule[1:]
                if cmp[-1] in ("-=≡"):
                    restLink += 1
                    molecule = molecule[:-1]
                elif cmp[-1] == "<":
                    restLink += 2
                    molecule = molecule[:-1]

                if grp in rest and rest[grp]:
                    rest[grp] -= 1
                else:
                    rest[grp] = restLink
                    group.append(atomic_decomposition(molecule))
        return group


@refDoc(__doi__, [21])
class Li(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Li-Xia-Xiang procedure (2016)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    Tb : float
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example A.1 in [21]_, critical temperature of hexanal

    >>> cmp = Li(Tb=401.45, group=[0, 1, 2, 10, 14, 17],
    ... contribution=[6, 12, 1, 5, 1, 12])
    >>> "%0.2f" % cmp.Tc
    '588.90'
    >>> cmp.formula
    'C6H12O'

    Example A.2 in [21]_, critical temperature of chlorotrimethylsilane

    >>> cmp = Li(Tb=330.75, group=[0, 1, 6, 9, 24, 25, 17],
    ... contribution=[3, 9, 1, 1, 3, 1, 9])
    >>> "%0.2f" % cmp.Tc
    '497.69'

    Example A.3 in [21]_, critical pressure of n-propyl formate

    >>> cmp = Li(Tb=353.97, group=[0, 1, 2, 10, 14, 17, 13],
    ... contribution=[4, 8, 2, 2, 1, 8, 2])
    >>> "%0.3f" % cmp.Pc.MPa
    '4.033'

    Example A.4 in [21]_, critical volume of acetic acid
    The paper has a bug, the calculated value is erroneous

    >>> cmp = Li(Tb=351.44, group=[0, 1, 2, 10, 14, 17, 16],
    ... contribution=[2, 4, 2, 1, 1, 3, 1])
    >>> "%0.2f" % (cmp.Vc.ccg*cmp.M)
    '180.65'
    """
    __title__ = "Li-Xia-Xiang (2016)"

    _coeff = {
        "tc": [-0.0003, -0.0016, -0.0472, -0.0533, 0.0179, -0.0478, -0.0380,
               0.0071, 0.0091, 0.0224, -0.0168, -0.0080, -0.0079, 0.0059,
               0.0200, 0, -0.0114, -0.0001, 0.0285, 0.0358, 0.0071, 0.0091,
               -0.0153, 0.0360, -0.0316, 0.0209, 0.0268, 0.0029, -0.0192,
               0.0152, -0.0024, 0, 0.0180, 0.0171, 0, 0.0528, -0.0061, -0.0055,
               0.0238, 0.0022, 0.0209, 0.0285, -0.0398],
        "Pc": [0.3304, 0.0118, 1.0330, 0.4774, 0.5733, 0.9346, 1.1961, 0.9630,
               1.4724, 0.0654, 0.3317, 0.2350, 0.0604, -0.1030, -0.6309,
               0.0000, -0.6777, 0.0179, -0.2329, -0.1194, 0.9630, 1.4724,
               0.2999, 0.2856, 0.6898, 0.1418, 0.1493, 0.0000, 0.0000, 0.0973,
               0.2519, 0.0000, -0.0060, -0.1233, 0.0000, -0.6638, 0.1793,
               0.0368, -0.2843, 0.1783, -0.0381, -0.1360, 0.4207],
        "vc": [60.0289, 5.6116, 60.5913, -58.5859, 58.2972, 51.3276, 90.5686,
               88.0064, 128.9134, 30.6559, 5.3217, 3.3178, -5.4347, -5.3629,
               -15.8297, 0.0000, -12.3706, -1.3769, -9.4243, 8.1548, 88.0064,
               128.9134, 29.1266, 29.0457, 27.0244, 46.3533, 49.2458, 0, 0,
               19.6420, 48.3117, 0, 116.9220, 39.4229, 0, 23.5185, -2.1851,
               -15.5277, -14.5939, 9.7344, 28.5488, 91.7583, -50.7082],

        "txt": [("C",),                          # 0
                ("H",),
                ("O",),
                ("N",),
                ("S",),
                ("F",),
                ("Cl",),
                ("Br",),
                ("I",),
                ("Si",),

                # 2nd Order term
                ("C-C",),                        # 10
                ("C=C",),
                ("C≡C",),
                ("C-O",),
                ("C=O",),
                ("O-O",),
                ("O-H",),
                ("C-H",),
                ("C-F",),
                ("C-Cl",),
                ("C-Br",),                       # 20
                ("C-I",),
                ("C-S",),
                ("C=S",),
                ("C-Si",),
                ("Si-Cl",),
                ("Si-O",),
                ("Si-H",),
                ("S-S",),
                ("S-H",),
                ("C-N",),                        # 30
                ("C=N",),
                ("C≡N",),
                ("N-H",),
                ("N-N",),
                ("N=O",),
                ("C-C [r]",),
                ("C=C [r]",),
                ("C-O [r]",),
                ("C-S [r]",),
                ("C-N [r]",),                    # 40
                ("C=N [r]",),
                ("benzene",)]}

    FirstOrder = 10
    SecondOrder = 43

    def isCalculable(self):
        """Procedure to define the status of input parameter"""
        if not self.kwargs["Tb"]:
            self.msg = QApplication.translate(
                    "pychemqt", "undefined boiling point")
            self.status = 0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.Tb = unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        tc, pc, vc = 0, 0, 0
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc += c*self._coeff["tc"][i]
            pc += c*self._coeff["Pc"][i]
            vc += c*self._coeff["vc"][i]

        self.Tc = unidades.Temperature(self.Tb*(1.5530+tc)+18.9999)
        self.Pc = unidades.Pressure(self.M/(1.2220+pc)**2, "MPa")
        self.Vc = unidades.SpecificVolume(
                (19.6531-1.2603*self.M+vc)/self.M, "ccg")

        GroupContribution.calculo(self)


@refDoc(__doi__, [22])
class MarreroGani(GroupContribution):
    """
    Group contribution for definition of unknown component using the
    Marrero-Gani procedure (2001)

    The resulting instance has all the necessary properties to use in PFD as a
    predefined compound.

    Parameters
    ----------
    group : array
        List with group index
    contribution : float
        List with group count ocurrences
    Tb : float, optional
        Normal boiling temperature, [K]
    M : float, optional
        Molecular weight, [-]
    SG : float, optional
        Specific gravity, [-]

    Notes
    -----
    Tb, M and SG are optional input, anyway know them improve the estimation

    Examples
    --------
    Example 1 from [22]_, Tb of N-phenyl-1,4-benzenodiamine

    >>> c1 = MarreroGani(group=[61, 62, 17, 14], contribution=[1, 1, 1, 9])
    >>> c2 = MarreroGani(group=[61, 62, 17, 14, 287],
    ... contribution=[1, 1, 1, 9, 1])
    >>> c3 = MarreroGani(group=[61, 62, 17, 14, 287, 346],
    ... contribution=[1, 1, 1, 9, 1, 1])
    >>> "%0.2f %0.2f %0.2f" % (c1.Tb, c2.Tb, c3.Tb)
    '614.62 616.03 623.94'

    Example 2 from [22]_, Tb of pyrene

    >>> c1 = MarreroGani(group=[15, 14], contribution=[6, 10])
    >>> c3 = MarreroGani(group=[15, 14, 361, 365], contribution=[6, 10, 2, 2])
    >>> "%0.2f %0.2f" % (c1.Tb, c3.Tb)
    '652.43 673.96'

    Example 3 from [22]_, Tb of 4-aminobutanol

    >>> c1 = MarreroGani(group=[28, 53, 1], contribution=[1, 1, 3])
    >>> c3 = MarreroGani(group=[28, 53, 1, 306], contribution=[1, 1, 3, 1])
    >>> "%0.2f %0.2f" % (c1.Tb, c3.Tb)
    '448.64 478.42'

    Example 4 from [22]_, Tm of 3,3'-methylenebis-4-hydroxycoumarin

    >>> c1 = MarreroGani(group=[28, 16, 14, 172, 179, 178, 1],
    ... contribution=[2, 4, 8, 2, 2, 2, 1])
    >>> c3 = MarreroGani(contribution=[2, 4, 8, 2, 2, 2, 1, 2, 2, 2, 1],
    ... group=[28, 16, 14, 172, 179, 178, 1, 354, 317, 353, 330])
    >>> "%0.2f %0.2f" % (c1.Tf, c3.Tf)
    '475.46 552.34'

    Example 5 from [22]_, Tm of flutemazepan
    The third order example has a undocumented contribution, using value from
    AROMFUSED[2]s2 contribution

    >>> c1 = MarreroGani(
    ... group=[28, 122, 123, 0, 17, 16, 14, 168, 179, 175, 177],
    ... contribution=[1, 1, 1, 1, 1, 2, 7, 1, 1, 1, 1])
    >>> c2 = MarreroGani(
    ... group=[28, 122, 123, 0, 17, 16, 14, 168, 179, 175, 177, 285, 283, 265],
    ... contribution=[1, 1, 1, 1, 1, 2, 7, 1, 1, 1, 1, 1, 1, 1])
    >>> c3 = MarreroGani(
    ... group=[28, 122, 123, 0, 17, 16, 14, 168, 179, 175, 177, 285, 283, 265,
    ... 356, 345, 347, 349],
    ... contribution=[1, 1, 1, 1, 1, 2, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    >>> "%0.2f %0.2f %0.2f" % (c1.Tf, c2.Tf, c3.Tf)
    '469.19 473.37 419.63'

    Example 6 from [22]_, Tm of 1,9-nonadiol

    >>> c1 = MarreroGani(group=[28, 1], contribution=[2, 9])
    >>> c3 = MarreroGani(group=[28, 1, 307], contribution=[2, 9, 1])
    >>> "%0.2f %0.2f" % (c1.Tf, c3.Tf)
    '303.66 315.70'
    """
    __title__ = "Marrero-Gani (2001)"

    # Table 6-8
    _coeff = {
        "tf": [0.6953, 0.2515, -0.3730, 0.0256, 1.1728, 0.9460, 0.7662, 0.1732,
               0.3928, 1.7036, 1.5453, 1.2850, 2.2276, 2.0516, 0.5860, 1.8955,
               1.2065, 0.9176, 2.0438, 1.0068, 0.1065, -0.5197, -0.1041,
               1.2832, 1.7744, 1.2612, 1.7495, 0, 2.7888, 5.1473, 7.4042,
               12.4296, 2.9588, 2.5232, 1.1565, 1.0638, 2.9157, 3.0186, 2.4744,
               2.1657, 1.6329, 1.0668, 0.3983, 2.0223, 1.3348, 0, 4.8044,
               1.5038, 1.3643, 0.8733, 0.2461, -0.4446, 1.3045, 3.2742,
               30.8394, 11.7400, 2.4034, 1.7746, 1.7577, 0.9607, 0.0442,
               3.9889, 1.4837, 1.7618, 3.3478, 8.8492, 1.4621, 2.5760, 2.1393,
               3.3807, 5.1346, 3.2747, 4.2256, 0, 9.1492, 2.2327, 3.2131,
               0.7812, 5.6280, 4.3531, 3.0376, 0, 2.5974, 0, 0, 13.2124,
               5.4720, 5.8825, 4.1720, 0, 0, 9.1763, 3.2657, 12.8071, 5.6631,
               3.3602, 6.5160, 9.8204, 7.2552, 9.3110, 14.2020, 13.0856,
               8.4447, 3.5041, 13.4695, 23.2570, 3.0882, 1.9253, 1.0224,
               1.8424, 2.5196, 3.6491, 4.4493, 1.5597, 1.1289, 2.5398,
               2.1689, 0.1312, 1.4828, 3.2035, 0, 1.7510, 1.7134, 0.9782,
               2.1905, 2.4741, 1.9444, 1.7641, 1.2308, 1.5454, 3.9813, 3.5484,
               10.5579, 2.3651, 0, 1.5791, 4.8181, 2.2992, 0.9704, 4.2329,
               2.8464, 0.9600, 1.7150, 1.0063, 0.7892, 1.1170, 0.9646, 5.3663,
               7.0778, 0, 5.8426, 3.6976, 3.9911, 5.2948, 0, 0, 1.0306, 0, 0,
               2.7461, 2.0330, -1.7840, 0.2337, 3.6593, 1.3135, 0, 0, 0.5699,
               0.0335, 0.1695, 1.1936, 0.4344, 0.3048, 0.2220, 3.4814, 0.6040,
               5.5779, 6.6382, 1.3828, 3.2119, 1.6023, 6.1006, 0.1175, -0.1214,
               0.2390, -0.3276, 3.3297, 0.7451, 0.0524, -0.1077, -0.2485,
               0.5715, -0.0968, -0.6024, -3.1734, 0.2114, -1.2441, -0.3489,
               0.3695, 0.9886, -1.1810, -0.1526, -0.0414, -0.5941, 0.3258,
               -1.8403, 11.5351, 12.3481, 0.9327, 7.5057, -0.4531, 14.1593,
               -2.3026, -2.1535, -2.7514, 4.0747, -0.9493, 0.2232, -0.5946,
               2.5962, -0.2509, 0.6304, -0.0811, -0.2568, -0.4329, 0, 0.0446,
               0.1027, 0.2117, -0.7191, 2.4103, 0.8092, -1.0802, 0.8607,
               0.8981, 0.1088, 1.9470, 1.2057, 0.3666, -0.2363, 0.4506, 0, 0,
               2.2421, -0.6997, -0.2636, -1.1057, 0.0642, 0.0790, -10.8058,
               -1.0516, -6.9427, 0.6572, 0.0416, -0.3151, 1.5819, -0.8604,
               -0.1326, -0.4669, -0.3548, -0.1727, 0.6817, -1.0631, 0.5124,
               2.8497, 1.3691, 1.5069, 0.0370, 0, -0.3312, 0, -2.0822, 0.7743,
               -0.8578, -0.8638, 0.5076, -0.3978, 0, 0, -0.4666, 0.1737,
               -1.9233, 0.7334, -0.0383, 1.0497, -0.6388, -0.6218, 0.9840,
               -0.2762, -0.3689, -0.3841, 1.7722, 0.4553, 2.0561, -0.5769,
               -0.2556, 1.6282, -0.1341, -1.6848, -0.9802, 0.3018, 0.1018,
               0.2811, -0.3189, -1.5257, 11.2271, 0.7732, 0.6674, -0.1073,
               -1.3891, -0.0781, 0, -0.0604, -1.1888, 0.6669, -0.3798, -2.6542,
               0.2479, 1.1395, 0.0570, -0.5640, 1.9902, 0, 0.5460, 0.4497,
               0.6647, 0.0792, 0.6457, 0.9608, 16.2235, 16.8558, -1.0394,
               -0.4486, -0.1376, 0.4361, 3.6847, 4.9038, -7.0038, 5.9653,
               1.5629, -9.1856, 0.2612, -1.8403, 0.0393, 0.9514, 3.4983,
               -0.3048, -1.3060, -4.9289, -10.1007, 1.0834, -0.4803, -3.2676,
               -0.3545, 0.2825, -1.2836, 0.3378, 1.8941, -2.7585, -3.0362,
               -3.2228, 1.6600, 7.0402, -3.3463, 6.8373, -1.5856, 2.0821,
               -4.4725, -2.5898, 1.0358],
        "tb": [0.8491, 0.7141, 0.2925, -0.0671, 1.5596, 1.5597, 1.3621, 1.2971,
               1.2739, 2.6840, 2.4014, 2.5400, 1.7618, 1.6767, 0.8365, 1.7324,
               1.1995, 1.5468, 1.3977, 1.5653, 1.4925, 0.8665, 0.5229, 2.4308,
               2.9262, 2.1472, 2.3057, 2.7341, 2.5670, 3.3205, 5.1108, 6.0677,
               3.1178, 2.6761, 2.1748, 1.7287, 3.4650, 2.5388, 3.5172, 3.1228,
               2.9850, 2.2869, 1.6918, 2.5972, 3.1952, 0.4621, 3.0854, 2.1903,
               1.7703, 1.3368, 0.8924, 0.4983, 1.8522, 2.7987, 2.0948, 1.6525,
               2.2514, 1.8750, 1.2317, 1.3841, 1.1222, 3.8298, 2.9230, 2.1918,
               2.0315, 1.5332, 1.4291, 4.5871, 3.9774, 2.8870, 4.1424, 3.0972,
               3.4891, 3.1220, 0, 3.1853, 4.5311, 3.8069, 3.3059, 4.5750,
               3.2069, 1.8896, 3.2656, 5.8779, 7.4566, 6.5652, 5.0724, 6.6810,
               6.0070, 0, 5.0664, 7.6172, 5.6487, 8.3775, 7.3497, 5.1373,
               7.5850, 7.4955, 0, 8.9406, 16.3539, 2.0796, 7.1529, 4.1459,
               5.7604, 1.1633, 0, 2.6364, 2.0246, 1.7049, 3.3420, 2.9609,
               3.9093, 1.5022, 1.3738, 1.0084, 2.2238, 0.5142, 1.1916, 2.5053,
               2.0542, 1.7227, 2.0669, 0.7945, 3.7739, 2.8414, 3.1778, 2.4231,
               0.8504, 1.5147, 4.5721, 4.0142, 0, 4.8721, 4.2329, 3.6653,
               3.1669, 3.1974, 2.5910, 2.0902, 3.2675, 2.3323, 2.9892, 2.6524,
               2.0965, 1.6412, 2.9731, 6.2796, 7.0976, 3.9199, 6.7785, 5.5627,
               6.1185, 8.4333, 2.0536, 1.0984, 2.7900, 5.6433, 4.5468, 5.1567,
               3.7657, 2.3522, 2.9272, 2.8847, 2.8451, 2.6124, 2.2036, 0.8234,
               0.5946, 0.0386, 1.5985, 1.2529, 1.1975, 1.5109, 2.1634, 1.6541,
               6.5230, 6.6710, 1.0245, 2.8793, 2.3256, 0, -0.0035, 0.0072,
               0.3160, 0.3976, 0.4487, 0.1097, 0.0369, -0.0537, -0.0093,
               -0.1286, -0.0215, -0.0803, -0.3203, -0.2066, -0.0500, -0.2825,
               -0.5325, -0.2987, 0.2981, -0.2310, 0.8854, 0.5082, -0.0064,
               0.2318, 0, 0, -0.1222, 0, -0.4625, 0, 0.9198, 0, 0, 1.8957,
               1.3434, 0.1815, 0.3401, 0.5794, 1.2171, 0.2427, 0.1399, 0.0591,
               -0.3192, -0.3486, -0.0268, 0.0653, -0.0430, 0.1102, 0.0667,
               0.4537, 0.2590, -0.0425, 0.1005, 1.0587, -0.0177, 0.1702,
               0.1584, 0.3094, 0.1030, 0.2238, 0.5390, -0.2197, 0.0886, 0.0352,
               0, 0.0196, 0.0494, -1.5974, 0.4267, 0.0879, 0.6115, 0.0173,
               -0.0504, -0.2474, -0.5736, -0.1210, -0.0148, 0.1395, 0.1829,
               -0.1192, -0.0455, 0.2667, -0.1899, -0.3179, -0.3576, -0.7458,
               0.1218, -0.0569, 0.4649, 0.1506, 0.1300, 0.6540, 0.0043,
               -0.2692, -0.2787, -0.2107, 0.0926, -0.4495, 0.0722, 0.0319,
               -0.6775, 0.0604, -0.3080, -0.1590, 0.0217, 0.1007, -0.1647,
               -0.1387, -0.1314, 0.2745, 0.1645, 0.0754, -0.1196, 0.0494,
               0.1344, 0.0032, -0.0817, -0.1564, -0.5176, 0.5477, 0.3533,
               -0.3888, 1.6498, 0, 1.0750, 0.7193, 1.1867, 0, 0.2991, -0.4605,
               0.0060, -0.1819, 0.4516, 1.3440, 0, -0.3741, -0.4961, -0.4574,
               -0.1736, 0.3138, 0.5928, 0.4387, 0.5632, 0.1415, 0, 0.2391,
               0.7192, 0, 0, 1.0171, 0.9674, 0.1126, 0.9317, 0.5031, 0, 0, 0,
               0, 0, 0.2242, 0.0185, 0, -0.0850, 1.1457, 0.5768, -0.5335,
               -5.2736, 0, 0.6571, -0.8252, 0.2790, -0.6848, 0.0441, -0.1666,
               -0.2692, -0.2807, -0.3294, -0.2931, -0.3360, 0.0402, 1.0466,
               -7.8521, 0, 0.9126, 0, -0.9432, -0.5844, 0.1733],
        "tc": [1.7506, 1.3327, 0.5960, 0.0306, 3.2295, 3.0741, 2.7717, 2.5666,
               2.6391, 5.4330, 4.8219, 0, 3.7897, 4.5870, 2.0337, 5.4979,
               3.1058, 4.5344, 4.0954, 3.4611, 2.9003, 1.9512, 0.8576, 5.7861,
               6.5062, 4.9967, 6.4572, 0, 5.2188, 9.3472, 14.6038, 15.4515,
               7.0058, 5.7157, 4.4743, 0, 9.4806, 5.8013, 9.4795, 6.3179,
               5.9619, 4.7558, 0, 5.6064, 6.7311, 0, 0, 4.7346, 3.4393, 2.4217,
               0.7889, 0.2511, 3.6588, 8.1745, 4.2847, 2.8546, 4.5529, 3.2422,
               2.0057, 3.0106, 2.1673, 10.2155, 8.4081, 5.8536, 4.7420, 0, 0,
               12.9827, 8.4309, 5.8829, 10.4124, 8.1381, 0, 0, 0, 6.5884,
               10.9507, 9.5487, 0, 12.1243, 0, 0, 0, 0, 0, 25.1184, 20.5590, 0,
               15.4603, 0, 0, 0, 0, 0, 19.8979, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 6.2561, 4.3756, 3.7063, 7.8956, 0, 8.8073, 3.3179, 2.6702,
               2.1633, 3.5702, 0.8543, 1.7737, 5.1653, 0, 3.0593, 5.7046,
               1.5491, 12.4470, 8.4199, 8.5775, 4.5036, 0.8976, 4.0947, 0, 0,
               0, 10.4579, 0, 0, 5.8307, 7.7300, 5.8527, 4.6431, 9.5115,
               7.7987, 6.9733, 6.4871, 0, 0, 0, 19.8953, 17.2586, 8.6910, 0,
               18.9366, 0, 0, 0, 0, 0, 0, 0, 0, 16.9914, 0, 38.6148, 6.6804,
               6.6418, 6.0159, 0, 1.8815, 1.1020, -0.2399, 3.6426, 3.5475, 0,
               4.4913, 5.9726, 4.3905, 0, 0, 2.7409, 12.6396, 5.5523, 24.3995,
               -0.0471, -0.1778, 0.5602, 0.8994, 1.5535, 0.4214, -0.0172,
               0.0262, -0.1526, -1.0434, -0.0338, -0.3658, -4.7275, -0.5537,
               -0.3576, -0.6768, -1.5224, -0.3940, 0.3414, 0, 1.9395, 1.2342,
               -3.3555, -1.1598, 0, 62.4740, 1.9595, 0.7686, 0, 0, 0.4750, 0,
               -0.2697, 1.9699, 0.2311, 2.1272, 1.5418, 0, 2.7051, 0.7502,
               0.2900, 0, 0, 0, -0.0188, -1.1249, -0.0880, 0, -1.7762, 2.2630,
               1.4069, 0.2698, -1.0107, 2.4950, 0, 0.8705, 0, 0, 0, 1.7860, 0,
               0, 1.1629, 0, 0, 0.1565, 0.8016, 0, 2.4070, 0, 0, -0.2509,
               -1.1019, 0, 0, -0.1233, 0.3816, 0.1093, 0, 0, -0.2832, 0, 0,
               0.8973, -0.9610, -2.0833, 0, -0.6447, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               -0.3450, 0.1607, 0.1090, -2.1303, -0.0003, 0, -0.3161, -0.0693,
               0.0803, 1.0088, 0.0908, -0.6412, 2.1116, 0.9353, 0.6241,
               -1.0256, 0.5784, 0.6595, 0, 0, 0, -2.2773, 0, 0, 0, -1.6986, 0,
               0.4950, 0.1725, 6.6872, 2.6769, 0, 0, -4.3195, 0, 0, 0.0834, 0,
               -0.0185, 6.1894, -0.2474, 0.5060, 3.0321, 0, 2.1761, 0, 0.4963,
               0, 0.1174, 0.7039, 0, 0, -0.2678, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               3.5541, 0, 0, 0, 3.5541, 0.9519, 0, 0, 0, 0, -0.9785, -1.4002,
               0, -1.0095, 0.1605, -0.6765, 0, 0, 0, 0, -1.0430, 3.3011, 0, 0,
               2.8885, 0, 1.1251, 3.9241, 7.7134],
        "pc": [0.018615, 0.013547, 0.007259, 0.001219, 0.025745, 0.023003,
               0.021137, 0.019609, 0.014114, 0.035483, 0.029678, 0, 0.014010,
               0.010888, 0.007260, 0.003564, 0.006512, 0.012859, -0.003339,
               0.020907, 0.018082, 0.011795, 0.011298, 0.030637, 0.026282,
               0.026371, 0.019507, 0, -0.005401, -0.008788, 0.009885, 0.017100,
               0.025227, 0.019619, 0.012487, 0, 0.011007, 0.010204, 0.019633,
               0.033812, 0.026983, 0.021990, 0, 0.015249, 0.018948, 0, 0,
               0.013087, 0.020084, 0.017954, 0.014487, 0.005613, 0.005115,
               0.011413, 0.013049, 0.010790, 0.015863, 0.020482, 0.005329,
               0.021186, 0.027454, 0.005335, -0.005596, -0.000838, 0.000571,
               0, 0, 0.036523, 0.029034, 0.024654, 0.020978, 0.024346, 0, 0, 0,
               0.025065, 0.021056, 0.014899, 0, 0.018311, 0, 0, 0, 0, 0,
               0.001467, 0.023455, 0, 0.043090, 0, 0, 0, 0, 0, 0.023447, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.021419, 0.015640, 0.009187,
               0.028236, 0, 0.036746, 0.023315, 0.020040, -0.010120, 0.031524,
               0.018572, 0.048565, 0.037948, 0, 0.041641, 0.016033, 0.014037,
               0.014403, 0.010199, -0.004637, -0.001460, 0.012034, 0.007923,
               0, 0, 0, 0.025986, 0, 0.018783, -0.002815, 0.017299, 0.008968,
               0.005118, 0.010086, 0.006399, 0.018013, 0.015254, 0, 0, 0,
               -0.005534, -0.000784, 0.004240, 0, -0.027208, 0, 0, 0, 0, 0, 0,
               0, 0, -0.029036, 0, -0.126108, 0.007235, 0.021238, 0.010678, 0,
               0.009884, 0.007596, 0.003268, 0.013815, 0.010576, 0, 0.019101,
               -0.003678, -0.001179, 0, 0, -0.000387, -0.000207, 0.001540,
               0.002487, 0.000473, 0.000340, -0.003207, -0.008733, -0.016852,
               0.000792, -0.000101, 0.000815, -0.000163, 0.005789, -0.000111,
               -0.001892, 0.006916, -0.000569, 0.001812, 0.000246, 0.003224,
               -0.002912, -0.000516, 0, -0.004712, 0.002581, 0.000726,
               0.000157, 0, -0.002696, -0.001479, 0.000090, 0, 0, -0.001445, 0,
               0.000655, 0.002330, -0.001022, 0.001321, -0.003385, 0,
               -0.001999, -0.000231, -0.000432, 0, -0.010021, 0, 0.000152,
               0.000893, 0.000044, 0, -0.000763, 0.002464, -0.000034,
               -0.000417, 0.002944, -0.000796, 0, 0.000183, 0, 0, 0, 0.004195,
               0, 0, -0.000384, 0, 0, -0.001446, -0.006495, 0, -0.002650, 0, 0,
               -0.000624, 0.003921, 0, 0, 0.000779, 0.001694, 0.000124, 0, 0,
               0.002114, 0, 0, 0.004640, 0.000039, -0.014535, 0, -0.000199, 0,
               0, 0, 0, 0, 0, 0, 0, 0, -0.000692, 0.001235, -0.000610,
               -0.004683, 0.000058, 0, 0.000522, 0.001790, 0.000467, -0.005598,
               0.000255, 0.004090, -0.007612, -0.001811, -0.000500, 0.007006,
               0.007006, 0.001283, 0, 0, 0, 0.008029, 0, 0, 0, 0.001544, 0,
               0.000728, -0.000327, 0.001937, 0.003792, 0.003254, 0, 0.006734,
               0, 0, -0.011090, 0, 0.000851, -0.040100, -0.005826, -0.003746,
               0.003007, 0, 0.002745, 0, -0.000985, 0, -0.002673, -0.004661, 0,
               0, -0.001837, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.004600, 0, 0, 0,
               0.017201, 0.008484, 0, 0, 0, 0, 0.001162, -0.004716, 0,
               -0.001332, -0.002030, -0.002436, 0, 0, 0, 0, 0.004695, 0.015244,
               0, 0, 0.007280, 0, -0.005369, -0.011207, -0.001275],
        "vc": [68.35, 56.28, 37.50, 16.01, 111.43, 98.43, 91.40, 83.89, 90.66,
               143.57, 146.36, 0, 84.60, 74.66, 42.39, 35.71, 34.65, 26.47,
               36.47, 97.33, 87.19, 73.51, 67.20, 134.69, 128.84, 110.74,
               112.08, 0, 30.61, 50.77, 90.66, 119.10, 127.99, 112.79, 97.16,
               0, 90.69, 71.08, 122.91, 148.91, 132.89, 125.52, 0, 93.29,
               105.53, 0, 0, 81.17, 88.20, 74.03, 60.06, 52.96, 47.27, 117.62,
               76.36, 80.01, 77.04, 95.15, 99.16, 94.94, 74.05, 81.40, 86.37,
               108.39, 63.39, 0, 0, 133.62, 134.73, 120.74, 119.08, 94.91, 0,
               0, 0, 141.24, 157.57, 143.36, 0, 133.06, 0, 0, 0, 0, 0, 138.71,
               190.71, 0, 244.71, 0, 0, 0, 0, 0, 162.08, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 112.12, 100.78, 87.01, 159.79, 0, 204.71, 87.71,
               78.08, 0, 102.71, 95.09, 108.85, 171.04, 0, 146.01, 92.67,
               54.36, 131.08, 104.12, 104.28, 77.99, 24.62, 57.77, 0, 0, 0,
               159.33, 0, 147.66, 58.01, 105.68, 109.36, 94.01, 95.08, 57.89,
               122.03, 106.60, 0, 0, 0, 82.36, 89.95, 115.80, 0, 144.58, 0,
               135.47, 0, 0, 0, 0, 0, 0, 85.59, 0, -142.79, 93.56, 125.43,
               194.36, 0, 49.24, 44.95, 33.32, 83.91, 70.98, 0, 83.96, 51.80,
               31.41, 0, 0, 17.69, 57.38, 45.45, 96.66, 1.71, 3.14, -3.75,
               -10.06, -8.70, -7.88, 0.50, 0.14, -2.67, 10.36, -4.08, 3.02,
               10.56, 4.28, 2.98, -3.04, 13.98, 5.17, 0.68, 0, 7.54, 5.58,
               20.82, -26.31, 0, 17.78, 12.46, 15.17, 0, 0, 7.91, 0, 20.43,
               24.82, 14.54, -10.31, -2.33, 0, -0.73, 1.69, -4.54, 0, 2.63, 0,
               2.80, 3.82, 0.21, 0, 4.36, -4.88, 2.50, -7.49, -0.25, -11.01, 0,
               2.00, 0, 0, 0, -3.40, 0, 0, -7.02, 0, 0, -2.04, -5.70, 0, 0.39,
               0, 0, 0.03, -4.43, 0, 0, 2.79, -2.95, 6.19, 0, 0, -16.97, 0, 0,
               -7.73, -2.50, -51.50, 0, -2.00, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               -12.03, 1.95, -5.17, -14.40, 0, 0, 2.86, 6.54, 3.70, -9.58,
               -2.05, -7.67, -7.04, -0.04, -0.04, 8.68, 8.68, 14.28, 0, 0, 0,
               -50.26, 0, 0, 0, -3.72, 0, -23.74, -0.84, 1.44, -1.62, -0.69, 0,
               6.69, 0, 0, -36.89, 0, -8.87, -26.26, -8.55, -11.56, 9.73, 0,
               7.72, 0, -3.33, 0, -4.67, 14.31, 0, 0, -7.05, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 12.60, 0, 0, 0, 0.44, 1.42, 0, 0, 0, 0, -2.63, 28.42,
               0, -6.88, -3.17, -3.85, 0, 0, 0, 0, 35.21, -6.96, 0, 0, -24.02,
               0, 63.29, -2.71, -12.04],
        "gf": [2.878, 8.064, 8.254, 16.413, 95.738, 92.656, 85.107, 88.691,
               93.119, 229.906, 226.710, 0, 230.029, 216.013, 26.732, 20.379,
               33.912, 23.331, 89.902, 24.919, 31.663, 30.393, 40.127, 114.531,
               111.216, 115.728, 263.205, 0, -144.051, -131.327, -337.090,
               -312.422, -120.667, -120.425, -116.799, 0, -91.812, -100.882,
               -80.222, -306.733, -298.332, -301.414, 0, -276.878, -291.662, 0,
               0, -299.803, -90.329, -105.579, -101.207, -92.804, -83.354,
               68.812, 61.452, 55.202, 88.512, 88.874, 73.101, 125.906,
               121.247, 66.470, 98.195, 143.280, 42.687, 0, 0, 134.997,
               142.475, 142.295, 162.175, 130.986, 0, 0, 0, 0, 25.783, 16.407,
               0, 57.352, 0, 0, 0, 0, 0, -127.512, -102.912, 0, -56.412, 0, 0,
               0, 0, 0, -44.595, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -19.484,
               -31.933, -37.848, -24.214, 0, -44.122, -180.212, -228.239, 0,
               -411.239, 0, -615.333, -249.020, 0, -396.814, 1.985, -141.306,
               91.505, 42.977, 43.910, 5.528, -182.973, -29.876, 0, 0, 0,
               -233.335, 0, -239.423, -75.568, 27.469, 27.030, 27.338, 48.905,
               15.818, 35.845, 42.684, 0, 0, 0, -52.231, -257.608, 0, 0,
               -519.853, 0, -314.643, 0, 0, 0, 0, 0, 0, 0, 0, 0, -447.186,
               11.149, 1.890, 0, 13.287, 6.107, -0.193, 86.493, 67.056, 0, 0,
               72.540, 83.779, 0, 0, -114.062, -156.672, 12.020, -241.601,
               -0.418, -2.776, 6.996, 8.938, 10.735, -6.562, -0.120, 1.006,
               3.857, -0.525, -1.543, 2.202, 3.920, -11.779, -16.075, -5.614,
               -25.382, 6.621, 4.833, 0, -1.051, -1.506, 0.344, 3.848, 0,
               3.145, -5.217, -4.281, 0, 0, -2.678, 0, -7.376, 18.974, 0.558,
               6.728, 1.306, 0, 0, 10.556, -10.098, 0, 14.470, 0, 8.207,
               -8.304, -12.085, 0, 10.194, -8.081, -2.044, 6.043, 0, 25.157,
               0, 16.725, 0, 0, 0, 3.020, 0, 0, 1.556, 0, 0, 1.238, 0.354, 0,
               -6.438, 0, 0, 28.972, -22.533, 0, 0, 4.178, 5.332, 6.084, 0, 0,
               6.768, 0, 0, -3.024, 2.046, -11.965, 0, -16.723, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 4.358, 0.107, 18.755, -18.970, 0, 0, 1.577, -1.037,
               -0.709, 7.731, -2.767, -2.148, 14.226, 4.926, -0.474, -9.713,
               -2.523, -4.703, 0, 0, 0, -16.570, 0, 0, 0, -4.708, 0, 3.079,
               7.536, -8.397, 10.194, 1.662, 0, 4.100, 0, 0, -7.035, 0, -1.601,
               -4.459, -5.267, -4.203, 1.318, 0, -67.517, 0, 0, 0, -0.729,
               -0.702, 0, 0, 11.125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.333, 0, 0, 0,
               -2.221, -0.596, 0, 0, 0, 0, 2.668, -4.229, 0, 1.993, -2.940,
               -1.873, 0, 0, 0, 0, 3.896, 13.843, 0, 0, -16.040, 0, 8.688,
               -5.112, 20.073],
        "hf": [-42.479, -20.829, -7.122, 8.928, 57.509, 69.664, 61.625, 81.835,
               95.710, 198.840, 208.490, 0, 224.902, 228.282, 12.861, 20.187,
               30.768, 24.701, 70.862, -19.258, 4.380, 18.440, 35.297, 77.863,
               88.084, 90.927, 257.448, 0, -178.360, -164.191, -389.931,
               -361.249, -180.604, -163.090, -139.909, 0, -106.965, -130.816,
               -107.159, -387.458, -364.204, -352.057, 0, -327.678, -307.727,
               0, 0, -331.397, -156.062, -152.239, -147.709, -121.608,
               -101.783, -10.703, 0.730, 2.019, 24.740, 23.610, 21.491, 55.024,
               65.331, 17.501, 53.274, 115.606, -8.556, 0, 0, 99.245, 151.390,
               124.770, 148.968, 124.917, 0, 0, 0, 0, -65.620, -60.750, 0,
               -22.931, 0, 0, 0, 0, 0, -201.369, -203.069, -183.613, -188.069,
               -48.210, 0, 0, 0, 0, -125.052, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, -65.056, -65.127, -62.881, -80.812, 0, -105.369, -227.469,
               -261.901, 0, -463.901, 0, -673.875, -306.765, 0, -458.074,
               -17.002, -160.965, 95.048, 38.917, 47.632, -1.703, -201.968,
               -46.963, 0, 0, 0, -343.903, 0, -333.385, -125.111, -8.021,
               3.510, 12.589, 41.648, 11.339, -3.337, 21.492, 0, 0, 0, -71.050,
               -305.498, -430.833, 0, -621.412, 0, -370.493, 0, 0, 0, 0, 0, 0,
               -1060.325, -1005.161, 72.339, -516.282, -52.241, -51.390, 0,
               -18.575, -12.464, -2.098, 59.841, 64.295, 0, 0, 23.138, 65.622,
               0, 0, -137.353, -180.166, 15.453, -283.839, -0.419, -1.967,
               6.065, 8.078, 10.535, -11.786, -0.048, 1.449, 3.964, 1.514,
               0.033, 4.994, 1.121, -12.295, -14.140, -4.422, -25.929, 8.244,
               0.000, 0, -0.592, -0.959, -1.443, 3.608, 0, 6.598, -6.058,
               -6.929, 0, 0, -1.727, 0, 7.292, 5.661, -3.906, 0.794, 4.025, 0,
               0, -7.261, -9.411, 0, 17.014, 0, 9.715, -16.903, -12.509, 0,
               9.090, -8.570, -3.447, 5.486, 0, 16.950, 0, 7.568, 0, 0, 0,
               4.145, 0, 0, 4.066, 0, 0, -0.751, -0.192, 0, -12.517, 0, 0,
               24.560, -12.044, 0, 0, 4.452, 4.428, -4.128, 0, 0, 10.390, 0, 0,
               -8.050, 3.446, 14.531, 0, -7.569, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               -15.751, 0.238, 21.498, -21.975, 0, 0, 1.486, 0.294, 0.384,
               5.743, -0.449, -7.538, 12.710, 5.220, -1.340, -9.644, -2.446,
               -6.466, 0, 0, 0, -17.778, 0, 0, 0, -6.572, 0, 4.171, 5.411,
               -8.651, 8.164, 1.753, 0, 0.371, 0, 0, 7.782, 0, 2.689, -4.558,
               -5.914, -4.863, 0.084, 0, -66.870, 0, 0, 0, 0.866, -2.291, 0, 0,
               7.108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.212, 0, 0, 0, -16.080,
               -1.994, 0, 0, 0, 0, -5.074, -2.303, 0, 1.904, -2.274, -1.316, 0,
               0, 0, 0, 5.819, 11.387, 0, 0, -19.089, 0, 13.586, -0.314,
               15.786],
        "hv": [0.217, 4.910, 7.962, 10.730, 4.031, 9.456, 8.602, 14.095,
               19.910, 11.310, 0, 0, 6.144, 12.540, 3.683, 6.631, 6.152, 6.824,
               9.420, 8.279, 11.981, 13.519, 16.912, 0, 0, 0, 0, 0, 24.214,
               34.099, 17.002, 0, 15.195, 19.392, 20.350, 0, 25.036, 12.370, 0,
               19.342, 21.100, 24.937, 23.739, 15.422, 25.206, 0, 0, 0, 5.783,
               9.997, 14.620, 13.850, 16.151, 15.432, 16.048, 17.257, 11.831,
               13.067, 14.048, 9.493, 12.636, 23.335, 23.026, 22.249, 13.425,
               0, 0, 21.923, 24.963, 24.967, 0, 16.639, 0, 0, 0, 0, 29.640,
               29.173, 0, 24.863, 0, 0, 0, 0, 46.490, 44.240, 0, 52.723,
               38.290, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               11.754, 12.048, 16.597, 17.251, 20.473, 20.550, 8.238, 0, 6.739,
               0, 1.621, 7.352, 8.630, 0, 8.086, 11.224, 3.965, 0, 14.393,
               14.171, 9.888, 0, 0, 0, 0, 0, 31.493, 0, 0, 0, 16.815, 17.098,
               18.397, 17.413, 9.813, 14.296, 16.965, 19.038, 19.996, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21.613, 0, 0, 0,
               3.341, 6.416, 7.017, 7.767, 7.171, 0, 0, 13.700, 0, 0, 0, 6.877,
               17.124, 12.262, 0, -0.399, -0.417, 0.532, 0.623, 5.086, 1.632,
               0.064, -0.060, 0.004, -0.550, -0.403, 0.723, 7.422, -1.871, 0,
               -0.206, -1.579, 0, 0, 0, -6.611, 0, 2.384, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, -0.683, 1.203, 0, 0, 0, 0.372, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 4.608, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.030,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0.096, -0.428, 0.153, 0, 0, 0, 0, 0,
               2.134, -4.607, 0, 0, 0, 0, 0, -0.616, 0, 0, 0, 0, 0, 0, 0,
               0.808, 0.585, 0, 0, 0, 1.164, -1.910, 0.331, 1.433, 0.313,
               -0.117, 0, 0, 0, -1.683, 0.277, 0.397, -0.939, -1.269, -1.719,
               -3.419, 1.742, 0.572, -2.744, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "hm": [1.660, 2.639, 0.134, -1.232, 1.268, 4.441, 2.451, 3.032, 2.616,
               7.076, 7.435, 6.000, -1.548, 6.128, 1.948, 0.845, 1.095, -0.531,
               2.555, 2.969, 0.948, -1.037, -2.856, 4.013, 8.274, 3.324, 2.514,
               0, 4.786, 8.427, 10.692, 14.649, 8.062, 8.826, 7.205, 0, 4.852,
               11.325, 7.273, 7.910, 9.479, 9.317, 0, 8.115, 8.149, 0, 5.875,
               10.573, 5.089, 4.891, 4.766, 2.458, -0.118, 13.482, 6.283, 0,
               4.490, 7.711, 2.561, 6.008, 1.756, 6.542, 0.624, -2.576, 6.158,
               0, 0, 7.303, 9.464, 4.166, 6.788, 6.867, 0, 0, 0, 0, 10.989, 0,
               -4.187, 7.572, 6.302, 0, 9.353, 0, 0, 16.840, 17.429, 0, 11.553,
               0, 0, 0, 0, 16.811, 8.658, 0, 10.959, 4.370, 0, 9.862, 12.845,
               10.958, 12.098, 9.557, 16.703, 18.460, 0, 6.353, 0, -0.082,
               6.781, 1.823, 3.492, 7.139, 3.917, 0, 7.011, 0, 2.526, 3.114, 0,
               2.156, 4.435, 2.003, 2.814, 5.734, 6.103, 4.826, 3.096, 5.181,
               0, 0, 0, 8.454, 0, 12.594, 0, 10.068, 4.266, -0.623, 4.513,
               5.829, 7.497, 4.096, 0, 0, 0, 13.403, 17.748, 0, 0, 0, 0, 3.281,
               0, 0, 0, 0, 0, 0, 0, 4.256, -5.654, 8.363, 0, 0, 0, 1.069,
               2.511, -0.921, 1.185, 2.559, 0, 5.351, 8.655, 0.269, 3.993, 0,
               3.806, 6.137, 5.170, 9.934, 0.396, 0.554, -1.766, 0.351, -1.089,
               1.408, 0.070, -0.632, -0.368, -0.369, 0.105, 1.005, 5.475,
               1.208, -2.666, -0.599, -0.459, 0, -0.149, 0, -0.306, -0.041,
               -1.575, 0, 0, 7.032, 4.264, 29.245, 0, 0, 0, 0, -3.623, -8.038,
               -4.371, -0.931, 0, 2.303, 1.100, 0, 3.169, 2.823, 2.212, 0,
               -0.480, -0.405, -0.014, 0, 1.291, 0, -0.639, 0.969, -2.754, 0,
               0, 0.890, -4.086, 0, 0, 0, 0, 0, 0, 0, 0, -0.270, -0.878, 0,
               -1.670, 0, 0, 2.235, 0.961, 0, 0, 0.033, -1.137, 2.421, 0, 0, 0,
               0, 0, 0, 0.328, 0.402, 0, -0.878, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               -1.237, 0, 0.235, 0, 0, -1.470, -1.059, 1.244, 0.473, -0.302,
               -2.530, -1.736, -2.246, 8.034, -0.786, 3.671, 5.975, 0, 0, 0,
               -1.487, 0, 0, 0, -7.583, 0, -4.840, -0.272, 1.661, -3.479,
               0.301, 0, 5.666, 0, 0, -0.607, 0, -2.703, -0.385, -0.442,
               -0.143, 5.377, 0, 0, 0, 0.223, 0, -0.958, 3.275, 0, 0, -4.091,
               0, 0, -3.687, 2.047, 7.327, 0, 2.510, 0.018, 0, -0.784, 0, 0,
               -2.485, 0.196, 1.934, 0, -0.599, 0, 0, 1.193, -3.971, -1.153,
               0.694, -3.699, 2.037, 2.150, 0, 0, 0, 1.176, 5.027, 0, 0,
               -3.417, 0, -4.967, -2.587, -1.365],

        "txt": [
            ("CH3-",),                      # 0
            ("-CH2-",),
            (">CH-",),
            (">C<",),
            ("CH2=CH-",),
            ("-CH=CH-",),
            ("CH2=C<",),
            ("-CH=C<",),
            (">C=C<",),
            ("CH2=C=CH-",),
            ("CH2=C=C<",),                  # 10
            ("-CH=C=CH-",),
            ("CH≡C-",),
            ("-C≡C-",),
            ("(a)CH",),
            ("(a)C (fused)",),
            ("(a)C (ring fused)",),
            ("(a)C (others)",),
            ("(a)N-",),
            ("(a)C-CH3",),
            ("(a)C-CH2-",),                 # 20
            ("(a)C-CH<",),
            ("(a)C-C<",),
            ("(a)C-CH=CH2",),
            ("(a)C-CH=CH-",),
            ("(a)C-C=CH2",),
            ("(a)C-C≡CH",),
            ("(a)C-C≡C-",),
            ("-OH",),
            ("(a)C-OH",),
            ("-COOH",),                     # 30
            ("(a)C-COOH",),
            ("CH3-CO-",),
            ("-CH2-CO-",),
            (">CH-CO-",),
            (">C-CO-",),
            ("(a)C-CO-",),
            ("-CHO",),
            ("(a)C-CHO",),
            ("CH3-COO-",),
            ("-CH2-COO-",),                 # 40
            (">CH-COO-",),
            (">C-COO-",),
            ("HCOO-",),
            ("(a)C-COO-",),
            ("(a)C-OOCH",),
            ("(a)C-OOC-",),
            ("COO- (others)",),
            ("CH3-O-",),
            ("-CH2-O-",),
            (">CH-O-",),                    # 50
            (">C-O-",),
            ("(a)C-O-",),
            ("-CH2-NH2",),
            (">CH-NH2",),
            (">C-NH2",),
            ("CH3-NH-",),
            ("-CH2-NH-",),
            (">CH-NH-",),
            ("CH3-N<",),
            ("-CH2-N<",),                   # 60
            ("(a)C-NH2",),
            ("(a)C-NH-",),
            ("(a)C-N<",),
            ("NH2 (others)",),
            ("-CH=N-",),
            (">C=N-",),
            ("-CH2-CN",),
            (">CH-CN",),
            (">C-CN",),
            ("(a)C-CN",),                   # 70
            ("CN (others)",),
            ("-CH2-NCO",),
            (">CH-NCO",),
            (">C-NCO",),
            ("(a)C-NCO",),
            ("-CH2-NO2",),
            (">CH-NO2",),
            (">C-NO2",),
            ("(a)C-NO2",),
            ("NO2 (others)",),              # 80
            ("-ONO",),
            ("-ONO2",),
            ("-HCON-CH2CH2",),
            ("-HCONH-CH2",),
            ("-CO-NH2",),
            ("-CO-NH-CH3",),
            ("-CO-NH-CH2-",),
            ("-CO-N-CH3CH3",),
            ("-CO-NCH3-CH2-",),
            ("-CO-N-CH2CH2<",),            # 90
            ("-CO-NH-CO-",),
            ("-CO-N-CO-",),
            ("(a)C-CO-NH2",),
            ("(a)C-NHCOH",),
            ("(a)C-N(CO)H-",),
            ("(a)C-CONH-",),
            ("(a)C-NHCO-",),
            ("(a)C-NCO-",),
            ("-NH-CO-NH-",),
            ("NH2-CO-NH-",),               # 100
            ("NH2-CO-N<",),
            ("-NH-CO-N<",),
            (">N-CO-N<",),
            ("(a)C-NHCONH2",),
            ("(a)C-NHCONH-",),
            ("-NHCO- (others)",),
            ("-CH2Cl",),
            (">CHCl",),
            (">CCl",),
            ("-CHCl2",),                   # 110
            (">CCl2",),
            ("-CCl3",),
            ("-CH2F",),
            (">CHF",),
            (">CF",),
            ("-CHF2",),
            (">CF2",),
            (">CF3",),
            ("-CCl2F",),
            ("-CHClF",),                   # 120
            ("-CClF2",),
            ("(a)C-Cl",),
            ("(a)C-F",),
            ("(a)C-I",),
            ("(a)C-Br",),
            ("I- (others)",),
            ("Br- (others)",),
            ("F- (others)",),
            ("Cl- (others)",),
            ("-CHNOH",),                   # 130
            (">CNOH",),
            ("(a)C-CHNOH",),
            ("-O-CH2-CH2-OH",),
            ("-O-CH-CH2-OH",),
            ("-O-CH2-CH-OH",),
            ("-O-OH",),
            ("-CH2-SH",),
            (">CH-SH",),
            (">C-SH",),
            ("(a)C-SH",),                  # 140
            ("-SH (others)",),
            ("CH3-S-",),
            ("-CH2-S-",),
            (">CH-S-",),
            (">C-S-",),
            ("(a)C-S-",),
            ("-SO-",),
            ("-SO2",),
            ("-SO3 (sulfite)",),
            ("-SO3 (sulfonate)",),         # 150
            ("-SO4 (sulfate)",),
            ("(a)C-SO-",),
            ("(a)C-SO2",),
            ("PH (phosphine)",),
            ("P (phosphine)",),
            ("PO3 (phosphite)",),
            ("PHO3 (phosphonate)",),
            ("PO3 (phosphonate)",),
            ("PHO4 (phosphate)",),
            ("PO4 (phosphate)",),          # 160
            ("(a)C-PO4",),
            ("(a)C-P",),
            ("CO3 (carbonate)",),
            ("C2H3O (oxirane)",),
            ("C2H2O (oxirane)",),
            ("C2O (oxirane)",),
            ("-CH2- (cyclic)",),
            ("-CH< (cyclic)",),
            (">C< (cyclic)",),
            ("-CH=CH- (cyclic)",),         # 170
            ("-CH=C< (cyclic)",),
            (">C=C< (cyclic)",),
            ("CH2=C< (cyclic)",),
            ("-NH- (cyclic)",),
            (">N- (cyclic)",),
            ("-CH=N- (cyclic)",),
            (">C=N- (cyclic)",),
            ("-O- (cyclic)",),
            ("-CO- (cyclic)",),
            ("-S- (cyclic)",),             # 180
            ("-SO2- (cyclic)",),

            # 2nd order term, index here are 181+, 1 -> 182
            ("CH(CH3)2", ),
            ("C(CH3)3", ),
            ("CHCH3CHCH3", ),
            ("CH(CH3)C(CH3)2", ),
            ("C(CH3)2C(CH3)2", ),
            ("CHn=CHm-CHp=CHk (m, p (0,1); k, n (0,2))", ),
            ("CH3-CHm=CHn (m (0,1); n (0,2))", ),
            ("CH2-CHm=CHn (m (0,1); n (0,2))", ),
            ("CHp-CHm=CHn (m, p (0,1); n (0,2))", ),          # 190
            ("CHCHO or CCHO", ),
            ("CH3COCH2", ),
            ("CH3COCH or CH3COC", ),
            ("CHCOOH or CCOOH", ),
            ("CH3COOCH or CH3COOC", ),
            ("CO-O-CO", ),
            ("CHOH", ),
            ("COH", ),
            ("CH3(OH)CHn(OH) (n (0,2))", ),
            ("NCCHOH or NCCOH", ),                            # 200
            ("OH-CHn-COO (n (0,2))", ),
            ("CHm(OH)CHn(OH) (0<m,n<2)", ),
            ("CHm(OH)CHn(-) (0<m,n,p<2)", ),
            ("CHm(NH2)CHn(NH2) (0<m,n<2)", ),
            ("CHm(NH)CHn(NH2) (1<m,n<2)", ),
            ("H2NCOCHnCHmCONH2 (1<m,n<2)", ),
            ("CHm(NHn)-COOH (0<m,n<2)", ),
            ("HOOC-CHn-COOH (1<n<2)", ),
            ("HOOC-CHn-CHm-COOH (1<n,m<2)", ),
            ("HO-CHn-COOH (1<n<2)", ),                        # 210
            ("NH2-CHn-CHm-COOH (1<n,m<2)", ),
            ("CH3-O-CHn-COOH (1<n<2)", ),
            ("SH-CH-COOH", ),
            ("SH-CHn-CHm-COOH (1<n,m<2)", ),
            ("CN-CHn-CHm-CN (1<n,m<2)", ),
            ("OH-CHn-CHm-CN (1<n,m<2)", ),
            ("SH-CHn-CHm-SH (1<n,m<2)", ),
            ("COO-CHn-CHm-OOC (1<n,m<2)", ),
            ("OOC-CHn-CHm-COO (1<n,m<2)", ),
            ("NC-CHn-COO (1<n<2)", ),                         # 220
            ("CO-CHn-COO (1<n<2)", ),
            ("CHm-O-CHn=CHp (0<n,m,p<3)", ),
            ("CHm=CHn-F (0<n,m<2)", ),
            ("CHm=CHn-Br (0<n,m<2)", ),
            ("CHm=CHn-I (0<n,m<2)", ),
            ("CHm=CHn-Cl (0<n,m<2)", ),
            ("CHm=CHn-CN (0<n,m<2)", ),
            ("CHn=CHm-COO-CHp (0<n,m,p<3)", ),
            ("CHm=CHn-CHO (0<n,m<2)", ),
            ("CHm=CHn-COOH (0<n,m<2)", ),                     # 230
            ("aC-CHn-X (1<n<2; X:halogen)", ),
            ("aC-CHn-NHm (1<n<2; 0<m<2)", ),
            ("aC-CHn-O- (1<n<2)", ),
            ("aC-CHn-OH (1<n<2)", ),
            ("aC-CHn-CN (1<n<2)", ),
            ("aC-CHn-CHO (1<n<2)", ),
            ("aC-CHn-SH (1<n<2)", ),
            ("aC-CHn-COOH (1<n<2)", ),
            ("aC-CHn-CO- (1<n<2)", ),
            ("aC-CHn-S- (1<n<2)", ),                          # 240
            ("aC-CHn-OOCH (1<n<2)", ),
            ("aC-CHn-NO2 (1<n<2)", ),
            ("aC-CHn-CONH2 (1<n<2)", ),
            ("aC-CHn-OOC (1<n<2)", ),
            ("aC-CHn-COO (1<n<2)", ),
            ("aC-SO2-OH", ),
            ("aC-CH(CH3)2", ),
            ("aC-C(CH3)3", ),
            ("aC-CF3", ),
            ("(CHn=C)(cyclic)-CHO (0<n<2)", ),                # 250
            ("(CHn=C)(cyclic)-COO-CHm (0<n,m<3)", ),
            ("(CHn=C)(cyclic)-CO- (0<n<2)", ),
            ("(CHn=C)(cyclic)-CH3 (0<n<2)", ),
            ("(CHn=C)(cyclic)-CH2- (0<n<2)", ),
            ("(CHn=C)(cyclic)-CN (0<n<2)", ),
            ("(CHn=C)(cyclic)-Cl (0<n<2)", ),
            ("CH(cyclic)-CH3", ),
            ("CH(cyclic)-CH2-", ),
            ("CH(cyclic)-CH<", ),
            ("CH(cyclic)-C<", ),                              # 260
            ("CH(cyclic)-CH=CHn (1<n<2)", ),
            ("CH(cyclic)-C=CHn (1<n<2)", ),
            ("CH(cyclic)-Cl", ),
            ("CH(cyclic)-F", ),
            ("CH(cyclic)-OH", ),
            ("CH(cyclic)-NH2", ),
            ("CH(cyclic)-NH-CHn (0<n<3)", ),
            ("CH(cyclic)-N-CHn (0<n<3)", ),
            ("CH(cyclic)-SH", ),
            ("CH(cyclic)-CN", ),                              # 270
            ("CH(cyclic)-COOH", ),
            ("CH(cyclic)-CO", ),
            ("CH(cyclic)-NO2", ),
            ("CH(cyclic)-S-", ),
            ("CH(cyclic)-CHO", ),
            ("CH(cyclic)-O-", ),
            ("CH(cyclic)-OOCH", ),
            ("CH(cyclic)-COO-", ),
            ("CH(cyclic)-OOC-", ),
            ("C(cyclic)-CH3", ),                              # 280
            ("C(cyclic)-CH2-", ),
            ("C(cyclic)-OH", ),
            (">N(cyclic)-CH3", ),
            (">N(cyclic)-CH2-", ),
            ("AROMRING s1s2", ),
            ("AROMRING s1s3", ),
            ("AROMRING s1s4", ),
            ("AROMRING s1s2s3", ),
            ("AROMRING s1s2s4", ),
            ("AROMRING s1s2s5", ),                            # 290
            ("AROMRING s1s2s3s4", ),
            ("AROMRING s1s2s3s5", ),
            ("AROMRING s1s2s4s5", ),
            ("PYRIDINE s2", ),
            ("PYRIDINE s3", ),
            ("PYRIDINE s4", ),
            ("PYRIDINE s2s3", ),
            ("PYRIDINE s2s4", ),
            ("PYRIDINE s2s5", ),
            ("PYRIDINE s2s6", ),                              # 300
            ("PYRIDINE s3s4", ),
            ("PYRIDINE s3s5", ),
            ("PYRIDINE s2s3s6", ),

            # 3rd order term, index here are 303+, 1 -> 304
            ("HOOC-(CHn)m-COOH (m>2; 0<n<2)", ),
            ("NHn-(CHn)m-COOH (m>2; 0<n<2)", ),
            ("NH2-(CHn)m-OH (m>2; 0<n<2)", ),
            ("OH-(CHn)m-OH (m>2; 0<n<2)", ),
            ("OH-(CHp)k-O-(CHn)m-OH (m,k>0; 0<n,p<2)", ),
            ("OH-(CHp)k-S-(CHn)m-OH (m,k>0; 0<n,p<2)", ),
            ("OH-(CHp)k-NHx-(CHn)m-OH (m,k>0; 0<n,p,x<2)", ),             # 310
            ("(CHp)-O-(CHn)m-OH (m>2; 0<n,p<2)", ),
            ("NH2-(CHn)m-NH2 (m>2; 0<n<2)", ),
            ("NHk-(CHn)m-NH2 (m>2; 0<k<1; 0<n<2)", ),
            ("SH-(CHn)m-SH (m>2; 0<n<2)", ),
            ("NC-(CHn)m-CN (m>2; 0<n<2)", ),
            ("COO-(CHn)m-OOC (m>2; 0<n<2)", ),
            ("aC-(CHn=CHm)(cyclic)(fused rings) (0<n,m<1)", ),
            ("aC-aC (different rings)", ),
            ("aC-CHn(cyclic) (different rings) (0<n<1)", ),
            ("aC-CHn(cyclic) (fused rings) (0<n<1)", ),                   # 320
            ("aC-(CHn)m-aC (different rings) (m>1; 0<n<2)", ),
            ("aC-(CHn)m-CH(cyclic) (different rings) (m>0; 0<n<2)", ),
            ("CH(cyclic)-CH(cyclic) (different rings)", ),
            ("CH(cyclic)-(CHn)m-CH(cyclic)(different rings) (m>0;0<n<2)",),
            ("CH (multiring)", ),
            ("C (multiring)", ),
            ("aC-CHn-aC (different rings) (0<n<2)", ),
            ("aC-(CHm=CHn)-aC (different rings) (0<n,m<2)", ),
            ("(CHm=C)(cyclic)-CH=CH-(C=CHn)(cyclic) (different rings)", ),
            ("(CHm=C)(cyclic)-CHp-(C=CHn)(cyclic) (different rings)", ),  # 330
            ("aC-CO-aC (different rings)", ),
            ("aC-CHm-CO-aC (different rings) (0<m<2)", ),
            ("aC-CO-(C=CHn)(cyclic) (different rings) (0<n<1)", ),
            ("aC-CO-CO-aC (different rings)", ),
            ("aC-CO(cyclic) (fused rings)", ),
            ("aC-CO-(CHn)m-CO-aC (different rings) (m>0; 0<n<2)", ),
            ("aC-CO-CHn(cyclic) (different rings) (0<n<1)", ),
            ("aC-CO-NHn-aC (different rings) (0<n<1)", ),
            ("aC-NHn-CO-NHm-aC (different rings) (0<n,m<1)", ),
            ("aC-CO-N(cyclic) (different rings)", ),                      # 340
            ("aC-S(cyclic) (fused rings)", ),
            ("aC-S-aC (different rings)", ),
            ("aC-POn-aC (different rings) (0<n<4)", ),
            ("aC-SOn-aC (different rings) (0<n<4)", ),
            ("aC-NHn(cyclic) (different rings) (0<n<1)", ),
            ("aC-NH-aC (different rings)", ),
            ("aC-(C=N)(cyclic) (different rings)", ),
            ("aC-(N=CHn)(cyclic) (fused rings) (0<n<1)", ),
            ("aC-(CHn=N)(cyclic) (fused rings) (0<n<1)", ),
            ("aC-O-CHn-aC (different rings) (0<n<2)", ),                  # 350
            ("aC-O-aC (different rings)", ),
            ("aC-CHn-O-CHm-aC (different rings) (0<n,m<2)", ),
            ("aC-O(cyclic) (fused rings)", ),
            ("AROMFUSED[2]", ),
            ("AROMFUSED[2] s1", ),
            ("AROMFUSED[2] s2", ),
            ("AROMFUSED[2] s2s3", ),
            ("AROMFUSED[2] s1s4", ),
            ("AROMFUSED[2] s1s2", ),
            ("AROMFUSED[2] s1s3", ),                                      # 360
            ("AROMFUSED[3]", ),
            ("AROMFUSED[4a]", ),
            ("AROMFUSED[4a] s1", ),
            ("AROMFUSED[4a] s1s4", ),
            ("AROMFUSED[4p]", ),
            ("AROMFUSED[4p] s3s4", ),
            ("PYRIDINEFUSED[2]", ),
            ("PYRIDINEFUSED[2-iso]", ),
            ("PYRIDINEFUSED[4]", )]}

    FirstOrder = 182
    SecondOrder = 304
    ThirdOrder = 370

    def calculo(self):
        if self.kwargs["M"]:
            self.M = self.kwargs["M"]
        else:
            self.M = self._M()

        tf, tb, tc, pc = 0, 0, 0, 0
        vc = 7.95
        gf = -34.967
        hf = 5.549
        hv = 11.733
        hm = -2.806
        for i, c in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tf += c*self._coeff["tf"][i]
            tb += c*self._coeff["tb"][i]
            tc += c*self._coeff["tc"][i]
            pc += c*self._coeff["pc"][i]
            vc += c*self._coeff["vc"][i]
            gf += c*self._coeff["gf"][i]
            hf += c*self._coeff["hf"][i]
            hv += c*self._coeff["hv"][i]
            hm += c*self._coeff["hm"][i]

        # Equations of Table 1
        self.Tf = unidades.Temperature(147.45*log(tf))
        self.Tb = unidades.Temperature(222.543*log(tb))
        self.Tc = unidades.Temperature(231.239*log(tc))
        self.Pc = unidades.Pressure(5.9827+(pc+0.108998)**-2, "bar")
        self.Vc = unidades.SpecificVolume(vc/self.M, "ccg")
        self.Gf = unidades.Enthalpy(gf/self.M, "kJg")
        self.Hf = unidades.Enthalpy(hf/self.M, "kJg")
        self.Hv = unidades.Enthalpy(hv/self.M, "kJg")
        self.Hm = unidades.Enthalpy(hm/self.M, "kJg")

        GroupContribution.calculo(self)


@refDoc(__doi__, [20])
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


_methods = [Joback, Constantinou, Wilson, Marrero, Elliott, Ambrose,
            Klincewicz, Lydersen, Valderrama, Nannoolal, Wen, Li, MarreroGani]


# Other methods not implement, poor accuracy

# Daubert, T.E., Bartakovits, R.
# Prediction of Critical Temperature and Pressure of Organic Compounds by Group
# Contribution
# Ind. Eng. Chem. Res., 28, 638 (1989)
# Jalowka, J.W., Daubert, T.E.
# Group Contribution Method to Predict Critical Temperature and Pressure of
# Hydrocarbons
# Ind. Eng. Process Des. Dev., 25, 139 (1986).

# Somayajulu, G.R.
# Estimation procedures for critical constants
# J. Chem. Eng. Data 1989, 34, 106–120

# Wang, Q., Ma, P., Jia, Q., Xia, S.
# Position Group Contribution Method for the Prediction of Critical
# Temperatures of Organic Compounds
# J. Chem. Eng. Data 2008, 53, 1103-1109
# Wang, Q., Jia, Q., Ma, P.
# Position Group Contribution Method for the Prediction of Critical Pressure of
# Organic Compounds
# J. Chem. Eng. Data 2008, 53, 1877-1885
# Jia, Q., Wang, Q., Ma, P.
# Position Group Contribution Method for the Prediction of Critical Volume of
# Organic Compounds
# J. Chem. Eng. Data 2008, 53, 2606-2612

# Fishtine, S.H.
# The Modified Lydersen Method for Predicting the Critical Constants of Pure
# Substanes
# Zeitschrift für Physikalische Chemie Neue Folge, 123 (1980), 39-49


# Other method not implement because calculate only a property, not general
# and don't appropiated to define a compound
# Aditional reference in Perry's, pag. 2-471

# Vetere, A.
# An Empirical Method for Evaluating Critical Molar Volumes
# AIChE I., 22,950 (1976); Errata, ibtd., 23,406 (1977)

# Fedors, R.F.
# A Method to Estimate Critical Volumes
# AIChE J., 25(1), 202 (1979).
# Fedors, R.F.
# A Relationship between Chemical Structure and the Critical Temperature
# Chem. Eng. Commun., 16, 149 (1982).

# Pailhes, F.
# Estimation of the Boiling Temperature at Normal Pressure for Organic
# Compounds from their Chemical Formula and a Known Boiling Temperature at Low
# Pressure
# Fluid Phase Equilib., 41 (1988) 97-107

# Dalmazzone, D., Salmon, A., Guella, S.
# A second order group contribution method for the prediction of critical
# temperatures and enthalpies of vaporization of organic compounds
# Fluid Phase Equilib. 242, (2006) 29–42.

# Thinh, T.P., Trong, T.K.
# Estimation of Standard Heats of Formation ΔHf, Standard Entropies of
# Formation ΔSf, Standard Free Energies of Formation ΔGf and Absolute Entropies
# S of Hydrocarbons from Group Contributions: An Accurate Approach
# Can. J.  Chern. Eng. 54 (1976), 344-357

# van Krevelen, S.W., Chermin, H.A.G.
# Estimation of the Free Enthalpy (Gibbs Free Energy) of Formation of Organic
# Compounds from Group Contributions
# Chem. Eng. Sci. 1 (1951) 66-80

# Ruzicka method for liquid specific heat
# Ruzicka, V., Domalski, E.S.
# J. Phys. Chem. Ref. Data, 22 (1993): 597, 619

# Benson method for ideal gas specific heat
# Benson, S.W., Cruickshank, F.R., Golden, n.M., Haugen, G.R., O'Neal, H.E.,
# Rogers, A.S., Shaw, R., Walsh, R.
# Additivity Rules for the Estimation of Thermochemical Properties
# Chem. Rev. 69 (1969) 279
# Benson, S.W., Buss, J.H.
# Additivity Rules for the Estimation of Molecular Properties. Thermodynamic
# Properties
# J. Chem. Phys. 29 (1958) 546-572
# Cohen, N., Benson, S.W.
# Estimation of Heats of Formation of Organic Compounds by Additivity Methods
# Chem. Rev. 93 (1993) 2419-2438
# O'Neal, H.E., Benson, S.W.
# Entropies and Heat Capacities of Cyclic and Polycyclic Compounds
# J. Chem. Eng. Data 15 (1970) 266-276

# Rihani, D.N., Doraiswamy, L.K.
# Estimation of Heat Capacity of Organic Compounds from Group Contributions
# Ind. Eng. Chern. Fundam. 4, (1965) 17-21

# Verma, K.K., Doraiswamy, L.K.
# Estimation of Heats of Formation of Organic Compounds
# Ind. Eng. Chern. Fundam. 4, (1965) 389-396

# van Velzen, D., Cardozo, R.L., Langenkarnp, H.
# A Liquid Viscosity-Temperature-Chemical Constitution Relation for Organic
# Compounds
# Ind. Eng. Chem. Fundam. 11 (1972) 20-25

# Franklin, I.L.
# Prediction of Heat and Free Energies of Organic Compounds
# Ind. Eng. Chem. 41, (1949) 1070-1076

# Chickos, J.S., Braton, C.M., Hesse, D.G.
# Estimating Entropies and Enthalpies of Fusion of Organic Compounds
# J. Org. Chem. 56 (1991) 927-938
# Chickos, J.S., Hosseini, S.
# A Group Additivity Approach for the Estimation of Vapor Pressures of Liquid
# Hydrocarbons from 298 to 500 K
# J. Org. Chem. 58 (1993) 5345-5350
# Chickos, J.S., Hesse, D.G., Liebman, J.F., Panshin, S.Y.
# Estimations of the Heats of Vaporization of Simple Hydrocarbon Derivatives at
# 298 K
# J. Org. Chem. 53 (1988) 3424-3429
# Chickos, J.S., Hesse, D.G., Liebman, J.F.
# Estimating Entropies and Enthalpies of Fusion of Hydrocarbons
# J. Org. Chem. 55 (1990) 3833-3840
# Chickos, J.S., Annunziata, R., Ladon, L.H., Liebman, J.F.
# Estimating Heats of Sublimation of Hydrocarbons. A Semiempirical Approach
# J. Org. Chem. 51 (1986) 4311-4314

# Domalski method for formation properties
# Domalski, E. S., and E. D. Hearing, J. Phys. Chem. Ref. Data, 22 (1993): 805

# Goodman method for solid specific heat
# Goodman, B. T., et al., J. Chem. Eng. Data, 49 (2004): 24.

# Reichenberg method for gas viscosity
# Reichenberg, D., AIChE J., 21 (1975): 181.

# Hsu method for liquid viscosity
# Hsu, H.-C., Y.-W. Sheu, and C.-H. Tu, Chem. Eng. J., 88 (2002): 27.

# Parachor method for surface tension
# Macleod, D. B., Trans. Faraday Soc., 19 (1923): 38
# Sugden, S. J., Chem. Soc., 125 (1924): 32
# Knotts, T. A., et al., J. Chem. Eng. Data, 46 (2001): 1007.
