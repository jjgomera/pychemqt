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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


.. include:: newComponent.rst


Base classes with common functionality
    * :class:`newComponente`: Base class for new component definition
    * :class:`GroupContribution`: Common group contribution functionality

List of classes with the implemented methods:
    * :class:`Joback`: Group contribution method of Joback
    * :class:`Constantinou`: Group contribution method of Constaninou and Gani
    * :class:`Wilson`: Group contribution method of Wilson and Jasperson
    * :class:`Marrero`: Group contribution method of Marrero and Pardillo
    * :class:`Elliott`: Group contribution method of Elliott
    * :class:`Ambrose`: Group contribution method of Ambrose
    * :class:`Klincewicz`: Group contribution method of Klincewicz
    * :class:`Lydersen`: Group contribution method of Lydersen
    * :class:`Valderrama`: Group contribution method of Valderrama
    * :class:`Nannoolal`: Group contribution method of Nannoolal
    * :class:`Wen`: Group contribution method of Wen-Qiang
    * :class:`Li`: Group contribution method of Li-Xia-Xiang
    * :class:`MarreroGani`: Group contribution method of Marrero-Gani

Other methods:
    * :func:`cpLS_Hurst`: Liquid-Solid heat capacities using the Hurst-Harrison
    method
'''


###############################################################################
# lib.newComponent module
# module with group contribution methods of definition of custom compound
###############################################################################


import glob
import os

from lib.newComponent._base import newComponente, GroupContribution


files = glob.glob(os.path.join("lib", "newComponent", "*.py"))

__all__ = []
for file in files:
    fname, ext = os.path.splitext(os.path.basename(file))
    if fname != "__init__" and fname not in ["_base"]:
        __all__.append(fname)
        __import__(f"lib.newComponent.{fname}")

_methods = GroupContribution.__subclasses__()


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
