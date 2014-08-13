#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Virial ecuation of state implementation
###############################################################################

from lib.EoS import EoS
from physics import R_atml

class Virial(EoS):
    """Virial equation of state"""
    pass

    def _B1(self):
        """First derivative first virial term versus T"""
        pass

    def _C1(self):
        """First derivative second virial term versus T"""
        pass

    def _B2(self):
        """Second derivative first virial term versus T"""
        pass

    def _C2(self):
        """Second derivative second virial term versus T"""
        pass

    def _Enthalpy(self):
        return R*T*((B-B1)/V+(2*C-C1)/2/V**2)

    def _B_Tsonopoulos(self, componente, T):
        # TODO: Tras añadir características quimicas a la base de datos poder calcular a y b
        a = b = 0
        Tr = T/componente.Tc
        f0 = 0.1445-0.33/Tr-0.1385/Tr**2-0.0121/Tr**3-0.000607/Tr**8
        f1 = 0.0637+0.331/Tr**2-0.423/Tr**3-0.008/Tr**8
        f2 = 1/Tr**6
        f3 = -1/Tr**8
        f = f0 + componente.f_acent*f1 + a*f2 + b*f3
        return f*R_atml*componente.Tc/componente.Pc.atm


    def _C_Orbey_vera(self, componente, T):
        Tr = T/componente.Tc
        g0 = 0.01407+0.02432/Tr**2.8-0.00313/Tr**10.5
        g1 = -0.02676+0.0177/Tr**2.8+0.04/Tr**3-0.003/Tr**6-0.00228/Tr**10.5
        g = g0+componente.f_acent*g1
        return g*R_atml**2*componente.Tc**2/componente.Pc.atm**2
