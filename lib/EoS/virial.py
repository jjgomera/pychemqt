#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Virial ecuation of state implementation
###############################################################################

from scipy import zeros

from lib import unidades
from lib.eos import EoS
from lib.physics import R_atml

Database = {
    98: [3.4162e1, -1.2087e4, -7-6702e5, -1.96e7],
    630: [9.1039e1, -5.9081e4, 1.0478e7, 3.0463e9],
    113: [-6.0611e3, 4.6389e6, 9.8352e8],
    48: [4.8826e1, -1.5614e4, -2.757e5, -4.7684e7],
    49: [5.74e1, -3.8829e4, 4.2899e5, -1.4661e9],
    102: "Carbon disulfide pag 35"
}


class Virial(EoS):
    """Class to model virial equation of state"""
    def __init__(self, T, P, mezcla):
        self.T = unidades.Temperature(T)
        self.P = unidades.Pressure(P, "atm")
        self.mezcla = mezcla
        self.componente = mezcla.componente
        self.xi = mezcla.fraccion

    def _Bi(self):
        B = []
        for comp in self.componente:
            if comp.indice in Database:
                Bi = 0
                for i, a in enumerate(Database[comp.indice]):
                    Bi += a/self.T**i
            else:
                # Use general correlation
                Bi = self._B_Tsonopoulos(comp.Tc, comp.Pc, comp.f_acent)
        B.append(Bi)
        return B

    def _Ci(self):
        C = []
        for comp in self.componente:
            C.append(self._C_Orbey_vera(comp.Tc, comp.Pc, comp.f_acent))
        return C

    def B(self):
        Bi = self._Bi()
        B = 0
        for i, xi in enumerate(self.xi):
            for j, xj in enumerate(self.xi):
                if i == j:
                    Bij = Bi[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/(ci.Vc**(1./3)+cj.Vc**(1./3))**3
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    Bij = self._B_Tsonopoulos(Tcij, Pcij, wij)
                B += xi*xj*Bij
        return B

    def C(self):
        Ci = self._Ci()
        Cij = zeros(len(Ci), len(Ci))
        for i, xi in enumerate(self.xi):
            for j, xj in enumerate(self.xi):
                if i == j:
                    Cij[i, j] = Ci[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/(ci.Vc**(1./3)+cj.Vc**(1./3))**3
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    Cij[i, j] = self._C_Orbey_vera(Tcij, Pcij, wij)

        C = 0
        for i, xi in enumerate(self.xi):
            for j, xj in enumerate(self.xi):
                for k, xk in enumerate(self.xi):
                    C += xi*xj*xk*(Cij[i, j]*Cij[j, k]*Cij[j, k])**(1./3)
        return C

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

    def _physics(self):
        self.U_exc = -R*T*(B1/V+C1/2/V**2)
        self.H_exc = R*T*((B-B1)/V+(2*C-C1)/2/V**2)
        self.Cv_exc = -R*((2*B1+B2)/V+(2*C1+C2)/2/V**2)
        self.Cp_exc = -R*(B2/V-((B-B1)**2-(C-C1)-C2/2)/V**2)
        self.S_exc = -R*(log(P)+B1/V+(B**2-C+C1)/2/V**2)
        self.A_exc = R*T*(log(P)+(B**2-C/2/V**2))
        self.G_exc = R*T*(log(P)+B/V+(B**2+C)/2/V**2)

        self.fug = P*exp(B/V+(C+B**2)/2/V**2)
        self.Z = 1+B*(P/R_atml/T)+(C-B**2)*(P/R/T)**2

    def _B_Tsonopoulos(self, Tc, Pc, w):
        # TODO: Tras añadir características quimicas a la base de datos poder calcular a y b
        a = b = 0
        Tr = self.T/Tc
        f0 = 0.1445-0.33/Tr-0.1385/Tr**2-0.0121/Tr**3-0.000607/Tr**8
        f1 = 0.0637+0.331/Tr**2-0.423/Tr**3-0.008/Tr**8
        f2 = 1/Tr**6
        f3 = -1/Tr**8
        f = f0 + w*f1 + a*f2 + b*f3
        return f*R_atml*Tc/Pc.atm

    def _C_Orbey_vera(self, Tc, Pc, w):
        Tr = self.T/Tc
        g0 = 0.01407+0.02432/Tr**2.8-0.00313/Tr**10.5
        g1 = -0.02676+0.0177/Tr**2.8+0.04/Tr**3-0.003/Tr**6-0.00228/Tr**10.5
        g = g0+w*g1
        return g*R_atml**2*Tc**2/Pc.atm**2


if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[62], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Virial(T, 1, mezcla)
        print eq.B()
