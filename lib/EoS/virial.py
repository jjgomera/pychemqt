#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Virial equation of state implementation
###############################################################################

from scipy import zeros, log, exp
from scipy.constants import R

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
    def __init__(self, *args, **kwargs):
        EoS.__init__(self, *args, **kwargs)
        self._physics(*args)

    def _Bi(self):
        B = []
        Bt = []
        Btt = []
        for comp in self.componente:
            if comp.indice in Database:
                Bi = 0
                Bit = 0
                Bitt = 0
                for i, a in enumerate(Database[comp.indice]):
                    Bi += a/self.T**i
                    Bit += -a*i*self.T**(i-1)
                    Bitt += a*i*(i-1)*self.T**(i-2)
            else:
                # Use general correlation
                Bi, Bit, Bitt = self._B_Tsonopoulos(comp.Tc, comp.Pc, comp.f_acent)
        B.append(Bi)
        Bt.append(Bit)
        Btt.append(Bitt)
        return B, Bt, Btt

    def _Ci(self):
        C = []
        Ct = []
        Ctt = []
        for comp in self.componente:
            if self.kwargs.get("C", 0):
                Ci, Cit, Citt = self._C_Orbey_Vera(comp.Tc, comp.Pc, comp.f_acent)
            else:
                Ci, Cit, Citt = self._C_Liu_Xiang(comp.Tc, comp.Pc, comp.f_acent, comp.Zc)
            C.append(Ci)
            Ct.append(Cit)
            Ctt.append(Citt)
        return C, Ct, Ctt

    def B(self):
        Bi, Bit, Bitt = self._Bi()
        B, Bt, Btt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                if i == j:
                    Bij = Bi[i]
                    Bijt = Bit[i]
                    Bijtt = Bitt[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/(ci.Vc**(1./3)+cj.Vc**(1./3))**3
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    Bij, Bijt, Bijtt = self._B_Tsonopoulos(Tcij, Pcij, wij)
                B += xi*xj*Bij
                Bt += xi*xj*Bijt
                Btt += xi*xj*Bijtt
        return B, Bt, Btt

    def C(self):
        Ci, Cit, Citt = self._Ci()
        Cij = zeros((len(Ci), len(Ci)))
        Cijt = zeros((len(Ci), len(Ci)))
        Cijtt = zeros((len(Ci), len(Ci)))
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                if i == j:
                    Cij[i, j] = Ci[i]
                    Cijt[i, j] = Cit[i]
                    Cijtt[i, j] = Citt[i]
                else:
                    ci = self.componente[i]
                    cj = self.componente[j]
                    Tcij = (ci.Tc*cj.Tc)**0.5
                    Vcij = (ci.Vc**(1./3)+cj.Vc**(1./3))**3
                    Pcij = 4*Tcij*(ci.Zc+cj.Zc)/Vcij
                    wij = 0.5*(ci.f_acent+cj.f_acent)
                    if self.kwargs.get("C", 0):
                        Cij[i, j] = self._C_Orbey_Vera(Tcij, Pcij, wij)
                    else:
                        Zcij = Pcij*Vcij/R_atml/Tcij
                        Cij[i, j] = self._C_Liu_Xiang(Tcij, Pcij, wij, Zcij)

        C, Ct, Ctt = 0, 0, 0
        for i, xi in enumerate(self.fraccion):
            for j, xj in enumerate(self.fraccion):
                for k, xk in enumerate(self.fraccion):
                    C += xi*xj*xk*(Cij[i, j]*Cij[j, k]*Cij[j, k])**(1./3)
        return C, Ct, Ctt

    def _physics(self, T, P, mezcla):
        B, B1, B2 = self.B()
        C, C1, C2 = self.C()

        self.Z = 1+B*(P/R_atml/T)+(C-B**2)*(P/R/T)**2
        V = self.Z*R*T/P
        self.U_exc = -R*T*(B1/V+C1/2/V**2)
        self.H_exc = R*T*((B-B1)/V+(2*C-C1)/2/V**2)
        self.Cv_exc = -R*((2*B1+B2)/V+(2*C1+C2)/2/V**2)
        self.Cp_exc = -R*(B2/V-((B-B1)**2-(C-C1)-C2/2)/V**2)
        self.S_exc = -R*(log(P)+B1/V+(B**2-C+C1)/2/V**2)
        self.A_exc = R*T*(log(P)+(B**2-C/2/V**2))
        self.G_exc = R*T*(log(P)+B/V+(B**2+C)/2/V**2)

        self.fug = P*exp(B/V+(C+B**2)/2/V**2)

    def _B_Tsonopoulos(self, Tc, Pc, w):
        # TODO: Tras añadir características quimicas a la base de datos poder calcular a y b
        a = b = 0
        Tr = self.T/Tc
        f0 = 0.1445-0.33/Tr-0.1385/Tr**2-0.0121/Tr**3-0.000607/Tr**8
        f1 = 0.0637+0.331/Tr**2-0.423/Tr**3-0.008/Tr**8
        f2 = 1/Tr**6
        f3 = -1/Tr**8
        f = f0 + w*f1 + a*f2 + b*f3
        return f*R_atml*Tc/Pc.atm, 0, 0

    def _C_Orbey_Vera(self, Tc, Pc, w):
        Tr = self.T/Tc
        g0 = 0.01407+0.02432/Tr**2.8-0.00313/Tr**10.5
        g1 = -0.02676+0.0177/Tr**2.8+0.04/Tr**3-0.003/Tr**6-0.00228/Tr**10.5
        g = g0+w*g1
        return g*R_atml**2*Tc**2/Pc.atm**2, 0, 0

    def _C_Liu_Xiang(self, Tc, Pc, w, Zc):
        """Liu, D.X., Xiang, H.W.: Corresponding-States Correlation and Prediction of Third Virial Coefficients for a Wide Range of Substances. International Journal of Thermophysics, November 2003, Volume 24, Issue 6, pp 1667-1680"""
        Tr = self.T/Tc
        X = (Zc-0.29)**2
        g0 = 0.1623538+0.3087440/Tr**3-0.01790184/Tr**6-0.02789157/Tr**11
        g1 = -0.5390344+1.783526/Tr**3-1.055391/Tr**6+0.09955867/Tr**11
        g2 = 34.22804-74.76559/Tr**3+279.9220/Tr**6-62.85431/Tr**11
        g = g0+w*g1+X*g2
        return g*Zc**2*R_atml**2*Tc**2/Pc.atm**2, 0, 0


if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Virial(T, 1, mezcla)
        print(eq.H_exc)
