#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to add EoS common functionality
###############################################################################

from scipy import exp, log, log10, tan, sinh, tanh, arctan, sqrt
from scipy import roots, r_
from scipy.constants import pi, Avogadro, R
from scipy.optimize import fsolve

from lib import unidades, config
from lib.physics import R_atml, factor_acentrico_octano

from lib.eos import EoS

class Grayson_Streed(EoS):
    """Ecuación de estado de Grayson Streed modificada por Chao-Seader
    Chao, K.C. and Seader, J.D.; A General Correlation of Vapor-Liquid Equilibria in Hydrocarbon Mixtures, AIChE Journal, 7, No 4 (December 1961"""
    __title__="Grayson Streed (1961)"
    __status__="GS"

    def __init__(self, T, P, mezcla):
        self.T=unidades.Temperature(T)
        self.P=unidades.Pressure(P, "atm")
        self.mezcla=mezcla
        self.componente=mezcla.componente
        self.fraccion=mezcla.fraccion
        p=unidades.Pressure(P, "atm")

        self.rk=RK(T, P, mezcla)
        self.Z=self.rk.Z
        self.V=self.rk.V
        self.x, self.xi, self.yi, self.Ki=self._Flash()

    def _k(self, xi, yi):
        fi=self.rk._fug(self.rk.Z[0], yi)
        suma1=0
        suma2=0
        Vi=[]
        for i in range(len(self.componente)):
            Vi.append(self.componente[i].Vc*self.componente[i].M)
            suma1+=xi[i]*Vi[i]*self.componente[i].parametro_solubilidad
            suma2+=xi[i]*Vi[i]
        dim=suma1/suma2

        gi=[]
        for i in range(len(self.componente)):
            gi.append(exp(Vi[i]/1e6*(self.componente[i].parametro_solubilidad-dim)**2/R/self.T))
        nio=[]
        for i in self.componente:
            tr=i.tr(self.T)
            pr=i.pr(self.P.atm)
            if i.indice==1:
                A=[1.50709, 2.74283, -0.02110, 0.00011, 0.0, 0.008585, 0., 0., 0., 0.]
            elif i.indice==2:
                A=[1.36822, -1.54831, 0., 0.02889, -0.01076, 0.10486, -0.02529, 0., 0., 0.]
            else:
                A=[2.05135, -2.10899, 0., -0.19396, 0.02282, 0.08852, 0., -0.00872, -0.00353, 0.00203]
            logn0=A[0]+A[1]/tr+A[2]*tr+A[3]*tr**2+A[4]*tr**3+(A[5]+A[6]*tr+A[7]*tr**2)*pr+(A[8]+A[9]*tr)*pr**2-log10(pr)
            logn1=-4.23893+8.65808*tr-1.2206/tr-3.15224*tr**3-0.025*(pr-0.6)
            nio.append(10**(logn0+i.f_acent*logn1))


        tital=[]
        for i in range(len(self.componente)):
            tital.append(nio[i]*gi[i])

        return tital, fi


    def _Flash(self):
        """Cálculo de los coeficientes de reparto entre fases, Ref Naji - Conventional and rapid flash claculations"""
        #Estimación inicial de K mediante correlación wilson Eq 19
        Ki=[]
        for i in self.componente:
            Ki.append(i.Pc/self.P*exp(5.37*(1.+i.f_acent)*(1.-i.Tc/self.T)))

        Rachford=lambda x: sum([zi*(ki-1.)/(1.-x+x*ki) for zi, ki in zip(self.fraccion, Ki)])

        if Rachford(0)>0 and Rachford(1)>0: #x>1, por tanto solo fase vapor
            xi=self.fraccion
            yi=self.fraccion
            x=1.
        elif Rachford(0)<0 and Rachford(1)<0: #x<0, por tanto solo fase líquida
            xi=self.fraccion
            yi=self.fraccion
            x=0.
        else:
            x=0.5
            while True:
                xo=x
                solucion=fsolve(Rachford, x, full_output=True)
                if solucion[2]!=1:
                    print solucion
                    break
                else:
                    x=solucion[0]
                    xi=[]
                    yi=[]
                    for zi, ki in zip(self.fraccion, Ki):
                        xi.append(zi/(1-x+x*ki))
                        yi.append(zi*ki/(1-x+x*ki))

                    titas=self._k(xi, yi)
                    tital=titas[0]
                    titav=titas[1]

                    fiv=[z*t*self.P for z, t in zip(yi, titav)]
                    fil=[z*t*self.P for z, t in zip(xi, tital)]
                    #criterio de convergencia Eq 21
                    if sum([abs(l/v-1) for l, v in zip(fil, fiv)])< 1e-12 and (x-xo)**2 < 1e-15:
                        break
                    else:
                        Ki=[l/v for l, v in zip(tital, titav)]

        return x, xi, yi, Ki





if __name__ == "__main__":
    from corriente import Mezcla

#    eq=SRK_API(340, 1, Mezcla([10, 38, 22, 61], [0.3, 0.5, 0.05, 0.15]))
#    print eq.x
#    eq=SRK(340, 1., Mezcla([10, 38, 22, 61], [0.3, 0.5, 0.05, 0.15]))
#    print eq.x
#    print eq._Dew_T()

#    p=unidades.Pressure(2000, "psi")
#    t=unidades.Temperature(100, "F")
#    eq=van_Waals(t.K, p.atm, Mezcla([1, 4], [.805, 0.195]))
#    print eq.Z

#    mezcla=SRK(350, 1, Mezcla([4, 5, 6, 7, 8, 10, 11, 12, 13], [0.02361538, 0.2923077, 0.3638462, 0.02769231, 0.01153846, 0.01769231, 0.03007692, 0.2093846, 0.02384615]))
#    mezcla=Grayson_Streed(350, 1, Mezcla([4, 5, 6, 7, 8, 10, 11, 12, 13], [0.02361538, 0.2923077, 0.3638462, 0.02769231, 0.01153846, 0.01769231, 0.03007692, 0.2093846, 0.02384615]))

#    print mezcla._Bubble_T(), mezcla._Dew_T()
#    print mezcla.T, mezcla.x

    mezcla=Mezcla(ids=[10, 38, 22, 61], fraccionMolar=[0.3, 0.5, 0.05, 0.15])
#    for t in range(300, 345, 1):
#        eq=SRK(t, 1., mezcla )
#        print t, eq.x
    eq=SRK(340, 1., mezcla )
    print eq.x
