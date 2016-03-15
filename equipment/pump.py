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
#   Library for pump equipment definition
###############################################################################

import os

from PyQt5.QtWidgets import QApplication
from scipy import log, exp, optimize, polyval, roots, r_
from scipy.constants import g
#import xlwt

from lib.unidades import (Pressure, Length, Power, Temperature, VolFlow,
                          Currency, Dimensionless, DeltaP, DeltaT)
from lib.corriente import Corriente
#from lib.datasheet import pdf
from .parents import equipment


class Pump(equipment):
    """Clase que modela un equipo de bombeo de líquidos

    Parámetros:
        entrada: instancia de la clase corriente que define la corriente de entrada en la bomba
        usarCurva:
            0   -   Funcionamiento fijo
            1   -   Funcionamiento siguiendo la curva característica
        incognita: Indica la variable a calcular en el caso de que se use la curva caracteristica
            0   -   Carga (incremento de presión)
            1   -   Caudal, en este caso sobreescribe el caudal de la corriente de entrada
        rendimiento: rendimiento de la bomba, no necesario si se ha definido la curva característica de la bomba
        deltaP: incremento de presión proporcionado por la bomba, no será necesario si se define la curva característica y es el caudal la variable a calcular, en atmosferas
        Pout: Presión a la salida de la bomba
        Carga: Cambio de presión basada en altura de columna de fluido
        curvaCaracteristica: array que define la curva característica, en la forma: [Dia, rpm, [Q1,..Qn], [h1,...,hn], [Pot1,...,Potn], [NPSH1,...NPSHn]].
        diametro: diametro nominal bomba
        velocidad: velocidad de giro de la bomba

    Coste
        tipo_bomba
            0   -   Centrifugal pumps
            1   -   Reciprocating pumps
            2   -   Gear pumps
            3   -   Vertical mixed flow
            4   -   Vertical axial flow
        tipo_centrifuga
            0   -   One stage, 3550 rpm, VSC
            1   -   One stage, 1750 rpm, VSC
            2   -   One stage, 3550 rpm, HSC
            3   -   One stage, 1750 rpm, HSC
            4   -   Two stage, 3550 rpm, HSC
            5   -   Multistage, 3550 rpm, HSC
        Material
            0   -   Cast iron
            1   -   Case steel
            2   -   304 or 316 fittings
            3   -   Stainless steel 304 or 316
            4   -   Case Gould's alloy no. 20
            5   -   Nickel
            6   -   Monel
            7   -   ISO B
            8   -   ISO B
            9   -   Titanium
            10  -   Hastelloy C
            11  -   Ductile iron
            12  -   Bronze
        motor: Tipo de motor
            0   -   Open drip-proof
            1   -   Totally enclosed, fan-cooled
            2   -   Explosion-proof
        rpm
            0   -   3600 rpm
            1   -   1800 rpm
            2   -   1200 rpm

    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    >>> bomba=Pump(entrada=agua, rendimiento=0.75, deltaP=20*101325, tipo_bomba=1)
    >>> print bomba.power.hp
    3.63596053358
    >>> print bomba.C_inst
    3493.24497823
    """
    title=QApplication.translate("pychemqt", "Pump")
    help=""
    kwargs={"entrada": None,
                    "usarCurva": 0,
                    "incognita": 0,
                    "rendimiento": 0.0,
                    "deltaP": 0.0,
                    "Pout": 0.0,
                    "Carga": 0.0,
                    "curvaCaracteristica": [],
                    "diametro": 0.0,
                    "velocidad": 0.0,

                    "f_install": 2.8,
                    "Base_index": 0.0,
                    "Current_index": 0.0,
                    "tipo_bomba": 0,
                    "tipo_centrifuga": 0,
                    "material": 0,
                    "motor": 0,
                    "rpm": 0}
    kwargsInput=("entrada", )
    kwargsCheck=("usarCurva", )
    kwargsValue=("Pout", "deltaP", "rendimiento", "Carga", "diametro", "velocidad")
    kwargsList=("incognita", "tipo_bomba", "tipo_centrifuga", "material", "motor", "rpm")
    calculateValue=("PoutCalculada", "power", "headCalculada", "volflow", "rendimientoCalculado")
    calculateCostos=("C_bomba", "C_motor", "C_adq", "C_inst")
    indiceCostos=7
    salida = [None]

    TEXT_BOMBA=[QApplication.translate("pychemqt", "Centrifugal"),
                            QApplication.translate("pychemqt", "Reciprocating"),
                            QApplication.translate("pychemqt", "Gear pump"),
                            QApplication.translate("pychemqt", "Vertical mixed flow"),
                            QApplication.translate("pychemqt", "Vertical axial flow")]
    TEXT_CENTRIFUGA=[QApplication.translate("pychemqt", "One stage, 3550 rpm, VSC"),
                                    QApplication.translate("pychemqt", "One stage, 1750 rpm, VSC"),
                                    QApplication.translate("pychemqt", "One stage, 3550 rpm, HSC"),
                                    QApplication.translate("pychemqt", "One stage, 1750 rpm, HSC"),
                                    QApplication.translate("pychemqt", "Two stage, 3550 rpm, HSC"),
                                    QApplication.translate("pychemqt", "Multistage, 3550 rpm, HSC")]
    TEXT_MATERIAL=[QApplication.translate("pychemqt", "Cast iron"),
                                QApplication.translate("pychemqt", "Case steel"),
                                QApplication.translate("pychemqt", "304 or 316 fittings"),
                                QApplication.translate("pychemqt", "Stainless steel 304 or 316"),
                                QApplication.translate("pychemqt", "Case Gould's alloy no. 20"),
                                QApplication.translate("pychemqt", "Nickel"),
                                QApplication.translate("pychemqt", "Monel (Ni-Cu)"),
                                QApplication.translate("pychemqt", "ISO B"),
                                QApplication.translate("pychemqt", "ISO C"),
                                QApplication.translate("pychemqt", "Titanium"),
                                QApplication.translate("pychemqt", "Hastelloy C (Ni-Fe-Mo)"),
                                QApplication.translate("pychemqt", "Ductile iron"),
                                QApplication.translate("pychemqt", "Bronze")]
    TEXT_MOTOR=[QApplication.translate("pychemqt", "Open drip-proof"),
                            QApplication.translate("pychemqt", "Totally enclosed, fan-cooled"),
                            QApplication.translate("pychemqt", "Explosion-proof")]
    TEXT_RPM=["3600 RPM", "1800 RPM", "1200 RPM"]


    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and self.kwargs["Current_index"]:
            self.statusCoste=True
        else:
            self.statusCoste=False

        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
        else:
            presion=self.kwargs["Pout"] or self.kwargs["deltaP"] or self.kwargs["Carga"]
            if self.kwargs["usarCurva"]:
                if self.kwargs["incognita"]:
                    if presion and self.kwargs["curvaCaracteristica"]:
                        self.msg=""
                        self.status=1
                        return True
                    elif presion:
                        self.msg=QApplication.translate("pychemqt", "undefined pump curve")
                        self.status=0
                    else:
                        self.msg=QApplication.translate("pychemqt", "undefined out pressure condition")
                        self.status=0
                elif self.kwargs["curvaCaracteristica"]:
                    self.msg=""
                    self.status=1
                    return True
                else:
                    self.msg=QApplication.translate("pychemqt", "undefined pump curve")
                    self.status=0
            else:
                if presion and self.kwargs["rendimiento"]:
                    self.msg=""
                    self.status=1
                    return True
                elif presion:
                    self.msg=QApplication.translate("pychemqt", "undefined efficiency")
                    self.status=0
                else:
                    self.msg=QApplication.translate("pychemqt", "undefined out pressure condition")
                    self.status=0


    def calculo(self):
        entrada=self.kwargs["entrada"]
        self.rendimientoCalculado=Dimensionless(self.kwargs["rendimiento"])

        if self.kwargs["Pout"]:
            DeltaP=Pressure(self.kwargs["Pout"]-entrada.P)
        elif self.kwargs["deltaP"]:
            DeltaP=Pressure(self.kwargs["deltaP"])
        elif self.kwargs["Carga"]:
            DeltaP=Pressure(self.kwargs["Carga"]*entrada.Liquido.rho*g)
        else:
            DeltaP=Pressure(0)

        if self.kwargs["usarCurva"]:
            if self.kwargs["diametro"]!=self.kwargs["curvaCaracteristica"][0] or self.kwargs["velocidad"]!=self.kwargs["curvaCaracteristica"][1]:
                self.curvaActual=self.calcularCurvaActual()
            else:
                self.curvaActual=self.kwargs["curvaCaracteristica"]
            self.Ajustar_Curvas_Caracteristicas()

        if not self.kwargs["usarCurva"]:
            head=Length(DeltaP/g/entrada.Liquido.rho)
            power=Power(head*g*entrada.Liquido.rho*entrada.Q/self.rendimientoCalculado)
            P_freno=Power(power*self.rendimientoCalculado)
        elif not self.kwargs["incognita"]:
            head=Length(polyval(self.CurvaHQ,entrada.Q))
            DeltaP=Pressure(head*g*entrada.Liquido.rho)
            power=Power(entrada.Q*DeltaP)
            P_freno=Power(polyval(self.CurvaPotQ,entrada.Q))
            self.rendimientoCalculado=Dimensionless(power/P_freno)
        else:
            head=Length(self.DeltaP/g/entrada.Liquido.rho)
            caudalvolumetrico=roots([self.CurvaHQ[0], self.CurvaHQ[1], self.CurvaHQ[2]-head])[0]
            power=Power(caudalvolumetrico*self.DeltaP)
            entrada=Corriente(entrada.T, entrada.P.atm, caudalvolumetrico*entrada.Liquido.rho*3600, entrada.mezcla, entrada.solido)
            P_freno=Power(polyval(self.CurvaPotQ,caudalvolumetrico))
            self.rendimientoCalculado=Dimensionless(power/P_freno)

        self.deltaP = DeltaP
        self.headCalculada=head
        self.power=power
        self.P_freno=P_freno
        self.salida=[entrada.clone(P=entrada.P+DeltaP)]
        self.Pin=entrada.P
        self.PoutCalculada=self.salida[0].P
        self.volflow=entrada.Q

    def Ajustar_Curvas_Caracteristicas(self):
        """Define la curva característica de la bomba a partir de los datos, todos ellos en forma de array de igual dimensión
        Q: caudal, m3/s
        h: carga, m
        mu: rendimiento
        NPSHr: carga neta de aspiración requerida por la bomba para no entrar en cavitación
        """
        Q=r_[self.curvaActual[2]]
        h=r_[self.curvaActual[3]]
        Pot=r_[self.curvaActual[4]]
        NPSH=r_[self.curvaActual[5]]

        funcion_h = lambda p, x: p[0]*x**2+p[1]*x+p[2] # Función a ajustar
        residuo_h = lambda p, x, y: funcion_h(p, x) - y # Residuo
        inicio=r_[0, 0, 0]
        ajuste_h, exito_h=optimize.leastsq(residuo_h,inicio,args=(Q, h))
        self.CurvaHQ=ajuste_h

        funcion_Pot = lambda p, x: p[0]*x**2+p[1]*x+p[2] # Función a ajustar
        residuo_Pot = lambda p, x, y: funcion_Pot(p, x) - y # Residuo
        inicio=r_[0, 0, 0]
        ajuste_Pot, exito_Pot=optimize.leastsq(residuo_Pot,inicio,args=(Q, Pot))
        self.CurvaPotQ=ajuste_Pot

        funcion_NPSH = lambda p, x: p[0]+p[1]*exp(p[2]*x) # Función a ajustar
        residuo_NPSH = lambda p, x, y: funcion_NPSH(p, x) - y # Residuo
        inicio=r_[0, 0, 0]
        ajuste_NPSH, exito_NPSH=optimize.leastsq(residuo_NPSH,inicio,args=(Q, NPSH))
        self.CurvaNPSHQ=ajuste_NPSH


    def calcularCurvaActual(self):
        """Método que define curvas caracteristica de la bomba a una velocidad de giro y diametro del propulsor diferente de la curva característica original haciendo uso de las leyes de afinidad, perry 10.25 Table 10.7"""

        D1=self.kwargs["curvaCaracteristica"][0]
        N1=self.kwargs["curvaCaracteristica"][1]
        D2=self.kwargs["diametro"]
        N2=self.kwargs["velocidad"]

        Q1=r_[self.kwargs["curvaCaracteristica"][2]]
        h1=r_[self.kwargs["curvaCaracteristica"][3]]
        Pot1=r_[self.kwargs["curvaCaracteristica"][4]]
        npsh1=r_[self.kwargs["curvaCaracteristica"][5]]
        Q2=Q1*D2/D1*N2/N1
        h2=h1*N2**2/N1**2*D2**2/D1**2
        Pot2=Pot1*N2**3/N1**3*D2**3/D1**3
        npsh2=npsh1*N2**2/N1**2*D2**2/D1**2 #Esta relación no es tan fiable como las anteriores

        return [D2, N2, Q2, h2, Pot2, npsh2]


    def coste(self):
        HP=self.power.hp
        LnHP = log(self.power.hp)
        Q=self.kwargs["entrada"].Q.galUSmin

        #Coste Bomba
        if self.kwargs["tipo_bomba"]==0:                       #Centrifugal pumps
            QH=log(Q*self.power.hp**0.5)
            Fm=[1., 1.35, 1.15, 2., 2., 3.5, 3.3, 4.95, 4.6, 9.7, 2.95, 1.15, 1.90]

            b1=[0., 5.1029, 0.0632, 2.0290, 13.7321, 9.8849][self.kwargs["tipo_centrifuga"]]
            b2=[0., -1.2217, 0.2744, -0.2371, -2.8304, -1.6164][self.kwargs["tipo_centrifuga"]]
            b3=[0., 0.0771, -0.0253, 0.0102, 0.1542, 0.0834][self.kwargs["tipo_centrifuga"]]

            Ft = exp(b1 + b2 * QH + b3 * QH**2)
            Cb=Fm[self.kwargs["material"]]*Ft*1.55*exp(8.833-0.6019*QH+0.0519*QH**2)

        elif self.kwargs["tipo_bomba"]==1:                     #Reciprocating pumps
            if self.kwargs["material"] == 0:                                #Case iron
                Cb=40.*Q**0.81
            elif self.kwargs["material"] == 3:                             #316 Staineless steel
                Cb = 410.*Q**0.52
            elif self.kwargs["material"] == 12:                           #Bronze
                Cb = 410.* 1.4 * Q**0.52
            elif self.kwargs["material"] == 5:                             #Nickel
                Cb = 410. * 1.86 * Q**0.52
            elif self.kwargs["material"] == 6:                             #Monel
                Cb = 410. * 2.20 * Q**0.52
            else:                                                  #Material not available. Assume case iron
                Cb = 40. * Q**0.81

        elif self.kwargs["tipo_bomba"]==2:                     #Gear pumps
            Cb=1000*exp(-0.0881+0.1986*log(Q)+0.0291*log(Q)**2)
        elif self.kwargs["tipo_bomba"]==3:                     #Vertical mixed flow
            Cb=0.036*Q**0.82*1000
        elif self.kwargs["tipo_bomba"]==4:                     #Vertical axial flow
            Cb=0.02*Q**0.78*1000

        C_bomba=Cb * self.kwargs["Current_index"] / self.kwargs["Base_index"]

        #Coste motor
        if self.kwargs["motor"] == 0:                                     #Open, drip-proof
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 4.8314, 0.0966, 0.10960
            elif self.kwargs["rpm"] == 0 and 7.5< HP <= 250.:
                a1, a2, a3 = 4.1514, 0.5347, 0.05252
            elif self.kwargs["rpm"] == 0 and HP > 250.:
                a1, a2, a3 = 4.2432, 1.03251, -0.03595
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 4.7075, -0.01511, 0.22888
            elif self.kwargs["rpm"] == 1 and 7.5< HP <= 250:
                a1, a2, a3 = 4.5212, 0.47242, 0.04820
            elif self.kwargs["rpm"] == 1 and HP > 250.:
                a1, a2, a3 = 7.4044, -0.06464, 0.05448
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 4.9298, 0.30118, 0.12630
            elif self.kwargs["rpm"] == 2 and 7.5< HP <= 250:
                a1, a2, a3 = 5.0999, 0.35861, 0.06052
            elif self.kwargs["rpm"] == 2 and HP > 250.:
                a1, a2, a3 = 4.6163, 0.88531, -0.02188
        elif self.kwargs["motor"] == 1:                                   #Totally enclosed, fan-cooled
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 5.1058, 0.03316, 0.15374
            elif self.kwargs["rpm"] == 0 and 7.5< HP <= 250.:
                a1, a2, a3 = 3.8544, 0.83311, 0.02399
            elif self.kwargs["rpm"] == 0 and HP > 250.:
                a1, a2, a3 = 5.3182, 1.08470, -0.05695
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 4.9687, -0.00930, 0.22616
            elif self.kwargs["rpm"] == 1 and HP > 7.5:
                a1, a2, a3 = 4.5347, 0.57065, 0.04609
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 5.1532, 0.28931, 0.14357
            elif self.kwargs["rpm"] == 2 and HP > 7.5:
                a1, a2, a3 = 5.3858, 0.31004, 0.07406
        elif self.kwargs["motor"] == 2:                                    #Explosion-proof
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 5.3934, -0.00333, 0.15475
            elif self.kwargs["rpm"] == 0 and HP > 7.5:
                a1, a2, a3 = 4.4442, 0.60820, 0.05202
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 5.2851, 0.00048, 0.19949
            elif self.kwargs["rpm"] == 1 and HP > 7.5:
                a1, a2, a3 = 4.8178, 0.51086, 0.05293
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 5.4166, 0.31216, 0.10573
            elif self.kwargs["rpm"] == 2 and HP > 7.5:
                a1, a2, a3 = 5.5655, 0.31284, 0.07212

        C_motor = 1.2 * exp(a1 + a2 * LnHP + a3 * LnHP**2) * self.kwargs["Current_index"] / self.kwargs["Base_index"]

        self.C_bomba=Currency(C_bomba)
        self.C_motor=Currency(C_motor)
        self.C_adq=Currency(C_bomba+C_motor)
        self.C_inst=Currency(self.C_adq*self.kwargs["f_install"])


    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.kwargs["entrada"].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Head"), self.headCalculada.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Brake horsepower"), self.P_freno.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Volumetric Flow"), self.volflow.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Power"), self.power.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimientoCalculado)+os.linesep
        txt+="%-25s\t %0.4f" %("Cp/Cv", self.kwargs["entrada"].Liquido.cp_cv)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Pump type"), self.TEXT_BOMBA[self.kwargs["tipo_bomba"]])
            if self.kwargs["tipo_bomba"]==0:
                txt+=", "+self.TEXT_CENTRIFUGA[self.kwargs["tipo_centrifuga"]]+os.linesep
            else:
                txt+=os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Material"), self.TEXT_MATERIAL[self.kwargs["material"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pump Cost"), self.C_bomba.str)+os.linesep
            txt+="%-25s\t %s, %s" %(QApplication.translate("pychemqt", "Motor type"), self.TEXT_MOTOR[self.kwargs["motor"]], self.TEXT_RPM[self.kwargs["rpm"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Motor Cost"), self.C_motor.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "PoutCalculada", Pressure),
                (QApplication.translate("pychemqt", "Head"), "headCalculada", Length),
                (QApplication.translate("pychemqt", "Brake horsepower"), "P_freno", Power),
                (QApplication.translate("pychemqt", "Volumetric Flow"), "volflow", VolFlow),
                (QApplication.translate("pychemqt", "Power"), "power", Power),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimientoCalculado", Dimensionless),
                (QApplication.translate("pychemqt", "Pump Type"), ("TEXT_BOMBA", "tipo_bomba"), str),
                (QApplication.translate("pychemqt", "Centrifuge Type"), ("TEXT_CENTRIFUGA", "tipo_centrifuga"), str),
                (QApplication.translate("pychemqt", "Material"), ("TEXT_MATERIAL", "material"), str),
                (QApplication.translate("pychemqt", "Motor Type"), ("TEXT_MOTOR", "motor"), str),
                (QApplication.translate("pychemqt", "Motor RPM"), ("TEXT_RPM", "rpm"), str),
                (QApplication.translate("pychemqt", "Pump Cost"), "C_bomba", Currency),
                (QApplication.translate("pychemqt", "Motor Cost"), "C_motor", Currency),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", Currency)]
        return list

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["deltaP"] = self.deltaP
        state["rendimientoCalculado"] = self.rendimientoCalculado
        state["headCalculada"] = self.headCalculada
        state["power"] = self.power
        state["P_freno"] = self.P_freno
        state["Pin"] = self.Pin
        state["PoutCalculada"] = self.PoutCalculada
        state["volflow"] = self.volflow
        state["statusCoste"] = self.statusCoste

        if self.statusCoste:
            state["C_bomba"] = self.C_bomba
            state["C_motor"] = self.C_motor
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

        if self.kwargs["usarCurva"]:
            pass

    def readStatefromJSON(self, state):
        self.deltaP = DeltaP(state["deltaP"])
        self.rendimientoCalculado = Dimensionless(state["rendimientoCalculado"])
        self.headCalculada = Length(state["headCalculada"])
        self.power = Power(state["power"])
        self.P_freno = Power(state["P_freno"])
        self.Pin = Pressure(state["Pin"])
        self.PoutCalculada = Pressure(state["PoutCalculada"])
        self.volflow = VolFlow(state["volflow"])
        self.statusCoste = state["statusCoste"]

        if self.statusCoste:
            self.C_bomba = Currency(state["C_bomba"])
            self.C_motor = Currency(state["C_motor"])
            self.C_adq = Currency(state["C_adq"])
            self.C_inst = Currency(state["C_inst"])

        if self.kwargs["usarCurva"]:
            pass
        self.salida = [None]

    def datamap2xls(self):
        datamap=(("PoutCalculada", "value", "H15"),
                ("PoutCalculada", "unit", "I15"),
                ("Pin", "value", "H16"),
                ("Pin", "unit", "I16"), )
        return datamap

    def export2pdf(self):
        bomba=pdf("Bomba")
        bomba.bomba(self)
        bomba.dibujar()
        os.system("evince datasheet.pdf")

    def export2xls(self):
        font0 = xlwt.Font()
        font0.bold = True
        font0.height = 300
        print((font0.height))


        style0 = xlwt.XFStyle()
        style0.font = font0

        style1 = xlwt.XFStyle()
        style1.num_format_str = 'D-MMM-YY'

        wb = xlwt.Workbook()
        ws = wb.add_sheet('A Test Sheet')

        ws.write(0, 0, 'Test', style0)
        ws.write(2, 0, 1)
        ws.write(2, 1, 1)
        ws.write(2, 2, xlwt.Formula("A3+B3"))

        wb.save('datasheet.xls')
        os.system("gnumeric datasheet.xls")

if __name__ == '__main__':
#    import doctest
#    doctest.testmod()


#    t=Temperature(27, "C")
#    DeltaP= Pressure(1, "atm")
#    agua=Corriente(t, 1, 3600, Mezcla([62], [1]))
#    bomba=Pump()
#    Q=r_[0, 2, 4, 6, 8, 10, 12, 14, 16, 18]/3600. #en m³/s
#    h=r_[15.5, 15.4, 15.3, 15.1, 14.8, 14.5, 14.1, 13.5, 12.7, 11.6]
#    Pot=r_[0.5, 0.53, 0.57, 0.61, 0.66, 0.71, 0.77, 0.83, 0.89, 0.97]*1000
#    NHPS=r_[0]*9
#    bomba(entrada=agua, usarCurva=1, calculo=1, DeltaP=DeltaP.atm, curvaCaracteristica=[3.5, 1500, Q, h, Pot, NHPS])
##    print bomba.head, "m"
##    print bomba.entrada.caudal_volumetrico.m3h, "m3h"
##    print bomba.rendimiento*100,  "%"
##    print bomba.power.hp, "HP, ",  bomba.power.kW, "kW"
##    print bomba.salida.P.atm,  "atm"
##    print agua.caudal_volumetrico.m3h, bomba.entrada.caudal_volumetrico.m3h
#    bomba.coste(2.8, 421.2, 600)
##    print bomba.C_bomba, bomba.C_motor, bomba.C_adq, bomba.C_inst
##    print bomba.C_inst
#    bomba.export2xls()

    agua=Corriente(T=275, P=101325., ids=[62], caudalMasico=1, fraccionMolar=[1.])
    bomba=Pump(entrada=agua, rendimiento=0.75, deltaP=20*101325, tipo_bomba=1)
#    print bomba.power.hp
#    print bomba.C_inst
#    print bool(bomba)
#    bomba.clear()
#    print bomba.__dict__
#    print bool(bomba)

#    bomba=Pump()
#    print bomba.__class__
