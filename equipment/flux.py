#!/usr/bin/python
# -*- coding: utf-8 -*-

###Modulo que define los equipos de control de flujo, mezcladores, divisores...

import os

from PyQt5.QtWidgets import  QApplication
from scipy.optimize import fsolve

from lib.corriente import Corriente
from lib import unidades
from .parents import equipment


class Divider(equipment):
    """Clase que define un divisor
    
    Parámetros
        entrada: instancia de clase corriente que define la corriente de entrada al divisor
        salidas: numero de corrientes de salida
        criterio: criterio de división
            0   -   fracciones, se indican las fracciones con respecto a la corriente de entrada que salen por cada corriente,
            1   -   caudales, kg/h. Se sobreescribirá el caudal de la corriente de entrada.
        fracciones: array con las fracciones del caudal que irán por cada una de las corrientes de salida, indicará
        deltaP: caída de presión en la unidad, opcional, en atm

    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMasica=[1, 0, 0, 0])
    >>> divisor=Divider(entrada=agua, criterio=0, salidas=3, split=[0.7, 0.2, 0.1])
    >>> print divisor.entrada.caudalmasico.kgh
    3600.0
    >>> for salida in divisor.salida: print salida.caudalmasico.kgh, 
    2520.0 720.0 360.0
    """
    title=QApplication.translate("pychemqt", "Divider")
    help=""
    kwargs={"entrada": None, 
                    "criterio": 0, 
                    "salidas": 0, 
                    "split": [], 
                    "deltaP": 0.0}
    kwargsInput=("entrada", )
    kwargsList=("criterio", )
    kwargsValue=("deltaP", )

    TEXT_CRITERIO=[QApplication.translate("pychemqt", "Flux ratio"), 
                            QApplication.translate("pychemqt", "Flowrate (overwrite input flow)")]

    @property
    def isCalculable(self):
        if self.kwargs["salidas"]==1:
            self.kwargs["split"]=[1.]
            
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
        elif not self.kwargs["split"]:
            self.msg=QApplication.translate("pychemqt", "undefined split fraction")
            self.status=0
        elif self.kwargs["salidas"]!=len(self.kwargs["split"]):
            self.msg=QApplication.translate("pychemqt", "incomplete split fraction")
            self.status=0
        else:
            self.status=1
            self.msg=""
            return True
        
        
    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.criterio=self.kwargs["criterio"]
        self.split=[unidades.Dimensionless(i) for i in self.kwargs["split"]]
        self.deltaP=unidades.DeltaP(self.kwargs["deltaP"])

        #Comprobamos que las fracciones están normalizadas:
        if self.criterio==0:
            if sum(self.split)!=1:
                fracciones_normalizadas=[]
                total=sum(self.split)
                for i in self.split:
                    fracciones_normalizadas.append(i/total)
                self.split=fracciones_normalizadas
        
        self.salida=[]
        if self.criterio==0:
            for i in self.split:
                self.salida.append(self.entrada.clone(P=self.entrada.P-self.deltaP, split=i))
        else:
            self.entrada=self.entrada.clone(caudalmasico=sum(self.split))
            for i in self.split:
                self.salida.append(self.entrada.clone(P=self.entrada.P-self.deltaP, split=i))
                
        self.inputMolarFlow=self.entrada.caudalmolar
        self.inputMassFlow=self.entrada.caudalmasico
        self.inputVolFlow=self.entrada.Q
        self.inputT=self.entrada.T
        self.inputP=self.entrada.P
        self.output=unidades.Dimensionless(self.kwargs["salidas"])

    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input stream"), self.kwargs["entrada"].caudalmasico.str)+os.linesep
        for i, salida in enumerate(self.salida):
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output stream")+str(i), salida.caudalmasico.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure loss"), self.deltaP.str)+os.linesep
        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Feed Molar Flow"), "inputMolarFlow", unidades.MolarFlow), 
                (QApplication.translate("pychemqt", "Feed Mass Flow"), "inputMassFlow", unidades.MassFlow), 
                (QApplication.translate("pychemqt", "Feed Volumetric Flow"), "inputVolFlow", unidades.VolFlow), 
                (QApplication.translate("pychemqt", "Feed Temperature"), "inputT", unidades.Temperature), 
                (QApplication.translate("pychemqt", "Feed Pressure"), "inputP", unidades.Pressure), 
                (QApplication.translate("pychemqt", "Number of Product Streams"), "output", unidades.Dimensionless), 
                (QApplication.translate("pychemqt", "Split Fractions"), "split", unidades.Dimensionless), 
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", unidades.DeltaP)] 
        return list

    def propertiesListTitle(self, index):
        lista=[]
        for i in range(self.kwargs["salidas"]):
            lista.append("Output %i" %(i+1))
        return lista
    
    def writeStatetoStream(self, stream):
        stream.writeInt32(self.criterio)
        stream.writeInt32(len(self.split))
        for value in self.split:
            stream.writeFloat(value)
        stream.writeFloat(self.deltaP)
        stream.writeFloat(self.inputMolarFlow)
        stream.writeFloat(self.inputMassFlow)
        stream.writeFloat(self.inputVolFlow)
        stream.writeFloat(self.inputT)
        stream.writeFloat(self.inputP)
        stream.writeInt32(self.output)

    def readStatefromStream(self, stream):
        self.criterio = stream.readInt32()
        self.split = []
        for i in range(stream.readInt32()):
            self.split.append(unidades.Dimensionless(stream.readFloat()))
        self.deltaP = unidades.DeltaP(stream.readFloat())
        self.inputMolarFlow = unidades.MolarFlow(stream.readFloat())
        self.inputMassFlow = unidades.MassFlow(stream.readFloat())
        self.inputVolFlow = unidades.VolFlow(stream.readFloat())
        self.inputT = unidades.Temperature(stream.readFloat())
        self.inputP = unidades.Pressure(stream.readFloat())
        self.output = unidades.Dimensionless(stream.readInt32())
        self.salida = [None]*self.kwargs["salidas"]


class Mixer(equipment):
    """Clase que define un mezclador

    Parámetros
        n_entradas: Numero de entradas
        entrada: instancia de clase corriente que define la corriente de entrada al mezclador, si se introduce en forma de lista indicará todas las entradas
        id_entrada: id de la entrada a definida
        criterio: criterio de definición de la presión de salida
            0   -   Presión mínima de las corrientes de entrada
            1   -   Presión media de las corrientes de entrada
            2   -   Presión especificada
        Pout: Presión de la corriente de salida

    >>> agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMasica=[1., 0, 0, 0])
    >>> agua2=Corriente(T=350, P=101325, caudalMasico=5, fraccionMasica=[1., 0, 0, 0])
    >>> mezclador=Mixer(entrada=[agua, agua2], criterio=0)
    >>> print mezclador.salida[0].T
    341.376858051
    """
    title=QApplication.translate("pychemqt", "Mixer")  
    help=""
    kwargs={"entrada": [], 
                    "id_entrada": 0, 
                    "criterio": 0, 
                    "Pout": 0.0}
    kwargs_forbidden=["entrada", "id_entrada"]
    kwargsValue=("Pout", )
    kwargsList=("criterio", )

    TEXT_METODO=[QApplication.translate("pychemqt", "Inputs minimum pressure"), 
                                QApplication.translate("pychemqt", "Inputs mean pressure"), 
                                QApplication.translate("pychemqt", "Custom")]

    @property
    def isCalculable(self):
        if self.kwargs["criterio"]==2 and not self.kwargs["Pout"]:
            self.msg=QApplication.translate("pychemqt", "pressure output not defined")
            self.status=0
        elif not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
        elif sum([s.status for s in self.kwargs["entrada"]])==0:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
        elif len(self.kwargs["entrada"])!=sum([s.status for s in self.kwargs["entrada"]]):
            self.msg=QApplication.translate("pychemqt", "some input stream isn't defined")
            self.status=3
            return True
        else: 
            self.msg=""
            self.status=1
            return True

    def cleanOldValues(self, **kwargs):
        if "entrada" in kwargs:
            if isinstance(kwargs["entrada"], list):
                kwargs["id_entrada"]=None
            else:
                corriente=kwargs["entrada"]
                kwargs["entrada"]=self.kwargs["entrada"][:]
                while len(kwargs["entrada"])<kwargs["id_entrada"]+1:
                    kwargs["entrada"].append(Corriente())
                kwargs["entrada"][kwargs["id_entrada"]]=corriente
                kwargs["id_entrada"]=None
        self.kwargs.update(kwargs)

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.criterio=self.kwargs["criterio"]
        Pout=self.kwargs["Pout"]

        if self.criterio==2:
            self.Pout=unidades.Pressure(Pout)
        else:
            lst=[]
            for entrada in self.entrada:
                if entrada.status:
                    lst.append(entrada.P)
            if self.criterio==0:
                self.Pout=unidades.Pressure(min(lst))
            else:
                self.Pout=unidades.Pressure(sum(lst, 0.0) / len(lst))

        h_in=0
        To=0
        caudalunitariomasico=[0]*len(self.entrada[0].fraccion)
        for entrada in self.entrada:
            if entrada.status:
                h_in+=entrada.h
                To+=entrada.T*entrada.caudalmasico
                for i, caudal in enumerate(entrada.caudalunitariomasico):
                    caudalunitariomasico[i]+=caudal
        To/=sum(caudalunitariomasico)
        
        if self.entrada[0].Config.get("Components", "Solids"):
            #TODO: Add solid mixer
            pass
        
        f=lambda T: Corriente(T=T, P=self.Pout, caudalUnitarioMasico=caudalunitariomasico).h-h_in
        T=fsolve(f, To)[0]
        salida=Corriente(T=T, P=self.Pout, caudalUnitarioMasico=caudalunitariomasico)
        self.salida=[salida]
        
        self.outT=salida.T
        self.outP=salida.P
        self.outX=salida.x
        self.outMolarFlow=salida.caudalmolar
        self.outMassFlow=salida.caudalmasico
        self.outVolFlow=salida.Q

    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Output Vapor Fraction"), self.salida[0].x)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Mass Flow"), self.salida[0].caudalmolar.str)+os.linesep
        txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Output Molar Composition")+os.linesep
        for componente, fraccion in zip(self.salida[0].componente, self.salida[0].fraccion):
            txt+="%-25s\t %0.4f" %(componente.nombre, fraccion)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Output Temperature"), "outT", unidades.Temperature), 
                (QApplication.translate("pychemqt", "Output Pressure"), "outP", unidades.Pressure), 
                (QApplication.translate("pychemqt", "Output vapor fraction"), "outX", unidades.Dimensionless), 
                (QApplication.translate("pychemqt", "Output Molar Flow"), "outMolarFlow", unidades.MolarFlow), 
                (QApplication.translate("pychemqt", "Output Mass Flow"), "outMassFlow", unidades.MassFlow), 
                (QApplication.translate("pychemqt", "Output Volumetric Flow"), "outVolFlow", unidades.VolFlow)]
        return list

    def writeStatetoStream(self, stream):
        stream.writeInt32(self.criterio)
        stream.writeFloat(self.Pout)
        stream.writeFloat(self.outT)
        stream.writeFloat(self.outP)
        stream.writeFloat(self.outX)
        stream.writeFloat(self.outMolarFlow)
        stream.writeFloat(self.outMassFlow)
        stream.writeInt32(self.outVolFlow)

    def readStatefromStream(self, stream):
        self.criterio = stream.readInt32()
        self.Pout = unidades.Pressure(stream.readFloat())
        self.outT = unidades.Temperature(stream.readFloat())
        self.outP = unidades.Pressure(stream.readFloat())
        self.outX = unidades.Dimensionless(stream.readFloat())
        self.outMolarFlow = unidades.MolarFlow(stream.readFloat())
        self.outMassFlow = unidades.MassFlow(stream.readFloat())
        self.outVolFlow = unidades.VolFlow(stream.readInt32())
        self.salida = [None]


class Valve(equipment):
    """Clase que define una válvula
    
    Parámetros:
        entrada: instancia de clase corriente que define la corriente de entrada a la vávula
        off: Estado de la válvula
            0   -   abierta
            1   -   semiabierta
            2   -   cerrada
        Pout: Presión a la salida, atm
        DeltaP: caída de presión en la unidad, atm
        Dew: Temperatura punto rocio, K
        Bubble: Temperatura punto burbuja, K

    >>> agua=Corriente(T=300, P=202650, caudalMasico=1, fraccionMasica=[1, 0, 0, 0])
    >>> valvula=Valve(entrada=agua, Pout=101325, off=1)
    >>> print agua.P.atm, valvula.salida[0].P.atm
    2.0 1.0
    >>> valvula(DeltaP=2650)
    >>> print valvula.salida[0].P.atm
    1.97384653343
    >>> valvula(off=0)
    >>> print valvula.salida[0].P.atm
    2.0
    """
    title=QApplication.translate("pychemqt", "Valve")  
    help=""
    kwargs={"entrada": None, 
                    "off": 0, 
                    "Pout": 0.0, 
                    "DeltaP": 0.0, 
                    "Dew": 0.0, 
                    "Bubble": 0.0}
    kwargsInput=("entrada", )
    kwargsValue=("Pout", "DeltaP", "Dew", "Bubble")
    kwargsList=("off", )

    TEXT_WORKING=[QApplication.translate("pychemqt", "Totally open"), 
                                  QApplication.translate("pychemqt", "Partially open"), 
                                  QApplication.translate("pychemqt", "Close")]
    
    @property
    def isCalculable(self):
        if self.kwargs["off"]==1:
            if not self.kwargs["entrada"]:
                self.msg=QApplication.translate("pychemqt", "undefined input")
                self.status=0
            elif self.kwargs["Pout"] or self.kwargs["DeltaP"] or self.kwargs["Dew"] or self.kwargs["Bubble"]:
                self.status=1
                self.msg=""
                return True
            else:
                self.msg=QApplication.translate("pychemqt", "undefined exit condition")
                self.status=0
        elif self.kwargs["off"]==2:
            self.msg=""
            self.status=1
            return True
        else: 
            if not self.kwargs["entrada"]:
                self.msg=QApplication.translate("pychemqt", "undefined input")
                self.status=0
            else:
                self.msg=""
                self.status=1
                return True

    def cleanOldValues(self, **kwargs):
        if "Pout" in kwargs:
            self.kwargs["DeltaP"]=0
            self.kwargs["Dew"]=0
            self.kwargs["Bubble"]=0
        elif "DeltaP" in kwargs:
            self.kwargs["Pout"]=0
            self.kwargs["Dew"]=0
            self.kwargs["Bubble"]=0
        elif "Dew" in kwargs:
            self.kwargs["Pout"]=0
            self.kwargs["DeltaP"]=0
            self.kwargs["Bubble"]=0
        elif "Bubble" in kwargs:
            self.kwargs["Pout"]=0
            self.kwargs["DeltaP"]=0
            self.kwargs["Dew"]=0
        self.kwargs.update(kwargs)

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        Pout=self.kwargs["Pout"]
        DeltaP=self.kwargs["DeltaP"]
        Dew=self.kwargs["Dew"]
        Bubble=self.kwargs["Bubble"]
        
        if self.kwargs["off"]==1:
            if Pout:
                self.Pout=unidades.Pressure(Pout)
            elif DeltaP:
                self.Pout=unidades.Pressure(self.entrada.P-DeltaP)
            elif Dew:
                corriente=self.entrada.clone(T=Dew)
                self.Pout=corriente.eos._Dew_P()
            elif Bubble:
                corriente=self.entrada.clone(T=Bubble)
                self.Pout=corriente.eos._Bubble_P()
            self.salida=[self.entrada.clone(P=self.Pout)]
            
        elif self.kwargs["off"]==2:
            self.entrada=Corriente()
            self.salida=[self.entrada]
    
        else:
            self.salida=[self.entrada]
    
        self.outT=self.salida[0].T
        self.outP=self.salida[0].P
        self.outX=self.salida[0].x


    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Output Vapor Fraction"), self.salida[0].x)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Working Condition"), self.TEXT_WORKING[self.kwargs["off"]])+os.linesep
        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Output Temperature"), "outT", unidades.Temperature), 
                (QApplication.translate("pychemqt", "Output Pressure"), "outP", unidades.Pressure), 
                (QApplication.translate("pychemqt", "Output vapor fraction"), "outX", unidades.Dimensionless), 
                (QApplication.translate("pychemqt", "Working Condition"), ("TEXT_WORKING", "off"), str)]
        return list

    def writeStatetoStream(self, stream):
        stream.writeFloat(self.Pout)
        stream.writeFloat(self.outT)
        stream.writeFloat(self.outP)
        stream.writeFloat(self.outX)

    def readStatefromStream(self, stream):
        self.Pout = unidades.Pressure(stream.readFloat())
        self.outT = unidades.Temperature(stream.readFloat())
        self.outP = unidades.Pressure(stream.readFloat())
        self.outX = unidades.Dimensionless(stream.readFloat())
        self.salida = [None]


if __name__ == '__main__':
#    import doctest
#    doctest.testmod()  
    
    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMasica=[1.])
#    agua2=Corriente(T=300, P=101325*2, caudalMasico=2, fraccionMasica=[1.])
    mezclador=Mixer(entrada=[agua, Corriente()], criterio=0)
#    print mezclador.status, mezclador.msg
    print((mezclador.salida[0].kwargs))
#    print mezclador.salida[0].caudalmasico, 
#    agua3=Corriente(T=300, P=101325, caudalMasico=4, fraccionMasica=[1., 0, 0, 0])
#    mezclador(id_entrada=2, entrada=agua3)
#    print mezclador.salida[0].caudalmasico

#    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[0.3, 0.2, 0.25, 0.25])
#    divisor=Divider(entrada=agua, salidas=3)
#    divisor(split=[0.3, 0.45, 0.25])
#    print divisor.status, divisor.kwargs["salidas"], divisor.msg
#
#    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMasica=[1., 0, 0, 0])
#    valvula=Valve(entrada=agua)
#    print valvula.status, valvula.msg, valvula.salida[0].P.atm
#    valvula(off=1, Pout=100000)
#    print valvula.status, valvula.msg, valvula.salida[0].P.atm
