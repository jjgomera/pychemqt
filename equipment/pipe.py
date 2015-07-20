#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################
#        librería de definición de tuberias       ###
###############################

import os

from PyQt5.QtWidgets import  QApplication
from scipy.constants import g, pi

from lib import unidades
from lib.utilities import representacion
from lib.corriente import Corriente
from lib.physics import f_friccion, Re
from equipment.parents import equipment
from equipment.heatExchanger import Heat_Exchanger


class Pipe(equipment):
    """Clase que modela un una tubería

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        metodo: Método de cálculo
            0   -   Flujo unifásico
            1   -   Flujo de agua (Eq. de Hazen-Williams)
            2   -   Flujo de vapor (Eq de Fritzsche)
            3   -   Flujo de gas isotérmico
            4   -   Flujo bifásico (Método Baker)
            5   -   Flujo bifásico (Método de Beggs and Brill)
        thermal: Método de funcionamiento adiabático:
            0   -   Adiabático
            1   -   Flujo de calor fíjo
            2   -   Calcular el intercambio de calor conocidas U y Tª externa
        Material: Array con los datos del material
            0   -   Nombre
            1   -   Clase
            2   -   Rugosidad, mm
            3   -   Catalogo
            4   -   Di, mm
            5   -   Espesor, mm
            6   -   De, mm
            7   -   Peso, kg/m
            8   -   Volumen, m³/100m
            9   -   Superficie, m²/100m
            10 -    Indice de la lista de materiales
            11 -    Indice de la lista de de dimensiones
        Accesorios: Array en el que cada entrada representa un equipo con perdida de carga:
            Indice tipo equipo
            Indice diametro
            K
            Numero de equipos
            Tipo
            Diametro en mm
            Diametro en pulgadas
            Texto explicativo

        l: Longitud en metros de la tubería
        h: diferencia de altura entre la entrada y la salida de la tubería, m
        K: Coeficiente del los accesorios
        C:Coeficiente para el método de Hazen-Williams
        T_ext: Temperatura en el exterior de la tubería
        U: Coeficiente global de transimisión de calor entre la tubería y el exterior
        Q: Calor transmitido a través de las paredes de la tubería

    Coste:
        solo disponible para las tuberías de acero, Ref Darby pag 217

    >>> agua=agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    >>> tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Cast Iron", "Class A", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06])
    >>> print tuberia.Di, tuberia.V,  tuberia.Re, tuberia.DeltaP
    0.1524 0.0626814744044 9427.99142792 1.8362005711
    """
    title=QApplication.translate("pychemqt", "Pipe")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "thermal": 0,
                    "material": [],
                    "accesorios": [],
                    "l": 0.0,
                    "h": 0.0,
                    "K": 0.0,
                    "C": 100.,
                    "T_ext": 0.0,
                    "U": 0.0,
                    "Q": 0.0,

                    "f_install": 2.8,
                    "Base_index": 0.0,
                    "Current_index": 0.0}
    kwargsInput=("entrada", )
    kwargsValue=("l", "h", "C")
    kwargsList=("metodo", "thermal")
    calculateValue=("DeltaP", "DeltaP_f", "DeltaP_ac", "DeltaP_h", "DeltaP_v", "DeltaP_100ft", "V", "f", "Re", "Tout")
    calculateCostos=("C_adq", "C_inst")
    indiceCostos=5

    TEXT_METODO=[QApplication.translate("pychemqt", "Single Phase flow"),
                                QApplication.translate("pychemqt", "Water (Hazen-Williams)"),
                                QApplication.translate("pychemqt", "Steam (Fritzsche)"),
                                QApplication.translate("pychemqt", "Isotermic gas flow"),
                                QApplication.translate("pychemqt", "Two Phase flow (Baker method)"),
                                QApplication.translate("pychemqt", "Two Phase flow (Beggs and Brill method)")]
    TEXT_THERMAL=[QApplication.translate("pychemqt", "Adiabatic"),
                                QApplication.translate("pychemqt", "Heat flux"),
                                QApplication.translate("pychemqt", "Heat transfer")]

    C_adq=unidades.Currency(None)
    C_inst=unidades.Currency(None)

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and self.kwargs["Current_index"] and self.kwargs["material"] and self.kwargs["material"][0] in ['Stainless Steel (ANSI)', 'Steel Galvanised (ANSI)', 'Steel (ANSI)']:
            self.statusCoste=True
        else:
            self.statusCoste=False
            self.C_adq=unidades.Currency(None)
            self.C_inst=unidades.Currency(None)


        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return
        if not self.kwargs["l"]:
            self.msg=QApplication.translate("pychemqt", "undefined pipe length")
            self.status=0
            return
        if not self.kwargs["material"]:
            self.msg=QApplication.translate("pychemqt", "undefined material")
            self.status=0
            return

        if self.kwargs["thermal"]==1 and not self.kwargs["Q"]:
            self.msg=QApplication.translate("pychemqt", "undefined heat flux")
            self.status=0
            return
        elif self.kwargs["thermal"]==2 and (not self.kwargs["T_ext"]  or  not self.kwargs["U"]):
            self.msg=QApplication.translate("pychemqt", "undefined heat transfer conditions")
            self.status=0
            return

        if self.kwargs["metodo"]==1 and not self.kwargs["C"]:
            self.msg=QApplication.translate("pychemqt", "undefined C William Factor")
            self.status=0
            return

        self.msg=""
        self.status=1
        return True


    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.L=unidades.Length(self.kwargs["l"])

        if self.entrada.x==0:
            self.rho=self.entrada.Liquido.rho
            self.mu=self.entrada.Liquido.mu
        else:
            self.rho=self.entrada.Gas.rho
            self.mu=self.entrada.Gas.mu

        self.material=self.kwargs["material"][0] + " " + self.kwargs["material"][1]
        self.Dn=self.kwargs["material"][3]
        self.rugosidad=unidades.Length(self.kwargs["material"][2], "mm")
        self.De=unidades.Length(self.kwargs["material"][6], "mm")
        self.w=unidades.Length(self.kwargs["material"][5], "mm")
        self.Di=unidades.Length((self.De-2*self.w))
        self.eD=unidades.Dimensionless(self.rugosidad/self.Di)
        self.seccion=unidades.Area(pi/4*self.Di**2)
        self.A=unidades.Area(pi*self.De*self.L)
        self.V=unidades.Speed(self.entrada.Q/self.seccion)
        self.Re=Re(self.Di, self.V, self.rho, self.mu)
        K=0
        for accesorio in self.kwargs["accesorios"]:
            K+=accesorio[2]*accesorio[3]
        self.K=unidades.Dimensionless(K)
        self.DeltaP_h=unidades.Pressure(g*self.kwargs["h"]*self.rho)
        self.DeltaP_ac=unidades.Pressure(self.K*self.V**2/2*self.rho)

        self.f=f_friccion(self.Re, self.eD)
        self.DeltaP_f=self.__DeltaP_friccion()
        #TODO:
        self.DeltaP_v=unidades.Pressure(0)

        self.DeltaP=unidades.Pressure(self.DeltaP_f+self.DeltaP_ac+self.DeltaP_h)
        self.DeltaP_100ft=self.DeltaP*100/self.L.ft
        self.Pout=unidades.Pressure(self.entrada.P-self.DeltaP)

        if self.kwargs["thermal"]==0:
            self.Tout=self.entrada.T
            self.Heat=unidades.Power(0)
        else:
            cambiador=Heat_Exchanger()
            cambiador.calculo(entrada=self.entrada, modo=self.kwargs["thermal"], Heat=self.kwargs["Q"], deltaP=self.DeltaP, A=self.A, U=self.kwargs["U"], Text=self.kwargs["Text"])
            self.Tout=cambiador.salida[0].T
            self.Heat=cambiador.Heat

        self.salida=[self.entrada.clone(T=self.Tout, P=self.Pout)]
        self.Pin=self.entrada.P
        self.Pout=self.salida[0].P


    def __DeltaP_friccion(self):
        """Método para el calculo de la perdida de presión"""
        if self.kwargs["metodo"]==0:
            delta=unidades.Pressure(self.L*self.V**2/self.Di*self.f*self.rho/2)
        elif self.kwargs["metodo"]==1:
            delta=unidades.Pressure((self.entrada.Q.galUSmin*self.L.ft**0.54/0.442/self.Di.inch**2.63/self.kwargs["C"])**(1./0.54), "psi")
        elif self.kwargs["metodo"]==2:
            delta=unidades.Pressure(2.1082*self.L.ft*self.entrada.caudalmasico.lbh**1.85/self.rho.lbft3/1e7/self.Di.inch**4.97, "psi")
        elif self.kwargs["metodo"]==3:
            pass

        elif self.kwargs["metodo"]==4:
            pass

        elif self.kwargs["metodo"]==5:
            pass

        return delta


    def coste(self):
        """
        Coste solo disponible para las tuberías de acero
        Ref Darby pag 217

        kwargs:
            schedule: Clase de acero
        """
        codigo=str(self.kwargs["material"][1])
        if codigo in ('Sch. 40', 'Sch. 5S'):
            a=30.
            p=1.31
        elif codigo in ('Sch. 80', 'Sch. 10S'):
            a=38.1
            p=1.35
        elif codigo in ('Sch. 160', 'Sch.  40S'):
            a=55.3
            p=1.39
        else:
            a=0
            p=1

        self.C_adq=unidades.Currency(a*self.Di.ft**p*self.L * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst=unidades.Currency(self.C_adq*self.kwargs["f_install"])


    def writeListtoStream(self, stream, key, value):
        """Personalizar en el caso de equipos con listas complejas"""
        if key=="material":
            stream.writeString(value[0].encode())
            stream.writeString(value[1].encode())
            stream.writeFloat(value[2])
            stream.writeString(value[3].encode())
            for val in value[4:-2]:
                stream.writeFloat(val)
            for val in value[-2:]:
                stream.writeInt32(val)
        elif key=="accesorios":
            stream.writeInt32(len(value))
            for accesorio in value:
                stream.writeInt32(accesorio[0])
                stream.writeInt32(accesorio[1])
                stream.writeFloat(accesorio[2])
                stream.writeInt32(accesorio[3])
                for cadena in accesorio[4:]:
                    stream.writeString(cadena.encode())


    def readListFromStream(self, stream, key):
        """Personalizar en el caso de equipos con listas complejas"""
        valor=[]
        if key=="material":
            valor.append(stream.readString().decode("utf-8"))
            valor.append(stream.readString().decode("utf-8"))
            valor.append(float(representacion(stream.readFloat())))
            valor.append(stream.readString().decode("utf-8"))
            for i in range(6):
                valor.append(float(representacion(stream.readFloat())))
            valor.append(stream.readInt32())
            valor.append(stream.readInt32())
        elif key=="accesorios":
            for i in range(stream.readInt32()):
                accesorio=[]
                accesorio.append(stream.readInt32())
                accesorio.append(stream.readInt32())
                accesorio.append(float(representacion(stream.readFloat())))
                accesorio.append(stream.readInt32())
                for j in range(4):
                    accesorio.append(stream.readString().decode("utf-8"))
                valor.append(accesorio)
        return valor


    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Catalog")+"-----------------#"+os.linesep
        txt+="%-25s\t %s %s" %(QApplication.translate("pychemqt", "Material"), self.kwargs["material"][0], self.kwargs["material"][1])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Nominal Diameter"), self.kwargs["material"][3])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Length"), self.L.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Roughness"), self.rugosidad.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Internal Diamter"), self.Di.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Diamter"), self.De.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thickness"), self.w.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Transversal section"), self.seccion.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Area"), self.A.str)+os.linesep

        if self.kwargs["accesorios"]:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Fittings")+"-----------------#"+os.linesep
            txt+="%-25s\t %s" %("K "+QApplication.translate("pychemqt", "Total"), self.K)+os.linesep
            for accesorio in self.kwargs["accesorios"]:
                txt+="%5i %-22s\t %s" %(accesorio[3], accesorio[7], accesorio[2])+os.linesep

        txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Method"), self.TEXT_METODO[self.kwargs["metodo"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP Total", None), self.DeltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP friction", None), self.DeltaP_f.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP fittings", None), self.DeltaP_ac.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP elevation", None), self.DeltaP_h.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP acceleration", None), self.DeltaP_v.str)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Thermal Condition"), self.TEXT_THERMAL[self.kwargs["thermal"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Fluid Speed"), self.V.str)+os.linesep
        txt+="%-25s\t %s" %("Reynolds", self.Re)+os.linesep
        txt+="%-25s\t %s" %("ε/D", self.eD)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Factor Friction"), self.f)+os.linesep

        if self.kwargs["thermal"]:
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.Tout.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Heat Transfer"), self.Heat.str)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt


    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Material"), "material", str),
                (QApplication.translate("pychemqt", "Nominal Diameter"), "Dn", str),
                (QApplication.translate("pychemqt", "Length"), "L", unidades.Length),
                (QApplication.translate("pychemqt", "Roughness"), "rugosidad", unidades.Length),
                (QApplication.translate("pychemqt", "Internal Diamter"), "Di", unidades.Length),
                (QApplication.translate("pychemqt", "External Diamter"), "De", unidades.Length),
                (QApplication.translate("pychemqt", "Thickness"), "w", unidades.Length),
                (QApplication.translate("pychemqt", "Transversal section"), "seccion", unidades.Area),
                (QApplication.translate("pychemqt", "External Area"), "A", unidades.Area),
                (QApplication.translate("pychemqt", "K total"), "K", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Fittings"), "accesorios", None),
                (QApplication.translate("pychemqt", "Method"), ("TEXT_METODO", "metodo"),  str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", unidades.Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", unidades.Pressure),
                (QApplication.translate("pychemqt", "ΔP Total", None), "DeltaP", unidades.DeltaP),
                (QApplication.translate("pychemqt", "ΔP friction", None), "DeltaP_f", unidades.DeltaP),
                (QApplication.translate("pychemqt", "ΔP fittings", None), "DeltaP_ac", unidades.DeltaP),
                (QApplication.translate("pychemqt", "ΔP elevation", None), "DeltaP_h", unidades.DeltaP),
                (QApplication.translate("pychemqt", "ΔP acceleration", None), "DeltaP_v", unidades.DeltaP),
                (QApplication.translate("pychemqt", "Thermal Condition"), ("TEXT_THERMAL", "thermal"),  str),
                (QApplication.translate("pychemqt", "Fluid Speed"), "V", unidades.Speed),
                (QApplication.translate("pychemqt", "Reynolds number"), "Re", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Relative roughness"), "eD", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Factor Friction"), "f", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Output Temperaturet"), "Tout", unidades.Temperature),
                (QApplication.translate("pychemqt", "Heat Transfer"), "Heat", unidades.Power),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", unidades.Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", unidades.Currency)]
        return list

    def propertiesListTitle(self, index):
        lista=[]
        for accesorio in self.kwargs["accesorios"]:
            lista.append("%3i %s" %(accesorio[3], accesorio[7]))
        return lista




if __name__ == '__main__':
#    import doctest
#    doctest.testmod()


#    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
#    tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Cast Iron", "Class A", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, 0, 2])
#    print tuberia.Di, tuberia.V,  tuberia.Re, tuberia.DeltaP

    agua=Corriente(T=300, P=101325, caudalMasico=1, fraccionMolar=[1.])
    tuberia=Pipe(entrada=agua, metodo=0, l=5, material=["Steel (ANSI)", "Sch. 40", 0.12192, '6"', 152.4, 11.43, 175.26, 42.07, 1.824, 55.06, -1, 2], notas="Tuberia")
    print((tuberia.DeltaP))
