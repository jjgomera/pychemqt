#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# library for definition of equipment with gas, solid and liquid interaction:
#     -Dryer
#     -Scrubber
###############################################################################


from PyQt4.QtGui import QApplication
from scipy import pi, exp, sqrt, log

from lib import unidades
from lib.physics import Cunningham
from lib.corriente import Corriente
from lib.psycrometry import Punto_Psicrometrico
from parents import equipment


class Dryer(equipment):
    """Clase que define un equipo de secado de sólidos

    Parámetros:
        entradaSolido: Instancia de clase Corriente que define la entrada de sólidos
        entradaAire: Instancia de clase Psychrometry que define le aire entrante
        entrada: array con ambas entradas al equipo
        tipoCalculo: integer que indica el tipo de cálculo
            0   -   Cálculo, se calculan la humedad del aire a la salida
            1   -   Diseño, se imponen las condiciones de salida del gas y se calcula el caudal necesario
        HR: Humedad relativa del aire a la salida, por defecto 100%
        TemperaturaSolid: Temperatura del sólido a la salida, por defecto se considerará la temperatura de salida del vapor
        HumedadResidual: Humedad residual del sólido, en kg/kg por defecto se considerará 0
        Heat: Calor intercambiado por el equipo, por defecto se considerará funcionamiento adiabático, heat=0
        deltaP: Perdida de presión del equipo

    """
    title=QApplication.translate("pychemqt", "Solid dryer")
    help=""
    kwargs={"entradaSolido": None,
            "entradaAire": None,
            "entrada": [],
            "tipoCalculo": 0,
            "HR": 0.0,
            "TemperaturaSolid": 0.0,
            "HumedadResidual": 0.0,
            "Heat": 0.0,
            "deltaP": 0.0}

    @property
    def isCalculable(self):
        if not self.kwargs["entradaSolido"]:
            self.msg=QApplication.translate("pychemqt", "undefined solid stream input")
            self.status=0
        elif not self.kwargs["entradaAire"]:
            self.msg=QApplication.translate("pychemqt", "undefined air stream input")
            self.status=0
        elif not self.kwargs["HR"]:
            self.msg=QApplication.translate("pychemqt", "using default air output relative humid 100%")
            self.status=3
            return True
        else:
            self.msg=""
            self.status=1
            return True

    def cleanOldValues(self, **kwargs):
        """Actualización de los kwargs con los nuevos introducidos si es necesario para cada equipo"""
        if kwargs.has_key("entrada"):
            kwargs["entradaSolido"]=kwargs["entrada"][0]
            kwargs["entradaSolido"]=kwargs["entrada"][1]
            del kwargs["entrada"]
        self.kwargs.update(kwargs)


    def calculo(self):
        #TODO: De momento, no se implementan las cuestiones de cinetica de intercambio de calor y materia que definirían las dimensiones necesarias del equipo
        HR=self.kwargs.get("HR", 100)
        self.Heat=unidades.Power(self.kwargs["Heat"])
        self.deltaP=unidades.Pressure(self.kwargs["deltaP"])
        self.entradaAire=self.kwargs["entradaAire"]

        Pout=min(self.kwargs["entradaSolido"].P.atm, self.kwargs["entradaAire"].P.atm)-self.deltaP.atm

        aguaSolidoSalida=self.kwargs["HumedadResidual"]*self.kwargs["entradaSolido"].solido.caudal.kgh
        aguaSolidoEntrada=self.kwargs["entradaSolido"].caudalmasico.kgh
        if self.kwargs["tipoCalculo"]==0:
            Caudal_aguaenAireSalida=aguaSolidoEntrada-aguaSolidoSalida+self.entradaAire.caudal.kgh*self.entradaAire.Xw
            Caudal_airesalida=self.entradaAire.caudal.kgh*self.entradaAire.Xa
            if self.entradaAire.Hs>Caudal_aguaenAireSalida/Caudal_airesalida:
                H=Caudal_aguaenAireSalida/Caudal_airesalida
            else:
                H=self.entradaAire.Hs
                aguaSolidoSalida+=Caudal_aguaenAireSalida/Caudal_airesalida-self.entradaAire.Hs
            self.SalidaAire=Punto_Psicrometrico(caudal=Caudal_aguaenAireSalida+Caudal_airesalida, tdb=self.entradaAire.Tdb, H=H)
            self.SalidaSolido=self.kwargs["entradaSolido"].clone(T=self.SalidaAire.Tdb, P=Pout, split=aguaSolidoSalida/aguaSolidoEntrada)
        else:
            pass


#        if self.HumedadResidual==0:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        else:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        self.SalidaAire=Corriente(self.entradaAire.T, self.entradaAire.P.atm-self.deltaP.atm, self.entradaAire.caudalmasico.kgh, self.entradaAire.mezcla)


class Scrubber(equipment):
    """ Model equipment for scrubber

    Parameters:
        entradaGas: Stream instance of Corriente for define gas input
        entradaLiquido: Stream instance of Corriente for define liquid input
        tipo_calculo:
            0   -   Rating, known dimensions
            1   -   Design, fixed efficiency, calculate dimensions
        diametro: scruber diameter
        rendimiento: Efficiency admisible
        modelo_rendimiento: Scrubber calculate model:
            0   -   Johnstone (1954)
            1   -   Calvert (1972)
        k: Johnstone Method - Empiric constant of scrubber, range 500-1000
        f: Calvert Method - Empiric parameter 0.2 (hydrophobic) to 0.7 (hydrophilic)
        modelo_DeltaP:
            0   -   Young (1977)
    """
    title = QApplication.translate("pychemqt", "Scrubber")
    help = ""
    kwargs = {"entradaGas": None,
              "entradaLiquido": None,
              "tipo_calculo": 0,
              "diametro": 0.0,
              "rendimiento": 0.0,
              "modelo_rendimiento": 0,
              "k": 0.0,
              "f": 0.0,
              "modelo_DeltaP": 0,
              "Lt": 0.0}
    kwargsInput = ("entradaGas", "entradaLiquido")
    kwargsValue = ("diametro", "rendimiento", "k", "Lt")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP")
    calculateValue = ("deltaP", "efficiencyCalc")

    TEXT_TIPO = [QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
                 QApplication.translate("pychemqt", "Design: Calculate diameter")]
    TEXT_MODEL = ["Johnstone (1954)",
                  "Calvert (1972)"]
    TEXT_MODEL_DELTAP = ["Calvert (1968)",
                         "Hesketh (1974)",
                         "Gleason (1971)",
                         "Volgin (1968)",
                         "Young (1977)"]

    __doi__ = [
            {"autor": "S. JENNINGS",
              "title": "The Mean Free Path in Air", 
              "ref": "Journal of Aerosol Science - J AEROSOL SCI 01/1988; 19(2):159-166",
              "doi":  "10.1016/0021-8502(88)90219-4"}, 
            {"autor": "Seymour Calvert, Dale Lundgren & Dilip S. Mehta",
              "title": "Venturi Scrubber Performance", 
              "ref": "Journal of the Air Pollution Control Association, 22:7, 529-532",
              "doi":  "10.1080/00022470.1972.10469674"}, 
            {"autor": "Shui-Chow Yung, Harry F. Barbarika & Seymour Calvert",
              "title": "Pressure Loss in Venturi Scrubbers", 
              "ref": "Journal of the Air Pollution Control Association, 27:4, 348-351",
              "doi":  "10.1080/00022470.1977.10470432"}, 
            {"autor": "Howard E. Hesketh",
              "title": "Fine Particle Collection Efficiency Related to Pressure Drop, Scrubbant and Particle Properties, and Contact Mechanism", 
              "ref": "Journal of the Air Pollution Control Association, 24:10, 939-942",
              "doi":  "10.1080/00022470.1974.10469992"}, 
              ]
            
    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""

        self.statusDeltaP = 1
        if self.kwargs["modelo_DeltaP"] in (3, 4) and not self.kwargs["Lt"]:
            self.statusDeltaP = 0

        if not self.kwargs["entradaGas"]:
            self.msg = QApplication.translate("pychemqt", "undefined gas stream")
            self.status = 0
            return
        if not self.kwargs["entradaLiquido"]:
            self.msg = QApplication.translate("pychemqt", "undefined liquid stream")
            self.status = 0
            return
        if self.kwargs["tipo_calculo"] == 0 and not self.kwargs["diametro"]:
            self.msg = QApplication.translate("pychemqt", "undefined diameter")
            self.status = 0
            return
        elif self.kwargs["tipo_calculo"] == 1 and not self.kwargs["rendimiento"]:
            self.msg = QApplication.translate("pychemqt", "undefined efficiency")
            self.status = 0
            return

        if self.kwargs["modelo_rendimiento"] == 0 and not self.kwargs["k"]:
            self.msg = QApplication.translate("pychemqt", "undefined venturi constant")
            self.status = 3
        elif self.kwargs["modelo_rendimiento"] == 1 and not self.kwargs["f"]:
            self.msg = QApplication.translate("pychemqt", "undefined calvert coefficient")
            self.status = 3

        return True

    def calculo(self):
        self.entradaGas = self.kwargs["entradaGas"]
        self.entradaLiquido = self.kwargs["entradaLiquido"]
        self.Dt = unidades.Length(self.kwargs["diametro"])
        self.Lt = unidades.Length(self.kwargs["Lt"])

        if self.kwargs["k"]:
            self.k = self.kwargs["k"]
        else:
            self.k = 1000.
        if self.kwargs["f"]:
            self.f = self.kwargs["f"]
        else:
            self.f = 0.5

        self.At = unidades.Area(pi/4*self.Dt**2)
        self.Vg = unidades.Speed(self.entradaGas.Q/self.At)
        self.R = self.entradaLiquido.Q/self.entradaGas.Q
        self.dd = 58600*(self.entradaLiquido.Liquido.epsilon/self.entradaLiquido.Liquido.rho)**0.5/self.Vg**2 +\
            597*(self.entradaLiquido.Liquido.mu/self.entradaLiquido.Liquido.epsilon**0.5/self.entradaLiquido.Liquido.rho**0.5)**0.45*(1000*self.R)**1.5

        self.efficiency_i = self._Efficiency()
        self.efficiencyCalc = self._GlobalEfficiency(self.efficiency_i)

        if self.statusDeltaP:
            self.deltaP=self._deltaP()
        else:
            self.deltaP=unidades.DeltaP(0)

        Solido_NoCapturado, Solido_Capturado=self.entradaGas.solido.Separar(self.efficiency_i)
        self.Pin=min(self.entradaGas.P, self.entradaLiquido.P)
        self.Pout=self.Pin-self.deltaP
        self.Min=self.entradaGas.solido.caudal
        self.Dmin=self.entradaGas.solido.diametro_medio
        self.Mr=Solido_NoCapturado.caudal
        self.Dmr=Solido_NoCapturado.diametro_medio
        self.Ms=Solido_Capturado.caudal
        self.Dms=Solido_Capturado.diametro_medio

        self.salida=[]
        self.salida.append(self.entradaGas.clone(solido=Solido_NoCapturado, P=self.Pin-self.deltaP))
        self.salida.append(self.entradaLiquido.clone(solido=Solido_Capturado, P=self.Pin-self.deltaP))

    def _Efficiency(self):
        rendimiento_fraccional=[]
        if self.kwargs["modelo_rendimiento"] == 0:  # Modelo de Johnstone (1954)
            l=sqrt(pi/8)*self.entradaGas.Gas.mu/0.4987445/sqrt(self.entradaGas.Gas.rho*self.entradaGas.P)
            for dp in self.entradaGas.solido.diametros:
                Kn=l/dp*2
                C=Cunningham(l, Kn)
                kp=C*self.entradaGas.solido.rho*dp**2*self.Vg/9/self.entradaGas.Gas.mu/self.dd
                penetration=exp(-self.k*self.R*kp**0.5)
                rendimiento_fraccional.append(1-penetration)

        elif self.kwargs["modelo_rendimiento"] == 1:  # Modelo de Calvert (1972)
            l=sqrt(pi/8)*self.entradaGas.Gas.mu/0.4987445/sqrt(self.entradaGas.Gas.rho*self.entradaGas.P)
            for dp in self.entradaGas.solido.diametros:
                Kn=l/dp*2
                C=Cunningham(l, Kn)
                kp=C*self.entradaGas.solido.rho*dp**2*self.Vg/9/self.entradaGas.Gas.mu/self.dd
                b=(-0.7-kp*self.f+1.4*log((kp*self.f+0.7)/0.7)+0.49/(0.7+kp*self.f))/kp
                penetration=exp(self.R*self.Vg*self.entradaLiquido.Liquido.rho*self.dd/55/self.entradaGas.Gas.mu*b)
                rendimiento_fraccional.append(1-penetration)

        return rendimiento_fraccional

    def _GlobalEfficiency(self, rendimientos):
        rendimiento_global = 0
        for i, fraccion in enumerate(self.entradaGas.solido.fracciones):
            rendimiento_global += rendimientos[i]*fraccion
        return unidades.Dimensionless(rendimiento_global)

#        DeltaP=Pressure(1.002*V**2*R, "kPa")
#
#        self.DeltaP=DeltaP
#
#        print V, DeltaP
        """        1. Calvert: 'dp = 5.4e-04 * (v^2) * rho_gas * (L/G)'
        2a. Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * (0.56 + 0.125*L/G + 0.0023*(L/G)^2) / 507'
        2b. Simplified Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * ((L/G)^0.78) / 1270'
                """
    def _deltaP(self):
        if self.kwargs["modelo_DeltaP"] == 0: #Calvert (1968)
            deltaP=0.85*self.entradaLiquido.Liquido.rho*self.Vg**2*self.R
        elif self.kwargs["modelo_DeltaP"] == 1: #Hesketh (1974)
            deltaP=1.36e-4*self.Vg.cms**2*self.entradaGas.Gas.rho.gcc*self.At.cm2**0.133*(0.56+935*self.R+1.29e-2*self.R**2)
        elif self.kwargs["modelo_DeltaP"] == 2: #Gleason (1971)
            deltaP=2.08e-5*self.Vg.cms**2*(0.264*self.entradaLiquido.Liquido.Q.ccs+73.8)
        elif self.kwargs["modelo_DeltaP"] == 3: #Volgin (1968)
            deltaP=3.32e-6*self.Vg.cms**2*self.R*0.26*self.Lt**1.43
        elif self.kwargs["modelo_DeltaP"] == 4: #Yung (1977)
            Re=self.dd*self.Vg+self.entradaGas.Gas.rho/self.entradaGas.Gas.mu
            Cd=0.22+(24/Re * (1+0.15*Re**0.6))
            X=3*self.Lt*Cd*self.entradaGas.Gas.rho/16/self.dd/self.entradaLiquido.Liquido.rho + 1
            deltaP=2*self.entradaLiquido.Liquido.rho*self.Vg**2*self.R*(1-X**2+(X**4-X**2)**0.5)

#        elif self.kwargs["modelo_DeltaP"] == 0: #Matrozov (1953)
#            deltaP=dPd+1.38e-3*self.Vg.cms**1.08*self.R**0.63
#        elif self.kwargs["modelo_DeltaP"] == 1: #Yoshida (1960)
#            pass
#        elif self.kwargs["modelo_DeltaP"] == 3: #Tohata (1964)
#            pass
#        elif self.kwargs["modelo_DeltaP"] == 4: #Geiseke (1968)
#            pass
#        elif self.kwargs["modelo_DeltaP"] == 8: #Boll (1973)
#            pass
#        elif self.kwargs["modelo_DeltaP"] == 9: #Behie & Beeckman (1973)
#            pass

        return unidades.DeltaP(deltaP)

    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["tipo_calculo"]].split(":")[0])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Model"), self.TEXT_MODEL[self.kwargs["modelo_rendimiento"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Pressure Loss Model"), self.TEXT_MODEL_DELTAP[self.kwargs["modelo_DeltaP"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Ciclon Model"), self.TEXT_MODEL_CICLON[self.kwargs["modelo_ciclon"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Critic Particle Diameter"), self.dc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Internal Cycles"), self.N.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Speed"), self.V.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep
        txt+="%-25s\t %i" %(QApplication.translate("pychemqt", "Cyclone number"), self.num_ciclones)+os.linesep

        txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Cyclone Dimensions")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Ciclon Diameter"), self.Dc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Inlet Height"), self.Hc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Inlet Width"), self.Bc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Diameter"), self.Jc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Cylinder Cyclone Section Length"), self.Lc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Conical Cyclone Section Length"), self.Zc.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Clean Gas Output Diameter"), self.De.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Clean Gas Inlet Orifice Length"), self.Sc.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
        if len(self.entrada.solido.diametros)>=1:
            txt+="%10s %10s %10s %10s %10s" %(u"Dp(µm)", "Xi", u"ηi", "Xis", "Xig")+os.linesep
            for i in range(len(self.rendimiento_parcial)):
                txt+="%10.1f %10.2f %10.3f %10.3f %10.3f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Model"), self.TEXT_COST[self.kwargs["tipo_costo"]])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt


    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "tipo_calculo"), str),
                (QApplication.translate("pychemqt", "Model"), ("TEXT_MODEL", "modelo_rendimiento"), str),
                (QApplication.translate("pychemqt", "Pressure Loss Model"), ("TEXT_MODEL_DELTAP", "modelo_DeltaP"), str),
                (QApplication.translate("pychemqt", "Ciclon Model"), ("TEXT_MODEL_CICLON", "modelo_ciclon"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Critic Particle Diameter"), "dc", Length),
                (QApplication.translate("pychemqt", "Gas Internal Cycles"), "N", Dimensionless),
                (QApplication.translate("pychemqt", "Gas Speed"), "V", Speed),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Cyclone number"), "num_ciclones", Dimensionless),
                (QApplication.translate("pychemqt", "Ciclon Diameter"), "Dc", Length),
                (QApplication.translate("pychemqt", "Inlet Height"), "Hc", Length),
                (QApplication.translate("pychemqt", "Inlet Width"), "Bc", Length),
                (QApplication.translate("pychemqt", "Solid Output Diameter"), "Jc", Length),
                (QApplication.translate("pychemqt", "Cylinder Cyclone Section Length"), "Lc", Length),
                (QApplication.translate("pychemqt", "Conical Cyclone Section Length"), "Zc", Length),
                (QApplication.translate("pychemqt", "Clean Gas Output Diameter"), "De", Length),
                (QApplication.translate("pychemqt", "Clean Gas Inlet Orifice Length"), "Sc", Length),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length),
                (QApplication.translate("pychemqt", "Cost Mode"), ("TEXT_COST", "tipo_costo"), str),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", Currency)]
        return list



if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

    from lib.corriente import Corriente, Solid
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[1/3600.], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    aire=Corriente(T=350, P=101325, caudalMasico=0.01, ids=[475], fraccionMolar=[1.], solido=solido)
    agua=Corriente(T=300, P=101325, caudalMasico=0.1, ids=[62], fraccionMolar=[1.])
    secador=Scrubber(entradaGas=aire, entradaLiquido=agua, modelo_rendimiento=1, diametro=0.25, f=0.5)
    print secador.deltaP.bar, secador.efficiencyCalc
