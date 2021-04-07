#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# library for definition of equipment with gas, solid and liquid interaction:
#     -Scrubber
#     -Dryer
###############################################################################


from math import pi, exp, sqrt, log
import os

from PyQt5.QtWidgets import QApplication

from lib.unidades import (Pressure, DeltaP, Area, Speed, Dimensionless,
                          Length, Power)
from lib.physics import Cunningham
from lib.corriente import Corriente
from lib.psycrometry import PsychroState
from equipment.parents import equipment
from equipment.gas_solid import Separador_SolidGas


class Scrubber(Separador_SolidGas):
    """Class to model a scrubber equipment

    Parameters:
        entradaGas: Corriente instance for define the gas input to equipment
        entradaLiquido: Corriente instance for define the liquid input
        tipo_calculo:
            0 - Rating, known dimensions
            1 - Design, fixed efficiency, calculate dimensions
        diametro: Scruber diameter
        rendimientoAdmisible: Efficiency admisible
        modelo_rendimiento: Scrubber calculate model:
            0 - Johnstone (1954)
            1 - Calvert (1972)
        k: Empiric constant of scrubber for Johnstone Method, range 500-1000
        f: Empiric constant of scrubber for Calvert Method, range from
            0.2 (hydrophobic) to 0.7 (hydrophilic)
        modelo_DeltaP:
            0 - Young (1977)
    """
    title = QApplication.translate("pychemqt", "Scrubber")
    help = ""
    kwargs = {"entradaGas": None,
              "entradaLiquido": None,
              "tipo_calculo": 0,
              "diametro": 0.0,
              "rendimientoAdmisible": 0.0,
              "modelo_rendimiento": 0,
              "k": 0.0,
              "f": 0.0,
              "modelo_DeltaP": 0,
              "Lt": 0.0}
    kwargsInput = ("entradaGas", "entradaLiquido")
    kwargsValue = ("diametro", "rendimientoAdmisible", "k", "Lt")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP")
    calculateValue = ("deltaP", "rendimiento")

    TEXT_TIPO = [
        QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
        QApplication.translate("pychemqt", "Design: Calculate diameter")]
    TEXT_MODEL = ["Johnstone (1954)",
                  "Calvert (1972)"]
    TEXT_MODEL_DELTAP = ["Calvert (1968)",
                         "Hesketh (1974)",
                         "Gleason (1971)",
                         "Volgin (1968)",
                         "Young (1977)"]

    __doi__ = [
        {"autor": "Jennings, S.G.",
         "title": "The Mean Free Path in Air",
         "ref": "J. Aerosol Sci. 19(2) (1988) 159-166"
                "19(2):159-166",
         "doi":  "10.1016/0021-8502(88)90219-4"},
        {"autor": "Calvert, S., Lundgren, D., Mehta, D.S.",
         "title": "Venturi Scrubber Performance",
         "ref": "J. Air Pollution Control Assoc., 22(7) (1972) 529-532",
         "doi":  "10.1080/00022470.1972.10469674"},
        {"autor": "Yung, S.-C., Barbarika, H.F., Calvert, S.",
         "title": "Pressure Loss in Venturi Scrubbers",
         "ref": "J. Air Pollution Control Assoc., 27(4) (1977) 348-351",
         "doi":  "10.1080/00022470.1977.10470432"},
        {"autor": "Hesketh, H.E.",
         "title": "Fine Particle Collection Efficiency Related to Pressure "
                  "Drop, Scrubbant and Particle Properties, and Contact "
                  "Mechanism",
         "ref": "J. Air Pollution Control Assoc., 24(10) (1974) 939-942",
         "doi":  "10.1080/00022470.1974.10469992"},
        {"autor": "Johnstone, H.F., Feild, R.B., Tassler, M.C.",
         "title": "Gas Absorption and Aerosol Collection in a Venturi "
                  "Atomizer",
         "ref": "Ind. Eng. Chemistry 46(8) (1954) 1601-1608",
         "doi":  "10.1021/ie50536a028"},
    ]

    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""

        self.statusDeltaP = 1
        if self.kwargs["modelo_DeltaP"] in (3, 4) and not self.kwargs["Lt"]:
            self.statusDeltaP = 0

        if not self.kwargs["entradaGas"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined gas stream")
            self.status = 0
            return
        if not self.kwargs["entradaLiquido"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined liquid stream")
            self.status = 0
            return
        if self.kwargs["tipo_calculo"] == 0 and not self.kwargs["diametro"]:
            self.msg = QApplication.translate("pychemqt", "undefined diameter")
            self.status = 0
            return
        elif self.kwargs["tipo_calculo"] == 1 and \
                not self.kwargs["rendimientoAdmisible"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined efficiency")
            self.status = 0
            return

        if self.kwargs["modelo_rendimiento"] == 0 and not self.kwargs["k"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined venturi constant")
            self.status = 3
        elif self.kwargs["modelo_rendimiento"] == 1 and not self.kwargs["f"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined calvert coefficient")
            self.status = 3

        return True

    def calculo(self):
        Gas = self.kwargs["entradaGas"]
        Liquido = self.kwargs["entradaLiquido"]
        sigma = Liquido.sigma
        rhoL = Liquido.Liquido.rho
        muL = Liquido.Liquido.mu

        self.Dt = Length(self.kwargs["diametro"])
        self.Lt = Length(self.kwargs["Lt"])

        if self.kwargs["k"]:
            self.k = Dimensionless(self.kwargs["k"])
        else:
            self.k = Dimensionless(1000.)
        if self.kwargs["f"]:
            self.f = Dimensionless(self.kwargs["f"])
        else:
            self.f = Dimensionless(0.5)

        self.At = Area(pi/4*self.Dt**2)
        self.Vg = Speed(Gas.Q/self.At)
        self.R = Liquido.Q/Gas.Q
        self.dd = Length(58600/self.Vg*(sigma/rhoL)**0.5+597*(
            muL/sigma**0.5/rhoL**0.5)**0.45*(1000*self.R)**1.5)

        self.rendimiento_parcial = self._Efficiency()
        self.rendimiento = self._GlobalEfficiency(self.rendimiento_parcial)

        if self.statusDeltaP:
            self.deltaP = self._deltaP()
        else:
            self.deltaP = DeltaP(0)

        self.CalcularSalidas(Gas)
        self.Pin = min(Gas.P, Liquido.P)

    def _Salidas(self, Gas):
        Liquido = self.kwargs["entradaLiquido"]
        unfiltered, filtered = Gas.solido.Separar(self.rendimiento_parcial)
        Pout = min(Gas.P, Liquido.P)-self.deltaP
        self.salida = []
        self.salida.append(Gas.clone(solido=unfiltered, P=Pout))
        self.salida.append(Liquido.clone(solido=filtered, P=Pout))

    def _Efficiency(self):
        Gas = self.kwargs["entradaGas"]
        Liquido = self.kwargs["entradaLiquido"]
        rhoS = Gas.solido.rho
        muG = Gas.Gas.mu
        rhoL = Liquido.Liquido.rho

        rendimiento_fraccional = []
        if self.kwargs["modelo_rendimiento"] == 0:
            # Modelo de Johnstone (1954)
            l = sqrt(pi/8)*Gas.Gas.mu/0.4987445/sqrt(Gas.Gas.rho*Gas.P)
            for dp in Gas.solido.diametros:
                Kn = l/dp*2
                C = Cunningham(l, Kn)
                kp = C*rhoS*dp**2*self.Vg/9/Gas.Gas.mu/self.dd
                penetration = exp(-self.k*self.R*kp**0.5)
                rendimiento_fraccional.append(1-penetration)

        elif self.kwargs["modelo_rendimiento"] == 1:
            # Modelo de Calvert (1972)
            l = sqrt(pi/8)*muG/0.4987445/sqrt(Gas.Gas.rho*Gas.P)
            for dp in Gas.solido.diametros:
                Kn = l/dp*2
                C = Cunningham(l, Kn)
                kp = C*rhoS*dp**2*self.Vg/9/muG/self.dd
                b = (-0.7-kp*self.f+1.4*log((kp*self.f+0.7)/0.7)+0.49 /
                     (0.7+kp*self.f))
                penetration = exp(self.R*self.Vg*rhoL*self.dd/55/muG*b/kp)
                if penetration > 1:
                    penetration = 1
                elif penetration < 0:
                    penetration = 0
                rendimiento_fraccional.append(1-penetration)

        return rendimiento_fraccional

    def _GlobalEfficiency(self, rendimientos):
        Gas = self.kwargs["entradaGas"]

        rendimiento_global = 0
        for i, fraccion in enumerate(Gas.solido.fracciones):
            rendimiento_global += rendimientos[i]*fraccion
        return Dimensionless(rendimiento_global)

#        DeltaP=Pressure(1.002*V**2*R, "kPa")
#
#        self.DeltaP=DeltaP
#
#        print V, DeltaP

#Calvert: 'dp = 5.4e-04 * (v^2) * rho_gas * (L/G)'
#Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * (0.56 + 0.125*L/G + 0.0023*(L/G)^2) / 507'
#Simplified Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * ((L/G)^0.78) / 1270'

    def _deltaP(self):
        Gas = self.kwargs["entradaGas"]
        Liquido = self.kwargs["entradaLiquido"]
        rhoL = Liquido.Liquido.rho

        if self.kwargs["modelo_DeltaP"] == 0:
            # Calvert (1968)
            deltaP = 0.85*rhoL*self.Vg**2*self.R
        elif self.kwargs["modelo_DeltaP"] == 1:
            # Hesketh (1974)
            deltaP = 1.36e-4*self.Vg.cms**2*rhoL.gcc*self.At.cm2**0.133 * \
                (0.56+935*self.R+1.29e-2*self.R**2)
        elif self.kwargs["modelo_DeltaP"] == 2:
            # Gleason (1971)
            deltaP = 2.08e-5*self.Vg.cms**2*(0.264*Liquido.Liquido.Q.ccs+73.8)
        elif self.kwargs["modelo_DeltaP"] == 3:
            # Volgin (1968)
            deltaP = 3.32e-6*self.Vg.cms**2*self.R*0.26*self.Lt**1.43
        elif self.kwargs["modelo_DeltaP"] == 4:
            # Yung (1977)
            Re = self.dd*self.Vg+Gas.Gas.rho/Gas.Gas.mu
            Cd = 0.22+(24/Re * (1+0.15*Re**0.6))
            X = 3*self.Lt*Cd*Gas.Gas.rho/16/self.dd/rhoL + 1
            deltaP = 2*rhoL*self.Vg**2*self.R*(1-X**2+(X**4-X**2)**0.5)

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

        return DeltaP(deltaP)

    def propTxt(self):
        Gas = self.kwargs["entradaGas"]

        txt = os.linesep + "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(range(13))
        txt += Separador_SolidGas.propTxt(self, 13, Gas)
        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Mode"),
              ("TEXT_TIPO", "tipo_calculo"), str),
             (QApplication.translate("pychemqt", "Model"),
              ("TEXT_MODEL", "modelo_rendimiento"), str),
             (QApplication.translate("pychemqt", "Pressure Loss Model"),
              ("TEXT_MODEL_DELTAP", "modelo_DeltaP"), str),
             (QApplication.translate("pychemqt", "Throat Diameter"), "Dt",
              Length),
             (QApplication.translate("pychemqt", "Throat Length"), "Lt",
              Length),
             (QApplication.translate(
                 "pychemqt", "Johnstone method scrubber constant"),
                 "k", Dimensionless),
             (QApplication.translate(
                 "pychemqt", "Calvert method scrubber constant"),
                 "f", Dimensionless),
             (QApplication.translate("pychemqt", "Drops Diameter"), "dd",
              Length),
             (QApplication.translate("pychemqt", "Throat Cross Area"), "At",
              Area),
             (QApplication.translate("pychemqt", "Throat Speed"), "Vg", Speed)]

        for prop in Separador_SolidGas.propertiesEquipment():
            l.append(prop)

        return l


class Dryer(equipment):
    """Clase que define un equipo de secado de sólidos

    Parámetros:
        entradaSolido: Instancia de clase Corriente que define la entrada de sólidos
        entradaAire: Instancia de clase Psychrometry que define el aire entrante
        mode: integer que indica el tipo de cálculo
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

            "mode": 0,
            "HR": 0.0,
            "TemperaturaSolid": 0.0,
            "HumedadResidual": 0.0,
            "Heat": 0.0,
            "deltaP": 0.0}

    kwargsInput = ("entradaAire", "entradaSolido")
    kwargsValue=("HR", "TemperaturaSolid", "HumedadResidual", "Heat", "deltaP")
    kwargsList=("mode", )
    calculateValue=("CombustibleRequerido", "Heat")

    TEXT_MODE=(
        QApplication.translate("pychemqt", "Rating, calculate output conditions"),
        QApplication.translate("pychemqt", "Design, calculate air flow necessary"))


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
        if "entrada" in kwargs:
            kwargs["entradaSolido"]=kwargs["entrada"][0]
            kwargs["entradaSolido"]=kwargs["entrada"][1]
            del kwargs["entrada"]
        self.kwargs.update(kwargs)


    def calculo(self):
        #TODO: De momento, no se implementan las cuestiones de cinetica de intercambio de calor y materia que definirían las dimensiones necesarias del equipo
        HR=self.kwargs.get("HR", 100)
        self.Heat=Power(self.kwargs["Heat"])
        self.deltaP=Pressure(self.kwargs["deltaP"])
        self.entradaAire=self.kwargs["entradaAire"]

        Pout=min(self.kwargs["entradaSolido"].P.atm, self.kwargs["entradaAire"].P.atm)-self.deltaP.atm

        aguaSolidoSalida=self.kwargs["HumedadResidual"]*self.kwargs["entradaSolido"].solido.caudal.kgh
        aguaSolidoEntrada=self.kwargs["entradaSolido"].caudalmasico.kgh
        if self.kwargs["mode"]==0:
            Caudal_aguaenAireSalida=aguaSolidoEntrada-aguaSolidoSalida+self.entradaAire.caudalMasico.kgh*self.entradaAire.Xw
            Caudal_airesalida=self.entradaAire.caudalMasico.kgh*self.entradaAire.Xa
            if self.entradaAire.Hs>Caudal_aguaenAireSalida/Caudal_airesalida:
                H=Caudal_aguaenAireSalida/Caudal_airesalida
            else:
                H=self.entradaAire.Hs
                aguaSolidoSalida+=Caudal_aguaenAireSalida/Caudal_airesalida-self.entradaAire.Hs
            self.SalidaAire=PsychroState(caudal=Caudal_aguaenAireSalida+Caudal_airesalida, tdb=self.entradaAire.Tdb, H=H)
            self.SalidaSolido=self.kwargs["entradaSolido"].clone(T=self.SalidaAire.Tdb, P=Pout, split=aguaSolidoSalida/aguaSolidoEntrada)
        else:
            pass


#        if self.HumedadResidual==0:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        else:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        self.SalidaAire=Corriente(self.entradaAire.T, self.entradaAire.P.atm-self.deltaP.atm, self.entradaAire.caudalmasico.kgh, self.entradaAire.mezcla)




if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

    from lib.solids import Solid
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(T=300, caudalSolido=[1/3600.], distribucion_diametro=diametros, distribucion_fraccion=fracciones, solids=[638])
    kw = {"fraccionMolar": [1.], "MEoS": True}
    aire=Corriente(T=350, P=101325, caudalMasico=0.01, ids=[475], solido=solido, **kw)
    agua=Corriente(T=300, P=101325, caudalMasico=0.1, ids=[62], **kw)
    secador=Scrubber(entradaGas=aire, entradaLiquido=agua, modelo_rendimiento=1, diametro=0.25, f=0.5)
    print(secador.propTxt())
