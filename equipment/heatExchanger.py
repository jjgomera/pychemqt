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
# library for heat exchanger equipment calculation
# - Heat_Exchanger
# - Heat_ExchangerDesign
# - Shell_Tube
# - Hairpin
# - Fired_Heater
###############################################################################


import os
from math import factorial

from PyQt5.QtWidgets import QApplication
from scipy import sqrt, exp, log, pi, arccos, sin, cos, tanh
from scipy.optimize import fsolve
from scipy.constants import g

from lib import unidades
from lib.corriente import Corriente
from lib.friction import f_friccion
from lib.adimensional import Re, Pr, Gr, Gz
from lib.heatTransfer import *
from .parents import equipment


# Equipment
class Heat_Exchanger(equipment):
    """Define a simple heat exchanger, only make energy balance

    Parameters:
        entrada: Corriente instance to define input stream
        Heat: global heat value exchanged
        Tout: Output temperature
        DeltaT: Increase temperature in stream
        A: area for heat exchange
        U: Global coefficient of heat transmision
        Text: Ambient external temperature
        DeltaP: Opcional pressure losses of equipment
    """
    title = QApplication.translate("pychemqt", "Heat Exchanger")
    help = ""
    kwargs = {
        "entrada": None,
        "Heat": 0.0,
        "Tout": 0.0,
        "DeltaT": 0.0,
        "A": 0.0,
        "U": 0.0,
        "Text": 0.0,
        "DeltaP": 0.0}

    kwargsInput = ("entrada", )
    kwargsValue = ("Heat", "Tout", "DeltaT", "A", "U", "Text", "DeltaP")
    calculateValue = ("HeatCalc", "ToutCalc")

    def cleanOldValues(self, **kwargs):
        """Clean kwargs values for old heat exchanger definition"""
        if kwargs.get("Tout", 0):
            self.kwargs["DeltaT"] = 0
            self.kwargs["Heat"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        elif kwargs.get("DeltaT", 0):
            self.kwargs["Tout"] = 0
            self.kwargs["Heat"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        elif kwargs.get("Heat", 0):
            self.kwargs["Tout"] = 0
            self.kwargs["DeltaT"] = 0
            self.kwargs["A"] = 0
            self.kwargs["U"] = 0
            self.kwargs["Text"] = 0
        self.kwargs.update(kwargs)

    @property
    def isCalculable(self):
        """
        modo: unknown variable to calculate
            1 - Known output temperature, calculate other variables
            2 - known heat exchange
            3 - known heat exchanger characteristic, calculate output stream
        """
        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            self.modo = 0
            return

        if self.kwargs["A"] and self.kwargs["U"] and self.kwargs["Text"]:
            self.modo = 3
        elif self.kwargs["Heat"]:
            self.modo = 2
        elif self.kwargs["Tout"] or self.kwargs["DeltaT"]:
            self.modo = 1
        else:
            self.msg = QApplication.translate(
                "pychemqt", "undefined output temperature specification")
            self.status = 0
            self.modo = 0

        if self.modo:
            self.msg = ""
            self.status = 1
            return True

    def calculo(self):
        entrada = self.kwargs["entrada"]
        self.DeltaP = unidades.DeltaP(self.kwargs["DeltaP"])
        self.HeatCalc = unidades.Power(self.kwargs["Heat"])
        if self.kwargs["Tout"]:
            Tout = unidades.Temperature(self.kwargs["Tout"])
        elif self.kwargs["DeltaT"]:
            Tout = unidades.Temperature(entrada.T+self.kwargs["DeltaT"])
        A = unidades.Area(self.kwargs["A"])
        U = unidades.HeatTransfCoef(self.kwargs["U"])
        Text = unidades.Temperature(self.kwargs["Text"])

        if self.modo == 1:
            self.salida = [entrada.clone(T=Tout, P=entrada.P-self.DeltaP)]
            self.HeatCalc = unidades.Power(self.salida[0].h-entrada.h)
        else:
            if self.modo == 2:
                self.HeatCalc = unidades.Power(0)
            else:
                self.HeatCalc = unidades.Power(A*U*(Text-entrada.T))

            def f():
                output = entrada.clone(T=T, P=entrada.P-self.DeltaP)
                return output.h-entrada.h-self.HeatCalc
            T = fsolve(f, entrada.T)[0]
            if T > max(Text, entrada.T) or T < min(Text, entrada.T):
                T = self.Text
            self.salida = [entrada.clone(T=T, P=entrada.P-self.DeltaP)]

        self.Tin = entrada.T
        self.ToutCalc = self.salida[0].T
        self.DeltaT = unidades.DeltaT(self.ToutCalc-entrada.T)

    def propTxt(self):
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Calculate properties")
        txt += "-----------------#" + os.linesep
        txt += "%-25s\t%s" % (QApplication.translate("pychemqt", "Input Temperature"), self.kwargs["entrada"].T.str)+os.linesep
        txt += "%-25s\t%s" % (QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt += "%-25s\t%s" % (QApplication.translate("pychemqt", "Temperature increase"), self.DeltaT.str)+os.linesep
        txt += "%-25s\t%s" % (QApplication.translate("pychemqt", "Pressure increase"), self.DeltaP.str)+os.linesep
        txt += "%-25s\t%s" % (QApplication.translate("pychemqt", "Heat"), self.HeatCalc.str)+os.linesep
        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Input Temperature"), "Tin", unidades.Temperature),
             (QApplication.translate("pychemqt", "Output Temperature"), "ToutCalc", unidades.Temperature),
             (QApplication.translate("pychemqt", "Temperature increase"), "DeltaT", unidades.DeltaT),
             (QApplication.translate("pychemqt", "Pressure increase"), "DeltaP", unidades.DeltaP),
             (QApplication.translate("pychemqt", "Heat"), "HeatCalc", unidades.Power)]
        return l


class Heat_ExchangerDesign(equipment):
    """Clase generica que define las caracteristicas comunes de los cambiadores de calor que requieren diseño"""

#Heat Exchanger design methods
    @classmethod
    def efectividad(cls, NTU, C_, flujo, **kwargs):
        """Calculo de la efectividad del cambiador
        Flujo vendra definido por su acronimo
        CF: Counter flow
        PF: Parallel flow
        CrFMix: Crossflow, both fluids mixed
        CrFSMix: Crossflow, one fluid mixed, other unmixed
        CrFunMix: Crossflow, both fluids unmixed
        1-2TEMAE: 1-2 pass shell and tube exchanger

        kwargs: Opciones adicionales:
            mixed: corriente mezclada para CrFSMix
                Cmin, Cmax
        """
        if C_ == 0:
            ep = 1-exp(-NTU)

        elif flujo == "PF":
            if C_ == 1:
                ep = (1-exp(-2*NTU))/2
            else:
                ep = (1-exp(-NTU*(1+C_)))/(1+C_)

        elif flujo == "CF":
            if C_ == 1:
                ep = NTU/(1+NTU)
            else:
                ep = (1-exp(-NTU*(1-C_)))/(1-C_*exp(-NTU*(1-C_)))

        elif flujo == "CrFunMix":
            def P(n, y):
                suma = 0
                for j in range(1, n+1):
                    suma += (n+1-j)/factorial(j)*y**(n+j)
                return suma/factorial(n+1)
            n = 1
            suma = 0
            while True:
                inc = C_**n*P(n, NTU)
                suma += inc
                n += 1
                if inc < 1e-12:
                    break
            ep = 1-exp(-NTU)-exp(-(1+C_)*NTU)*suma

        elif flujo == "CrFMix":
            if C_ == 1:
                ep = 1/(2/(1-exp(-NTU))-1/NTU)
            else:
                ep = 1/(1/(1-exp(-NTU))+C_/(1-exp(-NTU*C_))-1/NTU)

        elif flujo == "CrFSMix":
            if C_ == 1:
                ep = 1-exp(-(1-exp(-NTU)))
            else:
                if kwargs["mixed"] == "Cmin":
                    ep = 1-exp(-(1-exp(-NTU*C_))/C_)
                else:
                    ep = (1-exp(-C_*(1-exp(-NTU))))/C_

        elif flujo == "1-2TEMAE":
            if C_ == 1:
                ep = 2/(2+2**0.5/tanh(2**0.5*NTU/2))
            else:
                ep = 2/((1+C_)+(1+C_**2)**0.5/tanh(NTU*(1+C_**2)**0.5/2))

        return ep

    @classmethod
    def TemperatureEffectiveness(cls, NTU, R, flujo, **kwargs):
        """Calculo de la temperatura efectividad del cambiador
        Flujo vendra definido por su acronimo
        CF: Counter flow
        PF: Parallel flow
        CrFMix: Crossflow, both fluids mixed
        CrFSMix: Crossflow, one fluid mixed, other unmixed
        CrFunMix: Crossflow, both fluids unmixed
        1-2TEMAE: 1-2 TEMA E
        1-2TEMAE2: 1-2 TEMA E, shell fluid flow divided
        1-3TEMAE: 1-3 TEMA E
        1-4TEMAE: 1-4 TEMA E
        1-1TEMAG: 1-1 TEMA G
        1-2TEMAG: 1-2 TEMA G
        1-1TEMAH: 1-1 TEMA H
        1-2TEMAH: 1-2 TEMA H
        1-1TEMAJ: 1-1 TEMA J
        1-2TEMAJ: 1-2 TEMA J
        1-4TEMAJ: 1-4 TEMA J

        kwargs: Opciones adicionales:
            mixed: corriente mezclada para CrFSMix
                1, 2
        """

        if flujo == "PF":
            if R == 1:
                ep = NTU/(1+NTU)
            else:
                ep = (1-exp(-NTU*(1-R)))/(1-R*exp(-NTU*(1-R)))

        elif flujo == "CF":
            if R == 1:
                ep = (1-exp(-2*NTU))/2.
            else:
                ep = (1-exp(-NTU*(1+R)))/(1+R)

        elif flujo == "CrFunMix":
            ep = 1-exp(NTU**0.22/R*(exp(-R*NTU**0.78)-1))
            #            def P(n, y):
            #                suma=0
            #                for j in range(1, n+1):
            #                    suma+=(n+1-j)/factorial(j)*y**(n+j)
            #                return suma/factorial(n+1)
            #            n=1
            #            suma=0
            #            while True:
            #                inc=R**n*P(n, NTU)
            #                suma+=inc
            #                n+=1
            #                if inc<1e-12:
            #                    break
            #            ep=1-exp(-NTU)-exp(-(1+R)*NTU)*suma

        elif flujo == "CrFMix":
            K1 = 1-exp(-NTU)
            if R == 1:
                ep = 1/(2/K1-1/NTU)
            else:
                K2 = 1-exp(-R*NTU)
                ep = 1/(1/K1+R/K2-1/NTU)

        elif flujo == "CrFSMix":
            K = 1-exp(-NTU)
            if R == 1:
                ep = 1-exp(-K)
            else:
                if kwargs["mixed"] == "1":
                    ep = (1-exp(-R*K))/R
                else:
                    K = 1-exp(-R*NTU)
                    ep = 1-exp(-K/R)

        elif flujo == "1-2TEMAE":
            if R == 1:
                ep = 1/(1+1/tanh(NTU/2**0.5)/2**0.5)
            else:
                E = (1+R**2)**0.5
                ep = 2/(1+R+E/tanh(E*NTU/2))

        elif flujo == "1-2TEMAE2":
            E = exp(NTU)
            if R == 2:
                ep = 0.5*(1-(1+E**-2)/2/(1+NTU))
            else:
                B = exp(-NTU*R/2.)
                ep = 1/R*(1-(2-R)*(2.*E+R*B)/(2+R)/(2.*E-R/B))

        elif flujo == "1-3TEMAE":
            l1 = -3./2+(9./4+R*(R-1))**0.5
            l2 = -3./2-(9./4+R*(R-1))**0.5
            l3 = R
            d = l1-l2
            X1 = exp(l1*NTU/3.)/2/d
            X2 = exp(l2*NTU/3.)/2/d
            X3 = exp(l3*NTU/3.)/2/d
            if R == 1:
                A = -exp(-NTU)/18-exp(NTU/3)/2+(NTU+5)/9
            else:
                A = X1*(R+l1)*(R-l2)/2/l1-X3*d-X2*(R+l2)*(R-l1)/2/l2+1/(1-R)
            B = X1*(R-l2)-X2*(R-l1)+X3*d
            C = X2*(3*R+l1)-X1*(3*R+l2)+X3*d
            ep = 1/R*(1-C/(A*C+B**2))

        elif flujo == "1-4TEMAE":
            if R == 1:
                A = 1/tanh(5**0.5*NTU/4)
                B = tanh(NTU/4)
                ep = 4/(4+5**0.5*A+B)
            else:
                D = (4+R**2)**0.5
                A = 1/tanh(D*NTU/4)
                B = tanh(NTU*R/4)
                ep = 4/(2*(1+R)+D*A+R*B)

        elif flujo == "1-1TEMAG":
            if R == 1:
                B = NTU/(2+NTU)
            else:
                D = exp(-NTU*(1-R)/2)
                B = (1-D)/(1-R*D)
            A = 1/(1+R)*(1-exp(-NTU*(1+R)/2))
            ep = A+B-A*B*(1+R)+R*A*B**2

        elif flujo == "1-2TEMAG":
            if R == 2:
                alfa = exp(-NTU)
                ep = (1+2*NTU-alfa**2)/(4+4*NTU-(1-alfa)**2)
            else:
                alfa = exp(-NTU*(2+R)/4)
                beta = exp(-NTU*(2-R)/2)
                A = -2*R*(1-alfa)**2/(2+R)
                B = (4-beta*(2+R))/(2-R)
                ep = (B-alfa**2)/(A+2+R*B)

        elif flujo == "1-1TEMAH":
            A = 1/(1+R/2)*(1-exp(-NTU*(1+R/2)/2))
            if R == 2:
                B = NTU/(2+NTU)
            else:
                D = exp(-NTU*(1-R/2)/2)
                B = (1-D)/(1-R*D/2)
            E = (A+B-A*B*R/2)/2
            ep = E*(1+(1-B*R/2)*(1-A*R/2+A*B*R))-A*B*(1-B*R/2)

        elif flujo == "1-2TEMAH":
            if R == 4:
                H = NTU
                E = NTU/2
            else:
                beta = NTU*(4-R)/8
                H = (1-exp(-2*beta))/(4/R-1)
                E = (1-exp(-beta))/(4/R-1)
            alfa = NTU*(4-R)/8
            D = (1-exp(-alfa))/(4/R+1)
            G = (1-D)**2*(D**2+E**2)+D**2*(1+E)**2
            B = (1+H)*(1+E)**2
            ep = 1/R*(1-(1-D)**4/(B-4*G/R))

        elif flujo == "1-1TEMAJ":
            A = exp(NTU)
            if R == 2:
                ep = 0.5*(1-(1+1/A**2)/2/(1+NTU))
            else:
                B = exp(-NTU*R/2)
                ep = 1/R*(1-(2-R)*(2*A+R*B)/(2+R)/(2*A-R/B))

        elif flujo == "1-2TEMAJ":
            l = (1+R**2/4)**0.5
            A = exp(NTU)
            B = (A**l+1)/(A**l-1)
            C = A**((1+l)/2)/(l-1+(1+l)*A**l)
            D = 1+l*A**((l-1)/2)/(A**l-1)
            ep = 1/(1+R/2+l*B-2*l*C*D)

        elif flujo == "1-4TEMAJ":
            l = (1+R**2/16)**0.5
            A = exp(NTU)
            B = (A**l+1)/(A**l-1)
            C = A**((1+l)/2)/(l-1+(1+l)*A**l)
            D = 1+l*A**((l-1)/2)/(A**l-1)
            E = exp(R*NTU/2)
            ep = 1/(1+R/4*(1+3*E)/(1+E)+l*B-2*l*C*D)

        return ep

    @classmethod
    def CorrectionFactor(cls, P, R, flujo, **kwargs):
        """Calculo de la factor de correccion
        Flujo vendra definido por su acronimo
        CF: Counter flow
        PF: Parallel flow
        CrFMix: Crossflow, both fluids mixed
        CrFSMix: Crossflow, one fluid mixed, other unmixed
        CrFunMix: Crossflow, both fluids unmixed
        1-2TEMAE: 1-2 pass shell and tube exchanger

        kwargs: Opciones adicionales:
            mixed: corriente mezclada para CrFSMix
                Cmin, Cmax
        """
        if flujo == "PF" or flujo == "CF":
            f = 1

        elif flujo == "CrFSMix":
            if kwargs["mixed"] == "1":
                f = log((1-R*P)/(1-P))/(1-1/R)/log(1+R*log(1-P))
            else:
                f = log((1-R*P)/(1-P))/(R-1)/log(1+log(1-R*P)/R)

        elif flujo == "1-2TEMAE":
            if R == 1:
                if P*(2+2**0.5) >= 2:
                    f = 0
                else:
                    f = 2**0.5*P/(1-P)/log((2-P*(2-2**0.5))/(2-P*(2+2**0.5)))
            else:
                E = (1+R**2)**0.5
                if P*(1+R+E) >= 2:
                    f = 0
                else:
                    f = E*log((1-R*P)/(1-P))/(1-R)/log((2-P*(1+R-E))/(2-P*(1+R+E)))

        else:  # Para los ordenamientos de flujo sin solucion analitica
            NTU = cls.NTU_fPR(P, R, flujo, **kwargs)
            if R == 1:
                f = P/NTU/(1-P)
            else:
                f = log((1-R*P)/(1-P))/NTU/(1-R)

        return f

    @classmethod
    def NTU_fPR(cls, P, R, flujo, **kwargs):
        """Calculo de la factor de correccion
        Flujo vendra definido por su acronimo
        CF: Counter flow
        PF: Parallel flow
        CrFMix: Crossflow, both fluids mixed
        CrFSMix: Crossflow, one fluid mixed, other unmixed
        CrFunMix: Crossflow, both fluids unmixed
        1-2TEMAE: 1-2 pass shell and tube exchanger

        kwargs: Opciones adicionales:
            mixed: corriente mezclada para CrFSMix
                Cmin, Cmax
        """

        if flujo == "1-2TEMAE":
            if R == 1:
                NTU = log((1-P)/2-3*P)
            else:
                E = (1+R**2)**0.5
                NTU = log((2-P*(1+R-E))/(2-P*(1+R+E)))/E

        else:
            if R == 1:
                NTU = P/(1-P)
            else:
                NTU = log((1-R/P)/(1-P))/(1-R)

        return NTU

    @classmethod
    def Fi(cls, P, R, flujo, **kwargs):
        F = cls.CorrectionFactor(P, R, flujo, **kwargs)
        if R == 1:
            Fi = F*(1-P)
        else:
            Fi = F*P*(1-R)/log((1-R*P)/(1-P))
        return Fi


class Shell_Tube(Heat_ExchangerDesign):
    """Clase que define un cambiador de calor de carcasa y tubos

    Parámetros:
        entrada: Array con dos Instancia de clase corriente que define las corrientes que fluye por el equipo, en el orden [tubo, carcasa]
        entradaTubo: Instancia de clase corriente que define la corriente que pasa por los tubos
        entradaCarcasa: Instancia de calse corriente que define la corriente que pasa por la carcasa

    Standard:
        class_: Clase del standard TEMA:
            0   -   Clase R
            1   -   Clase B
            2   -   Clase C
        frontHead: Tipo de cabezal inicial
            0   -   Channel & Removable Cover
            1   -   Bonnet
            2   -   Removable Bundle
            3   -   Special High Pressure Closure
            4   -   Channel with Tubesheet & Removable Cover
        shell: Tipo de carcasa
            0   -   One Pass
            1   -   Two Pass
            2   -   Split Flow
            3   -   Double Split Flow
            4   -   Divided Flow
            5   -   Kettle Reboiler
            6   -   Cross Flow
        rearHead: Tipo de cabezal final
            0   -   Fixed Tubesheet (A head)
            1   -   Fixed Tubesheet (B head)
            2   -   Fixed Tubesheet (N head)
            3   -   Outside Packed Flt Head
            4   -   Flt Head with Backing Dev
            5   -   Pull Throught Flt Heat
            6   -   U-Tube Bundle
            7   -   Exit Sealed Flt Tubesheet
        orientation: Orientaction del cambiador
            0   -   Horizontal
            1   -   Vertical

    Métodos:
        tubesideLaminar: Método de cálculo de h en el lado del tubo en regimen laminar
            0   -   Eubank-Proctor
            1   -   VDI mean Nusselt
            2   -   Hausen
            3   -   Sieder-Tate
        tubesideTurbulent: Método de cálculo de h en el lado del tubo en regimen turbulento
            0   -   Sieder-Tate
            1   -   Colburn
            2   -   Dittus-Boelter
            3   -   ESDU
            4   -   Gnielinski
            5   -   VDI mean Nusselt
        shellsideSensible: Método de cálculo de h en el lado de la carcasa
            0   -   Stream analysis
            1   -   Bell-Delaware
            2   -   Kern

    Tubo:
        NTubes: Número de tubos
        NPases: Número de pasos de los tubos por la carcasa
        LTube: Lóngitud de tubos
        DeTube: Diametro externo
        wTube: Espesor de la tubería
        rTube: rugosidad interna de los tubos
        kTube: Conductividad térmica
        distribucionTube: Distribucion de tubos
            0   -   Triangular, 30º
            1   -   Diamante, 45º
            2   -   Rotated Triangular, 60º
            3   -   Square, 90º
        pitch: Espacio entre tuberias
        finned: boolean que indica que la tubería tiene alerones
            0   -   tubería lisa
            1   -   tubería con alerones
        Nfin: numero de aletas por metro de tubería
        heightFin:
        foulingTube: resistencia por depositos en la parte del tubo

    Carcasa:
        parallel: Número de intercambiadores en paralelo
        serie: Número de intercambiadores en serie
        DShell: Diematro de la carcasa
        foulingShell: resistencia por depositos en la parte de la carcasa

    Baffle:
        typeBaffle: Tipo de baffle
            0   -   Single segmental
            1   -   Double segmental
            2   -   Triple segmental
            3   -   No tubes in window
            4   -   Disk & donut
            5   -   Rod
        baffleSpacingIn: Espacio de separación anterior al primer bafle
        baffleSpacing: Espacio de separación entre baffles
        baffleSpacingOut: Espacio de separación posterior al último bafle
        baffleThickness: Espesor de los baffles
        BaffleOrientation: Orientación de los baffles
            0   -   Horizontal
            1   -   Vertical
        baffleCut: Porcentaje de corte de los baffles
        baffleCutBase: Base de cálculo del porcentaje de corte
            0   -   Diámetro
            1   -   Área

    Clearances:
        clearanceTubeBaffle
        clearanceShellBaffle
        clearanceShellBundle
        sealingStrips
    Coste:
        tipo: tipo de cambiador
            0   -   Fired head
            1   -   Kettle reboiler
            2   -   U-tubes
        material:
            0   -   Carbon Steel
            1   -   Stainless steel 316
            2   -   Stainless steel 304
            3   -   Stainless steel 347
            4   -   Nickel 200
            5   -   Monel 400
            6   -   Inconel 600
            7   -   Incoloy 825
            8   -   Titanium
            9   -   Hastelloy
        P_dis: Presión de diseño, si no se especifica se usará la máxima presión del las corrientes del proceso
    """

    title = QApplication.translate("pychemqt", "Shell and Tube Heat Exchanger")
    help = ""
    kwargs = {
        "entrada": [],
        "entradaTubo": None,
        "entradaCarcasa": None,

        "class_": 0,
        "frontHead": 0,
        "shell": 0,
        "rearHead": 0,
        "orientation": 0,

        "tubesideLaminar": 0,
        "tubesideTurbulent": 0,
        "shellsideSensible": 0,

        "NTube": 0,
        "NPases": 0,
        "LTube": 0.0,
        "DeTube": 0.0,
        "wTube": 0.0,
        "rTube": 0.0,
        "kTube": 0.0,
        "distribucionTube": 0,
        "pitch": 0,
        "finned": 0,
        "Nfin": 0,
        "heightFin": 0.0,
        "foulingTube": 0.0,

        "parallel": 0,
        "serie": 0,
        "DShell": 0.0,
        "foulingShell": 0.0,

        "baffleType": 0,
        "baffleSpacingIn": 0.0,
        "baffleSpacing": 0.0,
        "baffleSpacingOut": 0.0,
        "baffleThickness": 0.0,
        "BaffleOrientation": 0,
        "baffleCut": 0.0,
        "baffleCutBase": 0,

        "nozzleInTubesideDiameter": 0.0,
        "nozzleOutTubesideDiameter": 0.0,
        "nozzleInShellsideDiameter": 0.0,
        "nozzleOutShellsideDiameter": 0.0,

        "clearanceTubeBaffle": 0.0,
        "clearanceShellBaffle": 0.0,
        "clearanceShellBundle": 0.0,
        "sealingStrips": 0.0,

        "modo": 0,

        "f_install": 3.,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "tipoCoste": 0,
        "materialCoste": 0,
        "P_dis": 0.0}

    indiceCostos = 2

    TEXT_METHOD_TUBE_LAMINAR = ["Eubank-Proctor", "VDI mean Nusselt",
                                "Hausen", "Sieder-Tate"]
    TEXT_METHOD_TUBE_TURBULENT = ["Sieder-Tate", "Colburn", "Dittus-Boelter",
                                  "ESDU", "Gnielinski", "VDI mean Nusselt"]
    TEXT_METHOD_SHELL = ["Stream analysis", "Bell-Delaware", "Kern"]
    TEXT_CLASS = ["TEMA R", "TEMA B", "TEMA C"]
    TEXT_FRONTHEAD = [
        "A - "+QApplication.translate("pychemqt", "Channel & Removable Cover"),
        "B - "+QApplication.translate("pychemqt", "Bonnet"),
        "C - "+QApplication.translate("pychemqt", "Removable Bundle"),
        "D - "+QApplication.translate("pychemqt", "Special High Pressure Closure"),
        "N - "+QApplication.translate("pychemqt", "Channel with Tubesheet & Removable Cover")]
    TEXT_SHELL = [
        "E - "+QApplication.translate("pychemqt", "One Pass"),
        "F - "+QApplication.translate("pychemqt", "Two Pass"),
        "G - "+QApplication.translate("pychemqt", "Split Flow"),
        "H - "+QApplication.translate("pychemqt", "Double Split Flow"),
        "J - "+QApplication.translate("pychemqt", "Divided Flow"),
        "K - "+QApplication.translate("pychemqt", "Kettle Reboiler"),
        "X - "+QApplication.translate("pychemqt", "Cross Flow")]
    TEXT_REARHEAD = [
        "L - "+QApplication.translate("pychemqt", "Fixed Tubesheet (A head)"),
        "M - "+QApplication.translate("pychemqt", "Fixed Tubesheet (B head)"),
        "N - "+QApplication.translate("pychemqt", "Fixed Tubesheet (N head)"),
        "P - "+QApplication.translate("pychemqt", "Outside Packed Flt Head"),
        "S - "+QApplication.translate("pychemqt", "Flt Head with Backing Dev"),
        "T - "+QApplication.translate("pychemqt", "Pull Throught Flt Heat"),
        "U - "+QApplication.translate("pychemqt", "U-Tube Bundle"),
        "W - "+QApplication.translate("pychemqt", "Exit Sealed Flt Tubesheet")]
    TEXT_ORIENTATION = [
        QApplication.translate("pychemqt", "Horizontal"),
        QApplication.translate("pychemqt", "Vertical")]
    TEXT_DISTRIBUTION_TUBE = [
        QApplication.translate("pychemqt", "Triangular")+", 30º",
        QApplication.translate("pychemqt", "Diamond")+", 45º",
        QApplication.translate("pychemqt", "Rotated Triangular")+", 60º",
        QApplication.translate("pychemqt", "Square")+", 90º"]
    TEXT_BAFFLE_TYPE = [
        QApplication.translate("pychemqt", "Single segmental"),
        QApplication.translate("pychemqt", "Double segmental"),
        QApplication.translate("pychemqt", "Triple segmental"),
        QApplication.translate("pychemqt", "No tubes in window"),
        QApplication.translate("pychemqt", "Disk & donut"),
        QApplication.translate("pychemqt", "Rod")]
    TEXT_COST_TYPE = [
        QApplication.translate("pychemqt", "Fixed Head"),
        QApplication.translate("pychemqt", "Kettle Reboiler"),
        QApplication.translate("pychemqt", "U-Tube")]
    TEXT_COST_MATERIAL = [
        QApplication.translate("pychemqt", "Carbon Steel"),
        QApplication.translate("pychemqt", "Stainless Steel 316"),
        QApplication.translate("pychemqt", "Stainless Steel 304"),
        QApplication.translate("pychemqt", "Stainless Steel 347"),
        QApplication.translate("pychemqt", "Nickel 200"),
        QApplication.translate("pychemqt", "Monel 400"),
        QApplication.translate("pychemqt", "Inconel 600"),
        QApplication.translate("pychemqt", "Incoloy 825"),
        QApplication.translate("pychemqt", "Titanium"),
        QApplication.translate("pychemqt", "Hastelloy")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entradaTubo"]:
            self.msg = QApplication.translate("pychemqt", "undefined tubeside input")
            self.status = 0
            return
        if not self.kwargs["entradaCarcasa"]:
            self.msg = QApplication.translate("pychemqt", "undefined shellside input")
            self.status = 0
            return

        return True

    def calculo(self):
        if self.kwargs["modo"]:
            #Diseño
            pass

        else:  # Evaluación
            N = self.kwargs["NTube"]
            De = unidades.Length(self.kwargs["DeTube"])
            Di = unidades.Length(De-2*self.kwargs["wTube"])
            L = unidades.Length(self.kwargs["LTube"])

            #TubeSide
            rho = self.kwargs["entradaTubo"].Liquido.rho
            mu = self.kwargs["entradaTubo"].Liquido.mu
            cp = self.kwargs["entradaTubo"].Liquido.cp
            k = self.kwargs["entradaTubo"].Liquido.k
            w = self.kwargs["entradaTubo"].Liquido.caudalmasico
            beta = self.kwargs["entradaTubo"].Liquido.alfav
            v = self.kwargs["entradaTubo"].Q/N*4/pi/Di
            re = Re(D=Di, V=v, rho=rho, mu=mu)
            pr = self.kwargs["entradaTubo"].Liquido.Pr
            gz = Gz(w=w, cp=cp, k=k, L=L)
            gr = Gr(beta=beta, T=self.kwargs["entradaTubo"].T, To=self.kwargs["entradaTubo"].T, L=L, mu=mu)

            if re < 2300:
                if self.kwargs["tubesideLaminar"] == 0:
                    Nu = h_tubeside_laminar_Eubank_Proctor(Pr=pr, Gz=gz, Gr=gr, D=Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 1:
                    Nu = h_tubeside_laminar_VDI(Re=re, Pr=pr, D=Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 2:
                    Nu = h_tubeside_laminar_Hausen(Gz=gz)
                elif self.kwargs["tubesideLaminar"] == 3:
                    Nu = h_tubeside_laminar_Sieder_Tate(Gz=gz, Gr=gr)
            else:
                if self.kwargs["tubesideTurbulent"] == 0:
                    Nu = h_tubeside_turbulent_Sieder_Tate(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 1:
                    Nu = h_tubeside_turbulent_Colburn(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 2:
                    frio = self.kwargs["entradaCarcasa"] > self.kwargs["entradaTubo"]
                    Nu = h_tubeside_turbulent_Dittus_Boelter(Re=re, Pr=pr, calentamiento=frio)
                elif self.kwargs["tubesideTurbulent"] == 3:
                    Nu = h_tubeside_turbulent_ESDU(Re=re, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 4:
                    Nu = h_tubeside_turbulent_Gnielinski(Re=re, Pr=pr, D=Di, L=L)
                elif self.kwargs["tubesideTurbulent"] == 5:
                    line = self.kwargs["distribucionTube"] == 3
                    filas = self.kwargs["NTube"]**0.5
                    Nu = h_tubeside_turbulent_VDI(Re=re, Pr=pr, filas_tubos=filas, alineados=line)

            hi = unidades.HeatTransfCoef(Nu*k/Di)

            # ShellSide
            if self.kwargs["shellsideSensible"] == 0:
                h, DP = self.h_shellside_turbulent_Stream_Analysis()
            elif self.kwargs["shellsideSensible"] == 1:
                self.h_shellside_turbulent_Bell_Delaware()
            else:
                h = self.h_shelside_turbulent_Kern()

            # Fouling
            fi = self.kwargs["foulingShell"]
            fo = self.kwargs["foulingTube"]

        #F: Heat exchanger Design, Operation, Maintenance and Enhancement - Ali A. Rabah
#        U=1/((1/ho+fo)/Ef+fw+fi(Ao/Ai)+Ao/Ai/hi)
#        deltaT=(deltaT2-deltaT1)/log(deltaT2/deltaT1)
#        U=1/(Do/hi/Di+Do*log(Do/Di)/2/k+1/ho)
#        #TODO: añadir resistencias de depositos en la pared
#        """Serth - Process heat transfer_ principles and applications pag 102"""
#        q=U*Ao*deltaT
        self.area = unidades.Area(25)


    def fw(self):
        if self.kwargs["finned"]:
            fw=self.kwargs["wTube"]/self.kwargs["kTube"]*((self.kwargs["DeTube"]+2*self.kwargs["Nfin"]*self.kwargs["heightFin"]*(self.kwargs["DeTube"]+self.kwargs["heightFin"]))/(self.kwargs["DeTube"]-self.kwargs["wTube"]))
        else:
            fw=self.kwargs["DeTube"]/2/self.kwargs["kTube"]*log(self.kwargs["DeTube"]/(self.kwargs["DeTube"]-2*self.kwargs["wTube"]))
        return fw


    @staticmethod
    def h_shellside_turbulent_Stream_Analysis():
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        Serth - Process Heat Transfer - Principles and applications Cap. 7
        """


    def h_shellside_turbulent_Bell_Delaware(self):
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        Serth - Process Heat Transfer - Principles and applications Cap. 6
        """
        fi = 1

        P = self.kwargs["pitch"]/self.kwargs["DeTube"]
        mo = self.kwargs["entradaCarcasa"].caudalmasico
        if self.kwargs["distribucionTube"] == 2:
            Ptef = self.kwargs["pitch"]/2**0.5
        else:
            Ptef = self.kwargs["pitch"]

        if self.kwargs["distribucionTube"] < 3:
            tita_tp = unidades.Angle([30, 45, 60][self.kwargs["distribucionTube"]], "deg")
            Pt_ = self.kwargs["pitch"]*cos(tita_tp)
        else:
            Pt_ = self.kwargs["pitch"]
        Nc = self.kwargs["DShell"]*(1-2*self.kwargs["baffleCut"])/Pt_
        Ncw = 0.8*self.kwargs["baffleCut"]*self.kwargs["DShell"]/Pt_

        Dotl = self.kwargs["DShell"]-2*self.kwargs["clearanceShellBundle"]
        Sm = self.kwargs["baffleSpacing"]*(self.kwargs["DShell"]-Dotl+(Dotl-self.kwargs["DeTube"])/Ptef*(self.kwargs["pitch"]-self.kwargs["DeTube"]))

        G = mo/Sm
        Re = self.kwargs["DeTube"]*G/self.kwargs["entradaCarcasa"].Liquido.mu
        Pr = self.kwargs["entradaCarcasa"].Liquido.Pr
        rho = self.kwargs["entradaCarcasa"].Liquido.rho
        cp = self.kwargs["entradaCarcasa"].Liquido.cp

        Dctl = Dotl-self.kwargs["DeTube"]
        tita_ctl = 2*arccos(self.kwargs["DShell"]*(1-2*self.kwargs["baffleCut"])/Dctl)
        Fc = 1+1./pi*(sin(tita_ctl)-tita_ctl)
        Fw = 1./2/pi*(tita_ctl-sin(tita_ctl))
        Stb = pi/8*((self.kwargs["DeTube"]+2*self.kwargs["clearanceTubeBaffle"])**2-self.kwargs["DeTube"]**2)*self.kwargs["NTube"]*(1+Fc)

        tita_ds = 2*arccos(1-2*self.kwargs["baffleCut"])
        Ssb = self.kwargs["DShell"]*self.kwargs["clearanceShellBaffle"]*(pi-0.5*tita_ds)

        Sb = self.kwargs["baffleSpacing"]*(self.kwargs["DShell"]-Dotl)
        Sw = 1./8*self.kwargs["DShell"]**2*(tita_ds-sin(tita_ds))-1./4*self.kwargs["NTube"]*Fw*pi*self.kwargs["DeTube"]**2
        Dw = 4*Sw/(pi*self.kwargs["DeTube"]*self.kwargs["NTube"]*0.5*(1-Fc)+self.kwargs["DShell"]*tita_ds)

        Jc = 0.55+0.72*Fc

        rs = Ssb/(Ssb+Stb)
        rl = (Ssb+Stb)/Sm
        Jl = 0.44*(1-rs)+(1-0.44*(1-rs))*exp(-2.2*rl)
        Rl = exp(-1.33*(1+rs)*rl**(0.8-0.15*(1+rs)))

        rss = self.kwargs["sealingStrips"]/Nc
        if rss < 0.5:
            if Re < 100:
                Cj = 1.35
                Cr = 4.5
            else:
                Cj = 1.25
                Cr = 3.7
            Jb = exp(-Cj*Sb/Sm*(1-(2*rss)**(1./3)))
            Rb = exp(-Cr*Sb/Sm*(1-(2*rss)**(1./3)))
        else:
            Jb = 1.
            Rb = 1.

        if Re < 100:
            n1 = 1./3
            n2 = 1.
        else:
            n1 = 0.6
            n2 = 0.2
        nb = self.kwargs["LTube"]/self.kwargs["baffleSpacing"]
        if self.kwargs["baffleSpacingIn"] > self.kwargs["baffleSpacing"]:
            nb += 1
        if self.kwargs["baffleSpacingOut"] > self.kwargs["baffleSpacing"]:
            nb += 1
        Js = (nb-1+(self.kwargs["baffleSpacingIn"]/self.kwargs["baffleSpacing"])**(1-n1)+(self.kwargs["baffleSpacingOut"]/self.kwargs["baffleSpacing"])**(1-n1))/(nb-1+(self.kwargs["baffleSpacingIn"]/self.kwargs["baffleSpacing"])+self.kwargs["baffleSpacingOut"]/self.kwargs["baffleSpacing"])
        Rs = 0.5*((self.kwargs["baffleSpacing"]/self.kwargs["baffleSpacingIn"])**(2-n2)+(self.kwargs["baffleSpacing"]/self.kwargs["baffleSpacingOut"])**(2-n2))

        Nct = (nb+1)*(Nc+Ncw)
        if Re <= 20:
            Jr = (10/Nct)**0.18
        elif Re >= 100:
            Jr = 1.
        else:  # Interpolacion entre los valores de arriba
            Jr = 0.853379+0.0014662*Re

        if self.kwargs["distribucionTube"] == 1:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 1.55
                a2 = -0.667
            elif Re < 100:
                a1 = 0.498
                a2 = -0.656
            elif Re < 1000:
                a1 = 0.73
                a2 = -0.500
            elif Re < 10000:
                a1 = 0.37
                a2 = -0.396
            else:
                a1 = 0.37
                a2 = -0.396
                a3 = 1.93
                a4 = 0.5
        elif self.kwargs["distribucionTube"] == 3:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 0.97
                a2 = -0.667
            elif Re < 100:
                a1 = 0.9
                a2 = -0.631
            elif Re < 1000:
                a1 = 0.408
                a2 = -0.46
            elif Re < 10000:
                a1 = 0.107
                a2 = -0.266
            else:
                a1 = 0.37
                a2 = -0.395
                a3 = 1.187
                a4 = 0.37
        else:
            a3, a4 = 0, 0
            if Re < 10:
                a1 = 1.4
                a2 = -0.667
            elif Re < 100:
                a1 = 1.36
                a2 = -0.657
            elif Re < 1000:
                a1 = 0.593
                a2 = -0.477
            elif Re < 10000:
                a1 = 0.321
                a2 = -0.388
            else:
                a1 = 0.321
                a2 = -0.388
                a3 = 1.45
                a4 = 0.519

        if self.kwargs["distribucionTube"] == 1:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 32
                b2 = -1.
            elif Re < 100:
                b1 = 26.2
                b2 = -0.913
            elif Re < 1000:
                b1 = 3.5
                b2 = -0.476
            elif Re < 10000:
                b1 = 0.333
                b2 = -0.136
            else:
                b1 = 0.303
                b2 = -0.126
                b3 = 6.59
                b4 = 0.52
        elif self.kwargs["distribucionTube"] == 3:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 35.0
                b2 = -1.
            elif Re < 100:
                b1 = 32.1
                b2 = -0.963
            elif Re < 1000:
                b1 = 6.09
                b2 = -0.602
            elif Re < 10000:
                b1 = 0.0815
                b2 = 0.022
            else:
                b1 = 0.391
                b2 = -0.148
                b3 = 6.3
                b4 = 0.378
        else:
            b3, b4 = 0, 0
            if Re < 10:
                b1 = 48.0
                b2 = -1.
            elif Re < 100:
                b1 = 45.1
                b2 = -0.973
            elif Re < 1000:
                b1 = 4.570
                b2 = -0.476
            elif Re < 10000:
                b1 = 0.486
                b2 = -0.152
            else:
                b1 = 0.372
                b2 = -0.123
                b3 = 7.0
                b4 = 0.5

        a = a3/(1+0.14*Re**a4)
        b = b3/(1+0.14*Re**b4)

        j = a1*(1.33/P)**a*Re**a2
        f = b1*(1.33/P)**b*Re**b2
        hid = j*cp*G*fi/Pr**(2./3)
        h = unidades.HeatTransfCoef(hid*Jc*Jl*Jb*Jr*Js)

        DPideal = 2*f*Nc*G**2/g/rho/fi
        DPc = (nb-1)*DPideal*Rl*Rb
        if Re > 100:
            DPwideal = (2+0.6*Ncw)*mo**2/2/g/rho/Sm/Sw
        else:
            DPwideal = 26*nu*mo/g/(Sm*Sw)**0.5*(Ncw/P+self.kwargs["baffleCut"]*self.kwargs["DShell"]/Dw**2)+mo**2/g/rho/Sm/Sw
        DPw = nb*DPwideal*Rl
        DPe = 2*DPideal*(1+Ncw/Nc)*Rb*Rs

        mon = mo/self.kwargs["parallel"]

        DPn = 0
        for D in self.kwargs["nozzleInShellsideDiameter"], self.kwargs["nozzleOutShellsideDiameter"]:
            Ren = D*mon/self.kwargs["entradaCarcasa"].Liquido.mu
            Sn = pi/4*D**2
            if Ren > 100:
                DPn += 2e-13*self.kwargs["serie"]*mon**2/Sn
            else:
                DPn += 4e-13*self.kwargs["serie"]*mon**2/Sn
        DP = unidades.DeltaP(DPc+DPw+DPe+DPn)
        return h, DP

    @staticmethod
    def h_shelside_turbulent_Kern(Re, Pr):
        """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
        10<Re<1e6
        kern pag 137
        """
        return 0.36*Re**0.55*Pr**(1./3)*(mu/muw)**0.14


    def h_tubeside_laminar_condensation_Kern(self):
        return 0.815*(k**3*rho_l*(rho_l-rho_g)*g*l/(pi*mu_l*Do*(T-Tw)))**0.25


    def h_tubeside_laminar_condensation_Nusselt(self):
        return 0.72*eg**0.75*(k**3*rho_l*(rho_l-rho_g)*g*hlg/(mu_l*Do*(T-Tw)))**0.25

    def h_tubeSide_fined_Young(self):
        """Briggs, Katz, and Young, Chem.Eng. Prog., 59(11), 49–59 (1963)"""
        return 0.1378*Re**0.718*Pr**(1./3)*(finSpacing/finHeight)**0.296

    def coste(self):
        if self.kwargs["P_dis"]:
            Pd = unidades.Pressure(self.kwargs["P_dis"])
        else:
            Pd = unidades.Pressure(max(self.kwargs["entradaTubo"].P, self.kwargs["entradaCarcasa"].P))

        if self.kwargs["tipoCoste"] == 0:  # Fired head
            Fd = exp(-1.1156+0.09060*log(self.area.ft2))
        elif self.kwargs["tipoCoste"] == 1:  # Kettle reboiler
            Fd = 1.35
        else:  # U-tubes
            Fd = exp(-0.9816+0.0803*log(self.area.ft2))

        g1 = [0., 0.8603, 0.8193, 0.6116, 1.5092, 1.2989, 1.204, 1.1854, 1.5420, 0.1549][self.kwargs["materialCoste"]]
        g2 = [1., 0.23296, 0.15984, 0.22186, 0.60859, 0.43377, 0.50764, 0.49706, 0.42913, 0.51774][self.kwargs["materialCoste"]]
        Fm = g1+g2*log(self.area.ft2)

        if Pd.psi <= 300:
            Fp = 0.771+0.04981*log(self.area.ft2)
        elif Pd.psi <= 600:
            Fp = 1.0305+0.0714*log(self.area.ft2)
        else:
            Fp = 1.14+0.12088*log(self.area.ft2)

        C_base = exp(8.821-0.30863*log(self.area.ft2)+0.0681*log(self.area.ft2)**2)
        C = Fd*Fm*Fp*C_base

        self.C_adq = unidades.Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst = unidades.Currency(self.C_adq*self.kwargs["f_install"])


class Air_Cooler(equipment):
    """Clase que modela los enfriadores de aire"""
    title = "Enfriador de aire"
    help = ""

    def coste(self, *args):
        self._indicesCoste(*args)

        C = 24.6*self.area.ft2**0.4*1000

        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)


class Evaporator(equipment):
    """Clase que modela los evaporadores"""
    title = "Evaporador"
    help = ""

    def coste(self, *args, **kwargs):
        """
        tipo:
            0   -   Forced circulation
            1   -   Long tube
            2   -   Falling film
        material:
        forced circulation materials shell/tube
            0   -  Steel/copper
            1   -  Monel/cupronickel
            2   -  Nickel/nickel
        Long tube evaporator materials shell/tube
            0   -  Steel/copper
            1   -  Steel/steel
            2   -  Steel/aluminum
            3   -  Nickel/nickel
        """
        self._indicesCoste(*args)

        self.tipo = kwargs.get("tipo", 0)
        self.material = kwargs.get("material", 0)

        A = self.area.ft2

        if self.tipo == 0:
            C_base = exp(5.9785-0.6056*log(A)+0.08514*log(A)**2)*1000
            if self.material == 0:
                Fm = 1.
            elif self.material == 1:
                Fm = 1.35
            else:
                Fm = 1.8
        elif self.tipo == 1:
            C_base = 0.36*A**0.85*1000
            if self.material == 0:
                Fm = 1.
            elif self.material == 1:
                Fm = 0.6
            elif self.material == 2:
                Fm = 0.7
            else:
                Fm = 3.3
        else:
            C_base = exp(3.2362-0.0126*log(A)+0.0244*log(A)**2)*1000
            Fm = 1

        C = Fm*C_base
        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)


class Hairpin(Heat_ExchangerDesign):
    """Clase que modela los intercambiadores de calor de doble tubo

    Parámetros:
        entradaInterior: Instancia de clase corriente que define la corriente que pasa por la seccion interior
        entradaExterior: Instancia de calse corriente que define la corriente que pasa por la seccion exterior

        modo: Modo de cálculo
            0   -   Diseño
            1   -   Evaluacion
        flujo: Tipo de flujo
            0   -   Contracorriente
            1   -   Cocorriente
        orientacion: Orientacion de la tuberia
            0   -   Horizontal
            1   -   Vertical, interior descendente
            2   -   Vertical, interior ascendente
        metodo:
            0   -   Temperaturas medias
            1   -   Division del equipo en segmentos
        tubesideLaminar: Método de cálculo de h en el lado del tubo en regimen laminar
            0   -   Eubank-Proctor
            1   -   VDI mean Nusselt
            2   -   Hausen
            3   -   Sieder-Tate
        tubesideTurbulent: Método de cálculo de h en el lado del tubo en regimen turbulento
            0   -   Sieder-Tate
            1   -   Colburn
            2   -   Dittus-Boelter
            3   -   ESDU
            4   -   Gnielinski
            5   -   VDI mean Nusselt

        LTube: Longitud del tubo
        DeeTube: Diametro interno tuberia externa
        DeTube: Diametro externo pared tuberia interna
        DiTube: Diametro interno pared tuberia interna
        wTube: Espesor tuberia
        rTube: rugosidad interna de los tubos
        kTube: Conductividad térmica

        tubeFouling: Fouling en el lado del tubo
        annulliFouling: Fouling en el lado del anillo

        tubeFinned: Boolean que indica si el tubo tiene aletas
        hFin: Altura de la aleta
        thicknessBaseFin: Espesor en la base de la aleta
        thicknessTopFin: Espesor en lo alto de la aleta
        rootDoFin: Diametro externo en la base de las aletas
        kFin: Conductividad termica del material de la aleta
        nFin: Numero de aletas por metro de tuberia

        tubeTout: Temperatura de salida del fluido del tubo
        annulliTout: Temperatura de salida del fluido del anillo

    Coste:
        material:
            0   -  Carbon steel/carbon steel
            1   -  Carbon steel/304 stainless
            2   -  Carbon steel/316 stainless
        P_dis: Presion de diseño
    """
    title = QApplication.translate("pychemqt", "Hairpin Heat Exchanger")
    help = ""
    kwargs = {
        "entradaTubo": None,
        "entradaExterior": None,

        "modo": 0,
        "flujo": 0,
        "orientacion": 0,
        "tubesideLaminar": 0,
        "tubesideTurbulent": 0,
        "metodo": 0,
        "phase": 0,

        "DeeTube": 0.0,
        "DeTube": 0.0,
        "DiTube": 0.0,
        "wTube": 0.0,
        "rTube": 0.0,
        "kTube": 0.0,
        "LTube": 0.0,
        "nTube": 0.0,

        "tubeFouling": 0.0,
        "annulliFouling": 0.0,

        "tubeFinned": 0,
        "hFin": 0.0,
        "thicknessBaseFin": 0.0,
        "thicknessTopFin": 0.0,
        "rootDoFin": 0.0,
        "kFin": 0.0,
        "nFin": 0,

        "tubeTout": 0.0,
        "tubeXout": -1.0,
        "annulliTout": 0.0,
        "annulliXout": -1.0,

        "f_install": 3.,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "material": 0,
        "P_dis": 0}

    kwargsInput = ("entradaTubo", "entradaExterior")
    kwargsValue = ("DeTube", "DiTube", "wTube", "rTube", "kTube", "LTube",
                   "nTube", "tubeFouling", "annulliFouling", "P_dis",
                   "tubeTout", "annulliTout")
    kwargsList = ("modo", "flujo", "orientacion")
    kwargsCheck = ("tubeFinned", )
    calculateValue = ("Q", "ToutAnnulli", "ToutTube", "U", "A", "L",
                      "deltaPTube", "deltaPAnnulli", "CF")
    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 2

    TEXT_MODO = [
        QApplication.translate("pychemqt", "Design"),
        QApplication.translate("pychemqt", "Rating")]
    TEXT_FLUJO = [
        QApplication.translate("pychemqt", "Counterflow"),
        QApplication.translate("pychemqt", "Parallelflow")]
    TEXT_ORIENTACION = [
        QApplication.translate("pychemqt", "Horizontal"),
        QApplication.translate("pychemqt", "Vertical, (in down)"),
        QApplication.translate("pychemqt", "Vertical, (in up)")]
    TEXT_MATERIAL = [
        QApplication.translate("pychemqt", "Carbon steel/carbon steel"),
        QApplication.translate("pychemqt", "Carbon steel/304 stainless"),
        QApplication.translate("pychemqt", "Carbon steel/316 stainless")]
    CODE_FLUJO = ("CF", "PF")

    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entradaTubo"]:
            self.msg = QApplication.translate("pychemqt", "undefined internal stream input")
            self.status = 0
            return
        if not self.kwargs["entradaExterior"]:
            self.msg = QApplication.translate("pychemqt", "undefined external stream input")
            self.status = 0
            return

        if not self.kwargs["DeeTube"]:
            self.msg = QApplication.translate("pychemqt", "undefined pipe external diameter")
            self.status = 0
            return

        self.statusPipe = 0
        if self.kwargs["DeTube"] and self.kwargs["DiTube"]:
            self.statusPipe = 1
        elif self.kwargs["DeTube"] and self.kwargs["wTube"]:
            self.statusPipe = 2
        elif self.kwargs["DiTube"] and self.kwargs["wTube"]:
            self.statusPipe = 3
        else:
            self.msg = QApplication.translate("pychemqt", "undefined pipe diameters")
            self.status = 0
            return

        if not self.kwargs["kTube"]:
            self.msg = QApplication.translate("pychemqt", "undefined pipe material thermal conductivity")
            self.status = 0
            return

        self.statusFinned = 0
        self.tubefinned = QApplication.translate("pychemqt", "Bare Tube")
        if self.kwargs["tubeFinned"]:
            self.tubefinned = QApplication.translate("pychemqt", "Finned Tube")
            espesor = self.kwargs["thicknessBaseFin"] or self.kwargs["thicknessTopFin"]
            if self.kwargs["hFin"] and espesor:
                self.statusFinned = 1
                self.msg = ""
            else:
                self.msg = QApplication.translate("pychemqt", "fin not specified, using bare tube")
                self.status = 3

        if self.kwargs["modo"]:
            if not self.kwargs["LTube"]:
                self.msg = QApplication.translate("pychemqt", "undefined pipe length")
                self.status = 0
                return
        else:
            self.statusOut = 0
            o1 = self.kwargs["tubeTout"] or self.kwargs["tubeXout"] != -1.
            o2 = self.kwargs["annulliTout"] or self.kwargs["annulliXout"] != -1.
            if o1 and o2:
                self.statusOut = 1
            elif o1:
                self.statusOut = 2
            elif o2:
                self.statusOut = 3
            else:
                self.msg = QApplication.translate("pychemqt", "undefined output condition")
                self.status = 0
                return
        return True

    def calculo(self):
        # Define tipo de flujo
        if self.kwargs["flujo"]:
            self.flujo = "PF"
        else:
            self.flujo = "CF"

        # Calculate pipe dimension
        if self.statusPipe == 1:
            self.De = unidades.Length(self.kwargs["DeTube"])
            self.Di = unidades.Length(self.kwargs["DiTube"])
            self.w = unidades.Length((self.De-self.Di)/2)
            if self.kwargs["wTube"] and w != self.kwargs["wTube"]:
                self.msg = QApplication.translate("pychemqt", "Pipe thickness discard")
                self.status = 3
        elif self.statusPipe == 2:
            self.De = unidades.Length(self.kwargs["DeTube"])
            self.w = unidades.Length(self.kwargs["wTube"])
            self.Di = unidades.Length(self.De-w*2)
        else:
            self.Di = unidades.Length(self.kwargs["DiTube"])
            self.w = unidades.Length(self.kwargs["wTube"])
            self.De = unidades.Length(self.Di+w*2)
        self.Dee = unidades.Length(self.kwargs["DeeTube"])
        self.rugosidad = unidades.Length(self.kwargs["rTube"])
        self.k = unidades.ThermalConductivity(self.kwargs["kTube"])
        self.fi = unidades.Fouling(self.kwargs["tubeFouling"])
        self.fo = unidades.Fouling(self.kwargs["annulliFouling"])

        if self.kwargs["modo"]:
            self.rating()
        else:
            self.design()

        eD = unidades.Dimensionless(self.kwargs["rTube"]/self.Di)
        f = f_friccion(self.ReTube, eD)
        self.deltaPTube = unidades.DeltaP(self.L*self.VTube**2/self.Di*f*self.rhoTube/2)

        f_a = f_friccion(self.ReAnnulli, geometria=6)
        self.deltaPAnnulli = unidades.DeltaP(self.L*self.VAnnulli**2/self.De*f_a*self.rhoAnnulli/2)

        self.salida = [self.outTube.clone(P=self.outTube.P-self.deltaPTube)]
        self.salida.append(self.outAnnulli.clone(P=self.outAnnulli.P-self.deltaPAnnulli))
        self.ToutTube = self.salida[0].T
        self.ToutAnnulli = self.salida[1].T
        self.XoutTube = self.salida[0].x
        self.XoutAnnulli = self.salida[1].x
        self.TinTube = self.kwargs["entradaTubo"].T
        self.XinTube = self.kwargs["entradaTubo"].x
        self.TinAnnulli = self.kwargs["entradaExterior"].T
        self.XinAnnulli = self.kwargs["entradaExterior"].x

    def rating(self):
        """Evaluacion de una tuberia existente"""
        # Input stream
        inTube = self.kwargs["entradaTubo"]
        inAnnulli = self.kwargs["entradaExterior"]

        self.L = unidades.Length(self.kwargs["LTube"])

        # Metodo temperaturas medias
        if self.kwargs["metodo"] == 0:
            self.phaseTube = self.ThermalPhase(inTube, inTube)
            self.phaseAnnulli = self.ThermalPhase(inAnnulli, inAnnulli)

            Ci = inTube.Liquido.cp*inTube.caudalmasico
            Co = inAnnulli.Liquido.cp*inAnnulli.caudalmasico
            Cmin = min(Ci, Co)
            Cmax = max(Ci, Co)
            C_ = Cmin/Cmax

            self.A = unidades.Area(self.L*pi*self.De)
            hi = self.hTube(inTube)
            ho = self.hAnnulli(inAnnulli)
            ni, no = self.rendimientoAletas(hi, ho)
            self.Ug(hi, ni, ho, no)

            NTU = self.A*self.U/Cmin
            ep = self.efectividad(NTU, C_, self.CODE_FLUJO[self.kwargs["flujo"]])
            self.Q = unidades.Power(ep*Cmin*abs(inTube.T-inAnnulli.T))

            if inTube.T > inAnnulli.T:
                QTube = self.Q
                QAnnulli = -self.Q
            else:
                QTube = -self.Q
                QAnnulli = self.Q
            f = lambda T: inTube.clone(T=T).h-inTube.h+QTube
            T = fsolve(f, inTube.T)[0]
            self.outTube = inTube.clone(T=T)
            f = lambda T: inAnnulli.clone(T=T).h-inAnnulli.h+QAnnulli
            T = fsolve(f, inAnnulli.T)[0]
            self.outAnnulli = inAnnulli.clone(T=T)

    def design(self):
        """Diseño a partir de unas condiciones de salida indicadas"""
        # Input stream
        inTube = self.kwargs["entradaTubo"]
        inAnnulli = self.kwargs["entradaExterior"]

        # Metodo temperaturas medias
        if self.kwargs["metodo"] == 0:

            # Calculate output condition and sensible/latent thermal situation,
            # global thermal balance
            if self.statusOut == 1:
                if self.kwargs["tubeTout"]:
                    self.outTube = inTube.clone(T=self.kwargs["tubeTout"])
                else:
                    self.outTube = inTube.clone(x=self.kwargs["tubeXout"])
                if self.kwargs["annulliTout"]:
                    self.outAnnulli = inAnnulli.clone(T=self.kwargs["annulliTout"])
                else:
                    self.outAnnulli = inAnnulli.clone(x=self.kwargs["annulliXout"])

                DTi = abs(self.outTube.T-inTube.T)
                DTo = abs(self.outAnnulli.T-inAnnulli.T)
                Qo = abs(self.outAnnulli.h-inAnnulli.h)
                Qi = abs(self.outTube.h-inTube.h)
                self.Q = unidades.Power((Qo+Qi)/2.)

            elif self.statusOut == 2:
                if self.kwargs["tubeTout"]:
                    self.outTube = inTube.clone(T=self.kwargs["tubeTout"])
                else:
                    self.outTube = inTube.clone(x=self.kwargs["tubeXout"])

                DTi = abs(self.outTube.T-inTube.T)
                Qi = abs(self.outTube.h-inTube.h)
                self.Q = unidades.Power(Qi)
                f = lambda T: inAnnulli.clone(T=T).h-inAnnulli.h-Qi
                T = fsolve(f, inAnnulli.T)[0]
                DTo = abs(T-inAnnulli.T)
                self.outAnnulli = inAnnulli.clone(T=T)

            elif self.statusOut == 3:
                if self.kwargs["annulliTout"]:
                    self.outAnnulli = inAnnulli.clone(T=self.kwargs["annulliTout"])
                else:
                    self.outAnnulli = inAnnulli.clone(x=self.kwargs["annulliXout"])

                DTo = abs(self.outAnnulli.T-inAnnulli.T)
                Qo = abs(self.outAnnulli.h-inAnnulli.h)
                self.Q = unidades.Power(Qo)
                f = lambda T: inTube.clone(T=T).h-inTube.h-Qi
                T = fsolve(f, inTube.T)[0]
                DTi = abs(T-inTube.T)
                self.outTube = inTube.clone(T=T)

            self.phaseTube = self.ThermalPhase(inTube, self.outTube)
            self.phaseAnnulli = self.ThermalPhase(inAnnulli, self.outAnnulli)

            fluidTube = inTube.clone(T=(inTube.T+self.outTube.T)/2.)
            fluidAnnulli = inAnnulli.clone(T=(inAnnulli.T+self.outAnnulli.T)/2.)

            hi = self.hTube(fluidTube)
            ho = self.hAnnulli(fluidAnnulli)
            ni, no = self.rendimientoAletas(hi, ho)
            self.Ug(hi, ni, ho, no)

            DTi = abs(self.kwargs["tubeTout"]-inTube.T)
            DTo = abs(self.kwargs["annulliTout"]-inAnnulli.T)

            if self.kwargs["flujo"]:
                DTin = abs(inAnnulli.T-inTube.T)
                DTout = abs(self.kwargs["tubeTout"]-self.kwargs["annulliTout"])
            else:
                DTin = abs(self.kwargs["tubeTout"]-inAnnulli.T)
                DTout = abs(self.kwargs["annulliTout"]-inTube.T)
            if DTin == DTout:
                DTm = DTin
            else:
                DTm = (DTin-DTout)/log(DTin/DTout)

        self.A = unidades.Area(self.Q/self.U/DTm)
        self.L = unidades.Length(self.A/2/pi)


    def Ug(self, hi, ni, ho, no):
        U = 1/(self.De/self.Di/hi/ni+self.De*self.fi/self.Di/ni+self.De*log(self.De/self.Di)/2/self.k+self.fo/no+1/ho/no)
        Uc = 1/(self.De/self.Di/hi/ni+self.De*log(self.De/self.Di)/2/self.k+1/ho/no)
        self.U = unidades.HeatTransfCoef(U)
        self.CF = unidades.Dimensionless(U/Uc)
        self.OS = unidades.Dimensionless(Uc*(self.fi+self.fo))


    def ThermalPhase(self, input, output):
        # Calculate thermal fundamentals
        if input.x == output.x:
            if input.x == 0:
                phase = "Latent-Liquid"
            else:
                phase = "Latent-Vapor"
        elif input.x < output.x:
            phase = "Evaporator"
        else:
            phase = "Condenser"
        return phase

    def rendimientoAletas(self, hi, ho):
        """Metodo de calculo del rendimiento termico de aletas"""
        if self.kwargs["tubeFinned"]:
            if self.statusFinned:
                # de momento solo asumimos el caso de aleta circular
                do = self.kwargs["rootDoFin"]
                D = do+self.kwargs["hFin"]*2
                phi = (D/do-1)*(1+0.35*log(D/do))
                if self.kwargs["thicknessBaseFin"] and self.kwargs["thicknessTopFin"]:
                    delta = (self.kwargs["thicknessBaseFin"]+self.kwargs["thicknessTopFin"])/2
                else:
                    delta = self.kwargs["thicknessBaseFin"]+self.kwargs["thicknessTopFin"]
                if self.kwargs["kFin"]:
                    kf = self.kwargs["kFin"]
                else:
                    kf = self.kwargs["kTube"]

                X = phi*do/2*sqrt(2*ho/kf/delta)
                no = tanh(X)/X
            else:
                no = 1
            # Falta definir aletas interiores
            ni = 1
        else:
            ni = 1
            no = 1
        return ni, no


    def hTube(self, fluidTube):
        """Cálculo del coeficiente de transferencia de calor por conveccion en la parte del tubo"""
        if self.phaseTube[:6] == "Latent":
            if fluidTube.x == 0:
                fluido = fluidTube.Liquido
            else:
                fluido = fluidTube.Vapor

            rho_i = fluido.rho
            mu = fluido.mu
            k = fluido.k
            v_i = fluidTube.Q*4/pi/self.Di**2
            re_i = Re(D=self.Di, V=v_i, rho=rho_i, mu=mu)
            self.VTube = unidades.Speed(v_i)
            self.rhoTube = rho_i
            self.ReTube = unidades.Dimensionless(re_i)
            pr = fluido.Prandt

            if re_i < 2300:
                L = self.L
                cp = fluido.cp
                w = fluido.caudalmasico
                gz = Gz(w=w, cp=cp, k=k, L=L)
                beta = fluido.alfav
                gr = Gr(beta=beta, T=fluidTube.T, To=fluidTube.T, L=L, mu=mu)
                if self.kwargs["tubesideLaminar"] == 0:
                    Nu = h_tubeside_laminar_Eubank_Proctor(Pr=pr, Gz=gz, Gr=gr, D=self.Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 1:
                    Nu = h_tubeside_laminar_VDI(Re=re_i, Pr=pr, D=self.Di, L=L)
                elif self.kwargs["tubesideLaminar"] == 2:
                    Nu = h_tubeside_laminar_Hausen(Gz=gz)
                elif self.kwargs["tubesideLaminar"] == 3:
                    Nu = h_tubeside_laminar_Sieder_Tate(Gz=gz, Gr=gr)
            else:
                if self.kwargs["tubesideTurbulent"] == 0:
                    Nu = h_tubeside_turbulent_Sieder_Tate(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 1:
                    Nu = h_tubeside_turbulent_Colburn(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 2:
                    frio = self.kwargs["entradaCarcasa"].T > fluidTube.T
                    Nu = h_tubeside_turbulent_Dittus_Boelter(Re=re_i, Pr=pr, calentamiento=frio)
                elif self.kwargs["tubesideTurbulent"] == 3:
                    Nu = h_tubeside_turbulent_ESDU(Re=re_i, Pr=pr)
                elif self.kwargs["tubesideTurbulent"] == 4:
                    Nu = h_tubeside_turbulent_Gnielinski(Re=re_i, Pr=pr, D=self.Di, L=L)
                elif self.kwargs["tubesideTurbulent"] == 5:
                    line = self.kwargs["distribucionTube"] == 3
                    filas = self.kwargs["NTube"]**0.5
                    Nu = h_tubeside_turbulent_VDI(Re=re_i, Pr=pr, filas_tubos=filas, alineados=line)

        if self.phaseTube == "Condenser":
            if self.kwargs["orientation"] == 0:  # Condensacion en tubos horizontales
                if 0 < fluidTube.x < 1:
                    X_lockhart = ((1-fluidTube.x)/fluidTube.x)**0.9*(fluidTube.Vapor.rho/fluidTube.Liquido.rho)**0.5*(fluidTube.Liquido.mu/fluidTube.Vapor.mu)**0.1
                    G = fluidTube.caudalmasico*4/pi/self.Di**2
                    j = fluidTube.x*G/(g*self.Di*fluidTube.Vapor.rho*(fluidTube.Liquido.rho-fluidTube.Vapor.rho))**0.5
                    print((j, X_lockhart))

            else:   # condensacion veritcal
                pass

        return unidades.HeatTransfCoef(Nu*k/self.Di)

    def hAnnulli(self, fluidAnnulli):
        """Cálculo del coeficiente de transferencia de calor por conveccion en la parte del anillo"""
        a = self.Dee/self.De
        dh = self.Dee-self.De

        rho = fluidAnnulli.Liquido.rho
        mu = fluidAnnulli.Liquido.mu
        k = fluidAnnulli.Liquido.k
        v = fluidAnnulli.Q*4/pi/(self.Dee**2-self.De**2)
        re = Re(D=dh, V=v, rho=rho, mu=mu)
        self.VAnnulli = unidades.Speed(v)
        self.rhoAnnulli = rho
        self.ReAnnulli = unidades.Dimensionless(re)
        pr = fluidAnnulli.Liquido.Prandt

        if re <= 2300:
            Nu = h_anulli_Laminar(re, pr, a)
        elif re >= 1e4:
            Nu = h_anulli_Turbulent(re, pr, a)
        else:
            Nu = h_anulli_Transition(re, pr, a)

        return unidades.HeatTransfCoef(Nu*k/self.Di)

    def coste(self):
        self.material = self.kwargs["material"]

        if self.kwargs["P_dis"]:
            self.P_dis = unidades.Pressure(self.kwargs["P_dis"])
        else:
            self.P_dis = max(self.kwargs["entradaTubo"].P, self.kwargs["entradaExterior"].P)

        Pd = self.P_dis.psi
        Fm = [1., 1.9, 2.2][self.kwargs["material"]]

        if Pd < 4:
            Fp = 1.
        elif Pd < 6:
            Fp = 1.1
        else:
            Fp = 1.25

        C = Fm*Fp*900*self.A.ft2**0.18
        self.C_adq = unidades.Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst = unidades.Currency(self.C_adq * self.kwargs["f_install"])

    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Catalog")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Length"), self.L.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pipe Internal Diameter"), self.Di.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pipe External Diameter"), self.De.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli External Diameter"), self.Dee.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thickness"), self.w.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Roughness"), self.rugosidad.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Area"), self.A.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thermal Conductivity"), self.k.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Internal Fouling"), self.fi.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "External Fouling"), self.fo.str)+os.linesep

        if self.kwargs["tubeFinned"]:
            txt+="%s" %(QApplication.translate("pychemqt", "Finned Tube"))+os.linesep
#        hFin: Altura de la aleta
#        thicknessBaseFin: Espesor en la base de la aleta
#        thicknessTopFin: Espesor en lo alto de la aleta
#        rootDoFin: Diametro externo en la base de las aletas
#        kFin: Conductividad termica del material de la aleta
#        nFin: Numero de aletas por metro de tuberia

        else:
            txt+="%s" %(QApplication.translate("pychemqt", "Bare Tube"))+os.linesep

        txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_MODO[self.kwargs["modo"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Arrangement Flow"), self.TEXT_FLUJO[self.kwargs["flujo"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Layout"), self.TEXT_ORIENTACION[self.kwargs["orientacion"]])+os.linesep

        txt+=os.linesep+"%-25s\t %s" %(QApplication.translate("pychemqt", "Tube Mechanism"), self.phaseTube)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube Fluid Speed"), self.VTube.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube Reynolds"), self.ReTube.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube In Temperature"), self.kwargs["entradaTubo"].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube In Quality"), self.kwargs["entradaTubo"].x.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube Out Temperature"), self.ToutTube.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Tube Out Quality"), self.XoutTube.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP Tube", None), self.deltaPTube.str)+os.linesep

        txt+=os.linesep+"%-25s\t %s" %(QApplication.translate("pychemqt", "Annulli Mechanism"), self.phaseAnnulli)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli Fluid Speed"), self.VAnnulli.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli Reynolds"), self.ReAnnulli.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli In Temperature"), self.kwargs["entradaExterior"].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli In Quality"), self.kwargs["entradaExterior"].x.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli Out Temperature"), self.ToutAnnulli.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Annulli Out Quality"), self.XoutAnnulli.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "ΔP Annulli", None), self.deltaPAnnulli.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "U"), self.U.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Clean Factor"), self.CF.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Over Surface"), self.OS.str)+os.linesep

        if self.statusCoste:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Material"), self.TEXT_MATERIAL[self.kwargs["material"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Design Pressure"), self.P_dis.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt


    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Length"), "L", unidades.Length),
             (QApplication.translate("pychemqt", "Pipe Internal Diameter"), "Di", unidades.Length),
             (QApplication.translate("pychemqt", "Pipe External Diameter"), "De", unidades.Length),
             (QApplication.translate("pychemqt", "Annulli External Diameter"), "Dee", unidades.Length),
             (QApplication.translate("pychemqt", "Thickness"), "w", unidades.Length),
             (QApplication.translate("pychemqt", "Roughness"), "rugosidad", unidades.Length),
             (QApplication.translate("pychemqt", "External Area"), "A", unidades.Area),
             (QApplication.translate("pychemqt", "Thermal Conductivity"), "k", unidades.ThermalConductivity),
             (QApplication.translate("pychemqt", "Internal Fouling"), "fo", unidades.Fouling),
             (QApplication.translate("pychemqt", "External Fouling"), "fi", unidades.Fouling),
             (QApplication.translate("pychemqt", "Finned Tube"), "tubefinned", str),
             (QApplication.translate("pychemqt", "Mode"), ("TEXT_MODO", "modo"), str),
             (QApplication.translate("pychemqt", "Arrangement Flow"), ("TEXT_FLUJO", "flujo"), str),
             (QApplication.translate("pychemqt", "Layout"), ("TEXT_ORIENTACION", "orientacion"), str),
             (QApplication.translate("pychemqt", "Tube Mechanism"), "phaseTube", str),
             (QApplication.translate("pychemqt", "Tube Fluid Speed"), "VTube", unidades.Speed),
             (QApplication.translate("pychemqt", "Tube Reynolds"), "ReTube", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Tube In Temperature"), "TinTube", unidades.Temperature),
             (QApplication.translate("pychemqt", "Tube In Quality"), "XinTube", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Tube Out Temperature"), "ToutTube", unidades.Temperature),
             (QApplication.translate("pychemqt", "Tube Out Quality"), "XoutTube", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "ΔP Tube", None), "deltaPTube", unidades.DeltaP),
             (QApplication.translate("pychemqt", "Annulli Mechanism"), "phaseAnnulli", str),
             (QApplication.translate("pychemqt", "Annulli Fluid Speed"), "VAnnulli", unidades.Speed),
             (QApplication.translate("pychemqt", "Annulli Reynolds"), "ReAnnulli", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Annulli In Temperature"), "TinAnnulli", unidades.Temperature),
             (QApplication.translate("pychemqt", "Annulli In Quality"), "XinAnnulli", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Annulli Out Temperature"), "ToutAnnulli", unidades.Temperature),
             (QApplication.translate("pychemqt", "Annulli Out Quality"), "XoutAnnulli", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "ΔP Annulli", None), "deltaPAnnulli", unidades.DeltaP),
             (QApplication.translate("pychemqt", "U"), "U", unidades.HeatTransfCoef),
             (QApplication.translate("pychemqt", "Clean Factor"), "CF", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Over Surface"), "OS", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Material"), ("TEXT_MATERIAL", "material"), str),
             (QApplication.translate("pychemqt", "Design Pressure"), "P_dis", unidades.Pressure),
             (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", unidades.Currency),
             (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", unidades.Currency)]
        return l


class Refrigeration(equipment):
    """Clase que modela las unidades de refrigeración"""
    title = "Refrigerador"
    help = ""

    def Coste(self, *args, **kwargs):
        self._indicesCoste(*args)

        Tmax = self.Tmax.C
        Q = self.calor.MBtuh

        if Tmax >= 0:
            F = 1.
        elif Tmax >= -10:
            F = 1.55
        elif Tmax >= -20:
            F = 2.1
        elif Tmax >= -30:
            F = 2.65
        elif Tmax >= -3:
            F = 3.2
        else:
            F = 4.

        C = 146*F*Q**0.65*1000
        C_adq = C*self.Current_index/self.Base_index
        C_inst = C_adq*self.kwargs["f_install"]

        self.C_adq = C_adq
        self.C_inst = C_inst


class Fired_Heater(equipment):
    """Clase que modela los equipos de calentamiento por fuego directo

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        Tout: requerido, temperatura que queremos que alcance la coriente a calentar
        deltaP: Perdida de presión en el equipo
        Hmax: Flujo de calor de diseño del equipo, el valor del calor transmitido no podrá ser superior a este valor
        eficiencia: eficiencia térmica en la combustión, por defecto 75%
        poderCalorifico: poder calorífico del combustible, por defecto 900 Btu/stdft³

    Coste:
        tipo:
            0   -   Box type
            1   -   Cylindrical type
        subtipoBox:
            0   -   Process heater
            1   -   Pyrolysis
            2   -   Reformer without catalyst
        subtipoCylindrical:
            0   -   Cylindrical
            1   -   Dowtherm
        material:
            0   -   Carbon steel
            1   -   CrMo steel
            2   -   Stainless
        P_dis: Presión de diseño del equipo
    """
    title = QApplication.translate("pychemqt", "Fired Heater")
    help = os.environ["pychemqt"] + "doc/fireHeater.htm"
    kwargs = {
        "entrada": None,
        "Tout": 0.0,
        "deltaP": 0.0,
        "Hmax": 0.0,
        "eficiencia": 0.0,
        "poderCalorifico": 0.0,

        "f_install": 1.3,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "tipo": 0,
        "subtipoBox": 0,
        "subtipoCylindrical": 0,
        "material": 0,
        "P_dis": 0.0}

    kwargsInput = ("entrada", )
    kwargsValue = ("Tout", "deltaP", "Hmax", "eficiencia", "poderCalorifico",
                   "P_dis")
    kwargsList = ("tipo", "subtipoBox", "subtipoCylindrical", "material")
    calculateValue = ("CombustibleRequerido", "Heat")
    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 3
    salida = [None]

    TEXT_TIPO = [QApplication.translate("pychemqt", "Box Type"),
                 QApplication.translate("pychemqt", "Cylindrical type")]
    TEXT_SUBTIPOBOX = [
        QApplication.translate("pychemqt", "Process heater"),
        QApplication.translate("pychemqt", "Pyrolysis"),
        QApplication.translate("pychemqt", "Reforme without catalysis")]
    TEXT_SUBTIPOCYLINDRICAL = [
        QApplication.translate("pychemqt", "Cylindrical"),
        QApplication.translate("pychemqt", "Dowtherm")]
    TEXT_MATERIAL = [QApplication.translate("pychemqt", "Carbon steel"),
                     QApplication.translate("pychemqt", "Cr-Mo steel"),
                     QApplication.translate("pychemqt", "Stainless steel")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            return

        if not self.kwargs["Tout"]:
            self.msg = QApplication.translate("pychemqt", "undefined output temperature condition")
            self.status = 0
            return
        if self.kwargs["Tout"]<=self.kwargs["entrada"].T:
            self.msg = QApplication.translate("pychemqt", "bad output temperature condition")
            self.status = 0
            return

        if not self.kwargs["eficiencia"]:
            self.msg = QApplication.translate("pychemqt", "using default efficiency")
            self.status = 3
            return True
        if not self.kwargs["poderCalorifico"]:
            self.msg = QApplication.translate("pychemqt", "using default fuel calorific value")
            self.status = 3
            return True

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):
        entrada = self.kwargs["entrada"]
        self.Tout = unidades.Temperature(self.kwargs["Tout"])
        self.deltaP = unidades.DeltaP(self.kwargs["deltaP"])
        self.Hmax = unidades.Power(self.kwargs["Hmax"])
        if self.kwargs["eficiencia"]:
            eficiencia = self.kwargs["eficiencia"]
        else:
            eficiencia = 0.75

        if self.kwargs["poderCalorifico"]:
            poderCalorifico = self.kwargs["poderCalorifico"]
        else:
            poderCalorifico = 900

        salida = entrada.clone(T=self.Tout, P=entrada.P-self.deltaP)
        Ho = entrada.h
        H1 = salida.h
        Heat = unidades.Power(H1-Ho)

        if self.Hmax and Heat > self.Hmax:
            self.Heat = unidades.Power(self.Hmax)
            To = (entrada.T+self.Tout)/2
            T = fsolve(lambda T: entrada.clone(T=T, P=entrada.P-self.deltaP).h-Ho-self.Hmax, To)[0]
            self.salida = [entrada.clone(T=T, P=entrada.P-self.deltaP)]
        else:
            self.Heat = Heat
            self.salida = [salida]

        self.CombustibleRequerido = unidades.VolFlow(self.Heat.Btuh/poderCalorifico/eficiencia, "ft3h")
        self.deltaT = unidades.DeltaT(self.salida[0].T-entrada.T)
        self.eficiencia = unidades.Dimensionless(eficiencia)
        self.poderCalorifico = unidades.Dimensionless(poderCalorifico)
        self.Tin = entrada.T
        self.Tout = self.salida[0].T

    def coste(self):
        """
        tipo:
            0   -   Box type
            1   -   Cylindrical type
        subtipoBox:
            0   -   Process heater
            1   -   Pyrolysis
            2   -   Reformer without catalyst
        subtipoCylindrical:
            0   -   Cylindrical
            1   -   Dowtherm
        material:
            0   -   Carbon steel
            1   -   CrMo steel
            2   -   Stainless
        P_dis: Presión de diseño del equipo
        """
        if self.kwargs["P_dis"]:
            P_dis = unidades.Pressure(self.kwargs["P_dis"])
        else:
            P_dis = self.kwargs["entrada"].P

        if self.kwargs["tipo"] == 0:  # Boxtype
            k = [25.5, 33.8, 45.][self.kwargs["material"]]
            Fd = [0, 0.1, 0.35][self.kwargs["subtipoBox"]]

            if P_dis.psi <= 500:
                Fp = 0
            elif P_dis.psi <= 1000:
                Fp = 0.1
            elif P_dis.psi <= 1500:
                Fp = 0.15
            elif P_dis.psi <= 2000:
                Fp = 0.25
            elif P_dis.psi <= 2500:
                Fp = 0.4
            else:
                Fp = 0.6

            C = k*(1+Fd+Fp)*self.Heat.MBtuh**0.86*1000

        else:  # Cylindrical
            k = [27.3, 40.2, 42.][self.kwargs["material"]]
            Fd = [0, 0.33][self.kwargs["subtipoCylindrical"]]

            if P_dis.psi <= 500:
                Fp = 0
            elif P_dis.psi <= 1000:
                Fp = 0.15
            else:
                Fp = 0.2

            C = k*(1+Fd+Fp)*self.Heat.MBtuh**0.82*1000

        self.P_dis = P_dis
        self.C_adq = unidades.Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst = unidades.Currency(self.C_adq * self.kwargs["f_install"])

    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Temperature"), self.Tin.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Temperature increase"), self.deltaT.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure increase"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Maximum heat"), self.Hmax.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thermal Efficiency"), self.eficiencia.str)
        if self.kwargs["eficiencia"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        else:
            txt+=os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Fuel Heating Value"), self.poderCalorifico.str)
        if self.kwargs["poderCalorifico"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        else:
            txt+=os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Required Fuel"), self.CombustibleRequerido.str)+os.linesep

        if self.statusCoste:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep

            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "FireHeater type"), self.TEXT_TIPO[self.kwargs["tipo"]])+os.linesep
            if self.kwargs["tipo"]:
                txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Cylindrical type"), self.TEXT_SUBTIPOCYLINDRICAL[self.kwargs["subtipoCylindrical"]])+os.linesep
            else:
                txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Box type"), self.TEXT_SUBTIPOBOX[self.kwargs["subtipoBox"]])+os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Material"), self.TEXT_MATERIAL[self.kwargs["material"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Design Pressure"), self.P_dis.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Input Temperature"), "Tin", unidades.Temperature),
             (QApplication.translate("pychemqt", "Output Temperature"), "Tout", unidades.Temperature),
             (QApplication.translate("pychemqt", "Temperature increase"), "deltaT", unidades.DeltaT),
             (QApplication.translate("pychemqt", "Pressure increase"), "deltaP", unidades.DeltaP),
             (QApplication.translate("pychemqt", "Maximum heat"), "Hmax", unidades.Power),
             (QApplication.translate("pychemqt", "Thermal Efficiency"), "eficiencia", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Fuel Heating Value"), "poderCalorifico", unidades.Dimensionless),
             (QApplication.translate("pychemqt", "Required Fuel"), "CombustibleRequerido", unidades.VolFlow),
             (QApplication.translate("pychemqt", "FireHeater type"), ("TEXT_TIPO", "tipo"), str),
             (QApplication.translate("pychemqt", "Cylindrical type"), ("TEXT_SUBTIPOCYLINDRICAL", "subtipoCylindrical"), str),
             (QApplication.translate("pychemqt", "Box type"), ("TEXT_SUBTIPOBOX", "subtipoBox"), str),
             (QApplication.translate("pychemqt", "Material"), ("TEXT_MATERIAL", "material"), str),
             (QApplication.translate("pychemqt", "Design Pressure"), "P_dis", unidades.Pressure),
             (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", unidades.Currency),
             (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", unidades.Currency)]
        return l

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["Tout"] = self.Tout
        state["deltaP"] = self.deltaP
        state["Hmax"] = self.Hmax
        state["Heat"] = self.Heat
        state["CombustibleRequerido"] = self.CombustibleRequerido
        state["deltaT"] = self.deltaT
        state["eficiencia"] = self.eficiencia
        state["poderCalorifico"] = self.poderCalorifico
        state["Tin"] = self.Tin
        state["statusCoste"] = self.statusCoste
        if self.statusCoste:
            state["P_dis"] = self.P_dis
            state["C_adq"] = self.C_adq
            state["C_inst"] = self.C_inst

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.Tout = unidades.Temperature(state["Tout"])
        self.deltaP = unidades.DeltaP(state["deltaP"])
        self.Hmax = unidades.Power(state["Hmax"])
        self.Heat = unidades.Power(state["Heat"])
        self.CombustibleRequerido = unidades.VolFlow(state["CombustibleRequerido"])
        self.deltaT = unidades.DeltaT(state["deltaT"])
        self.eficiencia = unidades.Dimensionless(state["eficiencia"])
        self.poderCalorifico = unidades.Dimensionless(state["poderCalorifico"])
        self.Tin = unidades.Temperature(state["Tin"])
        self.statusCoste = state["statusCoste"]
        if self.statusCoste:
            self.P_dis = unidades.Pressure(state["P_dis"])
            self.C_adq = unidades.Currency(state["C_adq"])
            self.C_inst = unidades.Currency(state["C_inst"])
        self.salida = [None]



def unsteady():
#========================================================================
#Input constants
#========================================================================
    dtau = 0.001 # Set dimensionless time increments
    dx = 0.05 # Set dimensionless length increments
    Tmax = 0.95 # Set maximum dimensionless temperature
    M = 21 # Counter for length discretization
#========================================================================
#Calculate parameters
#========================================================================
    dx = 1.0/(M-1)
    dx_x = 1.0/(M-1)
    ratio = dtau/(dx**2)
    const = 1.0 - 2.0*ratio
#========================================================================
#Set counters to zero
#========================================================================
    i = 0
    tau = 0.0
#========================================================================
# Set up arrays for solution and print
#========================================================================
    Tnew = zeros(M, dtype = float)
    T = zeros(M, dtype = float)
    T[0] = 1.0
    T[-1] = 1.0
    print(("T initial = ", T))
#========================================================================
# I just pick 400 on trial and error for the total array
#========================================================================
    T_sol = zeros((400,M), dtype = float)
    T_sol[:,0] = 1.0
    T_sol[:,-1] = 1.0
#========================================================================
# While loop to iterate until mid-point temperature reaches Tmax
#========================================================================
    while T[10] < Tmax:
        i = i + 1
        tau = tau + dtau
#========================================================================
# Calculate new tempertures
#========================================================================
        for j in range(1,M-1):
            Tnew[j] = ratio*(T[j-1] + T[j+1]) + const*T[j]
#========================================================================
# Substitute new Temperatures in array for T
#========================================================================
        for k in range(1,M-1):
            T[k] = Tnew[k]
            T_sol[i,k] = T[k]

    print(("Tau and T_final =", tau, T_sol[i,:]))
#========================================================================
# Set up array for spatial values of x to plot
#========================================================================
    x = [i*dx_x for i in range(M-1)]
    x.append(1.0)
#========================================================================
# Plot the solutions
#=======================================================================
    plot(x,T_sol[50,:])
    plot(x,T_sol[100,:])
    plot(x,T_sol[150,:])
    plot(x,T_sol[250,:])
    plot(x,T_sol[i,:])
#legend(['Tau = 0.5','Tau = 0.1','Tau = 0.15','Tau = 0.25',
#'Tau = final time'])
    title('Normalized Slab Temperatures')
    plt.show()
    grid()

if __name__ == "__main__":
#    agua=Corriente(T=400, P=101325., caudalMasico=1., fraccionMolar=[1.])
#    cambiador=Fired_Heater(entrada=agua, Tout=450)
#    print cambiador.Heat.MJh

#    agua=Corriente(T=300, P=101325., caudalMasico=1., fraccionMolar=[1.])
#    Cambiador=Heat_Exchanger(entrada=agua, Tout=350)
#    print Cambiador.Heat.MJh

    aguaTubo=Corriente(T=283.15, P=101325., ids=[62], caudalMasico=10., fraccionMolar=[1.])
    aguaCarcasa=Corriente(T=370, P=101325., ids=[62], caudalMasico=1000., fraccionMolar=[1.])
    Ds=unidades.Length(19.25, "inch")
    Do=unidades.Length(1, "inch")
    L=unidades.Length(14, "ft")
    pt=unidades.Length(1.25, "inch")
    B=unidades.Length(3.85, "inch")
    dotl=unidades.Length(0.5*(19.25-17.91), "inch")
    cambiador = Shell_Tube(entradaTubo=aguaTubo, entradaCarcasa=aguaCarcasa, shellsideSensible=1, DShell=Ds, NTube=124, DeTube=Do, wTube=0.006, LTube=L, pitch=pt, distribucionTube=3, baffleSpacing=B, baffleSpacingIn=B, baffleSpacingOut=B, baffleCut=0.2, clearanceTubeBaffle=0.0004, clearanceShellBaffle=0.0025279, clearanceShellBundle=dotl, sealingStrips=0.1*9.24)

#    caliente=Corriente(T=140+273.15, P=361540., caudalMasico=1.36, ids=[62], fraccionMolar=[1.])
#    fria=Corriente(T=20+273.15, P=101325., caudalMasico=5000/3600., ids=[62], fraccionMolar=[1.])
#    Cambiador=Hairpin(modo=1)
#    print Cambiador.status, Cambiador.msg
#    Cambiador(entradaTubo=caliente)
#    print Cambiador.status, Cambiador.msg
#    Cambiador=Hairpin(entradaTubo=caliente, entradaExterior=fria, modo=1,
#                      DiTube=0.0525, DeTube=0.0603, DeeTube=0.0779, kTube=54, rTube=0.0459994e-3,
#                      annulliFouling= 0.000352, tubeFouling=0.000176, LTube=2.5)

#    from scipy import *
#    from pylab import *
#    import matplotlib.pyplot as plt
#    unsteady()
