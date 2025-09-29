#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


from numpy import logspace

from lib.unidades import Dimensionless
from tools.qt import translate
from equipment.parents import equipment


def h_shellside_turbulent_Bell_Delaware(**kw):
    """Coeficiente de transferencia de calor por calor sensible en el parte de la carcasa en regimen turbulento
    Serth - Process Heat Transfer - Principles and applications Cap. 6
    kw include
    DeTube (d) - Tube outside diameter
    tipo - Tube layout pattern
    Ds - Shell inside diameter
    Dotl - Tube bank outer tube limit diameter
    Ltl - Effective tube length
    Bc - Baffle cut as a percent of D
    Lbc - Central baffle spacing
    Nss - Number of sealing strips per side
    """

    # Bundle-to-shell clearance
    Lbb = 12+0.005*kw["Ds"]

    # Bundle diameter
    Dotl = kw["Ds"] - Lbb



    fi = 1

    PtD = kw["pitch"]/kw["DeTube"]
    mo = kw["entradaCarcasa"].caudalmasico
    if kw["distribucionTube"] == 2:
        Ptef = kw["pitch"]/2**0.5
    else:
        Ptef = kw["pitch"]

    if kw["distribucionTube"] < 3:
        tita_tp = unidades.Angle([30, 45, 60][kw["distribucionTube"]], "deg")
        Pt_ = kw["pitch"]*cos(tita_tp)
    else:
        Pt_ = kw["pitch"]
    Nc = kw["Ds"]*(1-2*kw["baffleCut"])/Pt_
    Ncw = 0.8*kw["baffleCut"]*kw["Ds"]/Pt_

    Dotl = kw["Ds"]-2*kw["clearanceShellBundle"]
    Sm = kw["baffleSpacing"]*(kw["Ds"]-Dotl+(Dotl-kw["DeTube"])/Ptef*(kw["pitch"]-kw["DeTube"]))

    G = mo/Sm
    Re = kw["DeTube"]*G/kw["entradaCarcasa"].Liquido.mu
    Pr = kw["entradaCarcasa"].Liquido.Pr
    rho = kw["entradaCarcasa"].Liquido.rho
    cp = kw["entradaCarcasa"].Liquido.cp

    Dctl = Dotl-kw["DeTube"]
    tita_ctl = 2*arccos(kw["Ds"]*(1-2*kw["baffleCut"])/Dctl)
    Fc = 1+1./pi*(sin(tita_ctl)-tita_ctl)
    Fw = 1./2/pi*(tita_ctl-sin(tita_ctl))
    Stb = pi/8*((kw["DeTube"]+2*kw["clearanceTubeBaffle"])**2-kw["DeTube"]**2)*kw["NTube"]*(1+Fc)

    tita_ds = 2*arccos(1-2*kw["baffleCut"])
    Ssb = kw["Ds"]*kw["clearanceShellBaffle"]*(pi-0.5*tita_ds)

    Sb = kw["baffleSpacing"]*(kw["Ds"]-Dotl)
    Sw = 1./8*kw["Ds"]**2*(tita_ds-sin(tita_ds))-1./4*kw["NTube"]*Fw*pi*kw["DeTube"]**2
    Dw = 4*Sw/(pi*kw["DeTube"]*kw["NTube"]*0.5*(1-Fc)+kw["Ds"]*tita_ds)

    # Jc: Correction factor for baffle configuration
    Jc = 0.55+0.72*Fc

    # Jl: Correction factor for baffle leakage effects, including both
    # tube-to-baffle and baffle-to-shell leakages.
    rs = Ssb/(Ssb+Stb)
    rl = (Ssb+Stb)/Sm
    Jl = 0.44*(1-rs)+(1-0.44*(1-rs))*exp(-2.2*rl)
    Rl = exp(-1.33*(1+rs)*rl**(0.8-0.15*(1+rs)))

    # Jb: Correction factor for bundle and pass partition bypass stream
    rss = kw["sealingStrips"]/Nc
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

    # Js: Correction factor for larger baffle spacing at the inlet and outlet
    # sections compared to the central baffle spacing.
    if Re < 100:
        n1 = 1./3
        n2 = 1.
    else:
        n1 = 0.6
        n2 = 0.2
    nb = kw["LTube"]/kw["baffleSpacing"]
    if kw["baffleSpacingIn"] > kw["baffleSpacing"]:
        nb += 1
    if kw["baffleSpacingOut"] > kw["baffleSpacing"]:
        nb += 1
    Js = (nb-1+(kw["baffleSpacingIn"]/kw["baffleSpacing"])**(1-n1)+(kw["baffleSpacingOut"]/kw["baffleSpacing"])**(1-n1))/(nb-1+(kw["baffleSpacingIn"]/kw["baffleSpacing"])+kw["baffleSpacingOut"]/kw["baffleSpacing"])
    Rs = 0.5*((kw["baffleSpacing"]/kw["baffleSpacingIn"])**(2-n2)+(kw["baffleSpacing"]/kw["baffleSpacingOut"])**(2-n2))

    # Jr: Correction factor for any adverse temperature gradient in laminar
    # flow. Only applies for Re < 100
    Nct = (nb+1)*(Nc+Ncw)
    if Re <= 20:
        Jr = (10/Nct)**0.18
    elif Re >= 100:
        Jr = 1.
    else:  # Interpolacion entre los valores de arriba
        Jr = 0.853379+0.0014662*Re

    j = jFactor(Re, PtD, kw["distribucionTube"])
    f = fFactor(Re, PtD, kw["distribucionTube"])
    hid = j*cp*G*fi/Pr**(2./3)
    h = unidades.HeatTransfCoef(hid*Jc*Jl*Jb*Jr*Js)

    DPideal = 2*f*Nc*G**2/g/rho/fi
    DPc = (nb-1)*DPideal*Rl*Rb
    if Re > 100:
        DPwideal = (2+0.6*Ncw)*mo**2/2/g/rho/Sm/Sw
    else:
        DPwideal = 26*nu*mo/g/(Sm*Sw)**0.5*(Ncw/PtD+kw["baffleCut"]*kw["Ds"]/Dw**2)+mo**2/g/rho/Sm/Sw
    DPw = nb*DPwideal*Rl

    # Eq 6.15
    DPe = 2*DPideal*(1+Ncw/Nc)*Rb*Rs

    mon = mo/kw["parallel"]

    DPn = 0
    for D in kw["nozzleInShellsideDiameter"], kw["nozzleOutShellsideDiameter"]:
        Ren = D*mon/kw["entradaCarcasa"].Liquido.mu
        Sn = pi/4*D**2
        if Ren > 100:
            DPn += 2e-13*kw["serie"]*mon**2/Sn
        else:
            DPn += 4e-13*kw["serie"]*mon**2/Sn
    DP = unidades.DeltaP(DPc+DPw+DPe+DPn)
    return h, DP


def jFactor(Re, PtD, tubeConf):
    """
    Parameters
    ----------
    Re : float
        Reynolds number
    PtD : float
        Ratio between tube pitch and tube OD
    tubeConf : str
        Code with tube bundle configuration, square|triangular|rotSquare

    Return
    ------
    j : float
        adimensional heat transfer coefficient
    """
    # Table 6.1
    if tubeConf == "rotSquare":
        a3, a4 = 1.93, 0.5
        if Re < 10:
            a1, a2 = 1.55, -0.667
        elif Re < 100:
            a1, a2 = 0.498, -0.656
        elif Re < 1000:
            a1, a2 = 0.73, -0.500
        else:
            a1, a2 = 0.37, -0.396

    elif tubeConf == "square":
        a3, a4 = 1.187, 0.37
        if Re < 10:
            a1, a2 = 0.97, -0.667
        elif Re < 100:
            a1, a2 = 0.9, -0.631
        elif Re < 1000:
            a1, a2 = 0.408, -0.46
        elif Re < 10000:
            a1, a2 = 0.107, -0.266
        else:
            a1, a2 = 0.37, -0.395

    else:
        a3, a4 = 1.45, 0.519
        if Re < 10:
            a1, a2 = 1.4, -0.667
        elif Re < 100:
            a1, a2 = 1.36, -0.657
        elif Re < 1000:
            a1, a2 = 0.593, -0.477
        else:
            a1, a2 = 0.321, -0.388

    a = a3/(1+0.14*Re**a4)                                             # Eq 6.3
    j = a1*(1.33/PtD)**a*Re**a2                                        # Eq 6.1

    return Dimensionless(j)


def fFactor(Re, PtD, tubeConf):
    """
    Parameters
    ----------
    Re : float
        Reynolds number
    PtD : float
        Ratio between tube pitch and tube OD
    tubeConf : str
        Code with tube bundle configuration, square|triangular|rotSquare

    Return
    ------
    f : float
        friction factor coefficient
    """
    # Table 6.1
    if tubeConf == "rotSquare":
        b3, b4 = 6.59, 0.52
        if Re < 10:
            b1, b2 = 32, -1
        elif Re < 100:
            b1, b2 = 26.2, -0.913
        elif Re < 1000:
            b1, b2 = 3.5, -0.476
        elif Re < 10000:
            b1, b2 = 0.333, -0.136
        else:
            b1, b2 = 0.303, -0.126

    elif tubeConf == "square":
        b3, b4 = 6.3, 0.378
        if Re < 10:
            b1, b2 = 35, -1
        elif Re < 100:
            b1, b2 = 32.1, -0.963
        elif Re < 1000:
            b1, b2 = 6.09, -0.602
        elif Re < 10000:
            b1, b2 = 0.0815, 0.022
        else:
            b1, b2 = 0.391, -0.148

    else:
        b3, b4 = 7, 0.5
        if Re < 10:
            b1, b2 = 48, -1
        elif Re < 100:
            b1, b2 = 45.1, -0.973
        elif Re < 1000:
            b1, b2 = 4.570, -0.476
        elif Re < 10000:
            b1, b2 = 0.486, -0.152
        else:
            b1, b2 = 0.372, -0.123

    b = b3/(1+0.14*Re**b4)                                             # Eq 6.4
    f = b1*(1.33/PtD)**b*Re**b2                                        # Eq 6.2
    return Dimensionless(f)


class Shell_Tube(equipment):
    """Class that defines a shell and tubes heat exchanger

    Parameters
    ----------
    entrada : Corriente
        Input stream to equipment

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
        Ds: Diematro de la carcasa
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

    title = translate("equipment", "Shell and Tube Heat Exchanger")
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
        "Ds": 0.0,
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
        "A - " + translate("equipment", "Channel & Removable Cover"),
        "B - " + translate("equipment", "Bonnet"),
        "C - " + translate("equipment", "Removable Bundle"),
        "D - " + translate("equipment", "Special High Pressure Closure"),
        "N - " + translate("equipment", "Channel with Tubesheet & Removable Cover")]
    TEXT_SHELL = [
        "E - " + translate("equipment", "One Pass"),
        "F - " + translate("equipment", "Two Pass"),
        "G - " + translate("equipment", "Split Flow"),
        "H - " + translate("equipment", "Double Split Flow"),
        "J - " + translate("equipment", "Divided Flow"),
        "K - " + translate("equipment", "Kettle Reboiler"),
        "X - " + translate("equipment", "Cross Flow")]
    TEXT_REARHEAD = [
        "L - " + translate("equipment", "Fixed Tubesheet (A head)"),
        "M - " + translate("equipment", "Fixed Tubesheet (B head)"),
        "N - " + translate("equipment", "Fixed Tubesheet (N head)"),
        "P - " + translate("equipment", "Outside Packed Flt Head"),
        "S - " + translate("equipment", "Flt Head with Backing Dev"),
        "T - " + translate("equipment", "Pull Throught Flt Heat"),
        "U - " + translate("equipment", "U-Tube Bundle"),
        "W - " + translate("equipment", "Exit Sealed Flt Tubesheet")]
    TEXT_ORIENTATION = [
        translate("equipment", "Horizontal"),
        translate("equipment", "Vertical")]
    TEXT_DISTRIBUTION_TUBE = [
        translate("equipment", "Triangular")+", 30º",
        translate("equipment", "Diamond")+", 45º",
        translate("equipment", "Rotated Triangular")+", 60º",
        translate("equipment", "Square")+", 90º"]
    TEXT_BAFFLE_TYPE = [
        translate("equipment", "Single segmental"),
        translate("equipment", "Double segmental"),
        translate("equipment", "Triple segmental"),
        translate("equipment", "No tubes in window"),
        translate("equipment", "Disk & donut"),
        translate("equipment", "Rod")]
    TEXT_COST_TYPE = [
        translate("equipment", "Fixed Head"),
        translate("equipment", "Kettle Reboiler"),
        translate("equipment", "U-Tube")]
    TEXT_COST_MATERIAL = [
        translate("equipment", "Carbon Steel"),
        translate("equipment", "Stainless Steel 316"),
        translate("equipment", "Stainless Steel 304"),
        translate("equipment", "Stainless Steel 347"),
        translate("equipment", "Nickel 200"),
        translate("equipment", "Monel 400"),
        translate("equipment", "Inconel 600"),
        translate("equipment", "Incoloy 825"),
        translate("equipment", "Titanium"),
        translate("equipment", "Hastelloy")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entradaTubo"]:
            self.msg = translate("equipment", "undefined tubeside input")
            self.status = 0
            return
        if not self.kwargs["entradaCarcasa"]:
            self.msg = translate("equipment", "undefined shellside input")
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
        Nc = self.kwargs["Ds"]*(1-2*self.kwargs["baffleCut"])/Pt_
        Ncw = 0.8*self.kwargs["baffleCut"]*self.kwargs["Ds"]/Pt_

        Dotl = self.kwargs["Ds"]-2*self.kwargs["clearanceShellBundle"]
        Sm = self.kwargs["baffleSpacing"]*(self.kwargs["Ds"]-Dotl+(Dotl-self.kwargs["DeTube"])/Ptef*(self.kwargs["pitch"]-self.kwargs["DeTube"]))

        G = mo/Sm
        Re = self.kwargs["DeTube"]*G/self.kwargs["entradaCarcasa"].Liquido.mu
        Pr = self.kwargs["entradaCarcasa"].Liquido.Pr
        rho = self.kwargs["entradaCarcasa"].Liquido.rho
        cp = self.kwargs["entradaCarcasa"].Liquido.cp

        Dctl = Dotl-self.kwargs["DeTube"]
        tita_ctl = 2*arccos(self.kwargs["Ds"]*(1-2*self.kwargs["baffleCut"])/Dctl)
        Fc = 1+1./pi*(sin(tita_ctl)-tita_ctl)
        Fw = 1./2/pi*(tita_ctl-sin(tita_ctl))
        Stb = pi/8*((self.kwargs["DeTube"]+2*self.kwargs["clearanceTubeBaffle"])**2-self.kwargs["DeTube"]**2)*self.kwargs["NTube"]*(1+Fc)

        tita_ds = 2*arccos(1-2*self.kwargs["baffleCut"])
        Ssb = self.kwargs["Ds"]*self.kwargs["clearanceShellBaffle"]*(pi-0.5*tita_ds)

        Sb = self.kwargs["baffleSpacing"]*(self.kwargs["Ds"]-Dotl)
        Sw = 1./8*self.kwargs["Ds"]**2*(tita_ds-sin(tita_ds))-1./4*self.kwargs["NTube"]*Fw*pi*self.kwargs["DeTube"]**2
        Dw = 4*Sw/(pi*self.kwargs["DeTube"]*self.kwargs["NTube"]*0.5*(1-Fc)+self.kwargs["Ds"]*tita_ds)

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

        # Table 6.1
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
                b2 = -0.12
                b3 = 7.0
                b4 = 0.5

        a = a3/(1+0.14*Re**a4)                                       # Eq 6.3
        b = b3/(1+0.14*Re**b4)                                       # Eq 6.4

        j = a1*(1.33/P)**a*Re**a2                                    # Eq 6.1
        f = b1*(1.33/P)**b*Re**b2                                    # Eq 6.2
        hid = j*cp*G*fi/Pr**(2./3)
        h = unidades.HeatTransfCoef(hid*Jc*Jl*Jb*Jr*Js)

        DPideal = 2*f*Nc*G**2/g/rho/fi
        DPc = (nb-1)*DPideal*Rl*Rb
        if Re > 100:
            DPwideal = (2+0.6*Ncw)*mo**2/2/g/rho/Sm/Sw
        else:
            DPwideal = 26*nu*mo/g/(Sm*Sw)**0.5*(Ncw/P+self.kwargs["baffleCut"]*self.kwargs["Ds"]/Dw**2)+mo**2/g/rho/Sm/Sw
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



if __name__ == "__main__":
    from matplotlib import pyplot
    Re = logspace(0, 5, 10000)
    print(Re)
    PtD = (1.25, 1.33, 1.5)
    for p in PtD:
        ji = [jFactor(Rei, p, "rotSquare") for Rei in Re]
        pyplot.plot(Re, ji)
        fi = [fFactor(Rei, p, "rotSquare") for Rei in Re]
        pyplot.plot(Re, fi)
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.show()

