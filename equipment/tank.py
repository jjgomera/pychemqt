#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''



###Modulo que define los equipos de almacenamiento

from math import log, exp, pi

from lib.adimensional import Re
from lib.unidades import Density, Length, Currency, Volume
from lib.corriente import Corriente
from .parents import equipment
from tools.qt import translate


class Tank(equipment):
    """Class to model tank

    Parameters:

    Cost:
        material:
            0   -   Carbon steel
            1   -   Stainless steel 316
            2   -   Stainless steel 304
            3   -   Stainless steel 347
            4   -   Nickel
            5   -   Monel
            6   -   Inconel
            7   -   Zirconium
            8    -  Titanium
            9    -   Brick and rubber or brick and polyester-lined steel
            10  -   Rubber or lead-lined steel
            11  -   Polyester, fiberglass-reinforced
            12  -   Aluminum
            13  -   Copper
            14  -   Concrete

        head : Type of head
            0   -   Ellipsoidal
            1   -   Hemispherical
            2   -   Bumped
    """

    title="Deposito de almacenamiento"
    help=""

    kwargs = {
        "entrada": None,
        "Di": 0,
        "L": 0,

        "hasHelicalCoil": False,
        "helicalCoil": None,

        "f_install": 1.7,
        "Base_index": 0.0,
        "Current_index": 0.0,
        "material": 0,
        "head": 0}

    kwargsInput = ("entrada",)
    kwargsValue = ("Di", "L")
    kwargsMandatory = ("helicalCoil", )

    calculateCostos = ("C_adq", "C_inst")
    indiceCostos = 3

    TEXT_MATERIAL = [
        translate("equipment", "Carbon steel"),
        translate("equipment", "Stainless steel 316"),
        translate("equipment", "Stainless steel 304"),
        translate("equipment", "Stainless steel 347"),
        translate("equipment", "Nickel"),
        translate("equipment", "Monel"),
        translate("equipment", "Inconel"),
        translate("equipment", "Zirconium"),
        translate("equipment", "Titanium"),
        translate("equipment", "Brick and rubber or brick and polyester-lined steel"),
        translate("equipment", "Rubber or lead-lined steel"),
        translate("equipment", "Polyester, fiberglass-reinforced"),
        translate("equipment", "Aluminum"),
        translate("equipment", "Copper"),
        translate("equipment", "Concrete")]

    TEXT_HEAD = [
        translate("equipment", "Ellipsoidal"),
        translate("equipment", "Hemispherical"),
        translate("equipment", "Bumped"),
        translate("equipment", "Flat")]

    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = translate("equipment", "undefined stream input")
            self.status = 0
            return

        if not self.kwargs["Di"]:
            self.msg = translate("equipment", "undefined internal diameter")
            self.status = 0
            return
        if not self.kwargs["L"]:
            self.msg = translate("equipment", "undefined length")
            self.status = 0
            return

        return True

    def calculo(self):
        self.V = self.volumen()

        if self.kwargs["hasHelicalCoil"] and self.kwargs["helicalCoil"]:
            fluid = self.kwargs["entrada"]
            if fluid.x == 0:
                fluido = fluid.Liquido
            else:
                fluido = fluid.Vapor

            rho = fluido.rho
            mu = fluido.mu
            k = fluido.k
            v = fluid.Q*4/pi/self.kwargs["Di"]**2
            re = Re(D=self.kwargs["Di"], V=v, rho=rho, mu=mu)
            pr = fluido.Prandt

            f = self.kwargs["helicalCoil"].f(re)
            Nu = self.kwargs["helicalCoil"].Nu(re, pr)
            print("f: ", f)
            print("Nu: ", Nu)

    def volumen(self):
        """Calculate volume of shell of equipment"""
        V_shell = pi/4*self.kwargs["Di"]**2*self.kwargs["L"]

        if self.kwargs["head"] == 0:
            V_head = 4/3*pi/8*self.kwargs["Di"]**3
        elif self.kwargs["head"] == 1:
            V_head = 4/3*pi/8/2*self.kwargs["Di"]**3
        elif self.kwargs["head"] == 2:
            V_head = 0.215483/2*self.kwargs["Di"]**3
        else:
            V_head = 0

        return Volume(V_shell+V_head)


    def coste(self):
        """Calculate cost parameters"""
        self.material = self.kwargs["material"]
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]

        V = self.V.galUS

        Fm = [1, 2.7, 2.4, 3.0, 3.5, 3.3, 3.8, 11.0, 11.0, 2.75, 1.9, 0.32,
              2.7, 2.3, 0.55][self.kwargs["material"]]

        if V <= 21000:
            C = Fm*exp(2.631+1.3673*log(V)-0.06309*log(V)**2)
        else:
            C = Fm*exp(11.662+0.6104*log(V)-0.04536*log(V)**2)

        self.C_adq = Currency(C * CI / BI)
        self.C_inst = Currency(self.C_adq * self.kwargs["f_install"])


if __name__ == '__main__':
    tanque=Tank()
    print((tanque.C_inst))
    print((tanque.V))
