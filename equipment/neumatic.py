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
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Neumatic conveying is a way that bulk solids can be pumped from one location to
another suspended in a flowing gas. Pneumatic conveying of bulk solids is
widely used in the process industry, such as food, pharmaceuticals, fine
chemikcal, paper pulp, fertilizer, cement or metal ores.

Its advantages includes:
    * Enclosure of bulk product avoiding mechanical incident, fire or
    explosion risk. Avoid too contamination of material
    * Saving space as the conveying line can be above ground.
    * Automation fairly easy

Its disadvantages includes:
    * Abrasive materials can be erosive pipes and other equipment coming into
    contact with the product
    * Particles must be dry
    * Friable materials have particle breakage
    * Material that oxidize easily need a inert gas as carrier fluid

Any pneumatic conveying systems have this basic components:
    * The gas mover, equipment like a blowers, a centrifugal fans or a
    reciprocating compressors
    * The transport pipeline
    * Product-Gas separador like a cyclone, filter or settling chamber
'''


from math import pi

from equipment.parents import equipment
from lib.drag import terminalVelocity
from lib.friction import f_friccion
from lib.neumatic import V_saltation, f_solid
from lib.unidades import Pressure, Speed
from tools.qt import translate


class Neumatic(equipment):
    """Class to define the pneumatic conveying equipment

    Parameters:
        entrada : Corriente instance to define the input stream to equipment
        D : Pipe diameter
        orientation : Geometric orientation of pipe
            0 - Horizontal
            1 - Vertical
        saltation : Method to calculate saltation velocity of gas:
            0 - Kizk (1973)
            1 - Matsumoto (1977)
            2 - Matsumoto (1975)
            3 - Matsumoto (1974)
            4 - Ochi (1991)
        fs : Method to calculate solid friction factor
            0 - Stemerding (1962)
            1 - Reddy-Pei (1969)
            2 - Swaaij-Buurman-Breugel (1970)
            3 - Capes-Nakamura (1973)
            4 - Konno-Saito (1969)
            5 - Yang (1978)
        eD : Pipe roughness, optional for Ochi saltation velocity method
        L : Pipe length
        K : Equivalent length of pipe fittings

    """
    title = translate("equipment", "Penumatic conveying")
    help = ""
    kwargs = {
        "entrada": None,
        "orientation": 0,
        "D": None,
        "eD": None,
        "L": None,
        "K": 0,
        "saltation": 0,
        "fs": 0}

    kwargsInput = ("entrada", )
    kwargsValue = ("D", "L", "K", "eD")
    kwargsList = ("saltation", "orientation", "fs")
    calculateValue = ("SLR", "DR", "VF", "vs", "V", "Re", "f", "DeltaP")

    TEXT_ORIENTATION = ("Horizontal", "Vertical")
    TEXT_SALTATION = ("Kizk", "Matsumoto 1977", "Matsumoto 1975",
                      "Matsumoto 1974", "Ochi")
    TEXT_FS = ("Stemerding", "Reddy-Pei", "Swaaij-Buurman-Breugel",
               "Capes-Nakamura", "Konno-Saito", "Yang")

    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg = translate("equipment", "undefined input")
            self.status = 0
            return

        if not self.kwargs["D"]:
            self.msg = translate("equipment", "undefined pipe diameter")
            self.status = 0
            return

        if not self.kwargs["L"]:
            self.msg = translate("equipment", "undefined pipe length")
            self.status = 0
            return

        self.msg = ""
        self.status = 1
        return True

    def calculo(self):

        entrada = self.kwargs["entrada"]

        # Solid loading ratio (SLR), ratio of the mass flow rate of solids to
        # the mass flow rate of gas
        self.SLR = entrada.solido.caudal/entrada.Gas.caudalmasico

        # Density ratio (DR), ratio of solids density to gas density, normally
        # in the order 500-1000
        self.DR = entrada.solido.rho/entrada.Gas.rho

        # Volume fraction (VF), the ratio of the volume occupied by solids to
        # the volume occupied by the gas
        self.VF = entrada.solido.Q/entrada.Q

        # Saltation velocity calculation
        # Set argument set as each method needs
        # All methods needs massflow, Dp, rhog and pipe diameter
        args = [entrada.solido.caudal,
                entrada.solido.diametro_medio,
                entrada.Gas.rho,
                self.kwargs["D"]]

        vt = terminalVelocity(entrada.solido.diametro_medio,
                              entrada.solido.rho,
                              entrada.Gas.rho, entrada.Gas.mu)

        area = pi/4*self.kwargs["D"]**2
        self.V = Speed(entrada.Gas.Q/area)
        self.Re = self.V*self.kwargs["D"]/entrada.Gas.nu*entrada.Gas.rho
        self.f = f_friccion(self.Re, eD=self.kwargs["eD"])

        # Matsumoto methods needs too solid density and terminal velocity
        if self.kwargs["saltation"] in (1, 2, 3):
            args.append(entrada.solido.rho)
            args.append(vt)

        # Ochi method need too terminal velocity and darcy friction factor
        if self.kwargs["saltation"] == 4:
            args.append(vt)
            args.append(self.f)

        self.vs = V_saltation(self.kwargs["saltation"], *args)

        # Solid friction factor calculation
        if self.kwargs["orientation"] == 1:
            eps = 0
            self.fs = fs_Yang_Horizontal(eps, self.V, self.kwargs["D"])

        args = []
        if self.kwargs["fs"] != 0:
            args.append(self.vs)
        if self.kwargs["fs"] == 4:
            args.append(self.kwargs["D"])
        if self.kwargs["fs"] == 5:
            args.append(0)  # eps
            args.append(vt)
            args.append(self.V)
        self.fs = f_solid(self.kwargs["fs"], *args)

        DeltaP_ac = self.kwargs["K"]*self.V**2/2*entrada.Gas.rho
        DeltaP_f = self.kwargs["L"]*self.V**2/self.kwargs["D"]*self.f*entrada.Gas.rho/2
        self.DeltaP = Pressure(DeltaP_ac + DeltaP_f)

        self.salida = [entrada.clone(P=entrada.P-self.DeltaP)]


if __name__ == "__main__":
    from lib.solids import Solid
    from lib.corriente import Corriente
    dm = [17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6,
          60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    sol = Solid(caudalSolido=[0.1], distribucion_diametro=dm,
                distribucion_fraccion=fracciones, solids=[638])
    kw = {"ids": [475], "fraccionMolar": [1.], "MEoS": True}
    entrada = Corriente(T=300, P=1e5, caudalMasico=1, solido=sol, **kw)
    neumatic = Neumatic(entrada=entrada, D=0.25, eD=0, saltation=4, L=50)
    print(neumatic.status, neumatic.msg)
