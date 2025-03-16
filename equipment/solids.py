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


#######################################################################
#   Library for solid treatment equipment
#   * Grinder: Simple divider equipment
#   * Tamices
#######################################################################


from numpy import log

from equipment.parents import equipment
from lib.corriente import Corriente
from lib.solids import Solid
from lib.unidades import Currency, MassFlow, Power
from tools.qt import translate


class Screen(equipment):
    """Clase que define los tamices"""
    title = "Screen"

    area = None
    C_adq = None
    C_inst = None

    def coste(self):
        C = 3.1*self.area.ft2**0.59*1000
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]

        self.C_adq = Currency(C*CI/BI)
        self.C_inst = Currency(self.C_adq*self.kwargs["f_install"])


class Grinder(equipment):
    """Class to model solid crusher grinder"""
    title = translate("equipment", "Grinder")

    # Table 21-8 from [1]
    BOND_INDEX = (
        (translate("equipment", "Andesite"), 22.13),
        (translate("equipment", "Barite"), 6.24),
        (translate("equipment", "Basalt"), 20.41),
        (translate("equipment", "Bauxite"), 9.45),
        (translate("equipment", "Cement clinker"), 13.49),
        (translate("equipment", "Cement raw material"), 10.57),
        (translate("equipment", "Chrome ore"), 9.6),
        (translate("equipment", "Clay"), 7.1),
        (translate("equipment", "Clay, calcined"), 1.43),
        (translate("equipment", "Coal"), 11.37),
        (translate("equipment", "Coke"), 20.7),
        (translate("equipment", "Coke, fluid petroleum"), 38.6),
        (translate("equipment", "Coke, petroleum"), 73.8),
        (translate("equipment", "Copper ore"), 13.13),
        (translate("equipment", "Coral"), 10.16),
        (translate("equipment", "Diorite"), 19.40),
        (translate("equipment", "Dolomite"), 11.31),
        (translate("equipment", "Emery"), 58.18),
        (translate("equipment", "Feldspar"), 11.67),
        (translate("equipment", "Ferrochrome"), 8.87),
        (translate("equipment", "Ferromanganese"), 7.77),
        (translate("equipment", "Ferrosilicon"), 12.83),
        (translate("equipment", "Flint"), 26.16),
        (translate("equipment", "Fluorspar"), 9.76),
        (translate("equipment", "Gabbro"), 18.45),
        (translate("equipment", "Galena"), 10.19),
        (translate("equipment", "Garnet"), 12.37),
        (translate("equipment", "Glass"), 3.08),
        (translate("equipment", "Gneiss"), 20.13),
        (translate("equipment", "Gold ore"), 14.83),
        (translate("equipment", "Granite"), 14.39),
        (translate("equipment", "Graphite"), 45.03),
        (translate("equipment", "Gravel"), 25.17),
        (translate("equipment", "Gypsum rock"), 8.16),
        (translate("equipment", "Ilmenite"), 13.11),
        (translate("equipment", "Iron ore"), 15.44),
        (translate("equipment", "Hematite"), 12.68),
        (translate("equipment", "Hematite—specular"), 15.4),
        (translate("equipment", "Oolitic"), 11.33),
        (translate("equipment", "Limanite"), 8.45),
        (translate("equipment", "Magnetite"), 10.21),
        (translate("equipment", "Taconite"), 14.87),
        (translate("equipment", "Kyanite"), 18.87),
        (translate("equipment", "Lead ore"), 11.40),
        (translate("equipment", "Lead-zinc ore"), 11.35),
        (translate("equipment", "Limestone"), 11.61),
        (translate("equipment", "Limestone for cement"), 10.18),
        (translate("equipment", "Manganese ore"), 12.46),
        (translate("equipment", "Magnesite, dead burned"), 16.8),
        (translate("equipment", "Mica"), 134.5),
        (translate("equipment", "Molybdenum"), 12.97),
        (translate("equipment", "Nickel ore"), 11.88),
        (translate("equipment", "Oil shale"), 18.1),
        (translate("equipment", "Phosphate fertilizer"), 13.03),
        (translate("equipment", "Phosphate rock"), 10.13),
        (translate("equipment", "Potash ore"), 8.88),
        (translate("equipment", "Potash salt"), 8.23),
        (translate("equipment", "Pumice"), 11.93),
        (translate("equipment", "Pyrite ore"), 8.9),
        (translate("equipment", "Pyrrhotite ore"), 9.57),
        (translate("equipment", "Quartzite"), 12.18),
        (translate("equipment", "Quartz"), 12.77),
        (translate("equipment", "Rutile ore"), 12.12),
        (translate("equipment", "Sandstone"), 11.53),
        (translate("equipment", "Shale"), 16.4),
        (translate("equipment", "Silica"), 13.53),
        (translate("equipment", "Silica sand"), 16.46),
        (translate("equipment", "Silicon carbide"), 26.17),
        (translate("equipment", "Silver ore"), 17.3),
        (translate("equipment", "Sinter"), 8.77),
        (translate("equipment", "Slag"), 15.76),
        (translate("equipment", "Slag, iron blast furnace"), 12.16),
        (translate("equipment", "Slate"), 13.83),
        (translate("equipment", "Sodium silicate"), 13),
        (translate("equipment", "Spodumene ore"), 13.7),
        (translate("equipment", "Syenite"), 14.9),
        (translate("equipment", "Tile"), 15.53),
        (translate("equipment", "Tin ore"), 10.81),
        (translate("equipment", "Titanium ore"), 11.88),
        (translate("equipment", "Trap rock"), 21.1),
        (translate("equipment", "Uranium ore"), 17.93),
        (translate("equipment", "Zinc ore"), 12.42))

    TEXT_TIPO_COSTOS = (
        translate("equipment", "Cone crusher"),
        translate("equipment", "Gyratory crusher"),
        translate("equipment", "Jaw crusher"),
        translate("equipment", "Hammer mill"),
        translate("equipment", "Ball mill"),
        translate("equipment", "Pulverizer"))

    calculateCostos = ("C_adq", "C_inst")

    indiceCostos = 7
    kwargs = {
        "entrada": None,
        "BondIndex": 0,
        "exponent": 1.5,
        "D80": 0.0,

        "tipoCoste": 0,
        "f_install": 1.8,
        "Base_index": 0.0,
        "Current_index": 0.0}
    kwargsInput = ("entrada", )
    kwargsValue = ("BondIndex", "exponent", "D80")
    calculateValue = ("power", "solidflow")

    power = None
    solidflow = None
    C_adq = None
    C_inst = None

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and \
                self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if self.kwargs["entrada"] and self.kwargs["BondIndex"] and \
                self.kwargs["exponent"] and self.kwargs["D80"]:
            self.status = 1
            return True

    def calculo(self):
        G = self.kwargs["entrada"].solido.caudal
        self.solidflow = MassFlow(G)
        Xf = self.kwargs["entrada"].solido.diametro_medio
        Xp = self.kwargs["D80"]

        if self.kwargs["exponent"] == 1:
            # Kick's law
            W = 100 * self.kwargs["BondIndex"] * log(Xf/Xp)*G*1e-3
        else:
            # for other options use the general equation 21.74 from [1]
            # This become for n = 2 the Rittinger's law, and for n = 1.5 the
            # Bond Law
            ex = self.kwargs["exponent"]
            W = 100 * self.kwargs["BondIndex"] * G * 1e3 * \
                (1/Xp**(ex-1) - 1/Xf**(ex-1))
        self.power = Power(W, "kW")
        sol_in = self.kwargs["entrada"].solido
        sol_out = Solid(caudalSolido=[sol_in.caudal])
        sol_out(distribucion_diametro=[self.kwargs["D80"]],
                distribucion_fraccion=[1])
        self.salida = [self.kwargs["entrada"].clone(solido=sol_out)]

    def coste(self):
        """
        tipo:
            0   -   Cone crusher
            1   -   Gyratory crusher
            2   -   Jaw crusher
            3   -   Hammer mill
            4   -   Ball mill
            5   -   Pulverizer
        """
        tipo = self.kwargs.get("tipoCoste", 0)
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]
        W = self.kwargs["entrada"].solido.caudal.Tonh

        if tipo == 0:
            C = 1.55 * W**1.05 * 1000
        elif tipo == 1:
            C = 8. * W**0.6 * 1000
        elif tipo == 2:
            C = 6.3 * W**0.57 * 1000
        elif tipo == 3:
            C = 2.44 * W**0.78 * 1000
        elif tipo == 4:
            C = 50 * W**0.69 * 1000
        else:
            C = 22.6 * W**0.39 * 1000

        self.C_adq = Currency(C*CI/BI)
        self.C_inst = Currency(self.C_adq*self.kwargs["f_install"])


if __name__ == "__main__":
    dm = [17.5, 22.4, 26.2, 31.8, 37, 42.4, 48, 54,
          60, 69, 81.3, 96.5, 109, 127]
    fracciones = [0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                  0.05, 0.03, 0.02]
    solid = Solid(caudalSolido=[1000], distribucion_diametro=dm,
                  distribucion_fraccion=fracciones, solids=[638])
    stream = Corriente(solido=solid)
    print(stream.solido)
    grinder = Grinder(entrada=Corriente(solido=solid), D80=1e-5, BondIndex=10)
    print(grinder.status)
