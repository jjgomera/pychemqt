#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


###Modulo que define los equipos de tratamiento de sólidos
#######################################################################
###librería de definición de equipos de tratamiento de sólidos:
###     -Tamices
###     -Molinos
#######################################################################


from numpy import log

from equipment.parents import equipment
from lib.corriente import Corriente
from lib.solids import Solid
from lib.unidades import Currency, MassFlow, Power
from tools.qt import tr


class Screen(equipment):
    """Clase que define los tamices"""
    title = "Screen"
    help = ""

    def coste(self):
        C = 3.1*self.area.ft2**0.59*1000
        CI = self.kwargs["Current_index"]
        BI = self.kwargs["Base_index"]

        self.C_adq = Currency(C*CI/BI)
        self.C_inst = Currency(self.C_adq*self.kwargs["f_install"])


class Grinder(equipment):
    """Clase que define los molinos de trituración de sólidos"""
    title = "Molino"
    help = ""

    # Table 21-8 from [1]
    BOND_INDEX = (
        (tr("pychemqt", "Andesite"), 22.13),
        (tr("pychemqt", "Barite"), 6.24),
        (tr("pychemqt", "Basalt"), 20.41),
        (tr("pychemqt", "Bauxite"), 9.45),
        (tr("pychemqt", "Cement clinker"), 13.49),
        (tr("pychemqt", "Cement raw material"), 10.57),
        (tr("pychemqt", "Chrome ore"), 9.6),
        (tr("pychemqt", "Clay"), 7.1),
        (tr("pychemqt", "Clay, calcined"), 1.43),
        (tr("pychemqt", "Coal"), 11.37),
        (tr("pychemqt", "Coke"), 20.7),
        (tr("pychemqt", "Coke, fluid petroleum"), 38.6),
        (tr("pychemqt", "Coke, petroleum"), 73.8),
        (tr("pychemqt", "Copper ore"), 13.13),
        (tr("pychemqt", "Coral"), 10.16),
        (tr("pychemqt", "Diorite"), 19.40),
        (tr("pychemqt", "Dolomite"), 11.31),
        (tr("pychemqt", "Emery"), 58.18),
        (tr("pychemqt", "Feldspar"), 11.67),
        (tr("pychemqt", "Ferrochrome"), 8.87),
        (tr("pychemqt", "Ferromanganese"), 7.77),
        (tr("pychemqt", "Ferrosilicon"), 12.83),
        (tr("pychemqt", "Flint"), 26.16),
        (tr("pychemqt", "Fluorspar"), 9.76),
        (tr("pychemqt", "Gabbro"), 18.45),
        (tr("pychemqt", "Galena"), 10.19),
        (tr("pychemqt", "Garnet"), 12.37),
        (tr("pychemqt", "Glass"), 3.08),
        (tr("pychemqt", "Gneiss"), 20.13),
        (tr("pychemqt", "Gold ore"), 14.83),
        (tr("pychemqt", "Granite"), 14.39),
        (tr("pychemqt", "Graphite"), 45.03),
        (tr("pychemqt", "Gravel"), 25.17),
        (tr("pychemqt", "Gypsum rock"), 8.16),
        (tr("pychemqt", "Ilmenite"), 13.11),
        (tr("pychemqt", "Iron ore"), 15.44),
        (tr("pychemqt", "Hematite"), 12.68),
        (tr("pychemqt", "Hematite—specular"), 15.4),
        (tr("pychemqt", "Oolitic"), 11.33),
        (tr("pychemqt", "Limanite"), 8.45),
        (tr("pychemqt", "Magnetite"), 10.21),
        (tr("pychemqt", "Taconite"), 14.87),
        (tr("pychemqt", "Kyanite"), 18.87),
        (tr("pychemqt", "Lead ore"), 11.40),
        (tr("pychemqt", "Lead-zinc ore"), 11.35),
        (tr("pychemqt", "Limestone"), 11.61),
        (tr("pychemqt", "Limestone for cement"), 10.18),
        (tr("pychemqt", "Manganese ore"), 12.46),
        (tr("pychemqt", "Magnesite, dead burned"), 16.8),
        (tr("pychemqt", "Mica"), 134.5),
        (tr("pychemqt", "Molybdenum"), 12.97),
        (tr("pychemqt", "Nickel ore"), 11.88),
        (tr("pychemqt", "Oil shale"), 18.1),
        (tr("pychemqt", "Phosphate fertilizer"), 13.03),
        (tr("pychemqt", "Phosphate rock"), 10.13),
        (tr("pychemqt", "Potash ore"), 8.88),
        (tr("pychemqt", "Potash salt"), 8.23),
        (tr("pychemqt", "Pumice"), 11.93),
        (tr("pychemqt", "Pyrite ore"), 8.9),
        (tr("pychemqt", "Pyrrhotite ore"), 9.57),
        (tr("pychemqt", "Quartzite"), 12.18),
        (tr("pychemqt", "Quartz"), 12.77),
        (tr("pychemqt", "Rutile ore"), 12.12),
        (tr("pychemqt", "Sandstone"), 11.53),
        (tr("pychemqt", "Shale"), 16.4),
        (tr("pychemqt", "Silica"), 13.53),
        (tr("pychemqt", "Silica sand"), 16.46),
        (tr("pychemqt", "Silicon carbide"), 26.17),
        (tr("pychemqt", "Silver ore"), 17.3),
        (tr("pychemqt", "Sinter"), 8.77),
        (tr("pychemqt", "Slag"), 15.76),
        (tr("pychemqt", "Slag, iron blast furnace"), 12.16),
        (tr("pychemqt", "Slate"), 13.83),
        (tr("pychemqt", "Sodium silicate"), 13),
        (tr("pychemqt", "Spodumene ore"), 13.7),
        (tr("pychemqt", "Syenite"), 14.9),
        (tr("pychemqt", "Tile"), 15.53),
        (tr("pychemqt", "Tin ore"), 10.81),
        (tr("pychemqt", "Titanium ore"), 11.88),
        (tr("pychemqt", "Trap rock"), 21.1),
        (tr("pychemqt", "Uranium ore"), 17.93),
        (tr("pychemqt", "Zinc ore"), 12.42))

    TEXT_TIPO_COSTOS = (
        tr("pychemqt", "Cone crusher"),
        tr("pychemqt", "Gyratory crusher"),
        tr("pychemqt", "Jaw crusher"),
        tr("pychemqt", "Hammer mill"),
        tr("pychemqt", "Ball mill"),
        tr("pychemqt", "Pulverizer"))

    calculateCostos = ("C_adq", "C_inst")
        # self.tipo.addItem(tr("equipment", "De cono"))
        # self.tipo.addItem(tr("equipment", "Giratorio"))
        # self.tipo.addItem(tr("equipment", "Dentado"))
        # self.tipo.addItem(tr("equipment", "De martillo"))
        # self.tipo.addItem(tr("equipment", "De bolas"))

        # self.tipo.addItem(tr("equipment", "Pulverizador"))

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
        sol_out(distribucion_diametro=[self.kwargs["D80"]], distribucion_fraccion=[1])
        self.salida = [self.kwargs["entrada"].clone(solido=sol_out)]
        # self.salida = [self.kwargs["entrada"].clone(distribucion_diametro=self.kwargs["D80"])]

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
