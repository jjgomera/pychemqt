#!/usr/bin/python
# -*- coding: utf-8 -*-

###Modulo que define los equipos de tratamiento de sólidos
#######################################################################
###librería de definición de equipos de tratamiento de sólidos:
###     -Tamices
###     -Molinos
#######################################################################

from lib import unidades
from lib.corriente import Corriente
from .parents import equipment

class Screen(equipment):
    """Clase que define los tamices"""
    title = "Screen"
    help = ""

    def coste(self, *args):
        self._indicesCoste(*args)

        C = 3.1*self.area.ft2**0.59*1000
        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)


class Grinder(equipment):
    """Clase que define los molinos de trituración de sólidos"""
    title = "Molino"
    help = ""

    def coste(self, *args, **kwargs):
        """
        tipo:
            0   -   Cone crusher
            1   -   Gyratory crusher
            2   -   Jaw crusher
            3   -   Hammer mill
            4   -   Ball mill
            5   -   Pulverizer
        """
        self._indicesCoste(*args)
        self.tipo = kwargs.get("tipo", 0)
        W = self.entrada.caudal.Tonh

        if self.tipo == 0:
            C = 1.55 * W**1.05 * 1000
        elif self.tipo == 1:
            C = 8. * W**0.6 * 1000
        elif self.tipo == 2:
            C = 6.3 * W**0.57 * 1000
        elif self.tipo == 3:
            C = 2.44 * W**0.78 * 1000
        elif self.tipo == 4:
            C = 50 * W**0.69 * 1000
        else:
            C = 22.6 * W**0.39 * 1000

        self.C_adq = unidades.Currency(C*self.Current_index/self.Base_index)
        self.C_inst = unidades.Currency(self.C_adq*self.f_install)
