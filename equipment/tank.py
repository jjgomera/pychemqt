#!/usr/bin/python
# -*- coding: utf-8 -*-

###Modulo que define los equipos de almacenamiento

from scipy import log, exp, pi
from PyQt5.QtWidgets import  QApplication

from lib.unidades import Density, Length, Currency, Volume
from lib.corriente import Corriente
from parents import equipment


class Tank(equipment):
    """Clase que define los tanques de almacenamiento"""
    title="Deposito de almacenamiento"
    help=""
    
    def __call__(self, **kwargs):
        
        self.Di=1.
        self.L=3.
        self.volumen(3)
        self.Coste(1.7, 0)
    
    
    def volumen(self, cabeza):
        """
        cabeza: tipo de cabeza del recipiente
            0   -   Ellipsoidal
            1   -   Hemispherical
            2   -   Bumped
            3   -   Flat
        """
        V_carcasa=pi/4*self.Di**2*self.L
        
        if cabeza==0:
            V_cabeza=4./3*pi/8*self.Di**3
        elif cabeza==1:
            V_cabeza=4./3*pi/8/2*self.Di**3
        elif cabeza==2:
            V_cabeza=0.215483/2*self.Di**3
        else:
            V_cabeza=0.
            
        self.Volumen=Volume(V_carcasa+V_cabeza)
        
    
    def coste(self, *args, **kwargs):
        """
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
        """
        self._indicesCoste(*args)
        
        self.material=kwargs["material"]
        
        V=self.Volumen.galUS

        Fm=[1., 2.7, 2.4, 3.0, 3.5, 3.3, 3.8, 11.0, 11.0, 2.75, 1.9, 0.32, 2.7, 2.3, 0.55][self.material]

        if V<=21000:
            C=Fm*exp(2.631+1.3673*log(V)-0.06309*log(V)**2)
        else:
            C=Fm*exp(11.662+0.6104*log(V)-0.04536*log(V)**2)
        
        self.C_adq=Currency(C*self.Current_index/self.Base_index)
        self.C_inst=Currency(self.C_adq*self.f_install)





if __name__ == '__main__':

    tanque=Tank()
    print tanque.C_inst
    print tanque.Volumen
    
    
    
#    flash.Coste(1.7, 0, 0, 3, 1, 2, 0.05, 0.05, 0)
