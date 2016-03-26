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
# library for heat transfer calculation
###############################################################################


from scipy import exp, log, pi, log10


# Pipe Laminar flow
def h_tubeside_laminar_Eubank_Proctor(Gz, Gr, Pr, D, L):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen laminar
    Eubank, D. C. and Proctor W. S. - Effect of natural convection on heat transfer with laminar flow in tubes, MS Thesis, Chemical Engineering Department, Massachusetts Institute of Technology, 1951."""
    # TODO: De momento se devuelve el valor de h*d/k*(mu_w/mu_a)**0.14
    return 1.8*(Gz+12.6*(Gr*Pr*D/L)**0.4)**(1./3)


def h_tubeside_laminar_VDI(Re, Pr, D, L):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen laminar
    VDI Heat Atlas G1 Pag 695"""
    Nu1 = 4.364
    Nu2 = 1.953*(Re*Pr*D/L)**(1./3)
    Nu3 = 0.924*Pr**(1./3)*(Re*D/L)**0.5
    return (Nu1**3+0.6**3+(Nu2-0.6)**3+Nu3**3)**(1./3)


def h_tubeside_laminar_Hausen(Gz):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen laminar
    Perry Capitulo 5 pag 15
    0.1<Gz<1e4"""
    return 3.66+0.19*Gz**0.8/(1+0.117*Gz**0.467)


def h_tubeside_laminar_Sieder_Tate(Gz, Gr):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen laminar
    Sieder and Tate - Heat Transfer and Pressure Drop of Liquids in Tubes, Industrial Engineering Chemistry, Vol. 28, p. 1429, 1936.
    Perry Capitulo 5 pag 15
    Gz>100"""
    return 1.86*Gz**(1./3)+0.87*(1+0.015*Gr**(1./3))


# Pipe Turbulent flow
def h_tubeside_turbulent_Sieder_Tate(Re, Pr):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento
    Sieder and Tate - Heat Transfer and Pressure Drop of Liquids in Tubes, Industrial Engineering Chemistry, Vol. 28, p. 1429, 1936.
    Re>10000
    0.7<Pr<16700
    L/D > 10"""
    return 0.027*Re**0.8*Pr**(1./3)#*(mu/mu_w)**0.14


def h_tubeside_turbulent_Colburn(Re, Pr):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento
    DeltaT pequeña
    Re>10000
    0.7<Pr<160
    L/D > 10"""
    return 0.023*Re**0.8*Pr**(1./3)


def h_tubeside_turbulent_Dittus_Boelter(Re, Pr, calentamiento):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento
    DeltaT pequeña
    Re>10000
    0.7<Pr<160
    L/D > 10"""
    if calentamiento:
        return 0.0243*Re**0.8*Pr**0.4
    else:
        return 0.0265*Re**0.8*Pr**0.3


def h_tubeside_turbulent_ESDU(Re, Pr):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento
    40000<Re<1e6
    0.3<Pr<300
    L/D>60"""
    return 0.0225*Re**0.795*Pr**0.495*exp(-0.0225*log(Pr)**2)


def h_tubeside_turbulent_Gnielinski(Re, Pr, D, L):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento y de transición
    3000<Re<5e6
    0.5<Pr<2000
    Serth - Process heat transfer_ principles and applications pag 63"""
    f = (0.782*log(Re-1.51))**-2
    return f/8*(Re-1000.)*Pr/(1+12.7*(f/8)**0.5*(Pr**(2./3)-1))*(1+(D/L)**(2./3))


def h_tubeside_turbulent_VDI(Re, Pr, filas_tubos, alineados):
    """Coeficiente de transferencia de calor por calor sensible en el interior de tubos horizontales en regimen turbulento
    Re>10
    Pr<600
    alineados: indica si los tubos estan colocados en linea
    filas_tubos: numeros de filas de tubos"""

    if alineados:
        if Re < 300:
            a = 0.742
            m = 0.431
        elif Re < 2e5:
            a = 0.211
            m = 0.651
        else:
            a = 0.116
            m = 0.7
    else:
        if Re < 300:
            a = 1.309
            m = 0.360
        elif Re < 2e5:
            a = 0.273
            m = 0.635
        else:
            a = 0.124
            m = 0.7

    F1 = (Pr/Pr_w)**0.26
    if filas_tubos > 10:
        F2 = 1
    else:
        F2 = 0.9

    return a*Re**m*Pr**0.34*F1*F2


# Double pipe
def h_anulli_Laminar(Re, Pr, a, dhL=0, boundary=0):
    """VDI Heat Atlas G2 Pag.702"""
    if boundary == 0:         #Inner surface heated
        Nu1 = 3.66+1.2*a**-0.8
        fg = 1.615*(1+0.14*a**-0.5)
    elif boundary == 1:       #Outer surface heated
        Nu1 = 3.66+1.2*a**0.5
        fg = 1.615*(1+0.14*a**(1./3))
    elif boundary == 2:       #Both surfaces heated
        Nu1 = 3.66+(4-0.102/(a+0.02))*a**0.04
        fg = 1.615*(1+0.14*a**0.1)

    Nu2 = fg*(Re*Pr*dhL)**(1./3)
    Nu3 = (2/(1+22*Pr))**(1./6)*(Re*Pr*dhL)**0.5

    return (Nu1**3+Nu2**3+Nu3**3)**(1./3.)


def h_anulli_Turbulent(Re, Pr, a, dhL=0, boundary=0):
    """VDI Heat Atlas G2 Pag.703"""
    if boundary == 0:         #Inner surface heated
        Fann = 0.75*a**-0.17
    elif boundary == 1:       #Outer surface heated
        Fann = (0.9-0.15*a**0.6)
    elif boundary == 2:       #Both surfaces heated
        Fann = (0.75*a**-0.17+(0.9-0.15*a**0.6))/(1+a)

    Re_ = Re*((1+a**2)*log(a)+(1-a**2))/((1-a)**2*log(a))
    Xann = (1.8*log10(Re_)-1.5)**-2
    k1 = 1.07+900/Re-0.63/(1+10*Pr)
    Nu = Xann/8*Re*Pr/(k1+12.7*(Xann/8)**0.5*(Pr**(2./3.)-1))*(1+dhL**(2./3.))*Fann
    return Nu


def h_anulli_Transition(Re, Pr, a, dhL=0, boundary=0):
    """VDI Heat Atlas G2 Pag.704"""
    g = (Re-2300.)/(1.e4-2300.)
    Nu_lam = h_anulli_Laminar(2300, Pr, a, dhL, boundary)
    Nu_turb = h_anulli_Turbulent(1.e4, Pr, a, dhL, boundary)

    Nu = (1-g)*Nu_lam+g*Nu_turb
    return Nu


#
# Unused
def Nu_Convection_Free_External_Horizontal_Plate(Pr, Ra):
    """Calculo del Nusselt en convección natural externa de una pared vertical"""
    f = (1+(0.492/Pr)**(9./16))**(-16./9)
    Nu = (0.825+0.387*(Ra*f)**(1./6))**2
    return Nu


def Nu_Convection_Free_External_Horizontal_Plate(Pr, Ra):
    """Calculo del Nusselt en convección natural externa de una pared horizontal"""
    f = (1+(0.322/Pr)**(11./20))**(-20./11)
    x = Ra*f
    if x < 7e4:
        Nu = 0.766*x**0.2
    else:
        Nu = 0.15*x**(1./3)
    return Nu


# Convection
def h_tube_Condensation_Akers(fluid, Di):
    """ref Pag 557 Kakac: Boiler..."""
    Ge = fluid.caudalmasico*4/pi/Di**2*((1-fluid.x)+fluid.x*(fluid.Liquido.rho/fluid.Vapor.rho)**0.5)
    Re = Di*Ge/fluid.Liquido.mu
    if Re < 5e4:
        C = 5.03
        n = 1./3
    else:
        C = 0.0265
        n = 0.8
    return C*Re**n*fluid.Liquido.Prandt**(1./3)


def h_tube_Condensation_Cavallini(fluid, Di):
    """ref Pag 557 Kakac: Boiler..."""
    Ge = fluid.caudalmasico*4/pi/Di**2*((1-fluid.x)+fluid.x*(fluid.Liquido.rho/fluid.Vapor.rho)**0.5)
    Re = Di*Ge/fluid.Liquido.mu
    return 0.05*Re**0.8*fluid.Liquido.Prandt**(1./3)


def h_tube_Condensation_Boyko(fluid, Di):
    """ref Pag 557 Kakac: Boiler..."""
    G = fluid.caudalmasico*4/pi/Di**2
    Re = Di*G/fluid.Liquido.mu
    return 0.021*Re**0.8*fluid.Liquido.Prandt**0.43*(1+fluid.x*(fluid.Liquido.rho/fluid.Vapor.rho-1))**0.5


def h_tube_Condensation_Shah(fluid, Di):
    """ref Pag 557 Kakac: Boiler..."""
    G = fluid.caudalmasico*4/pi/Di**2
    Re = Di*G/fluid.Liquido.mu
    Nul = 0.023*Re**0.8*fluid.Liquido.Prandt**0.4
    return Nul*((1-fluid.x)**0.8+3.8*fluid.x**0.76*(1-fluid.x)**0.04/fluid.Pr**0.38)


def h_tube_Condensation_Kosky(fluid, Di):
    """ref Pag 558 Kakac: Boiler..."""
    pass


def h_tube_Condensation_Traviss(fluid, Di, X):
    """ref Pag 558 Kakac: Boiler..."""
    G = fluid.caudalmasico*4/pi/Di**2
    Re = Di*G*(1-fluid.x)/fluid.Liquido.mu
    F1 = 0.15*(1/X+2.85*X**-0.476)
    if Re < 50:
        F2 = 0.707*fluid.Liquido.Prandt*Re
    elif Re < 1125:
        F2 = 5*fluid.Liquido.Prandt+5*log(1+fluid.Liquido.Prandt*(0.0964*Re**0.585-1))
    else:
        F2 = 5*fluid.Liquido.Prandt+5*log(1+5*fluid.Liquido.Prandt)+2.5*log(0.0031*Re**0.812)

    return fluid.Pr*Re**0.9*F1/F2

