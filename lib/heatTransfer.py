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


###############################################################################
# library for heat transfer calculation
###############################################################################


from math import exp, factorial, pi

from numpy import tanh
from numpy.lib.scimath import log, log10

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "VDI-Gesellschaft",
         "title": "VDI Heat Atlas 2nd Edition",
         "ref": "Berlin, New York. Springer 2010.",
         "doi": ""},
    2:
        {"autor": "Bergman, T.L., Lavine, A.S., Incropera, F.P., DeWitt, D.P.",
         "title": "Introduction to Heat Transfer. 6th Ed.",
         "ref": "Wiley, 2011.",
         "doi": ""},
    3:
        {"autor": "Churchill, S.W., Chu, H.H.S.",
         "title": "Correlating Equations for Laminar and Turbulent Free "
                  "Convection from a Vertical Plate",
         "ref": "Int. J. Heat Mass Transfer 18(11) (1975) 1323-29",
         "doi": "10.1016/0017-9310(75)90243-4"},
    4:
        {"autor": "Sieder, E.N., Tate, G.E.",
         "title": "Heat Transfer and Pressure Drop of Liquids in Tubes",
         "ref": "Ind. & Eng. Chemistry 28(12) (1936) 1929-1935",
         "doi": "10.1021/ie50324a027"},
    5:
        {"autor": "Triboix, A.",
         "title": "Exact and approximate formulas for cross flow heat "
                  "exchangers with unmixed fluids",
         "ref": "Int. Comm. Heat Mass Transfer 36(2) (2009)",
         "doi": "10.1016/j.icheatmasstransfer.2008.10.012"},

    6:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},

}


@refDoc(__doi__, [1, 2, 3])
def Nu_vertical_Churchill(Pr, Ra):
    r"""Calculates Nusselt number for laminar and turbulent flows near a
    vertical surface with Churchill-Chu correlation [3]_

    .. math::
        Nu^{1/2} = 0.825 + \frac{0.387 Ra^{1/6}}
        {\left(1+(0.492/Pr)^{9/16}\right)^{8/27}}

    For laminar flow (Ra ≤ 10⁹) use this correlation with a slightly better
    accuracy

    .. math::
        Nu = 0.68 + \frac{0.67 Ra^{1/4}}
        {\left(1+(0.492/Pr)^{9/16}\right)^{4/9}}

    The correlation originally developed for vertical flat surfaces can be too
    applied to vertical cylinders if this conditions is satisfied:

    .. math::
        \frac{D}{L} \ge \frac{35}{Gr_L^{1/4}}

    Parameters
    ----------
    Pr : float
        Prandtl number [-]
    Ra : float
        Rayleigh number [-]

    Returns
    -------
    Nu : float
        Nusselt number, [-]

    Notes
    -----
    This method does not exactly represent the behavior in the transition zone
    between laminar and turbulent flows (10⁸ ≤ Ra ≤ 20⁹), but its accuracy is
    enought for engineering applications in the entire range of Rayleigh
    numbers.

    Examples
    --------
    From [1]_, Example 1, pag. 667

    >>> print("%0.f" % Nu_vertical_Churchill(0.7, 9.26e8))
    120

    From [2]_, Example 9.2, pag. 574:

    >>> print("%0.f" % Nu_vertical_Churchill(0.69, 1.813E9))
    147
    """
    f = 1+(0.492/Pr)**(9/16)
    Nu = (0.825+0.387*Ra**(1/6)/f**(8/27))**2                        # Eq 9
    return Nu


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
    if boundary == 0:
        # Inner surface heated
        Fann = 0.75*a**-0.17                                            # Eq 18
    elif boundary == 1:
        # Outer surface heated
        Fann = 0.9-0.15*a**0.6                                          # Eq 19
    elif boundary == 2:
        # Both surfaces heated
        Fann = (0.75*a**-0.17+(0.9-0.15*a**0.6))/(1+a)                  # Eq 20

    Re_ = Re*((1+a**2)*log(a)+(1-a**2))/((1-a)**2*log(a))               # Eq 17
    Xann = (1.8*log10(Re_)-1.5)**-2                                     # Eq 16
    k1 = 1.07+900/Re-0.63/(1+10*Pr)                                     # Eq 15

    # Eq 14
    Nu = Xann/8*Re*Pr / (k1+12.7*(Xann/8)**0.5*(Pr**(2/3)-1)) \
        * (1+dhL**(2/3)) * Fann
    return Nu


def h_anulli_Transition(Re, Pr, a, dhL=0, boundary=0):
    """VDI Heat Atlas G2 Pag.704"""
    g = (Re-2300)/(1e4-2300)                                            # Eq 25
    Nu_lam = h_anulli_Laminar(2300, Pr, a, dhL, boundary)
    Nu_turb = h_anulli_Turbulent(1.e4, Pr, a, dhL, boundary)

    Nu = (1-g)*Nu_lam+g*Nu_turb                                         # Eq 24
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


# Heat Exchanger design methods
@refDoc(__doi__, [2, 4, 5])
def effectiveness(NTU, Cr, flux, mixed="Cmin", exact="True"):
    """Calculate heat exchanger efectiveness

    Parameters
    ----------
    NTU : float
        Number of transfer units, [-]
    Cr : float
        Heat capacity rates ratio, Cmin/Cmax, in 0-1 range, [-]
    flux : str
        The flux type of heat exchanger
    mixed : str, optional
        Mixed stream definition only necessary for CrFSMix model, Cmin or Cmax
    exact : boolean, optional
        Use the exact recursion solutión

    Notes
    -----
    Flux would be the key code of exchanger:
        * CF: Counter flow
        * PF: Parallel flow
        * CrFMix: Crossflow, both fluids mixed
        * CrFSMix: Crossflow, one fluid mixed, other unmixed
        * CrFunMix: Crossflow, both fluids unmixed
        * 1-2TEMAE: 1-2 pass shell and tube exchanger

    In case CrFunMix is possible calculate the exact recursion solution or the
    approximate correlation from Triboix [5_]

    Returns
    -------
    effectiveness : float
        Thermal effectiveness of heat exchanger, [-]

    Notes
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


def TemperatureEffectiveness(NTU, R, flujo, **kwargs):
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


def CorrectionFactor(P, R, flujo, **kwargs):
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
        NTU = NTU_fPR(P, R, flujo, **kwargs)
        if R == 1:
            f = P/NTU/(1-P)
        else:
            f = log((1-R*P)/(1-P))/NTU/(1-R)

    return f


def NTU_fPR(P, R, flujo, **kwargs):
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


def Fi(P, R, flujo, **kwargs):
    F = CorrectionFactor(P, R, flujo, **kwargs)
    if R == 1:
        Fi = F*(1-P)
    else:
        Fi = F*P*(1-R)/log((1-R*P)/(1-P))
    return Fi
