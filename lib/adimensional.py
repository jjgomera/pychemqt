#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


This module implements physics adimensional groups

    * :func:`Ar`: Archimedes number
    * :func:`Bi`: Biot number
    * :func:`Bo`: Bond number
    * :func:`Eu`: Euler number
    * :func:`Fo`: Fourier number
    * :func:`Fr`: Froude number
    * :func:`Ga`: Galilei number
    * :func:`Gr`: Grashof number
    * :func:`Gz`: Graetz number
    * :func:`Kn`: Knudsen number
    * :func:`Le`: Lewis number
    * :func:`Ma`: Mach number
    * :func:`Nu`: Nusselt number
    * :func:`Pe`: Péclet number
    * :func:`Pr`: Prandtl number
    * :func:`Ra`: Rayleigh number
    * :func:`Re`: Reynolds number
    * :func:`Sh`: Sherwood number
    * :func:`Sc`: Schmidt number
    * :func:`St`: Stanton number
    * :func:`We`: Weber number
'''


from scipy.constants import g

from lib.unidades import Dimensionless
from lib.utilities import refDoc

__doi__ = {
    1:
        {"autor": "VDI-Gesellschaft",
         "title": "VDI Heat Atlas 2nd Edition",
         "ref": "Berlin, New York. Springer 2010.",
         "doi": ""},
    2:
        {"autor": "Maloney, J.O.",
         "title": "Perry's Chemical Engineers' Handbook 8th Edition",
         "ref": "McGraw-Hill (2008)",
         "doi": ""}}


@refDoc(__doi__, [1])
def Ar(L, rho_p, rho, mu=None, nu=None, g=g):
    r"""Calculates Archimedes number `Ar` for two phases densities, a
    geometric length `L` or any viscosity definition.

    .. math::
        Ar=\frac{gL^3\delta\rho}{\rho v^2}

    Inputs either of any of the following sets:

    * L, density of both phases (rho_p, rho) and kinematic viscosity `nu`
    * L, density of both phases (rho_p, rho) and dynamic viscosity `mu`

    Parameters
    ------------
    L : float
        Characteristic length [m]
    rho_p : float
        Density of particle phase [kg/m³]
    rho : float
        Density of bulk phase [kg/m³]
    mu : float, optional
        Dynamic viscosity, [Pa*s]
    nu : float, optional
        Kinematic viscosity, [m²/s]
    g : float, optional
        Acceleration due to gravity, [m/s²]

    Returns
    -------
    Ar : float
        Archimedes number, [-]

    Notes
    -----
    It may be seen as a ratio of weight minus bouyancy and the inertial forces
    It is often used in equations describing the motion of particles (solid
    particles, drop, bubbles) in other fluid phase.

    Examples
    --------
    >>> print("%0.0f" % Ar(5e-4, 2610, 0.6072, nu=48.09e-6))
    2278
    """
    if rho and mu:
        nu = mu/rho
    elif not nu:
        raise Exception("undefined")
    deltarho = abs(rho_p-rho)
    return Dimensionless(g*L**3*deltarho/(rho*nu**2))


@refDoc(__doi__, [1])
def Bi(h, L, k):
    r"""Calculates Biot number `Bi` for heat transfer coefficient `h`,
    geometric length `L`, and thermal conductivity `k`.

    .. math::
        Bi=\frac{hL}{k}

    Parameters
    ----------
    h : float
        Heat transfer coefficient, [W/m^2/K]
    L : float
        Characteristic length, no typical definition [m]
    k : float
        Thermal conductivity, within the object [W/m/K]

    Returns
    -------
    Bi : float
        Biot number, [-]

    Notes
    -----
    It may be seen as a ratio of two heat transfer resistences in series.

    .. math::
        Bi = \frac{\text{Surface thermal resistance}}
        {\text{Internal thermal resistance}}

    It is useful in calculations of transient heating or cooling processes
    of solid bodies in liquid or gas flows.

    Examples
    --------
    >>> Bi(60, 0.02, 0.15)
    8.0
    """
    return Dimensionless(h*L/k)


@refDoc(__doi__, [2])
def Bo(rhol, rhog, sigma, L):
    r"""Calculates Bond number, `Bo` also known as Eotvos number.

    .. math::
        Bo = \frac{g(\rho_l-\rho_g)L^2}{\sigma}

    Parameters
    ----------
    rhol : float
        Density of liquid, [kg/m³]
    rhog : float
        Density of gas, [kg/m³]
    sigma : float
        Surface tension, [N/m]
    L : float
        Characteristic length, [m]

    Returns
    -------
    Bo : float
        Bond number, [-]
    """
    return (g*(rhol-rhog)*L**2/sigma)


@refDoc(__doi__, [2])
def Eu(dP, rho, V):
    r"""Calculates Euler number or `Eu` for a fluid of velocity `V` and
    density `rho` experiencing a pressure drop `dP`.

    .. math::
        Eu = \frac{\Delta P}{\rho V^2}

    Parameters
    ----------
    dP : float
        Pressure drop experience by the fluid, [Pa]
    rho : float
        Density of the fluid, [kg/m³]
    V : float
        Velocity of fluid, [m/s]

    Returns
    -------
    Eu : float
        Euler number, [-]

    Notes
    -----
    Used in pressure drop calculations.

    .. math::
        Eu = \frac{\text{Pressure drop}}{2\cdot \text{velocity head}}
    """
    Eu = dP/(rho*V**2)
    return Eu


@refDoc(__doi__, [1])
def Fo(k, L, t):
    r"""Calculates Fourier number `Fo`.

    .. math::
        Fo = \frac{kt}{l^2}

    Parameters
    ----------
    k : float
        Thermal diffusivity [m²/s]
    L : float
        Characteristic length [m]
    t : float
        Time of cooling or heating [s]

    Returns
    -------
    Fo : float
        Fourier number, [-]

    Notes
    -----
    Can be seen as a adimensional time.
    It is commonly used in transient conduction problems.

    Examples
    --------
    >>> print("%0.1f" % Fo(7e-6, 0.01, 60))
    4.2
    """
    return Dimensionless(k*t/L**2)


@refDoc(__doi__, [1])
def Fr(V, L, g=g):
    r"""Calculates Froude number `Fr` for velocity `V` and geometric length
    `L`. If desired, gravity can be specified as well. Normally the function
    returns the result of the equation below.

    .. math::
        Fr = \frac{V}{\sqrt{gL}}

    Parameters
    ----------
    V : float
        Velocity of the particle or fluid, [m/s]
    L : float
        Characteristic length, no typical definition [m]
    g : float, optional
        Acceleration due to gravity, [m/s^2]

    Returns
    -------
    Fr : float
        Froude number, [-]

    Notes
    -----
    Can be seen as a ratio of inertial force and gravity.

    .. math::
        Fr = \frac{\text{Inertial Force}}{\text{Gravity Force}}

    Appears in problems of forced motion when gravity has same influence,
    example in free liquid surfaces or multiphase flow.

    Examples
    --------
    >>> print("%0.0f" % Fr(5, L=0.025))
    102
    """
    return Dimensionless(V**2/(L*g))


@refDoc(__doi__, [1])
def Ga(L, rho=None, mu=None, nu=None, g=g):
    r"""Calculates Galilei number `Ga`.

    .. math::
        Ar=\frac{gL^3}{v^2}

    Inputs either of any of the following sets:

    * L and kinematic viscosity `nu`
    * L, density `rho` and dynamic viscosity `mu`

    Parameters
    ----------
    L : float
        Characteristic length [m]
    rho : float, optional
        Density of bulk phase [kg/m³]
    mu : float, optional
        Dynamic viscosity, [Pa*s]
    nu : float, optional
        Kinematic viscosity, [m^2/s]
    g : float, optional
        Acceleration due to gravity, [m/s^2]

    Returns
    -------
    Ga : float
        Galilei number, [-]

    Examples
    --------
    >>> print("%0.3f" % Ga(5e-4, nu=48.09e-6))
    0.530
    """
    if rho and mu:
        nu = mu/rho
    elif not nu:
        raise Exception("undefined")
    return Dimensionless(g*L**3/nu**2)


@refDoc(__doi__, [1])
def Gr(L, beta, T1, T2=0, rho=None, mu=None, nu=None, g=g):
    r"""Calculates Grashof number or `Gr` for a fluid with the given
    properties, temperature difference, and characteristic length.

    .. math::
        Gr = = \frac{g\beta (T_s-T_\infty)L^3}{\nu^2}
        = \frac{g\beta (T_s-T_\infty)L^3\rho^2}{\mu^2}

    Inputs either of any of the following sets:

    * L, beta, T1 and T2, and density `rho` and dynamic viscosity `mu`
    * L, beta, T1 and T2, and kinematic viscosity `nu`

    Parameters
    ----------
    L : float
        Characteristic length [m]
    beta : float
        Volumetric thermal expansion coefficient [1/K]
    T1 : float
        Temperature 1, usually a film temperature [K]
    T2 : float, optional
        Temperature 2, usually a bulk temperature (or 0 if only a difference
        is provided to the function) [K]
    rho : float, optional
        Density, [kg/m^3]
    mu : float, optional
        Dynamic viscosity, [Pa*s]
    nu : float, optional
        Kinematic viscosity, [m^2/s]
    g : float, optional
        Acceleration due to gravity, [m/s^2]

    Returns
    -------
    Gr : float
        Grashof number, [-]

    Notes
    -----
    .. math::
        Gr = \frac{\text{Buoyancy forces}}{\text{Viscous forces}}

    An error is raised if none of the required input sets are provided.
    It is a relative difference of densities within one phase only (liquid
    or gaseous), which occurs because of a spatial temperature diference.
    Important to describe heat transfer in natural convection flow problems.

    Examples
    --------
    >>> print("%0.2e" % Gr(L=0.6, beta=0.0031, T1=60, T2=20, nu=1.692e-5))
    9.17e+08
    """
    if rho and mu:
        nu = mu/rho
    elif not nu:
        raise Exception("undefined")
    Gr = g*beta*abs(T2-T1)*L**3/nu**2
    return Dimensionless(Gr)


@refDoc(__doi__, [1])
def Gz(k, D, t=None, L=None, V=None):
    r"""Calculates Graetz number `Gz`.

    .. math::
        Gz = \frac{D^2}{kt}

    Parameters
    ----------
    k : float
        Thermal diffusivity [m²/s]
    D : float
        Characteristic length [m]
    t : float, optional
        Time of cooling or heating [s]
    L : float, optional
        Second characteristic length, optional time definition [m]
    V : float, optional
        Mean flow velocity [m/s]

    Returns
    -------
    Gz : float
        Graetz number, [-]

    Notes
    -----
    It is the reciprocal of Fourier number.
    It is mainly used in calculations for steady flow, in which the time (the
    residence time of the fluid in a heated or cooled portion of a channel)
    can be expressed via the length L and the mean flow velocity V

    Examples
    --------
    >>> print("%0.0f" % Gz(1.48e-7, 0.018, V=1.5, L=3))
    1095
    """
    if V and L:
        t = L/V
    elif not t:
        raise Exception("undefined")
    return Dimensionless(D**2/k/t)


@refDoc(__doi__, [2])
def Kn(path, L):
    r"""Calculates Knudsen number or `Kn`.

    .. math::
        Kn = \frac{\lambda}{L}

    Parameters
    ----------
    path : float
        Mean free path between molecular collisions, [m]
    L : float
        Characteristic length, [m]

    Returns
    -------
    Kn : float
        Knudsen number, [-]

    Notes
    -----
    Used in mass transfer calculations.

    .. math::
        Kn = \frac{\text{Mean free path length}}{\text{Characteristic length}}
    """
    Kn = path/L
    return Kn


@refDoc(__doi__, [2])
def Le(D=None, alpha=None, Cp=None, k=None, rho=None):
    r"""Calculates Lewis number or `Le`.

    .. math::
        Le = \frac{k}{\rho C_p D} = \frac{\alpha}{D}

    Inputs can be either of the following sets:

    * Diffusivity and Thermal diffusivity
    * Diffusivity, heat capacity, thermal conductivity, and density

    Parameters
    ----------
    D : float
        Diffusivity of a species, [m²/s]
    alpha : float, optional
        Thermal diffusivity, [m²/s]
    Cp : float, optional
        Heat capacity, [J/kg/K]
    k : float, optional
        Thermal conductivity, [W/m/K]
    rho : float, optional
        Density, [kg/m³]

    Returns
    -------
    Le : float
        Lewis number, [-]

    Notes
    -----
    .. math::
        Le=\frac{\text{Thermal diffusivity}}{\text{Mass diffusivity}} =
        \frac{Sc}{Pr}

    An error is raised if none of the required input sets are provided.
    """
    if k and Cp and rho:
        alpha = k/(rho*Cp)
    elif alpha:
        pass
    else:
        raise Exception("undefined")
    Le = alpha/D
    return Le


@refDoc(__doi__, [2])
def Ma(V, c):
    r"""Calculates Mach number or `Ma` for a fluid.

    .. math::
        Ma = \frac{V}{c}

    Parameters
    ----------
    V : float
        Velocity of fluid, [m/s]
    c : float
        Speed of sound in fluid, [m/s]

    Returns
    -------
    Ma : float
        Mach number, [-]

    Notes
    -----
    Used in compressible flow calculations.

    .. math::
        Ma = \frac{\text{fluid velocity}}{\text{sonic velocity}}
    """
    return V/c


@refDoc(__doi__, [1])
def Nu(alfa, L, k):
    r"""Calculates Nusselt number or `Nu`.

    .. math::
        Nu = \frac{\alpha L}{k}

    Parameters
    ----------
    k : float
        Thermal conductivity, [W/m/K]
    L : float
        Characteristic length, [m]
    alpha : float
        Heat transfer coefficient, [W/m²/K]

    Returns
    -------
    Nu : float
        Nusselt number, [-]

    Notes
    -----
    Nusselt number is a dimensionless heat transfer coeffcient.
    Be careful with the characteristic lenght definition, it deppend of system
    configuration (internal diameter in a flow channel).

    Examples
    --------
    >>> print("%0.1f" % Nu(102.3, 0.03927, 0.03181))
    126.3
    """
    return Dimensionless(alfa*L/k)


@refDoc(__doi__, [2])
def Pe(V, L, rho=None, Cp=None, k=None, alpha=None):
    r"""Calculates heat transfer Peclet number or `Pe`

    .. math::
        Pe = \frac{VL\rho C_p}{k} = \frac{LV}{\alpha}

    Inputs either of any of the following sets:

    * V, L, density `rho`, heat capcity `Cp`, and thermal conductivity `k`
    * V, L, and thermal diffusivity `alpha`

    Parameters
    ----------
    V : float
        Velocity [m/s]
    L : float
        Characteristic length [m]
    rho : float, optional
        Density, [kg/m³]
    Cp : float, optional
        Heat capacity, [J/kg/K]
    k : float, optional
        Thermal conductivity, [W/m/K]
    alpha : float, optional
        Thermal diffusivity, [m²/s]

    Returns
    -------
    Pe : float
        Peclet number, [-]

    Notes
    -----
    .. math::
        Pe = \frac{\text{Bulk heat transfer}}{\text{Conduction heat transfer}}

    An error is raised if none of the required input sets are provided.
    """
    if rho and Cp and k:
        alpha = k/(rho*Cp)
    elif not alpha:
        raise Exception("undefined")
    Pe = V*L/alpha
    return Pe


@refDoc(__doi__, [1])
def Pr(cp=None, k=None, mu=None, nu=None, rho=None, alpha=None):
    r"""Calculates Prandtl number or `Pr` for a fluid with the given
    parameters.

    .. math::
        Pr = \frac{C_p \mu}{k} = \frac{\nu}{\alpha} = \frac{C_p \rho \nu}{k}

    Inputs can be any of the following sets:

    * Heat capacity, dynamic viscosity, and thermal conductivity
    * Thermal diffusivity and kinematic viscosity
    * Heat capacity, kinematic viscosity, thermal conductivity, and density

    Parameters
    ----------
    Cp : float
        Heat capacity, [J/kg/K]
    k : float
        Thermal conductivity, [W/m/K]
    mu : float, optional
        Dynamic viscosity, [Pa*s]
    nu : float, optional
        Kinematic viscosity, [m^2/s]
    rho : float
        Density, [kg/m^3]
    alpha : float
        Thermal diffusivity, [m^2/s]

    Returns
    -------
    Pr : float
        Prandtl number, [-]

    Notes
    -----
    .. math::
        Pr=\frac{\text{kinematic viscosity}}{\text{thermal diffusivity}} = \
            \frac{\text{momendum diffusivity}}{\text{thermal diffusivity}}

    An error is raised if none of the required input sets are provided.

    Examples
    --------
    >>> print("%0.2f" % Pr(cp=1821., k=0.134, mu=43.61e-5))
    5.93
    """
    if k and cp and mu:
        Pr = cp*mu/k
    elif nu and rho and cp and k:
        Pr = nu*rho*cp/k
    elif nu and alpha:
        Pr = nu/alpha
    else:
        raise Exception("undefined")
    return Dimensionless(Pr)


@refDoc(__doi__, [1])
def Ra(Pr, Gr):
    r"""Calculates Rayleigh number or `Ra` using Prandtl number `Pr` and
    Grashof number `Gr` for a fluid with the given
    properties, temperature difference, and characteristic length used
    to calculate `Gr` and `Pr`.

    .. math::
        Ra = PrGr

    Parameters
    ----------
    Pr : float
        Prandtl number, [-]
    Gr : float
        Grashof number, [-]

    Returns
    -------
    Ra : float
        Rayleigh number, [-]

    Notes
    -----
    Used in free convection problems only.
    """
    return Dimensionless(Pr*Gr)


@refDoc(__doi__, [1])
def Re(D, V, rho=None, mu=None, nu=None):
    r"""Calculates Reynolds number or `Re` for a fluid with the given
    properties for the specified velocity and diameter.

    .. math::
        Re = {D \cdot V}{\nu} = \frac{\rho V D}{\mu}

    Inputs either of any of the following sets:

    * V, D, density `rho` and kinematic viscosity `mu`
    * V, D, and dynamic viscosity `nu`

    Parameters
    ----------
    D : float
        Diameter [m]
    V : float
        Velocity [m/s]
    rho : float, optional
        Density, [kg/m^3]
    mu : float, optional
        Dynamic viscosity, [Pa*s]

    Returns
    -------
    Re : float
        Reynolds number, [-]

    Notes
    -----
    .. math::
        Re = \frac{\text{Momentum}}{\text{Viscosity}}

    Can be seen as a ratio of inertial forces to frictional forces.
    It's the crucial criterium to define the flow mode: laminar, turbulent.

    Examples
    --------
    >>> Re(0.052, 1.05, rho=999.8, nu=1.3e-6)
    42000.0
    """
    if rho and mu:
        nu = mu/rho
    elif not nu:
        raise Exception("undefined")
    return Dimensionless(V*D/nu)


@refDoc(__doi__, [2])
def Sh(K, L, D):
    r"""Calculates Sherwood number or `Sh`.

    .. math::
        Sh = \frac{KL}{D}

    Parameters
    ----------
    K : float
        Mass transfer coefficient, [m/s]
    L : float
        Characteristic length, no typical definition [m]
    D : float
        Diffusivity of a species [m/s²]

    Returns
    -------
    Sh : float
        Sherwood number, [-]

    Notes
    -----

    .. math::
        Sh = \frac{\text{Mass transfer by convection}}
        {\text{Mass transfer by diffusion}} = \frac{K}{D/L}
    """
    Sh = K*L/D
    return Sh


@refDoc(__doi__, [2])
def Sc(D, mu=None, nu=None, rho=None):
    r"""Calculates Schmidt number or `Sc`.

    .. math::
        Sc = \frac{\mu}{D\rho} = \frac{\nu}{D}

    Inputs can be any of the following sets:

    * Diffusivity, dynamic viscosity, and density
    * Diffusivity and kinematic viscosity

    Parameters
    ----------
    D : float
        Diffusivity of a species, [m²/s]
    mu : float, optional
        Dynamic viscosity, [Pa·s]
    nu : float, optional
        Kinematic viscosity, [m²/s]
    rho : float, optional
        Density, [kg/m³]

    Returns
    -------
    Sc : float
        Schmidt number, [-]

    Notes
    -----
    .. math::
        Sc =\frac{\text{kinematic viscosity}}{\text{molecular diffusivity}}
        = \frac{\text{viscous diffusivity}}{\text{species diffusivity}}

    An error is raised if none of the required input sets are provided.
    """
    if rho and mu:
        Sc = mu/(rho*D)
    elif nu:
        Sc = nu/D
    else:
        raise Exception("undefined")
    return Sc


@refDoc(__doi__, [2])
def St(Nu=None, Pe=None, alfa=None, rho=None, cp=None, V=None):
    r"""Calculate Stanton number `St`

    .. math::
        St = \frac{Nu}{Pe} = \frac{\alpha}{\rho c_pV}

    Inputs either of any of the following sets:

    * Péclet and Nusselt number
    * V, density `rho`, heat specific `cp` and heat transfer coefficient `alfa`

    Parameters
    ----------
    Nu : float, optional
        Nusselt number [-]
    Pe : float, optional
        Péclet number [-]
    V : float, optional
        Velocity of fluid [m/s]
    alfa : float, optional
        Heat transfer coefficient [W/m²/K]
    rho : float, optional
        Density [kg/m³]
    cp : float, optional
        Constant pressure specific heat [J/kg/K]

    Returns
    -------
    St : float
        Stanton number, [-]
    """
    if Nu and Pe:
        st = Nu/Pe
    elif alfa and rho and cp and V:
        st = alfa/(rho*cp*V)
    else:
        raise Exception("undefined")
    return Dimensionless(st)


@refDoc(__doi__, [1])
def We(V, L, rho, sigma):
    r"""Calculates Weber number, `We`, for a fluid with the given density,
    surface tension, velocity, and geometric parameter (usually diameter
    of bubble).

    .. math::
        We = \frac{V^2 L\rho}{\sigma}

    Parameters
    ----------
    V : float
        Velocity of fluid, [m/s]
    L : float
        Characteristic length, typically bubble diameter [m]
    rho : float
        Density of fluid, [kg/m^3]
    sigma : float
        Surface tension, [N/m]

    Returns
    -------
    We : float
        Weber number, [-]

    Notes
    -----
    Used in bubble calculations.

    .. math::
        We = \frac{\text{inertial force}}{\text{surface tension force}}

    Examples
    --------
    >>> print("%0.2f" % We(V=11, L=0.005, rho=1.188, sigma=0.07278))
    9.88
    """
    return Dimensionless(V**2*L*rho/sigma)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
