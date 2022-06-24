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
'''


from math import exp, log, log10, tanh

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Barati, R., Neyshabouri, S.A.A.S, Ahmadi, G.",
         "title": "Development of Empirical Models with High Accuracy for "
                  "Estimation of Drag Coefficient of Flow around a Smooth "
                  "Sphere: An Evolutionary Approach",
         "ref": "Powder Technology 257 (2014) 11-19",
         "doi": "10.1016/j.powtec.2014.02.045"},
    2:
        {"autor": "Clift, R., Grace, J.R., Weber, M.E.",
         "title": "Bubbles, Drops, and Particles",
         "ref": "Academic Press, 1978.",
         "doi": ""},
    3:
        {"autor": "Ceylan, K., Altunbaş, A., Kelbaliyev, G.",
         "title": "A New Model for Estimation of Drag Force in the Flow of "
                  "Newtonian Fluids around Rigid or Deformable Particles",
         "ref": "Powder Technology 119 (2001) 250-56",
         "doi": "10.1016/S0032-5910(01)00261-3"},
    4:
        {"autor": "Almedeij, J.",
         "title": "Drag Coefficient of Flow around a Sphere: Matching "
                  "Asymptotically the Wide Trend",
         "ref": "Powder Technology 186(3) (2008) 218-223",
         "doi": "10.1016/j.powtec.2007.12.006"},

    # 30:
        # {"autor": "",
         # "title": "",
         # "ref": "",
         # "doi": ""},
}


@refDoc(__doi__, [1])
def Barati(Re, extended=False):
    r'''Calculates drag coefficient of a smooth sphere using the method in
    [1]_.

    For Re < 2e5

    .. math::
        \begin{align*}
        C_d = 5.4856·10^9\tanh(4.3774\times10^{-9}/Re)
        + 0.0709\tanh(700.6574/Re) \\
        {} + 0.3894\tanh(74.1539/Re) - 0.1198\tanh(7429.0843/Re) \\
        {} + 1.7174\tanh[9.9851/(Re+2.3384)] + 0.4744
        \end{align*}

    For 2e5 <= Re < 1e6

    .. math::
        \begin{align*}
        C_d = 8·10^{-6}\left[(Re/6530)^2+\tanh(Re) - 8\ln(Re)/\ln(10)\right] \\
        {} - 0.4119\exp(-2.08x10^{43}/[Re + Re^2]^4) \\
        {} - 2.1344\exp(-{\left[\ln(Re^2 + 10.7563)/\ln(10)\right]^2
        + 9.9867}/Re) \\
        {} + 0.1357\exp(-[(Re/1620)^2 + 10370]/Re) \\
        {} - 8.5\times 10^{-3}\{2\ln[\tanh(\tanh(Re))]/\ln(10) - 2825.7162\}/Re
        \end{align*}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]
    extended : bool
        Use the extended version on all Re range of equation

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Examples
    --------
    Selected values from Table 6 in [1]_.

    >>> print("%0.2f" % Barati(0.002))
    12008.86
    >>> print("%0.2f" % Barati(0.002, extended=True))
    12034.71
    >>> print("%0.2f" % Barati(1000))
    0.47
    '''

    if not extended and Re < 2e5:
        # Eq 22
        Cd = 5.4856e9*tanh(4.3774e-9/Re) + 0.0709*tanh(700.6574/Re) \
            + 0.3894*tanh(74.1539/Re) - 0.1198*tanh(7429.0843/Re) \
            + 1.7174*tanh(9.9851/(Re+2.3384)) + 0.4744
    elif Re <= 1e6:
        # Eq 23
        Cd = 8e-6*((Re/6530)**2 + tanh(Re) - 8*log(Re)/log(10)) \
            - 0.4119*exp(-2.08e43/(Re+Re**2)**4) \
            - 2.1344*exp(-((log(Re**2 + 10.7563)/log(10))**2 + 9.9867)/Re) \
            + 0.1357*exp(-((Re/1620)**2 + 10370)/Re) \
            - 8.5e-3*(2*log(tanh(tanh(Re)))/log(10) - 2825.7162)/Re + 2.4795
    return Cd


@refDoc(__doi__, [2])
def Clift(Re):
    r'''Calculates drag coefficient of a smooth sphere using the method in
    [2]_.

    This method use different correlation for several ranges of Re, as describe
    in Table 5.2, pag 112

    .. math::
        C_d = \left\{ \begin{array}{ll}
        \frac{24}{Re} + \frac{3}{16} & \mbox{if $Re < 0.01$}\\
        \frac{24}{Re}(1 + 0.1315Re^{0.82 - 0.05w}) & \mbox{if $0.01 < Re < 20$}\\
        \frac{24}{Re}(1 + 0.1935Re^{0.6305}) & \mbox{if $20 < Re < 260$}\\
        10^{[1.6435 - 1.1242w + 0.1558w^2} & \mbox{if $260 < Re < 1500$}\\
        10^{[-2.4571 + 2.5558w - 0.9295w^2 + 0.1049w^3} &
        \mbox{if $1500 < Re < 1.2x10^4$}\\
        10^{[-1.9181 + 0.6370w - 0.0636w^2} &
        \mbox{if $1.2x10^4 < Re < 4.4x10^4$}\\
        10^{[-4.3390 + 1.5809w - 0.1546w^2} &
        \mbox{if $4.4x10^4 < Re < 3.38x10^5$}\\
        29.78 - 5.3w & \mbox{if $3.38x10^5 < Re < 4x10^5$}\\
        0.1w - 0.49 & \mbox{if $4x10^5 < Re < 10^6$}\\
        0.19w - \frac{8x10^4}{Re} & \mbox{if $10^6 < Re$}\end{array}\right.

    where :math:`w = \log_{10}{Re}`

    Parameters
    ----------
    Re : float
        Reynolds number, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Examples
    --------
    There isn´t testing values but checking a value similar to Barati
    correlation would be enough

    >>> print("%0.0f" % Clift(0.002))
    12000
    >>> print("%0.2f" % Clift(1000))
    0.47

    Notes
    -----
    This method define the total range of drag coefficiente, including above
    1e6.
    '''
    w = log10(Re)

    if Re < 0.01:
        Cd = 3/16 + 24/Re
    elif Re < 20:
        Cd = 24/Re*(1 + 0.1315*Re**(0.82 - 0.05*w))
    elif Re < 260:
        Cd = 24/Re*(1 + 0.1935*Re**(0.6305))
    elif Re < 1500:
        Cd = 10**(1.6435 - 1.1242*w + 0.1558*w**2)
    elif Re < 1.2e4:
        Cd = 10**(-2.4571 + 2.5558*w - 0.9295*w**2 + 0.1049*w**3)
    elif Re < 4.4e4:
        Cd = 10**(-1.9181 + 0.6370*w - 0.0636*w**2)
    elif Re < 3.38e5:
        Cd = 10**(-4.3390 + 1.5809*w - 0.1546*w**2)
    elif Re < 4e5:
        Cd = 29.78 - 5.3*w
    elif Re < 1e6:
        Cd = 0.1*w - 0.49
    else:
        Cd = 0.19 - 8e4/Re
    return Cd


@refDoc(__doi__, [3])
def Ceylan(Re):
    r'''Calculates drag coefficient of a smooth sphere using the method in
    [3]_.

    .. math::
        \begin{align*}
        C_d = 1 - 0.5e^{0.182} + 10.11Re^{-2/3}e^{0.952Re^{-1/4}}
        - 0.03859Re^{-4/3}e^{1.30Re^{-1/2}} \\
        {} + 0.037\times10^{-4}Re e^{-0.125\times10^{-4}Re}
        -0.116\times10^{-10}Re^2 e^{-0.444\times10^{-5}Re}
        \end{align*}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Examples
    --------
    Selected values from Table 2, pag 253

    # >>> print("%0.0f" % Ceylan(0.1))
    # 238
    # >>> print("%0.2f" % Ceylan(0.5))
    # 49.50
    # >>> print("%0.2f" % Ceylan(1e3))
    # 0.46

    >>> print("%0.2f" % Ceylan(1e6))
    0.26

    This correlation dont return the expected vaues in paper, possible a typo
    in paper
    '''

    # Eq 15
    K1 = 1 - 0.5*exp(0.182) + 10.11*Re**(-2/3)*exp(0.952*Re**-0.25) \
        - 0.03859*Re**(-4/3)*exp(1.3*Re**(-0.5))
    K2 = 0.037e-4*Re*exp(-0.125e-4*Re) - 0.116e-10*Re**2*exp(-0.444e-5*Re)

    # Eq 14
    Cd = K1 + K2

    return Cd


def Almedeij(Re):
    r'''Calculates drag coefficient of a smooth sphere using the method in
    [4]_.

    .. math::
        \begin{align*}
        C_d = \left(\frac{1}{(\phi_1 + \phi_2)^{-1} + (\phi_3)^{-1}} + \phi_4
        \right)^{0.1} \\
        {} \phi_1 = \left(\frac{24}{Re}\right)^{10} +
        \left(\frac{21}{Re^{0.67}}\right)^{10} +
        \left(\frac{4}{Re^{0.33}}\right)^{10} + 0.4^{10} \\
        {} \phi_2 = \frac{1}{\left[(0.148 Re^{0.11})^{-10} +
        (0.5)^{-10}\right]} \\
        {} \phi_3 = \left(\frac{1.57x10^8} {Re^{1.625}}\right)^{10} \\
        {} \phi_4 = \frac{1}{(6x10^{-17}Re^{2.63})^{-10} + (0.2)^{-10}}
        \end{align*}

    Parameters
    ----------
    Re : float
        Reynolds number, [-]

    Returns
    -------
    Cd : float
        Drag coefficient [-]

    Notes
    -----
    Range is Re <= 1E6.

    Examples
    --------
    There isn´t testing values but checking a value similar to Barati
    correlation would be enough

    >>> print("%0.0f" % Almedeij(0.002))
    12000
    >>> print("%0.2f" % Almedeij(1000))
    0.44
    '''
    f1 = (24/Re)**10 + (21/Re**0.67)**10 + (4/Re**0.33)**10 + 0.4**10
    f2 = 1/((0.148*Re**0.11)**-10 + 0.5**-10)
    f3 = (1.57e8*Re**-1.625)**10
    f4 = 1/((6e-17*Re**2.63)**-10 + 0.2**-10)
    Cd = (1/((f1 + f2)**-1 + f3**-1) + f4)**0.1
    return Cd



