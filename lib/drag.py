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


from math import exp, log, tanh

from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Barati, R., Neyshabouri, S.A.A.S, Ahmadi, G.",
         "title": "Development of Empirical Models with High Accuracy for "
                  "Estimation of Drag Coefficient of Flow around a Smooth "
                  "Sphere: An Evolutionary Approach",
         "ref": "Powder Technology 257 (2014) 11-19",
         "doi": "10.1016/j.powtec.2014.02.045"},

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
        Reynolds number of the sphere, [-]
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

