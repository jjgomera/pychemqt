#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from math import exp

from lib.EoS.Cubic.PR import PR


class PRGasem(PR):
    r"""Peng-Robinson cubic equation of state with a modified dependence of
    temperature by Gassem [1]_

    .. math::
        \alpha = \exp\left(\left(A+BT_r\right)
        \left(1-T_r^{C+D\omega+E\omega^2}\right)\right)

    where A, B, C, D and E are correlation parameters generalized
    """
    __title__ = "PR Gasem (2001)"
    __status__ = "PRGasem"
    __doi__ = {
        "autor": "Gasem, K.A.M., Gao, W., Pan, R.L., Robinson Jr, R.L.",
        "title": "A modified temperature dependence for the Peng-Robinson "
                 "equation of state",
        "ref": "Fluid Phase Equilibria 181 (2001) 113-125",
        "doi": "10.1016/s0378-3812(01)00488-5"},

    def _alfa(self, cmp, T):
        """Modified correlation"""
        Tr = T/cmp.Tc

        # Parameters, Table 5
        if cmp.id == 1:
            # Special case for hydrogen
            A = -4.3
            B = 10.4
        else:
            A = 2
            B = 0.836

        C = 0.134
        D = 0.508
        E = -0.0467

        # Using Case 3 model from table 2
        alfa = exp((A+B*Tr)*(1-Tr**(C+D*cmp.f_acent+E*cmp.f_acent**2)))
        return 0, alfa


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRGasem(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRGasem(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
