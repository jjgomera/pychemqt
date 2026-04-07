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


from lib.EoS.Cubic.PR78 import PR78


# Group interaction parameters, from Table 1 in [3]_, values in MPa
Aij = [[74.81, 261.5, 396.7, 32.94, 8.579, 90.25, 62.8, 40.38, 98.48],
       [51.47, 88.53, 36.72, 31.23, 29.78, 3.775, 12.78, -54.9],
       [-305.7, 145.2, 174.3, 103.3, 6.177, 101.9, -226.5],
       [263.9, 333.2, 158.9, 79.61, 177.1, 17.84],
       [13.04, 67.26, 139.3, 36.37, 40.15],
       [41.18, -3.088, 8.579, 10.29],
       [-13.38, 29.17, -26.42],
       [34.31, -105.7],
       [-50.1]]
Bij = [[165.7, 388.8, 804.3, -35, -29.51, 146.1, 41.86, 95.9, 231.6],
       [79.61, 315, 108.4, 84.76, 58.17, 144.8, 28.37, -319.5],
       [-250.8, 301.6, 352.1, 191.8, -33.97, -90.93, -51.47],
       [531.5, 203.8, 613.2, -326.0, 601.9, -109.5],
       [6.863, 167.5, 464.3, 26.42, 255.3],
       [50.79, 13.04, 76.86, -52.84],
       [20.25, 69.32, -789.2],
       [95.39, -286.5],
       [-891.1]]


# Map unifac group in compound database to group used in the PPR78 EoS
MAP_UNIFAC = {}


class PPR78(PR78):
    r"""Predictive Peng-Robinson cubic equation of state

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.37464 + 1.54226\omega-0.26992\omega^2 if \omega < 0.491\\
        m = 0.379642 + 1.48503\omega - 0.164423*\omega^2 + 0.016666*\omega^3\\
        k_{ij} = \frac{-\frac{1}{2}\sum_{k=1}^{Ng} \sum_{l=1}^{Ng}
        \left(\alpha_{ik}-\alpha_{jk}\right)\left(\alpha_{il}-\alpha_{jl}\right)
        A_{kl} \left(\frac{298.15}{T}\right)^{\left(\frac{B_{kl}}{A_{kl}}-1\right)}
        -\left(\frac{\sqrt{a_i}}{b_i}-\frac{\sqrt{a_j}}{b_j}\right)^2}
        {2\frac{\sqrt{a_i a_j}}{b_i b_j}}
        \end{array}

    Examples
    --------
    kij between propane and n-butane at T=303.15 from Appendix A of [1]_, tiny
    variation because differences in critical properties of components

    >>> mix = Mezcla(5, ids=[4, 6], caudalMolar=1, fraccionMolar=[0.5, 0.5])
    >>> eq = PPR78(303.15, 101325, mix)
    >>> print("%0.4f" % self.kij[0][1])
    0.0029
    """

    __title__ = "Predictive Peng-Robinson (1978)"
    __status__ = "PPR78"
    __doi__ = (
      {"autor": "Jaubert, J.-N., Mutelet, F.",
       "title": "VLE predictions with the Peng–Robinson equation of state and "
                "temperature dependent kij calculated through a group "
                "contribution method",
       "ref": "Fluid Phase Equilibria 224 (2004) 285-304",
       "doi": "10.1016/j.fluid.2004.06.059"},
      {"autor": "Jaubert, J.-N., Vitu, S., Mutelet, F., Corriou, J.-P.",
       "title": "Extension of the PPR78 model (predictive 1978, "
                "Peng-Robinson EOS with temperature dependent kij calculated "
                "through a group contribution method) to systems containing "
                "aromatic compounds",
       "ref": "Fluid Phase Equilibria 237 (2005) 193-211",
       "doi": "10.1016/j.fluid.2005.09.003"},
      {"autor": "Vitu, S., Jaubert, J.-N., Mutelet, F.",
       "title": "Extension of the PPR78 model (predictive 1978, "
                "Peng-Robinson EOS with temperature dependent kij calculated "
                "through a group contribution method) to systems containing "
                "naphtenic compounds",
       "ref": "Fluid Phase Equilibria 243 (2006) 9-28",
       "doi": "10.1016/j.fluid.2006.02.004"})

    def _Kij(self, eq=None):
        """Calculate binary interaction parameters"""
        group = self.mapUNIFAC()

        alfa = []
        for i, cmp in enumerate(self.componente):
            alfai = []
            for gi in group[i]:
                total = sum(group[i].values())
                alfai.append(group[i][gi]/total)
            alfa.append(alfai)

        # Eq 5 in [1]_
        kij = []
        for i, ci in enumerate(self.componente):
            kiji = []
            for j, cj in enumerate(self.componente):
                suma = 0
                for alfaik, alfajk, Aki, Bki in zip(alfa[i], alfa[j], Aij, Bij):
                    for alfail, alfajl, Akl, Bkl in zip(alfa[i], alfa[j], Aki, Bki):
                        suma += (alfaik-alfajk)*(alfail-alfajl)*Akl*1e6 \
                            * (298.15/self.T)**(Bkl/Akl-1)
                k = -0.5*suma
                k -= (self.ai[i]**0.5/self.bi[i]-self.ai[j]**0.5/self.bi[j])**2
                k /= 2*(self.ai[i]*self.ai[j])**0.5/self.bi[i]/self.bi[j]

                # Applying recomentation of unphysical result for kij > 1
                # Polishuk, I.
                # Comments on “VLE predictions with the Peng–Robinson equation
                # of state and temperature dependent kij calculated through a
                # group contribution method” by J.-N. Jaubert and F. Mutelet
                # [Fluid Phase Equilibria, 224 (2004) 285–304]
                # Fluid Phase Equilibria 249(1-2) (2006) 198-199
                # doi: 10.1016/j.fluid.2006.09.002
                if k > 1:
                    k = 1

                kiji.append(k)
            kij.append(kiji)
        return kij

    def mapUNIFAC(self):
        """Convert UNIFAC group saved in compound database to group
        contribution as used in PPR78"""
        group = []
        for cmp in self.componente:
            gi = {}
            gi[1] = cmp.UNIFAC.get(1, 0)
            gi[2] = cmp.UNIFAC.get(2, 0)
            gi[3] = cmp.UNIFAC.get(3, 0)
            gi[4] = cmp.UNIFAC.get(4, 0)
            # Special case for methane
            if cmp.id == 2:
                gi[5] = 1
            else:
                gi[5] = 0

            # Special case for ethane
            if cmp.id == 3:
                gi[6] = 1
            else:
                gi[6] = 0

            gi[7] = cmp.UNIFAC.get(10, 0)
            gi[8] = sum([cmp.UNIFAC.get(i, 0) for i in range(11, 15)])
            group.append(gi)

        return group


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    # mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    mix = Mezcla(5, ids=[4, 6, 40, 41], caudalMolar=1, fraccionMolar=[0.4, 0.4, 0.1, 0.1])
    # mix = Mezcla(5, ids=[4, 6], caudalMolar=1, fraccionMolar=[0.5, 0.5])
    eq = PPR78(303.15, 101325, mix)
    print(eq.kij)
    # print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    # eq = PR78(300, 42.477e5, mix)
    # print('%0.1f' % (eq.Vl.ccmol))
