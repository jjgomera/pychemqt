#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import log10
from scipy.constants import pi, Avogadro

from lib.meos import MEoS
from lib import unidades


class Ne(MEoS):
    """Multiparameter equation of state for neon"""
    name = "neon"
    CASNumber = "7440-01-9"
    formula = "Ne"
    synonym = "R-720"
    rhoc = unidades.Density(481.914888)
    Tc = unidades.Temperature(44.4918)
    Pc = unidades.Pressure(2678.6, "kPa")
    M = 20.179  # g/mol
    Tt = unidades.Temperature(24.556)
    Tb = unidades.Temperature(27.104)
    f_acent = -0.0387
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 107

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [], "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for neon of Katti et al. (1986).",
        "__doc__": u"""Katti, R.S., Jacobsen, R.T, Stewart, R.B., Jahangiri, M. Thermodynamic properties for neon for temperatures from the triple point to 700 K at pressures up to 700 MPa. Adv. Cryo. Eng. 31 (1986), 1189 – 1197.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 700000.0, "rhomax": 90.56, 
        "Pmin": 43.464, "rhomin": 62.059, 

        "nr1": [0.3532653449e1, -0.4513954384e1, -0.1524027959, 0.2188568609e1,
                -0.744299997e1, 0.7755627402e1, -0.3122553128e1, 0.1014206899e1,
                -0.5289214086e-1, 0.1566849239, -0.222852705, -0.1410150942e-1,
                0.7036229719e-1, -0.5882048367e-1, 0.1571172741e-1,
                0.1292202769e-2, 0.7902035603e-3, -0.3794403616e-3],
        "d1": [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 6, 6, 6],
        "t1": [0.5, 0.75, 3.5, 0.5, 0.75, 1, 1.5, 2.5, 0.25, 0.5, 2.5, 1, 3, 4,
               5, 1, 5, 6],

        "nr2": [0.4652799333e-1,  0.4524001818e-1, -0.2383421991, 0.629359013e-2,
                -0.1272313644e-2, -0.175235256e-6, 0.7188419232e-2,
                -0.5403006914e-1, 0.7578222187e-1, -0.3808588254e-1, 0.6034022431e-2],
        "d2": [1, 2, 2, 2, 2, 2, 4, 8, 8, 8, 8],
        "t2": [4, 1, 5, 8, 12, 32, 10, 6, 7, 8, 9],
        "c2": [3, 2, 2, 4, 6, 6, 2, 2, 2, 2, 2],
        "gamma2": [1]*11,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1,

    _surface = {"sigma": [0.0082803, 0.0173278, -0.0027336],
                "exp": [1.25, 2.25, 3.25]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [0.9969], "expt1": [0], "expd1": [1],
                   "a2": [-0.109, 0.0708, -2.88, -1.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 43.36814,
                "Tmin": Tt, "Tmax": 700.0,
                "a1": [1., 4437.], "exp1": [0, 1.33],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 43.464,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-10.65], "exp2": [1], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.55805e1, 0.68795e-1, 0.54840e1, -0.83760e1, 0.34276e1],
        "exp": [1.0, 1.5, 2.3, 2.8, 3.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.10601e1, 0.12076e3, -0.38553e3, 0.81655e3, -0.89907e3, 0.35466e3],
        "exp": [0.33, 1.4, 1.7, 2.2, 2.6, 3.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.23338e1, -0.36834e1, -0.85368e2, 0.22769e3, -0.17290e3],
        "exp": [0.444, 0.95, 3.5, 4.1, 4.5]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "Rabinovich (1988)"}

    _viscosity = visco0,

    def _visco0(self):
        """Rabinovich, V.A., Vasserman, A.A., Nedostup, V.I. and Veksler, L.S. "Thermophysical Properties of Neon, Argon, Krypton, and Xenon," Hemisphere Publishing Corp., 1988."""
        # FIXME: Da buenos resultados, pero los resultados difierente en la tercera cifra significativa.
        a = [17.67484, -2.78751, 311498.7, -48826500, 3938774000, -1.654629e11, 2.86561e12]
        Tr = self.T/0.29944
        y = 0.68321*(a[0]+a[1]*log10(Tr)+a[2]/Tr**2+a[3]/Tr**3+a[4]/Tr**4+a[5]/Tr**5+a[6]/Tr**6)
        nt = 266.93*(self.T*self.M)**0.5/y
        om = self.rho/1673.0
        c = [1.03010, -0.99175, 2.47127, -3.11864, 1.57066]
        b = [0.48148, -1.18732, 2.80277, -5.41058, 7.04779, -3.76608]
        sigma = 0.000000000305*(sum([ci*om**i for i, ci in enumerate(c)])-sum([bi*om**i for i, bi in enumerate(b)])*log10(self.T/122.1))
        br = 2.0/3.0*pi*Avogadro*sigma**3
        brho = self.rho/self.M*1000*br
        d = [1, 0.27676, 0.014355, 2.6480, -1.9643, 0.89161]
        nd = sum([di*brho**i for i, di in enumerate(d)])
        return unidades.Viscosity(nd*nt/100, "muPas")


if __name__ == "__main__":
    import doctest
    doctest.testmod()
