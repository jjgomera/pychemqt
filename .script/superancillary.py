#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Calculate of superancillary equation for mEoS fluid

import json

from numpy import arange, concatenate, linspace
from numpy.polynomial import chebyshev

from iapws import IAPWS95, D2O
from lib import mEoS


###############################################################################
# Configure fluid
# fluid = mEoS.F2
# fluid = IAPWS95
fluid = D2O

# Numerical tolerance
tol = 1e-7
###############################################################################


lst = []

t = concatenate([
    arange(fluid.Tt, fluid.Tc-0.1, 1),
    # arange(fluid.Tc-0.1, fluid.Tc-0.01, 0.001),
    # arange(fluid.Tc-0.01, fluid.Tc-0.001, 0.0001),
    # arange(fluid.Tc-0.001, fluid.Tc-0.0001, 0.00001),
    # arange(fluid.Tc-0.0001, fluid.Tc-0.00001, 0.000001),
    # arange(fluid.Tc-0.00001, fluid.Tc-0.000001, 0.0000001),
    # (fluid.Tc,)
    ])
states = fluid.from_list("x", 0.5, "T", t)

ps = [st.P for st in states]

coefP = chebyshev.chebfit(t, ps, 8, full=True)
cP = chebyshev.Chebyshev(coefP[0], domain=[fluid.Tt, fluid.Tc])
ajustP = [{"xmin": fluid.Tt, "xmax": fluid.Tc,
          "coef": coefP[0], "res": coefP[1][0][0]}]

rhoL = [st.Liquid.rho for st in states]
coefL = chebyshev.chebfit(t, rhoL, 8, full=True)
cL = chebyshev.Chebyshev(coefL[0], domain=[fluid.Tt, fluid.Tc])
ajustL = [{"xmin": fluid.Tt, "xmax": fluid.Tc,
          "coef": coefL[0], "res": coefL[1][0][0]}]

rhoG = [st.Gas.rho for st in states]
coefG = chebyshev.chebfit(t, rhoL, 8, full=True)
cG = chebyshev.Chebyshev(coefG[0], domain=[fluid.Tt, fluid.Tc])
ajustG = [{"xmin": fluid.Tt, "xmax": fluid.Tc,
          "coef": coefG[0], "res": coefG[1][0][0]}]


while True:
    stop = True

    print(len(ajustP))
    for i, (ajP, ajL, ajG) in reversed(list(enumerate(zip(ajustP, ajustL, ajustG)))):
        tmin = ajP["xmin"]
        tmax = ajP["xmax"]
        tmed = (tmin+tmax)/2
        print(i, tmin, tmax, ajP["res"], ajL["res"], ajG["res"])

        if ajP["res"]+ajG["res"]+ajL["res"] > tol:
            stop = False

            t1 = linspace(tmin, tmed, 1000)
            states1 = fluid.from_list("x", 0.5, "T", t1)
            ps1 = [st.P for st in states1]
            rhoL1 = [st.Liquid.rho for st in states1]
            rhoG1 = [st.Gas.rho for st in states1]

            t2 = linspace(tmed, tmax, 1000)
            states2 = fluid.from_list("x", 0.5, "T", t2)
            ps2 = [st.P for st in states2]
            rhoL2 = [st.Liquid.rho for st in states2]
            rhoG2 = [st.Gas.rho for st in states2]

            cP1 = cP.fit(t1, ps1, 12, full=True)
            # cP1 = chebyshev.chebfit(t1, ps1, 12, full=True)
            ajustP[i] = {"xmin": tmin, "xmax": tmed,
                         "coef": list(cP1[0]), "res": cP1[1][0][0]}
            cL1 = cL.fit(t1, rhoL1, 12, full=True)
            ajustL[i] = {"xmin": tmin, "xmax": tmed,
                         "coef": list(cL1[0].coef), "res": cL1[1][0][0]}
            cG1 = cG.fit(t1, rhoG1, 12, full=True)
            ajustG[i] = {"xmin": tmin, "xmax": tmed,
                         "coef": list(cG1[0].coef), "res": cG1[1][0][0]}

            cP2 = cP.fit(t2, ps2, 12, full=True)
            ajustP.insert(i+1, {"xmin": tmed, "xmax": tmax,
                                "coef": list(cP2[0].coef), "res": cP2[1][0][0]})
            cL2 = cL.fit(t2, rhoL2, 12, full=True)
            ajustL.insert(i+1, {"xmin": tmed, "xmax": tmax,
                                "coef": list(cL2[0].coef), "res": cL2[1][0][0]})
            cG2 = cG.fit(t2, rhoG2, 12, full=True)
            ajustG.insert(i+1, {"xmin": tmed, "xmax": tmax,
                                "coef": list(cG2[0].coef), "res": cG2[1][0][0]})

    if stop:
        break

dat = {}
dat["fluid"] = fluid.__name__
dat["Tt"] = fluid.Tt
dat["Tc"] = fluid.Tc
dat["jexpansions_p"] = ajustP
dat["jexpansions_rhoL"] = ajustL
dat["jexpansions_rhoV"] = ajustG
json.dump(dat, open("file.json", "w"), indent=2)
