#!/usr/bin/python3
# -*- coding: utf-8 -*-


# Calculate parameters for ancillary equation


from matplotlib import pyplot
import numpy as np
from scipy.optimize import curve_fit

from lib import mEoS


name = "R1336mzzZ"

fluid = mEoS.__getattribute__(name)


# T = []
# Pv = []
# rhol = []
# rhog = []
# ti = np.concatenate(
    # [np.arange(fluid.Tt, fluid.Tc-10, 1),
     # np.arange(fluid.Tc-10, fluid.Tc-1, 0.1),
     # np.arange(fluid.Tc-1, fluid.Tc-0.1, 0.01),
     # np.arange(fluid.Tc-0.1, fluid.Tc+1e-5, 0.001)])

# for t in ti:
    # st = fluid(T=t, x=0.5)
    # if st.status == 1:
        # T.append(st.T)
        # Pv.append(st.P)
        # rhol.append(st.Liquido.rho)
        # rhog.append(st.Gas.rho)

# # Remove unphysical points
# # Use inverse order index to avoid possible bug with unupdated index
# for ix in [19, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]:
    # del T[ix]
    # del Pv[ix]
    # del rhol[ix]
    # del rhog[ix]

# data = {}
# data["T"] = T
# data["Pv"] = Pv
# data["rhol"] = rhol
# data["rhog"] = rhog

# Option su save data and load to speed debug
import json
# json.dump(data, open("/home/jjgomera/Descargas/data.json", "w"), indent=2)
data = json.load(open("/home/jjgomera/Descargas/data.json", "r"))

T = np.r_[data["T"]]
Pv = np.r_[data["Pv"]]
rhol = np.r_[data["rhol"]]
rhog = np.r_[data["rhog"]]

def f_Pv(t, n1, n2, n3, n4, t1, t2, t3, t4):
    tita = 1-t/fluid.Tc
    suma = n1*tita**t1 + n2*tita**t2 + n3*tita**t3 + n4*tita**t4
    return fluid.Tc/t*suma

# bounds = ([-10, -10, -10, -10, 0, 1, 2, 3], [10, 10, 10, 10, 1, 2, 3, 5])
Pv_c, pcov = curve_fit(f_Pv, T, np.log(Pv/fluid.Pc))
print(Pv_c, pcov)

def f_rhol(t, n1, n2, n3, n4, t1, t2, t3, t4):
    tita = 1-t/fluid.Tc
    suma = 1 + n1*tita**t1 + n2*tita**t2 + n3*tita**t3 + n4*tita**t4
    return fluid.rhoc*suma


bounds = ([-10, -10, -10, -10, 0, 0.5, 1, 2], [10, 10, 10, 10, 0.5, 1, 2, 5])
rhol_c, pcov2 = curve_fit(f_rhol, T, rhol, bounds=bounds)
print(rhol_c, pcov2)

def f_rhog(t, n1, n2, n3, n4, t1, t2, t3, t4):
    tita = 1-t/fluid.Tc
    suma = n1*tita**t1 + n2*tita**t2 + n3*tita**t3 + n4*tita**t4
    return suma


rhog_c, pcov3 = curve_fit(f_rhog, T, np.log(rhog/fluid.rhoc))
print(rhog_c, pcov3)

print("Pv")
print("--------")
print("n", np.exp(Pv_c.tolist()[:4]))
print("t", Pv_c.tolist()[4:])
print("")
print("rhoLiq")
print("--------")
print("n", rhol_c.tolist()[:4])
print("t", rhol_c.tolist()[4:])
print("")
print("rhoGas")
print("--------")
print("n", rhog_c.tolist()[:4])
print("t", rhog_c.tolist()[4:])
print("")

ax1 = pyplot.subplot(311)
pyplot.plot(T, Pv, ls="", marker="x", mec="red", mfc="#00000000")
pyplot.plot(T, np.exp(f_Pv(T, *Pv_c))*fluid.Pc, color="black")
pyplot.plot(T, [fluid()._Vapor_Pressure(t) for t in T], color="blue")
ax2 = pyplot.subplot(312, sharex=ax1)
pyplot.plot(T, rhol, ls="", marker="x", mec="red", mfc="#00000000")
pyplot.plot(T, f_rhol(T, *rhol_c), color="black")
pyplot.plot(T, [fluid()._Liquid_Density(t) for t in T], color="blue")
ax3 = pyplot.subplot(313, sharex=ax1)
pyplot.plot(T, rhog, ls="", marker="x", mec="red", mfc="#00000000")
pyplot.plot(T, np.exp(f_rhog(T, *rhog_c))*fluid.rhoc, color="black")
pyplot.plot(T, [fluid()._Vapor_Density(t) for t in T], color="blue")
pyplot.show()
