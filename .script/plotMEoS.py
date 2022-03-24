#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Script to generate the image for mEoS compound


import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from scipy.optimize import fsolve

from lib import mEoS
from iapws import IAPWS95


###############################################################################
# Configuration section
###############################################################################

# Define standard to use in plot, IAPWS95 very slow!
fluid = mEoS.C3

# fluid = IAPWS95
# fluid._constants["Tmin"] = IAPWS95.Tt
# fluid._constants["Tmax"] = IAPWS95.Tc
# fluid._constants["Pmax"] = 1e8


Pt = fluid(T=fluid.Tt, x=0)
Pc = fluid(T=fluid.Tc*0.999, x=0)
Vt = fluid(T=fluid.Tt, x=1)
Tmin = fluid()._constants["Tmin"]
Tmax = fluid()._constants["Tmax"]
Pmax = fluid()._constants["Pmax"]*1e3
Pmin = Pt.P

smin = Pt.s.kJkgK
smax = max(Pc.s.kJkgK, Vt.s.kJkgK)
smax += abs(smax*0.1)
hmin = Pt.h.kJkg
hmax = max(Pc.h.kJkg, Vt.h.kJkg)
hmax += abs(hmax*0.1)

# smin = Pt.s
# smax = max(Pc.s, Vt.s)
# smax += abs(smax*0.1)
# hmin = Pt.h
# hmax = max(Pc.h, Vt.h)
# hmax += abs(hmax*0.1)

# Point count for line, high value get more definition but slow calculate time
points = 50

# Saturation line format
isosat_kw = {"ls": "-", "color": "black", "lw": 1}

# Melting line format
mel_kw = {"ls": ":", "color": "black", "lw": 1}

# Isoquality lines to plot
isoq = np.arange(0.1, 1, 0.1)
isoq_kw = {"ls": "--", "color": "black", "lw": 0.8}

# Ancillary saturation equations
anc_kw = {"ls": "--", "color": "green", "lw": 0, "marker": "*"}

# Isotherm lines to plot
isoT = np.linspace(Tmin, Tmax, 10, endpoint=True)
isoT_kw = {"ls": ":", "color": "red", "lw": 0.8}

# Isobar lines to plot
isoP = np.logspace(np.log10(Pt.P), np.log10(Pmax), 10)
isoP_kw = {"ls": ":", "color": "blue", "lw": 0.8}

# # Isoenthalpic lines to plot
isoh = np.linspace(hmin, hmax, 10)
isoh_kw = {"ls": ":", "color": "orange", "lw": 0.8}

# Isoentropic lines to plot
isos = np.linspace(smin, smax, 10)
isos_kw = {"ls": ":", "color": "brown", "lw": 0.8}

# Isochor lines to plot
isov = np.logspace(np.log10(1/1200), 8, 10)
isov_kw = {"ls": ":", "color": "darkgreen", "lw": 0.8}

# Validity region
validity_kw = {"facecolor": "#F2F2F2"}

# Ideal curves
ideal_kw = {"ls": "--", "color": "black", "lw": 0.8}
label_kw = {"size": "small", "ha": "center", "va": "bottom"}
###############################################################################


# Define plot
fig = plt.figure(figsize=(15, 15))

ax1 = fig.add_subplot(3, 2, 1)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_title("Ideal curves")
ax1.set_xlabel("Tr")
ax1.set_ylabel("Pr")

ax_Ph = fig.add_subplot(3, 2, 2)
ax_Ph.set_yscale("log")
ax_Ph.set_title("P-h Diagram")
ax_Ph.set_xlabel("h, kJ/kg")
ax_Ph.set_ylabel("P, MPa")

ax_Ts = fig.add_subplot(3, 2, 3)
ax_Ts.set_title("T-s Diagram")
ax_Ts.set_xlabel("s, kJ/kgK")
ax_Ts.set_ylabel("T, K")

ax_Tv = fig.add_subplot(3, 2, 4)
ax_Tv.set_xscale("log")
ax_Tv.set_title("T-v Diagram")
ax_Tv.set_xlabel("v, m³/kg")
ax_Tv.set_ylabel("T, K")

ax_hs = fig.add_subplot(3, 2, 5)
ax_hs.set_title("h-s Diagram")
ax_hs.set_xlabel("s, kJ/kgK")
ax_hs.set_ylabel("h, kJ/kg")

ax_vu = fig.add_subplot(3, 2, 6)
ax_vu.set_yscale("log")
ax_vu.set_title("v-u Diagram")
ax_vu.set_xlabel("u, kJ/kg")
ax_vu.set_ylabel("v, m³/kg")


Ts = list(np.concatenate([
    np.linspace(fluid.Tt, 0.9*fluid.Tc, points, endpoint=False),
    np.linspace(0.9*fluid.Tc, 0.99*fluid.Tc, points, endpoint=False),
    np.linspace(0.99*fluid.Tc, fluid.Tc, points, endpoint=True)]))
Pl = list(np.concatenate([
    np.logspace(
        np.log10(Pt.P), np.log10(0.9*fluid.Pc), points, endpoint=False),
    np.linspace(0.9*fluid.Pc, 0.99*fluid.Pc, points, endpoint=False),
    np.linspace(0.99*fluid.Pc, 1.1*fluid.Pc, points, endpoint=False),
    np.logspace(
        np.log10(1.1*fluid.Pc), np.log10(Pmax), points, endpoint=True)]))
Tl = list(np.concatenate([
    np.linspace(Tmin, 0.9*fluid.Tc, points, endpoint=False),
    np.linspace(0.9*fluid.Tc, 1.1*fluid.Tc, points, endpoint=False),
    np.linspace(1.1*fluid.Tc, Tmax, points, endpoint=True)]))


# Melting pressure
if fluid is not IAPWS95 and fluid._melting:
    print("Calculating melting line...")
    Tm = np.logspace(np.log10(fluid._melting["Tmin"]),
                     np.log10(fluid._melting["Tmax"]), points*10)
    Pm = [fluid._Melting_Pressure(t) for t in Tm]
    mel = [fluid(T=t, P=p) for t, p in zip(Tm, Pm)]

    P = []
    T = []
    h = []
    s = []
    v = []
    u = []
    for p in mel:
        if p.status == 1 and p.P < Pmax:
            P.append(p.P.MPa)
            T.append(p.T)
            h.append(p.h.kJkg)
            s.append(p.s.kJkgK)
            v.append(p.v)
            u.append(p.u.kJkg)

    ax1.plot(Tm/fluid.Tc, [p/fluid.Pc for p in Pm], **mel_kw)
    ax_Ph.plot(h, P, label="Melting line", **mel_kw)
    ax_Ts.plot(s, T, label="Melting line", **mel_kw)
    ax_Tv.plot(v, T, label="Melting line", **mel_kw)
    ax_hs.plot(s, h, label="Melting line", **mel_kw)
    ax_vu.plot(u, v, label="Melting line", **mel_kw)


# Calculate saturation line
print("Calculating saturation lines...")
liq = [fluid(T=t, x=0) for t in Ts]
vap = [fluid(T=t, x=1) for t in Ts]

xliq = [line.T/fluid.Tc for line in liq]
yliq = [line.P/fluid.Pc for line in liq]
ax1.plot(xliq, yliq, **isosat_kw)
xvap = [v.T/fluid.Tc for v in vap]
yvap = [v.P/fluid.Pc for v in vap]
ax1.plot(xvap, yvap, **isosat_kw)

hliq = [line.h.kJkg for line in liq]
Pliq = [line.P.MPa for line in liq]
sliq = [line.s.kJkgK for line in liq]
vliq = [line.v for line in liq]
uliq = [line.u.kJkg for line in liq]
ax_Ph.plot(hliq, Pliq, label="Saturated Liquid", **isosat_kw)
ax_Ts.plot(sliq, Ts, label="Saturated Liquid", **isosat_kw)
ax_Tv.plot(vliq, Ts, label="Saturated Liquid", **isosat_kw)
ax_hs.plot(sliq, hliq, label="Saturated Liquid", **isosat_kw)
ax_vu.plot(uliq, vliq, label="Saturated Liquid", **isosat_kw)

hvap = [v.h.kJkg for v in vap]
Pvap = [v.P.MPa for v in vap]
svap = [v.s.kJkgK for v in vap]
vvap = [v.v for v in vap]
uvap = [v.u.kJkg for v in vap]
ax_Ph.plot(hvap, Pvap, label="Saturaded Vapor", **isosat_kw)
ax_Ts.plot(svap, Ts, label="Saturaded Vapor", **isosat_kw)
ax_Tv.plot(vvap, Ts, label="Saturaded Vapor", **isosat_kw)
ax_hs.plot(svap, hvap, label="Saturaded Vapor", **isosat_kw)
ax_vu.plot(uvap, vvap, label="Saturaded Vapor", **isosat_kw)


# Set graphics limit to only saturation region
hmax = max(hvap)
hmax += abs(hmax*0.1)
umax = max(uvap)
umax += abs(umax*0.1)
smax = max(svap)
smax += abs(smax*0.1)

ax_Ph.set_xlim(Pt.h.kJkg, hmax)
ax_Ph.set_ylim(Pt.P.MPa, Pmax/1e6)
ax_Ts.set_xlim(Pt.s.kJkgK, smax)
ax_Ts.set_ylim(Pt.T, Pc.T*1.1)
ax_Tv.set_xlim(Pt.v, Vt.v)
ax_Tv.set_ylim(Pt.T, Pc.T*1.1)
ax_hs.set_xlim(Pt.s.kJkgK, smax)
ax_hs.set_ylim(Pt.h.kJkg, hmax)
ax_vu.set_xlim(Pt.u.kJkg, umax)
ax_vu.set_ylim(Pt.v, Vt.v)


# Ancillary equation
Tanc = np.concatenate([
    np.linspace(fluid.Tt, 0.9*fluid.Tc, 10),
    np.linspace(0.9*fluid.Tc, fluid.Tc, 10, endpoint=True)])
Panc = [fluid()._Vapor_Pressure(t)/fluid.Pc for t in Tanc]
ax1.plot(Tanc/fluid.Tc, Panc, label="Ancillary Vapor Pressure", **anc_kw)

lanc = [1/fluid()._Liquid_Density(t) for t in Tanc]
vanc = [1/fluid()._Vapor_Density(t) for t in Tanc]
ax_Tv.plot(lanc, Tanc, label="Ancillary Liquid Density", **anc_kw)
ax_Tv.plot(vanc, Tanc, label="Ancillary Vapor Density", **anc_kw)

# Calculate isoquality lines
# print("Calculating isoquality lines...")
# for q in isoq:
    # txt = "x: %s" % q
    # print("    %s" % txt)
    # pts = [fluid(T=t, x=q) for t in Ts]

    # x = [p.h.kJkg for p in pts]
    # y = [p.P.MPa for p in pts]
    # ax_Ph.plot(x, y, label=txt, **isoq_kw)

    # x = [p.s.kJkgK for p in pts]
    # y = [p.T for p in pts]
    # ax_Ts.plot(x, y, label=txt, **isoq_kw)

    # x = [p.v for p in pts]
    # y = [p.T for p in pts]
    # ax_Tv.plot(x, y, label=txt, **isoq_kw)

    # x = [p.s.kJkgK for p in pts]
    # y = [p.h.kJkg for p in pts]
    # ax_hs.plot(x, y, label=txt, **isoq_kw)

    # x = [p.u.kJkg for p in pts]
    # y = [p.v for p in pts]
    # ax_vu.plot(x, y, label=txt, **isoq_kw)


# Calculate isotherm lines
# print("Calculating isotherm lines...")
# for T in isoT:
    # txt = "T: %s K" % T
    # print("    %s" % txt)
    # # Calculate the saturation point if available
    # if T < fluid.Tc:
        # sat_pnt = []
        # for x in np.linspace(1, 0.1, 10, endpoint=False):
            # sat_pnt.append(fluid(T=T, x=x))
        # for x in np.linspace(0.1, 0, 11):
            # sat_pnt.append(fluid(T=T, x=x))
        # sat = True
    # else:
        # sat = False

    # pts = []
    # melting = True

    # for p in Pl:
        # point = fluid(P=p, T=T)
        # if point.status == 1:

            # # Discard point below the melting line
            # if melting and fluid._melting:
                # pm = fluid()._Melting_Pressure(point.T)
                # if pm < point.P:
                    # point = fluid(T=T, P=pm)
                    # melting = False

            # # Add saturation point if neccesary
            # if sat and T < fluid.Tc and sat_pnt and point.s < sat_pnt[-1].s:
                # for p in sat_pnt:
                    # pts.append(p)
                # sat = False
            # pts.append(point)

            # if not melting:
                # break

    # h = []
    # P = []
    # s = []
    # u = []
    # v = []
    # for p in pts:
        # if p.status:
            # h.append(p.h.kJkg)
            # P.append(p.P.MPa)
            # s.append(p.s.kJkgK)
            # u.append(p.u.kJkg)
            # v.append(p.v)
    # ax_Ph.plot(h, P, label=txt, **isoT_kw)
    # ax_hs.plot(s, h, label=txt, **isoT_kw)
    # ax_vu.plot(u, v, label=txt, **isoT_kw)


# Calculate isobar lines
# print("Calculating isobar lines...")
# for P in isoP:
    # txt = "P: %0.5g MPa" % (P/1e6)
    # print("    %s" % txt)
    # # Calculate the saturation point if available
    # if P < fluid.Pc:
        # sat_pnt = []
        # for x in np.linspace(1, 0.1, 10):
            # sat_pnt.append(fluid(P=P, x=x))
        # for x in np.linspace(0.1, 0, 11):
            # sat_pnt.append(fluid(P=P, x=x))
        # sat = True
    # else:
        # sat = False

    # pts = []
    # for t in Tl:
        # point = fluid(P=P, T=t)

        # if point.status == 1:

            # # Discard point below the melting line
            # if fluid._melting:
                # pm = fluid()._Melting_Pressure(point.T)
                # if pm < point.P:

                    # def f(t):
                        # pm = fluid()._Melting_Pressure(t)
                        # return pm - P

                    # T = fsolve(f, point.T)
                    # point = fluid(T=T, P=P)

            # # Add saturation point if neccesary
            # if sat and P < fluid.Pc and sat_pnt and point.s > sat_pnt[-1].s:
                # for p in sat_pnt[::-1]:
                    # pts.append(p)
                # sat = False

            # pts.append(point)

    # h = []
    # T = []
    # s = []
    # u = []
    # v = []
    # for p in pts:
        # if p.status == 1:
            # h.append(p.h.kJkg)
            # T.append(p.T)
            # s.append(p.s.kJkgK)
            # u.append(p.u.kJkg)
            # v.append(p.v)
    # ax_Ts.plot(s, T, label=txt, **isoP_kw)
    # ax_Tv.plot(v, T, label=txt, **isoP_kw)
    # ax_hs.plot(s, h, label=txt, **isoP_kw)
    # ax_vu.plot(u, v, label=txt, **isoP_kw)


# Calculate isoenthalpic lines
print("Calculating isoenthalpic lines...")
for h in isoh[-5:-2]:
    txt = "h: %0.5g kJ/kg" % (h/1000)
    print("    %s" % txt)

    pts = []
    for p in Pl:
        point = fluid(P=p, h=h*1000)
        if point.status == 1:

            # Discard point below the melting line
            if fluid._melting:
                pm = fluid()._Melting_Pressure(point.T)
                if pm < point.P:

                    def f(t):
                        pm = fluid()._Melting_Pressure(t)
                        return fluid(T=t, P=pm).h.kJkg-h

                    T = fsolve(f, point.T)
                    pm = fluid()._Melting_Pressure(T)
                    point = fluid(T=T, P=pm)

            pts.append(point)
            print(point.T, point.s, point.x)
        else:
            print("Error", p, point.status, point.msg)

    T = []
    s = []
    P = []
    u = []
    v = []
    h = []
    for p in pts:
        T.append(p.T)
        s.append(p.s.kJkgK)
        h.append(p.h.kJkg)
        u.append(p.u.kJkg)
        P.append(p.P.MPa)
        v.append(p.v)
    ax_Ph.plot(h, P, label=txt, **isoh_kw)
    ax_Ts.plot(s, T, label=txt, **isoh_kw)
    ax_Tv.plot(v, T, label=txt, **isoh_kw)
    # ax_hs.plot(s, h, label=txt, **isoh_kw)
    # ax_vu.plot(u, v, label=txt, **isoh_kw)

# Calculate isoentropic lines
# print("Calculating isoentropic lines...")
# for s in isos:
    # txt = "s: %0.5g kJ/kgK" % s
    # print("    %s" % txt)

    # phase = True
    # ts = Pt.Tt

    # # Find saturated liquid point
    # if Pt.s.kJkgK > s:
        # def f(ti):
            # if ti < fluid.Tt:
                # ti = fluid.Tt
            # return fluid(T=ti, x=0).s - s*1000
        # ts = fsolve(f, (fluid.Tc-fluid.Tt)/2)
    # else:
        # phase = False

    # pts = []
    # for t in Tl:
        # point = fluid(T=t, s=s*1000)
        # if point.status == 1:

            # # Discard point over the melting line
            # if fluid._melting:
                # pm = fluid()._Melting_Pressure(point.T)
                # if pm < point.P and point.x == 1:

                    # # Fuction to find the intersection between isoentropic line
                    # # and melting line
                    # def f(ti):
                        # pm = fluid()._Melting_Pressure(ti)
                        # return fluid(T=ti, P=pm).s.kJkgK-s

                    # T = fsolve(f, point.T)
                    # pm = fluid()._Melting_Pressure(T)
                    # point = fluid(T=T, P=pm)
                    # break

            # # Add saturation point and near to that
            # if phase and point.T > ts:
                # phase = False

                # if not pts:
                    # tmin = Tmin
                # else:
                    # tmin = pts[-1].T
                # ti = np.concatenate([
                    # np.linspace(Tmin, ts, points, endpoint=False),
                    # np.linspace(ts, point.T, points*10, endpoint=False)])

                # for ts in ti:
                    # ptss = fluid(T=ts, s=s*1000)
                    # if ptss.status == 1:
                        # pts.append(ptss)

            # pts.append(point)

    # T = []
    # h = []
    # P = []
    # u = []
    # v = []
    # for p in pts:
        # T.append(p.T)
        # h.append(p.h.kJkg)
        # u.append(p.u.kJkg)
        # P.append(p.P.MPa)
        # v.append(p.v)

    # ax_Ph.plot(h, P, label=txt, **isos_kw)
    # ax_Tv.plot(v, T, label=txt, **isos_kw)
    # ax_vu.plot(u, v, label=txt, **isos_kw)

# Calculate isochor lines
# print("Calculating isochor lines...")
# for v in isov:
    # txt = "v: %0.5g m³/kg" % v
    # print("    %s" % txt)

    # pts = [fluid(T=t, v=v) for t in Tl]

    # h = []
    # T = []
    # s = []
    # u = []
    # P = []
    # for p in pts:
        # if p.status == 1:

            # # Discard point over the melting line
            # if fluid._melting:
                # pm = fluid()._Melting_Pressure(p.T)
                # if pm < p.P:

                    # # Fuction to find the intersection between isochor line
                    # def f(ti):
                        # pm = fluid()._Melting_Pressure(ti)
                        # return fluid(T=ti, P=pm).v-v

                    # t = None
                    # for to in [p.Tt, p.T, p.Tc]:
                        # rinput = fsolve(f, p.T, full_output=True)
                        # if rinput[2] == 1 and abs(rinput[1]["fvec"]) < 1e-5:
                            # t = rinput[0]
                            # break

                    # if t:
                        # pm = fluid()._Melting_Pressure(t)
                        # p = fluid(T=t, P=pm)
                    # else:
                        # continue

            # h.append(p.h.kJkg)
            # T.append(p.T)
            # s.append(p.s.kJkgK)
            # u.append(p.u.kJkg)
            # P.append(p.P.MPa)

    # ax_Ph.plot(h, P, label=txt, **isov_kw)
    # ax_Ts.plot(s, T, label=txt, **isov_kw)
    # ax_hs.plot(s, h, label=txt, **isov_kw)

# Show the validity range of equation
Prmin = Pt.P/fluid.Pc
Trmin = Tmin/fluid.Tc
Prmax = Pmax/fluid.Pc
Trmax = Tmax/fluid.Tc

if fluid._melting:
    vert = [(Trmin, Prmin)]
    for t, p in zip(Tm, Pm):
        vert.append((t/fluid.Tc, min(p/fluid.Pc, Prmax)))
else:
    vert = [(Trmin, Prmin), (Trmin, Prmax)]

vert.append((Trmax, Prmax))
vert.append((Trmax, Prmin))

poly = Polygon(vert, **validity_kw)
ax1.add_patch(poly)


# Ideal Curves
# print("Calculating ideal curves...")
# T = np.logspace(2, 4, points*5, endpoint=True)

# for curva in ["ideal", "boyle", "joule-thomson", "joule"]:
    # print(curva)
    # T_line = []
    # P_line = []
    # pmin = Pmin
    # for t in T:
        # P = fluid()._IdealCurve(curva, t)
        # if fluid._melting:
            # Pm = fluid()._Melting_Pressure(t)
        # else:
            # Pm = 1e999

        # # Discard point below the saturaion line
        # if curva in ["boyle", "joule-thomson"]:
            # if t < fluid.Tc:
                # st = fluid(T=max(t, Tmin), x=0)
                # if st.P > P:
                    # P = st.P
                    # pmin = P

        # if 0 < P <= Pm:
            # T_line.append(t/fluid.Tc)
            # P_line.append(P/fluid.Pc)

    # # Append a point in a x axis
    # T_line.append(T_line[-1])
    # P_line.append(Pt.Pr)

    # ax1.plot(T_line, P_line, label=curva, **ideal_kw)

    # # Set ideal curve label
    # if curva == "joule-thomson":
        # txt = "j-t"
        # ax1.set_ylim(bottom=0.1*pmin/fluid.Pc)
    # else:
        # txt = curva

    # # Location of text, in boyle and j-t using the maximum, in other a position
    # # at left of plot
    # if curva in ["joule", "boyle", "joule-thomson"]:
        # i = P_line.index(max(P_line))
    # else:
        # i = int(len(T_line)*0.4)

    # xmin, xmax = ax1.get_xlim()
    # ymin, ymax = ax1.get_ylim()

    # fx = (np.log10(T_line[i+1])-np.log10(T_line[i-1]))/(
        # np.log10(xmax)-np.log10(xmin))
    # fy = (np.log10(P_line[i+1])-np.log10(P_line[i-1]))/(
        # np.log10(ymax)-np.log10(ymin))
    # rot = np.arctan(fy/fx)*360/2/np.pi
    # ax1.annotate(txt, (T_line[i], P_line[i]), rotation=rot, **label_kw)


ax1.annotate("  c", (1, 1), va="center", ha="left")


legend1 = [
    Patch(label='Validity range of EoS', **validity_kw),
    Line2D([0], [0], label='Saturarion', **isosat_kw)]
if fluid._melting:
    legend1.append(Line2D([0], [0], label='Melting pressure', **mel_kw))
ax1.legend(handles=legend1, loc='upper right')

legend_tv = [
    Line2D([0], [0], label='Saturarion', **isosat_kw),
    Line2D([0], [0], label='Ancillary equations', **anc_kw),
    Line2D([0], [0], label='Isoquality lines', **isoq_kw),
    Line2D([0], [0], label='Isoterm lines', **isoT_kw),
    Line2D([0], [0], label='Isobar lines', **isoP_kw),
    Line2D([0], [0], label='Isoenthalpic lines', **isoh_kw),
    Line2D([0], [0], label='Isoentropic lines', **isos_kw),
    Line2D([0], [0], label='Isochor lines', **isov_kw)]
ax_Tv.legend(handles=legend_tv, loc='upper right')

fig.set_tight_layout(True)

plt.show()

file = "/home/jjgomera/Descargas/ImagesPychemqt/%s.png" % fluid.__name__
fig.savefig(file)

subprocess.Popen(['gpicview', file])
