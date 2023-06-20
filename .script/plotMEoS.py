#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=protected-access


"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


# Script to generate the image for mEoS compound


import argparse
import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from scipy.optimize import fsolve

from iapws import IAPWS95
from lib import mEoS


desc = 'Generate several phase plot for a compound in the mEos library'
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("cmp", type=str, nargs="?", help='Name of compound')
parser.add_argument(
    "points", type=int, nargs="?",
    help='Definition of plot, as number of point for each line')
parser.add_argument("-m", "--matplotlib", help="Show matplotlib dialog",
                    action="store_true")
parser.add_argument("-s", "--show", help="Show image", action="store_true")

args = parser.parse_args()


###############################################################################
# Configuration section
###############################################################################

# Point count for line, high value get more definition but slow calculate time
points = 10

# Define standard to use in plot, IAPWS95 very slow!
fluid = mEoS.Ar

# fluid = IAPWS95
# fluid._constants["Tmin"] = IAPWS95.Tt
# fluid._constants["Tmax"] = IAPWS95.Tc
# fluid._constants["Pmax"] = 1e8


# Use the input parameter if given
if args.cmp:
    fluid = mEoS.__getattribute__(args.cmp)
if args.points:
    points = args.points

Tc = fluid.eq[0].get("Tc", fluid.Tc)

Pt = fluid(T=max(fluid.Tt, fluid()._constants["Tmin"]), x=0)
Pc = fluid(T=Tc*0.999, x=0)
Vt = fluid(T=max(fluid.Tt, fluid()._constants["Tmin"]), x=1)
Vmin = fluid(T=fluid()._constants["Tmax"], P=Pt.P)
Tmin = fluid()._constants["Tmin"]
Tmax = fluid()._constants["Tmax"]
Pmax = fluid()._constants["Pmax"] * 1e3
Pmin = Pt.P

smin = Pt.s.kJkgK
smax = max(Pc.s.kJkgK, Vt.s.kJkgK)
smax += abs(smax * 0.1)
hmin = Pt.h.kJkg
hmax = max(Pc.h.kJkg, Vt.h.kJkg)
hmax += abs(hmax * 0.1)

# smin = Pt.s
# smax = max(Pc.s, Vt.s)
# smax += abs(smax*0.1)
# hmin = Pt.h
# hmax = max(Pc.h, Vt.h)
# hmax += abs(hmax*0.1)

# Saturation line format
isosat_kw = {"ls": "-", "color": "black", "lw": 1}

# PIP lines
PIP_kw = {"ls": "-", "color": "black", "lw": 0.5}

# Melting line format
mel_kw = {"ls": ":", "color": "black", "lw": 1}

# Isoquality lines to plot
isoq = np.arange(0.1, 1, 0.1)
isoq_kw = {"ls": "--", "color": "black", "lw": 0.8}

# Ancillary saturation equations
anc_kw = {"ls": "", "color": "green", "marker": "*"}

# Isotherm lines to plot
if Tmax > 10*Tc:
    extended = np.logspace(np.log10(Tc), np.log10(Tmax), 10)
else:
    extended = np.linspace(Tc, Tmax, 5)
#     extended = np.r_[Tc, Tmax]
isoT = np.concatenate([
    np.linspace(Tmin, Tc, 5, endpoint=False), extended])
isoT_kw = {"ls": ":", "color": "red", "lw": 0.8}

# Isobar lines to plot
isoP = np.logspace(np.log10(Pt.P), np.log10(Pmax), 10)
isoP_PIP = np.r_[0.1, 0.2, 0.5, 1, 2, 5, fluid.Pc.MPa, 15, 20, 25, 30, 35, 50,
             100, 200, 500, 1000]*1e6
isoP_kw = {"ls": ":", "color": "blue", "lw": 0.8}

# # Isoenthalpic lines to plot
isoh = np.linspace(hmin, hmax, 10)
# isoh = np.r_[0, 250, 500, 750, 1000, 1250, 1500, 1750]
isoh_kw = {"ls": ":", "color": "orange", "lw": 0.8}

# Isoentropic lines to plot
isos = np.linspace(smin, smax, 10)
isos_kw = {"ls": ":", "color": "brown", "lw": 0.8}

# Isochor lines to plot
isov = np.logspace(
    np.log10(1/fluid()._constants["rhomax"]/fluid.M), np.log10(Vmin.v), 10)
# isov = np.r_[1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.2,
#              0.5, 1, 2, 5, 10, 20, 50]
isov_kw = {"ls": ":", "color": "darkgreen", "lw": 0.8}

# Validity region
validity_kw = {"facecolor": "#F2F2F2"}

# Ideal curves
ideal_kw = {"ls": "--", "color": "black", "lw": 0.8}
label_kw = {"size": "small", "ha": "center", "va": "bottom"}
###############################################################################


# Define plot
fig = plt.figure(figsize=(15, 15))

ax1 = fig.add_subplot(4, 2, 1)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_title("Ideal curves")
ax1.set_xlabel("Tr")
ax1.set_ylabel("Pr")

ax_PIP = fig.add_subplot(4, 2, 2)
ax_PIP.set_title("Phase Identification Parameter")
ax_PIP.set_xlabel("T, K")
ax_PIP.set_ylabel("$\Pi$")

ax_Ph = fig.add_subplot(4, 2, 3)
ax_Ph.set_yscale("log")
ax_Ph.set_title("P-h Diagram")
ax_Ph.set_xlabel("h, kJ/kg")
ax_Ph.set_ylabel("P, MPa")

ax_Ts = fig.add_subplot(4, 2, 4)
ax_Ts.set_title("T-s Diagram")
ax_Ts.set_xlabel("s, kJ/kgK")
ax_Ts.set_ylabel("T, K")

ax_Tv = fig.add_subplot(4, 2, 5)
ax_Tv.set_xscale("log")
ax_Tv.set_title("T-v Diagram")
ax_Tv.set_xlabel("v, m³/kg")
ax_Tv.set_ylabel("T, K")

ax_PT = fig.add_subplot(4, 2, 6)
ax_PT.set_yscale("log")
ax_PT.set_title("P-T diagram")
ax_PT.set_xlabel("T, K")
ax_PT.set_ylabel("P, MPa")

ax_hs = fig.add_subplot(4, 2, 7)
ax_hs.set_title("h-s Diagram")
ax_hs.set_xlabel("s, kJ/kgK")
ax_hs.set_ylabel("h, kJ/kg")

ax_vu = fig.add_subplot(4, 2, 8)
ax_vu.set_yscale("log")
ax_vu.set_title("v-u Diagram")
ax_vu.set_xlabel("u, kJ/kg")
ax_vu.set_ylabel("v, m³/kg")


Ts = list(np.concatenate([
    np.linspace(max(fluid.Tt, fluid()._constants["Tmin"]), 0.9 * Tc, points, endpoint=False),
    np.linspace(0.9 * Tc, 0.99 * Tc, points, endpoint=False),
    np.linspace(0.99 * Tc, 0.999 * Tc, points, endpoint=False),
    np.linspace(0.999 * Tc, Tc, points, endpoint=True)]))
Pl = list(np.concatenate([
    np.logspace(
        np.log10(Pt.P), np.log10(0.9 * fluid.Pc), points, endpoint=False),
    np.linspace(0.9 * fluid.Pc, 0.99 * fluid.Pc, points, endpoint=False),
    np.linspace(0.99 * fluid.Pc, 1.1 * fluid.Pc, points, endpoint=False),
    np.logspace(
        np.log10(1.1 * fluid.Pc), np.log10(Pmax), points, endpoint=True)]))
Tl = list(np.concatenate([
    np.linspace(Tmin, 0.9 * Tc, points, endpoint=False),
    np.linspace(0.9 * Tc, 1.1 * Tc, points, endpoint=False),
    np.logspace(np.log10(1.1 * Tc), np.log10(Tmax), points)]))


# Melting pressure
if fluid is not IAPWS95 and fluid._melting:
    print("Calculating melting line...")
    if fluid in (mEoS.NH3, mEoS.Xe):
        Tm = np.logspace(np.log10(fluid._melting1["Tmin"]),
                         np.log10(fluid._melting["Tmax"]), points * 10)
    else:
        Tm = np.logspace(np.log10(fluid._melting["Tmin"]),
                         np.log10(fluid._melting["Tmax"]), points * 10)
    Pm = [fluid._Melting_Pressure(t) for t in Tm]
    mel = []
    for t, p in zip(Tm, Pm):
        Ps = fluid()._Vapor_Pressure(t)
        if abs(p-Ps)/Ps < 1e-2:
            mel.append(fluid(T=t, x=0))
        else:
            mel.append(fluid(T=t, P=p))

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

    ax1.plot(Tm / Tc, [p / fluid.Pc for p in Pm], **mel_kw)
    ax_Ph.plot(h, P, label="Melting line", **mel_kw)
    ax_Ts.plot(s, T, label="Melting line", **mel_kw)
    ax_Tv.plot(v, T, label="Melting line", **mel_kw)
    ax_PT.plot(T, P, label="Melting line", **mel_kw)
    ax_hs.plot(s, h, label="Melting line", **mel_kw)
    ax_vu.plot(u, v, label="Melting line", **mel_kw)


# Calculate saturation line
print("Calculating saturation lines...")
liq = [fluid(T=t, x=0) for t in Ts]
vap = [fluid(T=t, x=1) for t in Ts]

xliq = [line.T / Tc for line in liq]
yliq = [line.P / fluid.Pc for line in liq]
ax1.plot(xliq, yliq, **isosat_kw)
xvap = [v.T / Tc for v in vap]
yvap = [v.P / fluid.Pc for v in vap]
ax1.plot(xvap, yvap, **isosat_kw)

PIPliq = [line.PIP for line in liq]
hliq = [line.h.kJkg for line in liq]
Pliq = [line.P.MPa for line in liq]
sliq = [line.s.kJkgK for line in liq]
vliq = [line.v for line in liq]
uliq = [line.u.kJkg for line in liq]
ax_PIP.plot(Ts, PIPliq, label="Saturated Liquid", **PIP_kw)
ax_Ph.plot(hliq, Pliq, label="Saturated Liquid", **isosat_kw)
ax_Ts.plot(sliq, Ts, label="Saturated Liquid", **isosat_kw)
ax_Tv.plot(vliq, Ts, label="Saturated Liquid", **isosat_kw)
ax_PT.plot(Ts, Pliq, label="Saturated Liquid", **isosat_kw)
ax_hs.plot(sliq, hliq, label="Saturated Liquid", **isosat_kw)
ax_vu.plot(uliq, vliq, label="Saturated Liquid", **isosat_kw)

PIPvap = [line.PIP for line in vap]
hvap = [v.h.kJkg for v in vap]
Pvap = [v.P.MPa for v in vap]
svap = [v.s.kJkgK for v in vap]
vvap = [v.v for v in vap]
uvap = [v.u.kJkg for v in vap]
ax_PIP.plot(Ts, PIPvap, label="Saturated Vapor", **PIP_kw)
ax_Ph.plot(hvap, Pvap, label="Saturaded Vapor", **isosat_kw)
ax_Ts.plot(svap, Ts, label="Saturaded Vapor", **isosat_kw)
ax_Tv.plot(vvap, Ts, label="Saturaded Vapor", **isosat_kw)
ax_hs.plot(svap, hvap, label="Saturaded Vapor", **isosat_kw)
ax_vu.plot(uvap, vvap, label="Saturaded Vapor", **isosat_kw)


# Set graphics limit to only saturation region
hmax = max(hvap)
hmax += abs(hmax * 0.1)
umax = max(uvap)
umax += abs(umax * 0.1)
smax = max(svap)
smax += abs(smax * 0.1)

ax_PIP.set_ylim(-10, 15)
ax_PIP.set_xlim(Tmin, 2*fluid.Tc)
ax_Ph.set_xlim(Pt.h.kJkg, hmax)
ax_Ph.set_ylim(Pt.P.MPa, Pc.P.MPa * 10)
ax_Ts.set_xlim(Pt.s.kJkgK, smax)
ax_Ts.set_ylim(Pt.T, Pc.T * 1.1)
ax_Tv.set_xlim(Pt.v, Vt.v)
ax_Tv.set_ylim(Pt.T, Pc.T * 1.1)
ax_PT.set_xlim(Tmin, Tmax)
ax_PT.set_ylim(Pmin/1e6, Pmax/1e6)
ax_hs.set_xlim(Pt.s.kJkgK, smax)
ax_hs.set_ylim(Pt.h.kJkg, hmax)
ax_vu.set_xlim(Pt.u.kJkg, umax)
ax_vu.set_ylim(Pt.v, Vt.v)

# Ancillary equation
Tanc = np.concatenate([
    np.linspace(fluid.Tt, 0.9 * Tc, 10),
    np.linspace(0.9 * Tc, Tc, 10, endpoint=True)])
Panc = [fluid()._Vapor_Pressure(t) / fluid.Pc for t in Tanc]
ax1.plot(Tanc / Tc, Panc, label="Ancillary Vapor Pressure", **anc_kw)

lanc = [1 / fluid()._Liquid_Density(t) for t in Tanc]
vanc = [1 / fluid()._Vapor_Density(t) for t in Tanc]
ax_Tv.plot(lanc, Tanc, label="Ancillary Liquid Density", **anc_kw)
ax_Tv.plot(vanc, Tanc, label="Ancillary Vapor Density", **anc_kw)

# Calculate isoquality lines
print("Calculating isoquality lines...")
for q in isoq:
    txt = f"x: {q:0.5g}"
    print("    %s" % txt)
    pts = [fluid(T=t, x=q) for t in Ts]

    x = [p.h.kJkg for p in pts]
    y = [p.P.MPa for p in pts]
    ax_Ph.plot(x, y, label=txt, **isoq_kw)

    x = [p.s.kJkgK for p in pts]
    y = [p.T for p in pts]
    ax_Ts.plot(x, y, label=txt, **isoq_kw)

    x = [p.v for p in pts]
    y = [p.T for p in pts]
    ax_Tv.plot(x, y, label=txt, **isoq_kw)

    x = [p.s.kJkgK for p in pts]
    y = [p.h.kJkg for p in pts]
    ax_hs.plot(x, y, label=txt, **isoq_kw)

    x = [p.u.kJkg for p in pts]
    y = [p.v for p in pts]
    ax_vu.plot(x, y, label=txt, **isoq_kw)


# Calculate isotherm lines
print("Calculating isotherm lines...")
for T in isoT:
    txt = "T: %s K" % T
    print("    %s" % txt)

    # Calculate the saturation point if available
    sat_pnt = []
    if T < Tc:
        for x in np.linspace(1, 0.1, 10, endpoint=False):
            sat_pnt.append(fluid(T=T, x=x))
        for x in np.linspace(0.1, 0, 11):
            sat_pnt.append(fluid(T=T, x=x))
        sat = True
    else:
        sat = False

    pts = []
    melting = True

    for p in Pl:
        point = fluid(P=p, T=T)
        if point.status in (1, 3):

            # Discard point below the melting line
            if melting and fluid._melting:
                pm = fluid()._Melting_Pressure(point.T)
                if pm < point.P:
                    point = fluid(T=T, P=pm)
                    melting = False

            # Add saturation point if neccesary
            if sat and T < Tc and sat_pnt and point.s < sat_pnt[-1].s:
                for pt in sat_pnt:
                    pts.append(pt)
                sat = False
            pts.append(point)

            if not melting:
                break

    h = []
    P = []
    s = []
    u = []
    v = []
    for p in pts:
        if p.status in (1, 3):
            h.append(p.h.kJkg)
            P.append(p.P.MPa)
            s.append(p.s.kJkgK)
            u.append(p.u.kJkg)
            v.append(p.v)
    ax_Ph.plot(h, P, label=txt, **isoT_kw)
    ax_hs.plot(s, h, label=txt, **isoT_kw)
    ax_vu.plot(u, v, label=txt, **isoT_kw)


# Calculate isobar lines
print("Calculating isobar lines...")
for P in np.concatenate([isoP, isoP_PIP]):
    txt = "P: %0.5g MPa" % (P/1e6)
    print("    %s" % txt)

    # Calculate the saturation point if available
    if P < fluid.Pc:
        sat_pnt = []
        for x in np.linspace(1, 0.1, 10):
            sat_pnt.append(fluid(P=P, x=x))
        for x in np.linspace(0.1, 0, 11):
            sat_pnt.append(fluid(P=P, x=x))
        sat = True
    else:
        sat = False

    pts = []
    for t in Tl:
        point = fluid(P=P, T=t)

        if point.status in (1, 3):

#             if point.rhoM > point._constants["rhomax"]:
#                 continue

            # Discard point below the melting line
            if fluid._melting:
                pm = fluid()._Melting_Pressure(point.T)
                if pm < point.P:

                    def f_p(t):
                        pm = fluid()._Melting_Pressure(t)
                        return pm - P

                    T = fsolve(f_p, point.T)[0]
                    point = fluid(T=T, P=P)

            # Add saturation point if neccesary
            if sat and P < fluid.Pc and sat_pnt and point.s > sat_pnt[-1].s:
                for p in sat_pnt[::-1]:
                    if p.status in (1, 3):
                        pts.append(p)
                sat = False

            pts.append(point)

    h = []
    T = []
    s = []
    u = []
    v = []
    PIP = []
    for p in pts:
        if p.status in (1, 3):
            h.append(p.h.kJkg)
            T.append(p.T)
            s.append(p.s.kJkgK)
            u.append(p.u.kJkg)
            v.append(p.v)
            if p.x in (0, 1):
                PIP.append(p.PIP)
            else:
                PIP.append(None)
    if P in isoP:
        ax_Ts.plot(s, T, label=txt, **isoP_kw)
        ax_Tv.plot(v, T, label=txt, **isoP_kw)
        ax_hs.plot(s, h, label=txt, **isoP_kw)
        ax_vu.plot(u, v, label=txt, **isoP_kw)
    ax_PIP.plot(T, PIP, label=txt, **PIP_kw)


# Calculate isoenthalpic lines
print("Calculating isoenthalpic lines...")
for h in isoh:
    txt = "h: %0.5g kJ/kg" % (h)
    print("    %s" % txt)

    pts = []
    for p in Pl:
        if pts:
            T0 = pts[-1].T
            rho0 = pts[-1].T
        else:
            T0 = None
            rho0 = None

        point = fluid(P=p, h=h*1000)

        if point.status in (1,3):

            # Discard point below the melting line
            if fluid in (mEoS.NH3, mEoS.Xe):
                Tmin = fluid._melting1["Tmin"]
            else:
                Tmin = fluid._melting["Tmin"]

            if fluid._melting and Tmin < point.T:
                pm = fluid()._Melting_Pressure(point.T)
                if pm < point.P:

                    def f_h(t):
                        pm = fluid()._Melting_Pressure(t)
                        return fluid(T=t, P=pm).h.kJkg - h

                    T = fsolve(f_h, point.T)[0]
                    pm = fluid()._Melting_Pressure(T)
                    point = fluid(T=T, P=pm)
                    pts.append(point)
                    break

            pts.append(point)

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
    ax_Ts.plot(s, T, label=txt, **isoh_kw)
    ax_Tv.plot(v, T, label=txt, **isoh_kw)
    ax_PT.plot(T, P, label=txt, **isoh_kw)

# Calculate isoentropic lines
print("Calculating isoentropic lines...")
for s in isos:
    txt = "s: %0.5g kJ/kgK" % s
    print("    %s" % txt)

    phase = True
    ts = Pt.Tt

    # Find saturated liquid point
    if Pt.s.kJkgK > s:
        def f_sm(ti):
            if ti < fluid.Tt:
                ti = fluid.Tt
            return fluid(T=ti, x=0).s - s*1000
        ts = fsolve(f_sm, (Tc-fluid.Tt)/2)
    else:
        phase = False

    pts = []
    for t in Tl:
        point = fluid(T=t, s=s*1000)

        if point.status in (1,3):
            # Discard point over the melting line
            if fluid in (mEoS.NH3, mEoS.Xe):
                Tmin = fluid._melting1["Tmin"]
            else:
                Tmin = fluid._melting["Tmin"]

            if fluid._melting and Tmin < point.T:
                pm = fluid()._Melting_Pressure(point.T)
                if pm < point.P:

                    # Fuction to find the intersection between isoentropic
                    # line and melting line
                    def f_s(ti):
                        ti = max(ti, Tmin)
                        pm = fluid()._Melting_Pressure(ti)
                        return fluid(T=ti, P=pm).s.kJkgK-s

                    T = fsolve(f_s, point.T)
                    pm = fluid()._Melting_Pressure(T)
                    point = fluid(T=T, P=pm)
                    if point.status in (1,3):
                        pts.append(point)

                    break

            # Add saturation point and near to that
            if phase and point.T > ts:
                phase = False

                if not pts:
                    tmin = Tmin
                else:
                    tmin = pts[-1].T
                ti = np.concatenate([
                    np.linspace(Tmin, ts, points, endpoint=False),
                    np.linspace(ts, point.T, points*10, endpoint=False)])

                for ts in ti:
                    ptss = fluid(T=ts, s=s*1000)
                    if ptss.status == 1:
                        pts.append(ptss)

            pts.append(point)

    T = []
    h = []
    P = []
    u = []
    v = []
    for p in pts:
        T.append(p.T)
        h.append(p.h.kJkg)
        u.append(p.u.kJkg)
        P.append(p.P.MPa)
        v.append(p.v)

    ax_Ph.plot(h, P, label=txt, **isos_kw)
    ax_Tv.plot(v, T, label=txt, **isos_kw)
    ax_vu.plot(u, v, label=txt, **isos_kw)
    ax_PT.plot(T, P, label=txt, **isos_kw)

# Calculate isochor lines
print("Calculating isochor lines...")
for v in isov:
    txt = "v: %0.5g m³/kg" % v
    print("    %s" % txt)

    pts = [fluid(T=t, v=v) for t in Tl]

    h = []
    T = []
    s = []
    u = []
    P = []
    for p in pts:
        if p.status == 1:

            # Discard point over the melting line
            if fluid._melting:
                pm = fluid()._Melting_Pressure(p.T)
                if pm < p.P:

                    # Function to find the intersection between isochor line
                    def f_v(ti):
                        ti = min(ti, fluid.eq[0]["Tmax"])
                        pm = fluid()._Melting_Pressure(ti)
                        return fluid(T=ti, P=pm).v-v

                    t = None
                    for to in [p.Tt, p.T, p.Tc]:
                        rinput = fsolve(f_v, p.T, full_output=True)
                        if rinput[2] == 1 and abs(rinput[1]["fvec"]) < 1e-5:
                            t = rinput[0]
                            break

                    if t:
                        pm = fluid()._Melting_Pressure(t)
                        p = fluid(T=t, P=pm)
                    else:
                        continue

            h.append(p.h.kJkg)
            T.append(p.T)
            s.append(p.s.kJkgK)
            u.append(p.u.kJkg)
            P.append(p.P.MPa)

    ax_Ph.plot(h, P, label=txt, **isov_kw)
    ax_PT.plot(T, P, label=txt, **isov_kw)
    ax_Ts.plot(s, T, label=txt, **isov_kw)
    ax_hs.plot(s, h, label=txt, **isov_kw)

# Show the validity range of equation
Prmin = Pt.P / fluid.Pc
Trmin = Tmin / Tc
Prmax = Pmax / fluid.Pc
Trmax = Tmax / Tc

if fluid._melting:
    vert = [(Trmin, Prmin)]
    for t, p in zip(Tm, Pm):
        vert.append((t / Tc, min(p / fluid.Pc, Prmax)))
else:
    vert = [(Trmin, Prmin), (Trmin, Prmax)]

vert.append((Trmax, Prmax))
vert.append((Trmax, Prmin))

poly = Polygon(vert, **validity_kw)
ax1.add_patch(poly)


# Ideal Curves
print("Calculating ideal curves...")
T = np.logspace(np.log10(Tmin), 4, points*5, endpoint=True)

for curva in ["ideal", "boyle", "joule-thomson", "joule"]:
    print("    %s" % curva)
    T_line = []
    P_line = []
    pmin = Pmin
    for t in T:
        P = fluid()._IdealCurve(curva, t)
        if P is None or P < 0:
            continue

        if fluid._melting:
            Pm = fluid()._Melting_Pressure(t)
        else:
            Pm = 1e999

        # Discard point below the saturaion line
        if curva in ["boyle", "joule-thomson"]:
            if t < Tc:
                st = fluid(T=max(t, Tmin), x=0)
                if st.P > P:
                    P = st.P
                    pmin = P

                # Discard any boyle point over saturation
                if curva == "boyle":
                    continue

        if 0 < P <= Pm:
            T_line.append(t/Tc)
            P_line.append(P/fluid.Pc)

    # Append a point in a x axis
    T_line.append(T_line[-1])
    P_line.append(Pt.Pr)

    ax1.plot(T_line, P_line, label=curva, **ideal_kw)

    # Set ideal curve label
    if curva == "joule-thomson":
        txt = "j-t"
        ax1.set_ylim(bottom=0.1*pmin/fluid.Pc)
    else:
        txt = curva

    # Location of text, in boyle and j-t using the maximum, in other a
    # position at left of plot
    if curva in ["joule", "boyle", "joule-thomson"]:
        i = P_line.index(max(P_line))
    else:
        i = int(len(T_line)*0.4)

    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()

    fx = (np.log10(T_line[i+1])-np.log10(T_line[i-1]))/(
        np.log10(xmax)-np.log10(xmin))
    fy = (np.log10(P_line[i+1])-np.log10(P_line[i-1]))/(
        np.log10(ymax)-np.log10(ymin))
    rot = np.arctan(fy/fx)*360/2/np.pi
    ax1.annotate(txt, (T_line[i], P_line[i]), rotation=rot, **label_kw)


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
ax_hs.legend(handles=legend_tv, loc='lower right')

fig.set_tight_layout(True)

if args.matplotlib:
    plt.show()

file = "/home/jjgomera/Descargas/ImagesPychemqt/%s.png" % fluid.__name__
fig.savefig(file)

if args.show:
    subprocess.call(["xdg-open", file])
else:
    print("Image generated at %s" % file)
