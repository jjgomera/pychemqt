#!/usr/bin/python3
# -*- coding: utf-8 -*-


from PyQt5.QtWidgets import QApplication
from scipy import exp, log10

from lib import unidades
from lib.compuestos import Componente
from lib.petro import Petroleo
from lib.physics import R_atml
from lib.sql import databank


class Crudo(Petroleo):
    """Class to model a hypotetical component from a defined crude oil
    Clase que define una fracción de petroleo a partir de la base de datos

    Parameters
    ------------
    index : index
        index of crude in crude databank
    Cplus : index
        Number of fraction to split the crude oil for heavier component than C6

    Notes
    -----
    The resultant instancecan be used as hypotetical component as Petroleo.
    """
    kwargs = Petroleo.kwargs.copy()
    kwarg = {"index": 0,
             "Cplus": 0,

             "Rgo": 0.0,
             "gas": None,
             "water": None}
    kwargs.update(kwarg)

    status = 0
    _bool = False
    msg = ""
    hasCurve = False
    hasSG = True
    hasRefraction = False

    def isCalculable(self):
        if self.kwargs["index"]:
            self.status = 1
            self.msg = ""
            return True
        else:
            self.status = 0
            self.msg = QApplication.translate("pychemqt", "Undefined petrol")

    def calculo(self):
        id = self.kwargs["index"]
        databank.execute("SELECT * FROM CrudeOil WHERE id=='%i'" % id)
        prop = databank.fetchone()

        API = prop[4]
        SG = 141.5/(API+131.5)
        PP = unidades.Temperature(prop[8], "F")
        v100 = prop[10]
        Tb = unidades.Temperature(
            (-753.-136*(1.-exp(-.15*v100))+572*SG-0.0512*v100+PP.R)/0.139, "R")

        self.definicion = 1
        self.kwargs["name"] = ", ".join(prop[1:3])
        self.kwargs["SG"] = SG
        self.kwargs["Tb"] = Tb
        self.kwargs["S"] = prop[5]
        self.kwargs["N"] = prop[6]
        self.kwargs["v100"] = prop[10]
        Petroleo.calculo(self)

        self.vanadium = prop[11]
        self.nickel = prop[12]
        self.carbonResid = prop[13]
        self.asphaltene = prop[14]
        self.S2 = prop[6]
        self.H2S = prop[18]
        self.nNeutralization = prop[19]
        self.ash = prop[21]
        self.salt = prop[22]
        self.water = prop[20]
        self.NPentane = prop[15]
        if prop[16]:
            self.reidVP = unidades.Pressure(prop[16], "psi")
        else:
            self.reidVP = None
        if prop[17]:
            self.FlashP = unidades.Temperature(prop[17], "F")
        self.PourP = PP

        self.C1 = prop[23]/100.
        self.C2 = prop[24]/100.
        self.C3 = prop[25]/100.
        self.iC4 = prop[26]/100.
        self.nC4 = prop[27]/100.
        self.iC5 = prop[28]/100.
        self.nC5 = prop[29]/100.

        # SGo = 0.7
        # SG_ = (SG-SGo)/SGo
        # B = 3.
        # A = SG_**3/0.619**3
        # Cplus = int(self.kwargs["Cplus"])
        # Tbi=[unidades.Temperature(1090-exp(6.9955-0.11193*Nc**(2./3))) for Nc in range(6, Cplus)]
        # SGi=[1.07-exp(3.65097-3.8864*Nc**0.1) for Nc in range(6, Cplus)]
        # x=[1-1/exp(A/B*(SG-SGo)**B/SGo**B) for SG in SGi]
        # Mi=[prop_Riazi_Alsahhaf(1, g, reverse=True) for g in SGi]
        # APIi=[141.5/SG-131.5 for SG in SGi]
        # Kwi=[Tbi[i].R**(1./3)/SGi[i] for i in range(len(Mi))]
        # di=[unidades.Density(prop_Riazi_Alsahhaf(2, M), "gcc") for M in Mi]
        # Ii=[prop_Riazi_Alsahhaf(3, M) for M in Mi]
        # Tci=[unidades.Temperature(Tbi[i]/prop_Riazi_Alsahhaf(4, M)) for i, M in enumerate(Mi)]
        # Pci=[unidades.Pressure(prop_Riazi_Alsahhaf(5, M), "bar") for M in Mi]
        # Vci=[unidades.SpecificVolume(1/prop_Riazi_Alsahhaf(6, M), "ccg") for M in Mi]
        # Wi=[prop_Riazi_Alsahhaf(7, M) for M in Mi]
        # Tensioni=[unidades.Tension(prop_Riazi_Alsahhaf(8, M), "dyncm") for M in Mi]
        # ParSoli=[unidades.SolubilityParameter(prop_Riazi_Alsahhaf(9, M), "calcc") for M in Mi]


    def pb_Standing(self, T):
        """Standing, M.B.: Volumetric and Phase Behavior of Oil Field Hydrocarbon Systems, SPE, Dallas (1977)"""
        t=unidades.Temperature(T)
        F=(self.Rgo.ft3bbl/self.gas.SG)**0.83*10**(0.00091*t.F-0.0125*self.API)
        return unidades.Pressure(18.2*(F-1.4), "psi")

    def pb_Lasater(self, T):
        """Lasater, J.A: "Bubble Point Pressure Correlation," Trans., AIME (1958) 213, 379-381"""
        t=unidades.Temperature(T)
        if self.API<=40:
            M=630-10.*self.API
        else:
            M=73110.*self.API**-1.562
        yg=self.Rgo.ft3bbl/379.3/(self.Rgo.ft3bbl/379.3+350*self.SG/M)
        if yg<=0.6:
            pb=0.679*exp(2.786*yg)-0.323
        else:
            pb=8.26*yg**3.56+1.95
        return unidades.Pressure(pb*t.R/self.gas.SG, "psi")

    def pb_Vazquez_Beggs(self, T, ts=350, ps=10):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3=0.0362, 1.0937, 25.724
        else:
            C1, C2, C3=0.0178, 1.187, 23.931

        gravity_corr=self.gas.SG*(1.+5.912e-5*self.API*ts.F*log10(ps.psi/114.7))
        return unidades.Pressure((self.Rgo.ft3bbl/C1/gravity_corr/exp(C3*self.API/T.R))**(1./C2), "psi")

    def pb_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        F=(self.Rgo.ft3bbl/self.gas.SG)**0.816*t.F**0.172/self.API**0.989
        return unidades.Pressure(10**(1.7669+1.7447*log10(F)-0.30218*log10(F)**2), "psi")

    def pb_Total(self, T):
        """TOTAL Compagnie Francaise Des Petroles: "Proyectos de Inyección de Fluidos - Correlaciones PVT para Crudos del Oriente de Venezuela," S.A. MENEVEN, Sept. 1983"""
        t=unidades.Temperature(T)
        if self.API<=10:
            C1, C2, C3, C4=12.847, 0.9636, 0.000993, 0.03417
        elif self.API<=35:
            C1, C2, C3, C4=25.2755, 0.7617, 0.000835, 0.011292
        else:
            C1, C2, C3, C4=216.4711, 0.6922, -0.000427, 0.02314
        return unidades.Pressure(C1*(self.Rgo.ft3bbl/self.gas.SG)**C2*10**(C3*t.F-C4*self.API), "psi")

    def pb_Al_Marhoun(self, T):
        """Al-Marhoun, M.A.: "PVT Correlation for Middle East Crude Oils," J. Pet. Tech (May 1988), 650-666"""
        t=unidades.Temperature(T)
        return unidades.Pressure(5.38088e-3*self.Rgo.ft3bbl**0.715082*self.gas.SG**-1.87784*self.SG**3.1437*t.R**1.32657, "psi")

    def pb_Dokla_Osman(self, T):
        """Dokla, M.E. and Osman, M.E.: "Correlation of PVT properties for UAE Crudes," Trans., AIME (1992) 293, 41-46"""
        t=unidades.Temperature(T)
        return unidades.Pressure(0.836386e4*self.Rgo.ft3bbl**0.724047*self.gas.SG**-1.01049*self.SG**0.107991*t.R**-0.952584, "psi")

    def pb_Petrosky_Farshad(self, T):
        """Petrosky, G.E., Jr. and Farshad, F.F.: "Pressure-Volume-Temperature Correlations for Gulf of Mexico Crude Oils," paper SPE 26644 presented at the 68th Annual Technical Conference and Exhibition, Houston, Texas, Oct. 3-6,1993."""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.5774/self.gas.SG**0.8439*10**(4.561e-5*t.F**1.3911-7.916e-4*self.API**1.541)
        return unidades.Pressure(112.727*(F-12.34), "psi")

    def pb_Kartoatmodjo_Schmidt(self, T, ts=350, ps=10):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3, C4=0.05958, 0.7972, 13.1405, 0.9986
        else:
            C1, C2, C3, C4=0.0315, 0.7587, 11.2895, 0.9143

        gravity_corr=self.gas.SG*(1.+0.1595*self.API**0.4078*ts.F**-0.2506*log10(ps.psi/114.7))
        return unidades.Pressure((self.Rgo.ft3bbl/C1/gravity_corr**C2/10**(C3*self.API/t.R))**C4, "psi")

    def pb(self, T):
        """Presión de burbujeo del crudo"""
        t=unidades.Temperature(T)
        metodo_pb=[self.pb_Standing, self.pb_Lasater, self.pb_Vazquez_Beggs, self.pb_Glaso, self.pb_Total, self.pb_Al_Marhoun, self.pb_Dokla_Osman, self.pb_Petrosky_Farshad, self.pb_Kartoatmodjo_Schmidt][Preferences.getint("petro", "pb")]
        pb=metodo_pb(t)
        CN2=1.+((-2.65e-4*self.API+5.5e-3)*t.F+(0.0931*self.API-0.895))*self.gas.N2
        CCO2=1.-693.8*self.gas.CO2*t.F**-1.553
        CH2S=1.-(0.9035+0.0015*self.API)*self.gas.H2S+0.019*(45-self.API)*self.gas.H2S**2
        return unidades.Pressure(CN2*CH2S*CCO2*pb)


    def B_Standing(self, T):
        """Standing, M.B.: Volumetric and Phase Behavior of Oil Field Hydrocarbon Systems, SPE, Dallas (1977)"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl*(self.gas.SG/self.SG)**0.5+1.25*t.F
        return 0.9759+12e-5*F**1.2

    def B_Vazquez_Beggs(self, T, ts=350, ps=10):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")
        if self.API<=30:
            C1, C2, C3=4.667e-4, 1.751e-5, -1.8106e-6
        else:
            C1, C2, C3=4.67e-4, 1.1e-5, 1.337e-9

        gravity_corr=self.gas.SG*(1.+5.912e-5*self.API*ts.F*log10(ps.psi/114.7))
        return 1.+C1*self.Rgo.ft3bbl+C2*(t.F-60)*self.API/gravity_corr+C3*self.Rgo.ft3bbl*(t.F-60)*self.API/gravity_corr

    def B_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl(self.gas.SG/self.SG)**0.526*+0.968*t.F
        return 1.+10**(-6.58511+2.91329*log10(F)-0.27683*log10(F)**2)

    def B_Total(self, T):
        """TOTAL Compagnie Francaise Des Petroles: "Proyectos de Inyección de Fluidos - Correlaciones PVT para Crudos del Oriente de Venezuela," S.A. MENEVEN, Sept. 1983"""
        t=unidades.Temperature(T)
        return 1.022+4.857e-4*self.Rgo.ft3bbl-2.009e-6*(t.F-60)*self.API/self.gas.SG+17.569e-9*self.Rgo.ft3bbl*(t.F-60)*self.API/self.gas.SG

    def B_Al_Marhoun(self, T):
        """Al-Marhoun, M.A.: "PVT Correlation for Middle East Crude Oils," J. Pet. Tech (May 1988), 650-666"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.74239*self.gas.SG**0.323294*self.SG**-1.20204
        return 0.497069+0.862963e-3*t.R+0.182594e-2*F+0.318099e-5*F**2

    def B_Dokla_Osman(self, T):
        """Dokla, M.E. and Osman, M.E.: "Correlation of PVT properties for UAE Crudes," Trans., AIME (1992) 293, 41-46"""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.773572*self.gas.SG**0.40402*self.SG**-0.882605
        return 0.431935e-1+0.156667e-2*t.R+0.139775e-2*F+0.380525e-5*F**2

    def B_Petrosky_Farshad(self, T):
        """Petrosky, G.E., Jr. and Farshad, F.F.: "Pressure-Volume-Temperature Correlations for Gulf of Mexico Crude Oils," paper SPE 26644 presented at the 68th Annual Technical Conference and Exhibition, Houston, Texas, Oct. 3-6,1993."""
        t=unidades.Temperature(T)
        F=self.Rgo.ft3bbl**0.3738*self.gas.SG**0.2914/self.SG**0.6265+0.24626*t.F**0.5371
        return 1.0113+7.2046e-5*F**3.0936

    def B_Kartoatmodjo_Schmidt(self, T, ts=350, ps=10):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        ts=unidades.Temperature(ts)
        ps=unidades.Pressure(ps, "atm")

        gravity_corr=self.gas.SG*(1.+0.1595*self.API**0.4078*ts.F**-0.2506*log10(ps.psi/114.7))
        F=self.Rgo.ft3bbl**0.755*gravity_corr**0.25*self.SG**-1.5+0.45*t.F
        return 0.98496+1e-4*F**1.5

    def B(self, T, P):
        """Factor volumétrico, relación entre el volumen a las condiciones del yacimiento y las condiciones normales a la presión de burbujeo"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        metodo_B=[self.pb_Standing, self.pb_Vazquez_Beggs, self.pb_Glaso, self.pb_Total, self.pb_Al_Marhoun, self.pb_Dokla_Osman, self.pb_Petrosky_Farshad, self.pb_Kartoatmodjo_Schmidt][Preferences.getint("petro", "Vol_factor")]
        Bpb=metodo_B(t)
        pb=self.pb(T)
        if p>pb:
            c=self.compressibilidad(T)
            return Bpb*exp(c*(pb-p))
        else:
            return Bpb


    def Mu_Gas(self, T):
        """Método de cáluclo de la viscosidad de vapores, API procedure 11B3.1, pag 1105"""
        #FIXME: no sale
        t=unidades.Temperature(T)
        return unidades.Viscosity(-0.0092696+t.R**0.5*(0.0010310+4.4507e-5*self.M**0.5)+1.1249e-5*self.M, "cP")
#        return unidades.Viscosity(-0.0092696+T**0.5*(0.001383-5.9712e-5*self.M**0.5)+1.1249e-5*self.M, "cP")


    def Mu_Beal(self, T):
        """Beal, C.: "The Viscosity of Air, Water, Natural Gas, crude Oil and its Associated Gases at Oil-Field Temperatures and Pressures," Trans., AIME (1946) 165, 94-115"""
        t=unidades.Temperature(T)
        a=10**(0.43+8.33/self.API)
        mu=(0.32+1.8e7/self.API**4.53)*(360/(t.F+200))**a
        return unidades.Viscosity(mu, "cP")

    def Mu_Beggs_Robinson(self, T):
        """Beggs, H.D. and Robinson, J.R.: "Estimating the Viscosity of Crude Oil Systems," J. Pet. Tech. Forum (Sept. 1975), 1140-1141"""
        t=unidades.Temperature(T)
        z=3.0324-0.02023*self.API
        y=10**z
        x=y*t.F**-1.163
        return unidades.Viscosity(10**x-1, "cP")

    def Mu_Glaso(self, T):
        """Glaso, O.: "Generalized Pressure-Volume-Temperature Correlations," J. Pet. Tech. (May 1980), 785-795"""
        t=unidades.Temperature(T)
        mu=3.141e10*t.F**-3.444*log10(self.API)**(10.313*log10(t.F)-36.447)
        return unidades.Viscosity(mu, "cP")

    def Mu_Egbogah(self, T):
        """Egbogah, E.O.: "An Improved Temperature-Viscosity Correlation for Crude Oil Systems," paper 83-34-32 presented at the 1983 Annual Technical Meeting of the Petroleum Society of CIM, Banff, Alberta, May 10-13, 1983"""
        t=unidades.Temperature(T)
        mu=1.8653-0.025086*self.API-0.5644*log10(t.F)
        return unidades.Viscosity(10**(10**mu)-1, "cP")

    def Mu_Kartoatmodjo_Schmidt(self, T):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        mu=16e8*t.F**-2.8177*log10(self.API)**(5.7526*log10(t.F)-26.9718)
        return unidades.Viscosity(mu, "cP")

    def Mu_Muerto(self, T):
        """Viscosidad de petroleos muertos (sin gas disuelto)"""
        metodos=[self.Mu_Beal, self.Mu_Beggs_Robinson, self.Mu_Glaso, self.Mu_Egbogah, self.Mu_Kartoatmodjo_Schmidt][Preferences.getint("petro", "mu_dead")]
        return metodos(T)

    def Mu_Chew_Connally(self, T, R):
        """Chew, J.N. and Connally, C.A. Jr.: "A Viscosity Correlation for Gas-Saturated Crude Oils," Trans, AIME (1959) 216, 23-25"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        A=10**(r.ft3bbl*(2.2e-7*r.ft3bbl-7.4e-4))
        b=0.68/10**(8.62e-5*r.ft3bbl)+0.25/10**(1.1e-3*r.ft3bbl)+0.062/10**(3.74e-3*r.ft3bbl)
        return unidades.Viscosity(A*muo.cP**b, "cP")

    def Mu_Beggs_Robinson_vivo(self, T, R):
        """Beggs, H.D. and Robinson, J.R.: "Estimating the Viscosity of Crude Oil Systems," J. Pet. Tech. Forum (Sept. 1975), 1140-1141"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        A=10.715*(r.ft3bbl+100)**-0.515
        b=5.44*(r.ft3bbl+150)**-0.338
        return unidades.Viscosity(A*muo.cP**b, "cP")

    def Mu_Kartoatmodjo_Schmidt_vivo(self, T, R):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        t=unidades.Temperature(T)
        r=unidades.V2V(R)
        muo=self.Mu_Muerto(T)
        b=10**(-0.00081*r.ft3bbl)
        A=(0.2001+0.8428*10**(-0.000845*r.ft3bbl))*muo.cP**(0.43+0.5165*b)
        return unidades.Viscosity(-0.06821+0.9824*A+40.34e-5*A**2, "cP")

    def Mu_Vivo(self, T, R):
        """Viscosidad de petroleos vivos (con gas disuelto)"""
        metodos=[self.Mu_Chew_Connally, self.Mu_Beggs_Robinson_vivo, self.Mu_Kartoatmodjo_Schmidt_vivo][Preferences.getint("petro", "mu_live")]
        return metodos(T, R)


    def Mu_Beal_presion(self, T, P, R):
        """Beal, C.: "The Viscosity of Air, Water, Natural Gas, crude Oil and its Associated Gases at Oil-Field Temperatures and Pressures," Trans., AIME (1946) 165, 94-115"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        mu=(0.024*muo.cP**1.6+0.038*muo.cP**0.56)*0.001*(p.psi-pb.psi)+muo.cP
        return unidades.Viscosity(mu, "cP")

    def Mu_Vazquez_Beggs(self, T, P, R):
        """Vázquez, M.E. and Beggs, H.D.: "Correlations for Fluid Physical Property Prediction," J.Pet. Tech. (June 1980), 968-970"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        m=2.6*p.psi**1.187*exp(-11.513-8.98e-5*p.psi)
        return unidades.Viscosity(muo.cP*(p.psi/pb.psi)**m, "cP")

    def Mu_Kartoatmodjo_Schmidt_presion(self, T, P, ):
        """Kartoatmodjo, T. and Schmidt, Z.: "Large Data Bank Improve Crude Physical Property Correlations," Oil and Gas J. (July 4, 1994) 51-55"""
        p=unidades.Pressure(P, "atm")
        pb=self.pb(T)
        muo=self.Mu_Vivo(T, R)
        mu=1.00081*muo.cP+1.127e-3*(p.psi-pb.psi)*(-65.17e-4*muo.cP**1.8148+0.038*muo.cP**1.59)
        return unidades.Viscosity(mu, "cP")

    def Mu_Presion(self, T, P):
        """Viscosidad de petroleos vivos (con gas disuelto)"""
        metodos=[self.Mu_Beal_presion, self.Mu_Vazquez_Beggs, self.Mu_Kartoatmodjo_Schmidt_presion][Preferences.getint("petro", "mu_live")]
        return metodos(T, P)



class Water(Componente):
    """Clase que define el agua específica que acompaña al petroleo, con las propiedades específicas"""
    def __init__(self):
        Componente.__init__(self, 62)

    def Solubilidad_Culberson_McKetta(self, T, P, S=0):
        """Culberson, O.L. and McKetta, J.J., Jr.: "Phase Equilibria in Hydrocarbon-Water Systems III - The solubility of Methane in Water at Pressures to 10,000 psia," Trans., AIME (1951) 192, 223-226
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=8.15839-6.12265e-2*t.F+1.91663e-4*t.F**2-2.1654e-7*t.F**3
        B=1.01021e-2-7.44241e-5*t.F+3.05553e-7*t.F**2-2.94883e-10*t.F**3
        C=(-9.02505+0.130237*t.F-8.53425e-4*t.F**2+2.34122e-6*t.F**3-2.37049e-9*t.F**4)*1e-7
        R=A+B*p.psi+C*p.psi**2
        Rs=R*10**(-0.0840655*S*t.F**-0.285854)
        return unidades.V2V(Rs, "ft3bbl")

    def Solubilidad_McCoy(self, T, P, S=0):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983).
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=2.12+3.45e-3*t.F-3.59e-5*t.F**2
        B=0.0107-5.26e-5*t.F+1.48e-7*t.F**2
        C=-8.75e-7+3.9e-9*t.F-1.02e-11*t.F**2
        R=A+B*p.psi+C*p.psi**2
        Rs=R*(1-(0.0753-1.73e-4*t.F)*S)
        return unidades.V2V(Rs, "ft3bbl")


    def Factor_Volumetrico_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990.
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        DVp=-1.0001e-2+1.33391e-4*t.F+5.50654e-7*t.F**2
        DVt=-1.95301e-9*p.psi*t.F-1.72834e-13*p.psi**2*t.F-3.58922e-7*p.psi-2.25341e-10*p.psi**2
        B=(1+DVp)*(1+DVt)
        return B*(1+S*(5.1e-8*p.psi+(5.47e-6-1.95e-10*p.psi)*(t.F-60)-(3.23e-8-8.5e-13*p.psi)*(t.F-60)**2))

    def Factor_Volumetrico_McCoy(self, T, P, S=0, R=1):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983).
        R: razon de gas disuelto
        S: salinidad en % en peso"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        if R==0:
            A=0.9947+5.8e-6*t.F+1.02e-6*t.F**2
            B=-4.228e-6+1.8276e-8*t.F-6.77e-11*t.F**2
            C=1.3e-10-1.3855e-12*t.F+4.285e-15*t.F**2
        else:
            A=0.9911+6.35e-5*t.F+8.5e-7*t.F**2
            B=-1.093e-6-3.497e-9*t.F+4.57e-12*t.F**2
            C=-5e-11+6.429e-13*t.F-1.43e-15*t.F**2
        B=A+B*p.psi+C*p.psi**2
        return B*(1+S*(5.1e-8*p.psi+(5.47e-6-1.95e-10*p.psi)*(t.F-60)-(3.23e-8-8.5e-13*p.psi)*(t.F-60)**2))

    def Rho(self, T, P, S=0):
        B=self.Factor_Volumetrico_McCain(T, P, S)
        s=S*1e7/58443
        g=1.+0.695e-6*s
        return unidades.Density(62.4*g/B, "lbft3")

    def Rho_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990."""
        B=self.Factor_Volumetrico_McCain(T, P, S)
        g=62.368+0.438603*S+1.60074e-3*S**2
        return unidades.Density(g/B, "lbft3")

    def Compresibilidad_Dodson_Standing(self, T, P, S=0, R=0):
        """Dodson, C.R. and Standing, M.B.: "Pressure-Volume-Temperature and Solubility RElations for Natural Gas-Water-Mixtures," Drill. and Prod. Prac., API (1944) 173-179"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        R=unidades.V2V(R)
        a=3.8546-1.34e-4*p.psi
        b=-0.01052+4.77e-7*p.psi
        c=3.9267e-5-8.8e-10*p.psi
        cw=(a+b*t.F+c*t.F**2)/1e6
        cr=1.+8.9e-3*R.ft3bbl
        cs=1+S**0.7*(-5.2e-2+2.7e-4*t.F-1.14e-6*t.F**2+1.121e-9*t.F**3)
        return cw*cr*cs

    def Compresibilidad_Osif(self, T, P, S=0):
        """Osif, T.L.: "The Effects of Salt, Gas, Temperature and Pressure on the Compressibility of Water," SPE Res.Eng. (Feb. 1988) 3, No.1. 175-181"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        S=S*1e4/58443
        return 1/(7.033*p.psi+541.5*S-537.*t.F+403300.)


    def Mu_Van_Wingen(self, T):
        """Van Wingen, N.: "Viscosity of Air, Water, Natural Gas, and Crude Oil at Varying Pressure and Temperatures," Secondary Recovery of Oil in the United States, API (1950) 127."""
        t=unidades.Temperature(T)
        return unidades.Viscosity(exp(1.003-1.479e-2*t.F+1.982e-5*t.F**2), "cP")

    def Mu_Mattews_Russel(self, T, P, S=0):
        """Mathews, C.S and Russel, D.G.: Pressure Buildup and Flow Text in Wells. Monograph Series. Society of Petroleum Engineers of AIME, Dallas (1967) 1, Appendix G."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=-0.04518+0.009313*S-0.000383*S**2
        B=70.634+0.09576*S**2
        mu=A+B/t.F
        f=1.+3.5e-12*p.psi**2*(t.F-40.)
        return unidades.Viscosity(mu*f, "cP")

    def Mu_McCain(self, T, P, S=0):
        """McCain, W.D., Jr: The Properties of Petroleum Fluids, 2nd ed. Tulsa, OK: PennWell Books, 1990."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=109.574-8.40564*S+0.313314*S**2+8.72213e-3*S**3
        B=-1.12166+2.63951e-2*S-6.79461e-4*S**2-5.47119e-5*S**3+1.55586e-6*S**4
        mu=A*t.F**B
        f=0.9994+4.0295e-5*p.psi+3.1062e-9*p.psi**2
        return unidades.Viscosity(mu*f, "cP")

    def Mu_McCoy(self, T, S=0):
        """McCoy, R.L.: Microcomputer Programs for Petroleum Engineers: Vol. 1, Reservoir Engineering and Formation Evaluation, Gulf Publishing Co., Houston (1983)."""
        t=unidades.Temperature(T)
        mu=0.02414*10**(247.8/(t-140))
        f=1.-1.87e-3*S**0.5+2.18e-4*S**2.5+(t.F**0.5-1.35e-2*t.F)*(2.76e-3*S-3.44e-4*S**1.5)
        return unidades.Viscosity(mu*f, "cP")

    def Tension_Jennings_Newman(self, T, P):
        """Jennings, H.Y., Jr. and Newman, G.I.L.:"The Effect of Temperature and Pressure on the Interfacial Tension of Water Against Mechane-Normal Decane Mixtures," Trans., AIME (1971) 251, 171-175."""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        A=79.1618-0.118978*t.F
        B=-5.28473e-3+9.87913e-6*t.F
        C=(2.33814-4.57194e-4*t.F-7.52678e-6*t.F**2)*1e-7
        return unidades.Tension(A+B*p.psi+C*p.psi**2, "dyncm")


class Natural_Gas(object):
    """Clase que define un gas natural como fracción desconocida de petroleo"""
    def __init__(self, SG=None, composicion=[], wet=False, CO2=0., H2S=0., N2=0.):
        """
        g: gravedad específica
        composicion: en la forma de array [[componentes],[fracciones molares]]
        wet: parámetro opcional que indica si se trata de un gas húmedo (con una pequeña fracción de fase líquida
        CO2: fracción molar de CO2 en el gas
        H2S: fracción molar de H2S en el gas
        """
        self.SG=SG
        self.wet=wet
        self.CO2=CO2
        self.H2S=H2S
        self.N2=N2
        if wet:
            self.tpc=unidades.Temperature(187.+330.*SG-72.5*SG**2, "R")
            self.ppc=unidades.Pressure(706-51.7*SG-11.1*SG**2, "psi")
        else:
            self.tpc=unidades.Temperature(168.+325.*SG-12.5*SG**2, "R")
            self.ppc=unidades.Pressure(677+15.*SG-37.5*SG**2, "psi")
        if N2!=0:
            self.tpc, self.ppc=self.Critical_Carr_Kobayashi_Burrows()
        if CO2+H2S!=0.:
            self.tpc, self.ppc=self.Critical_Wichert_Aziz()
        self.M=28.96*SG

    def Critical_Wichert_Aziz(self):
        """Wichert, E., and K. Aziz. “Calculation of Z’s for Sour Gases.” Hydrocarbon Processing 51, no. 5 (1972): 119–122."""
        """Wichert, E., Aziz, K., 1972. Calculation of Z’s for sour gases. Hydrocarb. Process. 51 (5), 119–122."""
        A=self.CO2+self.H2S
        e=120.*(A**0.9-A**1.6)+15*(self.H2S**0.5-self.H2S**4)
        tpc=self.tpc.R-e
        ppc=self.ppc.psi*tpc/(self.tpc.R+self.H2S*(1-self.H2S)*e)
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def Critical_Carr_Kobayashi_Burrows(self):
        """Carr, N., Kobayashi, R., Burrows, D., 1954. Viscosity of hydrocarbon gases under pres-
        sure. Trans. AIME 201, 270–275."""
        tpc=self.tpc.R-80*self.CO2+130*self.H2S-250*self.N2
        ppc=self.ppc.psi+440*self.CO2+600*self.H2S-170*self.N2
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def Critical_Whitson_Brule(self):
        """Whitson, C. H., and M. R. Brule. Phase Behavior. Richardson, TX: Society of Petroleum Engineers, 2000."""
        CO2=Componente(49)
        H2S=Componente(50)
        N2=Componente(46)
        g=(28.96*self.SG-(N2.M*self.N2+CO2.M*self.CO2+H2S.M*self.H2S))/28.96/(1-self.N2-self.CO2-self.H2S)
        tpcHC=168.+325.*g-12.5*g**2
        ppcHC=677+15.*g-37.5*g**2
        tpc=(1-self.N2-self.CO2-self.H2S)*tpcHC+N2.tc.R*self.N2+CO2.tc.R*self.CO2+H2S.tc.R*self.H2S
        ppc=(1-self.N2-self.CO2-self.H2S)*ppcHC+N2.pc.psi*self.N2+CO2.pc.psi*self.CO2+H2S.pc.psi*self.H2S
        return unidades.Temperature(tpc, "R"), unidades.Pressure(ppc, "psi")

    def tc_gas(self):
        """Cálculo de la temperatura crítica de un gas natural de composición conocida con metano como componente principal, API procedure 4C1.1 pag 329"""
        Tc1=exp(-5.624853)*exp(-0.0105852*self.MABP().R-1.4401126*self.SG+0.013200830*self.SG*self.MABP().R)
        Tc2=self.MABP().R**2.4289880
        Tc3=self.SG**-0.299808
        return unidades.Temperature(Tc1*Tc2*Tc3, "R")

    def Z_factor(self, T, P, Z_method=0):
        """Calculo del factor de compresibilidad del gas
        T: temperatura, en kelvin
        P: presión, en atm
        Z_method: método de cálculo del factor de compresibilidad:
            0   -   Hall_Yarborough
            1   -   Dranchuk_Abu_Kassem
            2   -   Dranchuk_Purvis_Robinson
            3   -   Shell-Oil Company
            4   -   Beggs_Brill
            5   -   Sarem
            6   -   Gopal
            7   -   Papay
        """
        Z=[Z_Hall_Yarborough, Z_Dranchuk_Abu_Kassem, Z_Dranchuk_Purvis_Robinson, Z_ShellOil, Z_Beggs_Brill, Z_Sarem, Z_Gopal, Z_Papay][Z_method]
        return Z(T/self.tpc, P/self.ppc.atm)

    def RhoG(self, T, P):
        Z=self.Z_factor(T, P)
        return unidades.Density(P*self.M/Z/R_atml/T, "gl")

    def Compressibility_Mattar_Brar_Aziz(self, T, P):
        """Mattar, L. G., S. Brar, and K. Aziz. “Compressibility of Natural Gases.” Journal of Canadian Petroleum Technology (October–November 1975): 77–80."""
        Tr=T/self.tpc
        Z=self.Z_factor(T, P)
        g=0.27*P/self.ppc.atm/Tr/Z
        T1=0.31506237-1.0467099/Tr-0.5783272/Tr**3
        T2=0.53530771-0.61232032/Tr
        T3=0.61232032*0.10488813/Tr
        T4=0.68157001/Tr**3
        T5=0.27*Pr/Tr
        dZ=T1+2*T2*g+5*T3*g**4+2*T4*g*(1+0.68446549*g**2-0.68446549**2*g**4)*exp(-0.68446549*g**2)-T5/g
        return self.ppc.atm/P-0.27/Z**2/Tr*dz/(1+g/Z*dz)

    def Gas_Formation_Volume_Factor(self, T, P):
        return Z*T/288.9/P

    def Viscosity_Carr_Kobayashi_Burrows(self, T):
        """Carr, N., R. Kobayashi, and D. Burrows. “Viscosity of Hydrocarbon Gases under Pressure.” Transactions of the AIME 201 (1954): 270–275."""
        muo=8.118e-3-6.15e-3*log10(self.SG)+(1.709e-5-2.062e-6*self.SG)*unidades.Temperature(T, "R").F
        muN2=self.N2*(8.49e-3*log10(self.SG)+9.59e-3)
        muCO2=self.CO2*(9.08e-3*log10(self.SG)+6.24e-3)
        muH2S=self.H2S*(8.49e-3*log10(self.SG)+3.73e-3)
        return unidades.Viscosity(muo+muN2+muCO2+muH2S, "cP")

    def Viscosity_Lee_Gonzalez_Eakin(self, T, P):
        """Lee, A. L., M. H. Gonzalez, and B. E. Eakin. “The Viscosity of Natural Gases.” Journal of Petroleum Technology (August 1966): 997–1000."""
        t=unidades.Temperature(T).R
        K=(9.4+0.02*self.M)*t**1.5/(209+19*self.M+t)
        X=3.5+986/t+0.01*self.M
        Y=2.4-0.2*X
        return unidades.Viscosity(1e-4*K*exp(X*(self.RhoG(T, P).lbft3/62.4)**Y), "cP")



if __name__ == '__main__':
    Crudo(index=1)
