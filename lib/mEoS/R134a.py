#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class R134a(MEoS):
    """Ecuación de estado de multiparametros para el R134a

    >>> r134A=R134a(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r134A.T, r134A.rho, r134A.u.kJkg, r134A.h.kJkg, r134A.s.kJkgK, r134A.cv.kJkgK, r134A.cp.kJkgK,r134A.w)
    300.0 4.8836 342.90 363.38 1.7177 0.72455 0.79910 149.16
    """
    name="1,1,1,2-tetrafluoroethane"
    CASNumber="811-97-2"
    formula="CF3CH2F"
    synonym="R134a"
    rhoc=unidades.Density(511.9)
    Tc=unidades.Temperature(374.21)
    Pc=unidades.Pressure(4059.28, "kPa")
    M=102.032      #g/mol
    Tt=unidades.Temperature(169.85)
    Tb=unidades.Temperature(247.076)
    f_acent=0.32684
    momentoDipolar=unidades.DipoleMoment(2.058, "Debye")
    id=236
#    id=1235

    CP1={  "ao": -0.629789,
                "an": [7.292937/374.18**0.5, 5.154411/374.18**0.75],
                "pow": [0.5, 0.75],
                "ao_exp": [], "exp": [],
                "ao_hyp": [],"hyp": []}

    CP2={  "ao": 19.4006,
                "an": [0.258531, -1.29665e-4],
                "pow": [1, 2],
                "ao_exp": [], "exp": [],
                "ao_hyp": [],"hyp": []}

    CP3={  "ao": 1.838736,
                "an": [3.01994e-2, -1.78455e-5, 4.42442e-9],
                "pow": [1, 2, 3],
                "ao_exp": [], "exp": [],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Tillner-Roth & Baehr (1994).",
        "__doc__":  u"""Tillner-Roth, R. and Baehr, H.D., "An international standard formulation of the thermodynamic properties of 1,1,1,2-tetrafluoroethane (HFC-134a) for temperatures from 170 K to 455 K at pressures up to 70 MPa," J. Phys. Chem. Ref. Data, 23:657-729, 1994.""",
        "R": 8.314471,
        "Tref": 374.18,
        "rhoref": 4.978830171*M,
        "cp": CP1,

        "nr1": [0.5586817000e-1, 0.4982230000, 0.2458698000e-1, 0.8570145000e-3, 0.4788584000e-3 , -0.1800808000e1 , 0.2671641000, -0.4781652000e-1],
        "d1": [2, 1, 3, 6, 6, 1, 1, 2],
        "t1": [-0.5, 0, 0, 0, 1.5, 1.5, 2, 2],

        "nr2": [0.1423987000e-1, 0.3324062000, -0.7485907000e-2, 0.1017263000e-3, -0.5184567000, -0.8692288000e-1, 0.2057144000, -0.5000457000e-2, 0.4603262000e-3, -0.3497836000e-2, 0.6995038000e-2, -0.1452184000e-1, -0.1285458000e-3 ],
        "d2": [5, 2, 2, 4, 1, 4, 1, 2, 4, 1, 5, 3, 10],
        "t2": [1, 3, 5, 1, 5, 5, 6, 10, 10, 10, 18, 22, 50],
        "c2": [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4],
        "gamma2": [1]*13}

    MBWR={
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-134a of Huber and McLinden (1992)",
        "__doc__":  u"""Huber, M.L. and McLinden, M.O., "Thermodynamic properties of R134a (1,1,1,2-tetrafluoroethane)," International Refrigeration Conference, West Lafayette, IN, July 14-17, 453-462, 1992.""",
        "R": 8.314471,
        "cp": CP2,

        "b": [None, 0.965209362217e-1, -0.401824768889e1, 0.395239532858e2, 0.134532868960e4, -0.139439741347e7, -0.309281355175e-2, 0.292381512283e1, -0.165146613555e4, 0.150706003118e7, 0.534973948313e-4, 0.543933317622, -0.211326049762e3, -0.268191203847e-1, -0.541067125950, -0.851731779398e3, 0.205188253646, -0.733050188093e-2, 0.380655963862e1, -0.105832087589, -0.679243084424e6, -0.126998378601e9, -0.426234431829e5, 0.101973338234e10, -0.186699526782e3, -0.933426323419e5, -0.571735208963e1, -0.176762738787e6, -0.397282752308e-1, 0.143016844796e2, 0.803085294260e-4, -0.171959073552, 0.226238385661e1]}

    helmholtz2={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-134a of Span and Wagner (2003).",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. III. Results for Polar Fluids," Int. J. Thermophys., 24(1):111-162, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [0.106631890000e1, -0.244959700000e1, 0.446457180000e-1, 0.756568840000e-1, 0.206520890000e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.42006912, 0.76739111, 0.17897427e-2, -0.36219746, -0.6780937e-1, -0.10616419, -0.18185791e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz3={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Astina and Sato (2004)",
        "__doc__":  u"""Astina, I.M. and Sato, H., "A Fundamental Equation of State for 1,1,1,2-Tetrafluoroethane with an Intermolecular Potential Energy Background and Relialbe Ideal-Gas Properties," Fluid Phase Equilib., 221:103-111, 2004.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [1.832124209, -2.940698861, 5.156071823e-1, 2.756965911e-1, 1.225264939, -6.486749497e-1, -9.286738053e-1, 3.920381291e-1, 1.056692108e-1],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4],
        "t1": [0.5, 1.125, 3.25, 0.5, 1.875, 2.75, 1.625, 2.125, 1.125],

        "nr2": [-7.586523371e-1, -1.198140136, -2.878260390e-1, -9.723032379e-2, 5.307113358e-2, -4.681610582e-2, -9.604697902e-3, 6.668035048e-3, 2.361266290e-3],
        "d2": [1, 2, 3, 2, 3, 4, 4, 5, 6],
        "t2": [3.75, 1.5, 0.75, 9, 8.5, 5.5, 32, 23, 31],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    eq=helmholtz1, MBWR, helmholtz2, helmholtz3
    _PR=0.001032
    
    _surface={"sigma": [0.06016], "exp": [1.26]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.77513e1, 0.29263e1, -0.26622e1, -0.39711e1], "exp": [1.0, 1.5, 1.9, 4.25]}
    _liquid_Density={ "eq": 1, "ao": [0.12449e2, -0.41023e2, 0.73641e2, -0.64635e2, 0.22551e2], "exp": [0.5, 0.7, 0.9, 1.1, 1.3]}
    _vapor_Density={ "eq": 3, "ao": [-0.29174e1, -0.72542e1, -0.23306e2, 0.59840e1, -0.71821e2], "exp": [0.383, 1.21, 3.3, 5.6, 7.0]}

    visco0={"eq": 1, "omega": 1,
                    "__name__": "Huber (2003)",
                    "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., "Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a", Ind. Eng. Chem. Res., 42:3163-3178, 2003.""",
                    "ek": 299.363, "sigma": 0.468932,
                    "collision": [0.355404, -0.464337, 0.257353e-1],
                    "Tref": 1., "rhoref": 1.*M,
                    "n_chapman": 0.215729/M**0.5,

                    "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4, 0.24710125e4, -0.33751717e4, 0.24916597e4, -0.78726086e3, 0.14085455e2, -0.34664158],
                    "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
                    "Tref_virial": 299.363, "etaref_virial": 0.0620984,

                    "Tref_res": 374.21, "rhoref_res": 5.0170613*M, "etaref_res": 1000,
                    "n_packed": [3.163695635587490, -0.8901733752064137e-1, 0.1000352946668359],
                    "t_packed": [0, 1, 2],
                    "n_poly": [-0.2069007192080741e-1, 0.3560295489828222e-3, 0.2111018162451597e-2, 0.1396014148308975e-1, -0.4564350196734897e-2, -0.3515932745836890e-2, -0.2147633195397038],
                    "t_poly": [0, -6, -2, -0.5, 2, 0, 0],
                    "d_poly": [1, 2, 2, 2, 2, 3, 0],
                    "g_poly": [0, 0, 0, 0, 0, 0, -1],
                    "c_poly": [0, 0, 0, 0, 0, 0, 0],
                    "n_num": [0.2147633195397038],
                    "t_num": [0],
                    "d_num": [0],
                    "g_num": [0],
                    "c_num": [0],
                    "n_den": [1, -1],
                    "t_den": [0, 0],
                    "d_den": [0, 1],
                    "g_den": [1, 0],
                    "c_den": [0, 0]}

    visco1={"eq": 4, "omega": 1,
                    "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
                    "__name__": "Quiñones-Cisneros (2006)",
                    "Tref": 374.21, "muref": 1.0,
                    "ek": 299.363, "sigma": 0.468932, "n_chapman": 0,
                    "n_ideal": [31.2515, -89.6122, 73.0823],
                    "t_ideal": [0, 0.25, 0.5],

                    "a": [1.07271318464787e-4, -4.41655360682255e-5, 0.0],
                    "b": [1.66457266522365e-4, -4.80292908400793e-5, 0.0],
                    "c": [8.08333416284215e-5, -4.90359549823121e-5, 0.0],
                    "A": [-2.12476175599662e-8, 2.81647242085073e-9, 0.0],
                    "B": [1.35593527573090e-8, 0.0, 3.17549774078234e-10],
                    "C": [0.0, 4.81768878752129e-7, -1.17148596093671e-7],
                    "D": [0.0, 0.0, 0.0 ]}

    visco2={"eq": 1, "omega": 1,
                    "__name__": "Laesecke (2003)",
                    "__doc__": """Laesecke, A., "Data reassessment and full surface correlation of the viscosity of HFC-134a (1,1,1,2-tetrafluoroethane),""",
                    "ek": 288.82, "sigma": 0.50647,
                    "collision": [0.355404, -0.464337, 0.257353e-1],
                    "Tref": 1., "rhoref": 1.*M,
                    "n_chapman": 0.215729/M**0.5,

                    "n_virial": [-0.17999496e1, 0.46692621e2, -0.53460794e3, 0.33604074e4, -0.13019164e5, 0.33414230e5, -0.58711743e5, 0.71426686e5, -0.59834012e5, 0.33652741e5, -0.12027350e5, 0.24348205e4, -0.20807957e3],
                    "t_virial": [0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5, -5.5, -6],
                    "Tref_virial": 288.82, "etaref_virial": 0.07823693,

                    "Tref_res": 374.18, "rhoref_res": 4.9788302*M, "etaref_res": 1000,
                    "n_packed": [3.07383, 0.482539055],
                    "t_packed": [0, 1],
                    "n_poly": [-0.331249e-1, -0.468509e-3, 0.306398],
                    "t_poly": [0, 0, 0],
                    "d_poly": [1, 2, 0],
                    "g_poly": [0, 0, -1],
                    "c_poly": [0, 0, 0],
                    "n_num": [-0.306398, 0.215221],
                    "t_num": [0, 0],
                    "d_num": [0, 1],
                    "g_num": [0, 0],
                    "c_num": [0, 0],
                    "n_den": [1, -1],
                    "t_den": [0, 0],
                    "d_den": [0, 1],
                    "g_den": [1, 0],
                    "c_den": [0, 0]}

    _viscosity=visco0, visco1, visco2

    thermo0={"eq": 1,
                "__name__": "Perkins (2000)",
                "__doc__": """Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L.,"Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)," NISTIR, 2000""",

                "Tref": 1., "kref": 1.,
                "no": [-1.05248e-2, 8.00982e-5],
                "co": [0, 1],

                "Trefb": 339.173, "rhorefb": 5.049886, "krefb": 2.055e-3,
                "nb": [1.836526, 5.126143, -1.436883, 6.261441e-1],
                "tb": [0]*4,
                "db": [1, 2, 3, 4],
                "cb": [0]*4,

                "critical": 3,
                "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
                "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.285356e-10, "Tcref": 561.411}

    _thermal=thermo0,

