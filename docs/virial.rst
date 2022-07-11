
Example code of usage, plot the correlations for R32 and compare with some
sources of experimental values

.. code-block:: python

    from matplotlib import pyplot
    from numpy import linspace, r_

    from lib.mEoS import R32
    from lib.EoS import virial

    D = R32.momentoDipolar.Debye
    Vc = R32.M/R32.rhoc

    # 2nd virial coefficient
    B1 = []
    B2 = []
    B3 = []
    B4 = []
    B5 = []
    B6 = []
    coef = virial.B_Database[R32.id]

    Ti = linspace(250, 450, 50)
    for T in Ti:
        B1.append(virial.B_Tsonopoulos(T, R32.Tc, R32.Pc, R32.f_acent)[0])
        B2.append(virial.B_Meng(T, R32.Tc, R32.Pc, R32.f_acent, D)[0])
        B3.append(virial.B_IglesiasSilva(T, R32.Tc, R32.Pc, Vc, R32.f_acent, D)[0])
        B4.append(virial.B_Orbey(T, R32.Tc, R32.Pc, R32.f_acent, R32.id)[0])
        B5.append(virial.B_Tarakad(T, R32.Tc, R32.Pc, R32.id)[0])
        B6.append(virial._B_Database(T, coef)[0])
    pyplot.plot(Ti, B1, label="Tsonopoulos", ls=":", c="k")
    pyplot.plot(Ti, B2, label="Meng-Duan", ls="--", c="k")
    pyplot.plot(Ti, B3, label="Iglesias-Silva", ls="-.", c="k")
    pyplot.plot(Ti, B4, label="Orbey", ls=":", c="g")
    pyplot.plot(Ti, B5, label="Tarakad", ls=":", c="b")
    pyplot.plot(Ti, B6, label="Database", ls="-", c="k")

    # Experimental date
    # Qian, Z.Y., Nishimura, A., Sato, H., Watanabe, K.
    # Compressibility Factors and Virial Coefficients of Difluoromethane
    # (HFC-32) Determined by Burnett Method
    # JSME Int. J. Ser. B. 36(4) (1993) 665-670
    T = [290, 300, 310, 320, 330, 340, 350, 360, 370]
    B = r_[-0.33975, -0.30666, -0.28011, -0.25594, -0.23379, -0.21422,
           -0.19777, -0.18327, -0.17231]
    pyplot.plot(T, B, ls='', marker="s", mec="k", mfc="k", label="93-Qia/Nis")

    # Sato, T., Sato, H., Watanabe, K.
    # PVT Property Measurements for Difluromethane
    # J. Chem. Eng. Data 39(4) (1994) 851-854
    T = [340, 350, 360, 370, 380, 390, 400, 410, 420]
    B = r_[-207.9, -191.4, -178.2, -166.2, -155.2, -144.7, -135.6, -128.1,
           -119.5]*1e-3
    pyplot.plot(T, B, ls='', marker="d", mec="k", mfc="w", label="94-Sat/Sat")

    pyplot.ylabel("B, l/mol")
    pyplot.xlabel("T, K")
    pyplot.legend()
    pyplot.show()


.. image:: images/Bvirial.png
    :alt: Second virial coefficient for R32


.. code-block:: python

    from matplotlib import pyplot
    from numpy import linspace, r_
    from scipy.constants import R

    from lib.EoS import virial
    from lib.mEoS import R32

    Zc = R32.Pc/R32.rhoc/R*R32.M/R32.Tc/1000
    D = R32.momentoDipolar.Debye

    # 3rd Virial coefficient
    C1 = []
    C2 = []
    C3 = []
    Ti = linspace(250, 450, 50)
    for T in Ti:
        C1.append(virial.C_OrbeyVera(T, R32.Tc, R32.Pc, R32.f_acent)[0]*1e6)
        C2.append(virial.C_LiuXiang(T, R32.Tc, R32.Pc, R32.f_acent, Zc)[0]*1e6)
        B = virial.B_Meng(T, R32.Tc, R32.Pc, R32.f_acent, D)
        C3.append(virial.C_Meng(T, R32.Tc, R32.Pc, D, B)[0]*1e6)

    pyplot.plot(Ti, C1, ls=":", c="k", label="Orbey-Vera")
    pyplot.plot(Ti, C2, ls="--", c="k", label="Liu-Xiang")
    pyplot.plot(Ti, C3, ls="-.", c="k", label="Meng-Duan")

    # Sato, T., Sato, H., Watanabe, K.
    # PVT Property Measurements for Difluromethane
    # J. Chem. Eng. Data 39(4) (1994) 851-854
    T = [340, 350, 360, 370, 380, 390, 400, 410, 420]
    C = r_[0.01625, 0.01431, 0.01325, 0.01226, 0.01133, 0.01034, 0.009646, 0.009418, 0.00848]
    pyplot.plot(T, C, ls='', marker="d", mec="k", mfc="w", label="94-Sat/Sat")

    # Defibaugh, D.R., Morrison, G., Weber, L.A.
    # Thermodynamic Properties of Difluoromethane
    # J. Chem. Eng. Data 39(2) (1994) 333-340
    T = [267, 273, 283, 293, 303, 313, 323, 333, 343, 353, 363, 373]
    C = r_[0.0263, 0.027, 0.0274, 0.0268, 0.0256, 0.0242, 0.0226, 0.0209, 0.0193, 0.0178, 0.0162, 0.0149]
    pyplot.plot(T, C, ls='', marker="s", mec="k", mfc="k", label="94-Def/Mor")

    # Zhang, H., Sato, H., Watanabe, K.
    # Gas Phase PVT Properties for the Difluoromethane + Pentafluoroethane
    # (R32+R125) System
    # J. Chem. Eng. Data 41(6) (1996) 1401-1408
    T = [290, 300, 310, 320, 330, 340, 350, 360, 370]
    C = r_[0.0341, 0.0305, 0.0275, 0.0248, 0.0224, 0.0203, 0.0185, 0.0168, 0.0153]
    pyplot.plot(T, C, ls='', marker="v", mec="k", mfc="w", label="96-Zha/Sat")

    pyplot.ylabel("C, cm⁶/mol²")
    pyplot.xlabel("T, K")
    pyplot.legend()
    pyplot.show()


.. image:: images/Cvirial.png
    :alt: Third virial coefficient for R32



Other correlations don't implemented
------------------------------------


Romero-Lielmezs (1989)
======================

Correlation for second virial coefficient using the Redlich-Kwong cubic EoS.
The correlation need the compound dependent coefficient x.

.. math::
  \begin{array}[t]{l}
    B_r = \frac{BP_c}{RT_c} \\
    B_r = \Omega_b - \Omega_a  \frac{\alpha}{T_r} \\
    \alpha(T_r, \omega, \mu_r^x) = \alpha^{(0)}(T_r) +
    \omega\alpha^{(1)}(T_r) + \mu_r^x\alpha^{(2)}(T_r) \\
    \alpha^{(0)} = -1.4524905 + 14.360017/T_r - 45.000285/T_r^2 +
    14.078304/T_r^6 + 1.7835426/T_r^7 \\
    \alpha^{(1)}(T_r) = -4.3816022 + 15.205023/T_r - 20.874489/T_r^2 +
    12.697209/T_r^3 - 2.5851848/T_r^4 \\
    \alpha^{(2)}(T_r) = -0.008019636/T_r^7 + 0.01092317/T_r^8 -
    0.003505639/T_r^9 \\
  \end{array}

Romero, A., Lielmezs, J. Correlation of Second Virial Coefficient for Polar
Fluids. Thermochimica Acta 145 (1989) 257-264,
http://dx.doi.org/10.1016/0040-6031(89)85145-7.


Besher-Lielmezs (1992)
======================

Correlation for third virial coefficient using the Peng-Robinson cubic EoS. The
correlation need the compound dependent coefficient x.

.. math::
  \begin{array}[t]{l}
    C_r = C \left(\frac{P_c}{RT_c}\right)^2 \\
    C_r = \Omega_b^2 + 2\Omega_a \Omega_b F(T_r, \omega, \mu_r^x) \\
    F(T_r, \omega, \mu_r^x) = F^{(0)}(T_r) + \omega F^{(1)}(T_r) +
    \mu_r^xF^{(2)}(T_r) \\
    F^{(0)}(T_r) = -0.01175 - 0.80483/T_r + 7.29366/T_r^2 -
    16.98304/T_r^3 + 16.86138/T_r^4 - 5.94613/T_r^5 \\
    F^{(1)}(T_r) = 9.25492 - 65.50763/T_r + 162.02620/T_r^2 -
    183.23773/T_r^3 + 96.85253/T_r^4 - 19.08807/T_r^5 \\
    F^{(2)}(T_r) = -1.30215 + 6.56985/T_r - 13.08120/T_r^2 +
    12.85166/T_r^3 - 6.24687/T_r^4 + 1.20561/T_r^5 \\
  \end{array}

Besher, E.M., Lielmezs, J. Correlation for the third virial coefficient for
non-polar and polar compounds using a cubic equation of state. Thermochimica
Acta 200 (1992) 1-13, http://dx.doi.org/10.1016/0040-6031(92)85101-z.


Chueh-Prausnitz (1967)
======================

Correlation for third virial coefficient. The correlation need the compound
dependent coefficient d, give in paper for few compound.

.. math::
    \frac{C}{V_c^2} = \left(0.232T_R^{-0.25} + 0.468T_R^{-5}\right)
    \left(1-e^{1-1.89T_R^2}\right) + d e^{-2.49+2.3T_R-2.7T_R^2}

Chueh, P.L., Prausnitz, J.M. Third Virial Coefficients of Nonpolar Gases and
Their Mixtures. AIChE J. 13(5) (1967) 896-902,
http://dx.doi.org/10.1002/aic.690130516.


de Santis-Grande (1979)
=======================

Correlation for third virial coefficient. The correlation need two aditional
molecular parameters as dipole polarizability and bondi molecular volume. The
paper give these parameters for several compounds but not very much.

.. math::
  \begin{array}[t]{l}
    \frac{C}{v_c^2} = C^0(T_R) + dC´(T_R) + d^2C"(T_R) \\
    C^0 = \frac{0.1961}{T_R^{0.25}} + \frac{0.3972}{T_R^5} +
    \left(0.06684T_R^4 - \frac{0.5428}{T_R^6}\right) e^{-T_R^2} \\
    C´ = \frac{64.5}{T_R^9} \left(1-2.085 e^{-T_R^2}\right) \\
    C" = \frac{801.7}{T_R^7} \\
    d = \frac{\omega \alpha N}{b}
  \end{array}

de Santis, R., Grande, B. An Equation for Predicting Third Virial Coeffcients
of Nonpolar Gases. AIChE J. 25(6) (1979) 931-938,
http://dx.doi.org/10.1002/aic.690250603.
