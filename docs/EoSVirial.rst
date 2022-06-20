


Beattie-Bridgeman (1928)
Benedict-Webb-Rubin (1940)
Benedict 1951: fugacity correlations
Cooper-Goldfrank (1967) BWR Coefficients
Orye (1969) BWRO
Morsy (1970) BWRM
Starling (1973) BWRS

Low T modification
Bloomer-Rao (1952)
Motard-Organick (1960)
Barner-Schreiner (1966)
Starling (1971)
Orye (1969)

Critical region
Eubank-Fort (1969)
Hirschfelder (1958)

More terms modifications
Strobridge (1962)
Bender (1970)
Morsy (1970)
Starling (1971)
Jacobsen-Stewart (1973)
Lee-Kesler (1975)
Nishiumi-Saito (1975)
Schmidt-Wagner (1985)
AGA Starling (1991)
Angus (1972, 1976, 1978, 1979 1980, 1985)
Angus-de Reuck (1976)

Starling-Han (1972)
Nishiumi (1980)

Lee-Kesler (1975)
Joffe (1976)
Plocker (1978)
Oellrich (1981)
Yu (1982)
Wu-Stiel (1985)



Yamada
^^^^^^

Yamada create a generalized version of BWR equation with 24 parameters (BWR24)
with the parameters dependence of acentric factor

.. math::
    \begin{align*}
    Z = 1 + \left(B_o-\frac{A_o}{T_R}-\frac{C_o}{T_R}\right)\frac{1}{V_R} \\
    {} + \left(a_\frac{b_R}{T_R}\right)\frac{1}{V_R^2} +
    \frac{\alfa_R}{V_R^5 T_R} \\
    {} + \frac{c_R}{T_R^3 V_R^2}\left(1+\frac{\gamma_R}{V_R^2}\right)
    \exp\left(-\frac{\gamma_R}{V_R^2}\right) \\
    \end{align*}

The equation don't use the critical volume as reducing value for volume, use

.. math::
   V_{sc} = \frac{V_{0.6}}{0.3862-0.0866\omega}

where V_0.6 is the the liquid density at a reduced temperature of 0.6.

The fitting of parameters use data points only for argon, methane, ethane,
n-butane, n-pentane, n-heptane and nitrogen, so the generalization is poor and
no include polar compounds, in fact the author recommend use the correlation
only for compounds with acentric factor not larger than 0.35.

The range of application of equation is only below reduced density of 1.8. The
paper give other correlation with 44 parameters to try to fit a wider range of
densities, but using same databank

.. math::
    \begin{align*}
    Z = 1 + \left(B_o+\frac{B_1}{T_R}+\frac{B_2}{T_R}+\frac{B_3}{T_R^3}\right)
    \frac{1}{V_R} \\
    {} + \left[C_o+\frac{C_1}{T_R}+\frac{C_2}{T_R^3}\left(1+\frac{C_4}{V_R^2}
    \right) \exp \left(-\frac{C_4}{V_R^2\right)\right] \frac{1}{V_R^2} \\
    {} + \left(D_o+\frac{D_1}{T_R}\right)\frac{1}{V_R^3} \\
    {} + \left[E_o+\frac{E_1}{T_R} \left(1+\frac{E_2}{V_R^4}\right)
    \exp \left(-\frac{E_2}{V_R^4}\right)\right] \frac{1}{V_R^4} \\
    {} + \left(F_o+\frac{F_1}{T_R}\right) \frac{1}{V_R^5} \\
    \end{align*}

Yamada, T. An Improved Generalized Equation of State. AIChE J. 19(2) (1973)
286-291, http://dx.doi.org/10.1002/aic.690190212.

Nishiumi (1980)
^^^^^^^^^^^^^^^

Trying to extend the BWRS equation of state to polar compounds Nishiumi
developed this EoS adding 3 additional polar parameters

    
.. math::
    \begin{align*}
    Z = 1 + \left(B_o - \frac{A_o}{T_r} - \frac{C_o}{T_r^3} - \frac{D_o}{Tr^4}
    \frac{E_o+\Psi_E}{T_r^5}\right) \rho_r \\
    {} + \left(b - \frac{a}{T_r} - \frac{d}{Tr_2} - \frac{e}{T_r^5} -
    \frac{f}{T_r^{24}}\right) \rho_r^2 \\
    {} + \alpha \left(\frac{a}{T_r} + \frac{d}{T_r^2} + \frac{e}{T_r^5} +
    \frac{f}{T_r^{24}}\right) \rho_r^5 \\
    {} + \left(\frac{c}{T_r^3} + \frac{g}{T_r^9} + \frac{h}{T_r^{18}} + y(T_r)
    \right) \rho_r^2 \left(1+\gamma\rho_r^2\right)\exp{-\gamma\rho_r^2} \\
    \end{align*}

The aditioanl parameter are:

.. math::
    \begin{array}[t]{l}
    y(T_r) = \frac{s_1}{T_r^{s_2}}\\
    \Psi_E = 2.83 \frac{\mu^2}{T_c V_c}
    \end{array}

This equation is only applicable to pure compound, don't define mixing rules
to apply to mixtures.


Nishiumi, H. An Improved Generalized BWR Equation of State with Three Polar
Parameters Applicable to Polar Substances. J. Chem Eng. Japan 13(3) (1980)
178-183, http://dx.doi.org/10.1252/jcej.13.178.


