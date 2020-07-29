
Other equations don't implemented


Clausius (1881)
^^^^^^^^^^^^^^^

Historical equation of state

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{T \left(V+c\right)^2}\\
    \Omega_a = 27/64\\
    \Omega_b = Z_c-0.25\\
    \Omega_c = 3/8-Z_c\\
    \end{array}

Clausius, R., “Ueber das Verhaiten der Kohlensaure in Bezug auf Druck,
Volumen und Temperatur”, Ann. Phys. Chem. 9 (1880) 337-359

Wilson (1964)
^^^^^^^^^^^^^

Modified RK temperature dependence

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)}\\
    \Omega_a = \Omega_{ac}\left(1+m\left(\frac{1}{T_r}-1\right)\right)T_r\\
    m = 1.57 + 1.62\omega\\
    \end{array}

Wilson, G.M. Vapor-Liquid Equilibria, Correlations by Means of a Modified
Redlich-Kwong Equation of State. Adv. Cryog. Eng. 9(D2) (1964) 168-176.


HPW (1976)
^^^^^^^^^^

Hederer, Peter, Wenzel equation of state, published like PR in 1976.

.. math::
    P = \frac{RT}{V-b}-\frac{aT^{\alpha}}{V\left(V+b\right)}\\

α, a and b are compound specific parameters, the paper give correlations for
parameter for several homologous series like alkanes, alkenes or alkynes, but
there isn't any general correlation.

Hederer, H., Peter, S., Wenzel, H. Calculation of Thermodynamic Properties from
a Modified Redlich-Kwong Equation of State. Chem. Eng. J., 11 (1976) 183-190,
http://dx.doi.org/10.1016/0300-9467(76)80039-1


vdW Adachi (1984)
^^^^^^^^^^^^^^^^^

Adachi modification to van der Waals original correlation using a logaritmic
temperature dependence for a

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V^2}\\
    a(T) = \frac{27}{64}10^{m\left(1-T_r\right)}\\
    m = 0.228165 + 0.791981\omega - 0.648552*\omega^2 + 0.654505*\omega^3\\
    b = 0.125\frac{RT_c}{P_c}\\
    \end{array}


Valid only by VLE calculation, not applicable in a general use

Adachi, Y., Lu, B.C.-Y. Simplest Equation of State for Vapor-Liquid Equilibrium
calculation: a Modification of the van der Walls Equation. AIChE J. 30(6)
(1984) 991-993, http://dx.doi.org/10.1002/aic.690300619


Usdin-McAuliffe (1976)
^^^^^^^^^^^^^^^^^^^^^^

Incomplete study to improve SRK72 EoS liquid density prediction

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V\left(V+d\right)}\\
    a = \frac{R^2T_c^2}{A_cP_c}\\
    b = \frac{RT_c}{B_cP_c}\\
    d = \frac{RT_c}{D_cP_c}\\
    \alpha^{0.5} = 1 + m\left(1-\sqrt{T_r}\right)\\
    \end{array}

If T_r ≤ 0.7:

.. math::
    m = 0.48049 + 4.516\omega Z_c^* + \left(0.67713\left(\omega-0.35\right)
    -0.02\right)x\left(T_r-0.7\right)

If 0.7 < T_r ≤ 1.0:

.. math::
    m = 0.48049 + 4.516\omega Z_c^* + \left(37.7846\omega{Z_c^*}^3+0.78662
    \right)\left(T_r-0.7\right)^2

:math:`D_c` is the most positive root of

.. math::
    D_c^3 + D_c^2\left(6Z_c-1\right) + D_c\left(4Z_c-1\right)3Z_c +
    \left(8Z_c-3\right)Z_c^2 = 0

:math:`B_c` and :math:`A_c` are calculated known :math:`D_c` from

.. math::
    \begin{array}[t]{l}
    B_c = D_c + 3Z_c - 1\\
    A_c = Z_c^3/B_c
    \end{array}

Usdin, E., McAuliffe, J.C. A One Parameter Family of Equations of State. Chem.
Eng. Sci., 31(ll) (1976) 1077-1084,
http://dx.doi.org/10.1016/0009-2509(76)87030-3


Toghiani-Vismanath (1986)
^^^^^^^^^^^^^^^^^^^^^^^^^

Equation intended to polar substances using the associating paramter of Halm
and Stiel, χ. Not implemented because this parameter is not available in
database.

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a\alpha\left(T_r\right)}{V^2+bV+cbV-cb^2}\\
    a_c = \Omega_a\frac{R^2T_c^2}{P_c}\\
    b_c = \Omega_b\frac{RT_c}{P_c}\\
    c_c = \frac{1-3\zeta_c}{\zeta_c \eta_c}\\
    \Omega_a = \left(1-\zeta_c\left(1-\eta_c\right)\right)^3\\
    \Omega_b = \eta_c \zeta_c\\
    \alpha^{0.5} = 1 + m\left(1-\sqrt{T_r}\right)\\
    \end{array}

The paper define the following generalization for parameters

.. math::
    \begin{array}[t]{l}
    m = A_1 + A_2\omega + A_3\omega^2 + A_4\chi + A_5\chi^2 + A_6\omega \chi +
    A_7\omega_3 + A_8\chi^3\\
    \zeta_c = B_1 + B_2\omega + B_3\omega^2 + B_4\chi + B_5\chi^2 +
    B_6\omega \chi + B_7\omega_3 + B_8\chi^3\\
    \end{array}

with the following parameters

+--------------------+-------------------+
|         m          |   :math:`\zeta_c` |
+====================+===================+
| A₁ = 0.441926073   | B₁ = 0.324020789  |
+--------------------+-------------------+
| A₂ = 1.342755128   | B₂ = -0.056675895 |
+--------------------+-------------------+
| A₃ = -0.328431972  | B₃ = -0.001268996 |
+--------------------+-------------------+
| A₄ = -15.020572758 | B₄ = -3.259131762 |
+--------------------+-------------------+
| A₅ = -2.226936391  | B₅ = 1.399475880  |
+--------------------+-------------------+
| A₆ = 21.112706213  | B₆ = 4.769592801  |
+--------------------+-------------------+
| A₇ = 0.015 079204  | B₇ = 0.003431898  |
+--------------------+-------------------+
| A₈ = 234.965900518 | B₈ = 59.450596850 |
+--------------------+-------------------+

:math:`\eta_c` may be determined using the relation

.. math::
    \eta_c^3 + \eta_c^2\left(\frac{2}{\zeta_c}-3\right) + 3\eta_c - 1 = 0


The polar factor of Halm and Stiel is defined using acentric factor as

.. math::
    \chi = \log{P_r|_{T_r=0.6}} + 1.7\omega + 1.552

so it could be calculated from compound with vapor pressure data available.

Toghiani, H., Viswanath, D.S. A Cubic Equation of State for Polar and Apolar
Fluids. Ind. Eng. Chem. Proc. Des. Dev. 25(2) (1986) 531-536,
http://dx.doi.org/10.1021/i200033a032


Harmens-Knapp (1980)
^^^^^^^^^^^^^^^^^^^^

Developed only with data of normal fluids, alkanes and several nonpolar gases
so not appropiate for general use

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a\alpha}{V^2+cbV-\left(c-1\right)b^2}\\
    a = \Omega_a\frac{R^2T_c^2}{P_c}\\
    b = \Omega_b\frac{RT_c}{P_c}\\
    c = 1 + \frac{1-3\zeta}{\beta \zeta}\\
    \Omega_a = 1-3\zeta+3\zeta^2+\beta\zeta\left(3-6\zeta+\beta\zeta\right)\\
    \Omega_b = \beta \zeta\\
    \beta = 0.10770 + 0.76405\zeta - 1.24282\zeta^2 + 0.96210\zeta^3\\
    \zeta = 0.3211 - 0.08\omega + 0.0384\omega^2\\
    \alpha = \left(1+A\left(1-\sqrt{T_r}\right)-B\left(1-\frac{1}{T_r}\right)
    \right)^2\\
    \end{array}

The paramters of α for ω ≤ 0.2

.. math::
    \begin{array}[t]{l}
    A = 0.5 + 0.27767\omega + 2.17225\omega^2\\
    B = -0.022 + 0.338\omega - 0.845\omega^2\\
    \end{array}

If ω > 0.2

.. math::
    \begin{array}[t]{l}
    A = 0.41311 + 1.14657\omega\\
    B = 0.0118\\
    \end{array}

For temperatures greater than critical the temperature dependence is

.. math::
    \alpha = 1 - \left(0.6258+1.5227\omega\right)\ln{T_r} +
    \left(0.1533+0.41\omega\right)\left(\ln{T_r}\right)^2

Harmens, A., Knapp, H. Three-Parameter Cubic Equation of State for Normal
Substances. Ind. Eng. Chem. Fundam. 19(3) (1980) 291-294,
http://dx.doi.org/10.1021/i160075a010


Fuller (1976)
^^^^^^^^^^^^^

Fuller modified the SRK to improve liquid densities prediction accuracy.

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{V\left(V+cb\right)}\\
    a = \frac{\Omega_a R^2T_c^2 \alpha}{P_c}\\
    b = \Omega_b \frac{RT_c}{P_c}\\
    c = \frac{1}{\beta}\\
    \left(\sqrt{\frac{1}{\beta}-\frac{3}{4}}-\frac{3}{2}\right)\\
    \Omega_a = \frac{\left(1+c\beta\right)^2\Omega_b}
    {\beta\left(1-\beta\right)^2\left(2+c\beta\right)}\\
    \Omega_b = \beta \frac{\left(1-\beta\right)\left(2+c\beta\right) -
    \left(1+c\beta\right)}{\left(2+c\beta\right)\left(1-\beta\right)^2}\\
    \alpha^{0.5} = 1 + q\left(1-Tr^{0.5}\right)\\
    q = \left(\frac{\beta}{0.26}\right)^{1/4}m\\
    m = 0.48 + 1.574\omega - 0.176\omega^2\\
    \beta = \beta_c + \left(\beta_o-\beta_c\right)\left(\frac{2}
    {1+e^{\theta\left(T_r-1\right)}}-1\right)\\
    \frac{\beta_o}{\beta_c} = 7.788 - 36.8316Z_c + 50.7061Z_c^2\\
    \theta = 10.9356+0.0285\bar{P}\\
    \end{array}

βc can be calculated from critical properties:

.. math::
    Z_c = \frac{P_cV_c}{RT_c} = \frac{\left(1-\beta_c\right)\left(2+c_c
    \beta_c\right)-\left(1+c_c\beta_c\right)}
    {\left(2+c_c\beta_c\right)\left(1-\beta_c\right)^2}

This equation isn't implemented because need the parachor parameter for earch
compound, not available in database and calculable by any group contribution
method.

Fuller, G.G. A Modified Redlich-Kwong-Soave Equation of State Capable of
Representing the Liquid State. Ind. Eng. Chem. Fundam. 15(4) (1976) 254-257,
http://dx.doi.org/10.1021/i160060a005


Mathias (1983)
^^^^^^^^^^^^^^

Modification of temperature dependence of α for SRK

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V\left(V+b\right)}\\
    a(T) = 0.42747\frac{R^2T_c^2}{P_c}\alpha\\
    b = 0.08664\frac{RT_c}{P_c}\\
    \alpha = 1 + m\left(1-T_r^0.5\right) - p\left(1-T_r\right)\left(0.7-T_r
    \right)\\
    m = 0.48508 + 1.55191\omega - 0.15613\omega^2\\
    \end{array}

This method need a compound specific parameter, p, a polar parameter.

Mathias, P.M. A Versatile Phase Equilibrium Equation of State. Ind. Eng. Chem.
Process Des. Dev. 22(3) (1983) 385-391, http://dx.doi.org/10.1021/i200022a008.


Martin (1979)
^^^^^^^^^^^^^

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{\left(V+c\right)^2}\\
    a(T) = \frac{27}{64} \frac{R^2T_c^2}{P_c}T_r^{-n}\\
    b = \left(0.857Z_c-0.1674\right)\frac{RT_c}{P_c}\\
    c = \left(-0.857Z_c+0.2924\right)\frac{RT_c}{P_c}\\
    \end{array}

This equation require a compound specific parameters, the temperature exponent
n, determined by equating the slope of the critical isochore to the slope of
the vapour pressure curve at the critical point.

Martin, J.J. Cubic Equations of State-Which?. Ind. Eng. Chem. Fundam. 18(2)
(1979) 81-97, http://dx.doi.org/10.1021/i160070a001.


Rogalski (1990)
---------------

Modified version of the volume-corrected Peng-Robinson equation state

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{\left(V+4.82843b\right)}\\
    \end{array}

This method has several alpha temperature dependence calculation by chemical
type of compounds, furthermore the calculation of pseudovolume is calculated
with a group contribution method.

Rogalski, M., Carrier, B., Solimando, R., Péneloux, A. Correlation and
Prediction of Physical Properties of Hydrocarbons with the Modified
Peng-Robinson Equation of State. 2. Representation of the Vapor Pressures and
of the Molar Volumes. Ind. Eng. Chem. Res. 29(4) (1990) 659-666,
http://dx.doi.org/10.1021/ie00100a026.


Raimondi (1980)
^^^^^^^^^^^^^^^

Modification of temperature dependence of α for SRK

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a(T)}{V\left(V+b\right)}\\
    a(T) = 0.42747\frac{R^2T_c^2}{P_c}\alpha\\
    b = 0.08664\frac{RT_c}{P_c}\\
    \alpha = 1 + \mu(\omega)\left(1-T_r^0.5\right) +
    \delta\left(T_r\right) \nu(\omega, \chi) \left(1-\sqrt{T_r/0.7}\right)^2\\
    \mu(\omega) = m_o + m_1\omega + m_2\omega^2 + m_3\omega^3 + m_4\omega^4 +
    m_5\omega^5\\
    \delta\left(T_r\right) = 1 for T_r < 0.7\\
    \delta\left(T_r\right) = 0 for T_r ≥ 0.7\\
    \nu(\omega, \chi) = \nu_0(\omega) + \chi(\omega)\\
    \nu_0(\omega) = r_0^0 + r_1^0\omega + r_2^0\omega^2 + r_3^0\omega^3\\ 
    \nu_1(\omega) = r_0^1 + r_1^1\omega + r_2^1\omega^2 + r_3^1\omega^3\\ 
    \end{array}

where mi and ri are fixed numerical coefficient originate from regression
analysis.

+------------------+------------------+-------------------+
|         m        |   :math:`r^0`    |   :math:`r^1`     |
+==================+==================+===================+
| m₀ = 0.99717930  | r₀ = 0.033445    | r₀ = 4.6970402    |
+------------------+------------------+-------------------+
| m₁ = 3.39879325  | r₁ = -1.2353988  | r₁ = -0.38872809  |
+------------------+------------------+-------------------+
| m₂ = -0.00727715 | r₂ = -0.36104233 | r₂ = 0.20135167   |
+------------------+------------------+-------------------+
| m₃ = -0.03785315 | r₃ = 0.16525837  | r₃ = -0.068628773 |
+------------------+------------------+-------------------+
| m₄ = -0.03426992 |                  |                   |
+------------------+------------------+-------------------+
| m₅ = 0.05631978  |                  |                   |
+------------------+------------------+-------------------+


This method isn't implemented, because it needs another compound specific
parameter, χ, it represents a second parameter to improve the experimental
vapour pressure fitting, like the acentric factor but defined at Tr=0.5.

Raimondi, L., A Modified Redlich-Kwong Equation of State for Vapour-Liquid
Equilibrium Calculations. Chem. Eng. Sci. 35(6) (1980) 1269-1275,
http://dx.doi.org/10.1016/0009-2509(80)85119-0


Gibbons-Laughton (1984)
^^^^^^^^^^^^^^^^^^^^^^^

Modification of temperature dependence of α for SRK

.. math::
    \alpha = 1 + X\left(T_r-1\right) + Y\left(\sqrt{T_r}-1\right)\\

X and Y are compound specific properties chosen by minimising the error in
the complete vapour pressure curve.

Gibbons, R.M., Laughton, A.P. An Equation of State for Polar and Non-polar
Substances and Mixtures. J. Chem. Soc., Faraday Trans. 2 80(9) (1984)
1019-1038, http://dx.doi.org/10.1039/F29848001019.


Ishikawa-Chung-Lu (1980)
^^^^^^^^^^^^^^^^^^^^^^^^

Mixture of hard sphere model with RK atractive term.

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V}\frac{\left(2V+b\right)}{\left(2V-b\right)} -
    \frac{a(T)}{T^{0.5}V\left(V+b\right)}\\
    a(T) = \Omega_a\frac{R^2T_c^2}{P_c}\alpha\\
    b = \Omega_b\frac{RT_c}{P_c}\\
    \Omega_a = \sum_{i} a_iT_r^i\\
    \Omega_b = \sum_{i} b_iT_r^i\\
    m = 0.48508 + 1.55191\omega - 0.15613\omega^2\\
    \end{array}

The pure component parameter are only available for 22 compound in paper, eight
for each compound for a accuracy not better than other equations with less
parameters, furthermore only with parameter available for alkanes and several
inorganic gases.

Ishikawa, T., Chung, W.K., Lu, B.C.-Y. A Cubic Perturbed, Hard Sphere Equation
of State for Thermodynamic Properties and Vapor-Liquid Equilibrium
Calculations. AIChE J. 26(3) (1980) 372-378,
https://doi.org/10.1002/aic.690260307.


vdW711 (1989)
^^^^^^^^^^^^^

van der Waals volume translation equation with a modified α temperature
dependence.

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V+t-b}-\frac{a_c\alpha}{\left(V+t\right)^2}\\
    a_c = \frac{27}{64}\frac{R^2T_c^2}{P_c}\\
    b = \frac{1}{8}\frac{RT_c}{P_c}\\
    t = t_o + \left(t_c-t_o\right)\exp{\beta \left(1-T_r\right)}\\
    t_o = \frac{RT_c}{P_c}
    \left(0.03901+0.0.04451\omega-0.02274\omega^2\right)\\
    t_c = \frac{RT_c}{P_c}\left(\frac{3}{8}-Z_c\right)\\
    Z_c = 0.2890 - 0.0701\omega - 0.0207\omega^2\\
    \beta = -7.35356-24.5176\omega+9.19829\omega^2\\
    \alpha = \left(1+m\left(1-\sqrt{T_r}\right)\right)^2\\
    m = 0.48553 + 1.62400\omega - 0.21884\omega^2\\
    \end{array}

Androulakis defined a enhanced α temperature dependence with three compound
specific parameters:

.. math::
    \alpha = 1 + d_1\left(1-T_r^{2/3}\right) + d_2\left(1-T_r^{2/3}\right)^2
    + d_3\left(1-T_r^{2/3}\right)^3

Watson, P., Cascella, M., May, D., Salerno, S., Tassios, D. Prediction of Vapor
Pressure and Saturated Molar Volumes with a Simple Cubic Equation of State:
Part II: The van der Waals- 711 EOS. Fluid Phase Equilibria, 27 (1986) 35-52,
http://doi.org/10.1016/0378-3812(86)87039-x

Androulakis, I.P., Kalospiros, N.S., Tassios, D.P. Thermophysical Properties of
Pure Polar and Nonpolar Compounds with a Modified vdW-711 Equation of State.
Fluid Phase Equilibria, 45 (1989) 135-163,
http://doi.org/10.1016/0378-3812(89)80254-7.


Schmidt-Wenzel (1980)
^^^^^^^^^^^^^^^^^^^^^

Generalized form of van der Waals cubic equation of state.


.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b}-\frac{a}{V^2+ubV+wb^2}\\
    u = 1-w\\
    w = -3\omega\\
    a = \Omega_a\alpha\frac{R^2T_c^2}{P_c}\\
    b = \Omega_b\frac{RT_c}{P_c}\\
    \Omega_a = \left(-\xi_c\left(1-\beta_c\right)\right)^3\\
    \Omega_b = \beta_c\xi_c\\
    \alpha = 1+k\left(1-T_r^{0.5}\right)\\
    k = k_o + \frac{\left(5T_r-3k_o-1\right)^2}{70}\\
    k_o = 0.465 + 1.347\omega - 0.528\omega^2\\
    \end{array}

βc can be calculated solving the cubic equation:

.. math::
   \left(6\omega+1\right)\beta_c^3+3\beta_c^2+3\beta_c-1=0
     

Schmidt, G., Wenzel, H. A Modified van der Waals Type Equation of State. Chem.
Eng. Sci. 35(7) (1980) 1503-1512, http://doi.org/10.1016/0009-2509(80)80044-3.


Lee-Erbar-Edmister (1973)
^^^^^^^^^^^^^^^^^^^^^^^^^

Modified version of Grayson-Streed-Chao-Seader equation with a modified
non-cubic equation for the gas phase, this EoS use interaction parameters for
mixtures.

.. math::
    P = \frac{RT}{V-b} - \frac{a}{V\left(V-b\right)} +
    \frac{bc}{V\left(V-b\right)\left(V+b\right)} {}\\

with all three parameters as functions of reduced temperature, critical
temperature, critical pressure and acentric factor:

.. math::
    \begin{array}[t]{l}
    a_i = \frac{R^2T_{ci}^2}{P_{ci}} \left[\left(0.246105+0.02869\omega_i
    \right) - \left(0.037472+0.149687\omega_i\right)T_{ri} \\
    {} + \frac{\left(0.16406+0.023727\omega_i\right)}{T_{ri}} +
    \frac{\left(0.04937+0.132433\omega_i\right)}{T_{ri}^2}\right]\\
    b_i = \frac{RT_{ci}}{P_{ci}} \left(0.086313+0.002\omega_i\right)\\
    c_i = \frac{R^2T_{ci}^2}{P_{ci}} \left[\frac{\left(0.451169+0.00948\omega_i
    \right)}{\sqrt{T_{ri}}}
    \frac{\left(0.387082+0.078842\omega_i\right)}{T_{ri}^2}\right]\\
    \end{array}

and a different correlation for liquid fugacity coefficient

.. math::
    \begin{align*}
    \ln \nu_i = A_1 + \frac{A_2}{T_r} + A_3\ln T_r + A_4T_r + A_5T_r^2 +
    A_6T_r^7\\
    {} + \left(A_7 + \frac{A_8}{T_r} + A_9\ln T_r + A_{10}T_r^2 +
    A_{11}T_r^7\right)P_r\\
    {} + A_{12}T_r^3P_r^2 + \left[\left(1-T_r\right)\left(A_{13}+\frac{A_{14}}
    {T_r} + A_{15}T_r\right)\\
    {} + A_{16}\frac{P_r}{T_r} + A_{17}T_rP_r^2\right]\omega - \ln P_r\\
    \end{align*}

17 coeeficient, with many set of values for different compounds.

The activity coefficient for liquid phase use a modified Scatchard-Hildebrand
versión with as many as 4 interaction parameters.

.. math::
    \begin{array}[t]{l}
    \ln{\gamma_i} = \frac{V_i^L}{RT}\left(\sum_j B_{ij}\Phi_j - \frac{1}{2}
    \sum_j \sum_m B_{jm}\Phi_j\Phi_m\right)\\
    B_{ij} = \left(\delta_i-\delta_j\right)^2 + 2l_{ij}\delta_i\delta_j\\
    \end{array}

Its little improved performance does not justify the increased complexity of the
correlation.

Lee, B-I, Erbar, J.H., Edmister, W.C. Prediction of Thermodynamnic Properties
for Low Temperature Hydrocarbon Process Calculations. AIChE J. 19(2) (1973)
349-356, http://doi.org/10.1002/aic.690190221.


Robinson-Chao (1971)
^^^^^^^^^^^^^^^^^^^^

Modified version of Grayson-Streed-Chao-Seader equation using a Redlich-Kwong
cubic equation modified by Chueh-Prausnitz

.. math::
    \begin{array}[t]{l}
    P = \frac{RT}{V-b} - \frac{a}{T^{0.5} V\left(V+b\right)}\\
    a_i = \frac{\Omega_a R^2 T_{ci}^{2.5}}{P_{ci}}\\
    b_i = \frac{\Omega_b R T_{ci}}{P_{ci}}\\
    a = \sum_i \sum_j y_i y_j a_{ij}\\
    b = \sum_i y_i b_i\\
    a_{ii} = \frac{\Omega_{ai} R^2 T_{ci}^{2.5}}{P_{ci}}\\
    a_{ij} = \frac{\left(\Omega_{ai}+\omega_{aj}\right) R^2 T_{cij}^{2.5}}
    {P_{cij}}\\
    P_{cij} = \frac{z_{cij}RT_{cij}}{v_{cij}}\\
    v_{cij}^{1/3} = \frac{1}{2}\left(v_{ci}^{1/3}+v_{cj}^{1/3}\right)\\
    z_{cij} = 0.291 - 0.08 \left(\frac{\omega_i+\omega_j}{2}\right)\\
    T_{cij} = \sqrt{T_{ci}T_{cj}}\left(1-k_{ij}\right)\\
    \end{array}

Forthemore using different mixing rules for a parameter, the Ω parameters are
compound specific.

The liquid fugacity has modified expresion

.. math::
    \begin{align*}
    \log \nu = \log nu^o + \omega \log nu^1\\
    \log \nu^o = B_o + B_1P_r + B_2P_r^2-\log P_r\\
    \log \nu^1 = \log \nu_{0.6}^1 + \left(P_r-0.6\right)
    \frac{\partial\log\nu^1}{\partial P_r}
    \end{align*}

B₀, B₁ and B₂ are function of Tr, with different dependences at different
values of Tr.

Robinson, R.L, Chao, K.-C. A Correlation of Vaporization Equilibrium Ratios
for Gas Processing Systems. Ind. Eng. Chem. Process Des. Develop. 10(2) (1971)
221-229, http://doi.org/10.1021/i260038a015.
Chueh, P.L., Prausnitz, J.M. Vapor-Liquid Equilibria at Hith Pressures. Vapor-
Phase Fugacity Coefficients in Nonpolar and Quantum-Gas Mixtures. Ind. Eng.
Chem. Fundam. 6(4) (1967) 492-498, http://doi.org/10.1021/i160024a003.

