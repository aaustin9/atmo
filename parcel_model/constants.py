""" Commonly used constants in microphysics and aerosol thermodynamics equations.

================= ============= ========== ==========        ======================
Symbol            Variable      Value      Units             Description
================= ============= ========== ==========        ======================
:math:`g`         ``g``         9.8        m s**-2           gravitational constant
:math:`C_p`       ``Cp``        1004.0     J/kg              specific heat of dry air
                                                             at constant pressure
:math:`\\rho_w`    ``rho_w``     1000.0     kg m**-3          density of water at STP
:math:`R_d`       ``Rd``        287.0      J/kg/K            gas constant for dry air
:math:`R_v`       ``Rv``        461.5      J/kg/K            gas constant for water vapor
:math:`R`         ``R``         8.314      J/mol/K           universal gas constant
:math:`M_w`       ``Mw``        0.0180153  kg/mol            molecular weight of water
:math:`M_a`       ``Ma``        0.0289     kg/mol            molecular weight of dry air
:math:`D_v`       ``Dv``        3e-5       m**2/s            diffusivity of water vapor
                                                             in air
:math:`\\alpha_c`  ``ac``        1.0        unitless          condensation coefficient
:math:`K_a`       ``Ka``        0.02       J/m/s/K           thermal conductivity of air
:math:`a_T`       ``at``        0.96       unitless          thermal accommodation
                                                             coefficient
:math:`\epsilon`  ``epsilon``   0.622      unitless          ratio of :math:`M_w/M_a`
================= ============= ========== ==========        ======================

Additionally, a reference table containing the
`1976 US Standard Atmosphere <http://www.pdas.com/atmos.html>`_ is implemented in the
constant ``std_atm``, which is a pandas DataFrame with the fields

- ``alt``, altitude in km
- ``sigma``, ratio of density to sea-level density
- ``delta``, ratio of pressure to sea-level pressure
- ``theta``, ratio of temperature to sea-level temperature
- ``temp``, temperature in K
- ``press``, pressure in Pa
- ``dens``, air density in kg/m**3
- ``k.visc``, air kinematic viscosity
- ``ratio``, ratio of speed of sound to kinematic viscosity in m**-1

Using default pandas functons, you can interpolate to any reference pressure or
height level.

"""

g = 9.81             #: Gravitational constant, m/s^2
Cp = 1004.0          #: Specific heat of dry air at constant pressure, J/kg
L = 2.5e6            #: Latent heat of condensation, J/kg
rho_w = 1e3          #: Density of water, kg/m^3
Rd = 287.0           #: Gas constant for dry air, J/(kg K)
Rv = 461.5           #: Gas constant for water vapor, J/(kg K)
R = 8.314            #: Universal gas constant, J/(mol K)
Mw = 18.0153/1e3     #: Molecular weight of water, kg/mol
Ma = 28.9/1e3        #: Molecular weight of dry air, kg/mol
Dv = 3.e-5           #: Diffusivity of water vapor in air, m^2/s
ac = 1.0             #: condensation constant
Ka = 2.e-2           #: Thermal conductivity of air, J/m/s/K
at = 0.96            #: thermal accomodation coefficient
epsilon = 0.622      #: molecular weight of water / molecular weight of dry air

import pandas as pd
from StringIO import StringIO
std_atm = """\
 alt  sigma  delta  theta  temp  press  dens   a    visc  k.visc ratio
-0.5 1.0489 1.0607 1.0113 291.4 107477 1.285 342.2 18.05 1.40E-5 24.36
 0.0 1.0000 1.0000 1.0000 288.1 101325 1.225 340.3 17.89 1.46E-5 23.30
 0.5 0.9529 0.9421 0.9887 284.9  95461 1.167 338.4 17.74 1.52E-5 22.27
 1.0 0.9075 0.8870 0.9774 281.7  89876 1.112 336.4 17.58 1.58E-5 21.28
 1.5 0.8638 0.8345 0.9662 278.4  84559 1.058 334.5 17.42 1.65E-5 20.32
 2.0 0.8217 0.7846 0.9549 275.2  79501 1.007 332.5 17.26 1.71E-5 19.39
 2.5 0.7812 0.7372 0.9436 271.9  74691 0.957 330.6 17.10 1.79E-5 18.50
 3.0 0.7422 0.6920 0.9324 268.7  70121 0.909 328.6 16.94 1.86E-5 17.64
 3.5 0.7048 0.6492 0.9211 265.4  65780 0.863 326.6 16.78 1.94E-5 16.81
 4.0 0.6689 0.6085 0.9098 262.2  61660 0.819 324.6 16.61 2.03E-5 16.01
 4.5 0.6343 0.5700 0.8986 258.9  57752 0.777 322.6 16.45 2.12E-5 15.24
 5.0 0.6012 0.5334 0.8873 255.7  54048 0.736 320.5 16.28 2.21E-5 14.50
 5.5 0.5694 0.4988 0.8760 252.4  50539 0.697 318.5 16.12 2.31E-5 13.78
 6.0 0.5389 0.4660 0.8648 249.2  47217 0.660 316.5 15.95 2.42E-5 13.10
 6.5 0.5096 0.4350 0.8535 245.9  44075 0.624 314.4 15.78 2.53E-5 12.44
 7.0 0.4816 0.4057 0.8423 242.7  41105 0.590 312.3 15.61 2.65E-5 11.80
 7.5 0.4548 0.3780 0.8310 239.5  38299 0.557 310.2 15.44 2.77E-5 11.19
 8.0 0.4292 0.3519 0.8198 236.2  35651 0.526 308.1 15.27 2.90E-5 10.61
 8.5 0.4047 0.3272 0.8085 233.0  33154 0.496 306.0 15.10 3.05E-5 10.05
 9.0 0.3813 0.3040 0.7973 229.7  30800 0.467 303.8 14.93 3.20E-5  9.51
 9.5 0.3589 0.2821 0.7860 226.5  28584 0.440 301.7 14.75 3.36E-5  8.99
10.0 0.3376 0.2615 0.7748 223.3  26499 0.414 299.5 14.58 3.53E-5  8.50
10.5 0.3172 0.2422 0.7635 220.0  24540 0.389 297.4 14.40 3.71E-5  8.02
11.0 0.2978 0.2240 0.7523 216.8  22699 0.365 295.2 14.22 3.90E-5  7.57
11.5 0.2755 0.2071 0.7519 216.6  20984 0.337 295.1 14.22 4.21E-5  7.00
12.0 0.2546 0.1915 0.7519 216.6  19399 0.312 295.1 14.22 4.56E-5  6.47
12.5 0.2354 0.1770 0.7519 216.6  17933 0.288 295.1 14.22 4.93E-5  5.99
13.0 0.2176 0.1636 0.7519 216.6  16579 0.267 295.1 14.22 5.33E-5  5.53
13.5 0.2012 0.1513 0.7519 216.6  15327 0.246 295.1 14.22 5.77E-5  5.12
14.0 0.1860 0.1398 0.7519 216.6  14170 0.228 295.1 14.22 6.24E-5  4.73
14.5 0.1720 0.1293 0.7519 216.6  13100 0.211 295.1 14.22 6.75E-5  4.37
15.0 0.1590 0.1195 0.7519 216.6  12111 0.195 295.1 14.22 7.30E-5  4.04
15.5 0.1470 0.1105 0.7519 216.6  11197 0.180 295.1 14.22 7.90E-5  3.74
16.0 0.1359 0.1022 0.7519 216.6  10352 0.166 295.1 14.22 8.54E-5  3.46
16.5 0.1256 0.0945 0.7519 216.6   9571 0.154 295.1 14.22 9.24E-5  3.19
17.0 0.1162 0.0873 0.7519 216.6   8849 0.142 295.1 14.22 9.99E-5  2.95
17.5 0.1074 0.0808 0.7519 216.6   8182 0.132 295.1 14.22 1.08E-4  2.73
18.0 0.0993 0.0747 0.7519 216.6   7565 0.122 295.1 14.22 1.17E-4  2.52
18.5 0.0918 0.0690 0.7519 216.6   6994 0.112 295.1 14.22 1.26E-4  2.33
19.0 0.0849 0.0638 0.7519 216.6   6467 0.104 295.1 14.22 1.37E-4  2.16
19.5 0.0785 0.0590 0.7519 216.6   5979 0.096 295.1 14.22 1.48E-4  2.00
20.0 0.0726 0.0546 0.7519 216.6   5529 0.089 295.1 14.22 1.60E-4  1.85
"""
std_atm = pd.read_csv(StringIO(std_atm), delim_whitespace=True)