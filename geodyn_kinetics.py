"""
Chemistry of aqueous solutions
Coefficients for disequilibrium
(c) Georg Kaufmann
"""
import numpy as np

R        = 8.314   # Gas constant, Pa m3 / K / mol
atm2Pa   = 101325. # 1 atm = 101325 Pa
liter2m3 = 1.e-3   # 1 l   = 1000 cm3 = 10^-3 m3 
year2sec = 365*24*60*60
day2sec  = 24*60*60
min2sec  = 60
tiny     = 1e-20


def k1pk1m(TC,S=0):
    """
    -----------------------------------------------------------------------
    k1-,k1-  - dissociation of CO2 in water
    CO2 + H2O <-> H+ + HCO3-
    -----------------------------------------------------------------------
    Johnsson (1982), after Schulz et al (2006): Determination of the rate 
    constants for the carbon dioxide to bicarbonate inter-conversion in 
    pH-buffered seawater systems, Mar. Chem., 100, 53-65
    input:
    TC [C]         - temperature
    S [permil]     - salinity
    output:
    k1+ [1/s]      - (CO2)+(H2O) -> (H+) + (HCO3-)
    k1- [mol/l/s]  - (CO2)+(H2O) <- (H+) + (HCO3-)
    k1- calculated from k1p/K1
    """
    TK = 273.16 + TC
    k1p = np.exp(1246.98 - 6.19e4/TK - 183.0*np.log(TK))
    k1m = k1p/K1K2(TC)[0]
    return k1p,k1m


def k4pk4m(TC,S=0):
    """
    -----------------------------------------------------------------------
    k4-,k4-  - dissolution of CO2 in OH-
    HCO3- <-> H+ + CO3--
    -----------------------------------------------------------------------
    Johnsson (1982), after Schulz et al (2006): Determination of the rate 
    constants for the carbon dioxide to bicarbonate inter-conversion in 
    pH-buffered seawater systems, Mar. Chem., 100, 53-65 
    input:
    TC [C]         - temperature
    S [permil]     - salinity
    output:
    k4+ [l/mol/s]  - (CO2)+(OH-) -> (HCO3-)
    k4- [1/s]      - (CO2)+(OH-) <- (HCO3-)
    k4- calculated from k1p/K1
    alternative for k4- from equilibrium: k4- = k4+*KW/K1
    """
    TK = 273.16 + TC
    A = -930.13; B = 0.110; D = 3.10e4; E = 140.9
    k4p = A + B*S*0.5 + D/TK + E*np.log(TK)
    k4p = np.exp(k4p) / KW(TC)
    A = -2225.22; B = -0.049; D = 8.91e4; E = 336.6
    k4m = A + B*S*0.5 + D/TK + E*np.log(TK)
    k4m = np.exp(k4m)
    return k4p,k4m


def k5pk5m(TC,S=0):
    """
    -----------------------------------------------------------------------
    k5-,k5-  - dissociation of HCO3- in water
    HCO3- <-> H+ + CO3--
    -----------------------------------------------------------------------
    Eigen (1964), after Schulz et al (2006): Determination of the rate 
    constants for the carbon dioxide to bicarbonate inter-conversion in 
    pH-buffered seawater systems, Mar. Chem., 100, 53-65 
    input:
    TC [C]         - temperature
    S [permil]     - salinity
    output:
    k5+ [l/mol/s]  - (H+)+(CO3--) -> (HCO3-)
    k5- [1/s]      - (H+)+(CO3--) <- (HCO3-)
    k5m = k5p*K2
    """
    TK = 273.16 + TC
    k5p = 5e10*TK/TK
    k5m = k5p*K1K2(TC)[1]
    return k5p,k5m


def k6pk6m(TC,S=0):
    """
    -----------------------------------------------------------------------
    k6-,k6-  - dissociation of HCO3- in OH-
    HCO3- + OH- <-> CO3-- + H2O
    -----------------------------------------------------------------------
    Eigen (1964), after Schulz et al (2006): Determination of the rate 
    constants for the carbon dioxide to bicarbonate inter-conversion in 
    pH-buffered seawater systems, Mar. Chem., 100, 53-65
    input:
    TC [C]         - temperature
    S [permil]     - salinity
    output:
    k6+ [l/mol/s]  - (HCO3-) + (OH-) -> (CO3--) + (H2O)
    k6- [1/s]      - (HCO3-) + (OH-) <- (CO3--) + (H2O)
    k6m calculated from k6p*KW/K2
    """
    TK = 273.16 + TC
    k6p = 6e9*TK/TK
    k6m = k6p*KW(TC)/K1K2(TC)[1]
    return k6p,k6m


def PWP(TC,pco2=0.00042):
    """
    -----------------------------------------------------------------------
    function calculates coefficients of PWP equation
    -----------------------------------------------------------------------
    version tabulated in Buhmann & Dreybrodt (1985)
    input:
    TC               - C
    pco2             - atm
    output:
    kappa1           - m/s
    kappa2           - m/s
    kappa3           - mol/m2/s
    kappa4           - m4/mol/s
    """
    TK = 273.16 + TC
    #-----------------------------------------------------------------------
    # CO_2 equilibrium coefficients
    #-----------------------------------------------------------------------
    kappa1 = 10.**(0.198 - 444./TK)
    kappa2 = 10.**(2.840 - 2177./TK)
    kappa3 = np.where(TC <= 25,
        10.**(-5.860 - 317./TK),
        10.**(-1.100 - 1737./TK))
    kappa4 = np.where(pco2 <= 0.05,
        10.**(-2.375 + 0.025*TC + 0.56*(-np.log10(pco2)-1.3)),
        10.**(-2.375 + 0.025*TC))
    #-----------------------------------------------------------------------
    # convert from original units
    #-----------------------------------------------------------------------
    kappa1 = 1.e-2 * kappa1          # cm/s -> m/s
    kappa2 = 1.e-2 * kappa2          # cm/s -> m/s
    kappa3 = 1.e-3 / 1.e-4 * kappa3  # mmol/cm2/s -> mol/m2/s
    kappa4 = 1.e-8 / 1.e-3 * kappa4  # cm4/mmol/s -> m4/mol/s
    return kappa1,kappa2,kappa3,kappa4
