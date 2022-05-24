"""
Chemistry of aqueous solutions
Coefficients for disequilibrium
(c) Georg Kaufmann
"""

import numpy as np
import sys
import geodyn_chem
# general constants
R = 83.1451 # gas constant [ml/bar/K/mol]

def KW(TC):
    """
    -----------------------------------------------------------------------
    KW - equilibrium constant dissociation of water
    KW: H2O <-> H+ + OH- 
    from:
    Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    refit data of Harned and Owen, The Physical Chemistry of
    Electrolyte Solutions, 1958
    this is on the SWS pH scale in (mol/kg-SW)^2
    input: 
    TC [C]:         temperature
    output:
    KW [mol^2/l^2]: H2O <-> H+ + OH- 
    -----------------------------------------------------------------------
    """
    TK = 273.16 + TC
    KW = 148.9802 - 13847.26/TK - 23.6521*np.log(TK)
    KW = np.exp(KW)
    return KW


def KH(TC):
    """
    -----------------------------------------------------------------------
    KH - Henry constant (solubility of CO2 in water)
    KH: CO2gas <-> CO2water
    from: 
    Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    input:
    TC [C]         - temperature
    output:
    KH [mol/l/atm] - CO2gas <-> CO2water
    -----------------------------------------------------------------------
    """
    TK = 273.16 + TC
    KH =  -60.2409 + (93.4517 / (TK/100)) + (23.3585 * np.log((TK/100)))
    KH = np.exp(KH)
    return KH


def K1K2(TC,S=0):
    """
    -----------------------------------------------------------------------
    K_1 and K_2 equilibrium constants dissociation of CO2 in water
    K1: CO2 + H2O <-> H+ + HCO3- 
    K2: HCO3-     <-> H+ + CO3--
    
    from:
    Millero et al (2006):
    Dissociation constants of carbonic acid in seawater as a function of 
    salinity and temperature, Marine Chemistry 100(1-2):80-94.
    input:
    TC [C]         - temperature
    S [permil]     - salinity
    output:
    K1 [mol/l]     - K1 = (H+)(HCO3-) / (CO2)
    K2 [mol/l]     - K2 = (H+)(CO3--) / (HCO3-)
    -----------------------------------------------------------------------
    """
    TK  = 273.16 + TC
    pK1 = -126.34048 + 6320.813/TK + 19.568224*np.log(TK)
    K1  = 10**(-pK1)
    pK2 = -90.18333 + 5143.692/TK + 14.613358*np.log(TK)
    K2  = 10**(-pK2)
    return K1,K2


def K5(TC):
    """
    -----------------------------------------------------------------------
    K_5 - equilibrium constant dissociation of carbonic acid
    K5: H2CO3 <-> H+ + HCO3-
    -----------------------------------------------------------------------
    input:
    TC [C]         - temperature
    output:
    K5 [mol/l]:    - H2CO3 <-> H+ + HCO3-
    """
    TK  = 273.16 + TC
    K5  = 1.707e-4 *TK/TK
    return K5


def KC(TC,S=0.):
    """
    -----------------------------------------------------------------------
    calcite solubility
    from:
    Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
      sd fit = .01 (for Sal part, not part independent of Sal)
      this is in (mol/kg-SW)^2
    -----------------------------------------------------------------------
    """
    TK  = 273.16 + TC
    KC = -171.9065 - 0.077993*TK + 2839.319/TK
    KC = KC + 71.595*np.log(TK)/np.log(10.0)
    KC = 10.0**KC
    return KC


def ion_debyehueckel(IS,TC):
    '''
    !-----------------------------------------------------------------------
    ! function calculates activity coefficients for different ions
    ! following the Debye-Hueckel model
    ! NOTE for Ca2+, Mg2+, Na+, Cl- the bdot extended model is used
    !      bdot value from phreeqc
    ! input:
    !  TC             - temperature [C]
    !  IS             - ionic strength [mol / l]
    ! output
    !  ion_ca2p       - activity [-]
    !  ion_hco3m      - "
    !  ion_mg2p       - "
    !  ion_hp         - "
    !  ion_co32m      - "
    !  ion_ohm        - "
    !  version using density and dielectric constant
    !  written by Georg Kaufmann 03/01/2008
    !-----------------------------------------------------------------------
    '''
    TK          = 273.160 + TC
    rho         = 1.0
    dielectric  = 87.720 - 0.397020 * TC + 0.000817840 * TC**2
    aa          = 1.82483e6 * np.sqrt(rho / (dielectric*TK)**3.)
    bb          = 50.29120 * np.sqrt(rho / (dielectric*TK))
    ion_hp      = 10.0**(-aa*1.0*np.sqrt(IS)/(1.0+bb*9.00*np.sqrt(IS)))
    ion_ca2p    = 10.0**(-aa*4.0*np.sqrt(IS)/(1.0+bb*5.00*np.sqrt(IS)) + 0.1650*IS)
    ion_mg2p    = 10.0**(-aa*4.0*np.sqrt(IS)/(1.0+bb*5.50*np.sqrt(IS)) + 0.2000*IS)
    ion_ohm     = 10.0**(-aa*1.0*np.sqrt(IS)/(1.0+bb*3.50*np.sqrt(IS)))
    ion_hco3m   = 10.0**(-aa*1.0*np.sqrt(IS)/(1.0+bb*5.40*np.sqrt(IS)))
    ion_co32m   = 10.0**(-aa*4.0*np.sqrt(IS)/(1.0+bb*5.40*np.sqrt(IS)))
    ion_so42m   = 10.0**(-aa*4.0*np.sqrt(IS)/(1.0+bb*5.00*np.sqrt(IS)))
    ion_nap     = 10.0**(-aa*1.0*np.sqrt(IS)/(1.0+bb*4.00*np.sqrt(IS)) + 0.040*IS)
    ion_clm     = 10.0**(-aa*1.0*np.sqrt(IS)/(1.0+bb*3.00*np.sqrt(IS)) + 0.040*IS)
    return ion_hp,ion_ca2p,ion_mg2p,ion_ohm,ion_hco3m,ion_co32m,ion_so42m,ion_nap,ion_clm


def CEQ_limestone_open (TC,pco2):
    '''
    !-----------------------------------------------------------------------
    !  function calculates calcium equilibrium concentration for limestone
    !  and the open system case
    ! input:
    !  K1                  - mol / l
    !  KC                  - mol^2 / l^2
    !  KH                  - mol / l atm
    !  K2                  - mol / l
    !  TC                  - C
    !  pco2                - atm
    ! output:
    !  chem_ceq_limestone_open - mol / l => mol / m^3
    !  written by Georg Kaufmann 03/01/2008
    !-----------------------------------------------------------------------
    '''
    # check for freezing conditions
    if (TC < 0.):
        sys.exit('chem_ceq_limestone_open: T<0')
    # calculate mass balance coeeficients
    K1 = geodyn_chem.K1K2(TC)[0]
    K2 = geodyn_chem.K1K2(TC)[1]
    KH = geodyn_chem.KH(TC)
    KC = geodyn_chem.KC(TC)
    # loop over ionis strength
    strength=1.e-4
    for i in range(1,7):
        [ion_hp,ion_ca2p,ion_mg2p,ion_ohm,ion_hco3m,ion_co32m,ion_so42m,ion_nap,ion_clm]=geodyn_chem.ion_debyehueckel(strength,TC)
        kk = K1*KC*KH / (4.0*K2*ion_ca2p*ion_hco3m**2)
        ceq = (pco2*kk)**(1.0/3.0)
        strength = 3.0*ceq
       #print (k0,k1,k2,k5,kc,ka,kh,kw)
       #print (ion_hp,ion_ca2p,ion_mg2p,ion_ohm,ion_hco3m,ion_co32m,ion_so42m,ion_nap,ion_clm)
    # rescale to mol/m^3
    ceq = 1000.*ceq
    return ceq


def CEQ_limestone_closed (TC,pco2):
    '''
    !-----------------------------------------------------------------------
    !  function calculates calcium equilibrium concentration for limestone
    !  and the closed system case
    ! input:
    !  K1                  - mol / l
    !  KC                  - mol^2 / l^2
    !  KH                  - mol / l atm
    !  K2                  - mol / l
    !  TC                  - C
    !  pco2                - atm
    ! output:
    !  chem_ceq_limestone_closed - mol / l => mol / m^3
    !  written by Georg Kaufmann 03/01/2008
    !-----------------------------------------------------------------------
    '''
    # check for freezing conditions
    if (TC < 0.):
        sys.exit('chem_ceq_limestone_closed: T<0')
    # calculate mass balance coeeficients
    K1 = geodyn_chem.K1K2(TC)[0]
    K2 = geodyn_chem.K1K2(TC)[1]
    KH = geodyn_chem.KH(TC)
    KC = geodyn_chem.KC(TC)
    # loop over ionis strength
    strength=1.e-4
    for i in range(1,7):
        [ion_hp,ion_ca2p,ion_mg2p,ion_ohm,ion_hco3m,ion_co32m,ion_so42m,ion_nap,ion_clm]=geodyn_chem.ion_debyehueckel(strength,TC)
        kk = K1*KC*KH / (4.0*K2*ion_ca2p*ion_hco3m**2)
        a1 = 1.0
        a2 = 0.0
        a3 = kk / KH
        a4 = -kk * pco2
        ceq = geodyn_chem.cubic_root (a1,a2,a3,a4)
        strength = 3.0*ceq
       #print (k0,k1,k2,k5,kc,ka,kh,kw)
       #print (ion_hp,ion_ca2p,ion_mg2p,ion_ohm,ion_hco3m,ion_co32m,ion_so42m,ion_nap,ion_clm)
    # rescale to mol/m^3
    ceq = 1000.*ceq
    return ceq


def cubic_root (a1,a2,a3,a4):
    '''
    !-----------------------------------------------------------------------
    ! find roots of a cubic polynomial
    ! a1*x**3 + a2*x**2 + a3*x + a4 = 0
    ! procedure follows Bronnstein, p. 131ff
    !-----------------------------------------------------------------------
    '''
    import numpy as np
    a  = a2 / a1
    b  = a3 / a1
    c  = a4 / a1
    q = (a**2 - 3.0*b) / 9.0
    r = (2.0*a**3 - 9.0*a*b + 27.0*c) / 54.0
    if (r**2 < q**3):
        phi = np.arccos(r / np.sqrt(q**3))
        x1  = -2.0*np.sqrt(q) * np.cos(phi/3.0) - a/3.0
        x2  = -2.0*np.sqrt(q) * np.cos((phi+2.0*np.pi)/3.0) - a/3.0
        x3  = -2.0*np.sqrt(q) * np.cos((phi-2.0*np.pi)/3.0) - a/3.0
        ceq = x1
    else:
        aa = -(r + np.sqrt(r**2 - q**3))**(1.0/3.0)
        if (aa == 0.):
            bb = 0.0
        elif (aa != 0.):
            bb = q / aa
        ceq = (aa + bb) - a/3.0
    return ceq
