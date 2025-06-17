import cantera as ct
import numpy as np

def compute_k_eff_cgs(T, P, M_conc_si):
    """
    Compute k_eff in CGS units (e.g., cm^3/mol/s or cm^6/mol^2/s).
    
    Parameters:
        T : Temperature in K
        P : Pressure in Pa
        M_conc_si : Third-body concentration [mol/m^3] (SI units)
    
    Returns:
        k_eff_cgs : Effective rate constant in CGS units
    """
    R = 8.3145  # J/(mol*K)
    
    # Arrhenius params
    A_high = 2.15e12        # 1/(mol*s)
    b_high = 0.292
    Ea_high = -592.0 * 4184  # cal/mol to J/mol
    
    A_low = 2.28e19         # 1/(mol^2*s)
    b_low = -0.977
    Ea_low = 0.0

    # Compute k_inf and k_0 in SI
    k_inf_si = A_high * T**b_high * np.exp(-Ea_high / (R * T))
    k_0_si = A_low * T**b_low * np.exp(-Ea_low / (R * T))
    
    # Reduced pressure
    Pr = (k_0_si * M_conc_si) / k_inf_si

    # Troe parameters
    a = 1.0
    T3 = 1.0e-10
    T1 = 1.0e30
    T2 = 1.0e30

    F_cent = (1 - a) * np.exp(-T / T3) + a * np.exp(-T / T1) + np.exp(-T2 / T)
    logF_cent = np.log10(F_cent)

    c = -0.4 - 0.67 * logF_cent
    n = 0.75 - 1.27 * logF_cent
    f = (np.log10(Pr) + c) / (n - 0.14 * (np.log10(Pr) + c))
    F = 10 ** (logF_cent / (1 + f**2))

    # Final rate in SI
    k_eff_si = (k_0_si * M_conc_si) / (1 + Pr) * F

    # Unit conversion to CGS:
    # mol/m^3 to mol/cm^3 = 1e-6
    # m^3 to cm^3 = 1e6
    # m^6 to cm^6 = 1e12

    # k_eff_si units ≈ m^3/mol/s (depends on regime)
    # So to convert m^3/mol/s to cm^3/mol/s: multiply by 1e6
    k_eff_cgs = k_eff_si * 1e6

    return k_eff_cgs

def CKYTCR(rho, T, y, imw):
    c = np.zeros(12)  # Initialize the concentration array with zeros

    for i in range(12):
        c[i] = rho * y[i] * imw[i]

    return c

gas = ct.Solution('mechanism.yaml')

T = 1600.1843011594888 
P = ct.one_atm*20.

d = {}
for i,name in enumerate(gas.species_names):
    d[name] = 0.0

d['H']    = 7.99153e-04
d['H2']   = 2.06370e-02
d['O']    = 1.45774e-03
d['OH']   = 2.63110e-03
d['H2O']  = 6.16204e-02
d['O2']   = 1.66889e-01
d['N2']   = 7.45220e-01
d['HO2']  = 7.12227e-04
d['H2O2'] = 3.39462e-05
d['OHV']  = 1.31581e-10
d['AR']   = 0.00000e+00
d['HE']   = 0.00000e+00

composition = ",".join(['{}:{}'.format(name,d["{}".format(name)]) for name in gas.species_names])
gas.TPY = T, P, composition

density = gas.density #[kg/m3]
# density *= 1.e-3 # [kg/m3] -> [g/cm3] CGS

massfrac = np.zeros(gas.n_species)
C        = np.zeros(gas.n_species)

for i,name in enumerate(d):
    massfrac[i] = d[name]

mw  = gas.molecular_weights
imw = 1.0/mw

C = CKYTCR(density, T, massfrac, imw)

k_eff = compute_k_eff_cgs(T, P, C[gas.species_index('H2O')])

wdot = k_eff * C[gas.species_index('H')] * C[gas.species_index('O2')] * C[gas.species_index('H2O')]
print("k_eff ", k_eff)
print("wdot ", wdot)
