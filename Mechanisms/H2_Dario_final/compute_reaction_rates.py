import numpy as np
import cantera as ct

def gibbs(species, T):
    T2 = T * T
    T3 = T * T * T
    T4 = T * T * T * T
    invT = 1.0 / T
    logT = np.log(T)

    # Species with midpoint at T=1000 Kelvin
    if T < 1000:
        # species 0: H
        species[0] = (
            +2.547366000000000e+04 * invT + 2.946682850000000e+00
            - 2.500000000000000e+00 * logT - 2.740033625000000e-16 * T
            + 2.665623750000000e-19 * T2 - 1.552516841666667e-22 * T3
            + 3.672898590000000e-26 * T4
        )
        # species 1: H2
        species[1] = (
            -9.511467020000000e+02 * invT + 3.559656969000000e+00
            - 2.721136050000000e+00 * logT - 2.396071225000000e-03 * T
            + 1.691370716666666e-06 * T2 - 7.474474574999999e-10 * T3
            + 1.353315510000000e-13 * T4
        )
        # species 2: O
        species[2] = (
            +2.912770090000000e+04 * invT + 8.052918099999999e-01
            - 3.106527720000000e+00 * logT + 1.378451980000000e-03 * T
            - 8.523930200000001e-07 * T2 + 3.579118225000000e-10 * T3
            - 6.737828450000000e-14 * T4
        )
        # species 3: OH
        species[3] = (
            +3.369554460000000e+03 * invT + 4.058480713800000e+00
            - 3.984540380000000e+00 * logT + 1.169039770000000e-03 * T
            - 7.387209883333333e-07 * T2 + 3.048454441666667e-10 * T3
            - 6.354742700000001e-14 * T4
        )
        # species 4: H2O
        species[4] = (
            -3.029732530000000e+04 * invT + 5.253373550000000e+00
            - 4.239470870000000e+00 * logT + 1.190968670000000e-03 * T
            - 1.255242963333333e-06 * T2 + 5.583659075000000e-10 * T3
            - 1.139007850000000e-13 * T4
        )
        # species 5: O2
        species[5] = (
            -1.052703670000000e+03 * invT - 5.176797700000000e-01
            - 3.654932980000000e+00 * logT + 9.588401700000001e-04 * T
            - 1.114957153333333e-06 * T2 + 4.912463450000000e-10 * T3
            - 8.317090400000000e-14 * T4
        )
        # species 6: N2
        species[6] = (
            -1.051964560000000e+03 * invT + 8.486600300000000e-01
            - 3.587600400000000e+00 * logT + 3.012735500000000e-04 * T
            - 1.497218761666667e-07 * T2 - 6.291014016666667e-11 * T3
            + 3.537338020000000e-14 * T4
        )
        # species 7: HO2
        species[7] = (
            +3.014058970000000e+02 * invT - 1.551889300000000e+00
            - 3.877615220000000e+00 * logT + 5.799244400000000e-04 * T
            - 1.775876083333333e-06 * T2 + 9.734864750000000e-10 * T3
            - 2.017821605000000e-13 * T4
        )
        # species 8: H2O2
        species[8] = (
            -1.766382310000000e+04 * invT - 1.411878630000000e+00
            - 3.828191160000000e+00 * logT - 1.636540200000000e-03 * T
            - 9.304980266666667e-07 * T2 + 6.848160891666667e-10 * T3
            - 1.527466825000000e-13 * T4
        )
        # species 9: OHV
        species[9] = (
            +5.002130060000000e+04 * invT + 2.278371300000000e+00
            - 3.637259270000000e+00 * logT - 9.257252050000000e-05 * T
            + 2.793854533333333e-07 * T2 - 1.989472975000000e-10 * T3
            + 4.216041795000000e-14 * T4
        )
        # species 10: AR
        species[10] = (
            -7.453750000000000e+02 * invT - 1.879674910000000e+00
            - 2.500000000000000e+00 * logT - 2.740033625000000e-16 * T
            + 2.665623750000000e-19 * T2 - 1.552516841666667e-22 * T3
            + 3.672898590000000e-26 * T4
        )
        # species 11: HE
        species[11] = (
            -7.453750000000000e+02 * invT + 1.571276026000000e+00
            - 2.500000000000000e+00 * logT - 2.740033625000000e-16 * T
            + 2.665623750000000e-19 * T2 - 1.552516841666667e-22 * T3
            + 3.672898590000000e-26 * T4
        )
    else:
        # species 0: H
        species[0] = (
            +2.547366000000000e+04 * invT + 2.946682850000000e+00
            - 2.500000000000000e+00 * logT
        )
        # species 1: H2
        species[1] = (
            -8.167484240000000e+02 * invT + 3.961172590000000e+00
            - 2.932865750000000e+00 * logT - 4.133040130000000e-04 * T
            + 2.440039400000000e-08 * T2 - 1.284170116666667e-12 * T3
            + 3.444024000000000e-17 * T4
        )
        # species 2: O
        species[2] = (
            +2.922661550000000e+04 * invT - 2.379309400000000e+00
            - 2.543636970000000e+00 * logT + 1.365812430000000e-05 * T
            + 6.983825333333333e-10 * T2 - 4.129015375000000e-13 * T3
            + 2.397768470000000e-17 * T4
        )
        # species 3: OH
        species[3] = (
            +3.697880840000000e+03 * invT - 3.006494770000000e+00
            - 2.838530330000000e+00 * logT - 5.537064450000000e-04 * T
            + 4.900003483333333e-08 * T2 - 3.505822741666666e-12 * T3
            + 1.211449450000000e-16 * T4
        )
        # species 4: H2O
        species[4] = (
            -2.988629310000000e+04 * invT - 4.205080219999999e+00
            - 2.677038900000000e+00 * logT - 1.486590800000000e-03 * T
            + 1.289614816666667e-07 * T2 - 7.869459500000001e-12 * T3
            + 2.134499550000000e-16 * T4
        )
        # species 5: O2
        species[5] = (
            -1.214730780000000e+03 * invT + 2.442514700000000e-01
            - 3.660960650000000e+00 * logT - 3.281829055000000e-04 * T
            + 2.352493783333333e-08 * T2 - 1.714982791666666e-12 * T3
            + 6.495671800000001e-17 * T4
        )
        # species 6: N2
        species[6] = (
            -9.245018400000000e+02 * invT - 2.918713730000000e+00
            - 2.952576370000000e+00 * logT - 6.984502000000001e-04 * T
            + 8.210526716666667e-08 * T2 - 6.550084958333333e-12 * T3
            + 2.303776020000000e-16 * T4
        )
        # species 7: HO2
        species[7] = (
            +3.516660660000000e+01 * invT + 1.210132220000000e+00
            - 4.172287410000000e+00 * logT - 9.405881350000000e-04 * T
            + 5.771288100000000e-08 * T2 - 1.622146241666667e-12 * T3
            - 8.812845250000000e-18 * T4
        )
        # species 8: H2O2
        species[8] = (
            -1.800241800000000e+04 * invT + 3.909661134000000e+00
            - 4.579773050000000e+00 * logT - 2.026630015000000e-03 * T
            + 2.164078833333333e-07 * T2 - 1.651761666666667e-11 * T3
            + 5.698439600000000e-16 * T4
        )
        # species 9: OHV
        species[9] = (
            +5.030140630000000e+04 * invT - 2.712986050000000e+00
            - 2.882730000000000e+00 * logT - 5.069871500000000e-04 * T
            + 3.794795000000000e-08 * T2 - 1.812235833333334e-12 * T3
            + 2.563152500000000e-17 * T4
        )
        # species 10: AR
        species[10] = (
            -7.453750000000000e+02 * invT - 1.879674910000000e+00
            - 2.500000000000000e+00 * logT
        )
        # species 11: HE
        species[11] = (
            -7.453750000000000e+02 * invT + 1.571276026000000e+00
            - 2.500000000000000e+00 * logT
        )

def production_rate(wdot, sc, g_RT, T, qdot_out):
    invT = 1.0 / T
    logT = np.log(T)

    # Reference concentration: P_atm / (RT) in inverse mol/m^3
    refC = 101325 / 8.31446 * invT
    refCinv = 1 / refC

    # Initialize wdot
    wdot[:] = 0.0
    qdot_out[:] = 0.0

    # Compute the mixture concentration
    mixture = np.sum(sc)

    # Compute the Gibbs free energy
    # g_RT = gibbs(T)  

    # Reaction 18: H2O2 (+M) <=> 2 OH (+M)
    k_f = 2000000000000 * np.exp((0.9) * logT - (24531.3092413143) * invT)
    Corr = mixture + (2.7) * sc[1] + (4.1) * sc[4] + (-0.21) * sc[5] + (4.2) * sc[8] + (-0.32) * sc[10] + (-0.56) * sc[11]
    redP = Corr / k_f * 2.35e+18 * np.exp(-2.293 * logT - 24528.2899413195 * invT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(0.57 * np.exp(-T * 1e+30) + 0.43 * np.exp(-T * 1e-30) + 0.0)
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    qf = Corr * k_f * (sc[8])
    qr = Corr * k_f * np.exp(-(-2.000000 * g_RT[3] + g_RT[8])) * (refCinv) * (sc[3] ** 2)
    qdot = qf - qr
    qdot_out[18] = qdot
    wdot[3] += 2.000000 * qdot
    wdot[8] -= qdot

    # H+O2(+M)=HO2(+M)   1.025E12  0.604  -241.1 
    # HE/0.35/ H2O/0.0/ AR/0.0/ H2/0.0/ HE/0.0/  
    # LOW/2.720e+19 -1.290  0.0/  
    # TROE/0.458 1.0E-10 1.0E30 1.0E30/

    # - equation: H + O2 (+M) <=> HO2 (+M)  # Reaction 31
    # type: falloff
    # low-P-rate-constant: {A: 2.72e+19, b: -1.29, Ea: 0.0}
    # high-P-rate-constant: {A: 1.025e+12, b: 0.604, Ea: -241.1}
    # Troe: {A: 0.458, T3: 1.0e-10, T1: 1.0e+30, T2: 1.0e+30}
    # efficiencies: {HE: 0.0, H2O: 0.0, AR: 0.0, H2: 0.0}

    # Reaction 30: H + O2 (+M) <=> HO2 (+M)
    k_f = 1025000 * np.exp((0.604) * logT - (-121.325538125518) * invT)
    Corr = mixture + (-1) * sc[1] + (-1) * sc[4] + (-1) * sc[10] + (-1) * sc[11]
    redP = Corr / k_f * 27200000 * np.exp(-1.29 * logT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(0.542 * np.exp(-T * 10000000000) + 0.458 * np.exp(-T * 1e-30) + np.exp(-1e+30 * invT))
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    qf = Corr * k_f * (sc[0] * sc[5])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7])
    qdot = qf - qr
    qdot_out[30] = qdot
    wdot[0] -= qdot
    wdot[5] -= qdot
    wdot[7] += qdot
    
    print("Reaction 30: H + O2 (+M) <=> HO2 (+M)")
    print("qf: ",qf/1000.0)
    print("qr: ",qr/1000.0)

    # Reaction 31: H + O2 (+AR) <=> HO2 (+AR)
    k_f = 564000 * np.exp((0.691) * logT - (-359.497986050892) * invT)
    Corr = mixture
    redP = Corr / k_f * 1190000 * np.exp(-1.04 * logT - 0.0553538332385192 * invT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(0.565 * np.exp(-T * 1e+10) + 0.435 * np.exp(-T * 1e-30) + np.exp(-1e+30 * invT))
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    qf = Corr * k_f * (sc[0] * sc[5] * sc[10])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7] * sc[10])
    qdot = qf - qr
    qdot_out[31] = qdot
    wdot[0] -= qdot
    wdot[5] -= qdot
    wdot[7] += qdot

    # Reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)

    # // - equation: H + O2 (+H2O) <=> HO2 (+H2O)  # Reaction 33
    # //   type: falloff
    # //   low-P-rate-constant: {A: 2.28e+19, b: -0.977, Ea: 0.0}
    # //   high-P-rate-constant: {A: 2.15e+12, b: 0.292, Ea: -592.0}
    # //   Troe: {A: 1.0, T3: 1.0e-10, T1: 1.0e+30, T2: 1.0e+30}

    # H+O2(+H2O)=HO2(+H2O)   2.15E12  0.292  -592.0  
    # LOW/2.28e+19 -0.977  0.0/  
    # TROE/1.0 1.0E-10 1.0E30 1.0E30/
    k_f = 2150000*(T**0.292)*np.exp(-297.904266156394/T)
    # k_f = 2150000 * np.exp((0.292) * logT - (-297.904266156394) * invT)
    Corr = sc[4] #mixture
    redP = Corr / k_f * 22800000 * np.exp(-0.977 * logT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(1 * np.exp(-T * 1e-30) + np.exp(-1e+30 * invT))
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    #correct up to here
    qf = Corr * k_f * (sc[0] * sc[5])
    #index 7 is HO2
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (sc[4] * sc[7])
    qdot = qf - qr
    qdot_out[32] = qdot
    wdot[0] -= qdot
    wdot[5] -= qdot
    wdot[7] += qdot
    print("Reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)")
    print("qf: ",qf/1000.0)
    print("qr: ",qr/1000.0)

    # Reaction 33: H + O2 (+H2) <=> HO2 (+H2)
    k_f = 4230000 * np.exp((0.646) * logT - (702.339500463647) * invT)
    Corr = sc[1] #mixture
    redP = Corr / k_f * 828000000 * np.exp(-1.72 * logT - -1.36371716433079 * invT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(0.647 * np.exp(-T * 1e+10) + 0.353 * np.exp(-T * 1e-30) + np.exp(-1e+30 * invT))
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    qf = Corr * k_f * (sc[0] * sc[5])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (sc[1] * sc[7])
    qdot = qf - qr
    qdot_out[33] = qdot
    wdot[0] -= qdot
    wdot[5] -= qdot
    wdot[7] += qdot
    print("Reaction 32: H + O2 (+H2O) <=> HO2 (+H2O)")
    print("qf: ",qf/1000.0)
    print("qr: ",qr/1000.0)

    # Reaction 34: H + O2 (+HE) <=> HO2 (+HE)
    k_f = 1025000 * np.exp((0.604) * logT - (-121.325538125518) * invT)
    Corr = mixture
    redP = Corr / k_f * 27200000 * np.exp(-1.29 * logT)
    F = redP / (1.0 + redP)
    logPred = np.log10(redP)
    logFcent = np.log10(0.542 * np.exp(-T * 1e+10) + 0.458 * np.exp(-T * 1e-30) + np.exp(-1e+30 * invT))
    troe_c = -0.4 - 0.67 * logFcent
    troe_n = 0.75 - 1.27 * logFcent
    troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred))
    F_troe = np.exp(np.log(10) * logFcent / (1.0 + troe * troe))
    Corr = F * F_troe
    qf = Corr * k_f * (sc[0] * sc[5] * sc[11])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[5] - g_RT[7])) * (refC) * (sc[7] * sc[11])
    qdot = qf - qr
    qdot_out[34] = qdot
    wdot[0] -= qdot
    wdot[5] -= qdot
    wdot[7] += qdot

    # Reaction 0: H2 + M <=> 2 H + M
    k_f = 45800000000000 * np.exp((-1.4) * logT - (52535.8199100127) * invT)
    Corr = mixture + (1.5) * sc[1] + (11) * sc[4] + (-0.17) * sc[11]
    qf = Corr * k_f * (sc[1])
    qr = Corr * k_f * np.exp(-(-2.000000 * g_RT[0] + g_RT[1])) * (refCinv) * ((sc[0] ** 2))
    qdot = qf - qr
    qdot_out[0] = qdot
    wdot[0] += 2.000000 * qdot
    wdot[1] -= qdot

    # Reaction 1: 2 O + M <=> O2 + M
    k_f = 6160 * np.exp((-0.5) * logT)
    Corr = mixture + (1.5) * sc[1] + (11) * sc[4] + (-0.17) * sc[10] + (-0.17) * sc[11]
    qf = Corr * k_f * ((sc[2] ** 2))
    qr = Corr * k_f * np.exp(-(2.000000 * g_RT[2] - g_RT[5])) * (refC) * (sc[5])
    qdot = qf - qr
    qdot_out[1] = qdot
    wdot[2] -= 2.000000 * qdot
    wdot[5] += qdot

    # Reaction 5: H + OH + M <=> H2O + M
    k_f = 35000000000 * np.exp((-2) * logT)
    Corr = mixture + (-0.27) * sc[1] + (2.65) * sc[4] + (-0.62) * sc[10] + (-1) * sc[11]
    qf = Corr * k_f * (sc[0] * sc[3])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC) * (sc[4])
    qdot = qf - qr
    qdot_out[5] = qdot
    wdot[0] -= qdot
    wdot[3] -= qdot
    wdot[4] += qdot

    # Reaction 6: H + OH + HE <=> H2O + HE
    k_f = 404000000000000 * np.exp((-3.234) * logT)
    qf = k_f * (sc[0] * sc[3] * sc[11])
    qr = k_f * np.exp(-(g_RT[0] + g_RT[3] - g_RT[4])) * (refC) * (sc[4] * sc[11])
    qdot = qf - qr
    qdot_out[6] = qdot
    wdot[0] -= qdot
    wdot[3] -= qdot
    wdot[4] += qdot

    # Reaction 8: H + O + M <=> OH + M
    k_f = 4710000 * np.exp((-1) * logT)
    Corr = mixture + (1.5) * sc[1] + (11) * sc[4] + (-0.25) * sc[10] + (-0.25) * sc[11]
    qf = Corr * k_f * (sc[0] * sc[2])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[2] - g_RT[3])) * (refC) * (sc[3])
    qdot = qf - qr
    qdot_out[8] = qdot
    wdot[0] -= qdot
    wdot[2] -= qdot
    wdot[3] += qdot

    # Reaction 9: H + O + M <=> OHV + M
    k_f = 15 * np.exp(-(3006.7195781832) * invT)
    Corr = mixture + (5.5) * sc[4] + (-0.6) * sc[5] + (-0.6) * sc[6] + (-0.65) * sc[10]
    qf = Corr * k_f * (sc[0] * sc[2])
    qr = Corr * k_f * np.exp(-(g_RT[0] + g_RT[2] - g_RT[9])) * (refC) * (sc[9])
    qdot = qf - qr
    qdot_out[9] = qdot
    wdot[0] -= qdot
    wdot[2] -= qdot
    wdot[9] += qdot

    # Reaction 2: H2 + O <=> H + OH
    k_f = 0.0508 * np.exp((2.7) * logT - (3150.13632793755) * invT)
    qf = k_f * (sc[1] * sc[2])
    qr = k_f * np.exp(-(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3])) * (sc[0] * sc[3])
    qdot = qf - qr
    qdot_out[2] = qdot
    wdot[0] += qdot
    wdot[1] -= qdot
    wdot[2] -= qdot
    wdot[3] += qdot

    # Reaction 3: H2 + OH <=> H + H2O
    k_f = 43800000 * np.exp(-(3517.48449397499) * invT)
    qf = k_f * (sc[1] * sc[3])
    qr = k_f * np.exp(-(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4])) * (sc[0] * sc[4])
    qdot = qf - qr
    qdot_out[3] = qdot
    wdot[0] += qdot
    wdot[1] -= qdot
    wdot[3] -= qdot
    wdot[4] += qdot

    # Reaction 4: H + O2 <=> O + OH
    k_f = 104000000 * np.exp(-(7803.88405329959) * invT)
    qf = k_f * (sc[0] * sc[5])
    qr = k_f * np.exp(-(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5])) * (sc[2] * sc[3])
    qdot = qf - qr
    qdot_out[4] = qdot
    wdot[0] -= qdot
    wdot[2] += qdot
    wdot[3] += qdot
    wdot[5] -= qdot

    # Reaction 7: H2O + O <=> 2 OH
    k_f = 2.97 * np.exp((2.02) * logT - (6743.10332178324) * invT)
    qf = k_f * (sc[2] * sc[4])
    qr = k_f * np.exp(-(g_RT[2] - 2.000000 * g_RT[3] + g_RT[4])) * ((sc[3] ** 2))
    qdot = qf - qr
    qdot_out[7] = qdot
    wdot[2] -= qdot
    wdot[3] += 2.000000 * qdot
    wdot[4] -= qdot

    # Reaction 10: H2O + OHV <=> H2O + OH
    k_f = 5930000 * np.exp((0.5) * logT - (-432.766332592059) * invT)
    qf = k_f * (sc[4] * sc[9])
    qr = k_f * np.exp(-(-g_RT[3] + g_RT[4] - g_RT[4] + g_RT[9])) * (sc[3] * sc[4])
    qdot = qf - qr
    qdot_out[10] = qdot
    wdot[3] += qdot
    wdot[4] -= qdot
    wdot[4] += qdot
    wdot[9] -= qdot

    # Reaction 11: H2 + OHV <=> H2 + OH
    k_f = 2950000 * np.exp((0.5) * logT - (-223.428199617296) * invT)
    qf = k_f * (sc[1] * sc[9])
    qr = k_f * np.exp(-(g_RT[1] - g_RT[1] - g_RT[3] + g_RT[9])) * (sc[1] * sc[3])
    qdot = qf - qr
    qdot_out[11] = qdot
    wdot[1] -= qdot
    wdot[1] += qdot
    wdot[3] += qdot
    wdot[9] -= qdot

    # Reaction 12: N2 + OHV <=> N2 + OH
    k_f = 108000 * np.exp((0.5) * logT - (-624.995098929462) * invT)
    qf = k_f * (sc[6] * sc[9])
    qr = k_f * np.exp(-(-g_RT[3] + g_RT[6] - g_RT[6] + g_RT[9])) * (sc[3] * sc[6])
    qdot = qf - qr
    qdot_out[12] = qdot
    wdot[3] += qdot
    wdot[6] -= qdot
    wdot[6] += qdot
    wdot[9] -= qdot

    # Reaction 13: OH + OHV <=> 2 OH
    k_f = 6010000 * np.exp((0.5) * logT - (-384.457532674806) * invT)
    qf = k_f * (sc[3] * sc[9])
    qr = k_f * np.exp(-(g_RT[3] - 2.000000 * g_RT[3] + g_RT[9])) * ((sc[3] ** 2))
    qdot = qf - qr
    qdot_out[13] = qdot
    wdot[3] -= qdot
    wdot[3] += 2.000000 * qdot
    wdot[9] -= qdot

    # Reaction 14: H + OHV <=> H + OH
    k_f = 1310000 * np.exp((0.5) * logT - (-84.0371831893882) * invT)
    qf = k_f * (sc[0] * sc[9])
    qr = k_f * np.exp(-(g_RT[0] - g_RT[0] - g_RT[3] + g_RT[9])) * (sc[0] * sc[3])
    qdot = qf - qr
    qdot_out[14] = qdot
    wdot[0] -= qdot
    wdot[0] += qdot

    # Reaction 15: AR + OHV <=> AR + OH
    k_f = 1690000 * np.exp(-(2080.80091310252) * invT)
    qf = k_f * (sc[9] * sc[10])
    qr = k_f * np.exp(-(-g_RT[3] + g_RT[9] + g_RT[10] - g_RT[10])) * (sc[3] * sc[10])
    qdot = qf - qr
    qdot_out[15] = qdot
    wdot[3] += qdot
    wdot[9] -= qdot
    wdot[10] -= qdot
    wdot[10] += qdot

    # Reaction 16: OHV <=> OH
    k_f = 1450000
    qf = k_f * (sc[9])
    qr = k_f * np.exp(-(-g_RT[3] + g_RT[9])) * (sc[3])
    qdot = qf - qr
    qdot_out[16] = qdot
    wdot[3] += qdot
    wdot[9] -= qdot

    # Reaction 17: O2 + OHV <=> O2 + OH
    k_f = 2100000 * np.exp((0.5) * logT - (-240.537566254656) * invT)
    qf = k_f * (sc[5] * sc[9])
    qr = k_f * np.exp(-(-g_RT[3] + g_RT[5] - g_RT[5] + g_RT[9])) * (sc[3] * sc[5])
    qdot = qf - qr
    qdot_out[17] = qdot
    wdot[3] += qdot
    wdot[5] -= qdot
    wdot[5] += qdot
    wdot[9] -= qdot

    # Reaction 19: H + H2O2 <=> H2O + OH
    k_f = 24100000 * np.exp(-(1997.77016324474) * invT)
    qf = k_f * (sc[0] * sc[8])
    qr = k_f * np.exp(-(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[8])) * (sc[3] * sc[4])
    qdot = qf - qr
    qdot_out[19] = qdot
    wdot[0] -= qdot
    wdot[3] += qdot
    wdot[4] += qdot
    wdot[8] -= qdot

    # Reaction 20: H + H2O2 <=> H2 + HO2
    k_f = 21500 * np.exp((1) * logT - (3019.29999482832) * invT)
    qf = k_f * (sc[0] * sc[8])
    qr = k_f * np.exp(-(g_RT[0] - g_RT[1] - g_RT[7] + g_RT[8])) * (sc[1] * sc[7])
    qdot = qf - qr
    qdot_out[20] = qdot
    wdot[0] -= qdot
    wdot[1] += qdot
    wdot[7] += qdot
    wdot[8] -= qdot

    # Reaction 21: H2O2 + O <=> HO2 + OH
    k_f = 9.55 * np.exp((2) * logT - (1997.77016324474) * invT)
    qf = k_f * (sc[2] * sc[8])
    qr = k_f * np.exp(-(g_RT[2] - g_RT[3] - g_RT[7] + g_RT[8])) * (sc[3] * sc[7])
    qdot = qf - qr
    qdot_out[21] = qdot
    wdot[2] -= qdot
    wdot[3] += qdot
    wdot[7] += qdot
    wdot[8] -= qdot

    # Reaction 22: H2O2 + OH <=> H2O + HO2
    k_f = 1740000 * np.exp(-(160.022899725901) * invT)
    qf = k_f * (sc[3] * sc[8])
    qr = k_f * np.exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4] * sc[7])
    qdot = qf - qr
    qdot_out[22] = qdot
    wdot[3] -= qdot
    wdot[4] += qdot
    wdot[7] += qdot
    wdot[8] -= qdot

    # Reaction 23: H2O2 + OH <=> H2O + HO2
    k_f = 75900000 * np.exp(-(3657.88194373451) * invT)
    qf = k_f * (sc[3] * sc[8])
    qr = k_f * np.exp(-(g_RT[3] - g_RT[4] - g_RT[7] + g_RT[8])) * (sc[4] * sc[7])
    qdot = qf - qr
    qdot_out[23] = qdot
    wdot[3] -= qdot
    wdot[4] += qdot
    wdot[7] += qdot
    wdot[8] -= qdot

    # Reaction 24: H + HO2 <=> 2 OH
    k_f = 70800000 * np.exp(-(148.448916412392) * invT)
    qf = k_f * (sc[0] * sc[7])
    qr = k_f * np.exp(-(g_RT[0] - 2.000000 * g_RT[3] + g_RT[7])) * ((sc[3] ** 2))
    qdot = qf - qr
    qdot_out[24] = qdot
    wdot[0] -= qdot
    wdot[3] += 2.000000 * qdot
    wdot[7] -= qdot

    # Reaction 25: H + HO2 <=> H2 + O2
    k_f = 11400 * np.exp((1.083) * logT - (278.782032855815) * invT)
    qf = k_f * (sc[0] * sc[7])
    qr = k_f * np.exp(-(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[7])) * (sc[1] * sc[5])
    qdot = qf - qr
    qdot_out[25] = qdot
    wdot[0] -= qdot
    wdot[1] += qdot
    wdot[5] += qdot
    wdot[7] -= qdot

    # Reaction 26: HO2 + O <=> O2 + OH
    k_f = 32500000
    qf = k_f * (sc[2] * sc[7])
    qr = k_f * np.exp(-(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[7])) * (sc[3] * sc[5])
    qdot = qf - qr
    qdot_out[26] = qdot
    wdot[2] -= qdot
    wdot[3] += qdot
    wdot[5] += qdot
    wdot[7] -= qdot

    # Reaction 27: HO2 + OH <=> H2O + O2
    k_f = 31648300 * np.exp((1.380154522) * logT - (15400.9555526519) * invT)
    qf = k_f * (sc[3] * sc[7])
    qr = k_f * np.exp(-(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[7])) * (sc[4] * sc[5])
    qdot = qf - qr
    qdot_out[27] = qdot
    wdot[3] -= qdot
    wdot[4] += qdot
    wdot[5] += qdot
    wdot[7] -= qdot

    # Reaction 28: 2 HO2 <=> H2O2 + O2
    k_f = 12100 * np.exp((0.422) * logT - (-745.26388205679) * invT)
    qf = k_f * ((sc[7] * sc[7]))
    qr = k_f * np.exp(-(-g_RT[5] + 2.000000 * g_RT[7] - g_RT[8])) * (sc[5] * sc[8])
    qdot = qf - qr
    qdot_out[28] = qdot
    wdot[5] += qdot
    wdot[7] -= 2.000000 * qdot
    wdot[8] += qdot

    # Reaction 29: 2 HO2 <=> H2O2 + O2
    k_f = 16900000000 * np.exp((-0.681) * logT - (6507.59792218663) * invT)
    qf = k_f * ((sc[7] * sc[7]))
    qr = k_f * np.exp(-(-g_RT[5] + 2.000000 * g_RT[7] - g_RT[8])) * (sc[5] * sc[8])
    qdot = qf - qr
    qdot_out[29] = qdot
    wdot[5] += qdot
    wdot[7] -= 2.000000 * qdot
    wdot[8] += qdot

def CKYTCR(rho, T, y, imw):
    c = np.zeros(12)  # Initialize the concentration array with zeros

    for i in range(12):
        c[i] = rho * y[i] * imw[i]
    #Returns [mole/cm3]
    return c

def CKWC(T, C, wdot, qdot_out):
    # Convert to SI (from mol/m³ to mol/L)
    for id in range(12):
        C[id] *= 1.0e6
    
    g_RT = np.zeros(len(C))
    gibbs(g_RT, T)

    production_rate(wdot, C, g_RT, T, qdot_out)  

    # Convert back to Chemkin units
    for id in range(12):
        C[id] *= 1.0e-6
        wdot[id] *= 1.0e-6

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

print("Index H2O: ", gas.species_index('H2O'))

wdot_cantera = gas.net_production_rates*gas.molecular_weights
net_reaction_rates = gas.net_rates_of_progress
forward_rates = gas.forward_rates_of_progress  # kmol/m^3/s
reverse_rates = gas.reverse_rates_of_progress


density = gas.density #[kg/m3]
density *= 1.e-3 # [kg/m3] -> [g/cm3] CGS

massfrac = np.zeros(gas.n_species)
C        = np.zeros(gas.n_species)
wdot     = np.zeros(gas.n_species)
qdot     = np.zeros(gas.n_reactions)

for i,name in enumerate(d):
    massfrac[i] = d[name]

mw  = gas.molecular_weights
imw = 1.0/mw

C = CKYTCR(density, T, massfrac, imw)
CKWC(T, C, wdot, qdot)

# [mole/cm3/s] -> [g/cm3/s]
for n in range(gas.n_species):
      wdot[n] *= mw[n]
wdot *= 1000.0 #CGS -> MKS
qdot /= 1000.0 #Not sure which units here but it matches Cantera

print(f"{'Species':>15} {'LMeX Net Production Rate [kg/m³/s]':>30} {'Cantera Net Production Rate [kg/m³/s]':>30}")
print("="*100)
for i,name in enumerate(gas.species_names):
    wdot_diff = wdot[i] - wdot_cantera[i]
    # if(abs(wdot[i]) > 1.e-15):
    print(f"{name:>15} {wdot[i]:30.4e}  {wdot_cantera[i]:30.4e}")
    # print(f"{name:>15} {wdot_diff:30.4e} ")

print("")

# Display net reaction rates for each reaction
print(f"{'Reaction':>40} {'LMeX Net Rate':>20} {'Cantera Net Rate':>20}")
print("="*110)
rxn_str = gas.reaction_equations()
for i in range(gas.n_reactions):
    print(f"{i} {rxn_str[i]:>40} {qdot[i]:20.4e} {net_reaction_rates[i]:20.4e}")
print("")

ratio_net_reactions = qdot/net_reaction_rates
deviation = 1.02 #percentage

print(f"{'Reaction':>40} {'Ratio LMeX/Cantera for Net Rate':>20}")
print("Reactions presented here have {:1.2f}% deviation".format((deviation-1.0)*100.))
print("="*110)
for i in range(gas.n_reactions):
    if(ratio_net_reactions[i] > deviation or ratio_net_reactions[i] < (1.0-deviation)):
        print(f" {i} {rxn_str[i]:>40} {ratio_net_reactions[i]:20.4e}")

reaction_numbers = [30,32,33]
for reaction_number in reaction_numbers:
    print(rxn_str[reaction_number])
    print(forward_rates[reaction_number])
    print(reverse_rates[reaction_number])


