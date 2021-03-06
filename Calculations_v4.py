import math
import pandas as pd
import numpy as np
import time

### Units Converter
class UnitsConverter:
    class Pressure:
        def bara_to_kgfpcm2g(pressure_bara: float):
            return pressure_bara * 1.0197162129779282 - 1.033227
        def kgfpcm2g_to_bara(pressure_kgfpcm2: float):
            return (pressure_kgfpcm2 + 1.033227) / 1.0197162129779282
        def psi_to_kgfpcm2(pressure_psi: float):
            return pressure_psi * 0.0703069579640175
        def bar_to_psi(pressure_bar: float):
            return pressure_bar * 14.503773773
        def bar_to_kPa(pressure_bar: float):
            return pressure_bar * 100
        def kPa_to_psi(pressure_kpa: float):
            return pressure_kpa * 0.1450377377
        def psi_to_kPa(pressure_psi: float):
            return pressure_psi / 0.1450377377
    class Flowrate:
        def sm3d_to_sm3y(flowrate_sm3pday: float):
            return flowrate_sm3pday * 365
    class Temperature:
        def C_to_R(temperature_C: float):
            return (temperature_C + 273.15) * 9 / 5
        def R_to_K(temperature_R: float):
            return temperature_R * 5 / 9
        def C_to_K(temperature_C: float):
            return temperature_C + 273.15


def get_zfactor(Aj: float,
                Bj: float,
                phase: str):
    ### Solving cubic Peng-Robinson Equation of State
    print('\nSolving Peng-Robinson EOS for {} phase compressibility factor...'.format(phase))
    roots = []
    a2 = -(1 - Bj)
    a1 = Aj - 2 * Bj - 3 * Bj**2
    a0 = -(Aj * Bj - Bj**2 - Bj**3)
    q = 1/3 * a1 - 1/9 * a2**2
    r = 1/6 * (a1 * a2 - 3 * a0) - 1/27 * a2**3
    if (q**3 + r**2) > 0:
        print('\tEquation has one real root and a pair of complex conjugate roots...')
    elif (q**3 + r**2) == 0:
        print('\tEquation has all real roots and at leas two of them are equal...')
    else:
        print('\tEquation has all real roots...')
    if abs(q**3 + r**2) <= 1e-5:  # this convergence criteria may affect equicomp results
        smallvar = 0  # introduced to avoid problems when expression is too small
    else:
        smallvar = (q ** 3 + r ** 2) ** (1 / 2)
    print('q={},\tr={}'.format(q, r))
    print('a0={},\ta1={},\ta2={}'.format(a0, a1, a2))
    print('smallvar={}'.format(smallvar))
    print('q**3+r**2='.format(abs(q**3 + r**2)))
    s1 = np.cbrt(r + smallvar)
    s2 = np.cbrt(r - smallvar)
    # print('-->', s1, s2)
    z1 = (s1 + s2) - a2 / 3
    z2 = complex(-1/2 * (s1 + s2) - a2 / 3, (3**(1/2)) / 2 * (s1 - s2))
    z3 = complex(-1/2 * (s1 + s2) - a2 / 3, -(3**(1/2)) / 2 * (s1 - s2))
    check1 = abs((z1 + z2 +z3) - complex(-1 * a2, 0)) < 0.001
    check2 = abs((z1 * z2 + z1 * z3 + z2 * z3) - complex(a1, 0)) < 0.001
    check3 = abs((z1 * z2 * z3) - complex(-1 * a0, 0)) < 0.001
    # print('-->', z1, z2, z3)
    if check1 and check2 and check3:
        print('\tRoots checked successfully!')
    else:
        print('\tCheck1: {} = {} --> {}'.format(z1 + z2 + z3, complex(-1 * a2, 0), check1))
        print('\tCheck2: {} = {} --> {}'.format(z1 * z2 + z1 * z3 + z2 * z3, complex(a1, 0), check2))
        print('\tCheck3: {} = {} --> {}'.format(z1 * z2 * z3, complex(-1 * a0, 0), check3))
        print('WARNING! Roots are NOT checked successfully!')
    for root in [z1, z2, z3]:
        if abs(root.imag) < 10**-6:
            root = root.real
        roots.append(root)
    ### Selecting proper root for compressibility factor
    zfactors = []
    for root in roots:
        if not(type(root) == complex) and root >= 0 and root >= Bj:
            zfactors.append(root)
    if len(zfactors) > 1:
        if phase == 'vapor':
            zfactor = min(zfactors)
        else:
            zfactor = max(zfactors)
    else:
        # zfactor = zfactors[0]
        try:
            zfactor = zfactors[0]
        except:  # ATTENTION!
            if Bj > 0:
                zfactor = Bj
            else:
                zfactor = 10**-5
    print()
    return zfactor


def Kvalues_comparison(Kvalues_df1: pd.DataFrame,  # For the time - comparison only for one interaction vapor-liquid
                          Kvalues_df2: pd.DataFrame):
    sum = np.sum((np.array(Kvalues_df1['Kign']) - np.array(Kvalues_df2['Kign'])) ** 2)
    return math.sqrt(sum / len(Kvalues_df1.index))


### Equlibrium composition calculation via K-values (as in GPSA M25) - works good for VLE
def get_equilibrium_composition_v1(streamcompostion: pd.DataFrame,
                                   Kvalues_df: pd.DataFrame,
                                   show_log: bool):
    locconvcrit = 1e-4
    start_time = time.perf_counter()
    print('\nCalculating equilibrium compositions...')
    df = pd.DataFrame({'xi': streamcompostion['Content [mol. fract.]'],
                       'Kign': Kvalues_df['Kign']}, index= streamcompostion.index)
    result_df = pd.DataFrame(columns= ['vapor', 'liquid'], index= df.index)
    gn_convergence_check = 0
    print('\tVapor-Liquid Equilibria...')
    def convergence_func(df: pd.DataFrame, L: float):
        Kign_arr = np.array(df['Kign'])
        res_arr = np.array(df['xi']) / (L + (1 - L) * Kign_arr)
        check_parameter = np.sum(res_arr) - 1
        result_df['liquid'] = res_arr
        result_df['vapor'] = res_arr * Kign_arr
        return check_parameter
    L_left = 0
    L_right = 1
    L_mid = (L_left + L_right) / 2
    for i in range(15):
        L_mid = (L_left + L_right) / 2
        check_l = convergence_func(df, L_left)
        check_mid = convergence_func(df, L_mid)
        if check_l * check_mid < 0:
            L_right = L_mid
        else:
            L_left = L_mid
        if show_log:
            print('\t\tvapor-liquid,\tstep-{:d},\tL = {:.4f},\tcheck = {:.3e}'.format(i + 1, L_mid, check_mid))
        if abs(check_mid) <= locconvcrit:
            gn_convergence_check = 1
            break
    if gn_convergence_check == 1:
        print('\tEquilibrium composition vapor/liquid converged!')
    else:
        if abs(check_mid) <= 0.05:
            print('\tEquilibrium composition vapor/liquid DID NOT converge in 15 interation, but still OK')
        else:
            print('\tWARNING! Equilibrium composition vapor/liquid DID NOT converge!')
    print('\tIteration time {:.3f} seconds'.format(time.perf_counter() - start_time))
    return result_df, L_mid


# Supposed to work with three phases but do not
def get_equilibrium_composition_v2(streamcompostion: pd.DataFrame,
                                   Kvalues_df: pd.DataFrame,
                                   show_log: bool):
    start_time = time.perf_counter()
    print('\nCalculating equilibrium compositions...')
    df = pd.DataFrame({'xi': streamcompostion['Content [mol. fract.]'],
                       'Kign': Kvalues_df['Kign'],
                       'Kigq': Kvalues_df['Kigq']}, index= streamcompostion.index)
    ANAE_result_df = pd.DataFrame(columns= ['non-aqueous', 'aqueous'], index= df.index)
    ANAE_convergence_check = 0
    print('\tVapor-non-Aqueous Equilibria...')
    def convergence_func_na(df: pd.DataFrame, W: float):
        for component in df.index:
            xiq = df.loc[component]['xi'] / (W + (1 - W) * df.loc[component]['Kigq'])
            xig = xiq * df.loc[component]['Kigq']
            ANAE_result_df.loc[component]['non-aqueous'] = xig
            ANAE_result_df.loc[component]['aqueous'] = xiq
        check_parameter = ANAE_result_df['aqueous'].sum() - 1
        return check_parameter
    W_left = 0
    W_right = 1
    W_mid = (W_left + W_right) / 2
    for i in range(50):
        W_mid = (W_left + W_right) / 2
        check_l = convergence_func_na(df, W_left)
        check_mid = convergence_func_na(df, W_mid)
        if check_l * check_mid < 0:
            W_right = W_mid
        else:
            W_left = W_mid
        if show_log:
            print('\t\tnon-aqueous/aqueous,\tstep-{:d},\tW = {:.4f},\tcheck = {:.3e}'.format(i + 1, W_mid, check_mid))
        if abs(check_mid) <= 10**-3:
            ANAE_convergence_check = 1
            break
    if ANAE_convergence_check == 1:
        print('\tEquilibrium composition vapor/liquid converged!')
    else:
        if abs(check_mid) <= 0.05:
            print('\tEquilibrium composition vapor/liquid DID NOT converge in 15 interation, but still OK')
        else:
            print('\tWARNING! Equilibrium composition vapor/liquid DID NOT converge!')
    print('\tIteration time {:.3} seconds'.format(time.perf_counter() - start_time))


    VLE_result_df = pd.DataFrame(columns=['vapor', 'liquid'], index=df.index)
    VLE_convergence_check = 0
    print('\tVapor-Liquid Equilibria...')
    def convergence_func_aq(df: pd.DataFrame, L: float):
        for component in df.index:
            xin = ANAE_result_df.loc[component]['non-aqueous'] / (L + (1 - L) * df.loc[component]['Kign'])
            xig = xin * df.loc[component]['Kign']
            VLE_result_df.loc[component]['vapor'] = xig
            VLE_result_df.loc[component]['liquid'] = xin
        check_parameter = VLE_result_df['liquid'].sum() - 1
        return check_parameter

    L_left = 0
    L_right = 1
    L_mid = (L_left + L_right) / 2
    for i in range(15):
        L_mid = (L_left + L_right) / 2
        check_l = convergence_func_aq(df, L_left)
        check_mid = convergence_func_aq(df, L_mid)
        if check_l * check_mid < 0:
            L_right = L_mid
        else:
            L_left = L_mid
        if show_log:
            print('\t\tvapor/liquid,\tstep-{:d},\tL = {:.4f},\tcheck = {:.3e}'.format(i + 1, L_mid, check_mid))
        if abs(check_mid) <= 10 ** -3:
            VLE_convergence_check = 1
            break
    if VLE_convergence_check == 1:
        print('\tEquilibrium composition vapor/liquid converged!')
    else:
        if abs(check_mid) <= 0.05:
            print('\tEquilibrium composition vapor/liquid DID NOT converge in 15 interation, but still OK')
        else:
            print('\tWARNING! Equilibrium composition vapor/liquid DID NOT converge!')
    print('\tIteration time {:.3} seconds'.format(time.perf_counter() - start_time))
    result_df = pd.DataFrame(columns= ['vapor', 'liquid', 'aqueous'], index= df.index)
    result_df['vapor'] = VLE_result_df['vapor']
    result_df['liquid'] = VLE_result_df['liquid']
    result_df['aqueous'] = ANAE_result_df['aqueous']
    return result_df, W_mid, L_mid


def get_initial_Kvalues(comppropDB: pd.DataFrame,
                        streamcomp: pd.DataFrame,
                        P: float,
                        T: float):  ### Pressure and Temperature in field units
    ### Calculation of initial guesses for K-values
    ### Basic properties
    Pc_arr = UnitsConverter.Pressure.kPa_to_psi(np.array(comppropDB['Pcrit [kPa]'])) # Field units
    Tc_arr = UnitsConverter.Temperature.C_to_R(comppropDB['Tcrit [C]']) # Field units
    w_arr = comppropDB['Acentricity']
    Pr_arr = P / Pc_arr
    Tr_arr = T / Tc_arr
    ### K-values initial guess
    Kvalues_df = pd.DataFrame(columns=['Kign', 'Kigq'], index=streamcomp.index)
    Kvalues_df['Kign'] = 1 / Pr_arr * math.e ** (5.37 * (1 + w_arr) * (1 - 1 / Tr_arr))
    Kvalues_df['Kigq'] = 10 ** 6 * (Pr_arr / Tr_arr)
    return Kvalues_df



def get_compdepvar(comppropDB: pd.DataFrame,
                   streamcomp: pd.DataFrame,
                   T: float): ### Temperature in field units
    ### Calculation of component-dependent variables for fugacity coefficients
    R_field = 10.731577089016  # [psi*ft3/(lbmol*R)] - Field
    Pc_arr = UnitsConverter.Pressure.kPa_to_psi(np.array(comppropDB['Pcrit [kPa]']))  # Field units
    Tc_arr = UnitsConverter.Temperature.C_to_R(comppropDB['Tcrit [C]'])  # Field units
    w_arr = comppropDB['Acentricity']
    Tr_arr = T / Tc_arr
    b_i_arr = 0.07780 * R_field * Tc_arr / Pc_arr
    kappa_arr = np.where(w_arr > 0.49,
                         0.379642 + 1.4853 * w_arr - 0.164423 * w_arr ** 2 + 0.01666 * w_arr ** 3, ### 1980 modification for heavy hydrocarbon components
                         0.37464 + 1.5422 * w_arr - 0.26992 * (w_arr ** 2))
    # print('new kappa\n', kappa_arr)
    alfa_arr = np.where(np.logical_and(Tc_arr == 374.149011230469, Tr_arr ** 0.5 < 0.85),
                        (1.0085677 + 0.82154 * (1 - Tr_arr ** 0.5)) ** 2, ### 1980 modification for water component
                        (1 + kappa_arr * (1 - Tr_arr ** 0.5)) ** 2)
    ac_arr = 0.45724 * (R_field ** 2) * (Tc_arr ** 2) / Pc_arr
    a_i_arr = ac_arr * alfa_arr
    compvar_df = pd.DataFrame(columns=['ai', 'bi'], index=streamcomp.index)
    compvar_df['ai'] = a_i_arr
    compvar_df['bi'] = b_i_arr
    return compvar_df


def get_phasedepvar(equicomp_df: pd.DataFrame,
                    compvar_df: pd.DataFrame,
                    binarycoefDB: pd.DataFrame,
                    P: float,
                    T: float):
    R_field = 10.731577089016  # [psi*ft3/(lbmol*R)] - Field
    phasevar_df = pd.DataFrame(columns=['aj', 'bj', 'Aj', 'Bj'], index=['vapor', 'liquid'])
    for phase in phasevar_df.index:
        # Mixing rules for each phase
        bj = (equicomp_df[phase] * compvar_df['bi']).sum()
        aj = 0
        for component_1 in equicomp_df.index:
            aj += (equicomp_df.loc[component_1][phase] * equicomp_df[phase] * (compvar_df.loc[component_1]['ai'] * compvar_df['ai']) ** 0.5\
                                                * (1 - binarycoefDB[component_1])).sum()
        phasevar_df.loc[phase]['aj'] = aj
        phasevar_df.loc[phase]['bj'] = bj
        ### Phase-dependent variables calculation
        phasevar_df.loc[phase]['Aj'] = phasevar_df.loc[phase]['aj'] * P / R_field ** 2 / T ** 2
        phasevar_df.loc[phase]['Bj'] = phasevar_df.loc[phase]['bj'] * P / R_field / T
    return phasevar_df


def get_phasecompdepvar(phasevar_df: pd.DataFrame,
                        compvar_df: pd.DataFrame,
                        equicomp_df: pd.DataFrame,
                        binarycoefDB: pd.DataFrame):
    ### Calculatin phase-component dependent variable for
    Aijprime_df = pd.DataFrame(columns=phasevar_df.index, index=compvar_df.index)
    Bijprime_df = pd.DataFrame(columns=phasevar_df.index, index=compvar_df.index)
    for phase in phasevar_df.index:
        Bijprime_df[phase] = compvar_df['bi'] / phasevar_df.loc[phase]['bj']
        for component in compvar_df.index:
            aux_var = (equicomp_df[phase] * (compvar_df['ai']) ** 0.5 * (1 - binarycoefDB[component])).sum()
            Aijprime_df.loc[component][phase] = 1 / phasevar_df.loc[phase]['aj'] * (
                        2 * math.sqrt(compvar_df.loc[component]['ai']) * aux_var)
    return Aijprime_df, Bijprime_df


def get_fugacities(streamcomp:pd.DataFrame,
                   phasevar_df: pd.DataFrame,
                   Aijprime_df: pd.DataFrame,
                   Bijprime_df: pd.DataFrame,
                   zfactors: pd.DataFrame):
    ### Fugacity coefficients calculation
    fugacit_df = pd.DataFrame(columns=phasevar_df.index, index=streamcomp.index)
    for phase in phasevar_df.index:
        Aj = phasevar_df.loc[phase]['Aj']
        Bj = phasevar_df.loc[phase]['Bj']
        for component in streamcomp.index:
            Aijprime = Aijprime_df.loc[component][phase]
            Bijprime = Bijprime_df.loc[component][phase]
            zj = zfactors.loc[phase]['zj']
            try:
                f = math.exp(-math.log(zj - Bj) + (zj - 1) * Bijprime \
                             - Aj / (2 * math.sqrt(2) * Bj) * (Aijprime - Bijprime) \
                             * math.log((zj + (math.sqrt(2) + 1) * Bj) / (zj - (math.sqrt(2) - 1) * Bj)))
            except:
                f = None
            fugacit_df.loc[component][phase] = f
    return fugacit_df


def get_Kvalues(fugacit_df: pd.DataFrame):
    ### Updated K-values calculation
    Kvalues_df = pd.DataFrame(columns=['Kign', 'Kigq'], index= fugacit_df.index)
    Kvalues_df['Kign'] = np.where(np.logical_or(fugacit_df['liquid'] is None, fugacit_df['vapor'] is None),
                                  1,
                                  fugacit_df['liquid'] / fugacit_df['vapor'])
    Kvalues_df['Kigq'] = np.array([None] * len(Kvalues_df['Kign']))
    return Kvalues_df


def redefine_equicomp(equicomp_df_in: pd.DataFrame,
                      equicomp_df_out: pd.DataFrame,
                      phase_fractions: pd.Series,
                      L: float,
                      phaseinteractnum: int,
                      phases_num: int):
    ### Arranging proper compositions and phase fractions to each phase
    if phases_num == 3:
        if phaseinteractnum == 0:
            equicomp_df_out['aqueous'] = equicomp_df_in['liquid']
            equicomp_df_out['vapor'] = equicomp_df_in['vapor']
            phase_fractions['Q'] = L
        elif phaseinteractnum == 1:
            equicomp_df_out['vapor'] = equicomp_df_in['vapor']
            equicomp_df_out['liquid'] = equicomp_df_in['liquid']
            phase_fractions['V'] = (1 - phase_fractions['Q']) * (1 - L)
            phase_fractions['L'] = (1 - phase_fractions['Q']) * L
    else:
        equicomp_df_out['vapor'] = equicomp_df_in['vapor']
        equicomp_df_out['liquid'] = equicomp_df_in['liquid']
        phase_fractions['V'] = 1 - L
        phase_fractions['L'] = L
    return equicomp_df_out, phase_fractions


### All functions combined
def flash_calc_PR_EOS(comppropDB: pd.DataFrame,
                      binarycoefDB: pd.DataFrame,
                      input_streamcomp: pd.DataFrame,
                      P_field: float,
                      T_field: float,
                      convcrit,
                      steps_limit):
    equicomp_df_sum = pd.DataFrame(columns=['vapor', 'liquid', 'aqueous'], index=input_streamcomp.index)
    streamcomp = pd.DataFrame(columns=['Content [mol. fract.'], index=input_streamcomp.index)
    phase_fractions = pd.Series([np.NaN, np.NaN, np.NaN], index=['V', 'L', 'Q'])
    if input_streamcomp.loc['H2O']['Content [mol. fract.]'] == 0:
        phases_num = 2
    else:
        phases_num = 3
    for interphase in range(phases_num - 1):
        if interphase == 0:
            streamcomp = input_streamcomp
        else:
            streamcomp['Content [mol. fract.]'] = equicomp_df_sum['vapor']

        ### STEP - 1: K's estimation using eq. (3-5) and (3-8)
        Kvalues_init_df = get_initial_Kvalues(comppropDB,
                                                   streamcomp,
                                                   P_field,
                                                   T_field)
        err_list = list()
        calc_err = 10 ** 6
        steps = 0
        while calc_err >= convcrit:
            steps += 1
            ### STEP - 2: Equlibrium compositions of phases calculation using K's from Step 1 (Methodology from GPSA)
            ### Basic properties
            compvar_df = get_compdepvar(comppropDB,
                                             streamcomp,
                                             T_field)
            ### Equlibrium copositions
            equicomp_df, L = get_equilibrium_composition_v1(streamcomp,
                                                                 Kvalues_init_df,
                                                                 True)
            ### STEP 3: - Fugacity coefficients calculation
            phasevar_df = get_phasedepvar(equicomp_df,
                                               compvar_df,
                                               binarycoefDB,
                                               P_field,
                                               T_field)
            ### Phase-component dependent variables calculation
            Aijprime_df, Bijprime_df = get_phasecompdepvar(phasevar_df,
                                                                compvar_df,
                                                                equicomp_df,
                                                                binarycoefDB)
            ### STEP 4: - Compressibility factors calculation
            zfactors = pd.DataFrame(columns=['zj'], index=phasevar_df.index)
            for phase in phasevar_df.index:
                zfactors.loc[phase]['zj'] = get_zfactor(phasevar_df.loc[phase]['Aj'],
                                                             phasevar_df.loc[phase]['Bj'],
                                                             phase)
            ### Fugacity coefficients factors calculation
            fugacit_df = get_fugacities(streamcomp,
                                             phasevar_df,
                                             Aijprime_df,
                                             Bijprime_df,
                                             zfactors)
            ### STEP 5: - New set of K-values calculation
            Kvalues_df = get_Kvalues(fugacit_df)
            calc_err = Kvalues_comparison(Kvalues_init_df, Kvalues_df)
            err_list.append(calc_err)
            Kvalues_init_df['Kign'] = Kvalues_df['Kign']
            print('K-values error at iteration: {:.3e}'.format(calc_err))
            if steps > steps_limit:
                print('WARNING: K-values did not converged!')
                break

        equicomp_df_final, L_final = get_equilibrium_composition_v1(streamcomp,
                                                                         Kvalues_init_df,
                                                                         True)
        if (abs(Kvalues_df['Kign'] - 1) < 10 ** -3).all():
            equicomp_df_sum['vapor'] = streamcomp['Content [mol. fract.]']
            phase_fractions['V'], phase_fractions['L'], phase_fractions['Q'] = 1, 0, 0
            break
        equicomp_df_sum, phase_fractions = redefine_equicomp(equicomp_df_final,
                                                                  equicomp_df_sum,
                                                                  phase_fractions,
                                                                  L_final,
                                                                  interphase,
                                                                  phases_num)
        if steps <= steps_limit:
            print('Converged in {} iterations\n'.format(steps - 1))
    return equicomp_df_sum, phase_fractions, zfactors


def get_phase_molar_weigh(comppropDB: pd.DataFrame,
                          equicomp_df: pd.DataFrame,
                          phase: str):
    MW = (np.array(equicomp_df[phase]) * np.array(comppropDB['MW [g/mole]'])).sum()
    return MW  # [g/mole]


def get_liquid_phase_density(comppropDB: pd.DataFrame,
                             equicomp_df: pd.DataFrame,
                             T_field: float,
                             phase: str):
    ### Mixing rules
    wSRKmix = (equicomp_df[phase] * comppropDB['SRK Acentricity']).sum()
    auxvar1 = (equicomp_df[phase] * comppropDB['Characteristic Volume [m3/kgmole]']).sum()
    auxvar2 = (equicomp_df[phase] * (comppropDB['Characteristic Volume [m3/kgmole]'] ** (1/3))).sum()
    auxvar3 = (equicomp_df[phase] * (comppropDB['Characteristic Volume [m3/kgmole]'] ** (2/3))).sum()
    Vasteriskmix = 0.25 * (auxvar1 + 3 * auxvar2 * auxvar3)
    auxvar4 = 0
    for component1 in equicomp_df.index:
        Tc1 = UnitsConverter.Temperature.C_to_R(comppropDB.loc[component1]['Tcrit [C]'])
        Tc2 = UnitsConverter.Temperature.C_to_R(comppropDB['Tcrit [C]'])
        auxvar4 += (equicomp_df.loc[component1][phase] * equicomp_df['liquid'] * \
                   np.sqrt(Tc1 * comppropDB.loc[component1]['Characteristic Volume [m3/kgmole]'] * \
                           Tc2 * comppropDB['Characteristic Volume [m3/kgmole]'])).sum()
    Tcmix = auxvar4 / Vasteriskmix
    ### Main variables calculation
    Tr = T_field / Tcmix
    if 0.25 < Tr < 1.0:
        VdeltaR = (-0.296123 + 0.386914 * Tr - 0.0427258 * Tr ** 2 - 0.0480645 * Tr ** 3) / (Tr - 1.00001)
        if Tr < 0.95:
            V0R = 1 - 1.52816 * (1 - Tr) ** (1 / 3) + 1.43907 * (1 - Tr) ** (2 / 3) - 0.81446 * (1 - Tr) + 0.190454 * (1 - Tr) ** (4 / 3)
        else:
            print('\tCostald density correlation ERROR! (Tr >= 0.95 for {} phase)'.format(phase))
            return np.NaN
    else:
        print('\tCostald density correlation ERROR! (Tr >= 1.0 for {} phase)'.format(phase))
        return np.NaN
    ### Density calculation
    Vs = V0R * (1 - wSRKmix * VdeltaR) * Vasteriskmix  # [m3/kgmol]
    MW = get_phase_molar_weigh(comppropDB, equicomp_df, phase)  # [kg/kgmole]
    density = MW / Vs  # [kg/m3]
    return density


def get_vapor_phase_density(comppropDB: pd.DataFrame,
                            equicomp_df: pd.DataFrame,
                            P_field: float,
                            T_field: float,
                            zfactor: float):
    P = UnitsConverter.Pressure.psi_to_kPa(P_field)
    T = UnitsConverter.Temperature.R_to_K(T_field)
    R = 8.31446261815324 # [J/(mole*K)]
    MW = get_phase_molar_weigh(comppropDB, equicomp_df, 'vapor')
    density = MW / (zfactor * R * T / P)
    return density

def get_mix_density(denisties: pd.Series,
                                 phase_fractions: pd.Series,
                                 phaseMW: pd.Series):
    mass = phase_fractions * phaseMW
    volume = mass / denisties
    density_mix = mass.sum() / volume.sum()
    return  density_mix