'''
_v4.3
'''

import math
import sys

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
                phase: str,
                show_log: bool):
    ### Solving cubic Peng-Robinson Equation of State
    if show_log:
        print('\nSolving Peng-Robinson EOS for {} phase compressibility factor...'.format(phase))
    roots = []
    a2 = -(1 - Bj)
    a1 = Aj - 2 * Bj - 3 * Bj**2
    a0 = -(Aj * Bj - Bj**2 - Bj**3)
    q = 1/3 * a1 - 1/9 * a2**2
    r = 1/6 * (a1 * a2 - 3 * a0) - 1/27 * a2**3
    if show_log:
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
    if np.isnan(smallvar):
        smallvar = 0
    # print('smallvar param =', (q ** 3 + r ** 2))
    # print('smallvar =',smallvar)
    # print('r =', r)
    # print('r + smallvar =', r +smallvar)
    s1 = np.cbrt(r + smallvar)
    s2 = np.cbrt(r - smallvar)
    # print('s1 = {}\ns2 = {}'.format(s1, s2))
    z1 = (s1 + s2) - a2 / 3
    z2 = complex(-1/2 * (s1 + s2) - a2 / 3, (3**(1/2)) / 2 * (s1 - s2))
    z3 = complex(-1/2 * (s1 + s2) - a2 / 3, -(3**(1/2)) / 2 * (s1 - s2))
    check1 = abs((z1 + z2 +z3) - complex(-1 * a2, 0)) < 0.001
    check2 = abs((z1 * z2 + z1 * z3 + z2 * z3) - complex(a1, 0)) < 0.001
    check3 = abs((z1 * z2 * z3) - complex(-1 * a0, 0)) < 0.001
    # print('-->', z1, z2, z3)
    if show_log:
        if check1 and check2 and check3:
            print('\tRoots checked successfully!\n')
        else:
            print('\tCheck1: {} = {} --> {}'.format(z1 + z2 + z3, complex(-1 * a2, 0), check1))
            print('\tCheck2: {} = {} --> {}'.format(z1 * z2 + z1 * z3 + z2 * z3, complex(a1, 0), check2))
            print('\tCheck3: {} = {} --> {}'.format(z1 * z2 * z3, complex(-1 * a0, 0), check3))
            print('WARNING! Roots are NOT checked successfully!\n')
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
                print('Bj =', Bj)
            else:
                zfactor = 10**-5
    return zfactor


def Kvalues_comparison(Kvalues_df1: pd.DataFrame,  # For the time - comparison only for one interaction vapor-liquid
                       Kvalues_df2: pd.DataFrame):
    sum = np.sum((np.array(Kvalues_df1['Kij']) - np.array(Kvalues_df2['Kij'])) ** 2)
    return math.sqrt(sum / len(Kvalues_df1.index))


### Equlibrium composition calculation via K-values (as in GPSA M25) - works good for VLE
def get_equilibrium_composition_v1(streamcompostion: pd.DataFrame,
                                   Kvalues_df: pd.DataFrame,
                                   show_log: bool):
    locconvcrit = 1e-10  ### This criteria directly affects stability tests
    start_time = time.perf_counter()
    if show_log:
        print('\nCalculating equilibrium compositions...')
    df = pd.DataFrame({'xi': streamcompostion['Content [mol. fract.]'],
                       'Kij': Kvalues_df['Kij']}, index= streamcompostion.index)
    result_df = pd.DataFrame(columns= ['vapor', 'liquid'], index= df.index)
    convergence_check = 0
    def convergence_func(df: pd.DataFrame, L: float):
        Kij_arr = np.array(df['Kij'])
        res_arr = np.array(df['xi']) / (L + (1 - L) * Kij_arr)
        check_parameter = np.sum(res_arr) - 1
        result_df['liquid'] = res_arr
        result_df['vapor'] = res_arr * Kij_arr
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
            convergence_check = 1
            break
    if convergence_check == 1:
        if show_log:
            print('\tEquilibrium composition vapor/liquid converged!')
    else:
        if abs(check_mid) <= 0.05:
            if True:
                print('\tEquilibrium composition vapor/liquid DID NOT converge in 15 interation, but still OK')
        else:
            print('\tWARNING! Equilibrium composition vapor/liquid DID NOT converge!')
    if show_log:
        print('\tIteration time {:.3f} seconds'.format(time.perf_counter() - start_time))
    return result_df, L_mid


def get_equilibrium_composition_v2(streamcompostion: pd.DataFrame,
                                   Kvalues_df: pd.DataFrame,
                                   locconvcrit: float,   ### This criteria directly affects stability tests
                                   show_log: bool):
    start_time = time.perf_counter()
    if show_log:
        print('\nCalculating equilibrium compositions...')
    df = pd.DataFrame({'xi': streamcompostion['Content [mol. fract.]'],
                       'Kij': Kvalues_df['Kij']}, index= streamcompostion.index)
    result_df = pd.DataFrame(columns= ['vapor', 'liquid'], index= df.index)
    def convergence_func(df: pd.DataFrame, L: float):
        Kij_arr = np.array(df['Kij'])
        res_arr = np.array(df['xi']) / (L + (1 - L) * Kij_arr)
        check_parameter = np.sum(res_arr) - 1
        result_df['liquid'] = res_arr
        result_df['vapor'] = res_arr * Kij_arr
        return check_parameter
    L_left = 0
    L_right = 1
    L_mid = (L_left + L_right) / 2
    steps_counter = 0
    check_mid = locconvcrit + 1
    check_l = check_mid + 1
    while abs(check_mid) > locconvcrit:
        if steps_counter > 150 or abs(check_mid - check_l) <= locconvcrit:
            print('ERROR! Equilibrium compositions did not converged in {} iterations! Residual error {:.3e}'.format(steps_counter,
                                                                                                                     (check_mid)))
            break
        steps_counter += 1
        L_mid = (L_left + L_right) / 2
        check_l = convergence_func(df, L_left)
        check_mid = convergence_func(df, L_mid)
        if check_l * check_mid < 0:
            L_right = L_mid
        else:
            L_left = L_mid
        if show_log:
            print('\t\tvapor-liquid,\tstep-{:d},\tL = {:.4f},\tcheck = {:.3e}'.format(steps_counter, L_mid, check_mid))
    if show_log:
        print('Converged in {} steps. Residual error: {:.2e}'.format(steps_counter, check_mid))
    return result_df, L_mid


# Supposed to work with three phases but do not
def get_equilibrium_composition_v3(streamcompostion: pd.DataFrame,
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
                        T: float):# Pressure and Temperature in field units
    ### Calculation of initial guesses for K-values
    ### Basic properties
    Pc_arr = UnitsConverter.Pressure.kPa_to_psi(np.array(comppropDB['Pcrit [kPa]'])) # Field units
    Tc_arr = UnitsConverter.Temperature.C_to_R(comppropDB['Tcrit [C]']) # Field units
    w_arr = comppropDB['Acentricity']
    Pr_arr = P / Pc_arr
    Tr_arr = T / Tc_arr
    ### K-values initial guess
    Kvalues_df = pd.DataFrame(columns=['Kij'], index=streamcomp.index)
    Kvalues_df['Kij'] = 1 / Pr_arr * math.e ** (5.37 * (1 + w_arr) * (1 - 1 / Tr_arr))
    '''
    Expression below for vapor - aqueous phase K-values initial guesses does not lead to adequate results
        Kvalues_df['Kij'] = 10 ** 6 * (Pr_arr / Tr_arr)
    '''
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
                   zfactors: dict):
    ### Fugacity coefficients calculation
    fugacit_df = pd.DataFrame(columns=phasevar_df.index, index=streamcomp.index)
    for phase in phasevar_df.index:
        Aj = phasevar_df.loc[phase]['Aj']
        Bj = phasevar_df.loc[phase]['Bj']
        for component in streamcomp.index:
            Aijprime = Aijprime_df.loc[component][phase]
            Bijprime = Bijprime_df.loc[component][phase]
            zj = zfactors[phase]
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
    Kvalues_df = pd.DataFrame(columns=['Kij'], index= fugacit_df.index)
    Kvalues_df['Kij'] = np.where(np.logical_or(fugacit_df['liquid'] is None, fugacit_df['vapor'] is None),
                                  1,
                                  fugacit_df['liquid'] / fugacit_df['vapor'])
    return Kvalues_df


def redefine_equicomp(equicomp_df_in: pd.DataFrame,
                      equicomp_df_out: pd.DataFrame,
                      phase_fractions: dict,
                      L: float,
                      Kvalues_single_df: pd.DataFrame,
                      Kvalues_df_out: pd.DataFrame,
                      phaseinteractnum: int,
                      phases_num: int):
    ### Arranging proper compositions and phase fractions to each phase
    if phases_num == 3:
        if phaseinteractnum == 0:
            equicomp_df_out['aqueous'] = equicomp_df_in['liquid']
            equicomp_df_out['vapor'] = equicomp_df_in['vapor']
            phase_fractions['aqueous'] = L
            Kvalues_df_out['Kigq'] = Kvalues_single_df['Kij']
        elif phaseinteractnum == 1:
            if abs(equicomp_df_in['vapor'].sum() - 1) > 0.01:
                equicomp_df_out['vapor'] = equicomp_df_in['vapor']
            elif abs(equicomp_df_in['liquid'].sum() - 1) <= 0.01:
                equicomp_df_out['vapor'] = equicomp_df_in['vapor']
            equicomp_df_out['liquid'] = equicomp_df_in['liquid']
            phase_fractions['vapor'] = (1 - phase_fractions['aqueous']) * (1 - L)
            phase_fractions['liquid'] = (1 - phase_fractions['aqueous']) * L
            Kvalues_df_out['Kign'] = Kvalues_single_df['Kij']
    else:
        equicomp_df_out['vapor'] = equicomp_df_in['vapor']
        equicomp_df_out['liquid'] = equicomp_df_in['liquid']
        equicomp_df_out['aqueous'] = pd.Series([np.nan] * len(equicomp_df_in.index))
        phase_fractions['vapor'] = 1 - L
        phase_fractions['liquid'] = L
        Kvalues_df_out['Kign'] = Kvalues_single_df['Kij']
    return equicomp_df_out, phase_fractions, Kvalues_df_out


### All functions combined
def flash_calc_PR_EOS_obsolete(comppropDB: pd.DataFrame,
                      binarycoefDB: pd.DataFrame,
                      input_streamcomp: pd.DataFrame,
                      P_field: float,
                      T_field: float,
                      convcrit,
                      steps_limit):
    equicomp_df_threephases = pd.DataFrame(columns=['vapor', 'liquid', 'aqueous'], index=input_streamcomp.index)
    streamcomp = pd.DataFrame(columns=['Content [mol. fract.'], index=input_streamcomp.index)
    phase_fractions = dict({'vapor' : 0., 'liquid' : 0., 'aqueous' : 0.})
    Kvalues_single_df = pd.DataFrame(columns=['Kij'], index=input_streamcomp.index) # df used to store K-values for single interaction inside iterations
    Kvalues_df = pd.DataFrame(columns=['Kign', 'Kigq'], index=input_streamcomp.index)

    if input_streamcomp.loc['H2O']['Content [mol. fract.]'] == 0:
        phases_num = 2
        interactions_name = ['vapor-liquid']
    else:
        phases_num = 3
        interactions_name = ['vapor-aqueous', 'vapor-liquid']

    for interphase in range(1):#range(phases_num - 1):  # runs calculations once for 'vapor-liquid' systems and twice for 'vapor-liquid-aqueous' systems
        if interphase == 0:
            streamcomp = input_streamcomp
        else:
            streamcomp['Content [mol. fract.]'] = equicomp_df['vapor']

        ### STEP - 1: K's estimation using eq. (3-5) and (3-8)
        Kvalues_init_df = get_initial_Kvalues(comppropDB, streamcomp, P_field, T_field)
        print('-->init K-values:', Kvalues_init_df)
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
            equicomp_df, L = get_equilibrium_composition_v2(streamcomp,
                                                            Kvalues_init_df,
                                                            1e-10,
                                                            False)
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
            zfactors = dict()
            for phase in phasevar_df.index:
                zfactors[phase] = get_zfactor(phasevar_df.loc[phase]['Aj'],
                                                        phasevar_df.loc[phase]['Bj'],
                                                        phase,
                                                        False)
            ### Fugacity coefficients factors calculation
            fugacit_df = get_fugacities(streamcomp,
                                        phasevar_df,
                                        Aijprime_df,
                                        Bijprime_df,
                                        zfactors)
            ### STEP 5: - New set of K-values calculation
            Kvalues_single_df = get_Kvalues(fugacit_df)
            calc_err = Kvalues_comparison(Kvalues_init_df, Kvalues_single_df)
            err_list.append(calc_err)
            Kvalues_init_df['Kij'] = Kvalues_single_df['Kij']
            print('K-values error at iteration {}: {:.3e}'.format(steps, calc_err))
            if steps > steps_limit:
                print('WARNING: K-values did not converged!')
                break
            print('K-values at iteration {}:\n'.format(steps), Kvalues_single_df)

        equicomp_df_final, L_final = get_equilibrium_composition_v2(streamcomp,
                                                                    Kvalues_single_df,
                                                                    1e-10,
                                                                    False)
        print('\n-->equicomp_df:',equicomp_df_final)
        print('-->L = {:.4f}'.format(L_final))

        if abs(equicomp_df_final['vapor'].sum() - 1) >= 1e-3:
            equicomp_df_final['vapor'] = equicomp_df_final['liquid']
            equicomp_df_final['liquid'] = pd.Series([0] * len(equicomp_df_final.index))
            L_final = 0
            if abs(equicomp_df_final['vapor'].sum() - 1) < 1e-3:
                # stop_flag = True
                print('flag switch')

        if abs(equicomp_df_final['vapor'].sum() - 1) >= 0.01:  # temporary measure to handle difficult steams (exmp. Stream6 @ 60 bara, 60 C
            equicomp_df_final['vapor'], equicomp_df_final['liquid'] = equicomp_df_final['liquid'], equicomp_df_final['vapor']
            L_final = 1 - L_final

        if (abs(Kvalues_single_df['Kij'] - 1) < 10 ** -3).all():
            equicomp_df_threephases['vapor'] = streamcomp['Content [mol. fract.]']
            phase_fractions['vapor'], phase_fractions['liquid'], phase_fractions['aqueous'] = 1, 0, 0
            break

        equicomp_df_threephases, phase_fractions, Kvalues_df = redefine_equicomp(equicomp_df_final,
                                                                     equicomp_df_threephases,
                                                                     phase_fractions,
                                                                     L_final,
                                                                     Kvalues_single_df,
                                                                     Kvalues_df,
                                                                     interphase,
                                                                     phases_num)
        if steps <= steps_limit:
            print('Converged in {} iterations for {}\n'.format(steps, interactions_name[interphase]))
    if abs(equicomp_df_threephases['liquid'].sum() - 1) > 0.01 and equicomp_df_threephases['liquid'].sum() > 0:
        equicomp_df_threephases['liquid'] = equicomp_df_threephases['vapor']
        phase_fractions['vapor'], phase_fractions['liquid'] = phase_fractions['liquid'], phase_fractions['vapor']
        equicomp_df_threephases['vapor'] = equicomp_df_threephases['aqueous']
        phase_fractions['vapor'], phase_fractions['aqueous'] = phase_fractions['aqueous'], phase_fractions['vapor']
        equicomp_df_threephases['aqueous'] = pd.Series([0] * len(equicomp_df_threephases.index),
                                                       index= equicomp_df_threephases.index)
        equicomp_df_threephases.loc['H2O', 'aqueous'] = 1
    # return equicomp_df_threephases, phase_fractions, zfactors
    return equicomp_df_final, phase_fractions, zfactors, Kvalues_df


def two_phase_VLE(comppropDB: pd.DataFrame,
                  binarycoefDB: pd.DataFrame,
                  streamcomp: pd.DataFrame,
                  P_field: float,
                  T_field: float,
                  convcrit_K: float,
                  convcrit_L: float,
                  steps_limit: int):
    ### Initial variables setup
    Kvalues_df = pd.DataFrame(columns=['Kij'], index=streamcomp.index) # df used to store K-values for single interaction inside iterations
    phase_fractions = dict()

    ### K-values estimation using Wilson's K-value equation
    Kvalues_init_df = get_initial_Kvalues(comppropDB, streamcomp, P_field, T_field)
    #print('\n-->initial Kvalues', Kvalues_init_df)

    calc_err = 1  # To initiate K-values cycle
    steps = 0
    while calc_err >= convcrit_K:
        if steps > steps_limit:
            print('WARNING: K-values did not converged!')
            break
        steps += 1
        #print('-----step-{}'.format(steps))
        ### Calculate component-dependent variables
        compvar_df = get_compdepvar(comppropDB, streamcomp, T_field)
        #print('\n-->compvar', compvar_df)

        ### Calculate equilibrium copositions
        equicomp_df, L = get_equilibrium_composition_v2(streamcomp, Kvalues_init_df, convcrit_L, False)
        #print('\n-->equicomp-1, L', equicomp_df, L)

        ### Calculate phase-dependent variables
        phasevar_df = get_phasedepvar(equicomp_df, compvar_df, binarycoefDB, P_field, T_field)
        #print('\n-->phasevar', phasevar_df)

        ### Calculate phase-component-dependent variables
        Aijprime_df, Bijprime_df = get_phasecompdepvar(phasevar_df, compvar_df, equicomp_df, binarycoefDB)
        #print('\n-->Aijprime, Bijprime', Aijprime_df, Bijprime_df)

        ### Calculate compressibility factors
        zfactors = dict()
        for phase in phasevar_df.index:
            zfactors[phase] = get_zfactor(phasevar_df.loc[phase]['Aj'], phasevar_df.loc[phase]['Bj'], phase, False)
        #print('\n-->zfactors', zfactors)

        ### Calculate fugacity coefficients factors
        fugacit_df = get_fugacities(streamcomp, phasevar_df, Aijprime_df, Bijprime_df, zfactors)
        #print('\n-->fugacities', fugacit_df)

        ### Calculate the new set of K-values
        Kvalues_df = get_Kvalues(fugacit_df)
        #print('\n-->Kvalues', Kvalues_df)

        ### Calculate error between K-values sets ande rearrange df's for next iteration
        calc_err = Kvalues_comparison(Kvalues_init_df, Kvalues_df)
        Kvalues_init_df['Kij'] = Kvalues_df['Kij']
        print('K-values error at iteration {}: {:.3e}'.format(steps, calc_err))
        #print('K-values at iteration {}:\n'.format(steps), Kvalues_df)

    ### Calculate final VLE equicomposition
    equicomp_df, L = get_equilibrium_composition_v2(streamcomp, Kvalues_df, convcrit_L, False)
    phase_fractions['vapor'] = 1 - L
    phase_fractions['liquid'] = L
    return equicomp_df, phase_fractions, zfactors, Kvalues_df


def flash_calc_PR_EOS(comppropDB: pd.DataFrame,
                      binarycoefDB: pd.DataFrame,
                      input_streamcomp: pd.DataFrame,
                      P_field: float,
                      T_field: float,
                      convcrit_K: float,
                      steps_limit: int):
    ### Initial variables setup
    equicomp_df_3ph = pd.DataFrame(columns=['vapor', 'liquid', 'aqueous'], index=input_streamcomp.index)
    phase_fractions_3ph = dict({'vapor' : 0., 'liquid' : 0., 'aqueous' : 0.})
    Kvalues_df_3ph = pd.DataFrame(columns=['Kign', 'Kigq'], index=input_streamcomp.index)
    zfactors_3ph = dict()

    ### First flash calcuation
    equicomp_df, phase_fractions, zfactors, Kvalues_df = two_phase_VLE(comppropDB,
                                                                       binarycoefDB,
                                                                       input_streamcomp,
                                                                       P_field,
                                                                       T_field,
                                                                       convcrit_K,
                                                                       1e-6,
                                                                       steps_limit)
    ### Fill 3-phase equicomposition df
    if zfactors['vapor'] > zfactors['liquid']:
        equicomp_df_3ph['vapor'] = equicomp_df['vapor']
        equicomp_df_3ph['liquid'] = equicomp_df['liquid']
        phase_fractions_3ph['vapor'] = phase_fractions['vapor']
        phase_fractions_3ph['liquid'] = phase_fractions['liquid']
        zfactors_3ph['vapor'] = zfactors['vapor']
        zfactors_3ph['liquid'] = zfactors['liquid']
        Kvalues_df_3ph['Kign'] = Kvalues_df['Kij']
    elif zfactors['vapor'] < zfactors['liquid']:
        equicomp_df_3ph['vapor'] = equicomp_df['liquid']
        equicomp_df_3ph['liquid'] = equicomp_df['vapor']
        phase_fractions_3ph['vapor'] = phase_fractions['liquid']
        phase_fractions_3ph['liquid'] = phase_fractions['vapor']
        zfactors_3ph['vapor'] = zfactors['liquid']
        zfactors_3ph['liquid'] = zfactors['vapor']
        Kvalues_df_3ph['Kign'] = Kvalues_df['Kij']
    else:
        print('Something went wrong with z-factors...')

    if phase_fractions_3ph['vapor'] <= 1e-5:
        phase_fractions_3ph['vapor'] = 0
        phase_fractions_3ph['liquid'] = 1
        equicomp_df_3ph['vapor'] = 0
    elif phase_fractions_3ph['liquid'] <= 1e-5:
        phase_fractions_3ph['liquid'] = 0
        phase_fractions_3ph['vapor'] = 1
        equicomp_df_3ph['liquid'] = 0
    else:
        ### Stability test for vapor phase and second flash calculation
        #print('---STABILITY TEST---')
        input_streamcomp['Content [mol. fract.]'] = equicomp_df_3ph['vapor']
        # print('-->input', input_streamcomp)
        try:
            equicomp_df, phase_fractions, zfactors, Kvalues_df = two_phase_VLE(comppropDB,
                                                                           binarycoefDB,
                                                                           input_streamcomp,
                                                                           P_field,
                                                                           T_field,
                                                                           convcrit_K,
                                                                           1e-10,
                                                                           steps_limit)

            if phase_fractions['vapor'] > 1e-5:
                ### Phase unstable and should be separated further
                equicomp_df_3ph['aqueous'] = equicomp_df_3ph['liquid']
                equicomp_df_3ph['vapor'] = equicomp_df['vapor']
                equicomp_df_3ph['liquid'] = equicomp_df['liquid']
                phase_fractions_3ph['aqueous'] = phase_fractions_3ph['liquid']
                phase_fractions_3ph['vapor'] = (1 - phase_fractions_3ph['aqueous']) * phase_fractions['vapor']
                phase_fractions_3ph['liquid'] = (1 - phase_fractions_3ph['aqueous']) * phase_fractions['liquid']
                zfactors_3ph['aqueous'] = zfactors_3ph['liquid']
                zfactors_3ph['vapor'] = zfactors['vapor']
                zfactors_3ph['liquid'] = zfactors['liquid']
                Kvalues_df_3ph['kigq'] = Kvalues_df_3ph['Kign']
                Kvalues_df_3ph['kign'] = Kvalues_df['Kij']
        except:
            print('Vapor phase is stable')

    return equicomp_df_3ph, phase_fractions_3ph, zfactors_3ph, Kvalues_df_3ph


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
    return density_mix


def get_phase_densities_actcond(comppropDB: pd.DataFrame,
                                equicomp_df: pd.DataFrame,
                                phase_fractions: dict,  # temporary (?) measure to detect missing phases
                                P_field: float,
                                T_field: float,
                                zfactors: dict):
    densities_dict = dict()
    ### Calculating vapor phase density in actual conditions
    P = UnitsConverter.Pressure.psi_to_kPa(P_field)
    T = UnitsConverter.Temperature.R_to_K(T_field)
    R = 8.31446261815324 # [J/(mole*K)]
    MW = get_phase_molar_weigh(comppropDB, equicomp_df, 'vapor')
    densities_dict['vapor'] = MW / (zfactors['vapor'] * R * T / P)

    ### Caclulation liquid and/or aqueous phase density in actual conditions
    # if equicomp_df['aqueous'].dtype == float:
    #     if np.isnan(equicomp_df['aqueous']).all():
    #         phases = ['liquid']
    #         densities_dict['aqueous'] = np.nan
    # else:
    #     phases = ['liquid', 'aqueous']
    phases = ['liquid', 'aqueous']
    if phase_fractions['liquid'] == 0 or phase_fractions['aqueous'] == 0:
        if phase_fractions['aqueous'] == 0:
            phases = ['liquid']
            densities_dict['aqueous'] = np.nan
        if phase_fractions['liquid'] == 0:
            phases = []
            densities_dict['liquid'] = np.nan

    for phase in phases:
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
                ### Density calculation
                Vs = V0R * (1 - wSRKmix * VdeltaR) * Vasteriskmix  # [m3/kgmol]
                MW = get_phase_molar_weigh(comppropDB, equicomp_df, phase)  # [kg/kgmole]
                densities_dict[phase] = MW / Vs  # [kg/m3]
            else:
                print('WARNING! \tCostald density correlation error. (Tr >= 0.95 for {} phase)'.format(phase))
                densities_dict[phase] = 0.
        else:
            print('\tWARNING! Costald density correlation error. (0.25 <= (Tr={:.2}) >= 1.0 for {} phase)'.format(Tr, phase))
            densities_dict[phase] = 0.
    return densities_dict