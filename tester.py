import pandas as pd
import Interfaces as intrf
import Calculations_v4 as calc
import time
import os
import sys

'''
Automated tests for different gas compositions at different conditions

Tests list:
-
'''

### Input files
cwd = os.getcwd()
comppropDB = pd.read_excel(intrf.get_comppropDB_names(cwd)[0], index_col= 'Name')
binarycoefDB = pd.read_excel(intrf.get_binarycoefDB_names(cwd)[0], index_col= 'index_col')
streamcomp_names_list = intrf.get_streamcomp_names(cwd)
P_list = [500., 2000., 4140., 6000., 10000.]
T_list = [-40., 20., 50., 100.]

### Solver settings
convcrit = 1e-3
steps_limit = 50

start_time = time.perf_counter()
result_df = pd.DataFrame()

tests_passed = 0
tests_failed = 0

### Solver loop
for P in P_list:
    P_field = calc.UnitsConverter.Pressure.kPa_to_psi(P)  # psi
    for T in T_list:
        T_field = calc.UnitsConverter.Temperature.C_to_R(T)  # R
        for streamcomp_name in streamcomp_names_list:
            print('\n\t{} @ {:.1f} bara, {:.1f} C'.format(streamcomp_name, P / 100, T))
            header_df = pd.DataFrame({('{} @ {:.1f} bara, {:.1f} C'.format(streamcomp_name, P / 100, T), '', '')},
                                     columns=(['vapor', 'liquid', 'aqueous']),
                                     index=(['STREAM']))
            input_streamcomp = pd.read_excel(streamcomp_name, index_col='Name')
            '''
            equicomp_df, phase_fract, zfactors, Kvalues = calc.flash_calc_PR_EOS(comppropDB,
                                                                                 binarycoefDB,
                                                                                 input_streamcomp,
                                                                                 P_field,
                                                                                 T_field,
                                                                                 convcrit,
                                                                                 steps_limit)
            print(equicomp_df)
            '''
            try:
                input_streamcomp = pd.read_excel(streamcomp_name, index_col='Name')
                equicomp_df, phase_fract, zfactors, Kvalues = calc.flash_calc_PR_EOS(comppropDB,
                                                                            binarycoefDB,
                                                                            input_streamcomp,
                                                                            P_field,
                                                                            T_field,
                                                                            convcrit,
                                                                            steps_limit)
                print(equicomp_df)
                try:
                    densities = calc.get_phase_densities_actcond(comppropDB,
                                                                 equicomp_df,
                                                                 phase_fract,
                                                                 P_field,
                                                                 T_field,
                                                                 zfactors)
                except:
                    densities = ({'vapor': 0, 'liquid': 0, 'aqueous': 0})

                phase_fract_df = pd.DataFrame({(phase_fract['vapor'], phase_fract['liquid'], phase_fract['aqueous'])},
                                              columns= (['vapor', 'liquid', 'aqueous']),
                                              index= (['Phase Fraction']))
                densities_df = pd.DataFrame({(densities['vapor'], densities['liquid'], densities['aqueous'])},
                                              columns= (['vapor', 'liquid', 'aqueous']),
                                              index= (['Density [kg/m3]']))
                tests_passed += 1
            except:
                equicomp_df = pd.DataFrame(columns= (['vapor', 'liquid', 'aqueous']), index= input_streamcomp.index)
                phase_fract_df = pd.DataFrame(columns= (['vapor', 'liquid', 'aqueous']), index= (['Phase Fraction']))
                densities_df = pd.DataFrame(columns= (['vapor', 'liquid', 'aqueous']), index= (['Densities']))
                tests_failed += 1
            result_df = result_df.append(header_df)
            result_df = result_df.append(equicomp_df)
            result_df = result_df.append(phase_fract_df)
            result_df = result_df.append(densities_df)
print(result_df)
print('Execution time: {:.3f} s'.format(time.perf_counter() - start_time))
print('\nTests performed: {}\n\tfailed: {}\n\tpassed: {}'.format(len(P_list) * len(T_list) * len(streamcomp_names_list),
                                                                 tests_failed, tests_passed))


try:
    writer = pd.ExcelWriter('tester_output.xlsx')
    result_df.to_excel(writer)
    writer.save()
except:
    print('Failed to write down results')

print('Open output file? ("y" - yes, "n" - no)')
confirm = str(input())
if confirm == 'y':
    try:
        os.system('start excel.exe "%s//tester_output.xlsx"' % (sys.path[0], ))
    except:
        print('Failed to open file')