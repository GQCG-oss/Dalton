import functionals_mu as mu_func
import functionals as func
import functionals_special as spec_func
from sympy import symbols, cse, fcode, Rational, sqrt, exp, ln, Max
import print_functional_to_DALTON as dalprint
import print_functional_to_DALTON_TPSSc_special_case as dalprint_spec
rho_a = symbols('rho_a', real=True, nonnegative=True)
gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
tau_a = symbols('tau_a', real=True)
rho_c = symbols('rho_c', real=True, nonnegative=True)
rho_s = symbols('rho_s', real=True)
gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
gamma_cs = symbols('gamma_cs', real=True)
tau_c, tau_s = symbols('tau_c tau_s', real=True)
mu = symbols('mu', real=True, nonnegative=True)


parameters = {}
parameterfile = open("constants.txt","r")
for line in parameterfile:
    parameters[line.split("=")[0].strip()] = str(line.split("=")[1])


#write_list = ["LDA_ERF_exchange","PBE_ERFGWS_exchange","TPSS_ERFGWS_exchange","PW92_ERF_correlation",
#              "PBE_ERFGWS_correlation","TPSS_ERFGWS_correlation","PBE_nomu_correlation","VWN5_ERF_correlation",
#              "wPBE_exchange","VWN5_nomu_correlation"]
write_list = ["VWN5_nomu_correlation"]

if "LDA_ERF_exchange" in write_list:
    out_file = open("../srfunctionals/LDA_ERF_exchange.F","w+")
    # ######################################################################
    # LDA exchange
    # ######################################################################
    # rho_c = 2*rho_a
    E_case_1 = mu_func.LDAx_mu_case_1(parameters)*2*rho_a
    d1E_rhoa_case_1 = E_case_1.diff(rho_a)
    d2E_rhoa2_case_1 = d1E_rhoa_case_1.diff(rho_a)
    
    E_case_2 = mu_func.LDAx_mu_case_2(parameters)*2*rho_a
    d1E_rhoa_case_2 = E_case_2.diff(rho_a)
    d2E_rhoa2_case_2 = d1E_rhoa_case_2.diff(rho_a)
    
    E_case_3 = mu_func.LDAx_mu_case_3(parameters)*2*rho_a
    d1E_rhoa_case_3 = E_case_3.diff(rho_a)
    d2E_rhoa2_case_3 = d1E_rhoa_case_3.diff(rho_a)		
    
    Kernel = [E_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_LDA_ERF_case_1", ["rho_a"],["Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_LDA_ERF_case_2", ["rho_a"],["Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_LDA_ERF_case_3", ["rho_a"],["Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1, d1E_rhoa_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_LDA_ERF_case_1", ["rho_a"],["Ea","d1Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2, d1E_rhoa_case_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_LDA_ERF_case_2", ["rho_a"],["Ea","d1Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3, d1E_rhoa_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_LDA_ERF_case_3", ["rho_a"],["Ea","d1Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1, d1E_rhoa_case_1, d2E_rhoa2_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_LDA_ERF_case_1", ["rho_a"],["Ea","d1Ea","d2Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2, d1E_rhoa_case_2, d2E_rhoa2_case_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_LDA_ERF_case_2", ["rho_a"],["Ea","d1Ea","d2Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3, d1E_rhoa_case_3, d2E_rhoa2_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_LDA_ERF_case_3", ["rho_a"],["Ea","d1Ea","d2Ea"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()


if "PBE_ERFGWS_exchange" in write_list:
    out_file = open("../srfunctionals/PBE_ERFGWS_exchange.F","w+")
    # ######################################################################
    # PBE exchange
    # ######################################################################
    # rho_c = 2*rho_a
    E_case_1 = mu_func.PBEx_mu_case_1(parameters)*2*rho_a
    d1E_rhoa_case_1 = E_case_1.diff(rho_a)
    d1E_gammaaa_case_1 = E_case_1.diff(gamma_aa)
    d2E_rhoa2_case_1 = d1E_rhoa_case_1.diff(rho_a)
    d2E_rhoagammaaa_case_1 = d1E_rhoa_case_1.diff(gamma_aa)
    d2E_gammaaa2_case_1 = d1E_gammaaa_case_1.diff(gamma_aa)
    
    E_case_2_1 = mu_func.PBEx_mu_case_2_1(parameters)*2*rho_a
    d1E_rhoa_case_2_1 = E_case_2_1.diff(rho_a)
    d1E_gammaaa_case_2_1 = E_case_2_1.diff(gamma_aa)
    d2E_rhoa2_case_2_1 = d1E_rhoa_case_2_1.diff(rho_a)
    d2E_rhoagammaaa_case_2_1 = d1E_rhoa_case_2_1.diff(gamma_aa)
    d2E_gammaaa2_case_2_1 = d1E_gammaaa_case_2_1.diff(gamma_aa)
    
    E_case_2_2 = mu_func.PBEx_mu_case_2_2(parameters)*2*rho_a
    d1E_rhoa_case_2_2 = E_case_2_2.diff(rho_a)
    d1E_gammaaa_case_2_2 = E_case_2_2.diff(gamma_aa)
    d2E_rhoa2_case_2_2 = d1E_rhoa_case_2_2.diff(rho_a)
    d2E_rhoagammaaa_case_2_2 = d1E_rhoa_case_2_2.diff(gamma_aa)
    d2E_gammaaa2_case_2_2 = d1E_gammaaa_case_2_2.diff(gamma_aa)
    
    E_case_2_3 = mu_func.PBEx_mu_case_2_3(parameters)*2*rho_a
    d1E_rhoa_case_2_3 = E_case_2_3.diff(rho_a)
    d1E_gammaaa_case_2_3 = E_case_2_3.diff(gamma_aa)
    d2E_rhoa2_case_2_3 = d1E_rhoa_case_2_3.diff(rho_a)
    d2E_rhoagammaaa_case_2_3 = d1E_rhoa_case_2_3.diff(gamma_aa)
    d2E_gammaaa2_case_2_3 = d1E_gammaaa_case_2_3.diff(gamma_aa)
    
    E_case_3 = mu_func.PBEx_mu_case_3(parameters)*2*rho_a
    d1E_rhoa_case_3 = E_case_3.diff(rho_a)
    d1E_gammaaa_case_3 = E_case_3.diff(gamma_aa)
    d2E_rhoa2_case_3 = d1E_rhoa_case_3.diff(rho_a)
    d2E_rhoagammaaa_case_3 = d1E_rhoa_case_3.diff(gamma_aa)
    d2E_gammaaa2_case_3 = d1E_gammaaa_case_3.diff(gamma_aa)
    
    Kernel = [E_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_PBE_ERF_case_1", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_PBE_ERF_case_2_1", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_PBE_ERF_case_2_2", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_PBE_ERF_case_2_3", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_PBE_ERF_case_3", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1, d1E_rhoa_case_1, d1E_gammaaa_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_PBE_ERF_case_1", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1, d1E_rhoa_case_2_1, d1E_gammaaa_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_PBE_ERF_case_2_1", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2, d1E_rhoa_case_2_2, d1E_gammaaa_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_PBE_ERF_case_2_2", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3, d1E_rhoa_case_2_3, d1E_gammaaa_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_PBE_ERF_case_2_3", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3, d1E_rhoa_case_3, d1E_gammaaa_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_PBE_ERF_case_3", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1, d1E_rhoa_case_1, d1E_gammaaa_case_1, d2E_rhoa2_case_1, d2E_rhoagammaaa_case_1, d2E_gammaaa2_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_PBE_ERF_case_1", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1, d1E_rhoa_case_2_1, d1E_gammaaa_case_2_1, d2E_rhoa2_case_2_1, d2E_rhoagammaaa_case_2_1, d2E_gammaaa2_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_PBE_ERF_case_2_1", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2, d1E_rhoa_case_2_2, d1E_gammaaa_case_2_2, d2E_rhoa2_case_2_2, d2E_rhoagammaaa_case_2_2, d2E_gammaaa2_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_PBE_ERF_case_2_2", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3, d1E_rhoa_case_2_3, d1E_gammaaa_case_2_3, d2E_rhoa2_case_2_3, d2E_rhoagammaaa_case_2_3, d2E_gammaaa2_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_PBE_ERF_case_2_3", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3, d1E_rhoa_case_3, d1E_gammaaa_case_3, d2E_rhoa2_case_3, d2E_rhoagammaaa_case_3, d2E_gammaaa2_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_PBE_ERF_case_3", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()


if "TPSS_ERFGWS_exchange" in write_list:
    out_file = open("../srfunctionals/TPSS_ERFGWS_exchange.F","w+")
    # ######################################################################
    # TPSS exchange
    # ######################################################################
    # rho_c = 2*rho_a
    E_case_1 = mu_func.TPSSx_mu_case_1(parameters)*2*rho_a
    d1E_rhoa_case_1 = E_case_1.diff(rho_a)
    d1E_gammaaa_case_1 = E_case_1.diff(gamma_aa)
    d1E_taua_case_1 = E_case_1.diff(tau_a)
    d2E_rhoa2_case_1 = d1E_rhoa_case_1.diff(rho_a)
    d2E_rhoagammaaa_case_1 = d1E_rhoa_case_1.diff(gamma_aa)
    d2E_gammaaa2_case_1 = d1E_gammaaa_case_1.diff(gamma_aa)
    d2E_rhoataua_case_1 = d1E_rhoa_case_1.diff(tau_a)
    d2E_gammaaataua_case_1 = d1E_gammaaa_case_1.diff(tau_a)
    d2E_taua2_case_1 = d1E_taua_case_1.diff(tau_a)
    
    E_case_2_1 = mu_func.TPSSx_mu_case_2_1(parameters)*2*rho_a
    d1E_rhoa_case_2_1 = E_case_2_1.diff(rho_a)
    d1E_gammaaa_case_2_1 = E_case_2_1.diff(gamma_aa)
    d1E_taua_case_2_1 = E_case_2_1.diff(tau_a)
    d2E_rhoa2_case_2_1 = d1E_rhoa_case_2_1.diff(rho_a)
    d2E_rhoagammaaa_case_2_1 = d1E_rhoa_case_2_1.diff(gamma_aa)
    d2E_gammaaa2_case_2_1 = d1E_gammaaa_case_2_1.diff(gamma_aa)
    d2E_rhoataua_case_2_1 = d1E_rhoa_case_2_1.diff(tau_a)
    d2E_gammaaataua_case_2_1 = d1E_gammaaa_case_2_1.diff(tau_a)
    d2E_taua2_case_2_1 = d1E_taua_case_2_1.diff(tau_a)
    
    E_case_2_2 = mu_func.TPSSx_mu_case_2_2(parameters)*2*rho_a
    d1E_rhoa_case_2_2 = E_case_2_2.diff(rho_a)
    d1E_gammaaa_case_2_2 = E_case_2_2.diff(gamma_aa)
    d1E_taua_case_2_2 = E_case_2_2.diff(tau_a)
    d2E_rhoa2_case_2_2 = d1E_rhoa_case_2_2.diff(rho_a)
    d2E_rhoagammaaa_case_2_2 = d1E_rhoa_case_2_2.diff(gamma_aa)
    d2E_gammaaa2_case_2_2 = d1E_gammaaa_case_2_2.diff(gamma_aa)
    d2E_rhoataua_case_2_2 = d1E_rhoa_case_2_2.diff(tau_a)
    d2E_gammaaataua_case_2_2 = d1E_gammaaa_case_2_2.diff(tau_a)
    d2E_taua2_case_2_2 = d1E_taua_case_2_2.diff(tau_a)
    
    E_case_2_3 = mu_func.TPSSx_mu_case_2_3(parameters)*2*rho_a
    d1E_rhoa_case_2_3 = E_case_2_3.diff(rho_a)
    d1E_gammaaa_case_2_3 = E_case_2_3.diff(gamma_aa)
    d1E_taua_case_2_3 = E_case_2_3.diff(tau_a)
    d2E_rhoa2_case_2_3 = d1E_rhoa_case_2_3.diff(rho_a)
    d2E_rhoagammaaa_case_2_3 = d1E_rhoa_case_2_3.diff(gamma_aa)
    d2E_gammaaa2_case_2_3 = d1E_gammaaa_case_2_3.diff(gamma_aa)
    d2E_rhoataua_case_2_3 = d1E_rhoa_case_2_3.diff(tau_a)
    d2E_gammaaataua_case_2_3 = d1E_gammaaa_case_2_3.diff(tau_a)
    d2E_taua2_case_2_3 = d1E_taua_case_2_3.diff(tau_a)
    
    E_case_3 = mu_func.TPSSx_mu_case_3(parameters)*2*rho_a
    d1E_rhoa_case_3 = E_case_3.diff(rho_a)
    d1E_gammaaa_case_3 = E_case_3.diff(gamma_aa)
    d1E_taua_case_3 = E_case_3.diff(tau_a)
    d2E_rhoa2_case_3 = d1E_rhoa_case_3.diff(rho_a)
    d2E_rhoagammaaa_case_3 = d1E_rhoa_case_3.diff(gamma_aa)
    d2E_gammaaa2_case_3 = d1E_gammaaa_case_3.diff(gamma_aa)
    d2E_rhoataua_case_3 = d1E_rhoa_case_3.diff(tau_a)
    d2E_gammaaataua_case_3 = d1E_gammaaa_case_3.diff(tau_a)
    d2E_taua2_case_3 = d1E_taua_case_3.diff(tau_a)
    
    Kernel = [E_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_TPSS_ERF_case_1", ["rho_a","gamma_aa","tau_a"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_TPSS_ERF_case_2_1", ["rho_a","gamma_aa","tau_a"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_TPSS_ERF_case_2_2", ["rho_a","gamma_aa","tau_a"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_TPSS_ERF_case_2_3", ["rho_a","gamma_aa","tau_a"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_TPSS_ERF_case_3", ["rho_a","gamma_aa","tau_a"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1,d1E_rhoa_case_1,d1E_gammaaa_case_1,d1E_taua_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_TPSS_ERF_case_1", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1,d1E_rhoa_case_2_1,d1E_gammaaa_case_2_1,d1E_taua_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_TPSS_ERF_case_2_1", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2,d1E_rhoa_case_2_2,d1E_gammaaa_case_2_2,d1E_taua_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_TPSS_ERF_case_2_2", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3,d1E_rhoa_case_2_3,d1E_gammaaa_case_2_3,d1E_taua_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_TPSS_ERF_case_2_3", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3,d1E_rhoa_case_3,d1E_gammaaa_case_3,d1E_taua_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_TPSS_ERF_case_3", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_1,d1E_rhoa_case_1,d1E_gammaaa_case_1,d1E_taua_case_1,
            d2E_rhoa2_case_1,d2E_rhoagammaaa_case_1,d2E_gammaaa2_case_1,d2E_rhoataua_case_1,
            d2E_gammaaataua_case_1,d2E_taua2_case_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_TPSS_ERF_case_1", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_1,d1E_rhoa_case_2_1,d1E_gammaaa_case_2_1,d1E_taua_case_2_1,
            d2E_rhoa2_case_2_1,d2E_rhoagammaaa_case_2_1,d2E_gammaaa2_case_2_1,d2E_rhoataua_case_2_1,
            d2E_gammaaataua_case_2_1,d2E_taua2_case_2_1]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_TPSS_ERF_case_2_1", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_2,d1E_rhoa_case_2_2,d1E_gammaaa_case_2_2,d1E_taua_case_2_2,
            d2E_rhoa2_case_2_2,d2E_rhoagammaaa_case_2_2,d2E_gammaaa2_case_2_2,d2E_rhoataua_case_2_2,
            d2E_gammaaataua_case_2_2,d2E_taua2_case_2_2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_TPSS_ERF_case_2_2", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_2_3,d1E_rhoa_case_2_3,d1E_gammaaa_case_2_3,d1E_taua_case_2_3,
            d2E_rhoa2_case_2_3,d2E_rhoagammaaa_case_2_3,d2E_gammaaa2_case_2_3,d2E_rhoataua_case_2_3,
            d2E_gammaaataua_case_2_3,d2E_taua2_case_2_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_TPSS_ERF_case_2_3", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E_case_3,d1E_rhoa_case_3,d1E_gammaaa_case_3,d1E_taua_case_3,
            d2E_rhoa2_case_3,d2E_rhoagammaaa_case_3,d2E_gammaaa2_case_3,d2E_rhoataua_case_3,
            d2E_gammaaataua_case_3,d2E_taua2_case_3]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_TPSS_ERF_case_3", ["rho_a","gamma_aa","tau_a"],
                                    ["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()


if "PW92_ERF_correlation" in write_list:
    out_file = open("../srfunctionals/PW92_ERF_correlation.F","w+")
    # ######################################################################
    # PW92 correlation, no-spin
    # ######################################################################
    E = mu_func.PW92c_mu(parameters).subs({rho_s: 0})*rho_c
    d1E_rhoc = E.diff(rho_c)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_PW92_ERF", ["rho_c"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_PW92_ERF", ["rho_c"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d2E_rhoc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_PW92_ERF", ["rho_c"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # PW92 correlation, spin
    # ######################################################################
    E = mu_func.PW92c_mu(parameters)*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_rhos = E.diff(rho_s)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocrhos = d1E_rhoc.diff(rho_s)
    d2E_rhos2 = d1E_rhos.diff(rho_s)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_SPIN_PW92_ERF", ["rho_c","rho_s"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_SPIN_PW92_ERF", ["rho_c","rho_s"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d2E_rhoc2,d2E_rhocrhos,d2E_rhos2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_SPIN_PW92_ERF", ["rho_c","rho_s"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()


if "VWN5_ERF_correlation" in write_list:
    out_file = open("../srfunctionals/VWN5_ERF_correlation.F","w+")
    # ######################################################################
    # VWN5 correlation, no-spin
    # ######################################################################
    E = mu_func.VWN5c_mu(parameters).subs({rho_s: 0})*rho_c
    d1E_rhoc = E.diff(rho_c)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_VWN5_ERF", ["rho_c"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_VWN5_ERF", ["rho_c"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d2E_rhoc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_VWN5_ERF", ["rho_c"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # VWN5 correlation, spin
    # ######################################################################
    E = mu_func.VWN5c_mu(parameters)*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_rhos = E.diff(rho_s)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocrhos = d1E_rhoc.diff(rho_s)
    d2E_rhos2 = d1E_rhos.diff(rho_s)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_SPIN_VWN5_ERF", ["rho_c","rho_s"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_SPIN_VWN5_ERF", ["rho_c","rho_s"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d2E_rhoc2,d2E_rhocrhos,d2E_rhos2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_SPIN_VWN5_ERF", ["rho_c","rho_s"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()



if "PBE_ERFGWS_correlation" in write_list:
    out_file = open("../srfunctionals/PBE_ERFGWS_correlation.F","w+")
    # ######################################################################
    # PBE correlatoin, no-spin
    # ######################################################################
    E = mu_func.PBEc_mu(parameters).subs({rho_s: 0})*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_gammacc = E.diff(gamma_cc)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
    d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_gammacc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,3],[1,4,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # PBE correlatoin, spin
    # ######################################################################
    E = mu_func.PBEc_mu(parameters)*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_rhos = E.diff(rho_s)
    d1E_gammacc = E.diff(gamma_cc)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocrhos = d1E_rhoc.diff(rho_s)
    d2E_rhos2 = d1E_rhos.diff(rho_s)
    d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
    d2E_rhosgammacc = d1E_rhos.diff(gamma_cc)
    d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_SPIN_PBE_GWS_ERF", ["rho_c","rho_s","gamma_cc"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d1E_gammacc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_SPIN_PBE_GWS_ERF", ["rho_c","rho_s","gamma_cc"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d1E_gammacc,d2E_rhoc2,d2E_rhocrhos,d2E_rhos2,d2E_rhocgammacc,d2E_rhosgammacc,d2E_gammacc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_SPIN_PBE_GWS_ERF", ["rho_c","rho_s","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()


if "TPSS_ERFGWS_correlation" in write_list:
    out_file = open("../srfunctionals/TPSS_ERFGWS_correlation.F","w+")
    # ######################################################################
    # TPSS correlation, no-spin
    # ######################################################################
    E_case_1 = mu_func.TPSSc_mu_case_1(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0})*rho_c
    d1E_rhoc_case_1 = E_case_1.diff(rho_c)
    d1E_gammacc_case_1 = E_case_1.diff(gamma_cc)
    d1E_tauc_case_1 = E_case_1.diff(tau_c)
    d2E_rhoc2_case_1 = d1E_rhoc_case_1.diff(rho_c)
    d2E_rhocgammacc_case_1 = d1E_rhoc_case_1.diff(gamma_cc)
    d2E_gammacc2_case_1 = d1E_gammacc_case_1.diff(gamma_cc)
    d2E_rhoctauc_case_1 = d1E_rhoc_case_1.diff(tau_c)
    d2E_gammacctauc_case_1 = d1E_gammacc_case_1.diff(tau_c)
    d2E_tauc2_case_1 = d1E_tauc_case_1.diff(tau_c)
    
    E_case_2 = mu_func.TPSSc_mu_case_2(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0})*rho_c
    d1E_rhoc_case_2 = E_case_2.diff(rho_c)
    d1E_gammacc_case_2 = E_case_2.diff(gamma_cc)
    d1E_tauc_case_2 = E_case_2.diff(tau_c)
    d2E_rhoc2_case_2 = d1E_rhoc_case_2.diff(rho_c)
    d2E_rhocgammacc_case_2 = d1E_rhoc_case_2.diff(gamma_cc)
    d2E_gammacc2_case_2 = d1E_gammacc_case_2.diff(gamma_cc)
    d2E_rhoctauc_case_2 = d1E_rhoc_case_2.diff(tau_c)
    d2E_gammacctauc_case_2 = d1E_gammacc_case_2.diff(tau_c)
    d2E_tauc2_case_2 = d1E_tauc_case_2.diff(tau_c)
    
    E_case_3 = mu_func.TPSSc_mu_case_3(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0})*rho_c
    d1E_rhoc_case_3 = E_case_3.diff(rho_c)
    d1E_gammacc_case_3 = E_case_3.diff(gamma_cc)
    d1E_tauc_case_3 = E_case_3.diff(tau_c)
    d2E_rhoc2_case_3 = d1E_rhoc_case_3.diff(rho_c)
    d2E_rhocgammacc_case_3 = d1E_rhoc_case_3.diff(gamma_cc)
    d2E_gammacc2_case_3 = d1E_gammacc_case_3.diff(gamma_cc)
    d2E_rhoctauc_case_3 = d1E_rhoc_case_3.diff(tau_c)
    d2E_gammacctauc_case_3 = d1E_gammacc_case_3.diff(tau_c)
    d2E_tauc2_case_3 = d1E_tauc_case_3.diff(tau_c)
    
    E_case_4 = mu_func.TPSSc_mu_case_4(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0})*rho_c
    d1E_rhoc_case_4 = E_case_4.diff(rho_c)
    d1E_gammacc_case_4 = E_case_4.diff(gamma_cc)
    d1E_tauc_case_4 = E_case_4.diff(tau_c)
    d2E_rhoc2_case_4 = d1E_rhoc_case_4.diff(rho_c)
    d2E_rhocgammacc_case_4 = d1E_rhoc_case_4.diff(gamma_cc)
    d2E_gammacc2_case_4 = d1E_gammacc_case_4.diff(gamma_cc)
    d2E_rhoctauc_case_4 = d1E_rhoc_case_4.diff(tau_c)
    d2E_gammacctauc_case_4 = d1E_gammacc_case_4.diff(tau_c)
    d2E_tauc2_case_4 = d1E_tauc_case_4.diff(tau_c)
    
    Kernel = [func.PBEc(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_alpha_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_beta_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}),
            E_case_1, 
            E_case_2, 
            E_case_3, 
            E_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,
                1,
                1,
                1]
    diff_idx = [[0],
                [0],
                [0],
                [0]]
    dalprint_spec.dalton_functional_printer(Kernel, "ESRC_TPSS_GWS_ERF", ["rho_c","gamma_cc","tau_c"],
                                    ["E",
                                        "E",
                                        "E",
                                        "E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [func.PBEc(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_alpha_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_beta_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}),
            E_case_1, d1E_rhoc_case_1, d1E_gammacc_case_1, d1E_tauc_case_1,
            E_case_2, d1E_rhoc_case_2, d1E_gammacc_case_2, d1E_tauc_case_2, 
            E_case_3, d1E_rhoc_case_3, d1E_gammacc_case_3, d1E_tauc_case_3, 
            E_case_4, d1E_rhoc_case_4, d1E_gammacc_case_4, d1E_tauc_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,
                1,3,
                1,3,
                1,3]
    diff_idx = [[0],[1,3,6],
                [0],[1,3,6],
                [0],[1,3,6],
                [0],[1,3,6]]
    dalprint_spec.dalton_functional_printer(Kernel, "D1ESRC_TPSS_GWS_ERF", ["rho_c","gamma_cc","tau_c"],
                                    ["E","d1E",
                                        "E","d1E",
                                        "E","d1E",
                                        "E","d1E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [func.PBEc(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_alpha_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}), 
            spec_func.PBEc_beta_replaced(parameters).subs({rho_s: 0, gamma_ss: 0, gamma_cs: 0}),
            E_case_1, d1E_rhoc_case_1, d1E_gammacc_case_1, d1E_tauc_case_1, d2E_rhoc2_case_1, d2E_rhocgammacc_case_1, d2E_gammacc2_case_1, d2E_rhoctauc_case_1, d2E_gammacctauc_case_1, d2E_tauc2_case_1,
            E_case_2, d1E_rhoc_case_2, d1E_gammacc_case_2, d1E_tauc_case_2, d2E_rhoc2_case_2, d2E_rhocgammacc_case_2, d2E_gammacc2_case_2, d2E_rhoctauc_case_2, d2E_gammacctauc_case_2, d2E_tauc2_case_2, 
            E_case_3, d1E_rhoc_case_3, d1E_gammacc_case_3, d1E_tauc_case_3, d2E_rhoc2_case_3, d2E_rhocgammacc_case_3, d2E_gammacc2_case_3, d2E_rhoctauc_case_3, d2E_gammacctauc_case_3, d2E_tauc2_case_3, 
            E_case_4, d1E_rhoc_case_4, d1E_gammacc_case_4, d1E_tauc_case_4, d2E_rhoc2_case_4, d2E_rhocgammacc_case_4, d2E_gammacc2_case_4, d2E_rhoctauc_case_4, d2E_gammacctauc_case_4, d2E_tauc2_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6,
                1,3,6,
                1,3,6,
                1,3,6]
    diff_idx = [[0],[1,3,6],[1,4,6,16,18,21],
                [0],[1,3,6],[1,4,6,16,18,21],
                [0],[1,3,6],[1,4,6,16,18,21],
                [0],[1,3,6],[1,4,6,16,18,21]]
    dalprint_spec.dalton_functional_printer(Kernel, "D2ESRC_TPSS_GWS_ERF", ["rho_c","gamma_cc","tau_c"],
                                    ["E","d1E","d2E",
                                        "E","d1E","d2E",
                                        "E","d1E","d2E",
                                        "E","d1E","d2E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # TPSS correlation, spin
    # ######################################################################
    E_case_1 = mu_func.TPSSc_mu_case_1(parameters)*rho_c
    d1E_rhoc_case_1 = E_case_1.diff(rho_c)
    d1E_rhos_case_1 = E_case_1.diff(rho_s)
    d1E_gammacc_case_1 = E_case_1.diff(gamma_cc)
    d1E_gammass_case_1 = E_case_1.diff(gamma_ss)
    d1E_gammacs_case_1 = E_case_1.diff(gamma_cs)
    d1E_tauc_case_1 = E_case_1.diff(tau_c)
    d2E_rhoc2_case_1 = d1E_rhoc_case_1.diff(rho_c)
    d2E_rhocrhos_case_1 = d1E_rhoc_case_1.diff(rho_s)
    d2E_rhos2_case_1 = d1E_rhos_case_1.diff(rho_s)
    d2E_rhocgammacc_case_1 = d1E_rhoc_case_1.diff(gamma_cc)
    d2E_rhosgammacc_case_1 = d1E_rhos_case_1.diff(gamma_cc)
    d2E_gammacc2_case_1 = d1E_gammacc_case_1.diff(gamma_cc)
    d2E_rhocgammass_case_1 = d1E_rhoc_case_1.diff(gamma_ss)
    d2E_rhosgammass_case_1 = d1E_rhos_case_1.diff(gamma_ss)
    d2E_gammaccgammass_case_1 = d1E_gammacc_case_1.diff(gamma_ss)
    d2E_gammass2_case_1 = d1E_gammass_case_1.diff(gamma_ss)
    d2E_rhocgammacs_case_1 = d1E_rhoc_case_1.diff(gamma_cs)
    d2E_rhosgammacs_case_1 = d1E_rhos_case_1.diff(gamma_cs)
    d2E_gammaccgammacs_case_1 = d1E_gammacc_case_1.diff(gamma_cs)
    d2E_gammassgammacs_case_1 = d1E_gammass_case_1.diff(gamma_cs)
    d2E_gammacs2_case_1 = d1E_gammacs_case_1.diff(gamma_cs)
    d2E_rhoctauc_case_1 = d1E_rhoc_case_1.diff(tau_c)
    d2E_rhostauc_case_1 = d1E_rhos_case_1.diff(tau_c)
    d2E_gammacctauc_case_1 = d1E_gammacc_case_1.diff(tau_c)
    d2E_gammasstauc_case_1 = d1E_gammass_case_1.diff(tau_c)
    d2E_gammacstauc_case_1 = d1E_gammacs_case_1.diff(tau_c)
    d2E_tauc2_case_1 = d1E_tauc_case_1.diff(tau_c)
    
    E_case_2 = mu_func.TPSSc_mu_case_2(parameters)*rho_c
    d1E_rhoc_case_2 = E_case_2.diff(rho_c)
    d1E_rhos_case_2 = E_case_2.diff(rho_s)
    d1E_gammacc_case_2 = E_case_2.diff(gamma_cc)
    d1E_gammass_case_2 = E_case_2.diff(gamma_ss)
    d1E_gammacs_case_2 = E_case_2.diff(gamma_cs)
    d1E_tauc_case_2 = E_case_2.diff(tau_c)
    d2E_rhoc2_case_2 = d1E_rhoc_case_2.diff(rho_c)
    d2E_rhocrhos_case_2 = d1E_rhoc_case_2.diff(rho_s)
    d2E_rhos2_case_2 = d1E_rhos_case_2.diff(rho_s)
    d2E_rhocgammacc_case_2 = d1E_rhoc_case_2.diff(gamma_cc)
    d2E_rhosgammacc_case_2 = d1E_rhos_case_2.diff(gamma_cc)
    d2E_gammacc2_case_2 = d1E_gammacc_case_2.diff(gamma_cc)
    d2E_rhocgammass_case_2 = d1E_rhoc_case_2.diff(gamma_ss)
    d2E_rhosgammass_case_2 = d1E_rhos_case_2.diff(gamma_ss)
    d2E_gammaccgammass_case_2 = d1E_gammacc_case_2.diff(gamma_ss)
    d2E_gammass2_case_2 = d1E_gammass_case_2.diff(gamma_ss)
    d2E_rhocgammacs_case_2 = d1E_rhoc_case_2.diff(gamma_cs)
    d2E_rhosgammacs_case_2 = d1E_rhos_case_2.diff(gamma_cs)
    d2E_gammaccgammacs_case_2 = d1E_gammacc_case_2.diff(gamma_cs)
    d2E_gammassgammacs_case_2 = d1E_gammass_case_2.diff(gamma_cs)
    d2E_gammacs2_case_2 = d1E_gammacs_case_2.diff(gamma_cs)
    d2E_rhoctauc_case_2 = d1E_rhoc_case_2.diff(tau_c)
    d2E_rhostauc_case_2 = d1E_rhos_case_2.diff(tau_c)
    d2E_gammacctauc_case_2 = d1E_gammacc_case_2.diff(tau_c)
    d2E_gammasstauc_case_2 = d1E_gammass_case_2.diff(tau_c)
    d2E_gammacstauc_case_2 = d1E_gammacs_case_2.diff(tau_c)
    d2E_tauc2_case_2 = d1E_tauc_case_2.diff(tau_c)
    
    E_case_3 = mu_func.TPSSc_mu_case_3(parameters)*rho_c
    d1E_rhoc_case_3 = E_case_3.diff(rho_c)
    d1E_rhos_case_3 = E_case_3.diff(rho_s)
    d1E_gammacc_case_3 = E_case_3.diff(gamma_cc)
    d1E_gammass_case_3 = E_case_3.diff(gamma_ss)
    d1E_gammacs_case_3 = E_case_3.diff(gamma_cs)
    d1E_tauc_case_3 = E_case_3.diff(tau_c)
    d2E_rhoc2_case_3 = d1E_rhoc_case_3.diff(rho_c)
    d2E_rhocrhos_case_3 = d1E_rhoc_case_3.diff(rho_s)
    d2E_rhos2_case_3 = d1E_rhos_case_3.diff(rho_s)
    d2E_rhocgammacc_case_3 = d1E_rhoc_case_3.diff(gamma_cc)
    d2E_rhosgammacc_case_3 = d1E_rhos_case_3.diff(gamma_cc)
    d2E_gammacc2_case_3 = d1E_gammacc_case_3.diff(gamma_cc)
    d2E_rhocgammass_case_3 = d1E_rhoc_case_3.diff(gamma_ss)
    d2E_rhosgammass_case_3 = d1E_rhos_case_3.diff(gamma_ss)
    d2E_gammaccgammass_case_3 = d1E_gammacc_case_3.diff(gamma_ss)
    d2E_gammass2_case_3 = d1E_gammass_case_3.diff(gamma_ss)
    d2E_rhocgammacs_case_3 = d1E_rhoc_case_3.diff(gamma_cs)
    d2E_rhosgammacs_case_3 = d1E_rhos_case_3.diff(gamma_cs)
    d2E_gammaccgammacs_case_3 = d1E_gammacc_case_3.diff(gamma_cs)
    d2E_gammassgammacs_case_3 = d1E_gammass_case_3.diff(gamma_cs)
    d2E_gammacs2_case_3 = d1E_gammacs_case_3.diff(gamma_cs)
    d2E_rhoctauc_case_3 = d1E_rhoc_case_3.diff(tau_c)
    d2E_rhostauc_case_3 = d1E_rhos_case_3.diff(tau_c)
    d2E_gammacctauc_case_3 = d1E_gammacc_case_3.diff(tau_c)
    d2E_gammasstauc_case_3 = d1E_gammass_case_3.diff(tau_c)
    d2E_gammacstauc_case_3 = d1E_gammacs_case_3.diff(tau_c)
    d2E_tauc2_case_3 = d1E_tauc_case_3.diff(tau_c)
    
    E_case_4 = mu_func.TPSSc_mu_case_4(parameters)*rho_c
    d1E_rhoc_case_4 = E_case_4.diff(rho_c)
    d1E_rhos_case_4 = E_case_4.diff(rho_s)
    d1E_gammacc_case_4 = E_case_4.diff(gamma_cc)
    d1E_gammass_case_4 = E_case_4.diff(gamma_ss)
    d1E_gammacs_case_4 = E_case_4.diff(gamma_cs)
    d1E_tauc_case_4 = E_case_4.diff(tau_c)
    d2E_rhoc2_case_4 = d1E_rhoc_case_4.diff(rho_c)
    d2E_rhocrhos_case_4 = d1E_rhoc_case_4.diff(rho_s)
    d2E_rhos2_case_4 = d1E_rhos_case_4.diff(rho_s)
    d2E_rhocgammacc_case_4 = d1E_rhoc_case_4.diff(gamma_cc)
    d2E_rhosgammacc_case_4 = d1E_rhos_case_4.diff(gamma_cc)
    d2E_gammacc2_case_4 = d1E_gammacc_case_4.diff(gamma_cc)
    d2E_rhocgammass_case_4 = d1E_rhoc_case_4.diff(gamma_ss)
    d2E_rhosgammass_case_4 = d1E_rhos_case_4.diff(gamma_ss)
    d2E_gammaccgammass_case_4 = d1E_gammacc_case_4.diff(gamma_ss)
    d2E_gammass2_case_4 = d1E_gammass_case_4.diff(gamma_ss)
    d2E_rhocgammacs_case_4 = d1E_rhoc_case_4.diff(gamma_cs)
    d2E_rhosgammacs_case_4 = d1E_rhos_case_4.diff(gamma_cs)
    d2E_gammaccgammacs_case_4 = d1E_gammacc_case_4.diff(gamma_cs)
    d2E_gammassgammacs_case_4 = d1E_gammass_case_4.diff(gamma_cs)
    d2E_gammacs2_case_4 = d1E_gammacs_case_4.diff(gamma_cs)
    d2E_rhoctauc_case_4 = d1E_rhoc_case_4.diff(tau_c)
    d2E_rhostauc_case_4 = d1E_rhos_case_4.diff(tau_c)
    d2E_gammacctauc_case_4 = d1E_gammacc_case_4.diff(tau_c)
    d2E_gammasstauc_case_4 = d1E_gammass_case_4.diff(tau_c)
    d2E_gammacstauc_case_4 = d1E_gammacs_case_4.diff(tau_c)
    d2E_tauc2_case_4 = d1E_tauc_case_4.diff(tau_c)
    
    Kernel = [func.PBEc(parameters), spec_func.PBEc_alpha_replaced(parameters), spec_func.PBEc_beta_replaced(parameters),
            E_case_1, 
            E_case_2, 
            E_case_3, 
            E_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,
                1,
                1,
                1]
    diff_idx = [[0],
                [0],
                [0],
                [0]]
    dalprint_spec.dalton_functional_printer(Kernel, "ESRC_SPIN_TPSS_GWS_ERF", ["rho_c","rho_s","gamma_cc","gamma_ss","gamma_cs","tau_c"],
                                    ["E",
                                        "E",
                                        "E",
                                        "E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [func.PBEc(parameters), spec_func.PBEc_alpha_replaced(parameters), spec_func.PBEc_beta_replaced(parameters),
            E_case_1, d1E_rhoc_case_1, d1E_rhos_case_1, d1E_gammacc_case_1, d1E_gammass_case_1, d1E_gammacs_case_1, d1E_tauc_case_1,
            E_case_2, d1E_rhoc_case_2, d1E_rhos_case_2, d1E_gammacc_case_2, d1E_gammass_case_2, d1E_gammacs_case_2, d1E_tauc_case_2, 
            E_case_3, d1E_rhoc_case_3, d1E_rhos_case_3, d1E_gammacc_case_3, d1E_gammass_case_3, d1E_gammacs_case_3, d1E_tauc_case_3, 
            E_case_4, d1E_rhoc_case_4, d1E_rhos_case_4, d1E_gammacc_case_4, d1E_gammass_case_4, d1E_gammacs_case_4, d1E_tauc_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,6,
                1,6,
                1,6,
                1,6]
    diff_idx = [[0],[1,2,3,4,5,6],
                [0],[1,2,3,4,5,6],
                [0],[1,2,3,4,5,6],
                [0],[1,2,3,4,5,6]]
    dalprint_spec.dalton_functional_printer(Kernel, "D1ESRC_SPIN_TPSS_GWS_ERF", ["rho_c","rho_s","gamma_cc","gamma_ss","gamma_cs","tau_c"],
                                    ["E","d1E",
                                        "E","d1E",
                                        "E","d1E",
                                        "E","d1E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [func.PBEc(parameters), spec_func.PBEc_alpha_replaced(parameters), spec_func.PBEc_beta_replaced(parameters),
            E_case_1, d1E_rhoc_case_1, d1E_rhos_case_1, d1E_gammacc_case_1, d1E_gammass_case_1, d1E_gammacs_case_1, d1E_tauc_case_1,d2E_rhoc2_case_1,d2E_rhocrhos_case_1,d2E_rhos2_case_1,d2E_rhocgammacc_case_1,d2E_rhosgammacc_case_1,d2E_gammacc2_case_1,d2E_rhocgammass_case_1,d2E_rhosgammass_case_1,d2E_gammaccgammass_case_1,d2E_gammass2_case_1,d2E_rhocgammacs_case_1,d2E_rhosgammacs_case_1,d2E_gammaccgammacs_case_1,d2E_gammassgammacs_case_1,d2E_gammacs2_case_1,d2E_rhoctauc_case_1,d2E_rhostauc_case_1,d2E_gammacctauc_case_1,d2E_gammasstauc_case_1,d2E_gammacstauc_case_1,d2E_tauc2_case_1,
            E_case_2, d1E_rhoc_case_2, d1E_rhos_case_2, d1E_gammacc_case_2, d1E_gammass_case_2, d1E_gammacs_case_2, d1E_tauc_case_2,d2E_rhoc2_case_2,d2E_rhocrhos_case_2,d2E_rhos2_case_2,d2E_rhocgammacc_case_2,d2E_rhosgammacc_case_2,d2E_gammacc2_case_2,d2E_rhocgammass_case_2,d2E_rhosgammass_case_2,d2E_gammaccgammass_case_2,d2E_gammass2_case_2,d2E_rhocgammacs_case_2,d2E_rhosgammacs_case_2,d2E_gammaccgammacs_case_2,d2E_gammassgammacs_case_2,d2E_gammacs2_case_2,d2E_rhoctauc_case_2,d2E_rhostauc_case_2,d2E_gammacctauc_case_2,d2E_gammasstauc_case_2,d2E_gammacstauc_case_2,d2E_tauc2_case_2, 
            E_case_3, d1E_rhoc_case_3, d1E_rhos_case_3, d1E_gammacc_case_3, d1E_gammass_case_3, d1E_gammacs_case_3, d1E_tauc_case_3,d2E_rhoc2_case_3,d2E_rhocrhos_case_3,d2E_rhos2_case_3,d2E_rhocgammacc_case_3,d2E_rhosgammacc_case_3,d2E_gammacc2_case_3,d2E_rhocgammass_case_3,d2E_rhosgammass_case_3,d2E_gammaccgammass_case_3,d2E_gammass2_case_3,d2E_rhocgammacs_case_3,d2E_rhosgammacs_case_3,d2E_gammaccgammacs_case_3,d2E_gammassgammacs_case_3,d2E_gammacs2_case_3,d2E_rhoctauc_case_3,d2E_rhostauc_case_3,d2E_gammacctauc_case_3,d2E_gammasstauc_case_3,d2E_gammacstauc_case_3,d2E_tauc2_case_3, 
            E_case_4, d1E_rhoc_case_4, d1E_rhos_case_4, d1E_gammacc_case_4, d1E_gammass_case_4, d1E_gammacs_case_4, d1E_tauc_case_4,d2E_rhoc2_case_4,d2E_rhocrhos_case_4,d2E_rhos2_case_4,d2E_rhocgammacc_case_4,d2E_rhosgammacc_case_4,d2E_gammacc2_case_4,d2E_rhocgammass_case_4,d2E_rhosgammass_case_4,d2E_gammaccgammass_case_4,d2E_gammass2_case_4,d2E_rhocgammacs_case_4,d2E_rhosgammacs_case_4,d2E_gammaccgammacs_case_4,d2E_gammassgammacs_case_4,d2E_gammacs2_case_4,d2E_rhoctauc_case_4,d2E_rhostauc_case_4,d2E_gammacctauc_case_4,d2E_gammasstauc_case_4,d2E_gammacstauc_case_4,d2E_tauc2_case_4]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,6,21,
                1,6,21,
                1,6,21,
                1,6,21]
    diff_idx = [[0],[1,2,3,4,5,6],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21],
                [0],[1,2,3,4,5,6],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21],
                [0],[1,2,3,4,5,6],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21],
                [0],[1,2,3,4,5,6],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]]
    dalprint_spec.dalton_functional_printer(Kernel, "D2ESRC_SPIN_TPSS_GWS_ERF", ["rho_c","rho_s","gamma_cc","gamma_ss","gamma_cs","tau_c"],
                                    ["E","d1E","d2E",
                                        "E","d1E","d2E",
                                        "E","d1E","d2E",
                                        "E","d1E","d2E"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()

    
 
if "PBE_nomu_correlation" in write_list:
    out_file = open("../srfunctionals/PBE_nomu_correlation.F","w+")
    # ######################################################################
    # PBE, mu=0.0, correlation
    # ######################################################################
    E = func.PBEc(parameters).subs({rho_s: 0})*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_gammacc = E.diff(gamma_cc)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
    d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)

    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_PBE", ["rho_c","gamma_cc"],["E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # # #
    Kernel = [E,d1E_rhoc,d1E_gammacc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_PBE", ["rho_c","gamma_cc"],["E","d1E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])

    # # #
    Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,3],[1,4,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # PBE correlatoin, spin
    # ######################################################################
    E = func.PBEc(parameters)*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_rhos = E.diff(rho_s)
    d1E_gammacc = E.diff(gamma_cc)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocrhos = d1E_rhoc.diff(rho_s)
    d2E_rhos2 = d1E_rhos.diff(rho_s)
    d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
    d2E_rhosgammacc = d1E_rhos.diff(gamma_cc)
    d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_SPIN_PBE", ["rho_c","rho_s","gamma_cc"],["E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d1E_gammacc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3]
    diff_idx = [[0],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_SPIN_PBE", ["rho_c","rho_s","gamma_cc"],["E","d1E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d1E_gammacc,d2E_rhoc2,d2E_rhocrhos,d2E_rhos2,d2E_rhocgammacc,d2E_rhosgammacc,d2E_gammacc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,3,6]
    diff_idx = [[0],[1,2,3],[1,2,3,4,5,6]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_SPIN_PBE", ["rho_c","rho_s","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=False,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()
    
    

if "wPBE_exchange" in write_list:
    out_file = open("../srfunctionals/wPBE_exchange.F","w+")
    # ######################################################################
    # PBE exchange
    # ######################################################################
    # rho_c = 2*rho_a
    E = mu_func.wPBEx(parameters)*2*rho_a
    d1E_rhoa = E.diff(rho_a)
    d1E_gammaaa = E.diff(gamma_aa)
    d2E_rhoa2 = d1E_rhoa.diff(rho_a)
    d2E_rhoagammaaa = d1E_rhoa.diff(gamma_aa)
    d2E_gammaaa2 = d1E_gammaaa.diff(gamma_aa)
    
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRX_wPBE", ["rho_a","gamma_aa"],["Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E, d1E_rhoa, d1E_gammaaa]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRX_wPBE", ["rho_a","gamma_aa"],["Ea","d1Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E, d1E_rhoa, d1E_gammaaa, d2E_rhoa2, d2E_rhoagammaaa, d2E_gammaaa2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRX_wPBE", ["rho_a","gamma_aa"],["Ea","d1Ea","d2Ea"],
                                    description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()
    
    

if "VWN5_nomu_correlation" in write_list:
    out_file = open("../srfunctionals/VWN5_nomu_correlation.F","w+")
    # ######################################################################
    # VWN5 correlation, no-spin
    # ######################################################################
    E = func.VWN5c(parameters).subs({rho_s: 0})*rho_c
    d1E_rhoc = E.diff(rho_c)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_VWN5", ["rho_c"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1]
    diff_idx = [[0],[1]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_VWN5", ["rho_c"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d2E_rhoc2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,1,1]
    diff_idx = [[0],[1],[1]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_VWN5", ["rho_c"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    
    # ######################################################################
    # VWN5 correlation, spin
    # ######################################################################
    E = func.VWN5c(parameters)*rho_c
    d1E_rhoc = E.diff(rho_c)
    d1E_rhos = E.diff(rho_s)
    d2E_rhoc2 = d1E_rhoc.diff(rho_c)
    d2E_rhocrhos = d1E_rhoc.diff(rho_s)
    d2E_rhos2 = d1E_rhos.diff(rho_s)
    
    Kernel = [E]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1]
    diff_idx = [[0]]
    dalprint.dalton_functional_printer(Kernel, "ESRC_SPIN_VWN5", ["rho_c","rho_s"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2]
    diff_idx = [[0],[1,2]]
    dalprint.dalton_functional_printer(Kernel, "D1ESRC_SPIN_VWN5", ["rho_c","rho_s"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    # # #
    Kernel = [E,d1E_rhoc,d1E_rhos,d2E_rhoc2,d2E_rhocrhos,d2E_rhos2]
    description = ["Implemented by E.R. Kjellgren.\n"]
    diff_order = [1,2,3]
    diff_idx = [[0],[1,2],[1,2,3]]
    dalprint.dalton_functional_printer(Kernel, "D2ESRC_SPIN_VWN5", ["rho_c","rho_s"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[out_file])
    out_file.close()
    
# ######################################################################
# ######################################################################
