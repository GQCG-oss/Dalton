from sympy import cse, atan, ln, fcode, sqrt, exp, erf
import sympy
import datetime

# TPSSc need a special case, because of the Max() inside the functional.
# This problem is solved by splitting TPSSc up in four cases, but to minimize
# terms that is needed to be calculated, since the terms inside the Max() are quite heavy,
# all four cases will be in the same subroutine, thus the "normal" print_functional_to_DALTON
# script cannot handle this.

def print_function(string_list, output_file):
    output_file.write("\n")
    output_file.write("\n")
    for line in string_list:
        if line[0] == "C":
            next_string = ""
            for i in line:
                next_string += i
            output_file.write(next_string)
        else:
            next_string = "      "
            for i in line:
                next_string += i
                if len(next_string) == 72:
                    if next_string[-1] != "\n":
                        # To prevent case were the last bit had a lenght of 72
                        #  if this is the case, priint it in the print statement after the loop
                        output_file.write(next_string+"\n")
                        next_string = "     &"
            output_file.write(next_string)
            if len(next_string) > 72:
                print("LINE TO LONG")
                break
    output_file.write("\n")

def dalton_functional_printer(kernel, routinename, input_variables, output_variables, shortrange=False, description=["No description.\n"], diff_order=[1], diff_idx=[[0]], output_files=[]):
    """
    shortrange = True, if mu is an input parameter
    diff_order, is a list that indicatices how many of of each derivative that is present
       diff_order[0] = Energy, diff_order[1] = number of first derivatives
    diff_idx, is a list of lists that containts the idx of the given derivatives.
    """
    if len(diff_order) != len (diff_idx):
        print("ERROR: order of derivatives does not match number of idx lists")
    for i in range(0, len(diff_order)):
        if diff_order[i] != len(diff_idx[i]):
            print("ERROR: Number of "+i+"'th derivatives does not match lenght of corrosponding idx list")
    # Format kernel expression
    #kernel = kernel.evalf() # evalf is super slow!
    E_cse = cse(kernel, optimizations="basic")
    E_clean = []
    E_clean_temp = []
    for i in range(len(E_cse[0])):
        E_clean.append((E_cse[0][i][0],E_cse[0][i][1].evalf()))
    for i in range(len(E_cse[1])):
        E_clean_temp.append(E_cse[1][i].evalf())
    E_clean = (E_clean, E_clean_temp)
    string_list = []
    subroutine_input = ", ".join(input_variables)
    if shortrange == True:
        subroutine_input += ", mu"
    subroutine_input += ", "+", ".join(output_variables[0:int(len(output_variables)/4)])
    string_list.append("C*****************************************************************************\n")
    string_list.append("pure subroutine "+routinename+"("+subroutine_input+")\n")
    string_list.append("C*****************************************************************************\n")
    for i in description:
        string_list.append("C   "+i)
    string_list.append("C\nC   Subroutine generated using Sympy "+str(sympy.__version__)+"\n")
    string_list.append("C   Generated: "+datetime.datetime.now().strftime("%B %d, %Y")+"\n")
    string_list.append("C*****************************************************************************\n")
    string_list.append("implicit none\n")
    string_list.append("real*8, intent(in) :: "+", ".join(input_variables))
    if shortrange == True:
        string_list[-1] += ", mu"
    string_list[-1] += "\n"
    string_list.append("real*8, intent(out) :: "+", ".join(output_variables[0:int(len(output_variables)/4)])+"\n")
    if "d1Ea" in output_variables:
        string_list[-1] = string_list[-1].replace("d1Ea","d1Ea(4)")
    else:
        string_list[-1] = string_list[-1].replace("d1E","d1E(9)")
    if "d2Ea" in output_variables:
        string_list[-1] = string_list[-1].replace("d2Ea","d2Ea(10)")
    else:
        string_list[-1] = string_list[-1].replace("d2E","d2E(45)")
    string_list.append("real*8 :: PBE_value, PBE_alpha_value, PBE_beta_value\n")
    if len(E_clean[0]) > 0:
        string_list.append("real*8 :: ")
        for term in E_clean[0]:
            for letter in str(fcode(term[0], source_format="free").replace('&\n      ',"")):
                string_list[-1] += letter
            if term != E_clean[0][-1]:
                string_list[-1] += ", "
        string_list[-1] += "\n"
    if len(diff_order)/4 >= 1 and "Ea" in output_variables:
        string_list.append("Ea = 0.0d0\n")
    elif len(diff_order)/4 >= 1:
        string_list.append("E = 0.0d0\n")
    if len(diff_order)/4 >= 2 and "d1Ea" in output_variables:
        string_list.append("d1Ea(:) = 0.0d0\n")
    elif len(diff_order)/4 >= 2:
        string_list.append("d1E(:) = 0.0d0\n")
    if len(diff_order)/4 >= 3 and "d2Ea" in output_variables:
        string_list.append("d2Ea(:) = 0.0d0\n")
    elif len(diff_order)/4 >= 3:
        string_list.append("d2E(:) = 0.0d0\n")
            
    for term in E_clean[0]:
        string_list.append(str(term[0])+" = ")
        for letter in str(fcode(term[1], source_format="free").replace('&\n      ',"").replace("parameter (pi = 3.1415926535897932d0)\n", "")):
            string_list[-1] += letter
        string_list[-1] += "\n"
    
    output_counter = 0
    diff_counter = 0
    term_counter = 0
    for term in E_clean[1]:
        if term_counter < 3:
            if term_counter == 0:
                string_list.append("PBE_value = ")
            elif term_counter == 2:
                string_list.append("PBE_alpha_value = ")
            elif term_counter == 1:
                string_list.append("PBE_beta_value = ")
            for letter in str(fcode(term, source_format="free").replace('&\n      ',"").replace("parameter (pi = 3.1415926535897932d0)\n", "")):
                string_list[-1] += letter
            string_list[-1] += "\n"
            term_counter += 1
        else:
            # divide by four because of 4 different cases
            if len(diff_order)/4 - len(diff_order)/4 == diff_counter:
                string_list.append("if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .ge. PBE_value)) then\n")
            elif len(diff_order)/4 * 2 - len(diff_order)/4 == diff_counter:
                string_list.append("end if\n")
                string_list.append("if ((PBE_alpha_value .ge. PBE_value) .and. (PBE_beta_value .lt. PBE_value)) then\n")
            elif len(diff_order)/4 * 3 - len(diff_order)/4 == diff_counter:
                string_list.append("end if\n")
                string_list.append("if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .ge. PBE_value)) then\n")
            elif len(diff_order)/4 * 4 - len(diff_order)/4 == diff_counter:
                string_list.append("end if\n")
                string_list.append("if ((PBE_alpha_value .lt. PBE_value) .and. (PBE_beta_value .lt. PBE_value)) then\n")
            
            if diff_idx[diff_counter][output_counter] == 0:
                string_list.append(output_variables[0]+" = ")
            else:
                string_list.append(output_variables[diff_counter]+"("+str(diff_idx[diff_counter][output_counter])+") = ")
            for letter in str(fcode(term, source_format="free").replace('&\n      ',"").replace("parameter (pi = 3.1415926535897932d0)\n", "")):
                string_list[-1] += letter
            string_list[-1] += "\n"
            output_counter += 1            
            if diff_order[diff_counter] == output_counter:
                output_counter = 0
                diff_counter += 1

    string_list.append("end if\n")
    
    string_list.append("end subroutine")
    for output_file in output_files:
        print_function(string_list, output_file)