
#Author:   Patrick Ettenhuber
#Purpose:  Check whether all arguments defined in a structure are communicated with the associated lsmpi_buffer routine, this script can easily be generalized even further

import os

def main():

   all_good = True
   good     = True
   compare_struct_and_its_bcast(good,"../lsutil/dec_typedef.F90","decsettings","../deccc/decmpi.F90","mpicopy_dec_settings",exceptions=["cc_models"])
   all_good = all_good and good
   #compare_struct_and_its_bcast(good,"../lsutil/TYPE-DEF.F90","LSSETTING","../lsutil/lsmpi-operations.F90","mpicopy_setting",exceptions=["MOLECULE","BASIS","FRAGMENT","SCHEME","IO","lst_dLHS","lst_dRHS","DmatLHS","DmatRHS","LST_GAB_LHS","LST_GAB_RHS","FRAGMENTS","OUTPUT","GGem","RedCS"])
   all_good = all_good and good

   if all_good:
      print "TESTSTATUS: GOOD"
   else:
      print "TESTSTATUS: FAIL"

#ALL OK:               returns the status of the test, True means all passed, False, something is wrong
#filename_with_struct: input the path to the file where the struct is defined relative to the path where the script is stored
#struct_name:          name of the type to be bcasted
#filename_with_bcast:  input the path to the file where the bcasting of the struct is defined relative to the path where the script is stored
#bcast_routine_name:   name of the subroutine that is used to bcast the struct
def compare_struct_and_its_bcast(ALL_OK,filename_with_struct,struct_name,filename_with_bcast,bcast_routine_name,exceptions=[]):

   swd = os.path.realpath(__file__)
   swd = swd[:swd.rfind("/")+1]

   print "checking: "+struct_name

   filename_with_struct = swd+filename_with_struct
   filename_with_bcast  = swd+filename_with_bcast

   dec_typedef = open(filename_with_struct,'r')
   found_in_bcast_list = []
   variable_name_list  = []
   
   dec_typedef_lines = dec_typedef.readlines()
   dec_typedef.close()
   
   #ADD SPECIAL VARIABLES TO IGNORE TO THIS LIST
   ignore_variables = exceptions
   
   for line_nr in range(len(dec_typedef_lines)):
      if "type "+struct_name.lower() in dec_typedef_lines[line_nr].lower():
         while not "end type "+struct_name.lower() in dec_typedef_lines[line_nr].lower():
            line_nr += 1
   
            if "end type "+struct_name.lower() in dec_typedef_lines[line_nr].lower() :
               break
   
            #ignore comment lines and empty lines
            line = dec_typedef_lines[line_nr].lower().strip()
            if( len(line) > 1 ):
               if( line[0] != "!" ):
                  variable_line = line.split(':',1)[-1].replace(':','').strip()
                  #truncate if there is an inline comment
                  if( "!" in variable_line):
                     variable_line = variable_line[0:variable_line.find('!')]
   
                  variables = variable_line.split(',')
   
                  for variable in variables:
                     #remove parenthesis from variable names, parenthesis are usually at the end of the variable names and begin with "("
                     if( "(" in variable):
                        variable = variable[0:variable.find('(')]
   
                     variable = variable.strip()
                     add_to_list = True
                     for i in range(len(ignore_variables)):
                        if variable == ignore_variables[i].lower() :
                           add_to_list = False
   
                     if add_to_list:
                        variable_name_list.append(variable)
                        found_in_bcast_list.append(False)
   
         break
   
   if( len(variable_name_list) != len(found_in_bcast_list) ):
      print "SOMETHING WRONG FINDING THE VARIABLES"
   
   decmpi_file = open(filename_with_bcast,'r')
   decmpi_lines = decmpi_file.readlines()
   decmpi_file.close()
   
   for line_nr in range(len(decmpi_lines)):
      if "subroutine "+bcast_routine_name.lower() in decmpi_lines[line_nr].lower().strip():
         while not "end subroutine "+bcast_routine_name.lower() in decmpi_lines[line_nr].lower().strip():
            line_nr += 1
            line = decmpi_lines[line_nr].lower().strip()
            if "call ls_mpi_buffer" in line and "%" in line:
               print line
               variable = line[line.find("(")+1:line.find(")")].split(",",1)[0].split("%")[1].strip()
               for var_nr in range(len(variable_name_list)):
                  if variable_name_list[var_nr] == variable:
                     found_in_bcast_list[var_nr] = True
         break
   
   ALL_FOUND = True
   for var_nr in range(len(found_in_bcast_list)):
      if not found_in_bcast_list[var_nr] :
         print "VARIABLE:",variable_name_list[var_nr]," not found in bcast"
         ALL_FOUND = False
   
   if ALL_FOUND:
      print "NO PROBLEMS DETECTED IN "+struct_name
   else:
      print "PROBLEMS DETECTED IN "+struct_name

   ALL_OK = ALL_FOUND



main()
