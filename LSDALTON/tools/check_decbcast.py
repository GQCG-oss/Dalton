import os

dec_typedef = open("../lsutil/dec_typedef.F90",'r')
found_in_bcast_list = []
variable_name_list  = []

dec_typedef_lines = dec_typedef.readlines()
dec_typedef.close()

for line_nr in range(len(dec_typedef_lines)):
   if "type decsettings" in dec_typedef_lines[line_nr].lower():
      while not "end type decsettings" in dec_typedef_lines[line_nr].lower():
         line_nr += 1

         if "end type decsettings" in dec_typedef_lines[line_nr].lower() :
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

                  variable_name_list.append(variable.strip())
                  found_in_bcast_list.append(False)

      break

decmpi_file = open("../deccc/decmpi.F90",'r')
decmpi_lines = decmpi_file.readlines()
decmpi_file.close()

for line_nr in range(len(decmpi_lines)):
   if "subroutine mpicopy_dec_settings(decitem)" in decmpi_lines[line_nr].lower().strip():
      while not "end subroutine mpicopy_dec_settings" in decmpi_lines[line_nr].lower().strip():
         line_nr += 1
         line = decmpi_lines[line_nr].lower().strip()
         if "call ls_mpi_buffer" in line:
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
   print "NO PROBLEMS DETECTED"
