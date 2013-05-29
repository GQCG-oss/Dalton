#!/usr/bin/env python2
#----------------------------------
# generate dft(gga) code for dirac
# Andreas Hesselmann <andreas@theochem.uni-duesseldorf.de>, 2003
# Adapted to new density property input and unrestricted functionals
# by Pawel Salek.
# Generalized for 4th-order derivatives by Olav Vahtras
#----------------------------------
# needs: maxima 
# usage: codegen.py 'file.max'
# 'file.max' is a  maxima input file of the following general form:
#      ....some variables....
#      ....some functions...
#      K(rhoa,ghroa,rhob,gradb):= 'expression of rhoa,grhoa,xa
#                                                rhob,grhob,xb'
# where xa=sqrt(grada*grada)/rhoa^(4/3)
# output: fun-'file'.c

import sys, string, os, re;

if __name__ == "__main__":
    funcfil=sys.argv[1];
    funcname=string.split(funcfil,'.')[0]

    #header
    try:
        f=open(funcfil,'r')
        func=f.readlines()
        f.close()
    except IOError:
        print "cannot open "+funcfil+" for reading. Program stops."
        sys.exit(1)
        
    f=open("fun-"+funcname+".c",'w')
    f.write("/* Automatically generated functional code: "+funcname+"\n")
    f.write("   Maxima input:\n")
    for line in func:
        line=string.replace(line,"/*","//")
        line=string.replace(line,"*/","")
        f.write(4*" "+">> "+line)
    f.write("*/\n\n")
    f.write("// add \"extern Functional "+funcname+"Functional;\" to 'functionals.h'\n")
    f.write("// add \"&"+funcname+"Functional,\" to 'functionals.c'\n")
    f.write("// add \"fun-"+funcname+".c\" to 'Makefile.in'\n\n")
    f.write("#include <math.h>\n")
    f.write("#include <stddef.h>\n\n")
    f.write("#define __CVERSION__\n\n")
    f.write("#include \"functionals.h\"\n")
    f.write("#define LOG log\n")
    f.write("#define ABS fabs\n")
    f.write("#define ASINH asinh\n")
    f.write("#define SQRT sqrt\n\n")
    f.write("/* INTERFACE PART */\n")
    f.write("static int "+funcname+"_isgga(void) {return 1;}\n")
    f.write("static int "+funcname+"_read(const char* conf_line);\n")
    f.write("static real "+funcname+"_energy(const FunDensProp* dp);\n")
    f.write("static void "+funcname+"_first(FunFirstFuncDrv *ds, real factor, \n\
                           const FunDensProp* dp);\n")
    f.write("static void "+funcname+"_second(FunSecondFuncDrv *ds, real factor,\n\
                            const FunDensProp* dp);\n")
    f.write("static void "+funcname+"_third(FunThirdFuncDrv *ds, real factor,\n\
                           const FunDensProp* dp);\n\n")
    f.write("static void "+funcname+"_fourth(FunFourthFuncDrv *ds, real factor,\n\
                           const FunDensProp* dp);\n\n")
    f.write("//static int fun_true(void) { return 1; }\n")
    f.write("Functional "+funcname+"Functional = {\n")
    f.write("  \""+funcname+"\",\n")
    f.write("  "+funcname+"_isgga,\n")
    f.write("  3,\n")
    f.write("  "+funcname+"_read,\n")
    f.write("  NULL,\n")
    f.write("  "+funcname+"_energy,\n")
    f.write("  "+funcname+"_first,\n")
    f.write("  "+funcname+"_second,\n")
    f.write("  "+funcname+"_third,\n")
    f.write("  "+funcname+"_fourth\n};\n\n")
    f.write("/* IMPLEMENTATION PART */\n")
    f.write("static int\n"+funcname+"_read(const char* conf_line)\n")
    f.write("{\n    fun_set_hf_weight(0);\n")
    f.write("    return 1;\n}\n\n")
    f.close()

    #maxima commands
    com1="string(float(subst(pow,\"^\",optimize("
    com2="))));\n"
    maxima_commands={
        "en": "zk=K(rhoa,grada,rhob,gradb,gradab)",
        "first": "[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1),\n"
                +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1),\n"
                +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1),\n"
                +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1),\n"
                +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,1)]\n",
        "second":"[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1),\n"
                +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1),\n"
                +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1),\n"
                +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1),\n"
                +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,1),\n"
                +" d2fdrara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2),\n"
                +" d2fdrarb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1),\n"
                +" d2fdraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1),\n"
                +" d2fdragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1),\n"
                +" d2fdraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,1),\n"
                +" d2fdrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2),\n"
                +" d2fdrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1),\n"
                +" d2fdrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1),\n"
                +" d2fdrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,1),\n"
                +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2),\n"
                +" d2fdgagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1),\n"
                +" d2fdgaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,1),\n"
                +" d2fdgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2),\n"
                +" d2fdgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,1),\n"
                +" d2fdabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,2)]\n",
        "third": "[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1),\n"
                +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1),\n"
                +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1),\n"
                +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1),\n"
                +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,1),\n"
                +" d2fdrara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2),\n"
                +" d2fdrarb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1),\n"
                +" d2fdraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1),\n"
                +" d2fdragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1),\n"
                +" d2fdraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,1),\n"
                +" d2fdrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2),\n"
                +" d2fdrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1),\n"
                +" d2fdrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1),\n"
                +" d2fdrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,1),\n"
                +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2),\n"
                +" d2fdgagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1),\n"
                +" d2fdgaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,1),\n"
                +" d2fdgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2),\n"
                +" d2fdgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,1),\n"
                +" d2fdabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,2),\n"
                +" d3fdrarara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3),\n"
                +" d3fdrararb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1),\n"
                +" d3fdraraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,1),\n"
                +" d3fdraragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradb,1),\n"
                +" d3fdraraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradab,1),\n"
                +" d3fdrarbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2),\n"
                +" d3fdrarbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,1),\n"
                +" d3fdrarbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradb,1),\n"
                +" d3fdrarbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradab,1),\n"
                +" d3fdragaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,2),\n"
                +" d3fdragagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradb,1),\n"
                +" d3fdragaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradab,1),\n"
                +" d3fdragbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,2),\n"
                +" d3fdragbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1,gradab,1),\n"
                +" d3fdraabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,2),\n"
                +" d3fdrbrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3),\n"
                +" d3fdrbrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,1),\n"
                +" d3fdrbrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradb,1),\n"
                +" d3fdrbrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradab,1),\n"
                +" d3fdrbgaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2),\n"
                +" d3fdrbgagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradb,1),\n"
                +" d3fdrbgaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradab,1),\n"
                +" d3fdrbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,2),\n"
                +" d3fdrbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1,gradab,1),\n"
                +" d3fdrbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,2),\n"
                +" d3fdgagaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,3),\n"
                +" d3fdgagagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradb,1),\n"
                +" d3fdgagaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradab,1),\n"
                +" d3fdgagbgb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,2),\n"
                +" d3fdgagbab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1,gradab,1),\n"
                +" d3fdgaabab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,2),\n"
                +" d3fdgbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,3),\n"
                +" d3fdgbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2,gradab,1),\n"
                +" d3fdgbabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,2),\n"
                +" d3fdababab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,3)]\n",
        "fourth":"[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1),\n"
                +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1),\n"
                +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1),\n"
                +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1),\n"
                +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,1),\n"
                +" d2fdrara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2),\n"
                +" d2fdrarb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1),\n"
                +" d2fdraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1),\n"
                +" d2fdragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1),\n"
                +" d2fdraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,1),\n"
                +" d2fdrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2),\n"
                +" d2fdrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1),\n"
                +" d2fdrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1),\n"
                +" d2fdrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,1),\n"
                +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2),\n"
                +" d2fdgagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1),\n"
                +" d2fdgaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,1),\n"
                +" d2fdgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2),\n"
                +" d2fdgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,1),\n"
                +" d2fdabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,2),\n"
                +" d3fdrarara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3),\n"
                +" d3fdrararb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1),\n"
                +" d3fdraraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,1),\n"
                +" d3fdraragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradb,1),\n"
                +" d3fdraraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradab,1),\n"
                +" d3fdrarbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2),\n"
                +" d3fdrarbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,1),\n"
                +" d3fdrarbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradb,1),\n"
                +" d3fdrarbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradab,1),\n"
                +" d3fdragaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,2),\n"
                +" d3fdragagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradb,1),\n"
                +" d3fdragaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradab,1),\n"
                +" d3fdragbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,2),\n"
                +" d3fdragbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1,gradab,1),\n"
                +" d3fdraabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,2),\n"
                +" d3fdrbrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3),\n"
                +" d3fdrbrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,1),\n"
                +" d3fdrbrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradb,1),\n"
                +" d3fdrbrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradab,1),\n"
                +" d3fdrbgaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2),\n"
                +" d3fdrbgagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradb,1),\n"
                +" d3fdrbgaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradab,1),\n"
                +" d3fdrbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,2),\n"
                +" d3fdrbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1,gradab,1),\n"
                +" d3fdrbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,2),\n"
                +" d3fdgagaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,3),\n"
                +" d3fdgagagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradb,1),\n"
                +" d3fdgagaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradab,1),\n"
                +" d3fdgagbgb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,2),\n"
                +" d3fdgagbab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1,gradab,1),\n"
                +" d3fdgaabab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,2),\n"
                +" d3fdgbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,3),\n"
                +" d3fdgbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2,gradab,1),\n"
                +" d3fdgbabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,2),\n"
                +" d3fdababab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,3),\n"
                +" d4fdrararara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,4),\n"
                +" d4fdrarararb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3,rhob,1),\n"
                +" d4fdrararaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3,grada,1),\n"
                +" d4fdrararagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3,gradb,1),\n"
                +" d4fdrararaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3,gradab,1),\n"
                +" d4fdrararbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,2),\n"
                +" d4fdrararbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1,grada,1),\n"
                +" d4fdrararbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1,gradb,1),\n"
                +" d4fdrararbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1,gradab,1),\n"
                +" d4fdraragaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,2),\n"
                +" d4fdraragagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,1,gradb,1),\n"
                +" d4fdraragaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,1,gradab,1),\n"
                +" d4fdraragbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradb,2),\n"
                +" d4fdraragbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradb,1,gradab,1),\n"
                +" d4fdraraabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradab,2),\n"
                +" d4fdrarbrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,3),\n"
                +" d4fdrarbrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2,grada,1),\n"
                +" d4fdrarbrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2,gradb,1),\n"
                +" d4fdrarbrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2,gradab,1),\n"
                +" d4fdrarbgaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,2),\n"
                +" d4fdrarbgagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,1,gradb,1),\n"
                +" d4fdrarbgaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,1,gradab,1),\n"
                +" d4fdrarbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradb,2),\n"
                +" d4fdrarbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradb,1,gradab,1),\n"
                +" d4fdrarbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradab,2),\n"
                +" d4fdragagaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,3),\n"
                +" d4fdragagagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,2,gradb,1),\n"
                +" d4fdragagaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,2,gradab,1),\n"
                +" d4fdragagbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradb,2),\n"
                +" d4fdragagbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradb,1,gradab,1),\n"
                +" d4fdragaabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1,gradab,2),\n"
                +" d4fdragbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,3),\n"
                +" d4fdragbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,2,gradab,1),\n"
                +" d4fdragbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1,gradab,2),\n"
                +" d4fdraababab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,3),\n"
                +" d4fdrbrbrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,4),\n"
                +" d4fdrbrbrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3,grada,1),\n"
                +" d4fdrbrbrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3,gradb,1),\n"
                +" d4fdrbrbrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3,gradab,1),\n"
                +" d4fdrbrbgaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,2),\n"
                +" d4fdrbrbgagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,1,gradb,1),\n"
                +" d4fdrbrbgaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,1,gradab,1),\n"
                +" d4fdrbrbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradb,2),\n"
                +" d4fdrbrbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradb,1,gradab,1),\n"
                +" d4fdrbrbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradab,2),\n"
                +" d4fdrbgagaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,3),\n"
                +" d4fdrbgagagb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2,gradb,1),\n"
                +" d4fdrbgagaab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2,gradab,1),\n"
                +" d4fdrbgagbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradb,2),\n"
                +" d4fdrbgagbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradb,1,gradab,1),\n"
                +" d4fdrbgaabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1,gradab,2),\n"
                +" d4fdrbgbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,3),\n"
                +" d4fdrbgbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,2,gradab,1),\n"
                +" d4fdrbgbabab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1,gradab,2),\n"
                +" d4fdrbababab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,3),\n"
                +" d4fdgagagaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,4),\n"
                +" d4fdgagagagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,3,gradb,1),\n"
                +" d4fdgagagaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,3,gradab,1),\n"
                +" d4fdgagagbgb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradb,2),\n"
                +" d4fdgagagbab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradb,1,gradab,1),\n"
                +" d4fdgagaabab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2,gradab,2),\n"
                +" d4fdgagbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,3),\n"
                +" d4fdgagbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,2,gradab,1),\n"
                +" d4fdgagbabab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1,gradab,2),\n"
                +" d4fdgaababab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,3),\n"
                +" d4fdgbgbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,4),\n"
                +" d4fdgbgbgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,3,gradab,1),\n"
                +" d4fdgbgbabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2,gradab,2),\n"
                +" d4fdgbababab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,3),\n"
                +" d4fdabababab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,4)]\n"
        }
                   

def make_maximafile(funcfil,com):
    try:
        f=open(funcfil,'r')
        func=f.readlines()
        f.close()
    except:
        print "Cannot open "+funcfil+"for reading."
        exit
        
    f=open('maxima.in','w')
    f.write("xa:sqrt(grada*grada)/rhoa^(4/3);\n");
    f.write("xb:sqrt(gradb*gradb)/rhob^(4/3);\n");
    for line in func:
        f.write(line)
    f.write(com1+maxima_commands[com]+com2)
    f.close()

def run_maxima(outfil):
    root, ext = os.path.splitext(outfil)
    if os.system("maxima < maxima.in >"+outfil)!=0:
        sys.exit("Error in run_maxima!")
    root, ext = os.path.splitext(outfil)
    infil = root + ".in"
    os.rename("maxima.in", infil)

def vars_one2field(s):
    f=0
    s=list(s)
    for i in range(len(s)):
        if s[i]=="%":
            s[i]="t["
            f=1
        elif f==1 and not s[i] in ["0","1","2","3","4","5","6","7","8","9"]:
            s[i]="]"+s[i]
            f=0
    if f==1:
        s.append("]")
    string=""
    for c in s:
        string=string+c
    return string

def f2c(expr):
    s=string.split(expr,",")
    nvars = 0
    for i in range(len(s)):
        if s[i][0]=="%" and s[i][-1]=="]":
            nvars=int(s[i][1:-1])
            s=s[i+1:]
            break
        if s[i][-4:]=="[%1]":
            nvars=1
            s=s[i+1:]
            break
    print "Temporary variables ", nvars
    s2=[]
    i=-1
    for s1 in s:
        if ":" in s1 or "=" in s1:
            s2.append(s1)
            i=i+1
        else:
            s2[i]=s2[i]+","+s1
    s=s2
    var={}
    for i in range(nvars):
        var["%"+str(i+1)]="t["+str(i+1)+"]"
    for i in range(len(s)):
        s[i]=vars_one2field(s[i]) #replace %... by t[...]
        if s[i][0:2]=="zk":
            s[i]=s[i][:-1]
        if s[i][0]=="[":
            s[i]=s[i][1:]
        s[i]=string.replace(s[i],":"," = ")+";"
    if s[-1][0:2] =="zk":
        s[-1]=s[-1][:-1]+";"
    else:
        s[-1]=s[-1][:-3]+";"
    return nvars,s
        
def read_maximaout(outfil):
    print outfil
    f=open(outfil,'r')
    expr=""
    block=0
    for line in f.readlines():
        data=string.split(line)
        if len(data)==0: continue
        if data[0]=="\n": continue
        if len(data)>=2:
            if len(data[1])>=5:
                if data[1][0:5]=="block":
                    block=1
                    expr=""
        if block==1:
            if data[0][0:3]!="(%i":
                line=string.strip(line)
                line=string.replace(line,"#","")
                line=string.replace(line,"\n","")
                line=string.replace(line,"\\","")
                expr=expr+line
            else:
                block=0
    nvars,cexpr=f2c(expr)
    return nvars,cexpr

def parse_maxima(outfil):
    lines = "".join([ l.strip().strip('\\') for l in open(outfil) ])
    has_block = re.match(".*(block\(.*\))\(%i\d+\)$", lines)
    has_zk = re.match(".* (\[?\w+ = .*)\(%i\d+\)$", lines)
    if has_block:
        expr =  has_block.groups()[0]
    elif has_zk:
        expr =  has_zk.groups()[0]
    else:
        expr = None
    return expr

def max2c(expr):
    """Transform maxima output to c-code

    Expected input: 
        alt 1: 'block([%1, %2, ...], %1:<def1>, %2:<def2>, ..., <code>)'
        alt 2: '<code>'
    """
    has_block = re.match(".*block\(\[(.*)\],(.*)\)$", expr)
    if has_block:
        # Separate the first argument to block from the rest
        block, rest  = has_block.groups()
        # first argument [%1, %2... ] gives the number of temporary variables
        nvars = len(block.split(','))
        # split remining into nvars definitions %1:<def> and actual source
        defs = re.split(",?%\d+:", rest)[1:]
        # not foolproof match - assumes no internal commas in definition
        lastdef, code = re.match('(.*?),(\s*\[?\w+\s*=.*)', defs[-1]).groups()
        defs[-1] = lastdef
        # Substitution of %n to t[n] in definitions and code string
        def p2t(match):
            pn = match.group()
            return 't[%s]' % pn[1:]
        p = re.compile(r'%\d+')
        c_defs = [ 't[%d] = %s;' % (i+1, p.sub(p2t, defs[i])) for i in range(nvars) ]
        c_code = p.sub(p2t, code) 
        #print "\nc_code\n",c_code
    else:
        nvars = 0
        c_defs = []
        c_code = expr

    #For derivatives split the codes 
    cmps = re.split('[,\[]?(\w+\s*=\s*)', c_code)[1:]
    # Remove potential final closing bracket
    cmps[-1] = re.sub('\]?$', '', cmps[-1])
    #print "\ncmps\n","\n", "\n".join(cmps)
    #Parwise concatenation
    for i in range(0,len(cmps),2):
        #print cmps[i] + cmps[i+1]
        c_defs.append(cmps[i] + cmps[i+1] + ";")

    return nvars, c_defs

def make_stdefs(f):
    f.write(4*" "+"real rhoa = dp->rhoa;\n"+
            4*" "+"real rhob = dp->rhob;\n"+
            4*" "+"real grada = dp->grada;\n"+
            4*" "+"real gradb = dp->gradb;\n"+
            4*" "+"real gradab = dp->gradab;\n\n")



def make_energycode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic real\n")
    f.write(funcname+"_energy(const FunDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"],zk;\n")
    make_stdefs(f)
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"return zk;\n")
    f.write("}\n")
    f.close()

def make_firstcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdrb, dfdga, dfdgb, dfdab;\n")
    make_stdefs(f)
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0100 += factor*dfdrb;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df0001 += factor*dfdgb;\n")
    f.write(4*" "+"ds->df00001 += factor*dfdab;\n")
    f.write("}\n")
    f.close()    

def make_secondcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdrb, dfdga, dfdgb, dfdab;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;\n")
    f.write(4*" "+"real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;\n")
    f.write(4*" "+"real d2fdraab, d2fdrbab;\n")
    f.write(4*" "+"real d2fdgaab, d2fdgbab, d2fdabab, d2fdgagb;\n")
    make_stdefs(f)
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0100 += factor*dfdrb;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df0001 += factor*dfdgb;\n")
    f.write(4*" "+"ds->df00001 += factor*dfdab;\n")
    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df1100 += factor*d2fdrarb;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df1001 += factor*d2fdragb;\n")
    f.write(4*" "+"ds->df10001 += factor*d2fdraab;\n")
    f.write(4*" "+"ds->df0200 += factor*d2fdrbrb;\n")
    f.write(4*" "+"ds->df0110 += factor*d2fdrbga;\n")
    f.write(4*" "+"ds->df0101 += factor*d2fdrbgb;\n")
    f.write(4*" "+"ds->df01001 += factor*d2fdrbab;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write(4*" "+"ds->df0011 += factor*d2fdgagb;\n")
    f.write(4*" "+"ds->df00101+= factor*d2fdgaab;\n")
    f.write(4*" "+"ds->df0002 += factor*d2fdgbgb;\n")
    f.write(4*" "+"ds->df00011+= factor*d2fdgbab;\n")
    f.write(4*" "+"ds->df00002+= factor*d2fdabab;\n")
    f.write("}\n")
    f.close()    

def make_thirdcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdrb, dfdga, dfdgb, dfdab;\n")

    f.write(4*" "+"real d2fdrara , d2fdrarb , d2fdraga , d2fdragb ;\n")
    f.write(4*" "+"real d2fdraab , d2fdrbrb , d2fdrbga , d2fdrbgb ;\n")
    f.write(4*" "+"real d2fdrbab , d2fdgaga , d2fdgagb , d2fdgaab ;\n")
    f.write(4*" "+"real d2fdgbgb , d2fdgbab , d2fdabab ;\n")

    f.write(4*" "+"real d3fdrarara , d3fdrararb , d3fdraraga , d3fdraragb ;\n")
    f.write(4*" "+"real d3fdraraab , d3fdrarbrb , d3fdrarbga , d3fdrarbgb ;\n")
    f.write(4*" "+"real d3fdrarbab , d3fdragaga , d3fdragagb , d3fdragaab ;\n")
    f.write(4*" "+"real d3fdragbgb , d3fdragbab , d3fdraabab , d3fdrbrbrb ;\n")
    f.write(4*" "+"real d3fdrbrbga , d3fdrbrbgb , d3fdrbrbab , d3fdrbgaga ;\n")
    f.write(4*" "+"real d3fdrbgagb , d3fdrbgaab , d3fdrbgbgb , d3fdrbgbab ;\n")
    f.write(4*" "+"real d3fdrbabab , d3fdgagaga , d3fdgagagb , d3fdgagaab ;\n")
    f.write(4*" "+"real d3fdgagbgb , d3fdgagbab , d3fdgaabab , d3fdgbgbgb ;\n")
    f.write(4*" "+"real d3fdgbgbab , d3fdgbabab , d3fdababab ;\n")

    
    make_stdefs(f)
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0100 += factor*dfdrb;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df0001 += factor*dfdgb;\n")
    f.write(4*" "+"ds->df00001 += factor*dfdab;\n")

    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df1100 += factor*d2fdrarb;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df1001 += factor*d2fdragb;\n")
    f.write(4*" "+"ds->df10001 += factor*d2fdraab;\n")
    f.write(4*" "+"ds->df0200 += factor*d2fdrbrb;\n")
    f.write(4*" "+"ds->df0110 += factor*d2fdrbga;\n")
    f.write(4*" "+"ds->df0101 += factor*d2fdrbgb;\n")
    f.write(4*" "+"ds->df01001 += factor*d2fdrbab;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write(4*" "+"ds->df0011 += factor*d2fdgagb;\n")
    f.write(4*" "+"ds->df00101+= factor*d2fdgaab;\n")
    f.write(4*" "+"ds->df0002 += factor*d2fdgbgb;\n")
    f.write(4*" "+"ds->df00011+= factor*d2fdgbab;\n")
    f.write(4*" "+"ds->df00002+= factor*d2fdabab;\n")

    f.write(4*" "+"ds->df3000 += factor*d3fdrarara;\n")
    f.write(4*" "+"ds->df2100 += factor*d3fdrararb;\n")
    f.write(4*" "+"ds->df2010 += factor*d3fdraraga;\n")
    f.write(4*" "+"ds->df2001 += factor*d3fdraragb;\n")
    f.write(4*" "+"ds->df20001 += factor*d3fdraraab;\n")
    f.write(4*" "+"ds->df1200 += factor*d3fdrarbrb;\n")
    f.write(4*" "+"ds->df1110 += factor*d3fdrarbga;\n")
    f.write(4*" "+"ds->df1101 += factor*d3fdrarbgb;\n")
    f.write(4*" "+"ds->df11001 += factor*d3fdrarbab;\n")
    f.write(4*" "+"ds->df1020 += factor*d3fdragaga;\n")
    f.write(4*" "+"ds->df1011 += factor*d3fdragagb;\n")
    f.write(4*" "+"ds->df10101+= factor*d3fdragaab;\n")
    f.write(4*" "+"ds->df1002 += factor*d3fdragbgb;\n")
    f.write(4*" "+"ds->df10011+= factor*d3fdragbab;\n")
    f.write(4*" "+"ds->df10002+= factor*d3fdraabab;\n")
    f.write(4*" "+"ds->df0300 += factor*d3fdrbrbrb;\n")
    f.write(4*" "+"ds->df0210 += factor*d3fdrbrbga;\n")
    f.write(4*" "+"ds->df0201 += factor*d3fdrbrbgb;\n")
    f.write(4*" "+"ds->df02001 += factor*d3fdrbrbab;\n")
    f.write(4*" "+"ds->df0120 += factor*d3fdrbgaga;\n")
    f.write(4*" "+"ds->df0111 += factor*d3fdrbgagb;\n")
    f.write(4*" "+"ds->df01101+= factor*d3fdrbgaab;\n")
    f.write(4*" "+"ds->df0102 += factor*d3fdrbgbgb;\n")
    f.write(4*" "+"ds->df01011+= factor*d3fdrbgbab;\n")
    f.write(4*" "+"ds->df01002+= factor*d3fdrbgbab;\n")
    f.write(4*" "+"ds->df0030 += factor*d3fdgagaga;\n")
    f.write(4*" "+"ds->df0021 += factor*d3fdgagagb;\n")
    f.write(4*" "+"ds->df00201+= factor*d3fdgagaab;\n")
    f.write(4*" "+"ds->df0012 += factor*d3fdgagbgb;\n")
    f.write(4*" "+"ds->df00111+= factor*d3fdgagbab;\n")
    f.write(4*" "+"ds->df00102+= factor*d3fdgaabab;\n")
    f.write(4*" "+"ds->df0003 += factor*d3fdgbgbgb;\n")
    f.write(4*" "+"ds->df00021+= factor*d3fdgbgbab;\n")
    f.write(4*" "+"ds->df00012+= factor*d3fdgbabab;\n")
    f.write(4*" "+"ds->df00003+= factor*d3fdababab;\n")

    f.write("}\n")
    f.close()    

def make_fourthcode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdrb, dfdga, dfdgb, dfdab;\n")

    f.write(4*" "+"real d2fdrara , d2fdrarb , d2fdraga , d2fdragb ;\n")
    f.write(4*" "+"real d2fdraab , d2fdrbrb , d2fdrbga , d2fdrbgb ;\n")
    f.write(4*" "+"real d2fdrbab , d2fdgaga , d2fdgagb , d2fdgaab ;\n")
    f.write(4*" "+"real d2fdgbgb , d2fdgbab , d2fdabab ;\n")

    f.write(4*" "+"real d3fdrarara , d3fdrararb , d3fdraraga , d3fdraragb ;\n")
    f.write(4*" "+"real d3fdraraab , d3fdrarbrb , d3fdrarbga , d3fdrarbgb ;\n")
    f.write(4*" "+"real d3fdrarbab , d3fdragaga , d3fdragagb , d3fdragaab ;\n")
    f.write(4*" "+"real d3fdragbgb , d3fdragbab , d3fdraabab , d3fdrbrbrb ;\n")
    f.write(4*" "+"real d3fdrbrbga , d3fdrbrbgb , d3fdrbrbab , d3fdrbgaga ;\n")
    f.write(4*" "+"real d3fdrbgagb , d3fdrbgaab , d3fdrbgbgb , d3fdrbgbab ;\n")
    f.write(4*" "+"real d3fdrbabab , d3fdgagaga , d3fdgagagb , d3fdgagaab ;\n")
    f.write(4*" "+"real d3fdgagbgb , d3fdgagbab , d3fdgaabab , d3fdgbgbgb ;\n")
    f.write(4*" "+"real d3fdgbgbab , d3fdgbabab , d3fdababab ;\n")

    f.write(4*" "+"real d4fdrararara , d4fdrarararb , d4fdrararaga , d4fdrararagb ;\n")
    f.write(4*" "+"real d4fdrararaab , d4fdrararbrb , d4fdrararbga , d4fdrararbgb ;\n")
    f.write(4*" "+"real d4fdrararbab , d4fdraragaga , d4fdraragagb , d4fdraragaab ;\n")
    f.write(4*" "+"real d4fdraragbgb , d4fdraragbab , d4fdraraabab , d4fdrarbrbrb ;\n")
    f.write(4*" "+"real d4fdrarbrbga , d4fdrarbrbgb , d4fdrarbrbab , d4fdrarbgaga ;\n")
    f.write(4*" "+"real d4fdrarbgagb , d4fdrarbgaab , d4fdrarbgbgb , d4fdrarbgbab ;\n")
    f.write(4*" "+"real d4fdrarbabab , d4fdragagaga , d4fdragagagb , d4fdragagaab ;\n")
    f.write(4*" "+"real d4fdragagbgb , d4fdragagbab , d4fdragaabab , d4fdragbgbgb ;\n")
    f.write(4*" "+"real d4fdragbgbab , d4fdragbabab , d4fdraababab , d4fdrbrbrbrb ;\n")
    f.write(4*" "+"real d4fdrbrbrbga , d4fdrbrbrbgb , d4fdrbrbrbab , d4fdrbrbgaga ;\n")
    f.write(4*" "+"real d4fdrbrbgagb , d4fdrbrbgaab , d4fdrbrbgbgb , d4fdrbrbgbab ;\n")
    f.write(4*" "+"real d4fdrbrbabab , d4fdrbgagaga , d4fdrbgagagb , d4fdrbgagaab ;\n")
    f.write(4*" "+"real d4fdrbgagbgb , d4fdrbgagbab , d4fdrbgaabab , d4fdrbgbgbgb ;\n")
    f.write(4*" "+"real d4fdrbgbgbab , d4fdrbgbabab , d4fdrbababab , d4fdgagagaga ;\n")
    f.write(4*" "+"real d4fdgagagagb , d4fdgagagaab , d4fdgagagbgb , d4fdgagagbab ;\n")
    f.write(4*" "+"real d4fdgagaabab , d4fdgagbgbgb , d4fdgagbgbab , d4fdgagbabab ;\n")
    f.write(4*" "+"real d4fdgaababab , d4fdgbgbgbgb , d4fdgbgbgbab , d4fdgbgbabab ;\n")
    f.write(4*" "+"real d4fdgbababab , d4fdabababab ;\n")

    make_stdefs(f)
    for i in range(len(code)):
        f.write(4*" "+code[i]+"\n")
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0100 += factor*dfdrb;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df0001 += factor*dfdgb;\n")
    f.write(4*" "+"ds->df00001 += factor*dfdab;\n")

    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df1100 += factor*d2fdrarb;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df1001 += factor*d2fdragb;\n")
    f.write(4*" "+"ds->df10001 += factor*d2fdraab;\n")
    f.write(4*" "+"ds->df0200 += factor*d2fdrbrb;\n")
    f.write(4*" "+"ds->df0110 += factor*d2fdrbga;\n")
    f.write(4*" "+"ds->df0101 += factor*d2fdrbgb;\n")
    f.write(4*" "+"ds->df01001 += factor*d2fdrbab;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write(4*" "+"ds->df0011 += factor*d2fdgagb;\n")
    f.write(4*" "+"ds->df00101 += factor*d2fdgaab;\n")
    f.write(4*" "+"ds->df0002 += factor*d2fdgbgb;\n")
    f.write(4*" "+"ds->df00011 += factor*d2fdgbab;\n")
    f.write(4*" "+"ds->df00002 += factor*d2fdabab;\n")

    f.write(4*" "+"ds->df3000 += factor*d3fdrarara;\n")
    f.write(4*" "+"ds->df2100 += factor*d3fdrararb;\n")
    f.write(4*" "+"ds->df2010 += factor*d3fdraraga;\n")
    f.write(4*" "+"ds->df2001 += factor*d3fdraragb;\n")
    f.write(4*" "+"ds->df20001 += factor*d3fdraraab;\n")
    f.write(4*" "+"ds->df1200 += factor*d3fdrarbrb;\n")
    f.write(4*" "+"ds->df1110 += factor*d3fdrarbga;\n")
    f.write(4*" "+"ds->df1101 += factor*d3fdrarbgb;\n")
    f.write(4*" "+"ds->df11001 += factor*d3fdrarbab;\n")
    f.write(4*" "+"ds->df1020 += factor*d3fdragaga;\n")
    f.write(4*" "+"ds->df1011 += factor*d3fdragagb;\n")
    f.write(4*" "+"ds->df10101 += factor*d3fdragaab;\n")
    f.write(4*" "+"ds->df1002 += factor*d3fdragbgb;\n")
    f.write(4*" "+"ds->df10011 += factor*d3fdragbab;\n")
    f.write(4*" "+"ds->df10002 += factor*d3fdraabab;\n")
    f.write(4*" "+"ds->df0300 += factor*d3fdrbrbrb;\n")
    f.write(4*" "+"ds->df0210 += factor*d3fdrbrbga;\n")
    f.write(4*" "+"ds->df0201 += factor*d3fdrbrbgb;\n")
    f.write(4*" "+"ds->df02001 += factor*d3fdrbrbab;\n")
    f.write(4*" "+"ds->df0120 += factor*d3fdrbgaga;\n")
    f.write(4*" "+"ds->df0111 += factor*d3fdrbgagb;\n")
    f.write(4*" "+"ds->df01101 += factor*d3fdrbgaab;\n")
    f.write(4*" "+"ds->df0102 += factor*d3fdrbgbgb;\n")
    f.write(4*" "+"ds->df01011 += factor*d3fdrbgbab;\n")
    f.write(4*" "+"ds->df01002 += factor*d3fdrbabab;\n")
    f.write(4*" "+"ds->df0030 += factor*d3fdgagaga;\n")
    f.write(4*" "+"ds->df0021 += factor*d3fdgagagb;\n")
    f.write(4*" "+"ds->df00201 += factor*d3fdgagaab;\n")
    f.write(4*" "+"ds->df0012 += factor*d3fdgagbgb;\n")
    f.write(4*" "+"ds->df00111 += factor*d3fdgagbab;\n")
    f.write(4*" "+"ds->df00102 += factor*d3fdgaabab;\n")
    f.write(4*" "+"ds->df0003 += factor*d3fdgbgbgb;\n")
    f.write(4*" "+"ds->df00021 += factor*d3fdgbgbab;\n")
    f.write(4*" "+"ds->df00012 += factor*d3fdgbabab;\n")
    f.write(4*" "+"ds->df00003 += factor*d3fdababab;\n")

    f.write(4*" "+"ds->df4000 += factor*d4fdrararara;\n")
    f.write(4*" "+"ds->df3100 += factor*d4fdrarararb;\n")
    f.write(4*" "+"ds->df3010 += factor*d4fdrararaga;\n")
    f.write(4*" "+"ds->df3001 += factor*d4fdrararagb;\n")
    f.write(4*" "+"ds->df30001 += factor*d4fdrararaab;\n")
    f.write(4*" "+"ds->df2200 += factor*d4fdrararbrb;\n")
    f.write(4*" "+"ds->df2110 += factor*d4fdrararbga;\n")
    f.write(4*" "+"ds->df2101 += factor*d4fdrararbgb;\n")
    f.write(4*" "+"ds->df21001 += factor*d4fdrararbab;\n")
    f.write(4*" "+"ds->df2020 += factor*d4fdraragaga;\n")
    f.write(4*" "+"ds->df2011 += factor*d4fdraragagb;\n")
    f.write(4*" "+"ds->df20101 += factor*d4fdraragaab;\n")
    f.write(4*" "+"ds->df2002 += factor*d4fdraragbgb;\n")
    f.write(4*" "+"ds->df20011 += factor*d4fdraragbab;\n")
    f.write(4*" "+"ds->df20002 += factor*d4fdraraabab;\n")
    f.write(4*" "+"ds->df1300 += factor*d4fdrarbrbrb;\n")
    f.write(4*" "+"ds->df1210 += factor*d4fdrarbrbga;\n")
    f.write(4*" "+"ds->df1201 += factor*d4fdrarbrbgb;\n")
    f.write(4*" "+"ds->df12001 += factor*d4fdrarbrbab;\n")
    f.write(4*" "+"ds->df1120 += factor*d4fdrarbgaga;\n")
    f.write(4*" "+"ds->df1111 += factor*d4fdrarbgagb;\n")
    f.write(4*" "+"ds->df11101 += factor*d4fdrarbgaab;\n")
    f.write(4*" "+"ds->df1102 += factor*d4fdrarbgbgb;\n")
    f.write(4*" "+"ds->df11011 += factor*d4fdrarbgbab;\n")
    f.write(4*" "+"ds->df11002 += factor*d4fdrarbabab;\n")
    f.write(4*" "+"ds->df1030 += factor*d4fdragagaga;\n")
    f.write(4*" "+"ds->df1021 += factor*d4fdragagagb;\n")
    f.write(4*" "+"ds->df10201 += factor*d4fdragagaab;\n")
    f.write(4*" "+"ds->df1012 += factor*d4fdragagbgb;\n")
    f.write(4*" "+"ds->df10111 += factor*d4fdragagbab;\n")
    f.write(4*" "+"ds->df10102 += factor*d4fdragaabab;\n")
    f.write(4*" "+"ds->df1003 += factor*d4fdragbgbgb;\n")
    f.write(4*" "+"ds->df10021 += factor*d4fdragbgbab;\n")
    f.write(4*" "+"ds->df10012 += factor*d4fdragbabab;\n")
    f.write(4*" "+"ds->df10003 += factor*d4fdraababab;\n")
    f.write(4*" "+"ds->df0400 += factor*d4fdrbrbrbrb;\n")
    f.write(4*" "+"ds->df0310 += factor*d4fdrbrbrbga;\n")
    f.write(4*" "+"ds->df0301 += factor*d4fdrbrbrbgb;\n")
    f.write(4*" "+"ds->df03001 += factor*d4fdrbrbrbab;\n")
    f.write(4*" "+"ds->df0220 += factor*d4fdrbrbgaga;\n")
    f.write(4*" "+"ds->df0211 += factor*d4fdrbrbgagb;\n")
    f.write(4*" "+"ds->df02101 += factor*d4fdrbrbgaab;\n")
    f.write(4*" "+"ds->df0202 += factor*d4fdrbrbgbgb;\n")
    f.write(4*" "+"ds->df02011 += factor*d4fdrbrbgbab;\n")
    f.write(4*" "+"ds->df02002 += factor*d4fdrbrbabab;\n")
    f.write(4*" "+"ds->df0130 += factor*d4fdrbgagaga;\n")
    f.write(4*" "+"ds->df0121 += factor*d4fdrbgagagb;\n")
    f.write(4*" "+"ds->df01201 += factor*d4fdrbgagaab;\n")
    f.write(4*" "+"ds->df0112 += factor*d4fdrbgagbgb;\n")
    f.write(4*" "+"ds->df01111 += factor*d4fdrbgagbab;\n")
    f.write(4*" "+"ds->df01102 += factor*d4fdrbgaabab;\n")
    f.write(4*" "+"ds->df0103 += factor*d4fdrbgbgbgb;\n")
    f.write(4*" "+"ds->df01021 += factor*d4fdrbgbgbab;\n")
    f.write(4*" "+"ds->df01012 += factor*d4fdrbgbabab;\n")
    f.write(4*" "+"ds->df01003 += factor*d4fdrbababab;\n")
    f.write(4*" "+"ds->df0040 += factor*d4fdgagagaga;\n")
    f.write(4*" "+"ds->df0031 += factor*d4fdgagagagb;\n")
    f.write(4*" "+"ds->df00301 += factor*d4fdgagagaab;\n")
    f.write(4*" "+"ds->df0022 += factor*d4fdgagagbgb;\n")
    f.write(4*" "+"ds->df00211 += factor*d4fdgagagbab;\n")
    f.write(4*" "+"ds->df00202 += factor*d4fdgagaabab;\n")
    f.write(4*" "+"ds->df0013 += factor*d4fdgagbgbgb;\n")
    f.write(4*" "+"ds->df00121 += factor*d4fdgagbgbab;\n")
    f.write(4*" "+"ds->df00112 += factor*d4fdgagbabab;\n")
    f.write(4*" "+"ds->df00103 += factor*d4fdgaababab;\n")
    f.write(4*" "+"ds->df0004 += factor*d4fdgbgbgbgb;\n")
    f.write(4*" "+"ds->df00031 += factor*d4fdgbgbgbab;\n")
    f.write(4*" "+"ds->df00022 += factor*d4fdgbgbabab;\n")
    f.write(4*" "+"ds->df00013 += factor*d4fdgbababab;\n")
    f.write(4*" "+"ds->df00004 += factor*d4fdabababab;\n")

    f.write("}\n")
    f.close()    
            

#---------------------------------------------------------------------
parts=[ ['en',     make_energycode], ['first', make_firstcode ],
        ['second', make_secondcode], ['third', make_thirdcode ],
        ['fourth', make_fourthcode]]

if __name__ == "__main__":
    for p in parts:
        make_maximafile(funcfil,p[0])
        run_maxima(p[0]+'.out')
        #nvars,code=read_maximaout(p[0]+'.out')
        nvars,code=max2c(parse_maxima(p[0]+'.out'))
        apply(p[1], [nvars, code])
        os.remove(p[0]+'.out')

