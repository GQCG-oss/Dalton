#!/usr/bin/env python2
#----------------------------------
# generate dft(gga) code for dirac
# Andreas Hesselmann <andreas@theochem.uni-duesseldorf.de>, 2003
# Adapted to new density property input and unrestricted functionals
# by Pawel Salek.
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
f.write("static int "+funcname+"_read(const char* conf_line);\n")
f.write("static real "+funcname+"_energy(const DftDensProp* dp);\n")
f.write("static void "+funcname+"_first(FirstFuncDrv *ds, real factor, \n\
                       const DftDensProp* dp);\n")
f.write("static void "+funcname+"_second(SecondFuncDrv *ds, real factor,\n\
                        const DftDensProp* dp);\n")
f.write("static void "+funcname+"_third(ThirdFuncDrv *ds, real factor,\n\
                       const DftDensProp* dp);\n\n")
f.write("//static int fun_true(void) { return 1; }\n")
f.write("Functional "+funcname+"Functional = {\n")
f.write("  \""+funcname+"\",\n")
f.write("  fun_true,\n")
f.write("  "+funcname+"_read,\n")
f.write("  NULL,\n")
f.write("  "+funcname+"_energy,\n")
f.write("  "+funcname+"_first,\n")
f.write("  "+funcname+"_second,\n")
f.write("  "+funcname+"_third\n};\n\n")
f.write("/* IMPLEMENTATION PART */\n")
f.write("static int\n"+funcname+"_read(const char* conf_line)\n")
f.write("{\n    dft_set_hf_weight(0);\n")
f.write("    return 1;\n}\n\n")
f.close()

#maxima commands
com1="string(float(subst(pow,\"^\",optimize("
com2="))));\n"
maxima_commands={
    "en": "zk=K(rhoa,grada,rhob,gradb,gradab)",
    "first": "[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa),\n"
            +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob),\n"
            +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada),\n"
            +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb),\n"
            +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab)]\n",
    "second":"[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa),\n"
            +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob),\n"
            +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada),\n"
            +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb),\n"
            +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab),\n"
            +" d2fdrara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2),"
            +" d2fdrarb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1),"
            +" d2fdraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1),"
            +" d2fdragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1),"
            +" d2fdraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,1),"
            +" d2fdrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2),"
            +" d2fdrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1),"
            +" d2fdrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1),"
            +" d2fdrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,1),"
            +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2),"
            +" d2fdgagb=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradb,1),"
            +" d2fdgaab=diff(K(rhoa,grada,rhob,gradb,gradab),grada,1,gradab,1),"
            +" d2fdgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2),"
            +" d2fdgbab=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,1,gradab,1),"
            +" d2fdabab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab,2)]",
    "third": "[dfdra=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa),"
            +" dfdrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob),\n"
            +" dfdga=diff(K(rhoa,grada,rhob,gradb,gradab),grada),\n"
            +" dfdgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb),\n"
            +" dfdab=diff(K(rhoa,grada,rhob,gradb,gradab),gradab),\n"
            +" d2fdrara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2),"
            +" d2fdrarb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1),"
            +" d2fdraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,1),"
            +" d2fdragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,1),"
            +" d2fdrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2),"
            +" d2fdraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradab,1),"
            +" d2fdrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradab,1),"
            +" d2fdgaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,2),"
            +" d2fdgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),gradb,2),"
            +" d2fdrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,1),"
            +" d2fdrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,1),"
            +" d3fdrararb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,rhob,1),"
            +" d3fdraraga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,grada,1),"
            +" d3fdraragb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradb,1),"
            +" d3fdrbrbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradab,1),"
            +" d3fdraraab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,2,gradab,1),"
            +" d3fdrarbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,2),"
            +" d3fdrarbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,grada,1),"
            +" d3fdrarbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradb,1),"
            +" d3fdrarbab=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,rhob,1,gradab,1),"
            +" d3fdragaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,grada,2),"
            +" d3fdragbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,1,gradb,2),"
           +" d3fdrarara=diff(K(rhoa,grada,rhob,gradb,gradab),rhoa,3),"
           +" d3fdrbrbrb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,3),"
           +" d3fdrbrbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,grada,1),"
           +" d3fdrbrbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,2,gradb,1),"
           +" d3fdrbgaga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2),"
           +" d3fdrbgbga=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,grada,2),"
           +" d3fdrbgbgb=diff(K(rhoa,grada,rhob,gradb,gradab),rhob,1,gradb,2),"
           +" d3fdgagaga=diff(K(rhoa,grada,rhob,gradb,gradab),grada,3)]"
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
    if os.system("maxima < maxima.in >"+outfil)!=0:
        sys.exit("Error in run_maxima!")

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
    f=open(outfil,'r')
    block=0
    for line in f.readlines():
        data=string.split(line)
        if len(data)==0: continue
        if data[0]=="\n": continue
        if len(data)>=2:
            if len(data[1])>=5:
                if data[1][0:5]=="BLOCK":
                    block=1
                    expr=""
        if block==1:
            if data[0][0:2]!="(C":
                line=string.strip(line)
                line=string.replace(line,"#","")
                line=string.replace(line,"\n","")
                expr=expr+line
            else:
                block=0
    nvars,cexpr=f2c(expr)
    return nvars,cexpr

def make_stdefs(f):
    f.write(4*" "+"real rhoa = dp->rhoa;\n"+
            4*" "+"real rhob = dp->rhob;\n"+
            4*" "+"real grada = dp->grada;\n"+
            4*" "+"real gradb = dp->gradb;\n"+
            4*" "+"real gradab = dp->gradab;\n\n")



def make_energycode(nvars,code):
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic real\n")
    f.write(funcname+"_energy(const DftDensProp* dp)\n")
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
    f.write(funcname+"_first(FirstFuncDrv *ds, real factor, const DftDensProp* dp)\n")
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
    f.write(funcname+"_second(SecondFuncDrv *ds, real factor, const DftDensProp* dp)\n")
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
    f.write(funcname+"_third(ThirdFuncDrv *ds, real factor, const DftDensProp* dp)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(nvars+1)+"];\n")
    f.write(4*" "+"real dfdra, dfdrb, dfdga, dfdgb, dfdab;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;\n")
    f.write(4*" "+"real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;\n")
    f.write(4*" "+"real d2fdraab, d2fdrbab;\n")
    f.write(4*" "+"real d3fdraraga, d3fdraragb, d3fdraraab, d3fdrbrbab;\n")  
    f.write(4*" "+"real d3fdrarara, d3fdrararb, d3fdragaga, d3fdrarbrb;\n")  
    f.write(4*" "+"real d3fdragbgb, d3fdrarbgb, d3fdrarbab, d3fdgagaga;\n")
    f.write(4*" "+"real d3fdrbrbrb, d3fdrbrbga, d3fdrbrbgb, d3fdrbgbgb;\n")
    f.write(4*" "+"real d3fdrbgbga, d3fdrarbga, d3fdrbgaga;\n")
    
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
    f.write(4*" "+"ds->df0002 += factor*d2fdgbgb;\n")
    f.write(4*" "+"ds->df2010 += factor*d3fdraraga;\n")
    f.write(4*" "+"ds->df2001 += factor*d3fdraragb;\n")
    f.write(4*" "+"ds->df1101 += factor*d3fdrarbgb;\n")
    f.write(4*" "+"ds->df11001 += factor*d3fdrarbab;\n")
    f.write(4*" "+"ds->df1020 += factor*d3fdragaga;\n")
    f.write(4*" "+"ds->df1002 += factor*d3fdragbgb;\n")
    f.write(4*" "+"ds->df3000 += factor*d3fdrarara;\n")
    f.write(4*" "+"ds->df2100 += factor*d3fdrararb;\n")
    f.write(4*" "+"ds->df20001 += factor*d3fdraraab;\n")
    f.write(4*" "+"ds->df02001 += factor*d3fdrbrbab;\n")
    f.write(4*" "+"ds->df1200 += factor*d3fdrarbrb;\n")
    f.write(4*" "+"ds->df1110 += factor*d3fdrarbga;\n")
    f.write(4*" "+"ds->df0300 += factor*d3fdrbrbrb;\n")
    f.write(4*" "+"ds->df0210 += factor*d3fdrbrbga;\n")
    f.write(4*" "+"ds->df0201 += factor*d3fdrbrbgb;\n")
    f.write(4*" "+"ds->df0120 += factor*d3fdrbgaga;\n")
    f.write(4*" "+"ds->df0102 += factor*d3fdrbgbgb;\n")
    f.write(4*" "+"ds->df0030 += factor*d3fdgagaga;\n")
    f.write("}\n")
    f.close()    
            
            

#---------------------------------------------------------------------
parts=[ ['en',     make_energycode], ['first', make_firstcode ],
        ['second', make_secondcode], ['third', make_thirdcode ]]
for p in parts:
    make_maximafile(funcfil,p[0])
    run_maxima(p[0]+'.out')
    nvars,code=read_maximaout(p[0]+'.out')
    apply(p[1], [nvars, code])
    os.remove(p[0]+'.out')

