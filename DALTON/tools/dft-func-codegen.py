#!/usr/bin/env python2
#----------------------------------
# generate dft(gga) code for dirac
# Andreas Hesselmann <andreas@theochem.uni-duesseldorf.de>, 2003
#----------------------------------
# needs: mupad 
# usage: codegen.py 'mupadfile'
# 'mupadfile' is a  mupad input file of the following general form:
#      ....some variables....
#      ....some functions...
#      K(rhoa,ghroa,rhob,gradb):= 'expression of rhoa,grhoa,xa
#                                                rhob,grhob,xb'
# where xa=sqrt(grada*grada)/rhoa^(4/3)
# output: fun-'mupadfile'.c

import sys, string, os, re;

funcfil=sys.argv[1];
funcname=string.split(funcfil,'.')[0]

#header
f=open(funcfil,'r')
func=f.readlines()
f.close()
f=open("fun-"+funcname+".c",'w')
f.write("/* Automatically generated functional code: "+funcname+"\n")
f.write("   Mupad input:\n")
for line in func:
    f.write(4*" "+">> "+line)
f.write("*/\n\n")
f.write("// add \"extern Functional "+funcname+"Functional;\" to 'functionals.h'\n")
f.write("// add \"&"+funcname+"Functional,\" to 'functionals.c'\n\n")
f.write("#include <math.h>\n")
f.write("#include <stddef.h>\n\n")
f.write("#define __CVERSION__\n\n")
f.write("#include \"functionals.h\"\n")
f.write("#include \"dcbham.h\"\n\n")
f.write("/* INTERFACE PART */\n")
f.write("static int "+funcname+"_read(const char* conf_line);\n")
f.write("static real "+funcname+"_energy(real rhoa, real rhob, real grada, real gradb);\n")
f.write("static void "+funcname+"_first(FirstDrv *ds, real factor, \n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n")
f.write("static void "+funcname+"_second(SecondFuncDrv *ds, real factor,\n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n")
f.write("static void "+funcname+"_third(ThirdFuncDrv *ds, real factor,\n\
                        real rhoa, real rhob, real ngrada, real ngradb);\n\n")
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
f.write("{\n#ifndef NO_BACKWARD_COMP\n  dcrham_.hfxfac = 0;\n")
f.write("#endif\n  return 1;\n}\n\n")
f.close()

#mupad commands
com1="print(Unquoted,generate::C(generate::optimize(["
com2="]))):\n"
mupad_commands={
    "en": com1+"zk=K(rhoa,grada,rhob,gradb)"+com2,
    "first": com1+"dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),\
                   dfdga=diff(K(rhoa,grada,rhob,gradb),grada)"+com2,
    "second":com1+"dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),\
                   dfdga=diff(K(rhoa,grada,rhob,gradb),grada),\
                   d2fdraga=diff(K(rhoa,grada,rhob,gradb),rhoa,grada),\
                   d2fdrara=diff(K(rhoa,grada,rhob,gradb),rhoa,rhoa),\
                   d2fdgaga=diff(K(rhoa,grada,rhob,gradb),grada,grada)"+com2,
    "third": com1+"dfdra=diff(K(rhoa,grada,rhob,gradb),rhoa),\
                   dfdga=diff(K(rhoa,grada,rhob,gradb),grada),\
                   d2fdraga=diff(K(rhoa,grada,rhob,gradb),rhoa,grada),\
                   d2fdrara=diff(K(rhoa,grada,rhob,gradb),rhoa,rhoa),\
                   d2fdgaga=diff(K(rhoa,grada,rhob,gradb),grada,grada),\
                   d3fdraraga=diff(K(rhoa,grada,rhob,gradb),rhoa,rhoa,grada),\
                   d3fdragaga=diff(K(rhoa,grada,rhob,gradb),rhoa,grada,grada),\
                   d3fdrarara=diff(K(rhoa,grada,rhob,gradb),rhoa,rhoa,rhoa),\
                   d3fdgagaga=diff(K(rhoa,grada,rhob,gradb),grada,grada,grada)"+com2
    }
                   

def make_mupadfile(funcfil,com):
    f=open(funcfil,'r')
    func=f.readlines()
    f.close()
    f=open('mupad.in','w')
    f.write("xa:=sqrt(grada*grada)/rhoa^(4/3):\n");
    f.write("xb:=sqrt(gradb*gradb)/rhob^(4/3):\n");
    for line in func:
        f.write(line)
    f.write(mupad_commands[com])
    f.close()

def run_mupad(outfil):
    if os.system("mupad -b "+outfil+" mupad.in")!=0:
        sys.exit("Error in run_mupad!")


def read_mupadout(outfil):
    f=open(outfil,'r')
    var=[]
    expression=[]
    i=-1
    for line in f.readlines():
        if line=="\n" : continue
        if string.find(line,'=')>-1:
            i=i+1
            var.append("")
            expression.append("")
            var[i],expression[i]=string.split(line,'=')
            var[i]=string.strip(var[i])
        elif line[-1]=="\\":
            expression[i]=expression[i]+line
        else:
            expression[i]=expression[i]+line
    f.close()

    # remove \\n 
    for i in range(len(expression)):
        expression[i]=string.replace(expression[i],'\\\n','')

    vardict={}
    i=0
    for v in var:
        if v[0]=='t':
            vardict[v]='t['+str(i)+']'
            i=i+1
        else:
            vardict[v]=v
    pat=re.compile(r'(t[\d]+)') #find all variables t...
    for i in range(len(expression)):
        l=re.findall(pat,expression[i])
        for j in range(len(l)):
            expression[i]=string.replace(expression[i],l[j],vardict[l[j]])
    return var,vardict,expression


def make_energycode(var,vardict,expression):
    l=var.index('zk')+1
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic real\n")
    f.write(funcname+"_energy(real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(l-1)+"];\n")
    for i in range(l):
        if var[i]=="zk":
            f.write(4*" "+"return "+expression[i])
            break
        f.write(4*" "+vardict[var[i]]+" = "+expression[i])
    f.write("}\n")
    f.close()

def make_firstcode(var,vardict,expression):
    l=var.index('dfdga')+1
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_first(FirstDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(l-2)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    for i in range(l):
        f.write(4*" "+vardict[var[i]]+" = "+expression[i])
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write("}\n")
    f.close()    

def make_secondcode(var,vardict,expression):
    l=var.index('d2fdgaga')+1
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_second(SecondFuncDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(l-4)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdgaga;\n")    
    for i in range(l):
        f.write(4*" "+vardict[var[i]]+" = "+expression[i])
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write("}\n")
    f.close()    

def make_thirdcode(var,vardict,expression):
    l=var.index('d3fdgagaga')+1
    f=open("fun-"+funcname+".c",'a')
    f.write("\nstatic void\n")
    f.write(funcname+"_third(ThirdFuncDrv *ds, real factor, real rhoa, real rhob, real grada, real gradb)\n")
    f.write("{\n")
    f.write(4*" "+"real t["+str(l-8)+"];\n")
    f.write(4*" "+"real dfdra, dfdga;\n")
    f.write(4*" "+"real d2fdraga, d2fdrara, d2fdgaga;\n")
    f.write(4*" "+"real d3fdraraga, d3fdragaga, d3fdrarara, d3fdgagaga;\n")  
    for i in range(l):
        f.write(4*" "+vardict[var[i]]+" = "+expression[i])
    f.write(4*" "+"ds->df1000 += factor*dfdra;\n")
    f.write(4*" "+"ds->df0010 += factor*dfdga;\n")
    f.write(4*" "+"ds->df1010 += factor*d2fdraga;\n")
    f.write(4*" "+"ds->df2000 += factor*d2fdrara;\n")
    f.write(4*" "+"ds->df0020 += factor*d2fdgaga;\n")
    f.write(4*" "+"ds->df2010 += factor*d3fdraraga;\n")
    f.write(4*" "+"ds->df1020 += factor*d3fdragaga;\n")
    f.write(4*" "+"ds->df3000 += factor*d3fdrarara;\n")
    f.write(4*" "+"ds->df0030 += factor*d3fdgagaga;\n")
    f.write("}\n")
    f.close()    
            
            

#---------------------------------------------------------------------
parts=['en','first','second','third']
for p in parts:
    make_mupadfile(funcfil,p)
    run_mupad(p+'.out')
    var,vardict,expression=read_mupadout(p+'.out')
    if p=='en':
        make_energycode(var,vardict,expression)
    elif p=='first':
        make_firstcode(var,vardict,expression)
    elif p=='second':
        make_secondcode(var,vardict,expression)
    elif p=='third':
        make_thirdcode(var,vardict,expression)
    os.remove(p+'.out')

