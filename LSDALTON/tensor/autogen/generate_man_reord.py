#GENERATING TENSOR PERMUTATIONS
#AUTHOR: PATRICK ETTENHUBER
#EMAIL : pett@chem.au.dk, pettenhuber@gmail.com
#DATE  : JUNE, 2013
import sys,datetime,os,time
from math import factorial
from random import randrange
"""This file is inteded for the automatic generation of a data sorting module in LSDALTON
As data sorting is under constant development and features many lines of code this script was
written for easy modification of the code. Furhtermore, one might note that the structure
of optimal code might depend on the architecture which can be at some point introduced here
as well

Currently the routine does not take any command-line arguments, but is dependet on the LSDALTON
folder structure, i.e. it needs LSDALTON/tensor and LSDALTON/tensor/autogen/reord_headers.F90
which can be found as long as LSDALTON is in the execution path, or one layer below execution
folder."""

##################################################################################################
##################################################################################################

def main():
  #SET ARGUMENT LIST, MAKE SURE THE WRITING HERE CORRESPONDS TO READING IT IN produce_file
  args = []
  args.append(sys.argv[0])
  args.append(False)
  args.append(False)
  args.append("")
  args.append(False)
  args.append(False)
# print sys.argv
  force_rewrite = False 
  acc_write = False
  real_sp_write = False

  lsdalton="/LSDALTON"
  tensor=lsdalton+"/tensor"


  for i in range(len(sys.argv)):
    if "debug_version" in sys.argv[i]:
      args[1] = True
    if "nocollapse" in sys.argv[i]:
      args[2] = True
    if "CMAKE_BUILD=" in sys.argv[i]:
      args[3] = sys.argv[i][sys.argv[i].find("=")+1:]
    if "acc" in sys.argv[i]:
      args[4] = True
      acc_write = True
    if "real_sp" in sys.argv[i]:
      args[5] = True
      real_sp_write = True

    if "FORCE_REWRITE" in sys.argv[i] :
      force_rewrite = True
    
  #GET THE FOLDER TO STORE THE manual_reorderings.F90
  cwd = os.getcwd()
  tensordir = args[0]
  if (tensor in tensordir ):

    tensordir = tensordir[0:tensordir.find(tensor)]+tensor+"/"

  else:

    for paths,dirs,files in os.walk(cwd+"/.."+lsdalton):
      if(tensor in paths):
        tensordir = paths[0:paths.find(tensor)]+tensor+"/"
        break

    if not os.path.exists(tensordir):
      print "COULD NOT FIND "+tensor+", exiting"
      sys.exit()


  installdir  = ""
  if args[3] != "" :
     installdir = args[3] +"/"
  else:
     installdir = tensordir

  # if output path does not exist, create it
  if not os.path.exists(installdir):
      os.mkdir(installdir)

  #DEFAULT OF WRITING THE FILES IS FALSE, CHECK THE PREVIOUS VARS IF A NEW PRODUCTION IS NECESSARY
  writenew = False

  #OVERRIDE OPTION
  if(force_rewrite):
    writenew = True

  #CHECK IF THE FILES EXIST AT ALL, IF NOT, WRITE
  if(not os.path.exists(installdir+"reorder_frontend.F90")):
    writenew = True
  else:
#    reordmod  = time.ctime(os.path.getmtime(installdir+"reorder_frontend.F90"))
    reordmod  = os.path.getmtime(installdir+"reorder_frontend.F90")
#    scriptmod = time.ctime(os.path.getmtime(sys.argv[0]))
    scriptmod = os.path.getmtime(sys.argv[0])
    if(os.path.getmtime(tensordir+"/autogen/reorder_header.F90")>scriptmod):
       scriptmod = os.path.getmtime(tensordir+"/autogen/reorder_header.F90")
  # radovan: this should not be done here but taken care of by (c)make
    if(scriptmod>reordmod or writenew):
      print "REORDER GENERATOR IS NEWER THAN REORDERING FILES - GENERATING NEW ONES"
      writenew = True


  #RANGE INDICES
  maxr = 4
  minr = 2

  #SET THE NAMES FOR THE FILES TO GENERATE
  num_types  = 6
  names      = [ [] for i in range(num_types)]
  names[0] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_reord.F90" for j in range(factorial(i + minr)) ] for i in range(maxr-minr+1) ]
  names[1] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_utils_f2t.F90" for j in range(factorial(i + minr))] for i in range(maxr-minr+1)]
  names[2] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_utils_t2f.F90" for j in range(factorial(i + minr))] for i in range(maxr-minr+1)]

  if(real_sp_write):
     names[3] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_reord_sp.F90" for j in range(factorial(i + minr)) ] for i in range(maxr-minr+1) ]
  if(acc_write):
     names[4] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_acc_reord.F90" for j in range(factorial(i + minr)) ] for i in range(maxr-minr+1) ]
  if(acc_write and real_sp_write):
     names[5] = [[ installdir+"reord"+str(i+minr)+"d_"+str(j+1)+"_acc_reord_sp.F90" for j in range(factorial(i + minr)) ] for i in range(maxr-minr+1) ]
  
  names_flat = [ path for sublist in names for subsublist in sublist for  path in subsublist]

  #check if all necessary files exist
  for path in names_flat:
     write_new_old = writenew
     if(not os.path.exists(path)):
        writenew = True
     if(writenew and not write_new_old):
        print "SOME FILES MISSING, GENERATING FROM SCRATCH"

  if(not writenew):
    c = open(installdir+"reorder_frontend.F90",'r')
    endvars_found = False
    #PARSE THE LINES TO CHECK WHETHER IT IS NECESSARY TO WRITE THE FILES FROM SCRATCH
    for line in c:
      line = line.strip()
 
      if "!ARG0:" in line:
        writenew = (not line.split()[-1] == str(args[0]))
        if writenew:
          break
 
      if "!ARG1:" in line:
        writenew = (not line.split()[-1] == str(args[1]))
        if writenew:
          break
 
      if "!ARG2:" in line:
        writenew = (not line.split()[-1] == str(args[2]))
        if writenew:
          break
 
      if "!ARG3:" in line:
        writenew = (not line.split()[-1] == str(args[3]))
        if writenew:
          break

      if "!ARG4:" in line:
        writenew = (not line.split()[-1] == str(args[4]))
        if writenew:
          break

      if "!ARG5:" in line:
        writenew = (not line.split()[-1] == str(args[5]))
        if writenew:
          break
 
      if "!END VARS" in line:
        endvars_found = True
        break
 
    c.close()
  
    #ASSUME THAT IF "!END VARS" IS NOT FOUND THE FILE IS AN OLD OR DAMAGED VERSION
    if not endvars_found:
      writenew = True


  if writenew:
    produce_files(installdir,tensordir,args,names,minr,maxr)

##################################################################################################
##################################################################################################
def produce_files(installdir,tensordir,args,names,minr,maxr):
   now         = datetime.datetime.now()

   #GET COMMAND LINE ARGUMENTS
   debug_loops = args[1]
   nocollapse  = args[2]
   acc         = args[4]
   real_sp     = args[5]

   interface_types = [ [] for i in range(len(names)) ]

   #SET THE NAMES FOR THE FILES TO GENERATE
   for iname,name in enumerate(names):

      #ONLY DO THE ITERATION ID THE NAME LIST IS NOT EMPTY
      if( len(name) > 0):

         #iname = 0 : normal full dense reorderings
         #iname = 1 : reoderings full to tiled
         #iname = 2 : reoderings tiled to full
         #iname = 3 : single precision reorderings
         #iname = 4 : acc dense reorderings
         #iname = 5 : acc single precision reorderings
         if( iname == 0 ):
            acc = False
            addition = "dp"
         elif (iname == 1):
            acc = False
            addition = "f2t"
         elif (iname == 2):
            acc = False
            addition = "t2f"
         elif (iname == 3):
            acc = False
            addition = "sp"
         elif (iname == 4):
            acc = True
            addition = "dp"
         elif (iname == 5):
            acc = True
            addition = "sp"
         else:
            print "MAKE SURE YOUR CASE IS IMPLEMENTED AND SET THE SWITCHES CORRECTLY"
            exit(0)

         #OPEN THE FILES TO GENERATE
         reord = [[ open(name[i][j],'w')        for j in range(factorial(i + minr)) ] for i in range(maxr-minr+1) ]

         #WRITE THE HEADERS OF ALL FILES
         for idx in range(maxr-minr+1):
           for i in range(factorial(minr+idx)):
             write_simple_module_header(reord[idx][i],idx+minr,i+1,now,args,acc)


         #SPECIFY THE ORDER OF REODERINGS
         for idx in range(maxr-minr+1):
           modes = idx + minr
           # GENERATE ORIGINAL ORDER AND STARTING POINT FOR NEW ORDER
           idxarr = [ i for i in range(modes) ]

           #GET ALL PERMUTATION
           all_permutations=permutations(idxarr)
          
           #LOOP OVER PERMUTATIONS AND WRITE SUBROUTINES
           for pnum,perm in enumerate(all_permutations):
             #CHECK IF THE ORDERING IS NECESSARY, ELSE JUST SKIP IT
             doreord = True
             for i in range(len(perm)-1):
               if(perm[i]+1==perm[i+1]):
                 doreord = False
             #always write the file for
             doreord = doreord or i

             if doreord :
               # WRITE THE CPU SUBROUTINE HEADER AND GET ITS NAME
               sub = write_subroutine_header(reord[idx][pnum],idxarr,perm,now,modes,addition,debug_loops,acc)
               #Write the CPU subroutine body
               write_subroutine_body(reord[idx][pnum],idxarr,perm,modes,args,addition,acc)
               #END THE CPU SUBROUTINE
               reord[idx][pnum].write("  end subroutine "+sub+"\n\n")
               write_simple_module_end_and_close(reord[idx][pnum])

               #skip interfaces for f2t and t2f
               if iname!=1 and iname!=2 :
                  interface_types[iname].append(sub)

             else:
               reord[idx][pnum].write("\n  subroutine dummy()\n  end subroutine dummy\n\n")
               write_simple_module_end_and_close(reord[idx][pnum])

   # WRITE REORDERINGS FRONTEND FILE
   f=open(installdir+"reorder_frontend.F90",'w')
   all_names = [get_namestub_from_path(item) for sublist in names for subsublist in sublist for  item in subsublist]
   write_main_header(f,now,all_names,args,tensordir,minr,maxr,interface_types)
   f.write("\nend module reorder_frontend_module")
   f.close()

   #TODO: MAKE TESTER INTERFACE TEST T2F, F2T, ACC AND SP VERSIONS
   write_testing_framework(installdir,minr,maxr)


def write_subroutine_body(f,idxarr,perm,modes,args,ad,acc):
  debug_loops = args[1]
  nocollapse  = args[2]

  #GENERAL CASE 
  if(not debug_loops or acc):
    if ad == "sp":
       prec = "tensor_sp"
    else:
       prec = "tensor_dp"

    cases = ["pre2 == 0.0E0_"+prec+" .and. pre1 == 1.0E0_"+prec]
    cases.append("pre2 == 0.0E0_"+prec+" .and. pre1 /= 1.0E0_"+prec)
    cases.append("pre2 == 1.0E0_"+prec+" .and. pre1 == 1.0E0_"+prec)
    cases.append("pre2 == 1.0E0_"+prec+" .and. pre1 /= 1.0E0_"+prec)
    cases.append("pre2 /= 1.0E0_"+prec+" .and. pre1 == 1.0E0_"+prec)
    cases.append("pre2 /= 1.0E0_"+prec+" .and. pre1 /= 1.0E0_"+prec)
  else:
    cases = [".true."]


  #defining acc directives and write acc parallel
  if acc:
     oaccloop_end = "!$acc end loop\n"
     oaccparallel_end = "!$acc end parallel\n"
     oaccparallel_async = " async(async_id1)"
     oaccparallel_wait  = "!$acc wait(async_id2)"
     oaccparallel_init = "!$acc parallel present(array_in,array_out)&\n"
     oaccparallel_init += "    !$acc& firstprivate(pre1,pre2,d"
     for j in range(modes):
        oaccparallel_init += abc[j]+",d"
     oaccparallel_init = oaccparallel_init[:-2] + ") private("
     for j in range(modes):
        oaccparallel_init += abc[j]+","
     oaccparallel_init = oaccparallel_init[:-1] + ")"

     if modes == 3:
       oaccloop_gang = "!$acc loop gang\n"
       oaccloop_worker = "!$acc loop worker\n"
     elif modes == 2 :
       oaccloop_gang = "!$acc loop gang, worker\n"
     else:
       oaccloop_gang = "!$acc loop gang collapse("+str(modes-2)+")\n"
       oaccloop_worker = "!$acc loop worker\n"

     oaccloop_vector = "!$acc loop vector\n"

     f.write("\n")
     f.write("    if(wait_arg)then\n\n")
     f.write("      "+oaccparallel_wait+oaccparallel_async+"\n\n")
     f.write("    endif\n")
     #f.write("    "+oaccparallel_wait+oaccparallel_async+"  if(wait_arg)\n\n")
     f.write("    "+oaccparallel_init+oaccparallel_async)
     f.write("\n")


  #write sfuff for all the cases above
  for cas in range(len(cases)):
    if(cas==0):
      f.write("\n    precase: if("+cases[cas]+")then\n")
    else:
      f.write("\n    elseif("+cases[cas]+")then\n")
  
    #WRITE OMP PARALLEL STATEMENT HERE
    omppar ="      !$OMP PARALLEL DEFAULT(NONE),PRIVATE("
    for j in range(modes):
      omppar += abc[j]+",b"+abc[j]+","
    if (ad == "f2t" or ad == "t2f"): 
      for j in range(modes):
        omppar += "b"+abc[j]+"f,"
    omppar = omppar[0:-1] + ")&\n      !$OMP SHARED(bcntr,pre1,pre2,array_in,array_out, &\n      !$OMP "
    for j in range(modes):
      omppar += "d"+abc[j] +",d"+abc[j]+"2,mod"+abc[j]+","
    if ((ad == "f2t") or (ad == "t2f")):
      omppar += "&\n      !$OMP "
      for j in range(modes):
        omppar += "f"+abc[j] +","
    omppar = omppar[0:-1]+")\n"
 
    if(not debug_loops and not acc): 
      f.write(omppar)
  
    #get the batched space
    casecounter = 1
    if(not debug_loops and not acc):
      fullrange = modes+1
    else:
      fullrange = 1

    for i in range(fullrange):

      #FIND THE RESTRICTED INDICES IN THE OLD ORDERING
      all_thingys = combinations(idxarr,i)

      for oldr in all_thingys:
  
        label1 = "r"+str(casecounter)
        label1 = ""
        casecounter += 1
        #get conditions for the reordering
        conditions = []
        oldu  = []
        newu  = []
        newr  = []
        for j in range(modes):
          oldu.append(idxarr[j])
          newu.append(perm[j])
          newr.append(perm[j])
          conditions.append("d"+abc[j]+"2>0")
  
        #modify the conditions accordingly
        for j in  range(len(oldr)):
          conditions[oldr[j]] = "mod"+abc[oldr[j]]
          for k in range(len(oldu)):
            if oldu[k] == oldr[j]:
              del oldu[k]
              break
          for k in range(len(newu)):
            if newu[k] == oldr[j]:
              del newu[k]
              break
        for j in  range(len(oldu)):
          for k in range(len(newr)):
            if newr[k] == oldu[j]:
              del newr[k]
              break
        
        #Build if-statement for current restrictions, or write 
        if(not debug_loops and not acc):
          conditionalstatement = "      "+label1+"if("
          for j in range(len(conditions)):
            conditionalstatement += conditions[j]
            if j != len(conditions)-1:
              conditionalstatement+=".and."
          conditionalstatement += ")then\n"
          f.write(conditionalstatement)
    
        #WRITE OMP DO STUFF HERE
        if(not debug_loops and not acc and modes-len(oldr)):
          if( nocollapse ):
             ncoll = min(modes-len(oldr),1)
          else:
             ncoll = modes-len(oldr)

          ompdo = "        !$OMP DO" 
          if ( ncoll > 1 ):
            ompdo += " COLLAPSE("+str(ncoll)+")"
          ompdo += "\n"
          f.write(ompdo)
 
        #ORDER THE LOOPS, this depends on the architecture and may be modified
        #THESE CONDITIONS ARE SET UP FOR INTEL, please adapt whenever a different compiler/architecture is used
        useold = False
        for j in range(len(oldu)):
          if(idxarr[oldu[j]]<perm[newu[j]]):
            useold = True
            break
          elif(idxarr[oldu[j]]>perm[newu[j]]):
            useold =False
            break
          
        #BUILD THE ORDER OF THE LOOPS ACCORDING TO THE PREVIOUS CONDITIONS
        if useold:
          outer = oldu
          inner = idxarr
        else:
          outer = newu
          inner = perm

        inneri = [j for j in reversed(inner)]
        outeri = [j for j in reversed(outer)]
  
        #WRITING THE OUTER FOR LOOPS HERE:
        if(not debug_loops and not acc):
          offsetstr="        "
          for j in  range(len(outeri)):
            f.write(offsetstr+"do b"+abc[outeri[j]]+"=1,d"+abc[outeri[j]]+"2,bs\n")
            offsetstr += "  "
          offsetstr = offsetstr[0:-2]
        else:
          offsetstr="    "

        f.write("\n")

        #wrtiting the first element counters for t2f and f2t versions
        if(not debug_loops):
          offsetstr2 = offsetstr + "  "
          for j in  range(len(outeri)):
            if ((ad == "t2f") or (ad == "f2t")):
              f.write(offsetstr2+"b"+abc[outeri[j]]+"f = f"+abc[outeri[j]]+" + b"+abc[outeri[j]]+"\n")          
        else:
          offsetstr2 = offsetstr + "  "        
       
        #WRITING THE INNER FOR LOOPS HERE: 
        if not debug_loops and not acc :
          for j in range(modes):
            if(inneri[j] in newr):
              f.write(offsetstr2+"do "+abc[inneri[j]]+"=d"+abc[inneri[j]]+"2+1,d"+abc[inneri[j]]+"\n")
            elif(inneri[j] in newu):
              f.write(offsetstr2+"do "+abc[inneri[j]]+"=0,bcntr\n")
            offsetstr2 += "  "
        else:
          #FOR OPENACC WRITE WORK SHARING DIRECTIVES
          if acc:
            f.write(offsetstr2+oaccloop_gang)
            for j in range(modes):
              if (modes == 4 and j == 2):
                f.write(offsetstr2+oaccloop_worker)
              if (modes == 3 and j == 1):
                f.write(offsetstr2+oaccloop_worker)
              if (modes == 4 and j == 3):
                f.write(offsetstr2+oaccloop_vector)
              if (modes == 3 and j == 2):
                f.write(offsetstr2+oaccloop_vector)
              if (modes == 2 and j == 1):
                f.write(offsetstr2+oaccloop_vector)
              f.write(offsetstr2+"do "+abc[inneri[j]]+"=1,d"+abc[inneri[j]]+"\n")
              offsetstr2 += "  "
          else:
             for j in range(modes):
               f.write(offsetstr2+"do "+abc[inneri[j]]+"=1,d"+abc[inneri[j]]+"\n")
               offsetstr2 += "  "
  
        offsetstr2 = offsetstr2[0:-2]
        
        #CENTRAL COPYING AND ADDITION STRING
        newidx = ""
        oldidx = ""
        if(not debug_loops and not acc):
          for j in range(modes):
            if(perm[j] in newu):
              if ad == "t2f":
                newidx += "b"+abc[perm[j]]+"f+"+abc[perm[j]]+","
              else:
                newidx += "b"+abc[perm[j]]+"+"+abc[perm[j]]+","
            elif(perm[j] in newr):
              if ad == "t2f":
                newidx += "f"+abc[perm[j]]+"+"+abc[perm[j]]+","
              else:
                newidx += abc[perm[j]]+","
            else:
              print "INVALID STUFF HAPPENING HERE"
            if(idxarr[j] in newu):
              if ad == "f2t":
                oldidx += "b"+abc[idxarr[j]]+"f+"+abc[idxarr[j]]+","
              else:
                oldidx += "b"+abc[idxarr[j]]+"+"+abc[idxarr[j]]+","
            elif(idxarr[j] in newr):
              if ad == "f2t":
                oldidx += "f"+abc[idxarr[j]]+"+"+abc[idxarr[j]]+","
              else:
                oldidx += abc[idxarr[j]]+","
            else:
              print "INVALID STUFF HAPPENING HERE"
        else:
          for j in range(modes):
            if ad == "t2f":
              newidx += "f"+abc[perm[j]]+"+"+abc[perm[j]]+","
              oldidx += abc[idxarr[j]]+","
            elif ad == "f2t":
              newidx += abc[perm[j]]+","
              oldidx += "f"+abc[idxarr[j]]+"+"+abc[idxarr[j]]+","
            else:
              newidx += abc[perm[j]]+","
              oldidx += abc[idxarr[j]]+","

        newidx = newidx[0:-1]
        oldidx = oldidx[0:-1]
  
        cpstr=offsetstr2+"  array_out("
        if(not debug_loops):
          if cas == 0:
            cpstr += newidx+")=array_in("+oldidx+")\n"
          elif cas==1:
            cpstr += newidx+")=pre1*array_in("+oldidx+")\n"
          elif cas==2:
            cpstr += newidx+")=array_out("+newidx+")&\n                                                  &+array_in("+oldidx+")\n"
          elif cas==3:
            cpstr += newidx+")=array_out("+newidx+")&\n                                                  &+pre1*array_in("+oldidx+")\n"
          elif cas==4:
            cpstr += newidx+")=pre2*array_out("+newidx+")&\n                                                  &+array_in("+oldidx+")\n"
          elif cas==5:
            cpstr += newidx+")=pre2*array_out("+newidx+")&\n                                                  &+pre1*array_in("+oldidx+")\n"
        else:
          cpstr += newidx+")=pre2*array_out("+newidx+")&\n                                                  &+pre1*array_in("+oldidx+")\n"
  
        f.write(cpstr)
        
 
        for j in  range(modes-1,-1,-1):
          f.write(offsetstr2+"enddo\n")
          if acc and (j==modes-1 or j==modes-2 or j==0):
             f.write(offsetstr2+oaccloop_end)
          offsetstr2 = offsetstr2[0:-2]
        f.write("\n")
        offsetstr2 = offsetstr2[0:-2]
  
        #WRITING THE OUTER ENDOFOR HERE:
        if(not debug_loops and not acc):
          for j in  range(modes-len(oldr)-1,-1,-1):
            f.write(offsetstr+"enddo\n")
            offsetstr = offsetstr[0:-2]

          if(modes-len(oldr)> 0):
            ompdo ="        !$OMP END DO NOWAIT\n"
            f.write(ompdo)

        if(not debug_loops and not acc):
         
          conditionalstatement="      endif "+label1+"\n"
          f.write(conditionalstatement)
        
      if(not debug_loops and not acc):
        if(i==modes-1):
          ompdo = "      !$OMP END PARALLEL\n"
          f.write(ompdo)
  
    if(cas==len(cases)-1):
      f.write("    endif precase\n")
      if acc:
         f.write("    "+oaccparallel_end+"\n")
   

#WRITE THE HEADER AND GET THE SUBROUTINE NAME
def write_subroutine_header(f,idxarr,perm,now,modes,ad,deb,acc):
  reordstr1 = ""
  reordstr2 = ""
  reordstr3 = ""

  if ad == "sp":
     prec = "tensor_sp"
  else:
     prec = "tensor_dp"


  for i in range(modes):
    reordstr1 += str(perm[i]+1)
    intermed = "dims("+str(idxarr[i]+1)+"),"
    if ad == "f2t":
      intermed = "f"+intermed
    reordstr2 += intermed
    intermed = "dims("+str(perm[i]+1)+"),"
    if ad == "t2f":
      intermed = "f"+intermed
    reordstr3 += intermed

  reordstr2 = reordstr2[0:-1]
  reordstr3 = reordstr3[0:-1]

  #BUILD SUBROUTINE NAME
  suppl = "_"+ad
  if acc:
    suppl += "_acc"

  sname =  "manual_"+reordstr1+"_reordering"+suppl

  #GET THE SUBROUTINE HEADER
  subheaderstr= "  !\> \\brief reorder a "+str(modes)+" diensional array  to get the indices\n"
  subheaderstr+= "  !   in the order "+reordstr1+" , this is a quite expensive reordering\n"
  subheaderstr+= "  !   and thus requires additional attention \n"
  subheaderstr+= "  !\> \\author Patrick Ettenhuber"
  if acc :
     subheaderstr+= " and Janus Juul Eriksen"
  subheaderstr+= "\n"
  subheaderstr+= "  !\> \date "+str(now.month)+", "+str(now.year)+"\n"
  subheaderstr+= "  subroutine "+sname+"(dims,"
  if ((ad == "f2t") or (ad == "t2f")):
    subheaderstr += "fdims,fels,"
  subheaderstr+= "pre1,array_in,pre2,array_out"
  if acc :
    subheaderstr += ",async_id1,async_id2,wait_arg"
  subheaderstr+=")\n"
  subheaderstr+= "    implicit none\n"
  if not acc:
    subheaderstr+= "    !> input for the block size in tiled reordering\n"
    subheaderstr+= "    integer, parameter :: bs=BS_"+str(modes)+"D\n"
  subheaderstr+= "    !>  the dimensions of the different modes in the original array\n"
  subheaderstr+= "    integer, intent(in) :: dims("+str(modes)+")"
  if ((ad == "f2t") or (ad == "t2f")):
    subheaderstr += ",fdims("+str(modes)+"),fels("+str(modes)+")"
  subheaderstr+= "\n"
  subheaderstr+= "    !> as this routine can be used for adding and scaling these are the prefactors\n"
  subheaderstr+= "    real("+prec+"),intent(in) :: pre1,pre2\n"
  subheaderstr+= "    !> array to be reordered\n"
  subheaderstr+= "    real("+prec+"),intent(in) :: array_in("+reordstr2+")\n"
  subheaderstr+= "    !> reordered array\n"
  subheaderstr+= "    real("+prec+"),intent(inout) :: array_out("+reordstr3+")\n"
  if acc :
    subheaderstr+= "    integer(acc_handle_kind),intent(in) :: async_id1\n"
    subheaderstr+= "    integer(acc_handle_kind),intent(in) :: async_id2\n"
    subheaderstr+= "    logical,intent(in) :: wait_arg\n"


  subheaderstr+= "    integer :: "
  if not acc:
    subheaderstr+= "bcntr,"
  for i in range(modes):
    subheaderstr+= abc[i]+",d"+abc[i]+","
    if not acc:
      subheaderstr+= "b"+abc[i]+",d"+abc[i]+"2,"
  subheaderstr = subheaderstr[0:-1]

  subheaderstr += "\n"
  if (ad == "f2t" or ad == "t2f"):
    subheaderstr+= "    integer :: "
    for i in range(modes):
      subheaderstr+= "b"+abc[i]+"f,f"+abc[i]+","
    subheaderstr = subheaderstr[0:-1] + "\n"

  if not acc:
    subheaderstr+= "    logical :: "
    for i in range(modes):
      subheaderstr+= "mod"+abc[i]+","
    subheaderstr = subheaderstr[0:-1] + "\n"

  subheaderstr += "\n"
  for i in range(modes):
    subheaderstr+= "    d"+abc[i]+"=dims("+str(i+1)+")\n"
  subheaderstr+= "\n"

  if (ad == "f2t" or ad == "t2f"):
    for i in range(modes):
      subheaderstr+= "    f"+abc[i]+"=fels("+str(i+1)+")-1\n"
  subheaderstr+= "\n"

  if not acc:
     for i in range(modes):
       subheaderstr+= "    d"+abc[i]+"2=(d"+abc[i]+"/bs)*bs\n"
     subheaderstr+= "\n"
     for i in range(modes):
       subheaderstr+= "    mod"+abc[i]+"=(mod(d"+abc[i]+",bs)>0)\n"
     subheaderstr+= "\n    bcntr=bs-1\n"

  if deb :
    if ad == "t2f":
      subheaderstr += "\n    if(pre2==0.0E0_"+prec+")array_out(&\n"
      for i in range(modes):
        subheaderstr += "                 &fels("+str(perm[i]+1)+"):fels("+str(perm[i]+1)+")+dims("+str(perm[i]+1)+")-1,&\n"
      subheaderstr = subheaderstr[0:-3]+") = 0.0E0_"+prec+"\n\n"
    else:
      subheaderstr += "\n    if(pre2==0.0E0_"+prec+")array_out = 0.0E0_"+prec+"\n\n"
  
  #WRITE THE HEADER TO FILE
  f.write(subheaderstr)

  #RETURN THE NAME OF THE SUBROUTINE
  return sname


def write_simple_module_header(f,idim,idx,now,args,acc):
   f.write("!\> \\brief this autogenerated module is inteded to contain high performance reorderings for\n!mutlidimensional arrays to tiles in a different distribution.\n!\> \\author Patrick Ettenhuber\n!\> \\date March 2013, file produced: "+str(now.month)+", "+str(now.year)+"\n")
   #WRITE THE VARIABLES WITH WHICH THE FILE WAS PRODUCED HERE 
   #--> FOR LATER READOUT AND SEE IF IT IS NECESSARY TO PRODUCE A NEW INSTANCE OF IT
   for i in range(len(args)):
     f.write("!ARG"+str(i)+": "+str(args[i])+"\n")

   f.write("!END VARS\n\n")

   namestub = get_namestub_from_path(f.name)

   f.write("module "+namestub+"_module\n")
   f.write("  use tensor_parameters_and_counters\n")
   if(acc):
     f.write("  use openacc\n")
   f.write("\n")
   f.write("  contains\n")

def get_namestub_from_path(path):
   return path[path.rfind("/")+1:].replace(".F90","") 

def write_simple_module_end_and_close(f):
   f.write("end module "+get_namestub_from_path(f.name)+"_module\n")
   f.close()


def write_main_header(f,now,names,args,tensordir,minr,maxr,interface_types):
   f.write("!\> \\brief this autogenerated module is inteded to contain high performance reorderings for\n!mutlidimensional arrays.\n!\> \\author Patrick Ettenhuber & Janus Juul Eriksen\n!\> \\date November 2012, file produced: "+str(now.month)+", "+str(now.year)+"\n")
   #WRITE THE VARIABLES WITH WHICH THE FILE WAS PRODUCED HERE 
   #--> FOR LATER READOUT AND SEE IF IT IS NECESSARY TO PRODUCE A NEW INSTANCE OF IT
   for i in range(len(args)):
     f.write("!ARG"+str(i)+": "+str(args[i])+"\n")

   f.write("!END VARS\n\n")

   f.write("module reorder_frontend_module\n")
   f.write("  use precision\n")
   f.write("  use tensor_parameters_and_counters\n")
   f.write("  use get_idx_mod\n")
   f.write("  use memory_handling\n")
   
   for name in names:
      f.write("  use "+name+"_module\n")

   f.write("  use LSTIMING\n\n")

   #write interfaces
   #SPECIFY THE ORDER OF REODERINGS
   for idx in range(maxr-minr+1):
     modes = idx + minr
     # GENERATE ORIGINAL ORDER AND STARTING POINT FOR NEW ORDER
     idxarr = [ i for i in range(modes) ]
     pdxarr = permutations(idxarr)
     for perm in pdxarr:
        pstr = ""
        for i in range(len(perm)):
           pstr += str(perm[i]+1)
        iname = "manual_"+pstr+"_reordering"
        interface = "  interface "+iname+"\n"
        interface += "     module procedure "
        write_interface = False
        for typ in interface_types:
           for sname in typ:
              if iname in sname:
                 if write_interface:
                    interface += "       &"
                 else:
                    write_interface = True
                 interface += sname +",&\n"
        interface = interface[:-3] + "\n  end interface "+iname+"\n\n"

        if write_interface:
           f.write(interface)

   #f.write("  contains\n")
   #Write the subroutines called by the user
   basic = open(tensordir+"/autogen/reorder_header.F90",'r')
   for line in basic:
     f.write(line)
   basic.close


def write_testing_framework(installdir,minr,maxr):
  f=open(installdir+"reorder_tester.F90",'w')
  header="\
module reorder_tester_module\n\
  use reorder_frontend_module\n\
  contains\n\
  subroutine test_array_reorderings(LUPRI)\n    implicit none\n\n\
    real(tensor_dp),pointer :: in1(:),sto(:)\n\
    real(tensor_dp),pointer :: res(:),til(:)\n\
    real(tensor_dp) :: ref(6),ref1s,ref2s,ref1,ref2\n\
    integer :: tile_idx,LUPRI\n\
    logical :: master,rigorous\n\
    integer :: "

  for i in range(maxr):
    header += abc[i]+",n"+abc[i]+","
  header = header[0:-1]+"\n    integer :: p1,p2\n\
    real(tensor_dp) :: pr1,pr2,begc1,begw1,endc1,endw1,begc2,begw2,endc2,endw2\n\
    character(len=7) :: teststatus\n\
    logical :: "
  for i in range(minr,maxr+1):
    header += "test"+str(i)+"d_normal,"
  header = header[0:-1] + "\n    logical :: "
  for i in range(minr,maxr+1):
    header += "test"+str(i)+"d_tiled,"
  header = header[0:-1] +"\n    master = .true.\n"
  for i in range(minr,maxr+1):
    header += "    test"+str(i)+"d_normal=.true.\n"
  for i in range(minr,maxr+1):
    #if i== 4:
    header += "    test"+str(i)+"d_tiled=.true.\n"
    #else:
    #  header += "    test"+str(i)+"d_tiled=.false.\n"

  header += "    rigorous=.true.\n"
 
  f.write(header)

  for mode in range(minr,maxr+1):
    #WRITING THE TESTING FRAMEWORK HERE
    words ="\
    if(test"+str(mode)+"d_normal)then\n\
      write (LUPRI,*)\"\"\n\
      write (LUPRI,*)\"TESTING "+str(mode)+"D REORDERINGS\"\n\
      write (LUPRI,*)\"**********************\"\n\
      write (LUPRI,*)\"\"\n"
    for i in range(mode):
      words += "      n"+abc[i]+"="+str(int(((8000.0*1000.0)/(8.0*2.0))**(1.0/(mode*1.0)))+1+randrange(10))+"\n"
    words +="\
      call mem_alloc(in1,"
    for i in range(mode):
      words += "n" + abc[i] + "*"
    words = words[0:-1]+")\n      call mem_alloc(res,"
    for i in range(mode):
      words += "n" + abc[i] + "*"
    words = words[0:-1]+")\n      call mem_alloc(sto,"
    for i in range(mode):
      words += "n" + abc[i] +"*"
    words = words[0:-1]+")\n      call random_number(in1)\n      call random_number(sto)\n\n" 

    words +="\
      do p1=1,2\n\
        do p2=0,2\n\
          pr1 = float(p1)\n\
          pr2 = float(p2)\n\
          if (p1==2) call random_number(pr1)\n\
          if (p2==2) call random_number(pr2)\n\n\
          write (LUPRI,'(A3,f4.1,A3,f4.1,A2)')\"B= \",pr1,\"*A+\",pr2,\"*B\"\n\n"
    f.write(words)
    pc = 0
    maxperms = factorial(mode)
    for perm in permutations(list(xrange(mode))):
      testcase="\
          call LSTIMER('START',begc1,begw1,LUPRI,.false.)\n\
          teststatus=\"SUCCESS\"\n\
          !res = sto\n\
          call dcopy("
      for i in range(mode):
        testcase += "n"+abc[i]+"*"
      testcase = testcase[0:-1] + ",sto,1,res,1)\n"
      testcase+="\
          call LSTIMER('START',begc2,begw2,LUPRI,.false.)\n\
          call array_reorder_"+str(mode)+"d(pr1,in1,"
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase +="["
      for i in range(mode):
        testcase += str(perm[i]+1)+","
      testcase = testcase[0:-1] +"],pr2,res)\n\
          call LSTIMER('START',endc2,endw2,LUPRI,.false.)\n\
          if(rigorous)then\n"
      ofstr="            "
      for i in range(mode):
        testcase += ofstr + "do "+abc[mode-1-i]+"=1,n"+abc[mode-1-i]+"\n"
        ofstr += "  "
      testcase += ofstr+"if(abs(pr1*in1("
      for i in range(mode):
        if i== 0:
          testcase += abc[i]
        else:
          testcase += "("+abc[i]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[j]
        testcase += "+"
      testcase =testcase[0:-1]+")+pr2*sto("
      for i in range(mode):
        if i== 0:
          testcase += abc[perm[i]]
        else:
          testcase += "("+abc[perm[i]]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[perm[j]]
        testcase += "+"
      testcase = testcase[0:-1]+")&\n"+ofstr+"&-res("
      for i in range(mode):
        if i== 0:
          testcase += abc[perm[i]]
        else:
          testcase += "("+abc[perm[i]]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[perm[j]]
        testcase += "+"
      testcase = testcase[0:-1]+"))>1.0E-11_tensor_dp)teststatus=\"FAILED \"\n"
      ofstr = ofstr[0:-2]
      for i in range(mode):
        testcase += ofstr+"enddo\n"
        ofstr = ofstr[0:-2]
      testcase+="          endif\n"
      testcase +="\
          call LSTIMER('START',endc1,endw1,LUPRI,.false.)\n\
          write (LUPRI,&\n\
          &'(I1,I1,\"-"
      for i in range(mode):
        testcase += str(perm[i]+1)
      testcase +=" C-T1: \",f9.4,\" W-T1: \",f9.4,\" C-T2: \",f9.4,\" W-T2: \",f9.4,\" STATUS=\",A7)')&\n          &p1,p2,endc1-begc1,endw1-begw1,endc2-begc2,endw2-begw2,teststatus\n\n"
    
      f.write(testcase)
      if(pc==maxperms-1):
        f.write("          write (LUPRI,*)\"\"\n\n")
      pc += 1
    f.write("        enddo\n      enddo\n")
    f.write("      call mem_dealloc(in1)\n      call mem_dealloc(res)\n      call mem_dealloc(sto)\n")
    f.write("    endif\n")


  #SPECIAL TESTS FOR TILED REORDERINGS
  for mode in range(minr,maxr+1):
    #WRITING THE TESTING FRAMEWORK HERE

    words ="\
    if(test"+str(mode)+"d_tiled)then\n\
      write (LUPRI,*)\"\"\n\
      write (LUPRI,*)\"TESTING "+str(mode)+"D TILED REORDERINGS\"\n\
      write (LUPRI,*)\"**********************\"\n\
      write (LUPRI,*)\"\"\n"


    for i in range(mode):
      sze=int(((8000.0*1000.0)/(8.0*2.0))**(1.0/(mode*1.0)))+1+randrange(10)
      #THIS IS INTRODUCED TO ENSURE WE HAVE EXACTLY 3 TILES
      if(i==mode-1):
        sze = 2*sze+1
      words += "      n"+abc[i]+"="+str(sze)+"\n"

    words +="      call mem_alloc(til,"
    for i in range(mode):
      if(i==mode-1):
        words+="(n"+abc[i]+"/2))\n"
      else:
        words+="n"+abc[i]+"*"
    words +="\
      call mem_alloc(in1,"
    for i in range(mode):
      words += "n" + abc[i] + "*"
    words = words[0:-1]+")\n      call mem_alloc(res,"
    for i in range(mode):
      words += "n" + abc[i] + "*"
    words = words[0:-1]+")\n      call mem_alloc(sto,"
    for i in range(mode):
      words += "n" + abc[i] +"*"
    words = words[0:-1]+")\n      call random_number(in1)\n      call random_number(sto)\n\n" 

    words +="\
      do p1=1,2\n\
        do p2=0,2\n\
          pr1 = float(p1)\n\
          pr2 = float(p2)\n\
          if (p1==2) call random_number(pr1)\n\
          if (p2==2) call random_number(pr2)\n\n\
          write (LUPRI,'(A3,f4.1,A3,f4.1,A2)')\"B= \",pr1,\"*A+\",pr2,\"*B\"\n\n"
    f.write(words)
    pc = 0
    maxperms = factorial(mode)
    for perm in permutations(list(xrange(mode))):

      #TEST THE FROM FORT REORDERINGS
      testcase="\
          call LSTIMER('START',begc1,begw1,LUPRI,.false.)\n\
          teststatus=\"SUCCESS\"\n\
          !res = sto\n\
          call dcopy("
      for i in range(mode):
        testcase += "n"+abc[i]+"*"
      testcase = testcase[0:-1] + ",sto,1,res,1)\n"
      testcase+="\
          call LSTIMER('START',begc2,begw2,LUPRI,.false.)\n\
          do tile_idx=1,3\n            call tile_from_fort(1.0E0_tensor_dp,in1,["
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase = testcase[0:-1] +"],"+str(mode)+",0.0E0_tensor_dp,til,tile_idx,["
      for i in range(mode):
        if(i==mode-1):
          testcase += "n"+abc[i]+"/2],["
        else:
          testcase += "n"+abc[i]+","
      for i in range(mode):
        testcase += str(i+1)+","
      testcase = testcase[0:-1] +"])\n"

      testcase +="            call tile_in_fort(pr1,til,tile_idx,["
      for i in range(mode):
        if(i==mode-1):
          testcase += "n"+abc[i]+"/2],pr2,res,["
        else:
          testcase += "n"+abc[i]+","
      for i in range(mode):
        testcase += "n"+abc[perm[i]]+","
      testcase = testcase [0:-1] + "],"+str(mode)+",["
      for i in range(mode):
        testcase += str(perm[i]+1)+","
      testcase = testcase [0:-1] + "])\n          enddo\n"

      testcase +="\
          call LSTIMER('START',endc2,endw2,LUPRI,.false.)\n\
          if(rigorous)then\n"
      ofstr="            "
      for i in range(mode):
        testcase += ofstr + "do "+abc[mode-1-i]+"=1,n"+abc[mode-1-i]+"\n"
        ofstr += "  "
      testcase += ofstr+"if(abs(pr1*in1("
      for i in range(mode):
        if i== 0:
          testcase += abc[i]
        else:
          testcase += "("+abc[i]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[j]
        testcase += "+"
      testcase =testcase[0:-1]+")+pr2*sto("
      for i in range(mode):
        if i== 0:
          testcase += abc[perm[i]]
        else:
          testcase += "("+abc[perm[i]]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[perm[j]]
        testcase += "+"
      testcase = testcase[0:-1]+")&\n"+ofstr+"&-res("
      for i in range(mode):
        if i== 0:
          testcase += abc[perm[i]]
        else:
          testcase += "("+abc[perm[i]]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[perm[j]]
        testcase += "+"
      testcase = testcase[0:-1]+"))>1.0E-11_tensor_dp)teststatus=\"FAILED \"\n"
      ofstr = ofstr[0:-2]
      for i in range(mode):
        testcase += ofstr+"enddo\n"
        ofstr = ofstr[0:-2]
      testcase+="          endif\n"
      testcase +="\
          call LSTIMER('START',endc1,endw1,LUPRI,.false.)\n\
          write (LUPRI,&\n\
          &'(I1,I1,\"-"
      for i in range(mode):
        testcase += str(perm[i]+1)
      testcase +=" C-T1: \",f9.4,\" W-T1: \",f9.4,\" C-T2: \",f9.4,\" W-T2: \",f9.4,\" STATUS=\",A7)')&\n          &p1,p2,endc1-begc1,endw1-begw1,endc2-begc2,endw2-begw2,teststatus\n\n"
    
      #TEST THE IN_FORT REORDERINGS
      testcase+="\
          call LSTIMER('START',begc1,begw1,LUPRI,.false.)\n\
          teststatus=\"SUCCESS\"\n\
          !res = sto\n\
          call dcopy("
      for i in range(mode):
        testcase += "n"+abc[i]+"*"
      testcase = testcase[0:-1] + ",sto,1,res,1)\n"
      testcase+="\
          call LSTIMER('START',begc2,begw2,LUPRI,.false.)\n\
          do tile_idx=1,3\n            call tile_from_fort(1.0E0_tensor_dp,sto,["
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase = testcase[0:-1] +"],"+str(mode)+",0.0E0_tensor_dp,til,tile_idx,["
      for i in range(mode):
        if(perm[i]==mode-1):
          testcase += "n"+abc[perm[i]]+"/2,"
        else:
          testcase += "n"+abc[perm[i]]+","
      testcase = testcase[0:-1] + "],["
      for i in range(mode):
        testcase += str(perm[i]+1)+","
      testcase = testcase[0:-1] + "])\n            call tile_from_fort(pr1,in1,["
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase = testcase[0:-1] +"],"+str(mode)+",pr2,til,tile_idx,["
      for i in range(mode):
        if(perm[i]==mode-1):
          testcase += "n"+abc[perm[i]]+"/2,"
        else:
          testcase += "n"+abc[perm[i]]+","
      testcase = testcase[0:-1] + "],["
      for i in range(mode):
        testcase += str(perm[i]+1)+","
      testcase = testcase[0:-1] + "])\n            call tile_in_fort(1.0E0_tensor_dp,til,tile_idx,["
      for i in range(mode):
        if(perm[i]==mode-1):
          testcase += "n"+abc[perm[i]]+"/2,"
        else:
          testcase += "n"+abc[perm[i]]+","
      testcase = testcase[0:-1] + "],0.0E0_tensor_dp,res,["
      for i in range(mode):
          testcase += "n"+abc[perm[i]]+","
      testcase = testcase[0:-1] + "],"+str(mode)+",["
      for i in range(mode):
          testcase += str(i+1)+","
      testcase = testcase[0:-1] + "])\n          enddo\n\
          call LSTIMER('START',endc2,endw2,LUPRI,.false.)\n\
          if(rigorous)then\n"
      ofstr="            "
      for i in range(mode):
        testcase += ofstr + "do "+abc[mode-1-i]+"=1,n"+abc[mode-1-i]+"\n"
        ofstr += "  "
      testcase += ofstr+"if(abs(pr1*in1("
      for i in range(mode):
        if i== 0:
          testcase += abc[i]
        else:
          testcase += "("+abc[i]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[j]
        testcase += "+"
      testcase =testcase[0:-1]+")+pr2*sto("
      for i in range(mode):
        if i== 0:
          testcase += abc[i]
        else:
          testcase += "("+abc[i]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[j]
        testcase += "+"
      testcase = testcase[0:-1]+")&\n"+ofstr+"&-res("
      for i in range(mode):
        if i== 0:
          testcase += abc[perm[i]]
        else:
          testcase += "("+abc[perm[i]]+"-1)"
          for j in range(i):
            testcase += "*n"+abc[perm[j]]
        testcase += "+"
      testcase = testcase[0:-1]+"))>1.0E-11_tensor_dp)teststatus=\"FAILED \"\n"
      ofstr = ofstr[0:-2]
      for i in range(mode):
        testcase += ofstr+"enddo\n"
        ofstr = ofstr[0:-2]
      testcase+="          endif\n"
      testcase +="\
          call LSTIMER('START',endc1,endw1,LUPRI,.false.)\n\
          write (LUPRI,&\n\
          &'(I1,I1,\"-"
      for i in range(mode):
        testcase += str(perm[i]+1)
      testcase +=" C-T1: \",f9.4,\" W-T1: \",f9.4,\" C-T2: \",f9.4,\" W-T2: \",f9.4,\" STATUS=\",A7)')&\n          &p1,p2,endc1-begc1,endw1-begw1,endc2-begc2,endw2-begw2,teststatus\n\n"



      f.write(testcase)
      if(pc==maxperms-1):
        f.write("          write (LUPRI,*)\"\"\n\n")
      pc += 1
    f.write("        enddo\n      enddo\n")
    f.write("      call mem_dealloc(til)\n      call mem_dealloc(in1)\n      call mem_dealloc(res)\n      call mem_dealloc(sto)\n")
    f.write("    endif\n")



  f.write("  end subroutine test_array_reorderings\nend module reorder_tester_module")
  f.close()

#ITERTOOLS COMBINATIONS AND PERMUTATIONS ARE NOT IMPLEMENTED BELOW PYTHON 2.6, THUS I COPIED THE
#SOURCE CODE OF THE HOMEPAGE DIRECTLY
##################################################################################################
##################################################################################################

#THIS IS A MODIFIED COPY OF http://docs.python.org/2/library/itertools.html#itertools.permutations
##################################################################################################
def permutations(iterable):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    r = n
    if r > n:
        return
    indices = range(n)
    cycles = range(n, n-r, -1)
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

#THIS IS A MODIFIED COPY OF http://docs.python.org/2/library/itertools.html#itertools.permutations
##################################################################################################
def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

##################################################################################################
##################################################################################################
def factorial(d):
   k=1
   for i in range(1,d+1):
     k *= i
   return k

abc = "abcdefghijklmnopqrstuvwxyz"
main()
