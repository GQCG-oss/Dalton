#GENERATING TENSOR PERMUTATIONS
#AUTHOR: PATRICK ETTENHUBER
#EMAIL : pett@chem.au.dk, pettenhuber@gmail.com
#DATE  : JUNE, 2013
import sys,datetime,os,time#,itertools,math
from random import randrange
"""This file is inteded for the automatic generation of a data sorting module in LSDALTON
As data sorting is under constant development and features many lines of code this script was
written for easy modification of the code. Furhtermore, one might note that the structure
of optimal code might depend on the architecture which can be at some point introduced here
as well

Currently the routine does not take any command-line arguments, but is dependet on the LSDALTON
folder structure, i.e. it needs LSDALTON/lsutil and LSDALTON/lsutil/autogen/reord_headers.F90
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
  for i in range(len(sys.argv)):
    if "VAR_LSDEBUG" in sys.argv[i]:
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
    
  
# print args
  #GET THE FOLDER TO STORE THE manual_reorderings.F90
  cwd = os.getcwd()
  lsutildir = args[0]
  if ("/LSDALTON/lsutil" in lsutildir ):

    lsutildir = lsutildir[0:lsutildir.find("/LSDALTON/lsutil")]+"/LSDALTON/lsutil/"

  else:

    for paths,dirs,files in os.walk(cwd+"/../LSDALTON"):
      if("/LSDALTON/lsutil" in paths):
        lsutildir = paths[0:paths.find("/LSDALTON/lsutil")]+"/LSDALTON/lsutil/"
        break

    if not os.path.exists(lsutildir):
      print "COULD NOT FIND LSDALTON/lsutil, exiting"
      sys.exit()


  installdir  = ""
  if args[3] != "" :
     installdir = args[3] +"/"
  else:
     installdir = lsutildir

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
  # radovan: this should not be done here but taken care of by (c)make
    if(scriptmod>reordmod or writenew):
      print "REORDER GENERATOR IS NEWER THAN REORDERING FILES - GENERATING NEW ONES"
      writenew = True
  if(not os.path.exists(installdir+"reord2d_2_reord.F90")):
    writenew = True
  if (real_sp_write):
    if(not os.path.exists(installdir+"reord2d_2_reord_sp.F90")):
      writenew = True
  if (acc_write):
    if(not os.path.exists(installdir+"reord2d_acc_reord.F90")):
      writenew = True
      if (real_sp_write):
        if(not os.path.exists(installdir+"reord2d_acc_reord_sp.F90")):
          writenew = True
  if(not os.path.exists(installdir+"reord3d_1_reord.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord3d_2_reord.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord3d_3_reord.F90")):
    writenew = True
  if (real_sp_write):
    if(not os.path.exists(installdir+"reord3d_1_reord_sp.F90")):
      writenew = True
    if(not os.path.exists(installdir+"reord3d_2_reord_sp.F90")):
      writenew = True
    if(not os.path.exists(installdir+"reord3d_3_reord_sp.F90")):
      writenew = True
  if (acc_write):
    if(not os.path.exists(installdir+"reord3d_acc_reord.F90")):
      writenew = True
      if (real_sp_write):
        if(not os.path.exists(installdir+"reord3d_acc_reord_sp.F90")):
          writenew = True
  if(not os.path.exists(installdir+"reord4d_1_reord.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_2_reord.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_3_reord.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_4_reord.F90")):
    writenew = True
  if (real_sp_write):
    if(not os.path.exists(installdir+"reord4d_1_reord_sp.F90")):
      writenew = True
    if(not os.path.exists(installdir+"reord4d_2_reord_sp.F90")):
      writenew = True
    if(not os.path.exists(installdir+"reord4d_3_reord_sp.F90")):
      writenew = True
    if(not os.path.exists(installdir+"reord4d_4_reord_sp.F90")):
      writenew = True
  if (acc_write):
    if(not os.path.exists(installdir+"reord4d_acc_reord.F90")):
      writenew = True
      if (real_sp_write):
        if(not os.path.exists(installdir+"reord4d_acc_reord_sp.F90")):
          writenew = True
  if(not os.path.exists(installdir+"reord4d_1_utils_f2t.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_2_utils_f2t.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_3_utils_f2t.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_4_utils_f2t.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_1_utils_t2f.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_2_utils_t2f.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_3_utils_t2f.F90")):
    writenew = True
  if(not os.path.exists(installdir+"reord4d_4_utils_t2f.F90")):
    writenew = True

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
    produce_files(installdir,lsutildir,args)

##################################################################################################
##################################################################################################
def produce_files(installdir,lsutildir,args):
   maxr = 4
   minr = 2
   hack_only_4d_for_utils = False
   f=open(installdir+"reorder_frontend.F90",'w')
   utils = []
   reord = []
   reord_sp = []
   acc_reord = []
   acc_reord_sp = []
   acc = args[4]
   real_sp = args[5]

   for idx in range(maxr-minr+1):
     #if(idx==maxr-minr):
     sho3 = []
     sho = []
     sho_sp = []
     sho_acc = []
     sho_acc_sp = []
     
     for i in range(idx+minr) :
       sho.append(open(installdir+"reord"+str(idx+minr)+"d_"+str(i+1)+"_reord.F90",'w'))

       if (real_sp):
         sho_sp.append(open(installdir+"reord"+str(idx+minr)+"d_"+str(i+1)+"_reord_sp.F90",'w'))

       sho2 = []
      
       if(idx+minr != 4 and hack_only_4d_for_utils):
         sho2.append(False)
         sho2.append(False)

       else:
         sho2.append(open(installdir+"reord"+str(idx+minr)+"d_"+str(i+1)+"_utils_f2t.F90",'w'))
         sho2.append(open(installdir+"reord"+str(idx+minr)+"d_"+str(i+1)+"_utils_t2f.F90",'w'))
       sho3.append(sho2)

     if (acc):
       sho_acc.append(open(installdir+"reord"+str(idx+minr)+"d_acc_reord.F90",'w'))       
       if (real_sp):
          sho_acc_sp.append(open(installdir+"reord"+str(idx+minr)+"d_acc_reord_sp.F90",'w'))

     reord.append(sho)
     reord_sp.append(sho_sp)
     utils.append(sho3)
     acc_reord.append(sho_acc)
     acc_reord_sp.append(sho_acc_sp)


   #GET COMMAND LINE ARGUMENTS
   debug_loops = args[1]
   nocollapse  = args[2]
    
    
   now = datetime.datetime.now()
   #FULL TO TILE OR TILE TO FULL
   forutils = ["f2t","t2f"]

   # WRITE REORDERINGS FRONTEND FILE
   write_main_header(f,now,args,lsutildir,minr,maxr,hack_only_4d_for_utils)
   write_testing_framework(f,minr,maxr)
   f.write("\nend module reorder_frontend_module")
   f.close()

   #WRITE THE HEADERS OF ALL FILES
   for idx in range(maxr-minr+1):
     for i in range(minr+idx):
       write_simple_module_header(reord[idx][i],idx+minr,i+1,now,args,"r")

       if (real_sp):
         write_simple_module_header(reord_sp[idx][i],idx+minr,i+1,now,args,"sp")

       if (idx+minr != 4 and hack_only_4d_for_utils):
         continue

       for k in range(len(forutils)):
         write_simple_module_header(utils[idx][i][k],idx+minr,i+1,now,args,forutils[k])

     if(acc):
       write_simple_module_header(acc_reord[idx][0],idx+minr,i+1,now,args,"acc")
       if (real_sp):
         write_simple_module_header(acc_reord_sp[idx][0],idx+minr,i+1,now,args,"acc_sp")

   #SPECIFY THE ORDER OF REODERINGS
   for idx in range(maxr-minr+1):
     modes = idx + minr
     # GENERATE ORIGINAL ORDER AND STARTING POINT FOR NEW ORDER
     idxarr = [0]*modes
     for i in range(modes):
       idxarr[i] = i
   
     #GET ALL PERMUTATION
     #all_permutations=itertools.permutations(idxarr)
     all_permutations=permutations(idxarr)
    
     #LOOP OVER PERMUTATIONS AND WRITE SUBROUTINES
     for perm in all_permutations:
       #CHECK IF THE ORDERING IS NECESSARY, ELSE JUST SKIP IT
       doreord = True
       for i in range(len(perm)-1):
         if(perm[i]+1==perm[i+1]):
           doreord = False

       if doreord :
         addition = "normal"
         # WRITE THE CPU SUBROUTINE HEADER AND GET ITS NAME
         sub = write_subroutine_header(reord[idx][perm[0]],idxarr,perm,now,modes,addition,debug_loops)
         #Write the CPU subroutine body
         write_subroutine_body(reord[idx][perm[0]],idxarr,perm,modes,args,addition)
         #END THE CPU SUBROUTINE
         reord[idx][perm[0]].write("  end subroutine "+sub+"\n\n")
         if (real_sp):
           addition = "sp"
           # WRITE THE CPU SUBROUTINE HEADER AND GET ITS NAME
           sub_sp = write_subroutine_header(reord_sp[idx][perm[0]],idxarr,perm,now,modes,addition,debug_loops)
           #Write the CPU subroutine body
           write_subroutine_body(reord_sp[idx][perm[0]],idxarr,perm,modes,args,addition)
           #END THE CPU SUBROUTINE
           reord_sp[idx][perm[0]].write("  end subroutine "+sub_sp+"\n\n")
         # WRITE THE ANALOGOUS GPU STUFF
         # we have six precases
         if (acc):
           addition = "normal"
           for acc_case in range(6):
             sub_acc = write_subroutine_header_acc(acc_reord[idx][0],idxarr,perm,now,modes,acc_case,addition)
             write_subroutine_body_acc(acc_reord[idx][0],idxarr,perm,modes,args,acc_case,addition)
             acc_reord[idx][0].write("  end subroutine "+sub_acc+"\n\n")
           if (real_sp):
             addition = "sp"
             for acc_case in range(6):
               sub_acc_sp = write_subroutine_header_acc(acc_reord_sp[idx][0],idxarr,perm,now,modes,acc_case,addition)
               write_subroutine_body_acc(acc_reord_sp[idx][0],idxarr,perm,modes,args,acc_case,addition)
               acc_reord_sp[idx][0].write("  end subroutine "+sub_acc_sp+"\n\n")

       if(idx+minr!=4 and hack_only_4d_for_utils):
         continue

       for  ad in range(len(forutils)):
         addition = forutils[ad]
         #if(hack_only_4d_for_utils and modes!=4):
         #  break
         # WRITE THE SUBROUTINE HEADER AND GET ITS NAME
         sub =  write_subroutine_header(utils[idx][perm[0]][ad],idxarr,perm,now,modes,addition,debug_loops)
         #Write the subroutine body
         write_subroutine_body(utils[idx][perm[0]][ad],idxarr,perm,modes,args,addition)
         #END THE SUBROUTINE
         #print "  end subroutine "+subroutinename+"\n\n"
         utils[idx][perm[0]][ad].write("  end subroutine "+sub+"\n\n")



   #END the file
   for idx in range(maxr-minr+1):
     for i in range(minr+idx):
       write_simple_module_end_and_close(reord[idx][i],idx+minr,i+1,now,args,'r')
       if (real_sp):
         write_simple_module_end_and_close(reord_sp[idx][i],idx+minr,i+1,now,args,'sp')
       if(idx+minr!=4 and hack_only_4d_for_utils):
         continue

       for ad in range(len(forutils)):
         write_simple_module_end_and_close(utils[idx][i][ad],idx+minr,i+1,now,args,forutils[ad])

     if (acc):
       write_simple_module_end_and_close(acc_reord[idx][0],idx+minr,i+1,now,args,'acc')
       if (real_sp):
         write_simple_module_end_and_close(acc_reord_sp[idx][0],idx+minr,i+1,now,args,'acc_sp')

   #remove empty file
   os.system("rm "+installdir+"reord2d_1_reord.F90")
   if (real_sp):
     os.system("rm "+installdir+"reord2d_1_reord_sp.F90")

def write_subroutine_body(f,idxarr,perm,modes,args,ad):
  debug_loops = args[1]
  nocollapse  = args[2]
  #GENERAL CASE pre1/=1 pre2/=0 or 1
  if(not debug_loops):
    if ad == "sp":
      cases = ["pre2 == 0.0E0_real_sp .and. pre1 == 1.0E0_real_sp"]
      cases.append("pre2 == 0.0E0_real_sp .and. pre1 /= 1.0E0_real_sp")
      cases.append("pre2 == 1.0E0_real_sp .and. pre1 == 1.0E0_real_sp")
      cases.append("pre2 == 1.0E0_real_sp .and. pre1 /= 1.0E0_real_sp")
      cases.append("pre2 /= 1.0E0_real_sp .and. pre1 == 1.0E0_real_sp")
      cases.append("pre2 /= 1.0E0_real_sp .and. pre1 /= 1.0E0_real_sp")
    else:
      cases = ["pre2 == 0.0E0_realk .and. pre1 == 1.0E0_realk"]
      cases.append("pre2 == 0.0E0_realk .and. pre1 /= 1.0E0_realk")
      cases.append("pre2 == 1.0E0_realk .and. pre1 == 1.0E0_realk")
      cases.append("pre2 == 1.0E0_realk .and. pre1 /= 1.0E0_realk")
      cases.append("pre2 /= 1.0E0_realk .and. pre1 == 1.0E0_realk")
      cases.append("pre2 /= 1.0E0_realk .and. pre1 /= 1.0E0_realk")
  else:
    cases = [".true."]

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
    if ((ad != "normal") and (ad != "sp")): 
      for j in range(modes):
        omppar += "b"+abc[j]+"f,"
    omppar = omppar[0:-1] + ")&\n      !$OMP SHARED(bcntr,pre1,pre2,array_in,array_out, &\n      !$OMP "
    for j in range(modes):
      omppar += "d"+abc[j] +",d"+abc[j]+"2,mod"+abc[j]+","
    if ((ad != "normal") and (ad != "sp")):
      omppar += "&\n      !$OMP "
      for j in range(modes):
        omppar += "f"+abc[j] +","
    omppar = omppar[0:-1]+")\n"
 
    if(not debug_loops): 
      f.write("#ifndef VAR_LSESSL\n")
      f.write(omppar)
      f.write("#endif\n")
  
    #get the batched space
    casecounter = 1
    if(not debug_loops):
      fullrange = modes+1
    else:
      fullrange = 1
    for i in range(fullrange):
      #FIND THE RESTRICTED INDICES IN THE OLD ORDERING
      all_thingys = combinations(idxarr,i)
      #all_thingys = itertools.combinations(idxarr,i)
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
        
        #Build if-statement for current restrictions
        if(not debug_loops):
          conditionalstatement = "      "+label1+"if("
          for j in range(len(conditions)):
            conditionalstatement += conditions[j]
            if j != len(conditions)-1:
              conditionalstatement+=".and."
          conditionalstatement += ")then\n"
          f.write(conditionalstatement)
    
        #WRITE OMP DO STUFF HERE
        if(not debug_loops):
          if(modes-len(oldr)>0):
            ompdo ="        !$OMP DO" 
            if (modes-len(oldr)>1 and (not nocollapse)):
              ompdo += " COLLAPSE("+str(modes-len(oldr))+")\n"
              f.write("#ifndef VAR_LSESSL\n")
              f.write(ompdo)
              f.write("#endif\n")
            else:
              ompdo += "\n"
              f.write("#ifndef VAR_LSESSL\n")
              f.write(ompdo)
              f.write("#endif\n")
 
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
          
        #f.write("!"+str(newu)+"    "+str(oldu)+"\n")
        #f.write("!useold "+str(useold)+"\n")
        #BUILD THE ORDER OF THE LOOPS ACCORDING TO THE PREVIOUS CONDITIONS
        outer =[]
        inner =[]
        for j in range(len(perm)):
          if useold:
            inner.append(idxarr[j])
          else:
            inner.append(perm[j])
        
        if useold:
          for j in range(len(oldu)):
            outer.append(oldu[j])
        else:
          for j in range(len(newu)):
            outer.append(newu[j])
        inneri = []
        outeri = []
        for j in reversed(inner):
          inneri.append(j)
        for j in reversed(outer):
          outeri.append(j)
      
  
        #WRITING THE OUTER FOR LOOPS HERE:
        if(not debug_loops):
          offsetstr="        "
          for j in  range(len(outeri)):
            f.write(offsetstr+"do b"+abc[outeri[j]]+"=1,d"+abc[outeri[j]]+"2,bs\n")
            offsetstr += "  "
          offsetstr = offsetstr[0:-2]
        else:
          offsetstr="    "

        f.write("\n")

        if(not debug_loops):
          offsetstr2 = offsetstr + "  "
          for j in  range(len(outeri)):
            if ((ad != "normal") and (ad != "sp")):
              f.write(offsetstr2+"b"+abc[outeri[j]]+"f = f"+abc[outeri[j]]+" + b"+abc[outeri[j]]+"\n")          
        else:
          offsetstr2 = offsetstr + "  "        
       
        #WRITING THE INNER FOR LOOPS HERE: 
        if(not debug_loops):
          for j in range(modes):
            if(inneri[j] in newr):
              f.write(offsetstr2+"do "+abc[inneri[j]]+"=d"+abc[inneri[j]]+"2+1,d"+abc[inneri[j]]+"\n")
            elif(inneri[j] in newu):
              f.write(offsetstr2+"do "+abc[inneri[j]]+"=0,bcntr\n")
            offsetstr2 += "  "
        else:
          for j in range(modes):
            f.write(offsetstr2+"do "+abc[inneri[j]]+"=1,d"+abc[inneri[j]]+"\n")
            offsetstr2 += "  "
  
        offsetstr2 = offsetstr2[0:-2]
        
        #CENTRAL COPYING AND ADDITION STRING
        newidx = ""
        oldidx = ""
        if(not debug_loops):
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
          offsetstr2 = offsetstr2[0:-2]
        f.write("\n")
  
        #WRITING THE OUTER ENDOFOR HERE:
        if(not debug_loops):
          for j in  range(modes-len(oldr)-1,-1,-1):
            f.write(offsetstr+"enddo\n")
            offsetstr = offsetstr[0:-2]

          if(modes-len(oldr)> 0):
            ompdo ="        !$OMP END DO NOWAIT\n"
            f.write("#ifndef VAR_LSESSL\n")
            f.write(ompdo)
            f.write("#endif\n")

        if(not debug_loops):
         
          conditionalstatement="      endif "+label1+"\n"
          f.write(conditionalstatement)
        
      if(not debug_loops):
        if(i==modes-1):
          ompdo = "      !$OMP END PARALLEL\n"
          f.write("#ifndef VAR_LSESSL\n")
          f.write(ompdo)
          f.write("#endif\n")
  
    if(cas==len(cases)-1):
      f.write("    endif precase\n")
   

def write_subroutine_body_acc(f,idxarr,perm,modes,args,acc_case,ad):
  
  x = 0
  while x < 2:
    x += 1
    #get the batched space
    casecounter = 1
    fullrange = modes+1
    #FIND THE RESTRICTED INDICES IN THE OLD ORDERING
    all_thingys = combinations(idxarr,0)
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
  
      #modify the conditions accordingly
      for j in  range(len(oldr)):
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
        
      #f.write("!"+str(newu)+"    "+str(oldu)+"\n")
      #f.write("!useold "+str(useold)+"\n")
      #BUILD THE ORDER OF THE LOOPS ACCORDING TO THE PREVIOUS CONDITIONS
      outer =[]
      inner =[]
      for j in range(len(perm)):
        if useold:
          inner.append(idxarr[j])
        else:
          inner.append(perm[j])
      
      if useold:
        for j in range(len(oldu)):
          outer.append(oldu[j])
      else:
        for j in range(len(newu)):
          outer.append(newu[j])
      inneri = []
      outeri = []
      for j in reversed(inner):
        inneri.append(j)
      for j in reversed(outer):
        outeri.append(j)
  
      offsetstr="    "
      offsetstr2 = offsetstr + "   "
  
      #WRITING OPENACC DIRECTIVES:
      oaccparallel_init1 = "!$acc parallel present(array_in,array_out)&\n"
      oaccparallel_init2 = "!$acc parallel present(array_in,array_out)&\n"
      if(modes == 4):
        oaccparallel_init1 += "!$acc& firstprivate(pre1,pre2,da,db,dc,dd) private(a,b,c,d) wait(async_id2) async(async_id1)\n"
        oaccparallel_init2 += "!$acc& firstprivate(pre1,pre2,da,db,dc,dd) private(a,b,c,d) async(async_id1)\n"
        oaccloop_gang = "!$acc loop gang collapse(2)\n"
        oaccloop_worker = "!$acc loop worker\n"
        oaccloop_vector = "!$acc loop vector\n"
      elif(modes == 3):
        oaccparallel_init1 += "!$acc& firstprivate(pre1,pre2,da,db,dc) private(a,b,c) wait(async_id2) async(async_id1)\n"
        oaccparallel_init2 += "!$acc& firstprivate(pre1,pre2,da,db,dc) private(a,b,c) async(async_id1)\n"
        oaccloop_gang = "!$acc loop gang\n"
        oaccloop_worker = "!$acc loop worker\n"
        oaccloop_vector = "!$acc loop vector\n"
      elif(modes == 2):
        oaccparallel_init1 += "!$acc& firstprivate(pre1,pre2,da,db) private(a,b) wait(async_id2) async(async_id1)\n"
        oaccparallel_init2 += "!$acc& firstprivate(pre1,pre2,da,db) private(a,b) async(async_id1)\n"
        oaccloop_gang = "!$acc loop gang, worker\n"
        oaccloop_vector = "!$acc loop vector\n"
  
      #WRITING THE INNER FOR LOOPS HERE:
      if (x == 1):
        f.write(offsetstr+"if (wait_arg) then\n\n")
        f.write(oaccparallel_init1)
      elif (x == 2):
        f.write(offsetstr+"else\n\n")
        f.write(oaccparallel_init2)
      f.write(oaccloop_gang)
      for j in range(modes):
        if (modes == 4 and j == 2):
          f.write(oaccloop_worker)
        if (modes == 3 and j == 1):
          f.write(oaccloop_worker)
        if (modes == 4 and j == 3):
          f.write(oaccloop_vector)
        if (modes == 3 and j == 2):
          f.write(oaccloop_vector)
        if (modes == 2 and j == 1):
          f.write(oaccloop_vector)
        f.write(offsetstr2+"do "+abc[inneri[j]]+"=1,d"+abc[inneri[j]]+"\n")
        offsetstr2 += "  "
  
      offsetstr2 = offsetstr2[0:-2]
      
      #CENTRAL COPYING AND ADDITION STRING
      newidx = ""
      oldidx = ""
      for j in range(modes):
        if(perm[j] in newu):
          newidx += abc[perm[j]]+","
        elif(perm[j] in newr):
          newidx += abc[perm[j]]+","
        else:
          print "INVALID STUFF HAPPENING HERE"
        if(idxarr[j] in newu):
          oldidx += abc[idxarr[j]]+","
        elif(idxarr[j] in newr):
          oldidx += abc[idxarr[j]]+","
        else:
          print "INVALID STUFF HAPPENING HERE"
  
      newidx = newidx[0:-1]
      oldidx = oldidx[0:-1]
  
      cpstr=offsetstr2+"  array_out("
  #    if cas == 0:
      if acc_case == 0:
        cpstr += newidx+") = array_in("+oldidx+")\n"
  #    elif cas==1:
      elif acc_case==1:
        cpstr += newidx+") = pre1*array_in("+oldidx+")\n"
  #    elif cas==2:
      elif acc_case==2:
        cpstr += newidx+") = array_out("+newidx+") + array_in("+oldidx+")\n"
  #    elif cas==3:
      elif acc_case==3:
        cpstr += newidx+") = array_out("+newidx+") + pre1*array_in("+oldidx+")\n"
  #    elif cas==4:
      elif acc_case==4:
        cpstr += newidx+") = pre2*array_out("+newidx+") + array_in("+oldidx+")\n"
  #    elif cas==5:
      elif acc_case==5:
        cpstr += newidx+") = pre2*array_out("+newidx+") + pre1*array_in("+oldidx+")\n"
  
      f.write(cpstr)
      
      oaccloop_end = "!$acc end loop\n"
      oaccparallel_end = "!$acc end parallel\n"
  
      for j in  range(modes-1,-1,-1):
        f.write(offsetstr2+"enddo\n")
        if (modes == 4 and j == 1):
          f.write("")
        else:
          f.write(oaccloop_end)
        offsetstr2 = offsetstr2[0:-2]
      f.write(oaccparallel_end)
      if (x == 1):
        f.write("\n")
      elif (x == 2):
        f.write("\n")
        f.write(offsetstr+"endif\n")
        f.write("\n")

#WRITE THE HEADER AND GET THE SUBROUTINE NAME
def write_subroutine_header(f,idxarr,perm,now,modes,ad,deb):
  reordstr1 = ""
  reordstr2 = ""
  reordstr3 = ""
  var_underscore = ""
  for i in range(modes):
    reordstr1 += str(perm[i]+1)
    intermed = "dims("+str(idxarr[i]+1)+"),"
    if ad == "f2t":
      intermed = "f"+intermed
      var_underscore = "_"
    if ad == "sp":
      var_underscore = "_"
    reordstr2 += intermed
    intermed = "dims("+str(perm[i]+1)+"),"
    if ad == "t2f":
      intermed = "f"+intermed
      var_underscore = "_"
    if ad == "sp":
      var_underscore = "_"
    reordstr3 += intermed
  reordstr2 = reordstr2[0:-1]
  reordstr3 = reordstr3[0:-1]
  if (ad == "normal"):
    sname =  "manual_"+reordstr1+"_reordering"
  else:
    sname =  "manual_"+reordstr1+"_reordering"+var_underscore+ad

  #GET THE SUBROUTINE HEADER
  subheaderstr= "  !\> \\brief reorder a "+str(modes)+" diensional array  to get the indices\n"
  subheaderstr+= "  !   in the order "+reordstr1+" , this is a quite expensive reordering\n"
  subheaderstr+= "  !   and thus requires additional attention \n"
  subheaderstr+= "  !\> \\author Patrick Ettenhuber\n"
  subheaderstr+= "  !\> \date "+str(now.month)+", "+str(now.year)+"\n"
  subheaderstr+= "  subroutine "+sname+"(dims,"
  if ((ad != "normal") and (ad != "sp")):
    subheaderstr += "fdims,fels,"
  subheaderstr+= "pre1,array_in,pre2,array_out)\n"
  subheaderstr+= "    implicit none\n"
  subheaderstr+= "    !> input for the block size in tiled reordering\n"
  subheaderstr+= "    integer, parameter :: bs=BS_"+str(modes)+"D\n"
  subheaderstr+= "    !>  the dimensions of the different modes in the original array\n"
  subheaderstr+= "    integer, intent(in) :: dims("+str(modes)+")"
  if ((ad != "normal") and (ad != "sp")):
    subheaderstr += ",fdims("+str(modes)+"),fels("+str(modes)+")"
  subheaderstr+= "\n"
  subheaderstr+= "    !> as this routine can be used for adding and scaling these are the prefactors\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(in) :: pre1,pre2\n"
  else:
    subheaderstr+= "    real(realk),intent(in) :: pre1,pre2\n"
  subheaderstr+= "    !> array to be reordered\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(in) :: array_in("+reordstr2+")\n"
  else:
    subheaderstr+= "    real(realk),intent(in) :: array_in("+reordstr2+")\n"
  subheaderstr+= "    !> reordered array\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(inout) :: array_out("+reordstr3+")\n"
  else:
    subheaderstr+= "    real(realk),intent(inout) :: array_out("+reordstr3+")\n"
  subheaderstr+= "    integer :: bcntr,"
  for i in range(modes):
    subheaderstr+= abc[i]+",b"+abc[i]+",d"+abc[i]+",d"+abc[i]+"2,"
  subheaderstr = subheaderstr[0:-1]
  subheaderstr += "\n"
  if ((ad != "normal") and (ad != "sp")):
    subheaderstr+= "    integer :: "
    for i in range(modes):
      subheaderstr+= "b"+abc[i]+"f,f"+abc[i]+","
    subheaderstr = subheaderstr[0:-1] + "\n"
  subheaderstr+= "    logical :: "
  for i in range(modes):
    subheaderstr+= "mod"+abc[i]+","
  subheaderstr = subheaderstr[0:-1] + "\n"
  subheaderstr += "\n"
  for i in range(modes):
    subheaderstr+= "    d"+abc[i]+"=dims("+str(i+1)+")\n"
  subheaderstr+= "\n"
  if ((ad != "normal") and (ad != "sp")):
    for i in range(modes):
      subheaderstr+= "    f"+abc[i]+"=fels("+str(i+1)+")-1\n"
  subheaderstr+= "\n"
  for i in range(modes):
    subheaderstr+= "    d"+abc[i]+"2=(d"+abc[i]+"/bs)*bs\n"
  subheaderstr+= "\n"
  for i in range(modes):
    subheaderstr+= "    mod"+abc[i]+"=(mod(d"+abc[i]+",bs)>0)\n\n"
  subheaderstr+= "    bcntr=bs-1\n"

  if deb :
    if ad == "t2f":
      subheaderstr += "\n    if(pre2==0.0E0_realk)array_out(&\n"
      for i in range(modes):
        subheaderstr += "                 &fels("+str(perm[i]+1)+"):fels("+str(perm[i]+1)+")+dims("+str(perm[i]+1)+")-1,&\n"
      subheaderstr = subheaderstr[0:-3]+") = 0.0E0_realk\n\n"
    else:
      if ad == "sp":
        subheaderstr += "\n    if(pre2==0.0E0_real_sp)array_out = 0.0E0_real_sp\n\n"
      else:
        subheaderstr += "\n    if(pre2==0.0E0_realk)array_out = 0.0E0_realk\n\n"
  
  f.write(subheaderstr)
  return sname


#WRITE THE HEADER AND GET THE SUBROUTINE NAME
def write_subroutine_header_acc(f,idxarr,perm,now,modes,acc_case,ad):
  reordstr1 = ""
  reordstr2 = ""
  reordstr3 = ""
  var_underscore = ""
  for i in range(modes):
    reordstr1 += str(perm[i]+1)
    intermed = "dims("+str(idxarr[i]+1)+"),"
    reordstr2 += intermed
    intermed = "dims("+str(perm[i]+1)+"),"
    reordstr3 += intermed
  reordstr2 = reordstr2[0:-1]
  reordstr3 = reordstr3[0:-1]
  if (ad == "normal"):
    sname =  "manual_acc_"+reordstr1+"_reordering_"+str(acc_case)
  else:
    var_underscore = '_'
    sname =  "manual_acc_"+reordstr1+"_reordering_"+str(acc_case)+var_underscore+ad

  #GET THE SUBROUTINE HEADER
  subheaderstr= "  !\> \\brief reorder a "+str(modes)+" diensional array (case "+str(acc_case)+") on a device\n"
  subheaderstr+= "  !   to get the indices in the order "+reordstr1+" , this is a quite expensive reordering\n"
  subheaderstr+= "  !   and thus requires additional attention \n"
  subheaderstr+= "  !\> \\author Janus Juul Eriksen & Patrick Ettenhuber\n"
  subheaderstr+= "  !\> \date "+str(now.month)+", "+str(now.year)+"\n"
  subheaderstr+= "  subroutine "+sname+"(dims,"
  subheaderstr+= "pre1,array_in,pre2,array_out,async_id1,async_id2,wait_arg)\n"
  subheaderstr+= "    implicit none\n"
  subheaderstr+= "    !>  the dimensions of the different modes in the original array\n"
  subheaderstr+= "    integer, intent(in) :: dims("+str(modes)+")\n"
  subheaderstr+= "    !> as this routine can be used for adding and scaling these are the prefactors\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(in) :: pre1,pre2\n"
  else:
    subheaderstr+= "    real(realk),intent(in) :: pre1,pre2\n"
  subheaderstr+= "    !> array to be reordered\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(in) :: array_in("+reordstr2+")\n"
  else:
    subheaderstr+= "    real(realk),intent(in) :: array_in("+reordstr2+")\n"
  subheaderstr+= "    !> reordered array\n"
  if ad == "sp":
    subheaderstr+= "    real(real_sp),intent(inout) :: array_out("+reordstr3+")\n"
  else:
    subheaderstr+= "    real(realk),intent(inout) :: array_out("+reordstr3+")\n"
  subheaderstr+= "    integer(acc_handle_kind),intent(in) :: async_id1\n"
  subheaderstr+= "    integer(acc_handle_kind),intent(in) :: async_id2\n"
  subheaderstr+= "    logical,intent(in) :: wait_arg\n"
  subheaderstr+= "    integer :: "
  for i in range(modes):
    subheaderstr+= abc[i]+",d"+abc[i]+","
  subheaderstr = subheaderstr[0:-1]
  subheaderstr += "\n"
  subheaderstr += "\n"
  for i in range(modes):
    subheaderstr+= "    d"+abc[i]+"=dims("+str(i+1)+")\n"
  subheaderstr += "\n"

  f.write(subheaderstr)
  return sname


def write_simple_module_header(f,idim,idx,now,args,kindof):
   f.write("!\> \\brief this autogenerated module is inteded to contain high performance reorderings for\n!mutlidimensional arrays to tiles in a different distribution.\n!\> \\author Patrick Ettenhuber\n!\> \\date March 2013, file produced: "+str(now.month)+", "+str(now.year)+"\n")
   #WRITE THE VARIABLES WITH WHICH THE FILE WAS PRODUCED HERE 
   #--> FOR LATER READOUT AND SEE IF IT IS NECESSARY TO PRODUCE A NEW INSTANCE OF IT
   for i in range(len(args)):
     f.write("!ARG"+str(i)+": "+str(args[i])+"\n")

   f.write("!END VARS\n\n")

   if("f2t" in kindof):
     f.write("module reord"+str(idim)+"d_"+str(idx)+"_utils_f2t_module\n")
   elif("t2f" in kindof):
     f.write("module reord"+str(idim)+"d_"+str(idx)+"_utils_t2f_module\n")
   elif(kindof=="r"):
     f.write("module reord"+str(idim)+"d_"+str(idx)+"_reord_module\n")
   elif(kindof=="sp"):
     f.write("module reord"+str(idim)+"d_"+str(idx)+"_reord_module_sp\n")
   elif(kindof=="acc"):
     f.write("module reord"+str(idim)+"d_acc_reord_module\n")
   elif(kindof=="acc_sp"):
     f.write("module reord"+str(idim)+"d_acc_reord_module_sp\n")
   else:
     print "NO VALID OPTION1 " + kindof
     sys.exit()

   f.write("  use precision\n")
   f.write("  use lsparameters\n")
   if(kindof=="acc"):
     f.write("  use openacc\n")
   if(kindof=="acc_sp"):
     f.write("  use openacc\n")
   f.write("\n")
   f.write("  contains\n")

def write_simple_module_end_and_close(f,idim,idx,now,args,kindof):

   if(kindof=='f2t'):
     f.write("end module reord"+str(idim)+"d_"+str(idx)+"_utils_f2t_module\n")
   elif(kindof=='t2f'):
     f.write("end module reord"+str(idim)+"d_"+str(idx)+"_utils_t2f_module\n")
   elif(kindof=='r'):
     f.write("end module reord"+str(idim)+"d_"+str(idx)+"_reord_module\n")
   elif(kindof=='sp'):
     f.write("end module reord"+str(idim)+"d_"+str(idx)+"_reord_module_sp\n")
   elif(kindof=='acc'):
     f.write("end module reord"+str(idim)+"d_acc_reord_module\n")
   elif(kindof=='acc_sp'):
     f.write("end module reord"+str(idim)+"d_acc_reord_module_sp\n")
   else:
     print "NO VALID OPTION2"+kindof
     sys.exit()
   f.close()







def write_main_header(f,now,args,lsutildir,minr,maxr,skip):
   f.write("!\> \\brief this autogenerated module is inteded to contain high performance reorderings for\n!mutlidimensional arrays.\n!\> \\author Patrick Ettenhuber & Janus Juul Eriksen\n!\> \\date November 2012, file produced: "+str(now.month)+", "+str(now.year)+"\n")
   #WRITE THE VARIABLES WITH WHICH THE FILE WAS PRODUCED HERE 
   #--> FOR LATER READOUT AND SEE IF IT IS NECESSARY TO PRODUCE A NEW INSTANCE OF IT
   for i in range(len(args)):
     f.write("!ARG"+str(i)+": "+str(args[i])+"\n")

   f.write("!END VARS\n\n")

   f.write("module reorder_frontend_module\n")
   f.write("  use precision\n")
   f.write("  use memory_handling\n")
   
   for mode in range(maxr,minr-1,-1):
     for i in range(mode):
       if (mode == 4 or ( mode == 3 or mode == 2 and not skip)):
         f.write("  use reord"+str(mode)+"d_"+str(i+1)+"_utils_t2f_module\n")
         f.write("  use reord"+str(mode)+"d_"+str(i+1)+"_utils_f2t_module\n")
       if(mode==2 and i == 0):
         continue
       f.write("  use ""reord"+str(mode)+"d_"+str(i+1)+"_reord_module\n")
   if(args[4]):
     f.write("#ifdef VAR_OPENACC\n")
     f.write("  use reord2d_acc_reord_module\n")
     f.write("  use reord3d_acc_reord_module\n")
     f.write("  use reord4d_acc_reord_module\n")
     f.write("#endif\n")
   f.write("#ifdef VAR_REAL_SP\n")
   for mode in range(maxr,minr-1,-1):
     for i in range(mode):
       if(mode==2 and i == 0):
         continue
       f.write("  use ""reord"+str(mode)+"d_"+str(i+1)+"_reord_module_sp\n")
   f.write("#endif\n")
   f.write("  use LSTIMING\n")
   #f.write("  contains\n")
   #Write the subroutines called by the user
   basic = open(lsutildir+"/autogen/reorder_header.F90",'r')
   for line in basic:
     f.write(line)
   basic.close


def write_testing_framework(f,minr,maxr):
  header="\
  subroutine test_array_reorderings(LUPRI)\n    implicit none\n\n\
    real(realk),pointer :: in1(:),sto(:)\n\
    real(realk),pointer :: res(:),til(:)\n\
    real(realk) :: ref(6),ref1s,ref2s,ref1,ref2\n\
    integer :: tile_idx,LUPRI\n\
    logical :: master,rigorous\n\
    integer :: "

  for i in range(maxr):
    header += abc[i]+",n"+abc[i]+","
  header = header[0:-1]+"\n    integer :: p1,p2\n\
    real(realk) :: pr1,pr2,begc1,begw1,endc1,endw1,begc2,begw2,endc2,endw2\n\
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
      testcase = testcase[0:-1]+"))>1.0E-11_realk)teststatus=\"FAILED \"\n"
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
          do tile_idx=1,3\n            call tile_from_fort(1.0E0_realk,in1,["
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase = testcase[0:-1] +"],"+str(mode)+",0.0E0_realk,til,tile_idx,["
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
      testcase = testcase[0:-1]+"))>1.0E-11_realk)teststatus=\"FAILED \"\n"
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
          do tile_idx=1,3\n            call tile_from_fort(1.0E0_realk,sto,["
      for i in range(mode):
        testcase += "n"+abc[i]+","
      testcase = testcase[0:-1] +"],"+str(mode)+",0.0E0_realk,til,tile_idx,["
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
      testcase = testcase[0:-1] + "])\n            call tile_in_fort(1.0E0_realk,til,tile_idx,["
      for i in range(mode):
        if(perm[i]==mode-1):
          testcase += "n"+abc[perm[i]]+"/2,"
        else:
          testcase += "n"+abc[perm[i]]+","
      testcase = testcase[0:-1] + "],0.0E0_realk,res,["
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
      testcase = testcase[0:-1]+"))>1.0E-11_realk)teststatus=\"FAILED \"\n"
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



  f.write("  end subroutine test_array_reorderings")

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
