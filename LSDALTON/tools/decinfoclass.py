#/usr/bin/python

# CLASS DEFINITIONS FOR DEC SPECIFIC INFO 
# AUTHOR: PATRICK ETTENHUBER
# YEAR: 2013
# EXPLANATION:
#    This file contains all classes for dec-specific information storage

from output_plot import *

# All fragment specific 
class fragment_class:
   def __init__(self):
      self.fragid    = 0
      self.fragpid   = 0
      self.dist      = 0.0
      self.ecorrocc  = []
      self.ecorrvirt = []
      self.ecorrlag  = []
      self.ecorrtype = []

   def print_frag_info(self,frag_out):
      for t, e in zip(self.ecorrtype, self.ecorrocc):
         frag_out.write("{0:7s}   {1:4d}   {2:4d}   {3:8.4f}     Eocc = {4:10.6e}\n".format(
            t, self.fragid, self.fragpid, self.dist, e))
      for t, e in zip(self.ecorrtype, self.ecorrvirt):
         frag_out.write("{0:7s}   {1:4d}   {2:4d}   {3:8.4f}     Evir = {4:10.6e}\n".format(
            t, self.fragid, self.fragpid, self.dist, e))
      for t, e in zip(self.ecorrtype, self.ecorrlag):
         frag_out.write("{0:7s}   {1:4d}   {2:4d}   {3:8.4f}     Elag = {4:10.6e}\n".format(
            t, self.fragid, self.fragpid, self.dist, e))


class decinfo_class:
   """A class used for storing all DEC related info"""
   #JUST INIT THE EMPTY CLASS AND FILL IT VIA OUTPUT
   def __init__(self):
      self.fotfloat   = 0
      self.nfragjobs  = 0
      self.sfragjobs  = 0
      self.pfragjobs  = 0
      self.nesti      = 0
      self.sfrags     = []
      self.pfrags     = []
      self.esti       = []
      self.jobs       = []
      self.ecorrocc   = []
      self.ecorrvirt  = []
      self.ecorrlag   = []
      self.ecorrtype  = []
      self.esterr     = 0.0
      self.OccSize    = [0.0,0.0,0.0]
      self.VirSize    = [0.0,0.0,0.0]
      self.enable_fragread = False

   ####################################
   #DEFINE DEC ANALYSIS OPERATIONS HERE
   ####################################

   plot_pair_energies = plot_pair_energies
   plot_SF_energy_errors = plot_SF_energy_errors

   #READ DEC SPECIFIC INFO
   def get_dec_info(self,filelines,fragtype,fromfrag):
      
      if(fragtype=="RPA"):
        self.ecorrtype.append("RPA")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="MP2"):
        self.ecorrtype.append("MP2")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="RIMP2"):
        self.ecorrtype.append("RI-MP2")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="CC2"):
        self.ecorrtype.append("CC2")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="CCD"):
        self.ecorrtype.append("CCD")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="CCSD"):
        self.ecorrtype.append("CCSD")
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
      elif(fragtype=="CCSD(T)"):
        self.ecorrtype.append("CCSD")
        self.ecorrtype.append("(T)")
        self.ecorrtype.append("CCSD(T)")
        self.ecorrocc.append(0.0)
        self.ecorrocc.append(0.0)
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
        self.ecorrlag.append(0.0)
        self.ecorrlag.append(0.0)
      else:
        print "ERROR(get_dec_info): fragtype,",fragtype,"not understood"
        exit()

      #FIRST ROUND GET BASIC INFO 
      found_fot = False
      found_sf = False
      found_pf = False
      found_nf = False
      found_Ec = False
      found_os = False
      found_vs = False
      esti     = False

      allfound = False

      for i in range(len(filelines)):
        line = filelines[i]
        if(fromfrag):
          #READ FOT:
          if ("FOT (Fragment Optimization Threshold)" in line):
            self.fotfloat = float(line.split()[-1])
            found_fot = True
          #READ NUMBER OF JOBS
          if ("DEC JOB SUMMARY: Number of single jobs =" in line):
            self.sfragjobs = int(line.split()[-1])
            found_sf = True
          if ("DEC JOB SUMMARY: Number of pair jobs   =" in line):
            self.pfragjobs = int(line.split()[-1])
            found_pf = True
          if ("DEC JOB SUMMARY: Total number of jobs  =" in line):
            self.nfragjobs = int(line.split()[-1])
            found_nf = True
          if ("SUMMARY FOR PAIR ESTIMATE ANALYSIS" in line):
            esti = True
          #READ DEC CORRELATION ENEGIES
          # RPA
          if("RPA occupied   correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("RPA virtual    correlation energy :"  in line):
            self.ecorrvirt[0] = float(line.split()[-1])
            found_Ec = True
          # MP2
          if("MP2 occupied   correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("MP2 virtual    correlation energy :"  in line):
            self.ecorrvirt[0] = float(line.split()[-1])
            found_Ec = True
          if("MP2 Lagrangian correlation energy :"  in line):
            self.ecorrlag[0] = float(line.split()[-1])
            found_Ec = True
          # RIMP2
          if("RI-MP2 occupied   correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("RI-MP2 virtual    correlation energy :"  in line):
            self.ecorrvirt[0] = float(line.split()[-1])
            found_Ec = True
          if("RI-MP2 Lagrangian correlation energy :"  in line):
            self.ecorrlag[0] = float(line.split()[-1])
            found_Ec = True
          # CCSD
          if("CCSD occupied correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("CCSD virtual  correlation energy :"  in line):
            self.ecorrvirt[0] = float(line.split()[-1])
            found_Ec = True
          # (T)
          if(" (T) occupied correlation energy :"   in line):
            self.ecorrocc[1]  = float(line.split()[-1])
            found_Ec = True
          if(" (T) virtual  correlation energy :"   in line):
            self.ecorrvirt[1] = float(line.split()[-1])
            found_Ec = True
          # CCSD(T)
          if("Total CCSD(T) occupied correlation energy :"   in line):
            self.ecorrocc[2]  = float(line.split()[-1])
            found_Ec = True
          if("Total CCSD(T) virtual  correlation energy :"   in line):
            self.ecorrvirt[2] = float(line.split()[-1])
            found_Ec = True
          # Estimated intrinsic error
          if ("*** Estimated intrinsic error :"    in line):
            self.esterr = float(line.split()[-1])
          #READ FRAGMENT SIZES:
          # if the estimates size are present they are first read
          # and then overwritten by the final sizes
          #read occ size
          if ("FRAGANALYSIS: Max/Ave/Min occ         :" in line):
            self.OccSize[0] = float(line.split()[-5]) # read max occ size
            self.OccSize[1] = float(line.split()[-3]) # read ave occ size
            self.OccSize[2] = float(line.split()[-1]) # read min occ size
            found_os = True
          #read vir size
          if ("FRAGANALYSIS: Max/Ave/Min unocc       :" in line):
            self.VirSize[0] = float(line.split()[-5]) # read max vir size
            self.VirSize[1] = float(line.split()[-3]) # read ave vir size
            self.VirSize[2] = float(line.split()[-1]) # read min vir size
            found_vs = True
        
          allfound = (found_fot and found_sf and found_pf and found_nf and found_Ec and found_os and found_vs)

        else:
          #READ NUMBER OF JOBS
          if ("FULL JOB SUMMARY: Number of single jobs =" in line):
            self.sfragjobs = int(line.split()[-1])
            found_sf = True
          if ("FULL JOB SUMMARY: Number of pair jobs   =" in line):
            self.pfragjobs = int(line.split()[-1])
            found_pf = True
          if ("FULL JOB SUMMARY: Total number of jobs  =" in line):
            self.nfragjobs = int(line.split()[-1])
            found_nf = True

          #READ REFERENCE FULL CALCULATION
          if("RPA correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("MP2 correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("RI-MP2 correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("CC2 correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("CCD correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if("CCSD correlation energy :"  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
            found_Ec = True
          if(" (T) correlation energy  :"  in line):
            self.ecorrocc[1]  = float(line.split()[-1])
            found_Ec = True
          if("CCSD(T) correlation energy :"  in line):
            self.ecorrocc[2]  = float(line.split()[-1])
            found_Ec = True

          allfound = (found_sf and found_pf and found_nf and found_Ec)

      if (not allfound):
        print "WARNING: Some basic DEC/CC information has not been found\n"


      if (esti):
        self.nesti = self.sfragjobs*(self.sfragjobs-1)/2 # tot num of pairs
        for i in range(self.nesti):
          self.esti.append(fragment_class())
          self.esti[i].ecorrocc.append(0.0)
      for i in range(self.nfragjobs):
        self.jobs.append(fragment_class())

      for i in range(self.pfragjobs):
        self.pfrags.append(fragment_class())
        for j in range(len(self.ecorrtype)):
          self.pfrags[i].ecorrocc.append(0.0)
          self.pfrags[i].ecorrvirt.append(0.0)
          self.pfrags[i].ecorrlag.append(0.0)
          self.pfrags[i].ecorrtype.append(self.ecorrtype[j])
      for i in range(self.sfragjobs):
        self.sfrags.append(fragment_class())
        for j in range(len(self.ecorrtype)):
          self.sfrags[i].ecorrocc.append(0.0)
          self.sfrags[i].ecorrvirt.append(0.0)
          self.sfrags[i].ecorrlag.append(0.0)
          self.sfrags[i].ecorrtype.append(self.ecorrtype[j])

      
      foundlags = False
      foundoccs = False
      foundvirts = False
      foundlagp = False
      foundoccp = False
      foundvirtp = False

      #SETTING THE OFFSETS FOR READING HERE
      skip       = 4
      elfragid   = 0
      elfragpid  = 1
      elensing   = 1
      eldist     = 2
      elenpair   = 3

      #SECOND ROUND GET THE FRAGMENT INFORMATION
      for i in range(len(filelines)):
        #EXCLUDE THE 4th and 5th ORDER ENERGIES
        exclude = False
        if("(fourth order)" in filelines[i] or "(fifth order)" in filelines[i]):
          exclude = True

        # Get estimated pairs info:
        if("Estimated occupied pair energies" in filelines[i] and not exclude):
          for j in range(self.nesti):
            self.esti[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
            self.esti[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
            self.esti[j].dist      = float(filelines[i+skip+j].split()[2])
            self.esti[j].ecorrocc[0]  = float(filelines[i+skip+j].split()[elenpair])
            self.esti[j].ecorrtype.append(self.ecorrtype[0])

        # Get job size (stored in dist, I know its ugly):
        if("DEC FRAGMENT JOB LIST" in filelines[i]):
           o = 5 # offset
           eljobsize = 1
           eljobid1  = 2
           eljobid2  = 3
           for j in range(self.nfragjobs):
              self.jobs[j].dist    = int(filelines[i+o+j].split()[eljobsize])
              self.jobs[j].fragid  = int(filelines[i+o+j].split()[eljobid1])
              # not valid for single fragments!!!
              self.jobs[j].fragpid = int(filelines[i+o+j].split()[eljobid2])

        for k in range(len(self.ecorrtype)):
          #if("pair energies" in filelines[i]):
          if(self.ecorrtype[k]+" Lagrangian single energies" in filelines[i]):
            foundlags = True
            for j in range(self.sfragjobs):
              self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
              self.sfrags[j].ecorrlag[k]  = float(filelines[i+skip+j].split()[elensing])
              self.sfrags[j].ecorrtype[k] = self.ecorrtype[k]
          if(self.ecorrtype[k]+" occupied single energies" in filelines[i] and not exclude):
            foundoccs = True
            for j in range(self.sfragjobs):
              self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
              self.sfrags[j].ecorrocc[k]  = float(filelines[i+skip+j].split()[elensing])
              self.sfrags[j].ecorrtype[k] = self.ecorrtype[k]
          if(self.ecorrtype[k]+" virtual single energies" in filelines[i] and not exclude):
            foundvirts = True
            for j in range(self.sfragjobs):
              self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
              self.sfrags[j].ecorrvirt[k] = float(filelines[i+skip+j].split()[elensing])
              self.sfrags[j].ecorrtype[k] = self.ecorrtype[k]
         
          if(self.ecorrtype[k]+" Lagrangian pair energies" in filelines[i] and not exclude):
            foundlagp = True
            for j in range(self.pfragjobs):
              self.pfrags[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
              self.pfrags[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
              self.pfrags[j].dist      = float(filelines[i+skip+j].split()[2])
              self.pfrags[j].ecorrlag[k]  = float(filelines[i+skip+j].split()[elenpair])
              self.pfrags[j].ecorrtype[k] = self.ecorrtype[k]
          if(self.ecorrtype[k]+" occupied pair energies" in filelines[i] and not exclude):
            foundoccp = True
            for j in range(self.pfragjobs):
              self.pfrags[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
              self.pfrags[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
              self.pfrags[j].dist      = float(filelines[i+skip+j].split()[2])
              self.pfrags[j].ecorrocc[k]  = float(filelines[i+skip+j].split()[elenpair])
              self.pfrags[j].ecorrtype[k] = self.ecorrtype[k]
          if(self.ecorrtype[k]+" virtual pair energies" in filelines[i] and not exclude):
            foundvirtp = True
            for j in range(self.pfragjobs):
              self.pfrags[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
              self.pfrags[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
              self.pfrags[j].dist      = float(filelines[i+skip+j].split()[2])
              self.pfrags[j].ecorrvirt[k] = float(filelines[i+skip+j].split()[elenpair])
              self.pfrags[j].ecorrtype[k] = self.ecorrtype[k]

      found_s = (foundlags or foundoccs or foundvirts)
      found_p = (foundlagp or foundoccp or foundvirtp)
      if (not (found_s or found_p)):
          print "WARNING: Single and/or pair fragment energies have not been found\n"

      #OUTSIDE SECOND LOOP    
        
     
   #PRINT DEC INFO TO FILE
   def print_dec_info(self,dec_out):
      for t, e in zip(self.ecorrtype, self.ecorrocc):
         dec_out.write("{0:7s} occupied correlation energy   = {1:10.6e} \n".format(t,e))
      for t, e in zip(self.ecorrtype, self.ecorrvirt):
         dec_out.write("{0:7s} virtual correlation energy    = {1:10.6e} \n".format(t,e))
      for t, e in zip(self.ecorrtype, self.ecorrlag):
         dec_out.write("{0:7s} lagrangian correlation energy = {1:10.6e} \n".format(t,e))

      dec_out.write("\nFOT = {0:6.2e} \n".format(self.fotfloat))

      dec_out.write("\nPair estimates: \n")
      dec_out.write("Model      Atom1  Atom2    dist        DEC energy\n")
      for e in self.esti:
         e.print_frag_info(dec_out)

      dec_out.write("\nSingle fragment: \n")
      dec_out.write("Model      Atom1  Atom2    dist        DEC energy\n")
      for s in self.sfrags:
         s.print_frag_info(dec_out)

      dec_out.write("\nPair fragment: \n")
      dec_out.write("Model      Atom1  Atom2    dist        DEC energy\n")
      for p in self.pfrags:
         p.print_frag_info(dec_out)

      


