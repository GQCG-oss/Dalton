from output_plot import *
#/usr/bin/python

# CLASS DEFINITIONS FOR DEC SPECIFIC INFO 
# AUTHOR: PATRICK ETTENHUBER
# YEAR: 2013
# EXPLANATION:
#    This file contains all classes for dec-specific information storage

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

class decinfo_class:
   """A class used for storing all DEC related info"""
   #JUST INIT THE EMPTY CLASS AND FILL IT VIA OUTPUT
   def __init__(self):
      self.fotint     = 0
      self.fotfloat   = 0
      self.nfragjobs  = 0
      self.sfragjobs  = 0
      self.pfragjobs  = 0
      self.sfrags     = []
      self.pfrags     = []
      self.frags      = []
      self.ecorrocc   = []
      self.ecorrvirt  = []
      self.ecorrlag   = []
      self.ecorrtype  = []
      self.esterr     = 0.0
      self.enable_fragread = False

   ####################################
   #DEFINE DEC ANALYSIS OPERATIONS HERE
   ####################################

   plot_pair_energies = plot_pair_energies

   #READ DEC SPECIFIC INFO
   def get_dec_info(self,filelines,fragtype,fromfrag):
      
      if(fragtype=="MP2"):
        self.ecorrtype.append("MP2")
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
        self.ecorrtype.append("CCSD(T)")
        self.ecorrtype.append("CCSD")
        self.ecorrtype.append("(T)")
        self.ecorrocc.append(0.0)
        self.ecorrocc.append(0.0)
        self.ecorrocc.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrvirt.append(0.0)
        self.ecorrlag.append(0.0)
        self.ecorrlag.append(0.0)
        self.ecorrlag.append(0.0)

      #FIRST ROUND GET BASIC INFO 
      found_sf = False
      found_pf = False
      found_nf = False

      allfound = False

      for i in range(len(filelines)):
        line = filelines[i]
        if(fromfrag):
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
          #READ DEC CORRELATION ENEGIES
          if ("Lagrangian scheme energy      :"    in line):
            self.ecorrlag[0]  = float(line.split()[-1])
          if ("Occupied scheme energy        :"    in line):
            self.ecorrocc[0]  = float(line.split()[-1])
          if ("Virtual scheme energy         :"    in line):
            self.ecorrvirt[0] = float(line.split()[-1])
          if("CCSD occupied correlation energy :"  in line):
            self.ecorrocc[1]  = float(line.split()[-1])
          if("CCSD virtual  correlation energy :"  in line):
            self.ecorrvirt[1] = float(line.split()[-1])
          if("(T) occupied correlation energy :"   in line):
            self.ecorrocc[2]  = float(line.split()[-1])
          if("(T) virtual  correlation energy :"   in line):
            self.ecorrvirt[2] = float(line.split()[-1])
          if ("*** Estimated intrinsic error :"    in line):
            self.esterr = float(line.split()[-1])
        else:
          #READ NUMBER OF JOBS
          if ("ORBITAL DISTRIBUTION INFORMATION" in line):
            j=4
            while("Total:" not in filelines[i+j]):
              if(int(filelines[i+j].split()[-2]) != 0 and int(filelines[i+j].split()[-2]) != 0):
                self.sfragjobs += 1
              j += 1
          #READ REFERENCE FULL CALCULATION
          if("Total CCSD(T) correlation energy        ="  in line):
            self.ecorrocc[0]  = float(line.split()[-1])
          if("Total CCSD correlation energy           ="  in line):
            self.ecorrocc[1]  = float(line.split()[-1])
          if("Total CCSD(T) energy contribution       ="  in line):
            self.ecorrocc[2]  = float(line.split()[-1])
          if("Total CCD correlation energy           ="  in line):
            self.ecorrocc[1]  = float(line.split()[-1])

        #allfound = (found_sf and found_pf and found_nf)
        #if(allfound):
        #  break

      if(not fromfrag):
        self.pfragjobs = self.sfragjobs*(self.sfragjobs - 1)/2
        self.nfragjobs = self.pfragjobs + self.sfragjobs
        nfithorderfragjobs = self.sfragjobs*self.sfragjobs

      #print fragtype,self.ecorrocc,self.ecorrvirt,self.ecorrlag,self.pfragjobs,self.sfragjobs


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
      
        for k in range(len(self.ecorrtype)):
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


      #OUTSIDE SECOND LOOP    
        
     

