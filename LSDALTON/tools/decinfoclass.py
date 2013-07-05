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

      #DUE TO DIFFERENT OUTPUTS A DIFFERENT NUMBER OF LINES HAS TO BE SKIPPED
      if(fromfrag):
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
      else:
        if(fragtype=="CCSD(T)"):
          skip       = 5
          skip2      = 6
          elfragid   = 1
          elfragpid  = 2
          elensing   = 2
          eldist     = 3
          elenpair   = 4
          elccsd     = 1
          elpt       = 2
          #SECOND ROUND GET THE FRAGMENT INFORMATION
          for i in range(len(filelines)):
            if("-- Atomic fragment energies (CCSD)" in filelines[i]):
              foundlags = True
              for j in range(self.sfragjobs):
                self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
                self.sfrags[j].ecorrocc[elccsd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrvirt[elccsd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrlag[elccsd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrtype[elccsd] = self.ecorrtype[elccsd]
            if("-- Atomic fragment energies (fourth--order E[4])" in filelines[i]):
              foundoccs = True
              for j in range(self.sfragjobs):
                self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
                self.sfrags[j].ecorrocc[elpt]  += float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrtype[elpt] = self.ecorrtype[elpt]
            
            if("-- Pair interaction energies (CCSD)" in filelines[i]):
              foundlagp = True
              for j in range(self.pfragjobs):
                self.pfrags[j].fragid    = int(filelines[i+skip2+j].split()[elfragid])
                self.pfrags[j].fragpid   = int(filelines[i+skip2+j].split()[elfragpid])
                self.pfrags[j].dist      = float(filelines[i+skip2+j].split()[eldist])
                self.pfrags[j].ecorrocc[elccsd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrvirt[elccsd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrlag[elccsd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrtype[elccsd] = self.ecorrtype[elccsd]
            if("-- Pair interaction energies (fourth--order E[4])" in filelines[i]):
              foundoccp = True
              for j in range(self.pfragjobs):
                self.pfrags[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
                self.pfrags[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
                self.pfrags[j].dist      = float(filelines[i+skip+j].split()[eldist])
                self.pfrags[j].ecorrocc[elpt]  += float(filelines[i+skip+j].split()[elenpair])
                self.pfrags[j].ecorrtype[elpt] = self.ecorrtype[elpt]
            #if("-- Pair fragment energies (fifth--order E[5])" in filelines[i]):
            #  foundvirtp = True
            #  for k in range(self.sfragjobs):
            #    for j in range(1,k):
            #      cb = j + k * self.sfragjobs
            #      tr = j + k * ((j*(j+1))/2)  
            #      #print j,k,cb,tr,self.pfragjobs
            #      #if(int(filelines[i+skip+cb].split()[elfragid])>int(filelines[i+skip+cb].split()[elfragpid])):
            #      #  self.pfrags[j].fragid    = int(filelines[i+skip+cb].split()[elfragid])
            #      #  self.pfrags[j].fragpid   = int(filelines[i+skip+cb].split()[elfragpid])
            #      #  self.pfrags[j].dist      = float(filelines[i+skip+cb].split()[eldist])
            #      #  self.pfrags[j].ecorrocc[elpt] += float(filelines[i+skip+cb].split()[elenpair])
            #      #  self.pfrags[j].ecorrtype[elpt] = self.ecorrtype[elpt]
            #    for j in range(k,self.sfragjobs):
            #      cb = j + k * self.sfragjobs
            #      #print j,k,cb,tr,self.pfragjobs
            #      #if(int(filelines[i+skip+cb].split()[elfragid])>int(filelines[i+skip+cb].split()[elfragpid])):
            #      #if(int(filelines[i+skip+cb].split()[elfragid])==int(filelines[i+skip+cb].split()[elfragpid])):
            #      #  self.sfrags[j].ecorrocc[elpt] += float(filelines[i+skip+j].split()[elenpair])
            #      #elif(int(filelines[i+skip+cb].split()[elfragid])<int(filelines[i+skip+cb].split()[elfragpid])):
            #      #  self.pfrags[j].fragid    = int(filelines[i+skip+j].split()[elfragid])
            #      #  self.pfrags[j].fragpid   = int(filelines[i+skip+j].split()[elfragpid])
            #      #  self.pfrags[j].dist      = float(filelines[i+skip+j].split()[eldist])
            #      #  self.pfrags[j].ecorrocc[elpt] += float(filelines[i+skip+j].split()[elenpair])
            #      #  self.pfrags[j].ecorrtype[elpt] = self.ecorrtype[elpt]
        if(fragtype=="CCD"):
          skip       = 5
          skip2      = 6
          elfragid   = 1
          elfragpid  = 2
          elensing   = 2
          eldist     = 3
          elenpair   = 4
          elccd      = 0
          elpt       = 2
          #SECOND ROUND GET THE FRAGMENT INFORMATION
          for i in range(len(filelines)):
            if("-- Atomic fragment energies (CCD)" in filelines[i]):
              foundlags = True
              for j in range(self.sfragjobs):
                self.sfrags[j].fragid       = int(filelines[i+skip+j].split()[elfragid])
                self.sfrags[j].ecorrocc[elccd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrvirt[elccd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrlag[elccd]  = float(filelines[i+skip+j].split()[elensing])
                self.sfrags[j].ecorrtype[elccd] = self.ecorrtype[elccd]
            
            if("-- Pair interaction energies (CCD)" in filelines[i]):
              foundlagp = True
              for j in range(self.pfragjobs):
                self.pfrags[j].fragid    = int(filelines[i+skip2+j].split()[elfragid])
                self.pfrags[j].fragpid   = int(filelines[i+skip2+j].split()[elfragpid])
                self.pfrags[j].dist      = float(filelines[i+skip2+j].split()[eldist])
                self.pfrags[j].ecorrocc[elccd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrvirt[elccd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrlag[elccd]  = float(filelines[i+skip2+j].split()[elenpair])
                self.pfrags[j].ecorrtype[elccd] = self.ecorrtype[elccd]
        elif(fragtype=="MP2"):
          skip       = 4
          skip2      = 4
          elfragid   = 1
          elfragpid  = 2
          elensing   = 2
          eldist     = 3
          elenpair   = 4
          elmp2      = 0
          foundsinglesocc  = False
          foundsinglesvirt = False
          foundpairsocc    = False
          foundpairsvirt   = False
          #SECOND ROUND GET THE FRAGMENT INFORMATION
          for i in range(len(filelines)):
            if("-- Atomic fragments - the four contributions"in filelines[i]):
              j=1
              while(not foundsinglesocc or not foundsinglesvirt):
                line = filelines[i+j].strip()
                if("-- Contribution 1" == line):
                  foundsinglesocc = True
                  for k in range(self.sfragjobs):
                    self.sfrags[k].fragid            =   int(filelines[i+skip+j+k].split()[elfragid])
                    self.sfrags[k].ecorrocc[elmp2]   = float(filelines[i+skip+j+k].split()[elensing])
                    self.sfrags[k].ecorrtype[elmp2]  = self.ecorrtype[elmp2]
                if("-- Contribution 3" == line):
                  foundsinglesvirt = True
                  for k in range(self.sfragjobs):
                    self.sfrags[k].fragid            =   int(filelines[i+skip+j+k].split()[elfragid])
                    self.sfrags[k].ecorrvirt[elmp2]  = float(filelines[i+skip+j+k].split()[elensing])
                    self.sfrags[k].ecorrtype[elmp2]  = self.ecorrtype[elmp2]
                j += 1

            if("-- Pair fragments - the four contributions"in filelines[i]):
              j=1
              while(not foundpairsocc or not foundpairsvirt):
                line = filelines[i+j].strip()
                if("-- Contribution 1" == line):
                  foundpairsocc = True
                  for k in range(self.pfragjobs):
                    self.pfrags[k].fragid    = int(filelines[i+skip2+j+k].split()[elfragid])
                    self.pfrags[k].fragpid   = int(filelines[i+skip2+j+k].split()[elfragpid])
                    self.pfrags[k].dist      = float(filelines[i+skip2+j+k].split()[eldist])
                    self.pfrags[k].ecorrocc[elmp2]  = float(filelines[i+skip2+j+k].split()[elenpair])
                    self.pfrags[k].ecorrtype[elmp2] = self.ecorrtype[elmp2]
                if("-- Contribution 3" == line):
                  foundpairsvirt = True
                  for k in range(self.pfragjobs):
                    self.pfrags[k].fragid    = int(filelines[i+skip2+j+k].split()[elfragid])
                    self.pfrags[k].fragpid   = int(filelines[i+skip2+j+k].split()[elfragpid])
                    self.pfrags[k].dist      = float(filelines[i+skip2+j+k].split()[eldist])
                    self.pfrags[k].ecorrvirt[elmp2] = float(filelines[i+skip2+j+k].split()[elenpair])
                    self.pfrags[k].ecorrtype[elmp2] = self.ecorrtype[elmp2]
                j += 1

    


      #OUTSIDE SECOND LOOP    
        
     

