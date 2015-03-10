#/usr/bin/python

# CLASS DEFINITIONS FOR CCC SPECIFIC INFO 
# AUTHOR: PATRICK ETTENHUBER
# YEAR: 2013
# EXPLANATION:
#    This file contains all classes for dec-specific information storage

class ccinfo_class:
   """A class used for storing all CCC related info"""
   #JUST INIT THE EMPTY CLASS AND FILL IT VIA OUTPUT
   def __init__(self):
      self.ecorrtye  = ""
      self.ecorr     = 0.0
      self.etot      = 0.0
      self.converged = False
      self.t2norm    = 0.0
      self.t1norm    = 0.0
      self.ttotnorm  = 0.0

   ####################################
   #DEFINE DEC ANALYSIS OPERATIONS HERE
   ####################################


   #READ CC SPECIFIC INFO
   def get_cc_info(self,filelines,cctype):
      
      self.ecorrtype = cctype
      for i in range(len(filelines)):
         line = filelines[i]
         if "E: Correlation energy" in line:
            self.ecorr = float(filelines[i].split()[-1])
            self.etot  = float(filelines[i+1].split()[-1])

         if "Hooray! CC equation is solved!" in line:
            self.converged = True
         
         if "Singles amplitudes norm" in line:
            self.t1norm = float(line.split()[-1])

         if "Doubles amplitudes norm" in line:
            self.t2norm = float(line.split()[-1])

         if "Total amplitudes norm" in line:
            self.ttotnorm = float(line.split()[-1])


   #PRINT CC INFO TO FILE
   def print_cc_info(self,cc_out):
      cc_out.write(self.ecorrtye)
      cc_out.write("Correlation energy     = {:10.6e} \n".format(self.ecorr))
      cc_out.write("Total energy           = {:10.6e} \n".format(self.etot))
      cc_out.write("CC equation converged  = "+str(self.converged)+"\n")
      cc_out.write("Double amplitudes norm = {:10.6e} \n".format(self.t2norm))
      cc_out.write("Single amplitudes norm = {:10.6e} \n".format(self.t1norm))
      cc_out.write("Total amplitudes norm  = {:10.6e} \n".format(self.ttotnorm))
      

