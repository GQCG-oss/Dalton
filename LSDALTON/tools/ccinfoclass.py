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


   #READ DEC SPECIFIC INFO
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
            self.t1norm = line.split()[-1]

         if "Doubles amplitudes norm" in line:
            self.t2norm = line.split()[-1]

         if "Total amplitudes norm" in line:
            self.ttotnorm = line.split()[-1]
