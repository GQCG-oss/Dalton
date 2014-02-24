#/usr/bin/python
from output_plot import *
from decinfoclass import *
from ccinfoclass import *

# EXTENDABLE LSDALTON.OUT PARSER FOR SIMPLE DATA COLLECTION
# AUTHOR: PATRICK ETTENHUBER
# YEAR: 2013
# EXPLANATION:
#    This is a parser class inteded for both developers and users for simple
#    data collection from an output file. This class can be used in any python
#    script by importing it at the beginning. 
#     - For python-newbes: write "from /path/to/this/filewithout.py import * "
#       at the beginning of your scritp, alternatively you may set $PYTHONPATH
#       to include the directory this scritpt is contained and just "form filewithout.py import *"
#
# USE OF THE CLASS:
#    If you included the class written above you can just declare an object of the class as
#    varname = lsoutput("/path/to/LSDALTON.OUT") 
#    and the parser should read all the info needed. Please extend the class if the info you
#    are looking for is not available. This is preferably done by adding new subclasses like
#    I have done for the decinfo subclass. So we can keep things separated. The Idea is to fill
#    the subclasses depending on the input file (which is read at the beginning from the LSDALTON.OUT).
#    I did not care too much about efficiency right now, both in memory and computation, so instead
#    of reading the file over and over, it is read once at the beginning and stored in the .lines variable.
#    Add specific functions at the end of the class. If it grows out of hand one probably has to declare the
#    functions outside of the class definition.

class lsoutput:
   """A class used for analyzing the LSDALTON.OUT files"""
   def __init__(self,pathtofile):
      # READ THE FILE AND FILL CLASS
      tmpstream = open(pathtofile,'r')
      self.path=pathtofile
      self.lines=tmpstream.readlines()
      tmpstream.close()
      
      ########################################### 
      #DECLARATION OF EMPTY ELEMENTS OF THE CLASS
      ########################################### 
      self.molinp = []
      self.dalinp = []

      #A CHARACTERISTIC STRING THAT CONTAINS THE ESSENTIALS
      self.calctype = [""]*2

      #DEC SPECIFIC CALULATION INFO AND OPERATIONS
      self.decinfo = decinfo_class()
      #CCC SPECIFIC CALULATION INFO AND OPERATIONS
      self.ccinfo  = ccinfo_class()

      #NUMBER OF BASIS FUNCTIONS
      self.nb = 0
      #NUMBER OF VIRTUAL ORBITALS
      self.nv = 0
      #NUMBER OF OCCUPIED ORBITALS
      self.no = 0

      ###########################################################
      #PARSE THE OUTPUT FILE LINE BY LINE AND EXTRACT INFORMATION
      ###########################################################
      found_molinp = False
      found_dalinp = False
      found_nb     = False
      found_nv     = False
      found_no     = False

      allfound     = False


      for i in range(len(self.lines)):
         line = self.lines[i].strip().upper()
       
         #SAVE THE MOLECULE.INP from the outputfile
         if("PRINTING THE MOLECULE.INP FILE" in line):
           found_molinp = True
           lineparser = ""
           j=i+2
           
           while ("PRINTING THE LSDALTON.INP FILE" not in lineparser):
             lineparser=self.lines[j].strip()
             self.molinp.append(lineparser)
             j+=1
           #TRIM THE MOLECULE.INP
           del self.molinp[-1]
           del self.molinp[-1]

         #SAVE LSDALTON.INP from the outputfile
         if("PRINTING THE LSDALTON.INP FILE" in line):
           found_dalinp = True
           lineparser = ""
           j=i+2
           
           while ("*END OF INPUT" not in lineparser):
             lineparser=self.lines[j].strip().upper()
             self.dalinp.append(lineparser)
            
             #SET DEC SPECIFIC CALCULTION INFO FROM INPUT FILE 
             if("**DEC" in lineparser):
               self.calctype[0] = "DEC"
               self.decinfo.enable_fragread = True
               lineparser2 = ""
               k=j+1
               found = False
               while("*" not in lineparser2):
                 lineparser2=self.lines[k].strip().upper()
                 if(".MP2" == lineparser2):
                   self.calctype[1] = "MP2"
                   found = True
                 elif(".MP2DEBUG" == lineparser2):
                   self.calctype.append("MP2DEBUG")
                   found = True
                 elif(".CCSD(T)" == lineparser2):
                   self.calctype[1] = "CCSD(T)"
                   found = True
                 elif(".CCSD" == lineparser2):
                   self.calctype[1] = "CCSD"
                   found = True
                 elif(".CCD" == lineparser2):
                   self.calctype[1] = "CCD"
                   found = True
                 elif(".CC2" == lineparser2):
                   self.calctype[1] = "CC2"
                   found = True
                 k+=1
               if(not found):
                 self.calctype[1] = "NONE"

               lineparser2 = ""
               k=j+1
               found = False
               while("*" not in lineparser2 and not found):
                 lineparser2=self.lines[k].strip().upper()
                 if(".FOT" in lineparser2):
                   self.calctype.append(" FOT="+self.lines[k+1].strip().upper())
                   self.decinfo.fotint   = int(self.lines[k+1].strip())
                   found = True
                 k+=1
               if(not found):
                 #SET THE DEFAULTS
                 self.calctype.append("FOT=4")
                 self.decinfo.fotint   = 4
               #USE THE INT TO CALCULATE THE 
               self.decinfo.fotfloat = 10**(-self.decinfo.fotint)

             #SET CC SPECIFIC CALCULTION INFO FROM INPUT FILE 
             if("**CC" in lineparser):
               self.calctype[0] = "CC"
               lineparser2 = ""
               k=j+1
               found = False
               while("*" not in lineparser2):
                 lineparser2=self.lines[k].strip().upper()
                 if(".MP2" in lineparser2):
                   self.calctype[1] = " MP2"
                   found = True
                 elif(".CCSD(T)" == lineparser2):
                   self.calctype[1] = "CCSD(T)"
                   found = True
                 elif(".CCSD" == lineparser2):
                   self.calctype[1] = "CCSD"
                   found = True
                 elif(".CCD" == lineparser2):
                   self.calctype[1] = "CCD"
                   found = True
                 elif(".CC2" == lineparser2):
                   self.calctype[1] = "CC2"
                   found = True
                 elif(".PRINTFRAGS"):
                   self.decinfo.enable_fragread = True
                 k+=1
               if(not found):
                 self.calctype[1] = " NONE"


             #ITER ONE LINE IN SCANNING THE LSDALTON.INP in the outputfile
             j+=1

         if("REGULAR BASISFUNCTIONS             :" in line):
           self.nb = int(line.strip().split()[-1])
           found_nb = True
  
         #RIGHT NOW TAKE THE OCCUPIED NUMBER FROM DEC PRINT, SHOULD BE CHANGED
         if("FULL: NUMBER OF ELECTRONS        :" in line):
           self.no = int(line.strip().split()[-1])
           found_no = True
         
         #THE PRELIMINARY FILLING IS FINISHED WHEN MOLINP AND DALINP ARE 
         #FILLED IN A SECOND LOOP THE DETAILS HAVE TO BE FILLED IN
         allfound = (found_molinp and found_dalinp and found_nb and found_no)
         if(allfound):
           break


      ##########################################
      #PARSING THE FILE IN A FIRST ROUND IS DONE
      ##########################################


      if(found_nb and found_no):
        print "ATTENTION: CALCULATION OF N_VIRTUAL = N_BASIS - N_OCCUPIED"
        self.nv = self.nb - self.no

      #REMOVE ALL TRAILING BLANK LINES FROM THE INPUT FILES
      while(self.molinp[0] == ''):
        del self.molinp[0]
      while(self.molinp[-1] == ''):
        del self.molinp[-1]
      while(self.dalinp[0] == ''):
        del self.dalinp[0]
      while(self.dalinp[-1] == ''):
        del self.dalinp[-1]

      #FIND MORE SPECIFIC INFORMATION ACCORDING TO THE JOB STRING
      
      #Read DEC fragments from DEC calculation
      if("DEC"==self.calctype[0] and not "MP2DEBUG" in self.calctype):
        self.decinfo.get_dec_info(self.lines,self.calctype[1],True)

      #Read CC information
      if("CC"==self.calctype[0]):
        self.ccinfo.get_cc_info(self.lines,self.calctype[1])

        #Read fragment info if availale
        if self.decinfo.enable_fragread:
           self.decinfo.get_dec_info(self.lines,self.calctype[1],False)


   ############################################################
   ###########         DEC SPECIFIC FUNCTIONS       ###########
   ############################################################
   #GET FRAGMENT ENERGIES FROM FULL CALCULATION
   def get_fraginfo_from_full(self):
      if(("DEC"==self.calctype[0] and "MP2DEBUG" in self.calctype)or self.decinfo.enable_fragread):
        print "reading fragment info"
        self.decinfo.get_dec_info(self.lines,self.calctype[1],False)
      else:
        print "ERROR(get_frag_from_full): cannot be performed for this type of calculation"
   #REFERENCE TO FULL CALCULATION
   def ref_to_full(self):
      print "NOT YET IMPLEMENTED"
   #WRITE FRAGMENT ENERGIES
   def write_sfrag_occ_energies(self,outstream,whichen,printabs):
      energ_id1 = 0
      founden = False
      for i in range(len(self.decinfo.ecorrtype)):
        if(self.decinfo.ecorrtype[i] == whichen):
          energ_id = i
          founden = True
      
      if(founden):
        for i in range(self.decinfo.sfragjobs):
          if(printabs):
            outstream.write(str(i) 
             +" "+ str('%.3e' % (abs(self.decinfo.sfrags[i].ecorrocc[energ_id])))
             +" "+ str(self.decinfo.sfrags[i].fragid) + "\n")
          else:
            outstream.write(str(i) 
             +" "+ str('%.3e' % (self.decinfo.sfrags[i].ecorrocc[energ_id]))
             +" "+ str(self.decinfo.sfrags[i].fragid) + "\n")
      else:
        print "ERROR(write_sfrag_occ_energies) called with the arguments"
        print self.decinfo.ecorrtype,"whichen=",whichen,"\tprintabs=",printabs
        print "THE SPECIFIED ENERGY COULD NOT BE FOUND IN THE CALCULATION"

   def write_pfrag_occ_energies(self,outstream,whichen,printabs):
      energ_id1 = 0
      founden = False
      for i in range(len(self.decinfo.ecorrtype)):
        if(self.decinfo.ecorrtype[i] == whichen):
          energ_id = i
          founden = True
    
      if(founden):
        for i in range(self.decinfo.pfragjobs):
          if(printabs):
            outstream.write(str(self.decinfo.pfrags[i].dist) 
            +" "+ str('%.3e' % (abs(self.decinfo.pfrags[i].ecorrocc[energ_id])))
            +" "+ str(self.decinfo.pfrags[i].fragid)+" "+ str(self.decinfo.pfrags[i].fragpid) + "\n")
          else:
            outstream.write(str(self.decinfo.pfrags[i].dist) 
            +" "+ str('%.3e' % (self.decinfo.pfrags[i].ecorrocc[energ_id]))
            +" "+ str(self.decinfo.pfrags[i].fragid)+" "+ str(self.decinfo.pfrags[i].fragpid) + "\n")
      else:
        print "ERROR(write_pfrag_occ_energies) called with the arguments"
        print self.decinfo.ecorrtype,"whichen=",whichen,"\tprintabs=",printabs
        print "THE SPECIFIED ENERGY COULD NOT BE FOUND IN THE CALCULATION"


   #COMPARE FRAGMENT ENERGIES
   def compare_sfrag_occ_energies(self,another,outstream,whichen,printabs,printstuff):
      energ_id1 = 0
      energ_id2 = 0
      founden1 = False
      founden2 = False
      accerr = 0.0

      for i in range(len(self.decinfo.ecorrtype)):
        if(self.decinfo.ecorrtype[i] == whichen):
          energ_id1 = i 
          founden1 = True
      for i in range(len(another.decinfo.ecorrtype)):
        if(another.decinfo.ecorrtype[i] == whichen):
          energ_id2 = i 
          founden2  = True

      if(founden1 and founden2):
        for i in range(self.decinfo.sfragjobs):
          if(printabs):
            outstream.write(str(i) 
             +" "+ str('%.3e' % (abs(self.decinfo.sfrags[i].ecorrocc[energ_id1] - another.decinfo.sfrags[i].ecorrocc[energ_id2])))
             +" "+ str(self.decinfo.sfrags[i].fragid) + "\n")
            accerr+= abs(self.decinfo.sfrags[i].ecorrocc[energ_id1] - another.decinfo.sfrags[i].ecorrocc[energ_id2])
 
          else:
            outstream.write(str(i) 
             +" "+ str('%.3e' % (self.decinfo.sfrags[i].ecorrocc[energ_id1] - another.decinfo.sfrags[i].ecorrocc[energ_id2]))
             +" "+ str(self.decinfo.sfrags[i].fragid) + "\n")
            accerr+= (self.decinfo.sfrags[i].ecorrocc[energ_id1] - another.decinfo.sfrags[i].ecorrocc[energ_id2])
        if(printstuff):
          if(printabs):
            print "ACCUMULATED ABSOLUTES "+whichen+" SINGLE FRAGMENT ERROR", accerr
          else:
            print "ACCUMULATED "+whichen+" SINGLE FRAGMENT ERROR", accerr
      else:
        print "ERROR(compare_sfrag_occ_energies) called with the arguments"
        print "self=",self.decinfo.ecorrtype,"\tanother=",another.decinfo.ecorrtype
        print "whichen=",whichen,"\tprintabs=",printabs,"\tprintstuff=",printstuff
        print "THE SPECIFIED ENERGY COULD NOT BE FOUND IN THE CALCULATION"
      return accerr

   def compare_pfrag_occ_energies(self,another,outstream,whichen,printabs,printstuff):
      energ_id1 = 0
      energ_id2 = 0
      founden1 = False
      founden2 = False
      accerr = 0.0

      for i in range(len(self.decinfo.ecorrtype)):
        if(self.decinfo.ecorrtype[i] == whichen):
          energ_id1 = i 
          founden1 = True
      for i in range(len(another.decinfo.ecorrtype)):
        if(another.decinfo.ecorrtype[i] == whichen):
          energ_id2 = i 
          founden2  = True

      if(founden1 and founden2):
        for i in range(self.decinfo.pfragjobs):
          if(printabs):
            outstream.write(str(self.decinfo.pfrags[i].dist) 
            +" "+ str('%.3e' % (abs(self.decinfo.pfrags[i].ecorrocc[energ_id1] - another.decinfo.pfrags[i].ecorrocc[energ_id2])))
            +" "+ str(self.decinfo.pfrags[i].fragid)+" "+ str(self.decinfo.pfrags[i].fragpid) + "\n")
            accerr+= abs(self.decinfo.pfrags[i].ecorrocc[energ_id1] - another.decinfo.pfrags[i].ecorrocc[energ_id2])
          else:
            outstream.write(str(self.decinfo.pfrags[i].dist) 
            +" "+ str('%.3e' % (self.decinfo.pfrags[i].ecorrocc[energ_id1] - another.decinfo.pfrags[i].ecorrocc[energ_id2]))
            +" "+ str(self.decinfo.pfrags[i].fragid)+" "+ str(self.decinfo.pfrags[i].fragpid) + "\n")
            accerr+= (self.decinfo.pfrags[i].ecorrocc[energ_id1] - another.decinfo.pfrags[i].ecorrocc[energ_id2])
        if(printstuff):
          if(printabs):
            print "ACCUMULATED ABSOLUTES "+whichen+" PAIR FRAGMENT ERROR", accerr
          else:
            print "ACCUMULATED "+whichen+" PAIR FRAGMENT ERROR", accerr
      else:
        print "ERROR(compare_pfrag_occ_energies) called with the arguments"
        print "self=",self.decinfo.ecorrtype,"\tanother=",another.decinfo.ecorrtype
        print "whichen=",whichen,"\tprintabs=",printabs,"\tprintstuff=",printstuff
        print "THE SPECIFIED ENERGY COULD NOT BE FOUND IN THE CALCULATION"
      return accerr

   def compare_pfrag_energies(self,another,outstream):
      print "NOT YET IMPLEMENTED"
      
