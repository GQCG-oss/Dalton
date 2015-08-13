import numpy as np
import matplotlib.ticker
try:
   import matplotlib.pyplot as plt
   imported = True
except:
   imported = False

def open_fig_container():
   if imported:
     fig = plt.figure()
     return fig
   else:
     print "This function is not available without a display"

#PLOT PAIR ENERGIES IS A FUNCTION OF THE decinfo STRUCTURE
#NECESSARY INPUTS ARE THE STRUCTURE WHICH CALLS THE FUNCTION AND FIG OBTAINED   
def plot_pair_energies(self,fig,ecorrtype="oMP2",title="DEFAULT TITLE",to_plot=0,color='b',marker="s",label="series",to_diff=0,xa=[0,0],ya=[0,0]):
   if(to_plot==0):
     ax1 = fig.add_subplot(111)
   else:
     ax1 = to_plot

   x=[]
   y=[]
   
   #DEFINING WHICH KIND OF DATA TO USE
   ov  = "o"
   ect = 0
   if ecorrtype=="oCCSD" or ecorrtype=="oMP2" or ecorrtype=="oRIMP2":
      ov  = "o"
      ect = 0
   elif ecorrtype=="vCCSD" or ecorrtype=="vMP2" or ecorrtype=="vRIMP2":
      ov  = "v"
      ect = 0
   elif ecorrtype=="o(T)":
      ov  = "o"
      ect = 1
   elif ecorrtype=="v(T)":
      ov  = "v"
      ect = 1
   elif ecorrtype=="oEsti":
      ov = "o"
      ect = 3
   else: 
      print "ERROR(plot_pair_energies):invalid choice of ecorrtype"
      exit()
   
   ymax = 0.0
   atoms = [0,0]
   #READING THE DATA FROM THE STRUCTURE
   if(to_diff==0):
      for i in range(self.pfragjobs):
        x.append(self.pfrags[i].dist)
        if ov == "o":
          if ecorrtype=="o(T)":
             if abs(self.pfrags[i].ecorrocc[ect]) > ymax:
                atoms = [self.pfrags[i].fragid,self.pfrags[i].fragpid]
          y.append(abs(self.pfrags[i].ecorrocc[ect]))
          ymax = max(ymax,y[-1])
        elif ov =="v":
          y.append(abs(self.pfrags[i].ecorrvirt[ect]))
   else:
      #READING THE DATA FROM THE STRUCTURE
      for i in range(self.pfragjobs):
        x.append(self.pfrags[i].dist)
        if ov == "o":
          y.append(abs(self.pfrags[i].ecorrocc[ect]-to_diff.pfrags[i].ecorrocc[ect]))
        elif ov =="v":
          y.append(abs(self.pfrags[i].ecorrvirt[ect]-to_diff.pfrags[i].ecorrvirt[ect]))


   #########################
   #SETTING AXIS PROPERTIES#
   #########################

   print "max pair was",ymax,atoms
   if(to_plot==0):
     xmin = min(x)
     xmax = max(x) 
     ymin = min(y) 
     ymax = max(y)
   else:
     xmin, xmax, ymin, ymax = ax1.axis()
     xmin = min(min(x),xmin)
     xmax = max(max(x),xmax)
     ymin = min(min(y),ymin)
     ymax = max(max(y),ymax)

   ax1.set_title(title)

   #y-axis properties
   ax1.set_autoscaley_on(False)
   ax1.set_yscale('log')
   ax1.set_ylim([min(ymin,0.9*ymin),max(ymax,1.1*ymax)])
   ax1.set_ylabel("Pair energy [Eh]")

   #x-axis properties
   ax1.set_autoscalex_on(False)
   ax1.set_xscale('log')
   xmin = min(xmin,0.9*xmin)
   xmax = max(xmax,1.1*xmax)
   ax1.set_xlim([xmin,xmax])
   ax1.set_xlabel("Pair distance [AA]")
   ax1.set_xticks(range(int(xmin)+1,int(xmax)+1))
   ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

   #plot data as scatter plot
   ax1.scatter(x, y, s=20, c=color, marker=marker, label = label)
   ax1.legend(loc=3)

   return ax1

# The plot_SF_energy_errors function plot the atomic fragment energy errors
# of a dec calculation compared to a (full) reference job.
#
# Necessary inputs are the structure from a reference job "self0",
# the structure from an actual (dec) calculation "self1" and fig obtained:
def plot_SF_energy_errors(self, info, fig, ecorrtype="oMP2", title="AF energy errors", 
    color='b', marker="o", label="series", to_plot=0):

   if(to_plot==0):
     ax1 = fig.add_subplot(111)
   else:
     ax1 = to_plot

   # DEFINING WHICH KIND OF DATA TO USE:
   ov  = "o"
   ect = 0
   if ecorrtype=="oCCSD" or ecorrtype=="oMP2" or ecorrtype=="oRIMP2":
      ov  = "o"
      ect = 0
   elif ecorrtype=="vCCSD" or ecorrtype=="vMP2" or ecorrtype=="vRIMP2":
      ov  = "v"
      ect = 0
   elif ecorrtype=="o(T)":
      ov  = "o"
      ect = 1
   elif ecorrtype=="v(T)":
      ov  = "v"
      ect = 1
   else: 
      print "ERROR(plot_SF_energy_error):invalid choice of ecorrtype"
      exit()
   
   ect = ect - 1
   # read reference AF energies:
   natoms=self.sfragjobs
   AF_full=np.zeros((natoms))
   xaxis=np.zeros((natoms),dtype=int)
   for i in range(natoms):
     xaxis[i] = i+1
     AF_full[i] = self.sfrags[i].ecorrocc[ect]


   # read actual AF energies and compute errors:
   natdec=info.sfragjobs
   # check:
   if (natdec!=natoms):
     print "natoms from ref:",natoms
     print "natoms from dec:",natdec
     exit("ERROR(plot_SF_energy_error):input structure not compatible")
    
   AF_dec=np.zeros((natoms))
   AF_err=np.zeros((natoms))
   AF_tot_err = 0.0
   for i in range(natoms):
     AF_dec[i] = info.sfrags[i].ecorrocc[ect]
     AF_err[i] = AF_dec[i] - AF_full[i]
     AF_tot_err += AF_err[i]
    
   # print total AF errors:
   print "Total Atomic Fragment Energy error = ",AF_tot_err

   #########################
   #SETTING AXIS PROPERTIES#
   #########################

   if (to_plot==0):
     xmin = min(xaxis)
     xmax = max(xaxis) 
     ymin = min(abs(AF_err)) 
     ymax = max(abs(AF_err))
   else:
     xmin, xmax, ymin, ymax = ax1.axis()

   ax1.set_title(title)

   #y-axis properties
   ax1.set_autoscaley_on(False)
   ax1.set_yscale('log')
   ax1.set_ylim([min(ymin,0.9*ymin),max(ymax,1.1*ymax)])
   ax1.set_ylabel("AF energy errors [Eh]")

   #x-axis properties
   ax1.set_autoscalex_on(False)
   #ax1.set_xscale('log')
   xmin = min(xmin,0.9*xmin)
   xmax = max(xmax,1.1*xmax)
   ax1.set_xlim([xmin,xmax])
   ax1.set_xlabel("Atomic center label")
   ax1.set_xticks(range(int(xmin)+1,int(xmax)+1))
   ax1.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

   # plot Atomic Fragment errors:
   ax1.scatter(xaxis, abs(AF_err), s=20,c=color, marker=marker, label = label)
   ax1.axhline(y=info.fotfloat, color=color)
   ax1.legend()

   return ax1

def show_plots():
   if imported:
     plt.show()
   else:
     print "This function is not available without a display"

def save_plots(name):
   if imported:
     plt.savefig(name)
   else:
     print "This function is not available without a display"

