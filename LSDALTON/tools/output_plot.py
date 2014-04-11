import matplotlib.pyplot as plt
import matplotlib.ticker
def open_fig_container():
   fig = plt.figure()
   return fig

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
   if ecorrtype=="oCCSD" or ecorrtype=="oMP2":
      ov  = "o"
      ect = 1
   elif ecorrtype=="vCCSD" or ecorrtype=="vMP2":
      ov  = "v"
      ect = 1
   elif ecorrtype=="o(T)":
      ov  = "o"
      ect = 2
   elif ecorrtype=="v(T)":
      ov  = "v"
      ect = 2
   else :
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


def show_plots():
   plt.show()
