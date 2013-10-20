import matplotlib.pyplot as plt
def open_fig_container():
   fig = plt.figure()
   return fig
   

def plot_pair_energies(self,fig,ecorrtype="oMP2",title="DEFAULT TITLE",to_plot=0,color='b',marker="s",to_diff=0):
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
      ect = 0
   elif ecorrtype=="vCCSD" or ecorrtype=="vMP2":
      ov  = "v"
      ect = 0
   elif ecorrtype=="o(T)":
      ov  = "o"
      ect = 1
   elif ecorrtype=="v(T)":
      ov  = "v"
      ect = 1
   else :
      print "ERROR(plot_pair_energies):invalid choice of ecorrtype"
      exit()
   
   #READING THE DATA FROM THE STRUCTURE
   if(to_diff==0):
      for i in range(self.pfragjobs):
        x.append(self.pfrags[i].dist)
        if ov == "o":
          y.append(abs(self.pfrags[i].ecorrocc[ect]))
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

   if(to_plot==0):
     xmin = min(x)
     xmax = max(x) 
     ymin = min(y) 
     ymax = max(y)
   else:
     xmin, xmax, ymin, ymax = ax1.axis()

   #SETTING AXIS PROPERTIES
   ax1.set_title(title)
   ax1.set_autoscaley_on(False)
   ax1.set_yscale('log')
   ax1.set_ylim([min(ymin,0.9*ymin),max(ymax,1.1*ymax)])
   ax1.set_autoscalex_on(False)
   ax1.set_xscale('log')
   ax1.set_xlim([min(xmin,0.9*xmin),max(xmax,1.1*xmax)])
   ax1.scatter(x, y, s=20, c=color, marker=marker)

   return ax1


def show_plots():
   plt.show()
