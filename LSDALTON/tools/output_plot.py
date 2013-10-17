import matplotlib.pyplot as plt
def open_fig_container():
   fig = plt.figure()
   return fig
   

def plot_pair_energies(fig,decinfo_struct,ecorrtype="oMP2",title="DEFAULT TITLE",to_plot=0):
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
   for i in range(decinfo_struct.pfragjobs):
     x.append(decinfo_struct.pfrags[i].dist)
     if ov == "o":
       y.append(abs(decinfo_struct.pfrags[i].ecorrocc[ect]))
     elif ov =="v":
       y.append(abs(decinfo_struct.pfrags[i].ecorrvirt[ect]))

   #SETTING AXIS PROPERTIES
   ax1.set_title(title)
   ax1.set_autoscaley_on(False)
   ax1.set_yscale('log')
   ax1.set_ylim([0.5*min(y),5*max(y)])
   ax1.set_autoscalex_on(False)
   ax1.set_xscale('log')
   ax1.set_xlim([0.9*min(x),1.1*max(x)])
   ax1.scatter(x, y, s=20, c='b', marker="s")

   return ax1

def plot_difference_pair_energies(fig,dis1,dis2,ecorrtype="oMP2",title="DEFAULT TITLE",to_plot=0):
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
   for i in range(dis1.pfragjobs):
     x.append(dis1.pfrags[i].dist)
     if ov == "o":
       y.append(abs(dis1.pfrags[i].ecorrocc[ect]-dis2.pfrags[i].ecorrocc[ect]))
     elif ov =="v":
       y.append(abs(dis1.pfrags[i].ecorrvirt[ect]-dis2.pfrags[i].ecorrvirt[ect]))

   #SETTING AXIS PROPERTIES
   ax1.set_title(title)
   ax1.set_autoscaley_on(False)
   ax1.set_yscale('log')
   ax1.set_ylim([0.5*min(y),5*max(y)])
   ax1.set_autoscalex_on(False)
   ax1.set_xscale('log')
   ax1.set_xlim([0.9*min(x),1.1*max(x)])
   ax1.scatter(x, y, s=20, c='r', marker="x")

   return ax1

def show_plots():
   plt.show()
