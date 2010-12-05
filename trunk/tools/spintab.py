#!/usr/bin/python
#
#
#Processing dalton.out-files for calculation of NMR spin-spin couplings.
#Written by Ola Berg Lutnaes 1999.
#
#Usage: spintab options file1 file2...
#
import sys, re,string,os
#
#DECLARATIONS
#
geometry=0
ac=0
ps=0
totalonly=0
equals=1
decimals=2
latex=1
latex2=0
plot=0
cont=0
fccont=0
dsocont=0
psocont=0
sdcont=0
inputfiles=[]
tot=[]
sd=[]
pso=[]
dso=[]
fc=[]
geom=[]
basissets=[]
couplingatoms=[]
numberofcouplings=0
outputfiles=[]
#
#READING COMMANDLINE
#
try:
	option=sys.argv[1]
except:
	print "Usage: spintab.py options inputfiles\nType \"spintab.py --help\" for help."
	sys.exit(1)
while len(sys.argv)>1:
	option=sys.argv[1]
	del sys.argv[1]
	if option=="--help":
		print """TAKES DALTON OUTPUT FOR CALCULATION OF NMR SPIN-SPIN COUPLING CONSTANTS
AND MAKES TABLE FILE FOR LATEX CONTAINING THE COUPLING CONSTANTS.

Latex file is named "table.tex" 

The program also checks for and removes equal coupling constants
(it is not checked if this is because of symmetry or not)
to reduce unnecessary output.

The program can read several outputfiles at once, collecting the
results in one table, provided that the outputfiles are from
similar calculations (for example only changing the basis set).
IMPORTANT: The program assumes in this case that the coupling
constants are ordered equally in the different outputfiles.

SPINTAB OPTIONS:
-d n	controls number of decimals in the latextable. 
	Value of n can be from 0 to 4. Default is 2.
-w	print latex table i "wide" format (can only
        process one file at a time with this option)
-e	don't remove equal coupling constants
-t	give only total coupling constants in output
-g	print more info on basisset, molecular geometry and 
	energy together with latex table
-b n n n n	the -b option makes the program sum contributions
		to the coupling constants from different 
		calculations. -b should be followed by four 
		numbers indicating from what dalton outputfile the 
		different contributions (fc, dso, pso and sd 
		respectively) should be collected.
		Example: spintab.py -b 1 2 2 2 file1.out file2.out
		sums fc-contribution from file1.out with dso- pso-
		and sd-contributions from file2.out.
***
The next options are more machine dependent, but some modification of
the script might make it work on your computer.
***
-f 	make postscript of latex table on the fly and show
	file in ghostview
-p	makes plot using gnuplot. Useful for basisset 
	convergence testing. Makes plots for total coupling
	constant and all contributions. Gnuplot data are 
	stored in file "gnuplot.data". Makes 
	psotscript(.eps-files. All files are merged in to 
	an overview-file "merged.eps", using epsmerge. For 
	printing, use "merged.pdf".
	The name of the dalton dal-file will be used to 
	identify the molecule and used in plot titles etc.
-pac	same as -p, but plots all contributions to the coupling 
	constant in the same diagram
-l	don't make latex table.
***

TROUBLESHOOTING
---------------
1. "spintab.py: Command not found." make sure the script look for
    python in the right place: The result of 'which python' on the
    commandline shall replace '/usr/bin/python' in the first line of
    spintab.py.
"""
		sys.exit(1)	
	if option=="-d":
		decimals=sys.argv[1]
		del sys.argv[1]
	elif option=="-l":
		latex=0
	elif option =="-p":
		plot=1
	elif option=="-e":
		equals=0
	elif option=="-t":
		totalonly=1
	elif option=="-f":
		ps=1
	elif option=="-g":
		geometry=1
	elif option=="-pac":
		plot=1
		ac=1
	elif option=="-w":
		latex=0
		latex2=1
	elif option=="-b":
		cont=1
		fccont=sys.argv[1]
		del sys.argv[1]
		dsocont=sys.argv[1]
		del sys.argv[1]
		psocont=sys.argv[1]
		del sys.argv[1]
		sdcont=sys.argv[1]
		del sys.argv[1]
	else:
		inputfiles.append(option)
		if not os.path.isfile(option):
			print "File "+option+" not found"
			sys.exit(1)
decimals=int(decimals)
if not inputfiles:
	print "No inputfiles specified."
	sys.exit()
##




#SUBROUTINES


#METHOD READING FILES ONE BY ONE AND STORING DATA

#Data stored:
#basis set used
#total coupling constants and induvidual contributions
#names of the two atoms involved in the coupling

def readfiles(filelist):
	numberofcouplings=0
	couplingatoms1=[]
	basissets1=[]
	tot1=[]
	fc1=[]
	pso1=[]
	dso1=[]
	sd1=[]
	geom=[]
	files=filelist[:]
	pattern1="  Basis set used\.*"				#SEARCHPATTERNS
	pattern2=r"                   ! ABACUS - Final spin-spin couplings !\.*"
	pattern3="          Indirect spin-spin coupling\.*"
	pattern4="  Isotropic coupling       \.*"
	pattern5="  Isotropic DSO contributio\.*"
	pattern6="  Isotropic PSO contributio\.*"
	pattern7="  Isotropic SD contribution\.*"
	pattern8="  Isotropic FC contributi\.*"
	pattern9=r"\.*  DALTON - An electronic structure program\.*"
	pattern10="  Threshold fo\.*"
	pattern11="   Interatomic separatio\.*"
	pattern12="  Symmetry Orb\.*"
	pattern13="\.*FINAL RESULTS FROM ABACUS\.*"
	pattern14="     Total energy     \.*"
	numberoffiles=len(files)
	while len(files)>0:
		totals=[]
		fcs=[]
		psos=[]
		dsos=[]
		sds=[]
		couplingatomses=[]
		infile=open(files[0],'r')
		del files[0]
		lines = infile.readlines()
		test=0
		for a in lines:				#TESTING INPUTFILES
			match = re.search(pattern9,a)
			if match:
				test=1
				break
		if test==0:
			print "One of the inputfiles is not a dalton output file"
			sys.exit(1)
		n=0
		while lines[n]:					#FINDING BASISSETS
			match=re.search(pattern1,lines[n])
			if match:
				parts=string.split(lines[n],"\"")
				basissets1.append(parts[1])
				while lines[n]:
					geom.append(lines[n])
					match=re.search(pattern10,lines[n])
					if match:
						break
					n=n+1
				break
			n=n+1
#		while lines[n]:					#STORING GEOMETRY DATA
#			n=n+1
#			match =re.search(pattern11,lines[n])
#			if match:
#				while lines[n]:
#					match=re.search(pattern12,lines[n])
#					if match:
#						break
#						break
#					geom.append(lines[n])
#					n=n+1
#				break
		while lines[n]:					#STORING ENERGY
			match=re.search(pattern14,lines[n])
			if match:
				geom.append(lines[n])
				break
			n=n+1
		geom.append("\n\n"+r"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"+"\n\n")
		while lines[n]:					#FAST FORWARD
			match=re.search(pattern2,lines[n])
			if match:
				break
			n=n+1
		while lines[n]:					#STORING COUPLING CONSTANTS
			match=re.search(pattern3,lines[n])
			if match:
				parts=string.split(lines[n],"between")
				parts=string.split(parts[1],":")
				couplingatomses.append(parts[0])
				numberofcouplings=numberofcouplings+1
			match=re.search(pattern4,lines[n])	
			if match:
				parts=string.split(lines[n],":")
				parts=string.split(parts[1],"Hz")
				parts=re.sub(" ","",parts[0])
				totals.append(parts)
			match=re.search(pattern5,lines[n])
			if match:
				parts=string.split(lines[n],":")
				parts=string.split(parts[1],"Hz")
				parts=re.sub(" ","",parts[0])
				dsos.append(parts)
			match=re.search(pattern6,lines[n])
			if match:
				parts=string.split(lines[n],":")
				parts=string.split(parts[1],"Hz")
				parts=re.sub(" ","",parts[0])
				psos.append(parts)
			match=re.search(pattern7,lines[n])
			if match:
				parts=string.split(lines[n],":")
				parts=string.split(parts[1],"Hz")
				parts=re.sub(" ","",parts[0])
				sds.append(parts)
			match=re.search(pattern8,lines[n])
			if match:
				parts=string.split(lines[n],":")
				parts=string.split(parts[1],"Hz")
				parts=re.sub(" ","",parts[0])
				fcs.append(parts)
			n=n+1
			try:
				lines[n]
			except:
				break
		tot1.append(totals)
		sd1.append(sds)
		dso1.append(dsos)
		pso1.append(psos)
		fc1.append(fcs)
		couplingatoms1.append(couplingatomses)	
		infile.close()
	numberofcouplings=numberofcouplings/numberoffiles
	return(tot1,fc1,dso1,pso1,sd1,couplingatoms1,basissets1,numberofcouplings,geom)
##

#TESTING INPUT TO BE ADEQUATE
#IF COUPLINGATOMS HAVE DIFFERENT NAME WARNING WILL BE PRINTED

def checkinput(couplingatoms):			
	if len(couplingatoms)>1:
		for c in couplingatoms:
			for d in couplingatoms:
				if c!=d:
					print "WARNING!! Your inputfiles may contain different sets of coupling constants."
##


#METHOD ADJUSTING DECIMALS IN OUTPUT

def decimalprint(couplinglist,decimals,n):
	tablestring=""
	for a in couplinglist:
		if decimals==0:
			tablestring=tablestring+"&"+"%.0f" % float(a[n])
		elif decimals==1:
			tablestring=tablestring+"&"+"%.1f" % float(a[n])
		elif decimals==2:
			tablestring=tablestring+"&"+"%.2f" % float(a[n])
		elif decimals==3:
			tablestring=tablestring+"&"+"%.3f" % float(a[n])
		elif decimals==4:
			tablestring=tablestring+"&"+"%.4f" % float(a[n])
		else:
			print "Decimal error"
			sys.exit(3)
	return tablestring
##


# METHOD THAT DELETES EQUAL COUPLING CONSTANTS


def sortcouplings(tot,fc,dso,pso,sd,couplingatoms,numberofcouplings):
	for a in tot[0]:
		for b in range (0,len(tot[0]),1):
			list = tot[0][:]
			ai=tot[0].index(a)
			if a==list[b] and ai != b:
				for n in tot:
					del n[b]
				for n in fc:
					del n[b]
				for n in pso:
					del n[b]
				for n in dso:
					del n[b]
				for n in sd:
					del n[b]
				for n in couplingatoms:
					del n[b]
				numberofcouplings=numberofcouplings-1
				break
	return(tot,fc,dso,pso,sd,couplingatoms,numberofcouplings)
##


#MERGING DIFFERENT BASISSET CALCULATIONS INTO ONE TABLE WHERE DIFFERENT CONTRIBUTIONS HAVE BEEN CALCULATED WITH DIFFERENT BASISSETS
totalitet=[]
def adjust(fccont,dsocont,psocont,sdcont,fc,dso,pso,sd):
	tot=[]
	total=[]
	a=[]
	b=[]
	c=[]
	d=[]
	fccont=int(fccont)-1
	dsocont=int(dsocont)-1
	psocont=int(psocont)-1
	sdcont=int(sdcont)-1
	for n in range (0,len(fc[fccont]),1):
		totale=float(fc[fccont][n])+float(dso[dsocont][n])+float(pso[psocont][n])+float(sd[sdcont][n])
		total.append(totale)
	tot.append(total)
	bas="?"
	a.append(fc[fccont])
	b.append(dso[dsocont])
	c.append(pso[psocont])
	d.append(sd[sdcont])
	return (tot,a,b,c,d,bas)
##

#WRITING TABLE FILE FOR LATEX

#Writing header

def table(totals,fcs,dsos,psos,sds,basissets,couplingatoms,decimals,inputfiles,totalonly,geom,geometry):
	ofile=open("table.tex",'w')
	header=r"Inputfiles: "
	for n in inputfiles:
		n=re.sub("_",r"\_",n)
		header=header+n+" "
	header=header+r"\\"+"\n"+r"\begin{tabular}{lc"
	for c in basissets:
		header=header+"c"
	header = header+r"}"+"\n"+r"Coupling&Contribution"
	for c in range(0,len(basissets),1):
		header=header+r"&"+basissets[c]
	header=header+r"\\ \hline"+"\n"
	ofile.write(header)

#Writing coupling constants

	marg="                  "
	n=0
	pairs=couplingatoms[0]
	while n<len(couplingatoms[0]):
		if totalonly==1:
			ofile.write(pairs[n]+r"&Total")
			number=decimalprint(totals,decimals,n)
			ofile.write(number)	
		else:
			ofile.write(marg+r"&Total")
			number=decimalprint(totals,decimals,n)
			ofile.write(number)
			ofile.write(r"\\"+"\n"+marg+"&FC   ")
			number=decimalprint(fcs,decimals,n)
			ofile.write(number)
			pairs=couplingatoms[0]
			ofile.write(r"\\"+"\n"+pairs[n]+"&DSO ")
			number=decimalprint(dsos,decimals,n)
			ofile.write(number)
			ofile.write(r"\\"+"\n"+marg+"&PSO   ")
			number=decimalprint(psos,decimals,n)
			ofile.write(number)
			ofile.write(r"\\"+"\n"+marg+"&SD   ")
			number=decimalprint(sds,decimals,n)
			ofile.write(number)	
		ofile.write(r"\\"+"\n")
		n=n+1
	ofile.write(r"\end{tabular}")
#Writing geometry and other data
	if geometry==1:
		ofile.write("\n"+r"\begin{verbatim}"+"\n")
		geometrytext=string.join(geom,"")
		ofile.write(geometrytext)
		ofile.write(r"\end{verbatim}")
	ofile.close()
	print "\'table.tex\' complete"
##


#WRITING TABLE FILE FOR LATEX IN WIDE FORMAT

#Writing header

def table2(totals,fcs,dsos,psos,sds,basissets,couplingatoms,decimals,inputfiles,totalonly,geom,geometry):
	couplingatoms=couplingatoms[0]
	ofile=open("table.tex",'w')
	header=r"Inputfiles: "
	for n in inputfiles:
		n=re.sub("_",r"\_",n)
		header=header+n+" "
	header=header+r"\\"+"\n"+r"\begin{tabular}{llrrrrr} \hline\hline"+"\n"
	header=header+r"Molekyl&Kobling&FC&DSO&PSO&SD&Total \\ \hline"+"\n"
	ofile.write(header)

#Writing coupling constants
	b=0
	while b<len(totals[0]):
		number=decimalprint(totals,decimals,b)
#		ofile.write(number)	
		line=r"&"+str(couplingatoms[b])+decimalprint(fcs,decimals,b)+decimalprint(dsos,decimals,b)+decimalprint(psos,decimals,b)+decimalprint(sds,decimals,b)+decimalprint(totals,decimals,b)+"\n"
#		line=str(couplingatoms[b])+r"&"+str(decimalprint(fcs,decimals,b))+r"&"+str(dsos[b])+r"&"+str(psos[b])+r"&"+str(sds[b])+r"&"+str(totals[b])+r"\\" +"\n"
		ofile.write(line)
		b=b+1
	ofile.write(r"\hline\hline"+"\n"+r"\end{tabular}")
#Writing geometry and other data
	if geometry==1:
		ofile.write("\n"+r"\begin{verbatim}"+"\n")
		geometrytext=string.join(geom,"")
		ofile.write(geometrytext)
		ofile.write(r"\end{verbatim}")
	ofile.close()
	print "\'table.tex\' complete"
##


#DATAFILE FOR GNUPLOT

#Header comment

def gnuplotdata(basissets,inputfiles,couplingatoms,tot,fc,dso,pso,sd):
	ofile=open("gnuplot.data",'w')
	header=r"#Datafile for GNUPLOT, spin-spin coupling constants"+"\n"+r"#basissets: "
	for n in basissets:
		header=header+n+" "
	header=header+"\n"+r"#inputfiles:"
	for n in inputfiles:
		header= header+n+"   "
	header=header+"\n"+r"#coupling constant:"
	for n in couplingatoms[0]:
		header=header+n+r"  ####  " 
	header=header+"\n"
	ofile.write(header)


#Printing data


	for n in range (0,len(inputfiles),1):
		linenumber=n+1
		ofile.write(str(linenumber)+"	")
		total1=tot[n]
		dso1=dso[n]
		pso1=pso[n]
		sd1=sd[n]
		fc1=fc[n]
		for a in range(0,len(total1),1):
			ofile.write(total1[a]+"	")
			ofile.write(fc1[a]+"	")
			ofile.write(pso1[a]+"	")
			ofile.write(dso1[a]+"	")
			ofile.write(sd1[a]+"	")
		ofile.write("\n")
	ofile.close()
	print "\'gnuplot.data\' complete"
##


#MAKING GNUPLOTCOMMAND FOR INDUVIDUAL PLOTS

def gnuplotcommand(inputfiles,couplingatoms,basissets):

	molecule=string.split(inputfiles[0],"_")
	molecule=molecule[0]
	outputfiles=[]
	gnuplotcommands=[]
	for n in range(0,len(couplingatoms[0]),1):
		a=1
		while a<6:
			plotcolumn=n*5+a+1
			coupling=couplingatoms[0]
			title="Spin-spin coupling "+coupling[n]+" in "+molecule+" using basissets"
			outputname=re.sub(" ","",coupling[n])+molecule+".eps"		
			for b in basissets:
				title= title+b+" "
			if a==1:
				title=title+"Total coupling"
				outputname="tot"+outputname
			if a==2:
				title=title+"FC contribution"
				outputname="fc"+outputname
			if a==3:
				title=title+"DSO contribution"
				outputname="pso"+outputname
			if a==4:
				outputname="dso"+outputname
				title=title+"PSO contribution"
			if a==5:
				outputname="sd"+outputname
				title=title+"SD contribution"
			outputfiles.append(outputname)
			gnuplotcmd="""set terminal postscript eps 14
set autoscale  
set nokey
set xlabel 'Økende basissett'
set ylabel 'Hz'
set title """
			gnuplotcmd=gnuplotcmd+"\'"+title+r"'"+"\nset output "+r"'"+outputname+r"'"+"\nplot "+r"'"+"gnuplot.data"+"'"+" using 1"+r":"+str(plotcolumn)+" w lp\nquit\n" 
			gnuplotcommands.append(gnuplotcmd)
			a=a+1
	return (gnuplotcommands,outputfiles)
##


#MAKING GNUPLOTCOMMAND FOR ALL CONTRIBUTIONS IN ONE PLOT
def gnuplotcommandac(inputfiles,couplingatoms,basissets):

	molecule=string.split(inputfiles[0],"_")
	molecule=molecule[0]
	outputfiles=[]
	gnuplotcommands=[]
	for n in range(0,len(couplingatoms[0]),1):
		theplotcolumns=[]
		a=1
		while a<6:
			plotcolumn=5*n+a+1
			theplotcolumns.append(plotcolumn)
			a=a+1
		coupling=couplingatoms[0]
		title="Spin-spin coupling "+coupling[n]+" in "+molecule+" using basissets"
		outputname=re.sub(" ","",coupling[n])+molecule+".eps"
		for b in basissets:
			title= title+b+" "
		title=title+"\\\nAll contributions"
		outputname="ac"+outputname
		outputfiles.append(outputname)
		gnuplotcmd="""set terminal postscript eps 14
set autoscale  
set xlabel 'Økende basissett'
set ylabel 'Hz'
set title """
		gnuplotcmd=gnuplotcmd+"\'"+title+r"'"+"\nset output "+r"'"+outputname+r"'"+"\nplot "
		for n in range(len(theplotcolumns)):
			newplot=r"'"+"gnuplot.data"+"'"+" using 1"+r":"+str(theplotcolumns[n])
			if n==0:
				newplot=newplot+" title \'Total\' w lp"
			if n==1:
				newplot=newplot+" title \'FC\' w lp"
			if n==2:
				newplot=newplot+" title \'DSO\' w lp"
			if n==3:
				newplot=newplot+" title \'PSO\' w lp"
			if n==4:
				newplot=newplot+" title \'SD\' w lp"
			if n<len(theplotcolumns)-1:
				newplot=newplot+","
			gnuplotcmd=gnuplotcmd+newplot
		gnuplotcmd=gnuplotcmd+"\nquit\n"
		gnuplotcommands.append(gnuplotcmd)
	return (gnuplotcommands,outputfiles)

#EXECUTING GNUPLOT
#
def gnuplotexecute(gnuplotcommands):
	for n in gnuplotcommands:
		ofilename=string.split(n,"set output \'")
		ofilename=string.split(ofilename[1],"\'")
		ofilename=ofilename[0]
		gnuplot=os.popen("gnuplot -persist",'w')
		gnuplot.write(n)
		gnuplot.close()
		print "\'"+ofilename+"\' complete"
		
##

#MERGING PLOTFILES INTO ONE USING SCRIPT EPSMERGE AND PS2PDF


def mergeplots(outputfiles):
	cmd="epsmerge -o merged.eps"
	for n in outputfiles:
		cmd=cmd+" "+n
	os.system(cmd)
	print "\'merged.eps\' complete"
	os.system("ps2pdf merged.eps merged.pdf")
	print "\'merged.pdf\' complete"
##

#MAKING POSTSCRIPTFILE FROM LATEXTABLE

def postscript():
	ofile=open("maintable.tex",'w')
	start="""\documentclass[12pt,norsk,a4paper]{report}
\usepackage{a4}
\usepackage{graphics}
\usepackage[norsk]{babel}
\usepackage{babel,varioref}
\usepackage{psfig}
\usepackage{epic}
\usepackage[latin1]{inputenc}
\usepackage{longtable}
\usepackage{dcolumn}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{supertabular}
\usepackage{multicol,ftnright} 
\usepackage{wrapfig}
\usepackage{cite}
\usepackage{changebar}
\usepackage{pstricks}
\usepackage{graphics}
\usepackage{tabularx}
\usepackage{multirow}
\usepackage{rotating}
\\begin{document}
\include{table}
\end{document}
"""
	ofile.write(start)
	ofile.close()
	cmd="""rm table.ps
latex maintable.tex
dvips -f maintable.dvi > table.ps
ghostview table.ps &
"""
	os.system(cmd)
	print"\'maintable.tex\' complete"
	print"\'table.ps\' complete"
##


###SCRIPT CONTROL

#Reading files
[tot,fc,dso,pso,sd,couplingatoms,basissets,numberofcouplings,geom]=readfiles(inputfiles)

print basissets
#Removing equal constants
if equals==1:
	[tot,fc,dso,pso,sd,couplingatoms,numberofcouplings]=sortcouplings(tot,fc,dso,pso,sd,couplingatoms,numberofcouplings)

#ADJUSTING TABLE DATA FOR  CONTRIBUTIONS TO SPINSPIN FROM DIFFERENT BASIS SETS


if cont==1:
#	inputfiles???
	[tot,fc,dso,pso,sd,basissets]=adjust(fccont,dsocont,psocont,sdcont,fc,dso,pso,sd)


#WRITING OUTPUT

#LATEX TABLE

if latex==1:
	table(tot,fc,dso,pso,sd,basissets,couplingatoms,decimals,inputfiles,totalonly,geom,geometry)		

#LATEX TABLE FORMAT  2 (FOR LARGE MOLECULES WITH MANY COUPLINGS)

if latex2==1:
	table2(tot,fc,dso,pso,sd,basissets,couplingatoms,decimals,inputfiles,totalonly,geom,geometry)



#PLOTTING

if plot==1:
	gnuplotdata(basissets,inputfiles,couplingatoms,tot,fc,dso,pso,sd)
	if ac==1:
		[gnuplotcommands,outputfiles]=gnuplotcommandac(inputfiles,couplingatoms,basissets)
	else:
		[gnuplotcommands,outputfiles]=gnuplotcommand(inputfiles,couplingatoms,basissets)
	gnuplotexecute(gnuplotcommands)

#TRANSFORMING FILE FORMATS

	mergeplots(outputfiles)
if ps==1:
	postscript()

