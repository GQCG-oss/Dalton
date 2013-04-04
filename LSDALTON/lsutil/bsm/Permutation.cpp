#include "Permutation.h"
#include <fstream>
#include <cmath>
#include <iostream>
namespace bsm
{
Permutation::Permutation(const double* xp,const double* yp,const double* zp,
			 const integer* atstart,integer nat,integer maxbsize)
  :natoms(nat),maxblocksize(maxbsize),nblocks(0)
{
try
  {
/* stinne: */
/* printf("Permutation - natoms = %d\n", nat);
printf("sizeof natoms = %d\n", sizeof(nat)); */
    if (nat<=0){throw InputF("Permutation: natoms<=0");}
    if (maxbsize<=0){throw InputF("Permutation: maxblocksize<=0");}
    xpos=new double[natoms];
    ypos=new double[natoms];
    zpos=new double[natoms];
    atomstart=new integer[natoms+1];
    bfperatom=new integer[natoms];
    atomorder=new integer[natoms];
  /* Since nblocks is unknown at this stage     */
  /* its upper bound natoms is used. We add one to deal nicely with
   * one atom systems. DivideSpace::sort adds two blocks at least. */
  natomsperblock=new integer[natoms+1];
  blockstart=new integer[natoms+2];
  blockstart[0]=0;
  integer tmpsize=0;
  integer n=0;
  for (integer i=0;i<natoms;i++)
    {
      xpos[i]=xp[i];
      ypos[i]=yp[i];
      zpos[i]=zp[i];
      atomstart[i]=atstart[i];
      bfperatom[i]=atstart[i+1]-atstart[i];
      if (bfperatom[i]<=0)
	{throw InputF("Permutation: atomstart claims negative number of "
                      "basis functions per atom");}
      if(atomstart[i]<0)
	{throw InputF("Permutation: atomstart contains negative index");}
      atomorder[i]=i;
      
      /* Naive blocking */
      tmpsize+=bfperatom[i];
      n++;
      if (i==natoms-1 || tmpsize+atstart[i+2]-atstart[i+1]>maxblocksize)
      {
	
	natomsperblock[nblocks]=n;
	nblocks++;
	blockstart[nblocks]=blockstart[nblocks-1]+n;
	tmpsize=0;
	n=0;
      }
    }
  atomstart[natoms]=atstart[natoms];  
  
  }
catch(...)
  {
    throw;
  }
}

const Permutation& Permutation::operator=(const Permutation& A)
{
  xpos=xpos;
  ypos=ypos;
  zpos=zpos;
  atomorder=atomorder;
  natoms=natoms;
  natomsperblock=natomsperblock;
  nblocks=nblocks;
  blockstart=blockstart;
  atomstart=atomstart;
  bfperatom=bfperatom;
  maxblocksize=maxblocksize;
  return *this;
}
double Permutation::distance(const integer* atomset1, const integer* atomset2,
		  const integer nset1,const integer nset2)
{
try
  {
    if (atomset1[0]<0 || atomset1[0]>=natoms || 
	atomset2[0]<0 || atomset2[0]>=natoms)
      {throw IndexF("Permutation::distance:atomset out of index bounds");}
    double dist=pow(xpos[atomset1[0]]-xpos[atomset2[0]],2)+
      pow(ypos[atomset1[0]]-ypos[atomset2[0]],2)+
      pow(zpos[atomset1[0]]-zpos[atomset2[0]],2);  
    double tmp;
    for(integer i=0;i<nset1;i++)
      {
	if (atomset1[i]<0 || atomset1[i]>=natoms)
	  {throw IndexF("Permutation::distance:atomset out of index bounds");}
	for(integer j=0;j<nset2;j++)
	  {
	    if (atomset2[j]<0 || atomset2[j]>=natoms)
	      {throw IndexF("Permutation::distance:atomset out of index bounds");}
	    tmp=pow(xpos[atomset1[i]]-xpos[atomset2[j]],2)+
	      pow(ypos[atomset1[i]]-ypos[atomset2[j]],2)+
	      pow(zpos[atomset1[i]]-zpos[atomset2[j]],2);
	    dist=(dist<=tmp ? dist : tmp);
	  }
      }
    return dist;
  }
catch(...)
  {
    throw;
  }
}

void Permutation::distmatrixtofile(char* file)
{
  try
    {
      std::ofstream output(file);
      if (!output)
	{throw IOF("Cannot open outputfile "+std::string(file));}
      integer* rowatomset=atomorder;
      integer* colatomset=atomorder;
      integer ncol,nrow;
      for (integer col=0;col<nblocks;col++)
	{
	  ncol=natomsperblock[col];
	  for (integer row=0;row<nblocks;row++)
	    {
	      nrow=natomsperblock[row];
	      output<<distance(colatomset,
			       rowatomset,
			       natomsperblock[col],
			       natomsperblock[row])<<' ';
	      rowatomset+=nrow;
	    }
	  output<<std::endl;
	  rowatomset=atomorder;
	  colatomset+=ncol;
	}
      
    }
  catch(...)
    {throw;}
}

void Permutation::permutationtofile(char* file)
{
try
  {
  std::ofstream output(file);
  if (!output)
    {throw IOF("Cannot open outputfile "+std::string(file));}
  for(integer i=0;i<natoms;i++)
    {
      output<<atomorder[i]<<' ';
    }
  output<<std::endl;
  for(integer i=0;i<nblocks;i++)
    {
      output<<natomsperblock[i]<<' ';
    }
  for (integer i=nblocks;i<natoms;i++)
    {
      output<<0<<' ';
    }
  }
  catch (...)
    {
      throw;
    }
}


Permutation::~Permutation()
    {
      delete[] xpos;
      delete[] ypos;
      delete[] zpos;
      delete[] atomorder;
      delete[] natomsperblock;
      delete[] blockstart;
      delete[] atomstart;
      delete[] bfperatom;
      

    }
}
