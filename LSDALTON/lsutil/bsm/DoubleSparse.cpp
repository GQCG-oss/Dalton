#include "DoubleSparse.h"
#include <fstream> 
#include <iostream>
#include <iomanip> /* For setprecision in fstream */
#include <cmath>

namespace bsm
{
DoubleSparse::DoubleSparse(char* file,integer rows,integer cols,double thres, integer cap)
  :Sparse<double>(rows,cols,cap)
{
  
  std::ifstream input(file);
  if (!input)
    {
      std::cerr<<"Cannot open inputfile "<<file<<std::endl;
      std::exit(1);
    }
  input>>std::setprecision(PREC);
  double tmpvalue;
  integer count=0;
  for (integer j=0;j<cols;j++)
  {
    colpoint[j]=nrofel;
    for (integer i=0;i<rows;i++)
      {
	if (input.eof())
	  {
	    std::cerr<<"Inputfile "<<file<<" does not match given matrix dimensions"<<std::endl;
	    std::exit(1);
	  }
	input>>tmpvalue;
	if (fabs(tmpvalue)>thres)
	  {
	    rowind[nrofel]=i;
	    elements[nrofel]=tmpvalue;
	    nrofel++;
	  }
      }
    
    /* Expand arrays if the total number of elements is unknown   */
    /* and the number of elements so far is near present capacity */
    if (cap==0 && nrofel>capacity-rows)
      {

	/* Assume the matrix has about the same number of       */
        /* elements in each column (totally : nrofel*cols/(j+1) */
        /* But be sure not to run out of storage (+rows)        */
	integer newcapacity=nrofel*cols/(j+1)+rows;
	requeststorage(newcapacity);
	count++;
      }
    
  }
  input.close();
  colpoint[cols]=nrofel;
  std::cout<<"nrofel: "<<nrofel<<std::endl;
  std::cout<<"reqstorage: "<<count<<std::endl;
  std::cout<<"Capacity: "<<capacity<<std::endl;
}

DoubleSparse::DoubleSparse(double* full,integer rows,integer cols,double thres, integer cap)
  :Sparse<double>(rows,cols,cap)
{
  try
    {
      double tmpvalue;
      integer count=0;
      for (integer j=0;j<cols;j++)
	{
	  colpoint[j]=nrofel;
	  for (integer i=0;i<rows;i++)
	    {
	      tmpvalue=full[j*rows+i];
	      if (fabs(tmpvalue)>thres)
		{
		  rowind[nrofel]=i;
		  elements[nrofel]=tmpvalue;
		  nrofel++;
		}
	    }
	  
	  /* Expand arrays if the total number of elements is unknown   */
	  /* and the number of elements so far is near present capacity */
	  if (cap==0 && nrofel>capacity-rows)
	    {
	      
	      /* Assume the matrix has about the same number of       */
	      /* elements in each column (totally : nrofel*cols/(j+1) */
	      /* But be sure not to run out of storage (+rows)        */
	      integer newcapacity=(integer)(nrofel*((double)cols/(j+1))+rows);
	      requeststorage(newcapacity);
	      count++;
	    }
	  
	}
      colpoint[cols]=nrofel;
      /*
	std::cout<<"nrofel: "<<nrofel<<std::endl;
	std::cout<<"reqstorage: "<<count<<std::endl;
	std::cout<<"Capacity: "<<capacity<<std::endl;
      */
    }
  catch(...)
    {throw;}
}


void DoubleSparse::storefullonfile(char* file)
{
  try
    {
      std::ofstream output(file);
      if (!output)
	{
	  throw IOF("DoubleSparse::storefullonfile:Cannot open outputfile "+
		    std::string(file));
	}
      output<<std::setprecision(PREC);
      
      const double zero=0;
      
      integer* tmpcol=new integer[nrofrows];
      for(integer v=0;v<nrofrows;v++)
	{
	  tmpcol[v]=-1;
	}
      
      for(integer j=0;j<nrofcols;j++)
	{
	  for(integer ind=colpoint[j];ind<colpoint[j+1];ind++)
	    {
	      tmpcol[rowind[ind]]=ind;
	    }
	  for(integer i=0;i<nrofrows;i++)
	    {
	      
	      if (tmpcol[i]>=colpoint[j])
		{
		  output<<elements[tmpcol[i]]<<' ';
		}
	      else
		{
		  output<<zero<<' ';
		}
	    }
	  output<<std::endl;
	}
      delete[] tmpcol;
    }
  catch(...)
    {throw;}
}
void DoubleSparse::truncate(double eps)
{
  integer newnrofel=0;
  integer colind=0;
  for(integer i=0;i<nrofel;i++)
    {
      if (colpoint[colind]==i)
	{
	  colpoint[colind]=newnrofel;
	  colind++;
	}
      if(elements[i]>eps)
	{
	  elements[newnrofel]=elements[i];
	  rowind[newnrofel]=rowind[i];
	  newnrofel++;	  
	}
     
    }
  colpoint[colind]=newnrofel;
  nrofel=newnrofel;
}
}
