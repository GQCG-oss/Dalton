#include "Symbolic.h"
#include "InputF.h"
#include "IOF.h"
#include <fstream>
#include <iostream>
namespace bsm
{
Symbolic::Symbolic(integer rows,integer cols,integer cap,integer el)
  :Matrix(rows,cols),capacity(cap),nrofel(el)
{
  try
    {
      if (cap<0 || el<0)
	{throw InputF("Symbolic:Capacity and nel must be at least zero");}
      colpoint=new integer[cols+1];
      /* Assume there is at least one element per row */
      if (capacity<rows)
	{
	  capacity=rows;
	}
      rowind=new integer[capacity];
      for(integer i=0; i<cols+1; i++) colpoint[i] = 0;
    }
  catch(...)
    {throw;}
} 

void Symbolic::spy(char* file)
{  
  try
    {
      std::ofstream output(file);
      if (!output)
	{
	  throw IOF("Symbolic::spy:Cannot open outputfile"+std::string(file));
	}
      const double zero=0;
      const double one=1;
      
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
		  output<<one<<' ';
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
Symbolic::~Symbolic()
{
  delete[] colpoint;
  delete[] rowind;
}   
}
