#include "DivideSpace.h"
namespace bsm
{
DivideSpace::DivideSpace(const double* xp,const double* yp,
			 const double* zp,const integer* atstart,
			 integer natoms,integer maxbsize)
  :Permutation(xp,yp,zp,atstart,natoms,maxbsize)
{  
  
  nblocks=0;
  sort(0,natoms-1);
}

void DivideSpace::sort(integer begin,integer end)
{
  /* Find direction with largest difference */
  double xmin=xpos[atomorder[begin]];
  double xmax=xpos[atomorder[begin]];
  double ymin=ypos[atomorder[begin]];
  double ymax=ypos[atomorder[begin]];
  double zmin=zpos[atomorder[begin]];
  double zmax=zpos[atomorder[begin]];
  double xdiff;
  double ydiff;
  double zdiff;
  integer ind;
  for(integer i=begin+1;i<=end;i++)
    {  
      ind=atomorder[i];
      xmin=(xmin<xpos[ind] ? xmin : xpos[ind]);
      xmax=(xmax>xpos[ind] ? xmax : xpos[ind]);
      ymin=(ymin<ypos[ind] ? ymin : ypos[ind]);
      ymax=(ymax>ypos[ind] ? ymax : ypos[ind]);
      zmin=(zmin<zpos[ind] ? zmin : zpos[ind]);
      zmax=(zmax>zpos[ind] ? zmax : zpos[ind]);
    }
  xdiff=xmax-xmin;
  ydiff=ymax-ymin;
  zdiff=zmax-zmin;

  /* Sort in direction with largest difference */
  if (xdiff>=ydiff && xdiff>=zdiff)
    {
      quicksort(xpos,begin,end);
    }
  else if (ydiff>zdiff)
    {
      quicksort(ypos,begin,end);
    }
  else
    {
      quicksort(zpos,begin,end);
    }


  /* Divide space in two boxes so that about the        */
  /* same number of basis functions ends up in each box */
  integer low=begin;
  integer high=end;
  integer bflow=bfperatom[atomorder[low]];
  integer bfhigh=bfperatom[atomorder[high]];
  while (high-low>1)
    {
      if (bflow<=bfhigh)
	{
	  low++;
	  bflow+=bfperatom[atomorder[low]];
	}
      else
	{
	  high--;
	  bfhigh+=bfperatom[atomorder[high]];
	}
    }

  /* bflow and bfhigh is the number of           */
  /* basis functions in the two boxes            */

  /* If the box contains more than one atom and the number */
  /* of basis functions is more than maximum blocksize     */
  /* divide again, otherwise add info of the new block     */
  if (bflow>maxblocksize && begin!=low)
    {
      sort(begin,low);
    }
  else
    {
      natomsperblock[nblocks]=low-begin+1;
      nblocks++;
      blockstart[nblocks]=blockstart[nblocks-1]+natomsperblock[nblocks-1];
    }
  if (bfhigh>maxblocksize && high!=end)
    {
      sort(high,end);
    }
  else
    {
      natomsperblock[nblocks]=end-high+1;
      nblocks++;
      blockstart[nblocks]=blockstart[nblocks-1]+natomsperblock[nblocks-1];
    }
}

void DivideSpace::quicksort(const double* pos,integer low,integer high)
{
  integer i=low;
  integer j=high;
  integer tmp;
  double pivot=pos[atomorder[(low+high)/2]];
  do
    {
      while (pos[atomorder[i]]<pivot) i++;
      while (pos[atomorder[j]]>pivot) j--;
      
      if (i<=j)
	{
	  tmp=atomorder[i];
	  atomorder[i]=atomorder[j];
	  atomorder[j]=tmp;
	  i++;
	  j--;
	}
    } while (i<=j);
  if(low<j) quicksort(pos,low,j);
  if(i<high) quicksort(pos,i,high);
}

}
