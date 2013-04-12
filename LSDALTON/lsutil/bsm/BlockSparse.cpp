#include "BlockSparse.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <assert.h>

namespace bsm
{
  BlockSparse::BlockSparse(const BlockSparse &A, bool transposed)
    :Sparse<Block >(A, transposed),perm(A.perm),
     totalnrows(transposed ? A.totalncols : A.totalnrows),
     totalncols(transposed ? A.totalnrows : A.totalncols)
  {}
  
  const BlockSparse& BlockSparse::operator=(const BlockSparse& A)
  {
    try
    {
      nrofrows=A.nrofrows;
      if(nrofcols!=A.nrofcols)
	{
	  delete[] colpoint;
	  nrofcols=A.nrofcols;
	  colpoint=new integer[nrofcols+1];
	}  
      if (capacity<A.nrofel)
	{
	  delete[] rowind;
	  delete[] elements;
	  capacity=A.nrofel;
	  rowind=new integer[capacity];
	  elements=new Block[capacity];
	}
      for(integer i=0;i<nrofcols+1;i++)
	{
	  colpoint[i]=A.colpoint[i];
	}
      nrofel=A.nrofel;
      for (integer i=0;i<nrofel;i++)
	{
	  rowind[i]=A.rowind[i];
	  elements[i]=A.elements[i];
	}
      perm=A.perm;
      totalnrows=A.totalnrows;
      totalncols=A.totalncols;
    }
  catch(...)
    {throw;}
    return *this;
  }

  BlockSparse::BlockSparse(Permutation& p,integer fullsize, const real* full,
			   const double cutoff)
    : Sparse<Block >(p.nblocks,p.nblocks), perm(p),
     totalnrows(fullsize),totalncols(fullsize)
  {
    try
      {
	
	if (cutoff<=0)
	  {throw InputF("BlockSparse:Cutoff is negative");}
	if (fullsize<p.natoms)
	  {throw InputF("BlockSparse:Total matrix dimension is smaller than the number of atoms");}
	integer ncol; /* Number of atoms in the current column atom set */
	integer nrow; /* Number of atoms in the current row atom set    */
	
	/* Two sets of atoms that defines a block */
	integer* colatomset=perm.atomorder; /* Column atom set */
	integer* rowatomset=perm.atomorder; /* Row atom set    */
	
	/* Loop over column blocks */
	for(integer col=0;col<nrofcols;col++)
	  {
	    colpoint[col]=nrofel;
	    ncol=perm.natomsperblock[col];
	    /* Loop over row blocks */
	    for(integer row=0;row<nrofrows;row++)
	      {
		nrow=perm.natomsperblock[row];
		
		if (cutoff>=1) /* Truncation using radius assumed */
		  {
		    /* Only if the current atomsets(row set and column set) */
		    /* lies close enough a new block is added               */
		    if (perm.distance(colatomset,rowatomset,ncol,nrow)
			<pow(cutoff,2)
			|| (row == col))
		      /* Keep submatrix if it lies on the diagonal */
		      {
			rowind[nrofel]=row;
			/* Each block findes his own elements in the full matrix */
			elements[nrofel]=Block(colatomset,ncol,rowatomset,nrow,
					       perm.atomstart,perm.bfperatom,
					       perm.natoms,full,
					       totalnrows,totalncols);
			nrofel++;
			
		      }
		    
		  }
		else /* Truncation using 1-norm of complete matrix assumed */
		  {
		    elements[nrofel]=Block(colatomset,ncol,rowatomset,nrow,
					       perm.atomstart,perm.bfperatom,
					       perm.natoms,full,
					       totalnrows,totalncols);
		    if (elements[nrofel].norm()>cutoff*elements[nrofel].getpart() 
			|| (row == col))
		      /* Keep submatrix if it lies on the diagonal */
		      {
			
			rowind[nrofel]=row;
			nrofel++;
		      }
		  }
		rowatomset+=nrow; /* Next atom set */
	      }
	    
	    /* Expand arrays if the number of elements    */
	    /* so far is near present capacity            */
	    if (nrofel>capacity-nrofrows)
	      {
		/* Assume the matrix has about the same number of       */
		/* elements in each column (totally : nrofel*cols/(j+1) */
		/* But be sure not to run out of storage (+rows)        */
		integer newcapacity=nrofel*nrofcols/(col+1)+nrofrows;
		requeststorage(newcapacity);
	      }
	    rowatomset=perm.atomorder;
	    colatomset+=ncol;
	  }
	colpoint[nrofcols]=nrofel;
	
      }
    catch(...)
      {throw;}
  }
  
  BlockSparse::BlockSparse(Permutation& p,integer fullsize)
    : Sparse<Block >(p.nblocks,p.nblocks), perm(p),
     totalnrows(fullsize),totalncols(fullsize)
  {}
  

  const BlockSparse& BlockSparse::operator=(const real* full)
  {
    try
      {
	integer ncol; /* Number of atoms in the current column atom set */
	integer nrow; /* Number of atoms in the current row atom set    */
	
	integer* colatomset;       /* Column atom set */
	integer* rowatomset;       /* Row atom set    */
	integer row;
	for (integer col=0;col<nrofcols;col++)
	  {
	    ncol=perm.natomsperblock[col];
	    colatomset=&perm.atomorder[perm.blockstart[col]];
	    for(integer ind=colpoint[col];ind<colpoint[col+1];ind++)
	      {
		row=rowind[ind];
		nrow=perm.natomsperblock[row];
		rowatomset=&perm.atomorder[perm.blockstart[row]];
		elements[ind].replace(colatomset,ncol,rowatomset,nrow,
				      perm.atomstart,perm.bfperatom,
				      perm.natoms,full,
				      totalnrows,totalncols);
	      }
	  }
	return *this;
      }
    catch(...)
      {throw;}
  }


  /*Fungerar endast med fock, dens struktur. Kan ha if-sats om 
    nrofrows==nrofcols och ropa på ny överlagrad fullmatrix i block
    Konstruktor bör kanske implementeras först */
  void BlockSparse::fullmatrix(real* full)
  {
    try
      {
	for(integer i=0;i<totalnrows*totalncols;i++)
	  {
	    full[i]=0;
	  }
	integer ncol; /* Number of atoms in the current column atom set */
	integer nrow; /* Number of atoms in the current row atom set    */
	
	integer* colatomset;//=perm.atomorder; /* Column atom set */
	integer* rowatomset;//=atomorder; /* Row atom set    */
	integer row;
	for (integer col=0;col<nrofcols;col++)
	  {
	    ncol=perm.natomsperblock[col];
	    colatomset=&perm.atomorder[perm.blockstart[col]];
	    for(integer ind=colpoint[col];ind<colpoint[col+1];ind++)
	      {
		row=rowind[ind];
		nrow=perm.natomsperblock[row];
		rowatomset=&perm.atomorder[perm.blockstart[row]];
		elements[ind].fullmatrix(colatomset,ncol,rowatomset,nrow,
					 perm.atomstart,perm.bfperatom,
					 perm.natoms,full,
					 totalnrows,totalncols);
	      }
	  }
      }
    catch(...)
      {throw;}
  }

  void BlockSparse::extract_diag(real* diag) const
  {
    integer ncol; /* Number of atoms in the current column atom set */
    integer nrow; /* Number of atoms in the current row atom set    */

    integer* colatomset;//=perm.atomorder; /* Column atom set */
    integer* rowatomset;//=atomorder; /* Row atom set    */
    integer row;
    for (integer col=0; col<nrofcols; col++) {
      ncol=perm.natomsperblock[col];
      colatomset=&perm.atomorder[perm.blockstart[col]];
      for(integer ind=colpoint[col]; ind<colpoint[col+1]; ind++) {
        row=rowind[ind];
        nrow=perm.natomsperblock[row];
        rowatomset=&perm.atomorder[perm.blockstart[row]];
        elements[ind].extract_diag(colatomset, ncol, rowatomset, nrow,
                                   perm.atomstart, perm.bfperatom,
                                   perm.natoms, diag,
                                   totalnrows,totalncols);
      }
    }
  }

  void BlockSparse::extract_diag_internal(real *diag)
  {
    integer offset = 0;
    for(integer col = 0; col<nrofcols; col++) {
      bool enc = false;
      for(integer rind = colpoint[col]; rind<colpoint[col+1]; rind++) {
        if(rowind[rind] == col) {
          elements[rind].extract_diag_intern(diag+offset);
          offset += elements[rind].getncols();
          enc = true;
          break;
        }
      }
      if(!enc)
        throw Failure("extract_diag_internal: No diagonal element "
                      "encountered");
    }
  }

  void BlockSparse::setup_rcdim_arrs(integer *rdim, integer *cdim) {
    for(integer i=0; i<nrofrows; i++) rdim[i] = -1;
    for(integer i=0; i<nrofcols; i++) cdim[i] = -1;
    for(integer col = 0; col<nrofcols; col++) {
      for(integer rind = colpoint[col]; rind<colpoint[col+1]; rind++) {
        integer row = rowind[rind];
        rdim[row] = elements[rind].getnrows();
        cdim[col] = elements[rind].getncols();
      }
    }
#if 1
    for(integer i=0; i<nrofrows; i++)
      if(rdim[i] == -1) throw Failure("Empty row!?");

    for(integer i=0; i<nrofcols; i++)
      if(cdim[i] == -1) throw("Empty column!?");
#endif
    integer off = 0;
    for(integer i=0; i<nrofrows; i++) {
      integer tmp = rdim[i]; rdim[i] = off;  off += tmp;
    }

    off = 0;
    for(integer i=0; i<nrofcols; i++) {
      integer tmp = cdim[i]; cdim[i] = off; off += tmp;
    }
  }
  void BlockSparse::precond_ao(integer symm, const real *fup, const real *fuq,
			       const real *du, real omega)
  {
    integer *rdim = new integer[nrofrows];
    integer *cdim = new integer[nrofcols];
    setup_rcdim_arrs(rdim, cdim);
    for(integer col = 0; col<nrofcols; col++) {
      for(integer rind = colpoint[col]; rind<colpoint[col+1]; rind++) {
        integer row = rowind[rind];
	if(symm == 1 || symm == 2) 
	  elements[rind].precond_ao_12(rdim[row], cdim[col],
				       fup, fuq, omega);
	else
	  elements[rind].precond_ao_other(rdim[row], cdim[col],
					  fup, fuq, du, omega);
	  
      }
    }  
    delete rdim; delete cdim;
  }

  /** Zero one of the triangles. 'U' upper half; 'L' - lower half with
     diagonal. */
  void BlockSparse::zerohalf(const char *triangle)
  {
    integer *rdim = new integer[nrofrows];
    integer *cdim = new integer[nrofcols];
    bool low_tr = *triangle == 'l' || *triangle == 'L';

    setup_rcdim_arrs(rdim, cdim);
    for(integer col = 0; col<nrofcols; col++) {
      for(integer rind = colpoint[col]; rind<colpoint[col+1]; rind++) {
        integer row = rowind[rind];
	if(low_tr) {
	  if(row>col) elements[rind] *= 0;
	  else if(row == col)
	    elements[rind].zero_lower_half();
	} else {
	  if(row<col) elements[rind] *= 0;
	  else if(row == col)
	    elements[rind].zero_upper_half();
	}
      }
    }
    delete rdim; delete cdim;
  }

 /** Multiply diagonal by given number */
  void BlockSparse::scal_dia(real alpha) {
    integer *rdim = new integer[nrofrows];
    integer *cdim = new integer[nrofcols];

    setup_rcdim_arrs(rdim, cdim);
    for(integer col = 0; col<nrofcols; col++) {
      for(integer rind = colpoint[col]; rind<colpoint[col+1]; rind++) {
        integer row = rowind[rind];
	if(col == row)
	  elements[rind].scal_dia(alpha);
      }
    }
    delete rdim; delete cdim;
  }

  float BlockSparse::elsparsity()
  {
    integer count=0;
    for (integer i=0;i<nrofel;i++)
      {
	count+=elements[i].getnrows()*elements[i].getncols();
      }
    
    return count/((float)(totalnrows*totalncols));
  }
  void BlockSparse::maxelementstofile(char* file)
  {
    try
    {
      std::ofstream output(file);
      if (!output)
        {
          throw IOF("BlockSparse::maxelementstofile:Cannot open outputfile "+
                    std::string(file));
        }
      output<<std::setprecision(PREC);
       
      
       
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
                  output<<elements[tmpcol[i]].maxabs()<<' ';
                }
              else
                {
                  output<<ZERO<<' ';
                }
            }
          output<<std::endl;
        }
      delete[] tmpcol;
    }
  catch(...)
    {throw;}

  }
  void BlockSparse::norm1ofblockstofile(char* file)
  {
    try
      {
	std::ofstream output(file);
	if (!output)
	  {
	    throw IOF("BlockSparse::norm1ofblockstofile:Cannot open outputfile "+
		      std::string(file));
	  }
	output<<std::setprecision(PREC);
	
	
	
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
		    output<<elements[tmpcol[i]].norm()/elements[tmpcol[i]].getpart()<<' ';
		  }
		else
		  {
		    output<<ZERO<<' ';
		  }
	      }
	    output<<std::endl;
	  }
	delete[] tmpcol;
      }
    catch(...)
      {throw;}
    
  }

bool BlockSparse::read(void *unit,
                       void (*read_int) (void *unit, integer *cnt, integer *where),
                       void (*read_real)(void *unit, integer *cnt, real *where))
{
  integer dt[5], sz;
  static integer FIVE = 5;

  read_int(unit, &FIVE, dt);
  if(dt[0] != totalnrows ||
     dt[1] != totalncols)
    throw Failure("BlockSparse::read: unexpected matrix size.");

  /* Symbolic part */
  nrofrows = dt[2]; nrofcols = dt[3];
  sz = nrofcols+1; read_int(unit, &sz, colpoint);
  if(dt[4]>capacity) {/* only extend */
    capacity = dt[4]; 
    delete[] rowind;
    delete[] elements;
    rowind   = new integer[capacity];
    elements = new Block[capacity];
  }
  nrofel = dt[4];
  read_int(unit, &dt[4], rowind);

  for(integer i=0; i<nrofel; i++)
    elements[i].read(unit, read_int, read_real);

  return true;
}

bool BlockSparse::write(void *unit,
                        void (*write_int) (void *unit, integer *cnt, integer *where),
                        void (*write_real)(void *unit, integer *cnt, real *where))
{
  integer dt[5];
  static integer FIVE = 5;
  dt[0] = totalnrows; 
  dt[1] = totalncols;
  dt[2] = nrofrows; 
  dt[3] = nrofcols;
  dt[4] = nrofel;
  write_int(unit, &FIVE, dt);
  dt[0] = nrofcols+1; write_int(unit, &dt[0], colpoint);
  write_int(unit, &nrofel, rowind);

  for(integer i=0; i<nrofel; i++)
    elements[i].write(unit, write_int, write_real);

  return true;
}

void BlockSparse::identity(void)
{
  if (capacity<nrofrows) { /* We need all the diagonal elements */
    delete[] rowind;
    delete[] elements;
    capacity = nrofrows;
    rowind   = new integer[capacity];
    elements = new Block[capacity];
    /* Observe that Blocks have wrong dimensions at this stage. */
  }

  integer ncol; /* Number of atoms in the current column atom set */
  integer nrow; /* Number of atoms in the current row atom set    */
	
  /* Two sets of atoms that defines a block */
  integer* colatomset = perm.atomorder; /* Column atom set */
  integer* rowatomset = perm.atomorder; /* Row atom set    */

  nrofel = 0;
  /* Loop over column blocks */
  for(integer col=0; col<nrofcols; col++) {
    colpoint[col]=nrofel;
    ncol=perm.natomsperblock[col];
    /* Loop over row blocks */
    for(integer row=0; row<nrofrows; row++) {
      nrow=perm.natomsperblock[row];
      if(col == row) {
        elements[nrofel]=Block(colatomset,ncol,rowatomset,nrow,
                               perm.atomstart,perm.bfperatom,
                               perm.natoms, NULL,
                               totalnrows,totalncols);
        elements[nrofel].identity();
        rowind[nrofel]=row;
        nrofel++;
      }
      rowatomset+=nrow; /* Next atom set */
    }
    rowatomset=perm.atomorder;
    colatomset+=ncol;
  }

  colpoint[nrofcols]=nrofel;
}

void BlockSparse::addscaleye(real scal)
{
  try
    {
      integer row;
      integer ind;
      for(integer col=0;col<nrofcols;col++)
        {
	  ind=colpoint[col];
	  row=rowind[ind];

	  /* Possible to do a faster search for the right row instead */
	  while(row!=col && ind<colpoint[col+1])
	    {
	      ind++;
	      row=rowind[ind];
	    }
	  if (row==col && ind<colpoint[col+1])
	    {
	      //std::cout<<"row="<<row<<"col="<<col<<std::endl;
	      elements[ind].addscaleye(scal);
	    }
	  else
	    {
	      throw Failure
		("BlockSparse::addscaleye: Nonzero diagonal blocks assumed, "
		 "probably something wrong with the matrix?\n"
		 "Maybe the constructor should always add blocks at "
		 "the diagonal, even if they're zero");
	    }
	}
    }
  catch(...)
    {throw;}
}

real BlockSparse::trace() const
{
  try
    {
      integer row;
      integer ind;
      real res=0;
      for(integer col=0;col<nrofcols;col++)
        {
	  ind=colpoint[col];
	  row=rowind[ind];
	  
	  /* Possible to do a search for the right row instead */
	  while(row!=col && ind<colpoint[col+1])
	    {
	      ind++;
	      row=rowind[ind];
	    }
	  if (row==col && ind<colpoint[col+1])
	    {
	      //std::cout<<"row="<<row<<"col="<<col<<std::endl;
	      res+=elements[ind].trace();
	    }
	}
      return res;
    }
  catch(...)
    {throw;}
}

/* Tr(A * B) = Tr(*this * B)*/
real BlockSparse::trace_ab(const BlockSparse& B) const {
  try {
    /* Check for correct dimensions */
    if (nrofcols == B.nrofrows){
      integer rowB;
      integer colA;
      integer rowA;
      integer indA;
      real tr = 0;
      /* For all columns in B */
      for(integer colB = 0; colB < B.nrofcols; colB++) {
	/* For all submatrices in column colB in B */
	for(integer indB = B.colpoint[colB]; indB < B.colpoint[colB + 1]; indB++) {
	  rowB = B.rowind[indB];
	  /* Find corresponding submatrix in A */
	  /* B_ij -> find submatrix A_ji if nonzero */
	  colA = rowB;
	  indA = colpoint[colA];
	  rowA = rowind[indA];
	  while(rowA != colB && indA < colpoint[colA + 1]) {
	    indA++;
	    rowA = rowind[indA];
	  }
	  /* If nonzero */
	  if (rowA == colB && indA < colpoint[colA + 1]) {
	    tr += elements[indA].trace_ab(B.elements[indB]);
	  }
	}
      }
      return tr;
    }
    else {
      throw DimensionF("BlockSparse:trace_ab: Incorrect matrix dimensions for multiplication");
    }
  }
  catch(...)
    {throw;}
}

/* Tr(A' * B) = Tr((*this)' * B)*/
real BlockSparse::trace_atransb(const BlockSparse& B) const {
  try {
    /* Check for correct dimensions */
    if (nrofrows == B.nrofrows){
      integer rowB;
      integer colA;
      integer rowA;
      integer indA;
      real tr = 0;
      /* For all columns in B */
      for(integer colB = 0; colB < B.nrofcols; colB++) {
	/* For all submatrices in column colB in B */
	colA = colB;
	for(integer indB = B.colpoint[colB]; indB < B.colpoint[colB + 1]; indB++) {
	  rowB = B.rowind[indB];
	  /* Find corresponding submatrix in A */
	  /* B_ij -> find submatrix A_ij if nonzero */
	  
	  indA = colpoint[colA];
	  rowA = rowind[indA];
	  while(rowA != rowB && indA < colpoint[colA + 1]) {
	    indA++;
	    rowA = rowind[indA];
	  }
	  /* If nonzero */
	  if (rowA == rowB && indA < colpoint[colA + 1]) {
	    tr += elements[indA].trace_atransb(B.elements[indB]);
	  }
	}
      }
      return tr;
    }
    else {
      throw DimensionF("BlockSparse:trace_atransb: Incorrect matrix dimensions for multiplication");
    }
  }
  catch(...)
    {throw;}
}



BlockSparse::BlockSparse(const BlockSparse& A, const BlockSparse& B,
                         const real alpha, const real beta) 
  :Sparse<Block >(A.perm.nblocks,A.perm.nblocks),perm(A.perm),
   totalnrows(A.totalnrows), totalncols(A.totalncols) {
  try {
    if (A.nrofcols==B.nrofcols && A.nrofrows==B.nrofrows) {
      integer indA;
      integer indB;
      integer rowA;
      integer rowB;
      nrofel = 0;
      /* For all columns */
      for(integer col = 0; col < nrofcols; col++) {
	colpoint[col]=nrofel;
	/* All submatrices in column col in B  */
	/* and common submatrices in A and B   */
	for(integer indb = B.colpoint[col]; indb < B.colpoint[col + 1]; indb++) {
	  rowB = B.rowind[indb];
	  indA = A.colpoint[col];
	  rowA = A.rowind[indA];	   
	  while(rowA != rowB && indA < A.colpoint[col + 1]) {
	    indA++;
	    rowA = A.rowind[indA];
	  }
	  if (rowA == rowB && indA < A.colpoint[col + 1]) {
	     elements[nrofel] = alpha * A.elements[indA] + beta * B.elements[indb];
	  }
	  else {
	    elements[nrofel] = beta * B.elements[indb];
	  }
	  rowind[nrofel]=rowB;
	  nrofel++;
	}
	/* The rest of A:s submatrices */
	for(integer inda = A.colpoint[col]; inda < A.colpoint[col + 1]; inda++) {
	  rowA = A.rowind[inda];
	  indB = B.colpoint[col];
	  rowB = B.rowind[indB];
	  while(rowA != rowB && indB < B.colpoint[col + 1]) {
	    indB++;
	    rowB = B.rowind[indB];
	  }
	  if (!(rowA == rowB && indB < B.colpoint[col + 1])) {
	    elements[nrofel] = alpha * A.elements[inda];
	    rowind[nrofel]=rowA;
	    nrofel++;	  
	  } 
	}
	/* Expand storage arrays if necessary */
	if (nrofel>capacity-nrofrows) {
	  /* Assume the matrix has about the same number of       */
	  /* elements in each column (totally : nrofel*cols/(j+1) */
	  /* But be sure not to run out of storage (+rows)        */
	  integer newcapacity=nrofel*nrofcols/(col+1)+nrofrows;
	  requeststorage(newcapacity);
	}
      }
      colpoint[nrofcols]=nrofel;
    }
    else
      {
	throw DimensionF("BlockSparse: Incorrect matrixdimensions for addition");
      }
  }
  catch(...)
    {throw;}
}

real BlockSparse::frob() const {
  real res = 0;
  real tmp;
  for (integer i = 0; i < nrofel; i++) {
    tmp = elements[i].frob_squared();
    res += tmp;
  }
  res = sqrt(res);
  return res;
} 

void BlockSparse::max_abs_diag(integer *pos, real *val) const
{
  real *v = new real[totalnrows];
  extract_diag(v);
  *pos = 0;
  real mav = fabs(v[*pos]);
  for(integer i=1; i<totalnrows; i++) {
    real av = fabs(v[i]);
    if(av>mav) {
      *pos = i;
      mav = av;
    }
  }
  *val = v[*pos];
  delete v;
}

real BlockSparse::max() const {
  assert(nrofel>0);
  real res = elements[0].max();
  for (integer i = 1; i < nrofel; i++) {
    real tmp = elements[i].max();
    if(fabs(tmp - elements[i].maxabs()) > 1e-7)
      /* printf("block %i: max: %g maxabs: %g\n", i, tmp, elements[i].maxabs()); */
    if(tmp>res)
      res = tmp;
  }
  return res;
} 

real BlockSparse::sum_outdia_sqnorm2() const {
  real res = 0;

  for(integer col = 0; col < nrofcols; col++) {
    /* For all submatrices in column col */
    for(integer ind = colpoint[col]; ind < colpoint[col + 1]; ind++) {
      integer row = rowind[ind];
      if(col == row)
        res += elements[ind].sum_outdia();
      else
        res += elements[ind].frob_squared();
    }
  }
  return res;
} 

void
BlockSparse::print_struct() const {
  puts("matrix colpoints");
  for(integer i=0; i<nrofcols; i++) printf("%d ", colpoint[i]);
  puts("matrix rowind");
  for(integer col = 0; col < nrofcols; col++)
    for(integer ind = colpoint[col]; ind < colpoint[col + 1]; ind++) {
      printf("Col: %d ind: %d rowind: %d ", col, ind, rowind[ind]);
    }
  for(integer i=0; i<nrofcols; i++) printf("%d ", colpoint[i]);
    puts("\nmatrix blocks:");
  for(integer col = 0; col < nrofcols; col++) {
    /* For all submatrices in column col */
    for(integer ind = colpoint[col]; ind < colpoint[col + 1]; ind++) {
      integer row = rowind[ind];
      printf("%2d: block in col: %d row: %d has size %d x %d | %p\n", ind, col, row,
             elements[ind].getnrows(), elements[ind].getncols(), &elements[ind]);
      elements[ind].print();
    }
  }
}
} /* namespace bsm */
