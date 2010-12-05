C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      dimension a(1), b(1), c(1)
      call sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      return
      end
      subroutine dgemv (trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      dimension a(1),x(1),y(1)
      call sgemv (trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      return
      end
      function ddot(n,x,nx,y,ny)
      dimension x(1),y(1)
      ddot=sdot(n,x,nx,y,ny)
      return
      end
      subroutine dswap(n,x,nx,y,ny)
      dimension x(1),y(1)
      call sswap(n,x,nx,y,ny)
      return
      end
      subroutine dcopy(n,x,nx,y,ny)
      dimension x(1),y(1)
      call scopy(n,x,nx,y,ny)
      return
      end
      subroutine dscal(n,a,x,nx)
      dimension x(1)
      call sscal(n,a,x,nx)
      return
      end
      subroutine daxpy(n,a,x,nx,y,ny)
      dimension x(1),y(1)
      call saxpy(n,a,x,nx,y,ny)
      return
      end
      integer function idamax(n,x,nx)
      dimension x(1)
      idamax=isamax(n,x,nx)
      return
      end
      function dnrm2(n,x,incx)
      dimension x(1)
      dnrm2=snrm2(n,x,incx)
      return
      end
      function dasum(n,x,incx)
      dimension x(1)
      dasum=sasum(n,x,incx)
      return
      end

