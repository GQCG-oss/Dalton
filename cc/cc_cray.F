!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
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

