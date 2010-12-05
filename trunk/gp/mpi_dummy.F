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
C     File: mpi_dummy.F
C
C     A set of dummy mpi routines for creating a sequential
C     dalton.x together with a parallel dalmpi.x
C     (See main Makefile) - hjaaj Mar 2005
C
C
#define LUPRI_MPI 0
C
C HJAaJ: these dummy MPI routines MUST return the necessary information
C        for sequential runs - they are only called for sequential dalton.x
C        linked with files where VAR_MPI was defined.
C        They must return MYNUM = 0 and NODTOT = 1. /061003
C
      SUBROUTINE MPI_INIT(IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_INIT called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_COMM_RANK(I1,MYNUM,IERR)
         MYNUM = 0
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_COMM_RANK called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_COMM_SIZE(I1,NODTOT,IERR)
         NODTOT = 1
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_COMM_SIZE called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_FINALIZE(IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_FINALIZE called, exiting ...'
         CALL QTRACE(LUPRI_MPI)
#endif
         STOP
C        ... normal exit of sequential calculation
      RETURN
      END
      SUBROUTINE MPI_ABORT()
         STOP 100
      RETURN
      END
      SUBROUTINE MPI_BCAST(i1,i2,i3,i4,i5,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_BCAST called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_RECV(BUFFER,COUNT,DATATYPE,MPI_ANY_SOURCE,TAG,
     &                    MPI_COMM_WORLD,ISTAT,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_RECV called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_SEND(BUFFER,COUNT,DATATYPE,DEST,TAG,
     &                    MPI_COMM_WORLD,IERR)
         IERR = 0
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_SEND called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_REDUCE()
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_REDUCE called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_BARRIER()
#ifdef MPI_DEBUG
         WRITE (LUPRI_MPI,*)
     &      'info: dummy MPI_BARRIER called, nothing done'
         CALL QTRACE(LUPRI_MPI)
#endif
      RETURN
      END
      SUBROUTINE MPI_ERROR_CLASS(IERR,IERRCL,IERR2)
        WRITE(LUPRI_MPI,*)
     &      'info: dummy MPI_ERROR_CLASS called, ierr =',IERR
        CALL QUIT('ERROR, dummy MPI_MYFAIL called')
      RETURN
      END
