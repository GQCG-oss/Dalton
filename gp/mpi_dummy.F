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
     &                    MPI_COMM_WORLD,my_STATUS,IERR)
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
