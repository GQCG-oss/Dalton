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
C-------------------------------------------------------------
C    Paal Dahle Jan-2001
C
C    Macros for MPI calls that make calls to MPI routines
C    easier and more transparent because:
C
C    a) Less arguments to MPIX... calls than to MPI_... calls
C    b) No need to specify the mpif.h include file or array
C       my_STATUS(MPI_SOURCE) in the calling subroutine.
C    c) Error handling is automatically taken care of.
C-------------------------------------------------------------
C
C  /* Deck mpixinit */
      SUBROUTINE MPIXINIT
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <infpar.h>
C
#if !defined (VAR_MPI)
      MYNUM  = 0
      MASTER = 0
      NODTOT = 0
#else
      INCLUDE 'mpif.h'
C
      CALL MPI_INIT(IERR)
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
C
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNUM,IERR)
#if defined(CRI_MEMTRACE)
      IF (MYNUM .LE. 1) CALL MEMTRACE(1)
#else
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
#endif
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NODTOT,IERR)
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
C
C     One process is lost as the master
C
      MASTER = 0
      NODTOT = NODTOT - 1
C
      IF (NODTOT .GT. MAXNOD) THEN
        WRITE(LUPRI,'(/5X,A/12X,A,2(/5X,A20,I8))')
     &  'ERROR: Number of MPI slave nodes allocated has exceeded',
     &  'hardcoded limit. Reduce number of nodes or recompile program.',
     &  'Allocated (NODTOT) :',NODTOT,'Limit (MAXNOD)     :',MAXNOD
        WRITE(LUPRI,'(//4(/5X,A))')
     &  'FOR RECOMPILATION:',
     &  'Increase parameter MAXNOD in infpar.h and execute:',
     &  '   make depend',
     &  '   make'
        CALL QUIT('MPI ERROR: Number of nodes requested is too large.')
      ENDIF
C
      DO I = 0, NODTOT
         NODEID(I) = I
      END DO
#endif
C
      RETURN
      END
C  /* Deck mpifinalize */
      SUBROUTINE MPIXFINALIZE
C
#include <implicit.h>
#include <priunit.h>
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
C
      CALL MPI_FINALIZE(IERR)
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
#endif
C
      RETURN
      END
C  /* Deck mpixbcast */
      SUBROUTINE MPIXBCAST(BUFFER,COUNT,TYPE,ROOT)
C
#include <implicit.h>
#include <priunit.h>
#if defined (VAR_MPI)
#include <maxorb.h>
#include <infpar.h>
      INCLUDE 'mpif.h'
      CHARACTER*6 TYPE
      INTEGER COUNT, DATATYPE, ROOT
      DIMENSION BUFFER(*)
C
      IF      (TYPE .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE .EQ. 'LOGICA') THEN
         DATATYPE = MPI_LOGICAL
      ELSE IF (TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE .EQ. 'STRING' .OR. TYPE .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' TYPE :',TYPE,' does not exist!'
         CALL QUIT('ERROR in MPIXBCAST: Nonexisting DATATYPE')
      END IF
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,ROOT,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
#endif
C
      RETURN
      END
C  /* Deck mpixrecv */
      SUBROUTINE MPIXRECV(BUFFER,COUNT,TYPE,SOURCE,TAG)
C
#include <implicit.h>
#include <priunit.h>
#if defined (VAR_MPI)
#include <maxorb.h>
#include <infpar.h>
      INCLUDE 'mpif.h'
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      CHARACTER*6 TYPE
      INTEGER COUNT, DATATYPE, SOURCE, TAG
      DIMENSION BUFFER(*)
C
      IF      (TYPE .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE .EQ. 'LOGICA') THEN
         DATATYPE = MPI_LOGICAL
      ELSE IF (TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE .EQ. 'STRING' .OR. TYPE .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' TYPE :',TYPE,' does not exist!'
         CALL QUIT('ERROR in MPIXRECV: Nonexisting DATATYPE')
      END IF
      IF (SOURCE.EQ.-1) THEN
         CALL MPI_RECV(BUFFER,COUNT,DATATYPE,MPI_ANY_SOURCE,TAG,
     &                 MPI_COMM_WORLD,my_STATUS,IERR)
         SOURCE = my_STATUS(MPI_SOURCE)
      ELSE
         CALL MPI_RECV(BUFFER,COUNT,DATATYPE,SOURCE,TAG,
     &                 MPI_COMM_WORLD,my_STATUS,IERR)
      END IF
      IF (IERR.GT.0) CALL MPI_MYFAIL(IERR)
#endif
C
      RETURN
      END
C  /* Deck mpixsend */
      SUBROUTINE MPIXSEND(BUFFER,COUNT,TYPE,DEST,TAG)
C
#include <implicit.h>
#include <priunit.h>
#if defined (VAR_MPI)
#include <maxorb.h>
#include <infpar.h>
      INCLUDE 'mpif.h'
      CHARACTER*6 TYPE
      INTEGER COUNT, DATATYPE, DEST, TAG
      DIMENSION BUFFER(*)
C
      IF      (TYPE .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE .EQ. 'LOGICA') THEN
         DATATYPE = MPI_LOGICAL
      ELSE IF (TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE .EQ. 'STRING' .OR. TYPE .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' MPIXSEND : TYPE "',TYPE,'" does not exist!'
         CALL QUIT('ERROR in MPIXBCAST: Nonexisting DATATYPE')
      END IF
      CALL MPI_SEND(BUFFER,COUNT,DATATYPE,DEST,TAG,MPI_COMM_WORLD,IERR)
      IF (IERR.GT.0) THEN
         write (lupri,*) 'ERROR: problem in MPIXSEND, ierr =',IERR
         write (lupri,*) '- count,dest,tag: ',count,dest,tag
         write (lupri,*) '- calling MPI_MYFAIL ...'
         CALL MPI_MYFAIL(IERR)
      END IF
#endif
C
      RETURN
      END

      subroutine mpixprocname(process_name,process_name_length)
C
      implicit none
      character*(*) process_name
      integer process_name_length, ierr
#ifdef VAR_MPI
      include 'mpif.h'
C
      call mpi_get_processor_name(process_name,process_name_length,ierr)
#endif
      END

#if defined (VAR_MPI)
C  /* Deck mpi_myfail */
      SUBROUTINE MPI_MYFAIL(IERR)
C
#include <implicit.h>
#include <priunit.h>
      INCLUDE 'mpif.h'
      CHARACTER ERRBUF*40
C
      CALL MPI_ERROR_CLASS(IERR,IERRCL,IERR2)
C
      IF (IERRCL.EQ.MPI_SUCCESS) THEN
         ERRBUF = 'No error'
      ELSE IF (IERRCL.EQ.MPI_ERR_BUFFER) THEN
         ERRBUF = 'Invalid buffer pointer'
      ELSE IF (IERRCL.EQ.MPI_ERR_COUNT) THEN
         ERRBUF = 'Invalid count argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TYPE) THEN
         ERRBUF = 'Invalid datatype argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TAG) THEN
         ERRBUF = 'Invalid tag argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_COMM) THEN
         ERRBUF = 'Invalid communicator'
      ELSE IF (IERRCL.EQ.MPI_ERR_RANK) THEN
         ERRBUF = 'Invalid rank'
      ELSE IF (IERRCL.EQ.MPI_ERR_REQUEST) THEN
         ERRBUF = 'Invalid request (handle)'
      ELSE IF (IERRCL.EQ.MPI_ERR_ROOT) THEN
         ERRBUF = 'Invalid root'
      ELSE IF (IERRCL.EQ.MPI_ERR_GROUP) THEN
         ERRBUF = 'Invalid group'
      ELSE IF (IERRCL.EQ.MPI_ERR_OP) THEN
         ERRBUF = 'Invalid operation'
      ELSE IF (IERRCL.EQ.MPI_ERR_TOPOLOGY) THEN
         ERRBUF = 'Invalid topology'
      ELSE IF (IERRCL.EQ.MPI_ERR_DIMS) THEN
         ERRBUF = 'Invalid dimension argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_ARG) THEN
         ERRBUF = 'Invalid argument of some other kind'
      ELSE IF (IERRCL.EQ.MPI_ERR_UNKNOWN) THEN
         ERRBUF = 'Unknown error'
      ELSE IF (IERRCL.EQ.MPI_ERR_TRUNCATE) THEN
         ERRBUF = 'Message truncated on receive'
      ELSE IF (IERRCL.EQ.MPI_ERR_OTHER) THEN
         ERRBUF = 'Known error not in this list'
      ELSE IF (IERRCL.EQ.MPI_ERR_INTERN) THEN
         ERRBUF = 'Internal MPI (implementation) error'
      ELSE IF (IERRCL.EQ.MPI_ERR_IN_STATUS) THEN
         ERRBUF = 'Error code is in status'
      ELSE IF (IERRCL.EQ.MPI_ERR_PENDING) THEN
         ERRBUF = 'Pending request'
      ELSE IF (IERRCL.EQ.MPI_ERR_LASTCODE) THEN
         ERRBUF = 'Last error code'
      ELSE
C        something we didn't know ...
         WRITE(ERRBUF,'(A,I4)') 'MPI error class',IERRCL
      END IF
C
      WRITE(LUPRI,'(/A)') ' ERROR in MPI : '//ERRBUF
C
      CALL QUIT('Error detected in MPI. Please consult dalton output!')
C
      RETURN
      END
#endif  /* ifdef VAR_MPI */
