!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
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
      implicit none
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"

      integer i
C
#if !defined (VAR_MPI)
      MYNUM  = 0
      MASTER = 0
      NODTOT = 0
#else
      INCLUDE 'mpif.h'
      integer(kind=MPI_INTEGER_KIND) :: ierr_mpi, value_mpi
C
#ifndef VAR_CHEMSHELL
      CALL MPI_INIT(ierr_mpi)
      IF (ierr_mpi.GT.0) CALL MPI_MYFAIL(ierr_mpi)
#endif
C
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,value_mpi,ierr_mpi)
      MYNUM = value_mpi
#if defined(CRI_MEMTRACE)
      IF (MYNUM .LE. 1) CALL MEMTRACE(1)
#else
      IF (ierr_mpi.GT.0) CALL MPI_MYFAIL(ierr_mpi)
#endif
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,value_mpi,ierr_mpi)
      NODTOT = value_mpi
      IF (ierr_mpi.GT.0) CALL MPI_MYFAIL(ierr_mpi)
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
        WRITE(LUPRI,'(//2(/5X,A))')
     &  'FOR RECOMPILATION:',
     &  'Increase parameter MAXNOD in infpar.h and rebuild with make'
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
#ifndef VAR_CHEMSHELL
#if defined (VAR_MPI)
      implicit none
      INCLUDE 'mpif.h'
C
      integer(kind=MPI_INTEGER_KIND) ierr_mpi
      CALL MPI_FINALIZE(ierr_mpi)
      IF (ierr_mpi.GT.0) CALL MPI_MYFAIL(ierr_mpi)
#endif
#endif
C
      RETURN
      END
C  /* Deck mpixbcast */
      SUBROUTINE MPIXBCAST(BUFFER,BUF_COUNT,TYPE_in,ROOT)
C
#if defined (VAR_MPI)
      implicit none
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
      INCLUDE 'mpif.h'
      CHARACTER*6 TYPE_in
      INTEGER BUF_COUNT,  ROOT
      INTEGER BUFFER(*) ! does not matter if wrong type here
      integer(kind=MPI_INTEGER_KIND) buf_count_mpi, root_mpi
      integer(kind=MPI_INTEGER_KIND) DATATYPE, ierr_mpi
C
      IF      (TYPE_in .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE_in .EQ. 'LOGICA') THEN
         DATATYPE = my_MPI_LOGICAL
      ELSE IF (TYPE_in .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE_in .EQ. 'STRING' .OR. TYPE_in .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*)
     &   'MPIXBCAST: DATA TYPE ',TYPE_in,' does not exist!'
         WRITE(*,*) 'MPIXBCAST: DATA TYPE ',TYPE_in,' does not exist!'
         WRITE(*,*) ' TYPE :',TYPE_in,' does not exist!'
         CALL QUIT('ERROR in MPIXBCAST: Nonexisting DATATYPE')
      END IF

      buf_count_mpi = BUF_COUNT
      root_mpi      = ROOT
      CALL MPI_BCAST(BUFFER,BUF_COUNT_mpi,DATATYPE,ROOT_mpi,
     &   MPI_COMM_WORLD,ierr_mpi)
      IF (ierr_mpi.GT.0) CALL MPI_MYFAIL(ierr_mpi)
#endif
C
      RETURN
      END
C  /* Deck mpixrecv */
      SUBROUTINE MPIXRECV(BUFFER,BUF_COUNT,TYPE_in,SOURCE,TAG)
C
#if defined (VAR_MPI)
      implicit none
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
      INCLUDE 'mpif.h'
      CHARACTER*6 TYPE_in
      INTEGER BUF_COUNT, SOURCE, TAG
      REAL*8  BUFFER(*)

      integer(kind=MPI_INTEGER_KIND) my_STATUS(MPI_STATUS_SIZE)
      integer(kind=MPI_INTEGER_KIND) buf_count_mpi,source_mpi,tag_mpi
      integer(kind=MPI_INTEGER_KIND) DATATYPE,ierr_mpi
C
      IF      (TYPE_in .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE_in .EQ. 'LOGICA') THEN
         DATATYPE = my_MPI_LOGICAL
      ELSE IF (TYPE_in .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE_in .EQ. 'STRING' .OR. TYPE_in .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' TYPE :',TYPE_in,' does not exist!'
         CALL QUIT('ERROR in MPIXRECV: Nonexisting DATATYPE')
      END IF
      buf_count_mpi = BUF_COUNT
      tag_mpi       = TAG
      IF (SOURCE.EQ.-1) THEN
         CALL MPI_RECV(BUFFER,BUF_COUNT_mpi,DATATYPE,MPI_ANY_SOURCE,
     &                 TAG_mpi,MPI_COMM_WORLD,my_STATUS,ierr_mpi)
         SOURCE = my_STATUS(MPI_SOURCE)
      ELSE
         SOURCE_mpi = SOURCE
         CALL MPI_RECV(BUFFER,BUF_COUNT_mpi,DATATYPE,SOURCE_mpi,
     &                 TAG_mpi,MPI_COMM_WORLD,my_STATUS,ierr_mpi)
      END IF
      IF (ierr_mpi.GT.0) THEN
         write (lupri,*) 'ERROR: problem in MPIXRECV, ierr =',ierr_mpi
         write (lupri,*) '- count,source,tag: ',BUF_count,source_mpi,tag
         write (lupri,*) '- calling MPI_MYFAIL ...'
         CALL MPI_MYFAIL(ierr_mpi)
      END IF
#endif
C
      RETURN
      END
C  /* Deck mpixsend */
      SUBROUTINE MPIXSEND(BUFFER,BUF_COUNT,TYPE_in,DEST,TAG)
C
#if defined (VAR_MPI)
      integer none
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
      INCLUDE 'mpif.h'
      CHARACTER*6 TYPE_in
      INTEGER BUF_COUNT, DEST, TAG
      REAL*8  BUFFER(*)

      integer(kind=MPI_INTEGER_KIND) my_STATUS(MPI_STATUS_SIZE)
      integer(kind=MPI_INTEGER_KIND) buf_count_mpi,dest_mpi,tag_mpi
      integer(kind=MPI_INTEGER_KIND) DATATYPE,ierr_mpi
C
      IF      (TYPE_in .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (TYPE_in .EQ. 'LOGICA') THEN
         DATATYPE = my_MPI_LOGICAL
      ELSE IF (TYPE_in .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (TYPE_in .EQ. 'STRING' .OR. TYPE_in .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' MPIXSEND : TYPE "',TYPE_in,'" does not exist!'
         CALL QUIT('ERROR in MPIXSEND: Nonexisting DATATYPE')
      END IF
      buf_count_mpi = buf_count
      dest_mpi      = DEST
      tag_mpi       = TAG
      CALL MPI_SEND(BUFFER,BUF_COUNT_mpi,DATATYPE,DEST_mpi,TAG_mpi,
     &   MPI_COMM_WORLD,ierr_mpi)
      IF (ierr_mpi.GT.0) THEN
         write (lupri,*) 'ERROR: problem in MPIXSEND, ierr =',ierr_mpi
         write (lupri,*) '- count,dest,tag: ',BUF_count,dest,tag
         write (lupri,*) '- calling MPI_MYFAIL ...'
         CALL MPI_MYFAIL(ierr_mpi)
      END IF
#endif
      END

      SUBROUTINE MPIX_ISEND(BUFFER,BUF_COUNT,BUF_TYPE,DEST,REQUEST,TAG)
C
#if defined (VAR_MPI)
      implicit none
#include "priunit.h"
      CHARACTER*6 BUF_TYPE
      INTEGER     BUFFER(*) ! does not matter if wrong type here
      INTEGER     BUF_COUNT, DEST, REQUEST, TAG
#include "maxorb.h"
#include "infpar.h"
      INCLUDE 'mpif.h'
      INTEGER(kind=MPI_INTEGER_KIND) :: buf_count_MPI
      INTEGER(kind=MPI_INTEGER_KIND) :: dest_MPI, request_MPI, tag_MPI
      INTEGER(kind=MPI_INTEGER_KIND) :: DATATYPE, ierr_mpi
C
      IF      (BUF_TYPE .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (BUF_TYPE .EQ. 'LOGICA') THEN
         DATATYPE = my_MPI_LOGICAL
      ELSE IF (BUF_TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (BUF_TYPE .EQ. 'STRING' .OR. BUF_TYPE .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' MPIX_ISEND : DATA TYPE "',BUF_TYPE,
     &      '" does not exist!'
         CALL QUIT('ERROR in MPIX_ISEND: Nonexisting DATATYPE')
      END IF
      buf_count_MPI = BUF_COUNT
      dest_MPI      = DEST
      request_MPI   = REQUEST
      tag_MPI       = TAG
      CALL MPI_ISEND(BUFFER,buf_count_MPI,DATATYPE,dest_MPI,tag_MPI,
     &   MPI_COMM_WORLD,request_MPI,ierr_mpi)
      IF (ierr_mpi.GT.0) THEN
         write (lupri,*) 'ERROR: problem in MPIX_ISEND, ierr =',ierr_mpi
         write (lupri,*) '- count,dest,tag: ',BUF_count,dest,tag
         write (lupri,*) '- calling MPI_MYFAIL ...'
         CALL MPI_MYFAIL(ierr_mpi)
      END IF
#endif
      END

      SUBROUTINE MPIX_IRECV(BUFFER,BUF_COUNT,BUF_TYPE,
     &                      SOURCE,REQUEST,TAG)
C
#if defined (VAR_MPI)
      implicit none
#include "priunit.h"
      CHARACTER*6 BUF_TYPE
      INTEGER     BUFFER(*) ! does not matter if wrong type here
      INTEGER     BUF_COUNT, SOURCE, REQUEST, TAG
#include "maxorb.h"
#include "infpar.h"
      INCLUDE 'mpif.h'
      INTEGER(kind=MPI_INTEGER_KIND) :: buf_count_MPI
      INTEGER(kind=MPI_INTEGER_KIND) :: source_MPI, request_MPI, tag_MPI
      INTEGER(kind=MPI_INTEGER_KIND) :: DATATYPE, ierr_mpi
C
      IF      (BUF_TYPE .EQ. 'INTEGE') THEN
         DATATYPE = my_MPI_INTEGER
      ELSE IF (BUF_TYPE .EQ. 'LOGICA') THEN
         DATATYPE = my_MPI_LOGICAL
      ELSE IF (BUF_TYPE .EQ. 'DOUBLE') THEN
         DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (BUF_TYPE .EQ. 'STRING' .OR. BUF_TYPE .EQ. 'CHARAC') THEN
         DATATYPE = MPI_CHARACTER
      ELSE
         WRITE(LUPRI,*) ' MPIX_IRECV : DATA TYPE "',BUF_TYPE,
     &      '" does not exist!'
         CALL QUIT('ERROR in MPIX_IRECV: Nonexisting DATATYPE')
      END IF
      buf_count_MPI = BUF_COUNT
      source_MPI    = SOURCE
      request_MPI   = REQUEST
      tag_MPI       = TAG
      CALL MPI_IRECV(BUFFER,buf_count_MPI,DATATYPE,source_MPI,tag_MPI,
     &   MPI_COMM_WORLD,request_MPI,ierr_mpi)
      IF (ierr_mpi.GT.0) THEN
         write (lupri,*) 'ERROR: problem in MPIX_IRECV, ierr =',ierr_mpi
         write (lupri,*) '- count,source,tag: ',BUF_count,source,tag
         write (lupri,*) '- calling MPI_MYFAIL ...'
         CALL MPI_MYFAIL(ierr_mpi)
      END IF
#endif
      END

      subroutine mpixprocname(process_name,process_name_length)
C
#ifdef VAR_MPI
      implicit none
      character*(*) process_name
      integer       process_name_length
      include 'mpif.h'
      integer(kind=MPI_INTEGER_KIND) name_length_mpi, ierr_mpi
C
      call mpi_get_processor_name(process_name,name_length_mpi,ierr_mpi)
      process_name_length = name_length_mpi
#endif
      END

#if defined (VAR_MPI)
C  /* Deck mpi_myfail */
      SUBROUTINE MPI_MYFAIL(IERR)
! only called from other routines in mpimacro.F
C
      implicit none
#include "priunit.h"
      INCLUDE 'mpif.h'

      integer(kind=MPI_INTEGER_KIND) ierr,ierrcl,ierr2
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
