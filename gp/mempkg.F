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
#ifdef UNDEF
!===========================================================================
!/* Comdeck memlog */
!060316-hjaaj: for DBG_LENID check memcb1.h
!910110-hjaaj MEMREL: removed non-standard char. concat. in CALL MEMCHK
!901004-hjaaj: NWNMEM introduced, no final check on NWNMEM (yet)
!===========================================================================
#endif
C  /* Deck memini */
      SUBROUTINE MEMINI(LUWME1,LUEME1)
C
C Copyright 8-Dec-1989,4-Oct-1990 Hans Joergen Aa. Jensen
C
C     Define output units for memory allocation routines
C     Initialize mempkg warning count NWNMEM to zero
C
#include <implicit.h>
#include "memcb1.h"
C
      LUWMEM = LUWME1
      LUEMEM = LUEME1
      IDENT8 = '** KFREE'
      READ (IDENT8,'(A8)') WIDENT_KFREE
      RETURN
      END
C  /* Deck memget */
      SUBROUTINE MEMGET(TYPE,KBASE,LENGTH,WORK,KFREE,LFREE)
C
C Copyright 9-Jan-1988/27-Jul-1993 Hans Joergen Aa. Jensen
C
C l.r. 980824-hjaaj: QENTER/QEXIT calls use too much time because MEMGET
C is called often. Now only call QENTER/QEXIT if warning or error.
C Restructured code (e.g. changed errors from IF () THEN to IF () GO TO)
C to get fewer logical tests and fewer jumps after logical tests.
C
C     Memory allocation
C
#include <implicit.h>
      CHARACTER*(*) TYPE
      DIMENSION     WORK(*)
      CALL MEMGET2(TYPE,'-??????-',KBASE,LENGTH,WORK,KFREE,LFREE)
      RETURN
      END
C  /* Deck memget */
      SUBROUTINE MEMGET2(TYPE,IDENT,KBASE,LENGTH,WORK,KFREE,LFREE)
C
C Copyright 9-Jan-1988/27-Jul-1993 Hans Joergen Aa. Jensen
C
C l.r. 980824-hjaaj: QENTER/QEXIT calls use too much time because MEMGET
C is called often. Now only call QENTER/QEXIT if warning or error.
C Restructured code (e.g. changed errors from IF () THEN to IF () GO TO)
C to get fewer logical tests and fewer jumps after logical tests.
C
C     Memory allocation
C
#include <implicit.h>
      CHARACTER*(*) TYPE, IDENT
      DIMENSION     WORK(*)
C
#include <iratdef.h>
#include "memcb1.h"
C
Chj1  CALL QENTER('MEMGET2 ')
      IF (TYPE(1:4) .EQ. 'REAL' .OR. TYPE(1:4) .EQ. 'INT8') THEN
         LREAL = LENGTH
      ELSE IF (TYPE(1:4) .EQ. 'COMP') THEN
         LREAL = 2*LENGTH
      ELSE IF (TYPE(1:4) .EQ. 'INTE') THEN
         LREAL = (LENGTH-1)/IRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'INT4') THEN
         LREAL = (LENGTH-1)/2 + 1
      ELSE IF (TYPE(1:4) .EQ. 'INT2') THEN
         LREAL = (LENGTH-1)/4 + 1
      ELSE IF (TYPE(1:4) .EQ. 'LOGI') THEN
         LREAL = (LENGTH-1)/LRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'WORK') THEN
         LENGTH = LFREE - 2*LENID
C                         ^- make sure OK for KFREE .eq. 1
         LREAL  = LENGTH
      ELSE
         CALL QENTER('MEMGET2 ')
         WRITE (LUWMEM,'(/2A/2A)')
     *      ' MEMGET2 ERROR, illegal type : ',TYPE,
     *      ' MEMGET2 ERROR, for ident    : ',IDENT
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,'(/2A/2A)')
     *      ' MEMGET2 ERROR, illegal type : ',TYPE,
     *      ' MEMGET2 ERROR, for ident    : ',IDENT
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET2: programming ERROR, illegal type')
      END IF
C
      IDENT8 = IDENT
C
      IF (KFREE .NE. 1) THEN

C        length error check:
         IF (LREAL+LENID .GT. LFREE) GO TO 8001

C        Error check, memory check of WORK(KFREE) OK ? :

         WMEMCK = WORK(KFREE-1)
         IF (MEMCK .NE. MEMID) GO TO 8002
C        ... comparing integers, to be sure no renormalization of float
#if DBG_LENID
         WMEMCK = WORK(KFREE-LENID+2)
         IF (MEMCK .NE. MEMID) GO TO 8002
#endif

C        Warning check, is length of WORK(KFREE) .ne. 0 ? :

         WMEMCK = WORK(KFREE-LENID+1)
         IF (MEMCK .NE. 0) GO TO 8003
C        ... go issue warning and continue at 7003
      ELSE
C        KFREE = 1 case, initialize memory vector info:
         IF (LREAL+2*LENID .GT. LFREE) GO TO 8001
#if DBG_LENID
         do i = 3,lenid
           WORK(i) = WMEMID
         end do
#else
         WORK(3) = WMEMID
#endif
         KFREE = KFREE + LENID
         LFREE = LFREE - LENID
      END IF

C     save ident and length for new allocation

 7003 CONTINUE
      READ (IDENT8,'(A8)') WORK(KFREE-LENID)
      WORK(KFREE-LENID+1) = WLREAL

      KBASE  = KFREE
      KFREE  = KBASE + (LREAL + LENID)
      LFREE  = LFREE - (LREAL + LENID)
      WORK(KFREE-LENID)   = WIDENT_KFREE
      WORK(KFREE-LENID+1) = WLZERO
#ifdef DBG_LENID
      do i = 2,lenid-1
         2ORK(KFREE-LENID+i) = WMEMID
      end do
#else
      WORK(KFREE-1) = WMEMID
#endif
Chj1  CALL QEXIT('MEMGET2 ') ! too expensive as MEMGET is called many times
      RETURN
C
C     error branches:
C
C        error branch  IF (LREAL+LENIDS .GT. LFREE) :
C
 8001    CALL QENTER('MEMGET2 ')
         IF (KFREE .GT. 1) THEN
            LENIDS = LENID
            WRITE (LUWMEM,1010) LREAL+LENIDS,LFREE
            CALL MEMCHK('MEMGET ERROR (insuff. memory)',WORK,1)
         ELSE
            LENIDS = 2*LENID
            WRITE (LUWMEM,1010) LREAL+LENIDS,LFREE
         END IF
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,1010) LREAL+LENIDS,LFREE
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET ERROR, insufficient work space in memory')
 1010 FORMAT(
     &   /' MEMGET ERROR, insufficient free space for next allocation',
     &   /T16,'( Need:',I10,', available (LFREE):',I10,' )')
C
C        error branch for invalid MEMID for WORK(KFREE)
C
 8002    CALL QENTER('MEMGET2 ')
         WRITE (LUWMEM,'(/A/A,I12,2(/A,I25))')
     *      ' MEMGET ERROR, not a valid memget id in work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory checks:',MEMCK,
     *      ' expected           :',MEMID
         CALL MEMCHK('MEMGET call of MEMCHK:',WORK,1)
         CALL QUIT('MEMGET ERROR, not a valid memget '//
     &             'id in work(kfree)')
C
C        warning branch for non-zero length for WORK(KFREE):
C
 8003    NWNMEM = NWNMEM + 1
         WRITE (LUWMEM,'(/A/A,I10,/A,I25)')
     *      ' MEMGET WARNING, nonzero allocation work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory allocation :',MEMCK
         CALL QENTER('MEMGET2 ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMGET WARNING, nonzero allocation work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QEXIT('MEMGET2 ')
C        Warning issued, continue and do the assignment now
         GO TO 7003
C
C     end of MEMGET2
C
      END
C  /* Deck memrel */
      SUBROUTINE MEMREL(TEXT,WORK,KFIRST,KREL,KFREE,LFREE)
C
C 14-Mar-1989 Hans Joergen Aa. Jensen
C
C Check memory allocation from WORK(KFIRST)
C then free memory from WORK(KREL)
C
#include <implicit.h>
      DIMENSION     WORK(*)
      CHARACTER*(*) TEXT
C
#include "memcb1.h"
C
      CHARACTER     TMPTXT*71
C
      CALL QENTER('MEMREL  ')
C
C     First check if memory allocation is intact
C     (check if anything has been out of bounds)
C
      LENTXT = LEN(TEXT)
      LENTXT = MIN(LENTXT,50)
      TMPTXT = TEXT(1:LENTXT) // ' (called from MEMREL)'
      LENTXT = LENTXT + 21
      CALL MEMCHK(TMPTXT(1:LENTXT),WORK,KFIRST)
C
C     Check if WORK(krel) has legal memid
C
      IF (KREL .GT. KFREE) GO TO 8001
      IF (KREL .NE. 1) THEN
         KFREE1 = KREL
      ELSE
         KFREE1 = 1 + LENID
      END IF

      WMEMCK = WORK(KFREE1-1)
      IF (MEMCK .NE. MEMID) GO TO 8002
C
C     release ...
C
      WORK(KFREE1-LENID  ) = WIDENT_KFREE
      WORK(KFREE1-LENID+1) = WLZERO
      LFREE = LFREE + (KFREE-KREL)
      KFREE = KREL
C
      CALL QEXIT('MEMREL  ')
      RETURN
C
C     error branch IF (KREL .GT. KFREE) THEN
C
 8001    WRITE (LUWMEM,'(/A/2A/A,3I10)')
     *      ' MEMREL ERROR, krel .gt. kfree',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KREL,KFREE =',KFIRST,KREL,KFREE
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMREL ERROR, krel .gt. kfree',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMREL ERROR, krel .gt. kfree')
C
C     error branch IF (MEMCK .NE. MEMID) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,2I10,2(/A,2I25))')
     *      ' MEMREL ERROR, not a valid memget id in work(krel)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KREL =',KFIRST,KREL,
     *      ' found memory check at work(krel):',MEMCK,
     *      ' expected                        :',MEMID
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMREL ERROR, not a valid memget id in work(krel)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMREL ERROR, not a valid memget id in work(krel)')
C
C     end of MEMREL
C
      END
C  /* Deck memchk */
      SUBROUTINE MEMCHK(TEXT,WORK,KFIRST)
C
C 17-Mar-1989 Hans Joergen Aa. Jensen
C
C Check if memory allocation from WORK(KFIRST) is intact
C (check if anything has been out of bounds).
C
#include <implicit.h>
      DIMENSION     WORK(*)
      CHARACTER*(*) TEXT
C
#include "memcb1.h"
C
Chj1  CALL QENTER('MEMCHK  ')
C
      KALLOC = KFIRST
      IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
C
      IALLOC = 0
  100 CONTINUE
         IALLOC = IALLOC + 1
#ifdef DBG_LENID
         do i = 2,lenid-1
            WMEMCK = WORK(KALLOC-LENID+i)
            IF (MEMCK .NE. MEMID) GO TO 8001
         end do
#else
         WMEMCK = WORK(KALLOC-1)
         IF (MEMCK .NE. MEMID) GO TO 8001
#endif
         WLREAL = WORK(KALLOC-LENID+1)
         IF (LREAL .LT. 0) GO TO 8002
         IF (LREAL .GT. 0) THEN
            KALLOC = KALLOC + LENID + LREAL
            GO TO 100
         END IF
C
Chj1  CALL QEXIT('MEMCHK  ')
      RETURN
C
C        error branch IF (MEMCK .NE. MEMID) THEN
C
 8001    WRITE (LUWMEM,'(/A/2A/A,3I10,2(/A,I25))')
     *      ' MEMCHK ERROR, not a valid memget id in work(kalloc-1)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KALLOC,IALLOC =',KFIRST,KALLOC,IALLOC,
     *      ' found memory checks:',MEMCK,
     *      ' expected           :',MEMID
         WRITE (LUWMEM,'(/A/A)')
     &      ' Dump of allocations before invalid memget id:',
     &      ' IALLOC #    IDENT     KBASE         Length'
         KALLOC = KFIRST
         IF (KALLOC .EQ. 1) KALLOC = 1 + LENID
         DO I = 1, IALLOC
            WRITE (IDENT8,'(A8)') WORK(KALLOC-LENID)
            WLREAL = WORK(KALLOC-LENID+1)
            WRITE (LUWMEM,'(I5,2X,A8,2I15)') I,IDENT8,KALLOC,LREAL
            KALLOC = KALLOC + LENID + LREAL
         END DO
#ifdef DBG_LENID
         do i = 2,lenid-1
            WMEMCK(1) = WORK(KALLOC-LENID+i)
            if (memck(1) .ne. memid(1)) write (luwmem,*)
     &         'memchk error in debug element: ',i+1,memck(1)
         end do
#endif
         CALL QENTER('MEMCHK  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMCHK ERROR, not a valid memget id in work(kalloc-1)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMCHK ERROR, not a valid memget '//
     &             'id in work(kalloc-1)')
C
C        error branch IF (LREAL .LT. 0) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,4I10)')
     *      ' MEMCHK ERROR, next allocation LENGTH is negative',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KALLOC,IALLOC,LENGTH =',KFIRST,KALLOC,IALLOC,LREAL
         CALL QENTER('MEMCHK  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMCHK ERROR, next allocation LENGTH is negative',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMCHK ERROR, next '//
     &             'allocation has negative length.')
C
C     end of MEMCHK
C
      END
