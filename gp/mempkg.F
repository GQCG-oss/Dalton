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
#ifdef UNDEF
!===========================================================================
!/* Comdeck memlog */
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
#include <memcb1.h>
C
      LUWMEM = LUWME1
      LUEMEM = LUEME1
      NWNMEM = 0
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
C
#include <iratdef.h>
#include <memcb1.h>
C
      PARAMETER ( LENID = 1/IRAT + 1)
      DIMENSION     WMEMID(2), WMEMCK(2)
      INTEGER       MEMID(2),  MEMCK(2)
      EQUIVALENCE  (MEMID, WMEMID), (MEMCK, WMEMCK)
      SAVE          MEMID
      DATA          MEMID(1) /1234567890/
C
Chj1  CALL QENTER('MEMGET  ')
      IF (TYPE(1:4) .EQ. 'REAL' .OR. TYPE(1:4) .EQ. 'INT8') THEN
         LREAL = LENGTH
      ELSE IF (TYPE(1:4) .EQ. 'COMP') THEN
         LREAL = 2*LENGTH
      ELSE IF (TYPE(1:4) .EQ. 'INTE' .OR. TYPE(1:4) .EQ. 'INT4') THEN
         LREAL = (LENGTH-1)/IRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'LOGI') THEN
         LREAL = (LENGTH-1)/LRAT + 1
      ELSE IF (TYPE(1:4) .EQ. 'INT2') THEN
         LREAL = (LENGTH-1)/IRAT2 + 1
      ELSE IF (TYPE(1:4) .EQ. 'WORK') THEN
         LENGTH = LFREE - 2*LENID
C                         ^- make sure OK for KFREE .eq. 1
         LREAL  = LENGTH
      ELSE
         CALL QENTER('MEMGET  ')
         WRITE (LUWMEM,'(/2A)')
     *      ' MEMGET ERROR, illegal type :',TYPE
         CALL QTRACE(LUWMEM)
         IF (LUEMEM .NE. LUWMEM) THEN
            WRITE (LUEMEM,'(/2A)')
     *         ' MEMGET ERROR, illegal type :',TYPE
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET ERROR, illegal type')
      END IF
C
C
      MEMID(2) = LREAL
      IF (KFREE .NE. 1) THEN
C
C        error and warning checks:
C
         WMEMCK(1) = WORK(KFREE-LENID)
         IF (LENID .EQ. 2) WMEMCK(2) = WORK(KFREE-1)
C        error check:
         IF (LREAL+LENID .GT. LFREE) GO TO 8001
         IF (MEMCK(1) .NE. MEMID(1)) GO TO 8002
C        warning check:
         IF (MEMCK(2) .NE. 0) GO TO 8003
C        ... go issue warning and continue at 8013
C
         WORK(KFREE-LENID) = WMEMID(1)
         IF (LENID .EQ. 2) WORK(KFREE-1) = WMEMID(2)
      ELSE
         IF (LREAL+2*LENID .GT. LFREE) GO TO 8001
         WORK(KFREE) = WMEMID(1)
         IF (LENID .EQ. 2) WORK(KFREE+1) = WMEMID(2)
         KFREE = KFREE + LENID
         LFREE = LFREE - LENID
      END IF
 8013 KBASE  = KFREE
      KFREE  = KBASE + (LREAL + LENID)
      LFREE  = LFREE - (LREAL + LENID)
      MEMID(2) = 0
      WORK(KFREE-LENID) = WMEMID(1)
      IF (LENID .EQ. 2) WORK(KFREE-1) = WMEMID(2)
Chj1  CALL QEXIT('MEMGET  ')
      RETURN
C
C     error branches:
C
C        error branch  IF (LREAL+LENIDS .GT. LFREE) :
C
 8001    CALL QENTER('MEMGET  ')
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
C        error branch IF (MEMCK(1) .NE. MEMID(1)) :
C
 8002    CALL QENTER('MEMGET  ')
         WRITE (LUWMEM,'(/A/A,I10,2(/A,I25))')
     *      ' MEMGET ERROR, not a valid memget id in work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory check :',MEMCK(1),
     *      ' expected           :',MEMID(1)
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMGET ERROR, not a valid memget id in work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMGET ERROR, not a valid memget '//
     &             'id in work(kfree)')
C
C        warning branch IF (MEMCK(2) .NE. 0) :
C
 8003    NWNMEM = NWNMEM + 1
         WRITE (LUWMEM,'(/A/A,I10,/A,I25)')
     *      ' MEMGET WARNING, nonzero allocation work(kfree)',
     *      ' KFREE =',KFREE,
     *      ' found memory allocation :',MEMCK(2)
         CALL QENTER('MEMGET  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMGET WARNING, nonzero allocation work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QEXIT('MEMGET  ')
C        Warning issued, do the assignment now and continue
         WORK(KFREE-LENID) = WMEMID(1)
         IF (LENID .EQ. 2) WORK(KFREE-1) = WMEMID(2)
         GO TO 8013
C
C     end of MEMGET
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
#include <iratdef.h>
C
#include <memcb1.h>
C
      CHARACTER     TMPTXT*71
      PARAMETER ( LENID = 1/IRAT + 1)
      DIMENSION     WMEMID(2), WMEMCK(2)
      INTEGER       MEMID(2),  MEMCK(2)
      EQUIVALENCE  (MEMID, WMEMID), (MEMCK, WMEMCK)
      SAVE          MEMID
      DATA          MEMID(1) /1234567890/
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
      WMEMCK(1) = WORK(KFREE1-LENID)
      IF (MEMCK(1) .NE. MEMID(1)) GO TO 8002
C
C     release ...
C
      MEMID(2) = 0
      WORK(KFREE1-LENID) = WMEMID(1)
      IF (LENID .EQ. 2) WORK(KFREE1-1) = WMEMID(2)
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
C     error branch IF (MEMCK(1) .NE. MEMID(1)) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,2I10,2(/A,I25))')
     *      ' MEMREL ERROR, not a valid memget id in work(krel)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KREL =',KFIRST,KREL,
     *      ' found memory check :',MEMCK(1),
     *      ' expected           :',MEMID(1)
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMREL ERROR, not a valid memget id in work(kfree)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMREL ERROR, not a valid memget id in work(kfree)')
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
#include <iratdef.h>
C
#include <memcb1.h>
C
      PARAMETER ( LENID = 1/IRAT + 1)
      DIMENSION     WMEMID(2), WMEMCK(2)
      INTEGER       MEMID(2),  MEMCK(2)
      EQUIVALENCE  (MEMID, WMEMID), (MEMCK, WMEMCK)
      SAVE          MEMID
      DATA          MEMID(1) /1234567890/
C
Chj1  CALL QENTER('MEMCHK  ')
C
      KFREE = KFIRST
      IF (KFREE .EQ. 1) KFREE = 1 + LENID
C
      IALLOC = 0
  100 CONTINUE
         IALLOC = IALLOC + 1
         WMEMCK(1) = WORK(KFREE-LENID)
         IF (LENID .EQ. 2) WMEMCK(2) = WORK(KFREE-1)
         IF (MEMCK(1) .NE. MEMID(1)) GO TO 8001
         LREAL = MEMCK(2)
         IF (LREAL .LT. 0) GO TO 8002
         IF (LREAL .GT. 0) THEN
            KFREE = KFREE + LENID + LREAL
            GO TO 100
         END IF
C
Chj1  CALL QEXIT('MEMCHK  ')
      RETURN
C
C        error branch IF (MEMCK(1) .NE. MEMID(1)) THEN
C
 8001    WRITE (LUWMEM,'(/A/2A/A,3I10,2(/A,I25))')
     *      ' MEMCHK ERROR, not a valid memget id in work(kfree-1)',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KFREE,IALLOC =',KFIRST,KFREE,IALLOC,
     *      ' found memory check :',MEMCK(1),
     *      ' expected           :',MEMID(1)
         WRITE (LUWMEM,'(/A/A)')
     &      ' Dump of allocations before invalid memget id:',
     &      ' IALLOC #    KBASE         Length'
         KFREE = KFIRST
         IF (KFREE .EQ. 1) KFREE = 1 + LENID
         DO I = 1, IALLOC-1
            WMEMCK(1) = WORK(KFREE-LENID)
            IF (LENID .EQ. 2) WMEMCK(2) = WORK(KFREE-1)
            LREAL = MEMCK(2)
            WRITE (LUWMEM,'(I5,2I15)') I,KFREE,LREAL
            KFREE = KFREE + LENID + LREAL
         END DO
         CALL QENTER('MEMCHK  ')
         CALL QTRACE(LUWMEM)
         IF (LUEMEM.NE.LUWMEM) THEN
            WRITE (LUEMEM,'(/A/A)')
     *      ' MEMCHK ERROR, not a valid memget id in work(kfree-1)',
     *      ' ---> see output file'
            CALL QTRACE(LUEMEM)
         END IF
         CALL QUIT('MEMCHK ERROR, not a valid memget '//
     &             'id in work(kfree-1)')
C
C        error branch IF (LREAL .LT. 0) THEN
C
 8002    WRITE (LUWMEM,'(/A/2A/A,4I10)')
     *      ' MEMCHK ERROR, next allocation LENGTH is negative',
     *      ' Text from calling routine : ',TEXT,
     *      ' KFIRST,KFREE,IALLOC,LENGTH =',KFIRST,KFREE,IALLOC,LREAL
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
