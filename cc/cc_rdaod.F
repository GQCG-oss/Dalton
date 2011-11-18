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
*======================================================================*
C  /* Deck cc_rdaod */
      SUBROUTINE CC_RDAOD(XINT,JSCOOR,JCOOR,JCORSY,IGAM,IDEL,NUMG,NUMD,
     *                    WORK,LWORK,IRECNR,NRECS,MXBUF,NFILES,
     *                    DIRINT,ITYPE,MXCOMP,LDERINT)
*----------------------------------------------------------------------*
C
C     Purpose: Read (**|GD) distribution of (derivative) AO integrals
C              for a given symmetry coordinate (JCOOR,JCORSY)
C
C     direct option (DIRINT=.TRUE.) disabled, because the derivative
C     integrals are only available non-direct (for the time beeing)
C
C     --> WORK is not used, IRECNR, LWORK are dummy
C
C     ITYPE = 0   --  two-electron integrals of unperturbed integrals
C     ITYPE = 1   --  geometric first derivatives of 2-el integrals
C     ITYPE = 5   --  magnetic first derivatives of 2-el integrals
C
C
C     Written by Christof Haettig 06-May-1998
C     based on Henrik Kochs CC_RDAO routine
C     magnetic derivatives introduced 20-Sep-1999
C
*----------------------------------------------------------------------*
#include <implicit.h>
C
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxash.h>
#include <iratdef.h>
#include <ccorb.h>
CCN #include <infind.h> replaced by #include <ccisao.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <eribuf.h>
#include <nuclei.h>
#include <chrnos.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      LOGICAL   DIRINT
      LOGICAL   LDERINT(8,*)
      INTEGER JSCOOR, JCOOR, JCORSY, MXCOMP, NFILES
      DIMENSION XINT(*),WORK(LWORK)
      INTEGER IGAM(NUMG),IDEL(NUMD)
      INTEGER IRECNR(MXBUF,0:NFILES-1), NRECS(0:JSCOOR)
C
      LOGICAL OLDDX, LONDON2
      INTEGER LENGTH(8), NDIMAB(8), NDIMDIS(8),IDSAO(8,8)
      INTEGER KSCOOR, KCOOR, KCORSY, MAXCMP, MCSYM, UNIAO2
      CHARACTER*8 NAME(16)
      CHARACTER*8 FNAME, FAODER
C
C
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8',
     *            'CCAODER1','CCAODER2','CCAODER3','CCAODER4',
     *            'CCAODER5','CCAODER6','CCAODER7','CCAODER8'/
      COMMON/SORTIO/LUAOIN(8)
C
      LOGICAL LDUMMY
 
*----------------------------------------------------------------------*
*     check ITYPE option and set some integers:
*       KCOOR  : required coordinate 
*       KCORSY : required symmetry class
*       KSCOOR : symmetry coordinate (used as file index)
*       MAXCMP : maximum # of coordinates
*       MCSYM  : maximum symmetry class 
*     for ITYPE = 0 the loops over the coordinates/symmetries are 
*     suppressed --> all for integers are set to 1
*----------------------------------------------------------------------*
      IF (ITYPE.EQ.0) THEN
        KCOOR  = 1  
        KCORSY = 1
        KSCOOR = 0  ! direct undiff. integrals are on file nb. 0
        MAXCMP = 1
        MCSYM  = 1
      ELSE IF (ITYPE.EQ.1) THEN
        KCOOR  = JCOOR
        KCORSY = JCORSY
        IF (DIRINT) THEN
           KSCOOR = JSCOOR 
        ELSE
           KSCOOR = JCOOR  ! symmetry coordinate = coordinate index
        END IF 
        MAXCMP = MXCOMP
        MCSYM  = NSYM
      ELSE IF (ITYPE.EQ.5) THEN
        KCOOR  = JCOOR
        KCORSY = JCORSY
        KSCOOR = JSCOOR
        MAXCMP = MXCOMP
        MCSYM  = NSYM
      ELSE
        WRITE (LUPRI,*) 'CC_RDAOD called with illegal value for ITYPE:',
     &        ITYPE
        CALL QUIT('CC_RDAOD called with illegal value for ITYPE.')
      END IF
C
C     set length of alpha,beta (2 index), and length and offsets of
C     alpha,beta,gamma (3 index) subblocks:
C
      IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
         DO ISYM = 1, NSYM
           NDIMAB(ISYM)  = NNBST(ISYM)
           NDIMDIS(ISYM) = NDISAO(ISYM)
           DO ISYMI = 1, NSYM
             IDSAO(ISYM,ISYMI) = IDSAOG(ISYM,ISYMI)
           END DO
         END DO
      ELSE IF (ITYPE.EQ.5) THEN
         DO ISYM = 1, NSYM
           NDIMAB(ISYM)  = N2BST(ISYM)
           NDIMDIS(ISYM) = NDISAOSQ(ISYM)
           DO ISYMI = 1, NSYM
             IDSAO(ISYM,ISYMI) = IDSAOGSQ(ISYM,ISYMI)
           END DO
         END DO
      ELSE
         CALL QUIT('Unknown ITYPE in CC_RDAOD.')
      END IF
C
C     precalculate length of distributions:
C
      DO ISYMD = 1, NSYM
        LENGTH(ISYMD) = 0
        DO ICOOR = 1, MAXCMP
          DO ICORSY = 1, MCSYM
            IF ( ITYPE.EQ.0 .OR. LDERINT(ICORSY,ICOOR) ) THEN
              ISYDIS = MULD2H(ISYMD,ICORSY)
              LENGTH(ISYMD) = LENGTH(ISYMD) + NDIMDIS(ISYDIS)
            END IF
          END DO
        END DO
      END DO

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'entered CC_RDAOD> DIRINT = ',DIRINT
        WRITE(LUPRI,*) 'entered CC_RDAOD> ITYPE,KSCOOR :',ITYPE,KSCOOR
        WRITE(LUPRI,*) 'entered CC_RDAOD> KCORSY,KCOOR :',KCORSY,KCOOR
      END IF

*----------------------------------------------------------------------*
*     Non-direct first derivative integrals.
*----------------------------------------------------------------------*

      IF (.NOT. DIRINT) THEN

         IF ( .NOT. (ITYPE.EQ.0 .OR. LDERINT(KCORSY,KCOOR)) ) THEN
C           WRITE (LUPRI,*) 
C    *            'CC_RDAOD> Warning: required integrals are zero...'
C           WRITE (LUPRI,*) 
C    *            'CC_RDAOD> ITYPE        = ',ITYPE
C           WRITE (LUPRI,*) 
C    *            'CC_RDAOD> KCORSY/KCOOR = ',KCORSY,KCOOR
C           WRITE (LUPRI,*) 
C    *            'CC_RDAOD> LDERINT      = ',LDERINT(KCORSY,KCOOR)

            KOFF = 1
            DO ID = 1, NUMD
            DO IG = 1, NUMG
              ISYMD  = ISAO(IDEL(ID))
              ISYMG  = ISAO(IGAM(IG))
              ISYDIS = MULD2H(ISYMD,KCORSY)
              ISYMAB = MULD2H(ISYDIS,ISYMG)                
              CALL DZERO(XINT(KOFF),NDIMAB(ISYMAB))
              KOFF = KOFF + NDIMAB(ISYMAB)
            END DO
            END DO
           
            RETURN
         END IF

         KOFF = 1
         DO ID = 1,NUMD
C
            ISYMD  = ISAO(IDEL(ID))
            D      = IDEL(ID) - IBAS(ISYMD)
            IOFF0  = LENGTH(ISYMD)*(D-1)

            IF (D.LT.0 .OR. D.GT.NBAST) THEN
              WRITE (LUPRI,*) 'ORBITAL INDEX OUT OF RANGE IN CC_RDAOD.'
              WRITE (LUPRI,*) 
     &             'CC_RDAOD: read distribution IDEL(ID) = ',IDEL(ID)
              WRITE (LUPRI,*) 'CC_RDAOD: ISYMD, D, IOFF0 = ',ISYMD,D,
     &             IOFF0
              CALL QUIT('ORBITAL INDEX OUT OF RANGE IN CC_RDAOD.')
            END IF 

            IF (ITYPE.EQ.0) THEN
              NFILE  = LUAOIN(ISYMD)
              FNAME  = NAME(ISYMD)
            ELSE IF (ITYPE.EQ.1 .OR. ITYPE.EQ.5) THEN
              NFILE  = 0
              FNAME  = NAME(8+ISYMD)
              CALL WOPEN2(NFILE,FNAME,64,0)
            END IF

            IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'CC_RDAOD: FNAME = ', FNAME
               WRITE (LUPRI,*) 'CC_RDAOD: NFILE = ', NFILE
               WRITE (LUPRI,*) 'CC_RDAOD: ITYPE = ', ITYPE
               WRITE (LUPRI,*) 'CC_RDAOD: LENGTH= ', LENGTH(ISYMD)
            END IF


            DO ICOOR  = 1, MXCOMP
            DO ICORSY = 1, NSYM

               ISYDIS = MULD2H(ISYMD,ICORSY)

               IF ( ITYPE.EQ.0 .OR. LDERINT(ICORSY,ICOOR) ) THEN
C
                  IF ( ICOOR.EQ.KCOOR .AND. ICORSY.EQ.KCORSY ) THEN

                     DO IG = 1,NUMG
                        ISYMG  = ISAO(IGAM(IG))
                        G      = IGAM(IG) - IBAS(ISYMG)
                        ISYMAB = MULD2H(ISYDIS,ISYMG)                
 
                        IOFF = IOFF0 + IDSAO(ISYMG,ISYDIS)
     *                               + NDIMAB(ISYMAB)*(G - 1) + 1
                        

                        CALL GETWA2 ( NFILE, FNAME, XINT(KOFF), IOFF,
     *                                NDIMAB(ISYMAB))
 
                        IF (LOCDBG) THEN
                          XNORM=DDOT(NDIMAB(ISYMAB),XINT(KOFF),1,
     *                                              XINT(KOFF),1)
                          WRITE (LUPRI,'(a,2i5,f10.5)') 'CC_RDAOD> ', 
     *                                       IDEL(ID),IGAM(IG),XNORM
                          CALL OUTPUT(XINT(KOFF),1,NDIMAB(ISYMAB),1,1,
     *                                     NDIMAB(ISYMAB),1,1,LUPRI)
                        END IF
 
                        KOFF = KOFF + NDIMAB(ISYMAB)
                    END DO

                  END IF 

                  IOFF0  = IOFF0 + NDIMDIS(ISYDIS)

               END IF 
            END DO
            END DO

            IF (ITYPE.EQ.1 .OR. ITYPE.EQ.5) THEN
              CALL WCLOSE2(NFILE,NAME(ISYM),'KEEP')
            END IF

         END DO
 
      ELSE
C
C        record length
C
C
         LBFINP = LBUF
C
#if defined (VAR_INT64)
         LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
         LRECL = LBFINP + LBFINP/IRAT + 1
#else
         IF (NBASIS .LE. 255) THEN
            LRECL = 2*LBFINP + LBFINP + 1
         ELSE
            LRECL = 2*LBFINP + 2*LBFINP + 1
         END IF
#endif
C
         IF (LOCDBG) THEN
           WRITE(LUPRI,*) 'LBFINP,NBASIS,LRECL:',LBFINP,NBASIS,LRECL
         END IF
C
C        Buffer allocation
C
         IBIT1 = 2**NBITS     - 1
         IBIT2 = 2**(2*NBITS) - 1
         LENINT = IRAT*LBUF + NIBUF*LBUF + 1
         LRINT  = IRAT*LBUF

         KAOAB = 1
         KDIST = KAOAB + (NBAST*NBAST+1)/IRAT + 1
         KIBUF = KDIST + (NUMD*NBAST+1)/IRAT + 1
         KRBUF = KIBUF + (NIBUF*LBUF-1)/IRAT + 1
         KEND1 = KRBUF + LBUF
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Insufficient work space in CC_RDAOD:'
            WRITE(LUPRI,*) 'Need:',KEND1
            WRITE(LUPRI,*) 'Available:',LWORK
            CALL QUIT('Insufficient work space in CC_RDAOD')
         ENDIF
       
         CALL CCSD_INIT2B(WORK(KAOAB),DUMMY,DUMMY,ITYPE,
     &                    .TRUE.,1,LDUMMY)

*----------------------------------------------------------------------*
*        open file, read integral distributions and close file again:
*----------------------------------------------------------------------*
         !SONIA, changed to -1, UNIAO2 = 0
         UNIAO2 = -1
         FAODER = 'AO2DIS'//CHRNOS(KSCOOR/10)
     &                    //CHRNOS(MOD(KSCOOR,10))
         CALL GPOPEN(UNIAO2,FAODER,'OLD','DIRECT',
     &               'UNFORMATTED',LRECL,OLDDX)
         IF (LOCDBG) WRITE(LUPRI,*) 'CC_RDAOD> opened file:',FAODER
C
C        read integrals and close file
C
         LONDON2  = .FALSE.
         CALL CCRDAOD1(XINT,WORK(KIBUF),WORK(KRBUF),UNIAO2,LONDON2,
     &                 KCORSY,IGAM,IDEL,NUMG,NUMD,ITYPE,
     &                 WORK(KAOAB),NDIMAB,WORK(KDIST),IBIT1,
     &                 IRECNR(1,KSCOOR),NRECS(KSCOOR))

         CALL GPCLOSE(UNIAO2,'KEEP')
 
*----------------------------------------------------------------------*
*        for London integrals read second half of the integral matrices:
*----------------------------------------------------------------------*
         IF (ITYPE.EQ.5) THEN
           LONDON2 = .TRUE.
           KSCOOR = KSCOOR + 3

           !SONIA, changed to -1, UNIAO2 = 0
           UNIAO2 = -1
           FAODER = 'AO2DIS'//CHRNOS(KSCOOR/10)
     &                      //CHRNOS(MOD(KSCOOR,10))
           CALL GPOPEN(UNIAO2,FAODER,'OLD','DIRECT',
     &                 'UNFORMATTED',LRECL,OLDDX)
 
           CALL CCRDAOD1(XINT,WORK(KIBUF),WORK(KRBUF),UNIAO2,LONDON2,
     &                   KCORSY,IGAM,IDEL,NUMG,NUMD,ITYPE,
     &                   WORK(KAOAB),NDIMAB,WORK(KDIST),IBIT1,
     &                   IRECNR(1,KSCOOR),NRECS(KSCOOR))

           CALL GPCLOSE(UNIAO2,'KEEP')
         END IF

      ENDIF

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'leaving CC_RDAOD> '
      END IF

      RETURN
      END
*======================================================================*
*======================================================================*
C  /* Deck ccrdaod1 */
      SUBROUTINE CCRDAOD1(XINT,IBUF,RBUF,UNIAO2,LONDON2,
     *                    ICORSY,IGAM,IDEL,NUMG,NUMD,ITYPE,
     *                    IADRPQ,NDIMPQ,IOFFRD,IBIT1,IRECNR,NRECORDS)
*----------------------------------------------------------------------*
C
C     Purpose: Read (**|GD) distribution of (derivative) AO integrals
C              for a given symmetry coordinate (JCOOR,JCORSY)
C
C     direct option (DIRINT=.TRUE.) disabled, because the derivative
C     integrals are only available non-direct (for the time beeing)
C
C     --> WORK is not used, IRECNR, LWORK are dummy
C
C     ITYPE = 0   --  two-electron integrals of unperturbed integrals
C     ITYPE = 1   --  geometric first derivatives of 2-el integrals
C     ITYPE = 5   --  magnetic first derivatives of 2-el integrals
C
C
C     Written by Christof Haettig 06-May-1998
C     based on Henrik Kochs CC_RDAO routine
C     magnetic derivatives introduced 20-Sep-1999
C
*----------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <eribuf.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ibtpar.h>
C
      LOGICAL LOCDBG
      PARAMETER ( LOCDBG = .FALSE. )
C
      DIMENSION XINT(*), RBUF(LBUF)
      LOGICAL LONDON2
      INTEGER IBUF(LBUF,NIBUF)
      INTEGER UNIAO2,ICORSY,NUMG,NUMD,ITYPE,IBIT1
      INTEGER IGAM(NUMG),IDEL(NUMD),IRECNR(*), NRECORDS
      INTEGER IADRPQ(NBAST,NBAST), IOFFRD(NBAST,NUMD)
      INTEGER NDIMPQ(NSYM)
C
C functions:
#include <ibtfun.h>
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'entered CCRDAOD1'
        WRITE(LUPRI,*) 'UNIAO2:',UNIAO2
        WRITE(LUPRI,*) 'NRECORDS:',NRECORDS
        CALL FLSHFO(LUPRI)
      END IF
C
      IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
        SIGN = 1.0d0
      ELSE IF (ITYPE.EQ.5) THEN
        SIGN = -1.0d0
      END IF

      IF (.NOT.LONDON2) THEN
C
C       initialize XINT and set up IOFFRD:
C       (set to -1 for those gammas we don't want to read)
C
        KOFF = 1
        DO ID = 1, NUMD
          DO IGAMMA = 1, NBAST
            IOFFRD(IGAMMA,ID) = -1
          END DO
          DO IG = 1, NUMG
            ISYMD = ISAO(IDEL(ID))
            ISYMG = ISAO(IGAM(IG))
            ISYPQ = MULD2H(MULD2H(ISYMD,ISYMG),ICORSY)
            IOFFRD(IGAM(IG),ID) = KOFF - 1
            CALL DZERO(XINT(KOFF),NDIMPQ(ISYPQ))
            KOFF = KOFF + NDIMPQ(ISYPQ)
            IF(LOCDBG) WRITE(LUPRI,'(A,5I5)')'ID,IDEL,IG,IGAM,IOFFRD:',
     &            ID,IDEL(ID),IG,IGAM(IG),IOFFRD(IGAM(IG),ID)
          END DO
        END DO
      
C
C       loop over all records and read all those with relevant
C       delta indeces
C
        DO NREC = 1, NRECORDS

          ! find delta index ID, if not found skip record...
          ID = 0
          DO JD = 1, NUMD
            IF (IRECNR(NREC).EQ.IDEL(JD)) ID = JD
            IF (LOCDBG)
     &      WRITE(LUPRI,*)'IRECNR,JD,IDEL,ID:',
     &                     IRECNR(NREC),JD,IDEL(JD),ID
          END DO

          IF (ID.GT.0) THEN
           READ(UNIAO2,ERR=2000,REC=NREC) RBUF,IBUF,NINT
           IF (NIBUF.EQ.1) THEN
             DO I = 1,NINT
                IP = IBTAND(       IBUF(I,1)         ,IBIT1)
                IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                IF (IOFFRD(IR,ID).GE.0) THEN
                  IADR = IOFFRD(IR,ID) + IADRPQ(IP,IQ)

C                 WRITE(LUPRI,*) 'IR,ID,IP,IQ,IADR:',IR,ID,IP,IQ,IADR
C                 IF (IP.LT.1 .OR. IP.GT.NBAST .OR.
C    &                IQ.LT.1 .OR. IQ.GT.NBAST .OR.
C    &                IR.LT.1 .OR. IR.GT.NBAST      ) THEN
C                   CALL QUIT('INDICES OUT OF RANGE IN CCRDAOD1.')
C                 END IF
C                 CALL FLSHFO(LUPRI)

                  XINT(IADR) = SIGN*RBUF(I)
                END IF
             END DO
           ELSE
             DO I = 1,NINT
                IP = IBTAND(       IBUF(I,2)       ,IBIT1)
                IQ = IBTAND(IBTSHR(IBUF(I,2),NBITS),IBIT1)
                IR = IBTAND(       IBUF(I,1)       ,IBIT1)
                IF (IOFFRD(IR,ID).GE.0) THEN
                  IADR = IOFFRD(IR,ID) + IADRPQ(IP,IQ)
                  XINT(IADR) = SIGN*RBUF(I)
                END IF
             END DO
           END IF

c          IF (LOCDBG .AND. NIBUF.EQ.1) THEN
c            WRITE(LUPRI,*) 'INTEGRALS READ FROM RECORD:',NREC
c            DO I = 1,NINT
c               IP = IBTAND(       IBUF(I,1)         ,IBIT1)
c               IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
c               IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
c               IF (IOFFRD(IR,ID).GE.0) THEN
c                 IADR = IOFFRD(IR,ID) + IADRPQ(IP,IQ)
c                 WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8,5X,I5)')
c    &               ' ## ', IP, IQ , IR, IDEL(ID), XINT(IADR), IADR
c               ELSE
c                 WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8,5X,A)')
c    &             ' ## ', IP, IQ , IR, IDEL(ID),SIGN*RBUF(I),'skipped'
c               END IF
c            END DO
c          END IF

          ENDIF

        END DO

      ELSE

C
C       loop over all records and read all those with relevant
C       delta indeces, exchange P and Q
C
        DO NREC = 1, NRECORDS

          ! find delta index ID, if not found skip record...
          ID = 0
          DO JD = 1, NUMD
            IF (IRECNR(NREC).EQ.IDEL(JD)) ID = JD
          END DO

          IF (ID.GT.0) THEN
           READ(UNIAO2,ERR=2000,REC=NREC) RBUF,IBUF,NINT
           IF (NIBUF.EQ.1) THEN
             DO I = 1,NINT
                IP = IBTAND(       IBUF(I,1)         ,IBIT1)
                IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                IF (IOFFRD(IR,ID).GE.0) THEN
                  IADR = IOFFRD(IR,ID) + IADRPQ(IQ,IP)

C                 WRITE(LUPRI,*) 'IR,ID,IP,IQ,IADR:',IR,ID,IP,IQ,IADR
C                 IF (IP.LT.1 .OR. IP.GT.NBAST .OR.
C    &                IQ.LT.1 .OR. IQ.GT.NBAST .OR.
C    &                IR.LT.1 .OR. IR.GT.NBAST      ) THEN
C                   CALL QUIT('INDICES OUT OF RANGE IN CCRDAOD1.')
C                 END IF
C                 CALL FLSHFO(LUPRI)

                  XINT(IADR) = SIGN*RBUF(I)
                END IF
             END DO
           ELSE
             DO I = 1,NINT
                IP = IBTAND(       IBUF(I,2)       ,IBIT1)
                IQ = IBTAND(IBTSHR(IBUF(I,2),NBITS),IBIT1)
                IR = IBTAND(       IBUF(I,1)       ,IBIT1)
                IF (IOFFRD(IR,ID).GE.0) THEN
                  IADR = IOFFRD(IR,ID) + IADRPQ(IQ,IP)
                  XINT(IADR) = SIGN*RBUF(I)
                END IF
             END DO
           END IF

           IF (LOCDBG .AND. NIBUF.EQ.1) THEN
             WRITE(LUPRI,*) 'INTEGRALS READ FROM RECORD:',NREC
             DO I = 1,NINT
                IP = IBTAND(       IBUF(I,1)         ,IBIT1)
                IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                IF (IOFFRD(IR,ID).GE.0) THEN
                  IADR = IOFFRD(IR,ID) + IADRPQ(IQ,IP)
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8,5X,I5)')
     &               ' ## ', IP, IQ , IR, IDEL(ID), XINT(IADR), IADR
                ELSE
                  WRITE (LUPRI,'(10X,A,2X,4I4,5X,1P,D16.8,5X,A)')
     &             ' ## ', IP, IQ , IR, IDEL(ID),SIGN*RBUF(I),'skipped'
                END IF
             END DO
           END IF

          ENDIF

        END DO

      END IF

      IF (LOCDBG) THEN
       IF (LONDON2) WRITE(LUPRI,*) 'CCRDAOD1> LONDON2!'
       DO ID = 1, NUMD
        DO IG = 1, NUMG
          ISYMD = ISAO(IDEL(ID))
          ISYMG = ISAO(IGAM(IG))
          ISYPQ = MULD2H(MULD2H(ISYMD,ISYMG),ICORSY)
          KOFF = IOFFRD(IGAM(IG),ID) + 1
          XNORM=DDOT(NDIMPQ(ISYPQ),XINT(KOFF),1,
     *                              XINT(KOFF),1)
          WRITE (LUPRI,'(a,2i5,f10.5)') 'CCRDAOD1> ', 
     *                           IDEL(ID),IGAM(IG),XNORM
          CALL OUTPUT(XINT(KOFF),1,NDIMPQ(ISYPQ),1,1,
     *                             NDIMPQ(ISYPQ),1,1,LUPRI)
        END DO
       END DO
      END IF


      RETURN
*----------------------------------------------------------------------*
C     I/O error handling:
*----------------------------------------------------------------------*
2000  CONTINUE
      CALL QUIT('I/O error in CCRDAOD...')

      END 
*======================================================================*
*======================================================================*
C  /* Deck RDERILBS */
      SUBROUTINE RDERILBS(IRECNR,NRECS,MXBUF,NFILES)
*----------------------------------------------------------------------*
C  read orbital (delta) indeces for all records on all files
C  into array IRECNR
*----------------------------------------------------------------------*
      IMPLICIT NONE

#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
#include <inftap.h>
#include <eribuf.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXBUF, NFILES, IDELTA, ISCOOR, ILINE, NBUFTOT
      INTEGER IRECNR(MXBUF,0:NFILES-1), NRECS(0:NFILES-1)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'entered RDERILBS...'
        WRITE(LUPRI,*) 'NBUFX,NFILES:',NBUFX,NFILES
        CALL FLSHFO(2)
      END IF

      IF (LUINTR .LE. 0) THEN
          IF (LOCDBG) WRITE(LUPRI,*) 'Open AOTWODIS file...'
          CALL GPOPEN(LUINTR,'AOTWODIS','OLD',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
      END IF

      REWIND (LUINTR)

      NBUFTOT = 0
      DO ISCOOR = 0, NFILES-1
        NRECS(ISCOOR) = 0
        NBUFTOT = NBUFTOT + NBUFX(ISCOOR)
      END DO

      DO ILINE = 1, NBUFTOT
        READ(LUINTR) IDELTA, ISCOOR
        NRECS(ISCOOR) = NRECS(ISCOOR) + 1
        IRECNR(NRECS(ISCOOR),ISCOOR) = IDELTA
        IF (LOCDBG) WRITE(LUPRI,*) 'ILINE,IDELTA,ISCOOR,NRECS:',
     &                              ILINE,IDELTA,ISCOOR,NRECS(ISCOOR)
      END DO

      IF (LOCDBG) WRITE(LUPRI,*) 'leaving RDERILBS...'
 
      RETURN
      END 
*======================================================================*
