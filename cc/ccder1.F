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
c /* deck ccder1 */
*=====================================================================*
       SUBROUTINE CCDER1(IATOM,LABELOP,LDERINT,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: get first derivative integrals from ABACUS
*             and sort them into the ordering used in CC
*             
*    Written by Christof Haettig, 06-May-1998
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <ccorb.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCDER1> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.) 

      LOGICAL LDERINT(8,3)
      CHARACTER*8 LABELOP
      INTEGER LWORK
      INTEGER IATOM, NDMAT, MAXDIF, MXCOMP
      INTEGER KDMAT, KEND, KFMAT, KISYMDM, KIFCTYP, LEND

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0
#else
      DOUBLE PRECISION WORK(LWORK)
#endif


      CALL QENTER('CCDER1')
*---------------------------------------------------------------------*
* check, if not the same integrals have been calculated the last time:
*---------------------------------------------------------------------*
C     IF (LABELOP .EQ. LAST_LABELOP) THEN
C       CALL QEXIT('CCDER1')
C       RETURN
C     END IF

*---------------------------------------------------------------------*
* some intializations
*---------------------------------------------------------------------*
      IF (LOCDBG) WRITE (LUPRI,*) MSGDBG, 'entered CCDER1.'
      
      KEND = 1          ! work space

      NDMAT = 0         ! # densityies for fock matrices

      MAXDIF = 1        ! order for derivatives --> first derivatives

      MXCOMP = 3        ! max. comp. per atom
 
      CALL RHSINI       ! initialize some ABACUS common blocks
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      KDMAT   = KEND
      KFMAT   = KDMAT + NDMAT*N2BASX
      KISYMDM = KFMAT + NDMAT*MXCOMP*NUCDEG(IATOM)*N2BASX
      KIFCTYP = KISYMDM + NDMAT*MXCOMP*NUCDEG(IATOM)
      KEND    = KIFCTYP + NDMAT*MXCOMP*NUCDEG(IATOM)
      LEND    = LWORK - KEND + 1
        
      CALL CCGETH2D(IATOM,MAXDIF,LDERINT,LABELOP,
     &              WORK(KDMAT),WORK(KFMAT), NDMAT,
     &              WORK(KISYMDM),WORK(KIFCTYP), MXCOMP,
     &              WORK(KEND),LEND)

C     LAST_LABELOP = LABELOP

*---------------------------------------------------------------------*
* return:
*---------------------------------------------------------------------*
      CALL FLSHFO(LUPRI)

      CALL QEXIT('CCDER1')
      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CCDER1                               *
*=====================================================================*
*======================================================================*
C  /* Deck ccsortderao */
      SUBROUTINE CCSORTDERAO(WORK,LWORK,MXCOMP,LDERINT,ITYPE,IPRINT)
*----------------------------------------------------------------------*
C
C     Purpose: sort derivative 2-el. AO integrals.
C
C     MXCOMP: maximum number of 'components'
C             --> three (x,y,z) for first derivative integrals 
C
C     ITYPE :  1 - geometric first derivatives
C              5 - magnetic first derivatives
C
C     Written by Christof Haettig 06-May-1998
C     based on Henrik Koch's CCSD_SORTAO  routine
C     magnetic derivatives added Sep-1999
C
*----------------------------------------------------------------------*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include <priunit.h>
#include <iratdef.h>
#include <inftap.h>
C
      PARAMETER (LBUF = 600)
C
      LOGICAL LOCDBG  ! local debug flag
      PARAMETER (LOCDBG = .false.)
C
      DIMENSION WORK(LWORK)
      LOGICAL LDERINT(8,MXCOMP)
C
      INTEGER LENGTH(8), LUAODR(8)
C
      CHARACTER*8 NAME(8)
C
      DATA NAME  /'CCAODER1','CCAODER2','CCAODER3','CCAODER4',
     *            'CCAODER5','CCAODER6','CCAODER7','CCAODER8'/
#include <ccorb.h>
#include <ccsdsym.h>

*----------------------------------------------------------------------*
* open files for sorted integrals:
*----------------------------------------------------------------------*
      DO ISYM = 1,NSYM
         LUAODR(ISYM) = -1
         CALL WOPEN2(LUAODR(ISYM),NAME(ISYM),64,0)
      END DO

*----------------------------------------------------------------------*
*     set buffer information.
*----------------------------------------------------------------------*
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBAST .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         CALL QUIT('CCSORTDERAO not implemented for more '//
     &        'than 255 functions.')
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1

      
*----------------------------------------------------------------------*
*     Buffer allocation.
*----------------------------------------------------------------------*
      KIBUF = 1
      KRBUF = KIBUF + (NIBUF*LBUF + 1)/IRAT + 1
      KAOAB = KRBUF + LBUF
      KAOG  = KAOAB + (N2BASX     + 1)/IRAT + 1
      KEND1 = KAOG  + (NBAST*NSYM*NSYM*MXCOMP + 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0)
     &     CALL QUIT('Insufficient work space in CCSORTDERAO.')
*----------------------------------------------------------------------*
*     Calculate in the index arrays needed in the sort.
*----------------------------------------------------------------------*

      CALL CCSD_INIT2B(WORK(KAOAB),WORK(KAOG),WORK(KAOG),ITYPE,
     &                 .FALSE.,MXCOMP,LDERINT)

*----------------------------------------------------------------------*
*     precalculate length of integral (* *|* del) distributions: 
*----------------------------------------------------------------------*
      DO ISYMD = 1, NSYM
         LENGTH(ISYMD) = 0
         DO ICOOR = 1, MXCOMP
           DO ICORSY = 1, NSYM
             IF ( LDERINT(ICORSY,ICOOR) ) THEN
               ISYDIS = MULD2H(ISYMD,ICORSY)

               IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
                 LENGTH(ISYMD) = LENGTH(ISYMD) + NDISAO(ISYDIS)
               ELSE IF (ITYPE.EQ.5) THEN
                 LENGTH(ISYMD) = LENGTH(ISYMD) + NDISAOSQ(ISYDIS)
               ELSE
                 CALL QUIT('Illegal ITYPE in CCSORTDERAO.')
               END IF

             END IF
           END DO
         END DO

      END DO

*----------------------------------------------------------------------*
*     Loop over batches of integrals.
*----------------------------------------------------------------------*
C
      DO 100 ISYMD = 1,NSYM
C
         IOFF2 = 1
C
         NTOTD  = NBAS(ISYMD)
	 IF (NTOTD .EQ. 0) GOTO 100

         NUMBAT = MIN(NTOTD,LWRK1/LENGTH(ISYMD))
C
         IF (NUMBAT .EQ. 0) THEN
            WRITE(LUPRI,*) 'In CCSD_SORTAO NUMBAT is zero'
            CALL QUIT('Insufiicient work space in CCRDAO')
         ENDIF
C
         ITOTBA = (NTOTD-1)/NUMBAT + 1
C
         ID1   = IBAS(ISYMD) + 1
         ID2   = IBAS(ISYMD)
         IOFF1 = IBAS(ISYMD)
C
         DO 200 I = 1,ITOTBA
C
            INUMBA = NUMBAT
            IF (NUMBAT*I .GT. NTOTD) THEN
               INUMBA = NTOTD - NUMBAT*(I-1)
            ENDIF
C
            ID2 = ID2 + INUMBA
C
            CALL DZERO(WORK(KEND1),LENGTH(ISYMD)*INUMBA)
C
            CALL CCSORTDER1(WORK(KEND1),WORK(KIBUF),WORK(KRBUF),
     *                      WORK(KAOAB),WORK(KAOG),ISYMD,LENGTH(ISYMD),
     *                      IOFF1,ID1,ID2,NIBUF,LBUF,NBITS,IBIT1,ITYPE)
C
c           CALL AROUND('Norm of integral distributions:')
c           IPRC = KEND1
c           DO IPRD = ID1,ID2
c              WRITE (LUPRI,*) 'D distribution',IPRD
c              DO ICOOR = 1, MXCOMP
c                DO ICORSY = 1, NSYM
c                  IF ( LDERINT(ICORSY,ICOOR) ) THEN
c                    WRITE (LUPRI,*) 'coordinate/symmetry',ICOOR,ICORSY
c                    ISYDIS = MULD2H(ISYMD,ICORSY)
c                    XNORM = 0.0d0
c                    DO IPSYMG = 1,NSYM
c                      ISYMAB = MULD2H(IPSYMG,ISYDIS)
c                      IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
c                        LEN = NNBST(ISYMAB) *  NBAS(IPSYMG)
c                        XNORM=XNORM+DDOT(LEN,WORK(IPRC),1,WORK(IPRC),1)
c                        IPRC = IPRC + NNBST(ISYMAB)*NBAS(IPSYMG)
c                      ELSE IF (ITYPE.EQ.5) THEN
c                        LEN = N2BST(ISYMAB) *  NBAS(IPSYMG)
c                        XNORM=XNORM+DDOT(LEN,WORK(IPRC),1,WORK(IPRC),1)
c                        DO G = 1, NBAS(IPSYMG)
c                         IGAM = G + IBAS(IPSYMG)
c                         KOFFX = IPRC + N2BST(ISYMAB)*(G-1)
c                         YNORM = DDOT(N2BST(ISYMAB),WORK(KOFFX),1,
c    *                                               WORK(KOFFX),1)
c                         WRITE (LUPRI,'(a,5i5,f10.5)'),'CCSORTDERAO> ',IPRD
c    *                      IGAM,ICOOR,ICORSY,IOFF2+KOFFX-KEND1,YNORM
c                        END DO
c                        IPRC = IPRC + N2BST(ISYMAB)*NBAS(IPSYMG)
c                      ELSE
c                        CALL QUIT('Unknown ITYPE in CCSORTDERAO.')
c                      END IF
c                    END DO
c                    WRITE (LUPRI,*) 'Norm of distribution:',XNORM
c      
c                  END IF
c                END DO
c              END DO
c           END DO
C
            CALL CCSORTDER2(WORK(KEND1),IOFF2,INUMBA,
     *                      LENGTH(ISYMD),ISYMD,NAME,LUAODR)
C
            IF ( (IPRINT.GT.50) .OR. LOCDBG ) THEN
               CALL AROUND('Integral distribution')
               IPRC = KEND1
               DO IPRD = ID1,ID2
                  WRITE (LUPRI,*) 'D distribution',IPRD
                  DO ICOOR = 1, MXCOMP
                    DO ICORSY = 1, NSYM
                      IF ( LDERINT(ICORSY,ICOOR) ) THEN
                        WRITE (LUPRI,*) 'coordinate/symmetry',
     &                        ICOOR,ICORSY
                        ISYDIS = MULD2H(ISYMD,ICORSY)
                        DO IPSYMG = 1,NSYM
                          WRITE(LUPRI,*) 'Gamma symmetry',IPSYMG
                          ISYMAB = MULD2H(IPSYMG,ISYDIS)
                          
                          IF (LOCDBG) THEN
                          
                            DO G = 1, NBAS(IPSYMG)
                              IG = IBAS(IPSYMG) + G
                              DO ISYMB = 1, NSYM
                               ISYMA = MULD2H(ISYMB,ISYMAB)
                               WRITE (LUPRI,*) 'symmet.:',
     *                                ISYMA,ISYMB,IPSYMG
                               IF (((ITYPE.EQ.0 .OR. ITYPE.EQ.1)
     *                               .AND. (ISYMB .GT. ISYMA)   )
     *                             .OR. (ITYPE.EQ.5)             )THEN
                                  DO B = 1, NBAS(ISYMB)
                                    IB = IBAS(ISYMB) + B
                                    DO A = 1, NBAS(ISYMA)
                                      IA = IBAS(ISYMA) + A
                                      WRITE (LUPRI,'(4I5,G20.10,I10)')
     *                                 IA,IB,IG,IPRD,WORK(IPRC),
     *                                 IPRC-KEND1+1
                                      IPRC = IPRC + 1
                                    END DO
                                  END DO
                               ELSE IF ((ITYPE.EQ.0 .OR. ITYPE.EQ.1)
     *                                  .AND. (ISYMB.EQ.ISYMA)) THEN
                                  DO B = 1, NBAS(ISYMB)
                                    IB = IBAS(ISYMB) + B
                                    DO A = 1, B
                                      IA = IBAS(ISYMA) + A
                                      WRITE (LUPRI,'(4I5,G20.10,I10)')
     *                                 IA,IB,IG,IPRD,WORK(IPRC),
     *                                 IPRC-KEND1+1
                                      IPRC = IPRC + 1
                                    END DO
                                  END DO
                               ELSE IF ((ITYPE.EQ.0 .OR. ITYPE.EQ.1)
     *                                  .AND. (ISYMB.LT.ISYMA)) THEN
                                  CONTINUE
                               ELSE
                                  CALL QUIT(
     *                                 'Unknown ITYPE in CCSORTDERAO.')
                               END IF
                              END DO
                            END DO
               
                          ELSE
                            IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
                              CALL OUTPUT(WORK(IPRC),1,NNBST(ISYMAB),1,
     *                                    NBAS(IPSYMG),NNBST(ISYMAB),
     *                                    NBAS(IPSYMG),1,LUPRI)
                              IPRC = IPRC + NNBST(ISYMAB)*NBAS(IPSYMG)
                            ELSE IF (ITYPE.EQ.5) THEN
                              CALL OUTPUT(WORK(IPRC),1,N2BST(ISYMAB),1,
     *                                    NBAS(IPSYMG),N2BST(ISYMAB),
     *                                    NBAS(IPSYMG),1,LUPRI)
                              IPRC = IPRC + N2BST(ISYMAB)*NBAS(IPSYMG)
                            ELSE
                                  CALL QUIT(
     *                                 'Unknown ITYPE in CCSORTDERAO.')
                            END IF
                          END IF
                        END DO
                      END IF
                    END DO
                  END DO
               END DO
            END IF
C
            ID1   = ID1   + INUMBA
            IOFF1 = IOFF1 + INUMBA
C
  200    CONTINUE
C
  100 CONTINUE
C
      DO ISYM = 1,NSYM
         CALL WCLOSE2(LUAODR(ISYM),NAME(ISYM),'KEEP')
      END DO
C
      CALL GPCLOSE(LU2DER,'DELETE')
C
      RETURN
      END
*======================================================================*
C  /* Deck ccsortder1 */
      SUBROUTINE CCSORTDER1(XINT,IBUF,RBUF,KAOAB,KAOGDER,ISYMD,LENGTH,
     *                      IOFF,ID1,ID2,NIBUF,LBUF,NBITS,IBIT1,ITYPE)
*----------------------------------------------------------------------*
C
C     Written by Christof Haettig 06-May-1998
C     based on Henrik Kochs CCSD_SORT1 routine
C     generalizations for London integrals Sep-1999
C
*----------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <ibtpar.h>
#include <ccorb.h>
      DIMENSION XINT(LENGTH,*),IBUF(LBUF,NIBUF),RBUF(LBUF)
      DIMENSION KAOAB(NBAST,NBAST),KAOGDER(NBAST,NSYM,NSYM,*)
      INTEGER ITYPE
C
#include <inftap.h>
#include <ibtfun.h>
C
      CALL REWSPL(LU2DER)
C
      IF (ITYPE.EQ.5) THEN
        CALL MOLLAB('AO2MGINT',LU2DER,LUPRI)
      END IF
C
      IF (ITYPE.EQ.0 .OR. ITYPE.EQ.1) THEN
        SIGN = +1.0D0
      ELSE IF (ITYPE.EQ.5) THEN
        SIGN = -1.0D0
      ELSE
        CALL QUIT('Unknown ITYPE in CCSORTDER1.')
      END IF
C
      IF (NIBUF .EQ. 1) THEN
C
   10    READ(LU2DER,ERR=2000) RBUF,IBUF,NINT
C
         IF (NINT .GT. 0) THEN

           DO I = 1, NINT

              LABEL = IBUF(I,1)
              IP = IBTAND(       LABEL         ,IBIT1)
 
              IF (IP .EQ. 0) THEN

*               ... FLAG : ICOOR,ICORSY...
                ICOOR  = IBTAND(IBTSHR(LABEL,NBITS),IBIT1)
                ICORSY = IBTAND(IBTSHR(LABEL,2*NBITS),IBIT1) + 1

              ELSE

*               ... INTEGRAL : (IP IQ | IR IS) ...
                IQ = IBTAND(IBTSHR(LABEL,  NBITS),IBIT1)
                IR = IBTAND(IBTSHR(LABEL,2*NBITS),IBIT1)
                IS = IBTAND(IBTSHR(LABEL,3*NBITS),IBIT1)
C
                IF ((IS .GE. ID1) .AND. (IS .LE. ID2)) THEN
                   ! store as (PQ|RS)
                   IADR = KAOGDER(IR,ISYMD,ICORSY,ICOOR) + KAOAB(IP,IQ)
                   XINT(IADR,IS-IOFF) = SIGN * RBUF(I)
                   !IF(ITYPE.EQ.5)WRITE (LUPRI,*)IP,IQ,IR,IS,XINT(IADR,IS-IOFF)
                ENDIF
                IF ((IR .GE. ID1) .AND. (IR .LE. ID2)) THEN
                   ! store as (QP|SR)
                   IADR = KAOGDER(IS,ISYMD,ICORSY,ICOOR) + KAOAB(IQ,IP)
                   XINT(IADR,IR-IOFF) = RBUF(I)
                   !IF(ITYPE.EQ.5)WRITE (LUPRI,*)IQ,IP,IS,IR,XINT(IADR,IR-IOFF)
                ENDIF
                IF ((IP .GE. ID1) .AND. (IP .LE. ID2)) THEN
                   ! store as (SR|QP)
                   IADR = KAOGDER(IQ,ISYMD,ICORSY,ICOOR) + KAOAB(IS,IR)
                   XINT(IADR,IP-IOFF) = RBUF(I)
                   !IF(ITYPE.EQ.5)WRITE (LUPRI,*)IS,IR,IQ,IP,XINT(IADR,IP-IOFF)
                ENDIF
                IF ((IQ .GE. ID1) .AND. (IQ .LE. ID2)) THEN
                   ! store as (RS|PQ)
                   IADR = KAOGDER(IP,ISYMD,ICORSY,ICOOR) + KAOAB(IR,IS)
                   XINT(IADR,IQ-IOFF) = SIGN * RBUF(I)
                   !IF(ITYPE.EQ.5)WRITE (LUPRI,*)IR,IS,IP,IQ,XINT(IADR,IQ-IOFF)
                ENDIF

              END IF
C
           END DO

         ELSE IF (NINT .LT. 0) THEN
           GOTO 100
         END IF
C
         GOTO 10
C
  100    CONTINUE
C
      ELSE
         CALL QUIT('Error: option NIBUF > 2 not '//
     &        'implemented in CCSORTDER1.')
      ENDIF
C
C
C
      RETURN
 
 2000 CALL QUIT('Error reading derivative integrals in CCSORTDER1')

      END
*======================================================================*
C  /* Deck ccsortder2 */
      SUBROUTINE CCSORTDER2(XINT,IOFF,INUMBA,LENGTH,ISYM,NAME,LUAODR)
*----------------------------------------------------------------------*
C
C     Written by Christof Haettig 06-May-1998
C     based on Henrik Kochs CCSD_SORT2 routine
C
*----------------------------------------------------------------------*
#include <implicit.h>
      DIMENSION XINT(*)
C
      INTEGER LUAODR(8)
      CHARACTER*8 NAME(8)
C
      CALL PUTWA2(LUAODR(ISYM),NAME(ISYM),XINT,IOFF,LENGTH*INUMBA)
      IOFF = IOFF + LENGTH*INUMBA
C
      RETURN
      END
*======================================================================*
C  /* Deck ccsd_init2b */
      SUBROUTINE CCSD_INIT2B(KAOAB,KAOG,KAOGDER,ITYPE,LDISTRIB,
     &                       MXCOMP,LDERINT)
*----------------------------------------------------------------------*
C
C     Henrik Koch and Alfredo Sanchez.       29-Jun-1994
C     derivative option by Christof Haettig  06-May-1998
C
C     Set up indexing arrays for integral sort
C
C     ITYPE = 0 : set up KAOG for undiff. 2-el integrals
C                 (KAOGDER ignored)
C
C     ITYPE = 1 : set up KAOGDER for geom. 1. derivative 2-el integrals
C                 (KAOG ignored)
C
C     ITYPE = 5 : set up KAOGDER for magn. 1. derivative 2-el integrals
C                 (KAOG ignored)
C
C     LDISTRIB : if true only KAOAB is set up 
C                (used to sort distributions)
C
C     LDERINT : array of logicals that tell if there are any
C               derivative integrals for a given component
C               and given symmetry 
C
*----------------------------------------------------------------------*
#include <implicit.h>
#include <ccorb.h>
      DIMENSION KAOAB(NBAST,NBAST)
      DIMENSION KAOG(NBAST,NSYM)
      DIMENSION KAOGDER(NBAST,NSYM,NSYM,MXCOMP)
      LOGICAL LDISTRIB, LDERINT(8,MXCOMP)
      DIMENSION NDIMAB(8)
#include <ccsdsym.h>
 
*----------------------------------------------------------------------*
* set up KAOAB index array for the leading indeces alpha,beta :
*----------------------------------------------------------------------*

      IF ( ITYPE.EQ.0   .OR.  ITYPE.EQ.1 ) THEN
C        --------------------------------------------------------
C        for undiff. Hamiltonian or geometric derivatives pack
C        fast running indeces alpha,beta in triangular form: 
C        --------------------------------------------------------
         DO ISYMAB = 1,NSYM
            NDIMAB(ISYMAB) = NNBST(ISYMAB)
            ICOUNT = 0
            DO ISYMB = 1,NSYM
               ISYMA = MULD2H(ISYMB,ISYMAB)
               IF (ISYMB .GT. ISYMA) THEN
                  DO B = 1,NBAS(ISYMB)
                     IB = IBAS(ISYMB) + B
                     DO A = 1,NBAS(ISYMA)
                        IA = IBAS(ISYMA) + A
                        ICOUNT = ICOUNT + 1
                        KAOAB(IA,IB) = ICOUNT
                        KAOAB(IB,IA) = ICOUNT
                     END DO
                  END DO
               ELSE IF (ISYMA .EQ. ISYMB) THEN
                  DO B = 1,NBAS(ISYMB)
                     IB = IBAS(ISYMB) + B
                     DO A = 1,B
                        IA = IBAS(ISYMA) + A
                        ICOUNT = ICOUNT + 1
                        KAOAB(IA,IB) = ICOUNT
                        KAOAB(IB,IA) = ICOUNT
                     END DO
                  END DO
               END IF
            END DO
         END DO
      ELSE
C        --------------------------------------------------------
C        if the Hamiltonian is not real (London orbitials) pack
C        fast running indeces alpha,beta in squared form: 
C        --------------------------------------------------------
         DO ISYMAB = 1,NSYM
            NDIMAB(ISYMAB) = N2BST(ISYMAB)
            DO ISYMB = 1,NSYM
               ISYMA = MULD2H(ISYMB,ISYMAB)
               ICOUNT = IAODIS(ISYMA,ISYMB)
               DO B = 1,NBAS(ISYMB)
                  IB = IBAS(ISYMB) + B
                  DO A = 1,NBAS(ISYMA)
                     IA = IBAS(ISYMA) + A
                     ICOUNT = ICOUNT + 1
                     KAOAB(IA,IB) = ICOUNT
                  END DO
               END DO
            END DO
         END DO
      END IF
 
 
*----------------------------------------------------------------------*
* set up KAOG/KAOGDER index array for the third index gamma:
*----------------------------------------------------------------------*
      IF (.NOT.LDISTRIB) THEN

       IF (ITYPE.EQ.0) THEN
C        --------------------------------------------------------
C        for undifferentiated integrals set up KAOG array:
C        --------------------------------------------------------
         DO ISYMD = 1,NSYM
            ISYDIS = MULD2H(ISYMD,ISYMOP)
            ICOUNT = 0
            DO ISYMG = 1,NSYM
               ISYMAB = MULD2H(ISYMG,ISYDIS)
               DO G = 1,NBAS(ISYMG)
                  IG = IBAS(ISYMG) + G
                  KAOG(IG,ISYMD) = ICOUNT
                  ICOUNT = ICOUNT + NDIMAB(ISYMAB)
               END DO
            END DO
         END DO
C
       ELSE IF (ITYPE.EQ.1 .OR. ITYPE.EQ.5) THEN
C        --------------------------------------------------------
C        for differentiated integrals set up KAOGDER array:
C        --------------------------------------------------------
         DO ISYMD = 1,NSYM
            ICOUNT = 0
            DO ICOOR = 1, MXCOMP
               DO ICORSY = 1, NSYM
                  IF ( LDERINT(ICORSY,ICOOR) ) THEN
                     ISYDIS = MULD2H(ISYMD,ICORSY)
                     DO ISYMG = 1,NSYM
                        ISYMAB = MULD2H(ISYMG,ISYDIS)
                        DO G = 1,NBAS(ISYMG)
                           IG = IBAS(ISYMG) + G
                           KAOGDER(IG,ISYMD,ICORSY,ICOOR) = ICOUNT
                           ICOUNT = ICOUNT + NDIMAB(ISYMAB)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
C
       ELSE
         CALL QUIT('Unknown ITYPE in CCSD_INIT2B.')
       END IF
C
      END IF

      RETURN
      END
*======================================================================*
c /* deck ccgeth2d */
*=====================================================================*
      SUBROUTINE CCGETH2D(IATOM, MAXDIF, LDERINT, LABELOP,
     &                    DMAT, FMAT, NDMAT, ISYMDM, IFCTYP, MXCOMP,
     &                    WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: call ABACUS to calculate derivatives of 2-el integrals 
*             and sort them to CC storage scheme
*
*    MAXDIF  : max. derivative order
*    LABELOP : operator label
*             
*    Christof Haettig, May 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <iratdef.h>
#include <aovec.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <nuclei.h>
#include <nodint.h>
#include <ccorb.h>
#include <cch2d.h>
#include <inftap.h>
#include <dummy.h>
#include <second.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCGETH2D> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.) 
      INTEGER NBUFS
      PARAMETER ( NBUFS = 600 )

      LOGICAL LDERINT(8,3)
      CHARACTER*8 LABELOP
      INTEGER IATOM, MAXDIF, NUMDIS, NDMAT
      INTEGER LWORK
      INTEGER ISYMDM(*), IFCTYP(*)


      LOGICAL RELCAL
      INTEGER MAXDIS
      INTEGER IBUF(NBUFS)
      INTEGER KJSTRS, KNPRIM, KNCONT, KIORBS, KJORBS, KKORBS, KLAST  
      INTEGER NFMAT, LFMAT, K2INT, L2INT, ITYPE, I2TYP
      INTEGER ICOOR, ICORSY, JR, JS, JP, JQ, JRS, JPQ, I
      INTEGER MXCOMP, ISYM, IPRHER

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL DMAT(*), FMAT(*)
      REAL TIM0, TIM1, TIM2, TIM3
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DMAT(*), FMAT(*)
      DOUBLE PRECISION TIM0, TIM1, TIM2, TIM3
#endif

* external functions
#if defined (SYS_CRAY)
      REAL BUF(NBUFS)
#else
      DOUBLE PRECISION BUF(NBUFS)
#endif

      CALL QENTER('CCGETH2D')
*---------------------------------------------------------------------*
* begin
*---------------------------------------------------------------------*
      IF (DEBUG.OR.LOCDBG) WRITE (LUPRI,*) MSGDBG, 'entered CCGETH2D.'
      
* initialize timing:
      TIM0 = SECOND()

* initialize address for next free element on work:
      KLAST = 1  

* set print level for integral part:
C     IPRHER = IPRINT / 5
      IPRHER = 0

*---------------------------------------------------------------------*
* allocate & initialize some vectors for TWOINT   
*---------------------------------------------------------------------*
      KJSTRS = KLAST
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
      IF (KLAST .GT. LWORK) THEN
        CALL QUIT('Insufficient work space in CCGETH2D.')
      END IF

      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &            WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),
     &            IATOM,.FALSE.,IPRHER)

      KLAST = KJORBS


*---------------------------------------------------------------------*
* initialize Fock matrices and related integer arrays ISYMDM,IFCTYP:
*---------------------------------------------------------------------*
      NFMAT = MXCOMP*NDMAT*NUCDEG(IATOM)
      LFMAT = NFMAT * N2BASX

      IF (NFMAT.GT.0) THEN
        IF (LFMAT.GT.0) CALL DZERO(FMAT,LFMAT)
        DO I = 1, NFMAT 
           ISYMDM(I) =  0
           IFCTYP(I) = 13
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate the derivatives of the twoelectron integrals and write
* them to file:
*---------------------------------------------------------------------*
      K2INT = KLAST                 
      L2INT = LWORK - K2INT + 1

      IF      ( LABELOP(1:5).EQ.'1DHAM' ) THEN
         ITYPE  = 1
         ! open file for derivative integrals here...
         LU2DER = -1
         CALL GPOPEN(LU2DER,'AO2DRINT','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LU2DER)
      ELSE IF ( LABELOP(1:5).EQ.'dh/dB' ) THEN
         ITYPE  = 5
         ! file for derivative integrals will be opened in DR2WRT
C        LU2DER = -1
C        CALL GPOPEN(LU2DER,'AO2MGINT','UNKNOWN','SEQUENTIAL',
C    &               'UNFORMATTED',IDUMMY,.FALSE.)
C        REWIND (LU2DER)
      ELSE
         WRITE (LUPRI,*) 'Unknown operator label in CCGETH2D.'
         CALL QUIT('Unknown operator label in CCGETH2D.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CCGETH2D> LABELOP,ITYPE,IATOM:',
     &                             LABELOP(1:5),ITYPE,IATOM
      END IF

      MAXDIS = 1
      I2TYP  = 0
      RELCAL = .FALSE.
      TKTIME = .FALSE.

      IF (LOCDBG) WRITE(LUPRI,*) 'CALL NOW TWOINT...'

      CALL TWOINT(WORK(K2INT),L2INT,DUMMY,
     &            FMAT, DMAT, NDMAT, ISYMDM, IFCTYP,
     &            DUMMY,IDUMMY,NUMDIS,MAXDIS,
     &            ITYPE,MAXDIF,IATOM,NODV,NOPV,NOCONT,
     &            TKTIME,IPRHER,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &            RETUR,IDUMMY,I2TYP,WORK(KJSTRS),
     &            WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &            DUMMY,RELCAL)

      TIM1 = SECOND()

      IF (LOCDBG) WRITE(LUPRI,*) 'RETURNED FROM TWOINT...'


*---------------------------------------------------------------------*
* sort the derivative integrals:
*---------------------------------------------------------------------*
      TIM2 = SECOND()

      IF (MAXDIF .EQ. 1) THEN
 
        IF (LOCDBG) WRITE (LUPRI,*) 'Setting up LDERINT array:'
        DO ICOOR = 1, 3
          DO ICORSY = 1, NSYM
            IF ( NDSINT(ICOOR,ICORSY-1) .GT. 0 ) THEN
              LDERINT(ICORSY,ICOOR) = .TRUE.
            ELSE
              LDERINT(ICORSY,ICOOR) = .FALSE.
            END IF
            IF (LOCDBG) WRITE (LUPRI,'(2X,3I5,L5)') ICOOR,ICORSY,
     &         NDSINT(ICOOR,ICORSY-1),LDERINT(ICORSY,ICOOR)
          END DO
        END DO

      ELSE
        CALL QUIT('MAXDIF <> 1 not implemented in CCGETH2D.')
      END IF

      IF (LOCDBG) CALL FLSHFO(LUPRI)

      CALL CCSORTDERAO(WORK,LWORK,3,LDERINT,ITYPE,IPRHER)

      TIM3 = SECOND()

*---------------------------------------------------------------------*
* print timing & return:
*---------------------------------------------------------------------*
      IF (IPRINT.GT.1 .OR. LOCDBG) THEN
         WRITE (LUPRI,'(/A,A,F12.2," seconds.")') 
     &        '>>> Time used in TWO',
     &   'INT:         ', TIM1 - TIM0
         WRITE (LUPRI,'(A,A,F12.2," seconds.")')
     &        '>>> Time used in CCS',
     &   'SORTDERAO:   ', TIM3 - TIM2
         WRITE (LUPRI,'(/A,A,F12.2," seconds.")')
     &        '>>> Total time used ',
     &   'in CCGETH2D :', SECOND() - TIM0
      END IF

      CALL FLSHFO(LUPRI)

      CALL QEXIT('CCGETH2D')
      RETURN
      END
*=====================================================================*
*              END OF SUBROUTINE CCGETH2D                             *
*=====================================================================*
