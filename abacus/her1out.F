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
===========================================================================
/* Comdeck log */
930209-hjaaj: add GSOLTT in SHDCHK
===========================================================================
#endif
C  /* Deck wrtund */
      SUBROUTINE WRTUND(STHMAT,NNBAST,NNBASX,IPRINT)
C
C     T. Helgaker Sep. 1986
C
C     This subroutine writes undifferentiated overlap, one-electron
C     Hamiltonian, and kinetic energy integrals on sequential file.
C
#include <implicit.h>
#include <maxaqn.h>
#include <priunit.h>
#include <drw2el.h>
#include <codata.h>
#include <inftap.h>
C
#include <ccom.h>
C
      PARAMETER (DP5 = 0.5D0)
      DIMENSION BUF(600), IBUF(600), STHMAT(NNBASX,*)
      CHARACTER*8 LAB123(3), TITLES, TITLEH, TITLET
      DATA LAB123/'********','        ','        '/
      DATA TITLES/'OVERLAP '/, TITLEH/'ONEHAMIL'/, TITLET/'KINETINT'/
      INDEX1(IJ) = INT(SQRT((2.0D0*IJ) + 0.25D0) + 0.4999D0)
      INDEX2(IJ) = IJ - INDEX1(IJ)*(INDEX1(IJ)-1)/2
C
      IF (IPRINT .GE. 4) CALL TITLER('Output from WRTUND','*',103)
C
      DO 100 ITYPE = 1, 3
         IF (ITYPE .EQ. 1) THEN
            IF (FINDPT) WRITE (LUPRI,'(/A)') 
     &      ' Metric is changed by direct perturbation theory (DPT) '
            IF (IPRINT .GE. 4) WRITE (LUPRI,'(/A)') ' Overlap matrix'
            WRITE (LUONEL) LAB123, TITLES
         ELSE IF (ITYPE .EQ. 2) THEN
            IF (FINDPT) WRITE (LUPRI,'(/A)') 
     &      ' One-electron Hamiltionian is changed by'//
     &      ' direct perturbation theory (DPT) '
            IF (IPRINT .GE. 4)
     &         WRITE (LUPRI,'(/A)')' One-el. Ham. matrix'
            WRITE (LUONEL) LAB123, TITLEH
         ELSE IF (ITYPE .EQ. 3) THEN
            IF (IPRINT .GE. 4) WRITE (LUPRI,'(/A)')' Kin. energy matrix'
            WRITE (LUONEL) LAB123, TITLET
         END IF
         IAB = 0
         ICOUNT = 0
         NBUF = 0
         DO IAB = 1, NNBAST
            AINT = STHMAT(IAB,ITYPE)
            IF (ITYPE .EQ. 1 .AND. FINDPT)
     *      AINT = AINT + DPTFAC * DP5 * ALPHAC**2 * STHMAT(IAB,3)
            IF (ABS(AINT) .GT. THRS) THEN
               ICOUNT = ICOUNT + 1
               BUF(ICOUNT) = AINT
               IBUF(ICOUNT) = IAB
               IF (ICOUNT .EQ. 600) THEN
                  WRITE (LUONEL) BUF, IBUF, ICOUNT
                  NBUF = NBUF + 1
                  IF (IPRINT .GE. 4) THEN
                     WRITE (LUPRI, '(//A,2I5,A,2I8)')
     *                  ' ICOUNT, NBUF  ', ICOUNT, NBUF,
     *                  ' NNBAST, NNBASX', NNBAST, NNBASX
                  END IF
                  IF (IPRINT .GE. 5) THEN
                     DO I = 1, ICOUNT
                        LABEL = IBUF(I)
                        IORBA = INDEX1(LABEL)
                        IORBB = INDEX2(LABEL)
                        WRITE (LUPRI,'(1X,3I5,5X,D24.12)')
     *                        LABEL, IORBA, IORBB, BUF(I)
                     END DO
                  END IF
                  ICOUNT = 0
               END IF
            END IF
         END DO
         NCOUNT = 600*NBUF + ICOUNT
         IF (ICOUNT .GT. 0) THEN
            WRITE (LUONEL) BUF, IBUF, ICOUNT
            NBUF = NBUF + 1
         END IF
         IF (IPRINT .GE. 4) THEN
            WRITE (LUPRI, '(//A,2I5,A,2I8)')
     *         ' ICOUNT, NBUF  ', ICOUNT, NBUF,
     *         ' NNBAST, NNBASX', NNBAST, NNBASX
         END IF
         IF (IPRINT .GE. 5) THEN
             DO 500 I = 1, ICOUNT
                LABEL = IBUF(I)
                IORBA = INDEX1(LABEL)
                IORBB = INDEX2(LABEL)
                WRITE (LUPRI,'(1X,3I5,5X,D24.12)')
     *                LABEL, IORBA, IORBB, BUF(I)
  500        CONTINUE
         END IF
         WRITE (LUONEL) BUF, IBUF, -1
         IF (IPRINT .GE. 3) THEN
            PERCNT = 100.0D0*NCOUNT
            PERCNT = PERCNT / IAB
            IF (ITYPE .EQ. 1) THEN
               WRITE (LUPRI,'(/I6,A,I4,A)') NCOUNT,
     *         ' atomic overlap integrals written in',NBUF,' buffers.'
            ELSE IF (ITYPE .EQ. 2) THEN
               WRITE (LUPRI,'(/I6,A,I4,A)') NCOUNT,
     *         ' one-el. Hamil. integrals written in',NBUF,' buffers.'
            ELSE
               WRITE (LUPRI,'(/I6,A,I4,A)') NCOUNT,
     *         ' kinetic energy integrals written in',NBUF,' buffers.'
            END IF
            WRITE (LUPRI,'(A,F7.2)') ' Percentage non-zero integrals:',
     *                              PERCNT
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck sorone */
      SUBROUTINE SORONE(SOINT,LSOINT,NMATS,INDMAX,IPRINT)
C
C     tuh Aug 1988
C
C     This subroutine sorts the differentiated overlap and one-electron
C     Hamiltonian elements from LUITMP on direct access unit LUDA1.
C     The elements read in on LUITMP are not the final SO integrals, but
C     AO integrals multiplied by symmetry factors and addresses of the
C     appropriate SO integrals.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
      PARAMETER (IBIT10 = 1023, IBUFMX = 600, NBUF = 1202)
      LOGICAL   OLDDX
      DIMENSION SOINT(LSOINT), BUF(600), IBUF(600)
C
#include <inftap.h>
#include <symmet.h>
#include <ccom.h>
#include <oneadr.h>
      DIMENSION BUFDA(NBUF), IBUFDA(IRAT*NBUF)
      EQUIVALENCE (BUFDA(1),IBUFDA(1))
#include <ibtfun.h>
C
C     **************************
C     ***** INITIALIZATION *****
C     **************************
C
C     Determine direct access buffer lengths
C
      LABUFD = NBUF
      LABUFI = IRAT*LABUFD
      MXABUF = (LABUFI - 2)/(IRAT + 1)
      IF (MXABUF .GT. INDMAX) THEN
         MXABUF = IRAT*((INDMAX + IRAT - 1)/IRAT)
         LABUFI = (IRAT + 1)*MXABUF + 2
         LABUFD = LABUFI/IRAT
      END IF
C
C     Make sure we have not been out of bounds for the labels. K.Ruud-06
C
      IF (INDMAX .GE. 2**22) THEN
         WRITE (LUPRI,*) ' Number of differentiated one-electron '//
     &        'integrals larger than available labels'
         WRITE (LUPRI,*) 'Please contact dalton-admin@kjemi.uio.no '//
     &        'for assistance'
         CALL QUIT('Differentiated integrals out of variable range')
      END IF
      MXABF2 = IRAT*MXABUF
      LABUF1 = LABUFI - 1
C
C     Open direct access file
C
      IF (LUDA1 .LE. 0) CALL GPOPEN(LUDA1,'ABACUS.DA1','UNKNOWN',
     &                              'DIRECT',' ',LABUFI,OLDDX)
      NMAT = LSOINT/INDMAX
      IF (NMAT .EQ. 0) THEN
         WRITE (LUPRI,'(//A/)') ' >>>> ERROR in SORONE: NMAT = 0.'
         CALL QUIT('NMAT zero in SORONE.')
      END IF
      IF (NMAT .GT. 2*NMATS) NMAT = 2*NMATS
      NPASS = (2*NMATS + NMAT - 1)/NMAT
      IF (IPRINT .GE. 10) THEN
         CALL HEADER(' <<<<< Subroutine SORONE >>>>>',-1)
         WRITE (LUPRI, '(2X,A,I10)')   'LSOINT:     ', LSOINT
         WRITE (LUPRI, '(2X,A,I5)')    'NPASS:      ', NPASS
         WRITE (LUPRI, '(2X,A,2I5)')   'LABUFD/I:   ', LABUFD, LABUFI
         WRITE (LUPRI, '(2X,A,I5)')    'MXABUF:     ', MXABUF
         WRITE (LUPRI, '(2X,A,D12.4)') 'THRS:       ', THRS
         WRITE (LUPRI, '(2X,A,I5)')    'INDMAX:     ', INDMAX
         WRITE (LUPRI, '(2X,A,I5)')    'NMATS:      ', NMATS
      END IF
C
C     ************************************
C     * CONSTRUCTED LOOP OVER LUITMP *****
C     ************************************
C
      IDISK = 1
      IFIRST = 1 - NMAT
      JDIR = 0
      ITYP = 1
      DO 200 IPASS = 1, NPASS
         IFIRST = IFIRST + NMAT
         ILAST  = IFIRST + NMAT - 1
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI, '(//2X,A,I5)') ' IPASS ', IPASS
            WRITE (LUPRI, '(2X,A,2I5)') ' IFIRST, ILAST ', IFIRST, ILAST
         END IF
C
C        Clear memory
C
         CALL DZERO(SOINT,NMAT*INDMAX)
C
C        Read through file and pick out elements in range.
C        Add elements. (Note that LUITMP does not contain final
C        SO integrals, only the AO integrals multiplied by symmetry
C        factors).
C
         REWIND LUITMP
  300    CONTINUE
            READ (LUITMP, END = 500) BUF, IBUF, LENGTH
            IF (IPRINT .GE. 10) WRITE (LUPRI,'(/,2X,A,I5)')
     *         'Buffer read from LUITMP, LENGTH =', LENGTH
            IF (LENGTH .EQ. 0) GO TO 300
C
C           Loop over elements in this buffer
C
            DO 400 I = 1, LENGTH
               LABEL = IBUF(I)
               IMAT  = IBTAND(LABEL,IBIT10)
               IF (IMAT.GE.IFIRST .AND. IMAT.LE.ILAST) THEN
                  IOFF = (IMAT - IFIRST)*INDMAX + IBTSHR(LABEL,10)
                  SOINT(IOFF) = SOINT(IOFF) + BUF(I)
C
C                 Print Section
C
                  IF (IPRINT .GE. 15) THEN
                     IF (MOD(I,20) .EQ. 1) WRITE (LUPRI,'(/5X,A/)')
     *                  'ITYP  ICAR IREPA INDEX       BUF         SOINT'
                     WRITE (LUPRI,'(1X,4I6,5X,1P,D12.4,5X,D12.4)')
     *                  (IMAT + NMATS - 1)/NMATS,
     *                  (MOD(IMAT-1,NMATS)+1+MAXREP)/(MAXREP+1),
     *                  MOD(IMAT - 1,MAXREP + 1),
     *                  IBTSHR(LABEL,10), BUF(I), SOINT(IOFF)
                  END IF
               END IF
  400       CONTINUE
C
C           Branch back and read next buffer
C
            GO TO 300
  500    CONTINUE
         IF (IPRINT .GE. 10) WRITE (LUPRI, '(/A/)')
     *      ' Last buffer in this pass processed. '
C
C        All elements in this pass have now been read.
C        Write matrices on direct access file.
C
         DO 600 IMAT = 1, NMAT
            JDIR = JDIR + 1
            IF (JDIR .EQ. NMATS + 1) THEN
               JDIR = 1
               ITYP = 2
            END IF
            ICAR  = (JDIR + MAXREP)/(MAXREP + 1)
            IREPA = MOD(JDIR - 1,MAXREP + 1)
C
            IF (IPRINT .GE. 10) WRITE (LUPRI, '(/A,4I5)')
     *         ' ITYP, JDIR, IREPA, ICAR ', ITYP, JDIR, IREPA, ICAR
C
C           Loop over elements for this perturbation (ICAR) and
C           symmetry of orbitals (IREPA).
C
            NOFF           = 0
            IBUFDA(LABUFI) = 0
            IOFF           = (IMAT - 1)*INDMAX
            LDISK          = 0
            DO 700 I = 1, INDMAX
               DINT = SOINT(IOFF + I)
               IF (ABS(DINT) .GE. THRS) THEN
                  NOFF                  = NOFF + 1
                  BUFDA(NOFF)           = DINT
                  IBUFDA(MXABF2 + NOFF) = I
               END IF
C
C              Write record when full or no more elements left
C
               IF ((NOFF.EQ.MXABUF).OR.(NOFF.GT.0.AND.I.EQ.INDMAX)) THEN
                  IBUFDA(LABUF1) = NOFF
                  CALL WRITDX (LUDA1,IDISK,LABUFI,IBUFDA(1))
                  LDISK = IDISK
                  IBUFDA(LABUFI) = IDISK
                  IF (IPRINT .GE. 10) THEN
                     WRITE (LUPRI, '(1X,3(A,I5))')
     *                'DA buffer for matrix type',ITYP,' and direction',
     *                JDIR, ' written:',IDISK
                     IF(IPRINT.GE.15)WRITE(LUPRI,'(20X,1P,D24.12,I10)')
     *                  (BUFDA(J), IBUFDA(MXABF2+J), J = 1, NOFF)
                  END IF
                  NOFF = 0
                  IDISK = IDISK + 1
               END IF
  700       CONTINUE
            LASTAD(ITYP,ICAR,IREPA) = LDISK
  600    CONTINUE
  200 CONTINUE
      IF (IPRINT .GE. 10) THEN
         NCRS = NMATS/(MAXREP + 1)
         CALL HEADER('LASTAD for overlap matrices',-1)
         DO 800 IREP = 0, MAXREP
            WRITE (LUPRI,'(/2X,10I5)') (LASTAD(1,I,IREP), I = 1, NCRS)
  800    CONTINUE
         CALL HEADER('LASTAD for one-el. Hamiltonian matrices',-1)
         DO 810 IREP = 0, MAXREP
            WRITE (LUPRI,'(/2X,10I5)') (LASTAD(2,I,IREP), I = 1, NCRS)
  810    CONTINUE
      END IF
      RETURN
      END
C  /* Deck getone */
      SUBROUTINE GETONE(ONEMAT,ICOOR,IREPA,NDIM,WORK,DNULL)
C
C     T. Helgaker Jan. 1985
C     Revised Jan. 28 1987 tuh
C     Revised Aug. 18 1988 for symmetry tuh
C
C     This subroutine retrieves differentiated overlap and one-electron
C     Hamiltonian matrices from direct access file (written in SORONE).
C     The input parameters are as follows:
C
C     ONEMAT - 'OMAT' for overlap matrix
C              'HMAT' for one-electron Hamiltonian matrix
C
C     ICOOR  - differentiation coordinate
C
C     IREPA  - symmetry of row orbitals (totally symmetric is 0)
C
C     NDIM   - dimension of matrix to be retrieved (triangular for
C              totally symmetric ICOOR, otherwise rectangular).
C
C     WORK   - work area of dimension NDIM
C
C     DNULL  - returns .TRUE. when the retrieved matrix is zero,
C              otherwise .FALSE.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <iratdef.h>
#include <nuclei.h>
      PARAMETER (IBIT22 = 4194303, NBUF = 1202)
      CHARACTER ONEMAT*(*)
      DIMENSION BUF(NBUF), IBUF(IRAT*NBUF), WORK(NDIM)
#include <oneadr.h>
#include <inftap.h>
      LOGICAL DNULL, NOTZER
      EQUIVALENCE (BUF(1),IBUF(1))
#include <ibtfun.h>
      IF (NDIM .EQ. 0) THEN
         WRITE(LUPRI,'(//A/)') ' >>>>> Error in GETONE: NDIM = 0 '
         CALL QUIT('NDIM zero in GETONE.')
      END IF
      IF (NBUF .LT. LABUFD) THEN
         WRITE(LUPRI,'(//A,I5,A,I5,/A)') ' NBUF = ', NBUF,
     *      ' is smaller than buffer length LABUFD =', LABUFD,
     *      ' Increase NBUF in GETONE.'
         CALL QUIT('Incorrect buffer length in GETONE.')
      END IF
      CALL DZERO(WORK,NDIM)
      IF (ONEMAT .EQ. 'OMAT') THEN
         IMAT = 1
      ELSE IF (ONEMAT .EQ. 'HMAT') THEN
         IMAT = 2
      ELSE
         CALL QUIT('Illegal matrix type in GETONE.')
      END IF
      NOTZER = .FALSE.
      IDISK = LASTAD(IMAT,ICOOR,IREPA)
      IF (IDISK .GT. 0) THEN
  200    CONTINUE
            CALL READDX(LUDA1,IDISK,LABUFI,IBUF(1))
            NOFF = IBUF(LABUF1)
            NOTZER = NOTZER .OR. (NOFF .GT. 0)
            DO 300 INT = 1, NOFF
               DINT = BUF(INT)
               LABEL = IBUF(MXABF2 + INT)
               IORBAB = IBTAND(LABEL,IBIT22)
               WORK(IORBAB) = DINT
  300       CONTINUE
            IDISK = IBUF(LABUFI)
         IF (IDISK .GT. 0) GO TO 200
      END IF
      DNULL = .NOT.NOTZER
      RETURN
      END
C  /* Deck shdpri */
      SUBROUTINE SHDPRI(WORK,LWORK)
C
C     tuh Aug 17 1988
C
C     This subroutine prints all differentiated overlap and one-electron
C     Hamiltonian matrices in the AO basis.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
#include <symmet.h>
      DIMENSION WORK(LWORK)
      LOGICAL DNULL
#include <ibtfun.h>
      CALL HEADER('Output from SHDPRI',-1)
      DO 100 IREPO = 0, MAXREP
         DO 200 ICENT = 1, NUCIND
            DO 300 ICOOR = 1, 3
               ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREPO,1)
               IF (ISCOOR .GT. 0) THEN
                  DO 400 IREPA = 0, MAXREP
                     IREPB = IBTXOR(IREPO,IREPA)
                     NDIMA = NAOS(IREPA + 1)
                     NDIMB = NAOS(IREPB + 1)
                     NDIM = NDIMA*NDIMB
                     IF (NDIM .GT. 0) THEN
                        IF (IREPO .EQ. 0) NDIM = NDIMA*(NDIMA + 1)/2
                        WRITE (LUPRI,'(//A,2X,A,I5)')
     *                         ' Coordinate and symmetry: ',
     *                          NAMEX(3*ICENT - 3 + ICOOR), IREPO
                        WRITE (LUPRI,'(/A,2I5)')
     *                     ' Symmetry of orbitals:', IREPA, IREPB
                        WRITE (LUPRI,'(A,2I5)')
     *                     ' Number of orbitals:  ', NDIMA, NDIMB
                        CALL GETONE('OMAT',ISCOOR,IREPA,NDIM,WORK,DNULL)
                        IF (DNULL) THEN
                           WRITE (LUPRI,'(/1X,A/)')
     *                       'Differentiated SO overlap matrix is zero.'
                        ELSE
                           CALL HEADER
     *                        ('Differentiated SO Overlap Matrix',-1)
                           IF (IREPO .EQ. 0) THEN
                              CALL OUTPAK(WORK,NDIMA,1,LUPRI)
                           ELSE
                              CALL OUTPUT(WORK,1,NDIMB,1,NDIMA,
     *                                    NDIMB,NDIMA,1,LUPRI)
                           END IF
                        END IF
                        CALL GETONE('HMAT',ISCOOR,IREPA,NDIM,WORK,DNULL)
                        IF (DNULL) THEN
                           WRITE (LUPRI,'(/1X,A/)')
     *                         'Differentiated SO one-electron '//
     *                         'Hamiltonian matrix is zero.'
                        ELSE
                           CALL HEADER('Differentiated SO '//
     *                        'One-Electron Hamiltonian Matrix',-1)
                           IF (IREPO .EQ. 0) THEN
                              CALL OUTPAK(WORK,NDIMA,1,LUPRI)
                           ELSE
                              CALL OUTPUT(WORK,1,NDIMB,1,NDIMA,
     *                                    NDIMB,NDIMA,1,LUPRI)
                           END IF
                        END IF
                     END IF
 400              CONTINUE
               END IF
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck shdchk */
      SUBROUTINE SHDCHK(WORK,LWORK,NODC,NODV,IPRINT)
C
C     tuh Aug 18 1988
C
C     This subroutine compares reorthonormalization and one-electron
C     Hamiltonian gradients calculated in AO and SO basis. Note that
C     only the totally symmetric derivatives can be checked this way.
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      LOGICAL NODC, NODV, DNULL
#include <abainf.h>
#include <energy.h>
#include <taysol.h>
#include <nuclei.h>
#include <symmet.h>
#include <inforb.h>
#include <ccom.h>
#include <rspprp.h>
#include <esg.h>
      PARAMETER (D0 = 0.0D0, D1000 = 1000.0D0)
      DIMENSION WORK(LWORK)
      IF (IPRINT .GE. 5) CALL HEADER('Output from SHDCHK',-1)
C
C     (1) Get the contravariant density and Fock matrices
C         Although these matrices have been calculated previously, have
C         have subsequently been transformed to the "symmetry distinct
C         AO basis" in DSYM1 and are therefore no longer available.
C
      KDER  = 1
      KDSO  = KDER  + NNBASX
      KFSO  = KDSO  + NNBAST
      KLAST = KFSO  + NNBAST
      LWRK  = LWORK - KLAST + 1
      IF (KLAST.GT.LWORK) CALL STOPIT('SHDCHK','DSOFSO',KLAST,LWORK)
      CALL DSOFSO(WORK(KDSO),WORK(KFSO),WORK(KLAST),LWRK,IPRINT,NODC,
     &            NODV)
C
C     (2) Loop over totally symmetric perturbations
C
      THRS10 = D1000*THRS
      DFMAX1 = D0
      DFMAX2 = D0
      DO 100 ICOOR = 1, 3*NUCDEP
         ISCOOR = IPTCNT(ICOOR,0,1)
         IF (ISCOOR .NE. 0) THEN
            IF (IPRINT .GE. 5) WRITE (LUPRI,'(/,A,I3/)')
     *         ' Symmetry coordinate ', ISCOOR
C
C           Loop over orbital symmetries
C
            TRACE1 = D0
            TRACE2 = D0
            DO 200 IREP = 0, MAXREP
               IREP1 = IREP + 1
               NDIM = NAOS(IREP1)*(NAOS(IREP1) + 1)/2
               IF (NDIM .GT. 0) THEN
                  CALL GETONE('OMAT',ISCOOR,IREP,NDIM,WORK(KDER),DNULL)
                  IF (.NOT.DNULL) TRACE1 = TRACE1
     *              - DDOT(NDIM,WORK(KDER),1,WORK(KFSO+NPARSU(IREP1)),1)
                  CALL GETONE('HMAT',ISCOOR,IREP,NDIM,WORK(KDER),DNULL)
                  IF (.NOT.DNULL) TRACE2 = TRACE2
     *              + DDOT(NDIM,WORK(KDER),1,WORK(KDSO+NPARSU(IREP1)),1)
               END IF
  200       CONTINUE
            DIFF1 = GRADFS(ISCOOR) - TRACE1
            DIFF2 = GRADKE(ISCOOR) + GRADNA(ISCOOR) + GSOLTT(ISCOOR)
     &            - TRACE2
            DFMAX1 = MAX(ABS(DIFF1),DFMAX1)
            DFMAX2 = MAX(ABS(DIFF2),DFMAX2)
            IF (ABS(DIFF1) .GT. THRS10 .AND. ( .NOT. ESG ) ) THEN
               WRITE(LUPRI,'(1P,//1X,A,I3,//3X,A/3X,3(D12.4,7X)/)')
     *           ' >>>>> WARNING - diff. between '//
     *           ' AO and SO reorth. gradient for coordinate',ISCOOR,
     *           ' AO gradient        SO gradient         Difference ',
     *           GRADFS(ISCOOR), TRACE1, DIFF1
               NWNABA = NWNABA + 1
            END IF
            IF (ABS(DIFF2) .GT. THRS10 .AND. ( .NOT. ESG ) ) THEN
               WRITE(LUPRI,'(1P,//1X,A,I3,//3X,A/3X,3(D12.4,7X)/)')
     *           ' >>>>> WARNING - diff. between '//
     *           ' AO and SO 1-el. int. gradient for coordinate',ISCOOR,
     *           ' AO gradient        SO gradient         Difference ',
     *           GRADKE(ISCOOR) + GRADNA(ISCOOR), TRACE2, DIFF2
               NWNABA = NWNABA + 1
            END IF
            IF (IPRINT .GE. 5) THEN
               WRITE (LUPRI,'(1P,1X,3(A,D12.4))')
     *           'Reorth. grad.:    AO = ',GRADFS(ISCOOR),
     *           ', SO = ',TRACE1, ', Diff. = ',DIFF1
               WRITE (LUPRI,'(1P,1X,3(A,D12.4))')
     *           '1-el. int. grad.: AO = ',
     *            GRADKE(ISCOOR) + GRADNA(ISCOOR),
     *           ', SO = ',TRACE2, ', Diff. = ',DIFF2
            END IF
         END IF
  100 CONTINUE
      IF (IPRINT .GT. 1) THEN
         WRITE (LUPRI,'(/A,1P,E10.2)')
     *      ' Greatest diff. between AO and SO'//
     *      ' reorthonormalization gradients:     ',DFMAX1
         WRITE (LUPRI,'(A,1P,E10.2)')
     *      ' Greatest diff. between AO and SO'//
     *      ' one-electron Hamiltonian gradients: ',DFMAX2
      END IF
      RETURN
      END
