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
C  /* Deck qpginp */
      SUBROUTINE QPGINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 7)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbiqpg.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NODC  ', '.NODV  ',
     *            'XXXXXXX', 'XXXXXXX', '.STOP  '/
C
      NEWDEF = (WORD .EQ. '*QPGCTL')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in QPGINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in QPGINP')
    1          CONTINUE
                  QSKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  QNODC = .TRUE.
               GO TO 100
    4          CONTINUE
                  QNODV = .TRUE.
               GO TO 100
    5             CONTINUE
               GO TO 100
    6             CONTINUE
               GO TO 100
    7             QCUT   = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in QPGINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in QPGINP')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for QPGCTL:',0)
         IF (QSKIP) THEN
            WRITE (LUPRI,'(A)') ' QPGCTL skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in QPGCTL:',IPRINT
            END IF
            IF (QNODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in QPGCTL.'
            IF (QNODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in QPGCTL.'
            IF (QTEST) WRITE (LUPRI,'(/,2A)') ' Test for quadrupole',
     *      ' moments and quadrupole reorthonormalization.'
            IF (QCUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after QPGCTL.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck qpgini */
      SUBROUTINE QPGINI
C
C     Initialize /CBIQPG/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbiqpg.h>
      IPRINT = IPRDEF
      QNODC   = .FALSE.
      QNODV   = .FALSE.
      QTEST   = .FALSE.
      QSKIP = .NOT. QPGRAD
      QCUT  = .FALSE.
C
      RETURN
      END
C  /* Deck qpgctl */
      SUBROUTINE QPGCTL(WORK,LWORK,PASS)
C
C     Based on DIPCTL by tuh 1985
C     Rewritten for symmetry January 1990, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL PASS
      DIMENSION WORK(LWORK)
#include <cbiqpg.h>
#include <nuclei.h>
#include <inforb.h>
#include <infdim.h>
C
      IF (QSKIP) RETURN
      CALL TIMER('START ',TIMSTR,TIMEND) 
      IF (IPRINT .GT. 0) CALL TITLER('Output from QPGCTL','*',103)
C
      KCMO    = 1
      KDVSP   = KCMO   + NCMOT
      KDV     = KDVSP  + NNASHX
      KQPSOP  = KDV    + N2ASHX
      KQPSO   = KQPSOP + 6*NNBASX
      KQPRHS  = KQPSO  + N2BASX
      KQPREO  = KQPRHS + NVARMA
      KWRK    = KQPREO + 3*NUCDEP
      LWRK    = LWORK  - KWRK + 1
      IF (KWRK .GE. LWORK) CALL STOPIT('QPGCTL',' ',KWRK,LWORK)
      CALL QPGCT1(WORK(KCMO),WORK(KDVSP),WORK(KDV),
     *            WORK(KQPSOP),WORK(KQPSO),
     *            WORK(KQPRHS),WORK(KQPREO),
     *            WORK(KWRK),LWRK)
      IF (IPRINT .GT. 1) CALL TIMER('QPGCTL',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (QCUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after QPGCTL as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in QPGCTL) *****')
      END IF
      RETURN
      END
C  /* Deck qpgct1 */
      SUBROUTINE QPGCT1(CMO,DVSP,DV,SOPACK,SOINT,SECRHS,SECREO,
     *                  WORK,LWORK)
C
C     December 1989, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <iratdef.h>
      PARAMETER (DM1 = -1.0D0)
C
      CHARACTER*8 LABINT(3*MXCOOR)
      DIMENSION CMO(NCMOT), DVSP(NNASHX), DV(NASHT,NASHT),
     *          SOPACK(NNBASX,6), SECRHS(NVARMA), SECREO(3*NUCDEP),
     *          WORK(LWORK), SOINT(NBAST,NBAST)
C
#include <dipole.h>
#include <moldip.h>
C
#include <difsec.h>
#include <cbiqpg.h>
#include <symmet.h>
#include <dorps.h>
#include <abainf.h>
#include <linaba.h>
#include <nuclei.h>
#include <inftap.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <inflin.h>
#include <gdvec.h>
#include <chrxyz.h>
      LOGICAL OLDDX
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from QPGCT1','*',103)
C
C     ***** Read orbitals *****
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
C
C     ***** Read one-electron density *****
C
      IF (NASHT .GT. 0) THEN
         CALL READI(LUSIFC,IRAT*NNASHX,DVSP)
         CALL DSPTSI(NASHT,DVSP,DV)
      ELSE
         READ (LUSIFC)
      END IF
      IF ((IPRINT .GT. 10) .AND. (NASHT .GT. 0)) THEN
         CALL HEADER('Active density matrix in QPGCT1',-1)
         CALL OUTPUT(DV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C     ***** Read SO integrals *****
C
      NCOMP  = 6
      NPATOM = 0
      KINT   = 1
      KREP   = KINT + (9*MXCENT + 1)/IRAT
      KLAST  = KREP + (9*MXCENT + 1)/IRAT
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('QPGCT1','GET1IN',KLAST,LWORK)
      CALL GET1IN(SOPACK,'SECMOM ',NCOMP,WORK(KLAST),LWRK,LABINT,
     &            WORK(KINT),WORK(KREP),IDUMMY,.FALSE.,NPATOM,.TRUE.,
     &            DUMMY,.FALSE.,DUMMY,IPRINT)

      IF (NNBASX .GT. LWORK) CALL STOPIT('QPGCT1','SYMUPK',NNBASX,LWORK)

      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Second moment (xx) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,1),NBAST,1,LUPRI)
         CALL HEADER('Second moment (xy) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,2),NBAST,1,LUPRI)
         CALL HEADER('Second moment (xz) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,3),NBAST,1,LUPRI)
         CALL HEADER('Second moment (yy) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,4),NBAST,1,LUPRI)
         CALL HEADER('Second moment (yz) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,5),NBAST,1,LUPRI)
         CALL HEADER('Second moment (zz) SO  matrix in QPGCT1',-1)
         CALL OUTPAK(SOPACK(1,6),NBAST,1,LUPRI)
      END IF
C
C     ***** Open file for right-hand side of response equations *****
C
      CALL GPOPEN(LUGDR,ABAGDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
C
C     ********************************
C     ***** Loop over components *****
C     ********************************
C
      J = 0
      DO IX = 1, 3
         DO IY = IX, 3
            J = J + 1
            IREP = IBTXOR(ISYMAX(IX,1),ISYMAX(IY,1))
            IXCOOR = IPTXYZ(IX,IREP,1)
            IYCOOR = IPTXYZ(IY,IREP,1)
            IF (IPRINT.GT.5) WRITE (LUPRI,'(1X,A,I5)')' Symmetry ',IREP
C
C           Number of orbital and configuration variables
C
            CALL ABAVAR(IREP+1,.FALSE.,IPRINT,WORK(KLAST),LWRK)
            CALL DSCAL(NNBASX,DM1,SOPACK(1,J),1) 
            CALL DSPTSI(NBAST,SOPACK(1,J),SOINT)
            IF (IPRINT .GT. 5) THEN
               CALL AROUND('Component of second moment:'
     &                           //CHRXYZ(IXCOOR)//CHRXYZ(IYCOOR))
ckr               IF (IPRINT .GT. 10) THEN
                  CALL HEADER('SO matrix in QPGCT1',-1)
                  CALL OUTPUT(SOINT,1,NBAST,1,NBAST,NBAST,NBAST,
     &                 1,LUPRI)
ckr               END IF
            END IF
C
C           ***** Reorthonormalization and right-hand side *****
            CALL OPGCTL(SECREO,SECRHS,CMO,DV,SOINT,SECACT,SECICT,
     &                  WORK,LWORK,IREP,QNODC,QNODV,IPRINT)
C
C           ***** Add reorthonormalization to DSECS *****
C
            IF (.NOT. HELFEY) THEN 
               CALL DCOPY(3*NUCDEP,
     &              SECREO,1,DSECS(IPTAX(IX,1),IPTAX(IY,1),1),9)
               DO ICOOR = 1, 3*NUCDEP
                  DSECS(IPTAX(IY,1),IPTAX(IX,1),ICOOR) = 
     &                 DSECS(IPTAX(IX,1),IPTAX(IY,1),ICOOR)
               END DO
            END IF
C
C           ***** Write right-hand side on file *****
C
            IDISK = 3*NUCDEP + 3 + J
            CALL WRITDX(LUGDR,IDISK,IRAT*NVARPT,SECRHS)
         END DO
      END DO
C
C     ***** Print static constribution to second moment gradient *****
C
      CALL GPCLOSE(LUGDR,'KEEP')
      IF (IPRINT .GT. 1) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         IF (KLAST .GT. LWORK)
     &        CALL STOPIT('QPGCT1','TRANUC',KLAST,LWORK)
         CALL HEADER('Reorthonorm. part of second moment gradient',-1)
         CALL PRISEC(DSECS,'SECDER',WORK(KCSTRA),WORK(KSCTRA))
         CALL HEADER('Static contribution to second moment gradient',-1)
         CALL DZERO(SEC1,27*NUCDEP)
         CALL SECADD(DSECN)
         CALL SECADD(DSECE) 
         CALL SECADD(DSECS) 
         CALL PRISEC(SEC1,'SECDER',WORK(KCSTRA),WORK(KSCTRA))
         CALL DZERO(SEC1,27*NUCDEP)
      END IF
      RETURN
      END
C  /* Deck secadd */
      SUBROUTINE SECADD(AMAT)
#include <implicit.h>
#include <mxcent.h>
      DIMENSION AMAT(3,3,MXCOOR)
#include <nuclei.h>
#include <difsec.h>
      NCOORD = 3*NUCDEP
      DO I = 1, 3
         DO J = 1, 3
            DO K = 1, NCOORD
               SEC1(I,J,K) = SEC1(I,J,K) + AMAT(I,J,K)
            END DO
         END DO
      END DO
      RETURN
      END
C  /* Deck atmdip */
      SUBROUTINE ATMDIP(CSTRA,SCTRA)
C
C     hs 130503
C
C Atomic dipole moments 
C
C according to J. Cioslowski, Phys. Rev. Lett. 62 (1989) 1469-1471
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (THIRD = 1.0D0/3.0D0)
      DIMENSION CSEC(3,3,MXCOOR), CDIP(3,MXCOOR)
      DIMENSION CSTRA(*), SCTRA(*)
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrxyz.h>
#include <difsec.h>
#include <moldip.h>
#include <dipole.h>
#include <inforb.h>
#include <infinp.h>
      DIMENSION DIPATM(3,NUCDEP), DIPEAT(3,NUCDEP)
      DIMENSION GEOM(3*NATOMS), NATTYP(NATOMS)
#include <ibtfun.h>
C
      NCOOR = 3*NUCDEP
C
C     Sum up electronic contributions to second moment gradient
C     and dipole gradient
C     
      CALL DZERO(SEC1,27*NUCDEP)
      CALL SECADD(DSECE) 
      CALL SECADD(DSECS) 
      CALL SECADD(DSECR)
C
      CALL DZERO(DIP1,9*NUCDEP)
      CALL DIPADD(DDIPE)
      CALL DIPADD(DDIPS)
      CALL DIPADD(DDIPR)
C
C     Transform to Cartesian coordinates
C
      CALL TRASEC(SEC1,CSEC,CSTRA,SCTRA,NCOOR,1,1)
      CALL TRADIP(DIP1,CDIP,CSTRA,SCTRA,NCOOR,1,1)
      JATOM = 0
      DO IATOM = 1, NUCIND
         DO ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               ICHARG = NINT(CHARGE(IATOM))
               NATTYP(JATOM + 1) = ICHARG
               DO KCOOR = 1, 3
                  GEOM(3*JATOM + KCOOR) =
     &                 PT(IBTAND(ISYMAX(KCOOR,1),ISYMOP))
     &                 *CORD(KCOOR,IATOM)
               END DO
               JATOM = JATOM + 1
            END IF
         END DO
      END DO
C
C     Add electronic contribution 
C
      MXCNT = NUCDEP
      NELEC = 2*NISHT + NACTEL
      CALL DZERO(DIPATM,3*NUCDEP)
      DO IATOM = 1, MXCNT
         ICOOR = 3*(IATOM - 1)
         DO IXCOOR = 1, 3
            DO IYCOOR = 1, 3
               DIPATM(IXCOOR,IATOM) = DIPATM(IXCOOR,IATOM) 
     &              + THIRD * CSEC(IXCOOR,IYCOOR,ICOOR+IYCOOR)
     &              + THIRD * DIPME(IPTAX(IYCOOR,1)) 
     &                      * CDIP(IXCOOR,ICOOR+IYCOOR) / DFLOAT(NELEC)
            END DO
         END DO
      END DO
C
C     Print total electronic contribution 
C
      CALL HEADER('Total electronic contribution to '//
     &     'atomic dipole moments',-1)
      WRITE (LUPRI,'(/,23X,3(A,13X),/)') 'Ex', 'Ey', 'Ez'
      DO IATOM = 1, MXCNT
         WRITE (LUPRI,'(8X,A6,3F15.8)') NAMDEP(IATOM),
     &        (DIPATM(K,IATOM),K=1,3)
      END DO
C
C     Add nuclear contribution
C      
      DO IATOM = 1, MXCNT
         DO IXCOOR = 1, 3
             DIPATM(IXCOOR,IATOM) = DIPATM(IXCOOR,IATOM)
     &           + NATTYP(IATOM) * GEOM(3*(IATOM-1)+IXCOOR)
         END DO
      END DO
C
C     Print total atomic dipole moments
C
      CALL HEADER('Total atomic dipole moments (au)',-1)
      WRITE (LUPRI,'(/,23X,3(A,13X),/)') 'Ex', 'Ey', 'Ez'
      DO IATOM = 1, MXCNT
         WRITE (LUPRI,'(8X,A6,3F15.8)') NAMDEP(IATOM),
     &        (DIPATM(K,IATOM),K=1,3)
      END DO
C
C     Contribution from atomic charges
C
      CALL HEADER('Contribution from atomic charges to '//
     &     'atomic dipole moments',-1)
      WRITE (LUPRI,'(/,23X,3(A,13X),/)') 'Ex', 'Ey', 'Ez'
      DO IATOM = 1, MXCNT
         WRITE (LUPRI,'(8X,A6,3F15.8)') NAMDEP(IATOM),
     &        (QAPT(IATOM) * GEOM(3*(IATOM-1)+K),K=1,3)
      END DO
C
C     Atomic dipoles without contribution from atomic charges
C
      CALL HEADER('Origin independent electronic contribution '//
     &     'to atomic dipole moments',-1)
      WRITE (LUPRI,'(/,23X,3(A,13X),/)') 'Ex', 'Ey', 'Ez'
      DO IATOM = 1, MXCNT
         DO K = 1, 3
            DIPEAT(K,IATOM) = DIPATM(K,IATOM)
     &           - QAPT(IATOM) * GEOM(3*(IATOM-1)+K)
         END DO
         WRITE (LUPRI,'(8X,A6,3F15.8)') NAMDEP(IATOM),
     &        (DIPEAT(K,IATOM),K=1,3)
      END DO
chs      
      call atmmul(dipatm,dipeat,geom)
chs
      CALL DZERO(SEC1,27*NUCDEP)
      CALL DZERO(DIP1,9*NUCDEP)
      RETURN 
      END
C  /* Deck atmmul */
      SUBROUTINE ATMMUL(DIPATM,DIPEAT,GEOM)
C
C     hs 270503
C
C Molecular multipole moments from atomic charges and dipole moments
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrxyz.h>
#include <difsec.h>
#include <moldip.h>
#include <dipole.h>
#include <inforb.h>
#include <infinp.h>
      DIMENSION DIPA(3), SECA(3,3), THIRDA(3,3,3)
      DIMENSION DIPATM(3,NUCDEP), DIPEAT(3,NUCDEP), GEOM(3*NATOMS)
#include <ibtfun.h>
C
C     Molecular dipole moment
C
      CALL DZERO(DIPA,3)
      DO IATOM = 1, NUCDEP
         DO ICOOR = 1, 3
            DIPA(ICOOR) = DIPA(ICOOR) 
     &           + QAPT(IATOM) * GEOM(3*(IATOM-1)+ICOOR)
         END DO
      END DO
      CALL HEADER('Molecular dipole moment obtained '//
     &     'from atomic charges',-1)
      CALL DP0PRI(DIPA)
C
      CALL DZERO(DIPA,3)
      DO IATOM = 1, NUCDEP
         DO I = 1, 3
            DIPA(I) = DIPA(I) + DIPATM(I,IATOM)
         END DO
      END DO
      CALL HEADER('Molecular dipole moment obtained '//
     &     'from atomic charges and dipoles',-1)
      CALL DP0PRI(DIPA)
C
C     Molecular second moment
C
      CALL DZERO(SECA,9)
      DO IATOM = 1, NUCDEP
         DO I = 1, 3
            DO J = I, 3
               SECA(IPTAX(I,1),IPTAX(J,1)) = 
     &              SECA(IPTAX(I,1),IPTAX(J,1)) 
     &              +QAPT(IATOM)*GEOM(3*(IATOM-1)+I)*GEOM(3*(IATOM-1)+J)
C
               SECA(IPTAX(J,1),IPTAX(I,1)) = 
     &              SECA(IPTAX(I,1),IPTAX(J,1))
            END DO
         END DO
      END DO
      CALL HEADER('Molecular second moment obtained '//
     &     'from atomic charges',-1)
      CALL POLPRI(SECA,'   ',1)
C
      CALL DZERO(SECA,9)
      DO IATOM = 1, NUCDEP
         DO I = 1, 3
            DO J = I, 3
               SECA(IPTAX(I,1),IPTAX(J,1)) = 
     &              SECA(IPTAX(I,1),IPTAX(J,1))
     &              +QAPT(IATOM)*GEOM(3*(IATOM-1)+I)*GEOM(3*(IATOM-1)+J)
     &              +DIPEAT(I,IATOM)*GEOM(3*(IATOM-1)+J)
     &              +DIPEAT(J,IATOM)*GEOM(3*(IATOM-1)+I)
C
               SECA(IPTAX(J,1),IPTAX(I,1)) = 
     &              SECA(IPTAX(I,1),IPTAX(J,1))
            END DO
         END DO
      END DO
      CALL HEADER('Molecular second moment obtained '//
     &     'from atomic charges and dipoles',-1)
      CALL POLPRI(SECA,'   ',1)
C     
C
C     Molecular third moment
C
      CALL DZERO(THIRDA,27)
      DO IATOM = 1, NUCDEP
         DO I = 1, 3
            DO J = I, 3
               DO K = J, 3
                  THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
     &                 + QAPT(IATOM)*GEOM(3*(IATOM-1)+I)
     &                   *GEOM(3*(IATOM-1)+J)*GEOM(3*(IATOM-1)+K) 
C
                  THIRDA(IPTAX(I,1),IPTAX(K,1),IPTAX(J,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(J,1),IPTAX(I,1),IPTAX(K,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(J,1),IPTAX(K,1),IPTAX(I,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(K,1),IPTAX(I,1),IPTAX(J,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(K,1),IPTAX(J,1),IPTAX(I,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               END DO
            END DO
         END DO
      END DO
      CALL HEADER('Molecular third moment obtained '//
     &     'from atomic charges',-1)
      CALL PRIOCT(THIRDA)
C
      CALL DZERO(THIRDA,27)
      DO IATOM = 1, NUCDEP
         DO I = 1, 3
            DO J = I, 3
               DO K = J, 3
                  THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
     &                 + QAPT(IATOM)*GEOM(3*(IATOM-1)+I)
     &                   *GEOM(3*(IATOM-1)+J)*GEOM(3*(IATOM-1)+K)
     &                 + DIPEAT(I,IATOM)*GEOM(3*(IATOM-1)+J)
     &                   *GEOM(3*(IATOM-1)+K)
     &                 + DIPEAT(J,IATOM)*GEOM(3*(IATOM-1)+I)
     &                   *GEOM(3*(IATOM-1)+K)
     &                 + DIPEAT(K,IATOM)*GEOM(3*(IATOM-1)+I)
     &                   *GEOM(3*(IATOM-1)+J)
C
                  THIRDA(IPTAX(I,1),IPTAX(K,1),IPTAX(J,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(J,1),IPTAX(I,1),IPTAX(K,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(J,1),IPTAX(K,1),IPTAX(I,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(K,1),IPTAX(I,1),IPTAX(J,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
                  THIRDA(IPTAX(K,1),IPTAX(J,1),IPTAX(I,1)) = 
     &                 THIRDA(IPTAX(I,1),IPTAX(J,1),IPTAX(K,1))
               END DO
            END DO
         END DO
      END DO
      CALL HEADER('Molecular third moment obtained '//
     &     'from atomic charges and dipoles',-1)
      CALL PRIOCT(THIRDA)
      RETURN
      END
