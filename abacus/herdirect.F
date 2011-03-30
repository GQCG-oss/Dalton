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
/* Comdeck herdirect_log */
941206-hjaaj
HERPAR: new routine for defining parallel hermit calc. of Fock matr.
HERFCK: call PARDRV if PARCAL true
        changed STATIS to KPRINT.gt.0 (STATIS missing in newest
        aba.cdk from oslo)
941205-from oslo:
HERFCK: also IPRINT = IPRFCK for PARCAL true (was = 0)
941103-hjaaj:
implemented ISYMDM(NDMAT), IFCTYP(NDMAT)
HERDI1,HERDI2: excluded; not used, are also incomplete
ABARUN: moved to here from abadrv.u (needed for HERMIT)
931026-hjaaj: sether now calls READIN; + consistency check
HERFCK: new parameter (print level in CBITWO)
931025-hjaaj: sirdirn04.u (i.e. HERFCK,SETHER) added to herdirect.u
GETFCK now calls HERFCK
AOLABE renamed to AOLAB4, generalized, and moved to GP.
GETFCK is now only internal test routine; HERFCK is called from the
outside (SIRIUS and RESPONSE).
921201-hjaaj:
   new Sirius input routine for Hermit parameters ? an extended sether?
   also direct calculation of solvent integrals ?
===========================================================================
#endif
C  /* Deck tstdir */
      SUBROUTINE TSTDIR(FMAT,DMAT,NDMAT,WORK,LWORK,IPRINT)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C     PURPOSE : Driver routine for the test calculation of the
C               two-electron part of the fock matrices.
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (MAXDIS=500)
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
#include <inforb.h>
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      KGMAT = 1
      KINDX = KGMAT + NBAST*NBAST*MAXDIS
      KEND  = KINDX + (MAXDIS + 1)/IRAT
      LFREE = LWORK - KEND
      IF (KEND .GT. LWORK) CALL STOPIT('TSTDIR',' ',KEND,LWORK)
C
      NUMDIS = -1
    1 CONTINUE
C
C     Get distributions
C
      CALL GETDIS(WORK(KGMAT),WORK(KINDX),NUMDIS,MAXDIS,WORK(KEND),
     &            LFREE)
      IF (NUMDIS .EQ. -1) GOTO 2
C
C     Add contributions to Fock matrix
C
      CALL DISFCK(FMAT,DMAT,NDMAT,WORK(KGMAT),WORK(KINDX),NUMDIS)
C
      GOTO 1
C
    2 CONTINUE
C
C     Fock matrix has now been constructed
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.2) THEN
         CALL HEADER('Density and Fock matrices from '//
     &               'distributions in TSTDIR',-1)
         DO 300 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
      CALL TIMER('GETDIS',TIMSTR,TIMEND)
      RETURN
      END
C  /* Deck disfck */
      SUBROUTINE DISFCK(FMAT,DMAT,NDMAT,GMAT,INDEX,NUMDIS)
C
C     Written by Henrik Koch and Trygve Helgaker 27-November-1991.
C
C
C     PURPOSE : Calculate fock matrices using distributions.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (ONE=1.0D00, HALF=0.5D00, FOURTH=0.25D00)
      INTEGER   P,Q,R,S
      DIMENSION FMAT(NBAST,NBAST,NDMAT), DMAT(NBAST,NBAST,NDMAT)
      DIMENSION GMAT(NBAST,NBAST,NUMDIS), INDEX(NUMDIS)
#include <inforb.h>
C
      DO 100 I = 1,NUMDIS
C
         CALL UNPKIJ(INDEX(I),R,S)
C
         IF (R .EQ. S) THEN
            FACCOU = HALF
            FACEXC = FOURTH
         ELSE
            FACCOU = ONE
            FACEXC = HALF
         END IF
C
         DO 110 J = 1,NDMAT
C
            VALUE = FACCOU*DDOT(NBAST*NBAST,DMAT(1,1,J),1,GMAT(1,1,I),1)
            FMAT(R,S,J) = FMAT(R,S,J) + VALUE
            FMAT(S,R,J) = FMAT(S,R,J) + VALUE
C
            DO 120 Q = 1,NBAST
C
               FMAT(R,Q,J) = FMAT(R,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,S,J),1,GMAT(1,Q,I),1)
               FMAT(S,Q,J) = FMAT(S,Q,J)
     &                  - FACEXC*DDOT(NBAST,DMAT(1,R,J),1,GMAT(1,Q,I),1)
C
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck getdis */
      SUBROUTINE GETDIS(GMAT,INDXAB,NUMDIS,MAXDIS,WORK,LWORK)
C
C     Written by Henrik Koch and Trygve Helgaker 26-November-1991.
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               distributions (**|cd).
C               The distributions are stored as full squares and without
C               symmetry reduction.
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
      LOGICAL ABA
      DIMENSION GMAT(NBAST,NBAST,MAXDIS), INDXAB(MAXDIS)
      DIMENSION WORK(LWORK)
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
C
      CALL ABARUN(ABA)
C
      I2TYP = 0
C
C-------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C-------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
         NODV = .FALSE.
         NOPV = .FALSE.
         IF (NASHT .EQ. 0) NODV = .TRUE.
         IF (NASHT .LT. 2) NOPV = .TRUE.
         IPRALL =  0
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      MAXDIF = 0
      ITYPE  = 4
      KJSTRS = 1
      KNPRIM = KJSTRS + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KNCONT = KNPRIM + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KIORBS = KNCONT + (MXSHEL*MXAOVC*2 + 1)/IRAT
      KJORBS = KIORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KKORBS = KJORBS + (MXSHEL*MXAOVC + 1)/IRAT
      KLAST  = KKORBS + (MXSHEL*MXAOVC + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('GETDIS','PAOVEC',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRALL)
      KLAST = KJORBS
      LWRK  = LWORK - KLAST + 1
C
C----------------------------
C     Calculate distributions
C----------------------------
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL TWOINT(WORK(KLAST),LWRK,DUMMY,DUMMY,DUMMY,NDMAT,
     &            IDUMMY,IDUMMY,GMAT,
     &            INDXAB,NUMDIS,MAXDIS,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &            TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &            IDUMMY,I2TYP,WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &            WORK(KIORBS),
     &            IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,.FALSE.)
      IF (IPRINT.GT.0) CALL TIMER('GETDIS/TWOINT',TIMSTR,TIMEND)
C
C--------------------------------------
C     Write out integral distributions
C--------------------------------------
C
      IF (IPRINT.GT.3) THEN
         CALL HEADER('Distribution matrices in GETDIS',-1)
         DO 300 I = 1, NUMDIS
            CALL UNPKIJ(INDXAB(I),IA,IB)
            WRITE(LUPRI,'(//A,2I5)')' Distribution matrix No.',I,NUMDIS
            WRITE(LUPRI,'(A,2I5)')  ' Orbital indices:       ',IA,IB
            CALL OUTPUT(GMAT(1,1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck fcktes */
#if !defined (VAR_GETFCK)
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
C     April 1995: FCKTES below must be changed to call new Fock matrix
C     routines instead of GETFCK.
      CALL QUIT('FCKTES not implemented in this version!')
      END
#else
      SUBROUTINE FCKTES(WORK,LWORK,MAXDIF,NODV,NOPV,NOCONT,TKTIME,
     &                  IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
      LOGICAL NODV, NOPV, NOCONT, TKTIME, RETUR, DIRECT
      DIMENSION WORK(LWORK)
#include <inforb.h>
C
      IF (NODV) THEN
         NDMAT = 1
      ELSE
         NDMAT = 2
      END IF
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),
     &            WORK(KDMAT+N2BASX),IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      DIRECT = .TRUE.
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL GETFCK(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,DIRECT,
     *          .FALSE.)
      CALL TIMER('FCKTES/GETFCK',TIMSTR,TIMEND)
      if (.true.) return
C
C     Test distributions
C
      KFMAT = 1
      KDMAT = KFMAT + N2BASX*NDMAT
      KCMO  = KDMAT + N2BASX*NDMAT
      KDV   = KCMO  + NCMOT
      KDTSO = KDV   + NNASHX
      KDASO = KDTSO + NNBAST
      KLAST = KDASO + NNBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FCKTES',' ',KLAST,LWORK)
      CALL ONEDSF(WORK(KCMO),WORK(KDV),WORK(KDMAT),WORK(KDMAT+N2BASX),
     &            IPRINT,NODC,NODV)
      CALL DAXPY(N2BASX,-D1,WORK(KDMAT+N2BASX),1,WORK(KDMAT),1)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KDTSO),NNBAST)
      CALL SQDENS(WORK(KDTSO),WORK(KFMAT),IPRINT)
      IF (.NOT.NODV) CALL DZERO(WORK(KFMAT+N2BASX),N2BASX)
C
      KLAST = KCMO
      LWRK  = LWORK - KLAST + 1
      CALL TSTDIR(WORK(KFMAT),WORK(KDMAT),NDMAT,WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
#endif
C  /* Deck sqdens */
      SUBROUTINE SQDENS(DTSO,DMAT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DTSO(*), DMAT(NBAST,NBAST)
#include <inforb.h>
C
      CALL DZERO(DMAT,NBAST*NBAST)
      ISOFF = 0
      DO 100 ISYM = 1, NSYM
         DO 200 I = 1, NBAS(ISYM)
            INDA = IBAS(ISYM) + I
            DO 300 J = 1, I
               INDB = IBAS(ISYM) + J
               DMAT(INDA,INDB) = DMAT(INDA,INDB) + DTSO(ISOFF+J)
               DMAT(INDB,INDA) = DMAT(INDA,INDB)
 300        CONTINUE
            ISOFF = ISOFF + I
 200     CONTINUE
 100  CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output of squared matrix from SQDENS','*',103)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck sether */
      SUBROUTINE SETHER(JPRINT,NEWGEO,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  SYMMET: MAXREP,NAOS(8))
C  INFORB: NSYM,NBAS(8)
C
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <cbirea.h>
#include <symmet.h>
#include <inforb.h>
C
      LOGICAL SFIRST, RELCAL, TSTINP, HERMIT, NEWGEO
      SAVE    SFIRST
      DATA    SFIRST /.TRUE./
C
      CALL QENTER('SETHER')
C
C     Initialize /CBIREA/ and call readin
C
      IF (SFIRST .OR. NEWGEO) THEN
chj   May 09: should we really call READIN here now,
chj   when they always are called in the beginning, also when RNHERM false??
Chj   June 09: maybe READIN may still be needed for NEWGEO, keep it for
Chj   now for safety, but REAINI should definitely NOT be called, it will
Chj   reiinitialize variables changed by user under *READIN
C
         IPREAD = JPRINT
         HERMIT = .FALSE.
         CALL READIN(WORK,LWORK,HERMIT)
C
C        consistency check
C
         IF (MAXREP+1 .NE. NSYM) THEN
            WRITE (LUPRI,'(2(/A,I5))')
     &      ' SETHER fatal error: NSYM  from SIRIUS is',NSYM,
     &      '                  MAXREP+1 from READIN is',MAXREP+1
            CALL QUIT('SETHER error: NSYM .ne. MAXREP+1')
         END IF
         NERR = 0
         DO 100 ISYM = 1,NSYM
            IF (NAOS(ISYM) .NE. NBAS(ISYM)) NERR = NERR + 1
  100    CONTINUE
         IF (NERR .GT. 0) THEN
            WRITE (LUPRI,'(/A/A/)')
     &      ' SETHER fatal error: NBAS(:) from LUONEL',
     &      '   is different from NAOS(:) from READIN.'
            WRITE (LUPRI,'(A,8I5)') ' NBAS(:) =',(NBAS(I),I=1,NSYM)
            WRITE (LUPRI,'(A,8I5)') ' NAOS(:) =',(NAOS(I),I=1,NSYM)
            CALL QUIT('SETHER error: NBAS(:) .ne. NAOS(:)')
         END IF
         SFIRST = .FALSE.
      END IF
      CALL QEXIT('SETHER')
      RETURN
      END
C  /* Deck herfck */
      SUBROUTINE HERFCK(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                  WORK,LWORK)
C
C     PURPOSE : Driver routine for direct calculation of the
C               two-electron part of the fock matrices.
C               We assume the densities and fock matrices are full
C               squares and without symmetry reduction .
C
C     10-Apr-1996-hjaaj: select if we want to use HRFCK1 or ERIFCK
C     (old integral routines or new ERI integral routines)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      LOGICAL USEERI, LDUMMY
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION FMAT(*), DMAT(*), ISYMDM(*), IFCTYP(*), 
     &          WORK(*)
C
C Used from common blocks:
C   gnrinf: SEGBAS, PARCAL
C   infpar: NODTOT
C
#include <infpar.h>
#include <inforb.h>
#include <gnrinf.h>
#include <abainf.h>
#include <blocks.h>
#include <cbihr2.h>
#include <cbieri.h>
#include <veclen.h>
C
#include <memint.h>
C
      CALL QENTER('HERFCK')
C
C     The new ERI routines have not been parallelized yet.
C     RUNERI is the flag from input processing, to tell that we should use eri,
C     USEERI is a temporary flag to check if we actually can in the present 
C     passage.
CNECgh980808: The construction of the Fock Matrix using the Temporary
C matrices leads to huge bank conflicts in case of an even NBASIS. 
C We therefore assure that the first dimension of FCKTMP is always odd.
C Changes for eri: HERFCK (allocation), HRFCK2 (zero out and summing)
C and FCKCON (dimensioning). Analog hermit. NODD added to nuclei.h & aobtch.h
C We need to include nuclei.h here, to have NBASIS and NODD.
C Initalization of NODD in herrdn.F and eriaob.F
C
      USEERI = RUNERI .AND.
     &        .NOT.SRINTS .AND. CHIVAL.EQ.-1.D0 .AND. PANAS.EQ.D0
      DO I = 1, NDMT
         USEERI = USEERI .AND. ((IFCTYP(I).EQ.13).OR.(IFCTYP(I).EQ.3))
      END DO
      IF (RUNERI .AND. .NOT.USEERI) WRITE(LUPRI,'(/A/A)')
     &   'INFO HERFCK: Using Hermit instead of ERI for 2-el. ints.',
     &   'because requested calculation is not implemented in ERI.'
#if defined (VAR_VECTOR)
      ICHUNK = MAX(IVECLN/NDMT,1)
#endif
C
      IF (.NOT. USEERI) THEN
         CALL HRFCK1(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &               ICEDIF,IFTHRS,WORK(KFREE),LFREE)
      ELSE
#if defined (VAR_MPI) || defined (VAR_PVM)
         IF (PARCAL) THEN
            if (iprfck .gt. 1) write (lupri,*) '- using parallel ERI'
            CALL MEMGET('INTE',KNSTAT,NODTOT,WORK,KFREE,LFREE)           
            IPRTYP = 2
            CALL PARDRV(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,WORK(KFREE),
     &                  WORK(KNSTAT),DUMMY,LFREE,IDUMMY,IDUMMY,IDUMMY,
     &                  LDUMMY,LDUMMY,LDUMMY,LDUMMY,LDUMMY,
     &                  IPRFCK,IPRTYP,
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY)

         ELSE
#endif
           CALL MEMGET('REAL',KCCFBT,MXPRIM*MXCONT  ,WORK,KFREE,LFREE)
           CALL MEMGET('INTE',KINDXB,8*MXSHEL*MXCONT,WORK,KFREE,LFREE)
#if defined (VAR_VECTOR)
           CALL MEMGET('REAL',KFCKTP,ICHUNK*NDMT*(NBASIS+NODD)*NBASIS,
     &                 WORK,KFREE,LFREE)
#else
           KFCKTP = KFREE
#endif
           if (iprfck .gt. 1) write (lupri,*) '- using sequential ERI'
           CALL ERIFCK(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                 WORK(KCCFBT),WORK(KINDXB),WORK(KFCKTP),
     &                 WORK(KFREE),LFREE)
#if defined (VAR_MPI) || defined (VAR_PVM)
         END IF
#endif
         CALL MEMREL('ERIFCK',WORK,KWORK,KWORK,KFREE,LFREE)
      END IF
      CALL QEXIT('HERFCK')
      RETURN
      END
C  /* Deck hrfck1 */
      SUBROUTINE HRFCK1(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,IPRFCK,
     &                  ICEDIF,IFTHRS,WORK,LWORK)
C
C     Written by Henrik Koch and Trygve Helgaker 22-November-1991.
C     March 1997 - tsaue Screening
C
C     PURPOSE : Driver routine for direct calculation of the
C               two-electron part of the fock matrices.
C               We assume the densities and fock matrices are full
C               squares and without symmetry reduction .
C
C               Common blocks needed to be initialized before this
C               subroutine is called:
C
C                      *CALL CBITWO   from TWOINP
C                      *CALL INFORB   from sirius
C                      *CALL GAMCOM   from READIN
C                      *CALL CCOM     from READIN
C                      *CALL PINCOM   from READIN
C                      *CALL PRIMIT   from READIN
C                      *CALL XYZPOW   from READIN
C                      *CALL SHELLS   from READIN
C                      *CALL NUCLEI   from READIN
C                      *CALL SYMMET   from READIN
C                      *CALL DORPS    from abatro.u
C                       ( only DOREPS and DOCOOR used)
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
      LOGICAL   ABA
      PARAMETER (D0 = 0.D0)
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT),
     &          ISYMDM(NDMAT), IFCTYP(NDMAT)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  INFORB: N2BASX,NBAST,?
C  CBITWO: IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,?
C  GNRINF: PARCAL
C  INFPAR: NODTOT
C
#include <abainf.h>
#include <dftcom.h>
#include <inforb.h>
#include <cbitwo.h>
#include <dorps.h>
#include <symmet.h>
#include <gnrinf.h>
#include <infpar.h>
#include <blocks.h>
C
#include <memint.h> 
      IF (NDMAT .LE. 0) RETURN  
      CALL QENTER('HRFCK1')
      I2TYP = 0
C
C     Inquire if this is an ABACUS run or not:
C
      CALL ABARUN(ABA)
C
C
C----------------------------------------------------------------
C     Setup information for the two-electron integral routines.
C----------------------------------------------------------------
C
      IF (.NOT.ABA) THEN
         NODV = (NASHT .EQ. 0)
         NOPV = (NASHT .LT. 2)
C
C        Remember to change this!!!!!!!
C
         IPRALL =  0
C
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
C
         IPRINT = IPRFCK
         IPRNTA = 0
         IPRNTB = 0
         IPRNTC = 0
         IPRNTD = 0
         RETUR  = .FALSE.
         NOCONT = .FALSE.
         TKTIME = .FALSE.
      END IF
C
      I2TYP  = 0
      ITYPE  = 3
      MAXDER = 0
      NPAO = MXSHEL*MXAOVC
      CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,.FALSE.,IPRINT)
      CALL MEMREL('HRFCK1.PAOVEC',WORK,KWORK,KJORBS,KFREE,LFREE)
C
C     Transform density to AO basis
C     Modified july-97, kr
C     FMAT still empty, and we use this array to store DMAT in AO-basis 
C     temporarily
C
      DO 200 IMAT = 1,NDMAT
         CALL DSOTAO(DMAT(1,IMAT),FMAT(1,1),NBAST,ISYMDM(IMAT),IPRINT)
         CALL DCOPY(N2BASX,FMAT(1,1),1,DMAT(1,IMAT),1)
  200 CONTINUE
      CALL DZERO(FMAT(1,1),NDMAT*N2BASX)
C
C     Prepare for screening
C     =====================
C
      N2DRAO = NSYMBL*NSYMBL*NDMAT
      N2GAB = NSYMBL*NSYMBL      
      CALL MEMGET('REAL',KGAB ,N2GAB, WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDRAO,N2DRAO,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDSKP,8     ,WORK,KFREE,LFREE)
C     
C     Get GAB-matrix (Cauchy-Schwarz)
C
      IF (SRINTS) THEN
         CALL GETGAB(WORK(KGAB),'GABSRXXX',N2GAB,
     &               WORK(KFREE),LFREE)
      ELSE
         CALL GETGAB(WORK(KGAB),'GABAOXXX',N2GAB,
     &               WORK(KFREE),LFREE)
      END IF
C
C     Make reduced density matrix for screening
C
      CALL MKDRAO(DMAT,WORK(KDRAO),NDMAT,
     &            WORK(KFREE),LFREE,IPRINT)
C
C--------------------------------
C     Calculate fock matrices.
C--------------------------------
C
      CALL GETTIM(CPU1,WALL1)
#if defined (VAR_MPI) || defined (VAR_PVM)
      IF (PARCAL) THEN
         IATOM = 0
         CALL MEMGET('INTE',KNSTAT,NODTOT,WORK,KFREE,LFREE)
C
         IPRTYP = 1
         HFXM0 = HFXMU
         HFXMU = D0
         CALL PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,
     &               WORK(KFREE),WORK(KNSTAT),DUMMY,LFREE,ITYPE,MAXDER,
     &               IATOM,NODV,NOPV,NOCONT,TKTIME,RETUR,IPRINT,IPRTYP,
     &               ICEDIF,IFTHRS,
     &               WORK(KGAB),WORK(KDRAO),DUMMY,WORK(KDSKP))
         IF (IPRINT .GT. 0) THEN
           CALL GETTIM(CPU2,WALL2)
           WALL   = WALL2 - WALL1
           CALL SCRSTA('HRFCK1/PARDRV',WORK(KDSKP),WALL)
         END IF

         HFXMU = HFXM0
C
C        Calculate attenuated exchange (hfxatt*erf(hfxmu*r12)/r12) 
C
         IF (HFXMU.NE.D0) THEN
            IPRTYP = 1
            IATOM  = 1
            HFXFC0 = HFXFAC
            HFXFAC = HFXATT 
            CALL MEMGET('REAL',KFCTMP,NDMAT*N2BASX,WORK,KFREE,LFREE)
c           We have to use a temporary array because PARDRV
c           overwrites the matrix from the phase one destroying 
c           entirely the Coulomb contribution.
c           This does not solve problems for LR phase of a QR
c           calculation.

            CALL DCOPY(NDMAT*N2BASX,FMAT,1,WORK(KFCTMP),1)
            CALL PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,
     &               WORK(KFREE),WORK(KNSTAT),DUMMY,LFREE,ITYPE,MAXDER,
     &               IATOM,NODV,NOPV,NOCONT,TKTIME,RETUR,IPRINT,IPRTYP,
     &               ICEDIF,IFTHRS,
     &               WORK(KGAB),WORK(KDRAO),DUMMY,WORK(KDSKP))
            IF (IPRINT .GT. 0) THEN
              CALL GETTIM(CPU1,WALL1)
              WALL   = WALL1 - WALL2
              CALL SCRSTA('HRFCK1/PARDRV-hfxatt',WORK(KDSKP),WALL)
            END IF
            CALL DAXPY(NDMAT*N2BASX,1D0,WORK(KFCTMP),1,FMAT,1)
            HFXFAC = HFXFC0
         END IF

         CALL MEMREL('HRFCK1.PARDRV',WORK,KWORK,KNSTAT,KFREE,LFREE)
      ELSE
#endif
C
C        Calculate full Coulomb and exchange weighted with hfxfac 
C
         CALL TIMER('START ',TIMSTR,TIMEND)
         HFXM0 = HFXMU
         HFXMU = D0
         CALL TWOINT(WORK(KFREE),LFREE,DUMMY,FMAT,DMAT,NDMAT,
     &               ISYMDM,
     &               IFCTYP,DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDER,0,
     &               NODV,NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,
     &               IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,I2TYP,
     &               WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &               WORK(KIORBS),ICEDIF,IFTHRS,
     &               WORK(KGAB),WORK(KDRAO),DUMMY,WORK(KDSKP),.FALSE.)
         HFXMU = HFXM0
         IF (IPRINT .GT. 0) THEN
           CALL GETTIM(CPU2,WALL2)
           WALL   = WALL2 - WALL1
           CALL SCRSTA('HRFCK1/TWOINT',WORK(KDSKP),WALL)
         END IF
C
C        Calculate attenuated exchange (hfxatt*erf(hfxmu*r12)/r12) 
C
         IF (HFXMU.NE.D0) THEN
            HFXFC0 = HFXFAC
            HFXFAC = HFXATT 
            CALL TWOINT(WORK(KFREE),LFREE,DUMMY,FMAT,DMAT,NDMAT,
     &               ISYMDM,
     &               IFCTYP,DUMMY,IDUMMY,NUMDIS,1,ITYPE,MAXDER,0,
     &               NODV,NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,
     &               IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,I2TYP,
     &               WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &               WORK(KIORBS),ICEDIF,IFTHRS,
     &               WORK(KGAB),WORK(KDRAO),DUMMY,WORK(KDSKP),.FALSE.)
            IF (IPRINT .GT. 0) THEN
              CALL GETTIM(CPU1,WALL1)
              WALL   = WALL1 - WALL2
              CALL SCRSTA('HRFCK1/TWOINT-hfxatt',WORK(KDSKP),WALL)
            END IF
            HFXFAC = HFXFC0
         END IF
         IF (IPRINT .GT. 0) CALL TIMER('HRFCK1/TWOINT',TIMSTR,TIMEND)
#if defined (VAR_MPI) || defined (VAR_PVM)
      END IF
#endif
      CALL SKLFCK(FMAT,DUMMY,WORK(KFREE),LFREE,IPRINT,.TRUE.,
     &            .FALSE.,.FALSE.,.FALSE.,NODV,MAXDER,.FALSE.,
     &            NDMAT,ISYMDM,IFCTYP,0,.FALSE.)
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.5) THEN
         CALL HEADER('Density and Fock matrices in HRFCK1',-1)
         DO 300 I = 1, NDMAT
            WRITE (LUPRI,'(//A,I3,A,I3)')
     &         ' Density matrix No.',I,' of',NDMAT
            CALL OUTPUTx(DMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A,I3,A,I3)')
     &         ' Fock matrix No.',I,' of',NDMAT
            CALL OUTPUTx(FMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
C
      CALL MEMREL('HRFCK1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('HRFCK1')
      RETURN
      END
C  /* Deck aindex */
      SUBROUTINE AINDEX(ISHELA,NAINTS,INDEXA,DOINDX,IORBSH,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
      LOGICAL DOINDX
      DIMENSION INDEXA(*), IORBSH(MXSHEL,MXAOVC)
#include <blocks.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IORB = 0
      MULA   = ISTBSH(ISHELA)
      DO 100 IORBA = 1, NORBSH(ISHELA)
         DO 200 ICOMPA = 1, KHKTSH(ISHELA)
            ITYNA = ISYMAO(NHKTSH(ISHELA),ICOMPA)
            IADR  = IORBSB(IORBSH(ISHELA,1))
     &            + KHKTSH(ISHELA)*(IORBA - 1) + ICOMPA
            DO 300 IREPA = 0, MAXREP
               IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                  IORB = IORB + 1
                  IF (DOINDX) INDEXA(IORB) = IPTSYM(IADR,IREPA)
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      NAINTS = IORB
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AINDEX',2)
         IF (.NOT.DOINDX) THEN
            WRITE (LUPRI,'(2X,A,I5)') ' Number of a orbitals:',NAINTS
         ELSE
            WRITE (LUPRI,'(2X,I3,A,8I5/,(40X,8I5))')
     &          NAINTS,' distributions in this TWOINT call:',
     &            (INDEXA(I),I=1,NAINTS)
         END IF
      END IF
      RETURN
      END
C  /* Deck abarun */
      SUBROUTINE ABARUN(RUNABA)
#include <implicit.h>
      LOGICAL ABA, RUNABA
      SAVE ABA
      DATA ABA /.FALSE./
      RUNABA = ABA
      RETURN
      ENTRY ABUSET
      ABA = .FALSE.
      RETURN
      ENTRY ABASET
      ABA = .TRUE.
      RETURN
      END
