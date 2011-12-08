C 
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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
C     BLKCLS(WORK(KLCIBT),NBLOCKS,WORK(KLBLKCLS),WORK(KLSPSPCL),
C    &            NOCTPA,NOCTPB)
      SUBROUTINE BLKCLS(IBLKS,NBLKS,IBLKCLS,ISPSPCL,
     &                  NCLS,LCLS,NOCTPA,NOCTPB)
*
* Class of each block, and dimension of each class
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLKS(8,NBLKS)
      INTEGER ISPSPCL(NOCTPA,NOCTPB)
*. Output
      INTEGER IBLKCLS(NBLKS),LCLS(NCLS)
*
C?    WRITE(6,*) ' ISPSPCL'
C?    CALL IWRTMA(ISPSPCL,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      IZERO = 0
      CALL ISETVC(LCLS,IZERO,NCLS)
      DO JBLK = 1, NBLKS
        IICLS = ISPSPCL(IBLKS(1,JBLK),IBLKS(2,JBLK))
        IBLKCLS(JBLK) = IICLS
        LCLS(IICLS) = LCLS(IICLS) + IBLKS(8,JBLK)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' BLKCLS Speaking '
        WRITE(6,*) ' ==============='
        WRITE(6,*)
        WRITE(6,*) ' Dimension of each class '
        CALL IWRTMA(LCLS,1,NCLS,1,NCLS)
        WRITE(6,*)
        WRITE(6,*) ' Class of each block : '
        CALL IWRTMA(IBLKCLS,1,NBLKS,1,NBLKS)
      END IF
*
      RETURN
      END
***********************************************************************
*
* Lucia.f : GAS implementing no pair relativistic Theory
*
* Version of May 97 , Jeppe Olsen
*
      SUBROUTINE CHK_ORBDIM(IGSFILL,NIRREP)
*
* Insert dimensions of orbital space IGSFILL
* Check number of shells in NGSSH with info from ENVIRONMENT
*
* Environment info must be available
*
* Jeppe Olsen, Feb. 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*
      IF(IGSFILL.NE.0) THEN
*. Fill GAS shell IGSFILL with remaining orbitals
       DO IRREP = 1, NIRREP
         LMO = 0
         DO IGAS = 1, NGAS
           IF(IGAS.NE.IGSFILL) LMO = LMO + NGSSH(IRREP,IGAS)
         END DO
         NGSSH(IRREP,IGSFILL) = NMOS_ENV(IRREP)-LMO
       END DO
      END IF
*. Make sure that no dimensions are negative
      LERROR = 0
      DO IGAS = 1, NGAS
       DO IRREP = 1, NIRREP
         IF(NGSSH(IGAS,IRREP).LT.0) THEN
           WRITE(6,*)
     &     ' Error : negative orbital dimension,IGAS,IRREP,N ',
     &     IGAS,IRREP,NGSSH(IGAS,IRREP)
           LERROR = LERROR + 1
         END IF
       END DO
      END DO
*. Make sure that all dimensions add correctly up
      DO IRREP = 1, NIRREP
        LMO = 0
        DO IGAS = 1, NGAS
          LMO = LMO + NGSSH(IRREP,IGAS)
        END DO
        IF(LMO.LT.NMOS_ENV(IRREP)) THEN
          WRITE(6,*)
     &    ' Warning : Number of orbitals in irrep reduced compared'
          WRITE(6,*)
     &    ' with information from environment, IRREP,NMO,NMO_ENV'
          WRITE(6,'(3I5)') IRREP,LMO,NMOS_ENV(IRREP)
        ELSE IF(LMO.NE.NMOS_ENV(IRREP)) THEN
          WRITE(6,*)
     &    ' Error : Number of orbitals in irrep not consistent'
          WRITE(6,*)
     &    ' with information from environment, IRREP,NMO,NMO_ENV'
          WRITE(6,'(3I5)') IRREP,LMO,NMOS_ENV(IRREP)
          LERROR = LERROR + 1
        END IF
      END DO
*
      IF(LERROR.NE.0) THEN
        WRITE(6,*) ' Problem with orbital dimensions'
        Call Abend2(       ' Problem with orbital dimensions' )
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE DXTYP2_GAS(NDXTP,ITP,JTP,KTP,LTP,
     &                     NOBTP,IL,IR,IPHGAS)
*
* Obtain types of I,J,K,L so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(J).ge.type(L)
* are included
*
* Intermediate occupations less than zero allowed for particle spaces
* (IPHGAS=2)
*
*
      INTEGER IL(NOBTP),IR(NOBTP),IPHGAS(NOBTP)
      INTEGER ITP(*),JTP(*),KTP(*),LTP(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DXTYP_GAS in action '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Occupation of left string '
        CALL IWRTMA(IL,1,NOBTP,1,NOBTP)
        WRITE(6,*) ' Occupation of right string '
        CALL IWRTMA(IR,1,NOBTP,1,NOBTP)
      END IF
*
*. Number of differing occupations
      NANNI = 0
      NCREA = 0
      NDIFT = 0
*
      ICREA1 = 0
      ICREA2 = 0
      IANNI1 = 0
      IANNI2 = 0
      DO IOBTP = 1, NOBTP
        NDIFT = NDIFT + ABS(IL(IOBTP)-IR(IOBTP))
        NDIF = IL(IOBTP)-IR(IOBTP)
        IF(NDIF.EQ.2) THEN
*. two electrons of type IOBTP must be created
          ICREA1 = IOBTP
          ICREA2 = IOBTP
          NCREA = NCREA + 2
        ELSE IF (NDIF .EQ. -2 ) THEN
*. Two electrons of type IOBTP must be annihilated
          IANNI1 = IOBTP
          IANNI2 = IOBTP
          NANNI = NANNI + 2
        ELSE IF (NDIF.EQ.1) THEN
*. one electron of type IOBTP must be created
          IF(NCREA.EQ.0) THEN
            ICREA1 = IOBTP
          ELSE
            ICREA2 = IOBTP
          END IF
          NCREA = NCREA + 1
        ELSE IF (NDIF.EQ.-1) THEN
* One electron of type IOBTP must be annihilated
          IF(NANNI.EQ.0) THEN
            IANNI1 = IOBTP
          ELSE
            IANNI2 = IOBTP
          END IF
          NANNI = NANNI + 1
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NCREA, NANNI ', NCREA, NANNI
        WRITE(6,*)  ' ICREA2, IANNI2 ', ICREA2,IANNI2
C       WRITE(6,*)  ' ICREA11,IANNI11 ', ICREA11,IANNI11
C       WRITE(6,*)  ' ICREA21,IANNI21 ', ICREA21,IANNI21
      END IF
*
      NDXTP = 0
      IF(NDIFT.GT.4) THEN
        NDXTP = 0
      ELSE
      IF(NCREA.EQ.0.AND.NANNI.EQ.0) THEN
*. strings identical, include diagonal excitions  itp = jtp, ktp=ltp
        DO IJTP = 1, NOBTP
          IF(IR(IJTP).GE.1.OR.IPHGAS(IJTP).EQ.2) THEN
            DO KLTP = 1, IJTP
             IF((IJTP.NE.KLTP.AND.(IR(KLTP).GE.1.OR.IPHGAS(KLTP).EQ.2))
     &      .OR.(IJTP.EQ.KLTP.AND.(IR(KLTP).GE.2.OR.IPHGAS(KLTP).EQ.2))
     &      ) THEN
                 NDXTP = NDXTP + 1
                 ITP(NDXTP) = IJTP
                 JTP(NDXTP) = IJTP
                 KTP(NDXTP) = KLTP
                 LTP(NDXTP) = KLTP
              END IF
            END DO
          END IF
        END DO
*. Strings differ by single excitation
      ELSE IF( NCREA.EQ.1.AND.NANNI.EQ.1) THEN
*. diagonal excitation plus creation in ICREA1 and annihilation in IANNI1
        DO IDIA = 1, NOBTP
          IF((IDIA.NE.IANNI1.AND.(IR(IDIA).GE.1.OR.IPHGAS(IDIA).EQ.2))
     &   .OR.(IDIA.EQ.IANNI1.AND.(IR(IDIA).GE.2.OR.IPHGAS(IDIA).EQ.2))
     &   )THEN
             NDXTP = NDXTP + 1
             ITP(NDXTP) = MAX(ICREA1,IDIA)
             KTP(NDXTP) = MIN(ICREA1,IDIA)
             JTP(NDXTP) = MAX(IANNI1,IDIA)
             LTP(NDXTP) = MIN(IANNI1,IDIA)
          END IF
        END DO
      ELSE IF(NCREA.EQ.2.AND.NANNI.EQ.2) THEN
*. Strings differ by double excitation
        NDXTP = 1
        ITP(1) = ICREA2
        KTP(1) = ICREA1
        JTP(1) = IANNI2
        LTP(1) = IANNI1
      END IF
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,I4)')
     &  ' Number of connecting double excitations ', NDXTP
        IF(NDXTP.NE.0) THEN
          WRITE(6,*) '  ITYP KTYP LTYP JTYP '
          WRITE(6,*) '  ===================='
          DO  IDX = 1,NDXTP
            WRITE(6,'(1X,4I5)')ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE EN_FROM_DENS(E,I12)
*
* Obtain energy from densities and integrals
*
*
* E = SUM(i,j) H(i,j) * RHO1(i,j)
*          + 0.5*SUM(i,j,K,L) (I J K L ) * RHO2( I J K L )
*
* Jeppe Olsen, Early 1997
*              Sept. 98    : I12 added
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLINT, KLDEN
      ! for addressing of WORK
      REAL*8 INPROD
*. Input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cgas.inc"
*
      COMMON/CECORE/ECORE
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
C?    WRITE(6,*) ' Energy under construction '
C?    WRITE(6,*) ' =========================='
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'EN_FRM')
      E1 = 0.0D0
      E2 = 0.0D0
*. Largest set of orbitals with given symmetry and type
CTF
* Using local MXTSOB_L (MXTSOB is now a common parameter!)
      MXTSOB_L = 0
      DO ISM = 1, NSMOB
      DO IGAS = 1, NGAS
        MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IGAS,ISM))
      END DO
      END DO
*. Allocate scratch space for 2-electron integrals and
*. two-electron densities
      MX4IBLK = MXTSOB_L ** 4
      CALL MEMMAN(KLINT,MX4IBLK,'ADDL  ',2,'KLINT')
      CALL MEMMAN(KLDEN,MX4IBLK,'ADDL  ',2,'KLDEN')
      ONE = 1.0D0
      DO ISM = 1, NSMOB
       DO JSM = 1, NSMOB
        CALL  SYMCOM(3,1,ISM,JSM,IJSM)
        DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
           NI = NOBPTS(IGAS,ISM)
           NJ = NOBPTS(JGAS,JSM)
           II = IOBPTS(IGAS,ISM)
           IJ = IOBPTS(JGAS,JSM)
           IF(ISM.EQ.JSM) THEN
*
* One-electron part
* =================
*
*. blocks of one-electron integrals and one-electron density

             CALL GETD1(WORK(KLDEN),ISM,IGAS,ISM,JGAS)
             CALL GETH1(WORK(KLINT),ISM,IGAS,ISM,JGAS)
C?           WRITE(6,*) ' Block of 1e integrals ISM,IGAS,JGAS',
C?   &                  ISM,IGAS,JGAS
C?           CALL WRTMT_LU(WORK(KLINT),NI,NJ,NI,NJ)
C?           WRITE(6,*) ' Block of 1e density ISM,IGAS,JGAS',
C?   &                  ISM,IGAS,JGAS
C?           CALL WRTMT_LU(WORK(KLDEN),NI,NJ,NI,NJ)
             E1 = E1 + INPROD(WORK(KLDEN),WORK(KLINT),NI*NJ)
           END IF
*
* Two-electron part
* =================
*
           IF(I12.EQ.2) THEN
           DO KSM = 1, NSMOB
*. Obtain LSM
             CALL  SYMCOM(3,1,IJSM,KSM,IJKSM)
             IJKLSM = 1
             CALL  SYMCOM(2,1,IJKSM,LSM,IJKLSM)
C?           WRITE(6,*) ' IJSM IJKSM LSM ',IJSM,IJKSM,IJKLSM
*
             DO KGAS = 1, NGAS
             DO LGAS = 1, NGAS
                NK = NOBPTS(KGAS,KSM)
                NL = NOBPTS(LGAS,LSM)
*. Blocks of density matrix and integrals
                IXCHNG = 0
                ICOUL  = 1
                CALL LGETINT(WORK(KLINT),
     &               IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &               IXCHNG,0,0,ICOUL)
                CALL GETD2 (WORK(KLDEN),
     &               ISM,IGAS,JSM,JGAS,KSM,KGAS,LSM,LGAS,ICOUL)
C?              write(6,*) ' Ism Jsm Ksm Lsm' , Ism,Jsm,Ksm,Lsm
C?              write(6,*)
C?   &          ' Igas Jgas Kgas Lgas' , Igas,Jgas,Kgas,Lgas
C?              WRITE(6,*) ' Integral block'
C?              CALL WRTMT_LU(WORK(KLINT),NI*NJ,NK*NL,NI*NJ,NK*NL)
C?              WRITE(6,*) ' Density block '
C?              CALL WRTMT_LU(WORK(KLDEN),NI*NJ,NK*NL,NI*NJ,NK*NL)
                NIJKL = NI*NJ*NK*NL
                E2 = E2 + 0.5D0*INPROD(WORK(KLDEN),WORK(KLINT),NIJKL)
C?              write(6,*) ' Updated 2e-energy ', E2
             END DO
             END DO
           END DO
           END IF
*
          END DO
         END DO
       END DO
      END DO
*
      E = ECORE + E1 + E2
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Output from EN_FROM_DEN'
        WRITE(6,*)
        WRITE(6,*) ' One-electron energy ', E1
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Two-electron energy ', E2
        END IF
        WRITE(6,*)
        WRITE(6,*) ' Total energy : ', E
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'EN_FRM')
*
      RETURN
      END
***********************************************************************
      SUBROUTINE EXCCLS2(NCLS,IACTIN,IACTOUT,IEXC,
     &                   IBASSPC_MX,IBASSPC)
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* If IBASSPC_MX .ne. 0, atmost basespaces belonging to this
*                       space is included
*
* Master routine
*
* Jeppe Olsen, Jan. 1999 - ved siden af ditte, paa MAS efter
*              hendes rygoperation
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IACTIN(NCLS)
*. Output
      DIMENSION IACTOUT(NCLS)
*. From the common blocks
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cgas.inc"
*
      CALL EXCCLS2_S(NGAS,WORK(KLOCCLS),NCLS,IACTIN,IACTOUT,IEXC,
     &               IBASSPC_MX,IBASSPC)
*
      RETURN
      END
***********************************************************************
      SUBROUTINE EXCCLS2_S(NGAS,ICLS,NCLS,IACTIN,IACTOUT,IEXC,
     &               IBASSPC_MX,IBASSPC)
*
* A set of classes ICLS are given with the active
* classes indicated by nonvanishing elements in IACTIN.
*
* Obtain classes that can be obtained by atmost IEXC excitations
*
* If IBASSPC_MX .ne. 0, atmost basespaces belonging to this
*                       space is included
* Slave routine
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS(NGAS,NCLS)
      DIMENSION IACTIN(NCLS)
      DIMENSION IBASSPC(*)
*. Output
      DIMENSION IACTOUT(NCLS)
*
C     write(6,*) ' ICLS in EXCCLS '
C     call iwrtma(icls,ngas,ncls,ngas,ncls)

      IZERO = 0
      CALL ISETVC(IACTOUT,IZERO,NCLS)
      DO ICLSIN = 1, NCLS
        IF(IACTIN(ICLSIN).NE.0) THEN
        DO ICLSOUT = 1, NCLS
*. Number of anihilations and creations required to connect classes
          NANNI = 0
          NCREA = 0
          DO IGAS = 1, NGAS
*
            IDIF = ICLS(IGAS,ICLSOUT)-ICLS(IGAS,ICLSIN)
            IF(IDIF .GT. 0 ) THEN
              NCREA = NCREA + IDIF
            ELSE IF (IDIF .LT. 0 ) THEN
              NANNI = NANNI - IDIF
            END IF
          END DO
*
          IF(NCREA.LE.IEXC) THEN
            IF(IBASSPC_MX.EQ.0.OR.IBASSPC(ICLSOUT).LE.IBASSPC_MX) THEN
              IACTOUT(ICLSOUT) = 1
            END IF
          END IF
*
        END DO
        END IF
      END DO
*
      NACTOUT = 0
      DO ICLSOUT = 1, NCLS
        IF(IACTOUT(ICLSOUT).GT.0) NACTOUT = NACTOUT + 1
      END DO
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Output from EXCCLS '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' Number of active output classes ',NACTOUT
         WRITE(6,*)
      END IF
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Active output classes '
         CALL IWRTMA(IACTOUT,NCLS,1,NCLS,1)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GASSPC
*
*
* Divide orbital spaces ( I_IAD ) into
*
*  Inactive spaces : Orbitals that are doubly occupied in all CI spaces
*  Active orbitals : Orbitals that have variable occ in atleast some spaces.
*  Secondary spaces: Orbitals that are unoccupied in all spaces
*
* Division based upon occupation in Compound CI spaces IGSOCC
*
* Jeppe Olsen, Summer of 98 ( not much of a summer !)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "strinp.inc"
#include "orbinp.inc"
*
      NEL_REF = NELEC(1) + NELEC(2)
*
      DO IGAS = 1, NGAS
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS)
       END IF
*
       IF(IGSOCC(IGAS,1) .EQ. NEL_MAX  .AND.
     &    IGSOCC(IGAS,2) .EQ. NEL_MAX       ) THEN
*. Inactive  space
          I_IAD(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCC(IGAS-1,1) .EQ. NEL_REF ) THEN
*. Delete space
          I_IAD(IGAS) = 3
       ELSE
*. Active space
          I_IAD(IGAS) = 2
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Division of orbitals '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' Inactive = 1, Active = 2, Delete = 3 '
        WRITE(6,*)
        CALL IWRTMA(I_IAD,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
***********************************************************************
*  Automatic splitting of GAS spaces if too many
*  orbitals per GAS:
*
      subroutine gassplit(NSHPGS,IGSSPLIT,MXPNGAS,
     &                    MXTSOB,NGASOUT,IPRNT)
*
      implicit real*8 (A-H,O-Z)
*
      dimension NSHPGS(MXPNGAS)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.5) then
        write(6,*) '++++++++++++++++++++++'
        write(6,*) ' GASSPLIT in action   '
        write(6,*) '++++++++++++++++++++++'
        write(6,*) 'Splitting space ...  with ',
     &             IGSSPLIT,NSHPGS(IGSSPLIT)
      end if
*
      ITER = 1
      NGASCT = 0
*
 10   ISAVE = NSHPGS(IGSSPLIT+ITER-1)
*
      if ((IGSSPLIT+ITER).gt.NGASOUT) then
        NGASOUT = NGASOUT + 1
        if (NTEST.ge.1) write(6,*) ' NGAS increased by 1.'
      end if
*
      NSHPGS(IGSSPLIT+ITER-1) = MXTSOB
      NSHPGS(IGSSPLIT+ITER) = ISAVE - MXTSOB
*
      if (NSHPGS(IGSSPLIT+ITER).gt.MXTSOB) then
        ITER = ITER + 1
        goto 10
      end if
*
      if (NTEST.ge.3) then
        write(6,*) 'New NSHPGS partition:'
        call iwrtma(NSHPGS,1,NGASOUT,1,MXPNGAS)
      end if
*
      return
      end
***********************************************************************
      SUBROUTINE GET_3BLKS(KVEC1,KVEC2,KC2)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* Jeppe Olsen, Jan 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KVEC1, KVEC2, KC2, KLCIOIO, KLCBLTP, KSVST
      ! for addressing of WORK
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
*
#include "cgas.inc"
#include "lucinp.inc"
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_3B')
*
      NTEST = 1000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*
      ISPC = ISSPC
      ISM  = ISSM
      CALL IAIBCM(ISPC,WORK(KLCIOIO))
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*. Largest number of strings of given symmetry and type
      MAXA = MXNSTR
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESCPH(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &              IPHGAS,NHLFSPGP,MNHL,IADVICE)
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'GET_3BL: MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'GET_3BL: MXADKBLK ', MXADKBLK
        WRITE(6,*) ' Space for resolution matrices ',LSCR2
      END IF
*
*. For Passive/active division
      LSCR12 = MAX(LBLOCK,4*LSCR2)
C     LSCR12 = MAX(LBLOCK,2*LSCR2)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_3B')
*
      CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      CALL MEMMAN(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
      CALL MEMMAN(KC2,LSCR12,'ADDS  ',2,'C2    ')
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_SPGP_INF(ISPGP,ITP,IGRP)
*
* Obtain groups defining supergroup ISPGP,ITP
*
* Jeppe Olsen, November 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Output
      DIMENSION IGRP(*)
*
      NTEST = 000
*. Absolute group number
      ISPGPABS = ISPGP + IBSPGPFTP(ITP) -1
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),IGRP,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' GET_SPGP_INF : ISPGP ITP ISPGPABS',
     &              ISPGP, ITP, ISPGPABS
        WRITE(6,*) ' Groups of supergroups'
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GTSPGP(IEL,ISPGP,IWAY)
*
*
* Relation between number of electrons in AS1, AS2 ... and
* supergoup number
*
* IWAY = 1 :
* Get ISPGP : Supergroup of strings
*             with IEL(*)  electrons in each AS
* IWAY = 2 :
* GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP
*
*
* Jeppe Olsen, Another lonely night in Lund
*               GAS version July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"
C     COMMON/GASSTR/MNGSOC(MXPNGAS),MXGSOC(MXPNGAS),NGPSTR(MXPNGAS),
C    &              IBGPSTR(MXPNGAS),NELFGP(MXPSTT),IGSFGP(MXPSTT),
C    &              NSTFGP(MXPSTT),MNELFGP(MXPNGAS),MXELFGP(MXPNGAS),
C    &              NELFTP(MXPSTT),NSPGPFTP(MXPSTT),IBSPGPFTP(MXPSTT),
C    &              ISPGPFTP(MXPNGAS,MXPSTT),NELFSPGP(MXPNGAS,MXPSTT),
C    &              NGRP,NSTTP,MXNSTR,NTSPGP
*
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,
C    &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC
*. input(IWAY = 2 ), output (IWAY = 1 )
      INTEGER IEL(*)
*
      IF(IWAY.EQ.1) THEN
*. Occupation => Number
        ISPGP = -1
        DO JSPGP = 1, NTSPGP
          IF(ISPGP.EQ.-1) THEN
            IEQUAL = 1
            DO IGAS = 1, NGAS
              IF(NELFSPGP(IGAS,JSPGP).NE.IEL(IGAS)) IEQUAL= 0
            END DO
            IF(IEQUAL.EQ.1) ISPGP = JSPGP
          END IF
        END DO
      ELSE IF (IWAY .EQ. 2 ) THEN
*. Number => Occupation
        DO IGAS = 1, NGAS
         IEL(IGAS) = NELFSPGP(IGAS,ISPGP)
        END DO
      END IF
*
      NTEST  = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from GTSPGP '
        WRITE(6,*)
     &   ' IWAY ISPGP IEL ', IWAY,ISPGP,(IEL(IGAS),IGAS = 1, NGAS)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE IAIBCM(ICISPC,IAIB)
*
* obtain allowed combinbation of alpha- and beta- supergroups
* for CI space ICISPC
*
* Master for IAIBCM_GAS
*
*      Jeppe Olsen, august 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "stinf.inc"
#include "strinp.inc"
*. Output
      INTEGER IAIB(*)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
C?    write(6,*) ' IAIB ::::::'
C?    write(6,*) ' LCMBSPC, ICISPC, ICMBSPC '
C?    WRITE(6,*) ICISPC,  LCMBSPC(ICISPC)
C?    WRITE(6,*) (ICMBSPC(II,ICISPC),II=1, LCMBSPC(ICISPC))

      CALL IAIBCM_GAS(LCMBSPC(ICISPC),ICMBSPC(1,ICISPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,IAIB,IPRDIA)
*
      RETURN
      END
***********************************************************************
      FUNCTION IBASSPC_FOR_CLS(ICLS)
*
*. Obtain base space for occupation class ICLS
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
*. Specific input
      INTEGER ICLS(NGAS)
C     COMMON/CGAS/IDOGAS,NGAS,NGSSH(MXPIRR,MXPNGAS),
C    &            NGSOB(MXPOBS,MXPNGAS),
C    &            NGSOBT(MXPNGAS),IGSOCC(MXPNGAS,2),IGSINA,IGSDEL,
C    &            IGSOCCX(MXPNGAS,2,MXPICI),NCISPC,
C    &            NCMBSPC, LCMBSPC(MXPICI),ICMBSPC(MXPSTT,MXPICI),
C    &            NMXOCCLS,IPHGAS(MXPNGAS),IHPVGAS(MXPNGAS)
*
      IBASE = 0
      DO ISPC = 1, NCMBSPC
        DO JJSPC = 1, LCMBSPC(ISPC)
          JSPC = ICMBSPC(JJSPC,ISPC)
*. Test for occupation constraints in CI space JSPC
          I_AM_OKAY = 1
          DO IGAS = 1, NGAS
            IF(IGAS.EQ.1) THEN
              NEL = ICLS(IGAS)
            ELSE
              NEL = NEL + ICLS(IGAS)
            END IF
*
            IF(NEL.LT.IGSOCCX(IGAS,1,JSPC).OR.
     &         NEL.GT.IGSOCCX(IGAS,2,JSPC)    ) THEN
                I_AM_OKAY = 0
            END IF
          END DO
*         ^ End of loop over gasspaces for given cispace
*
          IF(I_AM_OKAY.EQ.1.AND.IBASE.EQ.0) THEN
            IBASE = ISPC
          END IF
*
        END DO
*       ^ End of loop over cisspaces for given combination space
      END DO
*     ^ End of loop over combinations apaces
*
      IBASSPC_FOR_CLS = IBASE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation class and its basespace '
        CALL IWRTMA(ICLS,1,NGAS,1,NGAS)
        WRITE(6,*) IBASE
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE INTDIM(IPRNT)
*
* Number of integrals and storage mode
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "crun.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
* =======
*. Output
* =======
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*.1 : Number of one-electron integrals
      NINT1 =  NSXFSM(NSMOB,MXPOBS,NTOOBS,NTOOBS,ITSSX,ADSXA,1,IPRNT)
*.2 : Number of two-electron integrals
C     IF(INCORE.EQ.1.OR.EXTSPC.EQ.0) THEN
        IF(PNTGRP.EQ.1) THEN
*. Full eightfold symmetry can be used
          I12S = 1
          I34S = 1
          I1234S = 1
        ELSE
*. Only symmetry between 12 and 34
          I12S = 0
          I34S = 0
          I1234S = 1
        END IF
        NINT2 = NDXFSM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,
     &                  NTOOBS,ITSDX,ADSXA,SXDXSX,I12S,I34S,I1234S,
     &                  IPRNT)
C     END IF
*. Number of symmetry blocks of one- and two-electron integrals
      NBINT1 = NSMOB
      NBINT2 = NSMOB ** 3
      RETURN
      END
***********************************************************************
      SUBROUTINE INTPNT(IPNT1,ISL1,IPNT2,ISL2)
*
* Pointers to symmetry blocks of integrals
* IPNT1 : Pointer to given one-electron block, total symmetric
* ISL1  : Symmetry of last index for given first index, 1 e-
* IPNT2 : Pointer to given two-electron block
* ISL1  : Symmetry of last index for given first index, 1 e-
*
* In addition pointers to one-electron integrals with general
* symmetry is generated in WORK(KPGINT1(ISM))
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
C     INCLUDE
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*.CINTFO
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
* =======
*. Output
* =======
*
      INTEGER IPNT1(NSMOB),ISL1(NSMOB)
      INTEGER IPNT2(NSMOB,NSMOB,NSMOB),ISL2(NSMOB,NSMOB,NSMOB)
*.0 : Pointers to one-integrals, all symmetries, Lower half matrices
      DO ISM = 1, NSMOB
        CALL PNT2DM(1,NSMOB,NSMSX,ADSXA,NTOOBS,NTOOBS,
     &       ISM  ,ISL1,WORK(KPGINT1(ISM)),MXPOBS)
      END DO
*.0.5 : Pointers to one-electron integrals, all symmetries, complete form
      DO ISM = 1, NSMOB
        CALL PNT2DM(0,NSMOB,NSMSX,ADSXA,NTOOBS,NTOOBS,
     &       ISM  ,ISL1,WORK(KPGINT1A(ISM)),MXPOBS)
      END DO
*.1 : Number of one-electron integrals
      CALL PNT2DM(1,NSMOB,NSMSX,ADSXA,NTOOBS,NTOOBS,
     &            ITSSX,ISL1,IPNT1,MXPOBS)
*.2 : two-electron integrals
      CALL PNT4DM(NSMOB,NSMSX,MXPOBS,NTOOBS,NTOOBS,NTOOBS,NTOOBS,
     &            ITSDX,ADSXA,SXDXSX,I12S,I34S,I1234S,IPNT2,ISL2,
     &            ADASX)
C?    write(6,*) ' Memory check INTPNT 2 '
C?    CALL MEMCHK
      RETURN
      END
***********************************************************************
      SUBROUTINE LCISPC(IPRNT)
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*
* GAS VERSION
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLBLTP, KLCVST, KLIOIO
      ! for addressing of WORK
*
* ===================
*.Input common blocks
* ===================
*
*./LUCINP/
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "strbas.inc"
*./CSM/ : NSMST is used
#include "csm.inc"
#include "stinf.inc"
*./CSTATE/ : IDC is used
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.inc"
#include "parluci.h"
      CALL QENTER('LCISP')
*
*
*. Number of spaces
      NICISP = NCMBSPC
C?    write(luwrt,*) ' LCISPC : NICISP ', NICISP
*. Type of alpha- and beta strings
      IATP = 1
      IBTP = 2
*
      NOCTPA =  NOCTYP(IATP)
      NOCTPB =  NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*.Local memory
      CALL MEMMAN(KLBLTP,NSMST,'ADDL  ',2,'KLBLTP')
      IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
     &CALL MEMMAN(KLCVST,NSMST,'ADDL  ',2,'KLCVST')
      CALL MEMMAN(KLIOIO,NOCTPA*NOCTPB,   'ADDL  ',2,'KLIOIO')
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &CALL SIGVST(WORK(KLCVST),NSMST)

*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
      CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*
      MXSB = 0
      MXSOOB = 0
      DO 100 ICI = 1, NICISP
*. allowed combination of types
      CALL IAIBCM_GAS(LCMBSPC(ICI),ICMBSPC(1,ICI),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KLIOIO),IPRNT)
*
      DO  50 ISYM = 1, NSMCI
          if (IDC.eq.3.or.IDC.eq.4) then
            CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP),
     &                 WORK(KLCVST))
          else
            call zbltp_idc(ISMOST(1,ISYM),NSMST,IDC,WORK(KLBLTP))
          end if
          CALL NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),NGAS,ISYM,
     &                NSMST,NOCTPA,NOCTPB,WORK(KNSTSO(IATP)),
     &                WORK(KNSTSO(IBTP)),
     &                ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NELFGP,
     &                NCOMB,XNCOMB,MXS,MXSOO,WORK(KLBLTP),NTTSBL,
     &                LCOL,WORK(KLIOIO))

          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB           = MAX(MXSOOB,MXSOO)
          MXSB             = MAX(MXSB,MXS)
          NBLKIC(ISYM,ICI) = NTTSBL
          LCOLIC(ISYM,ICI) = LCOL
   50 CONTINUE
  100 CONTINUE
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      WRITE(luwrt,'(/a/a)')
     & '  Number of internal combinations per symmetry',
     & '  ============================================'
*
      DO ICI = 1, NCMBSPC
         WRITE(luwrt,*) ' CI space ', ICI
         WRITE(luwrt,'(4F20.0)')     (XISPSM(II,ICI),II=1,NSMCI)
         WRITE(luwrt,'(1P,4D20.10)') (XISPSM(II,ICI),II=1,NSMCI)
      END DO
      if(NTEST.ge.1) then
        WRITE(luwrt,*) ' Largest symmetry block           ',MXSB
        WRITE(luwrt,*) ' Largest Symmetry-type-type block ',MXSOOB
      end if
*
      IF(NTEST.GE.5) THEN
        WRITE(luwrt,'(/A/A)')
     &  '  Number of TTS subblocks per CI expansion',
     &  '  ========================================'
*
        DO ICI = 1,  NCMBSPC
          WRITE(luwrt,*) ' Internal CI space ', ICI
          CALL IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      END IF

*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
      if(NTEST.ge.1)then
        WRITE(luwrt,*) ' Largest number of blocks in a CI expansion',
     &                   MXNTTS
      end if
*
      IF(NTEST.GE.5) THEN
        WRITE(luwrt,'(/a/a)')
     & '  Number of columns per CI expansion',
     & '  =================================='
*
        DO  ICI = 1,  NCMBSPC
            WRITE(luwrt,*) ' Internal CI space ', ICI
            CALL IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      END IF
*
      CALL QEXIT('LCISP')
*
      END
***********************************************************************
*CADDB
      SUBROUTINE LCNHCN(LSCR)
*
* Amount of scratch Needed in the CNHCNM routine
*
* Jeppe Olsen, September 1993
*
* Amount of Memory required:2*NACTEL + MXCSFC**2 +
*                           6*MXDTFC+MXDTFC**2+MXCSFC*MXDTFC+
*                           MAX(MXDTFC*NACTEL+2*NACTEL,4*NACOB+2*NACTEL)
*
* Where NACTEL : Number of active electrons
*       NACOB  : Number of active orbitals
*       MXCSFC : Max number of CSF's for given COnfiguration
*       MXDTFC : Max number of Combs for given configuration
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "spinfo.inc"
*. NACTEL is obtained from lucinp
#include "lucinp.inc"
*. NACOB is obtained from orbinp
#include "orbinp.inc"
*./SPINFO/
C     COMMON/SPINFO/MULTSP,MS2P,
C    &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),
C    &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),
C    &              NCNASM(MXPCSM)
*
*. MXCSFC, MXSDFC
      MXCSFC = 0
      MXDTFC = 0
      DO 100 ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCPCNT(ITYP))
        MXDTFC = MAX(MXDTFC,NDPCNT(ITYP))
  100 CONTINUE
*
*
      LSCR  = 2*NACTEL + MXCSFC**2 +
     &        6*MXDTFC+MXDTFC**2+MXCSFC*MXDTFC+
     &        MAX(MXDTFC*NACTEL+2*NACTEL,4*NACOB+2*NACTEL)
*
C?    WRITE(6,*) ' LCNHCN : MXCSFC MXDTFC ',MXCSFC,MXDTFC
C?    WRITE(6,*) ' LCNHCN : LSCR ', LSCR
*
      RETURN
      END
***********************************************************************
      SUBROUTINE MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDIST)
*
* A combination of NIGRP groups are given (IGRP)
*. Find MIN and MAX for symmetry in each group
*
* Jeppe Olsen, September 1997
*              April 1998     From  MINMAX_SM_GP
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Include blocks
#include "mxpdim.inc"
#include "strbas.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "wrkspc.inc"
*. Input
      DIMENSION IGRP(NIGRP)
*.Output
      DIMENSION MNVAL(NIGRP),MXVAL(NIGRP)
*. Local scratch
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*
      NTEST = 0000
      IF(NTEST.GE.100) WRITE(6,*) ' >> Entering MINMAX_... <<'
*
      DO JGRP = 1, NIGRP
        MNVAL(JGRP) = MINMAX_SM_GP(1,IGRP(JGRP))
        MXVAL(JGRP) = MINMAX_SM_GP(2,IGRP(JGRP))
      END DO

*. Number of strings per sym and group
C     DO JGRP = 1, NIGRP
C       CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(JGRP)-1)*NSMST+1,
C    &               NSMST,LSMGP(1,JGRP))
C     END DO
C     IF(NTEST.GE.1000) THEN
C       WRITE(6,*) ' LSMGP '
C       CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
C     END IF
C. Max and min sym in each group
C     DO JGRP = 1, NIGRP
*
C       IMAX = 1
C       DO ISM=1, NSMST
C         IF(LSMGP(ISM,JGRP).GT.0) IMAX = ISM
C       END DO
C       MXVAL(JGRP) = IMAX
*
C       IMIN = NSMST
C       DO ISM = NSMST,1,-1
C        IF(LSMGP(ISM,JGRP).GT.0) IMIN = ISM
C       END DO
C       MNVAL(JGRP) = IMIN
C     END DO
*. Total number of symmetry distributions
      NDIST = 1
      DO JGRP = 1, NIGRP
        NDIST = NDIST*(MXVAL(JGRP)-MNVAL(JGRP)+1)
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Group combination : '
        WRITE(6,'(5X,10I3)') (IGRP(JGRP),JGRP=1, NIGRP)
        WRITE(6,*)
        WRITE(6,*) ' Group Minsym Maxsym'
        WRITE(6,*) ' ==================='
        DO JGRP = 1, NIGRP
          WRITE(6,'(3I6)') IGRP(JGRP),MNVAL(JGRP),MXVAL(JGRP)
        END DO
        WRITE(6,*)
        WRITE(6,*) ' Total number of distributions', NDIST
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving MINMAX_... <<'
*
      RETURN
      END
***********************************************************************
      SUBROUTINE MLSM(IML,IPARI,ISM,TYPE,IWAY)
*
* Transfer between ML,IPARI notation and compound notation ISM
*
* IWAY = 1 : IML,IPARI => Compound
* IWAY = 2 : IML,IPARI <= Compound
*
* TYPE : 'SX','OB','ST','DX','CI'
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
      CHARACTER*2 TYPE
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONAB/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*./CSM/
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
*.(Tired of warnings from 3090 Compiler so )
      NML = 0
      MXML= 0
      MNML= 0
*             )
      IF(TYPE.EQ.'OB') THEN
        NML = NMLOB
        MXML = MXMLOB
        MNML = MNMLOB
      ELSE IF(TYPE.EQ.'SX') THEN
        NML = NMLSX
        MXML = MXMLSX
        MNML = MNMLSX
      ELSE IF(TYPE.EQ.'DX') THEN
        NML = NMLDX
        MXML = MXMLDX
        MNML = MNMLDX
      ELSE IF(TYPE.EQ.'ST') THEN
        NML = NMLST
        MXML = MXMLST
        MNML = MNMLST
      ELSE IF(TYPE.EQ.'CI') THEN
        NML = NMLCI
        MXML = MXMLCI
        MNML = MNMLCI
      END IF
*
      IF(IWAY.EQ.1) THEN
C        ISM = (IPARI-1)*NML + MNML - 1
         ISM = (IPARI-1)*NML + IML - MNML + 1
      ELSE IF(IWAY.EQ.2) THEN
        IF(ISM.GT.NML) THEN
          IPARI = 2
          IML = ISM - NML + MNML - 1
        ELSE
          IPARI = 1
          IML = ISM       + MNML - 1
        END IF
      ELSE
        WRITE(6,*) ' Error in MLSM , IWAY = ' ,IWAY
        WRITE(6,*) ' MLSM stop !!! '
        Call Abend1( 20 )
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,A)') ' MLSM speaking ,type= ',TYPE
        WRITE(6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE LNATORB(RHO1,NSMOB,NTOPSM,NACPSM,NINPSM,
     &                   ISTOB,XNAT,RHO1SM,OCCNUM,
     &                   NACOB,SCR,IPRDEN)
*
* Obtain natural orbitals in symmetry blocks
*
* Jeppe Olsen, June 1994
*              Modification, Oct 94
*              Last modification, Feb. 1998 (reorder deg eigenvalues)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION RHO1(NACOB,NACOB)
      DIMENSION ISTOB(*)
      DIMENSION NTOPSM(NSMOB), NACPSM(NSMOB),NINPSM(NSMOB)
*. Output
      DIMENSION RHO1SM(*),OCCNUM(*),XNAT(*)
*. Scratch ( Largest symmetry block )
      DIMENSION SCR(*)
#include "parluci.h"
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRDEN)
*. IOBOFF : Offset for active orbitals in symmetry order
      DO ISMOB = 1, NSMOB
        IF(ISMOB.EQ.1) THEN
          IOBOFF = NINPSM(1)+1
          IMTOFF = 1
        ELSE
          IOBOFF =
     &    IOBOFF + NTOPSM(ISMOB-1)-NINPSM(ISMOB-1)+NINPSM(ISMOB)
          IMTOFF = IMTOFF + NACPSM(ISMOB-1)**2
        END IF
        LOB = NACPSM(ISMOB)
*
*. Extract symmetry block of density matrix
*
        DO IOB = IOBOFF,IOBOFF + LOB-1
           DO JOB = IOBOFF,IOBOFF + LOB-1
*. Corresponding type indeces
             IOBP = ISTOB(IOB)
             JOBP = ISTOB(JOB)
             RHO1SM(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1)
     &     = RHO1(IOBP,JOBP)
           END DO
        END DO
*
        IF (NTEST.GE.2) THEN
          WRITE(LUWRT,*)
          WRITE(LUWRT,*) ' Density matrix for symmetry  = ', ISMOB
          WRITE(LUWRT,*) ' ======================================='
          WRITE(LUWRT,*)
          CALL WRTMATMN(RHO1SM(IMTOFF),LOB,LOB,LOB,LOB,LUWRT)
        END IF
*. Pack and diagonalize
        SIGN = 1.0D0
        CALL TRIPAK_LUCI(RHO1SM(IMTOFF),SCR,1,LOB,LOB,SIGN)
        ONEM = -1.0D0
*. scale with -1 to get highest occupation numbers as first eigenvectors
        CALL SCALVE(SCR,ONEM,LOB*(LOB+1)/2)
        CALL EIGEN_LUCI(SCR,XNAT(IMTOFF),LOB,0,1)
*
        DO  I = 1, LOB
          OCCNUM(IOBOFF-1+I) = - SCR(I*(I+1)/2)
        END DO
*. Order the degenerate eigenvalues so diagonal terms are maximized
        TESTY = 1.0D-11
        DO IOB = 2, LOB
          IF(ABS(OCCNUM(IOBOFF-1+IOB)-OCCNUM(IOBOFF-2+IOB))
     &       .LE.TESTY) THEN
            XII   = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB  ))
            XI1I1 = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB-1))
            XII1  = ABS(XNAT(IMTOFF-1+(IOB-1-1)*LOB+IOB  ))
            XI1I  = ABS(XNAT(IMTOFF-1+(IOB-1)  *LOB+IOB-1))
*
            IF( XI1I.GT.XII.AND.XII1.GT.XI1I1 ) THEN
*. interchange orbital IOB and IOB -1
              CALL SWAPVE(XNAT(IMTOFF+(IOB-1)*LOB),
     &                    XNAT(IMTOFF+(IOB-1-1)*LOB),LOB)
              SS = OCCNUM(IOBOFF-1+IOB-1)
              OCCNUM(IOBOFF-1+IOB-1) = OCCNUM(IOBOFF-1+IOB)
              OCCNUM(IOBOFF-1+IOB)   = SS
              write(LUWRT,*) ' Orbitals interchanged ',
     &        IOBOFF-1+IOB,IOBOFF-2+IOB
            END IF
          END IF
        END DO
*
        IF (NTEST.GE.1) THEN
          WRITE(LUWRT,*) ' '
          WRITE(LUWRT,'(2X ,A,I3)')
     &    ' Natural occupation numbers for symmetry = ', ISMOB
          WRITE(LUWRT,'(2X ,A)')
     &    ' __________________________________________________'
          OCCSUM = DSUM(LOB,OCCNUM(IOBOFF),1)
          WRITE(LUWRT,9002) (OCCNUM(IOBOFF+J-1),J=1,LOB)
          WRITE(LUWRT,9003) OCCSUM
csk       CALL WRTMT_LU(OCCNUM(IOBOFF),1,LOB,1,LOB)
          IF (NTEST.GE.2) THEN
            WRITE(LUWRT,*)
            WRITE(LUWRT,*) ' Corresponding Eigenvectors '
            WRITE(LUWRT,*)
            CALL WRTMATMN(XNAT(IMTOFF),LOB,LOB,LOB,LOB,LUWRT)
          END IF
        END IF
      END DO
*. ( End of loop over orbital symmetries )
*
 9002 FORMAT(/,(5F14.9))
 9003 FORMAT(/'   Sum =',T15,F14.9)
C
      RETURN
      END
***********************************************************************
      SUBROUTINE NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
*
* an input  supergroup is given.
* apply an string of elementary operators to this supergroup and
* obtain supergroup mumber of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95
*
* ------
* Input
* ------
*
* INSPGP  : input super group ( given occupation in each AS )
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER ITPOP(*),IACOP(*)
*. Number of active spaces  (NGAS )
#include "mxpdim.inc"
#include "cgas.inc"
#include "strbas.inc"
#include "wrkspc.inc"
*. Local scratch
      DIMENSION IEL(MXPNGAS)
*. output
      INTEGER OUTSPGP
*
      INEW_OR_OLD = 1
      IF(INEW_OR_OLD.EQ.1) THEN
        CALL NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &       NGAS,WORK(KSPGPAN),WORK(KSPGPCR),OUTSPGP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
      ELSE IF(INEW_OR_OLD.EQ.2) THEN
*. Number of electrons in AS1,AS2, .... for input supergroup
        CALL GTSPGP(IEL,INSPGP,2)
*
        IDELTA = 0
        DO IOP = 1, NOP
*. Change in number of orbitals
          IF(IACOP(IOP).EQ.1) THEN
            IDELTA = IDELTA - 1
            IDELTAI = -1
          ELSE
            IDELTA = IDELTA + 1
            IDELTAI = +1
          END IF
          IEL(ITPOP(IOP)) = IEL(ITPOP(IOP)) + IDELTAI
        END DO
*. output supergroup
        CALL GTSPGP(IEL,OUTSPGP,1)
      END IF
*     /\ End  of NEW_OR_OLD switch
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP ,  OUTSPGP ', OUTSPGP
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE NEXT_SYM_DISTR(NGAS,MINVAL,MAXVAL,
     &           ISYM,ISYM_TOT,IFIRST,NONEW)
*
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
* Obtain next distribution of symmetries with given total
* Symmetry.
*
* Loop over first NGAS-1 spaces are performed, and the symmetry
* of the last space is then fixed by the required total sym
*
* Jeppe Olsen, Sept 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION MINVAL(NGAS),MAXVAL(NGAS)
*. Input and output
      DIMENSION ISYM(NGAS)
*
*. Symmetry of first NGAS -1 spaces
*
      IF(IFIRST.EQ.1) THEN
        DO IGAS = 1, NGAS-1
          ISYM(IGAS) = MINVAL(IGAS)
        END DO
        NONEW = 0
      END IF
 1001 CONTINUE
      IF(IFIRST.EQ.0) CALL NXTNUM3(ISYM,NGAS-1,MINVAL,MAXVAL,NONEW)
      IFIRST = 0
*
*. Symmetry of last space
*
      IF(NONEW.EQ.0) THEN
C       JSYM = 1
C       DO IGAS = 1, NGAS-1
C         CALL SYMCOM(3,0,JSYM,ISYM(IGAS),KSYM)
C         JSYM = KSYM
C       END DO
        JSYM = ISYMSTR(ISYM,NGAS-1)
        CALL SYMCOM(2,0,JSYM,ISYM(NGAS),ISYM_TOT)
*
        IF(MINVAL(NGAS).GT.ISYM(NGAS).OR.
     &     MAXVAL(NGAS).LT.ISYM(NGAS)    )GOTO 1001
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(NONEW.EQ.1) THEN
         WRITE(6,*) ' No new symmetry distributions '
        ELSE
         WRITE(6,*) ' Next symmetry distribution '
         CALL IWRTMA(ISYM,1,NGAS,1,NGAS)
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE NGASDT(IOCCMN,IOCCMX,NGAS,ITOTSM,
     &                  NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,
     &                  IAOCC,IBOCC,MXPNGAS,NELFGP,
     &                  NCOMB,XNCOMB,MXSB,MXSOOB,
     &                  IBLTP,NTTSBL,LCOL,IOCOC)
*
* Number of combinations with symmetry ITOTSM and
* occupation between IOCCMN and IOCCMX
*
* In view of the limited range of I*4, the number of dets
* is returned as integer and  real*8
*
* MXSB is largest UNPACKED symmetry block
* MXSOOB is largest UNPACKED symmetry-type-type block
* NTTSBL is number of TTS blocks in vector
* LCOL is the sum of the number of columns in each block
*
*
* Winter 94/95
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Allowed combinations of alpha and beta types
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Occupation constraints
      DIMENSION IOCCMN(NGAS),IOCCMX(NGAS)
*. Occupation of alpha and beta strings
      DIMENSION IAOCC(MXPNGAS,*),IBOCC(MXPNGAS,*)
*. Number of strings per supergroup and symmetry
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*),NELFGP(*)
*. block types
      DIMENSION IBLTP(*)
#include "parluci.h"
*
      CALL QENTER('NGASD')
      NTEST = 000
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' NGASDT speaking'
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' NGAS NOCTPA,NOCTPB ',NGAS,NOCTPA,NOCTPB
        WRITE(6,*) ' ITOTSM ', ITOTSM
        WRITE(6,*) ' Upper and lower occupation constraints'
        CALL IWRTMA(IOCCMN,1,NGAS,1,NGAS)
        CALL IWRTMA(IOCCMX,1,NGAS,1,NGAS)
        WRITE(6,*) ' IOCOC matrix '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        WRITE(6,*) ' Number of alpha and beta strings '
        CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
        CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
      END IF
*
      MXSB = 0
      MXSOOB = 0
      NCOMB = 0
      XNCOMB = 0.0D0
      NTTSBL = 0
      LCOL = 0
      DO 300 IASM = 1, NSMST
        IF(IBLTP(IASM).EQ.0) GOTO 300
        CALL SYMCOM(2,1,IASM,IBSM,ITOTSM)
C       write(6,*) ' Iasm ibsm ', iasm,ibsm
        LSB = 0
        IF(IBSM.NE.0) THEN
          IF(IBLTP(IASM).EQ.2) THEN
            ISYM = 1
          ELSE
            ISYM = 0
          END IF
          DO 200 IATP = 1, NOCTPA
           IF(ISYM.EQ.1) THEN
             MXBTP = IATP
           ELSE
             MXBTP = NOCTPB
           END IF
           LTSSA = NSSOA(IASM,IATP)
CSK           IF(LUCI_MYPROC.EQ.LUCI_MASTER) write(6,*)'LTSSA',LTSSA
           XLTSSA = FLOAT(LTSSA)
           DO 100 IBTP = 1, MXBTP
*
             IF(NTEST.GE.10) THEN
               WRITE(6,*) ' Alpha super group and beta super group'
               CALL IWRTMA(IAOCC(1,IATP),1,NGAS,1,NGAS)
               CALL IWRTMA(IBOCC(1,IBTP),1,NGAS,1,NGAS)
             END IF
*
             IF(IOCOC(IATP,IBTP).EQ.1) THEN
*. Size of unpacked block
               LTTSUP =  LTSSA*NSSOB(IBSM,IBTP)
CSK               IF(LUCI_MYPROC.EQ.LUCI_MASTER) write(6,*)'LTTSUP',LTTSUP
*. Size of packed block
               IF(ISYM.EQ.0.OR.IATP.NE.IBTP) THEN
                 LTTSBL = LTSSA*NSSOB(IBSM,IBTP)
CSK                 IF(LUCI_MYPROC.EQ.LUCI_MASTER) write(6,*)'LTTSBL',LTTSBL
                 XNCOMB = XNCOMB + XLTSSA*FLOAT(NSSOB(IBSM,IBTP))
               ELSE
                 LTTSBL = LTSSA *(LTSSA +1)/2
CSK                 IF(LUCI_MYPROC.EQ.LUCI_MASTER) write(6,*)'LTTSBL 2',LTTSBL
                 XNCOMB = XNCOMB + XLTSSA*(XLTSSA+1.0D0)*0.5
               END IF
               NCOMB = NCOMB + LTTSBL
               LSB = LSB + LTTSUP
               MXSOOB = MAX(MXSOOB,LTTSUP)
               NTTSBL = NTTSBL + 1
               LCOL = LCOL + NSSOB(IBSM,IBTP)
             END IF
  100      CONTINUE
  200     CONTINUE
          MXSB = MAX(MXSB,LSB)
        END IF
  300 CONTINUE
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' NGASDT : NCOMB XNCOMB ,NTTSBL',
     &               NCOMB,XNCOMB,NTTSBL
      END IF
*
      CALL QEXIT('NGASD')
*
      RETURN
      END
***********************************************************************
      SUBROUTINE OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX,
     &                  I_DO_BASSPC,IBASSPC)
*
* IWAY = 1 :
* obtain NOCCLS =
* Number of allowed ways of distributing the orbitals in the
* active spaces
*
* IWAY = 2 :
* OBTAIN NOCCLS and
* IOCCLS = allowed distributions of electrons
*
* Added Oct 98 : IBASSPC
* The basespace of
* a given class is the first space where this class occurs
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*. Input
      DIMENSION IGSMIN(NGAS),IGSMAX(NGAS)
*. Output
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION IBASSPC(*)
*. Local scratch
      DIMENSION IOCA(MXPNGAS),IOC(MXPNGAS)
*
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
         WRITE(6,*)  ' OCCLS in action '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' NGAS NEL ', NGAS,NEL
      END IF
*
      ISKIP = 1
      NOCCLS = 0
*. start with smallest allowed number
      DO IGAS = 1, NGAS
        IOCA(IGAS) = IGSMIN(IGAS)
      END DO
      NONEW = 0
      IFIRST = 1
*. Loop over possible occupations
 1000 CONTINUE
        IF(IFIRST.EQ.0) THEN
*. Next accumulated occupation
          CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
        END IF
        IF(NONEW.EQ.0) THEN
*. ensure that IOCA corresponds to an accumulating occupation,
*. i.e. a non-decreasing sequence
        IF(ISKIP.EQ.1) THEN
          KGAS = 0
          DO IGAS = 2, NGAS
            IF(IOCA(IGAS-1).GT.IOCA(IGAS)) KGAS = IGAS
          END DO
          IF(KGAS .NE. 0 ) THEN
            DO IGAS = 1, KGAS-1
              IOCA(IGAS) = IGSMIN(IGAS)
            END DO
            IOCA(KGAS) = IOCA(KGAS)+1
          END IF
        END IF
C?      WRITE(6,*) ' Another accumulated occupation: '
C?      CALL IWRTMA(IOCA,1,NGAS,1,NGAS)
*. corresponding occupation of each active space
        NEGA=0
        DO IGAS = 1, NGAS
          IF(IGAS.EQ.1) THEN
            IOC(IGAS) = IOCA(IGAS)
          ELSE
            IOC(IGAS) = IOCA(IGAS)-IOCA(IGAS-1)
            IF(IOC(IGAS).LT.0) NEGA = 1
          END IF
        END DO
C?      WRITE(6,*) ' Another occupation: '
C?      CALL IWRTMA(IOC,1,NGAS,1,NGAS)
        IFIRST = 0
*. Correct number of electrons
        IEL = IELSUM(IOC,NGAS)
        IF(IEL.EQ.NEL.AND.NEGA.EQ.0) THEN
          NOCCLS = NOCCLS + 1
          IF(IWAY.EQ.2) THEN
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Another allowed class : '
              CALL IWRTMA(IOC,1,NGAS,1,NGAS)
            END IF
            CALL ICOPVE(IOC,IOCCLS(1,NOCCLS),NGAS)
*
            IF(I_DO_BASSPC.EQ.1) THEN
              IBASSPC(NOCCLS) = IBASSPC_FOR_CLS(IOC)
            END IF
*
          END IF
        END IF
      END IF
      IF(NONEW.EQ.0) GOTO 1000
*
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' Number of Allowed occupation classes ', NOCCLS
         IF(IWAY.EQ.2.AND.NTEST.GE.20) THEN
           WRITE(6,*) ' Occupation classes '
           CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
         END IF
      END IF
*
      IF(I_DO_BASSPC.EQ.1) THEN
        WRITE(6,*) ' Base CI spaces for the classes '
        CALL IWRTMA(IBASSPC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE OCCLS_IN_CI(NOCCLS,IOCCLS,ICISPC,NINCCLS,INCCLS)
*
* A set of occupation classes are given.
* Find out the classes that are allowed for CI space ICISPC
*
*     Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
*. Specific input
      INTEGER IOCCLS(NGAS,NOCCLS)
*. Scratch
C     INTEGER IACOC(MXPNGAS)
*. Output
      INTEGER INCCLS(*)
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' OCCLS_IN_CI :  Input classes '
        CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCCLS)
      END IF
*
      NINCCLS = 0
      DO JOCCLS = 1, NOCCLS

*
CMI     INCLUDE = 0
        INCLUD = 0
        DO JJCMBSPC = 1, LCMBSPC(ICISPC)
          JCMBSPC = ICMBSPC(JJCMBSPC,ICISPC)
          I_AM_OK = 1
          DO IGAS = 1, NGAS
            IF(IGAS.EQ.1) THEN
              IEL = IOCCLS(1,JOCCLS)
            ELSE
              IEL = IEL + IOCCLS(IGAS,JOCCLS)
            END IF
            IF(IEL.LT.IGSOCCX(IGAS,1,JCMBSPC).OR.
     &      IEL.GT.IGSOCCX(IGAS,2,JCMBSPC)) I_AM_OK = 0
          END DO
CMI       IF(I_AM_OK .EQ. 1 ) INCLUDE = 1
          IF(I_AM_OK .EQ. 1 ) INCLUD = 1
        END DO
*
CMI     IF(INCLUDE.EQ.1) THEN
        IF(INCLUD.EQ.1) THEN
          NINCCLS = NINCCLS + 1
          INCCLS(JOCCLS) = 1
        ELSE
          INCCLS(JOCCLS) = 0
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output from OCCLS_IN_CI '
        WRITE(6,*) ' ========================'
        WRITE(6,*) ' CI space under study ', ICISPC
        WRITE(6,*) ' Number of occupation classes included ',NINCCLS
        WRITE(6,*) ' Active classes : 1 => active 0=> inactive '
        CALL IWRTMA(INCCLS,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE ORBINF(LUOUT,IPRNT)
*
* Obtain information about orbitals from shell information
*
* =====
* input
* =====
* Shell and symmetry information in /LUCINP/
*
* ======
* Output
* ======
* Orbital information in /ORBINP/
*
* Jeppe Olsen, Winter of 1991
*
#include "mxpdim.inc"
#include "lucinp.inc"
#include "cgas.inc"
*
#include "orbinp.inc"
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
************************************************
*                                              *
* Part 1 : From shell format to orbital format *
*                                              *
************************************************
      CALL OSPIR(NOSPIR,IOSPIR,PNTGRP,NIRREP,MXPIRR,MXPOBS,IPRNT)
*
* 2 : Shell information to orbital information for each group of orbital
*
*
* ===============
*     GAS case
* ===============
*
        DO IGAS = 1, NGAS
*. Shell => orbitals for each GAS space
          CALL SHTOOB(NGSSH(1,IGAS),NIRREP,MXPOBS,NSMOB,NOSPIR,
     &                IOSPIR,NGSOB(1,IGAS),NGSOBT(IGAS))
        END DO
*
*  ========================================================
*. Number of inactive, active, occupied , deleted orbitals
*  ========================================================
*
*
* current inactive and deleted orbitals are not identified so
       IGSINA = 0
       IGSDEL = 0
*
       call izero(ntoobs,nsmob)
       call izero(nocobs,nsmob)
       call izero(nacobs,nsmob)
*
       NTOOB = 0
       NACOB = 0
       NOCOB = 0
       DO IGAS = 1, NGAS
*. Inactive orbitals
         IF(IGAS.EQ.IGSINA) THEN
           CALL ICOPVE(NGSOB(1,IGAS),NINOBS,NSMOB)
           NINOB = NGSOBT(IGAS)
         END IF
*. Deleted orbitals
         IF(IGAS.EQ.IGSDEL) THEN
           CALL ICOPVE(NGSOB(1,IGAS),NDEOBS,NSMOB)
           NDEOB = NGSOBT(IGAS)
         END IF
*. Add to total number of orbitals
         CALL IVCSUM(NTOOBS,NTOOBS,NGSOB(1,IGAS),1,1,NSMOB)
         NTOOB = NTOOB + NGSOBT(IGAS)
*. Add to occupied orbitals
         IF(IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NOCOBS,NOCOBS,NGSOB(1,IGAS),1,1,NSMOB)
           NOCOB = NOCOB + NGSOBT(IGAS)
         END IF
*. Add to active orbitals
         IF(IGAS.NE.IGSINA.AND.IGAS.NE.IGSDEL) THEN
           CALL IVCSUM(NACOBS,NACOBS,NGSOB(1,IGAS),1,1,NSMOB)
           NACOB = NACOB + NGSOBT(IGAS)
         END IF
       END DO
* ===============================================
*. Well, report back
* ===============================================
        IF(NTEST.GT.0) THEN
          WRITE(LUOUT,*)
          WRITE(LUOUT,*) ' Number of orbitals per symmetry :'
          WRITE(LUOUT,*) ' ================================='
          WRITE(LUOUT,*)
          WRITE(LUOUT,'(1X,A,10I4,A)')
     &    '            Symmetry  ',(I,I = 1,NSMOB)
          WRITE(LUOUT,'(1X,A,2X,10A)')
     &    '           ========== ',('====',I = 1,NSMOB)
          DO IGAS = 1, NGAS
            WRITE(LUOUT,'(1X,A,I3,7X,A,10I4,8X,I3)')
     &      '   GAS',IGAS,'      ',(NGSOB(I,IGAS),I=1,NSMOB),
     &      NGSOBT(IGAS)
          END DO
*
          WRITE(LUOUT,*) ' Total number of orbitals ', NTOOB
          WRITE(LUOUT,*) ' Total number of occupied orbitals ', NOCOB
        END IF
*. Offsets for orbitals of given symmetry
        ITOOBS(1) = 1
        DO  ISMOB = 2, NSMOB
          ITOOBS(ISMOB) = ITOOBS(ISMOB-1)+NTOOBS(ISMOB-1)
        END DO
*
        IF(NTEST.GT.0) THEN
          WRITE(6,*) ' Offsets for orbital of given symmetry '
          CALL IWRTMA(ITOOBS,1,NSMOB,1,NSMOB)
        END IF

********************************************
*                                          *
* Part 2 : Reordering arrays for orbitals  *
*                                          *
********************************************
        CALL ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &       NOCOBS,NTOOBS,NTOOB,
     &       IREOST,IREOTS,ISMFTO,ITPFSO,
     &       IBSO,NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &       ISMFSO,ITPFTO,NOBPT,IPRNT)
*
      END
***********************************************************************
      SUBROUTINE ORBORD_GAS(NSMOB,MXPOBS,MXPNGAS,NGAS,NGSOB,NGSOBT,
     &                  NOCOBS,NTOOBS,NTOOB,
     &                  IREOST,IREOTS,ISFTO,ITFSO,IBSO,
     &                  NTSOB,IBTSOB,ITSOB,NOBPTS,IOBPTS,
     &                  ISFSO,ITFTO,NOBPT,IPRNT)
*
* Obtain Reordering arrays for orbitals
* ( See note below for assumed ordering )
*
*
* GAS version
*
* =====
* Input
* =====
*  NSMOB  : Number of orbital symmetries
*  MXPOBS : Max number of orbital symmetries allowed by program
*  MXPNGAS: Max number of GAS spaces allowed by program
*  NGAS   : Number of GAS spaces
*  NGSOB  : Number of GAS orbitals per symmetry and space
*  NGSOBT : Number of GAS orbitals per space
*  NOCOBS : Number of occupied orbitals per symmetry
*  NTOOBS : Number of orbitals per symmetry,all types
*
* ======
* Output
* ======
*  IREOST : Reordering array symmetry => type
*  IREOTS : Reordering array type     => symmetry
*  ISFTO  : Symmetry array for type ordered orbitals
*  ITFSO  : Type array for symmetry ordered orbitals( not activated )
*  IBSO   : First orbital of given symmetry ( symmetry ordered )
*  NOBPTS : Number of orbitals per subtype and symmetry
*  IOBPTS : Off sets for orbitals of given subtype and symmetry
*           ordered according to input integrals
*
* ISFSO  : Symmetry of orbitals, symmetry ordereing
* ITFTO  : Type of orbital, type ordering
*
* Jeppe Olsen, Winter 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NGSOB(MXPOBS,MXPNGAS),NOCOBS(*),NTOOBS(*)
      DIMENSION NGSOBT(MXPNGAS)
*. Output
      DIMENSION IREOST(*),IREOTS(*),ISFTO(*),ITFSO(*),IBSO(*)
      DIMENSION ISFSO(*),ITFTO(*)
      DIMENSION NOBPTS(MXPNGAS ,*),IOBPTS(MXPNGAS ,*)
      DIMENSION NOBPT(MXPNGAS )

* ==========================
* Note on order of orbitals
* ==========================
*
* The orbitals are supposed to be imported ordered symmetry-type
* ordered as
*
* Loop over symmetries of orbitals
*  Loop over GAS spaces
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of Loop over Gas spaces
* End of loop over symmetries
*
* Internally the orbitals are reordered to type symmetry order
* where the outer loop is over types and the inner loop is
* over symmetries, i.e.
*
* Loop over GAS spaces
*  Loop over symmetries of orbitals
*   Loop over orbitals of this sym and GAS
*   End of Loop over orbitals
*  End of loop over symmetries
* End of Loop over Gas spaces
*
*. 1:  Construct ISFTO, ITFTO, IREOST,IREOTS,NOBPTS,IOBPTS
*
      ITSOFF = 1
      DO IGAS = 1, NGAS
        DO ISYM = 1, NSMOB
          IF(ISYM.EQ.1) THEN
            IBSSM = 1
          ELSE
            IBSSM = IBSSM + NTOOBS(ISYM-1)
          END IF
          NPREV = 0
          DO JGAS = 1, IGAS-1
            NPREV = NPREV + NGSOB(ISYM,JGAS)
          END DO
          IADD = 0
          NOBPTS(IGAS,ISYM) = NGSOB(ISYM,IGAS)
          IOBPTS(IGAS,ISYM) = ITSOFF
C         NOBPTS(ISYM,IGAS) = NGSOB(ISYM,IGAS)
C         IOBPTS(ISYM,IGAS) = ITSOFF
          DO IORB = ITSOFF,ITSOFF+NGSOB(ISYM,IGAS)-1
            IADD = IADD + 1
            IREOTS(IORB) = IBSSM-1+NPREV+IADD
            IREOST(IBSSM-1+NPREV+IADD) = IORB
            ITFTO(IORB) = IGAS
            ISFTO(IORB) = ISYM
          END DO
          ITSOFF = ITSOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 2 : ISFSO,ITFSO
*
      ISTOFF = 1
      DO ISYM = 1, NSMOB
        DO IGAS = 1, NGAS
          DO IORB = ISTOFF,ISTOFF+NGSOB(ISYM,IGAS)-1
            ISFSO(IORB) = ISYM
            ITFSO(IORB) = IGAS
          END DO
          ISTOFF = ISTOFF + NGSOB(ISYM,IGAS)
        END DO
      END DO
*
* 3 IBSO, NOBPT
*
      IOFF = 1
      DO ISM = 1, NSMOB
       IBSO(ISM) = IOFF
       IOFF = IOFF + NTOOBS(ISM)
      END DO
      DO IGAS = 1, NGAS
        NOBPT(IGAS) = NGSOBT(IGAS)
      END DO
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' Output from ORBORD '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' Symmetry of orbitals , type ordered '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Symmetry => type reordering array '
        CALL IWRTMA(IREOST,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' Type => symmetry reordering array '
        CALL IWRTMA(IREOTS,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' IBSO array '
        CALL IWRTMA(IBSO,1,NSMOB,1,NSMOB)
*
        WRITE(6,*) ' NOBPTS '
        CALL IWRTMA(NOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
        WRITE(6,*) ' NOBPT '
        CALL IWRTMA(NOBPT,NGAS,1,MXPNGAS,1)
        WRITE(6,*) ' IOBPTS '
        CALL IWRTMA(IOBPTS,NGAS,NSMOB,MXPNGAS,MXPOBS)
*
        WRITE(6,*) ' ISFTO array : '
        CALL IWRTMA(ISFTO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFSO array : '
        CALL IWRTMA(ITFSO,1,NTOOB,1,NTOOB)
*
        WRITE(6,*) ' ISFSO array : '
        CALL IWRTMA(ISFSO,1,NTOOB,1,NTOOB)
        WRITE(6,*) ' ITFTO array : '
        CALL IWRTMA(ITFTO,1,NTOOB,1,NTOOB)
      END IF
*

      RETURN
      END
***********************************************************************
      SUBROUTINE OSPIR(NOSPIR,IOSPIR,PNTGRP,NIRREP,MXPIRR,MXPOBS,IPRNT)
*
* Number and symmetries of orbitals corresponding to a given shell
*
* =====
* Input
* =====
*
*   PNTGRP  : type of pointgroup
*         = 1 => D2h or a subgroup of D2H
*         = 2 => C inf v
*         = 3 => D inf h
*         = 4 => O 3
*   NIRREP : Number of irreducible representations per point group
*   MXPIRR : Largest allowed number of shell irreps
*   MXPOBS : Largest allowed number of orbital symmetries
*
* ======
* Output
* ======
*
*   NOSPIR : Number of orbital symmetries per irrep
*   IOSPIR : Orbital symmetries corresponding to a given irrep
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER PNTGRP
*. Output
      DIMENSION NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR)
*
      IF(PNTGRP.EQ.1) THEN
*=====
*.D2h
*=====
        NSMOB = 0
        DO 10 IRREP = 1, 8
          NOSPIR(IRREP) = 1
          IOSPIR(1,IRREP) = IRREP
   10   CONTINUE
      ELSE IF(PNTGRP.EQ.2) THEN
* =========
*. C inf V
* =========
* orbital symmetry is numbered as IML - MNMLOB + 1
        MNMLOB = -(NIRREP-1)
        DO 20 IRREP = 1, NIRREP
*.Irrep I contains orbitals with ML = -(IRREP-1),+(IRREP-1)
          IF(IRREP.EQ.1) THEN
            NOSPIR(IRREP) = 1
            IOSPIR(1,IRREP) = IRREP - 1 - MNMLOB + 1
          ELSE
            NOSPIR(IRREP) = 2
            IOSPIR(1,IRREP) = -(IRREP - 1) - MNMLOB + 1
            IOSPIR(2,IRREP) =  (IRREP - 1) - MNMLOB + 1
          END IF
   20   CONTINUE
      ELSE IF(PNTGRP.EQ.3) THEN
* ========
*. D inf H
* ========
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        IRREP = 0
        DO 35 IPARI = 1, 2
          IADD = (IPARI-1)*NMLOB
          DO 30 ML = 0,MXMLOB
            IRREP = IRREP + 1
            IF(ML.EQ.0) THEN
              NOSPIR(IRREP) = 1
              IOSPIR(1,IRREP) = IADD + ML - MNMLOB + 1
            ELSE
              NOSPIR(IRREP) = 2
              IOSPIR(1,IRREP) = IADD - ML - MNMLOB + 1
              IOSPIR(2,IRREP) = IADD + ML - MNMLOB + 1
            END IF
   30     CONTINUE
   35   CONTINUE

      ELSE IF(PNTGRP.EQ.4) THEN
* =====
*. O 3
* =====
* orbital symmetry is numbered as (PARITY-1) * NMLOB + IML - MNMLOB + 1
        MXMLOB =  NIRREP/2-1
        MNMLOB = -MXMLOB
        NMLOB =   NIRREP - 1
        DO 45 L = 0, NIRREP - 1
          IF(MOD(L,2).EQ.0) THEN
            IPARI = 1
          ELSE
            IPARI = 2
          END IF
          IADD = (IPARI-1)*NMLOB
          IRREP = L + 1
          NOSPIR(IRREP) = 2 * L + 1
          ICOMP = 0
          DO 40 ML = MNMLOB,MXMLOB
            ICOMP = ICOMP + 1
            IOSPIR(ICOMP,IRREP) = IADD + ML - MNMLOB + 1
   40     CONTINUE
   45   CONTINUE
      ELSE
        WRITE(6,*) ' Sorry  PNTGRP out of range , PNTGRP = ', PNTGRP
        WRITE(6,*) ' ORBIRR fatally wounded '
        Call Abend1( 5 )
      END IF
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' OSPIR speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Number of orbitals per irrep '
        CALL IWRTMA(NOSPIR,1,NIRREP,1,NIRREP)
        WRITE(6,*) ' Orbital symmetries per irrep '
        DO 100 IRREP = 1, NIRREP
          CALL IWRTMA(IOSPIR(1,IRREP),1,NOSPIR(IRREP),1,NOSPIR(IRREP))
  100   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE PERT_SUBSPACE(NPERT,H0,V,S,ECORE)
*
* Perturbation calculation has been performed
* giving a subspace of correction vectors.
*
* Analyze this subspace
*
* H0, V and S are matrices in subspace
*
*. Jeppe Olsen, July 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLMAT1, KLMAT2, KLMAT3, KLMAT4, KLMAT5, KLMAT6
      INTEGER*8 KLVEC1, KLVEC2, KLVEC3, KLC, KLEN
      ! for addressing of WORK
      DIMENSION S(*),H0(*),V(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ', IDUM,'PERT_S')
*
      NTEST = 10
      IF(NTEST.GE.1) THEN
         WRITE(6,*)
         WRITE(6,*) ' Welcome to pert_subspace '
         WRITE(6,*) ' ======================== '
         WRITE(6,*)
      END IF
*. Dimension of subspace : number of correction vectors + 1
      NDIM = NPERT + 1
*
*. A bit of local memory
      LENNY = NDIM ** 2
      IDUM = 0
      CALL MEMMAN(KLMAT1,LENNY,'ADDL  ',2,'MAT1  ')
      CALL MEMMAN(KLMAT2,LENNY,'ADDL  ',2,'MAT2  ')
      CALL MEMMAN(KLMAT3,3*LENNY,'ADDL  ',2,'MAT3  ')
      CALL MEMMAN(KLMAT4,LENNY,'ADDL  ',2,'MAT4  ')
      CALL MEMMAN(KLMAT5,LENNY,'ADDL  ',2,'MAT5  ')
      CALL MEMMAN(KLMAT6,LENNY,'ADDL  ',2,'MAT6  ')
      CALL MEMMAN(KLVEC1,NDIM ,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,NDIM ,'ADDL  ',2,'VEC2  ')
      CALL MEMMAN(KLVEC3,NDIM ,'ADDL  ',2,'VEC3  ')
*
*. Orthonormalize using symmetric orthogonalization
*
*. S**(-1/2)
       SIGN = 1.0D0
       CALL TRIPAK_LUCI(WORK(KLMAT5),S,2,NDIM,NDIM,SIGN)
       WRITE(6,*) ' Input S matrix '
       CALL WRTMT_LU(WORK(KLMAT5),NDIM,NDIM,NDIM,NDIM)
C          SQRTMT(A,NDIM,ITASK,ASQRT,AMSQRT,SCR)
      CALL SQRTMT(WORK(KLMAT5),NDIM,2,WORK(KLMAT2),
     &            WORK(KLMAT1),WORK(KLMAT3))
      IF(NTEST.GE.5) THEN
         WRITE(6,*) ' S-1/2 matrix '
         CALL WRTMT_LU(WORK(KLMAT1),NDIM,NDIM,NDIM,NDIM)
      END IF
*. Transform H0 and V to orthogonal basis
C          TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
      CALL TRAN_SYM_BLOC_MAT(H0,WORK(KLMAT1),1,NDIM,WORK(KLMAT2),
     &                       WORK(KLMAT3))
      CALL COPVEC(WORK(KLMAT2),H0,NDIM*(NDIM+1)/2)
      CALL TRAN_SYM_BLOC_MAT(V ,WORK(KLMAT1),1,NDIM,WORK(KLMAT2),
     &                       WORK(KLMAT3))
      CALL COPVEC(WORK(KLMAT2),V ,NDIM*(NDIM+1)/2)
*
      WRITE(6,*) ' H0 in symmetric orthogonal basis '
      CALL PRSYM (H0,NDIM)
      WRITE(6,*) ' V  in symmetric orthogonal basis '
      CALL PRSYM (V ,NDIM)
*. Find basis where H0 is diagonal and transform
*. Diagonalize H0, eigenvectors in MAT2
      CALL EIGEN_LUCI(H0,WORK(KLMAT2),NDIM,0,1)
      CALL COPDIA(H0,WORK(KLMAT3),NDIM,1)
*. And put back
      ZERO = 0.0D0
      CALL SETVEC(H0,ZERO,NDIM*(NDIM+1)/2)
      DO I = 1, NDIM
        H0(I*(I+1)/2) = WORK(KLMAT3-1+I)
      END DO
*. Transform  V to basis that diagonalizes H0
      CALL TRAN_SYM_BLOC_MAT(V ,WORK(KLMAT2),1,NDIM,WORK(KLMAT4),
     &                       WORK(KLMAT3))
      CALL COPVEC(WORK(KLMAT4),V,NDIM*(NDIM+1)/2)
*
      WRITE(6,*) ' H0 in basis of sub space zero order states '
      WRITE(6,*) ' ========================================== '
      WRITE(6,*)
      CALL PRSYM(H0,NDIM)
      WRITE(6,*)
      WRITE(6,*) ' V  in basis of sub space zero order states '
      WRITE(6,*) ' ========================================== '
      CALL PRSYM(V,NDIM)
*. Eigenvalues H = H0+V, eigenvectors in MAT4
      ONE = 1.0D0
      CALL VECSUM(WORK(KLMAT3),H0,V,ONE,ONE,NDIM*(NDIM+1)/2)
      CALL EIGEN_LUCI(WORK(KLMAT3),WORK(KLMAT4),NDIM,1,1)
      CALL COPDIA(WORK(KLMAT3),WORK(KLVEC1),NDIM,1)
      DO I = 1, NDIM
        WORK(KLVEC1-1+I) = WORK(KLVEC1-1+I) + ECORE
      END DO
*.
      WRITE(6,*)
      WRITE(6,*) ' Eigenvalues of H (with core-energy)in subspace '
      WRITE(6,*) ' =============================================== '
      WRITE(6,*)
      CALL WRTMT_LU(WORK(KLVEC1),NDIM,1,NDIM,1)
*
*. Perturbation expansion in subspace
*
*. Expand H0 and V  to complete matrices
C       TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
      SIGN = 1.0D0
      CALL TRIPAK_LUCI(WORK(KLMAT5),H0,2,NDIM,NDIM,SIGN)
      CALL TRIPAK_LUCI(WORK(KLMAT6),V ,2,NDIM,NDIM,SIGN)
*
      MAXORD = 500
      LEN = NDIM*(1+MAXORD)
      CALL MEMMAN(KLC,LEN,'ADDL  ',2,'KLC   ')
      CALL MEMMAN(KLEN,MAXORD+1,'ADDL  ',2,'KLC   ')
*. Zero order state
      ZERO = 0.0D0
      CALL SETVEC(WORK(KLC),ZERO,NDIM)
      WORK(KLC) = ONE
C          MATPERT(H0,V,NDIM,NORD,EN,C,VEC1,VEC2,VEC3)
      CALL MATPERT(WORK(KLMAT5),WORK(KLMAT6),NDIM,MAXORD,
     &             WORK(KLEN),WORK(KLC),WORK(KLVEC1),WORK(KLVEC2),
     &             WORK(KLVEC3),ECORE )
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PERTSB')
*
      RETURN
      END
***********************************************************************
      SUBROUTINE PNT2DM(I12SM,NSMOB,NSMSX,OSXO,IPSM,JPSM,
     &                  IJSM,ISM2,IPNTR,MXPOBS)
*
* Pointer to two dimensional array
*
* =====
* Input
* =====
* I12SM  : ne.0 => restrict to lower half
*          eq.0 => complete matrix
* NSMOB : Number of orbital symmetries
* NSMSX : Number of SX      symmetries
* OSXO  : Symmetry of orbital, SX => symmetry of other orbital
* IPSM : Number of orbitals per symmetry for index 1
* JPSM : Number of orbitals per symmetry for index 2
* IJSM  : Symmetry of two index array
*
* =======
* Output
* =======
* IPNTR : Pointer to block with first index of given symmetry
*         = 0 indicates forbidden block
* ISM2  : symmetry of second index for given first index
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER OSXO(MXPOBS,2*MXPOBS),IPSM(*),JPSM(*)
*.Output
      DIMENSION IPNTR(*),ISM2(*)
*
      CALL ISETVC(IPNTR,0,NSMOB)
      CALL ISETVC(ISM2 ,0,NSMOB)
      IOFF = 1
      DO 100 ISM = 1,NSMOB
        JSM = OSXO(ISM,IJSM)
        IF(JSM.EQ.0) GOTO 100
        IF(I12SM.EQ.0.OR.ISM.GE.JSM) THEN
*. Allowed block
          IPNTR(ISM) = IOFF
          ISM2(ISM) = JSM
          IF(I12SM.GT.0.AND.ISM.EQ.JSM) THEN
            IOFF = IOFF + IPSM(ISM)*(IPSM(ISM)+1)/2
          ELSE
            IOFF = IOFF + IPSM(ISM)*JPSM(JSM)
          END IF
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' dimension of two-dimensional array ',IOFF-1
      END IF
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Pointer '
        CALL IWRTMA(IPNTR,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' Symmetry of other array '
        CALL IWRTMA(ISM2,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE PNT4DM(NSMOB,NSMSX,MXPOBS,NO1PS,NO2PS,NO3PS,NO4PS,
     &           IDXSM,ADSXA,SXDXSX,IS12,IS34,IS1234,IPNTR,ISM4A,
     &           ADASX)
*
* Pointer for 4 dimensionl array with total symmetry IDXSM
* Pointer is given as 3 dimensional array corresponding
* to the first 3 indeces
* Symmetry of last index is give by ISM4
*
* IS12 (0,1,-1)   : Permutational symmetry between indeces 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between indeces 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between indeces 12 and 34
*
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER ADASX(MXPOBS,MXPOBS)
*. Specific input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*.Output
      INTEGER IPNTR(NSMOB,NSMOB,NSMOB),ISM4A(NSMOB,NSMOB,NSMOB)
*
      CALL ISETVC(IPNTR,0,NSMOB ** 3 )
      CALL ISETVC(ISM4A,0,NSMOB ** 3 )
*
C?    WRITE(6,*) 'NO1PS NO2PS NO3PS NO4PS '
C?    CALL IWRTMA(NO1PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO2PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO3PS,1,NSMOB,1,NSMOB)
C?    CALL IWRTMA(NO4PS,1,NSMOB,1,NSMOB)
      IOFF= 1
      DO 10 I1SM = 1, NSMOB
        DO 20 I2SM = 1, NSMOB
          I12SM = ADASX(I1SM,I2SM)
          I34SM = SXDXSX(I12SM,IDXSM)
          IF(I34SM.EQ.0) GOTO 20
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 20
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMOB+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
          DO 30 I3SM = 1, NSMOB
            I4SM = ADSXA(I3SM,I34SM)
            IF(I4SM.EQ.0) GOTO 30
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 30
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMOB+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 30
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF+ N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF= IOFF + N12*(N12+1)/2
            ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              IPNTR(I1SM,I2SM,I3SM) = IOFF
              ISM4A(I1SM,I2SM,I3SM) = I4SM
              IOFF=  IOFF+ N12*(N12-1)/2
            END IF
C?          WRITE(6,*) ' I1SM I2SM I3SM I4SM    IOFF'
C?          WRITE(6,'(1X,4I4,I9)')   I1SM,I2SM,I3SM,I4SM,IOFF
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*
*
C?    WRITE(6,*) ' PNT4DM , 64 elemets of IPNTR '
C?    call IWRTMA(IPNTR,1,64,1,64)
      NTEST = 0
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Length of 4 index array ', IOFF - 1
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE PTWFNRM(LU0,LUN,N,VEC1,VEC2,LUSCR1,LUSCR2,PTNORM)
*
* Norm of wavefunction obtained as perturbation expansion
* with intermediate normalization
*
* |0> = sum (i=0,n) |0(i)>
*
* Jeppe Olsen, Oct 97
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPRDD
*. Scratch blocks
      DIMENSION VEC1(*), VEC2(*)
*. Local scratch
      PARAMETER(MXORD = 100)
      DIMENSION SCR(MXORD)
*
      IF(N.GT.MXORD) THEN
        WRITE(6,*) ' PTWFNRM : N > MXORD : ',N,MXORD
        WRITE(6,*) ' Increase parameter MXORD '
        Call Abend2(       ' PTWFNRM : N > MXORD  ' )
      END IF
*
      LBLK = -1
      ONE = 1.0D0
*
      CALL COPVCD(LU0,LUSCR1,VEC1,1,LBLK)
      CALL SETVEC(SCR,ONE,N)
      CALL MVCSMD2(LUN,SCR,ONE,LUSCR1,LUSCR2,VEC1,VEC2,N,1,LBLK)
      PTNORM = INPRDD(VEC1,VEC2,LUSCR2,LUSCR2,1,LBLK)
      PTNORM = SQRT(PTNORM)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' PTWFNRM: Norm of summed vector ',PTNORM
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE REPART_CIV(IBATCH,NBATCH,LBATCH,LEBATCH,
     &                  MXLNG,ICLS_A,IBLK_TO_CLS,NCLS,NBLKS,
     &                  LBLOCK_VEC)
* A CI vector is defined through IBATCH (generated by PART_CIV)
*
*. Divide into batches, length atmost MXLNG, so only blocks
*  that are flagged active by ICLS_A are included
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : Inactive blocks are flagged by setting the first element
*          negative ( -1 * the original value )
*
* Input
*
* IBATCH :
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IBATCH(8,*)
      INTEGER ICLS_A(*)
*. General input
       INTEGER IBLK_TO_CLS(*)
*. Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
C     INTEGER I1BATCH(*)
      INTEGER LBLOCK_VEC(*)
*
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     REPART_CIV      '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*)
C       WRITE(6,*) ' Active classes '
C       CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
C       write(6,*) ' MXLNG = ',MXLNG
      END IF
*
      NBATCH = 0
      LENGTH = 0
      LENGTHP= 0
      IIBLOCK = 0
C     I1BATCH(1) = 1
*. Loop over blocks in batch
      DO 1000 IBLOCK = 1, NBLKS
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IF(ICLS_A(IBLK_TO_CLS(IBLOCK)).EQ.0) THEN
*. Block is inactive
          IBATCH(1,IBLOCK) = - ABS(IBATCH(1,IBLOCK))
          IIBLOCK = IIBLOCK + 1
          LBLOCK_VEC(IBLOCK) = -IBATCH(8,IBLOCK)
        ELSE
*. Block belongs to the active !
          IBATCH(1,IBLOCK) =   ABS(IBATCH(1,IBLOCK))
          LBLOCK =  IBATCH(7,IBLOCK)
          LBLOCKP =  IBATCH(8,IBLOCK)
          LBLOCK_VEC(IBLOCK) = LBLOCKP
          IF(LENGTH+LBLOCK.LE.MXLNG) THEN
            LENGTH = LENGTH + LBLOCK
            LENGTHP= LENGTHP+ LBLOCKP
            IIBLOCK = IIBLOCK + 1
          ELSE IF(LENGTH+LBLOCK.GT.MXLNG) THEN
*. This batch was finished by previous block, goto next batch
            NBATCH = NBATCH + 1
            LEBATCH(NBATCH) = LENGTHP
            LBATCH (NBATCH)  = IIBLOCK
C           I1BATCH(NBATCH+1) = IBLOCK
*. Current block is first block in new batch
            IIBLOCK = 1
            IBATCH(5,IBLOCK) = 1
            IBATCH(6,IBLOCK) = 1
            LENGTHP = LBLOCKP
            LENGTH  = LBLOCK
          END IF
        END IF
 1000 CONTINUE
*. Final batch
      IF( LENGTH .NE. 0 ) THEN
        NBATCH = NBATCH + 1
        LBATCH(NBATCH) = IIBLOCK
        LEBATCH(NBATCH) = LENGTHP
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' REPART.. : Number of batches ', NBATCH
      END IF
      IF(NTEST.GE.100) THEN
        IOFF = 1
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = IOFF, IOFF+LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
          IOFF = IOFF + LBATCH(JBATCH)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE RSMXMN_LUCI(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,
     &                       MIN1,MAX1,MIN3,MAX3,NTEST)
*
* Construct accumulated MAX and MIN arrays for a RAS set of strings
*
      IMPLICIT REAL*8           ( A-H,O-Z)
      DIMENSION  MINEL(*),MAXEL(*)
*
      NORB = NORB1 + NORB2 + NORB3
*. accumulated max and min in each of the three spaces
*. ( required max and min at final orbital in each space )
COLD  MIN1A = MIN1
      MIN1A = MAX(MIN1,NEL-MAX3-NORB2)
      MAX1A = MAX1
*
      MIN2A = NEL - MAX3
      MAX2A = NEL - MIN3
*
      MIN3A = NEL
      MAX3A = NEL
*
      DO 100 IORB = 1, NORB
        IF(IORB .LE. NORB1 ) THEN
          MINEL(IORB) = MAX(MIN1A+IORB-NORB1,0)
          MAXEL(IORB) = MIN(IORB,MAX1A)
        ELSE IF ( NORB1.LT.IORB .AND. IORB.LE.(NORB1+NORB2)) THEN
          MINEL(IORB) = MAX(MIN2A+IORB-NORB1-NORB2,0)
          IF(NORB1 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1))
          MAXEL(IORB) = MIN(IORB,MAX2A)
        ELSE IF ( IORB .GT. NORB1 + NORB2 ) THEN
          MINEL(IORB) = MAX(MIN3A+IORB-NORB,0)
          IF(NORB1+NORB2 .GT. 0 )
     &    MINEL(IORB) = MAX(MINEL(IORB),MINEL(NORB1+NORB2))
          MAXEL(IORB) = MIN(IORB,MAX3A)
        END IF
  100 CONTINUE
*
      IF( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' Output from RSMXMN '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SHTOOB(NSHPIR,NIRREP,MXPOBS,NSMOB,NOSPIR,IOSPIR,
     &           NOBPS,NOB)
*
* Number of shells per irrep => Number of orbitals per symmetry
*
* =====
* Input
* =====
*
*  NSHPIR : Number of shells per irrep
*  NIRREP : Number of irreps
*  MXPOBS : Largest allowed number of orbitals symmetries
*  NSMOB  : Number of orbital symmetries
*  NOSPIR : Number of orbital symmetries per irrep
*  IOSPIR : Orbital symmetries per irrep
*
* ======
* Output
* ======
*  NOBPS  : Number of orbitals per symmetry
*  NOB    : Number of orbitals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NSHPIR(*),NOSPIR(*),IOSPIR(MXPOBS,*)
*. Output
      DIMENSION NOBPS(*)
      CALL ISETVC(NOBPS,0,NSMOB)
      NOB = 0
      DO 100 IRREP = 1, NIRREP
        DO 90 ISM = 1, NOSPIR(IRREP)
          IISM = IOSPIR(ISM,IRREP)
          NOBPS(IISM) = NOBPS(IISM) + NSHPIR(IRREP)
          NOB = NOB + NSHPIR(IRREP)
   90   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' SHTOOB Speaking '
         WRITE(6,*) ' =============== '
         WRITE(6,*) ' Number of orbitals obtained ', NOB
         WRITE(6,*) ' Number of orbitals per symmetry '
         CALL IWRTMA(NOBPS,1,NSMOB,1,NSMOB)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SIGVST(ISGVST,NSMST)
*
* Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISGVST(*)
*
      DO 100 ISM = 1, NSMST
C            MLSM(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM(IML,IPARI,ISM,'ST',2)
        MIML = - IML
        CALL MLSM(MIML,IPARI,MISM,'ST',1)
        ISGVST(ISM) = MISM
  100 CONTINUE
*
      NTEST = 1
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' ISGVST array '
        WRITE(6,*) ' ============ '
        CALL IWRTMA(ISGVST,1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*
* ISMOST(ISYM,ITOTSM) : Symmetry of an internal state is ITOTSM
*                       if symmetry of 1 string is ISYM, the
*                       symmetry of the other string is
*                       ISMOST(ISYM,ITOTSM)
*
* Jeppe Olsen , Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISMOST(MXPCSM,MXPCSM)
*
      DO 1000 ITOTSM = 1, NSMCI
       DO 900 ISTSM  = 1, NSMST
C            SYMCOM(ITASK,IOBJ,I1,I2,I12)
        CALL SYMCOM(2,1,ISTSM,JSTSM,ITOTSM)
        ISMOST(ISTSM,ITOTSM) = JSTSM
  900  CONTINUE
 1000 CONTINUE
*
      NTEST = 0
      IF( NTEST.NE. 0 ) THEN
        WRITE(6,*) ' ==============='
        WRITE(6,*) ' Info from SMOST '
        WRITE(6,*) ' ==============='
        DO 1010 ITOTSM = 1, NSMCI
          WRITE(6,*) ' ISMOST array for ITOTSM = ', ITOTSM
          CALL IWRTMA(ISMOST(1,ITOTSM),1,NSMST,1,NSMST)
 1010   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE STSTSM(STSTSX,STSTDX,NSMST)
*
* construct  STSTSX and STSTDX giving
* symmetry of sx (dx) connecting two given string symmetries
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER STSTSX(NSMST,NSMST),STSTDX(NSMST,NSMST)
*
      DO 100 ILSTSM = 1, NSMST
        DO 50 IRSTSM = 1, NSMST
          CALL SYMCOM(1,5,ISXSM,IRSTSM,ILSTSM)
          CALL SYMCOM(1,6,IDXSM,IRSTSM,ILSTSM)
          STSTSX(ILSTSM,IRSTSM) = ISXSM
          STSTDX(ILSTSM,IRSTSM) = IDXSM
   50   CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' STSTSM : STSTSX, STSTDX '
        CALL IWRTMA(STSTSX,NSMST,NSMST,NSMST,NSMST)
        CALL IWRTMA(STSTDX,NSMST,NSMST,NSMST,NSMST)
      END IF
*
      RETURN
      END
***********************************************************************
*
      SUBROUTINE SXTYP_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP)
*
* Two supergroups are given. Find single excitations that connects
* these supergroups
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS)
*. Output
      DIMENSION ITP(*),JTP(*)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one sngle excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left  supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SYM_DIST_FOR_SPGRP(IGRP,NIGRP,ISM,NDIST,
     &           IDIST,LDIST,LENGTH,MXDIST)
*
* Symmetry distributions of a combination of groups
*
*
*. Input
* ======
*    IGRP  : The groups of the supergroup
*    NIGRP : Number of groups in the supergroup
*    ISM   : Required total symmetry
*    MXDIST: Largest allowed number of distributions
*
*. Output
* =======
*    NDIST : Number of symmetry distributions
*    IDIST : The symmetry distributions
*    LDIST : Length of each distribution
*    LENGTH: Total length of distributions with the given sym
*
* Jeppe Olsen, Sept. 97
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "csm.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "strbas.inc"
*. Specific Input
      DIMENSION IGRP(NIGRP)
*. Local scratch
      DIMENSION MXVAL(MXPNGAS),MNVAL(MXPNGAS),JDIST(MXPNGAS)
      DIMENSION LSMGP(MXPOBS,MXPNGAS)
*. Output
      DIMENSION IDIST(NIGRP,MXDIST),LDIST(MXDIST)
*. Max and Min arrays for symmetries
*
      NTEST = 00
      IF(NIGRP.EQ.0) THEN
*. Trivial zero supergroup, seperately handled to avoid infinite loops
       IF(ISM.EQ.1) THEN
         NDIST = 1
         LDIST(1) = 1
         LENGTH = 1
       ELSE
         NDIST = 0
         LENGTH = 0
       END IF
      ELSE
*. Nontrivial distributions
      CALL MINMAX_FOR_SYM_DIST(NIGRP,IGRP,MNVAL,MXVAL,NDISTX)
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Entering SYM_DIST ... <<'
*. Dimensions of given group and symmetry
      DO JGRP = 1, NIGRP
        CALL ICOPVE2(WORK(KNSTSGP(1)),(IGRP(JGRP)-1)*NSMST+1,
     &               NSMST,LSMGP(1,JGRP))
      END DO
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' LSMGP : '
        CALL IWRTMA(LSMGP,NSMST,NIGRP,MXPOBS,NIGRP)
      END IF
*. And generate symmetry distributions
      IFIRST = 1
      NDIST = 0
      LENGTH = 0
 1000 CONTINUE
C       NEXT_SYM_DISTR(NGAS,MINVAL,MAXVAL,ISYM,ISYM_TOT,IFIRST,NONEW)
        CALL NEXT_SYM_DISTR(NIGRP,MNVAL,MXVAL,
     &       JDIST,ISM,IFIRST,NONEW)
        IF(NONEW.EQ.0) THEN
           NDIST = NDIST + 1
           IF(NDIST.GT.MXDIST) THEN
             WRITE(6,*) 'SYM_DIST_FOR_SPGRP in problems '
             WRITE(6,*) ' NDIST .gt. MXDIST '
             WRITE(6,*) ' NDIST, MXDIST = ',NDIST,',',MXDIST
             Call Abend2( 'SYM_DIST_FOR_SPGRP Increase MXDIST' )
           END IF
           LDIM = 1
           DO JGRP = 1, NIGRP
            LDIM = LDIM*LSMGP(JDIST(JGRP),JGRP)
            IDIST(JGRP,NDIST) = JDIST(JGRP)
           END DO
           LDIST(NDIST) = LDIM
           LENGTH = LENGTH + LDIM

         END IF
      IF(NONEW.EQ.0) GOTO 1000
      END IF
*     ^ Switch for nontrivial/trivial distribution
*
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Symmetry distributions generated for : '
        WRITE(6,*) '    Total symmetry :  ', ISM
        WRITE(6,'(A,10I3,(10I3))')
     &             '     Groups         : ', (IGRP(JGRP),JGRP=1,NIGRP)
        WRITE(6,*)
        WRITE(6,*) '    Number of symmetry distributions ', NDIST
        WRITE(6,*) '    Total dimension                  ',LENGTH
        WRITE(6,*)
        WRITE(6,*) ' Number, Length, Symmetry distributions '
        WRITE(6,*) ' ======================================='
        DO KDIST = 1, NDIST
          WRITE(6,'(I5,I7,4X,10I3,(16X,10I3))')
     &    KDIST,LDIST(KDIST),(IDIST(JGRP,KDIST),JGRP=1,NIGRP)
        END DO
      END IF
*
      IF(NTEST.GE.1000) WRITE(6,*) ' >> Leaving SYM_DIST ... << '
*
      RETURN
      END
***********************************************************************
*
* lucia6.f
*
* ====================
* CI response routines lucia version
* ====================
*
* Arbitrary order perturbation theory for one or two perturbations
* Cauchy moments, effective spectra , Dispersion coefficients
*
* Jeppe Olsen
*
*
      SUBROUTINE SYM_FOR_OP(OPER,IXYZ,ISYM)
*
* Symmetry of perturbation operator OPER
* obtained from XYZ characters in start of OPER
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      CHARACTER*8 OPER
*. General input
      INTEGER IXYZ(3)
#include "multd2h.inc"
*
*. Number of cartesian components
      NCART=0
      DO ICHAR = 1, 8
        IF(OPER(ICHAR:ICHAR).EQ.'X'  .OR.
     &     OPER(ICHAR:ICHAR).EQ.'Y'  .OR.
     &     OPER(ICHAR:ICHAR).EQ.'Z'      ) THEN
             NCART = NCART + 1
        ELSE
*. End of cartesian part
          GOTO 1001
        END IF
      END DO
 1001 CONTINUE
*
      ISYM = 1
      DO ICART = 1, NCART
        IF(OPER(ICART:ICART).EQ.'X') ISYM = MULTD2H(ISYM,IXYZ(1))
        IF(OPER(ICART:ICART).EQ.'Y') ISYM = MULTD2H(ISYM,IXYZ(2))
        IF(OPER(ICART:ICART).EQ.'Z') ISYM = MULTD2H(ISYM,IXYZ(3))
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,'(A,A,A,I1)') ' Label ', OPER, ' has symmetry ', ISYM
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE SYMCOM(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1 I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* IOBJ = 1 : I1,I2 are strings I12 determinant
* ( Other things can follow )
* IOBJ = 2 : I1,I2,I3 are externals
* IOBJ = 3 : I1 is an external, I2,I3 are dets
* IOBJ = 4 : I1 is orbital, I2 is string,l, I12 is string
* IOBJ = 5 : I1 is single excitation, I2 is string,l, I12 is string
* IOBJ = 6 : I1 is orbital, I2 is Orbital I12 is single excitation
*
* If obtained symmetry I1 or I2 is outside bounds,
* zero is returned.
*
* Jeppe Olsen , Spring of 1991
*
* ================
*. Driver routine
* ================
*.LUCINP (PNTGRP is used )
C     INTEGER PNTGRP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*
      IF(PNTGRP.EQ.1) THEN
        CALL SYMCM1(ITASK,IOBJ,I1,I2,I12)
      ELSE IF(PNTGRP.GE.2.AND.PNTGRP.LE.4) THEN
        CALL SYMCM2(ITASK,IOBJ,I1,I2,I12)
      ELSE
        WRITE(6,*) ' PNTGRP parameter out of bounds ', PNTGRP
        WRITE(6,*) ' Enforced stop in SYMCOM '
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE SYMINF(IPRNT)
*
* Information about number of symmetries
*
* Input : /LUCINP/,/ORBINP
* Output : /CSM/,/CSMPRO/
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     INTEGER PNTGRP,EXTSPC
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*

*. Output
* NSMSX : number of symmetries of single excitations
* NSMDX : Number of symmetries of double excitations
* NSMST : Number of symmetries of strings
* NSMCI : NUmber of symmetries of CI spaces
* ITSSX : Total symmetrix single excitation
* ITSDX : Total symmetrix double excitation
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
* ADASX : symmetry of orbs i and i => symmetry of a+iaj
* ASXAD : symmetry of orb j and excit a+iaj => symmetry of i
* ADSXA : symmetry of orb i and excit a+iaj => symmetry of j
*
* SXSXDX : Symmetry of two single excitations
*          => symmetry of double  excitation
* SXDXSX : Symmetry of single excitation and double excitation
*          => symmetry of single  excitation

*.
      IF(PNTGRP.EQ.1) THEN
* =====
* D 2 h
* =====
        CALL ZSYM1(NIRREP,IPRNT)
      ELSE IF(PNTGRP.EQ.2) THEN
* ========
* C inf V
* ========
C            ZNONAB(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB(0,MAXML,NSMOB,IPRNT)
        CALL ZSYM2(IPRNT)
      ELSE IF(PNTGRP.EQ.3.OR.PNTGRP.EQ.4) THEN
* ===========
* D inf H O3
* ===========
C            ZNONAB(INVCEN,MAXMLO,NSMOB)
        CALL ZNONAB(1,MAXML,NSMOB,IPRNT)
        CALL ZSYM2(IPRNT)
      ELSE
        WRITE(6,*) ' You are to early , sorry '
        WRITE(6,*) ' Illegal PNTGRP in SYMINF ',PNTGRP
        Call Abend1( 11 )
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
     &                   KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK)
*
* Construct information about batch and block structure of CI space
* defined by ISPC,ISM,IATP,IBTP.
*
* Output is given in the form of pointers to vectors in WORK
* where the info is stored :
*
* KPCLBT : Length of each Batch ( in blocks)
* KPCLEBT : Length of each Batch ( in elements)
* KPCI1BT : Length of each block
* KPCIBT  : Info on each block
* KPCBLTP : BLock type for each symmetry
*
* NBATCH : Number of batches
* NBLOCK : Number of blocks
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KPCLBT, KPCLEBT, KPCI1BT, KPCIBT, KPCBLTP
      INTEGER*8 KLCIOIO, KSVST
      ! for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "cstate.inc"
#include "csm.inc"
#include "strbas.inc"
#include "crun.inc"
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================== '
        WRITE(6,*) ' Output from Z_BLKFO '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Pointers to output arrays
      NTTS = MXNTTS
      CALL MEMMAN(KPCLBT ,MXNTTS,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KPCLEBT,MXNTTS,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KPCI1BT,MXNTTS,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KPCIBT ,8*MXNTTS,'ADDL  ',1,'CIBT  ')
      CALL MEMMAN(KPCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*.    ^ These should be preserved after exit so put mark for flushing here
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'Z_BLKF')
*. Info needed for generation of block info
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ISPC,WORK(KLCIOIO))
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KPCBLTP),WORK(KSVST))
*
      IF(NTEST.GE.10) THEN
*. Allowed length of each batch
        WRITE(6,*) ' LBLOCK = ', LBLOCK
      END IF
*
*. Batches  of C vector
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KPCBLTP),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &              ISMOST(1,ISM),
     &              NBATCH,WORK(KPCLBT),WORK(KPCLEBT),
     &              WORK(KPCI1BT),WORK(KPCIBT),0,ITTSS_ORD)
*. Number of BLOCKS
      NBLOCK = IFRMR(WORK(KPCI1BT),1,NBATCH)
     &       + IFRMR(WORK(KPCLBT),1,NBATCH) - 1
      IF(NTEST.GE.1) THEN
         WRITE(6,*) ' Number of batches', NBATCH
         WRITE(6,*) ' Number of blocks ', NBLOCK
      END IF
*. Length of each block
      CALL EXTRROW(WORK(KPCIBT),8,8,NBLOCK,WORK(KPCI1BT))
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'Z_BLKF')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
* Codes for general symmetry handling
*
*                - ZSTINF : generate /STINF/ info on strings and mapping
*                - MEMSTR : allocates memory for string information
*                - WEIGHT : Weights for strings
*                - NSTRSO : Number of strings per sym and class
*                - ZBASE  : offset arrays for strings
*                - ZSMCL  : symmetry and class for each string
*                - GENSTR : Generate strings ordered by sym and class
*                - MEMEXT : Memory for external blocks
*
      SUBROUTINE ZBLTP(ISMOST,MAXSYM,IDC,ICBLTP,IMMLST)
*
      implicit real*8 (A-H,O-Z)
*
* Generate vector ICBLTP giving type of each block
*
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
      DIMENSION ISMOST(*),IMMLST(*)
*. Output
      DIMENSION ICBLTP(*)
*
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF (IBSYM.EQ.0) GOTO 100
        IF(((IDC.EQ.2.OR.IDC.EQ.4).AND.(IBSYM.GT.IASYM))
     &                    .OR.
     &       (IDC.EQ.3.AND.IMMLST(IASYM).GT.IASYM)) THEN
*.Symmetry block excluded
          ICBLTP(IASYM) = 0
        ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM)
     &                   .OR.
     &                IDC.EQ.1
     &                   .OR.
     &          (IDC.EQ.3.AND.IASYM.GE.IMMLST(IASYM))) THEN
*.Complete symmetry block included
          ICBLTP(IASYM) = 1
        ELSE
*.Lower half  symmetry block included
          ICBLTP(IASYM) = 2
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Block type of symmetry blocks '
         CALL IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
*
* Codes for general symmetry handling
*
*                - ZSTINF : generate /STINF/ info on strings and mapping
*                - MEMSTR : allocates memory for string information
*                - WEIGHT : Weights for strings
*                - NSTRSO : Number of strings per sym and class
*                - ZBASE  : offset arrays for strings
*                - ZSMCL  : symmetry and class for each string
*                - GENSTR : Generate strings ordered by sym and class
*                - MEMEXT : Memory for external blocks
*
      SUBROUTINE ZBLTP_IDC(ISMOST,MAXSYM,IDC,ICBLTP)
*
*  IMMLST removed. Only needed if IDC = 3 or 4
*   Timo Fleig
*
      implicit real*8 (A-H,O-Z)
*
* Generate vector ICBLTP giving type of each block
*
*
* ICBLTP gives type of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
*. Input
      DIMENSION ISMOST(*)
*. Output
      DIMENSION ICBLTP(*)
*
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF (IBSYM.EQ.0) GOTO 100
        IF(((IDC.EQ.2.OR.IDC.EQ.4).AND.(IBSYM.GT.IASYM))
     &                    .OR.
     &       (IDC.EQ.3)) THEN
*.Symmetry block excluded
          ICBLTP(IASYM) = 0
        ELSE IF((IDC.EQ.2.AND.IASYM.GT.IBSYM)
     &                   .OR.
     &                IDC.EQ.1
     &                   .OR.
     &          (IDC.EQ.3)) THEN
*.Complete symmetry block included
          ICBLTP(IASYM) = 1
        ELSE
*.Lower half  symmetry block included
          ICBLTP(IASYM) = 2
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Block type of symmetry blocks '
         CALL IWRTMA(ICBLTP,1,MAXSYM,1,MAXSYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE ZNELFSPGP(NTESTG)
*
* Generate for each supergroup the number of electrons in each active
* orbital space and store in NELFSPGP
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cgas.inc"
*. Input and Output ( NELFSPGP(MXPNGAS,MXPSTT) )
#include "gasstr.inc"
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      DO ITP = 1, NSTTP
        NSPGP = NSPGPFTP(ITP)
        IBSPGP = IBSPGPFTP(ITP)
        DO ISPGP = IBSPGP,IBSPGP + NSPGP - 1
          DO IGAS = 1, NGAS
            NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Distribution of electrons in Active spaces '
        DO ITP = 1, NSTTP
          WRITE(6,*) ' String type ', ITP
          WRITE(6,*) ' Row : active space, Column: supergroup '
          NSPGP = NSPGPFTP(ITP)
          IBSPGP = IBSPGPFTP(ITP)
          CALL IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE ZNONAB(INVCEN,MAXMLO,NSMOB,IPRNT)
*
*
* ============================
* Set up common block /NONAB/
* ============================
*
*========
* Input :
*========
*      INVCNT :inversion center is present(1), absent(0)
*      MAXMLO : Largest ML value of orbitals
*      NSMOB  : Number of symmetries of orbitals
*      Contents of common block /STRINP/,/ORBINP/
*=========
* output :
*=========
*======================
* Orbital Information
*======================
*      NORASM : Number of orbitals per abelian symmetry
*      MNMLOB : Smallest ML of orbitals
*      MXMLOB : largest ML of orbitals
*      NMLOB  : number of ML values for orbitals
*
*======================
* String Information
*======================
*      MNMLST : smallest ML of any strings
*      MXMLST : largest ML of any strings
*      NMLST  : Number of ML values of strings
*      NSMST  : Number of symmetries of strings
*
*==============================
* Single excitation Information
*==============================
*      MNMLSX : SMALLEST ML OF SINGLE EXCITATION
*      MXMLSX : LARGEST ML OF SINGLE EXCITATIONS
*      NMLSX  : NUMBER OF ML VALUES FOR SINGLE EXCITATIONS
*      NSMSX  : NUMBER OF SYMMETRIES FOR SINGLE EXCITATIONS
*=============================================
* External configurations (upto 4 electrons )
*=============================================
*      MNMLXT : SMALLEST ML OF External configurations
*      MXMLSX : LARGEST ML OF external configurations
*      NMLXT  : NUMBER OF ML VALUES FOR ext. configurations
*      NSMXT  : NUMBER OF SYMMETRIES FOR ext. configurations
*
* =============
* General input
* =============
*
*./ORBINP/
C     PARAMETER(MXPIRR = 20)
C     PARAMETER ( MXPOBS = 20 )
C     PARAMETER (MXPR4T = 10 )
C     PARAMETER(MXPORB=500)
C     INTEGER ORBSM
C     COMMON/ORBINP/NINOB,NACOB,NDEOB,NOCOB,NTOOB,
C    &              NORB0,NORB1,NORB2,NORB3,NORB4,
C    &              NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),
C    &              NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),
C    &              NR4OBS(MXPR4T,MXPOBS),NACOBS(MXPOBS),NOCOBS(MXPOBS),
C    &              NTOOBS(MXPOBS),NDEOBS(MXPOBS),NRS4TO(MXPR4T),
C    &              IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),
C    &              ITPFSO(MXPORB),IBSO(MXPOBS),
C    &              NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB)
#include "mxpdim.inc"
#include "orbinp.inc"
*./CSM/
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*./STRINP/
C     PARAMETER(MXPSTT=40)
C     COMMON/STRINP/NSTTYP,MNRS1(MXPSTT),MXRS1(MXPSTT),
C    &              MNRS3(MXPSTT),MXRS3(MXPSTT),NELEC(MXPSTT),
C    &              IZORR(MXPSTT),IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),
C    &              NZSTTP,NRSTTP,ISTTP(MXPSTT)
#include "strinp.inc"
* =======
*. Output
* =======
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONAB/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*
      NTEST = 0
      NTEST = MAX(IPRNT,NTEST)
*. Inversion symmetry
      IF( INVCEN .EQ. 0 ) THEN
        INVCNT = .FALSE.
        NIG = 1
      ELSE
        INVCNT = .TRUE.
        NIG = 2
      END IF
*
** 1 : Information about orbitals
*
      MXMLOB = MAXMLO
      MNMLOB =-MAXMLO
      NMLOB = MXMLOB - MNMLOB + 1
*. Number of orbitals per symmetry
      DO 10 ISYM = 1, NSMOB
        NORASM(ISYM) = IFREQ(ISMFTO,ISYM,NACOB)
   10 CONTINUE
      IF( NTEST.GE. 2 ) THEN
        WRITE(6,*) ' NORASM '
        CALL IWRTMA(NORASM,1,NSMOB,1,NSMOB)
        WRITE(6,*) ' MNMLOB,MXMLOB ',MNMLOB,MXMLOB
        WRITE(6,*) ' NMLOB, NSMOB ',NMLOB,NSMOB
      END IF
*
**  2. Information about strings
*
      MXMLST = 0
      MNMLST = 0
      DO 50 ITYPE = 1, NSTTYP
        IEL = NELEC(ITYPE)
*
        MXMLTP = 0
        DO 40 IML = MXMLOB,MNMLOB,-1
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MXMLTP = MXMLTP + IEL2*IML
          IEL = IEL - IEL2
   40   CONTINUE
        MXMLST = MAX(MXMLST,MXMLTP)
*
        MNMLTP = 0
        IEL = NELEC(ITYPE)
        DO 45 IML = MNMLOB,MXMLOB
          IORB = NORASM(IML-MNMLOB+1)
          IF(INVCNT) IORB = IORB + NORASM(NMLOB+IML-MNMLOB+1)
          IEL2 = MIN(IORB,IEL)
          MNMLTP = MNMLTP + IEL2*IML
          IEL = IEL - IEL2
   45   CONTINUE
        MNMLST = MIN(MNMLST,MNMLTP)
   50 CONTINUE
*
      NMLST  = MXMLST - MNMLST + 1
      NSMST  = NIG * NMLST
*
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' MXMLST,MNMLST,NSMST'
        WRITE(6,*)   MXMLST,MNMLST,NSMST
      END IF
*
** 3. Information about single excitations
*
      MNMLSX = MNMLOB - MXMLOB
      MXMLSX = MXMLOB - MNMLOB
      NMLSX  = MXMLSX - MNMLSX +1
      NSMSX  = NIG * NMLSX

      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLSX,NSMSX,MNMLSX ',NMLSX,NSMSX,MNMLSX
      END IF
*
** 4 : External configurations(double excitations)
*
      MXMLDX = 4*MAXMLO
      MNMLDX = -4*MAXMLO
      NMLDX  = MXMLDX - MNMLDX + 1
      NSMDX  = NIG * NMLDX
      IF( NTEST .GE.2 ) THEN
        WRITE(6,*) ' NMLDX,NSMDX,MNMLDX ',NMLDX,NSMDX,MNMLDX
      END IF
*
** 5 : Determinants
*
      MXMLCI =  2*MXMLST + MXMLDX
      MNMLCI = - MXMLCI
      NMLCI = 2 * MXMLCI + 1
      NSMCI = NIG * NMLCI
*
*.6 Total symmetrix single excitation and external
*
      ITSSX = 0 - MNMLSX + 1
      ITSDX = 0 - MNMLDX + 1

      IF ( NTEST .GE. 1 ) THEN
        WRITE(6,*)
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of orbitals     .. ', NSMOB
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of strings      .. ', NSMST
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of single excit. . ', NSMSX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of double excit. . ', NSMDX
        WRITE(6,'(A,I4)')
     &  '  Number of symmetries of determinants .. ', NSMCI
        WRITE(6,*)
*
        WRITE(6,*) ' Total symmetric single excitation .. ',ITSSX
        WRITE(6,*) ' Total symmetric double excitation .. ',ITSDX
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE ZSYM1(NIRREP,IPRNT)
*
* Number of symmetries for d2h
* Symmetry connecting arrays
* ( trivial, written for compatibility with higher point groups)
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
C     PARAMETER ( MXPOBS = 20 )
#include "mxpdim.inc"
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
      NSMSX = NIRREP
      NSMDX = NIRREP
      NSMST = NIRREP
      NSMCI = NIRREP
      NSMXT = NIRREP
      ITSSX = 1
      ITSDX = 1
      ITSXT = 1
*
C     COPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
      CALL ICPMT2(SYMPRO,ADASX,8,8,MXPOBS,MXPOBS,1)
      CALL ICPMT2(SYMPRO,ADSXA,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,ASXAD,8,8,MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXSXDX,8,8,2*MXPOBS,2*MXPOBS,1)
      CALL ICPMT2(SYMPRO,SXDXSX,8,8,2*MXPOBS,4*MXPOBS,1)
*
      RETURN
      END
***********************************************************************
*                                                                     *
      SUBROUTINE ZSYM2(IPRNT)
*
* Symmetry connecting arrays
*
* ======
*. Input
* ======
*
*./LUCINP/
C     INTEGER PNTGRP,CITYP,EXTSPC
C     PARAMETER(MXPIRR = 20)
C     PARAMETER(MXPOBS = 20 )
C     PARAMETE
C     COMMON/LUCINP/PNTGRP,NIRREP,NSMCMP,MAXML,MAXL,
C    &              INTSPC,EXTSPC,NRSSH(MXPIRR,3),
C    &              MNRS1R,MXRS1R,MNRS3R,MXRS3R,NACTEL,
C    &              NSMOB,NRS0SH(1,MXPIRR),NRS4SH(MXPR4T,MXPIRR),
C    &              MXR4TP, MXHR0,MXER4,
C    &              NINASH(MXPIRR),
C    &              INTXCI,NDELSH(MXPIRR),MNRS10,MXRS30
#include "mxpdim.inc"
#include "lucinp.inc"
*./NONAB/
      LOGICAL INVCNT
      COMMON/NONAB/ INVCNT,NIG,NORASM(MXPOBS),
     &              MNMLOB,MXMLOB,NMLOB,
     &              MXMLST,MNMLST,NMLST,
     &              NMLSX ,MNMLSX,MXMLSX,
     &              MNMLCI,MXMLCI,NMLCI,
     &              MXMLDX,MNMLDX,NMLDX
*./CSM
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.inc"
*
* ======
*.Output
* ======
*
*./CSMPRD/
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*
**. ADASX,ASXAD,ADSXA
      CALL ISETVC(ADASX,0,MXPOBS**2)
      CALL ISETVC(ASXAD,0,2*MXPOBS**2)
      CALL ISETVC(ADSXA,0,2*MXPOBS**2)
*
      DO 100 ISM = 1, NSMOB
C       MLSM(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM(IML,IPARI,ISM,'OB',2)
        DO 90 JSM = 1, NSMOB
          CALL MLSM(JML,JPARI,JSM,'OB',2)
*.a+ i a j symmetry
          IJML = IML - JML
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJSM = (IJPARI-1)*NMLSX + IJML - MNMLSX + 1
          ADASX(ISM,JSM) = IJSM
          ASXAD(JSM,IJSM) = ISM
          ADSXA(ISM,IJSM) = JSM
   90   CONTINUE
  100 CONTINUE
*.SXSXDX,SXDXSX
      DO 200 ISX = 1, NSMSX
C       MLSM(IML,IPARI,ISM,TYPE,IWAY)
        CALL MLSM(IML,IPARI,ISX,'SX',2)
        DO 190 JSX = 1, NSMSX
          CALL MLSM(JML,JPARI,JSX,'SX',2)
          IF((IPARI.EQ.1.AND.JPARI.EQ.1).OR.
     &       (IPARI.EQ.2.AND.JPARI.EQ.2)) THEN
            IJPARI = 1
          ELSE
            IJPARI = 2
          END IF
          IJML = IML + JML
          IJSM = (IJPARI-1)*NMLDX+IJML - MNMLDX + 1
          SXSXDX(ISX,JSX) = IJSM
          SXDXSX(ISX,IJSM) = JSX
  190   CONTINUE
  200 CONTINUE
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ADASX '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADASX,NSMOB,NSMOB,MXPOBS,MXPOBS)
         WRITE(6,*) ' ASXAD '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ASXAD,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' ADSXA '
         WRITE(6,*) ' ===== '
         CALL IWRTMA(ADSXA,NSMOB,NSMSX,MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXSXDX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXSXDX,NSMSX,NSMSX,2*MXPOBS,2*MXPOBS)
         WRITE(6,*) ' SXDXSX'
         WRITE(6,*) ' ======'
         CALL IWRTMA(SXDXSX,NSMSX,NSMDX,2*MXPOBS,4*MXPOBS)
      END IF
*
      RETURN
      END
