C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
C  /* Deck r12aux */  
      SUBROUTINE R12AUX(WORK,LWORK)
#include <implicit.h>
#include <priunit.h>   
C
      PARAMETER (D0 = 0D0, D1 = 1D0)
      DIMENSION WORK(LWORK)
#include <iratdef.h>
#include <dummy.h>
#include <thrldp.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <gnrinf.h>
#include <r12int.h>
      LOGICAL DELEMO, FOUND, LDUMMY
      INTEGER NBASF(8)
#include <memint.h>
C
      CALL QENTER('R12AUX')
      IF (.NOT. LMULBS) THEN
         DO ISYM = 1, NSYM
            MBAS1(ISYM) = NBAS(ISYM)
            MBAS2(ISYM) = 0
         END DO
      END IF
      NMOCCT = 2 * NOCCT * NOCCT
      CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAO,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFVAO,0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDCAO,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDVAO,0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAUX,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,MAX(NNBAST,NNBASX),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMO,MAX(7*NMOCCT,NNBAST),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAUX,NNBAST,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFQQ,NOCCT**4,WORK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK,LWRK,WORK,KFREE,LFREE)
C
C     ***** READ MOLECULAR ORBITALS *****
C
      JRDMO = 9
      CALL READMO(WORK(KCMO),JRDMO,WORK(KFREE),LFREE)
C
C     ***** INITIALIZE AUXILIARY BASIS SET *****               
C
      ISTC = KAUX
      DO ISYM = 1, NSYM
         NORBI = MBAS1(ISYM)
         NBASI = NBAS(ISYM)
         NAUXI = MBAS2(ISYM)
         NORB2(ISYM) = NAUXI
         IF (NAUXI .NE. 0) THEN
           IJ = ISTC
           DO J = 1, NAUXI
             DO I = 1, NBASI
               IF (I - NORBI .EQ. J) THEN
                 WORK(IJ) = D1
               ELSE
                 WORK(IJ) = D0
               END IF 
               IJ = IJ + 1
             END DO
           END DO
         END IF
         ISTC = ISTC + NAUXI*NBASI  
      END DO
C
C     ***** ORTHOGONALIZE AUXILIARY BASIS SET *****
C
      IF (R12ECO) CALL AUXORT(WORK(KCMO),WORK(KAUX),WORK(KWRK),LWRK)
      DELEMO = .FALSE.
      LAUXBS = .TRUE.
      CALL DELMO(WORK(KAUX),WORK(KWRK),LWRK,THROVL,DELEMO)
      LAUXBS = .FALSE.
C
C     **** CONSTRUCT FOCK MATRIX IN AUXILIARY BASIS ****
C
      CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDCAO),WORK(KDVAO),
     *               WORK(KCMO),DUMMY,WORK(KWRK),LWRK)
      R12TRA = .TRUE.
      NORXR = NORXR .OR. R12HYB
      IF (R12ECO) THEN
         MBSMAX = 6 
      ELSE
         IF (R12NOB .OR. NORXR) THEN
            MBSMAX = 5
         ELSE
            MBSMAX = 6
         END IF
      END IF
      CALL FCKMAO(.TRUE.,DUMMY,EMCMY,WORK(KFCAO),WORK(KFVAO),
     *     WORK(KDCAO),WORK(KDVAO),WORK(KCMO),WORK(KWRK),LWRK)
      WRITE(LUPRI,'(/A,I1,A)') ' Computation of exchange matrix'//
     *                         ' done [',MBSMAX,']'
      MBSMAX = 4
      R12TRA = .FALSE.
      CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KDCAO),1)
      CALL DGETSP(NBAST,WORK(KDCAO),WORK(KFCAO))
      IF (NSYM .GT. 1) CALL PKSYM1(WORK(KFCAO),WORK(KFCAO),NBAS,NSYM,2)
      CALL UTHUB(WORK(KFCAO),WORK(KFAUX),WORK(KAUX),WORK(KDCAO),NSYM,
     *           NBAS,NORB2)   
C
C     **** FORM CANONICAL AUXILIARY ORBITALS ****
C
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXBAS','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,111) NSYM
      ISTF = KFAUX
      ISTC = KAUX
      WRITE(LUPRI,'(/A/A)') ' AUXILIARY BASIS LOW EIGENVALUES',
     *                      ' ISYM         EIGENVALUE'
      DO ISYM = 1, NSYM
         NORBI  = NORB1(ISYM)
         NAUXI = NORB2(ISYM)
         NBASI = NBAS(ISYM)
         WRITE (LUMULB,111) ISYM,NORBI,NAUXI,NBASI
         IF (NAUXI .NE. 0) THEN         
            CALL JACO(WORK(ISTF),WORK(ISTC),NAUXI,NAUXI,NBASI,
     *                WORK(KWRK),WORK(KWRK+NAUXI))
            DO I=1,NAUXI
               WORK(KWRK + I - 1) = WORK(ISTF + I*(I+1)/2 - 1)
            END DO
            CALL ORDER(WORK(ISTC),WORK(KWRK),NAUXI,NBASI)
            WRITE(LUPRI,402) ISYM, WORK(KWRK)
  402       FORMAT(I4,F20.12)
            IJ = ISTC - 1
            DO I=1,NAUXI
               WRITE(LUMULB,111) I
               WRITE(LUMULB,112) WORK(KWRK + I - 1)
               WRITE(LUMULB,112) (WORK(IJ + J), J = 1, NBASI)  
               IJ = IJ + NBASI
            END DO
         END IF
         ISTF   = ISTF + NAUXI*(NAUXI+1)/2
         ISTC   = ISTC + NAUXI*NBASI  
      END DO
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM FOCK MATRIX IN ORTHOGONAL BASIS ****
C
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KSMAT),1)
      ISMO = KFCAO
      ISMO1 = KCMO
      ISMO2 = KAUX
      NNBASF = 0
      DO ISYM = 1, NSYM
         NORBI = NORB1(ISYM)
         NAUXI = NORB2(ISYM)
         NBASI = NBAS(ISYM) 
         CALL DCOPY(NORBI*NBASI,WORK(ISMO1),1,WORK(ISMO),1)
         ISMO = ISMO + NORBI*NBASI
         CALL DCOPY(NAUXI*NBASI,WORK(ISMO2),1,WORK(ISMO),1)
         ISMO = ISMO + NAUXI*NBASI
         ISMO1 = ISMO1 + NORBI*NBASI   
         ISMO2 = ISMO2 + NAUXI*NBASI
         NBASF(ISYM) = NORBI + NAUXI
         NNBASF = NNBASF + NBASF(ISYM)*(NBASF(ISYM)+1)/2
      END DO
      CALL UTHUB(WORK(KSMAT),WORK(KSMO),WORK(KFCAO),WORK(KDCAO),NSYM,
     *           NBAS,NBASF)    
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM OVERLAP IN ORTHOGONAL BASIS ****
C
      CALL GPOPEN(LUMULB,'AUXOVL','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      FOUND = .TRUE.
      CALL RDONEL('OVERLAP ',FOUND,WORK(KSMAT),NNBAST)
      CALL UTHUB(WORK(KSMAT),WORK(KSMO),WORK(KFCAO),WORK(KDCAO),NSYM,
     *           NBAS,NBASF)    
      WRITE(LUMULB,112) (WORK(KSMO + J), J = 0, NNBASF - 1) 
      CALL GPCLOSE(LUMULB,'KEEP')
C
C     **** FORM <ij|r12**2|kl> IN ORTHOGONAL BASIS ****
C
      ISMO = KFCAO
      ISMO1 = KCMO
      CALL DZERO(WORK(ISMO),NOCCT*NBAST)
      DO ISYM = 1, NSYM
         NOCCI = NOCC(ISYM)
         NORBI = NORB1(ISYM)
         NBASI = NBAS(ISYM) 
         DO INOC = 1, NOCCI 
           CALL DCOPY(NBASI,WORK(ISMO1),1,WORK(ISMO),1)
           ISMO = ISMO + NBAST
           ISMO1 = ISMO1 + NBASI
         END DO
         ISMO1 = ISMO1 + (NORBI - NOCCI) * NBASI   
         ISMO = ISMO + NBASI
      END DO 
      CALL R12QQ(WORK(KFQQ),WORK(KSMO+6*NMOCCT),WORK(KSMO+5*NMOCCT),
     *           WORK(KSMO+4*NMOCCT),WORK(KSMO+3*NMOCCT),
     *           WORK(KSMO+2*NMOCCT),WORK(KSMO+NMOCCT),WORK(KSMO), 
     *           WORK(KFCAO),WORK(KFCAO),WORK(KWRK),LWRK,
     *           NBAST,NOCCT)
      CALL GPOPEN(LUMULB,'AUXQ12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUMMY)
      WRITE(LUMULB,112) (WORK(KFQQ + J), J = 0, NOCCT**4 - 1)
      CALL GPCLOSE(LUMULB,'KEEP')
      CALL MEMREL('R12AUX',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('R12AUX')
      RETURN
  111 FORMAT(10I5)
  112 FORMAT(4E30.20)
      END
C  /* Deck dsori*/  
      SUBROUTINE DSORI(WORK,LWORK,UU,ZZ,NBAST,NOCCT)
#include <implicit.h>
#include <priunit.h>
#include <r12int.h>
      DIMENSION WORK(LWORK),UU(*),ZZ(*)                    
      LOGICAL FOUND
      KDSO = 1 + NBAST*NBAST
      KDRI = KDSO + NBAST*NBAST
      KEF1 = KDRI + NBAST*NBAST  
      KEF2 = KEF1 + NBAST*NBAST  
      KEND = KEF2 + NBAST*NBAST*2
      LWRK = LWORK - KEND
C
      CALL RDPROP('DSO 0603',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KDSO),WORK(KEND),LWRK,NBAST,NBAST)
      CALL AROUND('DSO 0603')
      CALL OUTPUT(WORK(KDSO),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
C
      CALL DZERO(WORK(KDRI),NBAST*NBAST)
C
      CALL RDPROP('NEF 05  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 02  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST)
C
      CALL RDPROP('NEF 04  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 01  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST) 
C
      CALL AROUND('DRI 0603')
      CALL OUTPUT(WORK(KDRI),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
C
      CALL RDPROP('DSO 0502',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KDSO),WORK(KEND),LWRK,NBAST,NBAST)
      CALL AROUND('DSO 0502')
      CALL OUTPUT(WORK(KDSO),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
C
      CALL DZERO(WORK(KDRI),NBAST*NBAST)
C
      CALL RDPROP('NEF 06  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 03  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST)
C
      CALL RDPROP('NEF 04  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 01  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST)
C
      CALL AROUND('DRI 0502')
      CALL OUTPUT(WORK(KDRI),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
C
      CALL RDPROP('DSO 0401',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KDSO),WORK(KEND),LWRK,NBAST,NBAST)
      CALL AROUND('DSO 0401')
      CALL OUTPUT(WORK(KDSO),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
C
      CALL DZERO(WORK(KDRI),NBAST*NBAST)
C
      CALL RDPROP('NEF 05  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 02  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST)
C
      CALL RDPROP('NEF 06  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF1),WORK(KEND),LWRK,NBAST,NBAST)
      CALL RDPROP('NEF 03  ',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,WORK(KEF2),WORK(KEND),LWRK,NBAST,NBAST)
      CALL DSORI2(WORK(KEF1),WORK(KEF2),WORK(KDRI),NBAST)
C
      CALL AROUND('DRI 0401')
      CALL OUTPUT(WORK(KDRI),1,MBAS1(1),1,MBAS1(1),NBAST,NBAST,1,LUPRI)
      RETURN
      END
C  /* Deck dsori2*/  
      SUBROUTINE DSORI2(X1,X2,XX,N)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (HALF = 0.5D0)
      DIMENSION X1(N,N), X2(N,N), XX(N,N)
#include <r12int.h>
      DO I = 1, MBAS1(1)
        DO J = 1, MBAS1(1)
          DO K = MBAS1(1) + 1, N
            XX(I,J) = XX(I,J) 
     *              + HALF*(X1(K,I)*X2(K,J) + X2(K,I)*X1(K,J))
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
C  /* Deck auxort */
      SUBROUTINE AUXORT(CMO,AUX,SCRA,LSCRA)
#include <implicit.h>
#include <priunit.h>
      DIMENSION CMO(*), AUX(*), SCRA(*)
      PARAMETER (D1 = 1.0D0)
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <r12int.h>
      KOVLP = 1
      KSMOS = KOVLP + NNBAST
      KSCR1 = KSMOS
      DO ISYM = 1,NSYM
        KSCR1 = KSCR1 + NBAS(ISYM) * (NORB1(ISYM) + NORB2(ISYM))
      END DO
      KSCR2 = KSCR1 + NBAST 
      IF (KSCR2 .GT. LSCRA) CALL ERRWRK('AUXORT',KSCR2+NBASMA,LSCRA)
      CALL RDONEL('OVERLAP ',.TRUE.,SCRA(KOVLP),NNBAST)
      ISMO = KSMOS
      ISMO1 = 1
      ISMO2 = 1
      DO ISYM = 1,NSYM
        ISSYM = KOVLP + IIBAS(ISYM)
        NBASI = NBAS(ISYM)
        NORBI = NORB1(ISYM)
        NAUXI = NORB2(ISYM)
        NORBN = NORBI + NAUXI
        IF (NORBN .NE. 0) THEN        
          CALL DCOPY(NORBI*NBASI,CMO(ISMO1),1,SCRA(ISMO),1)
          CALL DCOPY(NAUXI*NBASI,AUX(ISMO2),1,SCRA(ISMO+NORBI*NBASI),1)
          CALL AUXNRM(SCRA(ISSYM),SCRA(ISMO),NBASI,
     *                NORBI,NAUXI,SCRA(KSCR1))
          CALL DCOPY(NAUXI*NBASI,SCRA(ISMO+NORBI*NBASI),1,AUX(ISMO2),1)
          ISMO1 = ISMO1 + NORBI*NBASI
          ISMO2 = ISMO2 + NAUXI*NBASI  
          ISMO = ISMO + (NORBI + NAUXI) * NBASI
        END IF
      END DO
      RETURN
      END
C  /* Deck auxnrm */
      SUBROUTINE AUXNRM(S,VC,N,M1,M2,W)                
C
C     PROJECT OUT M1 VECTORS
C
#include <implicit.h>
      DIMENSION S(*), VC(*), W(*)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
C
      IF (M1. LE. 0 .OR. M2 .LE. 0) RETURN

      IVCI = M1 * N + 1
      DO I = 1, M2
        CALL MPAPV(N,S,VC(IVCI),W)
        IVCJ = 1
        DO J = 1, M1
          T = DDOT(N,VC(IVCJ),1,W,1)
          DO K = 0, N - 1
             VC(IVCI + K) = VC(IVCI + K) -  T * VC(IVCJ + K)  
          ENDDO
          IVCJ = IVCJ + N
        ENDDO
        IVCI = IVCI + N
      ENDDO   
      RETURN
      END
