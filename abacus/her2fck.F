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
Revision 1.2  2000/05/12 14:01:56  hjj
changed NFACTR to PTFACT to be sure integer truncation implied by NFACTR
will not give erroneous results (potential problem found with ftnchek).

961004 - tsaue INTFCK has undergone major surgery:
               - Inner loop now over integrals, prevents cache misses
               - Redefinition of IFCTYP
               - Possibility of screening on Coulomb/exchange
               - Indices unpacked in buffers
               - Deletetion of all "zero" integrals
941019-hjaaj
changed ISYMDM, IFCTYP to ISYMDM(NDMAT), IFCTYP(NDMAT)
implemented all six IFCTYPes (two was already implemented)
improved efficiency of INTFCK by moving GINT outside and
 moving common factors to SKLFC1
halved no. of  operations for IFCTYP=2 in INTFCK
===========================================================================
#endif
C  /* Deck intfck */
      SUBROUTINE INTFCK(FMAT,AOINT,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,SYMFAC,
     &                  IPRINT,NODV,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,WORK,LWORK)
C*****************************************************************************
C
C  Major surgery Aug 26 1996 by T.Saue
C
C  Written by Trygve Helgaker and Henrik Koch  93-93 ??
C
C  All IFCTYP's added and partly optimized Nov. 1994 by
C  Hans Joergen Aa. Jensen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER(D0 = 0.0D0)
      LOGICAL NODV, SUSCEP
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),
     &          ICORBD(NORBD),
     &          AOINT(NCCINT,NINTYP), FMAT(NBASIS,NBASIS,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), IFCTYP(NDMAT),
     &          DNSBUF(2,NDMAT),DINTSKP(2,4),WORK(LWORK)
#include <twocom.h>
#include <expcom.h>
C
      CALL QENTER('INTFCK')
#include <memint.h>
C
C     Scanning for zero integrals; this is done twice, but
C     saves quite a bit of memory
C
      IF (NINTYP .GT. 1) THEN
         CALL QUIT('Programming error please report. NINTYP .gt. 1')
      END IF
      NBUF = 0
      IF(DOSCRN) THEN
        FCKTOL = MAX((FCKTHR/DNSMAX),1.00D-15)
      ELSE
        FCKTOL = 1.00D-15
      ENDIF
      DO I = 1,NCCINT
        IF(ABS(AOINT(I,1)).GT.FCKTOL) NBUF = NBUF + 1
      ENDDO
      IF(NBUF.EQ.0) THEN
        DINTSKP(1,2) = DINTSKP(1,2) + NCCINT
        DINTSKP(2,2) = DINTSKP(2,2) + NCCINT
        GOTO 999
      ENDIF
      CALL MEMGET('INTE',KIND ,4*NBUF  ,WORK,KFREE,LFREE)
C
Cef begin
      CALL MEMGET('REAL',KAOINT,(NCCINT*NINTYP),WORK,KFREE,LFREE)
Cef end
      CALL INTFC1(FMAT,AOINT,WORK(KAOINT),DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,FCKTOL,SYMFAC,
     &                  IPRINT,NODV,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,WORK(KIND),NBUF)
C
      CALL MEMREL('INTFCK',WORK,KWORK,KWORK,KFREE,LFREE)
 999  CONTINUE
      CALL QEXIT('INTFCK')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck intfc1 */
Cef begin
C Added argument AOINT2 to subroutine INTFC1
Cef end
      SUBROUTINE INTFC1(FMAT,AOINT,AOINT2,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,FCKTOL,SYMFAC,
     &                  IPRINT,NODV,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,IND,NNBUF)
C*****************************************************************************
C
C  Revised Aug 26 1996 by T.Saue
C
C  Written by Trygve Helgaker and Henrik Koch  93-93 ??
C
C  All IFCTYP's added and partly optimized Nov. 1994 by
C  Hans Joergen Aa. Jensen
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <r12int.h>
      PARAMETER (DP25=0.25D0,DP5=0.5D0,D1=1.0D0,D2=2.0D0,D0=0.0D0)
      INTEGER A,B,C,D
      LOGICAL NODV, SUSCEP, NOTEST,LBIT,DOC,DOE
      DIMENSION ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          AOINT(NCCINT), AOINT2(NCCINT),
     &          FMAT(NBASIS,NBASIS,NDMAT,*),
     &          NINDAB(NORBA*NORBB,2), NINDCD(NORBC*NORBD,2),
     &          DMAT(NBASIS,NBASIS,NDMAT), IFCTYP(NDMAT),
     &          DNSBUF(2,NDMAT),IND(NNBUF,4),DINTSKP(2,4)
#include <twocom.h>
#include <expcom.h>
#include <dftcom.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 9) THEN
         CALL HEADER('Output from INTFC1',-1)
         WRITE (LUPRI, '(A,2L5)') ' NODV, SUSCEP',
     *                              NODV, SUSCEP
         WRITE (LUPRI, '(A,4I5)') ' NORB ', NORBA,NORBB,NORBC,NORBD
         WRITE (LUPRI, '(A,2L5)') ' DIAGAB/CD ', DIAGAB,DIAGCD
         WRITE (LUPRI, '(A,3L5)') ' SHAEQB, SHCEQD, SHABAB ',
     *                              SHAEQB, SHCEQD, SHABAB
         WRITE (LUPRI, '(A,4I5)') ' NHKTA', NHKTA,NHKTB,NHKTC,NHKTD
         WRITE (LUPRI, '(A,4I5)') ' KHKTA', KHKTA,KHKTB,KHKTC,KHKTD
         WRITE (LUPRI, '(A)') 'IFCTYP:'
         WRITE (LUPRI, '(10X,10I4)') (IFCTYP(I),I=1,NDMAT)
         WRITE (LUPRI, '(A,F12.6)') ' FCKTOL ', FCKTOL
         WRITE (LUPRI, '(A,F12.6)') ' SYMFAC ', SYMFAC
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals A '
         WRITE (LUPRI, '(20I5)') (ICORBA(I),I=1, NORBA)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals B '
         WRITE (LUPRI, '(20I5)') (ICORBB(I),I=1, NORBB)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals C '
         WRITE (LUPRI, '(20I5)') (ICORBC(I),I=1, NORBC)
         WRITE (LUPRI, '(/A/)') ' Start adresses of orbitals D '
         WRITE (LUPRI, '(20I5)') (ICORBD(I),I=1, NORBD)
      END IF
C
C     ****************************************************
C     ***** Generate indices and permutation factors *****
C     ****************************************************
C
      DOC = DOSCRN.AND.LBIT(ICEFLG,1)
      DOE = DOSCRN.AND.LBIT(ICEFLG,2)
      NOTEST = .NOT.(DIAGAB.OR.DIAGCD.OR.TCONAB.OR.TCONCD)
      SFAC = SYMFAC
      IF (.NOT.SHABAB) SFAC = SFAC+SFAC
      IF (.NOT.SHAEQB) SFAC = SFAC+SFAC
      IF (.NOT.SHCEQD) SFAC = SFAC+SFAC
      NBUF  = 0
      IOFF  = 0
      AOMAX = D0
      IF(SUSCEP) THEN
        IF(NOTEST) THEN
          DO ICOMPA = 1,KHKTA
                        KHKTBB = KHKTB
            DO ICOMPB = 1,KHKTBB
              DO ICOMPC = 1,KHKTC
                            KHKTDD = KHKTD
                DO ICOMPD = 1,KHKTDD
                  DO IORBAB = 1,NORBAB
                    IORBA = 1 + NINDAB(IORBAB,1)/KHKTA
                    IORBB = 1 + NINDAB(IORBAB,2)/KHKTB
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IORBC = 1 + NINDCD(IORBCD,1)/KHKTC
                        IORBD = 1 + NINDCD(IORBCD,2)/KHKTD
                        AOMAX       = MAX(AOMAX,ABS(AINT))
Cef begin                      
                        AOINT2(NBUF) = AINT
Cef end
                        IND(NBUF,1) = A
                        IND(NBUF,2) = B
                        IND(NBUF,3) = ICORBC(IORBC) + ICOMPC
                        IND(NBUF,4) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
Cef begin
          IF(SFAC.NE.D1) CALL DSCAL(NBUF,SFAC,AOINT2,1)
Cef end
        ELSE
          DO ICOMPA = 1,KHKTA
                        KHKTBB = KHKTB
            IF (DIAGAB) KHKTBB = ICOMPA
            DO ICOMPB = 1,KHKTBB
              FACAB = SFAC
              IF (DIAGAB.AND.ICOMPA.NE.ICOMPB) 
     &          FACAB = FACAB+FACAB
              DO ICOMPC = 1,KHKTC
                            KHKTDD = KHKTD
                IF (DIAGCD) KHKTDD = ICOMPC
                DO ICOMPD = 1,KHKTDD
                  FCABCD = FACAB
                  IF(DIAGCD.AND.ICOMPC.NE.ICOMPD) 
     &              FCABCD = FCABCD+FCABCD
                  DO IORBAB = 1,NORBAB
                    FAB = FCABCD
                    IORBA = 1 + NINDAB(IORBAB,1)/KHKTA
                    IORBB = 1 + NINDAB(IORBAB,2)/KHKTB
                    IF (TCONAB.AND.IORBA.NE.IORBB) 
     &                FAB = FAB+FAB
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IORBC = 1 + NINDCD(IORBCD,1)/KHKTC
                        IORBD = 1 + NINDCD(IORBCD,2)/KHKTD
                        FABCD = FAB
                        IF (TCONCD.AND.IORBC.NE.IORBD) 
     &                     FABCD = FABCD+FABCD
                        AOMAX       = MAX(AOMAX,ABS(AINT))
Cef begin                         
                        AOINT2(NBUF) = FABCD*AINT
Cef end
                        IND(NBUF,1) = A
                        IND(NBUF,2) = B
                        IND(NBUF,3) = ICORBC(IORBC) + ICOMPC
                        IND(NBUF,4) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE
        IF(NOTEST) THEN
          DO ICOMPA = 1,KHKTA
                        KHKTBB = KHKTB
            DO ICOMPB = 1,KHKTBB
              DO ICOMPC = 1,KHKTC
                            KHKTDD = KHKTD
                DO ICOMPD = 1,KHKTDD
                  DO IORBAB = 1,NORBAB
                    IORBA = NINDAB(IORBAB,1)
                    IORBB = NINDAB(IORBAB,2)
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IORBC = NINDCD(IORBCD,1)
                        IORBD = NINDCD(IORBCD,2)
                        AOMAX       = MAX(AOMAX,ABS(AINT))
Cef begin                   
                        AOINT2(NBUF) = AINT
Cef end
                        IND(NBUF,1) = A
                        IND(NBUF,2) = B
                        IND(NBUF,3) = ICORBC(IORBC) + ICOMPC
                        IND(NBUF,4) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
Cef end
          IF(SFAC.NE.D1) CALL DSCAL(NBUF,SFAC,AOINT2,1)
Cef end
        ELSE
          DO ICOMPA = 1,KHKTA
                        KHKTBB = KHKTB
            IF (DIAGAB) KHKTBB = ICOMPA
            DO ICOMPB = 1,KHKTBB
              FACAB = SFAC
              IF(DIAGAB.AND.ICOMPA.NE.ICOMPB) 
     &          FACAB = FACAB+FACAB
              DO ICOMPC = 1,KHKTC
                            KHKTDD = KHKTD
                IF (DIAGCD) KHKTDD = ICOMPC
                DO ICOMPD = 1,KHKTDD
                  FCABCD = FACAB
                  IF (DIAGCD.AND.ICOMPC.NE.ICOMPD) 
     &              FCABCD = FCABCD+FCABCD
                  DO IORBAB = 1,NORBAB
                    IORBA = NINDAB(IORBAB,1)
                    IORBB = NINDAB(IORBAB,2)
                    FAB   = FCABCD
                    IF(TCONAB.AND.IORBA.NE.IORBB) 
     &                FAB = FAB+FAB
                    A = ICORBA(IORBA) + ICOMPA
                    B = ICORBB(IORBB) + ICOMPB
                    DO IORBCD = 1, NORBCD
                      AINT  = AOINT(IOFF+IORBCD)
                      IF(ABS(AINT).GT.FCKTOL) THEN
                        NBUF  = NBUF + 1
                        IORBC = NINDCD(IORBCD,1)
                        IORBD = NINDCD(IORBCD,2)
                        FABCD = FAB
                        IF (TCONCD.AND.IORBC.NE.IORBD) 
     &                    FABCD = FABCD+FABCD
                        AOMAX       = MAX(AOMAX,ABS(AINT))
Cef begin                        
                        AOINT2(NBUF) = FABCD*AINT
Cef end
                        IND(NBUF,1) = A
                        IND(NBUF,2) = B
                        IND(NBUF,3) = ICORBC(IORBC) + ICOMPC
                        IND(NBUF,4) = ICORBD(IORBD) + ICOMPD
                      ENDIF
                    ENDDO
                    IOFF = IOFF + NORBCD
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     Statistics
C
      NDEL = NCCINT - NBUF
      DINTSKP(1,2) = DINTSKP(1,2) + NCCINT
      DINTSKP(2,2) = DINTSKP(2,2) + NDEL
      DINTSKP(1,3) = DINTSKP(1,3) + NCM*NBUF
      DINTSKP(1,4) = DINTSKP(1,4) + NEM*NBUF
C
C     **********************************
C     ***** Print section  *************
C     **********************************
C
      IF(IPRINT.GE.15) THEN
        DO INT = 1,NBUF
Cef begin 
          WRITE (LUPRI,'(1X,A,F16.12)') 'GINT,A,B,C,D',AOINT2(INT)
Cef end
        ENDDO
      ENDIF
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C           X = 3 Mixed symmetric and anti-symmetric symmetry blocks
C                 (for relativistic calculations)
C         Y indicates contributions
C           Y = 0 no contribution!
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
C We keep to outer loop over NDMAT, but the actual work is
C transfered to FCKCON, which is also used by eri and sirius.
C
      DO 400 I = 1,NDMAT
        FCM = DNSBUF(1,I)*AOMAX
        FEM = DNSBUF(2,I)*AOMAX
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF (IX.EQ.2) IC = 0
C       ... no Coulomb term for antisymmetric density matrix
        IF (HFXFAC.EQ.D0) IE = 0
C       Screening on Coulomb contributions
        IF(DOC.AND.IC.EQ.1.AND.FCM.LT.FCKTHR) THEN
          IC = 0
          DINTSKP(2,3) = DINTSKP(2,3) + NBUF
        ENDIF
C       Screening on exchange contributions
        IF(DOE.AND.IE.EQ.2.AND.FEM.LT.FCKTHR) THEN
          IE = 0
          DINTSKP(2,4) = DINTSKP(2,4) + NBUF
        ENDIF
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
C
C   FAC is a factor to account for different definiton of
C   integrals in twoint and eri. Here always 1.0
      IF (R12TRA) THEN
C        Compute exchange operator for auxiliary basis set
C        in R12 calculations (WK/UniKA/04-11-2002).
         FAC = -D1
         IY = 2
      ELSE
         FAC = D1
      END IF
Cef begin
      CALL FCKCON(FMAT,DMAT,NDMAT,I,AOINT2,IND,
     &            NCCINT,NBUF,NBUF,IX,IY,FAC)
Cef end
C
 400  CONTINUE
C
      RETURN
C
      END
C  /* Deck dsotao */
      SUBROUTINE DSOTAO(DSO,DAO,NBAST,IREPDM,IPRINT)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
C                                          880418  PRT
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      DIMENSION DSO(NBAST,NBAST), DAO(NBAST,NBAST)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <ibtfun.h>
      IF (IPRINT .GT. 10) CALL HEADER('Subroutine DSOTAO',-1)
C
C     Loop over all irreps in molecule
C
      ISTRA = 1
      CALL DZERO(DAO,NBAST*NBAST)
      DO 100 IREPA = 0, MAXREP
         NORBA = NAOS(IREPA+1)
         DO 200 I = ISTRA,ISTRA + NORBA - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            NHKTA  = NHKT(IA)
            KHKTA  = KHKT(IA)
            MULA   = ISTBAO(IA)
            INDA   = KSTRT(IA) + NA - KHKTA
            DO 300 ISYMA = 0, MAXOPR
            IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
               INDA = INDA + KHKTA
               FACA = PT(IBTAND(ISYMA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
               ISTRB = 1
               DO 400 IREPB = 0, MAXREP
                  NORBB = NAOS(IREPB+1)
                  IF (IBTXOR(IREPA,IREPB).EQ.IREPDM) THEN
                  DO 500 J = ISTRB,ISTRB + NORBB - 1
                     IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                     NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                     NHKTB  = NHKT(IB)
                     KHKTB  = KHKT(IB)
                     MULB   = ISTBAO(IB)
                     INDB   = KSTRT(IB) + NB - KHKTB
                     DO 600 ISYMB = 0, MAXOPR
                     IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                        INDB = INDB + KHKTB
                        FACB = PT(IBTAND(ISYMB,
     &                             IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                        DAO(INDA,INDB) = DAO(INDA,INDB)
     &                                 + FACA*FACB*DSO(I,J)
                     END IF
  600                CONTINUE
  500             CONTINUE
                  END IF
                  ISTRB = ISTRB + NORBB
  400          CONTINUE
            END IF
  300       CONTINUE
  200    CONTINUE
         ISTRA = ISTRA + NORBA
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Total density matrix in SO basis',-1)
         CALL OUTPUT(DSO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('Total density matrix in AO basis',-1)
         CALL OUTPUT(DAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck sklfck */
      SUBROUTINE SKLFCK(FMAT,HESSEE,WORK,LWORK,IPRINT,DIRFCK,DDFOCK,
     &                  EXPECT,PERTUR,NODV,MAXDER,LONDON,NDMAT,ISYMDM,
     &                  IFCTYP,IATOME,SPNORB)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <dummy.h>
#include <mxcent.h>
#include <iratdef.h>
#include <nuclei.h>
#include <inforb.h>
      LOGICAL DIRFCK, DDFOCK, NODV, EXPECT, LONDON, PERTUR, SPNORB
      DIMENSION FMAT(N2BASX,*), WORK(LWORK), HESSEE(*), ISYMDM(*),
     &          IFCTYP(*)
#include <abainf.h>
#include <ccom.h>
#include <shells.h>
#include <symmet.h>
#include <inftap.h>
C
      CALL QENTER('SKLFCK')
      KFREE  = 1
      LFREE  = LWORK
      KFRSAV = KFREE
C
C     Undifferentiated Fock matrix
C     ============================
C
      IF (DIRFCK) THEN
C
         CALL MEMGET('REAL',KSKLTN,N2BASX,WORK,KFREE,LFREE)
         LENMEM = NBAST*(MAXREP + 1)
         CALL MEMGET('REAL',KFAC  ,LENMEM,WORK,KFREE,LFREE)
         LENMEM = KMAX*(MAXREP + 1)*KHK(NHTYP)
         CALL MEMGET('INTE',KINDEX,LENMEM,WORK,KFREE,LFREE)
         LENMEM = NBAST*(MAXREP + 1)
         CALL MEMGET('INTE',KPOINT,LENMEM,WORK,KFREE,LFREE)
C
         IXYZ = 0
         CALL SKLFC1(FMAT,WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &               WORK(KPOINT),DUMMY,DUMMY,DUMMY,IDUMMY,IXYZ,
     &               NDMAT,ISYMDM,IFCTYP,IPRINT)
         CALL MEMREL('SKLFC1 for DIRFCK',WORK,1,KFRSAV,KFREE,LFREE)
      END IF
C
C     Derivative Fock matrices
C     ========================
C
      IF (DDFOCK .AND. .NOT.LONDON) THEN
         IF (PERTUR) THEN
            NMATMX = NUCDEG(IATOME)
         ELSE
            NMATMX = 0
            DO 100 IATOM = 1, NUCIND
               NMATMX = MAX(NMATMX,NUCDEG(IATOM))
  100       CONTINUE
         END IF
C
         KSKLTN = 1
         KFAC   = KSKLTN + 3*NMATMX*N2BASX
         KINDEX = KFAC   + NBAST*(MAXREP + 1)
         KPOINT = KINDEX + (KMAX*(MAXREP + 1)*KHK(NHTYP) + 1)/IRAT
         KNFAC  = KPOINT + (NBAST*(MAXREP + 1) + 1)/IRAT
         KNINDX = KNFAC  + NMATMX*(MAXREP + 1)
         KNPINT = KNINDX + ((MAXREP + 1) + 1)/IRAT
         KLAST  = KNPINT + (NMATMX*(MAXREP + 1) + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('SKLFCK',' ',KLAST,LWORK)
         IOFF = 0
         IFMAT = 0
         DO 200 ITYPE = 1, 2
         IF (ITYPE .EQ. 2 .AND. NODV) GO TO 200
            DO 300 IATOM = 1, NUCIND
            IF (.NOT.PERTUR .OR. IATOM.EQ.IATOME) THEN
               NMAT = NUCDEG(IATOM)
               CALL DCOPY(3*NMAT*N2BASX,FMAT(1,IOFF+1),1,WORK(KSKLTN),1)
               IOFFX = IOFF
               IOFFY = IOFF +   NMAT
               IOFFZ = IOFF + 2*NMAT
               DO 310 I = 1, NMAT
                  IADR = KSKLTN + 3*N2BASX*(I - 1)
                  CALL DCOPY(N2BASX,WORK(IADR),1,
     &                              FMAT(1,IOFFX + I),1)
                  CALL DCOPY(N2BASX,WORK(IADR+N2BASX),1,
     &                              FMAT(1,IOFFY + I),1)
                  CALL DCOPY(N2BASX,WORK(IADR+2*N2BASX),1,
     &                              FMAT(1,IOFFZ + I),1)
  310          CONTINUE
C     We ought also to reorder IFCTYP and ISYMDM, but as
C     all values are 1 and 0, resp., never mind. /941221-hjaaj
               DO 400 IXYZ = 1, 3
                  IADR = IOFF + NMAT*(IXYZ - 1) + 1
                  CALL SKLFC1(FMAT(1,IADR),WORK(KSKLTN),WORK(KFAC),
     &                        WORK(KINDEX),WORK(KPOINT),WORK(KNFAC),
     &                        WORK(KNINDX),WORK(KNPINT),IATOM,IXYZ,NMAT,
     &                        ISYMDM(IADR),IFCTYP(IADR),IPRINT)
                  CALL DSITSP(NBAST,FMAT(1,IADR),WORK(KSKLTN))
C
                  IF (EXPFCK) THEN
                     IF (NODV) THEN
                        IREC = 3*(IATOM - 1) + IXYZ
                     ELSE
                        IREC = 6*(IATOM - 1) + 3*(ITYPE - 1) + IXYZ
                     END IF
                     CALL WRITDX(LUDFCK,IREC,IRAT*NNBASX,WORK(KSKLTN))
                  ELSE
                     IFMAT = IFMAT + 1
                     CALL DCOPY(NNBASX,WORK(KSKLTN),1,FMAT(1,IFMAT),1)
                  END IF
  400          CONTINUE
               IOFF = IOFF + 3*NMAT
            END IF
  300       CONTINUE
  200    CONTINUE
      END IF
C
C     Expectation values
C     ==================
C
      IF (EXPECT) THEN
         KSKLTN = 1
         KFAC   = KSKLTN + MXCOOR**2
         KINDEX = KFAC   + 3*NUCDEP*(MAXREP + 1)
         KPOINT = KINDEX + (3*NUCIND*(MAXREP + 1) + 1)/IRAT
         KCSTRA = KPOINT + (3*NUCDEP*(MAXREP + 1) + 1)/IRAT
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         KLAST  = KSCTRA + MXCOOR*MXCOOR
         IF (KLAST .GT. LWORK) CALL STOPIT('SKLFCK','EXP',KLAST,LWORK)
         CALL SKLEXP(WORK(KSKLTN),HESSEE,WORK(KCSTRA),WORK(KSCTRA),
     &               WORK(KFAC),WORK(KINDEX),WORK(KPOINT),MAXDER,IPRINT)
      END IF
C
C     London Fock matrics
C     ===================
C
      IF ((DDFOCK .AND. LONDON) .OR. SPNORB) THEN
         KSKLTN = 1
         KFAC   = KSKLTN + 3*N2BASX
         KINDEX = KFAC   + NBAST*(MAXREP + 1)
         KPOINT = KINDEX + (KMAX*(MAXREP + 1)*KHK(NHTYP) + 1)/IRAT
         KLAST  = KPOINT + (NBAST*(MAXREP + 1) + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('SKLFCK','LND',KLAST,LWORK)
C
         CALL SKLOND(FMAT(1,1),WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &               WORK(KPOINT),NDMAT,IFCTYP,ISYMDM,SPNORB,IPRINT)
         IF (.NOT.NODV) THEN
            IV=4
            IF (SPNORB) IV=2
            CALL SKLOND(FMAT(1,IV),WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
     &                  WORK(KPOINT),NDMAT,IFCTYP,ISYMDM,SPNORB,IPRINT)
         END IF
C
         NTYPE = 3
         IF (.NOT.NODV) NTYPE = 6
         IF (.NOT.SPNORB) THEN
            DO 510 N = 1, NTYPE
               CALL DGETAP(NBAST,FMAT(1,N),WORK(KSKLTN))
               CALL DCOPY(NNBASX,WORK(KSKLTN),1,FMAT(1,N),1)
  510       CONTINUE
         END IF
      END IF
      CALL QEXIT('SKLFCK')
      RETURN
      END
C  /* Deck sklfc1 */
      SUBROUTINE SKLFC1(FMAT,SKLTON,FACTOR,INDEX,IPOINT,PTFACT,NINDEX,
     &                  NPOINT,IATOM,IXYZ,NMAT,ISYMDM,IFCTYP,IPRINT)
C
C Trygve Helgaker and Henrik Koch 93-94 ??
C 970304 - tsaue New scaling
C
C Implemented all IFCTYP's Oct-Nov 94 Hans Joergen Aa. Jensen
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
      PARAMETER (D1 = 1.0D0,D2I = 0.5D0)
      PARAMETER (D4I = 0.25D0, D8I = 0.125D0)
      INTEGER A, B, AP, BP
      DIMENSION FMAT(NBASIS,NBASIS,*), SKLTON(NBASIS,NBASIS,*),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP),
     &          NINDEX(0:MAXREP), NPOINT(NMAT,0:MAXREP),
     &          PTFACT(NMAT,0:MAXREP), ISYMDM(NMAT), IFCTYP(NMAT)
#include <pincom.h>
#include <shells.h>
#include <ibtfun.h>
C
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Subroutine SKLFC1',-1)
         CALL HEADER("Input skeleton Fock matrix in SKLFC1",1)
         DO 10 N = 1, NMAT
            CALL OUTPUT(FMAT(1,1,N),1,NBASIS,1,NBASIS,
     &                                NBASIS,NBASIS,1,LUPRI)
   10    CONTINUE
      END IF
C
C     Finish skeleton Fock matrices
C     ===================================
C     (a) symmetrize or antisymmetrize, if needed
C     (b) multiply with factors omitted in INTFCK
C         (see comments in INTFCK)
C     /941019-hjaaj
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Final skeleton Fock matrices in SKLFC1",1)
      END IF
      DO 100 IMAT = 1, NMAT
         IY = MOD(IFCTYP(IMAT),10)
         IX = (IFCTYP(IMAT)-IY)/10
         IF(IX.EQ.0) THEN
C        Non-symmetric matrix
           CALL DSCAL(NBASIS*NBASIS,D8I,FMAT(1,1,IMAT),1)
         ELSEIF(IX.EQ.1) THEN
C        Symmetric matrix
            DO 110 I = 1, NBASIS
            DO 110 J = 1, I
               FMATIJ = D4I*(FMAT(I,J,IMAT) + FMAT(J,I,IMAT))
               FMAT(I,J,IMAT) = FMATIJ
               FMAT(J,I,IMAT) = FMATIJ
  110       CONTINUE
         ELSEIF(IX.EQ.2) THEN
C        Anti-symmetric matrix
            DO 120 I = 1, NBASIS
            DO 120 J = 1, I
               FMATIJ = D4I*(FMAT(I,J,IMAT) - FMAT(J,I,IMAT))
               FMAT(I,J,IMAT) =   FMATIJ
               FMAT(J,I,IMAT) = - FMATIJ
  120       CONTINUE
         ENDIF
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,*) 'Skeleton Fock matrix no.',IMAT,
     &         ' out of',NMAT
            CALL OUTPUT(FMAT(1,1,IMAT),1,NBASIS,1,NBASIS,
     &                                NBASIS,NBASIS,1,LUPRI)
         END IF
  100 CONTINUE
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 200 ISHELL = 1, KMAX
            DO 210 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 220 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  220          CONTINUE
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBAO(ISHELL)))
               DO 230 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  230          CONTINUE
            END IF
  210       CONTINUE
  200    CONTINUE
C
         IORB = 0
         DO 300 ISHELL = 1, KMAX
            NHKTA = NHKT(ISHELL)
            DO 310 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 320 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 330 ISYM = 0, MAXREP
                     ISYMPR = IBTXOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  330             CONTINUE
  320          CONTINUE
            END IF
  310       CONTINUE
  300    CONTINUE
C
C        NPOINT and PTFACT
C
         IF (IXYZ .GT. 0) THEN
            IMAT = 0
            DO 215 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               IMAT = IMAT + 1
               NINDEX(ISYMOP) = IMAT
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBNU(IATOM)))
               NINDEX(ISYMOP) = NINDEX(ISYMF)
            END IF
  215       CONTINUE
C
            IMAT = 0
            DO 315 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               IMAT = IMAT + 1
               DO 335 ISYM = 0, MAXREP
                  ISYMPR = IBTXOR(ISYMOP,ISYM)
                  NPOINT(IMAT,ISYM) = NINDEX(ISYMPR)
                  PTFACT(IMAT,ISYM) = PT(IBTAND(ISYM,ISYMAX(IXYZ,1)))
  335          CONTINUE
            END IF
  315       CONTINUE
         END IF
C
C        Construct full matrix from skeleton
C
         IF (IXYZ. EQ. 0) THEN
C           one SKLTON matrix for IXYZ = 0 /961114-hjaaj
            DO 435 N = 1, NMAT
               CALL DCOPY(NBASIS**2,FMAT(1,1,N),1,SKLTON,1)
               DO 400 ISYMOP = 1, MAXREP
                  FP = PT(IBTAND(ISYMOP,ISYMDM(N)))
                  DO 420 B = 1, NBASIS
                     BP = IPOINT(B,ISYMOP)
                     FBFP = FACTOR(B,ISYMOP)*FP
                     DO 410 A = 1, NBASIS
                        AP = IPOINT(A,ISYMOP)
                        FA = FACTOR(A,ISYMOP)
                        FMAT(AP,BP,N) = FMAT(AP,BP,N)
     &                                + FA*FBFP*SKLTON(A,B,1)
  410                CONTINUE
  420             CONTINUE
  400          CONTINUE
  435      CONTINUE
         ELSE
C           NMAT SKLTON matrices for IXYZ .ne. 0  /961114-hjaaj
            CALL DCOPY(NMAT*(NBASIS**2),FMAT,1,SKLTON,1)
            DO 430 ISYMOP = 1, MAXREP
               DO 440 N = 1, NMAT
                  NP = NPOINT(N,ISYMOP)
                  FP = PTFACT(N,ISYMOP)
                  DO 450 A = 1, NBASIS
                     AP = IPOINT(A,ISYMOP)
                     FA = FACTOR(A,ISYMOP)
                     DO 460 B = 1, NBASIS
                        BP = IPOINT(B,ISYMOP)
                        FB = FACTOR(B,ISYMOP)
                        FMAT(AP,BP,NP) = FMAT(AP,BP,NP)
     &                                 + FA*FB*FP*SKLTON(A,B,N)
  460                CONTINUE
  450             CONTINUE
  440          CONTINUE
  430       CONTINUE
         END IF
C
         FAC = D1 / (1.0D0 + MAXREP)
         CALL DSCAL(NMAT*(NBASIS**2),FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrix in SKLFC1",1)
            DO 40 N = 1, NMAT
               CALL OUTPUT(FMAT(1,1,N),1,NBASIS,1,NBASIS,
     &                     NBASIS,NBASIS,1,LUPRI)
   40       CONTINUE
         END IF
C
C        Transformation to symmetry basis
C        ================================
C
         IF (IXYZ .EQ. 0) THEN
C           one SKLTON matrix for IXYZ = 0 /961114-hjaaj
            DO 550 N = 1, NMAT
              CALL DCOPY(NBASIS**2,FMAT(1,1,N),1,SKLTON,1)
              CALL DZERO(FMAT(1,1,N),NBASIS**2)
              ISTRA = 1
              DO 500 IREPA = 0, MAXREP
               NORBA = NAOS(IREPA+1)
               DO 510 I = ISTRA,ISTRA + NORBA - 1
                  IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
                  NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
                  NHKTA  = NHKT(IA)
                  KHKTA  = KHKT(IA)
                  MULA   = ISTBAO(IA)
                  INDA   = KSTRT(IA) + NA - KHKTA
                  DO 520 ISYMA = 0, MAXOPR
                  IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
                     INDA = INDA + KHKTA
                     FACA = PT(IBTAND(ISYMA,
     &                                IBTXOR(IREPA,ISYMAO(NHKTA,NA))))
                        ISTRB = 1
                        DO 525 IREPB = 0, MAXREP
                           NORBB = NAOS(IREPB+1)
                           IF (IBTXOR(IREPA,IREPB) .EQ. ISYMDM(N)) THEN
                           DO 530 J = ISTRB,ISTRB + NORBB - 1
                              IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                              NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                              NHKTB  = NHKT(IB)
                              KHKTB  = KHKT(IB)
                              MULB   = ISTBAO(IB)
                              INDB   = KSTRT(IB) + NB - KHKTB
#include <scalar.h>
                              DO 540 ISYMB = 0, MAXOPR
                              IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                                 INDB = INDB + KHKTB
                                 FACB = PT(IBTAND(ISYMB,
     &                                  IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
                                 FMAT(I,J,N) = FMAT(I,J,N)
     &                                 + FACA*FACB*SKLTON(INDA,INDB,1)
                              END IF
540                           CONTINUE
530                        CONTINUE
                           END IF
                           ISTRB = ISTRB + NORBB
  525                   CONTINUE
                  END IF
520               CONTINUE
510            CONTINUE
               ISTRA = ISTRA + NORBA
500           CONTINUE
  550       CONTINUE
         ELSE
C           NMAT SKLTON matrices for IXYZ .ne. 0  /961114-hjaaj
            CALL DCOPY(NMAT*NBASIS**2,FMAT,1,SKLTON,1)
            CALL DZERO(FMAT,NMAT*NBASIS*NBASIS)
            MULN = ISTBNU(IATOM)
            ISTRA = 1
            DO 600 IREPA = 0, MAXREP
               NORBA = NAOS(IREPA+1)
               DO 610 I = ISTRA,ISTRA + NORBA - 1
                  IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
                  NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
                  NHKTA = NHKT(IA)
                  KHKTA = KHKT(IA)
                  MULA  = ISTBAO(IA)
                  INDA  = KSTRT(IA) + NA - KHKTA
                  IVARA = IBTXOR(IREPA,ISYMAO(NHKTA,NA))
                  DO 620 ISYMA = 0, MAXOPR
                  IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
                     FA    = PT(IBTAND(ISYMA,IVARA))
                     INDA  = INDA + KHKTA
                     ISTRB = 1
                     DO 630 IREPB = 0, MAXREP
                        IVARN=IBTXOR(IBTXOR(IREPA,IREPB),ISYMAX(IXYZ,1))
                        NORBB = NAOS(IREPB+1)
                        IF (IBTAND(MULN,IVARN).EQ.0) THEN
                           DO 640 J = ISTRB,ISTRB + NORBB - 1
                              IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                              NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                              NHKTB = NHKT(IB)
                              KHKTB = KHKT(IB)
                              MULB  = ISTBAO(IB)
                              INDB  = KSTRT(IB) + NB - KHKTB
                              IVARB = IBTXOR(IREPB,ISYMAO(NHKTB,NB))
                              DO 650 ISYMB = 0, MAXOPR
                              IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                                 INDB = INDB + KHKTB
                                 FAB = FA*PT(IBTAND(ISYMB,IVARB))
                                 DO 660 ISYMN = 0, MAXOPR
                                 IF (IBTAND(ISYMN,MULN) .EQ. 0) THEN
                                    INDN = NINDEX(ISYMN)
                                    FAC  = FAB*PT(IBTAND(ISYMN,IVARN))
                                    FMAT(I,J,1) = FMAT(I,J,1)
     &                                      + FAC*SKLTON(INDA,INDB,INDN)
                                 END IF
  660                            CONTINUE
                              END IF
  650                         CONTINUE
  640                      CONTINUE
                           END IF
                        ISTRB = ISTRB + NORBB
  630                CONTINUE
                  END IF
  620             CONTINUE
  610          CONTINUE
               ISTRA = ISTRA + NORBA
  600       CONTINUE
         END IF
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrix in SKLFC1",1)
            DO 700 N = 1, NMAT
               CALL OUTPUT(FMAT(1,1,N),1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &                     1,LUPRI)
  700       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck sklexp */
      SUBROUTINE SKLEXP(SKLTON,HESSEE,CSTRA,SCTRA,FACTOR,INDEX,IPOINT,
     &                  MAXDER,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER A, B, AP, BP
      DIMENSION SKLTON(MXCOOR,MXCOOR), HESSEE(MXCOOR,MXCOOR),
     &          INDEX(NUCIND,3,0:MAXREP), FACTOR(3*NUCDEP,0:MAXREP),
     &          IPOINT(3*NUCDEP,0:MAXREP), CSTRA(*), SCTRA(*)
#include <energy.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL HEADER('Subroutine SKLEXP',-1)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton molecular gradient",1)
         CALL OUTPUT(GRADEE,1,1,1,3*NUCDEP,1,MXCOOR,1,LUPRI)
      END IF
      IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton molecular Hessian",1)
         CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,1,LUPRI)
      END IF
C
      IF (MAXDER .EQ. 2) THEN
         DO 10 I = 1, 3*NUCDEP
         DO 10 J = 1, I - 1
            HESSEE(I,J) = HESSEE(I,J) + HESSEE(J,I)
            HESSEE(J,I) = HESSEE(I,J)
   10    CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("Skeleton molecular Hessian symmetrized",1)
            CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,
     &                  1,LUPRI)
         END IF
      END IF
C
      IF (MAXREP .GT. 0) THEN
         ICOOR = 0
         DO 100 IATOM = 1, NUCIND
            DO 110 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               DO 120 ICART = 1, 3
                  ICOOR = ICOOR + 1
                  INDEX(IATOM,ICART,ISYMOP) = ICOOR
  120          CONTINUE
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBNU(IATOM)))
               DO 130 ICART = 1, 3
                  INDEX(IATOM,ICART,ISYMOP) = INDEX(IATOM,ICART,ISYMF)
  130          CONTINUE
            END IF
  110       CONTINUE
  100    CONTINUE
C
         ICOOR = 0
         DO 200 IATOM = 1, NUCIND
            DO 210 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
               DO 220 ICART = 1, 3
                  ICOOR = ICOOR + 1
                  DO 230 ISYM = 0, MAXREP
                     ISYMPR = IBTXOR(ISYMOP,ISYM)
                     IPOINT(ICOOR,ISYM)=INDEX(IATOM,ICART,ISYMPR)
                     FACTOR(ICOOR,ISYM)=PT(IBTAND(ISYM,ISYMAX(ICART,1)))
  230             CONTINUE
  220          CONTINUE
            END IF
  210       CONTINUE
  200    CONTINUE
C
C        Construct full gradient from skeleton
C
         CALL DCOPY(MXCOOR,GRADEE,1,SKLTON,1)
         DO 300 ISYMOP = 1, MAXREP
            DO 310 A = 1, 3*NUCDEP
               AP = IPOINT(A,ISYMOP)
               GRADEE(AP) = GRADEE(AP) + FACTOR(A,ISYMOP)*SKLTON(A,1)
  310       CONTINUE
  300    CONTINUE
         FAC = D1 / (MAXREP + 1.0D0)
         CALL DSCAL(MXCOOR,   FAC,GRADEE,1)
C
C        Construct full Hessian from skeleton
C
         IF (MAXDER. EQ. 2) THEN
            CALL DCOPY(MXCOOR**2,HESSEE,1,SKLTON,1)
            DO 400 ISYMOP = 1, MAXREP
               DO 410 A = 1, 3*NUCDEP
               DO 410 B = 1, 3*NUCDEP
                  AP  = IPOINT(A,ISYMOP)
                  BP  = IPOINT(B,ISYMOP)
                  FAC = FACTOR(A,ISYMOP)*FACTOR(B,ISYMOP)
                  HESSEE(AP,BP) = HESSEE(AP,BP) + FAC*SKLTON(A,B)
  410          CONTINUE
  400       CONTINUE
            FAC = D1/(MAXREP + D1)
            CALL DSCAL(MXCOOR**2,FAC,HESSEE,1)
         END IF
C
         IF (IPRINT .GT. 5) THEN
           CALL HEADER("Full molecular gradient (non-symmetry basis)",1)
           CALL OUTPUT(GRADEE,1,1,1,3*NUCDEP,1,MXCOOR,1,LUPRI)
         END IF
         IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
           CALL HEADER("Full molecular Hessian (non-symmetry basis)",1)
           CALL OUTPUT(HESSEE,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,
     &                 1,LUPRI)
         END IF
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(MXCOOR,GRADEE,1,SKLTON,1)
         CALL DZERO(GRADEE,MXCOOR)
         ICOOR = 0
         DO 500 IATOM = 1, NUCIND
         DO 500 ICART = 1, 3
            IF (IBTAND(ISTBNU(IATOM),ISYMAX(ICART,1)).EQ.0) THEN
               ICOOR = ICOOR + 1
               GRD = D0
               DO 510 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISTBNU(IATOM),ISYMOP).EQ.0) THEN
                  FAC = PT(IBTAND(ISYMOP,ISYMAX(ICART,1)))
                  GRD = GRD + FAC*SKLTON(INDEX(IATOM,ICART,ISYMOP),1)
               END IF
  510          CONTINUE
               GRADEE(ICOOR) = GRD
            END IF
  500    CONTINUE
C
         IF (MAXDER .EQ. 2) THEN
            CALL DCOPY(MXCOOR**2,HESSEE,1,SKLTON,1)
            CALL DZERO(HESSEE,MXCOOR**2)
            A = 0
            DO 600 IREPA  = 0, MAXREP
            DO 600 IATOMA = 1, NUCIND
            DO 600 ICARTA = 1, 3
               IVARA = IBTXOR(IREPA,ISYMAX(ICARTA,1))
               IF (IBTAND(ISTBNU(IATOMA),IVARA).EQ.0) THEN
                  A = A + 1
                  B = 0
                  DO 700 IREPB  = 0, MAXREP
                  DO 700 IATOMB = 1, NUCIND
                  DO 700 ICARTB = 1, 3
                     IVARB = IBTXOR(IREPB,ISYMAX(ICARTB,1))
                     IF (IBTAND(ISTBNU(IATOMB),IVARB).EQ.0) THEN
                        B = B + 1
                        HES = D0
                        DO 800 ISYMA = 0, MAXOPR
                        IF (IBTAND(ISTBNU(IATOMA),ISYMA).EQ.0) THEN
                           FA = PT(IBTAND(IVARA,ISYMA))
                           AP = INDEX(IATOMA,ICARTA,ISYMA)
                           DO 810 ISYMB = 0, MAXOPR
                           IF (IBTAND(ISTBNU(IATOMB),ISYMB).EQ.0) THEN
                              FB = PT(IBTAND(IVARB,ISYMB))
                              BP = INDEX(IATOMB,ICARTB,ISYMB)
                              HES = HES + FA*FB*SKLTON(AP,BP)
                           END IF
  810                      CONTINUE
                        END IF
  800                   CONTINUE
                        HESSEE(A,B) = HES
                     END IF
  700             CONTINUE
               END IF
  600       CONTINUE
         END IF
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Symmetrized molecular gradient",1)
         CALL PRIGRD(GRADEE,CSTRA,SCTRA)
      END IF
      IF (MAXDER .EQ. 2 .AND. IPRINT .GT. 5) THEN
         CALL HEADER("Symmetrized molecular Hessian",1)
         CALL PRIHES(HESSEE,'CENTERS',CSTRA,SCTRA)
      END IF
      RETURN
      END
C  /* Deck sklond */
      SUBROUTINE SKLOND(FMAT,SKLTON,FACTOR,INDEX,IPOINT,NDMAT,
     &                  IFCTYP,ISYMDM,SPNORB,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      LOGICAL SPNORB
      INTEGER A, B, AP, BP, X
      DIMENSION FMAT(NBASIS,NBASIS,3), SKLTON(NBASIS,NBASIS,3),
     &          INDEX(KMAX,0:MAXREP,*), FACTOR(NBASIS,0:MAXREP),
     &          IPOINT(NBASIS,0:MAXREP), IFCTYP(NDMAT), ISYMDM(NDMAT),
     &          ISYMFC(3)
#include <pincom.h>
#include <shells.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL HEADER('Subroutine SKLOND',-1)
C
      IF (SPNORB) THEN
         NFMAT=NDMAT
         DO I=1,NFMAT
            X=ABS(IFCTYP(I))
            ISYMFC(I)=IBTXOR(ISYMAX(X,2),ISYMDM(I))
         END DO
      ELSE
         NFMAT=3
         ISYMFC(1)=ISYMAX(1,2)
         ISYMFC(2)=ISYMAX(2,2)
         ISYMFC(3)=ISYMAX(3,2)
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Input skeleton Fock matrices in SKLOND",1)
         DO 100 X = 1, NFMAT
            WRITE(LUPRI,*) 'Input skeleton Fock matrix no.',X,NFMAT
            CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                                NBASIS,NBASIS,-1,LUPRI)
  100    CONTINUE
      END IF
C
C     Form symmetric skeleton Fock matrix
C     ===================================
C
      IF (.NOT.SPNORB) THEN
      DO 200 X = 1, NFMAT
         DO 210 I = 1, NBASIS
         DO 210 J = 1, I
            FMATIJ = DP5*(FMAT(I,J,X) - FMAT(J,I,X))
            FMAT(I,J,X) =   FMATIJ
            FMAT(J,I,X) = - FMATIJ
  210    CONTINUE
  200 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER("Symmetric skeleton Fock matrix in SKLOND",1)
         DO 300 X = 1, NFMAT
            WRITE(LUPRI,*) 'Symmetric skeleton Fock matrix no.',X,NFMAT
            CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,NBASIS,NBASIS,
     &                              1,LUPRI)
  300    CONTINUE
      END IF
      END IF
C
C     Form AO Fock matrix from symmetric skeleton matrix
C     ==================================================
C
      IF (MAXREP .GT. 0) THEN
C
C        First construct index and symmetry factor arrays
C
         IORB = 0
         DO 400 ISHELL = 1, KMAX
            DO 410 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 420 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  INDEX(ISHELL,ISYMOP,ICOMP) = IORB
  420          CONTINUE
            ELSE
               ISYMF = IBTXOR(ISYMOP,IBTAND(ISYMOP,ISTBAO(ISHELL)))
               DO 430 ICOMP  = 1, KHKT(ISHELL)
                  INDEX(ISHELL,ISYMOP,ICOMP) = INDEX(ISHELL,ISYMF,ICOMP)
  430          CONTINUE
            END IF
  410       CONTINUE
  400    CONTINUE
C
         IORB = 0
         DO 500 ISHELL = 1, KMAX
            NHKTA = NHKT(ISHELL)
            DO 510 ISYMOP = 0, MAXREP
            IF (IBTAND(ISYMOP,ISTBAO(ISHELL)) .EQ. 0) THEN
               DO 520 ICOMP  = 1, KHKT(ISHELL)
                  IORB = IORB + 1
                  DO 530 ISYM = 0, MAXREP
                     ISYMPR = IBTXOR(ISYMOP,ISYM)
                     IPOINT(IORB,ISYM) = INDEX(ISHELL,ISYMPR,ICOMP)
                     FACTOR(IORB,ISYM) =
     &                        PT(IBTAND(ISYM,ISYMAO(NHKTA,ICOMP)))
  530             CONTINUE
  520          CONTINUE
            END IF
  510       CONTINUE
  500    CONTINUE
C
C        Construct full matrix from skeleton
C
         CALL DCOPY(NFMAT*(NBASIS**2),FMAT,1,SKLTON,1)
         DO 600 ISYMOP = 1, MAXREP
            DO 610 X = 1, NFMAT
               FP = PT(IBTAND(ISYMOP,ISYMFC(X)))
               DO 620 A = 1, NBASIS
                  AP = IPOINT(A,ISYMOP)
                  FA = FACTOR(A,ISYMOP)
                  DO 630 B = 1, NBASIS
                     BP = IPOINT(B,ISYMOP)
                     FB = FACTOR(B,ISYMOP)
                     FMAT(AP,BP,X)=FMAT(AP,BP,X)+FA*FB*FP*SKLTON(A,B,X)
  630             CONTINUE
  620          CONTINUE
  610       CONTINUE
  600    CONTINUE
C
         FAC = D1/(MAXREP + D1)
         CALL DSCAL(NFMAT*(NBASIS**2),FAC,FMAT,1)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("AO Fock matrices in SKLOND",1)
            DO 700 X = 1, NFMAT
               WRITE(LUPRI,*) 'AO skeleton Fock matrix no.',X,NFMAT
               CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                                   NBASIS,NBASIS,1,LUPRI)
  700       CONTINUE
         END IF
C
C        Transformation to symmetry basis
C        ================================
C
         CALL DCOPY(NFMAT*(NBASIS**2),FMAT,1,SKLTON,1)
         CALL DZERO(FMAT,NFMAT*NBASIS*NBASIS)
         ISTRA = 1
         DO 800 IREPA = 0, MAXREP
            NORBA = NAOS(IREPA+1)
            DO 810 I = ISTRA,ISTRA + NORBA - 1
               IA    = IBTAND(IBTSHR(IPIND(I),16),65535)
               NA    = IBTAND(IBTSHR(IPIND(I), 8),  255)
               NHKTA = NHKT(IA)
               KHKTA = KHKT(IA)
               MULA  = ISTBAO(IA)
               INDA  = KSTRT(IA) + NA - KHKTA
               IVARA = IBTXOR(IREPA,ISYMAO(NHKTA,NA))
               DO 820 ISYMA = 0, MAXOPR
               IF (IBTAND(ISYMA,MULA) .EQ. 0) THEN
                  FA    = PT(IBTAND(ISYMA,IVARA))
                  INDA  = INDA + KHKTA
                  ISTRB = 1
                  DO 830 IREPB = 0, MAXREP
                     NORBB = NAOS(IREPB+1)
                     DO 840 J = ISTRB,ISTRB + NORBB - 1
                        IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
                        NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
                        NHKTB = NHKT(IB)
                        KHKTB = KHKT(IB)
                        MULB  = ISTBAO(IB)
                        INDB  = KSTRT(IB) + NB - KHKTB
                        IVARB = IBTXOR(IREPB,ISYMAO(NHKTB,NB))
                        DO 850 ISYMB = 0, MAXOPR
                        IF (IBTAND(ISYMB,MULB) .EQ. 0) THEN
                           INDB = INDB + KHKTB
                           FAB  = FA*PT(IBTAND(ISYMB,IVARB))
                           DO 860 X = 1, NFMAT
                           IF (IBTXOR(IREPA,IREPB).EQ.ISYMFC(X)) THEN
                              FMAT(I,J,X) = FMAT(I,J,X)
     &                                    + FAB*SKLTON(INDA,INDB,X)
                           END IF
  860                      CONTINUE
                        END IF
  850                   CONTINUE
  840                CONTINUE
                     ISTRB = ISTRB + NORBB
  830             CONTINUE
               END IF
  820          CONTINUE
  810       CONTINUE
            ISTRA = ISTRA + NORBA
  800    CONTINUE
         IF (IPRINT .GT. 5) THEN
            CALL HEADER("SO Fock matrix in SKLOND",1)
            DO 900 X = 1, NFMAT
               CALL OUTPUT(FMAT(1,1,X),1,NBASIS,1,NBASIS,
     &                                   NBASIS,NBASIS,1,LUPRI)
  900       CONTINUE
         END IF
      END IF
      RETURN
      END
C --- end of her2fck.F ---
