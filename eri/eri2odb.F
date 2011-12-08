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
C  /* Deck odbchs */
      SUBROUTINE ODBCHS(KODCL1,KODCL2,
     &                  KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,CCFBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <maxorb.h>
      LOGICAL SAMEOD, ODALL1, ODALL2
      DIMENSION CCFBT(*), WORK(LWORK)
#include <cbieri.h>
#include <aobtch.h>
#include <symmet.h>
#include <iratdef.h>
#include <odbtch.h>
#include <odclss.h>
#include <cbirea.h>
#include <r12int.h>
C
      CALL QENTER('ODBCHS')
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
C
C     Initialize
C
C     NITCL increased from 22 to 24 (WK/UniKA/04-11-2002).
      NITCL = 22
      NITBC = 4
      NRTBC = 1
      NITPP = 0
      NRTPA = 2
      NRTPB = MAXOPR + 1
      NRTPP = NRTPA + NRTPB
C
      NITTR = 4
C
      ODTRI1 = PMSAB
      ODTRI2 = PMSCD
      ODTR12 = PMS12
      IF (ODTRI1 .NEQV. ODTRI2) ODTR12 = .FALSE.
      IF (U12INT) ODTR12 = .FALSE.
C     No symmetry for [T1,r12] integrals (WK/UniKA/04-11-2002).
      ODALL1 = .TRUE.
      ODALL2 = .TRUE.
C
      DOSORT(1) = .TRUE.
      DOSORT(2) = .TRUE.
      DOSORT(3) = .NOT.OFFCNT
      DOSORT(4) =  DIASRT
      DOSORT(5) = MAXREP .GT. 0
      DOSORT(6) = MAXREP .GT. 0
      DOSORT(7) = .TRUE.
      DOSORT(8) = .TRUE.
C     Sort for basis-set identifiers (WK/UniKA/04-11-2002).
      DOSORT(9) = LMULBS
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(2X,A,  I5)') ' NPSORT in ODBCHS ',NPSORT
         WRITE (LUPRI,'(2X,A,10L5)') ' DOSORT in ODBCHS ',
     &                               (DOSORT(I),I=1,NPSORT)
      END IF
C
C     SAMEOD: ODs are the same for both electrons
C
      SAMEOD = ODTRI1 .EQV. ODTRI2 .AND. ODALL1 .EQV. ODALL2
C
C     Gab integrals for Cauchy-Schwarz inequality
C     ============================================
C
      CALL MEMGET('REAL',KGAB,NGAB*NGAB,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KGAB),NGAB*NGAB)
      KGABWRK = KFREE
C
C     CALL GETGAB(WORK(KGAB),LWRK,IPRINT)
C
C     **********************
C     ***** Electron 1 *****
C     **********************
C
C     OD batches
C     ==========
C
C     Maximum number of OD batches
C
      NAOBC2 = NAOBCH*NAOBCH
      IF (ODTRI1) NAOBC2 = NAOBCH*(NAOBCH + 1)/2
      NPRIM  = ISUM(NAOBCH,NPRFBT,1)
      NPRIM2 = NPRIM*NPRIM
C
      CALL MEMGET('INTE',KODBCH,NITBC*NAOBC2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRDBCH,NRTBC*NAOBC2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KODPPR,NITPP*NPRIM2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRDPPR,NRTPP*NPRIM2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KODTRA,NITTR*NAOBC2,WORK,KFREE,LFREE)
      CALL GODBCH(WORK(KODBCH),WORK(KODTRA),WORK(KRDBCH),NAOBC2,
     &            WORK(KODPPR),WORK(KRDPPR),NPRIM2,CCFBT,WORK(KFREE),
     &            LFREE,ODTRI1,IPRINT)
C
C     OD classes
C     ==========
C
      CALL MEMGET('INTE',KCLASS,NITCL*NODCLS,WORK,KFREE,LFREE)
      CALL CLSINF(WORK(KODBCH),WORK(KODTRA),WORK(KCLASS),
     &            WORK(KFREE),LFREE,IPRINT)
C
      CALL MEMCHK('ODBCHS after CLSINF',WORK,1)
C
C     Compress and save results
C     =========================
C
      KLAST0 = KFREE
      NODCL1 = NODCLS
      NODBC1 = NODBCH
      NODPP1 = NODPPR
C
      KODCL1 = KGABWRK
C              ------- We save GAB for next electron
      KODBC1 = KODCL1 + (NITCL*NODCL1 - 1)/IRAT + 1
      KRDBC1 = KODBC1 + (NITBC*NODBC1 - 1)/IRAT + 1
      KODPP1 = KRDBC1 +  NRTBC*NODBC1
      KRDPP1 = KODPP1 + (NITPP*NODPP1 - 1)/IRAT + 1
      KLAST1 = KRDPP1 +  NRTPP*NODPP1
C
      IF (KLAST1 .GT. KLAST0) THEN
         WRITE (LUPRI,'(/,1X,A,2(2X,A,I10))')
     &      ' Error in ODBCHS: KLAST1 greater than KLAST0.',
     &      ' KLAST1: ',KLAST1,' KLAST0: ',KLAST0
         CALL QUIT('Error in ODBCHS')
      END IF
C
      KODCLX = KODCL1 + (KLAST0 - KGABWRK)
      KODBCX = KODBC1 + (KLAST0 - KGABWRK)
      KRDBCX = KRDBC1 + (KLAST0 - KGABWRK)
      KODPPX = KODPP1 + (KLAST0 - KGABWRK)
      KRDPPX = KRDPP1 + (KLAST0 - KGABWRK)
      KRDPPY = KRDPPX +  NRTPP*NODPP1
      KLAST  = KRDPPY +  NRTPB*NODPP1
C
      IF (KLAST.GT.LWORK) CALL STOPIT('ODBCHS','COPY 1',KLAST,LWORK)
C
      CALL ITRMAT(NITCL,NODCL1,WORK(KCLASS),WORK(KODCLX))
      CALL ITRMAT(NITBC,NODBC1,WORK(KODBCH),WORK(KODBCX))
      CALL RTRMAT(NRTBC,NODBC1,WORK(KRDBCH),WORK(KRDBCX))
      CALL ITRMAT(NITPP,NODPP1,WORK(KODPPR),WORK(KODPPX))
      CALL RTRMAT(NRTPP,NODPP1,WORK(KRDPPR),WORK(KRDPPX))
C
      KSTART = KRDPPX + NRTPA*NODPP1
      CALL DCOPY(NRTPB*NODPP1,WORK(KSTART),1,WORK(KRDPPY),1)
      CALL RTRMAT(NODPP1,NRTPB,WORK(KRDPPY),WORK(KSTART))
C
      CALL DCOPY((KLAST1-KODCL1),WORK(KODCLX),1,WORK(KODCL1),1)
C
C     **********************
C     ***** Electron 2 *****
C     **********************
C
      IF (SAMEOD) THEN
         NODCL2 = NODCL1
         NODBC2 = NODBC1
         NODPP2 = NODPP1
         KODCL2 = KODCL1
         KODBC2 = KODBC1
         KRDBC2 = KRDBC1
         KODPP2 = KODPP1
         KRDPP2 = KRDPP1
         KLAST2 = KLAST1
      ELSE
C
C        Maximum number of OD batches
C
         NAOBC2 = NAOBCH*NAOBCH
         IF (ODTRI2) NAOBC2 = NAOBCH*(NAOBCH + 1)/2
         NPRIM  = ISUM(NAOBCH,NPRFBT,1)
         NPRIM2 = NPRIM*NPRIM
C
         KODBCH = KLAST1
         KRDBCH = KODBCH + (NITBC*NAOBC2 - 1)/IRAT + 1
         KODPPR = KRDBCH +  NRTBC*NAOBC2
         KRDPPR = KODPPR + (NITPP*NPRIM2 - 1)/IRAT + 1
         KODTRA = KRDPPR +  NRTPP*NPRIM2
         KLAST  = KODTRA + (NITTR*NAOBC2 - 1)/IRAT + 1
         LWRK   = LWORK  - KLAST + 1
         IF (KLAST.GT.LWORK)CALL STOPIT('ODBCHS','GODBCH 2',KLAST,LWORK)
         CALL GODBCH(WORK(KODBCH),WORK(KODTRA),WORK(KRDBCH),NAOBC2,
     &               WORK(KODPPR),WORK(KRDPPR),NPRIM2,CCFBT,WORK(KLAST),
     &               LWRK,ODTRI2,IPRINT)
C
C        OD classes
C        ==========
C
         KCLASS = KLAST
         KLAST0 = KCLASS + (NITCL*NODCLS - 1)/IRAT + 1
         LWRK   = LWORK  - KLAST0 + 1
         IF (KLAST.GT.LWORK)CALL STOPIT('ODBCHS','CLSINF 2',KLAST,LWORK)
         CALL CLSINF(WORK(KODBCH),WORK(KODTRA),WORK(KCLASS),
     &               WORK(KLAST0),LWRK,IPRINT)
C
C        Compress and save results
C        =========================
C
         NODCL2 = NODCLS
         NODBC2 = NODBCH
         NODPP2 = NODPPR
C
         KODCL2 = KLAST1
         KODBC2 = KODCL2 + (NITCL*NODCL2 - 1)/IRAT + 1
         KRDBC2 = KODBC2 + (NITBC*NODBC2 - 1)/IRAT + 1
         KODPP2 = KRDBC2 +  NRTBC*NODBC2
         KRDPP2 = KODPP2 + (NITPP*NODPP2 - 1)/IRAT + 1
         KLAST2 = KRDPP2 +  NRTPP*NODPP2
C
         IF (KLAST2 .GT. KLAST0) THEN
            WRITE (LUPRI,'(/,1X,A,2(2X,A,I10))')
     &         ' Error in ODBCHS: KLAST2 greater than KLAST0.',
     &         ' KLAST2: ',KLAST2,' KLAST0: ',KLAST0
            CALL QUIT('Error in ODBCHS')
         END IF
C
         KODCLX = KLAST0
         KODBCX = KODCLX + (NITCL*NODCL2 - 1)/IRAT + 1
         KRDBCX = KODBCX + (NITBC*NODBC2 - 1)/IRAT + 1
         KODPPX = KRDBCX +  NRTBC*NODBC2
         KRDPPX = KODPPX + (NITPP*NODPP2 - 1)/IRAT + 1
         KRDPPY = KRDPPX +  NRTPP*NODPP2
         KLAST  = KRDPPY +  NRTPB*NODPP2
C
         IF (KLAST.GT.LWORK) CALL STOPIT('ODBCHS','COPY 2',KLAST,LWORK)
C
         CALL ITRMAT(NITCL,NODCL2,WORK(KCLASS),WORK(KODCLX))
         CALL ITRMAT(NITBC,NODBC2,WORK(KODBCH),WORK(KODBCX))
         CALL RTRMAT(NRTBC,NODBC2,WORK(KRDBCH),WORK(KRDBCX))
         CALL ITRMAT(NITPP,NODPP2,WORK(KODPPR),WORK(KODPPX))
         CALL RTRMAT(NRTPP,NODPP2,WORK(KRDPPR),WORK(KRDPPX))
C
         KSTART = KRDPPX + NRTPA*NODPP2
         CALL DCOPY(NRTPB*NODPP2,WORK(KSTART),1,WORK(KRDPPY),1)
         CALL RTRMAT(NODPP2,NRTPB,WORK(KRDPPY),WORK(KSTART))
C
         CALL ICOPY(NITCL*NODCL2,WORK(KODCLX),1,WORK(KODCL2),1)
         CALL ICOPY(NITBC*NODBC2,WORK(KODBCX),1,WORK(KODBC2),1)
         CALL DCOPY(NRTBC*NODBC2,WORK(KRDBCX),1,WORK(KRDBC2),1)
         CALL ICOPY(NITPP*NODPP2,WORK(KODPPX),1,WORK(KODPP2),1)
         CALL DCOPY(NRTPP*NODPP2,WORK(KRDPPX),1,WORK(KRDPP2),1)
      END IF
C
      KFREE = KLAST2
      LFREE = LWORK - KFREE + 1
      IF (KFREE.GT.LWORK) CALL STOPIT('ODBCHS','FREE',KFREE,LWORK)
C
      CALL QEXIT('ODBCHS')
      RETURN
      END
C  /* Deck rtrmat */
      SUBROUTINE RTRMAT(NA,NB,AMAT,BMAT)
#include <implicit.h>
      DIMENSION AMAT(NA,NB), BMAT(NB,NA)
      DO 100 I = 1, NA
      DO 100 J = 1, NB
         BMAT(J,I) = AMAT(I,J)
  100 CONTINUE
      RETURN
      END
C  /* Deck itrmat */
      SUBROUTINE ITRMAT(NA,NB,IAMAT,IBMAT)
#include <implicit.h>
      DIMENSION IAMAT(NA,NB), IBMAT(NB,NA)
      DO 100 I = 1, NA
      DO 100 J = 1, NB
         IBMAT(J,I) = IAMAT(I,J)
  100 CONTINUE
      RETURN
      END
C  /* Deck godbch */
      SUBROUTINE GODBCH(IODBCH,IODTRA,RODBCH,NAOBC2,IODPPR,RODPPR,
     &                  NPRIM2,CCFBT,WORK,LWORK,ODTRI,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      LOGICAL ODTRI
      DIMENSION IODBCH(NITBC,NAOBC2), IODTRA(NITTR,NAOBC2),
     &          RODBCH(NRTBC,NAOBC2),
     &          IODPPR(NITPP,NPRIM2), RODPPR(NRTPP,NPRIM2),
     &          CCFBT(*), WORK(LWORK)
#include <aobtch.h>
#include <odbtch.h>
#include <odclss.h>
C
      IBATCH = 1
      ICLASS = 0
      IPRADR = 1
      DO 100 NHKTA = 1, MAXQN
         MAXB = MAXQN
         IF (ODTRI) MAXB = NHKTA
         DO 200 NHKTB = 1, MAXB
            CALL ODSORT(NHKTA,NHKTB,ODTRI,NBCHAB,IODBCH(1,IBATCH),
     &                  IODTRA(1,IBATCH),RODBCH(1,IBATCH),ICLASS,IODPPR,
     &                  RODPPR,IPRADR,CCFBT,WORK,LWORK,IPRINT)
            IBATCH = IBATCH + NBCHAB
  200    CONTINUE
  100 CONTINUE
      NODBCH = IBATCH - 1
      NODCLS = ICLASS
      NODPPR = IPRADR - 1
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from GODBCH',-1)
         WRITE(LUPRI,'(1X,A, I5)')' Number of OD batches NODBCH:',NODBCH
         WRITE(LUPRI,'(1X,A,3I5)')' Number of AO batches NAOBCH:',NAOBCH
         WRITE(LUPRI,'(1X,A,3I5)')' Number of OD classes NODCLS:',NODCLS
         WRITE(LUPRI,'(1X,A, L5)')' ODTRI:                      ',ODTRI
         CALL HEADER(
     &    '       #     class     A    B           Cauchy-Schwarz',1)
         DO 300 I = 1, NODBCH
            WRITE (LUPRI,'(6X,I5,5X,4I5,5X,F20.15)')
     &        I, (IODBCH(J,I),J=1,NITBC), RODBCH(1,I)
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck odsort */
      SUBROUTINE ODSORT(NHKTA,NHKTB,ODTRI,NBCHAB,IODBCH,IODTRA,RODBCH,
     &                  ICLASS,IODPPR,RODPPR,IPRADR,CCFBT,WORK,LWORK,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
      LOGICAL ODTRI, TRIBCH
      DIMENSION IODBCH(NITBC,*), IODTRA(NITTR,*), RODBCH(NRTBC,*),
     &          WORK(LWORK), IODPPR(NITPP,*), RODPPR(NRTPP,*),
     &          CCFBT(*)
#include <odclss.h>
#include <aobtch.h>
#include <odbtch.h>
C
      NQNBTA = NQNBT(NHKTA)
      NQNBTB = NQNBT(NHKTB)
      TRIBCH = ODTRI .AND. NHKTA .EQ. NHKTB
      NBCH   = NQNBTA*NQNBTB
      IF (TRIBCH) NBCH = NQNBTA*(NQNBTA + 1)/2
C
      KISORT = 1
      KDSORT = KISORT + (3*NBCH - 1)/IRAT + 1
      KSORT  = KDSORT +  NBCH
      KDSRT  = KSORT  + (NBCH*NPSORT - 1)/IRAT + 1
      KCLASS = KDSRT  + (NBCH - 1)/IRAT + 1
      KLAST  = KCLASS + (NBCH - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ODSORT',' ',KLAST,LWORK)
      LWRK   = KLAST - LWORK + 1
      CALL ODSOR1(NHKTA,NHKTB,TRIBCH,NBCH,NBCHAB,IODBCH,IODTRA,RODBCH,
     &            ICLASS,IODPPR,RODPPR,IPRADR,WORK(KISORT),WORK(KDSORT),
     &            WORK(KSORT),WORK(KDSRT),WORK(KCLASS),CCFBT,
     &            WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck odsor1 */
      SUBROUTINE ODSOR1(NHKTA,NHKTB,TRIBCH,NBCH,NBCHAB,IODBCH,IODTRA,
     &                  RODBCH,ICLASS,IODPPR,RODPPR,IPRADR,ISORT,DSORT,
     &                  NSORT,KODSRT,KCLASS,CCFBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
      LOGICAL TRIBCH, ODBTGT
      INTEGER A, B
#include <aobtch.h>
#include <odbtch.h>
#include <odclss.h>
      DIMENSION IODBCH(NITBC,*), IODTRA(NITTR,*), RODBCH(NRTBC,*),
     &          ISORT(3,NBCH), DSORT(1,NBCH), NSORT(NPSORT,NBCH),
     &          IODPPR(*), RODPPR(NRTPP,*),
     &          KODSRT(NBCH), KCLASS(NBCH), CCFBT(*), WORK(LWORK)
Cx    dimension of IODPPR was IODPPR(NITPP,*) but changed because NITPP = 0
Cx    in this version. /hjaaj Oct07 (gave out of bounds in CALL NODDIM below)
C
      CALL IZERO(NSORT,NPSORT*NBCH)
C
C     Collect ODs with non-zero elements
C
      IBATCH = 0
      DO 100 A = 1, NQNBT(NHKTA)
         MAXB = NQNBT(NHKTB)
         IF (TRIBCH) MAXB = A
         DO 200 B = 1, MAXB
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(/,2X,A,2I5)') 'A, B in ODSOR1 ',A,B
            END IF
            IA = KAOSRT(KQNBT(NHKTA) - 1 + A)
            IB = KAOSRT(KQNBT(NHKTB) - 1 + B)
            J_IODPPR = (IPRADR-1)*NITPP + 1
            CALL NODDIM(IA,IB,NPRFAB,NCTFAB,NCNTAB,NGTOAB,ISTBLA,ISTBLB,
     &                  ICMATA,ICMATB,PRXFAB,IODPPR(J_IODPPR), ! Cx IODPPR(1,IPRADR) not OK for NITPP=0
     &                  RODPPR(1,IPRADR),CCFBT)
            IF (NPRFAB .GT. 0) THEN
               IBATCH = IBATCH + 1
               ISORT(1,IBATCH) = IA
               ISORT(2,IBATCH) = IB
               ISORT(3,IBATCH) = IPRADR
               DSORT(1,IBATCH) = PRXFAB
               IF (DOSORT(1)) NSORT(1,IBATCH) = NPRFAB
               IF (DOSORT(2)) NSORT(2,IBATCH) = NCTFAB
               IF (DOSORT(3)) NSORT(3,IBATCH) = NCNTAB
               IF (DOSORT(4)) NSORT(4,IBATCH) = NGTOAB
               IF (DOSORT(5)) NSORT(5,IBATCH) = ISTBLA
               IF (DOSORT(6)) NSORT(6,IBATCH) = ISTBLB
               IF (DOSORT(7)) NSORT(7,IBATCH) = ICMATA
               IF (DOSORT(8)) NSORT(8,IBATCH) = ICMATB 
C              Sort for basis-set identifiers (WK/UniKA/04-11-2002).
               IF (DOSORT(9)) NSORT(9,IBATCH) = MBIDAB(IA,IB)
               IPRADR = IPRADR + NPRFAB
               IF (IPRINT .GT. 5) THEN
                  WRITE (LUPRI,'(2X,A, I5)') ' IBATCH ', IBATCH
                  WRITE (LUPRI,'(2X,A,2I5)') ' IA,IB  ', IA,IB
                  WRITE (LUPRI,'(2X,A, I5)') ' IPRADR ', IPRADR
                  WRITE (LUPRI,'(2X,A, F12.6)') ' PRXFAB ', PRXFAB
                  WRITE (LUPRI,'(2X,A,2I5)') ' NPRFAB,NCTFAB ',
     &                                         NPRFAB,NCTFAB
                  WRITE (LUPRI,'(2X,A,2I5)') ' NCNTAB,NGTOAB ',
     &                                         NCNTAB,NGTOAB
                  WRITE (LUPRI,'(2X,A,2I5)') ' ISTBLA,ISTBLB ',
     &                                         ISTBLA,ISTBLB
                  WRITE (LUPRI,'(2X,A,2I5)') ' ICMATA,ICMATB ',
     &                                         ICMATA,ICMATB
               END IF
            END IF
  200    CONTINUE
  100 CONTINUE
      NBCHAB = IBATCH
C
      IF (NBCHAB .GT. 0) THEN
C
C        Sort ODs
C
         DO 300 I = 1, NBCHAB
            KODSRT(I) = I
  300    CONTINUE
         DO 400 I = 2, NBCHAB
            IOD = KODSRT(I)
            DO 500 J = I - 1, 1, -1
               IF (ODBTGT(IOD,KODSRT(J),NSORT)) GO TO 550
               KODSRT(J+1) = KODSRT(J)
  500       CONTINUE
            J = 0
  550       KODSRT(J+1) = IOD
  400    CONTINUE
C
C        Classify ODs
C        ============
C
         CALL ODCLS(KCLASS,KODSRT,NSORT,NBCHAB,ICLASS,NCLASS,IPRINT)
C
C        Collect information on each OD batch
C        ====================================
C
         DO 600 I = 1, NBCHAB
C
C           IODBC1/IODBC2
C           =============
C
            IODBCH(1,I) = KCLASS(I)
            IODBCH(2,I) = ISORT(1,KODSRT(I))
            IODBCH(3,I) = ISORT(2,KODSRT(I))
            IODBCH(4,I) = ISORT(3,KODSRT(I))
C
C           RODBC1/RODBC2
C           =============
C
            RODBCH(1,I) = DSORT(1,KODSRT(I))
C
C           to be transferred to IODCL1/IODCL2 in CLSINF
C           ============================================
C
            IODTRA(1,I) = NSORT(1,KODSRT(I))
            IODTRA(2,I) = NSORT(2,KODSRT(I))
            IODTRA(3,I) = NSORT(3,KODSRT(I))
            IODTRA(4,I) = NSORT(4,KODSRT(I))
  600    CONTINUE
C
      END IF
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from ODSOR1',-1)
         WRITE (LUPRI,'(1X,A,2I5)') 'NHKTA, NHKTB: ',NHKTA,NHKTB
         WRITE (LUPRI,'(1X,A, I5)') 'NBCHAB:       ',NBCHAB
         WRITE (LUPRI,'(1X,A, I5)') 'NCLASS:       ',NCLASS
         WRITE (LUPRI,'(1X,A,10I5/,(15X,10I5))')
     &      'IODBCH(1,I):  ',(IODBCH(1,I),I=1,NBCHAB)
         WRITE (LUPRI,'(1X,A,10I5/,(15X,10I5))')
     &      'IODBCH(2,I):  ',(IODBCH(2,I),I=1,NBCHAB)
         WRITE (LUPRI,'(1X,A,10I5/,(15X,10I5))')
     &      'IODBCH(3,I):  ',(IODBCH(3,I),I=1,NBCHAB)
         WRITE (LUPRI,'(1X,A,10I5/,(15X,10I5))')
     &      'IODBCH(4,I):  ',(IODBCH(4,I),I=1,NBCHAB)
         WRITE (LUPRI,'(1X,A,9E8.2/,(8X,9E8.2))')
     &      'CSOD:  ',(RODBCH(1,I),I=1,NBCHAB)
      END IF
C
      RETURN
      END
C  /* Deck noddim */
      SUBROUTINE NODDIM(A,B,NPRFAB,NCTFAB,NCNTAB,NGTOAB,ISTBLA,ISTBLB,
     &                  ICMATA,ICMATB,PRXFAB,IODPPR,RODPPR,CCFBT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <pi.h>
      PARAMETER (D0 = 0.D0)
      LOGICAL GCNT
      INTEGER A, B, R, R1
      DIMENSION IODPPR(NITPP,*), RODPPR(NRTPP,*), PFAC(8)
      DIMENSION CCFBT(MXPRIM*MXCONT)
#include <cbieri.h>
#include <erithr.h>
#include <odclss.h>
#include <aobtch.h>
#include <symmet.h>
C
#include <ibtfun.h>
      XPAR(I) = PT(IBTAND(ISYMAX(1,1),I))
      YPAR(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZPAR(I) = PT(IBTAND(ISYMAX(3,1),I))
C
      ISTBLA = ISTBBT(A)
      ISTBLB = ISTBBT(B)
      ISTBLR = IBTOR(ISTBLA,ISTBLB)
      MLTPR  = MULT(ISTBLR)
C
      NCTFAB = NCTFBT(A)*NCTFBT(B)
      ICMATA = KCMTBT(A)
      ICMATB = KCMTBT(B)
C
      GCNT   =  GENCON_ERI .OR. NCTFAB .GT. 1
C
      NCNTAB = 2
      IF (NCNTBT(A) .EQ. NCNTBT(B) .AND. MLTPR.EQ.1) NCNTAB = 1
      NGTOAB = 2
      IF (A .EQ. B) NGTOAB = 1
C     IF (A .EQ. B .AND. MLTPR .EQ. 1) NGTOAB = 1
      IJ = 0
      PRXFAB = D0
      DO 100 I = 1, NPRFBT(A)
         DO 200 J = 1, NPRFBT(B)
C
            INDXPA = KEXPBT(A) - 1 + I
            INDXPB = KEXPBT(B) - 1 + J
            INDCFA = KCCFBT(A) - 1 + I
            INDCFB = KCCFBT(B) - 1 + J
C
            EXPA = EXPBT(INDXPA)
            EXPB = EXPBT(INDXPB)
            EXPR = EXPA*EXPB/(EXPA + EXPB)
C
            R1 = 0
            FACMAX = D0
            DO R = 0, MAXOPR
            IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
               R1 = R1 + 1
               DST2 = (CORXBT(A) - XPAR(R)*CORXBT(B))**2
     &              + (CORYBT(A) - YPAR(R)*CORYBT(B))**2
     &              + (CORZBT(A) - ZPAR(R)*CORZBT(B))**2
               FAC  = R2PI52*EXP(-EXPR*DST2)/(EXPA + EXPB)
               IF (.NOT.GCNT) FAC = CCFBT(INDCFA)*CCFBT(INDCFB)*FAC
               FACMAX = MAX(FACMAX,ABS(FAC))
               PFAC(R1) = FAC
            END IF
            END DO 
C
            IF (GCNT .OR. FACMAX.GT.THRSH) THEN
               IJ = IJ + 1
               RODPPR(1,IJ) = EXPA
               RODPPR(2,IJ) = EXPB
               DO R = 1, MLTPR
                  RODPPR(R+2,IJ) = PFAC(R)
               END DO
               PRXFAB = MAX(PRXFAB,FACMAX)
            END IF
C
  200    CONTINUE
  100 CONTINUE
      NPRFAB = IJ
      RETURN
      END
C  /* Deck odbtgt */
      LOGICAL FUNCTION ODBTGT(IOD1,IOD2,NSORT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <aobtch.h>
#include <odbtch.h>
      DIMENSION NSORT(NPSORT,*)
C
C     This routine sorts according to NPSORT criteria
C
      IPAR = 0
  100 CONTINUE
         IF (IPAR .EQ. NPSORT) GO TO 200
         IPAR = IPAR + 1
         ODBTGT = NSORT(IPAR,IOD1) .GT. NSORT(IPAR,IOD2)
         IF (NSORT(IPAR,IOD1) .EQ. NSORT(IPAR,IOD2)) GO TO 100
 200  CONTINUE
      RETURN
      END
C  /* Deck odcls */
      SUBROUTINE ODCLS(KCLASS,KODSRT,NSORT,NBCHAB,ICLASS,NCLASS,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <aobtch.h>
#include <odbtch.h>
      LOGICAL NEWCLS
      DIMENSION KCLASS(*), KODSRT(*), NSORT(NPSORT,*)
C
      ICLASS = ICLASS + 1
      NCLASS = 1
      KCLASS(1) = ICLASS
      DO 100 I = 2, NBCHAB
         M = KODSRT(I)
         N = KODSRT(I-1)
C
         NEWCLS = .FALSE.
         DO J = 1, NPSORT
            NEWCLS = NEWCLS .OR. NSORT(J,M) .NE. NSORT(J,N)
         END DO  
         IF (NEWCLS) THEN
            ICLASS = ICLASS + 1
            NCLASS = NCLASS + 1
         END IF
         KCLASS(I) = ICLASS
  100 CONTINUE
      RETURN
      END
C  /* Deck clsinf */
      SUBROUTINE CLSINF(IODBCH,IODTRA,IODCLS,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <shells.h>
#include <odbtch.h>
#include <aobtch.h>
      DIMENSION IODBCH(NITBC,*), IODTRA(NITTR,*), IODCLS(NITCL,NODCLS),
     &          WORK(LWORK)
#include <odclss.h>
C
      CALL IZERO(IODCLS,NITCL*NODCLS)
C
C     Number of batches in each class
C     ===============================
C
      DO I = 1, NODBCH
         IODCLS(1,IODBCH(1,I)) = IODCLS(1,IODBCH(1,I)) + 1
      END DO
C
C     Pointers to first element in each class
C     =======================================
C
      IJ = 1
      DO I = 1, NODCLS
         IODCLS(2,I) = IJ
         IJ = IJ + IODCLS(1,I)
      END DO  
C
C     Class info IODCL1/IODCL2
C     ========================
C
      DO 300 I = 1, NODCLS
         IPOINT = IODCLS(2,I)
         IA = IODBCH(2,IPOINT)
         IB = IODBCH(3,IPOINT)
         IODCLS( 3,I) = NHKTBT(IA)
         IODCLS( 4,I) = NHKTBT(IB)
         IODCLS( 5,I) = KHKTBT(IA)
         IODCLS( 6,I) = KHKTBT(IB)
         IODCLS( 7,I) = KCKTBT(IA)
         IODCLS( 8,I) = KCKTBT(IB)
         IODCLS( 9,I) = ISTBBT(IA)
         IODCLS(10,I) = ISTBBT(IB)
         IF (DOSORT(1)) IODCLS(11,I) = IODTRA(1,IPOINT)
         IF (DOSORT(2)) IODCLS(12,I) = IODTRA(2,IPOINT)
         IF (DOSORT(3)) IODCLS(13,I) = IODTRA(3,IPOINT)
         IF (DOSORT(4)) IODCLS(14,I) = IODTRA(4,IPOINT)
         IODCLS(15,I) = NPRFBT(IA)
         IODCLS(16,I) = NPRFBT(IB)
         IODCLS(17,I) = NCTFBT(IA)
         IODCLS(18,I) = NCTFBT(IB)
         IODCLS(19,I) = KCCFBT(IA)
         IODCLS(20,I) = KCCFBT(IB)
C        Basis-set identifiers (WK/UniKA/04-11-2002).
         IODCLS(21,I) = MBIDBT(IA)
         IODCLS(22,I) = MBIDBT(IB)
  300 CONTINUE
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from CLSINF',-1)
         WRITE (LUPRI,'(1X,A,I5)') ' Number of classes NODCLS:',NODCLS
         CALL HEADER('Array IODCLS in CLSINF',-1)
         DO I = 1, NODCLS
            WRITE (LUPRI,'(20I4)') I, (IODCLS(J,I), J=1,NITCL)
         END DO  
      END IF
C
      CALL ODSTAT(IODCLS,WORK,LWORK,IPRINT)
      RETURN
      END
C  /* Deck odstat */
      SUBROUTINE ODSTAT(IODCLS,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <odbtch.h>
      DIMENSION IODCLS(NITCL,NODCLS), WORK(LWORK)
#include <odclss.h>
C
C     IODCLS(1,*): # ods
C     IODCLS(2,*): pointer to first od
C     IODCLS(3,*): NHKTA
C     IODCLS(4,*): NHKTB
C     IODCLS(5,*): NPRFAB
C     IODCLS(6,*): NCTFAB
C     IODCLS(7,*): NCNTAB
C
      NPRFT  = 0
      NCTFT  = 0
      NONEC  = 0
      NTWOC  = 0
      NPRFT1 = 0
      NCTFT1 = 0
      NPRFT2 = 0
      NCTFT2 = 0
      DO 100 I = 1, NODCLS
         NODS   = IODCLS(1,I)
         NPRFAB = IODCLS(11,I)
         NCTFAB = IODCLS(12,I)
         NCNTAB = IODCLS(13,I)
         NPRFT  = NPRFT + NODS*NPRFAB
         NCTFT  = NCTFT + NODS*NCTFAB
         IF (I .EQ. 1) THEN
            MINPRB = NPRFAB
            MINCTB = NCTFAB
            MINPRC = NODS*NPRFAB
            MINCTC = NODS*NCTFAB
            MAXPRB = NPRFAB
            MAXCTB = NCTFAB
            MAXPRC = NODS*NPRFAB
            MAXCTC = NODS*NCTFAB
         ELSE
            MINPRB = MIN(MINPRB,NPRFAB)
            MINCTB = MIN(MINCTB,NCTFAB)
            MINPRC = MIN(MINPRC,NODS*NPRFAB)
            MINCTC = MIN(MINCTC,NODS*NCTFAB)
            MAXPRB = MAX(MAXPRB,NPRFAB)
            MAXCTB = MAX(MAXCTB,NCTFAB)
            MAXPRC = MAX(MAXPRC,NODS*NPRFAB)
            MAXCTC = MAX(MAXCTC,NODS*NCTFAB)
         END IF
         IF (NCNTAB .EQ. 1) THEN
            NONEC  = NONEC + 1
            NPRFT1 = NPRFT1 + NODS*NPRFAB
            NCTFT1 = NCTFT1 + NODS*NCTFAB
         ELSE
            NTWOC = NTWOC + 1
            NPRFT2 = NPRFT2 + NODS*NPRFAB
            NCTFT2 = NCTFT2 + NODS*NCTFAB
         END IF
  100 CONTINUE
      NPRFA = NINT(FLOAT(NPRFT)/FLOAT(NODCLS))
      NCTFA = NINT(FLOAT(NCTFT)/FLOAT(NODCLS))
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('OD statistics from ODSTAT',-1)
         WRITE (LUPRI,'(1X,3(A,I5))')
     &      ' Number of OD classes (one-center + two-center):',
     &      NONEC, ' +',NTWOC,' =',NODCLS
C
         WRITE (LUPRI,'(/,1X,A,I5)')
     &   ' Total number of primitive ODs:                   ',NPRFT
         WRITE (LUPRI,'(/,1X,A,2I5)')
     &   ' Number of primitive one- and two-center ODs:     ',NPRFT1,
     &                                                        NPRFT2
         WRITE (LUPRI,'(1X,A,I5)')
     &   ' Average number of primitive OD batches in class: ',NPRFA
         WRITE (LUPRI,'(1X,A,2I5)')
     &   ' Min/max number of primitives in OD batch:        ',MINPRB,
     &                                                        MAXPRB
         WRITE (LUPRI,'(1X,A,2I5)')
     &   ' Min/max number of primitives in OD class:        ',MINPRC,
     &                                                        MAXPRC
C
         WRITE (LUPRI,'(/1X,A,I5)')
     &   ' Total number of contracted ODs:                  ',NCTFT
         WRITE (LUPRI,'(/,1X,A,2I5)')
     &   ' Number of contracted one- and two-center ODs:    ',NCTFT1,
     &                                                        NCTFT2
         WRITE (LUPRI,'(1X,A,I5)')
     &   ' Average number of contracted OD batches in class:',NCTFA
         WRITE (LUPRI,'(1X,A,2I5)')
     &   ' Min/max number of contracted in OD batch:        ',MINCTB,
     &                                                        MAXCTB
         WRITE (LUPRI,'(1X,A,2I5)')
     &   ' Min/max number of contracted in OD class:        ',MINCTC,
     &                                                        MAXCTC
      END IF
C
      RETURN
      END
C  /* Deck MBIDAB */
      INTEGER FUNCTION MBIDAB(I,J)
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
#include <shells.h>
      MBIDAB = MBIDBT(I) + MBIDBT(J) * 2
      RETURN
      END
