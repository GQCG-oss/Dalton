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
C=====================================================================
C  /* Deck btrabidel */
      SUBROUTINE  BTRABIDEL(LUDENSTY,FNDENSTY,LUDEN4MO,FNDEN4MO,
     &                      ISYDEN,XCMO,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: backtranform the last index of the symmetrized (T) 
C              ABIC density to AO for later use in CCSD(T) FOPs 
C              and gradient
C     Loop idelta first, then C
C     S. Coriani, January 2002
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER LWORK
      integer LUDENSTY, LUDEN4MO
      integer IOABIC,IOABID,KCMOFF
      integer ISYMC,ISYABI,ISYDEL,ISYM, ISYDEN, ISYCMO
      integer IDEL, ID
      integer KABIDS,KABICS,KEND1,LWRK1,LENGTH,LENGTH1,iii

#if defined (SYS_CRAY)
      REAL WORK(LWORK),XCMO(*)
      REAL ZERO, ONE, TWO
      real ddot,xtest,xnorm,xnorm1
#else
      DOUBLE PRECISION WORK(LWORK),XCMO(*)
      DOUBLE PRECISION ZERO, ONE, TWO
      double precision ddot,xtest,xnorm,xnorm1
#endif


      CHARACTER*9 FNDENSTY, FNDEN4MO

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)' Backtransform symmetrized abic CCSD(T) density'
        WRITE(LUPRI,*)' LUDENSTY:', LUDENSTY,' FNDENSTY: ',FNDENSTY
        WRITE(LUPRI,*)' ISYDEN:', ISYDEN
        CALL FLSHFO(LUPRI)
      END IF
C
C-------------------------------
C     Symmetry
C-------------------------------
C
      ISYCMO = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      LENGTH = 0
      DO ISYM = 1, NSYM
        LENGTH = MAX(LENGTH,NCKASR(ISYM))
      END DO

      KABIDS = 1
      KABICS = KABIDS + LENGTH
      KEND1  = KABICS + LENGTH
      LWRK1  = LWORK - KEND1


      CALL DZERO(WORK(KABIDS),LENGTH)
      CALL DZERO(WORK(KABICS),LENGTH)

c---------------------------------------------------------
c Backtransform ds_abic (abi,c) = sum_c ds_abi,c CMO_del,c
c loop delta
c   loop c
c     d_abi^d = d_abi^d + d_abi^c * CMO^d,c (daxpy)
c   end loop c
c   dump on file for a given delta
c end loop delta
c---------------------------------------------------------

      xnorm1 = zero
      DO ISYDEL = 1, NSYM

         DO IDEL = 1, NBAS(ISYDEL)
         
            !ID = IDEL - IBAS(ISYDEL)

            ID = IDEL 

            ISYMC  = MULD2H(ISYCMO,ISYDEL)
            ISYABI = MULD2H(ISYDEN,ISYMC)
            
            xnorm = zero
            DO C = 1, NVIR(ISYMC)

               KCMOFF = IGLMVI(ISYDEL,ISYMC) 
     &                  + NBAS(ISYDEL)*(C-1) + ID

            !read in ABI;C block of 4MO density

               IOABIC = ICDKVI(ISYABI,ISYMC) + 
     &                  NCKASR(ISYABI)*(C-1) + 1

               CALL GETWA2(LUDEN4MO,FNDEN4MO,WORK(KABICS),
     &                     IOABIC,NCKASR(ISYABI))

               xnorm = xnorm +
     &           ddot(NCKASR(ISYABI),WORK(KABICS),1,WORK(KABICS),1)

!               write(lupri,*) '--------------------------'
!               write(lupri,*) 'Btrabidel: read in ds_abic'
!               write(lupri,*) 'isydel: ', isydel, ' id: ', id, ' c:', c
!               do iii = 1, NCKASR(ISYABI)
!                  write(lupri,*) 'ds_abi;c(',iii,')', work(kabics+iii-1)
!               end do
!               write(lupri,*) '--------------------------'

               CALL DAXPY(NCKASR(ISYABI),XCMO(KCMOFF),
     &                    WORK(KABICS),1,WORK(KABIDS),1)

            END DO
         if (locdbg) then
            write(lupri,*) '--------------------------'
            write(lupri,*)'isydel: ',isydel,' id: ',id,' isymc: ',isymc
            write(lupri,*) 'btrabic: norm of d_abic : ', xnorm

           write(lupri,*) '--------------------------'
           write(lupri,*) 'isydel: ', isydel, ' id: ', id, ' c:', c
          xtest = ddot(NCKASR(ISYABI),WORK(KABIDS),1,WORK(KABIDS),1)
           write(lupri,*)'Norm of abi;del of (T) dens', xtest
           xnorm1 = xnorm1 + xtest
         end if

            IOABID = ICDKAO(ISYABI,ISYDEL) + 
     &               NCKASR(ISYABI)*(ID-1) + 1

            CALL PUTWA2(LUDENSTY,FNDENSTY,WORK(KABIDS),
     &                  IOABID,NCKASR(ISYABI))

            CALL DZERO(WORK(KABIDS),NCKASR(ISYABI))

         END DO !IDEL

      END DO !ISYDEL
      IF (LOCDBG) THEN
         WRITE(LUPRI,*)'Total norm of abi;del of (T) dens', XNORM1
      END IF
C-----------------------------------------------------------------------------

      RETURN
      END

C=====================================================================
C  /* Deck btraibdel */
      SUBROUTINE  BTRAIBDEL(LUDENSTY,FNDENSTY,LUDEN4MO,FNDEN4MO,
     &                      ISYDEN,XCMO,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: backtranform the last index of the (T) AIBC density to AO 
C              for later use in CCSD(T) FOPs and gradient
C              Loop idelta first, then C
C     S. Coriani, February 2002
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>


      INTEGER LWORK
      INTEGER LUDENSTY, LUDEN4MO
      INTEGER IOFAIBC,IOFAIBD,KCMOFF
      INTEGER ISYCMO,ISYMC,ISYAIB,ISYDEL,ISYM, ISYDEN
      INTEGER IDEL, ID
      INTEGER KAIB_DE,KAIB_C,KEND1,LWRK1,LENGTH,III

      CHARACTER*9 FNDENSTY, FNDEN4MO

#if defined (SYS_CRAY)
      REAL WORK(LWORK),XCMO(*)
      REAL ZERO, ONE, TWO
      REAL DDOT,XTEST,XNORM,XNORM1
#else
      DOUBLE PRECISION WORK(LWORK),XCMO(*)
      DOUBLE PRECISION ZERO, ONE, TWO
      double precision ddot,xtest,xnorm,xnorm1
#endif


      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*)'                                              '
        WRITE(LUPRI,*)' LUDENSTY: ', LUDENSTY,' FNDENSTY: ', FNDENSTY
        WRITE(LUPRI,*)' ISYDEN:', ISYDEN
      END IF
C
C-------------------------------
C     Symmetry 
C-------------------------------
C
      ISYCMO = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      LENGTH = 0
      DO ISYM = 1, NSYM
        LENGTH = MAX(LENGTH,NCKATR(ISYM))
      END DO

      KAIB_DE = 1
      KAIB_C  = KAIB_DE + LENGTH
      KEND1   = KAIB_C  + LENGTH
      LWRK1   = LWORK - KEND1


      CALL DZERO(WORK(KAIB_DE),LENGTH)
      CALL DZERO(WORK(KAIB_C),LENGTH)

c---------------------------------------------------------
c Backtransform ds_aib;del (aib;del) = sum_c ds_aib;c CMO_del;c
c loop delta
c   loop c
c     d_aib^d = d_aib^d + d_aib^c * CMO^d,c (daxpy)
c   end loop c
c   dump on file for a given delta
c end loop delta
c---------------------------------------------------------

      xnorm1 = zero
      DO ISYDEL = 1, NSYM

         DO IDEL = 1, NBAS(ISYDEL)
         
            !ID = IDEL - IBAS(ISYDEL)

            ID = IDEL 
         
            ISYMC = MULD2H(ISYCMO,ISYDEL)
            ISYAIB = MULD2H(ISYDEN,ISYMC)
            
            xnorm = zero
            DO C = 1, NVIR(ISYMC)

               KCMOFF = IGLMVI(ISYDEL,ISYMC) 
     &                  + NBAS(ISYDEL)*(C-1) + ID

               !read in AIB;C block of 4MO density

               IOFAIBC = ICKBD(ISYAIB,ISYMC) 
     &                   + NCKATR(ISYAIB)*(C-1) + 1

               CALL GETWA2(LUDEN4MO,FNDEN4MO,WORK(KAIB_C),IOFAIBC,
     &                     NCKATR(ISYAIB))

               xnorm = xnorm + 
     &           ddot(NCKATR(ISYAIB),WORK(KAIB_C),1,WORK(KAIB_C),1)


!               write(lupri,*) '--------------------------'
!               write(lupri,*) 'Btraibdel: read in d_aibc or d_iabc'
!               write(lupri,*) 'isydel: ', isydel, ' id: ', id, ' c:', c
!               do iii = 1, NCKATR(ISYAIB)
!                  write(lupri,*) 'ds_abi;c(',iii,')', work(kaib_c+iii-1)
!               end do
!               write(lupri,*) '--------------------------'
  
!               write(lupri,*) 
!     &         'btraibdel: XCMO(',kcmoff,') = ', XCMO(kcmoff)

               CALL DAXPY(NCKATR(ISYAIB),XCMO(KCMOFF),WORK(KAIB_C),1,
     &                                  WORK(KAIB_DE),1)
            END DO
         if (locdbg) then
            write(lupri,*) '--------------------------'
            write(lupri,*)'isydel: ',isydel,' id: ',id,' isymc: ',isymc
            write(lupri,*) 'btraibc: norm of d_aibc : ', xnorm

        write(lupri,*) '--------------------------'
        write(lupri,*) 'isydel: ', isydel, ' id: ', id
        xtest = ddot(NCKATR(ISYAIB),WORK(KAIB_DE),1,WORK(KAIB_DE),1)
        write(lupri,*)'Norm of aib;del of (T) dens', xtest
        xnorm1 = xnorm1 + xtest
         end if

            IOFAIBD = ICKDAO(ISYAIB,ISYDEL) 
     &                + NCKATR(ISYAIB)*(ID-1)+1

            CALL PUTWA2(LUDENSTY,FNDENSTY,WORK(KAIB_DE),IOFAIBD,
     &                           NCKATR(ISYAIB))

            CALL DZERO(WORK(KAIB_DE),NCKATR(ISYAIB))

         END DO !IDEL

      END DO !ISYDEL
         if (locdbg) then
      write(lupri,*)'Total norm of aib;del of (T) dens', xnorm1
         end if
C-----------------------------------------------------------------------------

      RETURN
      END


C=====================================================================
C  /* Deck btraijdel */
      SUBROUTINE  BTRAIJDEL(XCMO,XDENS,ISYDEN,
     &                      LUDENSTY,FNDENSTY,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: backtransform the last index of the (T) densities d_aijk(aijk)
C              to AO basis. Dump result on file.
C     S. Coriani, February 2002
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER LWORK
      INTEGER ISYRES, KSTART, ISYDEN
      INTEGER LWRK1,KEND1
      INTEGER IOFF, ISYMA, ISYMI, ISYMJ, ISYMK, ISYDEL
      INTEGER KAIJK
      integer iii,nbasdel,kend2,lwrk2
      integer ISYAIJ,KCMOFFK
      integer ISYM
      integer ISYMAI
      integer isycmo, KDENSBT, koff, ISYMD
      integer kcmo,kcmoffj,koffden,koffres,ntotijk
      integer ntotaij
      integer LUDENSTY, LENAIJD

#if defined (SYS_CRAY)
      REAL WORK(LWORK),XCMO(*),XDENS(*)
      REAL ZERO, ONE, TWO
      real ddot,xnorm2,xnorm
#else
      DOUBLE PRECISION WORK(LWORK),XCMO(*),XDENS(*)
      DOUBLE PRECISION ZERO, ONE, TWO
      double precision ddot,xnorm2,xnorm
#endif


      CHARACTER FNDENSTY*9

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C
C#include <leinf.h>
C
      if (locdbg) then
        write(lupri,*) ' Backtransform (T) density aijk '
        write(lupri,*) ' LUDENSTY:',LUDENSTY,' FNDENSTY:',FNDENSTY
      end if

c---------------------------------------------------------
c Backtransform ds_aijk (aij,k) = sum_k d_aij,k CMO_del,k
c DEVO FARE UN LOOP SU ISYMD???????????????????????????
c Dipende da come li voglio prelevare dal file.....
c Backtransform: 
C d_aijdel (aij,del) = sum_k d_aijk(aij,k) CMO_del,k
c---------------------------------------------------------

      KEND1 = 1
      KDENSBT = KEND1
      KEND2   = KDENSBT + NTOTOC(ISYDEN)
      LWRK2   = LWORK - KEND2

      CALL DZERO(WORK(KDENSBT),NTOTOC(ISYDEN))

      ISYCMO = 1

      do isyaij = 1, nsym
         isymk  = muld2h(isyden,isyaij)
         isydel = muld2h(isycmo,isymk) 

         kcmoffk = IGLMRH(ISYDEL,ISYMK)  + 1
         koffden = ISAIKJ(ISYAIJ,ISYMK)  + 1
         koffres = ICKID(ISYAIJ,ISYDEL) + KDENSBT
 
         nbasdel = MAX(NBAS(ISYDEL),1)
         ntotaij = MAX(NCKI(ISYAIJ),1)

         call dgemm('N','T',ncki(isyaij),nbas(isydel),nrhf(isymk),
     &               ONE,XDENS(koffden),ntotaij,XCMO(kcmoffk),
     &               NBASDEL,ONE,WORK(KOFFRES),ntotaij)

      end do
C---------------------------------------------------------------------
      if (locdbg) then
         write(lupri,*) ' Backtransform aijdel after loop '
         call flshfo(lupri)
      end if
*
* Write on disk the COMPLETE (NTOTOC(ISYM0)) array AIJdelta
* (tutti i delta contemporaneamente e non uno alla volta)???
*
      KOFF   = KDENSBT
      DO ISYMD = 1, NSYM
         ISYAIJ  = ISYMD
         LENAIJD = NCKI(ISYAIJ)*NBAS(ISYMD)
         IOFF = ICKID(ISYAIJ,ISYMD) + 1
         IF (LENAIJD.GT.0) THEN
            CALL PUTWA2(LUDENSTY,FNDENSTY,WORK(KOFF),IOFF,LENAIJD)
            IF (LOCDBG) THEN
               XNORM = DDOT(LENAIJD,WORK(KOFF),1,WORK(KOFF),1)
               WRITE (LUPRI,*) 'ISYMD:',ISYMD, 'LENAIJD ', LENAIJD
               WRITE (LUPRI,*) 'IOFF', IOFF, 'KOFF ', KOFF
               WRITE (LUPRI,*) '(T) density AIJ;del ', XNORM
               call flshfo(lupri) 
            END IF
         END IF
         KOFF = KOFF + LENAIJD
      END DO
C-------------------------------------------------------------------- 
      RETURN
      END


C=====================================================================
C  /* Deck btriajdel */
      SUBROUTINE  BTRIAJDEL(IOPT,XCMO,XDENS1,XDENS2,ISYDEN,
     &                      LUDENSTY,FNDENSTY,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: backtransform the last index of the (T) densities d_iajk(aijk)
C              and d_iajb(aijb) to AO basis and add together to give the 
C              d_iajdelta density. Dump result on file.
C     S. Coriani, February 2002
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER LWORK, IOPT
      INTEGER ISYRES, KSTART, ISYDEN
      INTEGER LWRK1,KEND1
      INTEGER IOFF, ISYMA, ISYMI, ISYMJ, ISYMK, ISYDEL
      INTEGER ISYKJI,ISYMKJ,KIAJK, KIAKJ,KIAJKS
      integer iii,nbasdel,kend2,lwrk2
      integer KAIJK, KAIKJ, KAIJKS
      integer ISYMIJ, ISYIJK, KIJKA, KIJAK, KIJKAS
      integer ISYIAJ,KCMOFFK
      integer LENGTH, KAIB, KBIA, ISYM, KAIBS
      integer KDEN4, KDEN3, KCIBS1,KCIBS2
      integer ISYBIA, ISYAIB, KAIB_C, KBIA_C, KAIBS_C, ISYMB
      integer ISYBIC, ISYCIB, KCIB_A, KBIC_A, KCIBS_A, KCIBS
      integer KOFF2, KOFF3, KOFF1, ISYMAI, ISYMC, ISYMBI
      integer ISYMCI, isycmo, KDENSBT, koff, ISYMD
      integer kcmo,kcmoffj,kcmoffa,koffden,koffres,ntotijk
      integer kcmoffb, ntotiak
      integer ntotiaj
      integer LUDENSTY, LENIJKD, LENIAJD

#if defined (SYS_CRAY)
      REAL WORK(LWORK),XCMO(*),XDENS1(*),XDENS2(*)
      REAL ZERO, ONE, TWO
      real ddot,xnorm2,xnorm
#else
      DOUBLE PRECISION WORK(LWORK),XCMO(*),XDENS1(*),XDENS2(*)
      DOUBLE PRECISION ZERO, ONE, TWO
      double precision ddot,xnorm2,xnorm
#endif


      CHARACTER FNDENSTY*9

      logical locdbg
      parameter (locdbg = .false.)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C
C#include <leinf.h>
C
      if (locdbg) then
        write(lupri,*) ' Backtransform (T) densities iajk and iajb'
        write(lupri,*) ' LUDENSTY:',LUDENSTY,' FNDENSTY:', FNDENSTY
      end if

c---------------------------------------------------------
c Backtransform ds_iajk (aij,k) = sum_k ds_iaj,k CMO_del,k
c DEVO FARE UN LOOP SU ISYMD???????????????????????????
c Dipende da come li voglio prelevare dal file.....
c Backtransform: 
C d_iajdel (aij,del) = sum_k d_iajk(aij,k) CMO_del,k
C                     + sum_b d_iajb(aij,b) CMO_del,b
c---------------------------------------------------------

      KEND1 = 1
      KDENSBT = KEND1
      KEND2   = KDENSBT + NTOTOC(ISYDEN)
      LWRK2   = LWORK - KEND2

      ISYCMO = 1

      CALL DZERO(WORK(KDENSBT),NTOTOC(ISYDEN)) 

      do isyiaj = 1, nsym
         isymk  = muld2h(isyden,isyiaj)
         isydel = muld2h(isycmo,isymk) 

         kcmoffk = IGLMRH(ISYDEL,ISYMK)  + 1
         koffden = ISAIKJ(ISYIAJ,ISYMK)  + 1
         koffres = ICKID(ISYIAJ,ISYDEL)  + KDENSBT
 
         nbasdel = MAX(NBAS(ISYDEL),1)
         ntotiaj = MAX(NCKI(ISYIAJ),1)

         call dgemm('N','T',ncki(isyiaj),nbas(isydel),nrhf(isymk),
     &               ONE,XDENS1(koffden),ntotiaj,XCMO(kcmoffk),
     &               NBASDEL,ONE,WORK(KOFFRES),ntotiaj)

      end do
C---------------------------------------------------------------------
      IF (IOPT.EQ.2) THEN
         do isyiaj = 1, nsym
            isymb  = muld2h(isyden,isyiaj)
            isydel = muld2h(isycmo,isymb) 

            kcmoffb = IGLMVI(ISYDEL,ISYMB)  + 1
            koffden = IT2SP(ISYIAJ,ISYMB)  + 1
            koffres = ICKID(ISYIAJ,ISYDEL) + KDENSBT
 
            nbasdel = MAX(NBAS(ISYDEL),1)
            ntotiaj = MAX(NCKI(ISYIAJ),1)

            call dgemm('N','T',ncki(isyiaj),nbas(isydel),
     &               nvir(isymb),
     &               ONE,XDENS2(koffden),ntotiaj,XCMO(kcmoffb),
     &               NBASDEL,ONE,WORK(KOFFRES),ntotiaj)

         end do
      END IF
C---------------------------------------------------------------------
      if (locdbg) then
         write(lupri,*) ' Backtransform1 after loop '
         call flshfo(lupri)
      end if
!
! Write on disk the COMPLETE (NTOTOC(ISYM0)) array IAJdelta/AIJdelta
! (tutti i delta contemporaneamente e non uno alla volta)???
!
      KOFF   = KDENSBT
      DO ISYMD = 1, NSYM
         ISYIAJ  = ISYMD
         LENIAJD = NCKI(ISYIAJ)*NBAS(ISYMD)
         IOFF = ICKID(ISYIAJ,ISYMD) + 1
         IF (LENIAJD.GT.0) THEN
            CALL PUTWA2(LUDENSTY,FNDENSTY,WORK(KOFF),IOFF,LENIAJD)
            IF (LOCDBG) THEN
               XNORM = DDOT(LENIAJD,WORK(KOFF),1,WORK(KOFF),1)
               WRITE (LUPRI,*) 'ISYMD:',ISYMD, 'LENIAJD : ', LENIAJD
               WRITE (LUPRI,*) 'IOFF', IOFF, 'KOFF ', KOFF
               WRITE (LUPRI,*) '(T) density iaj,del : ', XNORM
            END IF
         END IF
         KOFF = KOFF + LENIAJD
      END DO
C-------------------------------------------------------------------- 
      RETURN
      END


C=====================================================================
C  /* Deck btrijkdel */
      SUBROUTINE  BTRIJKDEL(XCMO,ISYCMO,XDENS,ISYDEN,
     &                      LUDENSTY,FNDENSTY,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: backtranform the last index of the symmetrized (T) 
C              density occ.occ,occ;vir to AO for later use in CCSD(T) 
C              FOPs and gradient
C     S. Coriani, December 2001
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER LWORK,IOPT
      INTEGER ISYRES, KSTART, ISYDEN
      INTEGER KOCC1,KOCC2,KOCCS1,KOCCS2,LWRK1,KEND1
      INTEGER IOFF, ISYMA, ISYMI, ISYMJ, ISYMK, ISYDEL
      INTEGER ISYKJI,ISYMKJ,KIAJK, KIAKJ,KIAJKS
      integer iii,nbasdel,kend2,lwrk2
      integer KAIJK, KAIKJ, KAIJKS
      integer ISYMIJ, ISYIJK, KIJKA, KIJAK, KIJKAS
      integer ISYIAJ,KCMOFFK
      integer LENGTH, KAIB, KBIA, ISYM, KAIBS
      integer KDEN4, KDEN3, KCIBS1,KCIBS2
      integer ISYBIA, ISYAIB, KAIB_C, KBIA_C, KAIBS_C, ISYMB
      integer ISYBIC, ISYCIB, KCIB_A, KBIC_A, KCIBS_A, KCIBS
      integer KOFF2, KOFF3, KOFF1, ISYMAI, ISYMC, ISYMBI
      integer ISYMCI, isycmo, KDENSBT, koff, ISYMD
      integer kcmo,kcmoffj,kcmoffa,koffden,koffres,ntotijk
      integer ntotiaj
      integer LUDENSTY, LENIJKD, LENIAJD

#if defined (SYS_CRAY)
      REAL WORK(LWORK),XCMO(*),XDENS(*)
      REAL ZERO, ONE, TWO
      real ddot,xnorm2,xnorm
#else
      DOUBLE PRECISION WORK(LWORK),XCMO(*),XDENS(*)
      DOUBLE PRECISION ZERO, ONE, TWO
      double precision ddot,xnorm2,xnorm
#endif


      CHARACTER FNDENSTY*9

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      if (locdbg) then
        write(lupri,*) ' Backtransform CCSD(T) density ijk,del'
        write(lupri,*) ' LUDENSTY:',LUDENSTY,' FNDENSTY:', FNDENSTY
      end if

C
C--------------------------------------------------------------
C     Read in symmetrized densities from file and backfransform
C--------------------------------------------------------------
!
!ds_ijka is stored ijk,a --> dijk,delta = sum_a ds_ijk,a C_del,a  
!
c---------------------------------------------------------
c Backtransform ds_ijka (ijk,a) = sum_a ds_ijk,a CMO_del,a
c E SE FACESSI UN LOOP SU ISYMD???????????????????????????
c FORSE E' PROPRIO QUELLO CHE DEVO FARE..............
c---------------------------------------------------------

      KDENSBT = 1
      KEND1   = KDENSBT + N3ODEL(ISYDEN)
      LWRK1   = LWORK - KEND1

      CALL DZERO(WORK(KDENSBT),N3ODEL(ISYDEN))

      do isyijk = 1, nsym
         isyma  = muld2h(isyden,isyijk)
         isydel = muld2h(isycmo,isyma) 

         kcmoffa = IGLMVI(ISYDEL,ISYMA)  + 1
         koffden = I3OVIR(ISYIJK,ISYMA)  + 1
         koffres = I3ODEL(ISYIJK,ISYDEL) + KDENSBT
 
         nbasdel = MAX(NBAS(ISYDEL),1)
         ntotijk = MAX(NMAIJK(ISYIJK),1)

!         if (locdbg) then
!            write(lupri,*) ' Backtransform1 before DGEMM '
!            write(lupri,*) 'kcmoffa, koffden, koffres',
!     &                        kcmoffa, koffden, koffres
!            write(lupri,*) 'isyijk, isyma, isydel',
!     &                        isyijk, isyma, isydel
!            write(lupri,*) 'nbasdel, ntotijk',
!     &                        nbasdel, ntotijk
!            call flshfo(lupri)
!         end if
 
         CALL DGEMM('N','T',nmaijk(isyijk),nbas(isydel),nvir(isyma),
     &               ONE,XDENS(koffden),ntotijk,XCMO(kcmoffa),
     &               NBASDEL,ONE,WORK(KOFFRES),ntotijk)

      end do
!
! Write on disk the COMPLETE (N3ODEL(ISYM0)) array IJKdelta
! (tutti i delta contemporaneamente e non uno alla volta)
! E se avessi scritto direttamente N3ODEL(ISYRES?)
!
!

      xnorm = 
     &    ddot(n3odel(isyden),work(kdensbt),1,work(kdensbt),1)
      if (locdbg) then
        write(lupri,*) 'Total norm Ds(ijkdelta) ', xnorm
      end if

      KOFF   = KDENSBT
      DO ISYMD = 1, NSYM
         ISYIJK  = ISYMD
         LENIJKD = NMAIJK(ISYIJK)*NBAS(ISYMD)
         IOFF = I3ODEL(ISYIJK,ISYMD) + 1
         IF (LENIJKD.GT.0) THEN
            CALL PUTWA2(LUDENSTY,FNDENSTY,WORK(KOFF),IOFF,LENIJKD)
            IF (LOCDBG) THEN
               XNORM = DDOT(LENIJKD,WORK(KOFF),1,WORK(KOFF),1)
               WRITE (LUPRI,*) 'ISYMD:',ISYMD, 'Norm XIJKdel:',XNORM
            END IF
         END IF
         KOFF = KOFF + LENIJKD
      END DO
C
      RETURN
      END
