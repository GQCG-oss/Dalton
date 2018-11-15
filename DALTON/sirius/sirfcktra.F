#define FCKTRA_DEBUG -1
!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!

!  FILE: DALTON/sirius/sirfcktra.F
!  (c) Copyright Hans Joergen Aa. Jensen, hjj@sdu.dk (2016)

! ===================================================================
!     sirfcktra.F section 0: Sirius interface
! ===================================================================
      subroutine SIR_INTOPEN_FCKTRA()
!
!  Written by Hans Joergen Aa. Jensen October 2014
!
      implicit none
#include "priunit.h"

#if FCKTRA_DEBUG > 0
      write(LUPRI,*) 'SIR_INTOPEN_FCKTRA called'
#endif

      end

      subroutine SIR_FCKTRA_CTL(TYPE_TEXT, ITRSIR, THR_MOTWOINT, CMO,
     &   WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen Apr 2017
!
      implicit none
      integer ITRSIR, LWORK
      real*8  THR_MOTWOINT, CMO(*), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  gnrinf.h : PARCAL
!  inforb.h : NNORBX, NORBT, ???
!  inftra.h : IPRTRA, LBUF
!  inftap.h : LUINTM, LBINTM
!  ccom.h   : THRS
#include "iratdef.h"
#include "maxaqn.h"
#include "thrzer.h"
#include "gnrinf.h"
#include "inforb.h"
#include "inftra.h"
#include "inftap.h"
#include "cbihr2.h"
#include "ccom.h"

      character*4 TYPE_TEXT
      integer JTRSIR, NTRLVL, MX_NDMAT
      integer KFREE, LFREE, KFRSAV, KUCMO, KDMAT, KFMAT, I, J
      integer MISH_TEST(8), MASH_TEST(8), IDAMAX, L_H2CD
      integer, allocatable :: ICDTRA(:), ITRTYP(:)
      character*8 TABLE(3), LAB123(3), INT_LABEL
      data        TABLE /'MUL_2INT','DRC_2INT','END_2INT'/
      data        LAB123/'********','********','********'/
      logical     FEXIST, OLDDX, FNDLAB

C --- start of executable code

      call QENTER('SIR_FCKTRA_CTL')

      JTRSIR = abs(ITRSIR)

#if FCKTRA_DEBUG > 0
      IPRTRA = MAX(IPRTRA,FCKTRA_DEBUG)
      write(lupri,*) 'DEBUG: IPRTRA set to',IPRTRA
#endif

      IF (IPRTRA .GT. 0) THEN
         CALL TITLER('.FCKTRA AO->MO integral transformation','*',200)
         WRITE(LUPRI,'(A,I4,1P,D10.2)')
     &   ' Sirius integral transformation level & screening threshold:',
     &   JTRSIR,THR_MOTWOINT
         FLUSH(LUPRI)
      END IF

C     Test if MO integrals already available

      if (TYPE_TEXT .EQ. 'COUL') then
         INT_LABEL = TABLE(1)
         L_H2CD = NNORBT
      else if (TYPE_TEXT .EQ. 'EXCH') then
         INT_LABEL = TABLE(2)
         L_H2CD = N2ORBT
      else
         call QUIT('Unknown TYPE_TEXT : '//TYPE_TEXT)
      end if
      LBUF = IRAT*L_H2CD

      REWIND(LUINTM)
      IF ( FNDLAB(INT_LABEL,LUINTM) ) THEN
         REWIND(LUINTM)
         READ  (LUINTM)
         READ  (LUINTM) J,NTRLVL
         IF (NTRLVL .EQ. 3 .AND. JTRSIR .EQ. 2) JTRSIR = -1
         IF (NTRLVL .EQ. 2 .AND. JTRSIR .EQ. 3) JTRSIR =  3
         IF (NTRLVL .EQ. 5 .AND. JTRSIR .EQ. 6) JTRSIR =  6
         IF (NTRLVL .EQ. 6 .AND. JTRSIR .EQ. 4) JTRSIR = -1 ! (gg|ii) and (gi|gi) missing
         IF (NTRLVL .EQ. 6 .AND. JTRSIR .EQ. 5) JTRSIR = -1 ! (gg|gi) missing
         IF (NTRLVL .GE. JTRSIR) THEN
C           read CMO matrix from LUINTM
C           and subtract from input CMO matrix
            READ  (LUINTM) WORK(1:NCMOT)
            WORK(1:NCMOT) = WORK(1:NCMOT) - CMO(1:NCMOT)
            I = IDAMAX(NCMOT,WORK,1)
            READ (LUINTM) MISH_TEST(1:8), MASH_TEST(1:8)
            MISH_TEST(1:8) = MISH_TEST(1:8) - NISH(1:8)
            MASH_TEST(1:8) = MASH_TEST(1:8) - NASH(1:8)
            J = 0
            DO I = 1,8
               J = J + ABS(MISH_TEST(I)) + ABS(MASH_TEST(I))
            END DO
            IF (ABS(WORK(I)) .LE. THRZER .AND. J .EQ. 0) THEN
               IF (IPRTRA.GE.0) WRITE (LUPRI,'(/A)')
     &            ' SIR_FCKTRA_CTL: Integral transformation abandoned,'
     &          //' the required MO integrals are already available.'
               GO TO 9999
            END IF
         END IF
      END IF

      IF (THR_MOTWOINT .LT. THRS) THEN
         WRITE(LUPRI,'(/A,1P,D10.2,A,D10.2)')
     &   '.FCKTRA INFO: AO screening threshold raised from',
     &   THR_MOTWOINT,' to the general integral threshold',THRS
         THR_MOTWOINT = THRS
         FLUSH(LUPRI)
      END IF

      allocate(ICDTRA(NNORBX))

      KFREE  = 1
      LFREE  = LWORK
      KFRSAV = KFREE

      ! call SIR_INTOPEN_FCKTRA()

!     translate sirius integral level code to FCKTRA and NEWTRA integral level code
      call N_TRALVL(JTRSIR, NTRLVL) 
      IF (JTRSIR .GE. 1 .AND. JTRSIR .LE. 4) NTRLVL = NTRLVL + 1
      ! we need third order, because USEDRC not defined for FCKTRA
#if FCKTRA_DEBUG > 0
      write(LUPRI,*) 'DEBUG: SIR_FCKTRA_CTL called, debug level',
     &   FCKTRA_DEBUG
      write(LUPRI,*) 'DEBUG: ITRSIR, NTRLVL, LWORK:',ITRSIR,NTRLVL,LWORK
      flush(LUPRI)
#endif

!     call N_TRASET(sirntra.F) in order to be able to call N_TRALIM(sirntra.F)
!     which sets :
!        ICDTRA: index array for C,D distributions (**|CD) to calculate
!        ITRTYP(1:NORBT) = number of integral indices in which this orbital
!                          enters (i.e. 0,1,2,3, or 4)
      call N_TRASET(-1,LWORK)
#if FCKTRA_DEBUG > 0
      IPRTRA = MAX(IPRTRA,FCKTRA_DEBUG)
      write(lupri,*) 'DEBUG: IPRTRA set to',IPRTRA
#endif
      allocate(ITRTYP(NORBT))
      call N_TRALIM(NTRLVL,ICDTRA,ITRTYP)
      deallocate ( ITRTYP )

      REWIND(LUINTM)
      WRITE (LUINTM) (0, I=1,8) ! info not used for FCKTRA
      WRITE (LUINTM) LBINTM,JTRSIR,NSYM,NORB,NBAS
      WRITE (LUINTM) CMO(1:NCMOT)
      WRITE (LUINTM) NISH(1:8), NASH(1:8)

      call GETDAT(LAB123(2),LAB123(3))   ! place date in LAB123(2) and time in LAB123(3)
      WRITE (LUINTM) LAB123,INT_LABEL
      WRITE (LUINTM) NNORBX, ICDTRA
!
!...  Open direct access file for final MO integrals
!     (the ICDTRA array which links a C,D to a record has been saved on LUINTM)
!
      if (TYPE_TEXT .eq. 'COUL' ) then

         CALL GPINQ('MO2INT_FCKTRA','EXIST',FEXIST)
         IF (FEXIST) THEN
            IF (LUMINT .LE. 0) CALL GPOPEN(LUMINT,
     &         'MO2INT_FCKTRA','OLD','DIRECT',' ',LBUF,OLDDX)
            CALL GPCLOSE(LUMINT,'DELETE')
         END IF
         LUMINT = -1
         CALL GPOPEN(LUMINT,'MO2INT_FCKTRA','NEW',
     &      'DIRECT','UNFORMATTED',LBUF,OLDDX)
      else if (TYPE_TEXT .eq. 'EXCH') then

         CALL GPINQ('DRC2INT_FCKTRA','EXIST',FEXIST)
         IF (FEXIST) THEN
            IF (LUMINT .LE. 0) CALL GPOPEN(LUMINT,
     &         'DRC2INT_FCKTRA','OLD','DIRECT',' ',LBUF,OLDDX)
            CALL GPCLOSE(LUMINT,'DELETE')
         END IF
         LUMINT = -1
         CALL GPOPEN(LUMINT,'DRC2INT_FCKTRA','NEW',
     &      'DIRECT','UNFORMATTED',LBUF,OLDDX)
      else
         call QUIT('Unknown TYPE_TEXT : '//TYPE_TEXT)
      end if


      CALL FCKTRA_CTL(TYPE_TEXT, ICDTRA, THR_MOTWOINT, CMO, WORK, LWORK)


      WRITE (LUINTM) LAB123,TABLE(3)
      REWIND(LUINTM)
      CALL GPCLOSE(LUINTM, 'KEEP')

      deallocate ( ICDTRA )
 9999 call QEXIT('SIR_FCKTRA_CTL')
      return
      end

! ===================================================================
!     sirfcktra.F section 1: calculate MO integrals
! ===================================================================
      subroutine FCKTRA_CTL( TYPE_TEXT, ICDTRA, THR_MOTWOINT, CMO,
     &   WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      character*4 TYPE_TEXT
      integer   ICDTRA(1:NNORBX), LWORK
      real*8    THR_MOTWOINT, CMO(*), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  gnrinf.h : PARCAL
!  inforb.h : NNORBX, NORBT, ???
!  inftra.h : IPRTRA, LBUF
!  inftap.h : LUSUPM
!  cbihr2.h : IFTHRS
!  ccom.h   : THRS
#include "iratdef.h"
#include "maxaqn.h"
#include "thrzer.h"
#include "gnrinf.h"
#include "inforb.h"
#include "inftra.h"
#include "inftap.h"
#include "cbihr2.h"
#include "ccom.h"

      integer MX_NDMAT, LUSUPM_save, IFTHRS_save
      integer KFREE, LFREE, KFRSAV, KUCMO, KDMAT, KFMAT, I, J
      integer MISH_TEST(8), MASH_TEST(8), IDAMAX
      integer, allocatable :: ICD_NODE(:)
      character*8 TABLE(2),LAB123(3)
      data        TABLE /'MUL_2INT','END INTM'/
      data        LAB123/'********','********','********'/
      logical     FEXIST, OLDDX, FNDLAB

C --- start of executable code

      call QENTER('FCKTRA_CTL')

#if FCKTRA_DEBUG > 0
      IPRTRA = MAX(IPRTRA,FCKTRA_DEBUG)
      write(lupri,*) 'FCKTRA_DEBUG: IPRTRA set to',IPRTRA
      write(lupri,*) 'FCKTRA_DEBUG: THR_MOTWOINT ',THR_MOTWOINT
      flush(lupri)
#endif

      IF (IPRTRA .GT. 0) THEN
!        WRITE(LUPRI,'(A,1P,D10.2)')
!    &   ' Transformation screening threshold:', THR_MOTWOINT
         FLUSH(LUPRI)
      END IF

      IFTHRS_save = IFTHRS
      IFTHRS = NINT(-LOG10(THR_MOTWOINT))


      KFREE  = 1
      LFREE  = LWORK
      KFRSAV = KFREE


      if (FCKTRA_TYPE .eq. 1 .AND. PARCAL) then

#ifdef VAR_MPI
         CALL FCKTRA_DISTRIBUTED_MASTER( TYPE_TEXT, ICDTRA, IFTHRS,
     &      CMO, WORK, LWORK)

#else
         CALL QUIT('Parallel FCKTRA only possible if compiled with MPI')
#endif
      else ! not PARCAL

         call MEMGET2('REAL','UCMO',  KUCMO,NBAST*NORBT,
     &      WORK,KFREE,LFREE)

         MX_NDMAT = LFREE/(4*N2BASX) ! use half of free memory for DMAT and FMAT
         MX_NDMAT = MIN(NNORBX,MX_NDMAT)
         call MEMGET2('REAL','DMAT',  KDMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)
         call MEMGET2('REAL','FMAT',  KFMAT,MX_NDMAT*N2BASX,
     &      WORK,KFREE,LFREE)

         call UPKCMO(CMO,WORK(KUCMO))

         allocate(ICD_NODE(NNORBX))
         ICD_NODE(:) = 0

         LUSUPM_save = LUSUPM
         if (FCKTRA_TYPE .eq. 3) LUSUPM = -9999 ! use RDSUPM
         if (FCKTRA_TYPE .eq. 4) LUSUPM = -1    ! do not use RDSUPM

         call FCKTRA_FCK_DISTRIBUTED(TYPE_TEXT, ICDTRA, ICD_NODE,
     &      WORK(KUCMO), MX_NDMAT,
     &      WORK(KDMAT), WORK(KFMAT), WORK(KFREE), LFREE)

         LUSUPM = LUSUPM_save
         deallocate ( ICD_NODE )

         call MEMREL('after FCKTRA_FCK_DISTRIBUTED',WORK,
     &      1,KFRSAV,KFREE,LFREE)

      end if !  if (FCKTRA_TYPE .eq. 1 .AND. PARCAL) then ... else ...

      IFTHRS = IFTHRS_save
 9999 call QEXIT('FCKTRA_CTL')
      return
      end

      subroutine FCKTRA_FCK_DISTRIBUTED( TYPE_TEXT, ICDTRA,
     &   ICD_NODE, UCMO, MX_NDMAT, DMAT, FMAT, WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      character*4 TYPE_TEXT
      integer ICDTRA(NNORBX), ICD_NODE(NNORBX)
      integer MX_NDMAT, LWORK
      real*8  DMAT(N2BASX,MX_NDMAT), FMAT(N2BASX,MX_NDMAT)
      real*8  UCMO(NBAST,NORBT), WORK(LWORK)
#include "priunit.h"

! Used from include files:
!  infpri.h : IPRFCK
!  inftra.h : IPRTRA
!  infinp.h : DIRFCK, ADDSRI
!  inforb.h : N2BASX, NORBT,NNORBX, ...
!  infind.h : ISMO(:)
!  inftap.h : LUMINT
!  dftcom.h : HFXFAC
!  cbihrs.h : NOSSUP, ONLY_J
!  infpar.h : MYNUM, MASTER
!  mtags.h  : MTAG7
#include "maxorb.h"
#include "maxash.h"
#include "infpri.h"
#include "inftra.h"
#include "infinp.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "dftcom.h"
#include "cbihrs.h"
#include "infpar.h"
#include "mtags.h"

      integer ICD, IC, ID, NDMAT, L_H2CD
      integer, allocatable :: ICD_REC(:), IFCTYP(:), ISYMDM(:)
      logical FEXIST, OLDDX, FCKTRA_TIMING
      real*8  TCPU0,TCPU1,TCPU2,TWAL0,TWAL1,TWAL2
      real*8  TCPU_FCK,TCPU_WR,TWAL_FCK,TWAL_WR

      logical NOSSUP_SAVE, ONLY_J_SAVE, ADDSRI_save
      logical ABCD_EQ_BACD
      real*8  HFXFAC_SAVE

      allocate ( ICD_REC(MX_NDMAT) )
      allocate ( IFCTYP(MX_NDMAT) )
      allocate ( ISYMDM(MX_NDMAT) )

      FCKTRA_TIMING = IPRTRA .gt. 0  .or. FCKTRA_DEBUG .ge. 0
      IF (FCKTRA_TIMING) THEN
         call GETTIM(TCPU0,TWAL0)
         TCPU_FCK = 0.0D0
         TCPU_WR  = 0.0D0
         TWAL_FCK = 0.0D0
         TWAL_WR  = 0.0D0
      END IF
#if FCKTRA_DEBUG > 0
      write(LUPRI,*) 'FCKTRA_FCK_DISTRIBUTED called'
      write(LUPRI,*) 'MYNUM, MX_NDMAT, LWORK:',MYNUM,MX_NDMAT,LWORK
#if FCKTRA_DEBUG > 5
      write(LUPRI,*) 'Test output of UCMO:'
      call OUTPUT(UCMO,1,NBAST,1,NORBT,NBAST,NORBT,-1,LUPRI)
#endif
      flush(lupri)
#endif

!     if srdft, only long range integrals
      ADDSRI_save = ADDSRI
      ADDSRI = .FALSE.

!     for RDSUPM calls
      NOSSUP_SAVE = NOSSUP
      NOSSUP      = .true.
      ONLY_J_SAVE = ONLY_J

      HFXFAC_SAVE = HFXFAC

      if (TYPE_TEXT(1:4) .eq. 'COUL') then ! Coulomb/Mulliken integrals
         HFXFAC      = 0.0D0
         ONLY_J      = .true.
         ABCD_EQ_BACD= .true.
         IFCTYP(1:MX_NDMAT) = 11 ! symmetric density matrix, only Coulomb
      else if (TYPE_TEXT(1:4) .eq. 'EXCH') then ! Exchange/Dirac integrals
         HFXFAC      = 2.0D0
         ABCD_EQ_BACD= .false.
         IFCTYP(1:MX_NDMAT) =  2 ! non-symmetric density matrix, only exchange
         if (LUSUPM .ne. -1) then
            call quit(
     &      'ERROR: FCKTRA "EXCH" integrals not yet for "SUPMAT"')
         end if
      else if (TYPE_TEXT(1:4) .eq. 'ANTI') then ! antisymmetrized exchange integrals
         HFXFAC      = 2.0D0
         ABCD_EQ_BACD= .false.
         IFCTYP(1:MX_NDMAT) =  3 ! non-symmetric density matrix, Coulomb and exchange
         CALL QUIT('FCKTRA ERROR: '//
     &   'Antisymmetrized Dirac integrals not implemented yet.')
         if (LUSUPM .ne. -1) then
            call quit(
     &      'ERROR: FCKTRA "EXCH" integrals not yet for "SUPMAT"')
         end if
      else
         write(LUPRI,'(///2A)') 'FCKTRA ERROR: '//
     &   'Unknown MO integral type: ',TYPE_TEXT
         CALL QUIT('FCKTRA ERROR: '//
     &   'Requested MO integral type not implemented.')
      end if


      NDMAT = 0
      ICD   = 0
      do IC = 1, NORBT
         do ID = 1, IC
            ICD = ICD + 1
#if FCKTRA_DEBUG > 4
            write (lupri,*) 'mynum, icd, icdtra, icd_node',
     &      mynum, icd, icdtra(icd), icd_node(icd) ; flush(lupri)
#endif

            if (ICDTRA(ICD) .gt. 0 .and.
     &          ICD_NODE(ICD) .eq. MYNUM) then

               NDMAT = NDMAT + 1
               ISYMDM(NDMAT)  = MULD2H(ISMO(IC), ISMO(ID))
               ICD_REC(NDMAT) = ICDTRA(ICD)
               call FCKTRA_CD_DMAT(UCMO(1,IC),UCMO(1,ID),
     &            DMAT(1,NDMAT),ABCD_EQ_BACD)
#if FCKTRA_DEBUG > 2
               write(lupri,*)
     &           'mynum, NDMAT, IC, ID, ICD, ICDTRA, ISYMDM:',
     &            mynum,NDMAT,IC,ID,ICD,ICDTRA(ICD), ISYMDM(NDMAT)
#if FCKTRA_DEBUG > 5
               write(lupri,*) 'DMAT no. NDMAT'
               call output(DMAT(1,NDMAT),1,NBAST,1,NBAST,
     &            NBAST,NBAST,-1,LUPRI)
#endif
               flush(lupri)
#endif
            end if

            if (NDMAT .eq. MX_NDMAT .or. 
     &         (NDMAT .gt. 0 .and. ICD .eq. NNORBX) ) then
               call DZERO(FMAT,NDMAT*N2BASX)
               if (FCKTRA_TIMING) call GETTIM(TCPU1,TWAL1)

               if (LUSUPM .ne. -1) then
#if FCKTRA_DEBUG > 1
                  write (lupri,*) 'calling RDSUPM, NDMAT',NDMAT
                  flush(lupri)
#endif
                  call RDSUPM(NDMAT,FMAT,DMAT,ISYMDM,WORK,LWORK)
               else
#if FCKTRA_DEBUG > 2
                  IPRFCK = MAX(IPRFCK,FCKTRA_DEBUG)
                  write (lupri,*)'calling SIRFCK; NDMAT, IPRFCK, DIRFCK'
     &            ,NDMAT, IPRFCK, DIRFCK
                  flush(lupri)
#endif
                  call SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRFCK,
     &                     WORK,LWORK)
               end if
               if (FCKTRA_TIMING) then
                  call GETTIM(TCPU2,TWAL2)
#if FCKTRA_DEBUG > 2
                  write (lupri,*) 'CPU and wall time',
     &               TCPU2-TCPU1,TWAL2-TWAL1
                  flush(lupri)
#endif
                  TCPU_FCK = TCPU_FCK + TCPU2 - TCPU1
                  TWAL_FCK = TWAL_FCK + TWAL2 - TWAL1
               end if

               ! (**|CD) in AO basis now in FMAT
               if (ABCD_EQ_BACD) then
                  call FCKTRA_AB_TO_MO(NDMAT,UCMO,FMAT,DMAT)
                  ! (**|CD) in MO basis now in DMAT
                  CALL FCKTRA_AB_PACK(NDMAT,DMAT,FMAT,ISYMDM)
                  ! packed (**|CD) in MO basis now in FMAT
                  L_H2CD = NNORBT
               else
                  call FCKTRA_AB_TO_MO(NDMAT,UCMO,FMAT,FMAT)
                  ! (**|CD) in MO basis now in FMAT
                  call quit('packing abcd .ne. bacd not implemented')
                  L_H2CD = N2ORBT
               end if
#if FCKTRA_DEBUG > 2
               write (lupri,*) 'Back from FCKTRA_AB_TO_MO'; flush(lupri)
#endif
               if (MYNUM .ne. 0) then
                  call MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG7)
                  call MPIXSEND(NDMAT,1,'INTEGER',MASTER,MTAG7)
                  call MPIXSEND(ICD_REC,NDMAT,'INTEGER',MASTER,MTAG7)
                  call MPIXSEND(FMAT,NDMAT*L_H2CD,'DOUBLE',MASTER,MTAG7)
#if FCKTRA_DEBUG > 1
                  write(lupri,*)'done to master, mtag7',master,mtag7
                  flush(lupri)
#endif
               else
#if FCKTRA_DEBUG > 1
                  write(lupri,*) 'writing to MOTWOINT, ndmat',ndmat
                  flush(lupri)
#endif
                  call FCKTRA_WR_LUMINT(LUMINT,NDMAT,FMAT,
     &               L_H2CD,ICD_REC)
               end if
               if (FCKTRA_TIMING) then
                  call GETTIM(TCPU1,TWAL1)
                  TCPU_WR  = TCPU_WR  + TCPU1 - TCPU2
                  TWAL_WR  = TWAL_WR  + TWAL1 - TWAL2
               end if
               NDMAT = 0
            end if
         end do ! ID = 1, IC
      end do ! IC = 1, NORBT

      if (FCKTRA_TIMING) then
         call GETTIM(TCPU1,TWAL1)
         write(lupri,'(A,2F20.3)') 'CPU and WALL time (s) AO Fock',
     &      TCPU_FCK,TWAL_FCK
         write(lupri,'(A,2F20.3)') 'CPU and WALL time (s) AO->MO ',
     &      TCPU_WR,TWAL_WR
         TCPU2 = TCPU1 - TCPU0 - TCPU_FCK - TCPU_WR
         TWAL2 = TWAL1 - TWAL0 - TWAL_FCK - TWAL_WR
         write(lupri,'(A,2F20.3)') 'CPU and WALL time (s) rest   ',
     &      TCPU2,TWAL2
         flush(lupri)
      end if

      HFXFAC = HFXFAC_SAVE
      NOSSUP = NOSSUP_SAVE
      ONLY_J = ONLY_J_SAVE
      ADDSRI = ADDSRI_save

      deallocate ( ICD_REC )
      deallocate ( IFCTYP )
      deallocate ( ISYMDM )

      return
      end
! ===================================================================
      subroutine FCKTRA_CD_DMAT(CMO_C,CMO_D,DMAT_CD,SYM_CD)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      real*8  CMO_C(NBAST), CMO_D(NBAST), DMAT_CD(NBAST,NBAST)
      logical SYM_CD

! Used from include files:
!  inforb.h : NBAST
#include "inforb.h"

      integer JA, JB

      do JB = 1, NBAST
         do JA = 1,NBAST
            DMAT_CD(JA,JB) = CMO_C(JA)*CMO_D(JB)
         end do
      end do

      if (SYM_CD) then ! symmetrize
         call DGETSI(NBAST,DMAT_CD,DMAT_CD)
      end if

      end
! ===================================================================
      subroutine FCKTRA_AB_TO_MO(NDMAT,UCMO,H2CD_AO,H2CD_MO)
!
!  Written by Hans Joergen Aa. Jensen 2016-2017
!
!  Transform (ab|CD) to (AB|CD); ABCD MO indices, ab AO indices
!
!  H2CD_AO contains NDMAT "CD" records of (ab|CD) on input,
!  H2CD_MO contains NDMAT "CD" records of (AB|CD) on output.
!
!  Note: H2CD_AO and H2CD_MO may be equivalenced
!        (i.e. share memory), as N2ORBX .le. N2BASX
!
!
      implicit none
      integer NDMAT
      real*8  UCMO(NBAST,NORBT)
      real*8  H2CD_AO(N2BASX,NDMAT), H2CD_MO(N2ORBX,NDMAT)
      real*8, allocatable :: H2CD_TMP(:)
#include "priunit.h"

! Used from include files:
!  inforb.h : NBAST,NORBT,N2BASX,N2ORBX, ...
#include "inforb.h"

      integer IDMAT

      allocate (H2CD_TMP(N2BASX))

! TO DO: this should be done on GPU if available /hjaaj Aug. 2016
      do IDMAT = 1, NDMAT
         call DGEMM('T','N',NORBT,NBAST,NBAST, 1.0D0,
     &                 UCMO,NBAST,
     &                 H2CD_AO(1,IDMAT),NBAST, 0.0D0,
     &                 H2CD_TMP,NORBT)
         call DGEMM('N','N',NORBT,NORBT,NBAST, 1.0D0,
     &                 H2CD_TMP,NORBT,
     &                 UCMO,NBAST, 0.0D0,
     &                 H2CD_MO(1,IDMAT),NORBT)
#if FCKTRA_DEBUG > 5
            write (LUPRI,*) 'H2CD_MO, record',IDMAT
            call OUTPUT(H2CD_MO(1,IDMAT),1,NORBT,1,NORBT,
     &         NORBT,NORBT,-1,LUPRI)
#endif
      end do

      deallocate(H2CD_TMP)

      end
! ===================================================================
      subroutine FCKTRA_AB_PACK(NDMAT,H2CD_MO,H2CD_MO_PK,ICDSYM)
!
!  Written by Hans Joergen Aa. Jensen 2016-2017
!
!  Purpose: Pack Mulliken MO integrals before writing them to disk,
!           using (AB|CD) = (BA|CD) and (AB|CD)=0 if sym(AB) .ne. sym(CD).
!
!  H2CD_MO contains NDMAT "CD" records of (AB|CD) on input,
!  packed into H2CD_MO_PK on output.
!
!
      implicit none
      integer NDMAT, ICDSYM(NDMAT)
      real*8  H2CD_MO(N2ORBX,NDMAT), H2CD_MO_PK(NNORBT,NDMAT)

#include "priunit.h"

! Used from include files:
!  inforb.h : N2ORBX,NNORBT, NSYM, ...
#include "inforb.h"

      integer IDMAT
      real*8, allocatable :: H2CD_TMP(:)

      if (NSYM > 1) allocate (H2CD_TMP(NNORBX))

      do IDMAT = 1, NDMAT
         if (NSYM > 1) then
            call DGETSP(NORBT, H2CD_MO(1,IDMAT), H2CD_TMP)
            call TRDPAK( H2CD_TMP, H2CD_MO_PK(1,IDMAT),
     &         NORB,IORB,NORBT,ICDSYM(IDMAT),1)
         else
            call DGETSP(NORBT, H2CD_MO(1,IDMAT), H2CD_MO_PK(1,IDMAT))
         end if
#if FCKTRA_DEBUG > 5
         write (LUPRI,*) 'IDMAT, sym',
     &      IDMAT,ICDSYM(IDMAT)
         if (ICDSYM(idmat) == 1) then
            write (LUPRI,*) 'Packed H2CD_MO'
            call OUTPKB(H2CD_MO_PK(1,IDMAT),NORB,NSYM,-1,LUPRI)
         end if
#endif
      end do

      if (NSYM > 1) deallocate(H2CD_TMP)

      end

! ===================================================================
      subroutine FCKTRA_WR_LUMINT(LUMINT,NDMAT,H2CD,L_H2CD,ICD_REC)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      integer LUMINT, NDMAT, L_H2CD, ICD_REC(NDMAT), IOS
      real*8  H2CD(L_H2CD,NDMAT)

#include "priunit.h"

      integer IDMAT

      do IDMAT = 1, NDMAT
         ! write NDMAT (**|CD) MO integral distributions to file
         write(LUMINT, rec=ICD_REC(IDMAT), iostat=IOS)
     &      H2CD(1:L_H2CD,IDMAT)
         if (IOS .ne. 0) then
            call QENTER('FCKTRA_WR_LUMINT')
            write(lupri,*) 'IDMAT, rec, IOS',IDMAT,ICD_REC(IDMAT),IOS
            write(lupri,*) 'LUMINT, L_H2CD',LUMINT,L_H2CD
            call QUIT('IOS .ne. 0')
         end if
      end do

      end

! ===================================================================
!     sirfcktra.F section 2: FCKTRA type 1 routines (if VAR_MPI)
! ===================================================================
#ifdef VAR_MPI
      subroutine FCKTRA_DISTRIBUTED_MASTER( TYPE_TEXT, ICDTRA, IFTHRS,
     &   CMO, WORK, LWORK)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      character*4 TYPE_TEXT
      integer ICDTRA(NNORBX), IFTHRS, LWORK
      real*8  CMO(NCMOT), WORK(LWORK)

#include "priunit.h"

#include "mpif.h"

! Used from include files:
!  gnrinf.h : CHIVAL
!  blocks.h : NSYMBL
!  inforb.h : N2BASX, NORBT,NNORBX, ...
!  infind.h : ISMO(:)
!  inftra.h : IPRTRA
!  inftap.h : LUMINT
!  infpar.h : MYNUM
!  iprtyp.h : defined parallel calculation types
!  mtags.h  : MTAG7

#include "maxorb.h"
#include "maxash.h"
#include "gnrinf.h"
#include "blocks.h"
#include "inforb.h"
#include "infind.h"
#include "inftra.h"
#include "inftap.h"
#include "infpar.h"
#include "iprtyp.h"
#include "mtags.h"

      integer, allocatable :: ICD_NODE(:)
      integer IPRTYP, IPRINT_slaves, ICD, NCDTRA, IWHO, NWHO
      integer KFREE, LFREE
      integer NDMAT, IDIST, N2GAB, L_H2CD

      call qenter('FCKTRA_DISTRIBUTED_MASTER')

      if (TYPE_TEXT .eq. 'COUL') then
         L_H2CD = NNORBT
      else if (TYPE_TEXT .eq. 'EXCH') then
         L_H2CD = N2ORBT
      else
         call QUIT('Unknown TYPE_TEXT '//TYPE_TEXT)
      end if

!     start nodes on the distributed 2-electron integral transformation task

      IPRTYP = CALL_FCKTRA_DISTRIBUTED ! get code number from iprtyp.h
      call MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      IPRINT_slaves = IPRTRA ! TODO set print level in input
      call MPIXBCAST(IPRINT_slaves,1,'INTEGER',MASTER)

!     Now nodes are in subroutine FCKTRA_DISTRIBUTED_NODE.
!     First make sure basis set information is known to nodes

      KFREE = 1
      LFREE = LWORK
      CALL HER_sendreceive_molinfo(IPRINT_slaves,WORK,KFREE,LFREE)

!     Transfer integral transformation information to nodes

      allocate( ICD_NODE(NNORBX) )
      ICD_NODE(:) = 0
      NCDTRA = 0
      DO ICD = 1,NNORBX
         IF (ICDTRA(ICD) .EQ. 0) CYCLE
         ICD_NODE(ICD) = MOD(NCDTRA,NODTOT) + 1
         NCDTRA = NCDTRA + 1
      END DO

#if FCKTRA_DEBUG > 2
      write (*,*) 'FCKTRA_DISTRIBUTED_MASTER has started slaves'
      write (lupri,*) 'FCKTRA_DISTRIBUTED_MASTER has started slaves'
      write (lupri,*) 'IPRTYP,IPRINT_slaves=',IPRTYP,IPRINT_slaves
      write (lupri,*) 'ICD_NODE is'
      write (lupri,'(10I4,5X,10I4)') ICD_NODE(1:NNORBX)
      flush(lupri)
#endif
      call MPIXBCAST(TYPE_TEXT,4,'CHARACTER',MASTER)
      call MPIXBCAST(IFTHRS,1,'INTEGER',MASTER)
      call MPIXBCAST(CHIVAL,1,'DOUBLE',MASTER)
      call MPIXBCAST(ICDTRA,NNORBX,'INTEGER',MASTER)
      call MPIXBCAST(ISMO,NORBT,'INTEGER',MASTER)
      call MPIXBCAST(ICD_NODE,NNORBX,'INTEGER',MASTER)
      call MPIXBCAST(CMO,NCMOT,'DOUBLE',MASTER)

      IDIST = 0
  100 CONTINUE

      IWHO = -1
      call MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG7)
      call MPIXRECV(NDMAT,1,'INTEGER',NWHO,MTAG7)
      IF (NDMAT .LE. 0) THEN
         write(LUPRI,*)
     &   'ERROR: master received ndmat.le.0, ndmat, nwho',ndmat,nwho
         call quit('NDMAT .le. 0')
      END IF
      ! we reuse ICD_NODE to store ICD_REC
      ! CALL MPIXRECV(ICD_REC,NDMAT,'INTEGER',NWHO,MTAG7)
      CALL MPIXRECV(ICD_NODE,NDMAT,'INTEGER',NWHO,MTAG7)
      IF (NDMAT*L_H2CD .GT. LWORK) CALL QUIT('NDMAT*L_H2CD .gt. lwork')
      CALL MPIXRECV(WORK,NDMAT*L_H2CD,'DOUBLE',NWHO,MTAG7)

      call FCKTRA_WR_LUMINT(LUMINT,NDMAT,WORK,L_H2CD,ICD_NODE)

#if FCKTRA_DEBUG > 1
      write (lupri,*) 'Received from nwho =',nwho
      write (lupri,*) 'Received and saved ndmat =',ndmat
      write(lupri,'(A,(20I5))') 'The CD records:',ICD_NODE(1:NDMAT)
#endif

      IDIST = IDIST + NDMAT
      IF (IDIST .GT. NCDTRA) THEN
         write(lupri,*) 'ERROR: IDIST .gt. NCDTRA',IDIST,NCDTRA
         call quit('IDIST .gt. NCDTRA, uha')
      else if (IDIST .LT. NCDTRA) then
         go to 100
      end if

      deallocate( ICD_NODE )
      call qexit('FCKTRA_DISTRIBUTED_MASTER')
      return
      end subroutine FCKTRA_DISTRIBUTED_MASTER
! ===================================================================
      subroutine FCKTRA_DISTRIBUTED_NODE(WORK,LWORK,IPRINT)
!
!  Written by Hans Joergen Aa. Jensen 2016
!
      implicit none
      integer LWORK, IPRINT
      real*8  WORK(LWORK)

      integer, allocatable :: ICDTRA(:), ICD_NODE(:)
      integer ICD, NCDTRA, MX_NDMAT
      integer KFREE, LFREE, KDMAT, KFMAT, KUCMO
      logical PARHER_save, SLAVE_save, DIRFCK_save
      character*4 TYPE_TEXT
#include "priunit.h"

! Used from include files:
!  gnrinf.h : PARCAL, CHIVAL
!  infinp.h : DIRFCK
!  infind.h : ISMO(:)
!  inforb.h : NORBT, NNORBX, NCMOT
!  infpar.h : MYNUM,SLAVE,PARHER,?
!  inftap.h : LUSUPM
!  cbihr2.h : IFTHRS
#include "maxorb.h"
#include "maxash.h"
#include "gnrinf.h"
#include "infinp.h"
#include "infind.h"
#include "inforb.h"
#include "infpar.h"
#include "inftap.h"
#include "cbihr2.h"

!  mpif.h   : for MPI
#include "mpif.h"

!     First make sure basis set information is known to this node

      KFREE = 1
      LFREE = LWORK
      CALL HER_sendreceive_molinfo(IPRINT,WORK,KFREE,LFREE)

#if FCKTRA_DEBUG > 1
      write (lupri,*)
     &'FCKTRA_DISTRIBUTED_NODE has started, MYNYM =',MYNUM
      write (lupri,*) 'NORBT, NNORBX, NCMOT',
     &                 NORBT, NNORBX, NCMOT
      flush(lupri)
#endif

!    We check that orbital information has been transferred to slaves,
!    and then we call SETORB to define inforb.h and infind.h (we need ISMO(:))

      IF (NORBT .LE. 0) THEN
         CALL QUIT(
     &   'ERROR: Orbital information not transferred to slaves.')
      END IF
      CALL SETORB

!    Some set up for this task (AO to MO transformation)

      KFREE = 1
      LFREE = LWORK

      SLAVE_save  = SLAVE
      SLAVE  = .FALSE. ! Be sure the Fock matrices are NOT constructed in parallel,

      PARHER_save = PARHER
      PARHER = .FALSE. ! be sure the Fock matrices are NOT constructed in parallel
                       ! we are parallelizing on the distributions

      DIRFCK_save = DIRFCK
      DIRFCK = .TRUE.  ! NEVER read from AOTWOINT when parallel.

      LUSUPM = -1      ! Do not use supermatrix from file AO2SUPINT.


      allocate( ICDTRA(NNORBX) )
      allocate( ICD_NODE(NNORBX) )

      call MEMGET2('REAL','UCMO',  KUCMO,NBAST*NORBT,
     &      WORK,KFREE,LFREE)

      call MPIXBCAST(TYPE_TEXT,4,'CHARACTER',MASTER)
      call MPIXBCAST(IFTHRS,1,'INTEGER',MASTER)            ! accuracy of MO integrals (screening)
      call MPIXBCAST(CHIVAL,1,'DOUBLE',MASTER)
      call MPIXBCAST(ICDTRA,NNORBX,'INTEGER',MASTER)
      call MPIXBCAST(ISMO,NORBT,'INTEGER',MASTER)
      call MPIXBCAST(ICD_NODE,NNORBX,'INTEGER',MASTER)

      call MPIXBCAST(WORK(KFREE),NCMOT,'DOUBLE',MASTER)
      call UPKCMO(WORK(KFREE),WORK(KUCMO))

      NCDTRA = 0
      DO ICD = 1,NNORBX
         IF (ICD_NODE(ICD) .EQ. MYNUM) NCDTRA = NCDTRA + 1
      END DO

      IF (NCDTRA .EQ. 0) THEN
      ! nothing to do for me this time ...
         GO TO 1000
      END IF

      MX_NDMAT = LFREE/(4*N2BASX) ! use up to half of free memory for DMAT and FMAT
      MX_NDMAT = MIN(NCDTRA,MX_NDMAT)

      IF (IPRINT .GT. 0) THEN
         write(lupri,'(A,2I12)') 'FCKTRA # of distributions and'//
     &   ' # of distributions per load',NCDTRA, MX_NDMAT
         flush(lupri)
      END IF

      IF (MX_NDMAT .EQ. 0) THEN
         CALL QUIT('MX_NDMAT = 0, too little memory')
      END IF
      call MEMGET2('REAL','DMAT', KDMAT,MX_NDMAT*N2BASX,
     &     WORK,KFREE,LFREE)
      call MEMGET2('REAL','FMAT', KFMAT,MX_NDMAT*N2BASX,
     &     WORK,KFREE,LFREE)

#if FCKTRA_DEBUG > 2
      write(lupri,*) 'NCDTRA,MX_NDMAT,NNORBX =',NCDTRA,MX_NDMAT,NNORBX
      write (lupri,'(10I4,5X,10I4)') ICD_NODE(1:NNORBX)
      write(lupri,*) 'chival',chival
      write (lupri,*) mynum,'node; ICD_NODE is'
      write (lupri,*) MYNUM,' is calling fcktra_fck_distributed'
      flush(lupri)
#endif

      CALL  FCKTRA_FCK_DISTRIBUTED( TYPE_TEXT,
     &      ICDTRA, ICD_NODE, WORK(KUCMO), MX_NDMAT,
     &      WORK(KDMAT), WORK(KFMAT), WORK(KFREE), LFREE)


 1000 CONTINUE

      deallocate( ICDTRA )
      deallocate( ICD_NODE )

      PARHER = PARHER_save
      SLAVE  = SLAVE_save
      DIRFCK = DIRFCK_save

      return
      end subroutine FCKTRA_DISTRIBUTED_NODE
#endif   /* VAR_MPI */

! ===================================================================
!     sirfcktra.F section 3: read MO integrals
! ===================================================================
      subroutine FCKTRA_NXTH2M(IC,ID,H2CD,NEEDTP,WORK,IDIST)
!
!  Written by Hans Joergen Aa. Jensen 2016
!  This version is interface routine for .FCKTRA integral transformation.
!
! Purpose:
!    Read next Mulliken two-electron integral distribution (**|cd)
!    where |cd) distribution is needed according to NEEDTP(ITYPCD)
!
! Usage:
! 1) Set IDIST = 0 before first call.
!    Do **NOT** change IDIST in calling routine
!    until last distribution has been read (signalled by IDIST .eq. -1)
! 2) if IDIST<0 on input, then IC, ID must be set to correspond to
!    the distribution wanted in H2CD
!

      implicit none
#include "priunit.h"

      integer  :: IC, ID, NEEDTP(-4:6), IDIST
      real*8   :: H2CD(NORBT,NORBT), WORK(*)

! Used from include files:
!  inforb.h : NORBT, NNORBT, NNORBX, ???
!  infind.h : ISMO(:)
!  inftra.h : LBUF
!  inftap.h : LUINTM, LUMINT
#include "maxorb.h"
#include "maxash.h"
#include "dummy.h"
#include "iratdef.h"
#include "inforb.h"
#include "infind.h"
#include "inftra.h"
#include "inftap.h"
#include "orbtypdef.h"

      integer       :: MMORBX, ICD, ICD_1, ICDSYM
      integer       :: ITYPC, ITYPD, ITYPCD
      integer, save :: ICD_SAVE
      integer, save, allocatable :: ICDTRA(:)
      logical       :: OLDDX

#if FCKTRA_DEBUG > 3
      write(LUPRI,*) 'FCKTRA_NXTH2M called, IC, ID, IDIST =',
     &   IC, ID, IDIST
      flush(LUPRI)
#endif

      if (IDIST .le. 0) then
         if (allocated(ICDTRA)) deallocate(ICDTRA)
         allocate(ICDTRA(1:NNORBX))
         if (LUINTM .le. 0) THEN
            CALL GPOPEN(LUINTM,FNINTM,'OLD',' ',
     &                  'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         LBUF = IRAT*NNORBT
         IF (LUMINT .LE. 0) THEN
            CALL GPOPEN(LUMINT,'MO2INT_FCKTRA','OLD','DIRECT',' ',
     &         LBUF,OLDDX)
         END IF
         rewind (LUINTM)
         call MOLLAB('MUL_2INT',LUINTM,LUPRI)
         read (LUINTM) MMORBX, ICDTRA(1:MMORBX)
         if (NNORBX .ne. MMORBX) then
            call QENTER('FCKTRA_NXTH2M')
            write (LUPRI,'(//A,2I10)')  'FATAL ERROR '//
     &      'NNORBX on MO2INT_FCKTRA .ne. NNORBX',MMORBX,NNORBX
            call QUIT(
     &      'NNORBX on MO2INT_FCKTRA .ne. NNORBX in FCKTRA_NXTH2M')
         end if
         ICD_SAVE = 0
      end if

      if (IDIST .lt. 0) then

         if (IC .le. 0 .or. IC .gt. NORBT .and.
     &       ID .le. 0 .or. ID .gt. NORBT) then
            write (LUPRI,*) 'FCKTRA_NXTH2M: invalid input IC, ID',IC,ID
            call quit('FCKTRA_NXTH2M: invalid input IC, ID')
         end if
         if (IC .ge. ID) then
            ICD = IC*(IC-1)/2 + ID
         else
            ICD = ID*(ID-1)/2 + IC
         end if
         IDIST = ICDTRA(ICD)
         if (IDIST .le. 0) then
            call QUIT('IC,ID record no. .le. 0 in FCKTRA_NXTH2M')
         end if
         deallocate(ICDTRA)

      else ! IDIST .ge. 0

  100    continue
         IDIST = -1

         ICD_1 = ICD_SAVE + 1
         do ICD = ICD_1,NNORBX
            if (ICDTRA(ICD) .gt. 0) then
               IDIST    = ICDTRA(ICD)
               ICD_SAVE = ICD
               exit
            end if
         end do
         if (IDIST .eq. -1) then
            ! finished
            deallocate(ICDTRA)
            return
         end if

         do IC = 1, NORBT
            ICD = IC*(IC+1)/2
            if (ICD .lt. ICD_SAVE) cycle
            ID = ICD_SAVE - ICD + IC
            exit
         end do
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         ITYPCD = IDBTYP(ITYPC,ITYPD)
         if ( NEEDTP(ITYPCD) .eq. 0 ) go to 100

      end if

         ICDSYM = MULD2H( ISMO(IC), ISMO(ID) )
         read ( LUMINT, rec=IDIST ) WORK(1:NNORBT)
         call TRDPAK(WORK(1+NNORBT),WORK,NORB,IORB,NORBT,
     &      ICDSYM,-1)
         call DSPTGE(NORBT,WORK(1+NNORBT),H2CD)

#if FCKTRA_DEBUG > 9
         write (LUPRI,*) 'FCKTRA_NXTH2M, rec,IC,ID',IDIST,IC,ID
         write (LUPRI,*) 'Unpacked H2CD, ICDSYM',ICDSYM
         call output(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         if (ICDSYM.eq.1) then ! only print some of them, of ICDSYM.eq.1
            write (LUPRI,*) 'Packed H2CDPK'
            call outpkb(WORK,NORB,NSYM,-1,LUPRI)
         end if
#endif

      return
      end

      subroutine FCKTRA_NXTH2D(IC,ID,H2CD,NEEDTP,IDIST)
!
!  Written by Hans Joergen Aa. Jensen 2016
!  This version is interface routine for .FCKTRA integral transformation.
!
! Purpose:
!    Read next Dirac two-electron integral distribution <**|cd> = (*c|*d)
!    where (cd) distribution is needed according to NEEDTP(ITYPCD)
!
!

      implicit none
#include "priunit.h"

      integer  IC, ID, NEEDTP(*), IDIST
      real*8   H2CD(NORBT,NORBT)

! Used from include files:
!  inforb.h : NORBT, ???
#include "inforb.h"

      call quit('ERROR: FCKTRA_NXTH2D is not implemented')


      end
! -- end of DALTON/sirius/sirfcktra.F --
