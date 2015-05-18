!...
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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

module qmcmm

   implicit none

   public read_pot_qmnpmm
   public getdim_relmat
   public getdim_mmmat
   public comp_relmat
   public read_relmat
   public write_relmat
   public comp_mmrelmat
   public comp_dampvmat

   private

contains


      SUBROUTINE READ_POT_QMNPMM()
!
! Purpose:
!   reads NP and MM subsystems geometry, total charges, and
!   force field data from POTENTIAL.INP file
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

#include "codata.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      real(8), parameter :: xfact = 1.0d0/xtang
      integer :: i, istart, j, joff, luqmnp, idummy
!
      CHARACTER UNITS*2, NPWORD*2, FFWORD*4, NPLBL(MXNPATM)*2
      CHARACTER MMWORD*2, MMLBL(MXMMATM)*2
!
!     Open POTENTIAL.INP file
      LUQMNP=-1
      CALL GPOPEN(LUQMNP,'POTENTIAL.INP','OLD',' ',                     &
             'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUQMNP)
!     Read geometry input units, number of NP and MM
!     separate subsystems
      READ(LUQMNP,*) UNITS, TNPBLK, TMMBLK
      CALL UPCASE(UNITS)
!     Check input consistency
      IF ((UNITS.NE.'AA').AND.(UNITS.NE.'AU')) THEN
         CALL QUIT('Unknown units in POTENTIAL.INP')
      ENDIF
      IF ((DONPSUB).AND.(TNPBLK.LT.1)) THEN
         CALL QUIT('NP system is missing in POTENTIAL.INP')
      END IF
      IF ((DOMMSUB).AND.(TMMBLK.LT.1)) THEN
         CALL QUIT('MM system is missing in POTENTIAL.INP')
      END IF
      IF (TNPBLK.GT.MAXBLK) THEN
        WRITE(LUPRI,'(/2X,A)')                                          &
         'Maximum number of NP subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',               &
         TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
      IF (TMMBLK.GT.MAXBLK) THEN
        WRITE(LUPRI,'(/2X,A)')                                          &
         'Maximum number of MM subsystems exceeded!'
        WRITE(LUPRI,'(2X,A,I3,2X,A,I3)') 'Input TNPBLK=',               &
         TNPBLK, 'Maximum allowed:', MAXBLK
        CALL QUIT('Increase MAXBLK in qmnpmm.h')
      END IF
!     Read NP subsystems data
      IF (DONPSUB) THEN
         ISTART = 0
         DO I=1,TNPBLK
!           Read NP subsystem header
            READ(LUQMNP,*) NPWORD, NPCHRG(I), NPATOM(I)
!           Check input consistency
            CALL UPCASE(NPWORD)
            IF (NPWORD.NE.'NP') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')                           &
         'Incorrect NP subsystem I=', I,                                &
         'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (NPATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')                         &
         'Incorrect number of atoms ', NPATOM(I),                       &
         'in NP subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Read NP subsystem data
            DO J=1,NPATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) NPLBL(JOFF), NPCORD(1,JOFF),              &
               NPCORD(2,JOFF), NPCORD(3,JOFF),                          &
               NPFTYP(JOFF)
            END DO
            ISTART = ISTART + NPATOM(I)
            TNPATM = TNPATM + NPATOM(I)
            IF (TNPATM.GT.MXNPATM) THEN
               WRITE(LUPRI,'(/2X,A)')                                   &
          'Maximum number of NP atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                         &
          'Input TNPATM=', TNPATM,                                      &
          'Maximum allowed:', MXNPATM
               CALL QUIT('Increase MXNPATM in qmnpmm.h')
            END IF
         END DO
!        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TNPATM,XFACT,NPCORD,1)
         END IF
!        Read force field data
         READ(LUQMNP,*) FFWORD, TNPFF
         CALL UPCASE(FFWORD)
!        Check NP force field data consistency
         IF (FFWORD.NE.'NPFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted NP force field ',       &
            'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
         IF (TNPFF.GT.MXNPFF) THEN
            WRITE(LUPRI,'(/2X,A)')                                      &
             'Maximum number of NP force field types exceeded!'
            WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                            &
             'Input TNPFF=', TNPFF,                                     &
             'Maximum allowed:', MXNPFF
            CALL QUIT('Increase MXNPFF in qmnpmm.h')
         END IF
!        Read force field data
         DO I=1,TNPFF
            READ(LUQMNP,*) NPFPOL(I), NPFCAP(I), NPFOMG1(I),            &
                      NPFGAM1(I), NPFOMG2(I), NPFGAM2(I),               &
                      NPFFAC(I)
         END DO
!        Check force field definitions
         DO I=1,TNPATM
            IF (NPFTYP(I).GT.TNPFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')                           &
          'Unknown force field requested for atom', I,                  &
          'in NP region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_NPREGION(NPLBL)
         END IF
      END IF
!     MM subsystems input
      IF (DOMMSUB) THEN
         ISTART = 0
         DO I=1,TMMBLK
!           Read MM subsystem header
            READ(LUQMNP,*) MMWORD,  MMATOM(I)
!           Check input consistency
            CALL UPCASE(MMWORD)
            IF (MMWORD.NE.'MM') THEN
               WRITE(LUPRI,'(/2X,A,I2,1X,A)')                           &
         'Incorrect MM subsystem I=', I,                                &
         'header in POTENTIAL.INP file!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
            IF (MMATOM(I).LE.0) THEN
               WRITE(LUPRI,'(/2X,A,I5,A,I2,A)')                         &
         'Incorrect number of atoms ', MMATOM(I),                       &
         'in MM subsystem', I, '!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Read MM subsystem data
            DO J=1,MMATOM(I)
               JOFF = ISTART + J
               READ(LUQMNP,*) MMLBL(JOFF), MMMOL(JOFF),                 &
                         mm_cord(1,JOFF), mm_cord(2,JOFF),                &
                         mm_cord(3,JOFF), MMFTYP(JOFF)
            END DO
            ISTART = ISTART + MMATOM(I)
            TMMATM = TMMATM + MMATOM(I)
            IF (TMMATM.GT.MXMMATM) THEN
               WRITE(LUPRI,'(/2X,A)')                                   &
          'Maximum number of MM atoms exceeded!'
               WRITE(LUPRI,'(2X,A,I3,2X,A,I3)')                         &
          'Input TMMATM=', TMMATM,                                      &
          'Maximum allowed:', MXMMATM
               CALL QUIT('Increase MXMMATM in qmnpmm.h')
            END IF
         END DO
!        Convert to atomic units if neeeded
         IF (UNITS.EQ.'AA') THEN
            CALL DSCAL(3*TMMATM,XFACT,mm_cord,1)
         END IF
!        Read force field data
         READ(LUQMNP,*) FFWORD, TMMFF
         CALL UPCASE(FFWORD)
!        Check NP force field data consistency
         IF (FFWORD.NE.'MMFF') THEN
            WRITE(LUPRI,'(/2X,A,A)') 'Corrupted MM force field ',       &
            'header in POTENTIAL.INP file!'
            CALL QUIT('Corrupted POTENTIAL.INP')
         END IF
!        Read force field data
         DO I=1,TMMFF
            IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
               READ(LUQMNP,*) MMFM0(I)
            END IF
            IF ((.NOT.DOMMCAP).AND.DOMMPOL) THEN
               READ(LUQMNP,*) MMFM0(I),MMFPOL(I)
            END IF
         END DO
!        Check force field definitions
         DO I=1,TMMATM
            IF (MMFTYP(I).GT.TMMFF) THEN
               WRITE(LUPRI,'(/2X,A,I4,1X,A)')                           &
          'Unknown force field requested for atom', I,                  &
          'in MM region!'
               CALL QUIT('Corrupted POTENTIAL.INP')
            END IF
!           Set polarization centers
            IF (MMFPOL(MMFTYP(I)) .GT. 1.0D-6) THEN
               TPOLATM = TPOLATM+1
               MMSKIP(I) = 1
!              works only for one non-metallic MM region
            END IF
         END DO
         IF (IPRTLVL.GE.5) THEN
            CALL PRINT_MMREGION(MMLBL)
         END IF
      END IF
!     Close POTENTIAL.INP file
      CALL GPCLOSE(LUQMNP,'KEEP')
!
!
   end subroutine
      SUBROUTINE PRINT_NPREGION(ATMLBL)
!
! Purpose:
!     prints detailed information about nanoparticles
!     in QM/NP/MM embedding
!
! Input:
!  ATMLBL - List of atom labels for NP region
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

#include "priunit.h"
#include "qmnpmm.h"
!
      CHARACTER ATMLBL(MXNPATM)*2
      integer :: i, j, istart, joff
!
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: Nanoparticle(s) data','*',103)
      ELSE
        CALL TITLER('QM/NP Embedding: Nanoparticle(s) data','*',103)
      END IF
!
      IF (TNPBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'Nanoparticle region contains', &
          TNPBLK, 'separate nanoparticles.'
      ELSE
        WRITE(LUPRI,'(/2X,A,A)') 'Nanoparticle region contains single', &
         ' nanoparticle.'
      END IF
      ISTART = 0
      DO I=1,TNPBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A,F6.3,A)') 'Nanoparticle',I,      &
          'with charge', NPCHRG(I), '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')                                         &
          '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')                                         &
          'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         DO J=1,NPATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')      &
             ATMLBL(JOFF), NPCORD(1,JOFF), NPCORD(2,JOFF),              &
             NPCORD(3,JOFF), NPFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A,A)') 'Static force field(s) data for ',       &
       'nanoparticle region.'
      WRITE(LUPRI, '(2X,A)')                                            &
       '-------------------------------------'
      WRITE(LUPRI, '(2X,A)')                                            &
       'FF Type Polarizability    Capacitance'
      WRITE(LUPRI, '(2X,A)')                                            &
       '====================================='
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,5X,F11.5,5X,F11.5)') I, NPFPOL(I),         &
          NPFCAP(I)
      END DO
      WRITE(LUPRI, '(2X,A)')                                            &
       '====================================='
      WRITE(LUPRI,'(/2X,A,A)') 'Dynamic force field(s) data for ',      &
       'nanoparticle region.'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '---------------------------------------------------',           &
       '------'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       'FF Type Omega(1)  Gamma(1)  Omega(2)  Gamma(2)',                &
       '  Factor'
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '===================================================',           &
       '======'
      DO I=1,TNPFF
         WRITE(LUPRI,'(3X,I2,3X,5(1X,F9.5))') I, NPFOMG1(I), NPFGAM1(I),&
          NPFOMG2(I), NPFGAM2(I), NPFFAC(I)
      END DO
      WRITE(LUPRI, '(2X,A,A)')                                          &
       '===================================================',           &
       '======'
!
   end subroutine
      SUBROUTINE PRINT_MMREGION(ATMLBL)
!
! Purpose:
!     prints detailed information about MM region
!     in QM/NP/MM embedding
!
! Input:
!  ATMLBL - List of atom labels for MM region
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "priunit.h"
#include "qmnpmm.h"
!
      CHARACTER ATMLBL(MXNPATM)*2
      integer :: i, istart, j, joff
!
      IF (DOMMSUB) THEN
         CALL TITLER('QM/NP/MM Embedding: MM region(s) data','*',103)
      END IF
!
      IF (TMMBLK.GT.1) THEN
         WRITE(LUPRI,'(/2X,A,I3,1X,A)') 'MM region contains',           &
          TMMBLK, 'separate non-metallic subregions.'
      END IF
      ISTART = 0
      DO I=1,TMMBLK
         WRITE(LUPRI, '(/,2X,A,I3,1X,A)') 'MM subregion',I,             &
          '. Coordinates in au.'
         WRITE(LUPRI, '(2X,A)')                                         &
          '--------------------------------------------------'
         WRITE(LUPRI, '(2X,A)')                                         &
          'Atom    Coord. X    Coord. Y    Coord. Z   FF Type'
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
         DO J=1,MMATOM(I)
            JOFF = ISTART+J
            WRITE(LUPRI,'(3X,A,1X,F11.5,1X,F11.5,1X,F11.5,3X,I4)')      &
             ATMLBL(JOFF), mm_cord(1,JOFF), mm_cord(2,JOFF),              &
             mm_cord(3,JOFF), MMFTYP(JOFF)
         END DO
         WRITE(LUPRI, '(2X,A)')                                         &
          '=================================================='
!        works only for one MM region
         WRITE(LUPRI, '(2X,A,I4)') 'Number of polarizable centers   : ',&
          TPOLATM
         WRITE(LUPRI, '(2X,A,I4)') 'Number of nonpolarizable centers: ',&
          MMATOM(I)-TPOLATM
         ISTART = ISTART + NPATOM(I)
      END DO
      WRITE(LUPRI,'(/2X,A)') 'Force field(s) data for MM region'
!
      IF ((.NOT.DOMMCAP).AND.(.NOT.DOMMPOL)) THEN
         WRITE(LUPRI, '(2X,A)') '--------------------------------'
         WRITE(LUPRI, '(2X,A)') 'FF Type MM Charge Polarizability'
         WRITE(LUPRI, '(2X,A)') '================================'
         DO I=1,TMMFF
            WRITE(LUPRI,'(3X,I2,3X,1X,F9.5,4X,F9.5,1X,I2)') I, MMFM0(I),&
             MMFPOL(I), MMSKIP(I)
         END DO
         WRITE(LUPRI, '(2X,A)') '================================'
      END IF
   end subroutine


   pure subroutine getdim_relmat(imatdim, sqflag)
!
! Purpose:
!     determines Relay matrix dimensions.
!
! Input:
!   SQFLAG - Request total size of Relay matrix
! Output:
!  IMATDIM - Size of Relay matrix or it's dimension
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      integer, intent(out) :: imatdim
      logical, intent(in)  :: sqflag

#include "qmnpmm.h"
!
      IMATDIM = 0
!     Add nanoparticle contribution
      IF (DONPSUB) THEN
         IMATDIM = IMATDIM + 3*TNPATM
         IF (DONPCAP) THEN
            IMATDIM = IMATDIM + TNPATM
         END IF
      END IF
!     Add Lagrangian term
      IF (DONPCAP.OR.DOMMCAP) THEN
         IMATDIM = IMATDIM + 1
      END IF
!     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
!
   end subroutine


   pure subroutine getdim_mmmat(imatdim, sqflag)
!
! Purpose:
!     determines Relay matrix dimensions for MM region.
!
! Input:
!   SQFLAG - Request total size of Relay matrix
! Output:
!  IMATDIM - Size of Relay matrix or it's dimension
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      integer, intent(out) :: imatdim
      logical, intent(in)  :: sqflag

#include "qmnpmm.h"

      IMATDIM = 0
!     Add molecular contribution
      IF (DOMMSUB) THEN
         IMATDIM = IMATDIM + 3*TPOLATM
      END IF
!     Get requested dimmension
      IF ((.NOT.MQITER).AND.SQFLAG) THEN
         IMATDIM = IMATDIM*IMATDIM
      END IF
!
   end subroutine


      SUBROUTINE COMP_RELMAT(FMAT,WORK,LWORK)
!
! Purpose:
!     Computes Relay matrix and inverts it or estimates initial
!     MQ vector values for iterative MQ vector determination
!     algorithm.
!
! Input:
!   WORK  - Temporary memory array
!   LWORK - Size of temporary memory array
! Output:
! if .not. MQITER
!  FMAT   - Inverted Relay matrix
! else
!  FMAT   - Initial guess of real part of MQ vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!


#include "priunit.h"
#include "qmnpmm.h"

      real(8) :: fmat(:, :, :)
      integer :: lwork
      real(8) :: work(lwork)
      integer, allocatable :: ipiv(:)

      integer :: idimension, ierror
!
!     Initialize arrays
      CALL GETDIM_RELMAT(idimension,.TRUE.)
      fmat = 0.0d0
!     Reset matrix dimension parameter
      CALL GETDIM_RELMAT(idimension,.FALSE.)
      IF (.NOT.MQITER) THEN

!        compute polarizabilty dependent terms
         if (donppol .or. dommpol) then
            call get_amat(fmat)
         end if

!        compute capacitancy dependent term
         if (donpcap .or. dommcap) then
           !todo adapt these routines to handle complex case
           call get_cmat(fmat)
           call get_mmat(fmat)
           call get_qlag(fmat)
         end if

      ELSE
!       Conjugated gradient method via paradiso solver
      END IF
!     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay matrix ***'
            CALL OUTPUT(FMAT,1,idimension,1,idimension,idimension,idimension,1,LUPRI)
      END IF
!     Invert Relay matrix
      !todo adapt to complex case and invert that one
      IF (.NOT.MQITER) THEN
         allocate(ipiv(idimension))
         CALL DGETRF(idimension,idimension,FMAT,idimension,IPIV,IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_RELMAT')
         END IF
         CALL DGETRI(idimension,FMAT,idimension,IPIV,WORK,lwork,             &
                IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_RELMAT')
         END IF
         deallocate(ipiv)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay matrix ***'
            CALL OUTPUT(FMAT,1,idimension,1,idimension,idimension,idimension,1,LUPRI)
         END IF
      END IF
!
   end subroutine


   SUBROUTINE COMP_MMRELMAT(FMAT,WORK,LWORK)
!
! Purpose:
!     Computes Relay matrix and inverts it or estimates initial
!     MQ vector values for iterative MQ vector determination
!     algorithm. (MM region)
!
! Input:
!   WORK  - Temporary memory array
!   LWORK - Size of temporary memory array
! Output:
! if .not. MQITER
!  FMAT   - Inverted Relay matrix
! else
!  FMAT   - Initial guess of real part of MQ vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

#include "priunit.h"
#include "qmnpmm.h"
!
      real(8) :: FMAT(*), WORK(LWORK)
      integer, allocatable :: ipiv(:)
      integer :: idimension, ierror, lwork
!
!     Initialize arrays
      CALL GETDIM_MMMAT(idimension,.TRUE.)
      CALL DZERO(FMAT,idimension)
!     Reset matrix dimension parameter
      CALL GETDIM_MMMAT(idimension,.FALSE.)
      IF (.NOT.MQITER) THEN
!        Compute polarizabilty dependent terms
         CALL GET_MM_AMAT(FMAT,idimension)
      ELSE
!       Conjugated gradient method via paradiso solver
      END IF
!     Print Relay matrix
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed Relay (MM) matrix ***'
            CALL OUTPUT(FMAT,1,idimension,1,idimension,idimension,idimension,1,LUPRI)
      END IF
!     Invert Relay matrix
      IF (.NOT.MQITER) THEN
         allocate(ipiv(idimension))
         CALL DGETRF(idimension,idimension,FMAT,idimension,IPIV,IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('LU factorization failed in COMP_MMRELMAT')
         END IF
         CALL DGETRI(idimension,FMAT,idimension,IPIV,WORK,lwork,             &
                IERROR)
         IF (IERROR.NE.0) THEN
            CALL QUIT('Inversion failed in COMP_MMRELMAT')
         END IF
         deallocate(ipiv)
         IF (IPRTLVL.GE.15) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Inverted Relay(MM) matrix ***'
            CALL OUTPUT(FMAT,1,idimension,1,idimension,idimension,idimension,1,LUPRI)
         END IF
      END IF
!
   end subroutine


   subroutine get_amat(fmat)
      ! computes a matrix component of relay matrix for np and mm regions

      real(8), intent(inout) :: fmat(:, :, :) ! relay matrix with m matrix contribution added up

#include "qmnpmm.h"

      integer :: i, j, k, l, m
      integer :: joff, koff, loff
      real(8) :: rij(3), rpol, rad, rad2, rad51, rval
      real(8) :: facta, factb

!     Set diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            RPOL = 1.0d0/NPFPOL(NPFTYP(I))
            DO J=1,3
               JOFF = (I-1)*3+J
               FMAT(JOFF,JOFF,1) = RPOL
            END DO
         END DO
      END IF
!     Set off-diagonal components
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            DO J=I+1,TNPATM
!              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD = dsqrt(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD2 = RAD*RAD
               RAD51 = 1.0d0/(RAD2*RAD2*RAD)
               IF (NPMQGAU) CALL GET_GG_AFACT(FACTA,FACTB,I,J,RAD)
!              Distribute interaction tensor
               DO K=1,3
                  KOFF = (I-1)*3+K
                  DO L=1,3
                     LOFF = (J-1)*3+L
                     RVAL = 3.0d0*RIJ(K)*RIJ(L)
                     IF (K.EQ.L) RVAL = RVAL-RAD2
                     RVAL = RVAL*RAD51
                     IF (NPMQGAU) THEN
                        RVAL = RVAL*FACTA-FACTB*RIJ(K)*RIJ(L)
                     END IF
                     FMAT(KOFF,LOFF,1) = -RVAL
                     FMAT(LOFF,KOFF,1) = -RVAL
                  END DO
               END DO
            END DO
         END DO
      END IF
   end subroutine


   pure subroutine get_mm_amat(fmat, idimension)
      ! computes a matrix component of relay matrix for mm region


      real(8), intent(inout) :: fmat(idimension, idimension) ! relay matrix with m matrix contribution added up
      integer, intent(in)    :: idimension

#include "qmnpmm.h"

      integer :: istart, i, j, k, l, m
      integer :: ioff, joff, koff, loff
      real(8) :: rij(3), rad, rad3, rval, fact
      real(8) :: rpol, rad2, rad51

!     Set diagonal components
       IOFF = 1
       DO I=1,TMMATM
         IF (MMSKIP(I) .EQ. 0) CYCLE
         RPOL = 1.0d0/MMFPOL(MMFTYP(I))
         DO J=1,3
            JOFF = (IOFF-1)*3+J
            FMAT(JOFF,JOFF) = RPOL
         END DO
         IOFF = IOFF+1
      END DO
!     Set off-diagonal components
      IOFF = 1
      DO I=1,TMMATM
         IF (MMSKIP(I) .EQ. 0) CYCLE
         JOFF = IOFF+1
         DO J=I+1,TMMATM
            IF (MMSKIP(J) .EQ. 0) CYCLE
!           Compute distance dependent parameters
            RIJ(1) = mm_cord(1,I)-mm_cord(1,J)
            RIJ(2) = mm_cord(2,I)-mm_cord(2,J)
            RIJ(3) = mm_cord(3,I)-mm_cord(3,J)
            RAD = dsqrt(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
            RAD2 = RAD*RAD
            RAD51 = 1.0d0/(RAD2*RAD2*RAD)
!           Distribute interaction tensor
            DO K=1,3
               KOFF = (IOFF-1)*3+K
               DO L=1,3
                  LOFF = (JOFF-1)*3+L
                  RVAL = 3.0d0*RIJ(K)*RIJ(L)
                  IF (K.EQ.L) RVAL = RVAL-RAD2
                  RVAL = RVAL*RAD51
                  FMAT(KOFF,LOFF) = -RVAL
                  FMAT(LOFF,KOFF) = -RVAL
               END DO
            END DO
            JOFF = JOFF+1
         END DO
         IOFF = IOFF+1
      END DO
!
   end subroutine


   subroutine get_gg_afact(facta, factb, iatm, jatm, rad)
!
! Purpose:
!     Determines damping factors for T(2) operator for
!     Gaussian/Gaussian dipole model.
!
! Input:
!   IATM - I-th atom in NP region
!   JATM - J-th atom in NP region
!   RAD  - Distance between I and J atoms
! Output:
!   FACTA  - Scalling factor
!   FACTB  - Additional factor
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      real(8), intent(out) :: facta
      real(8), intent(out) :: factb
      integer, intent(in)  :: iatm
      integer, intent(in)  :: jatm
      real(8), intent(in)  :: rad

#include "pi.h"
#include "qmnpmm.h"

      real(8) :: ripol
      real(8) :: rjpol
      real(8) :: rdim
      real(8) :: rim
      real(8) :: rjm
      real(8) :: rijm
      real(8) :: rval
      real(8),external :: erf
!     Get polarizabilities
      RIPOL = NPFPOL(NPFTYP(IATM))/3.0d0
      RJPOL = NPFPOL(NPFTYP(JATM))/3.0d0
!     Get damping radius
      RDIM = dsqrt(2.0d0)/SQRTPI
      RIM  = (RDIM*RIPOL)**(1.0d0/3.0d0)
      RJM  = (RDIM*RJPOL)**(1.0d0/3.0d0)
      RIJM = dsqrt(RIM*RIM+RJM*RJM)
      RVAL = RAD/RIJM
!     Compute factors
      FACTA = ERF(RVAL)-2.0d0*RVAL*DEXP(-RVAL*RVAL)/SQRTPI
      FACTB = 4.0d0*DEXP(-RVAL*RVAL)
      FACTB = FACTB/(RAD*RAD*RIJM*RIJM*RIJM*SQRTPI)
!
   end subroutine


   subroutine get_cmat(fmat)
      ! computes c matrix component of relay matrix for np and mm regions

      real(8), intent(inout) :: fmat(:, :, :) ! relay matrix with m matrix contribution added up

#include "qmnpmm.h"

      integer :: istart, i, j, m
      integer :: ioff, joff
      real(8) :: rij(3), rad, rad3, rval, fact

!     Set diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
!        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            FMAT(IOFF,IOFF,1) = -1.0d0/NPFCAP(NPFTYP(I))
         END DO
      END IF
!     Set off-diagonal components
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
!        Fix me: MM shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=I+1,TNPATM
               JOFF = ISTART+J
!              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD  = dsqrt(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RVAL = 1.0d0/RAD
               IF (NPMQGAU) THEN
                  CALL GET_GG_CFACT(FACT,I,J,RAD)
                  RVAL = FACT*RVAL
               END IF
               FMAT(IOFF,JOFF,1) = -RVAL
               FMAT(JOFF,IOFF,1) = -RVAL
            END DO
         END DO
      END IF
   end subroutine


   subroutine get_gg_cfact(fact, iatm, jatm, rad)
!
! Purpose:
!     Determines damping factors for T(0) operator for
!     Gaussian/Gaussian dipole model.
!
! Input:
!   IATM - I-th atom in NP region
!   JATM - J-th atom in NP region
!   RAD  - Distance between I and J atoms
! Output:
!   FACT  - Scalling factor
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

      real(8), intent(out) :: fact
      integer, intent(in)  :: iatm
      integer, intent(in)  :: jatm
      real(8), intent(in)  :: rad

#include "pi.h"
#include "qmnpmm.h"

      real(8) :: ricap
      real(8) :: rjcap
      real(8) :: rijm
      real(8),external :: erf
      
!     Get capacitancies
      RICAP = 1.41421356237309504880D0*NPFCAP(NPFTYP(IATM))/SQRTPI
      RJCAP = 1.41421356237309504880D0*NPFCAP(NPFTYP(JATM))/SQRTPI

!     Get damping radius & scalling factor
      RIJM = dsqrt(RICAP*RICAP+RJCAP*RJCAP)
      FACT = ERF(RAD/RIJM)

   end subroutine


   subroutine get_mmat(fmat)
      ! computes m matrix component of relay matrix for np and mm regions

      real(8), intent(inout) :: fmat(:, :, :) ! relay matrix with m matrix contribution added up

#include "qmnpmm.h"

      integer :: istart, i, j, m
      integer :: ioff, joff
      real(8) :: rij(3), rad, rad3, rval, fact

!     Set off-diagonal components
      IF (DONPPOL.AND.DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = ISTART+3*TNPATM
         DO I=1,TNPATM
            DO J=1,TNPATM
               IF (I.EQ.J) CYCLE
!              Compute distance dependent parameters
               RIJ(1) = NPCORD(1,I)-NPCORD(1,J)
               RIJ(2) = NPCORD(2,I)-NPCORD(2,J)
               RIJ(3) = NPCORD(3,I)-NPCORD(3,J)
               RAD  = dsqrt(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
               RAD3 = 1.0d0/(RAD*RAD*RAD)
!              Get damping factor
               IF (NPMQGAU) CALL GET_GG_MFACT(FACT,I,J,RAD)
!              Distribute contributions
               DO M=1,3
                  IOFF = (I-1)*3+M
                  JOFF = ISTART+J
                  RVAL = RIJ(M)*RAD3
                  IF (NPMQGAU) RVAL = RVAL*FACT
                  FMAT(IOFF,JOFF,1) = -RVAL
                  FMAT(JOFF,IOFF,1) = -RVAL
               END DO
            END DO
         END DO
      END IF

   end subroutine


   subroutine get_gg_mfact(fact, iatm, jatm, distance_i_j)
      ! determines damping factors for t(1) operator for
      ! gaussian/gaussian dipole model


      real(8), intent(inout) :: fact ! scaling factor
      integer, intent(in)    :: iatm ! i-th atom in np region
      integer, intent(in)    :: jatm ! j-th atom in np region
      real(8), intent(in)    :: distance_i_j

! sqrtpi
#include "pi.h"

! npftyp, npfpol, npfcap
#include "qmnpmm.h"

      real(8) :: ripol
      real(8) :: rdim
      real(8) :: rim
      real(8) :: rjcap
      real(8) :: rijm
      real(8) :: radx
      real(8),external :: erf

!     get polarizabilities
      ripol = npfpol(npftyp(iatm))/3.0d0
!     get damping radius
      rdim = 1.41421356237309504880d0/sqrtpi
      rim  = (rdim*ripol)**(1.0d0/3.0d0)
!     get j-th capacitancy
      rjcap = rdim*npfcap(npftyp(jatm))
!     get damping radius & scalling factor
      rijm = dsqrt(rim*rim+rjcap*rjcap)
      radx = distance_i_j/rijm
      fact = erf(radx)-2.0d0*radx*dexp(-radx*radx)/sqrtpi

   end subroutine


   pure subroutine get_qlag(fmat)
      ! sets charge constrain in relay matrix for np and mm regions

      real(8), intent(inout) :: fmat(:, :, :)

! donpcap, dommcap, donppol
#include "qmnpmm.h"

      integer :: i, istart, ioff, idimension
      idimension = size(fmat, 1)

      ! set diagonal components
      if (donpcap .or. dommcap) then
         istart = 0
         if (donppol) istart = istart + 3*tnpatm
!        fixme: mm shift
         if (donpcap) then
            do i = 1, tnpatm
               ioff = istart + i
               fmat(ioff, idimension, 1) = 1.0d0
               fmat(idimension, ioff, 1) = 1.0d0
            end do
         end if
      end if

   end subroutine



      SUBROUTINE WRITE_RELMAT(FMAT)
!
! Purpose:
!     Stores Relay matrix in binary file.
!
! Input:
!  FMAT   - Inverted Relay matrix
!
! Last updated: 16/08/2013 by Z. Rinkevicius.
!

#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
!
      real(8) :: fmat(*)
      integer :: idimension
      integer :: luqmnp

!     determine dimensions
      CALL GETDIM_RELMAT(idimension,.TRUE.)
!     write inverted Relay matrix
      LUQMNP = -1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED', &
             IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      CALL WRTIEF(FMAT, idimension, 'QQMNPMAT', LUQMNP)
      CALL GPCLOSE(LUQMNP,'KEEP')
!
      end subroutine
      SUBROUTINE READ_RELMAT(FMAT)
!
! Purpose:
!     Reads Relay matrix in binary file.
!
! Input:
!  FMAT   - Inverted Relay matrix
!
! Last updated: 16/08/2013 by Z. Rinkevicius.
!

#include "qmnpmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "inftap.h"
!
      real(8) :: fmat(*)
      integer :: idimension
      integer :: luqmnp
      logical :: fndlab

!     determine dimensions
      CALL GETDIM_RELMAT(idimension,.TRUE.)

!     read inverted Relay matrix
      LUQMNP=-1
      CALL GPOPEN(LUQMNP,'QMMMNP','UNKNOWN','SEQUENTIAL','UNFORMATTED', &
             IDUMMY,.FALSE.)
      REWIND(LUQMNP)
      IF (FNDLAB('QQMNPMAT',LUQMNP)) THEN
        CALL READT(LUQMNP,idimension,FMAT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMNP file.')
      ENDIF
      CALL GPCLOSE(LUQMNP,'KEEP')
!
      end subroutine




   subroutine comp_dampvmat(fmat, mqvec)
      ! computes damped potential matrix in ao basis


      real(8), intent(in)    :: fmat(*) ! inverted relay matrix
      real(8), intent(inout) :: mqvec(*)

#include "priunit.h"
#include "qmnpmm.h"
#include "maxorb.h"
#include "aovec.h"
#include "shells.h"
#include "primit.h"

#ifdef VAR_MPI
! qmcmm_work
#include "iprtyp.h"
#endif

      real(8), allocatable :: rdvec(:)
      real(8), allocatable :: rqvec(:)
      integer              :: idimension, i
      integer              :: iprint

      call getdim_relmat(idimension, .false.)

      ! allocate and set gaussian broadening paramenters

      allocate(rdvec(tnpatm))
      allocate(rqvec(tnpatm))

      call set_damparam(rdvec, rqvec)

      if (iprtlvl > 14) then
         write(lupri, '(/,2x,a)') '*** Computed MQ vector start ***'
         do i = 1, idimension
            write(lupri, '(i8, f18.8)') i, mqvec(i)
         end do
         write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
      end if

#ifdef VAR_MPI
      call mpixbcast(QMCMM_WORK, 1, 'INTEGER', 0)
      iprint = 0
      call mpixbcast(iprint, 1, 'INTEGER', 0)
#endif

#ifdef ENABLE_VPOTDAMP
      call vpotdamped(kmax, nhkt, nuco, nrco, jstrt, cent, priccf, priexp, &
                      fmat,                                                &
                      npcord,                                              &
                      mqvec(3*tnpatm+1), rqvec,                            &
                      mqvec            , rdvec,                            &
                      tnpatm)
#else
      call quit('VPOTDAMP not compiled in this version')
#endif

      deallocate(rdvec)
      deallocate(rqvec)

   end subroutine


   pure subroutine set_damparam(rdvec, rqvec)
      ! sets damping parameters vectors for dipoles and charges
      ! see Eqs 14 and 15 in J. Phys. Chem. C 112, 40 (2008)


      real(8), intent(inout) :: rdvec(*)
      real(8), intent(inout) :: rqvec(*)

! tnpatm, npfpol, npfcap
#include "qmnpmm.h"

! sqrtpi
#include "pi.h"

      real(8) :: rdim
      real(8) :: ripol
      integer :: i

      rdim = dsqrt(2.0d0)/sqrtpi
      do i = 1, tnpatm
         ripol = npfpol(npftyp(i))/3.0d0
         rdvec(i) = (rdim*ripol)**(1.0d0/3.0d0)
         rqvec(i) = rdim*npfcap(npftyp(i))
      end do

   end subroutine

end module
