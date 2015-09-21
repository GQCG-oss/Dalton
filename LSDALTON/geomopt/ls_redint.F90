!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program, Release 2.0
!...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
!...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
!...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
!...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
!...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
!...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
!...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
!...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
!...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
!...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
!...   For questions concerning this copyright write to:
!...      dalton-admin@kjemi.uio.no
!...
!...   For information on how to get a licence see:
!...      http://www.kjemi.uio.no/software/dalton/dalton.html
!
!
! File:  abacus/abarint.F
!
! New module taking care of redundant internal coordinates.
! The procedures are called from abaopt.F /960818-vebjornb
!
!  /* Deck inired */
      SUBROUTINE LS_INIRED(MXRCRD,MX2CRD,WILBMT,BMTRAN,BMTINV,PJINMT, &
     &                  TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,TMPMT6, &
     &                  Molecule,NAtoms,lupri,optinfo)
use ls_util
use optimization_input 
use files
use molecule_type
use molecule_typetype
use precision
!
!     Initializes everything that has to do with redundant internal
!     coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri,NAtoms
      Type(opt_setting) :: optinfo
      Type(MOLECULEINFO) :: Molecule
      Real(realk) WILBMT(MXRCRD*MX2CRD), BMTRAN(MXRCRD*MXRCRD)
      Real(realk) BMTINV(MXRCRD*MX2CRD), PJINMT(MXRCRD*MXRCRD)
      Real(realk) TMPMAT(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      Real(realk) TMPMT3(MX2CRD,MX2CRD), TMPMT4(MX2CRD,MX2CRD)
      Real(realk) TMPMT5(MX2CRD,MX2CRD), TMPMT6(MX2CRD,MX2CRD)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Integer IBOND(MXCENT,MXCENT)
!
!     First we have to determine the redundant internal coordinates.
!
      call ls_FNDRED(TMPMAT,IBOND,NAtoms,Molecule,lupri,optinfo,.FALSE.)
      IF (optinfo%RedRed) THEN
         call ls_RREDUN(lupri,optinfo)
         call lsquit('ICRD: no value assigned to this variable',-1)
!         IF (optinfo%IPrint .GE. IPRMIN) WRITE(LUPRI,'(/A,I6)') &
!     &        ' Number of internal coordinates after removals :',ICRD
      END IF
!
!     Then we construct Wilsons B matrix
!
      call ls_GETWIL(Molecule,MXRCRD,MX2CRD,TMPMAT,WILBMT, &
      & BMTRAN,TMPMT2,lupri,NAtoms,optinfo)
!
!     ... and the derivative of the B matrix
!
      IF (optinfo%IPrint .GE. IPRDBG) &
     &     call ls_GETDWL(Molecule,MXRCRD,TMPMAT,TMPMT2, &
     &     TMPMT3,WILBMT,lupri,NAtoms,optinfo)
!
!     We also need the inverse of the rectangular B matrix.
!
      call ls_GTBINV(MXRCRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,WILBMT, &
     &     BMTRAN,BMTINV,PJINMT,TMPMT5,TMPMT6,lupri,optinfo)
      RETURN
      END

!  /* Deck fndred */
      SUBROUTINE LS_FNDRED(ATMARR,IBOND,IATOM,Molecule,lupri,optinfo,FIND)
use ls_util   
use precision
use optimization_input 
use files
use molecule_type
use molecule_typetype
use memory_handling
!
!     Finds natural redundant internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Type(MOLECULEINFO) :: Molecule
      Integer :: lupri,IATOM
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8)
      Integer IBOND(MXCENT,MXCENT)
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      Logical :: FIND
!
!
      IPRSAV = optinfo%IPrint
      IF (optinfo%FindRe) optinfo%IPrint = MAX(optinfo%IPrint,IPRMIN)
      call ls_DZERO(ATMARR,8*MXCENT)
      optinfo%INTCRD = 0
!      call ls_IZERO(optinfo%INTCRD,48*MXCENT)
!      call ls_IZERO(IBOND,MXCENT*MXCENT)
      IBOND = 0
      optinfo%NIntCoord = 0
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Output from FNDRED')
      END IF
!
!     First we initialize the ATMARR array.
!
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.FALSE.,lupri)
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Expanded array of atoms (Aangstrom)')
         call output(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      optinfo%ICartCoord = 3*IATOM
!
!     We find all bonds.
!
      call ls_FNDBND(ATMARR,IBOND,IATOM,lupri,optinfo)
!
!     Then all angles
!
      call ls_FNDANG(ATMARR,IBOND,IATOM,lupri,optinfo)
!
!     Finally all dihedral angles
!
      IF (.NOT. optinfo%NoDihe) call ls_FNDIHD(ATMARR,IBOND,IATOM,lupri,optinfo)
!
!     Determine value of the various coordinates
!
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.TRUE.,lupri)
      ITMP = optinfo%NIntCoord
!
!!!!! Vladimir: if we just want to find redundant internals
!!!!! we allocate optinfo%CoordInt here
!!!!! Probably, temporary here
      If (FIND) Call mem_alloc(optinfo%CoordInt,optinfo%NIntCoord)

      call ls_GETINT(IATOM,optinfo%NIntCoord,ATMARR,optinfo%CoordInt,lupri,optinfo)
      optinfo%NIntCoord = ITMP
!
!     Output
!
      IF (optinfo%IPrint .GE. IPRMIN) call ls_PRNRED(optinfo%NIntCoord,optinfo%CoordInt, &
      & lupri,optinfo)
!
!     Check that the number of internal coordinates do
!     not exceed the limit of 8*MXCENT.
!     If it does, try giving intelligent suggestions.
!
      IF ( ((optinfo%NIntCoord .GT. 3*8*MXCENT) .AND. FIND) .OR. &
     & ((optinfo%RatFun .OR. optinfo%GDIIS) .AND. (optinfo%NIntCoord+1 .GT. 3*8*MXCENT))) THEN
         WRITE(LUPRI,'(//A//)') '***** CRITICAL ERROR IN OPTIMIZE  &
     &        MODULE *****'
         WRITE(LUPRI,'(A,I6,A)') 'The specified geometry requires', &
     &        optinfo%NIntCoord,' internal coordinates.'
         WRITE(LUPRI,'(A,I6,A)') 'This exceeds the current limit of', &
     &        8*MXCENT,' coordinates.'
         IF (optinfo%NIntCoord .GT. 15*optinfo%ICartCoord) THEN
            WRITE(LUPRI,'(/A)') 'WARNING: The molecular geometry has &
     &            probably been specified in units of Aangstrom'
            WRITE(LUPRI,'(A)') 'without setting the Aangstrom mark &
     &           in the input file (see manual).'
         ELSE IF (optinfo%CartCoord .AND. optinfo%InmdHess) THEN
            WRITE(LUPRI,'(/A)') 'Run the calculation with an &
     &           initial diagonal Hessian, for example:'
            WRITE(LUPRI,'(A)') '          .INITEV'
            WRITE(LUPRI,'(A/)') '          0.6E0_realk'
         ELSE
            WRITE(LUPRI,'(/A)') 'Try running the calculation in  &
     &           Cartesian coordinates.'
         END IF
         call lsquit('Too many internal ccordinates!',lupri)
      END IF
      IC = 1
 20   CONTINUE
      IF (IC .LE. optinfo%NIntCoord) THEN
         IF (optinfo%IConstr(IC) .EQ. 2) THEN
            DO 30 I = IC+1, optinfo%NIntCoord
               DO 32 J = 1, 6
                  optinfo%INTCRD(I-1,J) = optinfo%INTCRD(I,J)
 32            CONTINUE
               optinfo%IConstr(I-1) = optinfo%IConstr(I)
 30         CONTINUE
            optinfo%NIntCoord = optinfo%NIntCoord - 1
         ELSE
            IC = IC + 1
         END IF
         GOTO 20
      END IF
      optinfo%IPrint = IPRSAV
!!!!! Vladimir: we then deallocate optinfo%CoordInt if we just find 
!!!!! the number of internals
      IF (FIND) call mem_dealloc(optinfo%CoordInt)
!
      RETURN
      END
!==================!
! Write_topology   !
!==================!
Subroutine Write_topology(optinfo)
use ls_util
use precision
use optimization_input
use files
Implicit none
Type(opt_setting) :: optinfo
Integer :: i, ii_file
! First we open a formatted file
ii_file = -1 
Write(*,*) 'TOPOLOGY'
call lsopen(ii_file,'internal_interface.dat','NEW','FORMATTED')
! Write information
Write(ii_file,'(A21)')'Topology of internals'
Do i = 1,optinfo%NIntCoord
   ! First all kinds of bonds
   If (optinfo%INTCRD(i,1) .LT. 10) then
      WRITE(ii_file,'(I8,I5,I5,I5)') i, &
      1,optinfo%INTCRD(i,2),optinfo%INTCRD(i,3)
   Endif
   ! Second, angles
   If ((optinfo%INTCRD(i,1) .GT. 10) .AND. (optinfo%INTCRD(i,1) .LT. 20)) then
      WRITE(ii_file,'(I8,I5,I5,I5,I5)') i, &
      2,optinfo%INTCRD(i,2),optinfo%INTCRD(i,3),optinfo%INTCRD(i,4)
   Endif
   ! Third, dihedrals
   If (optinfo%INTCRD(i,1) .GT. 20) then
      WRITE(ii_file,'(I8,I5,I5,I5,I5,I5)') i, &
      3,optinfo%INTCRD(i,2),optinfo%INTCRD(i,3),optinfo%INTCRD(i,4),optinfo%INTCRD(i,5)
   Endif
Enddo
! Write Cartesian coordinates
Write(ii_file,'(A21)')'Cartesian coordinates'
Do i = 1, optinfo%ICartCoord/3
   Write(ii_file,'(F15.5,F15.5,F15.5)') optinfo%Coordinates(1,i), &
   optinfo%Coordinates(2,i), optinfo%Coordinates(3,i)
Enddo
! Write the step in internals
Write(ii_file,'(A17)')'Step in internals'
Do i = 1, optinfo%NIntCoord
   Write(ii_file,'(F15.5)') optinfo%StpInt(i)
Enddo
! Close the file
call lsclose(ii_file,'KEEP')
!
End subroutine write_topology

!  /* Deck prnred */
      SUBROUTINE LS_PRNRED(ICRD,VALINT,lupri,optinfo)
use ls_util 
use precision
use optimization_input 
use files
use Fundamental
!
!     Prints out information about redundant internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) VALINT(ICRD)
      CHARACTER CRDTXT*26, MARK*3, UNTTXT*5
!
      call lsheader(lupri,'Redundant internal coordinates')
      IF (optinfo%IPrint .LT. 3) WRITE(LUPRI,'(A/)') &
     &   ' (Auxiliary bonds are not printed at this print level.)'
      optinfo%RemCoord = .FALSE.
      DO 10 I = 1, ICRD
         MARK = '   '
         IF (optinfo%ConOpt .AND. (optinfo%IConstr(I) .EQ. 1)) MARK = '*  '
         IF (optinfo%IConstr(I) .EQ. 2) THEN
            MARK = '#  '
            optinfo%RemCoord = .TRUE.
         END IF
         CRDTXT = '                          '
         IF (optinfo%INTCRD(I,1) .LT. 10) THEN
            UNTTXT = ' Ang.'
            IF (optinfo%INTCRD(I,1) .EQ. 1) CRDTXT='Regular bond'
            IF (optinfo%INTCRD(I,1) .EQ. 3) CRDTXT='Interfragment bond'
            IF (optinfo%INTCRD(I,1) .EQ. 4) CRDTXT='Additional interfrag. bond'
            IF (optinfo%INTCRD(I,1) .EQ. 5) CRDTXT='Hydrogen bond'
            IF (optinfo%INTCRD(I,1) .EQ. 7) CRDTXT='Auxiliary bond'
            IF (optinfo%INTCRD(I,1) .EQ. 9) CRDTXT='Manual bond'
            IF ((optinfo%INTCRD(I,1) .NE. 7) .OR. (optinfo%IPrint .GE. 3)) &
     &           WRITE(LUPRI,'(I8,2A,2I4,F22.5,A5)') I, &
     &           MARK,CRDTXT,optinfo%INTCRD(I,2),optinfo%INTCRD(I,3), &
     &           bohr_to_angstrom*VALINT(I),UNTTXT
         ELSE IF (optinfo%INTCRD(I,1) .GT. 20) THEN
            UNTTXT = ' deg.'
            IF (optinfo%INTCRD(I,1) .EQ. 21) CRDTXT='Dihedral angle'
            IF (optinfo%INTCRD(I,1) .EQ. 23) CRDTXT='Extra dihedral angle'
            IF (optinfo%INTCRD(I,1) .EQ. 25) &
     &           CRDTXT='Interfrag. dihedral angle'
            IF (optinfo%INTCRD(I,1) .EQ. 29) CRDTXT='Manual dihedral angle'
            WRITE(LUPRI,'(I8,2A,4I4,F14.3,A5)') I, &
     &           MARK,CRDTXT,optinfo%INTCRD(I,2),optinfo%INTCRD(I,3), &
     &           optinfo%INTCRD(I,4),optinfo%INTCRD(I,5),1.8E2_realk*VALINT(I)/PI,UNTTXT
         ELSE
            UNTTXT = ' deg.'
            IF (optinfo%INTCRD(I,1) .EQ. 11) CRDTXT='Regular angle'
            IF (optinfo%INTCRD(I,1) .EQ. 12) CRDTXT='Additional angle'
            IF (optinfo%INTCRD(I,1) .EQ. 19) CRDTXT='Manual angle'
            WRITE(LUPRI,'(I8,2A,3I4,F18.3,A5)') I, &
     &           MARK,CRDTXT,optinfo%INTCRD(I,2),optinfo%INTCRD(I,3), &
     &           optinfo%INTCRD(I,4),180E0_realk*VALINT(I)/PI,UNTTXT
         END IF
 10   CONTINUE
      IF (optinfo%ConOpt .OR. optinfo%RemCoord) THEN
         WRITE(LUPRI,*)
         IF (optinfo%ConOpt) WRITE(LUPRI,*) &
     &        '*) Constrained (fixed) coordinate.'
         IF (optinfo%RemCoord) WRITE(LUPRI,*) &
     &        '#) Coordinate will be removed.'
      END IF
      WRITE(LUPRI,'(/A,I6)') &
     &     ' Total number of redundant internal coordinates:',ICRD
      RETURN
      END
!!!!!! Vladimir :: temporary disabled, since no exzct Hessian is available 
!!!!!!in LSDALTON 
!  /* Deck redvib */
!      SUBROUTINE LS_REDVIB(NCORD,NMODE,MXRCRD,MX2CRD,EVEC,IFRQCM,IMAGIN, &
!     &     ATMARR,CRDRIN,WILBMT,BMTRAN,TMPMAT,TMPMT2,lupri,optinfo)
!use ls_util 
!use optimization_input 
!use files
!use fundamental
!
!     Converts normal coordinates to redundant internal coordinates.
!
!Implicit Real(realk) (A-H,O-Z)
!#include "mxcent.h"
!#include "nuclei.h"
!      Integer :: lupri
!      Type(opt_setting) :: optinfo
!      Real(realk) EVEC(NCORD,NCORD), IFRQCM(NCORD)
!      Real(realk) IMAGIN(NCORD), ATMARR(MXCENT,8)
!      Real(realk) CRDRIN(MXRCRD), WILBMT(MXRCRD,MXCOOR)
!      Real(realk) BMTRAN(MXRCRD,MXRCRD), TMPMAT(MX2CRD,MX2CRD)
!      Real(realk) TMPMT2(MX2CRD,MX2CRD)
!      LOGICAL LOGDEL
!      CHARACTER*4 CRDTXT
!      CHARACTER*1 CHRIMG(0:1)
!      DATA CHRIMG /' ','i'/
!
!      SXFAMU = SQRT(XFAMU)
!      IF (.NOT. (optinfo%RedInt .OR. optinfo%DelInt)) RETURN
!      ITMP = optinfo%NIntCoord
!      LOGDEL = optinfo%DelInt
!      IC = optinfo%NIntCoord
!      IF (optinfo%DelInt) IC = optinfo%NIntCoord
!      optinfo%DelInt = .FALSE.
!      call ls_DZERO(TMPMT2,MX2CRD*MX2CRD)
!      call lsheader(lupri,'Analysis in redundant internal coordinates')
!      call ls_ATMINI(ATMARR,IATOM,.TRUE.,lupri)
!!
!!     We have to subtract the components of the recently
!!     determined step.
!
!!!!!!!!!!! Vladimir: LS_WLKCOR removed!     
!!      TMPMT2 = optinfo%STPSYM
!!      call ls_WLKCOR(TMPMT2,optinfo%ICartCoord,MXCOOR,-1,lupri,optinfo)
!      DO 10 ICENT = 1, NUCIND
!         DO 15 J = 1, 3
!            CORD(J,ICENT) = CORD(J,ICENT) - TMPMT2(3*(ICENT-1)+J,1)
! 15      CONTINUE
! 10   CONTINUE
!!
!!     We have to calculate the Wilson B matrix and the
!!     values of all redundant internal coordinates.
!!
!!      call ls_GETINT(IATOM,MXRCRD,ATMARR,CRDRIN,lupri,optinfo)
!      IP = optinfo%IPrint
!      optinfo%IPrint = -1
!      call ls_GETWIL(MXRCRD,MX2CRD,ATMARR,WILBMT,BMTRAN,TMPMAT,lupri,optinfo)
!      optinfo%IPrint = IP
!      call lsheader(lupri,'Expanded array of atoms (Bohr)')
!      call output(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
!      call ls_PRNRED(IC,CRDRIN,lupri,optinfo)
!
!     The molecular coordinates are restored.
!
!      DO 20 ICENT = 1, NUCIND
!         DO 25 J = 1, 3
!            CORD(J,ICENT) = CORD(J,ICENT) + TMPMT2(3*(ICENT-1)+J,1)
! 25      CONTINUE
! 20   CONTINUE
!!
!!     Finally the eigenvectors are transformed
!!
!      call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
!      call ls_DZERO(TMPMT2,MX2CRD*MX2CRD)
!      DO 40 IMODE = 1, NMODE
!         DO 45 II = 1, IC
!            DO 47 ICOOR = 1, 3*IATOM
!               TMPMAT(II,IMODE) = TMPMAT(II,IMODE) + &
!     &              SXFAMU*EVEC(ICOOR,IMODE)*WILBMT(II,ICOOR)
! 47         CONTINUE
! 45      CONTINUE
! 40   CONTINUE
!      call lsheader(lupri,'Normal coordinates (internal)')
!
!      ISTR = 1
!      NBATCH = (NMODE + 4)/5
!      DO 100 IBATCH = 1, NBATCH
!         IEND = MIN(ISTR + 4,NMODE)
!         NUMB = IEND - ISTR + 1
!         WRITE (LUPRI,'(/A12,5(I5,A2,I4,A1))') '            ', &
!     &      (II,'  ',IFRQCM(II),CHRIMG(IMAGIN(II)), II = ISTR,IEND)
!         LENH = 10 + NUMB*12
!         WRITE (LUPRI,'(2X,70A1)') ('-', II = 1,LENH)
!         WRITE (LUPRI,'()')
!         DO 110 ICOOR = 1, optinfo%NIntCoord
!            CRDTXT = 'bend'
!            IF (optinfo%INTCRD(ICOOR,1) .LT. 10) CRDTXT = 'stre'
!            IF (optinfo%INTCRD(ICOOR,1) .GE. 20) CRDTXT = 'tors'
!            IF (optinfo%INTCRD(ICOOR,1) .NE. 7) &
!     &           WRITE (LUPRI,'(I5,1X,A4,(T13,5F12.6))') ICOOR,CRDTXT, &
!     &           (TMPMAT(ICOOR,II),II=ISTR,IEND)
! 110     CONTINUE
!         ISTR = ISTR + 5
! 100  CONTINUE
!      optinfo%NIntCoord = ITMP
!      optinfo%DelInt = LOGDEL
!      RETURN
!      END

!  /* Deck fndbnd */
      SUBROUTINE LS_FNDBND(ATMARR,IBOND,IATOM,lupri,optinfo)
use precision
use ls_util 
use files
use optimization_input 
use Fundamental
!
!     Finds all bonds.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8), VEC1(3), VEC2(3) 
      Integer IBOND(MXCENT,MXCENT)
      LOGICAL POSSIB
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
!     We have to find all regular bonds.
!
      IBEFOR = optinfo%NIntCoord
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Regular bonds')
      END IF
      DO 10 I = 1, IATOM - 1
         DO 15 J = I + 1, IATOM
            RADI = ATMARR(I,5)
            RADJ = ATMARR(J,5)            
            VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
            VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
            VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
            DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
!
!     Regular bonds have dist. less than 1.3 times sum of cov. radii
!
            IF (DIST .LE. 1.3E0_realk*(RADI+RADJ)) THEN
               IBOND(I,J) = 1
               IBOND(J,I) = 1
               optinfo%NIntCoord = optinfo%NIntCoord + 1
               optinfo%INTCRD(optinfo%NIntCoord,1) = 1
               optinfo%INTCRD(optinfo%NIntCoord,2) = I
               optinfo%INTCRD(optinfo%NIntCoord,3) = J
               IF (optinfo%IPrint .GE. IPRMAX) THEN
                  WRITE(LUPRI,'(A,2I4)') &
     &                 ' Bond between atoms      : ',i,j
               END IF
            END IF
 15      CONTINUE
 10   CONTINUE
      IF ((optinfo%IPRINT .GE. IPRMAX) .AND. (IBEFOR .EQ. optinfo%NIntCoord)) THEN
         WRITE(LUPRI,'(A)') ' None were found.'
      END IF
!
!     Next we have to check if single atoms or fragments are too far
!     away to be bonded by normal bonds, in which case interfragment
!     bond(s) must be assigned.
!
!     The first step is to assign fragment numbers to all atoms, this
!     process is done iteratively.
!
      IBEFOR = optinfo%NIntCoord
      IFRAG = 1
      IAT = 1
 20   CONTINUE
      ATMARR(IAT,6) = IFRAG*1.0E0_realk
 25   CONTINUE
      ICHANG = 0
      DO 30 II = 1,IATOM
         IF (NINT(ATMARR(II,6)) .EQ. IFRAG) THEN
            DO 35 JJ = 1, IATOM
               IF ((IBOND(II,JJ) .GE. 1) .AND. &
     &              (NINT(ATMARR(JJ,6)) .EQ. 0)) THEN
                  ATMARR(JJ,6) = IFRAG*1.0E0_realk
                  ICHANG = 1
               END IF
 35         CONTINUE
         END IF
 30   CONTINUE
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Assigning atoms to fragments')
         call output(ATMARR,1,IATOM,6,6,MXCENT,8,1,LUPRI)
      END IF
      IF (ICHANG .GT. 0) GOTO 25
!
      II = 1
 32   CONTINUE
      IF (II .LE. IATOM) THEN
         IF (ATMARR(II,6) .GT. 0.5E0_realk) THEN
           II = II + 1
           GOTO 32
         ELSE
           IFRAG = IFRAG + 1
           IAT = II
           ATMARR(IAT,6) = IFRAG*1.0E0_realk
           GOTO 20
        END IF
      END IF
      IF (optinfo%IPRINT .GE. IPRDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,I4)')  &
     &        ' Number of separate fragments found: ', IFRAG
      END IF
!
!     Two fragments get connected through the shortest interfragment distance.
!     All distances shorter than the minimum of 2.0 Ang. and this distance
!     are also marked as interfragment bonds.
!
      IF (IFRAG .GT. 1) THEN
         DSTMIN = 1.0E10_realk
         DO 40 II = 1, IATOM
            DO 45 JJ = 1, IATOM
               IF (NINT(ATMARR(II,6)) .NE. NINT(ATMARR(JJ,6))) THEN
                  VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                  VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                  VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                  DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                  IF (DIST .LT. DSTMIN) THEN
                     DSTMIN = DIST
                     IATM1 = II
                     IFRG1 = NINT(ATMARR(II,6))
                     IATM2 = JJ
                     IFRG2 = NINT(ATMARR(JJ,6))
                  END IF
               END IF
 45         CONTINUE
 40      CONTINUE
         IBOND(IATM1,IATM2) = 3
         IBOND(IATM2,IATM1) = 3
         optinfo%NIntCoord = optinfo%NIntCoord + 1
         optinfo%INTCRD(optinfo%NIntCoord,1) = 3
         optinfo%INTCRD(optinfo%NIntCoord,2) = IATM1
         optinfo%INTCRD(optinfo%NIntCoord,3) = IATM2
         BNDLIM = MIN(1.3E0_realk*DSTMIN, 2.0E0_realk)
!
!     Additional interfragment bonds
!
         DO 50 II = 1, IATOM
            IF (NINT(ATMARR(II,6)) .EQ. IFRG1) THEN
               DO 55 JJ = 1, IATOM
                  IF ((IBOND(II,JJ) .EQ. 0) .AND. &
     &                 (NINT(ATMARR(JJ,6)) .EQ. IFRG2)) THEN
                     VEC1(1) = ATMARR(II,2)-ATMARR(JJ,2)
                     VEC1(2) = ATMARR(II,3)-ATMARR(JJ,3)
                     VEC1(3) = ATMARR(II,4)-ATMARR(JJ,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     IF (DIST .LE. BNDLIM) THEN
                        BNDTYP = 4
!
!     We make a distinction between interfragment bonds that are
!     approximately as short as the shortest one and additional
!     longer interfragment bonds.
!
                        IF (DIST .LE. DSTMIN*1.05E0_realk) BNDTYP = 3
                        IBOND(II,JJ) = BNDTYP
                        IBOND(JJ,II) = BNDTYP
                        optinfo%NIntCoord = optinfo%NIntCoord + 1
                        optinfo%INTCRD(optinfo%NIntCoord,1) = BNDTYP
                        IF (II .LT. JJ) THEN
                           optinfo%INTCRD(optinfo%NIntCoord,2) = II
                           optinfo%INTCRD(optinfo%NIntCoord,3) = JJ
                        ELSE
                           optinfo%INTCRD(optinfo%NIntCoord,2) = JJ
                           optinfo%INTCRD(optinfo%NIntCoord,3) = II
                        END IF
                     END IF
                  END IF
 55            CONTINUE
            END IF
 50      CONTINUE
!
!     The process of connecting fragments runs iteratively, until all
!     atoms are connected, i.e. the number of fragments is one.
!
         IF (IFRAG .GT. 1) THEN
            DO 57 I = 1, IATOM
               ATMARR(I,6) = D0
 57         CONTINUE
            IAT = 1
            IFRAG = 1
            GOTO 20
         END IF
      END IF
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPRINT .GE. IPRMAX)) THEN
         call lsheader(lupri,'Interfragment bonds')
         DO 60 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,2I4)') &
     &           ' Interfragment bond between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3)
 60      CONTINUE
      END IF
!
!     Then we have to find all hydrogen bonds in the system.
!     X --- H ... Y (X,Y = N, O, F, P, S, Cl)
!
      IBEFOR = optinfo%NIntCoord
      DO 70 I = 1, IATOM
         IF (NINT(ATMARR(I,1)) .EQ. 1) THEN
            POSSIB = .FALSE.
            DO 72 J = 1, IATOM
               ID = NINT(ATMARR(J,1))
               IF ((IBOND(I,J) .EQ. 1) .AND. ((ID .EQ. 7) .OR. &
     &              (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15) &
     &              .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                  POSSIB = .TRUE.
                  VEC2(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC2(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC2(3) = ATMARR(J,4)-ATMARR(I,4)
               END IF
 72         CONTINUE
            IF (POSSIB) THEN
               VDWR1 = vdwrad_ls(1,lupri)
               DO 75 J = 1, IATOM
                  ID = NINT(ATMARR(J,1))
                  IF ((IBOND(I,J) .EQ. 0) .AND. ((ID .EQ. 7) .OR. &
     &                 (ID .EQ. 8) .OR. (ID .EQ. 9) .OR. (ID .EQ. 15) &
     &                  .OR. (ID .EQ. 16) .OR. (ID .EQ. 17))) THEN
                     VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                     VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                     VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                     DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     ICHG = NINT(ATMARR(J,1))
!
!     We require the distance to be less than 0.9 times the sum of the
!     van der Waals radii and larger than the covalent radii. Finally
!     we say that the angle (X-H-Y) has to be larger than 90 degrees.
!
                     VDWRI = vdwrad_ls(ICHG,lupri)
                     IF ((DIST .LE. 0.9E0_realk*(VDWR1+VDWRI)) &
     &                .AND. (DIST .GE. (CovRad(1,lupri)+CovRad(ICHG,lupri)))) THEN
                        BNDANG = vecang_ls(VEC1,VEC2)
                        IF (ABS(BNDANG) .GE. 0.49E0_realk*PI) THEN
                           IBOND(I,J) = 5
                           IBOND(J,I) = 5
                           optinfo%NIntCoord = optinfo%NIntCoord + 1
                           optinfo%INTCRD(optinfo%NIntCoord,1) = 5
                           IF (I .LT. J) THEN
                              optinfo%INTCRD(optinfo%NIntCoord,2) = I
                              optinfo%INTCRD(optinfo%NIntCoord,3) = J
                           ELSE
                              optinfo%INTCRD(optinfo%NIntCoord,2) = J
                              optinfo%INTCRD(optinfo%NIntCoord,3) = I
                           END IF
                        END IF
                     END IF
                  END IF
 75            CONTINUE
            END IF
         END IF
 70   CONTINUE
!
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Hydrogen bonds')
         DO 77 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,2I4)') &
     &           ' Hydrogen bond between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3)
 77      CONTINUE
      END IF
!
!     Bond coordinates that have been manually specified are added
!     (if not already present)
!
      IBEFOR = optinfo%NIntCoord
      IF (optinfo%AddCoord .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Manually added bonds')
      END IF
    If (optinfo%AddCoord) then
      DO 700 I = 1, optinfo%NAdd
         IF ((optinfo%AddCoordArray(I,3) .EQ. 0) .AND. (optinfo%AddCoordArray(I,4) .EQ. 0)) THEN
            IF ((optinfo%AddCoordArray(I,1) .GT. IATOM) &
     &           .OR. (optinfo%AddCoordArray(I,2) .GT. IATOM)) THEN
               WRITE(LUPRI,'(A,I4)') &
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number &
     &              is larger than the number of atoms.'
               call lsquit(' Error detected in the specification of &
     &              additional bond coordinate(s)!',lupri)
            END IF
            IF (IBOND(optinfo%AddCoordArray(I,1),optinfo%AddCoordArray(I,2)) .EQ. 0) THEN
               optinfo%NIntCoord = optinfo%NIntCoord + 1
               optinfo%INTCRD(optinfo%NIntCoord,1) = 9
               optinfo%INTCRD(optinfo%NIntCoord,2) = optinfo%AddCoordArray(I,1)
               optinfo%INTCRD(optinfo%NIntCoord,3) = optinfo%AddCoordArray(I,2)
            END IF
         END IF
 700  CONTINUE
    Endif
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         DO 710 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,2I4)') &
     &           ' Manual bond between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3)
 710     CONTINUE
      END IF
!
!     We also determine "extra bonds" between atoms. These will
!     not generate angles and dihedral angles, but they seem to
!     make geometry optimization in redundant coordinates more
!     efficient. Only used for approx. Hessians.
!
      IBEFOR = optinfo%NIntCoord
!      IF (.NOT. optinfo%NoAux) THEN
      IF (.NOT. optinfo%NoAux .AND. (optinfo%ModHes .OR. optinfo%CMBMod .OR. optinfo%InmdHess)) THEN
         IF (optinfo%IPrint .GE. IPRMAX)  call lsheader(lupri,'Auxiliary bonds')
         DO 80 I = 1, IATOM - 1
            DO 85 J = I + 1, IATOM
               RADI = ATMARR(I,5)
               RADJ = ATMARR(J,5)            
               VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
               VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
               VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
               DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
!
!     Extra bonds have dist. less than 2.5 times sum of cov. radii
!
               IF ((DIST .LE. 2.5E0_realk*(RADI+RADJ)) .AND. &
     &              (IBOND(I,J) .EQ. 0)) THEN
                  optinfo%NIntCoord = optinfo%NIntCoord + 1
                  optinfo%INTCRD(optinfo%NIntCoord,1) = 7
                  optinfo%INTCRD(optinfo%NIntCoord,2) = I
                  optinfo%INTCRD(optinfo%NIntCoord,3) = J
                  IF (optinfo%IPrint .GE. IPRMAX) WRITE(LUPRI,'(A,2I4)') &
     &                 ' Auxiliary bond between atoms: ',i,j
               END IF
 85         CONTINUE
 80      CONTINUE
         IF ((optinfo%IPrint .GE. IPRMAX) .AND. (IBEFOR .EQ. optinfo%NIntCoord))  &
             WRITE(LUPRI,'(A)') ' None were added.'
      END IF

      RETURN
      END

!  /* Deck fndang */
      SUBROUTINE LS_FNDANG(ATMARR,IBOND,IATOM,lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use Fundamental
!
!     Finds all angles.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8), VEC1(3), VEC2(3)
      Integer  IBOND(MXCENT,MXCENT)
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
!     We have to find all angles between bonded atoms
!
      IBEFOR = optinfo%NIntCoord
      DO 10 I = 1, IATOM
         NUMBND = 0
         DO 15 J = 1, IATOM
            IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 9)) &
     &           NUMBND = NUMBND + 1
 15      CONTINUE
         IF (NUMBND .GT. 1) THEN
            DO 20 J = 1, IATOM - 1
               IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 4)  &
     &              .AND. (IBOND(I,J) .NE. 9)) THEN
                  VEC1(1) = ATMARR(J,2)-ATMARR(I,2)
                  VEC1(2) = ATMARR(J,3)-ATMARR(I,3)
                  VEC1(3) = ATMARR(J,4)-ATMARR(I,4)
                  DO 25 K = J + 1, IATOM
                     IF ((IBOND(I,K) .GT. 0) .AND. &
     &                    (IBOND(I,K) .NE. 4) .AND. &
     &                    (IBOND(I,K) .NE. 9)) THEN
                        optinfo%NIntCoord = optinfo%NIntCoord + 1
                        optinfo%INTCRD(optinfo%NIntCoord,1) = 11
                        optinfo%INTCRD(optinfo%NIntCoord,2) = J
                        optinfo%INTCRD(optinfo%NIntCoord,3) = I
                        optinfo%INTCRD(optinfo%NIntCoord,4) = K
!
!     We have to check for linear or near linear bonds (these get two
!     internal orthogonal bond coordinates). The "critical" bond angle
!     is 175 degrees.
!
                        VEC2(1) = ATMARR(K,2)-ATMARR(I,2)
                        VEC2(2) = ATMARR(K,3)-ATMARR(I,3)
                        VEC2(3) = ATMARR(K,4)-ATMARR(I,4)
                        BNDANG = vecang_ls(VEC1,VEC2)
                        IF ((ABS(BNDANG) .GE. 175E0_realk*PI/180E0_realk) .AND. &
     &                       (.NOT. optinfo%NoAdda)) THEN
                           optinfo%NIntCoord = optinfo%NIntCoord + 1
                           optinfo%INTCRD(optinfo%NIntCoord,1) = 12
                           optinfo%INTCRD(optinfo%NIntCoord,2) = J
                           optinfo%INTCRD(optinfo%NIntCoord,3) = I
                           optinfo%INTCRD(optinfo%NIntCoord,4) = K
                        END IF
                     END IF
 25               CONTINUE
               END IF
 20         CONTINUE
         END IF
 10   CONTINUE
!
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Angles')
         DO 30 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,3I4)') &
     &           ' Angles between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3),optinfo%INTCRD(I,4)
 30      CONTINUE
      END IF
!
!     Angle coordinates that have been manually specified are added
!     (if not already present)
!
      IBEFOR = optinfo%NIntCoord
      IF (optinfo%AddCoord .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Manually added angles')
      END IF
   If (optinfo%AddCoord) then
      DO 50 I = 1, optinfo%NAdd
         IF ((optinfo%AddCoordArray(I,3) .GT. 0) .AND. (optinfo%AddCoordArray(I,4) .EQ. 0)) THEN
            IF ((optinfo%AddCoordArray(I,1) .GT. IATOM) .OR. (optinfo%AddCoordArray(I,2) .GT. IATOM)  &
     &           .OR. (optinfo%AddCoordArray(I,3) .GT. IATOM))  THEN
               WRITE(LUPRI,'(A,I4)') &
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number ' &
     &              // 'is larger than the number of atoms.'
               call lsquit(' Error detected in the specification of  &
     &              additional angle coordinate(s)!',lupri)
            END IF
            IBND1 = IBOND(optinfo%AddCoordArray(I,1),optinfo%AddCoordArray(I,2))
            IBND2 = IBOND(optinfo%AddCoordArray(I,2),optinfo%AddCoordArray(I,3))
!     Check if the three atoms should already have produced an angle coordinate
            IF (.NOT.(((IBND1.GT. 0).AND.(IBND1.NE. 4).AND.(IBND1.NE. 9)) &
     &      .AND.((IBND2.GT. 0).AND.(IBND2.NE. 4).AND.(IBND2.NE. 9)))) THEN
               optinfo%NIntCoord = optinfo%NIntCoord + 1
               optinfo%INTCRD(optinfo%NIntCoord,1) = 19
               optinfo%INTCRD(optinfo%NIntCoord,2) = optinfo%AddCoordArray(I,1)
               optinfo%INTCRD(optinfo%NIntCoord,3) = optinfo%AddCoordArray(I,2)
               optinfo%INTCRD(optinfo%NIntCoord,4) = optinfo%AddCoordArray(I,3)
            END IF
         END IF
 50   CONTINUE
   Endif
!
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         DO 710 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,3I4)') &
     &           ' Manual angle between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3),optinfo%INTCRD(I,4)
 710     CONTINUE
      END IF
      RETURN
      END

!  /* Deck fndihd */
      SUBROUTINE LS_FNDIHD(ATMARR,IBOND,IATOM,lupri,optinfo)
use precision
use ls_util 
use files
use optimization_input 
use Fundamental
!
!     Finds all dihedral angles.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8)
      Integer  IBOND(MXCENT,MXCENT)
      Real(realk) VEC1(3), VEC2(3), VEC3(3), VEC4(3)
      Integer     :: IBNDCN(16)
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175E0_realk*PI/180E0_realk)
      LOGICAL FRGBND, FOUND
!
!     We have to find all dihedral angles between bonded atoms
!
      IF (IATOM .LE. 3) RETURN
      FRGBND = .FALSE.
      INODIH = optinfo%NIntCoord
      IBEFOR = optinfo%NIntCoord
      DO 10 I = 1, IATOM
         DO 15 J = 1, IATOM
            IF (IBOND(I,J) .EQ. 3) FRGBND = .TRUE.
            IF ((IBOND(I,J) .GT. 0) .AND. (IBOND(I,J) .NE. 4)) THEN
               DO 20 K = 1, IATOM
                  IF ((IBOND(J,K) .GT. 0) .AND. (K .NE. I) .AND. &
     &                 (IBOND(J,K) .NE. 4)) THEN 
                     DO 25 L = I, IATOM
                        IF ((IBOND(K,L) .GT. 0) .AND. (L .NE. J) .AND. &
     &                       (L .NE. I) .AND. (IBOND(K,L) .NE. 4)) THEN
!
!     Check for linearity.
!
                           VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
                           VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
                           VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
                           VEC2(1) = ATMARR(K,2)-ATMARR(J,2)
                           VEC2(2) = ATMARR(K,3)-ATMARR(J,3)
                           VEC2(3) = ATMARR(K,4)-ATMARR(J,4)
                           BNDAN1  = vecang_ls(VEC1,VEC2)
                           VEC1(1) = -VEC2(1)
                           VEC1(2) = -VEC2(2)
                           VEC1(3) = -VEC2(3)
                           VEC2(1) = ATMARR(L,2)-ATMARR(K,2)
                           VEC2(2) = ATMARR(L,3)-ATMARR(K,3)
                           VEC2(3) = ATMARR(L,4)-ATMARR(K,4)
                           BNDAN2  = vecang_ls(VEC1,VEC2)
                           IF ((ABS(BNDAN1) .LT. DEG175) &
     &                          .AND. (ABS(BNDAN2) .LT. DEG175)) THEN
                              optinfo%NIntCoord = optinfo%NIntCoord + 1
                              optinfo%INTCRD(optinfo%NIntCoord,1) = 21
                              optinfo%INTCRD(optinfo%NIntCoord,2) = I
                              optinfo%INTCRD(optinfo%NIntCoord,3) = J
                              optinfo%INTCRD(optinfo%NIntCoord,4) = K
                              optinfo%INTCRD(optinfo%NIntCoord,5) = L
                           END IF
                        END IF
 25                  CONTINUE
                  END IF
 20            CONTINUE
            END IF
 15      CONTINUE
 10   CONTINUE
!
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Dihedral angles')
         DO 30 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,4I4)') &
     &           ' Dihedral angles between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3),optinfo%INTCRD(I,4),optinfo%INTCRD(I,5)
 30      CONTINUE
      END IF
!
!     For weakly bound systems, we want to make sure that we add
!     dihedrals across the interfragment bonds (if not already present)
!
      IF (FRGBND) THEN
         DO 40 IFBCRD = 1, optinfo%NIntCoord
            IF (optinfo%INTCRD(IFBCRD,1) .EQ. 3) THEN
               I = optinfo%INTCRD(IFBCRD,2)
               J = optinfo%INTCRD(IFBCRD,3)
               FOUND = .FALSE.
!
!     Check if dihedral(s) already span this interfragment bond
!
               DO 41 ICRD = 1, optinfo%NIntCoord
                  IF (optinfo%INTCRD(ICRD,1) .EQ. 21) THEN
                     IF ((optinfo%INTCRD(ICRD,3) .EQ. I) .AND. &
     &                    (optinfo%INTCRD(ICRD,4) .EQ. J)) &
     &                    FOUND = .TRUE.
                     IF ((optinfo%INTCRD(ICRD,3) .EQ. J) .AND. &
     &                    (optinfo%INTCRD(ICRD,4) .EQ. I)) &
     &                    FOUND = .TRUE.
                  END IF
 41            CONTINUE
!
!     If not, we have to do a wider search. What we are looking for are
!     dihedrals A-B-C-D, where the actual bonding situation is A--B...C--X--D,
!     i.e. BC is the interfragment bond while C and D are both bonded
!     to the same atom X
!
               IF (.NOT. FOUND) THEN
                  DO 42 ICRD = 1, optinfo%NIntCoord
                     IF (optinfo%INTCRD(ICRD,1) .EQ. 1) THEN
                        II1 = 0
                        IF (optinfo%INTCRD(ICRD,2) .EQ. I) THEN
                           II1 = optinfo%INTCRD(ICRD,3)
                           II2 = I
                           II3 = J
                        ELSE IF (optinfo%INTCRD(ICRD,3) .EQ. I) THEN
                           II1 = optinfo%INTCRD(ICRD,2)
                           II2 = I
                           II3 = J
                        ELSE IF (optinfo%INTCRD(ICRD,2) .EQ. J) THEN
                           II1 = optinfo%INTCRD(ICRD,3)
                           II2 = J
                           II3 = I
                        ELSE IF (optinfo%INTCRD(ICRD,3) .EQ. J) THEN
                           II1 = optinfo%INTCRD(ICRD,2)
                           II2 = J
                           II3 = I
                        END IF
                        IF (II1 .GT. 0) THEN
                           DO 43 IAT1 = 1, IATOM
                              IF (IBOND(IAT1,II3) .EQ. 1) THEN
                                 DO 44 IAT2 = 1, IATOM
                                    IF ((IBOND(IAT1,IAT2) .EQ. 1) .AND. &
     &                                   (IAT2 .NE. II3)) THEN
                                       II4 = IAT2
!     Check for linearity.
                                       VEC1(1) = ATMARR(II1,2) &
     &                                      - ATMARR(II2,2)
                                       VEC1(2) = ATMARR(II1,3) &
     &                                      - ATMARR(II2,3)
                                       VEC1(3) = ATMARR(II1,4) &
     &                                      - ATMARR(II2,4)
                                       VEC2(1) = ATMARR(II3,2) &
     &                                      - ATMARR(II2,2)
                                       VEC2(2) = ATMARR(II3,3) &
     &                                      - ATMARR(II2,3)
                                       VEC2(3) = ATMARR(II3,4) &
     &                                      - ATMARR(II2,4)
                                       BNDAN1  = vecang_ls(VEC1,VEC2)
                                       VEC1(1) = -VEC2(1)
                                       VEC1(2) = -VEC2(2)
                                       VEC1(3) = -VEC2(3)
                                       VEC2(1) = ATMARR(II4,2) &
     &                                      - ATMARR(II3,2)
                                       VEC2(2) = ATMARR(II4,3) &
     &                                      - ATMARR(II3,3)
                                       VEC2(3) = ATMARR(II4,4) &
     &                                      - ATMARR(II3,4)
                                       BNDAN2  = vecang_ls(VEC1,VEC2)
                                       IF ((ABS(BNDAN1).LT.DEG175).AND. &
     &                                    (ABS(BNDAN2).LT.DEG175)) THEN
                                          optinfo%NIntCoord = optinfo%NIntCoord + 1
                                          optinfo%INTCRD(optinfo%NIntCoord,1) = 25
                                          optinfo%INTCRD(optinfo%NIntCoord,2) = II1
                                          optinfo%INTCRD(optinfo%NIntCoord,3) = II2
                                          optinfo%INTCRD(optinfo%NIntCoord,4) = II3
                                          optinfo%INTCRD(optinfo%NIntCoord,5) = II4
                                       END IF
                                    END IF
 44                              CONTINUE
                              END IF
 43                        CONTINUE
                        END IF
                     END IF
 42               CONTINUE
               END IF
            END IF
 40      CONTINUE
      END IF
!
!     Systems with more than three atoms where no regular dihedral
!     angles has been found, need addition of dihedral angles to take
!     care of out of plane bending. First we find the central atom.
!
      ICNT = 0
      IBNDNR = 0
      IF ((IATOM .GE. 4) .AND. (IBEFOR .EQ. optinfo%NIntCoord)) THEN
         DO 50 I = 1, IATOM
            INUMBN = 0
            DO 55 J = 1, IATOM
               IF (IBOND(I,J) .GT. 0) INUMBN = INUMBN + 1
 55         CONTINUE
            IF (INUMBN .GE. 3) THEN
               ICNT   = I
               IBNDNR = INUMBN
            END IF
 50      CONTINUE
      END IF
      IF ((ICNT .GT. 0) .AND. (IBNDNR .LE. 16)) THEN
         call ls_IZERO(IBNDCN,16)
         IBEFOR = optinfo%NIntCoord
         II = 1
         DO 60 I = 1, IATOM
            IF (IBOND(I,ICNT) .GT. 0) THEN
               IBNDCN(II) = I
               II = II + 1
            END IF
 60      CONTINUE
         DO 65 II = 1, IBNDNR
            VEC1(1) = ATMARR(IBNDCN(II),2)-ATMARR(ICNT,2)
            VEC1(2) = ATMARR(IBNDCN(II),3)-ATMARR(ICNT,3)
            VEC1(3) = ATMARR(IBNDCN(II),4)-ATMARR(ICNT,4)
            DO 67 JJ = 1, IBNDNR
               IF (JJ .NE. II) THEN
                  VEC2(1) = ATMARR(IBNDCN(JJ),2)-ATMARR(ICNT,2)
                  VEC2(2) = ATMARR(IBNDCN(JJ),3)-ATMARR(ICNT,3)
                  VEC2(3) = ATMARR(IBNDCN(JJ),4)-ATMARR(ICNT,4)
                  BNDAN1  = vecang_ls(VEC1,VEC2)
                  IF (ABS(BNDAN1) .LT. DEG175) THEN
                     VEC3(1) = -VEC2(1)
                     VEC3(2) = -VEC2(2)
                     VEC3(3) = -VEC2(3)
                     DO 69 KK = 1, IBNDNR
                        IF ((KK .NE. II) .AND. (KK .NE. JJ) &
     &                       .AND. (optinfo%NIntCoord .EQ.IBEFOR)) THEN
                           VEC4(1) = ATMARR(IBNDCN(KK),2) &
     &                          -ATMARR(IBNDCN(JJ),2)
                           VEC4(2) = ATMARR(IBNDCN(KK),3) &
     &                          -ATMARR(IBNDCN(JJ),3)
                           VEC4(3) = ATMARR(IBNDCN(KK),4) &
     &                          -ATMARR(IBNDCN(JJ),4)
                           BNDAN2  = vecang_ls(VEC3,VEC4)
                           IF (ABS(BNDAN2) .LT. DEG175) THEN
!
!     We need well-balanced dihedral angles, and once we have found one
!     suitable angle, we rely on permutational symmetry.
!
!     First we have to do some sorting.
!
                              I = II
                              J = JJ
                              K = KK
                              ITMP = IBNDCN(1)
                              IBNDCN(1) = IBNDCN(I)
                              IBNDCN(I) = ITMP
                              IF (J .EQ. 1) J = I
                              IF (K .EQ. 1) K = I
                              I = 1
                              ITMP = IBNDCN(2)
                              IBNDCN(2) = IBNDCN(J)
                              IBNDCN(J) = ITMP
                              IF (K .EQ. 2) K = J
                              J = 2
                              ITMP = IBNDCN(3)
                              IBNDCN(3) = IBNDCN(K)
                              IBNDCN(K) = ITMP
                              K = 3
!
!     Then we find all the dihedral angles
!
                              DO 100 LL = 0, IBNDNR - 1
                                 optinfo%NIntCoord = optinfo%NIntCoord + 1
                                 optinfo%INTCRD(optinfo%NIntCoord,1) = 23
                                 optinfo%INTCRD(optinfo%NIntCoord,2) = &
     &                                IBNDCN(MOD(LL,IBNDNR)+1)
                                 optinfo%INTCRD(optinfo%NIntCoord,3) = ICNT
                                 optinfo%INTCRD(optinfo%NIntCoord,4) = &
     &                                IBNDCN(MOD(LL+1,IBNDNR)+1)
                                 optinfo%INTCRD(optinfo%NIntCoord,5) = &
     &                                IBNDCN(MOD(LL+2,IBNDNR)+1)
 100                          CONTINUE

 123                          CONTINUE

                           END IF
                        END IF
 69                  CONTINUE
                  END IF
               END IF
 67         CONTINUE
 65      CONTINUE
         IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
            call lsheader(lupri,'Extra dihedral angles')
            DO 70 I = IBEFOR + 1, optinfo%NIntCoord
               WRITE(LUPRI,'(A,4I4)') &
     &              ' Dihedral angles between atoms: ', &
     &              optinfo%INTCRD(I,2),optinfo%INTCRD(I,3),optinfo%INTCRD(I,4),optinfo%INTCRD(I,5)
 70         CONTINUE
         END IF
      END IF
!
!     Dihedral angle coordinates that have been manually specified are added
!     (if not already present)
!
      IBEFOR = optinfo%NIntCoord
      IF (optinfo%AddCoord .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         call lsheader(lupri,'Manually added dihedral angles')
      END IF
   If (optinfo%AddCoord) then
      DO 90 I = 1, optinfo%NAdd
         IF ((optinfo%AddCoordArray(I,3) .GT. 0) .AND. (optinfo%AddCoordArray(I,4) .GT. 0)) THEN
            IF ((optinfo%AddCoordArray(I,1) .GT. IATOM) .OR. (optinfo%AddCoordArray(I,2) .GT. IATOM) &
     &           .OR. (optinfo%AddCoordArray(I,3) .GT. IATOM) &
     &           .OR. (optinfo%AddCoordArray(I,4) .GT. IATOM)) THEN
               WRITE(LUPRI,'(A,I4)') &
     &              ' Critical error in additional coordinate #', I
               WRITE(LUPRI,'(A)') ' Chosen atom number  &
     &               is larger than the number of atoms.'
               call lsquit(' Error detected in the specification of  &
     &               additional dihedral coordinate(s)!',lupri)
            END IF
            IBND1 = IBOND(optinfo%AddCoordArray(I,1),optinfo%AddCoordArray(I,2))
            IBND2 = IBOND(optinfo%AddCoordArray(I,2),optinfo%AddCoordArray(I,3))
            IBND3 = IBOND(optinfo%AddCoordArray(I,3),optinfo%AddCoordArray(I,4))
!     Check if the four atoms should already have produced a dihedral coordinate
            IF (.NOT.(((IBND1.GT. 0).AND.(IBND1.NE. 4).AND.(IBND1.NE. 9)) &
     &      .AND.((IBND2.GT. 0).AND.(IBND2.NE. 4).AND.(IBND2.NE. 9)) &
     &      .AND.((IBND3.GT. 0).AND.(IBND3.NE. 4).AND.(IBND3.NE. 9)))) THEN
               optinfo%NIntCoord = optinfo%NIntCoord + 1
               optinfo%INTCRD(optinfo%NIntCoord,1) = 29
               optinfo%INTCRD(optinfo%NIntCoord,2) = optinfo%AddCoordArray(I,1)
               optinfo%INTCRD(optinfo%NIntCoord,3) = optinfo%AddCoordArray(I,2)
               optinfo%INTCRD(optinfo%NIntCoord,4) = optinfo%AddCoordArray(I,3)
               optinfo%INTCRD(optinfo%NIntCoord,5) = optinfo%AddCoordArray(I,4)
            END IF
         END IF
 90   CONTINUE
   Endif
!
      IF ((IBEFOR .NE. optinfo%NIntCoord) .AND. (optinfo%IPrint .GE. IPRMAX)) THEN
         DO 95 I = IBEFOR + 1, optinfo%NIntCoord
            WRITE(LUPRI,'(A,4I4)') &
     &           ' Manual dihedral between atoms: ', &
     &           optinfo%INTCRD(I,2),optinfo%INTCRD(I,3),optinfo%INTCRD(I,4),optinfo%INTCRD(I,5)
 95      CONTINUE
      END IF
      RETURN
      END

!  /* Deck rredun */
      SUBROUTINE LS_RREDUN(lupri,optinfo)
use precision
use ls_util 
use files
use optimization_input 
use Fundamental
!
!     Reduces the redundancy of the redundant internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      NDIHD = 0
      IBEFOR = optinfo%NIntCoord
      DO 10 I = 1, optinfo%NIntCoord
         IF (optinfo%INTCRD(I,1) .GE. 20) NDIHD = NDIHD + 1
 10   CONTINUE
!
!     The number of dihedral angles should be larger than 10 and
!     comprise more than 1/3 of the total number of coordinates.
!
      IF ((NDIHD .LT. 10) .OR. (NDIHD .LT. optinfo%NIntCoord/3)) RETURN
      IFIRST = optinfo%NIntCoord - NDIHD + 1
      IC = IFIRST
!
!     We pick out every fifth dihedral angle.
!
 20   CONTINUE
      IF (MOD(IC*1.0E0_realk, 10.0E0_realk) .GT. 1.0E-3_realk) THEN
         IC = IC + 1
         GOTO 20
      END IF
!
!     Then we mark all the internal coordinates that have the same value.
!
 25   CONTINUE
      VAL = ABS(optinfo%CoordInt(IC))
      IF (VAL .LT. 1.0E6_realk) THEN
         DO 30 I = IFIRST, optinfo%NIntCoord
            IF (ABS(ABS(optinfo%CoordInt(I))-VAL) .LT. 1.0E-8_realk) THEN
               optinfo%CoordInt(I) = 1.1E6_realk
            END IF
 30      CONTINUE
      END IF
      IF (IC+5 .LE. optinfo%NIntCoord) THEN
         IC = IC + 10
         GOTO 25
      END IF
!
!     Then all marked coordinates are removed.
!
      IC = optinfo%NIntCoord - NDIHD + 1
 35   CONTINUE
      IF (IC .LE. optinfo%NIntCoord) THEN
         IF (optinfo%CoordInt(IC) .GT. 1.0E6_realk) THEN
            DO 40 I = IC, optinfo%NIntCoord-1
               DO 42 J = 1, 6
                  optinfo%INTCRD(I,J) = optinfo%INTCRD(I+1,J)
 42            CONTINUE
               optinfo%CoordInt(I) = optinfo%CoordInt(I+1)
 40         CONTINUE
            optinfo%NIntCoord = optinfo%NIntCoord - 1
            GOTO 35
         ELSE
            IC = IC + 1
            GOTO 35
         END IF
      END IF
      NDIHD = NDIHD - (IBEFOR-optinfo%NIntCoord)
      IF ((NDIHD .GT. 10) .AND. (NDIHD .GT. optinfo%NIntCoord/5)) THEN
         IC = optinfo%NIntCoord - NDIHD + 1
         GOTO 20
      END IF
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) IBEFOR-optinfo%NIntCoord, &
     &        ' dihedral angles were removed to reduce redundancy!'
         WRITE(LUPRI,*)
      END IF
      RETURN
      END

!  /* Deck getwil */
      SUBROUTINE LS_GETWIL(Molecule,MXRCRD,MX2CRD,ATMARR,WILBMT,BMTRAN,& 
      & TMPMAT,lupri,IATOM,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use Fundamental
use molecule_type
use molecule_typetype
!
!     Constructs the Wilson B matrix used for transformations between
!     Cartesian and redundant natural internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri,IATOM
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8), WILBMT(MXRCRD,MXCOOR)
      Real(realk) BMTRAN(MXRCRD,MXRCRD), TMPMAT(MX2CRD,MX2CRD)
      Real(realk) VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (DEG175 = 175E0_realk*PI/180E0_realk)
!
      call ls_DZERO(WILBMT,MXRCRD*MXCOOR)
      IF (MXCENT*8.GT.MXCOOR*MXCOOR) CALL lsquit('Rethink ATMARR in LS_GETWIL',lupri)
      call ls_DZERO(ATMARR,8*MXCENT)
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Output from GETWIL')
      END IF
!
!     First we initialize the ATMARR array.
!
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.TRUE.,lupri)
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Expanded array of atoms (Bohr)')
         call output(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      optinfo%ICartCoord = 3*IATOM
!
      IF (optinfo%Rebid) THEN
         optinfo%NIntCoord = -1
      END IF
      NRIC = optinfo%NIntCoord
      IF (optinfo%DelInt .AND. (optinfo%NIntCoord .GT. 0)) NRIC = optinfo%NIntCoord
      DO 20 IC = 1, NRIC
!
!     The connection between Cartesian coordinates and bonds:
!     -------------------------------------------------------
!
         IF (optinfo%INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
            call ls_NRMLVC(VEC1)
            DO 25 II = 1, 3
               WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+II) =  VEC1(II)
               WILBMT(IC,(optinfo%INTCRD(IC,3)-1)*3+II) = -VEC1(II)
 25         CONTINUE
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
               VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
               VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
               BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               DO 27 II = 1, IATOM
                  DO 28 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0E-6_realk
                     VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     BNDL2 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     optinfo%Change = (BNDL2-BNDL1)/1.0E-6_realk
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0E-6_realk
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,optinfo%Change
 28               CONTINUE
 27            CONTINUE
            END IF
!
!     The connection between Cartesian coordinates and angles:
!     --------------------------------------------------------
!
         ELSE IF (optinfo%INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
            VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) - ATMARR(optinfo%INTCRD(IC,3),2)
            VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) - ATMARR(optinfo%INTCRD(IC,3),3)
            VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) - ATMARR(optinfo%INTCRD(IC,3),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            call ls_NRMLVC(VEC1)
            call ls_NRMLVC(VEC2)
!
!     All regular angles
!
            IF (optinfo%INTCRD(IC,1) .EQ. 11) THEN
               call ls_VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0E-8_realk) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0E-8_realk) THEN
                  VEC3(1) =  VEC1(2)-VEC1(3)
                  VEC3(2) = -VEC1(1)-VEC1(3)
                  VEC3(3) =  VEC1(1)+VEC1(2)
               END IF
               call ls_NRMLVC(VEC3)
!
!     Second coordinate of angles larger than 175 degrees.
!
            ELSE
               call ls_VECPRD(VEC1,VEC3,VEC4)
               call ls_NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
!
            call ls_VECPRD(VEC1,VEC3,VEC4)
            call ls_VECPRD(VEC3,VEC2,VEC5)
            call ls_NRMLVC(VEC4)
            call ls_NRMLVC(VEC5)
            DO 30 II = 1, 3
               WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+II) = VEC4(II)/BNDL1
               WILBMT(IC,(optinfo%INTCRD(IC,3)-1)*3+II) = &
     &                          -VEC4(II)/BNDL1 - VEC5(II)/BNDL2
               WILBMT(IC,(optinfo%INTCRD(IC,4)-1)*3+II) = VEC5(II)/BNDL2
 30         CONTINUE
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
               VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
               VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
               VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) - ATMARR(optinfo%INTCRD(IC,3),2)
               VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) - ATMARR(optinfo%INTCRD(IC,3),3)
               VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) - ATMARR(optinfo%INTCRD(IC,3),4)
               call ls_VECPRD(VEC3,VEC1,VEC4)
               call ls_NRMLVC(VEC4)
               ANGLE = vecang_ls(VEC1,VEC2)
               IF ((optinfo%INTCRD(IC,1) .EQ. 12) .OR. (optinfo%INTCRD(IC+1,1) .EQ. 12)) &
     &              ANGLE = vecang_ls(VEC1,VEC4) + vecang_ls(VEC4,VEC2)
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               DO 47 II = 1, IATOM
                  DO 48 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0E-6_realk
                     VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     ANGLE2 = vecang_ls(VEC1,VEC2)
                     IF ((optinfo%INTCRD(IC,1) .EQ. 12) .OR. &
     &                  (optinfo%INTCRD(IC+1,1) .EQ. 12)) &
     &                    ANGLE2 = vecang_ls(VEC1,VEC4) + vecang_ls(VEC4,VEC2)
                     optinfo%Change = (ANGLE2-ANGLE)/1.0E-6_realk
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0E-6_realk
                     WRITE(LUPRI,'(2I5,F16.6)') II,JJ,optinfo%Change
 48               CONTINUE
 47            CONTINUE
            END IF
!
!     The connection between Cartesian coordinates and dihedral angles:
!     -----------------------------------------------------------------
!
         ELSE
            VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
            VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) - ATMARR(optinfo%INTCRD(IC,3),2)
            VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) - ATMARR(optinfo%INTCRD(IC,3),3)
            VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) - ATMARR(optinfo%INTCRD(IC,3),4)
            VEC3(1) = ATMARR(optinfo%INTCRD(IC,5),2) - ATMARR(optinfo%INTCRD(IC,4),2)
            VEC3(2) = ATMARR(optinfo%INTCRD(IC,5),3) - ATMARR(optinfo%INTCRD(IC,4),3)
            VEC3(3) = ATMARR(optinfo%INTCRD(IC,5),4) - ATMARR(optinfo%INTCRD(IC,4),4)
            BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
            BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
            BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
            call ls_NRMLVC(VEC1)
            call ls_NRMLVC(VEC2)
            call ls_NRMLVC(VEC3)
            UW =  DDOT(3,VEC1,1,VEC2,1)
            VW =  DDOT(3,VEC3,1,VEC2,1)
            call ls_VECPRD(VEC1,VEC2,VEC4)
            call ls_VECPRD(VEC3,VEC2,VEC5)
!
!     Check for linearity. The dihedral angle is undefined if three
!     or four of the atoms are linear.
!
            IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0E-16_realk) .AND. &
     &           (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0E-16_realk)) THEN
               DO 50 II = 1, 3
!
!     We construct a few "building blocks".
!
                  CMPU1 = VEC4(II)/(BNDL1*(1.0E0_realk-UW*UW))
                  CMPV1 = VEC5(II)/(BNDL3*(1.0E0_realk-VW*VW))
                  CMPU2 = VEC4(II)*UW/(BNDL2*(1.0E0_realk-UW*UW))
                  CMPV2 = VEC5(II)*VW/(BNDL2*(1.0E0_realk-VW*VW))
                  WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+II) =  CMPU1
                  WILBMT(IC,(optinfo%INTCRD(IC,3)-1)*3+II) = -CMPU1+CMPU2-CMPV2
                  WILBMT(IC,(optinfo%INTCRD(IC,4)-1)*3+II) =  CMPV1+CMPV2-CMPU2
                  WILBMT(IC,(optinfo%INTCRD(IC,5)-1)*3+II) = -CMPV1
 50            CONTINUE
            END IF
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
               VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
               VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
               VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) - ATMARR(optinfo%INTCRD(IC,3),2)
               VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) - ATMARR(optinfo%INTCRD(IC,3),3)
               VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) - ATMARR(optinfo%INTCRD(IC,3),4)
               VEC3(1) = ATMARR(optinfo%INTCRD(IC,5),2) - ATMARR(optinfo%INTCRD(IC,4),2)
               VEC3(2) = ATMARR(optinfo%INTCRD(IC,5),3) - ATMARR(optinfo%INTCRD(IC,4),3)
               VEC3(3) = ATMARR(optinfo%INTCRD(IC,5),4) - ATMARR(optinfo%INTCRD(IC,4),4)
               call ls_NRMLVC(VEC2)
               CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
               VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
               VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
               VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
               call ls_NRMLVC(VEC1)
               CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
               VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
               VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
               VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
               call ls_NRMLVC(VEC3)
               DIHED = vecang_ls(VEC1,VEC3)
               IF (ABS(DIHED) .GT. DEG175) THEN
                  call ls_VECPRD(VEC1,VEC2,VEC4)
                  call ls_NRMLVC(VEC4)
                  DIHED = vecang_ls(VEC1,VEC4)+vecang_ls(VEC4,VEC3)
               END IF
               VEC5(1) = VEC1(1)
               VEC5(2) = VEC1(2)
               VEC5(3) = VEC1(3)
               VEC6(1) = VEC3(1)
               VEC6(2) = VEC3(2)
               VEC6(3) = VEC3(3)
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Numerical results:'
               WRITE(LUPRI,*)
               FAC = -D1
               DO 97 II = 1, IATOM
                  DO 98 JJ = 1, 3
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) + 1.0E-6_realk
                     VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(optinfo%INTCRD(IC,5),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(optinfo%INTCRD(IC,5),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(optinfo%INTCRD(IC,5),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),4)
                     call ls_NRMLVC(VEC2)
                     CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
                     VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
                     VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
                     VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
                     call ls_NRMLVC(VEC1)
                     CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
                     VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
                     VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
                     VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
                     call ls_NRMLVC(VEC3)
                     IF (FAC .LT. D0) THEN
                        FAC = D1
                        call ls_VECPRD(VEC1,VEC3,VEC4)
                        IF (DDOT(3,VEC4,1,VEC2,1) .LT. D0) FAC = -D1
                     END IF
                     DIHED2 = vecang_ls(VEC1,VEC3)
                     CH1 = FAC*(vecang_ls(VEC1,VEC6) &
     &                    -vecang_ls(VEC5,VEC6))/1.0E-6_realk
                     CH2 = FAC*(vecang_ls(VEC3,VEC5) &
     &                    -vecang_ls(VEC6,VEC5))/1.0E-6_realk
                     IF (ABS(DIHED2) .GT. DEG175) THEN
                        call ls_VECPRD(VEC1,VEC2,VEC4)
                        call ls_NRMLVC(VEC4)
                        DIHED2 = vecang_ls(VEC1,VEC4)+vecang_ls(VEC4,VEC3)
                     END IF
                     optinfo%Change = FAC*(DIHED2-DIHED)/1.0E-6_realk
                     ATMARR(II,JJ+1) = ATMARR(II,JJ+1) - 1.0E-6_realk
                     WRITE(LUPRI,'(2I5,3F16.6)') II,JJ,optinfo%Change,CH1,CH2
 98               CONTINUE
 97            CONTINUE
            END IF
         END IF
!
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            TXT = 'Internal coordinate number:     '
            WRITE(TXT(28:32),'(I5)') IC
            call lsheader(lupri,TXT)
            DO 99 K = 1, IATOM
               WRITE(LUPRI,'(I25,A,F16.6)') &
     &              K,'x',WILBMT(IC,(K-1)*3+1)
               WRITE(LUPRI,'(I25,A,F16.6)') &
     &              K,'y',WILBMT(IC,(K-1)*3+2)
               WRITE(LUPRI,'(I25,A,F16.6)') &
     &              K,'z',WILBMT(IC,(K-1)*3+3)
 99         CONTINUE
         END IF
 20   CONTINUE
!
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Wilson B matrix [B(ij) = dq(i)/dx(j)]')
         call output(WILBMT,1,NRIC,1,optinfo%ICartCoord,MXRCRD,MXCOOR,1,LUPRI)
      END IF
!
!     If we're using delocalized internal coordinates, we make
!     the transformation to the active space.
!
      IF (optinfo%DelInt .AND. (optinfo%NIntCoord .GT. 0)) THEN
         call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
         DO 100 I = 1, NRIC
            DO 102 J = 1, optinfo%ICartCoord
               TMPMAT(I,J) = WILBMT(I,J)
 102        CONTINUE
 100     CONTINUE
         call ls_DZERO(WILBMT,MXRCRD*MXCOOR)
         DO 120 I = 1, optinfo%NIntCoord
            DO 122 J = 1, optinfo%ICartCoord
               DO 124 K = 1, optinfo%NIntCoord
                  WILBMT(I,J) = WILBMT(I,J) + BMTRAN(K,I)*TMPMAT(K,J)
 124           CONTINUE
 122        CONTINUE
 120     CONTINUE
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            call lsheader(lupri,'Active space B matrix')
            call output(WILBMT,1,optinfo%NIntCoord,1,optinfo%ICartCoord, &
     &           MXRCRD,MXCOOR,1,LUPRI)
         END IF
      END IF
      RETURN
      END

!  /* Deck getdwl */
      SUBROUTINE LS_GETDWL(Molecule,MXRCRD,DWILBM,ATMARR,TMPMAT, & 
      & WILBMT,lupri,IATOM,optinfo)
use precision
use ls_util 
use files
use optimization_input 
use molecule_type
use molecule_typetype
!
!     Constructs the derivative of the Wilson B matrix used for transformations
!     between Cartesian and redundant natural internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) DWILBM(MXCOOR,MXCOOR), WILBMT(MXRCRD,MXCOOR)
      Real(realk) ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      Real(realk) VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3), VEC6(3)
      CHARACTER TXT*32
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      call ls_DZERO(ATMARR,8*MXCENT)
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Output from GETDWL')
      END IF
!
!     First we initialize the ATMARR array.
!
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.TRUE.,lupri)
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Expanded array of atoms (Bohr)')
         call output(ATMARR,1,IATOM,1,4,MXCENT,8,1,LUPRI)
      END IF
      optinfo%ICartCoord = 3*IATOM
!
      DO 20 IC = 1, optinfo%NIntCoord
         call ls_DZERO(TMPMAT,MXCOOR*MXCOOR)
         call ls_DZERO(DWILBM,MXCOOR*MXCOOR)
         IF (optinfo%INTCRD(IC,1) .LT. 10) THEN
!
!     Numerically:
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               DO 30 II = 2, 3
                  DO 35 JJ = 1, 3
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) + 1.0E-6_realk
                     VEC1(1)= ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2)= ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3)= ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     call ls_NRMLVC(VEC1)
                     call ls_DZERO(VEC2,3)
                     DO 40 J = 1, 3
                        VEC2(J) = 1.0E0_realk
                        COMPNT = DDOT(3,VEC1,1,VEC2,1)
                        TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                       (optinfo%INTCRD(IC,3)-1)*3+J) = -COMPNT
                        TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                       (optinfo%INTCRD(IC,2)-1)*3+J) =  COMPNT
                        VEC2(J) = 0.0E0_realk
 40                  CONTINUE
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) - 1.0E-6_realk
 35               CONTINUE
 30            CONTINUE
               DO 50 II = 2, 3
                  DO 51 JJ = 1, 3
                     DO 53 KK = 2, 3
                        DO 54 LL = 1, 3
                           DWILBM((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) = &
     &                          (TMPMAT((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) - &
     &                          WILBMT(IC,(optinfo%INTCRD(IC,II)-1)*3+JJ)) &
     &                          /1.0E-6_realk
 54                     CONTINUE
 53                  CONTINUE
 51               CONTINUE
 50            CONTINUE
               DO 56 II = 1, optinfo%ICartCoord
                  DO 58 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0E0_realk
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0E-6_realk) &
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 58               CONTINUE
 56            CONTINUE
               call lsheader(lupri,'Numerically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
!
!     Analytically:
!
            call ls_GTDWLM(Molecule,MXRCRD,IC,DWILBM,ATMARR,WILBMT,lupri,optinfo)
            IF (optinfo%IPrint .GE. IPRMAX) THEN
               call lsheader(lupri,'Analytically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
!
         ELSE IF (optinfo%INTCRD(IC,1) .LT. 20) THEN
!
!     Numerically:
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               DO 60 II = 2, 4
                  DO 65 JJ = 1, 3
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) + 1.0E-6_realk
                     VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     call ls_NRMLVC(VEC1)
                     call ls_NRMLVC(VEC2)
!
!     All regular angles
!
                     IF (optinfo%INTCRD(IC,1) .EQ. 11) THEN
                        call ls_VECPRD(VEC1,VEC2,VEC3)
                        VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
                        IF (VNRM .LE. 1.0E-8_realk) THEN
                           VEC3(1) =  VEC1(2)+VEC1(3)
                           VEC3(2) = -VEC1(1)+VEC1(3)
                           VEC3(3) = -VEC1(1)-VEC1(2)
                        END IF
                        call ls_NRMLVC(VEC3)
!
!     Second coordinate of angles larger than 175 degrees.
!
                     ELSE
                        call ls_VECPRD(VEC1,VEC3,VEC4)
                        call ls_NRMLVC(VEC4)
                        VEC3(1) = VEC4(1)
                        VEC3(2) = VEC4(2)
                        VEC3(3) = VEC4(3)
                     END IF
!
                     call ls_VECPRD(VEC1,VEC3,VEC4)
                     call ls_NRMLVC(VEC4)
                     call ls_VECPRD(VEC3,VEC2,VEC5)
                     call ls_NRMLVC(VEC5)
!
                     DO 70 J = 1, 3
                        TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                       (optinfo%INTCRD(IC,2)-1)*3+J) =  (VEC4(J)/BNDL1)
                        TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                       (optinfo%INTCRD(IC,3)-1)*3+J) = -(VEC4(J)/BNDL1) &
     &                       - (VEC5(J)/BNDL2)
                        TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                       (optinfo%INTCRD(IC,4)-1)*3+J) =  (VEC5(J)/BNDL2)
 70                  CONTINUE
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) - 1.0E-6_realk
 65               CONTINUE
 60            CONTINUE
!
               DO 90 II = 2, 4
                  DO 91 JJ = 1, 3
                     DO 93 KK = 2, 4
                        DO 94 LL = 1, 3
                           DWILBM((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) = &
     &                          (TMPMAT((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) - &
     &                          WILBMT(IC,(optinfo%INTCRD(IC,II)-1)*3+JJ)) &
     &                          /1.0E-6_realk
 94                     CONTINUE
 93                  CONTINUE
 91               CONTINUE
 90            CONTINUE
               DO 96 II = 1, optinfo%ICartCoord
                  DO 98 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0E0_realk
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0E-6_realk) &
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 98               CONTINUE
 96            CONTINUE
               call lsheader(lupri,'Numerically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
!
!     Analytically:
!
            call ls_GTDWLM(Molecule,MXRCRD,IC,DWILBM,ATMARR,WILBMT,lupri,optinfo)
            IF (optinfo%IPrint .GE. IPRMAX) THEN
               call lsheader(lupri,'Analytically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
!
         ELSE
!
!     Numerically:
!
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               DO 100 II = 2, 5
                  DO 105 JJ = 1, 3
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) + 1.0E-6_realk
                     VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),2)
                     VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),3)
                     VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,3),4)
                     VEC3(1) = ATMARR(optinfo%INTCRD(IC,5),2) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),2)
                     VEC3(2) = ATMARR(optinfo%INTCRD(IC,5),3) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),3)
                     VEC3(3) = ATMARR(optinfo%INTCRD(IC,5),4) &
     &                    - ATMARR(optinfo%INTCRD(IC,4),4)
                     VNRM1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
                     VNRM2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
                     VNRM3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
                     call ls_NRMLVC(VEC1)
                     call ls_NRMLVC(VEC2)
                     call ls_NRMLVC(VEC3)
                     UW =  DDOT(3,VEC1,1,VEC2,1)
                     VW =  DDOT(3,VEC3,1,VEC2,1)
                     call ls_VECPRD(VEC1,VEC2,VEC4)
                     call ls_VECPRD(VEC2,VEC3,VEC5)
                     call ls_VECPRD(VEC3,VEC2,VEC6)
!
!     Check for linearity. The dihedral angle is undefined if three
!     or four of the atoms are linear.
!
                     IF ((DDOT(3,VEC4,1,VEC4,1) .GT. 1.0E-16_realk) .AND. &
     &                    (DDOT(3,VEC5,1,VEC5,1) .GT. 1.0E-16_realk)) THEN
                        DO 120 J = 1, 3
                           TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                          (optinfo%INTCRD(IC,2)-1)*3+J) &
     &                          = VEC4(J)/(VNRM1*(1.0E0_realk-UW*UW))
                           TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                          (optinfo%INTCRD(IC,3)-1)*3+J) = &
     &                          (VEC4(J)*(VNRM1*UW-VNRM2))/ &
     &                          (VNRM1*VNRM2*(1.0E0_realk-UW*UW)) &
     &                          - (VW*VEC6(J))/(VNRM2*(1.0E0_realk-VW*VW))
                           TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                          (optinfo%INTCRD(IC,4)-1)*3+J) = &
     &                          (VEC6(J)*(VNRM3*VW+VNRM2))/ &
     &                          (VNRM3*VNRM2*(1.0E0_realk-VW*VW)) &
     &                          - (UW*VEC4(J))/(VNRM2*(1.0E0_realk-UW*UW))
                           TMPMAT((optinfo%INTCRD(IC,II)-1)*3+JJ, &
     &                          (optinfo%INTCRD(IC,5)-1)*3+J) = &
     &                          VEC5(J)/(VNRM3*(1.0E0_realk-VW*VW))
 120                    CONTINUE
                     END IF
                     ATMARR(optinfo%INTCRD(IC,II),JJ+1) = &
     &                    ATMARR(optinfo%INTCRD(IC,II),JJ+1) - 1.0E-6_realk
 105              CONTINUE
 100           CONTINUE
!
               DO 200 II = 2, 5
                  DO 210 JJ = 1, 3
                     DO 220 KK = 2, 5
                        DO 230 LL = 1, 3
                           DWILBM((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) = &
     &                          (TMPMAT((optinfo%INTCRD(IC,KK)-1)*3+LL, &
     &                          (optinfo%INTCRD(IC,II)-1)*3+JJ) - &
     &                          WILBMT(IC,(optinfo%INTCRD(IC,II)-1)*3+JJ)) &
     &                          /1.0E-6_realk
 230                    CONTINUE
 220                 CONTINUE
 210              CONTINUE
 200           CONTINUE
               DO 250 II = 1, optinfo%ICartCoord
                  DO 255 JJ = 1, II
                     DWILBM(II,JJ) = (DWILBM(II,JJ)+DWILBM(JJ,II))/2.0E0_realk
                     IF (ABS(DWILBM(II,JJ)) .LT. 1.0E-6_realk) &
     &                    DWILBM(II,JJ) = D0
                     DWILBM(JJ,II) = DWILBM(II,JJ)
 255              CONTINUE
 250           CONTINUE
               call lsheader(lupri,'Numerically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
!
!     Analytically:
!
            call ls_GTDWLM(Molecule,MXRCRD,IC,DWILBM,ATMARR,WILBMT,lupri,optinfo)
            IF (optinfo%IPrint .GE. IPRMAX) THEN
               call lsheader(lupri,'Analytically diff. Wilson B matrix')
               call output(DWILBM,1,optinfo%ICartCoord,1,optinfo%ICartCoord,MXCOOR, &
     &              MXCOOR,1,LUPRI)
            END IF
         END IF
 20   CONTINUE
      RETURN
      END

!  /* Deck gtdwl0 */
      SUBROUTINE LS_GTDWL0(Molecule,MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN, &
     &     lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use memory_handling
use molecule_type
use molecule_typetype
use precision
!
!     Ths subroutine determines the coordinate system in work, and
!     allocates memory if needed.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      Real(realk), pointer :: DEL(:),RED(:)
!     Allocate some memory
      Call mem_alloc(DEL,optinfo%NIntCoord*optinfo%ICartCoord*optinfo%ICartCoord) 
      Call mem_alloc(RED,optinfo%NIntCoord*optinfo%ICartCoord*optinfo%ICartCoord) 
!
      IF (.NOT. optinfo%DelInt) THEN
         call ls_GTDWLM(Molecule,MXRCRD,IC,DWILBM,ATMARR,WILBMT,lupri,optinfo)
      ELSE
         call ls_GTDWL1(Molecule,optinfo,MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN, &
     &   RED,DEL,optinfo%NIntCoord,optinfo%NIntCoord, &
     &   optinfo%ICartCoord,lupri)
      END IF
!     Deallocate memory
      Call mem_dealloc(DEL) 
      Call mem_dealloc(RED) 
!
      RETURN
      END

!  /* Deck gtdwl1 */
      SUBROUTINE LS_GTDWL1(Molecule,optinfo,MXRCRD,IC,DWILBM,ATMARR,WILBMT,BMTRAN,& 
          REDMAT,DELMAT,IRED,IDEL,ICRT,lupri)
use precision
use ls_util 
use files
use optimization_input 
use molecule_type
use molecule_typetype
!
!     Ths subroutine acts as a buffer for the subroutine GTDWLM when
!     delocalized internal coordinates are used.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) REDMAT(IRED,ICRT,ICRT)
      Real(realk) DELMAT(IDEL,ICRT,ICRT)
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk)
!
!     We only calculate the redundant derivatives once.
!
      IF (IC .EQ. 1) THEN
         call ls_DZERO(REDMAT,IRED*ICRT*ICRT)
         call ls_DZERO(DELMAT,IDEL*ICRT*ICRT)
         DO 10 IR = 1, IRED
            call ls_GTDWLM(Molecule,MXRCRD,IR,DWILBM,ATMARR, &
            & WILBMT,lupri,optinfo)
            DO 20 I = 1, ICRT
               DO 22 J = 1, ICRT
                  REDMAT(IR,I,J) = DWILBM(I,J)
 22            CONTINUE
 20         CONTINUE
 10      CONTINUE
!
!     Then we transform it to delocalized coordinates.
!
         DO 30 I1 = 1, IDEL
            DO 32 J1 = 1, IRED
               DO 34 K1 = 1, ICRT
                  DO 36 K2 = 1, ICRT
                     DELMAT(I1,K1,K2) = DELMAT(I1,K1,K2) &
     &                    + BMTRAN(J1,I1)*REDMAT(J1,K1,K2)
 36               CONTINUE
 34            CONTINUE
 32         CONTINUE
 30      CONTINUE
      END IF
      call ls_DZERO(DWILBM,MXCOOR*MXCOOR)
      DO 50 I = 1, ICRT
         DO 52 J = 1, ICRT
            DWILBM(I,J) = DELMAT(IC,I,J)
 52      CONTINUE
 50   CONTINUE
      RETURN
      END

!  /* Deck gtdwlm */
      SUBROUTINE LS_GTDWLM(Molecule,MXRCRD,IC,DWILBM,ATMARR, &
      & WILBMT,lupri,optinfo)
use precision
use ls_util 
use files
use optimization_input 
use molecule_type
use molecule_typetype
!
!     Returns the derivative of the Wilson B matrix for _one_
!     redundant internal coordinate (IC).
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) DWILBM(MXCOOR,MXCOOR), ATMARR(MXCENT,8)
      Real(realk) WILBMT(MXRCRD,MXCOOR)
      Real(realk) VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(3)
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk)
!
      call ls_DZERO(DWILBM,MXCOOR*MXCOOR)
      IATOM = optinfo%ICartCoord / 3
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.TRUE.,lupri)
!
!     The derivative for bonds:
!     -------------------------
!
      IF (optinfo%INTCRD(IC,1) .LT. 10) THEN
         VEC1(1)=ATMARR(optinfo%INTCRD(IC,2),2)-ATMARR(optinfo%INTCRD(IC,3),2)
         VEC1(2)=ATMARR(optinfo%INTCRD(IC,2),3)-ATMARR(optinfo%INTCRD(IC,3),3)
         VEC1(3)=ATMARR(optinfo%INTCRD(IC,2),4)-ATMARR(optinfo%INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         call ls_NRMLVC(VEC1)
         DO 10 II = 1, 3
            DO 15 JJ = 1, 3
!
!     The variable FACIJ makes the distinction between "diagonal" and
!     "off-diagonal" elements, where "diagonal" refers to Cartesian
!     displacements (II = JJ) and not the atoms.
!
               FACIJ = D0
               IF (II .EQ. JJ) FACIJ = D1
               CMPUU = (FACIJ - VEC1(II)*VEC1(JJ))/BNDL1
               DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,2)-1)*3+JJ) &
     &              =  CMPUU
               DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &              = -CMPUU
               DWILBM((optinfo%INTCRD(IC,3)-1)*3+JJ,(optinfo%INTCRD(IC,2)-1)*3+II) &
     &          = DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ)
               DWILBM((optinfo%INTCRD(IC,3)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &              =  CMPUU
 15         CONTINUE
 10      CONTINUE
!
!     The derivative for angles:
!     --------------------------
!
      ELSE IF (optinfo%INTCRD(IC,1) .LT. 20) THEN
         VEC1(1)=ATMARR(optinfo%INTCRD(IC,2),2)-ATMARR(optinfo%INTCRD(IC,3),2)
         VEC1(2)=ATMARR(optinfo%INTCRD(IC,2),3)-ATMARR(optinfo%INTCRD(IC,3),3)
         VEC1(3)=ATMARR(optinfo%INTCRD(IC,2),4)-ATMARR(optinfo%INTCRD(IC,3),4)
         VEC2(1)=ATMARR(optinfo%INTCRD(IC,4),2)-ATMARR(optinfo%INTCRD(IC,3),2)
         VEC2(2)=ATMARR(optinfo%INTCRD(IC,4),3)-ATMARR(optinfo%INTCRD(IC,3),3)
         VEC2(3)=ATMARR(optinfo%INTCRD(IC,4),4)-ATMARR(optinfo%INTCRD(IC,3),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         call ls_NRMLVC(VEC1)
         call ls_NRMLVC(VEC2)
         COSA = DDOT(3,VEC1,1,VEC2,1)
         SINA = SQRT(1-COSA*COSA)
!
!     We can only calculate the derivative, if sin(q) is non-zero, as
!     all expressions contains a division by sin(q).
!
         IF (ABS(SINA) .GT. 1.0E-4_realk) THEN
            DO 25 II = 1, 3
               DADMI = WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(optinfo%INTCRD(IC,4)-1)*3+II)
               DADOI = WILBMT(IC,(optinfo%INTCRD(IC,3)-1)*3+II)
               DO 27 JJ = 1, 3
                  DADMJ = WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(optinfo%INTCRD(IC,4)-1)*3+JJ)
                  DADOJ = WILBMT(IC,(optinfo%INTCRD(IC,3)-1)*3+JJ)
!
!     The variable FACIJ makes the distinction between "diagonal" and
!     "off-diagonal" elements, where "diagonal" refers to Cartesian
!     displacements (II = JJ) and not the atoms.
!
                  FACIJ = D0
                  IF (II .EQ. JJ) FACIJ = D1
!
!     We only need a few "building blocks" to make all expressions.
!
                  CMPUU  = (VEC1(II)*VEC2(JJ) + VEC1(JJ)*VEC2(II) &
     &                 - 3.0E0_realk*COSA*VEC1(II)*VEC1(JJ) &
     &                 + FACIJ*COSA)/(SINA*BNDL1*BNDL1)
                  CMPVV  = (VEC2(II)*VEC1(JJ) + VEC2(JJ)*VEC1(II) &
     &                 - 3.0E0_realk*COSA*VEC2(II)*VEC2(JJ) &
     &                 + FACIJ*COSA)/(SINA*BNDL2*BNDL2)
                  CMPWW1 = (VEC1(II)*VEC1(JJ) + VEC2(II)*VEC2(JJ) &
     &                 - COSA*VEC1(II)*VEC2(JJ) &
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
                  CMPWW2 = (VEC2(II)*VEC2(JJ) + VEC1(II)*VEC1(JJ) &
     &                 - COSA*VEC2(II)*VEC1(JJ) &
     &                 - FACIJ)/(SINA*BNDL1*BNDL2)
!
!     All the different elements are constructed.
!
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,2)-1)*3+JJ) &
     &             = CMPUU - DADMI*DADMJ*COSA/SINA
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &             = -CMPUU - CMPWW1 - DADMI*DADOJ*COSA/SINA
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+JJ,(optinfo%INTCRD(IC,2)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &             = CMPWW1 - DADMI*DADNJ*COSA/SINA
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+JJ,(optinfo%INTCRD(IC,2)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &             = CMPUU + CMPVV + CMPWW1 + CMPWW2 &
     &             - DADOI*DADOJ*COSA/SINA
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &             = -CMPVV - CMPWW2 - DADNI*DADOJ*COSA/SINA
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+JJ,(optinfo%INTCRD(IC,4)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,4)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &             = CMPVV - DADNI*DADNJ*COSA/SINA
 27            CONTINUE
 25         CONTINUE
!
!     If the second derivative is undefined (sin(q) ~ 0), we set the
!     matrix equal to the zero matrix.
!
         ELSE
            DO 30 I = 1, optinfo%ICartCoord
               DO 32 J = 1, optinfo%ICartCoord
                  DWILBM(I,J) = D0
 32            CONTINUE
 30         CONTINUE
         END IF
!
!     The derivative for dihedral angles:
!     -----------------------------------
!
      ELSE 
         VEC1(1) = ATMARR(optinfo%INTCRD(IC,2),2) - ATMARR(optinfo%INTCRD(IC,3),2)
         VEC1(2) = ATMARR(optinfo%INTCRD(IC,2),3) - ATMARR(optinfo%INTCRD(IC,3),3)
         VEC1(3) = ATMARR(optinfo%INTCRD(IC,2),4) - ATMARR(optinfo%INTCRD(IC,3),4)
         VEC2(1) = ATMARR(optinfo%INTCRD(IC,4),2) - ATMARR(optinfo%INTCRD(IC,3),2)
         VEC2(2) = ATMARR(optinfo%INTCRD(IC,4),3) - ATMARR(optinfo%INTCRD(IC,3),3)
         VEC2(3) = ATMARR(optinfo%INTCRD(IC,4),4) - ATMARR(optinfo%INTCRD(IC,3),4)
         VEC3(1) = ATMARR(optinfo%INTCRD(IC,5),2) - ATMARR(optinfo%INTCRD(IC,4),2)
         VEC3(2) = ATMARR(optinfo%INTCRD(IC,5),3) - ATMARR(optinfo%INTCRD(IC,4),3)
         VEC3(3) = ATMARR(optinfo%INTCRD(IC,5),4) - ATMARR(optinfo%INTCRD(IC,4),4)
         BNDL1 = SQRT(DDOT(3,VEC1,1,VEC1,1))
         BNDL2 = SQRT(DDOT(3,VEC2,1,VEC2,1))
         BNDL3 = SQRT(DDOT(3,VEC3,1,VEC3,1))
         call ls_NRMLVC(VEC1)
         call ls_NRMLVC(VEC2)
         call ls_NRMLVC(VEC3)
         call ls_VECPRD(VEC1,VEC2,VEC4)
         call ls_VECPRD(VEC3,VEC2,VEC5)
         UW   = DDOT(3,VEC1,1,VEC2,1)
         VW   = DDOT(3,VEC3,1,VEC2,1)
         UWVW = DDOT(3,VEC4,1,VEC5,1)
         UVW  = DDOT(3,VEC1,1,VEC5,1)
!
!     We can only calculate the derivative, if sin(q) (or (u * (v x w)))
!     is non-zero, as all expressions contains a division by sin(q).
!
         IF (ABS(UVW) .GT. 1.0E-6_realk) THEN
            DO 40 II = 1, 3
               DADMI = WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+II)
               DADNI = WILBMT(IC,(optinfo%INTCRD(IC,5)-1)*3+II)
               DO 45 JJ = 1, 3
                  DADMJ = WILBMT(IC,(optinfo%INTCRD(IC,2)-1)*3+JJ)
                  DADNJ = WILBMT(IC,(optinfo%INTCRD(IC,5)-1)*3+JJ)
!
!     The variable FACIJ makes the distinction between "diagonal" and
!     "off-diagonal" elements, where "diagonal" refers to Cartesian
!     displacements (II = JJ) and not the atoms. The sign of this
!     factor follows this table:
!
!                                 i
!                              1  2  3
!                             ---------
!                           1| 0  +  -
!                         j 2| -  0  +
!                           3| +  -  0
!
                  FACIJ = (JJ-II)*(-0.5E0_realk)**(ABS(JJ-II))
!
!     For some of the elements we need the third Cartesian direction,
!     that is k =/ i,j.
!
                  KK = MAX(1,MIN(6-II-JJ,3))
!
!     With the following "building blocks" all expressions can be constructed:
!
                  CMPUU  = (VEC4(II)*(UW*VEC2(JJ)-VEC1(JJ)) &
     &                   + VEC4(JJ)*(UW*VEC2(II)-VEC1(II)))/ &
     &                     (BNDL1*BNDL1*(1.0E0_realk-UW*UW)**2)
                  CMPVV  = (VEC5(II)*(VW*VEC2(JJ)-VEC3(JJ)) &
     &                   + VEC5(JJ)*(VW*VEC2(II)-VEC3(II)))/ &
     &                     (BNDL3*BNDL3*(1.0E0_realk-VW*VW)**2)
                  CMPUW  = 0.5E0_realk*(VEC4(II)*(VEC2(JJ)-2.0E0_realk*UW*VEC1(JJ) &
     &                   + UW*UW*VEC2(JJ))+VEC4(JJ)*(VEC2(II) &
     &                   - 2.0E0_realk*UW*VEC1(II)+UW*UW*VEC2(II)))/ &
     &                     (BNDL1*BNDL2*(1.0E0_realk-UW*UW)**2)
                  CMPVW  = 0.5E0_realk*(VEC5(II)*(VEC2(JJ)-2.0E0_realk*VW*VEC3(JJ) &
     &                   + VW*VW*VEC2(JJ))+VEC5(JJ)*(VEC2(II) &
     &                   - 2.0E0_realk*VW*VEC3(II)+VW*VW*VEC2(II)))/ &
     &                     (BNDL3*BNDL2*(1.0E0_realk-VW*VW)**2)
                  CMPWW1 = 0.5E0_realk*(VEC4(II)*(VEC1(JJ)+UW*UW*VEC1(JJ) &
     &                   - 3.0E0_realk*UW*VEC2(JJ)+UW*UW*UW*VEC2(JJ)) &
     &                   + VEC4(JJ)*(VEC1(II)+UW*UW*VEC1(II) &
     &                   - 3.0E0_realk*UW*VEC2(II)+UW*UW*UW*VEC2(II)))/ &
     &                     (BNDL2*BNDL2*(1.0E0_realk-UW*UW)**2)
                  CMPWW2 = 0.5E0_realk*(VEC5(II)*(VEC3(JJ)+VW*VW*VEC3(JJ) &
     &                   - 3.0E0_realk*VW*VEC2(JJ)+VW*VW*VW*VEC2(JJ)) &
     &                   + VEC5(JJ)*(VEC3(II)+VW*VW*VEC3(II) &
     &                   - 3.0E0_realk*VW*VEC2(II)+VW*VW*VW*VEC2(II)))/ &
     &                     (BNDL2*BNDL2*(1.0E0_realk-VW*VW)**2)
                  CMPKK1 = (UW*VEC2(KK)-VEC1(KK))/ &
     &                     (BNDL1*BNDL2*(1.0E0_realk-UW*UW))
                  CMPKK2 = (VW*VEC2(KK)-VEC3(KK))/ &
     &                     (BNDL3*BNDL2*(1.0E0_realk-VW*VW))
!
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,2)-1)*3+JJ) &
     &                 = CMPUU
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &                 = -CMPUU + CMPUW + FACIJ*CMPKK1
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+JJ,(optinfo%INTCRD(IC,2)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &                 = -CMPUW - FACIJ*CMPKK1
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+JJ,(optinfo%INTCRD(IC,2)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,2)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,5)-1)*3+II,(optinfo%INTCRD(IC,5)-1)*3+JJ) &
     &                 = -CMPVV
                  DWILBM((optinfo%INTCRD(IC,5)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &                 = -CMPVW - FACIJ*CMPKK2
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+JJ,(optinfo%INTCRD(IC,5)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,5)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,5)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &                 = CMPVV + CMPVW + FACIJ*CMPKK2
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+JJ,(optinfo%INTCRD(IC,5)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,5)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+II,(optinfo%INTCRD(IC,3)-1)*3+JJ) &
     &                 = CMPUU - 2.0E0_realk*CMPUW - CMPWW1 + CMPWW2
                  DWILBM((optinfo%INTCRD(IC,3)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &                 = CMPUW + CMPVW + CMPWW1 - CMPWW2 &
     &                 + FACIJ*(CMPKK1-CMPKK2)
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+JJ,(optinfo%INTCRD(IC,3)-1)*3+II) &
     &             = DWILBM((optinfo%INTCRD(IC,3)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ)
                  DWILBM((optinfo%INTCRD(IC,4)-1)*3+II,(optinfo%INTCRD(IC,4)-1)*3+JJ) &
     &                 = -CMPVV - 2.0E0_realk*CMPVW + CMPWW2 - CMPWW1
 45            CONTINUE
 40         CONTINUE
!
!     If the second derivative is undefined (sin(q) ~ 0), we set the
!     matrix equal to the identity matrix.
!
         ELSE
            DO 50 I = 1, optinfo%ICartCoord
               DWILBM(I,I) = D1
 50         CONTINUE
         END IF
      END IF
      RETURN
      END

!  /* Deck gtbinv */
      SUBROUTINE LS_GTBINV(MXRCRD,UMAT,WMAT,VMAT,VAL,WILBMT,BMTRAN, &
     &     BMTINV,PJINMT,VMATRED,TMPMAT,lupri,optinfo)
use precision
use ls_util 
use files
use optimization_input 
!
!     Decomposes the transpose of Wilson's B matrix, then one can
!     construct the inverse of the rectangulat matrix.
!
!     7/5-2009 Optimization of code, mostly by replacing matrix multiplies 
!              with dgemm (and by introducing VMATRED). /Simen Reine
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXCOOR)
      Real(realk) VMAT(MXRCRD,MXRCRD), VAL(MXRCRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) BMTINV(MXRCRD,MXCOOR), PJINMT(MXRCRD,MXRCRD)
      Real(realk) VMATRED(MXRCRD,MXRCRD), TMPMAT(MXRCRD,MXRCRD)
      PARAMETER (D0 = 0.0E0_realk)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      IF (optinfo%IPrint .GE. IPRMED) call lsheader(lupri,'Output from GTBINV')
      IDIM1 = optinfo%ICartCoord
      IDIM2 = optinfo%NIntCoord
!     New, more efficient verison of DECOMP
!     VMATRED is used to store VMATTR in DECOMP temporarily
      call ls_DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT,TMPMAT, &
     &     VMATRED,VAL,lupri)
!
!     Some output
!
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'U matrix in GTBINV')
         call output(UMAT,1,IDIM1,1,IDIM2,MXCOOR,MXRCRD,1,LUPRI)
         call lsheader(lupri,'W matrix in GTBINV')
         call output(WMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
         call lsheader(lupri,'V matrix in GTBINV')
         call output(VMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
      END IF
!
!     Based on the singular value decomposition of the B matrix,
!     we can construct a non-redundant B matrix mixing the various
!     internal primitives. The eigenvectors is stored in BMTRAN.
!
!     Not optimized time-wise! /SR
!
      IF (optinfo%DelInt .AND. (optinfo%NIntCoord .LT. 0)) THEN
!
!     The redundant B matrix is copied to BMTRAN.
!
         call ls_DZERO(BMTRAN,MXRCRD*MXRCRD)
         DO 540 I = 1, IDIM2
            DO 545 J = 1, IDIM1
               BMTRAN(I,J) = WILBMT(I,J)
 545        CONTINUE
 540     CONTINUE
!
!     We remove all vectors in V corresponding to singular values.
!
!hjaaj   WMAT(IDIM2+1,IDIM2+1) = 1.1E6_realk
!hjaaj-Oct07: bugfix: this goes out of bounds!!!!!
         IDIM3 = IDIM2
         IC = 1
 550     CONTINUE
         IF (ABS(WMAT(IC,IC)) .LT. 1.0E-8_realk) THEN
            DO 560 I = IC, IDIM3 - 1
               DO 565 J = 1, IDIM2
                  VMAT(J,I) = VMAT(J,I+1)
 565           CONTINUE
               DO 567 J = 1, IDIM1
                  UMAT(J,I) = UMAT(J,I+1)
 567           CONTINUE
               WMAT(I,I) = WMAT(I+1,I+1)
 560        CONTINUE
!hjaaj-Oct07: new code follows next 5 lines
            IF (IDIM3 .EQ. IDIM2) THEN
               WMAT(IDIM3,IDIM3) = 1.1E6_realk
            ELSE
               WMAT(IDIM3,IDIM3) = WMAT(IDIM3+1,IDIM3+1)
            END IF
            IDIM3 = IDIM3 - 1
            GOTO 550
         ELSE IF (ABS(WMAT(IC,IC)) .LT. 1.0E6_realk) THEN
            IC = IC + 1
            IF (IC .LE. IDIM2) GOTO 550
         END IF
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Reduced U matrix')
            call output(UMAT,1,IDIM1,1,IDIM3,MXCOOR,MXRCRD,1,LUPRI)
            call lsheader(lupri,'Reduced W matrix')
            call output(WMAT,1,IDIM3,1,IDIM3,MXRCRD,MXRCRD,1,LUPRI)
            call lsheader(lupri,'Reduced V matrix (non-redund. eigenvectors)' )
            call output(VMAT,1,IDIM2,1,IDIM3,MXRCRD,MXRCRD,1,LUPRI)
         END IF
!
!     We calculate the weights of the different primitives
!
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Primitive weights in final space')
            DO 600 IPRIM = 1, IDIM2
               WG = 0.0E0_realk
               DO 610 I = 1, IDIM3
                  WG = WG + VMAT(IPRIM,I)*VMAT(IPRIM,I)
 610           CONTINUE
               WRITE(LUPRI,'(A,I4,F24.6)') '        ',IPRIM,WG
 600        CONTINUE
         END IF
!
!     We construct the B matrix for the delocalized
!     internal coordinates.
!
         call ls_DZERO(WILBMT,MXRCRD*MXCOOR)
         DO 620 I = 1, IDIM3
            DO 622 J = 1, IDIM1
               DO 624 K = 1, IDIM2
                  WILBMT(I,J) = WILBMT(I,J) + VMAT(K,I)*BMTRAN(K,J)
 624           CONTINUE
 622        CONTINUE
 620     CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'B matrix of active space')
            call output(WILBMT,1,IDIM3,1,IDIM1,MXRCRD,MXCOOR,1,LUPRI)
         END IF
!
!     The eigenvectors (V matrix) is stored in BMTRAN
!
         call ls_DZERO(BMTRAN,MXRCRD*MXRCRD)
         DO 630 I = 1, IDIM2
            DO 635 J = 1, IDIM3
               BMTRAN(I,J) = VMAT(I,J)
 635        CONTINUE
 630     CONTINUE
         optinfo%NIntCoord = optinfo%NIntCoord
         optinfo%NIntCoord = IDIM3
         IDIM2 = IDIM3
!
!     Then we start the decomposition again with this matrix
!     (which now should have no singular values).
!
         call ls_DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT,TMPMAT, &
     &        VMATRED,VAL,lupri)
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Non-redundant U matrix in GTBINV')
            call output(UMAT,1,IDIM1,1,IDIM2,MXCOOR,MXRCRD,1,LUPRI)
            call lsheader(lupri,'Non-redundant W matrix in GTBINV')
            call output(WMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
            call lsheader(lupri,'Non-redundant V matrix in GTBINV')
            call output(VMAT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
         END IF
      END IF
!
!     The inverse is constructed:
!
!           (B^t)^-1 = V (1/w_ii) U^t
!
!     All non-negative values of W are inverted, values close to zero
!     are set to zero.
!
      optinfo%nProjected = 0
      DO 130 I = 1, IDIM2
         IF (ABS(WMAT(I,I)) .LE. 1.0E-10_realk) THEN
            WMAT(I,I) = D0
         ELSE
            WMAT(I,I) = 1.0E0_realk/WMAT(I,I)
            optinfo%nProjected = optinfo%nProjected + 1
         END IF
 130  CONTINUE
!
!     A projection operator that will be used later,
!     is also constructed: P = V W W^-1 V^t
!
!     Set up reduced V
      NK = 0
      DO K=1,IDIM2
        IF (ABS(WMAT(K,K)) .GE. 1.0E-10_realk) THEN
          NK = NK + 1
          DO I=1,IDIM2
            VMATRED(I,NK) = VMAT(I,K)
          ENDDO
        ENDIF
      ENDDO
!     Calculate P
      call ls_DZERO(PJINMT,MXRCRD*MXRCRD)
      call dgemm('N','T',IDIM2,IDIM2,NK,1.0E0_realk,VMATRED,MXRCRD, &
     &           VMATRED,MXRCRD,0.0E0_realk,PJINMT,MXRCRD)

      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Projection matrix')
         call output(PJINMT,1,IDIM2,1,IDIM2,MXRCRD,MXRCRD,1,LUPRI)
      END IF
!
!     We construct (1/w_ii) U^t
!
      call dgemm('N','T',IDIM2,IDIM1,IDIM2,1.0E0_realk,WMAT,MXRCRD, &
     &           UMAT,MXCOOR,0.0E0_realk,BMTINV,MXRCRD)
!
!     Then V (1/w_ii) U^t
!
      call dgemm('N','N',IDIM2,IDIM1,IDIM2,1.0E0_realk,VMAT,MXRCRD, &
     &           BMTINV,MXRCRD,0.0E0_realk,WMAT,MXRCRD)

      call ls_DZERO(BMTINV,MXRCRD*MXCOOR)
      DO 180 J = 1, IDIM1
         DO 182 I = 1, IDIM2
            BMTINV(I,J) = WMAT(I,J)
 182     CONTINUE
 180  CONTINUE
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Inverse of B^t')
         call output(BMTINV,1,IDIM2,1,IDIM1,MXRCRD,MXCOOR,1,LUPRI)
      END IF
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Total number of internal coordinates: ',optinfo%NIntCoord
         WRITE(LUPRI,*) 'Number of non-redundant coordinates : ',optinfo%nProjected
         WRITE(LUPRI,*)
      END IF
!
!     Put number of coordinates to be projected away in optinfo%nProjected
!
      optinfo%nProjected = optinfo%NIntCoord - optinfo%nProjected
      RETURN
      END

!  /* Deck decomp */
      SUBROUTINE LS_DECOMP(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,WILBMT, &
     &                  WBMTTR,VMATTR,SVAL,lupri)
use precision
use ls_util 
use files
use optimization_input 
use memory_handling
!
!     Decomposes the transpose of Wilson's B matrix.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Real(realk) UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXRCRD)
      Real(realk) VMAT(MXRCRD,MXRCRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR)
      Real(realk) WBMTTR(MXCOOR,MXRCRD)
      Real(realk) VMATTR(MXRCRD,MXRCRD)
      Real(realk) SVAL(MXCOOR)
      Real(realk), pointer :: Work(:)
      Integer :: LWork
      Real(realk), PARAMETER ::D0 = 0.0E0_realk
      INFO=0
!
!     The transpose of Wilsons B matrix is decomposed to:
!
!           B^t = U W V^t  ;   U,V orthogonal, W diagonal
!
      call ls_DZERO(UMAT,MXRCRD*MXCOOR)
      call ls_DZERO(WMAT,MXRCRD*MXRCRD)
      call ls_DZERO(VMAT,MXRCRD*MXRCRD)
      call ls_DZERO(WBMTTR,MXCOOR*MXRCRD)
      call ls_DZERO(SVAL,MXCOOR)
!
!     The transpose is placed in WBMTTR
!
      DO I = 1, IDIM2
         DO J = 1, IDIM1
            WBMTTR(J,I) = WILBMT(I,J)
         ENDDO
      ENDDO
!
!     The transpose of Wilsons B matrix is decomposed to:
!
!           B^t = U W V^t  ;   U,V orthogonal, W diagonal
!
 
!     Allocate some memory for single-value decomposition
      LWork = MAX(1,3*MIN(IDim1,IDim2)+MAX(IDim1,IDim2),5*MIN(IDim1,IDim2))
!Simen Hack for ACML. The above Lwork isnot sufficient, and performance is better if LWork
!Simen is larger
      LWork = MAX(LWork,10000000)
      Call mem_alloc(Work,LWork)
      call DGESVD('A','A',IDIM1,IDIM2,WBMTTR,MXCOOR,SVAL,UMAT,MXCOOR, &
     &            VMATTR,MXRCRD,Work,LWork,INFO)
!     Deallocate memory
      Call mem_dealloc(Work)
!
      IF (INFO.NE. 0) THEN
        WRITE(*,*) 'Error in DECOMP! Info = ',INFO
        call lsquit('Error in DECOMP!',lupri)
      ENDIF
      DO I=1,IDIM2
        WMAT(I,I) = SVAL(I)
      ENDDO
      DO I=1,IDIM2
        DO J=1,IDIM2
          VMAT(J,I) = VMATTR(I,J)
        ENDDO
      ENDDO
      RETURN
      END

!  /* Deck dcmpol */
      SUBROUTINE LS_DCMPOL(MXRCRD,IDIM1,IDIM2,UMAT,WMAT,VMAT,VAL, &
     &     WILBMT,BMTINV,lupri)
use precision
use ls_util 
use optimization_input 
use files
!
!     Old version of subroutine to decompose the transpose of Wilson's B matrix.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Real(realk) UMAT(MXCOOR,MXRCRD), WMAT(MXRCRD,MXRCRD)
      Real(realk) VMAT(MXRCRD,MXRCRD), VAL(MXRCRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR)
      Real(realk) BMTINV(MXRCRD,MXCOOR)
      PARAMETER (D0 = 0.0E0_realk)
!
      call ls_DZERO(BMTINV,MXRCRD*MXCOOR)
      call ls_DZERO(UMAT,MXRCRD*MXCOOR)
      call ls_DZERO(WMAT,MXRCRD*MXRCRD)
      call ls_DZERO(VMAT,MXRCRD*MXRCRD)
      call ls_DZERO(VAL,MXRCRD)
!
!     The transpose is placed in UMAT
!
      DO 10 I = 1, IDIM2
         DO 12 J = 1, IDIM1
            UMAT(J,I) = WILBMT(I,J)
 12      CONTINUE
 10   CONTINUE
!
!     The transpose of Wilsons B matrix is decomposed to:
!
!           B^t = U W V^t  ;   U,V orthogonal, W diagonal
!
      BMTNRM = D0
      SCL    = D0
      SGNFAC = D0
      DO 20 I = 1, IDIM2
         VAL(I) = SCL*SGNFAC
         SCL    = D0
         SGNFAC = D0
         SUM    = D0
         II = I + 1
         IF (I .LE. IDIM1) THEN
            DO 25 J = I, IDIM1
               SCL = SCL + ABS(UMAT(J,I))
 25         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 30 J = I, IDIM1
                  UMAT(J,I) = UMAT(J,I)/SCL
                  SUM = SUM + UMAT(J,I)*UMAT(J,I)
 30            CONTINUE
               DIAG = UMAT(I,I)
               SGNFAC = -SIGN(SQRT(SUM),DIAG)
               FAC = SGNFAC*DIAG - SUM
               UMAT(I,I) = DIAG - SGNFAC
               DO 35 J = II, IDIM2
                  SUM = D0
                  DO 37 K = I, IDIM1
                     SUM = SUM + UMAT(K,I)*UMAT(K,J)
 37               CONTINUE
                  FAC2 = SUM/FAC
                  DO 39 K = I, IDIM1
                     UMAT(K,J) = UMAT(K,J)+FAC2*UMAT(K,I)
 39               CONTINUE
 35            CONTINUE
               DO 40 J = I, IDIM1
                  UMAT(J,I) = SCL*UMAT(J,I)
 40            CONTINUE
            END IF
         END IF
         WMAT(I,I) = SCL*SGNFAC
         SUM    = D0
         SGNFAC = D0
         SCL    = D0
         IF ((I .LE. IDIM1) .AND. (I .NE. IDIM2)) THEN
            DO 45 J = II, IDIM2
               SCL = SCL + ABS(UMAT(I,J))
 45         CONTINUE
            IF (SCL .NE. D0) THEN
               DO 50 J = II, IDIM2
                  UMAT(I,J) = UMAT(I,J)/SCL
                  SUM = SUM + UMAT(I,J)*UMAT(I,J)
 50            CONTINUE
               FAC = UMAT(I,II)
               SGNFAC = -SIGN(SQRT(SUM),FAC)
               FAC2 = SGNFAC*FAC - SUM
               UMAT(I,II) = FAC - SGNFAC
               DO 52 J = II, IDIM2
                  VAL(J) = UMAT(I,J)/FAC2
 52            CONTINUE
               DO 54 J = II, IDIM1
                  SUM = D0
                  DO 56 K = I+1, IDIM2
                     SUM = SUM + UMAT(J,K)*UMAT(I,K)
 56               CONTINUE
                  DO 58 K = II, IDIM2
                     UMAT(J,K) = UMAT(J,K) + SUM*VAL(K)
 58               CONTINUE
 54            CONTINUE
               DO 60 J = II, IDIM2
                  UMAT(I,J) = SCL*UMAT(I,J)
 60            CONTINUE
            END IF
         END IF
         BMTNRM = MAX(BMTNRM,(ABS(WMAT(I,I)) + ABS(VAL(I))))
 20   CONTINUE
      DO 70 I = IDIM2, 1, -1
         IF (I .LT. IDIM2) THEN
            IF (SGNFAC .NE. D0) THEN
               DO 72 J = II, IDIM2
                  VMAT(J,I) = UMAT(I,J)/(SGNFAC*UMAT(I,II))
 72            CONTINUE
               DO 74 J = II, IDIM2
                  SUM = D0
                  DO 76 K = II, IDIM2
                     SUM = SUM + UMAT(I,K)*VMAT(K,J)
 76               CONTINUE
                  DO 78 K = II, IDIM2
                     VMAT(K,J) = VMAT(K,J) + SUM*VMAT(K,I)
 78               CONTINUE
 74            CONTINUE
            END IF
            DO 79 J = II, IDIM2
               VMAT(I,J) = D0
               VMAT(J,I) = D0
 79         CONTINUE
         END IF
         VMAT(I,I) = 1.0E0_realk
         SGNFAC = VAL(I)
         II = I
 70   CONTINUE
      DO 80 I = MIN(IDIM1,IDIM2), 1, -1
         II = I + 1
         SGNFAC = WMAT(I,I)
         DO 82 J = II, IDIM2
            UMAT(I,J) = D0
 82      CONTINUE
         IF (SGNFAC .NE. D0) THEN
            SGNFAC = 1.0E0_realk/SGNFAC
            DO 84 J = II, IDIM2
               SUM = D0
               DO 86 K = II, IDIM1
                  SUM = SUM + UMAT(K,I)*UMAT(K,J)
 86            CONTINUE
               FAC = SGNFAC*(SUM/UMAT(I,I))
               DO 88 K = I, IDIM1
                  UMAT(K,J) = UMAT(K,J) + FAC*UMAT(K,I)
 88            CONTINUE
 84         CONTINUE
            DO 90 J = I, IDIM1
               UMAT(J,I) = SGNFAC*UMAT(J,I)
 90         CONTINUE
         ELSE
            DO 92 J = I, IDIM1
               UMAT(J,I) = D0
 92         CONTINUE
         END IF
         UMAT(I,I) = UMAT(I,I) + 1.0E0_realk
 80   CONTINUE
      DO 100 J = IDIM2, 1, -1
         DO 105 ITRS = 1, 50
            DO 110 K = J, 1, -1
               KK = K - 1
               IF ((ABS(VAL(K))+BMTNRM) .EQ. BMTNRM) GOTO 902
               IF ((ABS(WMAT(KK,KK))+BMTNRM) .EQ. BMTNRM) GOTO 901
 110        CONTINUE
 901        CONTINUE
            FC   = D0
            SUM = 1.0E0_realk
            DO 112  I = K, J
               FAC = SUM*VAL(I)
               VAL(I) = FC*VAL(I)
               IF ((ABS(FAC)+BMTNRM) .EQ. BMTNRM) GOTO 902
               SGNFAC = WMAT(I,I)
               FAC2 = SQRT(FAC*FAC + SGNFAC*SGNFAC)
               WMAT(I,I) = FAC2
               FAC2 = 1.0E0_realk/FAC2
               FC  = SGNFAC*FAC2
               SUM =   -FAC*FAC2
               DO 114 L = 1, IDIM1
                  FY = UMAT(L,KK)
                  FZ = UMAT(L,I)
                  UMAT(L,KK) = FY*FC + FZ*SUM
                  UMAT(L,I)  = FZ*FC - FY*SUM
 114           CONTINUE
 112        CONTINUE
 902        CONTINUE
            FZ = WMAT(J,J)
            IF (K .EQ. J) THEN
               IF (FZ .LT. D0) THEN
                  WMAT(J,J) = -FZ
                  DO 116 L = 1, IDIM2
                     VMAT(L,J) = -VMAT(L,J)
 116              CONTINUE
               END IF
               GOTO 903
            END IF
            IF (ITRS .EQ. 50) call lsquit('No convergence in GTBINV.',lupri)
            FX = WMAT(K,K)
            JJ = J - 1
            FY = WMAT(JJ,JJ)
            SGNFAC = VAL(JJ)
            FAC2 = VAL(J)
            FAC = ((FY-FZ)*(FY+FZ) + (SGNFAC-FAC2)*(SGNFAC+FAC2))/ &
     &           (2.0E0_realk*FAC2*FY)
            SGNFAC = SQRT(FAC*FAC + 1.0E0_realk)
            FAC = ((FX-FZ)*(FX+FZ)+FAC2*((FY/(FAC+SIGN(SGNFAC,FAC))) &
     &           -FAC2))/FX
            FC  = 1.0E0_realk
            SUM = 1.0E0_realk
            DO 120 L = K, JJ
               I = L + 1
               SGNFAC = VAL(I)
               FY = WMAT(I,I)
               FAC2   = SGNFAC*SUM
               SGNFAC = SGNFAC*FC
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               VAL(L) = FZ
               FC  = FAC/FZ
               SUM = FAC2/FZ
               FAC = (FX*FC) + (SGNFAC*SUM)
               SGNFAC = -(FX*SUM) + (SGNFAC*FC)
               FAC2 = FY*SUM
               FY = FY*FC
               DO 122 M = 1, IDIM2
                  FX = VMAT(M,L)
                  FZ = VMAT(M,I)
                  VMAT(M,L) = (FX*FC) + (FZ*SUM)
                  VMAT(M,I) =-(FX*SUM) + (FZ*FC)
 122           CONTINUE
               FZ = SQRT(FAC*FAC + FAC2*FAC2)
               WMAT(L,L) = FZ
               IF (FZ .NE. D0) THEN
                  FZ = 1.0E0_realk/FZ
                  FC = FAC*FZ
                  SUM = FAC2*FZ
               END IF
               FAC = (FC*SGNFAC) + (SUM*FY)
               FX = -(SUM*SGNFAC) + (FC*FY)
               DO 124 M = 1, IDIM1
                  FY = UMAT(M,L)
                  FZ = UMAT(M,I)
                  UMAT(M,L) = (FY*FC) + (FZ*SUM)
                  UMAT(M,I) = -(FY*SUM) + (FZ*FC)
 124           CONTINUE
 120        CONTINUE
            VAL(K) = D0
            VAL(J) = FAC
            WMAT(J,J) = FX
 105     CONTINUE
 903     CONTINUE
 100  CONTINUE
      RETURN
      END

!  /* Deck gq2gx */
      SUBROUTINE LS_GQ2GX(MXRCRD,GRADQ,GRADX,WILBMT,optinfo)
use precision
use ls_util 
use optimization_input 
use files
!
!     Transforms the gradient in redundant internal coordinates to
!     Cartesian coordinates:
!                              g_x = B^t g_q
!
Implicit Real(realk) (A-H,O-Z)
      Type(opt_setting) :: optinfo
      Real(realk) GRADQ(MXRCRD), GRADX(MXCOOR), WILBMT(MXRCRD,MXCOOR)
      call ls_DZERO(GRADX,MXCOOR)
      DO 10 I = 1, optinfo%ICartCoord
         DO 12 J = 1, optinfo%NIntCoord
            GRADX(I) = GRADX(I) + WILBMT(J,I)*GRADQ(J)
 12      CONTINUE
 10   CONTINUE
      RETURN
      END

!  /* Deck gx2gq */
      SUBROUTINE LS_GX2GQ(MXRCRD,GRADX,GRADQ,BMTINV,optinfo)
use ls_util 
use precision
use files
use optimization_input
!
!     Transforms the gradient in Cartesian coordinates to redundant
!     internal coordinates:
!                              g_q = (B^t)^-1 g_x
!
Implicit Real(realk) (A-H,O-Z)
      Type(opt_setting) :: optinfo
      Real(realk) GRADX(MXCOOR), GRADQ(MXRCRD), BMTINV(MXRCRD,MXCOOR)
      call ls_DZERO(GRADQ,MXRCRD)
      DO 10 I = 1, optinfo%NIntCoord
         DO 12 J = 1, optinfo%ICartCoord
            GRADQ(I) = GRADQ(I) + BMTINV(I,J)*GRADX(J)
 12      CONTINUE
 10   CONTINUE 
      RETURN
      END

!  /* Deck hq2hx */
      SUBROUTINE LS_HQ2HX(Molecule,MXRCRD,MX2CRD,ATMARR,TMPMAT, &
                TMPMT2,HESSQ,GRADQ,HESSX,WILBMT,BMTRAN,optinfo)
use ls_util 
use precision
use optimization_input 
use files
use molecule_type
use molecule_typetype
!
!     Transforms the Hessian in redundant internal coordinates to
!     Cartesian coordinates:
!                              H_x = B^t H_q B + B'^t g_q
!
Implicit Real(realk) (A-H,O-Z)
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8), TMPMAT(MX2CRD,MX2CRD)
      Real(realk) TMPMT2(MX2CRD,MX2CRD), HESSQ(MXRCRD,MXRCRD)
      Real(realk) GRADQ(MXRCRD), HESSX(optinfo%ICartCoord,optinfo%ICartCoord)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
      DO 10 IC = 1, optinfo%NIntCoord
         call ls_GTDWL0(Molecule,MXRCRD,IC,TMPMT2,ATMARR,WILBMT,BMTRAN,& 
         & lupri,optinfo)
         DO 12 I = 1, optinfo%ICartCoord
            DO 14 J = 1, optinfo%ICartCoord
               TMPMAT(I,J) = TMPMAT(I,J) + TMPMT2(I,J)*GRADQ(IC)
 14         CONTINUE
 12      CONTINUE
 10   CONTINUE
      call ls_DZERO(TMPMT2,MX2CRD*MX2CRD)
      DO 20 I = 1, optinfo%NIntCoord
         DO 22 J = 1, optinfo%ICartCoord
            DO 24 K = 1, optinfo%NIntCoord
               TMPMT2(I,J) = TMPMT2(I,J) + HESSQ(I,K)*WILBMT(K,J)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      call ls_DZERO(HESSX,optinfo%ICartCoord*optinfo%ICartCoord)
      DO 30 I = 1, optinfo%ICartCoord
         DO 32 J = 1, optinfo%ICartCoord
            DO 34 K = 1, optinfo%NIntCoord
               HESSX(I,J) = HESSX(I,J) + WILBMT(K,I)*TMPMT2(K,J)
 34         CONTINUE
            HESSX(I,J) = HESSX(I,J) + TMPMAT(I,J)
 32      CONTINUE
 30   CONTINUE
!
!     We make sure the Cartesian Hessian is symmetric.
!
      DO 40 I = 1, optinfo%ICartCoord
         DO 42 J = 1, I
            HESSX(J,I) = HESSX(I,J)
 42      CONTINUE
 40   CONTINUE
      RETURN
      END

!  /* Deck hx2hq */
      SUBROUTINE LS_HX2HQ(Molecule,MXRCRD,MX2CRD,ATMARR,TMPMAT, &
      TMPMT2,TMPMT3,HESSX,GRADQ,HESSQ,WILBMT,BMTINV,BMTRAN,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use molecule_type
use molecule_typetype
!
!     Transforms the Hessian in Cartesian coordinates to redundant
!     internal coordinates:
!                         H_q = (B^t)^-1 (H_x - B'^t g_q) B^-1
!
Implicit Real(realk) (A-H,O-Z)
      Type(opt_setting) :: optinfo
      Type(MOLECULEINFO) :: Molecule
      Real(realk) ATMARR(MXCENT,8), TMPMAT(MXCOOR,MXCOOR)
      Real(realk) TMPMT2(MX2CRD,MX2CRD), TMPMT3(MX2CRD,MX2CRD)
      Real(realk) HESSX(optinfo%ICartCoord,optinfo%ICartCoord), GRADQ(MXRCRD)
      Real(realk) HESSQ(MXRCRD,MXRCRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTINV(MXRCRD,MXCOOR)
      Real(realk) BMTRAN(MXRCRD,MXRCRD)
      LOGICAL INCLGT
      call ls_DZERO(HESSQ,MXRCRD*MXRCRD)
      call ls_DZERO(TMPMAT,MXCOOR*MXCOOR)
      call ls_DZERO(TMPMT2,MX2CRD*MX2CRD)
      call ls_DZERO(TMPMT3,MX2CRD*MX2CRD)
!
!     INCLGT determines wether the gradient term in the transformation
!     is included. There seems to be no reason NOT to include it!!
!     (in some articles this term has been accused of ruining the
!     eigenvalue structure of the Hessian).
!
      INCLGT = .TRUE.
!     
      IF (INCLGT) THEN
         DO 10 IC = 1, optinfo%NIntCoord
            call ls_GTDWL0(Molecule,MXRCRD,IC,TMPMAT,ATMARR,WILBMT,BMTRAN, &
     &           lupri,optinfo)
            DO 12 I = 1, optinfo%ICartCoord
               DO 14 J = 1, optinfo%ICartCoord
                  TMPMT2(I,J) = TMPMT2(I,J) + TMPMAT(I,J)*GRADQ(IC)
 14            CONTINUE
 12         CONTINUE
 10      CONTINUE
      END IF
      DO 16 J = 1, optinfo%ICartCoord
         DO 18 I = 1, optinfo%ICartCoord
            TMPMT2(I,J) = HESSX(I,J) - TMPMT2(I,J)
 18      CONTINUE
 16   CONTINUE
      DO 20 I = 1, optinfo%ICartCoord
         DO 22 J = 1, optinfo%NIntCoord
            DO 24 K = 1, optinfo%ICartCoord
               TMPMT3(I,J) = TMPMT3(I,J) + TMPMT2(I,K)*BMTINV(J,K)
 24         CONTINUE
 22      CONTINUE
 20   CONTINUE
      DO 30 I = 1, optinfo%NIntCoord
         DO 32 J = 1, optinfo%NIntCoord
            DO 34 K = 1, optinfo%ICartCoord
               HESSQ(I,J) = HESSQ(I,J) + BMTINV(I,K)*TMPMT3(K,J)
 34         CONTINUE
 32      CONTINUE
 30   CONTINUE
      RETURN
      END

!  /* Deck cghint */
      SUBROUTINE LS_CGHINT(Molecule,MXRCRD,MX2CRD,HESTMP,TMPMAT,TMPMT2, &
     &     TMPMT3,TMPMT4,WILBMT,BMTINV,BMTRAN,HESINT,lupri,optinfo)
use ls_util 
use files
use optimization_input 
use precision
use molecule_type
use molecule_typetype
!
!     Transforms Cartesian gradient and Hessian to redundant internal
!     coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) TMPMAT(MX2CRD,MX2CRD), TMPMT2(MX2CRD,MX2CRD)
      Real(realk) TMPMT3(MX2CRD,MX2CRD), TMPMT4(MX2CRD,MX2CRD)
      Real(realk) WILBMT(MXRCRD,MXCOOR), BMTINV(MXRCRD,MXCOOR)
      Real(realk) BMTRAN(MXRCRD,MXRCRD)
      Real(realk) HESINT(MXRCRD,MXRCRD), HESTMP(MXRCRD,MXRCRD)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Output from CGHINT')
      END IF
      Write(*,*)'CGHINT,MXRCRD=',MXRCRD
      call ls_GX2GQ(optinfo%NIntCoord,optinfo%GradMol,optinfo%GRDINT,BMTINV,optinfo)
      GRADNM = SQRT(DDOT(optinfo%ICartCoord,optinfo%GradMol,1,optinfo%GradMol,1))
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Cartesian gradient')
         call output(optinfo%GradMol,1,1,1,optinfo%ICartCoord,1,MXCOOR,1,LUPRI)
      END IF
      call ls_HX2HQ(Molecule,MXRCRD,MX2CRD,TMPMAT,TMPMT2,TMPMT3,TMPMT4,&
       optinfo%HessMol,optinfo%GRDINT,HESINT,WILBMT,BMTINV,BMTRAN,optinfo)
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Cartesian Hessian')
         call output(optinfo%HessMol,1,optinfo%ICartCoord,1, &
  &      optinfo%ICartCoord,MXCOOR,MXCOOR,1,LUPRI)
      END IF
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Internal gradient')
         call output(optinfo%GRDINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Internal Hessian')
         Write(*,*)'HESINT',HESINT
         call output(HESINT,1,optinfo%NIntCoord,1,optinfo%NIntCoord,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck prjint */
      SUBROUTINE LS_PRJINT(MXRCRD,NDIM,PJINMT,CONMAT,HESINT,TMPMT1,TMPMT2, &
     &     TMPMT3,TMPMT4,lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use precision
use memory_handling
!
!     Because of possible redundancies, we have to project
!     both the gradient and Hessian. If a constrained optimization
!     has been requested, we have to construct a matrix of constraints
!     and modify the projection matrix.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) PJINMT(MXRCRD,MXRCRD), CONMAT(MXRCRD,MXRCRD)
      Real(realk) HESINT(MXRCRD,MXRCRD), TMPMT1(NDIM,NDIM)
      Real(realk) TMPMT2(NDIM,NDIM), TMPMT3(NDIM,NDIM)
      Real(realk) TMPMT4(NDIM,NDIM)
      Real(realk), pointer :: Work(:)
      PARAMETER (IPRMIN = 0, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      NDIM2 = NDIM*NDIM
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Output from PRJINT')
      END IF
!

!
!     The matrix of constraints is constructed if necessary.
!
      IF (optinfo%ConOpt) THEN
         call ls_DZERO(CONMAT,MXRCRD*MXRCRD)
         call ls_DZERO(TMPMT1,NDIM2)
         call ls_DZERO(TMPMT2,NDIM2)
         call ls_DZERO(TMPMT3,NDIM2)
         call ls_DZERO(TMPMT4,NDIM2)
         DO 10 I = 1, NDIM
            IF (optinfo%IConstr(I) .EQ. 1) THEN
               CONMAT(I,I) = 1.0E0_realk
               TMPMT1(I,I) = PJINMT(I,I)
               DO 15 J = I+1, NDIM
!
!     At the same time we construct the product CP'C,
!     Where C is the matrix of constraints and P' the original
!     projection matrix.
!
                  TMPMT1(I,I) = PJINMT(I,I)
                  IF (optinfo%IConstr(J) .EQ. 1) THEN
                     TMPMT1(I,J) = PJINMT(I,J)
                     TMPMT1(J,I) = PJINMT(I,J)
                  END IF
 15            CONTINUE
            END IF
 10      CONTINUE
!
!     We find the inverse of CP'C
!
   
!     Allocate Work first
         Call mem_alloc(Work,NDIM)
!     Call LAPACK routine to get inverse     
         call DGEINV(NDIM,TMPMT1,TMPMT2,Work,INFO)
!     Deallocate Work
         Call mem_dealloc(Work)
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Matrix of constraints')
            call output(CONMAT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
            call lsheader(lupri,'CP''C')
            call output(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            call lsheader(lupri,'Inverse of CP''C')
            call output(TMPMT2,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
!
!     We find CP' and P'C
!
         call ls_DZERO(TMPMT1,NDIM2)
         call ls_DZERO(TMPMT3,NDIM2)
         DO 20 I = 1, NDIM
            IF (optinfo%IConstr(I) .EQ. 1) THEN
               DO 22 J = 1, NDIM
                  DO 24 K = 1, NDIM
                     TMPMT1(I,J) = TMPMT1(I,J) + CONMAT(I,K)*PJINMT(K,J)
                     TMPMT3(J,I) = TMPMT3(J,I) + PJINMT(J,K)*CONMAT(K,I)
 24               CONTINUE
 22            CONTINUE
            END IF
 20      CONTINUE
!
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'CP''')
            call output(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            call lsheader(lupri,'P''C')
            call output(TMPMT3,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
!
!     We calculate (CP'C)^-1 CP'
!
         DO 30 I = 1, NDIM
            IF (optinfo%IConstr(I) .EQ. 1) THEN
               DO 32 J = 1, NDIM
                  DO 34 K = 1, NDIM
                     TMPMT4(I,J) = TMPMT4(I,J) + TMPMT2(I,K)*TMPMT1(K,J)
 34               CONTINUE
 32            CONTINUE
            END IF
 30      CONTINUE
!
!     then we find P'C (CP'C)^-1 CP'
!
         call ls_DZERO(TMPMT1,NDIM2)
         DO 40 I = 1, NDIM
            IF (optinfo%IConstr(I) .EQ. 1) THEN
               DO 42 J = 1, NDIM
                  DO 44 K = 1, NDIM
                     TMPMT1(I,J) = TMPMT1(I,J) + TMPMT3(I,K)*TMPMT4(K,J)
 44               CONTINUE
 42            CONTINUE
            END IF
 40      CONTINUE
!
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'TMPMT4')
            call output(TMPMT4,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            call lsheader(lupri,'TMPMT1')
            call output(TMPMT1,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF
!
!     The projection matrix is then modified to include the constraints.
!
         DO 50 I = 1, NDIM
            DO 52 J = 1, NDIM
               PJINMT(I,J) = PJINMT(I,J) - TMPMT1(I,J)
 52         CONTINUE
 50      CONTINUE
      END IF
!
      IF (optinfo%IPrint .GE. IPRMED) THEN
         IF (optinfo%ConOpt) THEN
            call lsheader(lupri,'Projection matrix w/constraints')
         ELSE
            call lsheader(lupri,'Projection matrix')
         END IF
         call output(PJINMT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      IF (optinfo%IPrint .GE. IPRMAX) THEN
         call lsheader(lupri,'Unprojected gradient')
         call output(optinfo%GRDINT,1,1,1,NDIM,1,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Unprojected Hessian')
         call output(HESINT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
!
!     First we project the gradient
!
      call dgemm('N','N',NDIM,1,NDIM,1.0E0_realk,PJINMT,MXRCRD, &
     &           optinfo%GRDINT,NDIM,0.0E0_realk,TMPMT1,NDIM)
      call DCOPY(NDIM,TMPMT1,1,optinfo%GRDINT,1)

      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Projected gradient')
         call output(optinfo%GRDINT,1,1,1,NDIM,1,MXRCRD,1,LUPRI)
      END IF
!
!     Then the Hessian
!
      call dgemm('N','N',NDIM,NDIM,NDIM,1.0E0_realk,PJINMT,MXRCRD, &
     &           HESINT,MXRCRD,0.0E0_realk,TMPMT1,NDIM)
      call dgemm('N','N',NDIM,NDIM,NDIM,1.0E0_realk,TMPMT1,NDIM, &
     &           PJINMT,MXRCRD,0.0E0_realk,HESINT,MXRCRD)
!
!     The projected Hessian is "stabilized".
!     1.0E4_realk is an arbitrary (high) value.
!     This section has been commented out because this
!     stabilization seems unnecessary, it just messes up some
!     of the eigenvalues.
!
!      DO 90 I = 1, NDIM
!         DO 95 J = 1, NDIM
!            HESINT(I,J) = HESINT(I,J) - 1.0E4_realk*PJINMT(I,J)
! 95      CONTINUE
!         HESINT(I,I) = HESINT(I,I) + 1.0E4_realk
! 90   CONTINUE
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Projected Hessian')
         call output(HESINT,1,NDIM,1,NDIM,MXRCRD,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck diaint */
      SUBROUTINE LS_DIAINT(MXRCRD,MX2CRD,NCRDHS,EVEC,EVCTMP,HESPCK, &
     &     TMPHES,TMPMAT,THRIND,HESINT,DIDENT,lupri,optinfo)
use ls_util 
use optimization_input 
use files
use precision
use memory_handling
!
!     The Hessian in redundant internal coordinates is diagonalized.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri, LWork
      Type(opt_setting) :: optinfo
      Real(realk) EVEC(MX2CRD,MX2CRD), EVCTMP(NCRDHS,NCRDHS)
      Real(realk) HESPCK(NCRDHS*NCRDHS), TMPHES(NCRDHS,NCRDHS)
      Real(realk) TMPMAT(MX2CRD*MX2CRD), HESINT(MXRCRD,MXRCRD)
      Real(realk) DIDENT(MXRCRD,MXRCRD)
      Real(realk) :: Work(3*NCRDHS)
!      Real(realk), allocatable :: Work(:)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (D0 = 0.0E0_realk)
      info=0
!#ifdef OLD_DIAINT
!      DO 10 J = 1, optinfo%NIntCoord
!         DO 12 I = 1, optinfo%NIntCoord
!            TMPHES(I,J) = HESINT(I,J)
! 12      CONTINUE
! 10   CONTINUE
!      IF (optinfo%RatFun .AND. (.NOT. optinfo%Saddle)) THEN
!         DO 15 I = 1, optinfo%NIntCoord
!            TMPHES(I,NCRDHS) = optinfo%GRDINT(I)
!            TMPHES(NCRDHS,I) = optinfo%GRDINT(I)
! 15      CONTINUE
!         TMPHES(NCRDHS,NCRDHS) = D0
!      END IF
!      call ls_DZERO(HESPCK,NCRDHS*NCRDHS)
!      call ls_DSITSP(NCRDHS,TMPHES,HESPCK)
!
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Output from DIAINT')
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Internal gradient')
            call output(optinfo%GRDINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
               call lsheader(lupri,'Internal Hessian')
               call output(HESINT,1,optinfo%NIntCoord,1,optinfo%NIntCoord, &
     &              MXRCRD,MXRCRD,1,LUPRI)
            IF (optinfo%RatFun .AND. (.NOT. optinfo%Saddle)) THEN
               call lsheader(lupri,'Augmented internal Hessian')
               call output(TMPHES,1,NCRDHS,1,NCRDHS, &
     &              NCRDHS,NCRDHS,1,LUPRI)
            END IF
            !call lsheader(lupri,'Packed Hessian')
            !call output(HESPCK,1,1,1,(NCRDHS*(NCRDHS+1))/2, &
            !     1,NCRDHS*NCRDHS,1,LUPRI)
         END IF
      END IF
!
!     call ls_DUNIT(EVCTMP,NCRDHS)
!     call ls_JACO(HESPCK,EVCTMP,NCRDHS,NCRDHS,NCRDHS, &
!    &     TMPMAT(1),TMPMAT(MX2CRD+1))
!      optinfo%IndHes = 0
!      DO 20 J = 1, NCRDHS
!         optinfo%EVAL(J) = HESPCK(J*(J+1)/2)
!         optinfo%GRDDIA(J) = DDOT(NCRDHS,optinfo%GRDINT,1,EVCTMP(1,J),1)
!         DO 22 I = 1, NCRDHS
!            EVEC(I,J) = EVCTMP(I,J)
! 22      CONTINUE
!         IF (optinfo%EVAL(J) .LT. -THRIND) optinfo%IndHes = optinfo%IndHes + 1
! 20   CONTINUE
!#else

!
!     Code for circumventing the time-demanding call ls_to JACO. /SR 4/5-2009
!
      DO J = 1, optinfo%NIntCoord
         DO I = 1, optinfo%NIntCoord
           TMPHES(I,J) = HESINT(I,J)
         ENDDO
      ENDDO
      IF (optinfo%RatFun .AND. (.NOT. optinfo%Saddle)) THEN
         DO I = 1, optinfo%NIntCoord
            TMPHES(I,NCRDHS) = optinfo%GRDINT(I)
            TMPHES(NCRDHS,I) = optinfo%GRDINT(I)
         ENDDO
         TMPHES(NCRDHS,NCRDHS) = D0
      END IF

      call ls_DUNIT(DIdent,NCRDHS)
!
      LWORK = 3*NCRDHS
      Call ls_dzero(Work,LWork)
      Call ls_dzero(optinfo%EVAL,MXCOOR)
      Call dsyev('V','U',NCRDHS,TMPHES,NCRDHS,optinfo%Eval,Work,LWork,Info)
!
      optinfo%IndHes = 0
      DO J = 1, NCRDHS
         optinfo%GRDDIA(J) = DDOT(NCRDHS,optinfo%GRDINT,1,TMPHES(:,J),1)
         DO I = 1, NCRDHS
            EVEC(I,J) = TMPHES(I,J)
         ENDDO
         IF (optinfo%Eval(J) .LT. -THRIND) optinfo%IndHes = optinfo%IndHes + 1
      ENDDO
!#endif
      IF (optinfo%IPrint .GE. IPRMED) THEN
         call lsheader(lupri,'Eigenvalues in DIAINT')
         call output(optinfo%Eval,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Eigenvectors in DIAINT')
         call output(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
         call lsheader(lupri,'Gradient (diagonal rep.) in DIAINT')
         call output(optinfo%GRDDIA,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
      END IF
!

      DO 25 I = 1, NCRDHS
         IF (ABS(optinfo%EVAL(I)) .LT. 1.0E-6_realk) optinfo%EVAL(I) = optinfo%EVAL(I) +  1.0E5_realk
 25   CONTINUE
      DO 30 I = 1, NCRDHS
         JMIN = I
         EMIN = optinfo%EVAL(I)
         DO 35 J = (I + 1), NCRDHS
            IF (optinfo%EVAL(J) .LT. EMIN) THEN
               EMIN = optinfo%EVAL(J)
               JMIN = J
            END IF
 35      CONTINUE
         IF (JMIN .NE. I) THEN
            call dswap(1,optinfo%EVAL(I),1,optinfo%EVAL(JMIN),1)
            call dswap(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
            call dswap(1,optinfo%GRDDIA(I),1,optinfo%GRDDIA(JMIN),1)
         END IF
 30   CONTINUE
      IF (optinfo%RatFun .AND. (.NOT. optinfo%Saddle) .AND. (optinfo%EVAL(1) .LT. -THRIND) &
     &     .AND. optinfo%IndHes .GT. 0) optinfo%IndHes = optinfo%IndHes - 1
      DO 40 I = 1, NCRDHS
         IF (ABS((ABS(optinfo%EVAL(I)) - 1.0E5_realk)) .LT. 1.0E-3_realk) &
     &        optinfo%EVAL(I) = optinfo%EVAL(I) -  1.0E5_realk
 40   CONTINUE
      IF (optinfo%IPrint .GE. IPRMIN) THEN
         call lsheader(lupri,'Sorted eigenvalues in WLKEI1')
         call output(optinfo%EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
         call lsheader(lupri,'Sorted eigenvectors in WLKEI1')
         call output(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
         call lsheader(lupri,'Gradient (sorted) in WLKEI1')
         call output(optinfo%GRDDIA,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
      END IF
      RETURN
      END

!  /* Deck fnstin */
      SUBROUTINE LS_FNSTIN(Molecule,MXRCRD,MX2CRD,NCRDHS,HESINT,EVEC, &
     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,CSTEP,WILBMT,BMTRAN, &
     &     BMTINV,GRDARR,STPARR,ACTIVE,EMOD,VECMOD,STPLIN,lupri,optinfo)
use ls_util 
use optimization_input 
use files
use Fundamental
use molecule_type
use molecule_typetype
use memory_handling
use precision
!
!     We determine the step in redundant internal coordinates,
!     then we find the corresponding step vector in
!     Cartesian coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) HESINT(MXRCRD,MXRCRD), EVEC(MX2CRD,MX2CRD)
      Real(realk) TMPMAT(MX2CRD*MX2CRD), TMPMT2(MX2CRD*MX2CRD)
      Real(realk) TMPMT3(MX2CRD*MX2CRD), TMPMT4(MX2CRD,MX2CRD)
      Real(realk) TMPMT5(MX2CRD,MX2CRD), WILBMT(MX2CRD,MX2CRD)
      Real(realk) CSTEP(MXCOOR), BMTRAN(MXRCRD,MXRCRD)
      Real(realk) BMTINV(MXRCRD,MXCOOR)
      Real(realk) GRDARR(MXRCRD,25), STPARR(MXRCRD,25)
      Real(realk) STPLIN(MXRCRD), VECMOD(MXCOOR)
      Real(realk) CMPLIM
      LOGICAL   INSIDE, ACTIVE, DOSCAL
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, DP5 = 0.5E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!
      Real(realk), pointer :: VALORG(:),VALOLD(:),VALNEW(:),DIFF(:)
      Real(realk), pointer :: TMPVEC(:),CSTRA(:,:),&
      & SCTRA(:,:),CRDORG(:,:),CRDOLD(:,:),CRDNEW(:,:)
      Real(realk), external :: dv3dot_ls     
!
      NVEC = optinfo%NIntCoord - optinfo%nProjected
      IF (optinfo%LnSearch .AND. (.NOT. optinfo%RatFun) .AND. (optinfo%ItrNmr .GT. 0)) &
     &     call ls_LINSRC(optinfo%NIntCoord,MXRCRD,optinfo%GRDINT,GRDARR,STPLIN, &
     &     STPARR,TMPMAT,TMPMT2,ACTIVE,EMOD,lupri,optinfo)
      IF (ACTIVE) THEN
         DO 5 J = 1, optinfo%NIntCoord
            DO 7 I = 1, optinfo%KEPTIT
               STPARR(J,I) = STPARR(J,I) - STPLIN(J)
 7          CONTINUE
            IF (.NOT. optinfo%RatFun) &
     &           optinfo%GRDDIA(J) = DDOT(optinfo%NIntCoord,optinfo%GRDINT,1,EVEC(1,J),1)
 5       CONTINUE
         IF (.NOT. optinfo%RatFun) THEN
            IF (optinfo%IPrint .GT. 5) THEN
               call lsheader(lupri,'Diagonal interpolated gradient')
               call output(optinfo%GRDDIA,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
            END IF
         END IF
      END IF
!
      IF (optinfo%IPrint .GE. IPRMED) call lsheader(lupri,'Output from FNSTIN')
!
!     First comes the trust region method
!
      IF (optinfo%TrustRg .OR. (optinfo%GDIIS .AND. (optinfo%KEPTIT .LT. 3))) THEN
!     For saddle point optimizations, we construct the image function.
!
         IF (optinfo%Saddle) THEN
!            IMODE = optinfo%NSPMod
!
!     We can follow a specific eigenvector if needed...
!
!            IF (optinfo%NSPMod .GT. 0) THEN
!               call ls_FNDMOD(.TRUE.,MXRCRD,EVEC,WILBMT,VECMOD, &
!     &              TMPMT2,TMPMT3,IMODE,lupri,optinfo)
!            ELSE
!
!     ... otherwise we have to find the first mode with a non-zero
!     gradient-element.
!
               IMODE = 1
               Do j = 1, NVEC
                  If (ABS(optinfo%GRDDIA(IMODE)) .LE. 1.0E-10_realk) then
                     IMODE = IMODE + 1
                  Endif
               Enddo
               If (IMODE .GT. NVEC) IMODE = 1
!            END IF
            IF (optinfo%IPrint .GE. IPRMED) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,*) 'Making image function by changing  &
     &              the sign of mode ',IMODE
               WRITE(LUPRI,*)
            END IF
            call MAKIMG(optinfo%NProjected,optinfo%NIntCoord,&
               & optinfo%EVAL,optinfo%GRDDIA,optinfo%STPDIA,IMODE,.FALSE.)
         END IF
!

!
!     First we find the internal step in diagonal representation.
!
         call ls_DZERO(optinfo%STPDIA,MXRCRD)
         DO 10 I = 1, optinfo%NIntCoord
!
!     The eigenvalue threshold was changed from 1.0E-8_realk due to
!     problems with optimization in delocalized internals.
!
            IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-6_realk) THEN
               optinfo%STPDIA(I) = D0
            ELSE
               optinfo%STPDIA(I) = -optinfo%GRDDIA(I)/optinfo%EVAL(I)
            END IF
 10      CONTINUE
         optinfo%stepNorm = SQRT(DDOT(optinfo%NIntCoord,optinfo%STPDIA,1,optinfo%STPDIA,1))
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Newton step')
            call output(optinfo%STPDIA,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         END IF
!
!     If Newton step is larger that trust radius, we take a step
!     to the boundary. If the Hessian index is larger than zero,
!     the level-shifted step will also be employed, provided the
!     Newton step is larger than 0.5E-3_realk. For saddle points we
!     employ the level-shift when the index is different from 1.
!
         IF (((optinfo%stepNorm .GT. optinfo%TrustRad) .AND. (.NOT. optinfo%NoTrust)) .OR. &
     &        ((.NOT. optinfo%Saddle) .AND. (optinfo%IndHes .GT. 0) .AND. &
     &        (optinfo%stepNorm .GE. 0.5E-3_realk)) .OR. (optinfo%Saddle .AND. &
     &        (optinfo%IndHes .NE. 1))) THEN
            IF (optinfo%IPrint .GE. IPRMED) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm of Newton step:', optinfo%stepNorm
               WRITE(LUPRI,'(A,F15.10/)')' Trust radius       :', optinfo%TrustRad
            END IF
            INSIDE = .FALSE.
            IF (optinfo%stepNorm .LT. optinfo%TrustRad) INSIDE = .TRUE.
            call ls_LSHFT0(optinfo%NIntCoord,optinfo%NIntCoord,&
     &  MIN(optinfo%TrustRad,optinfo%stepNorm),RNU,optinfo%ZeroGrad,INSIDE,lupri,optinfo)
            BNDNRM = optinfo%TrustRad
            optinfo%stepNorm = SQRT(DDOT(optinfo%NIntCoord,optinfo%STPDIA,1,optinfo%STPDIA,1))
            IF (optinfo%IPrint .GE. IPRMED) THEN
               WRITE(LUPRI,'(/A,F15.10)')' Norm, boundary step:', optinfo%stepNorm
            END IF
         END IF
!     For saddle point optimizations, we restore the original function.
!
         IF (optinfo%Saddle) THEN
            DOSCAL = (.NOT. optinfo%Newton)
            IF (optinfo%InitHess .AND. (optinfo%ItrNmr .EQ. 0)) DOSCAL = .FALSE.
            IF (DOSCAL) THEN
               CMPLIM = MAX(optinfo%TrustRad*0.67E0_realk, 0.5E0_realk)
               DO I = 1, optinfo%NIntCoord
                  IF (ABS(optinfo%STPDIA(I)) .GT. CMPLIM) &
     &                 optinfo%STPDIA(I) = SIGN(CMPLIM,optinfo%STPDIA(I))
               ENDDO
            END IF
            call MAKIMG(optinfo%NProjected,optinfo%NIntCoord,&
               & optinfo%EVAL,optinfo%GRDDIA,optinfo%STPDIA,IMODE,.TRUE.)
         END IF


!
!     Energy is predicted, will be used later to update trust radius.
!
         optinfo%predictedChange = DDOT(optinfo%NIntCoord,optinfo%GRDDIA,1,optinfo%STPDIA,1) &
     &        + 0.5E0_realk*dv3dot_ls(optinfo%NIntCoord,optinfo%STPDIA,optinfo%EVAL,optinfo%STPDIA)
!
!     If the predicted energy is positive, it means the Newton step is
!     towards a maximum/saddle point. We then simply reverse the
!     step direction (a bit dirty, but seems to work).
!
         IF ((.NOT. optinfo%Saddle) .AND. (optinfo%predictedChange .GT. 0.0E0_realk)) THEN
            WRITE(LUPRI,*) 'Reversing step!'
            DO 40 I = 1, optinfo%NIntCoord
               optinfo%STPDIA(I) = -optinfo%STPDIA(I)
 40         CONTINUE
            optinfo%predictedChange = DDOT(optinfo%NIntCoord,optinfo%GRDDIA,1,optinfo%STPDIA,1) &
     &           + 0.5E0_realk*dv3dot_ls(optinfo%NIntCoord,optinfo%STPDIA,optinfo%EVAL,optinfo%STPDIA)
            IF (optinfo%IPrint .GT. 2) THEN
               WRITE(LUPRI,'(A,F25.15)') &
     &              ' New pred. energy change',optinfo%predictedChange
            END IF
         END IF
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            IF (optinfo%IPrint .GE. IPRDBG) THEN
               call lsheader(lupri,'Internal diagonal gradient')
               call output(optinfo%GRDDIA,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
               call lsheader(lupri,'Eigenvalues')
               call output(optinfo%EVAL,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
            END IF
            IF (optinfo%IPrint .GE. IPRMED) THEN
               call lsheader(lupri,'Internal diagonal step')
               call output(optinfo%STPDIA,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
            END IF
            WRITE (LUPRI,'(/A,F25.15)') ' Predicted energy change', &
     &           optinfo%predictedChange
         END IF
!
!     The diagonal step is converted to ordinary internal
!     coordinates.
!
         call ls_DZERO(optinfo%STPINT,MXRCRD)
         DO 20 I = 1, optinfo%NIntCoord
            DO 22 J = 1, optinfo%NIntCoord
               optinfo%STPINT(I) = optinfo%STPINT(I) + optinfo%STPDIA(J)*EVEC(I,J)
 22         CONTINUE
!
!     For angles and dihedral angles we have to avoid step components
!     giving multiples of 2*pi.
!
            IF (optinfo%RedInt .AND. optinfo%INTCRD(I,1) .GT. 10) &
     &           optinfo%STPINT(I) = MOD(optinfo%STPINT(I),2.0E0_realk*PI)
            IF (ABS(optinfo%STPINT(I)) .LE. 1.0E-8_realk) optinfo%STPINT(I) = D0
 20      CONTINUE
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            call lsheader(lupri,'Step in internal coordinates')
            call output(optinfo%STPINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         END IF
!
!     The rational function method
!
      ELSE IF (optinfo%RatFun) THEN
!!!!!! Vladimir: saddle point code disabled in LSDALTON
!!!!!!
!         IF (optinfo%Saddle) THEN
!            call ls_PRFSTI(MXRCRD,MX2CRD,NCRDHS,HESINT,EVEC, &
!     &           TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,VECMOD,lupri,optinfo)
!         ELSE
            call ls_RFSTP(MX2CRD,NCRDHS,MXRCRD,optinfo%NIntCoord,EVEC, &
     &  optinfo%STPINT,optinfo%GRDINT,TMPMAT,HESINT,lupri,optinfo)
!         END IF
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)') &
     &           ' Predicted energy change', optinfo%predictedChange
         END IF
!
!     The Geometrical DIIS method
!
      ELSE IF (optinfo%GDIIS) THEN
         call ls_DZERO(TMPMAT,MX2CRD*MX2CRD)
         call ls_DZERO(TMPMT4,MX2CRD*MX2CRD)
!
!     First we have to construct the inverse Hessian.
!
         DO 210 I = 1, optinfo%NIntCoord
            IF (ABS(optinfo%EVAL(I)) .GE. 1.0E-6_realk) THEN
               DO 212 J = 1, optinfo%NIntCoord
                  TMPMT4(I,J) = EVEC(J,I)/optinfo%EVAL(I)
 212           CONTINUE
            END IF
 210     CONTINUE
         DO 215 I = 1, optinfo%NIntCoord
            DO 217 J = 1, optinfo%NIntCoord
               DO 219 K = 1, optinfo%NIntCoord
                  TMPMAT(I+(J-1)*optinfo%NIntCoord) = TMPMAT(I+(J-1)*optinfo%NIntCoord) + &
     &                 EVEC(I,K)*TMPMT4(K,J)
 219           CONTINUE
 217        CONTINUE
 215     CONTINUE
!
!     Then the DIIS-step is determined
!
         call ls_GDISTP(MXRCRD,optinfo%NIntCoord,MXRCRD,&
     &        MX2CRD,optinfo%STPDIA,optinfo%GRDINT,HESINT, &
     &        TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,GRDARR,&
     &        STPARR,lupri,optinfo)
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            WRITE (LUPRI,'(/A,F25.15)') &
     &           ' Predicted energy change', optinfo%predictedChange
         END IF
         DO 250 I = 1, optinfo%NIntCoord
            IF (ABS(optinfo%STPDIA(I)) .GE. 1.0E-6_realk) THEN
               optinfo%STPINT(I) = optinfo%STPDIA(I)
            ELSE
               optinfo%STPINT(I) = D0
            END IF
 250     CONTINUE
      END IF
      IF (ACTIVE) THEN
         DO 300 I = 1, optinfo%NIntCoord
            optinfo%STPINT(I) = optinfo%STPINT(I) + STPLIN(I)
 300     CONTINUE
         optinfo%predictedChange = optinfo%predictedChange + (EMOD-optinfo%energy)
         WRITE (LUPRI,'(/A,F25.15)') &
     &        ' Modified energy prediction due to line search', optinfo%predictedChange
      END IF
!
      optinfo%stepNorm = SQRT(DDOT(optinfo%NIntCoord,optinfo%STPINT,1,optinfo%STPINT,1))
!
!     Find Cartesian step vector
!

!     Allocate temporary memory
      Call mem_alloc(CRDORG,MXCENT,8)
      Call mem_alloc(CRDOLD,MXCENT,8)
      Call mem_alloc(CRDNEW,MXCENT,8)
      Call mem_alloc(VALORG,MXCOOR)
      Call mem_alloc(VALOLD,MXCOOR)
      Call mem_alloc(VALNEW,MXCOOR)
      Call mem_alloc(DIFF,MXCOOR)
      Call mem_alloc(TMPVEC,MXCOOR)
      Call mem_alloc(CSTRA,MXCOOR,MXCOOR)
      Call mem_alloc(SCTRA,MXCOOR,MXCOOR)
!
      CRDORG = 0E0_realk
      CRDOLD = 0E0_realk
      CRDNEW = 0E0_realk
      VALORG = 0E0_realk
      VALOLD = 0E0_realk
      VALNEW = 0E0_realk
      DIFF = 0E0_realk
      TMPVEC = 0E0_realk
      CSTRA = 0E0_realk
      SCTRA = 0E0_realk
      call ls_STPI2C(Molecule,MXRCRD,CRDORG, &
     &     CRDOLD,CRDNEW ,VALORG, &
     &     VALOLD,VALNEW,DIFF,TMPVEC, &
     &     CSTEP,BMTRAN,BMTINV,CSTRA, &
     &     SCTRA,.FALSE.,lupri,optinfo)
!
!     If the optimization is constrained, we take a small step to
!     reimpose the constraints.
!
      IF (optinfo%ConOpt) THEN
      call ls_STPI2C(Molecule,MXRCRD,CRDORG, &
     &     CRDOLD,CRDNEW,VALORG, &
     &     VALOLD,VALNEW,DIFF,TMPVEC, &
     &     CSTEP,BMTRAN,BMTINV,CSTRA, &
     &     CSTRA,.FALSE.,lupri,optinfo)
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Internal values before step')
            call output(TMPMT2(1),1,1,1,optinfo%NIntCoord,1,MX2CRD,1,LUPRI)
         END IF
         IF (optinfo%IPrint .GE. IPRMED) THEN
            call lsheader(lupri,'Internal values after step')
            call output(TMPMT2(2*MX2CRD+1),1,1,1,optinfo%NIntCoord,1, &
     &           MX2CRD,1,LUPRI)
         END IF
      ELSE
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            NRIC = optinfo%NIntCoord
            call lsheader(lupri,'Internal values before step')
            !call output(TMPMT2(MX2CRD+1),1,1,1,optinfo%NIntCoord,1, &
            call output(VALOLD,1,1,1,optinfo%NIntCoord,1, &
     &           MX2CRD,1,LUPRI)
         END IF
         IF (optinfo%IPrint .GE. IPRMED) THEN
            call lsheader(lupri,'Internal values after step')
            !call output(TMPMT2(2*MX2CRD+1),1,1,1,optinfo%NIntCoord,1, &
            call output(VALNEW,1,1,1,optinfo%NIntCoord,1, &
     &           MX2CRD,1,LUPRI)
         END IF
      END IF
!
!     After the step is corrected, energy is predicted again, will be used later to update trust radius.
!     The same is done for the step norm
!
      optinfo%predictedChange = DDOT(optinfo%NIntCoord,optinfo%GRDDIA,1,optinfo%STPDIA,1) &
     &        + 0.5E0_realk*dv3dot_ls(optinfo%NIntCoord,optinfo%STPDIA,optinfo%EVAL,optinfo%STPDIA)
!      optinfo%stepNorm = SQRT(DDOT(optinfo%NIntCoord,optinfo%STPINT,1,optinfo%STPINT,1))
! Deallocate memory
      Call mem_dealloc(CRDORG)
      Call mem_dealloc(CRDOLD)
      Call mem_dealloc(CRDNEW)
      Call mem_dealloc(VALORG)
      Call mem_dealloc(VALOLD)
      Call mem_dealloc(VALNEW)
      Call mem_dealloc(DIFF)
      Call mem_dealloc(TMPVEC)
      Call mem_dealloc(CSTRA)
      Call mem_dealloc(SCTRA)
      RETURN
      END
!=====================!
! Merge_geometry      !
!=====================!
!  Finds how to displace atoms to make the chahge
!  in one internal (during scan) not affect the others
Subroutine Merge_Geometry(optinfo,Scan_Coord)
use precision
use ls_util 
use optimization_input 
use files
use Fundamental
use memory_handling
Implicit none
Type(opt_setting) :: optinfo
Integer :: Scan_Coord, i,j,k, atom,NAtoms,New
Logical :: GoOn
!
NAtoms = optinfo%ICartCoord/3
! Allocate and initialize Atoms_to_move 
Call mem_alloc(optinfo%Atoms_to_move,2,NAtoms)
optinfo%Atoms_to_move = 0
! Set "active" atoms, bonds for now
If (optinfo%INTCRD(Scan_Coord,1) .EQ. 1) then 
   optinfo%Atoms_to_move(1,1) = optinfo%INTCRD(Scan_Coord,2)
   optinfo%Atoms_to_move(2,1) = optinfo%INTCRD(Scan_Coord,3)
Endif
! Loop over internals to find all atoms bound to the "active" ones
Do atom = 1,2
   optinfo%N_to_move(atom) = 1
   Do i = 1, NAtoms   ! Loop over atoms to move
      Do j = 1, optinfo%NIntCoord ! Loop over internal coordinates
       If (optinfo%INTCRD(j,1) .EQ. 1) then ! A coordinate should be a bond
         If (j .NE. Scan_Coord) then 
            ! Find whether a bond is formed by a known atom to move 
            If ( (optinfo%INTCRD(j,2) .EQ. optinfo%Atoms_to_move(atom,i)) .OR. &
               & (optinfo%INTCRD(j,3) .EQ. optinfo%Atoms_to_move(atom,i)) ) then
               If ( optinfo%INTCRD(j,2) .EQ. optinfo%Atoms_to_move(atom,i)  ) then
                  New = optinfo%INTCRD(j,3)
               Else
                  New = optinfo%INTCRD(j,2)
               Endif
               ! Check whether the new atom_to move is not yet referenced
               Do k = 1,NAtoms
                  If (New .NE. optinfo%Atoms_to_move(atom,k)) then
                     GoOn = .TRUE.
                  Else
                     GoOn = .FALSE.
                  Endif
                  If (GoOn .EQV. .FALSE.) Exit   ! Atom referenced: it's not a new atom_to_move
               Enddo
               ! If new is not referenced than we add it to atoms_to_move
               Write(*,*) New,optinfo%N_to_move
               If (GoOn .EQV. .TRUE.) then
                   optinfo%N_to_move(atom) = optinfo%N_to_move(atom) + 1
                   optinfo%Atoms_to_move(atom,optinfo%N_to_move(atom)) = New
               Endif ! If not referenced
            Endif ! If atom contributes to bond
         Endif ! If not Scan_coord
       Endif  ! If bond
      Enddo
   Enddo
Enddo
!
End subroutine Merge_geometry
      
!  /* Deck getint */
      SUBROUTINE LS_GETINT(IATOM,MXRCRD,ATMCRD,VALINT,lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use Fundamental
!
!     Determines the value of all redundant internal coordinates.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(opt_setting) :: optinfo
      Real(realk) ATMCRD(MXCENT,8), VALINT(MXRCRD)
      Real(realk) VEC1(3), VEC2(3), VEC3(3), VEC4(3), VEC5(3)
      PARAMETER (D0 = 0.0E0_realk, DEG179 = 179E0_realk*PI/180E0_realk)
!     Remove comment mark after inserting to the LSDALTON framework
!      call ls_DZERO(VALINT,MXRCRD)
      NRIC = optinfo%NIntCoord
      IF (optinfo%DelInt) NRIC = optinfo%NIntCoord
      DO 10 IC = 1, NRIC
         IF (optinfo%INTCRD(IC,1) .LT. 10) THEN
            VEC1(1) = ATMCRD(optinfo%INTCRD(IC,2),2) - ATMCRD(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(optinfo%INTCRD(IC,2),3) - ATMCRD(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(optinfo%INTCRD(IC,2),4) - ATMCRD(optinfo%INTCRD(IC,3),4)
            VALINT(IC) = SQRT(DDOT(3,VEC1,1,VEC1,1))
         ELSE IF (optinfo%INTCRD(IC,1) .LT. 20) THEN
            VEC1(1) = ATMCRD(optinfo%INTCRD(IC,2),2) - ATMCRD(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(optinfo%INTCRD(IC,2),3) - ATMCRD(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(optinfo%INTCRD(IC,2),4) - ATMCRD(optinfo%INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(optinfo%INTCRD(IC,4),2) - ATMCRD(optinfo%INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(optinfo%INTCRD(IC,4),3) - ATMCRD(optinfo%INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(optinfo%INTCRD(IC,4),4) - ATMCRD(optinfo%INTCRD(IC,3),4)
!
!     Regular angles
!
            IF (optinfo%INTCRD(IC,1) .EQ. 11) THEN
               call ls_VECPRD(VEC1,VEC2,VEC3)
               VNRM = SQRT(DDOT(3,VEC3,1,VEC3,1))
               IF (VNRM .LE. 1.0E-8_realk) THEN
                  VEC3(1) =  VEC1(2)+VEC1(3)
                  VEC3(2) = -VEC1(1)+VEC1(3)
                  VEC3(3) = -VEC1(1)-VEC1(2)
               END IF
               call ls_NRMLVC(VEC3)
!
!     Second coordinate of angles larger than 175 degrees.
!
            ELSE
               call ls_VECPRD(VEC1,VEC3,VEC4)
               call ls_NRMLVC(VEC4)
               VEC3(1) = VEC4(1)
               VEC3(2) = VEC4(2)
               VEC3(3) = VEC4(3)
            END IF
            call ls_VECPRD(VEC3,VEC1,VEC4)
            call ls_NRMLVC(VEC4)
            VALINT(IC) = vecang_ls(VEC1,VEC2)
            IF ((optinfo%INTCRD(IC,1) .EQ. 12) .OR. (optinfo%INTCRD(IC+1,1) .EQ. 12)) &
     &           VALINT(IC) = vecang_ls(VEC1,VEC4) + vecang_ls(VEC4,VEC2)
!
!     Dihedral angles
!
         ELSE
            VEC1(1) = ATMCRD(optinfo%INTCRD(IC,2),2) - ATMCRD(optinfo%INTCRD(IC,3),2)
            VEC1(2) = ATMCRD(optinfo%INTCRD(IC,2),3) - ATMCRD(optinfo%INTCRD(IC,3),3)
            VEC1(3) = ATMCRD(optinfo%INTCRD(IC,2),4) - ATMCRD(optinfo%INTCRD(IC,3),4)
            VEC2(1) = ATMCRD(optinfo%INTCRD(IC,4),2) - ATMCRD(optinfo%INTCRD(IC,3),2)
            VEC2(2) = ATMCRD(optinfo%INTCRD(IC,4),3) - ATMCRD(optinfo%INTCRD(IC,3),3)
            VEC2(3) = ATMCRD(optinfo%INTCRD(IC,4),4) - ATMCRD(optinfo%INTCRD(IC,3),4)
            VEC3(1) = ATMCRD(optinfo%INTCRD(IC,5),2) - ATMCRD(optinfo%INTCRD(IC,4),2)
            VEC3(2) = ATMCRD(optinfo%INTCRD(IC,5),3) - ATMCRD(optinfo%INTCRD(IC,4),3)
            VEC3(3) = ATMCRD(optinfo%INTCRD(IC,5),4) - ATMCRD(optinfo%INTCRD(IC,4),4)
            call ls_NRMLVC(VEC2)
            CMPNT1 = DDOT(3,VEC1,1,VEC2,1)
            VEC1(1) = VEC1(1) - CMPNT1*VEC2(1)
            VEC1(2) = VEC1(2) - CMPNT1*VEC2(2)
            VEC1(3) = VEC1(3) - CMPNT1*VEC2(3)
            call ls_NRMLVC(VEC1)
            CMPNT2 = DDOT(3,VEC3,1,VEC2,1)
            VEC3(1) = VEC3(1) - CMPNT2*VEC2(1)
            VEC3(2) = VEC3(2) - CMPNT2*VEC2(2)
            VEC3(3) = VEC3(3) - CMPNT2*VEC2(3)
            call ls_NRMLVC(VEC3)
            IF ((DDOT(3,VEC1,1,VEC1,1) .GT. 1.0E-16_realk) .AND. &
     &           (DDOT(3,VEC3,1,VEC3,1) .GT. 1.0E-16_realk)) THEN
               VALINT(IC) = vecang_ls(VEC1,VEC3)
               IF (ABS(VALINT(IC)) .GT. DEG179) THEN
                  call ls_VECPRD(VEC1,VEC2,VEC4)
                  call ls_NRMLVC(VEC4)
                  DCMP2 = vecang_ls(VEC4,VEC3)
                  IF (DCMP2 .GE. PI/2E0_realk) DCMP2 = PI - DCMP2
                  VALINT(IC) = vecang_ls(VEC1,VEC4) + DCMP2
               END IF
               call ls_VECPRD(VEC1,VEC3,VEC4)
               IF (DDOT(3,VEC2,1,VEC4,1) .LT. 0.0E0_realk) &
     &              VALINT(IC) = -VALINT(IC)
            ELSE
               VALINT(IC) = D0
            END IF
         END IF
 10   CONTINUE
      RETURN
      END

!  /* Deck stpi2c */
      SUBROUTINE LS_STPI2C(Molecule,MXRCRD,CRDORG,CRDOLD,CRDNEW, &
     &   VALORG,VALOLD,VALNEW,DIFF,TMPVEC,CSTEP,BMTRAN,BMTINV, &
     &     CSTRA,SCTRA,CORREC,lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use Fundamental
use q_to_x_mod
use molecule_type
use molecule_typetype
!
!     Transforms internal step to Cartesian step iteratively.
!     The logical parameter CORREC indicates if the internal step
!     is a small correctional step to reimpose constraints in a
!     constrained optimization.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) CRDORG(MXCENT,8),CRDOLD(MXCENT,8)
      Real(realk) CRDNEW(MXCENT,8), VALOLD(MXRCRD), VALORG(MXRCRD)
      Real(realk) VALNEW(MXRCRD), DIFF(MXRCRD)
      Real(realk) TMPVEC(MXCOOR), CSTEP(MXCOOR)
      Real(realk) BMTRAN(MXRCRD,MXRCRD), BMTINV(MXRCRD,MXCOOR)
      Real(realk) CSTRA(MXCOOR*MXCOOR), SCTRA(MXCOOR*MXCOOR)
      LOGICAL CORREC, ADJUST
      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk)
      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
      PARAMETER (ITRLIM = 25)
      CHARACTER STPTXT*16
!     Remove by merging!
      Real(realk) STEPINT(MXCOOR)
!
!     Apply the new polynomial step-converter (back-transformation)
If (optinfo%New_stepping .OR. optinfo%IBT .OR. optinfo%OldIBT) then 
    Call Back_transform(CSTEP,optinfo%ICartCoord,optinfo%NIntCoord,lupri,optinfo)
Else
! Write the topology file if needed
      Call write_topology(optinfo)
!
      ADJUST = .FALSE.
      NRIC = optinfo%NIntCoord
      NNNRIC = 0
      IF (optinfo%DelInt) THEN
         NRIC = optinfo%NIntCoord
         NNNRIC = optinfo%NIntCoord
      END IF
!
!     To be able to properly control the step we must use
!     the primitive (redundant) internal coordinates. If we use
!     delocalized internals, we transform the step to redundant
!     internals. VALORG is used to store the non-redundant step.
!
      IF (optinfo%DelInt) THEN
         call ls_DZERO(VALORG,MXRCRD)
         DO 100 I = 1, NNNRIC
            VALORG(I) = optinfo%STPINT(I)
 100     CONTINUE
         call ls_DZERO(optinfo%STPINT,MXRCRD)
         DO 105 I = 1, NRIC
            DO 107 J = 1, NNNRIC
               optinfo%STPINT(I) = optinfo%STPINT(I) + BMTRAN(I,J)*VALORG(J)
 107        CONTINUE
 105     CONTINUE
      END IF
!
!     We have to find the value of all internal coordinates for the
!     old geometry. If a correctional step is requested, we copy
!     the values of the internal coordinates after the major step.
!
      IATOM = optinfo%ICartCoord / 3
      IF (.NOT. CORREC) THEN
      call Atom_Ini(CRDOLD,Molecule,optinfo,IATOM,.TRUE.,lupri)
      call ls_GETINT(IATOM,MXRCRD,CRDOLD,VALOLD,lupri,optinfo)
         RMSLIM = 1.0E-6_realk
      ELSE
!        ... BUGFIX hjaaj Oct 07, IATOM was not defined.
!            Is used for TMPVEC below, which has Real(realk) optinfo%ICartCoord,
!            based on this I divided by three to get IATOM (see DO 15 loop)
         DO 3 I = 1, NRIC
            VALORG(I) = VALOLD(I)
            VALOLD(I) = VALNEW(I)
 3       CONTINUE
         DO 5 I = 1, MXCENT
            DO 7 J = 1, 4
               CRDORG(I,J) = CRDOLD(I,J)
               CRDOLD(I,J) = CRDNEW(I,J)
 7          CONTINUE
 5       CONTINUE
         RMSLIM = 1.0E-9_realk
      END IF
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Previous coordinates')
         call output(CRDOLD,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         IF (optinfo%DelInt) THEN
            call lsheader(lupri,'Previous primitive internal values')
         ELSE
            call lsheader(lupri,'Previous internal values')
         END IF
         call output(VALOLD,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         IF (optinfo%DelInt) THEN
            call lsheader(lupri,'Step to take (primitives)')
         ELSE
            call lsheader(lupri,'Step to take')
         END IF
         call output(optinfo%STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
!
!     First estimate of Cartesian step.
!
      call ls_DZERO(TMPVEC,MXCOOR)
      IF (.NOT. CORREC) THEN
         IF (optinfo%DelInt) THEN
            DO 110 I = 1, optinfo%ICartCoord
               DO 111 J = 1, NNNRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*VALORG(J)
 111           CONTINUE
 110         CONTINUE
          ELSE
             DO 10 I = 1, optinfo%ICartCoord
               DO 11 J = 1, NRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*optinfo%STPINT(J)
 11            CONTINUE
 10         CONTINUE
         END IF
!
      ELSE
!
!     SCTRA is used for temporary storage
!
         DO 12 IC = 1, NRIC
            IF (optinfo%IConstr(IC) .GT. 0) THEN
               SCTRA(IC) = optinfo%CoordInt(IC) - VALOLD(IC)
               IF (optinfo%INTCRD(IC,1) .GT. 10) &
     &              SCTRA(IC) = MOD(SCTRA(IC),2.0E0_realk*PI)
               IF ((optinfo%INTCRD(IC,1) .GT. 20) .AND. &
     &              (ABS(SCTRA(IC)) .GT. PI)) THEN
                  IF (SCTRA(IC) .GT. 0.0E0_realk) THEN
                     SCTRA(IC) = SCTRA(IC) - 2.0E0_realk*PI
                  ELSE
                     SCTRA(IC) = SCTRA(IC) + 2.0E0_realk*PI
                  END IF
               END IF
            ELSE
               SCTRA(IC) = 0.0E0_realk
            END IF
 12      CONTINUE
         DO 13 I = 1, optinfo%ICartCoord
            DO 14 J = 1, NRIC
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*SCTRA(J)
 14         CONTINUE
 13      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Projected step')
            call output(SCTRA,1,1,1,NRIC,1,MXCOOR,1,LUPRI)
         END IF
      END IF
      ! 
      ! We combine the HOBIT and the iterative back transformation
      !
      If (optinfo%New_stepping) then
         Do i = 1,optinfo%ICartCoord
            TMPVEC(i) = CSTEP(i)
         Enddo
      Endif
      !
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Inverse of B^t')
         call output(BMTINV,1,optinfo%NIntCoord,1,optinfo%ICartCoord,MXRCRD,MXCOOR,1,LUPRI)
         call lsheader(lupri,'First estimate of Cartesian step')
         call output(TMPVEC,1,1,1,optinfo%ICartCoord,1,MXCOOR,1,LUPRI)
      END IF
      ITRCRD = 1
      RMSOLD = D0
      RMS1ST = -D1
 123  CONTINUE
      DO 15 I = 1, IATOM
         CRDNEW(I,1) = CRDOLD(I,1)
         DO 17 J = 1, 3
            CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 17      CONTINUE
 15   CONTINUE
      call ls_GETINT(IATOM,MXRCRD,CRDNEW,VALNEW,lupri,optinfo)
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'New coordinates')
         call output(CRDNEW,1,IATOM,1,4,MXCENT,8,1,LUPRI)
         IF (optinfo%DelInt) THEN
            call lsheader(lupri,'New primitive internal values')
         ELSE
            call lsheader(lupri,'New internal values')
         END IF
         call output(VALNEW,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
      call ls_DZERO(DIFF,MXRCRD)
      RMSINT = D0
      IF (.NOT. CORREC) THEN
         DO 18 I = 1, NRIC
            DIFF(I) = optinfo%STPINT(I) - (VALNEW(I) - VALOLD(I))
            IF (optinfo%INTCRD(I,1) .GT. 10) &
     &           DIFF(I) = MOD(DIFF(I),2.0E0_realk*PI)
            IF ((optinfo%INTCRD(I,1) .GT. 20) .AND. &
     &           (ABS(DIFF(I)) .GT. PI)) THEN
               IF (DIFF(I) .GT. 0.0E0_realk) THEN
                  DIFF(I) = DIFF(I) - 2.0E0_realk*PI
               ELSE
                  DIFF(I) = DIFF(I) + 2.0E0_realk*PI
               END IF
            END IF
            IF (ABS(DIFF(I)) .LT. 1.0E-14_realk) DIFF(I) = D0
 18      CONTINUE
         RMSINT = DDOT(NRIC,DIFF,1,DIFF,1)
      ELSE
         DO 19 I = 1, NRIC
            DIFF(I) = (optinfo%CoordInt(I) - VALNEW(I))*(optinfo%IConstr(I)*1.0E0_realk)
            IF (optinfo%INTCRD(I,1) .GT. 10) DIFF(I) = MOD(DIFF(I),2.0E0_realk*PI)
            IF ((optinfo%INTCRD(I,1) .GT. 20) .AND. (ABS(DIFF(I)) .GT. PI)) THEN
               IF (DIFF(I) .GT. 0.0E0_realk) THEN
                  DIFF(I) = DIFF(I) - 2.0E0_realk*PI
               ELSE
                  DIFF(I) = DIFF(I) + 2.0E0_realk*PI
               END IF
            END IF
            IF (ABS(DIFF(I)) .LT. 1.0E-14_realk) DIFF(I) = D0
            RMSINT = RMSINT + DIFF(I)*DIFF(I)
 19      CONTINUE
      END IF
!
!     The difference wich will be used for the next iteration is
!     transformed back to delocalized internals.
!
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         IF (optinfo%DelInt) THEN
            call lsheader(lupri,'Difference to wanted step (primitives)')
         ELSE
            call lsheader(lupri,'Difference to wanted step')
         END IF
         call output(DIFF,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
      END IF
      IF (optinfo%DelInt) THEN
         call ls_DZERO(SCTRA,MXCOOR*MXCOOR)
         RMSINT = 0.0E0_realk
         DO 200 I = 1, NRIC
            SCTRA(I) = DIFF(I)
 200     CONTINUE
         call ls_DZERO(DIFF,MXRCRD)
         DO 205 I = 1, NNNRIC
            DO 207 J = 1, NRIC
               DIFF(I) = DIFF(I) + BMTRAN(J,I)*SCTRA(J)
 207        CONTINUE
            RMSINT = RMSINT + DIFF(I)*DIFF(I)
 205     CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            call lsheader(lupri,'Difference to wanted step')
            call output(DIFF,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         END IF
      END IF
      RMSINT = SQRT(RMSINT/(1.0E0_realk*NRIC))
      IF (RMS1ST .LT. D0) RMS1ST = RMSINT
      DIFRMS = ABS(RMSOLD-RMSINT)
      RMSOLD = RMSINT
      IF (optinfo%IPrint .GE. IPRDBG) THEN
         call lsheader(lupri,'Root-mean-square of difference')
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',RMSINT
         call lsheader(lupri,'Change in root-mean-square of difference')
         WRITE(LUPRI,'(A,G16.6)') 'Value:    ',DIFRMS
      END IF
      IF ((RMSINT .GE. RMSLIM) .AND. (RMSINT .LE. 1.0E2_realk) &
     &     .AND. (DIFRMS .GE. 1.0E-12_realk) .AND. (ITRCRD .LE. ITRLIM)) THEN
         DO 20 I = 1, optinfo%ICartCoord
            DO 22 J = 1, optinfo%NIntCoord
               TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*DIFF(J)
 22         CONTINUE
 20      CONTINUE
         IF (optinfo%IPrint .GE. IPRDBG) THEN
            STPTXT = 'Updated step #XX'
            WRITE(STPTXT(15:16),'(I2)') ITRCRD
            call lsheader(lupri,STPTXT)
            call output(TMPVEC,1,1,1,optinfo%ICartCoord,1,MXCOOR,1,LUPRI)
         END IF
         ITRCRD = ITRCRD + 1
         GOTO 123
!      ELSE IF (((ITRCRD .GE. ITRLIM) .AND. (RMS1ST .LE. 1.0E2_realk*RMSINT))
!     &        .OR. ((ITRCRD .LT. ITRLIM) .AND.
!     &        (RMS1ST .LE. 1.0E2_realk*RMSINT) .AND.
!     &        (RMSINT .GE. 1.0E2_realk*RMSLIM))
!     &        .OR. (RMSINT .GE. 1.0E0_realk)) THEN
      ELSE IF ((ITRCRD .GE. ITRLIM) .OR. ((RMS1ST .LE. 1.0E2_realk*RMSINT) &
     &        .AND. (RMSINT .GE. 1.1E0_realk*RMSLIM)) .OR. &
     &        (RMSINT .GE. 1.0E0_realk)) THEN
!      ELSE IF ((ITRCRD .GE. ITRLIM) .OR. (RMSINT .GE. 1.0E2_realk)) THEN
         IF (optinfo%IPrint .GE. IPRMIN) THEN
            WRITE(LUPRI,*) &
     &           'Step does not converge, reverting to first estimate.'
            WRITE(LUPRI,*)
         END IF
         call ls_DZERO(TMPVEC,MXCOOR)
         IF (.NOT. CORREC) THEN
            IF (optinfo%DelInt) THEN
               DO 130 I = 1, optinfo%ICartCoord
                  DO 132 J = 1, NNNRIC
                     TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*VALORG(J)
 132              CONTINUE
 130           CONTINUE
            ELSE
               DO 30 I = 1, optinfo%ICartCoord
                  DO 32 J = 1, NRIC
                     TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*optinfo%STPINT(J)
 32               CONTINUE
 30            CONTINUE
            END IF
         ELSE
            DO 35 I = 1, optinfo%ICartCoord
               DO 37 J = 1, NRIC
                  TMPVEC(I) = TMPVEC(I) + BMTINV(J,I)*SCTRA(J)
 37            CONTINUE
 35         CONTINUE
         END IF
         DO 40 I = 1, IATOM
            DO 42 J = 1, 3
               CRDNEW(I,J+1) = CRDOLD(I,J+1) + TMPVEC((I-1)*3+J)
 42         CONTINUE
 40      CONTINUE
         call ls_GETINT(IATOM,MXRCRD,CRDNEW,VALNEW,lupri,optinfo)
         DO 45 I = 1, NRIC
            optinfo%STPINT(I) = VALNEW(I) - VALOLD(I)
            IF (optinfo%INTCRD(I,1) .GT. 10) optinfo%STPINT(I) = MOD(optinfo%STPINT(I),2.0E0_realk*PI)
            IF ((optinfo%INTCRD(I,1) .GT. 20) &
     &           .AND. (ABS(optinfo%STPINT(I)) .GT. PI)) THEN
               IF (optinfo%STPINT(I) .GT. 0.0E0_realk) THEN
                  optinfo%STPINT(I) = optinfo%STPINT(I) - 2.0E0_realk*PI
               ELSE
                  optinfo%STPINT(I) = optinfo%STPINT(I) + 2.0E0_realk*PI
               END IF
            END IF
 45      CONTINUE
!
!     If we use delocalized internal coordinates, we have to
!     transform the values. SCTRA is used for temporary storage.
!
         IF (optinfo%DelInt) THEN
            ADJUST = .TRUE.
            call ls_DZERO(SCTRA,MXRCRD)
            DO 300 I = 1, NRIC
               SCTRA(I) = optinfo%STPINT(I)
 300        CONTINUE
            call ls_DZERO(optinfo%STPINT,MXRCRD)
            DO 305 I = 1, NNNRIC
               DO 307 J = 1, NRIC
                  optinfo%STPINT(I) = optinfo%STPINT(I) + BMTRAN(J,I)*SCTRA(J)
 307           CONTINUE
 305        CONTINUE
         END IF
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Adjusted step')
            call output(optinfo%STPINT,1,1,1,optinfo%NIntCoord,1,MXRCRD,1,LUPRI)
         END IF
      END IF
      IF (CORREC) THEN
         call ls_GETINT(IATOM,MXRCRD,CRDNEW,VALNEW,lupri,optinfo)
         DO 50 I = 1, NRIC
            optinfo%STPINT(I) = VALNEW(I) - VALORG(I)
            IF (optinfo%INTCRD(I,1) .GT. 10) optinfo%STPINT(I) = MOD(optinfo%STPINT(I),2.0E0_realk*PI)
            IF ((optinfo%INTCRD(I,1) .GT. 20) &
     &           .AND. (ABS(optinfo%STPINT(I)) .GT. PI)) THEN
               IF (optinfo%STPINT(I) .GT. 0.0E0_realk) THEN
                  optinfo%STPINT(I) = optinfo%STPINT(I) - 2.0E0_realk*PI
               ELSE
                  optinfo%STPINT(I) = optinfo%STPINT(I) + 2.0E0_realk*PI
               END IF
            END IF
 50      CONTINUE
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Adjusted step')
            call output(optinfo%STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         END IF
         DO 60 I = 1, IATOM
            DO 62 J = 1, 3
               TMPVEC((I-1)*3+J) = CRDNEW(I,J+1) - CRDORG(I,J+1)
 62         CONTINUE
 60      CONTINUE
      END IF
!
!     The step in delocs is moved back to optinfo%STPINT
!
      IF (optinfo%DelInt .AND. (.NOT. ADJUST)) THEN
         call ls_DZERO(optinfo%STPINT,MXRCRD)
         DO 800 I = 1, optinfo%NIntCoord
            optinfo%STPINT(I) = VALORG(I)
 800     CONTINUE
      END IF
!     Remove the line at merging!!!!!
      STEPINT = optinfo%STPINT
!     Find the cartesian step
      call ls_DZERO(optinfo%STPSYM,MXRCRD)
!
      optinfo%STPSYM(1:MXCOOR) = TMPVEC(1:MXCOOR)
      DO 90 I = 1, optinfo%ICartCoord
         optinfo%STPSYM(I) = optinfo%STPSYM(I)
 90   CONTINUE
!     Remove the 2 lines at merging!!!!!
      optinfo%STPINT = STEPINT
!      call ls_WLKCOR(CSTEP,optinfo%ICartCoord,MXCOOR,optinfo%IPrint,lupri,optinfo)
      Do i = 1,MXCOOR
         CSTEP(i)=optinfo%STPSYM(i)
      Enddo
!
!     Experimental extra stabilization by freezing atoms.
!     All Cartesian step vector components are simply set to zero
!     for the selected frozen atoms. Note that this is only applied
!     to the correctional step
!
      IF (CORREC .AND. (optinfo%NFreeze .GT. 0)) THEN
         DO 510 I = 1, optinfo%NFreeze
            DO 515 J = 1, 3
               CSTEP(3*(optinfo%FreezeArray(I)-1)+J) = 0E0_realk
 515        CONTINUE
 510     CONTINUE
!     Calculate new internal coordinates
         DO 530 I = 1, IATOM
            DO 535 J = 1, 3
               CRDNEW(I,J+1) = CRDORG(I,J+1)+ CSTEP(3*(I-1)+J)
 535        CONTINUE
 530     CONTINUE
         call ls_GETINT(IATOM,MXRCRD,CRDNEW,VALNEW,lupri,optinfo)
         DO 550 I = 1, NRIC
            optinfo%STPINT(I) = VALNEW(I) - VALORG(I)
            IF (optinfo%INTCRD(I,1) .GT. 10) optinfo%STPINT(I) = MOD(optinfo%STPINT(I),2.0E0_realk*PI)
            IF ((optinfo%INTCRD(I,1) .GT. 20) &
     &           .AND. (ABS(optinfo%STPINT(I)) .GT. PI)) THEN
               IF (optinfo%STPINT(I) .GT. 0.0E0_realk) THEN
                  optinfo%STPINT(I) = optinfo%STPINT(I) - 2.0E0_realk*PI
               ELSE
                  optinfo%STPINT(I) = optinfo%STPINT(I) + 2.0E0_realk*PI
               END IF
            END IF
 550     CONTINUE
         IF (optinfo%IPrint .GE. IPRMAX) THEN
            call lsheader(lupri,'Adjusted step after freezing')
            call output(optinfo%STPINT,1,1,1,NRIC,1,MXRCRD,1,LUPRI)
         END IF
!     If .FRZITR has been specified, the freezing is turned off after
!     the requested number of iterations
         IF ((optinfo%IterFreeze .GT. 0) .AND. ((optinfo%ItrNmr+1) .GE. optinfo%IterFreeze)) &
     &        optinfo%NFreeze = 0
      END IF
! Polynomial stepping
Endif
      RETURN 
      END

!!!!!! Vladimir: saddle point code disabled in LSDALTON
!!!!!!
!  /* Deck prfsti */
!      SUBROUTINE LS_PRFSTI(MXRCRD,MX2CRD,NCRDHS,HESINT,EVEC, &
!     &     TMPMAT,TMPMT2,TMPMT3,TMPMT4,TMPMT5,VECMOD,lupri,optinfo)
!use ls_util 
!use optimization_input 
!use files
!
!     Controls saddle point optimization in redundant internal
!     coordinates using the partitioned rational function approach.
!
!Implicit Real(realk) (A-H,O-Z)
!#include "mxcent.h"
!      Integer :: lupri
!      Type(opt_setting) :: optinfo
!      Real(realk) HESINT(MXRCRD,MXRCRD),EVEC(MX2CRD,MX2CRD)
!      Real(realk) TMPMAT(MX2CRD*MX2CRD),TMPMT2(MX2CRD*MX2CRD)
!      Real(realk) TMPMT3(MX2CRD*MX2CRD),TMPMT4(MX2CRD,MX2CRD)
!      Real(realk) TMPMT5(MX2CRD,MX2CRD)
!      Real(realk) VECMOD(MXCOOR)
!      PARAMETER (D0 = 0.0E0_realk, D1 = 1.0E0_realk, DP5 = 0.5E0_realk)
!      PARAMETER (IPRMIN = 1, IPRMED = 3, IPRMAX = 5, IPRDBG = 12)
!     For saddle point optimizations, we can follow a specific eigenvector.
!     Due to the fact that we are separating one mode for maximization,
!     NCRDHS is temporarily reduced by one.
!
!      IMODE = optinfo%NSPMod
!      NCRDHS = NCRDHS-1
!      IF (optinfo%NSPMod .GT. 0) THEN
!         call ls_FNDMOD(.FALSE.,MXRCRD,EVEC,TMPMAT,VECMOD, &
!     &        TMPMT2,TMPMT3,IMODE,lupri,optinfo)
!
!
!     If the lowest mode has a gradient element of zero, we have to pick
!     another mode for maximization (or we will end up in a minimum!).
!
!      ELSE
!         IMODE = 1
! 50      CONTINUE
!         IF ((ABS(optinfo%GRDDIA(IMODE)) .LT. 1.0E-10_realk) .AND. &
!     &        (IMODE .LT. NCRDHS)) THEN
!            IMODE = IMODE + 1
!            GOTO 50
!     
!     If we find no such mode, we just set IMODE = 1, because we must
!     be at a stationary point.
!     
!         ELSE IF (ABS(optinfo%GRDDIA(IMODE)) .LT. 1.0E-10_realk) THEN
!            IMODE = 1
!         END IF
!      END IF
!      IF (optinfo%IPrint .GE. IPRMAX) THEN
!         WRITE(LUPRI,*)
!         WRITE(LUPRI,*) 'Mode ',IMODE,' will be partitioned ' // &
!     &        'out and maximized.'
!         WRITE(LUPRI,*)
!      END IF!
!     The selected mode is placed at the very end.
!
!      call ls_DZERO(TMPMT5,MX2CRD*MX2CRD)
!      call ls_DZERO(TMPMT4,MX2CRD*MX2CRD)
!      DO 400 I = 1, NCRDHS
!         DO 402 J = 1, IMODE-1
!            TMPMT5(I,J) = EVEC(I,J)
! 402     CONTINUE
!         DO 403 J = IMODE, NCRDHS-1
!            TMPMT5(I,J) = EVEC(I,J+1)
! 403     CONTINUE
!         TMPMT4(1,I) = optinfo%EVAL(I)
!         TMPMT4(2,I) = optinfo%GRDDIA(I)
!         TMPMT5(I,NCRDHS) = EVEC(I,IMODE)
! 400  CONTINUE
!      TMPVAL = TMPMT4(1,IMODE)
!      DO 406 I = IMODE, NCRDHS-1
!         TMPMT4(1,I) = TMPMT4(1,I+1)
!         TMPMT4(2,I) = TMPMT4(2,I+1)
! 406  CONTINUE
!      TMPMT4(1,NCRDHS) = TMPVAL
!      TMPMT4(2,NCRDHS) = optinfo%GRDDIA(IMODE)
!
!     We then make the augmented Hessian that will be minimized.
!
!      call ls_DZERO(TMPMT2,NCRDHS*NCRDHS)
!      DO 410 I = 1, NCRDHS-1
!         TMPMT2(I+(I-1)*NCRDHS) = TMPMT4(1,I)
!         TMPMT2(I+(NCRDHS-1)*NCRDHS) = TMPMT4(2,I)
!         TMPMT2(NCRDHS+(I-1)*NCRDHS) = TMPMT4(2,I)
! 410  CONTINUE
!      TMPMT2(NCRDHS+(NCRDHS-1)*NCRDHS) = D0
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Augmented Hessian')
!         call output(TMPMT2,1,NCRDHS,1,NCRDHS,NCRDHS, &
!     &        NCRDHS,1,LUPRI)
!      END IF
!      call ls_DZERO(TMPMT3,MX2CRD*MX2CRD)
!      call ls_DSITSP(NCRDHS,TMPMT2,TMPMT3)
!      call ls_DUNIT(TMPMT2,NCRDHS)
!      call ls_JACO(TMPMT3,TMPMT2,NCRDHS,NCRDHS,NCRDHS, &
!     &     TMPMAT(1),TMPMAT(MX2CRD+1))
!      DO 420 J = 1, NCRDHS
!         optinfo%EVAL(J) = TMPMT3(J*(J+1)/2)
!         DO 425 I = 1, NCRDHS
!            EVEC(I,J) = TMPMT2(I+(J-1)*NCRDHS)
! 425     CONTINUE
! 420  CONTINUE
!
!     We add 1.0E5_realk to all eigenvalues that are essentially zero
!     for the sorting.
!
!      DO 427 I = 1, NCRDHS
!         IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-8_realk) optinfo%EVAL(I) = optinfo%EVAL(I) + 1.0E5_realk
! 427  CONTINUE
!      DO 430 I = 1, NCRDHS
!         JMIN = I
!         EMIN = optinfo%EVAL(I)
!         DO 435 J = (I + 1), NCRDHS
!            IF (optinfo%EVAL(J) .LT. EMIN) THEN
!               EMIN = optinfo%EVAL(J)
!               JMIN = J
!            END IF
! 435     CONTINUE
!         IF (JMIN .NE. I) THEN
!            call dswap(1,  optinfo%EVAL  (I),1,optinfo%EVAL  (JMIN),1)
!            call dswap(MX2CRD,EVEC(1,I),1,EVEC(1,JMIN),1)
!     call dswap(1,optinfo%GRDDIA(I),1,optinfo%GRDDIA(JMIN),1)
!         END IF
! 430  CONTINUE
!      DO 440 I = 1, NCRDHS
!         IF (ABS(ABS(optinfo%EVAL(I))-1.0E5_realk) .LT. 1.0E-3_realk) &
!     &        optinfo%EVAL(I) = optinfo%EVAL(I) - 1.0E5_realk
! 440  CONTINUE
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'RF-eigenvalues')
!         call output(optinfo%EVAL,1,1,1,NCRDHS,1,MXRCRD,1,LUPRI)
!         call lsheader(lupri,'RF-eigenvectors')
!         call output(EVEC,1,NCRDHS,1,NCRDHS,MX2CRD,MX2CRD,1,LUPRI)
!      END IF
!      call ls_PRFSTP(MX2CRD,NCRDHS,MXRCRD,EVEC,optinfo%STPINT,optinfo%GRDINT, &
!     &     TMPMAT,HESINT,1,lupri,optinfo)
!
!     In the case of saddle point optimization, we also need
!     to take care of the second partition and combine the two.
!
!      CMPLIM = MAX(optinfo%TrustRad*0.67E0_realk, 0.30E0_realk)
!      DO 500 I = 1, NCRDHS-1
!         IF (ABS(optinfo%STPINT(I)) .GT. CMPLIM) &
!     &        optinfo%STPINT(I) = SIGN(CMPLIM,optinfo%STPINT(I))
!         TMPMT4(3,I) = optinfo%STPINT(I)
! 500  CONTINUE
!      NCRDHS = NCRDHS + 1
!
!     We then make the augmented Hessian that will be maximized.
!
!      call ls_DZERO(TMPMT2,NCRDHS*NCRDHS)
!      TMPMT2(1) = TMPMT4(1,NCRDHS-1)
!      TMPMT2(2) = TMPMT4(2,NCRDHS-1)
!      TMPMT2(3) = TMPMT4(2,NCRDHS-1)
!      TMPMT2(4) = D0
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Augmented Hessian')
!         call output(TMPMT2,1,2,1,2,2,2,1,LUPRI)
!      END IF
!      call ls_DZERO(TMPMT3,MX2CRD*MX2CRD)
!      call ls_DSITSP(2,TMPMT2,TMPMT3)
!      call ls_DUNIT(TMPMT2,2)
!      call ls_JACO(TMPMT3,TMPMT2,2,2,2,TMPMAT(1),TMPMAT(1+MX2CRD))
!      DO 510 J = 1, 2
!         optinfo%EVAL(J) = TMPMT3(J*(J+1)/2)
!         DO 515 I = 1, 2
!            EVEC(I,J) = TMPMT2(I+(J-1)*2)
! 515     CONTINUE
! 510  CONTINUE
!      DO 517 I = 1, 2
!         IF (ABS(optinfo%EVAL(I)) .LE. 1.0E-8_realk) optinfo%EVAL(I) = optinfo%EVAL(I) + 1.0E5_realk
! 517  CONTINUE
!
!     The eigenvalues are sorted
!
!      IF (optinfo%EVAL(1) .GT. optinfo%EVAL(2)) THEN
!         call dswap(1,  optinfo%EVAL  (1),1,optinfo%EVAL  (2),1)
!         call dswap(MX2CRD,EVEC(1,1),1,EVEC(1,2),1)
!      END IF
!      DO 520 I = 1, 2
!         IF (ABS(ABS(optinfo%EVAL(I))-1.0E5_realk) .LT. 1.0E-3_realk) &
!     &        optinfo%EVAL(I) = optinfo%EVAL(I) - 1.0E5_realk
! 520  CONTINUE
!      call ls_PRFSTP(MX2CRD,2,optinfo%NIntCoord,EVEC,optinfo%STPINT,optinfo%GRDINT, &
!     &     TMPMAT,HESINT,2,lupri,optinfo)
!      TMPVAL = optinfo%STPINT(1)
!      IF (ABS(TMPVAL) .GT. CMPLIM) &
!     &     TMPVAL = SIGN(CMPLIM,TMPVAL)
!      call ls_DZERO(optinfo%STPSYM,MXCOOR)
!      DO 530 I = 1, NCRDHS-2
!         optinfo%STPSYM(I) = TMPMT4(3,I)
! 530  CONTINUE
!      optinfo%STPSYM(NCRDHS-1) = TMPVAL
!      IF (optinfo%IPrint .GE. IPRDBG) THEN
!         call lsheader(lupri,'Diagonal RF-step')
!         call output(optinfo%STPSYM,1,1,1,NCRDHS-1,1,MXRCRD,1,LUPRI)
!      END IF
!
!     The symmetry step is constructed from the original eigenvectors
!     (of the normal Hessian) and the diagonal RF-step.
!
!      call ls_DZERO(optinfo%STPINT,MXRCRD)
!      call ls_DZERO(EVEC,MX2CRD*MX2CRD)
!      DO 540 I = 1, NCRDHS-1
!         DO 545 J = 1, NCRDHS-1
!            EVEC(I,J) = TMPMT5(I,J)
! 545     CONTINUE
!         optinfo%EVAL(I) = TMPMT4(1,I)
! 540  CONTINUE
!      DO 550 I = 1, NCRDHS-1
!         DO 555 J = 1, NCRDHS-1
!            optinfo%STPINT(I) = optinfo%STPINT(I) + EVEC(I,J)*optinfo%STPSYM(J)
! 555     CONTINUE
! 550  CONTINUE
!
!     The final call ls_to RFSTP to do scaling with respect to the
!     trust radius.
!
!      call ls_PRFSTP(MX2CRD,NCRDHS,MXRCRD,EVEC,optinfo%STPINT,optinfo%GRDINT, &
!     &     TMPMAT,HESINT,3,lupri,optinfo)
!      RETURN
!      END
!
!  /* Deck ls_rhoij */
      FUNCTION rhoij_ls(ATMARR,I,J,ALFA,DIST,ORIG_LINDH)
!
!     Calculates the function
!        rho_ij = exp[alfa_ij (r_ref,ij^2 - r_ij^2)]
!     for the model Hessian by Roland Lindh
!
!     Revised Nov 2002 hjaaj for using covalent/metallic radii
!     instead of simplified values in DIST
!     (the problem is that they do NOT work for heavier
!     nuclei - as Roland also notes in his paper;
!     with covalent/metallic radii we should get reasonable
!     r_ref for all bonds)
!
!     for atom no. I:
!     ATMARR(I,1) = charge
!     ATMARR(I,2:4) = (x,y,z) in Bohr
!     ATMARR(I,5) = covalent/metallic radius in Angstrom
!
use precision
use Fundamental
use optimization_input
Implicit Real(realk) (A-H,O-Z)
      Real(realk) ATMARR(MXCENT,8)
      Real(realk) ALFA(3,3), DIST(3,3), VEC(3)
      LOGICAL   ORIG_LINDH
      IROWI = 2
      IROWJ = 2
      IF (ATMARR(I,1) .LE. 2) IROWI = 1
      IF (ATMARR(J,1) .LE. 2) IROWJ = 1
      IF (ATMARR(I,1) .GE. 11) IROWI = 3
      IF (ATMARR(J,1) .GE. 11) IROWJ = 3
      VEC(1) = ATMARR(I,2) - ATMARR(J,2)
      VEC(2) = ATMARR(I,3) - ATMARR(J,3)
      VEC(3) = ATMARR(I,4) - ATMARR(J,4)
      V2 = VEC(1)*VEC(1) + VEC(2)*VEC(2) + VEC(3)*VEC(3)

      IF (ORIG_LINDH) THEN
         RREFIJ = DIST(IROWI,IROWJ)
      ELSE
         RREFIJ = (ATMARR(I,5) + ATMARR(J,5))/bohr_to_angstrom
      END IF
      rhoij_ls = EXP( ALFA(IROWI,IROWJ)*(RREFIJ*RREFIJ - V2) )
!     write (lupri,*) 'i,j,rrefij,v2,ls_rhoij',i,j,rrefij,v2,ls_rhoij
      RETURN
      END

!  /* Deck bldhes */
      SUBROUTINE LS_BLDHES(Molecule,MXRCRD,ATMARR,HESINT,lupri,optinfo)
use precision
use ls_util 
use optimization_input 
use files
use molecule_type
use molecule_typetype
!
!     Builds a simple model Hessian in redundant internal coordinates.
!     As described by Roland Lindh.
!
Implicit Real(realk) (A-H,O-Z)
      Integer :: lupri
      Type(MOLECULEINFO) :: Molecule
      Type(opt_setting) :: optinfo
      Real(realk) ATMARR(MXCENT,8), HESINT(MXRCRD,MXRCRD)
!
!     Here we assign all the necessary parameters
!
      PARAMETER (STRET = 0.450E0_realk, ROTAT = 0.150E0_realk, TORSN = 0.005E0_realk)
      Real(realk) ALFA(3,3), DIST(3,3)
      SAVE      ALFA, DIST
      DATA ALFA / 1.0000E0_realk,  0.3949E0_realk,  0.3949E0_realk, &
     &            0.3949E0_realk,  0.2800E0_realk,  0.2800E0_realk, &
     &            0.3949E0_realk,  0.2800E0_realk,  0.2800E0_realk/
      DATA DIST / 1.35E0_realk, 2.10E0_realk, 2.53E0_realk, &
     &            2.10E0_realk, 2.87E0_realk, 3.40E0_realk, &
     &            2.53E0_realk, 3.40E0_realk, 3.40E0_realk/
!
      call ls_DZERO(HESINT,MXRCRD*MXRCRD)
      iATOM = optinfo%IcartCoord/3
      call Atom_Ini(ATMARR,Molecule,optinfo,IATOM,.TRUE.,lupri)
!
!     We loop over all the internal coordinates and build up the
!     diagonal Hessian.
!
      NRIC = optinfo%NIntCoord
      IF (optinfo%DelInt) NRIC = optinfo%NIntCoord
      DO 10 IC = 1, NRIC
        IF (optinfo%INTCRD(IC,1) .LT. 10) THEN
          FAC = STRET
          RIJ = rhoij_ls(ATMARR,optinfo%INTCRD(IC,2),optinfo%INTCRD(IC,3),ALFA,DIST,optinfo%LINDHD)
          RJK = 1.0E0_realk
          RKL = 1.0E0_realk
        ELSE IF (optinfo%INTCRD(IC,1) .LT. 20) THEN
          FAC = ROTAT
          RIJ = rhoij_ls(ATMARR,optinfo%INTCRD(IC,2),optinfo%INTCRD(IC,3),ALFA,DIST,optinfo%LINDHD)
          RJK = rhoij_ls(ATMARR,optinfo%INTCRD(IC,3),optinfo%INTCRD(IC,4),ALFA,DIST,optinfo%LINDHD)
          RKL = 1.0E0_realk
        ELSE
          FAC = TORSN
          RIJ = rhoij_ls(ATMARR,optinfo%INTCRD(IC,2),optinfo%INTCRD(IC,3),ALFA,DIST,optinfo%LINDHD)
          RJK = rhoij_ls(ATMARR,optinfo%INTCRD(IC,3),optinfo%INTCRD(IC,4),ALFA,DIST,optinfo%LINDHD)
          RKL = rhoij_ls(ATMARR,optinfo%INTCRD(IC,4),optinfo%INTCRD(IC,5),ALFA,DIST,optinfo%LINDHD)
        END IF
!
!     Feb 2009 - VB:
!     Attempt to scale interfragment bonds, no definite improvement, thus turned off
!
!        IF ((optinfo%INTCRD(IC,1) .EQ. 3) .OR.(optinfo%INTCRD(IC,1) .EQ. 4)) THEN
!           FAC = FAC*0.5E0_realk
!        END IF
        HESINT(IC,IC) = FAC*RIJ*RJK*RKL
 10   CONTINUE
      RETURN
      END

!  /* Deck nrmlvc */
      SUBROUTINE LS_NRMLVC(VEC)
use optimization_input
use precision
use ls_util 
use files
!
!     Normalizes a three-Real(realk)al vector.
!
Implicit Real(realk) (A-H,O-Z)
      Real(realk) VEC(3)
      VECNRM = SQRT(DDOT(3,VEC,1,VEC,1))
      VEC(1) = VEC(1)/VECNRM
      VEC(2) = VEC(2)/VECNRM
      VEC(3) = VEC(3)/VECNRM
      RETURN
      END

!  /* Deck nrmlvx */
      SUBROUTINE LS_NRMLVX(ICRD,VEC)
use optimization_input
use ls_util 
use files
use precision
!
!     Normalizes any vector.
!
Implicit Real(realk) (A-H,O-Z)
      Real(realk) VEC(ICRD)
      VECNRM = SQRT(DDOT(ICRD,VEC,1,VEC,1))
      IF (VECNRM .LE. 1.0E-12_realk) VECNRM = 1.0E0_realk
      DO 10 I = 1, ICRD
         VEC(I) = VEC(I)/VECNRM
 10   CONTINUE
      RETURN
      END

!  /* Deck vecprd */
      SUBROUTINE LS_VECPRD(VEC1,VEC2,VPRD)
use ls_util 
use files
use precision
!
!     Finds the vector product of two three-Real(realk)al vectors.
!
Implicit Real(realk) (A-H,O-Z)
      Real(realk) VEC1(3),VEC2(3),VPRD(3)
      VPRD(1) =  VEC1(2)*VEC2(3) - VEC2(2)*VEC1(3)
      VPRD(2) = -VEC1(1)*VEC2(3) + VEC2(1)*VEC1(3)
      VPRD(3) =  VEC1(1)*VEC2(2) - VEC2(1)*VEC1(2)
      RETURN
      END
!==============!
! vecang_ls    !
!==============!
      FUNCTION vecang_ls(V1, V2)
! :::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::                                                 ::
! ::       Finds the angle between two vectors       ::
! ::                                                 ::
! :::::::::::::::::::::::::::::::::::::::::::::::::::::
Use Fundamental
Use precision
Implicit Real(realk) (A-H,O-Z)
!
Real(realk) V1(3), V2(3)
Real(realk) :: ZerTol
!
ZerTol = 1.0E-12_realk
TEMP = (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))*&
     &           (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
         IF (TEMP .GT. ZERTOL) THEN
            TEMP = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))/SQRT(TEMP)
! ::: Parallell vectors might yield a value here that is slightly greater :::
! ::: than one. ACOS is undefined for these values, so we have to round  :::
! ::: them off, removing the "excess", accumulated numerical error.      :::
! ::: Also, ACOS is extremely sensitive around +/-1, so we round off     :::
! ::: numbers close to these values
            IF ((ABS(ABS(TEMP)-1.0E0_realk) .LT. ZERTOL) &
     &           .OR. (ABS(TEMP).GT. 1.0E0_realk)) THEN
               TEMP = ANINT(TEMP)
            END IF
            vecang_ls = ACOS(TEMP)
         ELSE
            vecang_ls = 0.5E0_realk*PI
         END IF
         RETURN
END FUNCTION vecang_ls
!================!
!   vdwrad_ls    !
!================!
      FUNCTION vdwrad_ls(NCHARGE,lupri)
use precision
use ls_util
Implicit Real(realk) (A-H,O-Z)
!     Based on van der Waals radii in Angstrom.
!     Returns -1 where data is inavailable
      Real(realk) RAD(100)
      DATA (RAD(I), I = 1, 100)/&
     &       110.,  220.,  122.,   63.,&
     &155.,  155.,  140.,  135.,  130.,&
     &154.,  190.,  160.,  140.,  110.,&
     &202.,  220.,  150.,  150.,  220.,&
     &188.,  181.,  175.,  277.,  239.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100.,   -100.,   -100.,   -100.,   -100.,&
     & -100./
!
 IF (NCHARGE .EQ. 0) THEN
!hj  solvent cavity center or floating orbital
     vdwrad_ls = -1.0E0_realk
 ELSE IF (NCHARGE .LT. 1 .OR. NCHARGE .GT. 100) THEN
    WRITE (LUPRI,*) 'ERROR, vdwrad_ls called with NCHARGE =',NCHARGE
    CALL LSQUIT('vdwrad_ls called with illegal NCHARGE',lupri)
 ELSE
    vdwrad_ls = 0.01E0_realk * RAD(NCHARGE)
 END IF
 RETURN
END FUNCTION vdwrad_ls
!=============!
!  dv3dot_ls  !
!=============! 
      FUNCTION dv3dot_ls(N,V1,V2,V3)
!
!     7-Aug-1986 hjaaj
!
use precision
Implicit Real(realk) (A-H,O-Z)
      Real(realk) V1(*), V2(*), V3(*)
      Real(realk), PARAMETER :: D0 = 0.0E0_realk 
      T = D0
      DO 100 K = 1,N
         T = T + V1(K) * V2(K) * V3(K)
  100 CONTINUE
      dv3dot_ls = T
      RETURN
END FUNCTION dv3dot_ls
!==========!
!  DGEINV  !
!==========!
SUBROUTINE DGEINV(N,A,AINV,WRK,INFO)
Use precision
      INTEGER N, N2
      Real(realk)  A(*), AINV(*), WRK(*), DET(2)
      INTEGER IPVT(N)
!
      N2 = N*N
      CALL DCOPY(N2,A,1,AINV,1)
      CALL DGEFA(AINV,N,N,IPVT,INFO)
      IF (INFO .EQ. 0) CALL DGEDI(AINV,N,N,IPVT,DET,WRK,01)
      RETURN
END SUBROUTINE
