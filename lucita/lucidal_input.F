!***********************************************************************
      subroutine LUCITAINP(word,inperr,allopt)
!***********************************************************************
C
C       Input section for LUCITA module
C
C       Written by Timo Fleig for Dirac, December 2001
C       adapted for Dalton by Stefan Knecht, November 2010
C
C***********************************************************************
      use lucita_orbital_spaces
#include <implicit.h>
#include "priunit.h"
C
#include "luctinp.h"
#include "inforb.h"
C
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
      LOGICAL   ALLOPT,NEWDEF
      CHARACTER WORD*7,WORD1*7,PROMPT*1,TABLE(NTABLE)*7
!
!     temporary scratch arrays for input processing
      integer, allocatable :: tmp_inactive_shells(:)
      integer, allocatable :: tmp_active_shells(:,:)
      integer              :: max_number_of_ptg_irreps
      integer              :: max_number_of_gas_spaces
C
C
      DATA TABLE /'.TITLE ','.INIWFC','.CITYPE','.NROOTS','.SYMMET',
     &            '.NACTEL','.MULTIP','.PRINTG','.PRINTL','.SZCALC',
     &            '.INACTI','.GASSHE','.GASSPC','.FROORB','.RAS1  ',
     &            '.RAS2  ','.RAS3  ','.DENSI ','.RSTRCI','.MXCIVE',
     &            'xxxxxxx','xxxxxxx','.MAXITR','.LBLKSZ','xxxxxxx',
     &            '.DISTRT','xxxxxxx','.TRUNCF','.MEMFAC','XXXXXXX'/

      if(allopt) call prtab(ntable,table,word//' input keywords',
     &                      lupri)

!     Initialize keyword status flags driving interface.
!     (Therefore local keywords do not have to be initialized!)
      call izero(imokw,ntable)

!     allocate scratch space for input processing
!     note: we assume that
!           - we run in D2h or a subgroup of it, e.g. no super symmetry.
!           - the number of GAS spaces does not exceed 6 (see MXPNGAS in lucita/mxpdim.inc)
      if(nsym.gt.8) call quit(' *** error in LUCITAINP: only D2h 
     & and its sub-point groups are implemented. ***')

      max_number_of_ptg_irreps = 8
      max_number_of_gas_spaces = 6

      allocate(tmp_inactive_shells(max_number_of_ptg_irreps))
      allocate(tmp_active_shells(max_number_of_gas_spaces,
     &         max_number_of_ptg_irreps))

      call izero(tmp_inactive_shells,8)
      call izero(tmp_active_shells,6*8)

!     initialize common block luctinfi in luctinp.h
      NROOTD           = 0
      ISSYMD           = 0
      NACTED           = 0
      IMULTD           = 0
      IPRNGD           = 0
      IPRNLD           = 0
      IDENSD           = 0
      MXHL1D           = 0
      MXEL3D           = 0
      INGASD           = 0
      NSEQCD           = 0
      IRSTLT           = 0
      MXCIVE           = 0
      IPARMODEL        = 0
      ICIMAXITER       = 0
      IMAXBLKSIZE      = 0
      I_USE_DIST_ROUTE = 0
      IN_MEMFAC        = 0

!     Process input for LUCITA
!     ========================
      NEWDEF = (WORD .EQ. '*LUCITA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &            '" not recognized in LUCITA.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LUCITAINP.')
C
C&&& TITLE : One title line for CI specification
    1          CONTINUE
                  READ(LUCMD,'(A)') TITLUC
                  IMOKW(1) = 1
               GO TO 100
    2          CONTINUE
C&&& INIWFC: Initial DHF wave function (closed- or open-shell)
                  READ(LUCMD,'(A)') WAFFCD
                  IMOKW(2) = 1
               GO TO 100
    3          CONTINUE
C&&& CITYPE: Type of CI calculation (FCI, SD, SDTQ, RAS, GAS)
                  READ(LUCMD,'(A)') CALCTP
                  IMOKW(3) = 1
               GO TO 100
    4          CONTINUE
C&&& NROOTS: Number of eigenvalues
                  READ(LUCMD,*) NROOTD
                  IMOKW(4) = 1
               GO TO 100
    5          CONTINUE
C&&& SYMMET: State symmetry
                  READ(LUCMD,*) ISSYMD
                  IMOKW(5) = 1
               GO TO 100
    6          CONTINUE
C&&& NACTEL: Number of active electrons
                  READ(LUCMD,*) NACTED
                  IMOKW(6) = 1
               GO TO 100
    7          CONTINUE
C&&& MULTIP: State spin multiplicity
                  READ(LUCMD,*) IMULTD
                  IMOKW(7) = 1
               GO TO 100
    8          CONTINUE
C&&& PRINTG: Global LUCIA print flag
                  READ(LUCMD,*) IPRNGD
                  IMOKW(8) = 1
               GO TO 100
    9          CONTINUE
C&&& PRINTL: Local print flag DIRLUC
                  READ(LUCMD,*) IPRNLD
                  IMOKW(9) = 1
               GO TO 100
   10          CONTINUE
C&&& SZCALC: Rough size of calculation
                  READ(LUCMD,'(A)') SZCALD
                  IMOKW(10) = 1
               GO TO 100
   11          CONTINUE
C&&& INACTI: Inactive function space
                  read(lucmd,*) (tmp_inactive_shells(i), i=1,nsym)
                  IMOKW(11) = 1
               GO TO 100
   12          CONTINUE
C&&& GASSHE: Number of GAS and distribution
                  READ(LUCMD,*) INGASD
                  if(INGASD .gt. max_number_of_gas_spaces)
     &            call quit(' *** error in LUCITAINP:  number of GAS
     & spaces exceeds the current (efficiency) limit of 6. ***')
                  do i = 1,ingasd
                    read(lucmd,*) (tmp_active_shells(i,j), j=1,nsym)
                  end do
                  IMOKW(12) = 1
               GO TO 100
   13          CONTINUE
C&&& GASSPC: Cumulated max. and min. GAS occupations
                  READ(LUCMD,*) NSEQCD
                  DO IR = 1,INGASD,1
                    READ(LUCMD,'(A)') CRDGOC(IR)
                  END DO
                  IMOKW(13) = 1
               GO TO 100
   14          CONTINUE
C&&& FROORB: Frozen orbital function space
                  READ(LUCMD,'(A)') CRDFRO
                  write(6,'(6A,1A72)') 'CRDFRO',CRDFRO
                  IMOKW(14) = 1
               GO TO 100
   15          CONTINUE
C&&& RAS1  : RAS1 function space
                  READ(LUCMD,*) (tmp_active_shells(1,j), j=1,nsym)
                  READ(LUCMD,*) MXHL1D
                  IMOKW(15) = 1
               GO TO 100
   16          CONTINUE
C&&& RAS2  : RAS2 function space
                  READ(LUCMD,*) (tmp_active_shells(2,j), j=1,nsym)
                  IMOKW(16) = 1
               GO TO 100
   17          CONTINUE
C&&& RAS3  : RAS3 function space
                  READ(LUCMD,*) (tmp_active_shells(3,j), j=1,nsym)
                  READ(LUCMD,*) MXEL3D
                  IMOKW(17) = 1
               GO TO 100
   18          CONTINUE
C&&& DENSI : Level of density matrix calculation (1-,1- and 2-particle)
                  READ(LUCMD,*) IDENSD
                  IMOKW(18) = 1
               GO TO 100
   19          CONTINUE
C&&& RSTRCI: Restart option from CI vector
                  READ(LUCMD,*) IRSTLT
                  IMOKW(19) = 1
               GO TO 100
   20          CONTINUE
C&&& MXCIVE: Maximum number of subspace CI vectors
                  READ(LUCMD,*) MXCIVE
                  IMOKW(20) = 1
               GO TO 100
   21          CONTINUE
C&&& xxxxxxxxx:
               GO TO 100
   22          CONTINUE
C&&& xxxxxxxxx:
               GO TO 100
   23          CONTINUE
C&&& ICIMAXITER: Maximum number of CI Iterations
                  READ(LUCMD,*) ICIMAXITER
                  IMOKW(23) = 1
               GO TO 100
   24          CONTINUE
C&&& IMAXBLKSIZE: Maximum batchsize for c-/sigma-vector
                  READ(LUCMD,*) IMAXBLKSIZE
                  IMOKW(24) = 1
               GO TO 100
   25          CONTINUE
C&&& xxxxxxxxxx:
               GO TO 100
   26          CONTINUE
C&&& I_USE_DIST_ROUTE: Determination of distribution routine to use
                  READ(LUCMD,*) I_USE_DIST_ROUTE
                  IMOKW(26) = 1
               GO TO 100
   27          CONTINUE
C&&& xxxxxxxxxx:
               GO TO 100
   28          CONTINUE
C&&& TRUNCFactor: truncation factor for residual vectors for next trial
C                 vector
                  READ(LUCMD,*) CTRUNC_FAC
                  IMOKW(28) = 1
               GO TO 100
   29          CONTINUE
C&&& MEMFAC     : memory multiplier                                    
                  READ(LUCMD,*) IN_MEMFAC
                  IMOKW(29) = 1
               GO TO 100
   30          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in LUCITAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in LUCITAINP.')
            END IF
      END IF
  300 CONTINUE

!     lucita input has been read - proceed with compulsory settings
      if(ichang.gt.0)then

!       a. orbital spaces

!       init_input_type:
!       1: ci input
!       2: mcscf input
        init_input_type = 1

!       init_wave_f_type:
!       1: GAS input
!       2: RAS input
        init_wave_f_type = 1
        if(ingasd.eq.0) init_wave_f_type = 2

        call define_lucita_orb_spaces(nsym,
     &                                ingasd,
     &                                max_number_of_ptg_irreps,
     &                                max_number_of_gas_spaces,
     &                                tmp_inactive_shells,
     &                                tmp_active_shells,        
     &                                norb,
     &                                init_input_type,
     &                                init_wave_f_type)

!       b. more to come...???

      end if

!     reset inactive space keyword for GAS calculations (not in use there)
      if((ingasd .gt. 0) .and. (imokw(11) .eq. 1)) imokw(11) = 0

      deallocate(tmp_inactive_shells)
      deallocate(tmp_active_shells)

      end
