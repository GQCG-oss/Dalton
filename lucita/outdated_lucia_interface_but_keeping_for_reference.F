      subroutine readin_internal_lucia_input(luin,lucmol)
*
*
* File is supposed to be positioned at first record of input
* The end of the input stream is identified by END OF INPUT
* Unless MOLCS is specified,
* All keywords are initiated by a point ., while comments are
* initiated by a *.
*
* The keywords can broadly be divided into two types
*  1 : Keywords describing CI calculation to be carried out
*  2 : Keywords describing how CI optimization should be performed
*
*
* All input parameter concerning CI space are saved in /LUCIN1/
* All input concerning actual CI vectors are save in /CSTATE/
* All input paramters concerning run are saved in /CRUN/
*
* Since the keywords are read in from one pass over input file,
* the keywords must be in logical order.For example, the number
* of irreducible representations (irreps) must be give before
* the number of shells per irrep
*
* Jeppe Olsen,Spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      PARAMETER(MXPLNC = 72 )
      CHARACTER*72 TITLEC
      CHARACTER*72 CARD
      CHARACTER*72 CARD1
      CHARACTER*72 LASTCARD
      COMMON/CTITLE/TITLEC(3)
      CHARACTER*6 KEYWOR
      PARAMETER(MXPKW = 97)
      DIMENSION KEYWOR(MXPKW)
      DIMENSION ISETKW(MXPKW)
*. Local scratch for decoding multi-item lines, atmost 32 items per line
      PARAMETER(MX_ITEM = 32)
      CHARACTER*72 ITEM(MX_ITEM), CARDX, ITEMX
*
      DATA KEYWOR/'TITLE ','PNTGRP','NIRREP','INTSPC','EXTSPC',
     &     'NACTEL','INACT ','CORE  ','RAS1  ','RAS2  ',
     &     'RAS3  ','MXSCTP','SECOND','REFSPC','INTSEL',
     &     'MS2   ','MULTS ','IREFSM','ROOTS ','IDIAG ',
     &     'MAXIT ','EXPHAM','RESTRT','INTIMP','INCORE',
     &     'DELETE','MSCOMB','MLCOMB','IPRSTR','IPRCIX',
     &     'IPRORB','IPRDIA','MXCIV ','CISTOR','NOCSF ',
     &     'IPRXT ','NOINT ','DMPINT','RESDIM','CJKAIB',
     &     'INIREF','RESTRF','IPROCC','MOCAA ','MOCAB ',
     &     'ECORE ','PERTU ','APRREF','APRZER','GASSH ',
     &     'GASSPC','CMBSPC','CICONV','SEQUEN','EXTKOP',
     &     'MACHIN','C1DSC ','H0SPC ','H0FRM ','RFROOT',
     &     'H0EX  ','INIDEG','LAMBDA','LCSBLK','IPRDEN',
     &     'xxxxxx','ECHO  ','FINORB','E_THRE','C_THRE',
     &     'E_CONV','C_CONV','CLSSEL','DENSI ','PTEKT ',
     &     'H0ROOT','NORST2','SKIPEI','XYZSYM','PROPER',
     &     'TRAPRP','RESPON','MXITLE','IPRRSP','RTHOME',
     &     'USE_PH','ADVICE','TRACI ','USE_PS','PTFOCK',
     &     'PRNCIV','RES_CC','xxxxxx','xxxxxx','xxxxxx',
     &     'RELAX ','EXPERT'/

*
      INTEGER CITYP
*.Largest allowed number of allowed irreps for orbs

#include "mxpdim.inc"
#include "lucinp.inc"
#include "cstate.inc"
#include "crun.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "oper.inc"
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*./CECORE/
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
      CHARACTER IN_MPINN*7
      CHARACTER*11 NNIN_MPI
      INTEGER IN_MPILL
*
*. Flag for compatibility with normal MOLCAS input format
      MOLCS = 1
*
      NERROR = 0
      NWARN  = 0
      EXTSPC = 0
      IECHO  = 0
!     IECHO  = 1 ! hjaaj DEBUG , default 0
* No cc as default
      I_DO_CC = 0
*
*  Open file LUCIAIN from MOLUC conversion program
*
      IN_MPINN ="LUCIAIN"
      IF (LUCI_MYPROC .LT. 10) THEN    ! MPI ID has one digit
         WRITE (NNIN_MPI,'(A7,A1,I1)') IN_MPINN,'.',LUCI_MYPROC
         IN_MPILL=9
      ELSEIF (LUCI_MYPROC .LT. 100) THEN  ! MPI ID has two digits
         WRITE (NNIN_MPI,'(A7,A1,I2)') IN_MPINN,'.',LUCI_MYPROC
         IN_MPILL=10
      ELSE                                ! MPI ID has three digits
         WRITE (NNIN_MPI,'(A7,A1,I3)') IN_MPINN,'.',LUCI_MYPROC
         IN_MPILL=11
      ENDIF

      open(unit=LUCMOL,file=NNIN_MPI(1:IN_MPILL),
     &     status='OLD',form='FORMATTED')
*
*****************************************************************
*                                                               *
* Part 1 : Read in Keywords and perform some preliminary checks *
*                                                               *
*****************************************************************
*
*. Defaults for pointgroup and number of irreps must be set here
*. Default point group D2H
      PNTGRP = 1
*. Default number of irreps
      NIRREP = 8
*
      CALL IZERO(ISETKW,MXPKW)
 1000 CONTINUE
*. Next potential keyword
      READ(LUCMOL,'(A)') CARD
*. Left-position nonblank characters in CARD
      CALL LFTPOS(CARD,MXPLNC)
*. Change to upper case
C     UPPCAS(LINE,LENGTH)
      CALL UPPCAS(CARD,MXPLNC)
      IF(CARD(1:1).EQ.'*'.OR.CARD(1:1).EQ.'!'.OR.
     &     CARD(1:1).EQ.'#'                        )THEN
*. Skip comment cards
         GOTO 999
*. End of input card
      ELSE IF(CARD(1:5).EQ.'ENDOF') THEN
         GOTO 1001
      ELSE IF(MOLCS.EQ.0.AND.CARD(1:1).NE.'.') THEN
*. Line out of context
         WRITE(luwrt,'(1X ,A)') ' Warning, card out of context : '
         WRITE(luwrt,'(1X ,A)') CARD
         NWARN = NWARN + 1
      ELSE IF(MOLCS.EQ.1.OR.CARD(1:1).EQ.'.') THEN
         IF(MOLCS.EQ.1) THEN
*. Move characters one place to right
            DO 1286 ICHAR = 7,2,-1
               CARD(ICHAR:ICHAR) = CARD(ICHAR-1:ICHAR-1)
 1286       CONTINUE
            CARD(1:1) = ' '
         END IF
*. A keyword has been identified, match with possible keywords
         IF(CARD(2:6).EQ.'TITLE' ) THEN
*
* =========================
*.Keyword 1 :  TITLE cards
* =========================
*
*. Three title cards
            ISETKW(1) = 1
            DO 20 IC = 1, 1
               READ(LUCMOL,'(A)') TITLEC(IC)
 20         CONTINUE
            GOTO 999
         END IF
*
*
*================================================
*. Keyword 2 : <POINTG> : Point group of orbitals
*================================================
*
* Possible point groups : D2H,CINFV,DINFH,O3
         IF(CARD(2:4).EQ.'D2H'   .OR.
     &        CARD(2:6).EQ.'CINFV' .OR.
     &        CARD(2:6).EQ.'DINFH' .OR.
     &        CARD(2:3).EQ.'O3'    ) THEN
*
            ISETKW(2) = 1
            IF(CARD(2:4).EQ.'D2H') THEN
               PNTGRP = 1
            ELSE IF(CARD(2:6).EQ.'CINFV') THEN
               PNTGRP = 2
            ELSE IF(CARD(2:6).EQ.'DINFH') THEN
               PNTGRP = 3
            ELSE IF(CARD(2:3).EQ.'O3') THEN
               PNTGRP = 4
            END IF
            GOTO 999
         END IF
*
         IF(CARD(2:7).EQ.'NIRREP') THEN
*
*=====================================================
*. Keyword 3 : <NIRREP> : Number of irreps of orbitals
*=====================================================
*
* Number of irreducible representations in point group
* D2h              : 1,2,4,8
* C inf H, D inf H : largest ML
* O3               : Largest L
*.D2h or subgroup
* ===============
            IF(PNTGRP.EQ.1) THEN
               READ(LUCMOL,*) NIRREP
               NSMCMP = NIRREP
               NSMOB  = NIRREP
               ISETKW(3) = 1
*.Dimensions 3,5,6,7,8 are not allowed
               IF(NIRREP.EQ.3.OR.(NIRREP.GT.4.AND.NIRREP.LT.8)) THEN
                  WRITE(luwrt,*) ' Input error : NIRREP = ', NIRREP
                  WRITE(luwrt,*) ' Allowed values of NIRREP :1,2,4,8'
                  NERROR = NERROR + 1
                  ISETKW(3) = -1
               END IF
*. Zero values used for other pointgroups
               MAXML  = -1
               MAXL   = -1
               INVCNT = -1
            ELSE IF (PNTGRP.EQ.2) THEN
*. Cinf V
* =======
               READ(LUCMOL,*) MAXML
               ISETKW (3) = 1
               IF(MAXML.LT.0) THEN
                  WRITE(luwrt,*)
     &           ' Largest ML values of shells must be atleast be zero '
                  WRITE(luwrt,*) ' MAXML from input :' ,MAXML
                  NERROR = NERROR + 1
                  ISETKW(3) = -1
               END IF
               NIRREP =  MAXML + 1
               NSMCMP = 2 * MAXML + 1
               NSMOB = NSMCMP
               INVCNT = 0
               MAXL = -1
            ELSE IF (PNTGRP.EQ.3) THEN
*. Dinf H
* =======
               READ(LUCMOL,*) MAXML
               ISETKW (3) = 1
               IF(MAXML.LT.0) THEN
                  WRITE(luwrt,*)
     &           ' Largest ML values of shells must at least be zero '
                  WRITE(luwrt,*) ' MAXML from input : ',MAXML
                  NERROR = NERROR + 1
                  ISETKW(3) = -1
               END IF
               NIRREP = 2 * ( MAXML + 1)
               NSMCMP = 2 * ( 2*MAXML + 1 )
               NSMOB = NSMCMP
               INVCNT = 1
               MAXL  = -1
            ELSE IF (PNTGRP.EQ.4) THEN
*. O 3
* =======
               READ(LUCMOL,*) MAXL
               ISETKW (3) = 1
               IF(MAXL.LT.0) THEN
                  WRITE(luwrt,*)
     &            ' Largest L values of shells must be atleast be zero '
                  WRITE(luwrt,*) ' MAXL from input : ' , MAXL
                  NERROR = NERROR + 1
                  ISETKW(3) = -1
               END IF
               MAXML = MAXL
               NIRREP = MAXL + 1
               NSMCOM = 2 * (2 *MAXML + 1 )
               NSMOB = NSMCMP
               INVCNT = 1
            END IF
            IF(ISETKW(3).EQ.-1)
     &           WRITE(luwrt,*) ' .NIRREP input incorrect !! . '
            GOTO 999
         END IF
*
* ================================================
*. Keyword 4 : INTSPC : Type of internal CI space
* ================================================
*
         IF(CARD(2:4).EQ.'CAS'.OR.CARD(2:4).EQ.'FCI'.OR.
     &        (CARD(2:4).EQ.'RAS'.AND.CARD(5:5).EQ.' ')) THEN
            ISETKW(4) = 1
            IF(CARD(2:4).EQ.'CAS'.OR. CARD(2:4).EQ.'FCI' ) THEN
               INTSPC = 1
            ELSE IF (CARD(2:4).EQ.'RAS') THEN
               INTSPC = 2
*. Limits on allowed number of electrons in RASI and RAS III
               READ(LUCMOL,*) MNRS10,MXRS30
            END IF
            GOTO 999
         END IF
*
*===========================
* Keyword 5 : External space
*===========================
*
         IF(CARD(2:7).EQ.'EXTSPC') THEN
            ISETKW(5) = 1
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
            IF(CARD1(1:4).EQ.'NONE') THEN
               EXTSPC = 0
               MXER4 = 0
               MXHR0 = 0
            ELSE IF
     &          (CARD1(1:4).EQ.'CORE'.AND.CARD1(5:10).EQ.'SECOND') THEN
               EXTSPC = 3
               READ(LUCMOL,*) MXHR0,MXER4
            ELSE IF
     &          (CARD1(1:6).EQ.'SECOND'.AND.CARD1(7:10).EQ.'CORE') THEN
               EXTSPC = 3
               READ(LUCMOL,*) MXER4,MXHR0
            ELSE IF(CARD1(1:4).EQ.'CORE') THEN
               EXTSPC = 1
               READ(LUCMOL,*) MXHR0
            ELSE IF(CARD1(1:6).EQ.'SECOND') THEN
               EXTSPC = 2
               READ(LUCMOL,*) MXER4
            ELSE
               ISETKW(5) = - 1
               WRITE(luwrt,*) ' Illegal card for EXTSPC : '
               WRITE(luwrt,'(1X ,A)') CARD1
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
*
* =============================================
* Keyword 6 NACTEL : Number of active electrons
* =============================================
*
         IF(CARD(2:7).EQ.'NACTEL') THEN
            READ(LUCMOL,*)LUCI_NACTEL
            ISETKW(6) = 1
            IF(LUCI_NACTEL.LT.0) THEN
               WRITE(luwrt,*)
     &           ' ERROR : Illegal number of active electrons ', 
     &                     LUCI_NACTEL
               ISETKW(6) = -1
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
*==================
* 7 : Inactive shells
*==================
         IF(CARD(2:7).EQ.'INACTI') THEN
            READ(LUCMOL,*) (NINASH(IRREP),IRREP = 1, NIRREP)
            ISETKW(7) = 1
            GOTO 999
         END IF
*=================================
* 8 : Core shells ( = RAS0 shells)
*==================================
         IF(CARD(2:5).EQ.'CORE') THEN
            READ(LUCMOL,*) (NRS0SH(1,IRREP),IRREP = 1, NIRREP)
            ISETKW(8) = 1
            EXTSPC = EXTSPC + 1
            GOTO 999
         END IF
*===========
* 9 : RAS 1
*===========
         IF(CARD(2:5).EQ.'RAS1') THEN
*.Number of RAS 1 shells per irrep
            READ(LUCMOL,*) (NRSSH(IRREP,1),IRREP = 1, NIRREP)
*.Smallest allowed number of electrons in RAS 1
C!         READ(LUCMOL,*) MNER10
            ISETKW(9) = 1
            GOTO 999
         END IF
*===========
* 10 : RAS 2
*===========
         IF(CARD(2:5).EQ.'RAS2'.OR.CARD(2:7).EQ.'ACTIVE') THEN
            READ(LUCMOL,*) (NRSSH(IRREP,2),IRREP = 1, NIRREP)
            ISETKW(10) = 1
            GOTO 999
         END IF
*===========
* 11 : RAS 3
*===========
         IF(CARD(2:5).EQ.'RAS3') THEN
            ISETKW(11) = 1
*.Number of RAS 3 shells per irrep
            READ(LUCMOL,*) (NRSSH(IRREP,3),IRREP = 1, NIRREP)
*.Largest allowed number of electrons in RAS III
C!         READ(LUCMOL,*) MXER30
            GOTO 999
         END IF
* ==================================================
* 12 : RAS 4, subdivision of shells of given irrep
* ==================================================
*. Transferred to read in of RAS1
* =================================================
* 13 : Number of shells in secondary space per type
* =================================================
         IF(CARD(2:7).EQ.'SECOND') THEN
*. Number of subdivisions of RAS4
            READ(LUCMOL,*) MXR4TP
*. Number of secondary irreps per type
            DO 30 ITYPE = 1, MXR4TP
              READ(LUCMOL,*) (NRS4SH(IRREP,ITYPE),IRREP = 1, NIRREP)
 30         CONTINUE
            ISETKW(13) = 1
*.Largest allowed number of electron in secondary space
            GOTO 999
         END IF
* =========================
* 14 : Reference space
* =========================
*
* For RAS reference spaces,the reference space can be reduced compared
* to the zero order space by introducing different constraints on
* the number of allowed electrons in RAS1 and RAS3. The reference space
* must however be a subset of the zero order space.
*
         IF(CARD(2:7).EQ.'REFSPC') THEN
            READ(LUCMOL,*) MNRS1R,MXRS3R
            ISETKW(14) = 1
            IF(MNRS1R.LT.MNRS10.OR.MXRS3R.GT.MXRS30) THEN
               WRITE(luwrt,*)
     &              ' Reference space larger than zero order space'
               NERROR = NERROR + 1
               ISETKW(14) = -1
            END IF
            GOTO 999
         END IF
*
* =========================================================
* 15 : selection of active configurations in internal space
* =========================================================
*
         IF(CARD(2:7).EQ.'INTSEL' ) THEN
            ISETKW(15) = 1
*
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:4).EQ.'NONE') THEN
*. All internals are included
               INTSEL = 0
            ELSE IF(CARD1(1:6).EQ.'INDTST') THEN
*. Include coeffcients larger than CTHRES or having energy contributions
*. larger than ETHREA
               INTSEL = 1
               READ(LUCMOL,*)  CTHRES,ETHRES
            ELSE IF(CARD1(1:6).EQ.'TOTTST') THEN
*. Obtain CTHRES of the total wavefinction and ETHRES of the total
*. energy
               INTSEL = 2
               READ(LUCMOL,*)  CTHRES,ETHRES
            ELSE IF(CARD1(1:6).EQ.'INDWCN') THEN
*. Include all configutations with reference weights larger than
*. a given threshold in reference CI
               INTSEL = 3
               READ(LUCMOL,*) XWCNF
            ELSE IF(CARD1(1:6).EQ.'TOTWCN') THEN
*. Include the largest configurations so all a given fraction
*. of the Zero order reference is included in the CI
               INTSEL = 4
               READ(LUCMOL,*) XWCNF
            ELSE
               WRITE(luwrt,*) ' Error : Illegal card in INTSEL :'
               WRITE(luwrt,'(1X ,A)') CARD1
               ISETKW(15) = - 1
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
*===============================
* 16 : Two times spin projection
* ==============================
         IF(CARD(2:4).EQ.'MS2') THEN
            ISETKW(16) = 1
            READ(LUCMOL,*) MS2
            GOTO 999
         END IF
*========================
* 17 : spin multiplicity
* =======================
         IF(CARD(2:6).EQ.'MULTS') THEN
            ISETKW(17) = 1
            READ(LUCMOL,*) MULTS
            GOTO 999
         END IF
*========================
* 18 : Reference symmetry
* =======================
         IF(CARD(2:7).EQ.'IREFSM') THEN
            ISETKW(18) = 1
            IF(PNTGRP.EQ.1) THEN
               READ(LUCMOL,*) IREFSM
            ELSE IF(PNTGRP.EQ.2) THEN
               READ(LUCMOL,*) IREFML
            ELSE IF (PNTGRP.EQ.3) THEN
               READ(LUCMOL,*) IREFML,IREFPA
               IF(IREFPA.EQ.-1) IREFPA = 2
            ELSE IF (PNTGRP.EQ.3) THEN
               READ(LUCMOL,*) IREFL,IREFML,IREFPA
               IF(IREFPA.EQ.-1) IREFPA = 2
            END IF
            GOTO 999
         END IF
*==========================
* 19 : Roots to be obtained
* =========================
         IF(CARD(2:6).EQ.'ROOTS') THEN
            ISETKW(19) = 1
            READ(LUCMOL,*) NROOT
            DO I = 1, NROOT
               IROOT(I) = I
            END DO
C           READ(LUCMOL,*) (IROOT(I),I=1,NROOT)
            GOTO 999
         END IF
*===============================
* 20 : Diagonalization algorithm   : .MEGACI , .TERACI
*===============================
         IF(CARD(2:7).EQ.'MEGACI') THEN
            ISETKW(20) = 1
            IDIAG = 1
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'TERACI') THEN
            ISETKW(20) = 1
            IDIAG = 2
            GOTO 999
         END IF
*==================================
* 21 : Explicit hamilton matrix   : MXP1,MXP2,MXQ
*==================================
         IF(CARD(2:7).EQ.'EXPHAM') THEN
            READ(LUCMOL,*) MXP1,MXP2,MXQ
            ISETKW(21) = 1
            GOTO 999
         END IF
*===================================================
* 22 : Largest allowed number of Iterations per root : MAXIT
*===================================================
         IF(CARD(2:6).EQ.'MAXIT') THEN
            ISETKW(22) = 1
            READ(LUCMOL,*) MAXIT
            GOTO 999
         END IF
*====================
* 23 : Restart option
*====================
         IF(CARD(2:7).EQ.'RESTRT') THEN
            ISETKW(23) = 1
            IRESTR = 1
            GOTO 999
         END IF
*========================================
* 24 Import of integrals and environment
*========================================
         IF(CARD(2:7).EQ.'MOLCAS'.OR.CARD(2:7).EQ.'ENV=MO') THEN
*. Integrals imported from MOLCAS
            INTIMP = 1
            ISETKW(24) = 1
            ENVIRO(1:6) = 'MOLCAS'
            GOTO 999
         ELSE IF(CARD(2:6).EQ.'LUCAS')THEN
*. Integrals imported from LUCAS
            INTIMP = 2
            ISETKW(24) = 1
            ENVIRO(1:6) = 'LUCAS '
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'FMINSM'.OR.CARD(2:7).EQ.'ENV=FM'
     &           .OR.CARD(2:7).EQ.'ENV=LU') THEN
*. Internal LUCIA environment as generated by a previous LUCIA run.
*
*. Integrals read formatted in, only integrals differing from
*. zero by symmetry are  included
            INTIMP = 3
            ISETKW(24) = 1
            ENVIRO(1:6) = 'LUCIA '
            GOTO 999
         ELSE IF(CARD(2:7).EQ.'SIRIUS'.OR.CARD(2:7).EQ.'DALTON'
     &           .OR.CARD(2:7).EQ.'ENV=DA') THEN
*. Integrals imported from SIRIUS/DALTON
            INTIMP = 5
C         write(6,*) ' Sirius Flag activated '
            ISETKW(24) = 1
            ENVIRO(1:6) = 'DALTON'
            GOTO 999
         else if (CARD(2:7).eq.'DIRAC ') then
           read(LUCMOL,*) IDBGRP
           if (IDBGRP.eq.1) DOUGRP = 'D2h'
           if (IDBGRP.eq.2) DOUGRP = 'D2 '
           if (IDBGRP.eq.3) DOUGRP = 'C2v'
           if (IDBGRP.eq.4) DOUGRP = 'C2h'
           if (IDBGRP.eq.5) DOUGRP = 'C2 '
           if (IDBGRP.eq.6) DOUGRP = 'Cs '
           if (IDBGRP.eq.7) DOUGRP = 'Ci '
           if (IDBGRP.eq.8) DOUGRP = 'C1 '
           INTIMP = 6
           ISETKW(24) = 1
           ENVIRO(1:6) = 'DIRAC '
           goto 999
         ELSE IF(CARD(2:7).EQ.'ENV=NO' ) THEN
*. No program environment, integrals, coefs will just be set to zero
            ENVIRO(1:6) = 'NONE  '
            ISETKW(24) = 1
            GOTO 999
         END IF
* ===============================
* 25 :INCORE option for integrals
* ==============================
         IF(CARD(2:7).EQ.'INCORE') THEN
            ISETKW(25) = 1
            INCORE = 1
            GOTO 999
         END IF
* ===================
* 26 : Deleted shells
* ===================
         IF(CARD(2:7).EQ.'DELETE') THEN
            ISETKW(26) = 1
            READ(LUCMOL,*) (NDELSH(IRREP),IRREP= 1, NIRREP)
            GOTO 999
         END IF
* ===================
* 27: Ms combinations
* ===================
         IF(CARD(2:7).EQ.'MSCOMB') THEN
            ISETKW(27) = 1
            READ(LUCMOL,*) PSSIGN
            IF(.NOT.(PSSIGN.EQ.-1.0D0.OR.PSSIGN.EQ.1.0D0)) THEN
              WRITE(luwrt,*)' Illegal Spin combination factor ',PSSIGN
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
* ===================
* 28: Ml combinations
* ===================
         IF(CARD(2:7).EQ.'MLCOMB') THEN
            ISETKW(28) = 1
            READ(LUCMOL,*) PLSIGN
            IF(.NOT.(PLSIGN.EQ.-1.0D0.OR.PLSIGN.EQ.1.0D0)) THEN
               WRITE(luwrt,*)' Illegal ml combination factor ',PLSIGN
               NERROR = NERROR + 1
            END IF
            GOTO 999
         END IF
* ======================================
* 29 : Print flag for string information
* ======================================
         IF(CARD(2:7).EQ.'IPRSTR') THEN
            ISETKW(29) = 1
            READ(LUCMOL,*) IPRSTR
            GOTO 999
         END IF
* ========================================
* 30 : Print flag for CI space information
* ========================================
         IF(CARD(2:7).EQ.'IPRCIX') THEN
            ISETKW(30) = 1
            READ(LUCMOL,*) IPRCIX
            GOTO 999
         END IF
* =======================================
* 31 : Print flag for orbital information
* =======================================
         IF(CARD(2:7).EQ.'IPRORB') THEN
            ISETKW(31) = 1
            READ(LUCMOL,*) IPRORB
            GOTO 999
         END IF
* ===============================================
* 32 : Print flag for diagonalization information
* ===============================================
         IF(CARD(2:7).EQ.'IPRDIA') THEN
            ISETKW(32) = 1
            READ(LUCMOL,*) IPRDIA
            GOTO 999
         END IF
* ===============================================
* 36 : Print flag for Externals
* ===============================================
         IF(CARD(2:6).EQ.'IPRXT') THEN
            ISETKW(36) = 1
            READ(LUCMOL,*) IPRXT
            GOTO 999
         END IF
* =====================================
* 43 : Print occupation of lowest Dets
* =====================================
         IF(CARD(2:7).EQ.'IPROCC') THEN
            ISETKW(43) = 1
            READ(LUCMOL,*) IPROCC
            GOTO 999
         END IF
* ====================================
* 65 : Print level for density matrices
* ====================================
         IF(CARD(2:7).EQ.'IPRDEN') THEN
            ISETKW(65) = 1
            READ(LUCMOL,*) IPRDEN
            GOTO 999
         END IF
* ===========================================
* 84 : Print level for Response calculations
* ===========================================
         IF(CARD(2:7).EQ.'IPRRSP') THEN
            ISETKW(84) = 1
            READ(LUCMOL,*) IPRRSP
            GOTO 999
         END IF
*
*=========================================================
* 33 : Largest allowed number of Vectors in diagonalization
*=========================================================
         IF(CARD(2:6).EQ.'MXCIV') THEN
            ISETKW(33) = 1
            READ(LUCMOL,*) MXCIV
            GOTO 999
         END IF
* =============================
* 34 : Storage mode for vectors
* =============================
         IF(CARD(2:7).EQ.'CISTOR')THEN
            ISETKW(34) = 1
            READ(LUCMOL,*) ICISTR
            GOTO 999
         END IF
* ================================
* 35 : Do not employ CSF expansion
* ================================
         IF(CARD(2:6).EQ.'NOCSF') THEN
            ISETKW(35) = 1
            NOCSF = 1
            GOTO 999
         END IF
* ================================
* 37 : Do not read in integrals
* ================================
         IF(CARD(2:6).EQ.'NOINT') THEN
            ISETKW(37) = 1
            NOINT = 1
            GOTO 999
         END IF
* ================================
* 38 : Dump integrals in formatted form on file 90
* ================================
         IF(CARD(2:7).EQ.'DMPINT') THEN
            ISETKW(38) = 1
            IDMPIN = 1
            GOTO 999
         END IF
* ================================
* 39 : Define dimension of resolution matrices
* ================================
         IF(CARD(2:7).EQ.'RESDIM') THEN
            ISETKW(39) = 1
            READ(LUCMOL,*) MXINKA
            GOTO 999
         END IF
* ====================================================================
* 40 : Use CJKAIB matrices as intermediate matrices in alpha-beta-loop
* ====================================================================
         IF(CARD(2:7).EQ.'CJKAIB') THEN
            ISETKW(40) = 1
            ICJKAIB = 1
            GOTO 999
         END IF
* ====================================================================
* 44 : Use Minimal operatioon count method for alpha-alpha and beta-beta
* ====================================================================
         IF(CARD(2:6).EQ.'MOCAA') THEN
            ISETKW(44) = 1
            MOCAA = 1
            GOTO 999
         END IF
* ====================================================================
* 45 : Use Minimal operatioon count method for alpha-beta
* ====================================================================
         IF(CARD(2:6).EQ.'MOCAB') THEN
            ISETKW(45) = 1
            MOCAB = 1
            GOTO 999
         END IF
* ====================================================================
* 41 : Initial CI in reference space
* ====================================================================
         IF(CARD(2:7).EQ.'INIREF') THEN
            ISETKW(41) = 1
            INIREF  = 1
            GOTO 999
         END IF
* ====================================================================
* 42 : Restart from reference CI expansion
* ====================================================================
         IF(CARD(2:7).EQ.'RESTRF') THEN
            ISETKW(42) = 1
            IRESTRF = 1
*. Flag that restart will be used for zero space calculation
            ISETKW(23) = 1
            IRESTR = 1
            GOTO 999
         END IF
* ====================================================================
* 46 : Read in of core energy
* ====================================================================
         IF(CARD(2:6).EQ.'ECORE') THEN
            ISETKW(46) = 1
            READ(LUCMOL,*) ECORE
            GOTO 999
         END IF
*
* =====================================================================
* 47 : Use Perturbation theory for zero order space
* =====================================================================
*
         IF(CARD(2:6).EQ.'PERTU') THEN
*
*. Perturbation theory : Three parameters to be specified :
*
*      1 : Max order of correction vectors required
*      2 : Type of partitioning ( H0 )
*          Current choices : MP, EN, H0READ
*      3 : zero order energy : E0=EX ( use exact energy of reference state )
*                              E0=AV ( Use expectation value of H0 )
*                              E0=RE ; Readin zero order energy in
            ISETKW(47) = 1
            IPERT = 1
*. Number of correction vectors
            READ(LUCMOL,*) NPERT
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
C?       WRITE(6,'(A)') CARD1
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENP = 1
               IPART = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENP = 2
               IPART = 2
            ELSE IF(CARD1(1:6).EQ.'H0READ' ) THEN
*. Read in one body hamiltonian
               MPORENP = 0
               IPART = 3
*.
            ELSE
               WRITE(luwrt,*) ' Unknown partitioning '
               WRITE(luwrt,'(1X ,A)') CARD1
               NERROR = NERROR + 1
            END IF
* Zero order energy :
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
C?       WRITE(6,'(A)') CARD1
*
            IF(CARD1(1:5).EQ.'E0=AV') THEN
               IE0AVEX = 1
            ELSE IF(CARD1(1:5).EQ.'E0=EX') THEN
               IE0AVEX = 2
            ELSE IF(CARD1(1:5).EQ.'E0=RE') THEN
               IE0AVEX = 3
               READ(5,*) E0READ
               WRITE(luwrt,*) ' Zero order energy =',E0READ
            ELSE
               WRITE(luwrt,*) ' Unknown form of zero order energy '
               WRITE(luwrt,'(1X ,A)') CARD1
               NERROR = NERROR + 1
            END IF
*
            GOTO 999
         END IF
*
* =====================================================================
* 48 : Approximate Hamiltonian in reference space
* =====================================================================
*
         IF(CARD(2:7).EQ.'APRREF') THEN
            ISETKW(48) = 1
            READ(LUCMOL,*)  MNRS1RE,MXRS3RE
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENR = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENR = 2
            ELSE
               WRITE(luwrt) ' Unknown partitioning '
               WRITE(luwrt,'(1X ,A)') CARD1
               NERROR = NERROR + 1
            END IF
*
            IAPRREF = 1
            GOTO 999
         END IF
*
* =====================================================================
* 49 : Approximate Hamiltonian in zero order space
* =====================================================================
*
         IF(CARD(2:7).EQ.'APRZER') THEN
            ISETKW(49) = 1
            READ(LUCMOL,*)  MNRS1ZE,MXRS3ZE
*. Moeller-Plesset or Epstein-Nesbet partitioning
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
            CALL UPPCAS(CARD1,MXPLNC)
*
            IF(CARD1(1:2) .EQ. 'MP' ) THEN
               MPORENZ = 1
            ELSE  IF(CARD1(1:2) .EQ. 'EN' ) THEN
               MPORENZ = 2
            ELSE
               WRITE(luwrt) ' Unknown partitioning '
               WRITE(luwrt,'(1X ,A)') CARD1
               NERROR = NERROR + 1
            END IF
            IAPRZER = 1
            GOTO 999
         END IF
*
* =====================================================================
* 50 : Generalized active space concept invoked, orbital spaces
* =====================================================================
*
         IF(CARD(2:6).EQ.'GASSH') THEN
*. Generalized active space in use
            ISETKW(50) = 1
            IDOGAS = 1
            READ(LUCMOL,*) NGAS
            IGSFILL = 0
            DO IGAS = 1, NGAS
               READ(LUCMOL,'(A)') CARD1
               CALL LFTPOS(CARD1,MXPLNC)
               CALL UPPCAS(CARD1,MXPLNC)
*. A line can be one of the following
*  NIRREP numbers giving dim of each irrep for this space
* A character entry :
*                     NONE => No orbitals in this space
*                     ALL  => All remaining orbitals  in this space
*                     REST => All remaining orbitals  in this space
*. Note : Only a single space must be defined by ALL or REST
               CALL DECODE_LINE(CARD1,MXPLNC,NITEM,ITEM,MX_ITEM)
               ITEMX = ITEM(1)
               IF(ITEMX(1:4).EQ.'NONE') THEN
                  DO IRREP = 1, NIRREP
                     NGSSH(IRREP,IGAS) = 0
                  END DO
               ELSE IF(ITEMX(1:3).EQ.'ALL'.OR.ITEMX(1:4).EQ.'REST') THEN
*. Only a single space must be defined in this way
                  IF(IGSFILL.NE.0) THEN
                     WRITE(luwrt,*)
     &                   ' Several shell spaces defined by ALL or REST'
                     WRITE(luwrt,*)
     &                   ' This confuses and upsets me '
                     WRITE(luwrt,*)
     &                   '                                / Lucita '
                     ISETKW(50) = -1
                     NERROR = NERROR + 1
                  END IF
                  IGSFILL = IGAS
               ELSE
*. I expect that NIRREP integers are given
                  IF(NITEM.NE.NIRREP) THEN
                     WRITE(luwrt,*) ' Erroneous input to GASSH : '
                     WRITE(luwrt,'(72A)') CARD1
                     WRITE(luwrt,*) ' Specify either :   NONE '
                     WRITE(luwrt,*) '                     ALL'
                     WRITE(luwrt,*) '                    REST'
                     WRITE(luwrt,*) ' Or NIRREP integers '
                     NERROR = NERROR + 1
                     ISETKW(50) = -1
                  END IF
*. Well assume NIRREP integers
                  DO IRREP = 1, NIRREP
                     CALL CHAR_TO_INTEGER_MOLUC(ITEM(IRREP),
     &                                          NGSSH(IRREP,IGAS))
                  END DO
               END IF
*. Number of irreps per GAS
C        READ(LUCMOL,*) (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
            END DO
CTF
*  Parameter check!!
            do IGS=1,NGAS,1
              do IRR=1,NIRREP,1
                if (NGSSH(IRR,IGS).gt.MXTSOB) then
                write(luwrt,*) 'Too many orbitals per GAS and symmetry!'
                  write(luwrt,*) 'My maximum is  ',MXTSOB
                  write(luwrt,*) 'Please redefine your GAS spaces.'
                  Call Abend2( 'Quitting.' )
                end if
              end do
            end do
*
            GOTO 999
         END IF
*
* =====================================================================
* 51 : Generalized active space occupation restrictions
* =====================================================================
*
         IF(CARD(2:7).EQ.'GASSPC') THEN
*. Orbital constraints in gas spaces
*. GASSH must have been defined vefore, check this
            IF(ISETKW(50).EQ.0) THEN
               WRITE(luwrt,*) ' Dear User'
               WRITE(luwrt,*)
               WRITE(luwrt,*) ' GASSH must be specified before GASSPC'
               WRITE(luwrt,*)
     &          ' Else I do not know about the number of orbital spaces'
               WRITE(luwrt,*) ' So I will stop '
               Call Abend2( 'READIN: put GASSH before GASSPC' )
            END IF
            IDOGAS = 1
            ISETKW(51) = 1
*. Number of oribtal spaces
            READ(LUCMOL,*) NCISPC
            DO ISPC = 1, NCISPC
*. Upper and lower limits for each orbital space
               READ(LUCMOL,*)
     &          (IGSOCCX(IGAS,1,ISPC),IGSOCCX(IGAS,2,ISPC),IGAS=1,NGAS)
            END DO
            GOTO 999
         END IF
* ==========================================================================
* 52 : Calculations will be performed in combination of different GAS spaces
* ==========================================================================
*
         IF(CARD(2:7).EQ.'CMBSPC') THEN
         IDOGAS = 1
         ISETKW(52) = 1
*. Check if SEQUEN have been specified.
         IF(ISETKW(54).EQ.1) THEN
           WRITE(luwrt,*) ' Dear user '
           WRITE(luwrt,*)
           WRITE(luwrt,*)' SEQUEN flag has been specified before CMBSPC'
           WRITE(luwrt,*)' This confuses me and makes me wonder 
     & what the'
           WRITE(luwrt,*)' meaning of everything is. '
           WRITE(luwrt,*)' Please ensure that CMPSPC is given 
     &before SEQUEN'
           WRITE(luwrt,*)
           WRITE(luwrt,*)'                                  Lucita  '
           WRITE(luwrt,*)
           Call Abend2( 'READIN: Specify CMBSPC before SEQUEN' )
         END IF
*. Number of combination spaces
         READ(LUCMOL,*) NCMBSPC
         DO JCMBSPC  = 1, NCMBSPC
*. Number of gas spaces in this space
            READ(LUCMOL,*) LLCMBSPC
            LCMBSPC(JCMBSPC) = LLCMBSPC
*. Gasspaces included
           READ(LUCMOL,*) (ICMBSPC(IGASSPC,JCMBSPC),IGASSPC=1,LLCMBSPC)
         END DO
         GOTO 999
      END IF
*
* ===============================================
* 53 : Energy convergence of CI
* ===============================================
*
      IF(CARD(2:7).EQ.'CICONV') THEN
         READ(LUCMOL,*) THRES_E
         ISETKW(53) = 1
         GOTO 999
      END IF
*
      IF(CARD(2:7).EQ.'SEQUEN') THEN
         ISETKW(54) = 1
*
* =========================================
* 54 : SEQUEN KEYWORD
* =========================================
*
* Form of input is
*
* Loop over CI spaces
*  READ NCALC <= Number of calculations in this space
*  Loop Over the NCALC calculations
*    READ type_of_calculation, further info ( remember the comma)(see below)
*  End of loop over NCALC calulation
* End of loop over CI spaces
*
*. Is total number of CI spaces defined ?
         IF(ISETKW(52).EQ.0) THEN
*. Combination spaces were not explicitly defined,
*. assume each gas space is a conb space
            NCMBSPC = NCISPC
         END IF
*
         DO JCMBSPC = 1, NCMBSPC
            READ(LUCMOL,*) NSEQCI(JCMBSPC)
            DO ICI = 1, NSEQCI(JCMBSPC)
*. Read in as character line, and decode
*. Format : Type of calc, further info
*. Possible types of calculations :
* =================================
*    CI : Normal  CI
*    APR-CI : CI with approximate Hamiltonian
*    PERTU  : Perturbation theory, high order version with vectors on
*             disc
*    VECFREE: Various vector free calculations
*    CC     : Coupled Cluster calculation  (CCSD , Extremely inefficient !)
*    ICCI   : Internal contracted CI
*    ICPT   : Internal contracted    PT

*
               READ(LUCMOL,'(A)') CARD1
               CALL LFTPOS(CARD1,MXPLNC)
               CALL UPPCAS(CARD1,MXPLNC)
               CALL DECODE_LINE(CARD1,MXPLNC,NITEM,ITEM,MX_ITEM)
*. Type of calc :
               CARDX=ITEM(1)
               CSEQCI(ICI,JCMBSPC) = ITEM(1)
*
* CI or CI with approximate hamiltonian
*
               IF(CARDX(1:2).EQ.'CI'     .OR.
     &              CARDX(1:6).EQ.'APR-CI'     ) THEN
*. CI calculation, second item in line will be max number of its'
                  IF(NITEM.EQ.1) THEN
*. No second item, use default number of iterations: maybe not
*. defined presently, so flag by a minus and insert later
*
* At the moment : I want MAXIT as the second entry
                     WRITE(luwrt,*)
     &                   ' ERROR :  Number_of_iterations not specified'
                     WRITE(luwrt,*)
     &                   ' Required form of CI card is : '
                     WRITE(luwrt,*) ' CI , Number_of_iterations'
                     ISEQCI(ICI,JCMBSPC) = -1
                     NERROR = NERROR + 1
                     ISETKW(54) = -1
                  ELSE
                     CALL CHAR_TO_INTEGER_MOLUC(ITEM(2),
     &                                          ISEQCI(ICI,JCMBSPC))
                  END IF
               ELSE IF(CARDX(1:5).EQ.'PERTU') THEN
*. Perturbation calculation, following items are
* Maxord, Ipart, E0 with
* 1) Maxord : order to which perturbation vectors will be solved
* 2) Ipart  :  Partitioning of zero order Hamiltonian,
*              MP-DIAG : Diagonal Moller-Plesset operator
*              MP-FULL : Full nondiagonal Moller-Plesset operator
*              EN      : Epstein-Nesbet : Hamiltonian diagonal
*              GENH0   : General H0, specified by separate keyword
* 3) E0     :  Definition of zero order energy
*              E0=EX : Use exact energy of zero order state
*              E0=AV : Use average Zero order energy
*              E0=RE : Read in exact zero .
*
* First time around : No extra info, use normal perturbation keyword
* PERTU to specify perturbation calculation
*
               ELSE IF(CARDX(1:7).EQ.'VECFREE') THEN
*
* ========================
*. Vector free calculation
* ========================
*
*. Second entry is level of calculation
*
*              LEVEL = 1 => second order perturbation calc
*              LEVEL = 2 => + 1 CI it + third order calc
*              LEVEL = 3 => 1 MP4 in current CI space
*              LEVEL = 4 => Level 2 + MP4 in next space
*
                  IF(NITEM.EQ.1) THEN
                     WRITE(luwrt,*)
     &                    ' ERROR :  Level parameter not specified'
                     WRITE(luwrt,*)
     &                    ' Required form of VECFREE card is : '
                     WRITE(luwrt,*) ' VECFREE , LEVEL'
                     ISEQCI(ICI,JCMBSPC) = -1
                     NERROR = NERROR + 1
                  ELSE
                     CALL CHAR_TO_INTEGER_MOLUC(ITEM(2),
     &                                          ISEQCI(ICI,JCMBSPC))
*. Level parameter is traditionally specified by negative number,
                     ISEQCI(ICI,JCMBSPC) = -ISEQCI(ICI,JCMBSPC)
                  END IF
               ELSE IF(CARDX(1:2).EQ.'CC') THEN
*
* ==============================
*. Coupled Cluster calculation
* ==============================
*
                  WRITE(luwrt,*) ' CC routines will be called '
                  I_DO_CC = 1
*. Last calculation which is CC
                  LAST_CC_SPC = JCMBSPC
                  LAST_CC_RUN = ICI
                  IF(NITEM.EQ.1) THEN
* At the moment : I want MAXIT as the second entry
                     WRITE(luwrt,*)
     &                    ' ERROR :  Number_of_iterations not specified'
                     WRITE(luwrt,*)
     &                    ' Required form of CC card is : '
                     WRITE(luwrt,*) ' CC , Number_of_iterations'
                     ISEQCI(ICI,JCMBSPC) = -1
                     NERROR = NERROR + 1
                     ISETKW(54) = -1
                  ELSE
                     CALL CHAR_TO_INTEGER_MOLUC(ITEM(2),
     &                                          ISEQCI(ICI,JCMBSPC))
                  END IF
               ELSE IF(CARDX(1:4).EQ.'ICCI' ) THEN
*
* ==============================
*. Internal contracted CI calculation
* ==============================
*
                  WRITE(luwrt,*) ' ICCI routines will be called '
               ELSE IF(CARDX(1:4).EQ.'ICPT' ) THEN
*
* ==============================
*. Internal contracted PT calculation
* ==============================
*
                  WRITE(luwrt,*) ' Internal contracted PT '
               ELSE
*
                  WRITE(luwrt,'(A,A)')
     &              ' Unknown type of calculation specified in SEQUEN',
     &                 CARDX
                  WRITE(luwrt,*) ' Allowed ENTRIES : '
                  WRITE(luwrt,*) ' ================='
                  WRITE(luwrt,*) '     CI'
                  WRITE(luwrt,*) '     APR_CI'
                  WRITE(luwrt,*) '     PERTU '
                  WRITE(luwrt,*) '     VECFREE'
                  WRITE(luwrt,*) '     CC     '
                  WRITE(luwrt,*) '     ICCI   '
                  WRITE(luwrt,*) '     ICPT   '
                  NERROR = NERROR + 1
                  ISETKW(54) = -1
               END IF
            END DO
*          ^ End of loop over calculations for given CI space
         END DO
*        ^ End of loop over CI spaces
*. The old input for the SEQUEN : Short and numeric !:
C          IF(NSEQCI(JCMBSPC).GT.0)
C    &     READ(LUCMOL,*) (ISEQCI(ICI,JCMBSPC),ICI = 1, NSEQCI(JCMBSPC))
C        END DO
         GOTO 999
      END IF
*
* =====================================================================
* 55 : Call EXTENDED KOOPMANS' THEOREM ROUTINE
* =====================================================================
*
      IF(CARD(2:7).EQ.'EXTKOP') THEN
*. Oh yes, we will do it !
         IEXTKOP = 1
         ISETKW(55) = 1
         GOTO 999
      END IF
*
* ==========================
* 56 : What's your engine ?
* ==========================
*
      IF(CARD(2:7).EQ.'MACHIN') THEN
         ISETKW(56) = 1
         READ(LUCMOL,'(A)') CARD1
         CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
         CALL UPPCAS(CARD1,MXPLNC)
         !MACHINE(1:6) = CARD1(1:6)
         call quit('MACHINE keyword not active')
C?       WRITE(6,'(A,A)') ' Machine = ', MACHINE
         GOTO 999
      END IF
*
* ==========================================================
* 57 : Save first order correction to wavefunction on DISC?
* ==========================================================
*
* ( For vector free calculations )
*
      IF(CARD(2:6).EQ.'C1DSC') THEN
         ISETKW(57) = 1
         IC1DSC = 1
         GOTO 999
      END IF
*
* ==========================================================
*.58 :  Specify subspaces in which perturbation is nonvanishing
* ==========================================================
*
      IF(CARD(2:6).EQ.'H0SPC') THEN
*. Ensure that number of GASSPACES have been defined
         IF(ISETKW(50).EQ.0) THEN
            WRITE(luwrt,*) ' Dear User'
            WRITE(luwrt,*)
            WRITE(luwrt,*) ' GASSH must be specified before H0SPC'
            WRITE(luwrt,*)
     &         ' Else I do not know about the number of orbital spaces'
            WRITE(luwrt,*) ' So I will stop '
            Call Abend2( 'READIN: put GASSH before H0SPC ' )
         END IF
         READ(LUCMOL,*) NPTSPC
         IF(NPTSPC.GT.MXPPTSPC) THEN
*
            WRITE(luwrt,*) ' To many perturbation spaces '
            WRITE(luwrt,*)
     &           ' raise MXPPTSPC from ', MXPPTSPC ,' to ',NPTSPC
            Call Abend2( 'NPTSPC>MXPPTSPC in READIN ' )
         END IF
*
         IH0SPC = 1
         DO JPTSPC = 1, NPTSPC
*. Number of occupation spaces in this subspace
C          DO JGAS = 1, NGAS
            READ(LUCMOL,*)
     &           (IOCPTSPC(1,JGAS,JPTSPC),IOCPTSPC(2,JGAS,JPTSPC),
     &           JGAS = 1, NGAS)
C          END DO
         END DO
         ISETKW(58) = 2
         GOTO 999
      END IF
*
* ============================================
*.59 :  Specify Type of H0 for each subspace
* ============================================
*
      IF(CARD(2:6).EQ.'H0FRM') THEN
*. Ensure that number of Perturbation subspaces have been defined
         IF(ISETKW(58).EQ.0) THEN
            WRITE(luwrt,*) ' Dear User'
            WRITE(luwrt,*)
            WRITE(luwrt,*) ' H0SPC must be specified before H0FRM'
            WRITE(luwrt,*)
     &           ' Else I do not know about the number of spaces'
            WRITE(luwrt,*) ' So I will stop '
            Call Abend2( 'READIN: put H0SPC before H0FRM ' )
         END IF
*. Type of perturbation in this subspace
*
* 1 => Diagonal MP
* 2 => EN
* 3 => Nondiagonal MP
* 4 => Exact Hamiltonian
* 5 => Nondiagonal FI+FA + exact in orbital subspaces
*
         DO JPTSPC = 1, NPTSPC
            READ(LUCMOL,*) IH0INSPC(JPTSPC)
         END DO
         ISETKW(59) = 2
         GOTO 999
      END IF
*
* =============================================
* 60 : Reference root for Perturbation theory
* =============================================
*
      IF(CARD(2:7).EQ.'RFROOT') THEN
         ISETKW(60) = 1
         READ(LUCMOL,*) IRFROOT
C        WRITE(6,*) ' Reference Root = ',IRFROOT
         GOTO 999
      END IF
*
* ======================================================
* 61 : Orbital spaces in which Exact Hamiltonian is used
* ======================================================
*
      IF(CARD(2:5).EQ.'H0EX') THEN
         ISETKW(61) = 1
         GOTO 999
      END IF
*
* ================================================
* 62 : Treatment of degenerencies of initial guess
* ================================================
*
      IF(CARD(2:7).EQ.'INIDEG') THEN
         ISETKW(62) = 1
         READ(LUCMOL,*) INIDEG
         GOTO 999
      END IF
*
* ========================================================
* 63 : Use modified Hamilton operator in CI optimization
* ========================================================
*
      IF(CARD(2:7).EQ.'LAMBDA') THEN
         ISETKW(63) = 1
         READ(LUCMOL,*) XLAMBDA
         GOTO 999
      END IF
*
* =============================================================
* 64 : Length of smallest block for batch of C an Sigma vectors
* =============================================================
*
      IF(CARD(2:7).EQ.'LCSBLK') THEN
         ISETKW(64) = 1
         READ(LUCMOL,*) LCSBLK
         GOTO 999
      END IF
*
*
* =============================================================
* 66 : xxx
* =============================================================
*
* =============================================================
* 67 : ECHO the following keywords
* =============================================================
*
      IF(CARD(2:5).EQ.'ECHO') THEN
         IECHO = 1
         ISETKW(67) = 1
         GOTO 999
      END IF
*
*
* ====================
* 68 : Final orbitals
* ====================
*
*. Should be specified after NIRREP, I have not added the
* test!!
      IF(CARD(2:7).EQ.'FINORB') THEN
*. Type of final orbitals
         READ(LUCMOL,'(A)') CARD1
         CALL LFTPOS(CARD1,MXPLNC)
*. Change to upper case
         CALL UPPCAS(CARD1,MXPLNC)
*
         WRITE(luwrt,'(A,A)')
     &        ' Type of final orbitals ',CARD1
         ISETKW(68) = 1
*
         IF(CARD1(1:5).EQ.'NATUR') THEN
*. Natural orbitals
            IFINMO = 1
         ELSE IF ( CARD(1:5).EQ.'CANON' ) THEN
*. Canonical orbitals
            IFINMO = 2
         ELSE IF ( CARD1(1:6).EQ.'PS_NAT') THEN
*. Pseudo natural orbitals
            IFINMO = 3
         ELSE IF ( CARD1(1:6) .EQ. 'PS_CAN') THEN
*. Pseudo canonical orbitals
            IFINMO = 4
         ELSE IF (CARD1(1:6) .EQ. 'PS_NC') THEN
*. Pseudo natural-canonical orbitals
            IFINMO = 5
*. requires input of subshells in which to define
*. Pseudo-natural orbitals
            READ(LUCMOL,*) NPSSPC
            DO IPSSPC = 1, NPSSPC
               READ(LUCMOL,*) (NPSSH(IRREP,IPSSPC),IRREP=1,NIRREP)
            END DO
         ELSE
*. Unidentified type of final orbitals
            WRITE(luwrt,*) ' Unidentified type of final orbitals'
            WRITE(luwrt,'(A,A)') '  you suggested : ', CARD1
            WRITE(luwrt,*)
            WRITE(luwrt,*) ' Allowed types of final orbitals'
            WRITE(luwrt,*) ' ==============================='
            WRITE(luwrt,*)
            WRITE(luwrt,*) '     NATUR'
            WRITE(luwrt,*) '     CANON'
            WRITE(luwrt,*) '     PS_NAT'
            WRITE(luwrt,*) '     PS_CAN'
            WRITE(luwrt,*) '     PS_NC'
            NERROR = NERROR + 1
            ISETKW(68) = - 1
         END IF
         GOTO 999
*
      END IF
*
*
* ===================================================================
* 69 : Threshold on second order energy corrections, individual coefs
* ===================================================================
*
      IF(CARD(2:7).EQ.'E_THRE') THEN
         READ(LUCMOL,*) E_THRE
         ISETKW(69) = 1
         GOTO 999
      END IF
*
*
* =======================================================================
* 70 : Threshold on first order wavefunction corrections,individual coefs
* =======================================================================
*
      IF(CARD(2:7).EQ.'C_THRE') THEN
         READ(LUCMOL,*) C_THRE
         ISETKW(70) = 1
         GOTO 999
      END IF
*
* ===================================================================
* 71 : Threshold on second order energy corrections, Total Threshold
* ===================================================================
*
      IF(CARD(2:7).EQ.'E_CONV') THEN
         READ(LUCMOL,*) E_CONV
         ISETKW(71) = 1
         GOTO 999
      END IF
*
*
* =======================================================================
* 72 : Threshold on first order wavefunction corrections,Total Threshold
* =======================================================================
*
      IF(CARD(2:7).EQ.'C_CONV') THEN
         READ(LUCMOL,*) C_CONV
         ISETKW(72) = 1
         GOTO 999
      END IF
*
*
* ===============================
* 73 : Selection of classes
* ===============================
*
      IF(CARD(2:7).EQ.'CLSSEL') THEN
         ICLSSEL = 1
         ISETKW(73) = 1
         GOTO 999
      END IF
*
*
* =====================================
* 74 : Calculation of density matrices
* ======================================
*
      IF(CARD(2:6).EQ.'DENSI') THEN
         READ(LUCMOL,*) IDENSI
*. IDENSI = 0 => No calculation of density matrices
*  IDENSI = 1 =>  Calculation of one- body density matrix
*  IDENSI = 2 =>  Calculation of one- and two-body density matrices
         ISETKW(74) = 1
         GOTO 999
      END IF
*
*
*
* =====================================
* 75 : Perturbation expansion of EKT
* ======================================
*
      IF(CARD(2:6).EQ.'PTEKT') THEN
         IPTEKT = 1
*. Number of EKT to be analyzed, atmost 20
         READ(5,*)  NPTEKT
         IF(NPTEKT.GT.20) THEN
            WRITE(luwrt,*) ' Atmost 20 perturbation expansions'
            Call Abend2( ' NPTEKT in .PTEKT to Large ' )
         END IF
*. orbital and symmetry for zero order solution
         DO JEKT = 1, NPTEKT
            READ(5,*) LPTEKT(1,JEKT),LPTEKT(2,JEKT)
         END DO
         ISETKW(75) = 1
C?       WRITE(6,*) ' NPTEKT = ', NPTEKT
C?       WRITE(6,*) ' LPTEKT = ',LPTEKT(1,1),LPTEKT(2,1)
         GOTO 999
      END IF
*
* =================================================
* 76 : Root used to define Zero order Hamiltonian
* =================================================
*
      IF(CARD(2:7).EQ.'H0ROOT') THEN
         ISETKW(76) = 1
         READ(LUCMOL,*) IH0ROOT
C        WRITE(6,*) ' Reference Root = ',IH0ROOT
         GOTO 999
      END IF
*
* ======================================
* 77 : No restart in CI calculation 2
* =====================================
*
      IF(CARD(2:7).EQ.'NORST2') THEN
         ISETKW(77) = 1
         IRST2 =  0
         WRITE(luwrt,*) ' NORST2 flag read '
         GOTO 999
      END IF
*
* =====================================================
* 78 : Skip initial evaluation of energy from CI calc 2
* ====================================================
*
      IF(CARD(2:7).EQ.'SKIPEI') THEN
         ISETKW(78) = 1
         ISKIPEI =  1
         WRITE(luwrt,*) ' SKIPEI flag set  '
         GOTO 999
      END IF
*
* =================================================================
* 79 : Symmetry of X, Y and Z - Yes it could be obtained from files
* ================================================================
*
      IF(CARD(2:7).EQ.'XYZSYM') THEN
         ISETKW(79) = 1
         READ(LUCMOL,*) (IXYZSYM(I),I=1,3)
C?       WRITE(6,*) 'IXYZSYM', (IXYZSYM(I),I=1,3)
         GOTO 999
      END IF
*
* ==============================================
* 80 : One-electron properties to be calculated
* ==============================================
*
      IF(CARD(2:7).EQ.'PROPER') THEN
         ISETKW(80) = 1
         READ(LUCMOL,*) NPROP
         DO IPROP = 1, NPROP
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
            PROPER(IPROP)=CARD1(1:6)
            WRITE(6,'(A,A)') ' Property to be calculated ',
     &           PROPER(IPROP)
         END DO
         GOTO 999
      END IF
*
* ==============================================
* 81 : Transition properties
* ==============================================
*
      IF(CARD(2:7).EQ.'TRAPRP') THEN
         ISETKW(81) = 1
*. Number and symmetry of additional states
         READ(LUCMOL,*) IEXCSYM, NEXCSTATE
C        READ(LUCMOL,*) NEXCSTATE
         ITRAPRP = 1
         GOTO 999
      END IF
*
* ================================
* 82 : CI response calculations
* ================================
*
*. Input goes as
*
* Labels for operators for which average values will be calculated ( A-ops)
* Number of response calculations
* Loop over calculations
* Label for pertop1, Label for pertop1, order for op1, order for op2, freq
* End of loop over calculations
* The first per operator is static, the second can be dynamic ( freq.ne.0)
*
* Example
*
*  XDIPLEN, ZDIPLEN
*  1
*  XDIPLEN, YDIPLEN, 2, 2, 0.0D0
*  Labels of oper
      IF(CARD(2:7).EQ.'RESPON') THEN
         ISETKW(82) = 1
*. Yes I will do respons
         IRESPONS = 1
*. Labels for operators whose expectation values will be expanded
         MXNRESP =20
         READ(LUCMOL,'(A)') CARD1
         CALL LFTPOS(CARD1,MXPLNC)
         CALL UPPCAS(CARD1,MXPLNC)
         CALL DECODE_LINE(CARD1,MXPLNC,NITEM,ITEM,MX_ITEM)
         N_AVE_OP = NITEM
         IF(N_AVE_OP.GT. MXNRESP) THEN
            WRITE(luwrt,*) ' READIN : Error for keyword RESPON'
            WRITE(luwrt,*) ' Specified number of A ops = ', N_AVE_OP
            WRITE(luwrt,*) ' Larger than MAX = ', MXNRESP
            WRITE(luwrt,*) ' PLEASE reduce NAVE_OP and RETURN '
            Call Abend2( ' READIN, KEYWORD RESPON : NAVE_OP .gt. 20 ' )
         END IF
         DO JITEM = 1, NITEM
            AVE_OP(JITEM) = ITEM(JITEM)
         END DO
*. Number of respons calculations to be performed
         READ(LUCMOL,*) NRESP
         IF(NRESP.GT. MXNRESP) THEN
            WRITE(luwrt,*) ' READIN : Error for keyword RESPON'
            WRITE(luwrt,*) ' Specified number of calcs = ', NRESP
            WRITE(luwrt,*) ' Larger than MAX = ', MXNRESP
            WRITE(luwrt,*) ' PLEASE reduce NRESP and RETURN '
            Call Abend2( ' READIN, KEYWORD RESPON : NRESP .gt. 20 ' )
         END IF
         DO IRESP = 1, NRESP
*. Operator1, Operator 2, Maxord for op1, Maxord for op2, freq
* ( Remember commas in betweeen !!)
*. Read in as character line, and decode
            READ(LUCMOL,'(A)') CARD1
            CALL LFTPOS(CARD1,MXPLNC)
            CALL UPPCAS(CARD1,MXPLNC)
            CALL DECODE_LINE(CARD1,MXPLNC,NITEM,ITEM,MX_ITEM)
*. Entries 1 and 2 : the operators in character form
            RESP_OP(1,IRESP) = ITEM(1)
            RESP_OP(2,IRESP) = ITEM(2)
C?         WRITE(6,'(A,A,A)') ' RESP( ,1),RESP( ,2)=  ',
C?   &     RESP_OP(1,IRESP) , RESP_OP(2,IRESP)
*. Entries 3 and 4 : integers, maxord
            CALL CHAR_TO_INTEGER_MOLUC(ITEM(3),MAXORD_OP(1,IRESP))
            CALL CHAR_TO_INTEGER_MOLUC(ITEM(4),MAXORD_OP(2,IRESP))
            IF(NITEM.EQ.4) THEN
*. No frequency
               RESP_W(IRESP) = 0.0
            ELSE
               CALL CHAR_TO_REAL(ITEM(5),RESP_W(IRESP),MXPLNC)
            END IF
         END DO
*
         GOTO 999
      END IF
*
* ==============================================
* 83 : Max number of iterations in lin.eq
* ==============================================
*
      IF(CARD(2:7).EQ.'MXITLE') THEN
         ISETKW(83) = 1
*. Number and symmetry of additional states
         READ(LUCMOL,*) MXITLE
         GOTO 999
      END IF
*
* ==============================================
* 85 : Root homing
* ==============================================
*
      IF(CARD(2:7).EQ.'RTHOME') THEN
         ISETKW(85) = 1
         IROOTHOMING = 1
         GOTO 999
      END IF
*
* ==============================================
* 86 : Allow Particle-hole simplifications
* ==============================================
*
      IF(CARD(2:7).EQ.'USE_PH') THEN
         ISETKW(86) = 1
         IUSE_PH = 1
         GOTO 999
      END IF
*
* ==============================================
* 87 : Allow the sigma routine to take advice
* ==============================================
*
      IF(CARD(2:7).EQ.'ADVICE') THEN
         ISETKW(87) = 1
         IADVICE = 1
         GOTO 999
      END IF
*
* ================================================================
* 88 : Transform CI vectors to alternative orbital representation
* ================================================================
*
      IF(CARD(2:6).EQ.'TRACI') THEN
         ITRACI = 1
         ISETKW(88) = 1
*. Read Form or orbitals to which expansion should be formed
*
* Two pieces of info required :
*  1 : Complete rotations or just rotations internal rotations on GAS space
*      Keywords : Restrict or complete
*  2 : Form of final orbitals
*      Keywords : Canonical or Natural
*      As usual the input is written as keyword1, keyword2
*
         READ(LUCMOL,'(A)') CARD1
         CALL LFTPOS(CARD1,MXPLNC)
         CALL UPPCAS(CARD1,MXPLNC)
         CALL DECODE_LINE(CARD1,MXPLNC,NITEM,ITEM,MX_ITEM)
         IF(NITEM.LT. 2) THEN
            WRITE(luwrt,*) ' READIN : Error for keyword TRACI'
            WRITE(luwrt,*) ' Number of items read ', NITEM
            WRITE(luwrt,*)
     &      ' Form of line should be : complete/restrict, fock/natural'
         END IF
*
         ITRACI_CR=ITEM(1)
         ITRACI_CN=ITEM(2)
         IF(    ITRACI_CR(1:4).NE.'REST'
     &        .AND.ITRACI_CR(1:4).NE.'COMP') THEN
            WRITE(luwrt,*) ' Illegal entry under keyword TRACI '
            WRITE(luwrt,*) ' Your suggestion : ', ITRACI_CR
            WRITE(luwrt,*) ' Allowed entries : '
            WRITE(luwrt,*) ' =================='
            WRITE(luwrt,*)    ' COMPlete '
            WRITE(luwrt,*)    ' RESTrict'
            NERROR = NERROR + 1
            ISETKW(88) = -1
         END IF
         IF(    ITRACI_CN(1:4).NE.'CANO'
     &        .AND.ITRACI_CN(1:4).NE.'NATU') THEN
            WRITE(luwrt,*) ' Illegal entry under keyword TRACI '
            WRITE(luwrt,*) ' Your suggestion : ', ITRACI_CN
            WRITE(luwrt,*) ' Allowed entries : '
            WRITE(luwrt,*) ' =================='
            WRITE(luwrt,*)    ' CANOnica'
            WRITE(luwrt,*)    ' NATUral '
            NERROR = NERROR + 1
            ISETKW(88) = -1
         END IF
         GOTO 999
      END IF
*
* ====================================================
* 89 : Separate strings into active and passive parts
* ====================================================
*
      IF(CARD(2:7).EQ.'USE_PA') THEN
         ISETKW(89) = 1
         IUSE_PA = 1
         GOTO 999
      END IF
*
* ==========================================
* 90 : Perturbation expansion of Fock matrix
* ===========================================
*
      IF(CARD(2:7).EQ.'PTFOCK') THEN
         ISETKW(90) = 1
         IPTFOCK = 1
         GOTO 999
      END IF
*
* ==============================
* 91 : Print final CI vectors
* ==============================
*
      IF(CARD(2:7).EQ.'PRNCIV') THEN
         ISETKW(91) = 1
         READ(LUCMOL,*) IPRNCIV
         GOTO 999
      END IF
*
* =====================================================
* 92 : Restart CC calculation (with coefs on LU_CCAMP)
* =====================================================
*
      IF(CARD(2:7).EQ.'RES_CC') THEN
         ISETKW(92) = 1
         I_RESTRT_CC = 1
         GOTO 999
      END IF
*
* =====================================================
* 93 : xxx
* =====================================================
*
* =========================================================
* 94 : xxx
* =========================================================
*
* =========================================================
* 95 : xxx
* =========================================================
*
* =========================================================
* 96 : Use Relaxed densities for properties
* =========================================================
*
      IF(CARD(2:6).EQ.'RELAX') THEN
         ISETKW(96) = 1
         IRELAX= 1
         GOTO 999
      END IF
*
* =========================================================
* 97 : Expert mode : Input errors neglected
* =========================================================
*
      IF(CARD(2:7).EQ.'EXPERT') THEN
         ISETKW(97) = 1
         IEXPERT= 1
         GOTO 999
      END IF
*
*. KEYWORD was not identified
*
      WRITE(luwrt,*)
     &     '  ****  Error, unidentified KEYWORD in READIN   **** '
      WRITE(luwrt,*)
      WRITE(luwrt,*) ' Last line read  '
      WRITE(luwrt,*) ' ================'
      WRITE(luwrt,'(10X,A)') CARD
      WRITE(luwrt,*)
      WRITE(luwrt,*) ' Preceeding KEYWORD'
      WRITE(luwrt,*) ' ==================='
      WRITE(luwrt,'(10X,A)') LASTCARD
      NERROR = NERROR + 1
*
      END IF
 999  CONTINUE
      IF(IECHO.EQ.1)
     &     WRITE(luwrt,'(A,A)') ' processed KEYWORD/COMMENT  : ', CARD
*. Save previous keyword
      LASTCARD(1:72) = CARD(1:72)
      GOTO 1000
*.End of loop over KEYWORDS
 1001 CONTINUE
*

      IF(NERROR.NE.0) THEN
        WRITE(luwrt,'(A)')
     &  ' Run will be aborted due to input errors '
        WRITE(luwrt,'(A,I9)')
     &  ' Number of input errors detected in READIN ', NERROR
*
        WRITE(luwrt,*) 
     &  ' The following keywords were correctly identified'
        WRITE(luwrt,*) 
     &  ' ================================================'
        DO  IENTRY = 1, MXPKW
          IF(ISETKW(IENTRY).EQ.+1)
     &    WRITE(luwrt,'(10X,A)') KEYWOR(IENTRY)
        END DO
        WRITE(luwrt,*)
*
        WRITE(luwrt,*) 'ERRORS were detected for the following KEYWORDS'
        DO IENTRY = 1, MXPKW
           IF(ISETKW(IENTRY).EQ.-1) WRITE(luwrt,'(A)') KEYWOR(IENTRY)
        END DO
        WRITE(luwrt,*)
        WRITE(luwrt,*)
        WRITE(luwrt,*)
        WRITE(luwrt,*)
     &  '     An expert is a man who has made all the mistakes,'
        WRITE(luwrt,*)
     &  '     which can be made, in a very narrow field        '
        WRITE(luwrt,*)
     &  '                                                      '
        WRITE(luwrt,*)
     &  '                                      Niels Bohr      '

        IF(IEXPERT.EQ.0) THEN
          Call Abend2( ' Error in input' )
        ELSE
          WRITE(6,*) ' Program continues (EXPERT mode )'
        END IF
      END IF
*
*  Close LUCIAIN conversion file
*
      close (unit=LUCMOL)
*
**********************************************************************
*                                                                    *
* Part 2 : Insert defaults for missing optional keywords             *
*          and print error messages for missing mandatory keywords   *
*                                                                    *
**********************************************************************
*
      NMISS = 0
*
*.1 : Default title
*
      IF(ISETKW(1).EQ.0) THEN
        TITLEC(1) =
     &  ' Some molecule or some atom                                  '
        TITLEC(2) =
     &  ' Some type of CI expansion                                  '
        TITLEC(3) =
     &  ' Some user who is too lazy to supply a TITLE                 '
        ISETKW(1) = 2
      END IF
*
*.2  Missing pointgroup ( has actually been defaulted )
*
      IF(ISETKW(2).EQ.0) THEN
        PNTGRP = 1
        ISETKW(2) = 2
      END IF
*
*.3 Missing number of irreps, allowed for D2h, illegal else
*
      IF(ISETKW(3).EQ.0) THEN
        IF(PNTGRP .EQ.1 ) THEN
*. Repeat default
          NIRREP = 8
          NSMCMP = NIRREP
          NSMOB  = NIRREP
          ISETKW(3) = 2
        ELSE
*. Number of irreps is mandatory for CINV,DINFH,O3
          NMISS = NMISS + 1
          WRITE(luwrt,*)
     &    '  Input error ! .NIRREP must be specified for CinV,DinH,O3'
        END IF
      END IF
*
* 4 : Internal CI expansion
*
*.Default is CAS
      IF(ISETKW(50) .EQ. 0 ) THEN
      IF(ISETKW(4).EQ.0) THEN
        INTSPC = 1
        ISETKW(4) = 2
*. If a RAS1 or a RAS 3 space has been defined, RAS must have
*  been specified
        IF(ISETKW(9).EQ.1.OR.ISETKW(11).EQ.1) THEN
         ISETKW(4) = 0
         NMISS = NMISS + 1
          WRITE(luwrt,*)
     &    '  Input error ! .RAS must be specified when .RAS1 or .RAS3'
          WRITE(luwrt,*)
     &    '                 has been activated '

        END IF
      END IF
      ELSE IF (ISETKW(50) .EQ. 0 ) THEN
*. FCI expansion
        INTSPC = 3
      END IF
*
* 6 : Number of active electrons
*
*. Mandatory
      IF(ISETKW(6).EQ.0) THEN
        NMISS = NMISS + 1
          WRITE(luwrt,*)
     &    '  Input error ! .NACTEL must be specified '
      END IF
*
* 7 : Inactive orbitals
*
      IF(ISETKW(7).EQ.0) THEN
        CALL ISETVC(NINASH,0,NIRREP)
        ISETKW(7) = 0
      END IF
*
* 8 : Core orbitals, only of interest if EXTSPC .ne. 0
*
      IF(ISETKW(8).EQ.0) THEN
        CALL ISETVC(NRS0SH,0,NIRREP)
        MNHR0 = 0
        IF(EXTSPC.EQ.0) THEN
          ISETKW(8) = 3
        ELSE
          ISETKW(8) = 2
        END IF
      END IF
*
* 9 : RAS 1 orbitals
*
      IF(ISETKW(9).EQ.0) THEN
        CALL ISETVC(NRSSH(1,1),0,NIRREP)
        IF(INTSPC.EQ.1) THEN
          ISETKW(9) = 3
        ELSE
          ISETKW(9) = 2
        END IF
      END IF
*
* 10 : RAS 2 orbitals
*
      IF(ISETKW(10).EQ.0) THEN
        CALL ISETVC(NRSSH(1,2),0,NIRREP)
        ISETKW(10) = 2
      END IF
*
* 11 : RAS 3 orbitals
*
      IMLCR3 = 0
      IF(ISETKW(11).EQ.0) THEN
        CALL ISETVC(NRSSH(1,3),0,NIRREP)
        IF(MOLCS.EQ.1.AND.INTSPC.EQ.2) THEN
*. Use information from one-electron integral file to obtain
* default
          IMLCR3 = 1
        END IF
        IF(INTSPC.EQ.1) THEN
          ISETKW(11) = 3
        ELSE
          ISETKW(11) = 2
        END IF
      END IF
*
* 12 : Partitioning of secondary space ( default 1 set in SECOND)
*
C     IF(ISETKW(12).EQ.0.OR.ISETKW(12).EQ.2) THEN
C       MXR4TP = 1
C       IF(EXTSPC.EQ.0) THEN
C         ISETKW(12) = 3
C       ELSE
C         ISETKW(12) = 2
C       END IF
C     END IF
*
* 13 : Secondary space
*
      IF(ISETKW(13).EQ.0) THEN
        MXR4TP = 1
        DO 210 IR4TP = 1, MXR4TP
          CALL ISETVC(NRS4SH(1,IR4TP),0,NIRREP)
  210   CONTINUE
        MXER4 = 0
        IF(EXTSPC.EQ.0) THEN
          ISETKW(13) = 3
        ELSE
          ISETKW(13) = 2
        END IF
      END IF
*
* 14 : occupation restrictions for Reference space
*
      IF(ISETKW(14).EQ.0) THEN
        MNRS1R = MNRS10
        MXRS3R = MXRS30
        IF(EXTSPC.EQ.0.OR.INTSPC.EQ.1) THEN
          ISETKW(14) = 3
        ELSE
          ISETKW(14) = 2
        END IF
      END IF
*
* 15 : Selection of active configurations
*
      IF(ISETKW(15).EQ.0) THEN
*. Standard is no selection
        INTSEL = 0
      END IF
*
* 16 : Two times spin projection
*
      IF(ISETKW(16).EQ.0) THEN
        WRITE(luwrt,*)
     &  '  Input error ! .MS2 must be specified '
        NMISS = NMISS + 1
      END IF
*
* 17 : Spin multiplicity
*
*. Spin multiplicities : Not active at the moment
CO    IF(ISETKW(17).EQ.0) THEN
CO      WRITE(luwrt,*)
CO   &  '  Input error ! .MULTS must be specified '
CO      NMISS = NMISS + 1
CO    END IF
*
* 18 : Reference symmetry
*
      IF(ISETKW(18).EQ.0) THEN
        WRITE(luwrt,*)
     &  '  Input error ! .IREFSM must be specified '
        NMISS = NMISS + 1
      END IF
*
* 19 : Roots to be optimized
*
      IF(ISETKW(19).EQ.0) THEN
        WRITE(luwrt,*)
     &  '  Input error ! .ROOTS must be specified '
        NMISS = NMISS + 1
      END IF
*
* 20 : Diagonalization routine
*
      IF(ISETKW(20).EQ.0) THEN
*. Standard is currently MICDV*
        IDIAG = 1
        ISETKW(20) = 2
      END IF
*
* 21 : Explicit Hamiltonian
*
      IF(ISETKW(21).EQ.0) THEN
*. Default is no explicit Hamiltonian
        MXP1 = 0
        MXP2 = 0
        MXQ  = 0
        ISETKW(21) = 2
      END IF
*
* 22 : Largest allowed number of CI iterations per root
*
      IF(ISETKW(22).EQ.0) THEN
*. Default is 20 ( not active I expect )
        MAXIT = 20
        ISETKW(22) = 2
      END IF
*
* 23 : Restart option
*
      IF(ISETKW(23).EQ.0) THEN
*. Default is no explicit Hamiltonian
        IRESTR = 0
        ISETKW(23) = 3
      END IF
*
* 24 : Integral import
*
      IF(ISETKW(24).EQ.0) THEN
*. Default is - from NOV26 : Dalton
        INTIMP = 5
        ENVIRO(1:6) = 'DALTON'
        ISETKW(24) = 2
      END IF
*
* 25 : INCORE option for integrals
*
      IF(ISETKW(25).EQ.0) THEN
        IF(EXTSPC.EQ.0 ) THEN
          INCORE = 1
        ELSE
          INCORE = 0
        END IF
        ISETKW(25) = 2
C
C       IF(INTEXP.EQ.0) THEN
C         ISETKW(25) = 3
C       ELSE
C         ISETKW(25) = 3
C       END IF

      END IF
*
* 26 : DELETEd shells
*
      IF(ISETKW(26) .EQ. 0 ) THEN
*. If CAS + Active have been set or RAS + Ras3 have been set,
*. obtain for MOLCAS Interface from number of basis functions
        IF(INTSPC.EQ.1.OR.
     &    (INTSPC.EQ.2.AND.ISETKW(11).EQ.1)) THEN
          IMLCR3 = 2
        ELSE
          CALL ISETVC(NDELSH,0,NIRREP)
        END IF
        ISETKW(26) = 2
      END IF
*
* 27 : Ms combinations
*
      IF(ISETKW(27).EQ.0) THEN
        PSSIGN = 0.0D0
        ISETKW(27) = 2
      ELSE IF(MS2.NE.0) THEN
        WRITE(luwrt,*) ' Spin combinations only allowed with MS2 = 0'
        WRITE(luwrt,*) ' Your value of MS2 = ',MS2, ' differs from zero'
        WRITE(luwrt,*) ' LUCIA will neglect your nice suggestion '
        WRITE(luwrt,*)  ' to use spin combinations '
        PSSIGN = 0.0D0
        ISETKW(27) = 2
      END IF
*
* 28 : Ml combinations
*
      IF(ISETKW(28).EQ.0) THEN
        PLSIGN = 0.0D0
        ISETKW(28) = 2
      ELSE IF(PNTGRP.EQ.1) THEN
        WRITE(luwrt,*) ' Ml combinations not allowed with d2h '
        WRITE(luwrt,*) ' LUCIA will neglect your nice suggestion '
        WRITE(luwrt,*)  ' to use ML combinations '
        PLSIGN = 0.0D0
        ISETKW(28) = 2
      ELSE IF(IREFML.NE.0) THEN
        WRITE(luwrt,*) ' ML combinations only allowed with ML = 0'
        WRITE(luwrt,*)
     &  ' Your value of IREFML = ',IREFML, ' differs from zero'
        WRITE(luwrt,*) ' LUCIA will neglect your nice suggestion '
        WRITE(luwrt,*)  ' to use ML combinations '
        PLSIGN = 0.0D0
        ISETKW(28) = 2
      END IF
      IF(PSSIGN.EQ.0.0D0.AND.PLSIGN.EQ.0.0D0) THEN
        IDC = 1
      ELSE IF(PSSIGN.NE.0.0D0.AND.PLSIGN.EQ.0.0D0) THEN
        IDC = 2
      ELSE IF(PSSIGN.EQ.0.0D0.AND.PLSIGN.NE.0.0D0) THEN
        IDC = 3
      ELSE IF(PSSIGN.NE.0.0D0.AND.PLSIGN.NE.0.0D0) THEN
        IDC = 4
      END IF
C?    WRITE(6,* ) ' TEST readin IDC = ', IDC
*
* 29 : print flag for string information
*
      IF(ISETKW(29).EQ.0) THEN
        IPRSTR = 0
        ISETKW(29) = 2
      END IF
*
* 30 : print flag for string information
*
      IF(ISETKW(30).EQ.0) THEN
        IPRCIX = 0
        ISETKW(30) = 2
      END IF
*
* 31 : print flag for orbital information
*
      IF(ISETKW(31).EQ.0) THEN
        IPRORB = 0
        ISETKW(31) = 2
      END IF
*
* 32 : print flag for diagonalization information
*
      IF(ISETKW(32).EQ.0) THEN
        IPRDIA = 3
        ISETKW(32) = 2
      END IF
*
* 36 : print flag for External blocks
*
      IF(ISETKW(36).EQ.0) THEN
        IPRXT  = 0
        ISETKW(36) = 2
      END IF
*
* 43 : Print occupation of lowest SD's / configurations
*
      IF(ISETKW(43).EQ.0) THEN
        IPROCC = 0
        ISETKW(43) = 2
      END IF
*
* 65 : Print level for densities, default is to print
*      natural occupation numbers only
*
      IF(ISETKW(65).EQ.0) THEN
        IPRDEN = 1
        ISETKW(65) = 2
      END IF
*
* 84 : Print level for response, default is to print
*      final response functions as well as contributions
*
      IF(ISETKW(84).EQ.0) THEN
        IPRRSP = 3
        ISETKW(84) = 2
      END IF
*
* 33 : Number of Ci vectors in subspace
*
      IF(ISETKW(33).EQ.0) THEN
* default is 3/2 vectors per root
        IF(IDIAG.EQ.1) THEN
          MXCIV = 3 * NROOT
        ELSE
          MXCIV = 2 * NROOT
        END IF
        ISETKW(33) = 2
      END IF
*
      IF(ISETKW(33).EQ.1.AND.MXCIV .LT.2*NROOT) THEN
        WRITE(luwrt,*)
     &  '   The number of vectors is increased to 2*NROOT = ',2*NROOT
        MXCIV = 2 * NROOT
      END IF
*
      IF(IDIAG.EQ.2 .AND. MXCIV.GT.2 ) THEN
        MXCIV = 2
        NWARN = NWARN + 1
        WRITE(luwrt,*) ' Warning : You have specified TERACI '
        WRITE(luwrt,*) '           I allow myself to set MXCIV = 2 '
        WRITE(luwrt,*)
        WRITE(luwrt,*) '                   Best Wishes    '
        WRITE(luwrt,*) '                      Lucia       '
      END IF

*
* 34 : CI storage mode
*
      IF(ISETKW(34).EQ.0) THEN
*. Default is three type-type-symmetry blocks
        ICISTR = 3
        ISETKW(34) = 2
      END IF
      IF(ICISTR.EQ.1) THEN
*. complete vectors should not be used together with PICO
      WRITE(luwrt,*)
     &'    You have suggested the use of two complete vectors in core'
      WRITE(luwrt,*)
     &'    Although this could be an interesting suggestion '
      WRITE(luwrt,*)
     &'    I allow myself to reduce the storage mode to 3 sym. blocks '
      ICISTR = 2
      END IF
*
* 35 : Employ CSF expansion ?
*
*. Default is no ( only possibility at the moment )
      IF(ISETKW(35).EQ.0) THEN
        NOCSF = 1
        ISETKW(35) = 2
      END IF
* CSF expansion must only be used when two vectors are stored in CORE
      IF(NOCSF.EQ.0.AND.ICISTR.NE.1) THEN
        WRITE(luwrt,*)
     &  '   Do not use CSF expansions for blocked CI calculations ! '
        WRITE(luwrt,*)
     &  '   Employ the NOCSF option or set CISTOR to 1 ! '
        WRITE(luwrt,*)
     &  '   I will stop and wait for your decision  ! '
        Call Abend1( 11 )
      END IF
*
* 37 : Avoid any readin of integrals ( useful for obtaining
*      size of CI expansion etc.
*
      IF(ISETKW(37).EQ.0 ) THEN
        NOINT = 0
        ISETKW(37) = 2
      END IF
*
* 38 : Dump integrals in formatted form : Default is No
*
      IF(ISETKW(38).EQ.0) THEN
        IDMPIN = 0
        ISETKW(38) = 2
      END IF
*. If import is from LUCIA, dumping of integrals is disabled
*. Disabling is disabled (sic) : To allow for final integraltrans
C     IF(IDMPIN.EQ.1.AND.ENVIRO(1:5).EQ.'LUCIA') THEN
C       IDMPIN = 0
C       WRITE(6,*) 'Warning : Dump to LUCIA format disabled'
C       WRITE(6,*) '(input format is LUCIA !              )'
C       WRITE(6,*)
C       WRITE(6,*) '                     /Lucia            '
C     END IF

*
* 39 : Explicitly dimension of dimension of block of resolution strings
*
      IF(ISETKW(39).EQ.0) THEN
        MXINKA = 100
        ISETKW(39) = 2
      END IF
*
* 40 : Use CJKAKB intermediate matrices in alpha-beta loop,
*      Default is  YES !!!!!
*
      IF(ISETKW(40).EQ.0) THEN
        ICJKAIB = 1
        ISETKW(40) = 2
      END IF
*
*  41 : Initial CI in reference space, default is : No
*
      IF(ISETKW(41).EQ.0) THEN
         ISETKW(41) = 2
         INIREF = 0
      END IF
*
*  42 : Restart with CI in reference space
*
      IF(ISETKW(42).EQ.0) THEN
         ISETKW(42) = 2
         IRESTRF = 0
      END IF
*
*  44 : Use MOC method for alpha-alpha loop, default is NO !
*
      IF(ISETKW(44).EQ.0) THEN
         ISETKW(44) = 2
         MOCAA = 0
      END IF
*
*  45 : Use MOC method for alpha-beta loop, default is NO !
*
      IF(ISETKW(45).EQ.0) THEN
         ISETKW(45) = 2
         MOCAB = 0
      END IF
*
*  46 : Core energy : Default is 0 / MOLCAS : Value read in !
*
      IF(ISETKW(46).EQ.0) THEN
         ISETKW(46) = 2
         ECORE = 0.0D0
      END IF
*
*  47 : Use perturbation theory for zero order space . def is no !
*
      IF(ISETKW(47).EQ.0) THEN
        IPERT = 0
        NPERT = 0
        ISETKW(47) = 2
*. Else ensure that a CI in reference space is performed
      ELSE
        INIREF = 1
      END IF
*
*
*  48 : Approximate Hamiltonian in reference space : NO !!
*
      IF(ISETKW(48).EQ.0) THEN
        IAPRREF = 0
        MNRS1RE = MNRS1R
        MXRS3RE = MXRS3R
        ISETKW(48) = 2
      END IF
*
*  49 : Approximate Hamiltonian in zero order space : NO !!
*
      IF(ISETKW(49).EQ.0) THEN
        IAPRZER = 0
        MNRS1ZE = MNRS10
        MXRS3ZE = MXRS30
        ISETKW(49) = 2
      END IF
*
* 50 : GAS shells must be defined
*
      IF(ISETKW(50).EQ.0) THEN
        WRITE(luwrt,*) ' GASSH must be defined '
        NERROR = NERROR + 1
        IGSFILL = 0
        ISETKW(50) = -1
      END IF
*
* 52 : Combination of gasspaces : Default is just to take each  space
*      By itself
*
      IF(ISETKW(52).EQ.0) THEN
        NCMBSPC = NCISPC
        DO ICISPC = 1, NCISPC
          LCMBSPC(ICISPC) = 1
          ICMBSPC(1,ICISPC) = ICISPC
        END DO
        ISETKW(52) = 2
      END IF
*
* 53 : Convergence threshold for CI
*
      IF(ISETKW(53).EQ.0) THEN
!       THRES_E = 1.0D-10
        THRES_E = 1.0D-7
        ISETKW(53) = 2
      END IF
*
* 54 : General sequencer : default is just straight sequence
*      of CI with default number of iterations
      IF(ISETKW(54).EQ.0) THEN
        DO JCMBSPC = 1, NCMBSPC
          NSEQCI(JCMBSPC) = 1
          CSEQCI(1,JCMBSPC) = 'CI'
          ISEQCI(1,JCMBSPC) = MAXIT
        END DO
        ISETKW(54) = 2
      END IF
*
* 55 : EKT calculation : Default is no
*
      IF(ISETKW(55).EQ.0) THEN
        IEXTKOP = 0
        ISETKW(55) = 2
      END IF
*
* 57 : Allow first order correction to be saved on DISC
*     (For vector free calculations )
*     Default is : NO !!
      IF(ISETKW(57).EQ.0) THEN
        IC1DSC = 0
        ISETKW(57) = 2
      END IF
*
* 58 : Restrictions on interactions of perturbation
*
*. Default is : no
      IF(ISETKW(58).EQ.0) THEN
        IH0SPC = 0
        ISETKW(58) = 2
      END IF
*
* 59 : Type of perturbation in subspaces spaces
*
* Default is specified by IPART from keyword PERTU
      IF(ISETKW(59).EQ.0) THEN
       ISETKW(59) = 2
       IF(IH0SPC.NE.0) THEN
         DO JPTSPC = 1, NPTSPC
           IH0INSPC(JPTSPC) = IPART
         END DO
       END IF
      END IF
*
* 60 : Reference Root, default is NROOT
*
*. Should be less or equal to NROOT
      IF(ISETKW(60).EQ.1) THEN
        IF(IRFROOT.GT.NROOT) THEN
          WRITE(luwrt,*) ' Reference root (RFROOT) larger '
          WRITE(luwrt,*) ' than total number of roots (NROOT) '
          WRITE(luwrt,*) ' CHANGE NROOT or RFROOT '
          NMISS = NMISS + 1
        END IF
      END IF

      IF(ISETKW(60).EQ.0) THEN
       ISETKW(60) = 2
       IRFROOT = NROOT
      END IF
*
* 61 : H0EX : Orbital spaces in which exaxt Hamiltonian is used
*      No default
*.
*. Is H0EX required ( Has H0FRM = 5 been used )
      IUSED = 0
      IF(ISETKW(59).EQ.1) THEN
         IUSED = 0
         DO JPTSPC = 1, NPTSPC
           IF( IH0INSPC(JPTSPC) .EQ. 5 ) IUSED = 1
         END DO
       END IF
       IF(IUSED.EQ.0.AND.ISETKW(61).EQ.0) THEN
*. No exact spaces included and none have been defined !
       END IF
       IF(IUSED.EQ.1.AND.ISETKW(61).EQ.0) THEN
*. Needed, but not supplied
          WRITE(luwrt,*) ' You have specified that zero order operator'
          WRITE(luwrt,*) ' Include exact Hamilton operator in subspace'
          WRITE(luwrt,*) ' You should then also supply Keyword H0EX '
          NMISS = NMISS + 1
       END IF
*
*. If perturbation theory will be invoked be sure that the
*. form of perturbation theory has been specified through
* KEYWORD PERTU ( number 47 as you maybe know )
      IDOPERT = 0
      DO JCMBSPC = 1, NCMBSPC
        DO JSEQCI = 1, NSEQCI(JCMBSPC)
          IF(ISEQCI(JSEQCI,JCMBSPC).EQ.-5 ) IDOPERT = 1
        END DO
      END DO
*
      IF(IDOPERT.EQ.1 .AND. IPERT.EQ.0) THEN
        WRITE(luwrt,*) ' Perturbation theory will be used '
        WRITE(luwrt,*) ' Please specify form through PERTU keyword '
        NMISS = NMISS + 1
      END IF
*
*. 62 : Default Handling of degenrences of initial CI vectors
*.      Default is : No action
*
      IF(ISETKW(62).EQ.0) THEN
        INIDEG = 0
        ISETKW(62) = 2
      END IF
*
*. 63 : Use F + Lambda(H-F) as operator instead of H
*.      Default is : No i.e Lambda = 1
*
      IF(ISETKW(63).EQ.0) THEN
        XLAMBDA = 1.0D0
        ISETKW(63) = 2
      END IF
*
*. 64 : Smallest block in batch of C and sigma
*.      Default is zero
*
      IF(ISETKW(64).EQ.0) THEN
        LCSBLK = 0
        ISETKW(64) = 2
      END IF
*
*. 66 : xxx
*
*. 68 : Type of natural orbitals, default is natural orbitals
*
      IF(ISETKW(68).EQ.0) THEN
        ISETKW(68) = 2
        IFINMO = 1
      END IF
*
*. 69 : Default Threshold for individual energy correction = 0.0
*
      IF(ISETKW(69).EQ.0) THEN
        E_THRE = 0.0D0
        ISETKW(69) = 2
      END IF
*
*. 70 : Default Threshold for wave individual function corrections = 0.0
*
      IF(ISETKW(70).EQ.0) THEN
        C_THRE = 0.0D0
        ISETKW(70) = 2
      END IF
*
*. 71 : Default Threshold for total energy corrections = 0.0
*
      IF(ISETKW(71).EQ.0) THEN
        E_CONV = 0.0D0
        ISETKW(71) = 2
      END IF
*
*. 72 : Default Threshold for total wave function correction = 0.0
*
      IF(ISETKW(72).EQ.0) THEN
        C_CONV = 0.0D0
        ISETKW(72) = 2
      END IF
*
*. 73 : Perform Class selection : Default yes if TERACI is used
*
      IF(ISETKW(73).EQ.0) THEN
        IF(IDIAG.EQ.1) THEN
          ICLSSEL = 0
        ELSE IF (IDIAG.EQ.2) THEN
          ICLSSEL = 1
        END IF
        ISETKW(73) = 2
      END IF
*
* 74 : Calculation of density matrices : Default is
*       calculation of one-body density
*
      IF(ISETKW(74).EQ.0) THEN
        IDENSI = 1
        ISETKW(74) = 2
      END IF
*. If IDENSI was set to zero and properties were requested
*  overwrite input to obtain 1-el matrix
      IF(IDENSI.EQ.0.AND.ISETKW(80).EQ.1) THEN
        WRITE(luwrt,*) ' You have specified calculation of'
        WRITE(luwrt,*) ' one-electron properties, and this'
        WRITE(luwrt,*) ' requires the calculation of the '
        WRITE(luwrt,*) ' one-electron density. '
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' You have, however, specified IDENSI=0'
        WRITE(luwrt,*) ' corresponding  to no densities'
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' I will allow myself to modify your'
        WRITE(luwrt,*) ' input to allow calculation of the '
        WRITE(luwrt,*) ' one-electron densities, so property'
        WRITE(luwrt,*) ' calculation can proceed as planned '
        WRITE(luwrt,*)
        WRITE(luwrt,*)                        ' Lucia '
*. and do it
        IDENSI = 1
      END IF
*. If CC is performed, one- and two- particle densities are
*  used in current simple-minded implementation.
      IF(I_DO_CC .EQ. 1 .AND. IDENSI.LE.1 ) THEN
        IDENSI = 2
        WRITE(luwrt,*) ' IDENSI flag raised to two for CC calculation'
      END IF
*
* 75 : Perturbation expansion of EKT, default is no
*
      IF(ISETKW(75).EQ.0) THEN
        IPTEKT = 0
        ISETKW(75) = 2
      END IF
*
* 76 : Root for zero order operator , default is NROOT
*
*. Should be less or equal to NROOT
      IF(ISETKW(76).EQ.1) THEN
        IF(IH0ROOT.GT.NROOT) THEN
          WRITE(luwrt,*) ' Zero order operator root (H0ROOT) larger '
          WRITE(luwrt,*) ' than total number of roots (NROOT) '
          WRITE(luwrt,*) ' CHANGE NROOT or H0ROOT '
          NMISS = NMISS + 1
        END IF
      END IF
      IF(ISETKW(76).EQ.0) THEN
       ISETKW(76) = 2
       IH0ROOT = NROOT
      END IF
*
* 77 : NO restart from previous vectors in calculation 2
*      Deafault is NO NO, ie. restart in calc 2
*
      IF(ISETKW(77).EQ.0) THEN
        IRST2 = 1
        ISETKW(77) = 2
      END IF
*
* 78 : skip initial energy evaluations - if possible
*
      IF(ISETKW(78).EQ.0) THEN
        ISKIPEI = 1
        ISETKW(78) = 2
      END IF
*
* 79 : Symmetry of x,y,z - needed for property calculations
*
      IF(ISETKW(79).EQ.0) THEN
*. Problematic if Properties should be calculated
       IF(ISETKW(80).EQ.1.OR.ISETKW(81).EQ.1.OR.ISETKW(82).EQ.1)
     & THEN
         WRITE(luwrt,*) ' Symmetry of X,Y,Z has not been given'
         WRITE(luwrt,*) ' You have to specify this for property calc'
         WRITE(luwrt,*) ' Please add KEYWORD XYZSYM '
         NMISS = NMISS + 1
         ISETKW(79) = -1
       ELSE
*. Is not needed, just supply zeroes
         DO ICOMP = 1, 3
           IXYZSYM(ICOMP) = 0
         END DO
         ISETKW(79) = 2
       END IF
      END IF
*
* 80 : Property calculation, default is no
*
      IF(ISETKW(80).EQ.0) THEN
        NPROP = 0
        ISETKW(80) = 2
      END IF
*
* 81 : Transition properties , default is no
*
      IF(ISETKW(81).EQ.0) THEN
        ITRAPRP = 0
        ISETKW(81) = 2
      END IF
*
* 82 : Response properties , default is no
*
      IF(ISETKW(82).EQ.0) THEN
        IRESPONS = 0
        ISETKW(82) = 2
        NRESP = 0
        N_AVE_OP = 0
      END IF
*. Properties should be defined if transition properties are
*. invoked
      IF(ITRAPRP.NE.0.AND.NPROP.EQ.0) THEN
        WRITE(luwrt,*)
     &  ' You have specified transition property calculation'
        WRITE(luwrt,*)
     &  ' (keyword TRAPRP) but no property labels have been supplied'
        WRITE(luwrt,*)
     &  '(Keyword PROPER). Transition densities will be obtained '
      END IF
*
* 83 : Max number of iterations in linear equations
*
      IF(ISETKW(83).EQ.0) THEN
        MXITLE = 20
        ISETKW(83) = 2
      END IF
*
* 85 : Root homing, default is no
*
      IF(ISETKW(85).EQ.0) THEN
        IROOTHOMING = 0
        ISETKW(85) = 2
      END IF
*
* 86 : Particle hole simplifications, default is no
*
      IF(ISETKW(86).EQ.0) THEN
       IUSE_PH = 0
       ISETKW(86) = 2
      END IF
*
* 87 : Ask advice for route in sigma blocks, default is no
*      (It is said that programs reflects the minds of their creators)
*
      IF(ISETKW(87).EQ.0) THEN
       IADVICE = 0
       ISETKW(87) = 2
      END IF
*
* 88 : Transform CI vectors at end of each calculation
*      default is no
*
      IF(ISETKW(88).EQ.0) THEN
       ITRACI = 0
       ISETKW(88) = 2
       ITRACI_CR = 'undefine'
       ITRACI_CN = 'undefine'
      END IF
*
* 89 : Divide strings into active and passive parts
*      default is no
*
      IF(ISETKW(89).EQ.0) THEN
       IUSE_PA = 0
       ISETKW(89) = 2
      END IF
*
* 90 : Perturbation expansion of Fock matrix : default is no
*
      IF(ISETKW(90).EQ.0) THEN
       IPTFOCK = 0
       ISETKW(90) = 2
      END IF
*
* 91 : Print final CI vectors : default is no
*
      IF(ISETKW(91).EQ.0) THEN
       IPRNCIV = 0
       ISETKW(91) = 2
      END IF
*
* 92 : Restart CC calculation with coefs on LU_CCAMP
*
      IF(ISETKW(92).EQ.0) THEN
       I_RESTRT_CC = 0
       ISETKW(92) = 2
      END IF
*
* 93 : xxx
*
* 94 : xxx
*
* 95 : xxx
*
* Use relaxed densities for properties : default is no
*
      IF(ISETKW(96).EQ.0) THEN
        IRELAX = 0
        ISETKW(96) = 2
      END IF
*.
      IF(IRELAX.EQ.1) THEN
*. To obtain relaxed densities two-elec density must be calc, so
        IF(IDENSI.LT.2) THEN
          WRITE(luwrt,*) ' Density matrix flag (IDENSI) raised '
          WRITE(luwrt,*) ' to allow calculation of 2-elec densities'
          IDENSI = 2
        END IF
      END IF
*
* Expert mode ( neglect mistyped keywords ) : default is no expert
*
      IF(ISETKW(97).EQ.0) THEN
        IEXPERT = 0
        ISETKW(97) = 2
      END IF
*
*. Thresholds only active in connection with IDIAG = 2,
*. Check and maybe issue a warning
      IF(IDIAG.EQ.2) THEN
*. Check to ensure that zero or two thresholds were  set,
        IF(ISETKW(69).NE.ISETKW(70)) THEN
          WRITE(luwrt,*)
     &    ' Only a single threshold (E_THRE or C_THRE) '
          WRITE(luwrt,*)
     &    ' on individual determinants given. '
          WRITE(luwrt,*)
     &    ' One of the thresholds vanishes therefore and '
          WRITE(luwrt,*)
     &    ' all determinants will therefore be included  '
          WRITE(luwrt,*)
          WRITE(luwrt,*) '                   Warns '
          WRITE(luwrt,*)
          WRITE(luwrt,*) '                   LUCIA  '
        END IF
      ELSE
*. Good old diagonalization, thrsholds not active
        IF(ISETKW(69).EQ.1.OR.ISETKW(70).EQ.1) THEN
          WRITE(luwrt,*)
     &    ' Thresholds on selection of individual coefficients '
          WRITE(luwrt,*)
     &    ' are only active in connection with keyword TERACI  '
          WRITE(luwrt,*)
          WRITE(luwrt,*) '                   Warns '
          WRITE(luwrt,*)
          WRITE(luwrt,*) '                   LUCIA  '
        END IF
      END IF
*
      IF(ISETKW(50).EQ.1.AND. ISETKW(51).EQ.0) THEN
* Number of GAS shells given but no occupations !!
        WRITE(luwrt,*) ' GAS calculation (GASSH specified)'
        WRITE(luwrt,*) ' But no Occupation constraints (GASSPC) '
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' Please add GASSPC '
        NMISS = NMISS + 1
      END IF
*
*
*
      IF(NMISS.NE.0) THEN
        WRITE(luwrt,'(1X ,A,I9)')
     &  ' Number of missing required keyword ', NMISS
        WRITE(luwrt,*)
     &  ' You have wounded me I give up '
        WRITE(luwrt,*)
        WRITE(luwrt,*)
        WRITE(luwrt,*)
        WRITE(luwrt,*)
     & '     An expert is a man who has made all the mistakes,'
        WRITE(luwrt,*)
     &  '     which can be made, in a very narrow field        '
        WRITE(luwrt,*)
     &  '                                                      '
        WRITE(luwrt,*)
     &  '                                      Niels Bohr      '
        IF(IEXPERT.EQ.0) THEN
          Call Abend()
        ELSE
          WRITE(luwrt,*) ' Processing continues (EXPERT mode )'
        END IF
      END IF
*. Open one-electron file to obtain core energy and
*. Number of MO's and AO's
      IF(NOINT.EQ.0) THEN
!       CALL GET_ORB_DIM_ENV(ECORE_ENV)
!       IF(ISETKW(46).EQ.2) ECORE = ECORE_ENV
#ifdef PRG_DIRAC
        if(ENVIRO(1:5).ne.'DIRAC')
     &    CALL CHK_ORBDIM(IGSFILL,NIRREP)
#else
        if(ENVIRO(1:5).ne.'DALTO')
     &    CALL CHK_ORBDIM(IGSFILL,NIRREP)
#endif
      ELSE
        WRITE(luwrt,*) ' GETOBS not called since NOINT = 1'
        ECORE = 0.0D0
      END IF
*. Check number of orbitals and insert occupations for ALL/REST

************************************************************
*                                                          *
* Part 3 : Print input                                     *
*                                                          *
************************************************************
*
c     WRITE(luwrt,*)
c     WRITE(luwrt,*) '******************'
c     WRITE(luwrt,*) '*  Title of run  *'
c     WRITE(luwrt,*) '******************'
c     WRITE(luwrt,*)
c     CALL PRTITL(TITLEC)
c     WRITE(luwrt,*)
*
*. Type of reference state
      WRITE(luwrt,*)
      WRITE(luwrt,*) '*************************************'
      WRITE(luwrt,*) '*  Symmetry and spin of CI vectors  *'
      WRITE(luwrt,*) '*************************************'
      WRITE(luwrt,*)
*. Point group
      IF(PNTGRP.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Point group ............ D2H'
        if (ENVIRO(1:6).eq.'DIRAC ') then
          WRITE(luwrt,'(1X ,A,A)')
     &  '       Using subgroup ......... ',DOUGRP
          WRITE(luwrt,'(1X ,A,A)')
     &  '         (D2 equiv. C2v, C2 eq. Cs) '
        end if
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Point group ............ C inf v'
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Point group ............ D inf h'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Point group ............ O3'
      END IF
*.Spatial symmetry
      IF(PNTGRP.EQ.1) THEN
        WRITE(luwrt,'(1X ,A,I1)')
     &  '     Spatial symmetry ....... ', IREFSM
      ELSE IF(PNTGRP.EQ.2) THEN
        WRITE(luwrt,'(1X ,A,I1)')
     &  '     ML value ............... ', IREFML
      ELSE IF(PNTGRP.EQ.3) THEN
        WRITE(luwrt,'(1X ,A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(luwrt,'(1X ,A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(luwrt,'(1X ,A)')
     &  '     Parity   ..............  Ungerade'
      ELSE IF(PNTGRP.EQ.4) THEN
        WRITE(luwrt,'(1X ,A,I1)')
     &  '     L  value ............... ', IREFL
        WRITE(luwrt,'(1X ,A,I1)')
     &  '     ML value ............... ', IREFML
        IF(IREFPA.EQ.1) WRITE(luwrt,'(1X ,A)')
     &  '     Parity   ..............  Gerade'
        IF(IREFPA.EQ.2) WRITE(luwrt,'(1X ,A)')
     &  '     Parity   ..............  Ungerade'
      END IF
*.Spin
      WRITE(luwrt,'(1X ,A,I4)')
     &  '     2 times spinprojection', MS2
      WRITE(luwrt,'(1X ,A,I4)')
     &  '     Spin multiplicity ....', MULTS
*.Number of active electrons
      WRITE(luwrt,'(1X ,A,I4)')
     &  '     Active electrons .....', LUCI_NACTEL
      WRITE(luwrt,*)
      WRITE(luwrt,*) '*********************************************'
      WRITE(luwrt,*) '*  Shell spaces and occupation constraints  *'
      WRITE(luwrt,*) '********************************************* '
      WRITE(luwrt,*)
*
      IF(IDOGAS.EQ.0) THEN
*. Kept because outpu can lated be used for GAS
*
*. NOT a GAS expansion
*
*
      WRITE(luwrt,'(1X ,A,10I4)')
     &  '                Irrep ',(I,I = 1,NIRREP)
      WRITE(luwrt,'(1X ,A,2X,10A)')
     &  '                ===== ',('====',I = 1,NIRREP)
*
*. Inactive
      IF(ISETKW(7).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Inactive         ',(NINASH(I),I=1,NIRREP)
      END IF
*. Core
      IF(ISETKW(8).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Core             ',(NRS0SH(I,1),I=1,NIRREP)
      END IF
*. RAS1
      IF(ISETKW(9).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Ras1             ',(NRSSH(I,1),I=1,NIRREP)
      END IF
*. RAS2/ACTIVE
      IF(ISETKW(10).EQ.1) THEN
        IF(INTSPC.EQ.1) THEN
          WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Active           ',(NRSSH(I,2),I=1,NIRREP)
        ELSE IF(INTSPC.EQ.2) THEN
          WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Ras2             ',(NRSSH(I,2),I=1,NIRREP)
        END IF
      END IF
*. RAS3
      IF(ISETKW(11).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,10I4)')
     &  '     Ras3             ',(NRSSH(I,3),I=1,NIRREP)
      END IF
      IF(INTSPC.EQ.2.AND.IMLCR3.EQ.1) WRITE(luwrt,'(1X ,A)')
     &  '     ( RAS 3 space supplied by courtesy of TRAONE )'
*. Secondary space
      IF(ISETKW(13).EQ.1) THEN
        DO 310 ITP = 1,MXR4TP
          WRITE(luwrt,'(A,I2,A,10I4)')
     &  '      Secondary',ITP,'      ',(NRS4SH(I,ITP),I=1,NIRREP)
  310   CONTINUE
      END IF
*. Deleted space
      IF(ISETKW(26).EQ.1) THEN
        WRITE(luwrt,'(A,10I4)')
     &  '      Deleted          ',(NDELSH(I),I=1,NIRREP)
      END IF
      IF(IMLCR3.EQ.2) WRITE(luwrt,'(A)')
     &  '      ( Deleted shells supplied by courtesy of TRAONE )'
*.Core space
      WRITE(luwrt,*)
      IF(ISETKW(8).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,I4)')
     &  '     Largest number of excitations out of core ..... ',MXHR0
      END IF
*.Secondary space
      IF(ISETKW(13).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,I4)')
     &  '     Largest number of excitations to secondary space',MXER4
      END IF
      ELSE
*
*. GAS space
*
      WRITE(luwrt,*)
      WRITE(luwrt,*) ' *************************'
      WRITE(luwrt,*) ' Generalized active space '
      WRITE(luwrt,*) ' *************************'
      WRITE(luwrt,*)
      WRITE(luwrt,'(A)') ' Orbital subspaces:'
      WRITE(luwrt,'(A)') ' ================== '
      WRITE(luwrt,*)
      WRITE(luwrt,'(A,20I4)')
     &  '                 Irrep ',(I,I = 1,NIRREP)
      WRITE(luwrt,'(A,2X,20A)')
     &  '                 ===== ',('====',I = 1,NIRREP)
      DO IGAS = 1, NGAS
        WRITE(luwrt,'(A,I2,A,10I4,6X,2I6)')
     &  '        GAS',IGAS,'          ',
     &  (NGSSH(IRREP,IGAS),IRREP = 1, NIRREP)
      END DO
      IF(IGSFILL.NE.0) WRITE(luwrt,'(7X,A,I3)')
     &' Gas space provided by courtesy of LUCIA :',  IGSFILL
*
      WRITE(luwrt,*)
      WRITE(luwrt,*)  '*******************'
      WRITE(luwrt,*)  ' Occupation spaces '
      WRITE(luwrt,*)  '*******************'
      WRITE(luwrt,*)
      WRITE(luwrt,'(A,I3)')
     & ' Number of Occupation spaces : ',NCISPC
      WRITE(luwrt,*)
      DO ICISPC = 1, NCISPC
      WRITE(luwrt,'(A,I3)')
     &' Bounds on accumulated occupations for space : ',ICISPC
      WRITE(luwrt,'(A)')
     & ' ====================================================== '
      WRITE(luwrt,'(A)')
      WRITE(luwrt,'(A)') '         Min. occ    Max. occ '
      WRITE(luwrt,'(A)') '         ========    ======== '
      DO IGAS = 1, NGAS
        WRITE(luwrt,'(A,I2,3X,I3,9X,I3)')
     &  '   GAS',IGAS,IGSOCCX(IGAS,1,ICISPC),IGSOCCX(IGAS,2,ICISPC)
      END DO
      END DO
*
      WRITE(luwrt,*)
      WRITE(luwrt,*)
     &' **************************************************'
      WRITE(luwrt,*)
     &' Specification of CI Spaces (combinations of above)'
      WRITE(luwrt,*)
     &' **************************************************'
      WRITE(luwrt,*)

      WRITE(luwrt,*)
      WRITE(luwrt,'(A,I3)')
     &' Number of CI spaces included : ', NCMBSPC
      WRITE(luwrt,*)
      DO JCMBSPC = 1, NCMBSPC
        WRITE(luwrt,*)
        WRITE(luwrt,'(A,I3)') ' Information about CI space ', JCMBSPC
        WRITE(luwrt,'(A)')    ' =================================='
        WRITE(luwrt,'(1X ,3X,A,I3)')
     &  'Number of occupation spaces included  ',LCMBSPC(JCMBSPC)
        WRITE(luwrt,'(A,10I3)') '    Occupation spaces included ',
     &  (ICMBSPC(II,JCMBSPC),II=1,LCMBSPC(JCMBSPC))
*
      END DO
*
      WRITE(luwrt,*)
      WRITE(luwrt,*) ' ******************************************'
      WRITE(luwrt,*) ' Specification of Sequence of calculations '
      WRITE(luwrt,*) ' ******************************************'
      WRITE(luwrt,*)
      DO JCMBSPC = 1, NCMBSPC
        WRITE(luwrt,*)
        WRITE(luwrt,'(7X,A,I3)') ' CI space ', JCMBSPC
        WRITE(luwrt,'(7X,A)')    ' =============='
        WRITE(luwrt,*)
*
C       WRITE(6,'(A,I3)') ' Number of calculations in this CI space ',
C    &  NSEQCI(JCMBSPC)
C       WRITE(6,'(A)')   '  Calculations in this subspace '
C       WRITE(6,'(A)')   '  =============================='
        DO JSEQ = 1, NSEQCI(JCMBSPC)
          CARDX = CSEQCI(JSEQ,JCMBSPC)
          IF(CARDX(1:7).EQ.'VECFREE') THEN
            WRITE(luwrt,'(10X,A,I3)')
     &      '       Vector free calculation at level ',
     &      -ISEQCI(JSEQ,JCMBSPC)
          ELSE IF(CARDX(1:2).EQ.'CI') THEN
            WRITE(luwrt,'(10X,A,I3)')
     &      '       Normal CI with max. iterations = ',
     &      ISEQCI(JSEQ,JCMBSPC)
          ELSE IF(CARDX(1:6).EQ.'APR-CI') THEN
            WRITE(luwrt,'(10X,A,I3)')
     &      '       CI using approximate H with max. iterations = ',
     &      ISEQCI(JSEQ,JCMBSPC)
          ELSE IF(CARDX(1:5).EQ.'PERTU') THEN
            WRITE(luwrt,'(10X,A,I3)')
     &      '       Perturbation calculation '
          ELSE IF(CARDX(1:2).EQ.'CC'   ) THEN
            WRITE(luwrt,'(10X,A,I3)')
     &      '       Coupled Cluster Calculation, max. iterations =',
     &      ISEQCI(JSEQ,JCMBSPC)
            IF(JCMBSPC.EQ.LAST_CC_SPC.AND.JSEQ.EQ.LAST_CC_RUN) THEN
              WRITE(luwrt,'(10X,A)')
     &      '       (Expanded cc wf will be transferred to LUC ) '
            END IF
          ELSE IF(CARDX(1:4).EQ.'ICCI' ) THEN
            WRITE(luwrt,'(10X,A)')
     &      '       Internal Contracted calculation '

          END IF
        END DO
*       ^ End of loop over spaces in given CI space
      END DO
*     ^ End of loop over CI spaces
      WRITE(luwrt,*)
      END IF
*     ^ End of GAS/NOGAS switch
*
      IF(XLAMBDA.NE.1.0D0) THEN
        WRITE(luwrt,*)
        WRITE(luwrt,'(A,F13.8)')
     &  ' Modifed operator H(l) = l*F + l*(H-F) used with l =',XLAMBDA
        IF(IRESTR.EQ.0) THEN
        WRITE(luwrt,*)
     &  ' Notice : This madness starts  in second calculation'
        ELSE
         WRITE(luwrt,*) 
     & ' You have specified a calculation with modified '
         WRITE(luwrt,*) ' Hamiltonian (the LAMBDA option) and RESTART '
         WRITE(luwrt,*) ' so this is what I will do '
         WRITE(luwrt,*)
         WRITE(luwrt,*) 
     & '   1:) Perform CI in space 1 to obtain Hamiltonian'
         WRITE(luwrt,*) '       (no RESTART in this space )'
         WRITE(luwrt,*) '   2:) CI calculation in space 2  with '
         WRITE(luwrt,*) 
     & '       modified Hamiltonian and RESTART from LU21'
         WRITE(luwrt,*) ' Space 2 should therefore correspond to the'
         WRITE(luwrt,*) ' restarted calculation '
       END IF
      END IF
*
      WRITE(luwrt,*)
      WRITE(luwrt,*) '***********'
      WRITE(luwrt,*) '*  Roots  *'
      WRITE(luwrt,*) '*********** '
      WRITE(luwrt,*)
      WRITE(luwrt,'(1X ,A,I3)')
     &  '     Number of roots to be obtained ', NROOT
      WRITE(luwrt,'(1X ,A,(20I3))')
     &  '     Roots to be obtained ', (IROOT(I),I=1, NROOT )
*
      WRITE(luwrt,*)
      WRITE(luwrt,*) '**************************'
      WRITE(luwrt,*) '*  Run time definitions  *'
      WRITE(luwrt,*) '************************** '
      WRITE(luwrt,*)
*. Program environment
      WRITE(luwrt,'(A,A6)') '      Program environment... ', ENVIRO
*. Integral import
      IF(NOINT.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     No integrals will be read in       '
      ELSE IF(NOINT.EQ.0) THEN
      IF(INTIMP.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Integrals read in in MOLCAS format '
      ELSE IF(INTIMP.EQ.5) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Integrals read in in SIRIUS format '
      ELSE IF(INTIMP.EQ.2) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Integrals read in in LUCAS format '
      ELSE IF(INTIMP.EQ.3) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Integrals read in in formatted form (E22.15) ',
     &  '      From unit 13'
        WRITE(luwrt,'(1X ,A)')
     &  '     All integrals of correct symmetry combination read in'
      END IF
*. Integral storage
      IF(INCORE.EQ.1) WRITE(luwrt,'(1X ,A)')
     &  '     All integrals stored in core'
      END IF
      WRITE(luwrt,*)
* ( END IF for NOINT
*. CSF or SD expansion
      IF(NOCSF.EQ.0) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     CI optimization performed with CSF''s '
      ELSE
        WRITE(luwrt,'(1X ,A)')
     &  '     CI optimization performed with SD''s '
      END IF
*. Ms,Ml combinations
      IF(ISETKW(27).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,F8.3)')
     &  '     Spin combinations used with sign ',PSSIGN
      END IF
      IF(ISETKW(28).EQ.1) THEN
        WRITE(luwrt,'(1X ,A,F8.3)')
     &  '     ML   combinations used with sign ',PLSIGN
      END IF
*. Initial approximation to vectors
      WRITE(luwrt,*)
      IF(IRESTR.EQ.1.AND.IRESTRF.EQ.0) THEN
         WRITE(luwrt,'(1X ,A)')
     &  '     Restarted calculation '
      ELSE IF(IRESTRF.EQ.1) THEN
         WRITE(luwrt,'(1X ,A)')
     &  '     Restarted calculation from REFERENCE space expansion'
      ELSE
         IF(MXP1.NE.0) THEN
           WRITE(luwrt,'(1X ,A)')
     &  '     Initial vectors obtained from explicit Hamiltonian'
         ELSE IF(MXP1.EQ.0) THEN
           WRITE(luwrt,'(1X ,A)')
     &  '     Initial vectors obtained from diagonal'
         END IF
      END IF
      IF(I_RESTRT_CC.EQ.1) THEN
           WRITE(luwrt,'(1X ,A)') '     CC calculation restarted '
      END IF
*. Handling of degenerencies of initial vectors
      IF(INIDEG.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Symmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.-1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Antiymmetric combination of degenerate initial vectors'
      ELSE IF (INIDEG.EQ.0) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     No combination of degenerate initial vectors'
      END IF
*. Ms,Ml combinations
C     IF(ISETKW(27).EQ.1) THEN
C       WRITE(luwrt,'(1X ,A,F8.3)')
C    &  '     Spin combinations used with sign ',PSSIGN
C     END IF
C     IF(ISETKW(28).EQ.1) THEN
C       WRITE(luwrt,'(1X ,A,F8.3)')
C    &  '     ML   combinations used with sign ',PLSIGN
C     END IF
*. CI storage mode
      WRITE(luwrt,*)
      IF(ICISTR.EQ.1) THEN
        WRITE(luwrt,*)
     &  '     3 symmetry blocks and two vectors will be held in core '
      ELSE IF( ICISTR.EQ.2) THEN
        WRITE(luwrt,*)
     &  '     3 symmetry blocks will be held in core '
      ELSE IF( ICISTR.EQ.3) THEN
        WRITE(luwrt,*)
     &  '     3 symmetry-occ-occ blocks will be held in core '
      END IF
      IF(LCSBLK.NE.0) WRITE(luwrt,'(A,I10)')
     &  '      Smallest allowed size of sigma- and C-batch:',LCSBLK
      WRITE(luwrt,'(A,I10)')
     &  '      Dimension of block of resolution strings   :',MXINKA
      IF(IUSE_PH.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Particle-hole separation used '
      ELSE
        WRITE(luwrt,'(1X ,A)')
     & '      Particle-hole separation not used '
      END IF
*
      IF(IADVICE.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Advice routine call to optimize sigma generation'
      END IF
*
      IF(IUSE_PA.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '     Strings divided into active and passive parts'
      ELSE
        WRITE(luwrt,'(1X ,A)')
     &  '     Strings not divided into active and passive parts'
      END IF
*
      IF(IDENSI.EQ.0) THEN
        WRITE(luwrt,'(1X ,/A)')
     &  '      No calculation of density matrices'
      ELSE IF(IDENSI.EQ.1) THEN
        WRITE(luwrt,'(1X ,/A)')
     &  '     One-body density matrix calculated'
      ELSE IF(IDENSI.EQ.2) THEN
        WRITE(luwrt,'(1X ,/A)')
     &  '     One- and two-body density matrices  calculated'
      END IF
      WRITE(luwrt,*)
*
*. Diagonalization information
      WRITE(luwrt,'(1X ,A)')
     &  '     CI diagonalization : '
      WRITE(luwrt,'(1X ,A)')
     &  '     ==================== '
*. Subspace Hamiltinian
      IF(MXP1+MXP2+MXQ .EQ.0) THEN
        WRITE(luwrt,'(1X ,A)')
     &  '        No subspace Hamiltonian '
      ELSE
        WRITE(luwrt,'(1X ,A,3I4)')
     &  '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
      END IF
*. Diagonalizer
      IF(IDIAG.EQ.1.AND.ICISTR.EQ.1) THEN
        WRITE(luwrt,'(1X ,A)')
     &    '        Diagonalizer : MINDV4 '
      ELSE IF(IDIAG.EQ.1.AND.ICISTR.GE.2) THEN
        WRITE(luwrt,'(1X ,A)')
     &    '        Diagonalizer : MICDV* '
      ELSE IF(IDIAG.EQ.2) THEN
      WRITE(luwrt,'(1X ,A)')
     &    '        Diagonalizer : PICO*  '
      END IF
*. Root homing
      IF(IROOTHOMING.EQ.1) THEN
      WRITE(luwrt,'(1X ,A)')
     &  '        Root homing will be used '
      ELSE
      WRITE(luwrt,'(1X ,A)')
     &  '        No root homing '
      END IF
*. No restart in CI calc 2
      IF(IRST2.EQ.0) THEN
      WRITE(luwrt,'(1X ,A)')
     &  '        No restart from previous vectors in second calc '
      END IF
      IF(ISKIPEI.EQ.1) THEN
      WRITE(luwrt,'(1X ,A)')
     &  '        Initial energy evaluations skipped after first calc'
      WRITE(luwrt,'(1X ,A)')
     &  '        (Only active in connection with TERACI )'
      END IF
*. Number of iterations
C     WRITE(luwrt,'(A,I5)')
C    &  '         Allowed number of iterations    ',MAXIT
*. Number of CI vectors in subspace
      WRITE(luwrt,'(A,I5)')
     &  '         Allowed Dimension of CI subspace',MXCIV
*
      WRITE(luwrt,'(A,1P,E11.2)')
     &  '         Convergence threshold for energy',THRES_E
*
      IF(IDIAG.EQ.2) THEN
        WRITE(luwrt,'(1X ,/A,E11.5)')
     &   '        Individual second order energy threshold ',E_THRE
        WRITE(luwrt,'(1X ,A,E11.5)')
     &   '        Individual first order wavefunction threshold ',C_THRE
        IF(ICLSSEL.EQ.1) THEN
         WRITE(luwrt,*)
         WRITE(luwrt,'(1X ,A)')
     &   '         Class selection will be performed : '
         WRITE(luwrt,'(1X ,A)')
     &   '         =================================== '
         WRITE(luwrt,'(1X ,A,E11.5)')
     &    '          Total second order energy threshold ',E_CONV
         WRITE(luwrt,'(1X ,A,E11.5)')
     &    '          Total first order wavefunction threshold ',C_CONV
        ELSE
         WRITE(luwrt,'(1X ,A)')
     &'            No class selection in iterative procedure '
        END IF
      END IF
*
      IF(IPERT.NE.0) THEN
        WRITE(luwrt,'(1X ,A)')
     &    '     Perturbation calculation'
        WRITE(luwrt,'(1X ,A)')
     &  '     ======================= '
        WRITE(luwrt,*)
     &  '        Root Choosen as zero order state ', IRFROOT
        WRITE(luwrt,*)
     &  '        Root used for zero order operator ', IH0ROOT

        IF(IE0AVEX.EQ.1) THEN
          WRITE(luwrt,*)
     &  '        Expectation value of H0 used as zero order energy '
        ELSE IF( IE0AVEX.EQ.2) THEN
          WRITE(luwrt,*)
     &  '        Exact energy of reference used as zero order energy'
        END IF
        WRITE(luwrt,*)
     &  '        Correction vectors obtained through  order ', NPERT
        IF(IH0SPC.EQ.0) THEN
        WRITE(luwrt,*)
     &  '        No restrictions on perturbation interactions '
        ELSE
        WRITE(luwrt,*)
     &  '        Perturbation restricted to interactions in subspaces'
        END IF
*
        IF(IH0SPC.NE.0) THEN
        WRITE(luwrt,*)
     &  '        Number of perturbation subspaces ', NPTSPC
        WRITE(luwrt,*)
        WRITE(luwrt,*)
     &  '        ======================== '
        WRITE(luwrt,*)
     &  '        Perturbation subspaces : '
        WRITE(luwrt,*)
     &  '        ======================== '
        DO JPTSPC = 1, NPTSPC
          WRITE(luwrt,'(7X,/A)') '         Min. occ    Max. occ '
          WRITE(luwrt,'(7X,A)')  '         ========    ======== '
          DO IGAS = 1, NGAS
            WRITE(luwrt,'(7X,A,I2,3X,I3,9X,I3)')
     &      '   GAS',IGAS,IOCPTSPC(1,IGAS,JPTSPC)
     &                   ,IOCPTSPC(2,IGAS,JPTSPC)
          END DO
        END DO
*
       WRITE(luwrt,'(7X,/A)') '========================================'
       WRITE(luwrt,'(7X,A)') 'Approximate Hamiltonian in CI subspaces '
       WRITE(luwrt,'(7X,A)') '========================================'
       WRITE(luwrt,'(7X,A)')
       WRITE(luwrt,'(7X,A)') '    Subspace          H(apr)   '
       WRITE(luwrt,'(7X,A)') '  ============================='
       WRITE(luwrt,'(7X,A)')
       DO JPTSPC = 1, NPTSPC
         IF(IH0INSPC(JPTSPC).EQ.1) THEN
           WRITE(luwrt,'(12X,I3,8X,A)')
     &     JPTSPC, ' Diagonal Fock operator '
         ELSE IF(IH0INSPC(JPTSPC).EQ.2) THEN
           WRITE(luwrt,'(12X,I3,8X,A)')
     &     JPTSPC, ' Epstein-Nesbet operator'
         ELSE IF(IH0INSPC(JPTSPC).EQ.3) THEN
           WRITE(luwrt,'(12X,I3,8X,A)')
     &     JPTSPC, ' Nondiagonal Fock operator '
         ELSE IF(IH0INSPC(JPTSPC).EQ.4) THEN
           WRITE(luwrt,'(12X,I3,8X,A)')
     &     JPTSPC, ' Complete Hamiltonian  '
         ELSE IF(IH0INSPC(JPTSPC).EQ.5) THEN
           WRITE(luwrt,'(12X,I3,8X,A)')
     &     JPTSPC, ' Mix of Fock and Exact operator '
         END IF
        END DO
*
       END IF
       END IF
*
       IF(NPROP.ne.0) THEN
         WRITE(luwrt,'(7X,A,I3)')
     &   ' Number of properties to be calculated', NPROP
         WRITE(luwrt,*)
         WRITE(luwrt,'(9X,A)')    ' Properties : '
         WRITE(luwrt,'(9X,A)')   ' ============='
         DO IPROP = 1, NPROP
           WRITE(luwrt,'(16X,A)') PROPER(IPROP)
         END DO
*
         IF(IRELAX.EQ.0) THEN
           WRITE(luwrt,'(7X,A)') ' No use of relaxed densities '
         ELSE
           WRITE(luwrt,'(7X,A)') ' Relaxed densities used for '
           WRITE(luwrt,'(7X,A)') ' (implemented only for pert) '
         END IF
       END IF
*
       IF(IEXTKOP.EQ.0.AND.IPTEKT.EQ.0) THEN
C        WRITE(6,'(5X,A)') ' No extended Koopmans'' calculations '
       ELSE IF(IEXTKOP.NE.0) THEN
         WRITE(luwrt,'(5X,A)') ' Extended Koopmans'' calculations '
       ELSE IF(IPTEKT.NE.0) THEN
         WRITE(luwrt,'(5X,A)') 
     & ' Perturbation expansion of EKT equations'
       END IF
*
       IF(IPTFOCK.EQ.1) THEN
         WRITE(luwrt,*) ' Perturbation expansion of Fock matrix '
       ELSE
C        WRITE(6,*) 'No  Perturbation expansion of Fock matrix '
       END IF
*
      IF(ITRAPRP.EQ.0) THEN
C       WRITE(6,*)
C       WRITE(6,'(5X,A)')
C    &  ' No transition properties will be calculated'
      ELSE
        WRITE(luwrt,*)
        WRITE(luwrt,'(5X,A)')
     &  ' Transition properties will be calculated '
        WRITE(luwrt,*)  ' Symmetry of additional states :', IEXCSYM
        WRITE(luwrt,*)  ' Number   of additional states :', NEXCSTATE
        WRITE(luwrt,*)
      END IF
      IF(IRESPONS.NE.0) THEN
      WRITE(luwrt,*)
      WRITE(luwrt,*) '**************************'
      WRITE(luwrt,*) '*  Response Calculation  *'
      WRITE(luwrt,*) '************************** '
      WRITE(luwrt,*)
        WRITE(luwrt,*)
     &  ' CI-Response will be called after each CI calculation'
        WRITE(luwrt,*)
     &  ' Root used for response calculations (RFROOT) ',IRFROOT
        WRITE(luwrt,*)
C       WRITE(6,*) ' Number of A-operators : ', N_AVE_OP
        WRITE(luwrt,*) ' Labels of A-operators '
        WRITE(luwrt,*) ' ======================='
        WRITE(luwrt,*)
        DO IAVE = 1, N_AVE_OP
          WRITE(luwrt,'(1X , 6X,A)') AVE_OP(IAVE)
        END DO
        WRITE(luwrt,*)
C       WRITE(6,*) ' Number of response calculations ', NRESP
        WRITE(luwrt,*) ' Perturbations : '
        WRITE(luwrt,*) ' ================'
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' Calc  Op1    Op2    Mxord1     Mxord2    Freq '
        DO IRESP = 1, NRESP
          WRITE(luwrt,'(1X ,I2,2X,A,A,3X,I4,3X,I4,2X,F12.7)' )
     &    IRESP,RESP_OP(1,IRESP),RESP_OP(2,IRESP),MAXORD_OP(1,IRESP),
     &    MAXORD_OP(2,IRESP),RESP_W(IRESP)
        END DO
      END IF
*
        WRITE(luwrt,*)
        WRITE(luwrt,'(5X,A)') ' Final orbitals :'
        WRITE(luwrt,'(5X,A)') ' ================'
        IF(IFINMO.EQ.1) THEN
          WRITE(luwrt,'(8X,A)') ' Natural orbitals'
        ELSE IF(IFINMO.EQ.2) THEN
          WRITE(luwrt,'(8X,A)') ' Canonical orbitals'
        ELSE IF(IFINMO.EQ.3) THEN
          WRITE(luwrt,'(8X,A)') ' Pseudo-natural orbitals'
          WRITE(luwrt,'(8X,A)')
     &   ' (Density matrix diagonalized in orbital subspaces )'
        ELSE IF(IFINMO.EQ.4) THEN
          WRITE(luwrt,'(8X,A)') ' Pseudo-canonical orbitals'
          WRITE(luwrt,'(8X,A)')
     &   ' (FI+FA  diagonalized in orbital subspaces )'
         ELSE IF (IFINMO .EQ. 5 ) THEN
          WRITE(luwrt,'(8X,A)')
     &   ' Pseudo-natural-canonical orbitals (sic)'
          WRITE(luwrt,'(8X,A)')
     &   ' (Pseudo natural orbitals are first obtained'
          WRITE(luwrt,'(8X,A)')
     &   '  by diagonalizing density matrix in orbital subpspaces.'
          WRITE(luwrt,'(8X,A)')
     &   '  FI+FA is transformed to this basis, and the transformed'
          WRITE(luwrt,'(8X,A/)')
     &   '  matrix is block diagonalized) '
          WRITE(luwrt,'(8X,A)')
     &   ' Orbital spaces in which transformed FIFA is diagonalized'
          WRITE(luwrt,'(8X,A)')
     &   ' ========================================================'
          DO IPSSPC = 1, NPSSPC
            WRITE(luwrt,'(A,I2,A,10I4,6X,2I6)')
     &      '     SPACE',IPSSPC,'          ',
     &     (NPSSH(IRREP,IPSSPC),IRREP = 1, NIRREP)
          END DO
        END IF
*. Transformation of CI vectors
      IF(ITRACI.EQ.0) THEN
C       WRITE(6,'(5X,A)')  ' No transformation of CI vectors'
      ELSE
        WRITE(luwrt,'(5X,A)')   ' CI vectors transformed in each run'
        WRITE(luwrt,'(7X,A,A)')
     &        ' Complete or restricted rotations :',ITRACI_CR
        WRITE(luwrt,'(7X,A,A)')
     &        ' Type of Final orbitals           :',ITRACI_CN
      END IF
*
*.Print levels
*
      WRITE(luwrt,'(/A/A)')  '      Print levels : ',
     &                       '      =============='
      IF(ISETKW(29).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for string    information = ', IPRSTR
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for string    information = ', IPRSTR
      END IF
      IF(ISETKW(30).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for CI space  information = ', IPRCIX
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for CI space  information = ', IPRCIX
      END IF
      IF(ISETKW(31).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for orbital   information = ', IPRORB
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for orbital   information = ', IPRORB
      END IF
      IF(ISETKW(65).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for density matrix        = ', IPRDEN
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for density matrix        = ', IPRDEN
      END IF
      IF(ISETKW(32).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for iterative information = ', IPRDIA
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for iterative information = ', IPRDIA
      END IF
      IF(ISETKW(36).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for External blocks       = ', IPRXT
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for External blocks       = ', IPRXT
      END IF
*
      IF(IRESPONS.NE.0) THEN
      IF(ISETKW(84).EQ.2) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Default print level for response section      = ', IPRRSP
      ELSE
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Raised  print level for response section      = ', IPRRSP
      END IF
      END IF
*
      IF(IPROCC.NE.0) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Occupation of SD''s/ configurations           = ',IPROCC
      END IF
*
      IF(IPRNCIV.EQ.1 ) THEN
        WRITE(luwrt,'(3X ,A,I3)')
     &  '      Final CI vectors will be printed '
      END IF
*
      IF(IDMPIN.EQ.1) THEN
        WRITE(luwrt,'(3X ,A)')
        WRITE(luwrt,*) 
     & '      Integrals written in formatted form (E22.15)'
        WRITE(luwrt,*) '      on file 90 '
      END IF
      END
