#define SRDFT_DEBUG -1
#define LDAembed
!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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


C*****************************************************************************
      SUBROUTINE SRFUN_INPUT(SRFTMP)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Sep. 2003
C
C     Purpose :
C     Process input from .SRFUN keyword in *SCF INPUT
C     which defines the Exchange and Correlation functionals
C     to be used for the DFT-hybrid model in question.
C
C*****************************************************************************
      implicit none
#include "priunit.h"
#include "dftcom.h"
C     variables for SRCPBERI
      real*8 :: mu_0, mu_0_m6, mu_0_m7, mu_0_m4
      common /cb_srcpbe_RI/ mu_0, mu_0_m6, mu_0_m7, mu_0_m4
C     general variables
      character, intent(in) :: SRFTMP*80
      logical :: SETXFUN, SETCFUN
      integer, parameter :: NXTABLE=21, NCTABLE=36, NXCTABLE=6, maxfld=3
      integer :: lenfld(maxfld), I, J, IXWINT, ICWINT, nfld
      character(len=23) :: XTABLE(NXTABLE), CTABLE(NCTABLE),
     &                     XCTABLE(NXCTABLE), SRXCFUN, fields(maxfld)
      DATA XCTABLE/'SRLDA                  ','SRPBEGWS               ', ! 1:2
     &             'SRPBERI                ','SRPBE0GWS              ', ! 3:4
     &             'LRCLDA                 ','LRCPBEGWS              '/ ! 5:6
      DATA XTABLE /'NULL                   ','SRXLDA                 ', ! 1:2
     &             'SRXGGA                 ','SRBCK                  ', ! 3:4
     &             'HFEXCH                 ','SRGGA2                 ', ! 5:6
     &             'SRXWIB                 ','SRXWI2                 ', ! 7:8
     &             'SRXWI3                 ','SRXWI4                 ', ! 9:10
     &             'SRXPBEHSE              ','SRXPBETCS              ', ! 11:12
     &             'SRXPBERI               ','SRXPBEGWS              ', ! 13:14
     &             'LXLDAS                 ','LXLDA                  ', ! 15:16
     &             'LXPBEGWS               ','LXBCK                  ', ! 17:18
     &             'SRXLDA_C               ','SRXTPSS                ', ! 19:20
     &             'LRCWPBE                '/                           ! 21:22
      DATA CTABLE /'NULL                   ','SRCLDA                 ', ! 1:2
     &             'SRCGGA                 ','SRCMULO_A              ', ! 3:4
     &             'SRLYPT                 ','SRCWIB                 ', ! 5:6
     &             'SRCWI2                 ','SRCWI3                 ', ! 7:8
     &             'SRCWI4                 ','SRCPBETCS              ', ! 9:10
     &             'SRCPBETCSJ             ','SRCPBERI               ', ! 11:12
     &             'SRCPBELO               ','SRCPBEWI               ', ! 13:14
     &             'SRCPBEGWSC             ','SRCLDA_S               ', ! 15:16
     &             'SRCMULO_B              ','SRCMULO_C              ', ! 17:18
     &             'SRCMULO_D              ','SRCMULO_E              ', ! 19:20
     &             'SRCMDLDA               ','LNSCLDAS               ', ! 21:22
     &             'LNSCLDA                ','LSCLDA                 ', ! 23:24
     &             'LNSCPBEGWS             ','LSCPBEGWS              ', ! 25:26
     &             'LNCLYP                 ','LSCLYP                 ', ! 27:28
     &             'SRCLYPRI               ','SRCPBEGWS              ', ! 29:30
     &             'SRCLDA_C               ','SRCTPSS                ', ! 31:32
     &             'CPBE                   ','SRCVWN5                ', ! 33:34
     &             'SRCPW92                ','CVWN5                  '/ ! 35:36

C K. Sharkas and J. Toulouse beg
C List of lambda-dependent complement density functionals for linear separation of the electron-electron interaction.
C K. Sharkas, A. Savin, H. J. Aa. Jensen, J. Toulouse, J. Chem. Phys. 137, 044104 (2012)
C
C Complement exchange functionals:
C LXLDA: lambda-dependent exchange LDA
C LXLDAS: lambda-dependent exchange LDA with spin density
C LXPBEGWS: lambda-dependent exchange PBE (from srPBEGWS at mu=0)
C LXBCK: lambda-dependent exchange B88
C
C Complement correlation functionals:
C LSCLDA: lambda-dependent correlation LDA with density scaling
C LNSCLDA: lambda-dependent correlation LDA without density scaling
C LNSCLDAS: lambda-dependent correlation LDA without density scaling with spin density
C LSCPBEGWS: lambda-dependent correlation PBE (from srPBEGWS at mu=0) with density scaling
C LNSCPBEGWS: lambda-dependent correlation PBE (from srPBEGWS at mu=0) without density scaling
C LSCLYP: lambda-dependent correlation LYP with density scaling
C LNCLYP: lambda-dependent correlation LYP without density scaling
C
C See in herdrv.F for how to activate the scaling of the two-electron integrals by lambda.
C
C C K. Sharkas and J. Toulouse end

      DATA SETXFUN,SETCFUN /2*.FALSE./
      CALL QENTER('SRFUN_INPUT')
      IXWINT = 0
      ICWINT = 0
      mu_0   = -1.0D0

C Cut string SRFTMP
      call CUT_INTO_FIELDS(SRFTMP, ' ', maxfld, fields, lenfld, nfld)
      IF (nfld .EQ. 1) THEN
         SRXCFUN = fields(1)
         SRXFUN = ' '
         SRCFUN = ' '
         SRLOCALSPIN = ' '
         DFT_LOCALSPIN = .FALSE.
         DFT_SPINDNS = .TRUE.
      ELSE IF (nfld .GE. 2) THEN
         SRXCFUN = ' '
         SRXFUN = fields(1)
         SRCFUN = fields(2)
         SRLOCALSPIN = ' '
         DFT_LOCALSPIN = .FALSE.
         DFT_SPINDNS = .TRUE.
      END IF
      IF (nfld .GE. 3) THEN
         SRLOCALSPIN = fields(3)
         IF (INDEX(SRLOCALSPIN,'SRLOCALSP') .GT. 0) THEN
            DFT_LOCALSPIN = .TRUE.
            DFT_SPINDNS   = .TRUE.
            write (lupri,*) 'DFT_LOCALSPIN activated! ',SRLOCALSPIN
            ! Do not quit if field 3 not used here.
            ! field 3 might be mu_0 for SRCPBERI below 
            ! 160523-hjaaj
         ELSE IF (SRLOCALSPIN .EQ. 'NO_SPINDENSITY') THEN
            DFT_SPINDNS = .FALSE.
         END IF
      END IF
      
      DO I = 1, NXCTABLE
        IF (XCTABLE(I) .EQ. SRXCFUN) THEN
          GO TO (301, 302, 303, 304, 305, 306),I
        END IF
      END DO
      GO TO 10
C Short-range LDA functional (erk 30 jun. 20)
301   SRXFUN = 'SRXLDA                 '
      SRCFUN = 'SRCVWN5                '
      GO TO 10
C Short-range Goll-Werner-Stoll PBE functional (erk 30 jun. 20)
302   SRXFUN = 'SRXPBEGWS              '
      SRCFUN = 'SRCPBEGWS              '
      GO TO 10
C Short-range PBERI functional (erk 30 jun. 20)
303   SRXFUN = 'SRXPBEHSE              '
      SRCFUN = 'SRCPBERI               '
      GO TO 10
C Short-range Goll-Werner-Stoll PBE0 functional (erk 30 jun. 20)
304   SRXFUN = 'SRXPBEGWS              '
      SRCFUN = 'SRCPBEGWS              '
      HFXFAC = 0.25D0
      HFXSET = .TRUE.
      GO TO 10
C Long-range correction LDA functional (erk 30 jun. 20)
305   SRXFUN = 'SRXLDA                 '
      SRCFUN = 'CVWN5                  '
      GO TO 10
C Long-range correction PBE functional with GWS exchange (erk 30 jun. 20)
306   SRXFUN = 'SRXPBEGWS              '
      SRCFUN = 'CPBE                   '
      GO TO 10

  10  CONTINUE
      SETXFUN = .TRUE.
      DO I = 1, NXTABLE
        IF (XTABLE(I) .EQ. SRXFUN) THEN
          GO TO (101,102,103,104,105,106,107,108,109,110,111,112,113,
     &           114,115,116,117,118,119,120,121),I
        END IF
      END DO
      SETXFUN = .FALSE.
      GO TO 15
C     No Exchange functional
101   ISJT      = .TRUE.
      GO TO 15 ! NULL
C     No description
102   DOSRX_LDA = .TRUE.
      ISJT      = .TRUE.
      GO TO 15
C     No description
103   DOSRX_GGA = .TRUE.
      GO TO 15
C     No description
104   DOSRBCK  = .TRUE.
      GO TO 15
C     100 % Hartree-Fock like exchange
105   DOHFEXCH = .TRUE.
      IF (.NOT.HFXSET) THEN
         HFXFAC = 1.0D0
         HFXSET = .TRUE.
      END IF
      GO TO 15
C     No description
106   DOSRGGA2 = .TRUE.
      GO TO 15
C     No description
107   DOSRX_WIB = .TRUE.
      IXWINT    = 1
      GO TO 15
C     No description
108   DOSRX_WIB = .TRUE.
      IXWINT    = 2
      GO TO 15
C     No description
109   DOSRX_WIB = .TRUE.
      IXWINT    = 3
      GO TO 15
C     No description
110   DOSRX_WIB = .TRUE.
      IXWINT    = 12
      GO TO 15
C     Short-range exchange PBE functional of HSE
111   DOSRX_PBEHSE = .TRUE.
      ISJT        = .TRUE.
      GO TO 15
C     Short-range exchange PBE functional of TCS
112   DOSRX_PBETCS = .TRUE.
      ISJT        = .TRUE.
      GO TO 15  ! Manu 01-02-2006
C Short-range exchange rational interpolation functional
113   DOSRX_PBERI = .TRUE.
      ISJT    = .TRUE.
      GO TO 15
C Short-range exchange PBE functional of GWS (JT 11-08-09)
114   DOSRX_PBEGWS = .TRUE.
      ISJT    = .TRUE.
      GO TO 15
C Complement spin-dependent exchange LDA functional (KS 08-08-11)
C  based on DOSRX_LDA_S of JT for mu=0
115   DOLAX_LDAS = .TRUE.
      ISJT    = .TRUE.
      GO TO 15
C  based on standard dalton DOSRX_LDA for mu=0 (KS 2011) 
116   DOLAX_LDA = .TRUE.
      GO TO 15
C  Complement exchange PBE functional (KS 08-08-11)
117   DOLAX_PBEGWS = .TRUE.
      ISJT    = .TRUE.
      GO TO 15
C Complement Beck exchange functional (KS 08-08-11)
118   DOLAX_GGABCK = .TRUE.
      GO TO 15
C Short-range spin-dependent exchange LDA functional to fit with LDA (edh 10-02-16)
119   DOSRX_LDA = .TRUE.
      ISJT    = .TRUE.
      GO TO 15
C Short-range TPSS exchange functional (erk 27 aug. 18)
120   DOSRX_TPSS   = .TRUE.
      ISJT         = .TRUE.
      GO TO 15
C Short-range wPBE exchange functional (erk 27 feb. 19)
121   DOSRX_wPBE   = .TRUE.
      ISJT         = .TRUE.
      GO TO 15

  15  CONTINUE

      SETCFUN = .TRUE.
      DO I = 1, NCTABLE
        IF (CTABLE(I) .EQ. SRCFUN) THEN
          GO TO (201,202,203,204,205,206,207,208,209,210,211,212,213,
     &           214,215,216,217,218,219,220,221,222,223,224,225,226,
     &           227,228,229,230,231,232,233,234,235,236),I
        END IF
      END DO
      SETCFUN = .FALSE.
      GO TO 20

C     No correlation functional
201   ISJT      = .TRUE. 
      GO TO 20 ! NULL
C     No description
202   DOSRC_LDA = .TRUE.
      ISJT      = .TRUE. 
      GO TO 20
C     No description
203   DOSRC_GGA = .TRUE.
      GO TO 20
C     SRCMULO_A : AMUL = A/r_s
204   DOSRC_MULOCAL(0) = .TRUE.
      DOSRC_MULOCAL(1) = .TRUE.
      GO TO 20
C     No description
205   DOSRLYPT = .TRUE.
      GO TO 20
C     No description
206   DOSRC_WIB = .TRUE.
      ICWINT    = 1
      GO TO 20
C     No description
207   DOSRC_WIB = .TRUE.
      ICWINT    = 2
      GO TO 20
C     No description
208   DOSRC_WIB = .TRUE.
      ICWINT    = 3
      GO TO 20
C     No description
209   DOSRC_WIB = .TRUE.
      ICWINT    = 12
      GO TO 20
C     Short-range correlation PBE TCS functional
210   DOSRC_PBETCS = .TRUE.
      ISJT     = .TRUE.
      GO TO 20
C     Short-range correlation PBE TCS J functional
211   DOSRC_PBETCSJ = .TRUE.
      ISJT     = .TRUE.
      GO TO 20
C     Rational interpolation
212   DOSRC_PBERI = .TRUE.
      ISJT    = .TRUE.
      IF (nfld .ge. 3 .AND. .NOT.DFT_LOCALSPIN .AND. 
     &   (.NOT. SRLOCALSPIN .EQ. 'NO_SPINDENSITY')) THEN
         read (fields(3),*) MU_0
         mu_0_m4 = mu_0 ** (-4)
         mu_0_m6 = mu_0 ** (-6)
         mu_0_m7 = mu_0 ** (-7)
         write (lupri,*) 'SRC_PBE_RI mu_0 = ',mu_0,
     &      mu_0_m4,mu_0_m6,mu_0_m7
      END IF
      GO TO 20
C     No description
213   DOSRC_PBELO   = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C     Weighted interpolation
214   DOSRC_PBEWI = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Short-range correlation PBE function of GWS (JT 11-08-09)
215   DOSRC_PBEGWS = .TRUE.
      ISJT         = .TRUE. 
      GO TO 20
C Short-range spin-dependent correlation LDA function (JT 18-08-09)
216   DOSRC_LDA = .TRUE.
      ISJT    = .TRUE. 
      GO TO 20
C     SRCMULO_B :
C Short-range LDA type correlation functional with lokal mu value
C of mu=2*sqrt(alpha/pi)/sqrt(r_s)
217   DOSRC_MULOCAL(0) = .TRUE.   
      DOSRC_MULOCAL(2) = .TRUE.
      GO TO 20
C     SRCMULO_C :
C Short-range LDA type correlation functional with lokal mu value
C of mu=abs(grad n(r))/4*n(r)
218   DOSRC_MULOC_GGA = .TRUE.
      ISJT = .TRUE.
      GO TO 20
C     SRCMULO_D :
219   DOSRC_MULOD_GGA = .TRUE.
      ISJT = .TRUE.
      GO TO 20
C     SRCMULO_E :
C Extension of the SRCMULO_A functional: XMULFAC depends on grd/rho
220   DOSRC_MULOE_GGA = .TRUE.
      ISJT = .TRUE.
      GO TO 20
C Short-range correlation LDA function with multideterminant reference (for OEP) (JT 26-08-11)
C ... or for RSDHf calculations (Manu 06-12-2012) 
221   DOSRC_MD_LDA = .TRUE.
!     IF(DOSRX_OEP) THEN
!      DOSRXC_NOLRSYS=.TRUE.
!     END IF
      ISJT    = .TRUE.
      GO TO 20
C Complement spin-dependent (NON-Scaled) correlation LDA functional (KS 08-08-11)
222   DOLANSC_LDAS = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Complement (NON-Scaled) correlation LDA functional (KS 08-08-11)
223   DOLANSC_LDA = .TRUE.
      GO TO 20
C Complement (Scaled) correlation LDA functional (KS 08-08-11)
224   DOLASC_LDA = .TRUE.
      GO TO 20
C Complemnt Non-scaled correlation PBE function (KS 08-08-11)
225   DOLANSC_PBEGWS = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Complemnt (scaled) correlation PBE function (KS 08-08-11)
226   DOLASC_PBEGWS = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Non-scaled Complement correlation LYP functional (KS 08-08-11)
227   DOLANC_GGALYP = .TRUE.
      GO TO 20
C Scaled complement correlation LYP functional (KS 08-08-11)
228   DOLASC_GGALYP = .TRUE.
      GO TO 20
C Rational interpolation LYP functional                             
229   DOSRC_LYPRI = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Modified short-range spin-dependent correlation PBE functional of GWS (edh 21-03-15)
230   DOSRC_PBEGWS = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C Short-range spin-dependent correlation LDA functional w. PW92 param. (edh 29-03-15)
231   DOSRC_LDA_PW92 = .TRUE.
      ISJT           = .TRUE.
      GO TO 20
C Short-range TPSS spin-unpolarized correlation functional (erk 03 sep. 18)
232   DOSRC_TPSS   = .TRUE.
      ISJT         = .TRUE.
      GO TO 20
C mu independt PBE spin-polarized correlation functional (erk 27 feb. 19)
233   DOC_PBE_nomu = .TRUE.
      ISJT         = .TRUE.
      GO TO 20
C Short-range correlation LDA functional w. VWN5c param. (erk 27 feb. 19)
234   DOSRC_VWN5  = .TRUE.
      ISJT        = .TRUE.
      GO TO 20
C Short-range spin-dependent correlation LDA functional w. PW92 param. (edh 29-03-15)
235   DOSRC_LDA_PW92 = .TRUE.
      ISJT    = .TRUE.
      GO TO 20
C mu independt VWN5 spin-unpolarized correlation functional (erk 01 mar. 19)
236   DOC_VWN5_nomu = .TRUE.
      ISJT          = .TRUE.
      GO TO 20

  20  CONTINUE

      IF ((nfld .EQ. 1) .AND. (.NOT.SETXFUN) .AND. (.NOT.SETCFUN)) THEN
          WRITE(LUPRI,'(/2A/,A/)')
     & '  SRFUN_INPUT : Invalid or no SR exchange func. specified : ',
     &    SRXCFUN,'             Currently implemented functionals are :'
          DO J = 1,NXCTABLE
              WRITE(LUPRI,'(13X,A)') XCTABLE(J)
          ENDDO
      ENDIF

      IF (.NOT.SETXFUN) THEN
          WRITE(LUPRI,'(/2A/,A/)')
     & '  SRFUN_INPUT : Invalid or no SR exchange func. specified : ',
     &    SRXFUN,'             Currently implemented functionals are :'
          DO J = 1,NXTABLE
              WRITE(LUPRI,'(13X,A)') XTABLE(J)
          ENDDO
      ENDIF

      IF (.NOT.SETCFUN) THEN
          WRITE(LUPRI,'(/2A/A/)')
     & '  SRFUN_INPUT : Invalid or no SR correlation func. specified : '
     &   ,SRCFUN,'             Currently implemented functionals are :'
          DO J = 1,NCTABLE
              WRITE(LUPRI,'(13X,A)') CTABLE(J)
          ENDDO
      ENDIF
      IF (.NOT.(SETXFUN.AND.SETCFUN))
     &   CALL QUIT('Error in specfication of SR functionals')
      IF (IXWINT.EQ.ICWINT) THEN
         IWINT = IXWINT
      ELSE
         CALL QUIT ('Different extrapolation schemes specified for'//
     &              ' exchange and correlation')
      END IF

      IF (.NOT.HFXSET) THEN
C        if user has not specified HFXFAC with .HFXFAC or implicitly with HFEXCH
C        then set HF exchange factor to zero.
         HFXFAC = 0.0D0
         HFXSET = .TRUE.
      END IF
      CALL QEXIT('SRFUN_INPUT')
      RETURN
      END


C****************************************************************************
      SUBROUTINE SRDFT(ND_SIM,EXCMAT,DMAT,ESRDFT,DOERG,
     &                 DO_MOLGRAD,DOATR,TRIPLET,WORK,LWORK,IPRINT_in)
C*****************************************************************************
C    Purpose : Calculate SR DFT-energy and potential contribution for
C               SR-DFT hybrid methods.
C
C    Input:        ND_SIM : number of input density matrices and output potential matrices
C                  DOERG   true: xc energy and Vxc potential matrix
C                  DO_MOLGRAD   true: molecular gradient contribution
C                  DOATR   true: electronic Hessian, non-linear Exc second derivative terms
C                  TRIPLET true: triplet operators
C                  DMAT(1) = charge density matrix
C                  DMAT(2) = 1-index transformed charge density matrix (DOATR)
C                  DMAT(3) = spin density matrix (DFT_SPINDNS)
C                  DMAT(4) = 1-index transformed spin density matrix (DFT_SPINDNS and DOATR)
C            (if SRHYBR then different ordering of DMAT matrices)
C            (Note that the order of DMATs also fits for gradient calculations, DOATR false:
C             for energy and gradien calculations DMAT(2) will be the DVAO matrix needed
C             for Coulomb and exchange Fock matrices and not needed for SRDFT.)
C
C    Output:       EXCMAT(1) = charge density potential matrix
C                  EXCMAT(2) = spin density potential matrix (DFT_SPINDNS)
C                         or = special potential matrix (SRHYBR)
C                  ESRDFT(1) = srDFT energy
C                  ESRDFT(2) = srDFT exchange energy
C                  ESRDFT(3) = srDFT correlation energy
C
C    Created: Jesper Kielberg Pedersen, Mar. 2003
C    Modified: 17-08-09, J. Toulouse, add spin density matrix
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"

#ifdef VAR_MPI
!  mpif.h   : for MPI
#include "mpif.h"
      integer(kind=MPI_INTEGER_KIND) :: my_MPI_REAL8 = MPI_REAL8
      integer(kind=MPI_INTEGER_KIND) :: my_MPI_SUM   = MPI_SUM  
      integer(kind=MPI_INTEGER_KIND) :: ierr_mpi, len_mpi
#endif

! infpar.h : MYNUM, NODTOT
#include "maxorb.h"
#include "gnrinf.h"
#include "infpar.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "dftinf.h"
#include "dftcom.h"
      logical, intent(in) :: DO_MOLGRAD, DOERG, DOATR, TRIPLET
      logical :: DOLND, DOGGA, do_metagga
      logical :: first_call, LSRHYBR
      DATA       first_call /.true./
      real*8 :: EXCMAT(NBAST,NBAST,*),DMAT(N2BASX,*),ESRDFT(3)
      real*8 :: WORK(LWORK)
      real*8 :: RESULTS(11), RESULTS_tmp(11)
      real*8, allocatable :: excmat_tmp(:)
      integer*8 :: i,j,loc

      CALL QENTER('SRDFT')
      CALL TIMER('START ',TIMSTR,TIMEND)
      RESULTS(:) = 0.0D0

      IF (TRIPLET) THEN
#if SRDFT_DEBUG > 4
         WRITE(lupri,*) 'srdft.F: Calling srDFT with triplet'
         flush(lupri)
#endif
      END IF

      IPRINT  = max(IPRINT_in, SRDFT_DEBUG)
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      IF (DFT_SPINDNS .AND. LSRHYBR) THEN
         CALL QUIT('SRDFT: Spin density AND SRHYBR not implemented')
      END IF

#if SRDFT_DEBUG > 4
      write(lupri,*) 'Entered SRDFT, SDRDT_DEBUG defined'
      write(lupri,*) 'NDER,NTYPSO =',NDER,NTYPSO
#endif
#if SRDFT_DEBUG > 10
      WRITE(LUPRI,'(//A//)') ' --- SRDFT: Total EXCMAT on entry ---'
      CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      WRITE(LUPRI,'(//A//)') ' --- SRDFT: Total density matrix ---'
      CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      IF (LSRHYBR) THEN
        WRITE(LUPRI,'(//A//)') ' --- SRDFT: Valence density matrix ---'
        CALL OUTPUT(DMAT(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        IF (DOATR) THEN
          WRITE(LUPRI,'(//A//)') ' --- SRDFT: Total DTRMAT ---'
          CALL OUTPUT(DMAT(1,3),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
          WRITE(LUPRI,'(//A//)') ' --- SRDFT: Valence DTRMAT ---'
          CALL OUTPUT(DMAT(1,4),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        ENDIF
      ELSE IF (DOATR) THEN
        WRITE(LUPRI,'(//A//)') ' --- SRDFT: Total DTRMAT ---'
        CALL OUTPUT(DMAT(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      IF (DFT_SPINDNS) THEN
        WRITE(LUPRI,'(//A/)') ' --- SRDFT: Spin density matrix ---'
        CALL OUTPUT(DMAT(1,3),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        IF (DOATR) THEN
          WRITE(LUPRI,'(//A/)') ' --- SRDFT: Spin DTRMAT ---'
          CALL OUTPUT(DMAT(1,4),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
        WRITE(LUPRI,'(//A//)')' --- SRDFT: Total spin EXCMAT on entry'
        CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
#endif
!     Shall we go in MPI parallel mode for srdft ?

#ifdef VAR_MPI
      IF (MYNUM .EQ. 0 .AND. NODTOT .GT. 0) THEN
#if SRDFT_DEBUG > 4
         WRITE(lupri,*) 'calling srdft_par_master, iprint',iprint
         flush(lupri)
#endif
         CALL SRDFT_PAR_MASTER(
     &        ND_SIM,DMAT,ESRDFT,
     &        DOERG,DO_MOLGRAD,DOATR,TRIPLET,
     &        WORK,LWORK,IPRINT)
#if SRDFT_DEBUG > 4
         WRITE(lupri,*) 'returned from srdft_par_master'
         flush(lupri)
#endif
         GO TO 8000
      END IF
#endif

!     Enter here for sequential runs (NODTOT.eq.0) or for slaves
C     =================================
C     Initializations :
C     Number of AOs and their addresses
C     =================================
      DFTHR0 = 1.0D-8
      DFTHRL = 1.0D-10
      DFTHRI = 1.0D-13
      DOLND  = .FALSE.

C     ... erk : For setting meta-GGAs
      do_metagga = DOSRC_TPSS_S .OR. DOSRX_TPSS_S .OR. DOSRC_TPSS .OR.
     &             DOSRX_TPSS
C     ... jkp : For gradient correction to SR-LDA
      DOGGA = DOSRX_GGA.OR.DOSRC_GGA.OR.DOSRBCK.OR.DOSRGGA2.OR.
     &        DOSRX_WIB.OR.DOSRC_WIB.OR.DOSRX_PBEHSE.OR.DOSRX_PBETCS.OR.
     &        DOSRC_PBETCS.OR.DOSRX_PBERI.OR.DOSRC_PBERI.OR.DOSRX_PBEGWS
     &        .OR.DOSRC_PBEGWS.OR.DOSRC_MULOC_GGA.OR.DOSRC_MULOD_GGA
     &        .OR.DOSRC_MULOE_GGA.OR.DOSRC_PBELO.OR.DOLAX_PBEGWS.OR.
     &        DOLANSC_PBEGWS.OR.DOLASC_PBEGWS.OR.DOLAX_GGABCK.OR.
     &        DOLANC_GGALYP.OR.DOLASC_GGALYP.OR.DOSRC_LYPRI .OR.
     &        DOSRX_PBEGWS_S.OR.DOSRC_PBEGWS_S.OR.DOSRC_PBEGWS_VWN5
     &       .OR.DOC_PBE_nomu_S.OR.DOC_PBE_nomu.OR.DOSRX_wPBE_S
     &       .OR.DOSRX_wPBE
C     ... erk : Sets DOGGA active when using meta-GGAs
     &        .OR. do_metagga
C

      NDER  = 0
      IF (DO_MOLGRAD) NDER = 1
      IF (LSRHYBR) THEN
         NDER = NDER + 2
C        ... we use DFTPOT for SRHYBR
      ELSE IF (DOGGA) THEN
         NDER = NDER + 1
         IF (.NOT.ISJT) NDER = NDER + 1
C        ... we always use DFTPOT when not ISJT
CC       IF (DFTPOT) NDER = NDER + 1
      END IF
C
C
      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      IF (NDER.GT.2) CALL QUIT('NDER.gt.2 not implemented')
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1= (NSOB1-1)*NBAST + 1
C     =============================
C     Allocations
C     =============================
      NBUF   = 1 000 000
      KX     = 1
      KY     = KX    + NBUF
      KZ     = KY    + NBUF
      KW     = KZ    + NBUF
      KGSO   = KW    + NBUF
      KCNT   = KGSO  + NTYPSO*NBAST
      KDGASR = KCNT  + NBAST
      IF (DFT_SPINDNS) THEN
         KDGAT = KDGASR + 2*NBAST
      ELSE
         KDGAT = KDGASR + NBAST
      END IF
      IF (LSRHYBR) THEN
         KWRK  = KDGAT  + NBAST
      ELSE
         KWRK  = KDGAT
      END IF
      LWRK  = LWORK + 1 - KWRK
      IF (KWRK .GT. LWORK) CALL ERRWRK('SRDFT',-KWRK,LWORK)

      IF(ISJT) THEN
        CALL SRDFT1JT(ND_SIM,EXCMAT,DMAT,DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGASR),DOGGA,DOERG,DO_MOLGRAD,DOATR,TRIPLET,
     &            do_metagga,
     &            RESULTS,WORK(KWRK),LWRK,IPRINT)
      ELSE
        CALL SRDFT1(ND_SIM,EXCMAT,DMAT,DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGASR),WORK(KDGAT),DOGGA,DOERG,DO_MOLGRAD,DOATR,
     &            TRIPLET,RESULTS,WORK(KWRK),LWRK,IPRINT)
      ENDIF

C
 8000 CONTINUE

      IF (DFT_SPINDNS .OR. LSRHYBR) THEN
        NEXCMAT = 2
      ELSE
        NEXCMAT = 1
      ENDIF
!
!     parallel stuff
!
#if defined  (VAR_MPI)
      IF (NODTOT .GT. 0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
C
         IF (MYNUM .EQ. 0) THEN
            len_mpi = 11
            RESULTS_TMP(:) = 0.0d0 ! RESULTS(:) is already zeroed
            CALL MPI_REDUCE(RESULTS_TMP,RESULTS,len_mpi,
     &        my_MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr_mpi)

            length  = N2BASX*NEXCMAT
            allocate(excmat_tmp(length))
            excmat_tmp(1:length) = 0.0d0
            call dzero(EXCMAT,length)
            len_mpi = length
            CALL MPI_REDUCE(excmat_tmp,EXCMAT,len_mpi,
     &         my_MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr_mpi)
            deallocate(excmat_tmp)
         ELSE
            len_mpi = 11
            CALL MPI_REDUCE(RESULTS,dummy,len_mpi,
     &         my_MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr_mpi)

            len_mpi = N2BASX*NEXCMAT
            CALL MPI_REDUCE(EXCMAT,dummy,len_mpi,
     &         my_MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr_mpi)
         END IF
         IF (IPRINT .GE. 1) THEN
            write(lupri,*) 'srdft results on node, DOATR',mynum,DOATR
            write(lupri,'(6F15.10)') RESULTS(1:11)
            write(lupri,*) 'srdft excmat on node, DOATR',mynum,DOATR
            call output(excmat,1,nbast,1,nbast,nbast,nbast,-1,lupri)
            flush(lupri)
         END IF
      END IF
#endif

!     end of parallel stuff, now only master continues

      IF (MYNUM .GT. 0) GO TO 9000

      ELCTRN    = RESULTS(1)
      ELCTRSKIP = RESULTS(2)
      SPINDEN   = RESULTS(3)
      ExInt     = RESULTS(4)
      EcInt     = RESULTS(5)
      NRSGTR    = NINT(RESULTS(6))
      ! LSRHYBR
      ELCTSR    = RESULTS(7)
      EXSR      = RESULTS(8)
      ECSR      = RESULTS(9)
      EXTOT     = RESULTS(10)
      ECTOT     = RESULTS(11)
C
      IF (DOERG) THEN
         ESRDFT(1) = EXint + ECint
         ESRDFT(2) = Exint
         ESRDFT(3) = Ecint
         IF (LSRHYBR) THEN
            ESRDFT(1) = ESRDFT(1) + EXTOT - EXSR + ECTOT - ECSR
         END IF
         FAC = 0.5D0
      ELSE IF (DOATR) THEN
         FAC = 1.0D0
      ELSE
         CALL QUIT('both DOERG and DOATR are false')
      END IF

      DO K = 1, NEXCMAT
      DO I = 1, NBAST
      DO J = 1, I - 1
         AVERAG = FAC*(EXCMAT(I,J,K) + EXCMAT(J,I,K))
         EXCMAT(I,J,K) = AVERAG
         EXCMAT(J,I,K) = AVERAG
      END DO
      END DO
      END DO
C
C     Print section
C
      IF (IPRINT .GE. 1) THEN
      IF (DOERG) THEN
         ELCTOT = ELCTRN + ELCTRSKIP
         IF (LSRHYBR) THEN
            WRITE (LUPRI,'(3(/A,F18.10))')
     &' No. of active electrons from num. int.         :',ELCTRN,
     &' No. of active electrons skipped in num. int.   :',ELCTRSKIP,
     &' Sum of active and skipped electrons            :',ELCTOT
            WRITE (LUPRI,'(A,F18.10)')
     &' No. of valence electrons  (SRHYBR)             :',ELCTSR
            WRITE (LUPRI,'(/A,F18.10)')
     &' Exchange energy (SRHYBR: SR valence part)      :',EXint
            WRITE (LUPRI,'(A,F18.10)')
     &' Correlation energy  (SRHYBR: SR valence part)  :',ECint,
     &' Exchange energy (SRHYBR: BLYP valence part)    :',EXSR,
     &' Correlation energy (SRHYBR: BLYP valence part) :',ECSR,
     &' Exchange energy (SRHYBR: BLYP total)           :',EXTOT,
     &' Correlation energy (SRHYBR: BLYP total)        :',ECTOT
            EXint = EXint + EXTOT - EXSR
            ECint = ECint + ECTOT - ECSR
         ELSE
            WRITE (LUPRI,'(4(/A,F18.10))')
     & ' No. of non-skipped electrons from num. int.    :',ELCTRN,
     & ' No. of     skipped electrons   in num. int.    :',ELCTRSKIP,
     & ' Sum of non-skipped and skipped electrons       :',ELCTOT,
     & ' Integrated spin density (N_alpha-N_beta)       :',SPINDEN
         END IF

         WRITE (LUPRI,'(3(/A,F18.10),/A,I18)')
     & ' Exchange energy                                :',ExInt,
     & ' Correlation energy                             :',EcInt,
     & ' DFT exchange-correlation energy                :',ESRDFT(1),
     & ' Number of times RHO(2) > RHO(1)                :',NRSGTR

      END IF ! DOERG
      IF (DOATR) THEN
         WRITE(LUPRI,'(/A)') ' srDFT linear transformation finished'
      END IF

      IF (IPRINT.GE.8) THEN
         WRITE(LUPRI,'(//A/)') ' --- SRDFT: Total EXCMAT on exit'
         CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         IF (DFT_SPINDNS) THEN
            WRITE(LUPRI,'(//A/)')' --- SRDFT: Total spin EXCMAT on exit'
            CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
         END IF
         IF (LSRHYBR) THEN
            WRITE(LUPRI,'(//A/)')
     &         ' --- SRDFT: SRHYBR valence EXCMAT on exit'
            CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
         END IF
      END IF  ! IPRINT .GE. 5

      END IF  ! IPRINT .GE. 1

 9000 IF (IPRINT .GT. 0) CALL TIMER('SRDFT ',TIMSTR,TIMEND)
      CALL QEXIT('SRDFT')
      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRDFT1(ND_SIM,EXCMAT,DMAT,DOLND,CORX,CORY,CORZ,
     &                  WEIGHT,NBUF,GSO,NCNT,DMATGAOSR,DMATGAOT,DOGGA,
     &                  DOERG,DO_MOLGRAD,DOATR,TRIPLET,RESULTS,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate DFT-energy and potential contribution for
C               SR-DFT hybrid methods.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
! infpar.h : MYNUM, NODTOT
#include "maxorb.h"
#include "gnrinf.h"
#include "infpar.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
#include "chrnos.h"

      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
      PARAMETER (DP3 = 1.0D0/3.0D0)

      LOGICAL   DOLND,DOGGA,DOERG,DO_MOLGRAD,DOATR,TRIPLET
      LOGICAL   LDOGGA,LSRHYBR,RTEST,ldebug
      REAL*8    EXCMAT(NBAST,NBAST,*), DMAT(N2BASX,*),
     &          CORX(NBUF),CORY(NBUF),CORZ(NBUF),WEIGHT(NBUF),
     &          GSO(NBAST*NTYPSO),NCNT(NBAST),
     &          DMATGAOSR(NBAST),DMATGAOT(NBAST),RESULTS(11),
     &          WORK(LWORK),
     &          RHGT(3),RHGSR(3),ZKSR(3,3),XKSR(3,3), COR(3)

      CHARACTER*17 QUADNAME

      INTEGER, SAVE :: NTOT_DFTGRID_save = -1

#if SRDFT_DEBUG > 0
      write(lupri,*) 'Entered SRDFT1, SRDFT_DEBUG defined =',SRDFT_DEBUG
      WRITE(LUPRI,*) 'SRDFT1: ND_SIM =',ND_SIM
      WRITE(LUPRI,*) 'SRDFT1: CHIVAL =',CHIVAL
      WRITE(LUPRI,*) 'SRDFT1: HFXFAC =',HFXFAC
      WRITE(LUPRI,*) 'SRDFT1: DOERG,DOATR,DO_MOLGRAD,DOGGA,TRIPLET  =',
     &   DOERG,DOATR,DO_MOLGRAD,DOGGA,TRIPLET
      WRITE(LUPRI,*) 'SRDFT1: ERFEXP(0:2) =',ERFEXP(0:2)
#endif

      IF (ND_SIM .NE. 1) THEN
         WRITE(LUPRI,*) 'SRDFT1: illegal ND_SIM ',ND_SIM
         CALL QUIT('SRDFT1: illegal ND_SIM value')
      END IF

      IF (DO_MOLGRAD) THEN
         CALL QUIT('SRDFT1: DO_MOLGRAD not implemented yet')
      END IF

      IF (DFT_SPINDNS) THEN
         CALL QUIT('SRDFT1: spin density not implemented yet')
      END IF
      IF (TRIPLET) THEN
         CALL QUIT('SRDFT1: triplet not implemented yet')
      ENDIF


      ldebug = .false.
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      IF (LSRHYBR) THEN
         NEXCMAT = 2
      ELSE
         NEXCMAT = 1
      END IF
#if SRDFT_DEBUG > 10
!     debug print of quadrature points between IPNTstart and IPNTend
      !IPNTstart = 8168
      IPNTstart = 0
C     IPNTend   = 10210
      IPNTend   = SRDFT_DEBUG
#endif
C     =============================
C     Calculate DFT grid
C     (NTOT_DFTGRID is the total number of grid points)
C     =============================

      IF (.NOT.DFTGRID_DONE_OLD .OR. NTOT_DFTGRID_save .LE. 0) THEN
         CALL MAKE_DFTGRID(WORK,LWORK,NTOT_DFTGRID,1,.FALSE.)
         CALL CONDFT
         NTOT_DFTGRID_save = NTOT_DFTGRID
         DFTGRID_DONE_OLD = .TRUE.
      ELSE
         NTOT_DFTGRID = NTOT_DFTGRID_save
         ! because wrong NTOT_DFTGRID has been transferred from master in SRDFT_PAR_NODE
      END IF

      LDOGGA = DOGGA .OR. LSRHYBR

      ELCTRN = D0
      ELCTSR = D0
      ELCTRSKIP = D0
C     ... Initialize energies to zero. They should come out
C         zero if NULL functional is defined.
      EX     = D0
      EC     = D0
      EXC    = D0
      ETX    = D0
      ETC    = D0
      VXC    = D0
      VTX    = D0
      VTC    = D0
      VXBT   = D0
      VXBV   = D0
      VTOT   = D0

C     For SRHYBR:

      EXSR   = D0
      EXTOT  = D0
      ECSR   = D0
      ECTOT  = D0

C     HJAaJ: THRELCTRN makes sure that error in
C            no. of electrons from grid integration less than DFTHR0
C            (NTOT_DFTGRID is the total number of grid points)

C    THRELCTRN = 1.0D-6 / NTOT_DFTGRID
      THRELCTOT = 1.0D-7
      THRELCTRN = THRELCTOT/100
C     THRDRHO   = DFTHRI
      THRDRHO   = 1.0D-10
C Debug values: include everything
#ifdef DFT_NOSCREEN
      THRELCTRN = D0
      THRDRHO   = D0
      DFTHR0    = D0
      DFTHRI    = D0
      DFTHRL    = D0
#endif

C     HJAaJ: skip eval. of neighbour points with same RHO

      RHOOLD = D0
      NSKIP1 = 0
      NSKIP2 = 0

C     Make quadname: Can take 9999 procs

      IF (MYNUM .EQ. 0) THEN
         QUADNAME = 'DALTON.QUAD'
      ELSE
         QUADNAME = 'DALTON.QUAD.n'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
      END IF

      LUQUAD = -1
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)

      IF (MYNUM .EQ. 0) THEN
         my_first_point = 1
         my_last_point  = NTOT_DFTGRID
      ELSE
         NPNTS_PER_NODE = (NTOT_DFTGRID-1)/NODTOT + 1
         my_first_point = NPNTS_PER_NODE*(MYNUM-1) + 1
         my_last_point  = MIN(NPNTS_PER_NODE*MYNUM, NTOT_DFTGRID)
         if (iprint .gt. 5) write(lupri,*)
     &      'MYNUM,NPNTS_PER_NODE,my_first_point,my_last_point',
     &       MYNUM,NPNTS_PER_NODE,my_first_point,my_last_point
      END IF
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
CDBG  kpoint = npoint/10

      IF (NPOINT.GT.0) THEN
         NPNT1 = NPNTS + 1       ! first integration point from this record
         NPNTS = NPNTS + NPOINT  ! last  integration point from this record
         IPNT_ST  = max(NPNT1,my_first_point) + 1 - NPNT1
         IPNT_END = min(NPNTS,my_last_point) + 1 - NPNT1
         if (iprint .gt. 5) write(lupri,*)
     &      'MYNUM,IPNT_ST,IPNT_END', MYNUM,IPNT_ST,IPNT_END
         IF (IPNT_ST .GT. IPNT_END) THEN
            READ (LUQUAD) ! skip this record
            GO TO 200
         END IF
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
C   Loop over grid points -----------------------------------------------
         DO 300 IPNT = IPNT_ST, IPNT_END
            IF (IPRINT .GT. 100) THEN
C              Print grid information
               WRITE (LUPRI,'(I8,4F12.6)')
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF

            WGHT = WEIGHT(IPNT)

C           AOs
C           ===

            THRINT = DFTHRI/WGHT
            COR(1) = CORX(IPNT)
            COR(2) = CORY(IPNT)
            COR(3) = CORZ(IPNT)
            CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
     &                  NBAST,DOLND,LDOGGA,THRINT,IPRINT)
C           =============================
C           Density
C           =============================
            IF (LSRHYBR) THEN
               CALL GETRHO_OLD(DMAT(1,2),GSO(KSO0),RHOSR,
     &                     DMATGAOSR,THRINT,IPRINT)
               RHOSR13 = RHOSR ** DP3
               CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHOTOT,
     &                     DMATGAOT,THRINT,IPRINT)
               RHOT13 = RHOTOT ** DP3
            ELSE
               CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHOSR,
     &                     DMATGAOSR,THRINT,IPRINT)
               RHOSR13 = RHOSR ** DP3
               RHOTOT = RHOSR
               RHOT13 = RHOSR13
               CALL DCOPY(NBAST,DMATGAOSR,1,DMATGAOT,1)
            END IF
            DELCTRN = WGHT*RHOTOT
            DELCTSR = WGHT*RHOSR
C            IF (RHOTOT.LE.DFTHR0) THEN
            IF (ABS(DELCTRN).LE.THRELCTRN) THEN
               NSKIP1 = NSKIP1 + 1
               ELCTRSKIP = ELCTRSKIP + DELCTRN
               THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
            ELSE
C
               RHOTST = RHOTOT
               IF (LDOGGA) THEN
C             ... For the gradient correction to SR-LDA we need
C                 the gradient of the density for the energy and
C                 the hessian of the density for the potential.

C                 Gradient of density
C                 ===================
                  CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMATGAOSR,
     &                         1,D0,RHGSR,1)
                  RHOSRGRD = SQRT(RHGSR(1)**2 + RHGSR(2)**2
     &                     + RHGSR(3)**2)
                  IF (LSRHYBR) THEN
                     CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMATGAOT,
     &                          1,D0,RHGT,1)
                     RHOTGRD = SQRT(RHGT(1)**2 + RHGT(2)**2
     &                       + RHGT(3)**2)
                  ELSE
                     RHGT(1) = RHGSR(1)
                     RHGT(2) = RHGSR(2)
                     RHGT(3) = RHGSR(3)
                     RHOTGRD = RHOSRGRD
                  ENDIF

C                 Hessian of density
C                 ===================

                  IF (LSRHYBR) THEN
                     CALL DFTRHH(DMAT(1,2),DMATGAOSR,
     &                           GSO(KSO0),GSO(KSO1),GSO(KSO2),
     &                           RHGSR,RHOSRLAP,RHOSRGHG)
                     CALL DFTRHH(DMAT,DMATGAOT,GSO(KSO0),GSO(KSO1),
     &                           GSO(KSO2),RHGT,RHOTLAP,RHOTGHG)
                  ELSE
                     CALL DFTRHH(DMAT,DMATGAOSR,GSO(KSO0),GSO(KSO1),
     &                           GSO(KSO2),RHGSR,RHOSRLAP,RHOSRGHG)
                     RHOTLAP = RHOSRLAP
                     RHOTGHG = RHOSRGHG
                  END IF
                  RHOTST = RHOTST + RHOTGRD + RHOTLAP + RHOTGHG
               ENDIF
C              =============================
C              Number of electrons
C              =============================
               ELCTRN = ELCTRN + DELCTRN
               ELCTSR = ELCTSR + DELCTSR
C              =============================
C              ETX,ETC: Energy
C              VTX,VTC: Exchange-correlation potential
C              ==============================
               IF (DOERG) THEN
C              RTEST = (ABS(RHOTST-RHOOLD) .GT. THRDRHO)
C              IF (RTEST) THEN
               IF (.true.) THEN
C              ... hjaaj 10. July 2003: we could make this better by using that
C                  EXC(RHO) = EXC(RHOOLD+DRHO) ~ EXC(RHOOLD) + VXC(RHOOLD)*DRHO
C                  = ETX+ETC + (VTX+VTC)*DRHO (as EXC and VXC are calculated for RHOOLD)
                  RHOOLD = RHOTST
C                 ... only reevaluate if new RHOTOT /hjaaj July 2003

C                 -------------------------------------------------------
C                 Exchange Correlation energy and potential
C                 -------------------------------------------------------

#if SRDFT_DEBUG > 10
                  if (IPNT.GE.IPNTstart.AND.IPNT.LE.IPNTend) then
                     ldebug=.true.
                  else
                     ldebug=.false.
                  end if
#endif
                  CALL SRDFTEXC(XKSR,ZKSR,RHOSR,RHOSR13,RHOSRGRD,
     &                          RHOSRLAP,RHOSRGHG,CHIVAL,ERFEXP,
     &                         .TRUE.,VXB,ldebug,VLAMBDA)
               ELSE
C                 we can use ETX, ETC, VTX, VTC from last IPNT
                  NSKIP2 = NSKIP2 + 1
               END IF
               EX   = EX + WGHT*XKSR(1,1)
               EC   = EC + WGHT*ZKSR(1,1)
               VXC  = WGHT*(XKSR(2,1) + ZKSR(2,1))

C             ... SRHYBR (RHOTOT) contribution
C              IF (RTEST) THEN
               IF (.true.) THEN
                  IF (LSRHYBR) THEN
                     VXB   = -WGHT*VXB
                     EXSR  = EXSR + WGHT*XKSR(1,2)
                     ECSR  = ECSR + WGHT*ZKSR(1,2)
                     VXC   = VXC - WGHT*(XKSR(2,2)+ZKSR(2,2))
                     CALL SRDFTEXC(XKSR,ZKSR,RHOTOT,RHOT13,RHOTGRD,
     &                             RHOTLAP,RHOTGHG,CHIVAL,ERFEXP,
     &                            .FALSE.,VXB,ldebug,VLAMBDA)
                     EXTOT = EXTOT + WGHT*XKSR(1,2)
                     ECTOT = ECTOT + WGHT*ZKSR(1,2)
                     VXCT  = WGHT*(XKSR(2,2)+ZKSR(2,2))
                     VXBT  = WGHT*VXB
                  ENDIF
               ENDIF
#if SRDFT_DEBUG > 10
C               if (mod(IPNT,kpoint) .eq. 0) then
               if (IPNT.GE.IPNTSTART.AND.IPNT.LE.IPNTEND) then
                  write (lupri,*) 'IPNT,ex,ec,vx,vc',IPNT,ex,ec,vx,vc
                  write (lupri,*)
     &            'wght,xksr,zksr',wght,xksr(1,1),zksr(1,1)
                  write (lupri,*)
     &            'rhotot,rhot13,rhotgrd,rhotlap,rhotghg : ',
     &            rhotot,rhot13,rhotgrd,rhotlap,rhotghg
                  write (lupri,*)
     &            'rhosr,rhosr13,rhosrgrd,rhosrlap,rhosrghg : ',
     &            rhosr,rhosr13,rhosrgrd,rhosrlap,rhosrghg
                  write (lupri,*) 'nskip1,nskip2',nskip1,nskip2
                  write (lupri,*) 'excmat(1,1,1),gso(1) before: ',
     &                             excmat(1,1,1),gso(1)
               end if
#endif
C              =====================================================
C              Exchange-correlation contribution to Kohn-Sham matrix
C              =====================================================
C -hjaaj: implement screening>             IF (ABS(VXC) .GT. SCRTHR)
C  --- NOTE ! This uses FROMVX.
!      WRITE(LUPRI,*) 'SRDFT1JT: IPNT=',IPNT,' VXC=',VXC !JT

               VXB = 0.25d0*VXB ! 23-Apr-2018 hjaaj: new definition of VXB
               CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHGSR,VXC,VXB,
     &                     LDOGGA,.true.,DFTHRL,IPNT)
               IF (LSRHYBR) THEN
C             ... SRHYBR (RHOTOT) contribution
                  VXBT = 0.25d0*VXBT ! 23-Apr-2018 hjaaj: new definition of VXBT
                  CALL DFTKSM(EXCMAT(1,1,2),GSO(KSO0),GSO(KSO1),RHGT,
     &                        VXCT,VXBT,LDOGGA,.true.,DFTHRL,IPNT)
               END IF

               ENDIF

C              Hessian transformation
C              ======================

               IF (DOATR) THEN
                  IF (LSRHYBR) THEN
C                    matrix 1 & 2 are DTAO and DVAO
C                    matrix 3 & 4 are DXTAO and DXVAO
                     JDXAO = 4
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.TRUE.,RHOSR,RHOSR13,
     &                             RHOSRGRD,RHOSRLAP, RHOSRGHG,RHGSR,
     &                             DMATGAOSR,VLAMBDA)
                     JDXAO = 3
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT(1,1,2),
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.FALSE.,RHOTOT,RHOT13,RHOTGRD,
     &                             RHOTLAP,RHOTGHG,RHGT,
     &                             DMATGAOT,VLAMBDA)
                  ELSE
C                    matrix 1 is DTAO
C                    matrix 2 is DXTAO
                     JDXAO = 2
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.TRUE.,RHOTOT,RHOT13,
     &                             RHOTGRD,RHOTLAP,RHOTGHG,RHGT,
     &                             DMATGAOT,VLAMBDA)
                  END IF
               END IF
            ENDIF

  300    CONTINUE
C End of loop over points ------------------------------------------------
C
         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF

 1000 CALL GPCLOSE(LUQUAD,'KEEP')

      RESULTS(1) = ELCTRN
      RESULTS(2) = ELCTRSKIP
      RESULTS(3) = 0.0D0 ! SPINDEN
      RESULTS(4) = EX
      RESULTS(5) = EC
      RESULTS(6) = 0 ! NRSGTR
      ! LSRHYBR
      RESULTS(7) = ELCTSR
      RESULTS(8) = EXSR
      RESULTS(9) = ECSR
      RESULTS(10)= EXTOT
      RESULTS(11)= ECTOT

      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRDFTEXC(XKSR,ZKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                    CHIVAL,ERFEXP,LSRLDA,VXB,ldebug,VLAMBDA)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Dec. 2003
C
C     Purpose : Driver for Short-range DFT exchange and correlation
C               contributions to energy and potential.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "dftinf.h"
#include "dftcom.h"
      real*8   VLAMBDA
      LOGICAL ERFEXP(0:2),LSRLDA,ldebug
      DIMENSION XKSR(3,2), ZKSR(3,2), ZMU(3)
C
      VB    = 0.0D0
      VL    = 0.0D0
      VXB   = 0.0D0
      CALL DZERO(XKSR,6)
      CALL DZERO(ZKSR,6)
      CALL DZERO(ZMU ,3)

C     -------------------------------------------------------
C     Exchange energy and potential
C     -------------------------------------------------------

      IF (SRHYBR.AND.LSRLDA) THEN
          CALL VXSRLDA(XKSR,RHO,RHO13,CHIVAL,1,ERFEXP)
      END IF

      IF (DOLAX_LDA) THEN
          CALL VXLALDA (XKSR,RHO,RHO13,VLAMBDA)
      END IF

      IF (DOSRX_GGA) THEN
          CALL VXSRGGA(XKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                 CHIVAL,ERFEXP)
      END IF 

      IF (DOLAX_GGABCK) THEN
          CALL VLAXBEK(XKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                 VLAMBDA)
      END IF

      IF (DOSRGGA2) THEN
          CALL VXSRGGA2(XKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                  RHOGHG,CHIVAL,ERFEXP)
      END IF
      IF (DOSRX_WIB) THEN
          CALL EBCK(XKSR(1,1),RHO,RHO13,RHOGRD)
          CALL VBCK(XKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
          IF (IWINT.EQ.1) THEN
             CALL VXSRLDA(ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
             CALL VSRWINT(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3) THEN
             CALL VXSRLDA(ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
             CALL EDRC(ED,RHO,RHO13)
             CALL VDRC(VD,RHO13)
             XKSR(1,1) = XKSR(1,1) + ED
             XKSR(2,1) = XKSR(2,1) + VD
             CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GT.10) THEN
             CALL EDRC(ED,RHO,RHO13)
             CALL VDRC(VD,RHO13)
             XKSR(1,1) = XKSR(1,1) + ED
             XKSR(2,1) = XKSR(2,1) + VD
             CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE
            CALL QUIT ('Wrong dft functional interpolation scheme'//
     &                 'in DFTEXC (IWINT).')
          END IF
      END IF
      IF (DOSRBCK ) THEN
          CALL EDRC(ED,RHO,RHO13)
          CALL EBCK(EB,RHO,RHO13,RHOGRD)
          XKSR(1,1) = ED + EB
          CALL VDRC(VD,RHO13)
          CALL VBCK(VB,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
          XKSR(2,1) = VD + VB
      ENDIF
      IF (SRHYBR) THEN
          CALL EDRC(ED,RHO,RHO13)
          CALL VDRC(VD,RHO13)
#ifdef LDAembed
          EB  = 0.0D0
          VB  = 0.0D0
          VXB = 0.0D0
#else
          CALL EBCK(EB,RHO,RHO13,RHOGRD)
          CALL GBCK(VB,VXB,RHO,RHO13,RHOGRD)
#endif
          XKSR(1,2) = ED + EB
          XKSR(2,2) = VD + VB
      END IF
C
C     -------------------------------------------------------
C     Correlation energy and potential
C     -------------------------------------------------------
C
      IF (DOSRC_MULOCAL(0).OR.(SRHYBR.AND.LSRLDA)) THEN
          CALL VCSRLDA(ZKSR,RHO,RHO13,CHIVAL,DOSRC_MULOCAL,1,ERFEXP)
      END IF

      IF (DOLANSC_LDA) THEN
          CALL VCLALDA (ZKSR,RHO,RHO13,VLAMBDA)
      END IF
      IF (DOLASC_LDA) THEN
          CALL VCLASLDA (ZKSR,RHO,RHO13,VLAMBDA)
      END IF

      IF (DOSRC_GGA) THEN
          CALL VCSRGGA(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                 CHIVAL,ERFEXP)
      END IF

C  NON-SCALED
      IF (DOLANC_GGALYP) THEN
          CALL VNCLALYP(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,VLAMBDA)
      END IF

C SCALED 

      IF (DOLASC_GGALYP) THEN
          CALL VSCLALYP(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,VLAMBDA)
      ENDIF

      IF (DOSRC_WIB) THEN
          CALL ELYP(ZKSR(1,1),RHO,RHO13,RHOGRD)
          CALL VLYP(ZKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP)
          IF (IWINT.EQ.1) THEN
             CALL VCSRLDA(ZMU,RHO,RHO13,CHIVAL,DOSRC_MULOCAL,1,ERFEXP)
  
             CALL WVWN(VLDA,RHO,RHO13,ELDA,.TRUE.,.TRUE.)
             ZKSR(1,1) = ZKSR(1,1) - ELDA*RHO
             ZKSR(2,1) = ZKSR(2,1) - VLDA
             CALL VSRWINT(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GE.2.OR.IWINT.EQ.3) THEN
             CALL VCSRLDA(ZMU,RHO,RHO13,CHIVAL,DOSRC_MULOCAL,1,ERFEXP)
             CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GT.10) THEN
             CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE
             CALL QUIT ('Wrong dft functional interpolation scheme'//
     &                  'in DFTEXC (IWINT).')
          END IF
      END IF
      IF (DOSRLYPT) THEN
          CALL ELYP(ZKSR,RHO,RHO13,RHOGRD)
          CALL GLYP(ZKSR(2,1),VL2,RHO,RHO13,RHOGRD)
          VXBT = VL2
      ENDIF
      IF (SRHYBR) THEN
#ifdef LDAembed
          CALL WVWN(VL,RHO,RHO13,EL,.TRUE.,.TRUE.)
          EL  = EL*RHO
          VL2 = 0.0D0
#else
          CALL ELYP(EL,RHO,RHO13,RHOGRD)
          CALL GLYP(VL,VL2,RHO,RHO13,RHOGRD)
#endif
          VXB = VXB + VL2
          ZKSR(1,2) = EL
          ZKSR(2,2) = VL
      ENDIF
      
      RETURN
      END


#ifdef JKP_DEBUG
C*****************************************************************************
      SUBROUTINE MCRHO(CMO,INDXCI,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate active density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "inforb.h"
#include "infvar.h"
#include "dftinf.h"
#include "dftcom.h"

      LOGICAL DOGGA,DOBCK,DOLYP,DOLND,DO_MOLGRAD,DODCAO
      DIMENSION CMO(*),WORK(LWORK)
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C     =================================
C     Initializations :
C     Number of AOs and their addresses
C     =================================
      DOLND = .FALSE.
      DOBCK = .FALSE.
      DOLYP = .FALSE.
      DO_MOLGRAD = .FALSE.
      DOGGA = DOBCK .OR. DOLYP
      
      NDER  = 0
      IF (DO_MOLGRAD) NDER = 1
      IF (DOGGA) THEN
         NDER = NDER + 1
         IF (DFTPOT) NDER = NDER + 1
      END IF

      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1= (NSOB1-1)*NBAST + 1
C     =============================
C     Allocations
C     =============================
      DODCAO = (NISHT.GT.0)
      NBUF   = 1000000
      KCREF  = 1
      KX     = KCREF + NCONF
      KY     = KX    + NBUF
      KZ     = KY    + NBUF
      KW     = KZ    + NBUF
      KGSO   = KW    + NBUF
      KCNT   = KGSO  + NTYPSO*NBAST
      KDGA   = KCNT  + NBAST
      KDV    = KDGA  + NBAST
      KDVAO  = KDV   + NNASHX
      KDCAO  = KDVAO + N2BASX
      IF (DODCAO) THEN
          KWRK  = KDCAO  + N2BASX
      ELSE
          KWRK  = KDCAO
      ENDIF
      LWRK  = LWORK + 1 - KWRK
      CALL MCRHO1(CMO,WORK(KDV),DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGA),DOGGA,WORK(KCREF),INDXCI,DODCAO,
     &            WORK(KDCAO),WORK(KDVAO),WORK(KWRK),LWRK,IPRINT)

      CALL TIMER('MCRHO',TIMSTR,TIMEND)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE MCRHO1(CMO,DV,DOLND,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,
     &                  NCNT,DMATGAO,DOGGA,CREF,INDXCI,DODCAO,DCAO,DVAO,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "infpar.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
      LOGICAL   DOLND,DOGGA,DODCAO
      DIMENSION CMO(*),DV(NNASHX),CORX(NBUF),CORY(NBUF),
     &          CORZ(NBUF),WEIGHT(NBUF),GSO(NBAST*NTYPSO),
     &          NCNT(NBAST),DMATGAO(NBAST),CREF(NCONF),
     &          DCAO(NBAST,NBAST),DVAO(NBAST,NBAST),WORK(LWORK), COR(3)
      CHARACTER*17 QUADNAME
#include "chrnos.h"

C     =============================
C     Open file for dumping density
C     =============================
      CALL GPOPEN(LUDDUMP,'DENSITY.GRID','UNKNOWN',
     &                ' ','FORMATTED',IDUMMY,.FALSE.)

C     =============================
C     Retrieve CREF
C     =============================
C
      REWIND LUIT1
      CALL MOLLAB('STARTVEC ',LUIT1,LUERR)
      DO 100 I = 1, (ISTACI-1)
         READ (LUIT1)
  100 CONTINUE
      CALL READT(LUIT1,NCONF,CREF)

C     =============================
C     Calculate DV
C     =============================
C
      CALL MAKDV(CREF,DV,INDXCI,WORK,LWORK)

C     =============================
C     Calculate inactive and active
C     density matrises in AO-basis
C     =============================
C
      CALL FCKDEN(DODCAO,.TRUE.,DCAO,DVAO,CMO,DV,WORK,LWORK)
C
C     =============================
C     Calculate grid
C     =============================

      IF (.NOT.DFTGRID_DONE_OLD) THEN
         CALL MAKE_DFTGRID(WORK,LWORK,NTOT_DFTGRID,1,.FALSE.)
         CALL CONDFT
         DFTGRID_DONE_OLD = .TRUE.
      END IF

      CELCTRN = D0
      VELCTRN = D0

      LUQUAD = -1

C     Make quadname: Can take 10000 procs

      IF (MYNUM .GT. 0) THEN
         QUADNAME = 'DALTON.QUAD.n'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C
         CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)
      ELSE
         CALL GPOPEN(LUQUAD,'DALTON.QUAD','OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      END IF

      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
#if defined (VAR_SRDFT_MPI)
         DO 300 IPNT = 1+MYNUM, NPOINT, NODTOT+1
#else
         DO 300 IPNT = 1, NPOINT
#endif
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F12.6)')
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF

            WGHT = WEIGHT(IPNT)

C           AOs
C           ===

            THRINT = DFTHRI/WGHT
            COR(1) = CORX(IPNT)
            COR(2) = CORY(IPNT)
            COR(3) = CORZ(IPNT)
            CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,
     &                  NBAST,DOLND,DOGGA,THRINT,IPRINT)

C           Density & Number of electrons
C           =============================

C           ... inactive part
            IF (DODCAO) THEN
               CALL GETRHO_OLD(DCAO,GSO,RHO,DMATGAO,THRINT,IPRINT)
               CELCTRN = CELCTRN + WGHT*RHO
            ENDIF
C           ... active part
            CALL GETRHO_OLD(DVAO,GSO,RHO,DMATGAO,THRINT,IPRINT)
            VELCTRN = VELCTRN + WGHT*RHO

C           Dump active density to file
C           =============================

            WRITE(LUDDUMP,'(4F25.13)') CORX(IPNT),CORY(IPNT),
     &           CORZ(IPNT),RHO

  300    CONTINUE

         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF

      CALL GPCLOSE(LUQUAD,'KEEP')

C     Test on the number of electrons

      ELCTRX = FLOAT(2*NRHFT)
      ELCTRN = CELCTRN + VELCTRN
      ERROR  = ELCTRN - ELCTRX

C     Print section

      WRITE (LUPRI,'(3(/2X,A,F14.6))')
     &' Nr. of inactive electrons from numerical integration :',CELCTRN,
     &' Nr. of active electrons from numerical integration   :',VELCTRN,
     &' Total nr. of electrons from numerical integration    :',ELCTRN
      WRITE (LUPRI,'(2(/2X,A,F14.6))')
     &' Number of electrons from orbial occupations          :',ELCTRX,
     &' Error in the number of electrons                     :',ERROR
      IF (ABS(ERROR) .GT. DFTELS) THEN
         WRITE (LUPRI,'(/2X,A,F14.6,/2X,A)')
     &' Error larger than DFTELS (set input)                 :',DFTELS,
     &   ' Calculation aborted.'
         CALL QUIT
     &    ('Wrong number of electrons in DFTDRV. Calculation aborted.')
      END IF

      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(/2X,A,F14.7,6X,D8.2,I14)')
     &      ' Number of electrons/abscissas:  ',
     &        ELCTRN,ELCTRN-ELCTRX,NPNTS
      END IF
      RETURN
      END
#endif  /* JKP_DEBUG */


C*****************************************************************************
      SUBROUTINE SRDFTLTR(KSYMOP,DTRMAT,EXCMAT,WGHT,GAO,GAO1,CHIVAL,
     &                    ERFEXP,LSRLDA,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                    RHG,DTGAO,VLAMBDA)
C*****************************************************************************
C
C     Jesper K. Pedersen. Nov. 2003
C     Purpose : Generate the 1-index transformed dft Hessian
C               needed for short-range DFT hybrids.
C     Based on DFTLTR used in regular DFT.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0,
     &           D13=1.0D0/3.0D0,
     &          RSFAC = 0.62035 04908 99400 08660 D0)
C               RSFAC = (D3/(D4*PI))**THIRD
      LOGICAL INVGGA,DOGGA,ERFEXP(0:2),LSRLDA
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
      DIMENSION DTRMAT(NBAST,NBAST),GAO(NBAST),GAO1(NBAST,3),
     &          EXCMAT(NBAST,NBAST),RHG(3),DTGAO(NBAST),
     &          DTGTAO(NBAST),XKSR(3,3),ZKSR(3,3),ZMU(3),BT3(3),BV3(3)
      DIMENSION B3(3)

      FRZTMP = 0D0
      FZZTMP = 0D0

C     Check for valid functional

      INVGGA = DOSRX_GGA.OR.DOSRBCK.OR.DOSRGGA2.OR.DOSRC_GGA.OR.DOSRLYPT
      DOGGA  = DOSRX_WIB.OR.DOSRC_WIB.OR.SRHYBR

      CALL DGEMV('N',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
      B0 = DDOT(NBAST,DTGAO,1,GAO,1)

      IF (INVGGA) THEN
         CALL QUIT('SRDFTLTR called with GGA type functional. '//
     &             'This is not implemented!')
      END IF

C     ---------------------------------------------------------
C     SRLDA,MULOCAL,SRHYBR or Interpolated GGA/SRLDA functional
C     ---------------------------------------------------------

      IF (DOGGA) THEN
C        B3 = GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
C        DTGAO= DTRMAT'*GAO
         CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
C        B3 = B3 + GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
         BMAX = MAX(abs(B0),abs(B3(1)),abs(B3(2)),abs(B3(3)))

         IF (IWINT .GT. 0) THEN ! used for DOSRX_WIB and DOSRC_WIB
            RS   = RSFAC/RHO13
            X = CHIVAL*RS*COPFAC
            DERFCX = DERFC(X)
            GAUX=(D2/SQRTPI)*X*EXP(-D13*X**2) ! for IWINT .eq. 3
         END IF
      ELSE
         BMAX = abs(B0)
         CALL DZERO(B3,3)
      END IF

      IF (BMAX.GT.DFTHRL) THEN

C        *****************************
C        ** Calculate exchange part **
C        *****************************

         IF (.NOT.DOSRC_MULOCAL(0)) THEN

            IF (LSRLDA) THEN
C        ...  Do not do SRLDA if SRDFTLTR is called with LSRLDA false
C             (regular Kohm-Sham part with total density for SRHYBR).
              CALL VXSRLDA(XKSR,RHO,RHO13,CHIVAL,2,ERFEXP)
            END IF

C           Below are various corrections to the SRLDA energy/potential

            IF (DOSRX_WIB) THEN
C         ... Special case for Weighted Interpolation functionals
              CALL DCOPY(3,XKSR,1,ZMU,1)
              CALL DZERO(XKSR,3)
              CALL EBCK(XKSR(1,1),RHO,RHO13,RHOGRD)
              CALL VBCK(XKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
              CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,
     &                   FZZBCK,RHO,RHOGRD)
              XKSR(3,1) = FRRBCK
              IF (IWINT.EQ.1) THEN
                 CALL VSRWINT(XKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
              ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3.OR.IWINT.GT.10) THEN
                 IF (IWINT.GT.10) CALL DZERO(ZMU,3)
                 CALL EDRC(ED,RHO,RHO13)
                 CALL VDRC(VD,RHO13)
                 CALL V1DRC(VDRC0,VDRC1,RHO,RHO13)
                 XKSR(1,1) = XKSR(1,1) + ED
                 XKSR(2,1) = XKSR(2,1) + VD
                 XKSR(3,1) = XKSR(3,1) + VDRC1
                 CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
              ELSE
                 write(lupri,*) 'IWINT : ',IWINT
                 CALL QUIT ('Wrong dft functional interpolation '//
     &                      'scheme in SRDFTLTR (IWINT).')
              END IF
C             --- Remember mixed terms from Becke
              IF (IWINT.EQ.3) THEN
                 FRZTMP = FRZBCK*GAUX
                 FZZTMP = FZZBCK*GAUX
              ENDIF
              FRZBCK = FRZBCK*DERFCX + FRZTMP
              FZZBCK = FZZBCK*DERFCX + FZZTMP

            ELSE IF (SRHYBR) THEN
C          ... Special case for SRDFT in DFT embedding.
              CALL V1DRC(VDRC0,FRRDRC,RHO,RHO13)
#ifndef LDAembed
              CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,FZZBCK,
     &                   RHO,RHOGRD)
              IF (LSRLDA) THEN
                 XKSR(3,1) =  XKSR(3,1) - FRRDRC - FRRBCK
                 FR0BCK    = -FR0BCK
                 FZ0BCK    = -FZ0BCK
                 FRRBCK    = -FRRBCK
                 FRZBCK    = -FRZBCK
                 FZZBCK    = -FZZBCK
              ELSE
                  XKSR(3,1) = FRRDRC + FRRBCK
              END IF
#else
              FR0BCK = 0D0
              FZ0BCK = 0D0
              FRRBCK = 0D0
              FRZBCK = 0D0
              FZZBCK = 0D0
#endif
            ELSE
              FR0BCK = 0D0
              FZ0BCK = 0D0
              FRRBCK = 0D0
              FRZBCK = 0D0
              FZZBCK = 0D0
            ENDIF
         ENDIF  ! IF (.NOT.DOSRC_MULOCAL(0)) THEN

C        ********************************
C        ** Calculate correlation part **
C        ********************************

         IF (LSRLDA) THEN
C     ... Don't do SRLDA if SRDFTLTR is called with LSRHYBR true
C          (regular Kohm-Sham part with total density for SRHYBR).
           CALL VCSRLDA(ZKSR,RHO,RHO13,CHIVAL,DOSRC_MULOCAL,2,ERFEXP)
         END IF
C
C     ... Below are various corrections to the SRLDA energy/potential
C
         IF (DOSRC_WIB) THEN
C        ... Special case for Weighted Interpolation functionals
           CALL DCOPY(3,ZKSR,1,ZMU,1)
           CALL DZERO(XKSR,3)
C          -- lyp
           CALL ELYP(ZKSR(1,1),RHO,RHO13,RHOGRD)
           CALL VLYP(ZKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP)
           RHOA = DP5*RHO
           RHGA = (DP5*RHOGRD)**2
           CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
     &                 DF00001,RHO,RHO13,RHOGRD,.TRUE.)
           CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
     &                 DF0101,DF1001,DF0110,DF10001,
     &                 DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
           FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
           FRRLYP = DP5*(DF2000 + DF1100)
           FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
           FZZLYP = FZ0LYP/RHOGRD
           ZKSR(3,1) = FRRLYP
C
           IF (IWINT.EQ.1) THEN
              CALL WVWN(VLDA,RHO,RHO13,ELDA,.TRUE.,.TRUE.)
              CALL V1VWN(VVWN0,FRRVWN,RHO,RHO13)
              ZKSR(1,1) = ZKSR(1,1) - ELDA*RHO
              ZKSR(2,1) = ZKSR(2,1) - VLDA
              ZKSR(3,1) = ZKSR(3,1) - FRRVWN
              CALL VSRWINT(ZKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
           ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3 .OR. IWINT.GT.10) THEN
              IF (IWINT.GT.10) CALL DZERO(ZMU,3)
              CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
           ELSE
              CALL QUIT ('Wrong dft functional interpolation '//
     &       'scheme in correlation part of SRDFTLTR (IWINT).')
           END IF
C          --- Remember mixed terms from lyp
           IF (IWINT.EQ.3) THEN
              FRZTMP = FRZLYP*GAUX
              FZZTMP = FZZLYP*GAUX
           ENDIF
           FRZLYP = FRZLYP*DERFCX + FRZTMP
           FZZLYP = FZZLYP*DERFCX + FZZTMP

         ELSE IF (SRHYBR) THEN
C        ... Special case for SRDFT in DFT embedding.
            RHOA = DP5*RHO
            RHGA = (DP5*RHOGRD)**2

#ifdef LDAembed
            CALL V1VWN(VVWN0,FRRLYP,RHO,RHO13)
            FZ0LYP = 0.0D0
            FRZLYP = 0.0D0
            FZZLYP = 0.0D0
#else
            CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
     &                  DF00001,RHO,RHO13,RHOGRD,.TRUE.)
            CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
     &                  DF0101,DF1001,DF0110,DF10001,
     &                  DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
            FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
            FRRLYP = DP5*(DF2000 + DF1100)
            FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
            FZZLYP = FZ0LYP/RHOGRD

#endif
            IF (LSRLDA) THEN
               ZKSR(3,1) =  ZKSR(3,1) - FRRLYP
               FZ0LYP = -FZ0LYP
               FRRLYP = -FRRLYP
               FRZLYP = -FRZLYP
               FZZLYP = -FZZLYP
            ELSE
               ZKSR(3,1) =  FRRLYP
            END IF

         ELSE
            FRRVWN = 0D0
            FZ0LYP = 0D0
            FRRLYP = 0D0
            FRZLYP = 0D0
            FZZLYP = 0D0
         END IF

C        ************************
C        ** Total contribution **
C        ************************

#ifdef LDAembed
         ZNV     = 0.d0
#else
         ZNV     = D1/RHOGRD
#endif
         FRR = WGHT*(ZKSR(3,1) + XKSR(3,1))
         FZ0 = ZNV*WGHT*(FZ0BCK + FZ0LYP)
         FRZ = WGHT*(FRZBCK + FRZLYP)
         FZZ = WGHT*(FZZBCK + FZZLYP)

         RX = ZNV*RHG(1)
         RY = ZNV*RHG(2)
         RZ = ZNV*RHG(3)

C        ***************************
C        ** Linear transformation **
C        ***************************

         IF (.NOT.SRHYBR .OR. DOHFSRDFT) THEN
!        i.e. for normal MCSRDFT (with SRHYBR .false.) and always for
!        HFSRDFT
            IF (DOGGA) THEN
               BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
               FAC0 = FRR*B0 + FRZ*BR
               FACR = FRZ*B0 + FZZ*BR
            ELSE
               FAC0 = FRR*B0
            END IF
            IF (NSYM.EQ.1) THEN
               DO I = 1, NBAST
                  G0 = GAO(I)
                  IF (DOGGA) THEN
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                  END IF
                  DO J = 1, I
                     EXCMAT(J,I) = EXCMAT(J,I)+FAC0*G0*GAO(J)
                     IF (DOGGA) THEN
                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR = AX*RX + AY*RY + AZ*RZ
                        AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
                        EXCMAT(J,I) = EXCMAT(J,I)+FACR*AR+FZ0*AB
                     END IF
                  END DO
               END DO
            ELSE
               DO ISYM = 1, NSYM
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  JSYM = MULD2H(ISYM,KSYMOP)
                  IF (ISYM.GE.JSYM) THEN
                     JSTR = IBAS(JSYM) + 1
                     JEND = IBAS(JSYM) + NBAS(JSYM)
                     DO I = ISTR, IEND
                        G0 = GAO(I)
                        IF (DOGGA) THEN
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                        END IF
                        DO J = JSTR, MIN(I,JEND)
                           EXCMAT(J,I) = EXCMAT(J,I) + FAC0*G0*GAO(J)
                           IF (DOGGA) THEN
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR = AX*RX + AY*RY + AZ*RZ
                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -
     &                             AR*BR
                              EXCMAT(J,I) = EXCMAT(J,I)
     &                                     + FACR*AR + FZ0*AB
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END IF
         END IF   ! IF (.NOT.SRHYBR .OR. DOHFSRDFT) THEN
!        i.e. for normal MCSRDFT (with SRHYBR .false.) and always for
      END IF   ! IF (BMAX.GT.DFTHRL) THEN

      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRFMAT(FMAT,CMO,DVREF,EJCSR,EJVSR,EDSR,ESRDFT,
     &                  EMYDFTAUX,UEJCVSR,WRK,LWRK,IPRINT)
C*****************************************************************************
C
C     20-Mar-2003 Jesper K. Pedersen
C
C     Driver for Fock matrix construction of short-range terms
C     needed for CI-DFT (J_ee,E_xc)
C
C     EJCSR =  0.5 * Tr (D^c . J^sr . D^c)
C     EJVSR = -0.5 * Tr (D^{v,ref} . J^sr . D^{v,ref})
C     EDSR  = -1.0 * Tr (V^{sr,ref}_{xc} . D^{v,ref})
C     ESRDFT  = E^{sr}_{xc} [D^c+D^{v,ref}] 
C
C     23-Aug-2012 modified by Manu 
C   
C     * added a new argument, EMYDFTAUX, that contains the inactive
C       short-range contributions to auxiliary energy:
C
C       EMYDFTAUX = 1.0 * Tr (D^c . J^sr . D^c)
C                  +1.0 * Tr (D^{v,ref} . J^sr . D^c)
C                  +1.0 * Tr (V^{sr,ref}_{xc} . D^c) 
C       UEJCVSR   = 1.0 * Tr (D^{v,updated} . J^sr . D^c)
C
C     * Note that for ensemble CI-srDFT, D^{v,ref} should be, here, equal to the
C       ensemble density matrix
C*****************************************************************************
#include "implicit.h"
      REAL*8    FMAT(*), CMO(*), DVREF(*), WRK(*), ESRDFT(3)
      INTEGER   ISYMDM(2), IFCTYP(2)
      LOGICAL   LSRHYBR, DFT_SPINDNS_SAVE
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0, DM1 = -1.0D0)
      SAVE ISYMDM, IFCTYP
      DATA ISYMDM/1,1/, IFCTYP/999,999/

#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infvar.h"
#include "inftap.h"
#include "inforb.h"
#include "infinp.h"
#include "infpri.h"
#include "dftcom.h"

      CALL QENTER('SRFMAT')

      DFT_SPINDNS_SAVE = DFT_SPINDNS
      DFT_SPINDNS = .FALSE.
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KDCAO,2*N2BASX,WRK,KFREE,LFREE)
      KDVAO = KDCAO + N2BASX
      IF (LSRHYBR) THEN
         NF = 3
      ELSE
         NF = 2
      END IF
      CALL MEMGET('REAL',KFCAO,NF*N2BASX,WRK,KFREE,LFREE)
      KFVAO = KFCAO + N2BASX

      CALL FCKDEN(.TRUE.,.TRUE.,
     &            WRK(KDCAO),WRK(KDVAO),CMO,DVREF,WRK(KFREE),LFREE)

      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A/)')
     &   '----  Valence AO density matrix ----'
         CALL OUTPUT(WRK(KDVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C ... Get J^c,sr and J^v,sr
C
      IF(HFXFAC .NE. D0) THEN
C        Coulomb + Exchange
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      ELSE
C        only Coulomb (Hartree) term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ENDIF
      CALL DZERO(WRK(KFCAO),NF*N2BASX)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFCAO),WRK(KDCAO),2,
     &             ISYMDM,IFCTYP,WRK(KFREE),LFREE)
      EJCSR   =  DP5*DDOT(N2BASX,WRK(KDCAO),1,WRK(KFCAO),1)
      EJVSR   = -DP5*DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
      UEJCVSR =      DDOT(N2BASX,WRK(KDCAO),1,WRK(KFVAO),1)
     
C Manu: compute first and second terms in EMYDFTAUX

      EMYDFTAUX = 2.0D0*EJCSR 
      EMYDFTAUX = EMYDFTAUX 
     &         +DDOT(N2BASX,WRK(KDCAO),1,WRK(KFVAO),1)

      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A/)')
     &   '----  Short-range core Hartree matrix ----'
         CALL OUTPUT(WRK(KFCAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range valence Hartree matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF

C ... Add up in KFCAO

      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)

C ... Get V^sr in WRK(KFVAO) and E_xc^sr
C     using total density matrix (stored in KDCAO)
      CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDCAO),1)

      CALL DZERO(WRK(KFVAO),N2BASX)

C     SRDFT(ND_SIM,EXCMAT,DMAT,ESRDFT(1:3),DOERG,DO_MOLGRAD,DOATR,WORK,LWORK,IPRINT)
      CALL SRDFT(1,WRK(KFVAO),WRK(KDCAO),ESRDFT,
     &           .TRUE.,.FALSE.,.FALSE.,.FALSE.,
     &           WRK(KFREE),LFREE,IPRFCK)
      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range Exchange-Correlation Fock-matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF
      EDSR = -DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)

C Manu: compute the third term in EMYDFTAUX.
C       Note that WRK(KDCAO) actually contains the TOTAL density matrix
C       and we only need the core contribution.

      EMYDFTAUX = EMYDFTAUX
     &         +DDOT(N2BASX,WRK(KDCAO),1,WRK(KFVAO),1)
     &         +EDSR

      IF (LSRHYBR) THEN
         write (lupri,*) 'WARNING, srdft-hybr exc matrix not used'
C TODO
      END IF

C ... Add the Kohn-Sham Fock matrix

      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)
C
      IF (DOSRBCK) THEN
C     ... Get K^lr using total density and subtract it!
          IFCTYP(1) = 12
C         ... only exchange
          CALL DZERO(WRK(KFVAO),N2BASX)
          CALL SIRFCK2(LUINTA,'AOTWOINT',WRK(KFVAO),WRK(KDCAO),1,
     &                 ISYMDM,IFCTYP,WRK(KFREE),LFREE)
          EKVLR =     DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
          CALL DAXPY(N2BASX,DM1,WRK(KFVAO),1,WRK(KFCAO),1)
          EKTLR = DP5*DDOT(N2BASX,WRK(KDCAO),1,WRK(KFVAO),1)
          EDSR = EDSR + EKVLR
          ESRDFT(1) = ESRDFT(1) - EKTLR
      ENDIF

C     ===============================================
C      Transform short-range Fock-matrix to MO basis
C     ===============================================

      CALL DCOPY(N2BASX,WRK(KFCAO),1,WRK(KDVAO),1)
      CALL DGETSP(NBAST,WRK(KDVAO),WRK(KFCAO))
      IF (NSYM.GT.1) CALL PKSYM1(WRK(KFCAO),WRK(KFCAO),NBAS,NSYM,2)
      CALL UTHUB(WRK(KFCAO),FMAT,CMO,WRK(KDVAO),NSYM,NBAS,NORB)

      IF (P6FLAG(14)) THEN
        WRITE(LUPRI,1200)
        CALL OUTPKB(FMAT,NORB,NSYM,1,LUPRI)
      END IF
 1200 FORMAT(/' Short-range Fock matrix (MO-basis)')

      CALL MEMREL('SRFMAT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      DFT_SPINDNS = DFT_SPINDNS_SAVE
      CALL QEXIT('SRFMAT')
      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRFMAT_vensemble(FMAT,CMO,UDVREF,DVREF,DIFFDVREF,
     &                            diff_ejvsr,diff_edsr,uEDFT,
     &                            WRK,LWRK,IPRINT)
C*****************************************************************************
C
C     !> brief: calculate various terms needed to make a variational ensemble energy
C     !> author: S. Knecht and E. Fromager based on the routine by J. K. Pedersen
C     
C     EJCSR =  0.5 * Tr (D^c . J^sr . D^c)
C     EJVSR = -0.5 * Tr (D^{v,ref} . J^sr . D^{v,ref})
C     EDSR  = -1.0 * Tr (V^{sr,ref}_{xc} . D^{v,ref})
C     uEDFT(1)  = E^{sr}_{xc} [D^c+D^{v,ref}]
C     uEDFT(2:3) = exchange and correlation contributions, respectively
C
C     !> input: ensemble density matrices D^{v,ref}, D^{v,upd} and D^{v,diff} = D^{v,upd}-D^{v,ref}
C*****************************************************************************
 
#include "implicit.h"
      REAL*8    FMAT(*),CMO(*),UDVREF(*),DVREF(*),DIFFDVREF(*),WRK(*)
      REAL*8    uEDFT(3)
      LOGICAL   LSRHYBR, DFT_SPINDNS_SAVE
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0, DM1 = -1.0D0)
      INTEGER   ISYMDM(2), IFCTYP(2)
      SAVE ISYMDM, IFCTYP
      DATA ISYMDM/1,1/, IFCTYP/999,999/

#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infvar.h"
#include "inftap.h"
#include "inforb.h"
#include "infinp.h"
#include "infpri.h"
#include "dftcom.h"

      CALL QENTER('SRFMAT_vensemble')

      DFT_SPINDNS_SAVE = DFT_SPINDNS
      DFT_SPINDNS = .FALSE.
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KDCAO,2*N2BASX,WRK,KFREE,LFREE)
      KDVAO = KDCAO + N2BASX
      CALL MEMGET('REAL',KDCAO_backup,2*N2BASX,WRK,KFREE,LFREE)
      KDVAO_backup = KDCAO_backup + N2BASX
      IF (LSRHYBR) THEN
         NF = 3
      ELSE
         NF = 2
      END IF
      CALL MEMGET('REAL',KFCAO,NF*N2BASX,WRK,KFREE,LFREE)
      KFVAO = KFCAO + N2BASX

      !> start with computing uedft - requires u_rho_ensemble

      CALL FCKDEN(.TRUE.,.TRUE.,
     &            WRK(KDCAO),WRK(KDVAO),CMO,UDVREF,WRK(KFREE),LFREE)

      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Valence AO density matrix ----'
         CALL OUTPUT(WRK(KDVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF

C ... Get J^c,sr adn J^v,sr

      IF(HFXFAC .NE. D0) THEN
C        Coulomb + Exchange
C     ... Warning : assuming DIM(IFCTYP) = 2
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      ELSE
C        only Coulomb (Hartree) term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ENDIF
      CALL DZERO(WRK(KFCAO),NF*N2BASX)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFCAO),WRK(KDCAO),2,
     &             ISYMDM,IFCTYP,WRK(KFREE),LFREE)
      EJVSR   = -DP5*DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
C
      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range core Hartree matrix ----'
         CALL OUTPUT(WRK(KFCAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range valence Hartree matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF

c ... Add up in KFCAO

      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)

C ... Get V^sr in WRK(KFVAO) and E_xc^sr
C     using total density matrix (stored in KDCAO)
      CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDCAO),1)

      CALL DZERO(WRK(KFVAO),N2BASX)

C     SRDFT(ND_SIM,EXCMAT,DMAT,ESRDFT(1:3),DOERG,DO_MOLGRAD,DOATR,WORK,LWORK,IPRINT)
      CALL SRDFT(1,WRK(KFVAO),WRK(KDCAO),uEDFT,
     &           .TRUE.,.FALSE.,.FALSE.,.FALSE.,
     &           WRK(KFREE),LFREE,IPRFCK)
      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range Exchange-Correlation Fock-matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF
      EDSR = -DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)

      IF (LSRHYBR) THEN
         write (lupri,*) 'WARNING, srdft-hybr exc matrix not used'
C TODO
      END IF

      !> part 2: compute diff_edsr / requires DVREF and DIFFDVREF

      call dzero(WRK(KDCAO),N2BASX)
      call dzero(WRK(KDVAO),N2BASX)

      CALL FCKDEN(.TRUE.,.TRUE.,
     &            WRK(KDCAO),WRK(KDVAO),CMO,DVREF,WRK(KFREE),LFREE)

      call dzero(WRK(KDCAO_backup),N2BASX)
      call dzero(WRK(KDVAO_backup),N2BASX)

      CALL FCKDEN(.TRUE.,.TRUE.,
     &            WRK(KDCAO_backup),WRK(KDVAO_backup),CMO,DIFFDVREF,
     &            WRK(KFREE),LFREE)

C ... Get J^c,sr adn J^v,sr

      IF(HFXFAC .NE. D0) THEN
C        Coulomb + Exchange
C     ... Warning : assuming DIM(IFCTYP) = 2
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      ELSE
C        only Coulomb (Hartree) term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ENDIF
      CALL DZERO(WRK(KFCAO),NF*N2BASX)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFCAO),WRK(KDCAO),2,
     &             ISYMDM,IFCTYP,WRK(KFREE),LFREE)

      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range core Hartree matrix ----'
         CALL OUTPUT(WRK(KFCAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range valence Hartree matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF

C ... Add up in KFCAO

      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)

C ... Get V^sr in WRK(KFVAO) and E_xc^sr
C     using total density matrix (stored in KDCAO)
      CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDCAO),1)

      CALL DZERO(WRK(KFVAO),N2BASX)

C     SRDFT(ND_SIM,EXCMAT,DMAT,ESRDFT(1:3),DOERG,DO_MOLGRAD,DOATR,WORK,LWORK,IPRINT)
      CALL SRDFT(1,WRK(KFVAO),WRK(KDCAO),uEDFT,
     &           .TRUE.,.FALSE.,.FALSE.,.FALSE.,
     &           WRK(KFREE),LFREE,IPRFCK)
      IF (IPRINT.GT.9) THEN
         WRITE(LUPRI,'(//5X,A//)')
     &   '----  Short-range Exchange-Correlation Fock-matrix ----'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      ENDIF
      EDSR = -DDOT(N2BASX,WRK(KDVAO_backup),1,WRK(KFVAO),1)

      IF (LSRHYBR) THEN
         write (lupri,*) 'WARNING, srdft-hybr exc matrix not used'
C TODO
      END IF

      diff_edsr = EDSR

      !> step 3: compute diff_ejvsr / requires DIFFDVREF (sitting in _backup)


C ... Get J^c,sr and J^v,sr

      IF(HFXFAC .NE. D0) THEN
C        Coulomb + Exchange
C     ... Warning : assuming DIM(IFCTYP) = 2
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      ELSE
C        only Coulomb (Hartree) term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ENDIF
      CALL DZERO(WRK(KFCAO),NF*N2BASX)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFCAO),WRK(KDCAO_backup),2,
     &             ISYMDM,IFCTYP,WRK(KFREE),LFREE)
      EJVSR   = -DP5*DDOT(N2BASX,WRK(KDVAO_backup),1,WRK(KFVAO),1)

      IF (LSRHYBR) THEN
         write (lupri,*) 'WARNING, srdft-hybr exc matrix not used'
C TODO
      END IF

      diff_ejvsr = EJVSR

      CALL MEMREL('SRFMAT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      DFT_SPINDNS = DFT_SPINDNS_SAVE
      CALL QEXIT('SRFMAT_vensemble')
      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRDFTSO(DCVAL,UDV,V2C_SR,EOVEC)
C*****************************************************************************
C
C  Copyright Dec-2003 Hans Joergen Aa. Jensen
C  for MCSRDFT (based on SOLGO code).
C
C  Purpose: Add V2C_SR orbital sigma vector contribution to EOVEC.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION UDV(NASHT,*), V2C_SR(NORBT,*), EOVEC(*)
C  Used from common blocks:
C    INFORB: NISHT,NASHT
C    INFVAR: NWOPH,JWOP(2,*)
C    INFIND: IOBTYP(*),ISX(*)
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infvar.h"
#include "infind.h"
#include "dftcom.h"
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )

      CALL QENTER('SRDFTSO ')
      DCFAC = D2 * DCVAL
C     ... DCVAL is the value of DC = <0| Eii |B>,
C         = 2 <0|0> = 2 for singlet density matrix,
C         = 2 <0|B> = 0 for singlet transition density matrix,
C         = 0           for triplet density matrix
      DO 300 IG = 1,NWOPH
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF (ITYPK .EQ. JTINAC) THEN
C           first index inactive
            EOVEC(IG) = EOVEC(IG) + DCFAC * V2C_SR(K,L)
         ELSE
C           first index active
            NK   = ICH(K)
            TEMP = D0
            DO 100 NX = 1,NASHT
               DVKX = UDV(NX,NK)
               IX = ISX(NISHT + NX)
               TEMP = TEMP + DVKX*V2C_SR(IX,L)
  100       CONTINUE
            EOVEC(IG) = EOVEC(IG) + D2 * TEMP
         END IF
         IF (ITYPL .EQ. JTACT) THEN
C           second index active
            NL   = ICH(L)
            TEMP = D0
            DO 200 NX = 1,NASHT
               DVLX = UDV(NX,NL)
               IX = ISX(NISHT + NX)
               TEMP = TEMP + DVLX*V2C_SR(IX,K)
  200       CONTINUE
            EOVEC(IG) = EOVEC(IG) - D2 * TEMP
         END IF
  300 CONTINUE

      CALL QEXIT('SRDFTSO ')
      RETURN
C     end of srdftso.
      END


C*****************************************************************************
      SUBROUTINE RSP_SRDFTSO(UDV,V2C_SR,EVEC)
C*****************************************************************************
C
C  Copyright June 2010 Hans Joergen Aa. Jensen
C
C  for MCSRDFT response
C  (based on SRDFTSO and SLVSOR(rspsol.F) oneind=.true.  code).
C
C  Purpose: Add V2C_SR orbital sigma vector contribution to EVEC.
C
C        <0| [q_KL, V^([2c],xc-SR) ] |0>
C
C  UDV    : active density matrix of reference state |0>
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION UDV(NASHDI,*), V2C_SR(NORBT,*), EVEC(*)
C  Used from common blocks:
C    INFORB: NISHT,NASHT
C    INFVAR: NWOPH,JWOP(2,*)
C    INFIND: IOBTYP(*),ISX(*)
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infvar.h"
#include "infind.h"
#include "infdim.h"
#include "dftcom.h"
#include "infrsp.h"
#include "wrkrsp.h"
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )

      CALL QENTER('RSP_SRDFTSO ')

#if SRDFT_DEBUG > 15
      write (lupri,*) 'rsp_srdftso: input UDV:'
      call output(udv,1,nasht,1,nasht,nasht,nasht,-1,lupri)
      write (lupri,*) 'rsp_srdftso: input V2C_SR:'
      call output(v2c_sr,1,norbt,1,norbt,norbt,norbt,-1,lupri)
      write (lupri,*) 'rsp_srdftso:  input EVEC orbital part:'
      call output(evec(kzconf+1),1,kzwopt,1,2,kzvar,2,1,lupri)
#endif

      KYCONF = KZCONF + KZVAR

C     Distribute srDFT Fock-type matrix contribution to orbital
C     part of EVECS

      KSYM1 = -1
      DO IG = 1,KZWOPT
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         KSYM  = ISMO(K)
         LSYM  = ISMO(L)
         IF ( KSYM.NE.KSYM1 ) THEN
            KSYM1 = KSYM
            IORBK = IORB(KSYM)
            NASHK = NASH(KSYM)
            NISHK = NISH(KSYM)
            IASHK = IASH(KSYM)
            IORBL = IORB(LSYM)
            NASHL = NASH(LSYM)
            NISHL = NISH(LSYM)
            IASHL = IASH(LSYM)
         END IF
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF (ITYPK .EQ. JTINAC) THEN ! first index inactive
            EVEC(KZCONF+IG) = EVEC(KZCONF+IG) + D2*V2C_SR(L,K)
           IF (.NOT. (TDA .OR. CISRPA)) THEN
            EVEC(KYCONF+IG) = EVEC(KYCONF+IG) - D2*V2C_SR(K,L)
           END IF
         ELSE                        ! first index active
            NK    = ICH(K)
            TEMPZ = D0
            TEMPY = D0
            DO NX = 1,NASHK
               DVKX = UDV(NK,IASHK+NX)
               DVXK = UDV(IASHK+NX,NK)
               IX = (IORBK+NISHK)+NX
               TEMPZ = TEMPZ + DVKX*V2C_SR(L,IX)
               TEMPY = TEMPY + DVXK*V2C_SR(IX,L)
            END DO
            EVEC(KZCONF+IG) = EVEC(KZCONF+IG) + TEMPZ
           IF (.NOT. (TDA .OR. CISRPA)) THEN
            EVEC(KYCONF+IG) = EVEC(KYCONF+IG) - TEMPY
           END IF
         END IF
         IF (ITYPL .EQ. JTACT) THEN ! second index active
            NL    = ICH(L)
            TEMPZ = D0
            TEMPY = D0
            DO NX = 1,NASHL
               DVLX = UDV(NL,IASHL+NX)
               DVXL = UDV(IASHL+NX,NL)
               IX = (IORBL+NISHL)+NX
               TEMPZ = TEMPZ + DVXL*V2C_SR(IX,K)
               TEMPY = TEMPY + DVLX*V2C_SR(K,IX)
            END DO
            EVEC(KZCONF+IG) = EVEC(KZCONF+IG) - TEMPZ
           IF (.NOT. (TDA .OR. CISRPA)) THEN
            EVEC(KYCONF+IG) = EVEC(KYCONF+IG) + TEMPY
           END IF
         END IF
      END DO ! DO IG = 1, KZWOPT

#if SRDFT_DEBUG > 15
      write (lupri,*) 'rsp_srdftso: output EVEC orbital part:'
      call output(evec(kzconf+1),1,kzwopt,1,2,kzvar,2,1,lupri)
#endif
      CALL QEXIT('RSP_SRDFTSO ')
      RETURN
C     end of rsp_srdftso.
      END
C --- end of srdft/srdft.F ---
