*=====================================================================*
      SUBROUTINE CCFTRANSC(LISTL, IDLSTL, LISTR, IDLSTR, WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: shortcut to F matrix transformation for a single
*             transformation (should be avoided and maybe completely
*             eliminated)
*        
*     Written by Christof Haettig, Januar 1999, based on CC_FMATRIX.
*
*     Trash: Debug copy for FBTA. Sonia Coriani, 22/09-1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include "maxorb_dec.h"
#  include "ccorb_dec.h"
#  include "ccsdsym_dec.h"
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include "maxorb.h"
#include "ccorb.h"
#include "ccsdsym.h"

      INTEGER LUFMAT

      CHARACTER*(*) LISTL, LISTR
      CHARACTER*(8) FILFMA
      INTEGER IOPTRES, IDLSTL, IDLSTR, IDUM, NFTRAN, LWORK
      INTEGER IFTRAN(3,1)
      INTEGER KIADRBFD, KIADRZ0, KIADRPQ0, KIADRPQR, KIADRPQMO   
      INTEGER KIADRBFI, KRHO1, KRHO2, KFIM, KMGDL, KRIM, KZDEN, KGIM 
      INTEGER KLAMPA, KLAMHA, KCTR2, KLAMP, KLAMH, KDENS, KFOCK
      INTEGER KEND, LEND, ISYML, ISYMR, ISYRES, LENALL
      
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL WORK(LWORK), RDUM
#else 
      DOUBLE PRECISION WORK(LWORK), RDUM
#endif 

      NFTRAN      = 1
      IFTRAN(1,1) = IDLSTL
      IFTRAN(2,1) = IDLSTR
      IOPTRES     = 0
      FILFMA      = 'CC__FMAT'

      KIADRBFD    = 1
      KIADRZ0     = KIADRBFD    + MAXORB*NFTRAN
      KIADRPQ0    = KIADRZ0     + MAXORB*NFTRAN
      KIADRPQR    = KIADRPQ0    + MAXORB*NFTRAN
      KIADRPQMO   = KIADRPQR    + MAXORB*NFTRAN
      KIADRBFI    = KIADRPQMO   + MAXORB*NFTRAN
      KRHO1       = KIADRBFI    + MAXORB*NFTRAN
      KRHO2       = KRHO1       + NFTRAN
      KFIM        = KRHO2       + NFTRAN
      KMGDL       = KFIM        + NFTRAN
      KRIM        = KMGDL       + NFTRAN
      KZDEN       = KRIM        + NFTRAN
      KGIM        = KZDEN       + NFTRAN
      KLAMPA      = KGIM        + NFTRAN
      KLAMHA      = KLAMPA      + NFTRAN
      KCTR2       = KLAMHA      + NFTRAN
      KLAMP       = KCTR2       + NFTRAN
      KLAMH       = KLAMP       + NFTRAN
      KDENS       = KLAMH       + NFTRAN
      KFOCK       = KDENS       + NFTRAN
      KEND        = KFOCK       + NFTRAN
      LEND        = LWORK       - KEND

      IF (LEND .LT. 0) THEN
         WRITE (LUPRI,*) 'Insufficient work space in CCFTRANSC.'
         WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
         WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
         CALL QUIT('Insufficient work space in CCFTRANSC.')
      END IF

      CALL CCFMATNEWSC( IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES, FILFMA, 
     &                 IDUM, RDUM, 
     &                 WORK(KIADRBFD), WORK(KIADRZ0), WORK(KIADRPQ0),
     &                 WORK(KIADRPQR),WORK(KIADRPQMO),WORK(KIADRBFI),
     &                 WORK(KRHO1),WORK(KRHO2),WORK(KFIM),WORK(KMGDL),
     &                 WORK(KRIM),WORK(KZDEN),WORK(KGIM),WORK(KLAMPA),
     &                 WORK(KLAMHA), WORK(KCTR2), WORK(KLAMP),
     &                 WORK(KLAMH), WORK(KDENS), WORK(KFOCK),
     &                 0, WORK(KEND), LEND )

      ISYML = ILSTSYM(LISTL,IDLSTL)
      ISYMR = ILSTSYM(LISTR,IDLSTR)
      ISYRES = MULD2H(ISYML,ISYMR)

      LENALL = NT1AM(ISYRES) + NT2AM(ISYRES)
      IF (LENALL.GT.LWORK) THEN
         WRITE (LUPRI,*) 'Insufficient work space in CCFTRANSC.'
         WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
         WRITE (LUPRI,*) 'Need at least:',LENALL, ' words.'
         CALL QUIT('Insufficient work space in CCFTRANSC.')
      END IF

      LUFMAT = 0
      CALL WOPEN2(LUFMAT, FILFMA, 64, 0)
      CALL GETWA2(LUFMAT,FILFMA,WORK(1),1,LENALL) 
      CALL WCLOSE2(LUFMAT, FILFMA, 'DELETE')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCFTRANSC                            *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
c/* Deck CC_FMATNEW */
*=====================================================================*
      SUBROUTINE CCFMATNEWSC(IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES,
     &                      FILFMA, IFDOTS, FCONS, 
     &                      IADRBFD, IADRZ0, IADRPQ0, IADRPQR,
     &                      IADRPQMO, IADRBFI, KRHO1, KRHO2, KFIM, 
     &                      KMGDL, KRIM, KZDEN, KGIM, KLAMPA,
     &                      KLAMHA, KCTR2, KLAMP, KLAMH,
     &                      KDENS, KFOCK, MXVEC, WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: AO-direct calculation of a linear transformation of two
*             CC amplitude vectors, L and R, with the CC F matrix
*             (or B matrix for L different from zeroth-order multipl.)
*
*                L.EQ.'L0'  --->      F * R       
*                L.NE.'L0'  --->  L * B * R       
*          
*             the linear transformations are calculated for a list
*             of L vectors and a list of R vectors: 
*
*                LISTL       -- type of L vectors
*                LISTR       -- type of R vectors
*                IFTRAN(1,*) -- indeces of L vectors
*                IFTRAN(2,*) -- indeces of R vectors
*                IFTRAN(3,*) -- indeces or addresses of result vectors
*                NFTRAN      -- number of requested transformations
*                FILFMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IFDOTS      -- indeces of vectors to be dotted on
*                FCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILFMA is used as file name
*                          the start addresses of the vectors are
*                          returned in IFTRAN(3,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IFTRAN(3,*). N.B.: if WORK is not large
*                          enough IOPTRES is automatically reset to 0!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILFMA is used
*                          as list type and IFTRAN(3,*) as list index
*                          NOTE that IFTRAN(3,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, FILFMA is used
*                          as list type and IFTRAN(3,*) as list index
*                          NOTE that IFTRAN(3,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILFMA and the indeces from IFDOTS
*                          the result of the dot products is returned
*                          in the FCONS array
*
*     Written by Christof Haettig, November 1998, based on CC_FMAT.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include "priunit_dec.h"
#  include "ccsdsym_dec.h"
#  include "maxorb_dec.h"
#  include "mxorb_dec.h"
#  include "mxcent_dec.h"
#  include "ccorb_dec.h"
#  include "cbieri_dec.h"
#  include "distcl_dec.h"
#  include "iratdef_dec.h"
#  include "eribuf_dec.h"
#  include "ccisao_dec.h"
#  include "aovec_dec.h"
#  include "blocks_dec.h"
#else
#  include <implicit.h>
      INTEGER MAXSHL
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "mxorb.h"
#include "mxcent.h"
#include "ccsdio.h"
#include "ccorb.h"
#include "cciccset.h"
#include "cbieri.h"
#include "distcl.h"
#include "iratdef.h"
#include "eribuf.h"
#include "ccisao.h"
#include "ccfield.h"
#include "aovec.h"
#include "blocks.h"

* local parameters:
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CC_FMAT> ')

      LOGICAL LOCDBG
csonia      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL APPEND, NOAPPEND
      PARAMETER (APPEND = .TRUE., NOAPPEND = .FALSE.)

      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
      INTEGER ISYM0
      PARAMETER( ISYM0 = 1 ) ! symmetry of the reference state
      
      INTEGER LUBF,   LUBFD,  LUC,    LUD,    LUFK,   LURIM, LUBFI
      INTEGER LURHO,  LUFIM,  LUMGD,  LUXIM,  LUYIM,  LUBDZ0, LUO3
      INTEGER LUCBAR, LUDBAR, LUPQRR, LUPQR0, LUPQMO, LUFMAT, LUZDN
      INTEGER LUGIM

      CHARACTER*(8) BFFIL,  FNBFD,  CTFIL, DTFIL, FILFCK, FILRIM, FILO3
      CHARACTER*(8) RHOFIL,FIMFIL,FILMGD,FILXIM,FILYIM,FNBFI,FILZDN
      CHARACTER*(8) CBAFIL, DBAFIL, FNPQRR, FNPQR0, FNPQMO, FNBDZ0
      CHARACTER*(8) GIMFIL
      PARAMETER (BFFIL ='CCFM_BFI', FNBFD ='CCBFDENS',
     &           CBAFIL='CCFM_CBA', DBAFIL='CCFM_DBA',
     &           CTFIL ='CCFM_CIM', DTFIL ='CCFM_DIM', 
     &           FIMFIL='CCFM_FIM', RHOFIL='CCFM_RHO',
     &           FILRIM='CCFM_RIM', FILMGD='CCFM_MGD',
     &           FILXIM='CCFM_XIM', FILYIM='CCFM_YIM',
     &           FILFCK='CCFM_FCK', FNPQRR='CCFMPQRR',
     &           FNPQR0='CCFMPQR0', FNPQMO='CCFMPQMO',
     &           FNBDZ0='CCFMBDZ0', FNBFI ='CCFM_BZI',
     &           FILO3 ='CCFM_O3I', FILZDN='CCFM_ZDN',
     &           GIMFIL='CCFM_GIM'                     )


      CHARACTER*(*) LISTL, LISTR, FILFMA
      INTEGER IOPTRES
      INTEGER NFTRAN, MXVEC, LWORK
      INTEGER IFTRAN(3,NFTRAN)
      INTEGER IFDOTS(MXVEC,NFTRAN)


#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FCONS(MXVEC,NFTRAN)
      REAL DUM, XNORM, FF, DUMMY
      REAL DTIME, CONVRT, TIMALL
      REAL TIMFCK, TIMBF, TIMC, TIMD, TIMIMR, TIMXYM, TIMIO
      REAL TIMPQ, TIMIML, TIMIM2, TIMPRE, TIMINT, TIMRDAO
      REAL TIMTRBT, TIMIM0, TIMTRN
      REAL ZERO, ONE, TWO, HALF, FACT
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FCONS(MXVEC,NFTRAN)
      DOUBLE PRECISION XNORM, FF
      DOUBLE PRECISION DTIME, CONVRT, TIMALL
      DOUBLE PRECISION TIMFCK, TIMBF, TIMC, TIMD, TIMIMR, TIMXYM, TIMIO
      DOUBLE PRECISION TIMPQ, TIMIML, TIMIM2, TIMPRE, TIMINT, TIMRDAO
      DOUBLE PRECISION TIMTRBT, TIMIM0, TIMTRN
      DOUBLE PRECISION ZERO, ONE, TWO, HALF, FACT
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)

      CHARACTER*(10) MODEL, MODELW, CDUMMY
      LOGICAL LLSAME
      INTEGER INDEXA(MXCORB)
      INTEGER INTMEDR(2,MAXSIM), NINTR, IINTR
      INTEGER INTMEDL(4,MAXSIM), NINTL, IINTL
      INTEGER INTMED2(4,MAXSIM), NINT2, IINT2
      INTEGER IRHGH(0:MAXSIM), I2HGH(0:MAXSIM), IOFFCD(0:MAXSIM+1)
      INTEGER IADRBFD(MAXORB,NFTRAN), IADRZ0(MAXORB,NFTRAN)
      INTEGER IADRPQ0(MAXORB,NFTRAN), IADRPQR(MAXORB,NFTRAN)
      INTEGER IADRPQMO(MAXORB,NFTRAN),IADRBFI(MAXORB,NFTRAN)
      INTEGER KRHO1(NFTRAN),KFIM(NFTRAN),KMGDL(NFTRAN),KRIM(NFTRAN)
      INTEGER KZDEN(NFTRAN),KGIM(NFTRAN)
      INTEGER KLAMPA(NFTRAN),KLAMHA(NFTRAN),KRHO2(NFTRAN),KCTR2(NFTRAN)
      INTEGER KLAMP(NFTRAN),KLAMH(NFTRAN),KDENS(NFTRAN),KFOCK(NFTRAN)
      INTEGER LENX,LENY,LENMGD,LENR,LENFIM,LENRHO,LENFK,LENBF,LEN
      INTEGER MT2BGD,MDISAO,MDSRHF,MSCRATCH,MEMAVAIL,NWORK,NSECMAR
      INTEGER NNWORK, IDUMMY, ICDEL2, IBATCH, KEND2, JEND2, LWRK2
      INTEGER IOPTW, ITRAN, IOPT, ICORE, IF, IOPTRSP, ISTARTBFI, IVEC
      INTEGER ISYM, ICOUNT, ISYMAK, ISYBET, IBSRHF(8,8), NBSRHF(8)
      INTEGER ISYML, ISYMR, ISYRES, IDLSTL, IDLSTR, KEND1, JEND1
      INTEGER KFOCK0, KDENS0, KT1AMP0, KDSRHFA, ISYMXA, KFCKC0, KDNSC0
      INTEGER KLAMP0, KLAMH0, KTHETA1, KTHETA2, KZETA1, KT1AMPA
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2, NBATCH
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2, KFREE
      INTEGER LFREE, KEND, LWRK, KENDSV, LWRKSV, KEND0, LWRK0, LWRK1
      INTEGER NTOSYM, NTOT, KRECNR, ISYMD1, ILLL, NUMDIS, IDEL2, IDEL
      INTEGER ISYDEL, KXINT, KEND3, K3OINT, KBFRHF, KDCRHF, KDSRHF
      INTEGER KFINT, KRIMA, KFOCKA, KTHETA0, KYINT, KYBARA, LENZDN
      INTEGER LWRK3, IADRTH, KLAMDPA, KLAMDHA, LENALL, KT2AMPA
      INTEGER ISYCTR, ISYAMP, KLIAJB, KA2IM, KXIDJL, KXJLID
      INTEGER KXIAJB, NDBLE, LENGIM, KGZETA, KFOCK0AO, KXBARA

* external functions:
      INTEGER ICCSET1
      INTEGER ICCSET2
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL DDOT, SECOND
#else 
      DOUBLE PRECISION DDOT, SECOND
#endif 
  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_FMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        WRITE (LUPRI,*) MSGDBG, 'LISTL : ',LISTL
        WRITE (LUPRI,*) MSGDBG, 'LISTR : ',LISTR
        WRITE (LUPRI,*) MSGDBG, 'FILFMA: ',FILFMA
        WRITE (LUPRI,*) MSGDBG, 'NFTRAN: ',NFTRAN
        WRITE (LUPRI,*) MSGDBG, 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
        WRITE (LUPRI,'(/1x,2a)') 'CC_FMAT called for a Coupled Cluster'
     &          ,' method not implemented in CC_FMAT...'
        CALL QUIT('Unknown CC method in CC_FMAT.')
      END IF

      IF (.NOT. DUMPCD) THEN
        WRITE (LUPRI,*) 'DUMPCD = ',DUMPCD
        WRITE (LUPRI,*) 'CC_FMAT requires DUMPCD=.TRUE.'
        CALL QUIT('DUMPCD=.FALSE. , CC_FMAT requires DUMPCD=.TRUE.')
      END IF

      IF (NFTRAN .GT. MAXSIM) THEN
        WRITE (LUPRI,*) 
     *            'Error in CC_FMAT: NFTRAN is larger than MAXSIM.'
        CALL QUIT('Error in CC_FMAT: NFTRAN is larger than MAXSIM.')
      END IF

      IF (IPRINT.GT.0) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)') '+','+'

         WRITE (LUPRI,'(1x,A52)')
     &         '|        F MATRIX TRANSFORMATION SECTION           |'

         IF (IOPTRES.EQ.3) THEN
            WRITE (LUPRI,'(1X,A52)') 
     &         '|          (result is written to file)             |'
         ELSE IF (IOPTRES.EQ.4) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|     (result is added to a vector on file)        |'
         ELSE IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A52)') 
     &         '|    (result used to calculate dot products)       |'
         END IF
        
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

      END IF

* initialize timings:
      TIMALL  = SECOND()
      TIMIM0  = ZERO
      TIMIML  = ZERO
      TIMIMR  = ZERO
      TIMIM2  = ZERO
      TIMPRE  = ZERO
      TIMXYM  = ZERO
      TIMPQ   = ZERO
      TIMFCK  = ZERO
      TIMBF   = ZERO
      TIMC    = ZERO
      TIMD    = ZERO
      TIMINT  = ZERO
      TIMRDAO = ZERO
      TIMTRBT = ZERO
      TIMIO   = ZERO

* set option and model to write vectors to file:
      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IOPTW  = 3
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IOPTW  = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_FMAT.')
      END IF

* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
         LUFMAT = 0
         CALL WOPEN2(LUFMAT, FILFMA, 64, 0)
      ELSE IF (IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4) THEN
         CONTINUE
      ELSE IF (IOPTRES .EQ. 5) THEN
         IF (MXVEC*NFTRAN.NE.0) CALL DZERO(FCONS,MXVEC*NFTRAN)
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_FMAT.')
      END IF
 
* precalculate symmetry array for BSRHF:
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO

*=====================================================================*
* build nonredundant arrays of response vectors and pairs of them
* for which intermediates have to be calculated
*=====================================================================*
      DTIME = SECOND()

* array for intermediates which depend on the left vectors:
      NINTL = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET2(INTMEDL,LISTL,IFTRAN(1,ITRAN),
     &                    'R0 ',0,              NINTL,MAXSIM,APPEND)
      END DO 

* array for intermediates which depend on the right vectors:
      NINTR = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET1(INTMEDR,LISTR,IFTRAN(2,ITRAN),NINTR,MAXSIM,APPEND)
      END DO 

* array for intermediates that depend on left and right vectors: 
      NINT2 = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET2(INTMED2,LISTL,IFTRAN(1,ITRAN), 
     &                    LISTR,IFTRAN(2,ITRAN),NINT2,MAXSIM,APPEND)
      END DO 


      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A)') 
     *         'List of response vector for left intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMEDL(I,J),I=1,2),J=1,NINTL)
        WRITE (LUPRI,'(/A)') 
     *         'List of response vector for right intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMEDR(I,J),I=1,2),J=1,NINTR)
        WRITE (LUPRI,'(/A)') 
     *         'List of vector pairs for 2. order intermediates:'
        WRITE (LUPRI,'((/5X,4I5))') ((INTMED2(I,J),I=1,4),J=1,NINT2)
      END IF
 
      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* estimate scratch space requirements
*---------------------------------------------------------------------*
      DTIME = SECOND()

      MT2BGD = 0
      MDISAO = 0
      MDSRHF = 0
      DO ISYM = 1, NSYM
        MT2BGD = MAX(MT2BGD,NT2BGD(ISYM))
        MDISAO = MAX(MDISAO,NDISAO(ISYM))
        MDSRHF = MAX(MDSRHF,NDSRHF(ISYM))
      END DO

*     5 x a NT2BGD type intermediate  
*     + integral arrays + some reserve

      MSCRATCH = 5*MT2BGD + MDISAO + MDSRHF + 10*N2BASX  
      IF (CC2)               MSCRATCH = MSCRATCH+MDISAO+2*MDSRHF
      IF (.NOT.(CCS.OR.CC2)) MSCRATCH = MSCRATCH+MDISAO+4*MDSRHF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 
     *         'scratch space estimate MSCRATCH:',MSCRATCH
        CALL FLSHFO(LUPRI)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* estimate memory for 'in core' version and batched versions:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      MEMAVAIL = LWORK - MSCRATCH

      NWORK  = 0
      NBATCH = 1
      IF (CCS) THEN
        NSECMAR = 10 * N2BASX
      ELSE IF (CC2 .OR. CCSD) THEN
        NSECMAR = 10 * MT2BGD
      ELSE
          CALL QUIT('Unknown CC model in CC_FMAT.')
      END IF

      IRHGH(0) = 0
      I2HGH(0) = 0

* intermediates that dependent on right response vectors:
* (see routine CCBPRE1 for details)
      DO IINTR = 1, NINTR
        IDLSTR = INTMEDR(1,IINTR)
        ISYM   = ILSTSYM(LISTR,IDLSTR)
 
        NNWORK = 2*NGLMDT(ISYM) + 2*N2BST(ISYM)
        IF (CCSD) NNWORK = 2*NGLMDT(ISYM)+NT2AOIJ(ISYM)+NEMAT1(ISYM)

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          IRHGH(NBATCH) = IINTR - 1
          I2HGH(NBATCH) = 0
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_FMAT. (01)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK, ' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_FMAT. (01)')
        END IF
      END DO

* intermediates that dependent on left and right vectors:
* (see routine CCFPRE2 for details)
      DO IINT2 = 1, NINT2
        IDLSTL = INTMED2(1,IINT2)
        IDLSTR = INTMED2(3,IINT2)
        ISYML  = ILSTSYM(LISTL,IDLSTL)
        ISYMR  = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYML,ISYMR)
 
        IF (CCS) THEN
          NNWORK = 2*N2BST(ISYRES)
        ELSE IF (CC2) THEN
          NNWORK = 2*N2BST(ISYRES) + 2*NGLMDT(ISYMR) +
     &               NT2SQ(ISYML)  + NT1AM(ISYRES)
        ELSE IF (CCSD) THEN
          NNWORK = 2*NGLMDT(ISYMR) + 
     &               NT1AO(ISYRES) + NT1AM(ISYRES) + NT2AOIJ(ISYRES)
        ELSE
          CALL QUIT('Unknown CC model in CC_FMAT.')
        END IF

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          IRHGH(NBATCH) = NINTR
          I2HGH(NBATCH) = IINT2 - 1
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_FMAT. (02)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK,' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_FMAT. (02)')
        END IF
      END DO
 
      IRHGH(NBATCH) = NINTR
      I2HGH(NBATCH) = NINT2

      IF   (LOCDBG .AND. (NBATCH.EQ.1)) THEN
        WRITE (LUPRI,*) MSGDBG, 
     *         'one batch only... will be done in core.'
        WRITE (LUPRI,*) MSGDBG, 
     *         'memory for intermediates: ', NWORK
        WRITE (LUPRI,*) MSGDBG, 
     *         'remaining scratch space: ', LWORK-NWORK
        CALL FLSHFO(LUPRI)
      ELSE IF (LOCDBG .AND. (NBATCH.GT.1)) THEN
        WRITE (LUPRI,*) MSGDBG, 
     *         'more than one batch... choose I/O algorithm.'
        WRITE (LUPRI,*) MSGDBG, 
     *         'max. memory for intermediates: ',MEMAVAIL
        WRITE (LUPRI,*) MSGDBG, 
     *         'number of batches: ',NBATCH
        CALL FLSHFO(LUPRI)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* read zeroth-order singles amplitudes, allocate space for Fock matrix,
* and prepare zeroth-order lambda matrices and density: 
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KFOCK0AO = 1
      KFOCK0   = KFOCK0AO + N2BAST

      KDENS0   = KFOCK0   + N2BAST
      KT1AMP0  = KDENS0   + N2BAST
      KLAMP0   = KT1AMP0  + NT1AMX
      KLAMH0   = KLAMP0   + NLAMDT
      KEND0    = KLAMH0   + NLAMDT

      IF (FROIMP.OR.FROEXP) THEN
        KFCKC0   = KEND0
        KDNSC0   = KFCKC0 + N2BAST
        KEND0    = KDNSC0 + N2BAST
      END IF

      LWRK0    = LWORK - KEND0

      IF (LWRK0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_FMAT. (0)')
      END IF

* read zeroth order amplitudes:
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get unperturbed Lambda matrices:
      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND0),LWRK0)

* calculate the density matrix:
      ICORE = 1 ! include core contribution
      CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDENS0),
     &               ISYM0,ICORE, WORK(KEND0),LWRK0)

* calculate pure core contribution to the density matrix,
* and initialize core contribution to Fock matrix with zeros
      IF (FROIMP.OR.FROEXP) THEN
        ICORE = 0 ! exclude core contribution
        CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDNSC0),
     &                 ISYM0,ICORE, WORK(KEND0),LWRK0)
        CALL DSCAL(N2BAST,-ONE,WORK(KDNSC0),1)
        CALL DAXPY(N2BAST,+ONE,WORK(KDENS0),1,WORK(KDNSC0),1)
        CALL DZERO(WORK(KFCKC0),N2BAST)
      END IF

* initialize Fock matrix with the one-electron integrals:
      CALL CCRHS_ONEAO(WORK(KFOCK0),WORK(KEND0),LWRK0)
      DO IF= 1, NFIELD
        FF = EFIELD(IF)
        CALL CC_ONEP(WORK(KFOCK0),WORK(KEND0),LWRK0,FF,1,LFIELD(IF) )
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'norm of T1AMP0:',
     &        DDOT(NT1AMX,WORK(KT1AMP0),1,WORK(KT1AMP0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of XLAMP0:',
     &        DDOT(NLAMDT,WORK(KLAMP0),1,WORK(KLAMP0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of DENS0:',
     &        DDOT(N2BAST,WORK(KDENS0),1,WORK(KDENS0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of FOCK0:',
     &        DDOT(N2BAST,WORK(KFOCK0),1,WORK(KFOCK0),1)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend only on the left 
* vectors: X, Y, M, Chi, Yps, P, Q, and the BFZeta density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IOPTRSP = 1
      CALL CCFPREINT(INTMEDL,NINTL,WORK(KLAMP0),WORK(KLAMH0),
     &             CDUMMY,LENX,CDUMMY,LENY,CDUMMY,LENMGD,CDUMMY,LENZDN,
     &             FNPQR0,IADRPQ0,FNPQMO,IADRPQMO,FNBDZ0,IADRZ0,
     &             TIMXYM,TIMBF,TIMIO,TIMPQ,IOPTRSP,WORK(KEND0),LWRK0)

      TIMIML = TIMIML + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend only on the right
* vectors: CBAR, DBAR, and the BF density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL CCBPRE1INT(INTMEDR,NINTR,IOFFCD,IADRBFD,
     &                CBAFIL,DBAFIL,FNBFD,
     &                WORK(KLAMP0),WORK(KLAMH0),
     &                TIMIO,TIMC,TIMD,TIMBF,WORK(KEND0),LWRK0)

      TIMIMR = TIMIMR + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend on left and right
* vectors: response X, Y, M, Chi, Yps, P, Q, and BFZeta density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IOPTRSP = 2
      CALL CCFPREINT(INTMED2,NINT2,WORK(KLAMP0),WORK(KLAMH0),
     &             FILXIM,LENX,FILYIM,LENY,FILMGD,LENMGD,FILZDN,LENZDN,
     &             FNPQRR,IADRPQR,CDUMMY,IDUMMY,CDUMMY,IDUMMY,
     &             TIMXYM,TIMBF,TIMIO,TIMPQ,IOPTRSP,WORK(KEND0),LWRK0)

      TIMIM2 = TIMIM2 + SECOND() - DTIME

*---------------------------------------------------------------------*
* open files for local intermediates generated in the loop over the
* AO integrals and which need to be initialized here:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL CCFOPEN(LUBF,  LUFK,   LURIM,  LUFIM,  LUGIM,  LURHO,
     &             BFFIL, FILFCK, FILRIM, FIMFIL, GIMFIL, RHOFIL,
     &             LENBF, LENFK,  LENR,   LENFIM, LENGIM, LENRHO,
     &             NINTR, NINT2,  WORK(KEND0),    LWRK0           )

* open some other files needed in the loop over integrals, but need
* not to be initialized:
      IF (CCS.OR.CC2) THEN
         LUZDN = 0
         CALL WOPEN2(LUZDN, FILZDN, 64,0)
      ELSE IF (.NOT.(CCS.OR.CC2)) THEN
         LUC     = 0
         LUD     = 0
         LUMGD   = 0
         LUBFI   = 0
         LUBFD   = 0
         LUPQR0  = 0
         LUPQRR  = 0
         LUBDZ0  = 0
         LUO3    = 0
         CALL WOPEN2(LUC,   CTFIL, 64,0)
         CALL WOPEN2(LUD,   DTFIL, 64,0)
         CALL WOPEN2(LUMGD, FILMGD,64,0)
         CALL WOPEN2(LUBFI, FNBFI, 64,0)
         CALL WOPEN2(LUBFD, FNBFD, 64,0)
         CALL WOPEN2(LUPQR0,FNPQR0,64,0)
         CALL WOPEN2(LUPQRR,FNPQRR,64,0)
         CALL WOPEN2(LUBDZ0,FNBDZ0,64,0)
         CALL WOPEN2(LUO3 , FILO3, 64,0)
      END IF

* initialize offsets C and D intermediates:
      ICDEL2    = 0
      
* initialize start address for BFI intermediates:
      ISTARTBFI = 1

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* if all vectors and intermediates fit into the memory, read all
* response vectors before the loop over AO integral shells:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (NBATCH .EQ. 1) THEN

            CALL CCBPRE1(INTMEDR, 1, NINTR,
     &                   KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                   LUBF,BFFIL,LENBF,LUFK,FILFCK,LENFK,
     &                   LURIM,FILRIM,LENR,
     &                   WORK(KLAMP0), WORK(KLAMH0),
     &                   WORK, LWORK, KEND0, JEND1 )
            KEND1 = JEND1

            CALL CCFPRE2(INTMED2,1,NINT2,
     &                   KFIM, LUFIM,FIMFIL,LENFIM,
     &                   KGIM, LUGIM,GIMFIL,LENGIM,
     &                   KRHO1,LURHO,RHOFIL,LENRHO,
     &                   KMGDL,LUMGD,FILMGD,LENMGD,
     &                   KZDEN,LUZDN,FILZDN,LENZDN,
     &                   KCTR2,KLAMPA,KLAMHA,
     &                   WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK, LWORK, KEND1, JEND1  )
            KEND1 = JEND1

            IF (LOCDBG) THEN
             WRITE (LUPRI,*) MSGDBG,
     *              'allocated work space for intermediates:'
             WRITE (LUPRI,*) MSGDBG,'KRHO2 :',(KRHO2(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KLAMP :',(KLAMP(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KLAMH :',(KLAMH(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KDENS :',(KDENS(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KFOCK :',(KFOCK(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KRIM  :',(KRIM(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KRHO1 :',(KRHO1(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KZDEN :',(KZDEN(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KFIM  :',(KFIM(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KGIM  :',(KGIM(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KMGDL :',(KMGDL(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KLAMPA:',(KLAMPA(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KLAMHA:',(KLAMHA(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KEND1:',KEND1
             CALL FLSHFO(LUPRI)
            END IF

      ELSE
        KEND1 = KEND0
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_FMAT. (1)')
      END IF
 
      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* initialize integral calculation
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KEND = KEND1
      LWRK = LWRK1

      IF (DIRECT) THEN
         NTOSYM = 1

         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND),LWRK,IPRERI)
         ELSE
           KCCFB1 = KEND
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK   = LWORK  - KEND
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND),LWRK,IPRERI)
 
           KEND = KFREE
           LWRK = LFREE
         END IF

         KENDSV = KEND
         LWRKSV = LWRK
      ELSE
         NTOSYM = NSYM
      END IF

      TIMINT = TIMINT + SECOND() - DTIME
*---------------------------------------------------------------------*
* start loop over AO integral shells:
*---------------------------------------------------------------------*
      DO ISYMD1 = 1, NTOSYM
       
        IF (DIRECT) THEN
          IF (HERDIR) THEN
             NTOT = MAXSHL
          ELSE
             NTOT = MXCALL        
          ENDIF
        ELSE
          NTOT = NBAS(ISYMD1)
        END IF
 
        DO ILLL = 1, NTOT

          DTIME = SECOND()
 
          IF (DIRECT) THEN
            KEND = KENDSV
            LWRK = LWRKSV
 
            IF (HERDIR) THEN
               CALL HERDI2(WORK(KEND),LWRK,INDEXA,ILLL,NUMDIS,
     &                     IPRINT)
            ELSE
               CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                     WORK(KODCL1),WORK(KODCL2),
     *                     WORK(KODBC1),WORK(KODBC2),
     *                     WORK(KRDBC1),WORK(KRDBC2),
     *                     WORK(KODPP1),WORK(KODPP2),
     *                     WORK(KRDPP1),WORK(KRDPP2),
     *                     WORK(KCCFB1),WORK(KINDXB),
     *                     WORK(KEND),LWRK,IPRERI)
            END IF

            KRECNR = KEND
            KEND   = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LWRK   = LWORK - KEND
  
            IF (LWRK .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_FMAT. (1a)')
            END IF
 
          ELSE
            NUMDIS = 1
          END IF
 
          TIMINT = TIMINT + SECOND() - DTIME
*---------------------------------------------------------------------*
*        if out of core: allocate memory and get response vectors:
*---------------------------------------------------------------------*
          DO IBATCH = 1, NBATCH
             KEND2 = KEND ! reset memory for each batch
     
             IF (LOCDBG) THEN
               WRITE (LUPRI,*) MSGDBG, IBATCH,
     &                                  '-th. batch out of ',NBATCH
               WRITE (LUPRI,*) MSGDBG, 'IR:',IRHGH(IBATCH-1)+1,' -- ',
     &                                  IRHGH(IBATCH)
               WRITE (LUPRI,*) MSGDBG, 'I2:',I2HGH(IBATCH-1)+1,' -- ',
     &                                  I2HGH(IBATCH)
             END IF
 
             IF (NBATCH.GT.1) THEN
 
               DTIME = SECOND()

               CALL CCBPRE1(INTMEDR,IRHGH(IBATCH-1)+1,IRHGH(IBATCH),
     &                      KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                      LUBF,BFFIL,LENBF,LUFK,FILFCK,LENFK,
     &                      LURIM,FILRIM,LENR,
     &                      WORK(KLAMP0), WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2 )
               KEND2 = JEND2
 
               CALL CCFPRE2(INTMED2,I2HGH(IBATCH-1)+1,I2HGH(IBATCH),
     &                      KFIM, LUFIM,FIMFIL,LENFIM,
     &                      KGIM, LUGIM,GIMFIL,LENGIM,
     &                      KRHO1,LURHO,RHOFIL,LENRHO,
     &                      KMGDL,LUMGD,FILMGD,LENMGD,
     &                      KZDEN,LUZDN,FILZDN,LENZDN,
     &                      KCTR2,KLAMPA,KLAMHA,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2 )
               KEND2 = JEND2

               TIMPRE = TIMPRE + SECOND() - DTIME
 
             END IF

             LWRK2 = LWORK - KEND2
             IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_FMAT. (2)')
             END IF
 
*---------------------------------------------------------------------*
*        loop over number of distributions on the disk:
*---------------------------------------------------------------------*
          DO IDEL2  = 1, NUMDIS
    
            IF (DIRECT) THEN
              IDEL   = INDEXA(IDEL2)
              ISYDEL = ISAO(IDEL)
            ELSE
              IDEL   = IBAS(ISYMD1) + ILLL
              ISYDEL = ISYMD1
            END IF
 
*           read AO integral distribution and calculate integrals with
*           one index transformed to occupied MO (particle):
            
            KXINT  = KEND2
            KEND3  = KXINT  + NDISAO(ISYDEL)

            IF (CC2) THEN
               KDSRHF  = KEND3
               KDSRHFA = KDSRHF  + NDSRHF(ISYDEL)
               KEND3   = KDSRHFA + MDSRHF
            ELSE IF (CCSD) THEN
               K3OINT  = KEND3
               KBFRHF  = K3OINT + NMAIJK(ISYDEL)
               KDCRHF  = KBFRHF + NBSRHF(ISYDEL)
               KDSRHF  = KDCRHF + NBSRHF(ISYDEL)
               KEND3   = KDSRHF + NDSRHF(ISYDEL)
            END IF

            LWRK3  = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_FMAT. (3)')
            END IF
 
            DTIME = SECOND()
            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     &                  WORK(KRECNR),DIRECT)
            TIMRDAO = TIMRDAO + SECOND() - DTIME

C           XNORM = DDOT(NDISAO(ISYDEL),WORK(KXINT),1,WORK(KXINT),1)
C           WRITE (LUPRI,*) 'IDEL,XNORM:',IDEL,XNORM


            IF (CCSD) THEN
C              -----------------------------------------------------
C              some integral transformations and presortings for
C              BF, C, and D intermediate and 21G term:
C                   DSRHF  --  standard (**|kdel) integrals
C                   BFRHF  --  (**|kdel) presorted for B & D interm.
C                   DCRHF  --  (**|kdel) presorted for C & D interm.
C                   3OINT  --  (ij|kdel) integrals for 21G term
C              -----------------------------------------------------
               DTIME = SECOND()

               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMP0),ISYM0,
     &                     WORK(KEND3),LWRK3,ISYDEL)

               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMH0),
     &                       ISYM0,WORK(KLAMP0),WORK(KEND3),
     &                       LWRK3,IDEL,ISYDEL,LUO3,FILO3)

               CALL CC_BFBSORT(WORK(KDSRHF),WORK(KBFRHF),ISYDEL)
               
               CALL CCB_CDSORT(WORK(KXINT),ISYDEL,WORK(KDCRHF),
     &                         WORK(KLAMP0),ISYM0,WORK(KEND3),LWRK3)

               TIMTRBT = TIMTRBT + SECOND() - DTIME
               TIMIM0  = TIMIM0  + SECOND() - DTIME

               KEND3  = KDSRHF
               LWRK3  = LWORK  - KEND3

               KDSRHF = KDUM

            END IF


C           ----------------------------------------------------------
C           calculate zeroth-order Fock matrix (Fhat), as well as the
C           pure core contribution to the zeroth-order Fock matrix
C           ----------------------------------------------------------
            IF (IBATCH .EQ. 1) THEN
              DTIME = SECOND()

              CALL CC_AOFOCK(WORK(KXINT), WORK(KDENS0),
     *                       WORK(KFOCK0),WORK(KEND3),
     *                       LWRK3,IDEL,ISYDEL,ISYM0)

              IF (FROIMP.OR.FROEXP) THEN
                 CALL CC_AOFOCK(WORK(KXINT), WORK(KDNSC0),
     *                          WORK(KFCKC0),WORK(KEND3),
     *                          LWRK3,IDEL,ISYDEL,ISYM0)
              END IF

              TIMFCK = TIMFCK + SECOND() - DTIME
              TIMIM0 = TIMIM0 + SECOND() - DTIME
            END IF
 

C           ----------------------------------------------------------
C           calculate intermediates that depend only on right vectors:
C           ----------------------------------------------------------
            DO IINTR = IRHGH(IBATCH-1)+1, IRHGH(IBATCH)
              IDLSTR = INTMEDR(1,IINTR)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)

*             calculate addresses for C & D intermediates:
              IT2DLR(IDEL,IINTR) = ICDEL2
              ICDEL2 = ICDEL2 + NT2BCD(MULD2H(ISYDEL,ISYMR))

              DTIME = SECOND()
              CALL CCBINT1(WORK(KXINT), WORK(KBFRHF), WORK(KDCRHF),
     &                     IDEL, ISYDEL,  WORK(KRHO2(IINTR)),  
     &                     WORK(KLAMP0),  WORK(KLAMH0), 
     &                     WORK(KLAMP(IINTR)),  WORK(KLAMH(IINTR)),
     &                     ISYMR, IINTR,
     &                     WORK(KDENS(IINTR)),  WORK(KFOCK(IINTR)),  
     &                     WORK(KRIM(IINTR)),
     &                     LUC, CTFIL, LUD, DTFIL, 
     &                     LUBFD, FNBFD, IADRBFD(1,IINTR), 
     &                     WORK(KEND3),   LWRK3,          
     &                     TIMFCK, TIMBF, TIMC, TIMD  )
              TIMIMR = TIMIMR + SECOND() - DTIME

            END DO
 
 
C           ----------------------------------------------------------
C           calculate intermediates that depend on both the left and
C           the right vectors:
C           ----------------------------------------------------------
            DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
              IDLSTL = INTMED2(1,IINT2)
              IDLSTR = INTMED2(3,IINT2)
              ISYML  = ILSTSYM(LISTL,IDLSTL)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)
              ISYRES = MULD2H(ISYML,ISYMR)

              IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                         'R0 ',0,NINTL,MAXSIM,NOAPPEND)

              DTIME = SECOND()
              CALL CCFINT2(IDEL,ISYDEL,ISYML,WORK(KXINT),
     &                   WORK(KBFRHF),WORK(K3OINT),
     &                   WORK(KMGDL(IINT2)),WORK(KZDEN(IINT2)),
     &                   WORK(KFIM(IINT2)), WORK(KGIM(IINT2)),
     &                   WORK(KRHO1(IINT2)),
     &                   LUBFI,FNBFI,IADRBFI(1,IINT2),ISTARTBFI,
     &                   LUPQR0,FNPQR0,IADRPQ0(1,IINTL),
     &                   LUPQRR,FNPQRR,IADRPQR(1,IINT2),
     &                   LUBDZ0,FNBDZ0,IADRZ0(1,IINTL),
     &                   WORK(KLAMPA(IINT2)),WORK(KLAMHA(IINT2)),ISYMR,
     &                   WORK(KLAMP0), WORK(KLAMH0), 
     &                   WORK(KEND3),  LWRK3                          )
              TIMIM2 = TIMIM2 + SECOND() - DTIME

            END DO
 
C           ----------------------------------------------------------
C           for CC2 calculate the 21CD term contributions.
C           ----------------------------------------------------------
            IF (CC2) THEN
C            ........................................................
C            DSRHF  --  (**|kdel) integrals, transformed with XLAMH0:
C            ........................................................
             DTIME = SECOND()
             CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMH0),ISYM0,
     &                   WORK(KEND3),LWRK3,ISYDEL)
             TIMTRBT = TIMTRBT + SECOND() - DTIME
             TIMIM0  = TIMIM0  + SECOND() - DTIME

             DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
              IDLSTL = INTMED2(1,IINT2)
              IDLSTR = INTMED2(3,IINT2)
              ISYML  = ILSTSYM(LISTL,IDLSTL)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)
              ISYRES = MULD2H(ISYML,ISYMR)

              IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                         'R0 ',0,NINTL,MAXSIM,NOAPPEND)

C             .............................................
C             calculate the first term depending on DSRHF0:
C             .............................................
              LLSAME = .FALSE.
              FACT   = ONE
              DTIME = SECOND()
              CALL CCG_21CD( WORK(KRHO1(IINT2)),ISYRES,
     &                   WORK(KCTR2(IINT2)),ISYML,
     &                   WORK(KLAMH0), WORK(KLAMP0), 
     &                   WORK(KLAMHA(IINT2)),WORK(KLAMPA(IINT2)),ISYMR,
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KDSRHF), ISYDEL, IDEL, ISYDEL,
     &                   LLSAME, FACT, WORK(KEND3),  LWRK3 )

C             ..............................................
C             calculate the second term depending on DSRHFA:
C             ..............................................
              DTIME = SECOND()
              CALL CCTRBT(WORK(KXINT),WORK(KDSRHFA),
     &                    WORK(KLAMHA(IINT2)),ISYMR,
     &                    WORK(KEND3),LWRK3,ISYDEL)
              TIMTRBT = TIMTRBT + SECOND() - DTIME


              LLSAME = .TRUE.
              FACT   = HALF
              ISYMXA = MULD2H(ISYDEL,ISYMR)
              DTIME = SECOND()
              CALL CCG_21CD( WORK(KRHO1(IINT2)),ISYRES,
     &                   WORK(KCTR2(IINT2)),ISYML,
     &                   WORK(KLAMH0), WORK(KLAMP0), 
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KDSRHFA), ISYMXA, IDEL, ISYDEL,
     &                   LLSAME, FACT, WORK(KEND3),  LWRK3 )
              TIMIM2 = TIMIM2 + SECOND() - DTIME

             END DO
            END IF

          END DO ! IDEL2
*---------------------------------------------------------------------*
*         end of the loop over integral distributions:
*         if batched I/O algorithm used, save result on disc:
*---------------------------------------------------------------------*
          IF (NBATCH.GT.1) THEN
            DTIME = SECOND()
            CALL CCFSAVE(IBATCH,  IRHGH, I2HGH, INTMEDR, INTMED2,
     &                   KRHO2,   LUBF,  BFFIL, LENBF,
     &                   KFOCK,   LUFK,  FILFCK,LENFK,
     &                   KRIM,    LURIM, FILRIM,LENR,
     &                   KFIM,    LUFIM, FIMFIL,LENFIM,
     &                   KGIM,    LUGIM, GIMFIL,LENGIM,
     &                   KRHO1,   LURHO, RHOFIL,LENRHO,
     &                   NINTR,   NINT2, WORK,  LWORK )
            TIMIO = TIMIO + SECOND() - DTIME
          END IF


        END DO ! IBATCH
       END DO ! ILLL
      END DO ! ISYMD1
*=====================================================================*
* End of Loop over AO-integrals
*=====================================================================*

*---------------------------------------------------------------------*
* if in-core algorithm used, save results now on disc:
*---------------------------------------------------------------------*
      IF (NBATCH.EQ.1) THEN
        DTIME = SECOND()
        CALL CCFSAVE(1,       IRHGH, I2HGH, INTMEDR, INTMED2,
     &               KRHO2,   LUBF,  BFFIL, LENBF,
     &               KFOCK,   LUFK,  FILFCK,LENFK,
     &               KRIM,    LURIM, FILRIM,LENR,
     &               KFIM,    LUFIM, FIMFIL,LENFIM,
     &               KGIM,    LUGIM, GIMFIL,LENGIM,
     &               KRHO1,   LURHO, RHOFIL,LENRHO,
     &               NINTR,   NINT2, WORK,  LWORK )
        TIMIO = TIMIO + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* close & delete some files which are no longer needed:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (CCS.OR.CC2) THEN
         CALL WCLOSE2(LUZDN, FILZDN,'DELETE')
      ELSE IF (.NOT.(CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUMGD, FILMGD,'DELETE')
         CALL WCLOSE2(LUPQR0,FNPQR0,'DELETE')
         CALL WCLOSE2(LUPQRR,FNPQRR,'DELETE')
         CALL WCLOSE2(LUBDZ0,FNBDZ0,'DELETE')
         CALL WCLOSE2(LUO3 , FILO3, 'DELETE')
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME

*---------------------------------------------------------------------*
* recover workspace:
*---------------------------------------------------------------------*
      KEND1 = KEND0
      LWRK1 = LWRK0


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'Loop over AO-integrals completed ',
     &           ' & AO intermediates saved on file.'
        WRITE (LUPRI,*) MSGDBG,
     &           'recover work space: KEND1,LWRK1=',KEND1,LWRK1
        WRITE (LUPRI,*) MSGDBG,'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* open files for CBAR, DBAR, X & Y intermediates:
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN
         LUCBAR = 0
         LUDBAR = 0
         LUPQMO = 0
         CALL WOPEN2(LUCBAR,CBAFIL, 64,0)
         CALL WOPEN2(LUDBAR,DBAFIL, 64,0)
         CALL WOPEN2(LUPQMO,FNPQMO, 64,0)
      END IF

      IF (.NOT.CCS) THEN
         LUXIM = 0
         LUYIM = 0
         CALL WOPEN2(LUXIM, FILXIM, 64,0)
         CALL WOPEN2(LUYIM, FILYIM, 64,0)
      END IF
*---------------------------------------------------------------------*
* save a copy of the zeroth-order AO Fock matrix:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL DCOPY(N2BAST,WORK(KFOCK0),1,WORK(KFOCK0AO),1)

      TIMFCK = TIMFCK + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

*---------------------------------------------------------------------*
* do some printout:
*---------------------------------------------------------------------*
      IF (IPRINT.GT.0) THEN

         WRITE (LUPRI,'(1X,A,F10.2,A)')
     &    '| time for zero order intermediat.:',TIMIM0 ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINTL,' sets of left interm.:',TIMIML ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINTR,' sets of right inter.:',TIMIMR ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINT2,' sets of 2. ord. int.:',TIMIM2 ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A)')
     &    '| intermediates calculated in ',NBATCH,' batches          |'

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

         IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A)')
     &         '| L vector | R vector |  # products  |             |'
            WRITE (LUPRI,'(1X,3(A,A3),A)') '|  ',LISTL(1:3), ' No. |  '
     &       ,LISTR(1:3),' No. |   with ', FILFMA(1:3),
     &       '   |  time/secs  |'
         ELSE 
            WRITE (LUPRI,'(1X,A)') 
     &         '| L vector | R vector |    result    |             |'
            WRITE (LUPRI,'(1X,A2,A3,2(A,A3),A)')  '|  ',LISTL(1:3), 
     &        '  No. |  ', LISTR(1:3),' No. |   ', FILFMA(1:3),
     &        '  No.   |  time/secs  |'
         END IF

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
      END IF

*=====================================================================*
* calculate F matrix transformations:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NFTRAN

        IDLSTL = IFTRAN(1,ITRAN)
        IDLSTR = IFTRAN(2,ITRAN)

        ISYCTR = ILSTSYM(LISTL,IDLSTL)
        ISYAMP = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYAMP,ISYCTR)

        KTHETA1 = KEND1
        KEND2   = KTHETA1 + NT1AM(ISYRES)
        KTHETA2 = KDUM

        CALL DZERO(WORK(KTHETA1),NT1AM(ISYRES))

        IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                          'R0 ',0,      NINTL,MAXSIM,NOAPPEND)
        IINTR  = ICCSET1(INTMEDR,LISTR,IDLSTR,NINTR,MAXSIM,NOAPPEND)
        IINT2  = ICCSET2(INTMED2,LISTL,IDLSTL,
     &                           LISTR,IDLSTR,NINT2,MAXSIM,NOAPPEND)

        TIMTRN = SECOND()

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,
     *            'F matrix transformation for ITRAN,',ITRAN
           WRITE (LUPRI,*) MSGDBG,
     *            'IADRTH:',IADRTH
           WRITE (LUPRI,*) MSGDBG,
     *            'LISTL,IDLSTL:',LISTL,IDLSTL
           WRITE (LUPRI,*) MSGDBG,
     *            'LISTR,IDLSTR:',LISTR,IDLSTR
           WRITE (LUPRI,*) MSGDBG,
     *            'ISYCTR,ISYAMP,ISYRES:',ISYCTR,ISYAMP,ISYRES
           WRITE (LUPRI,*) MSGDBG,
     *            'IINTL,IINTR,IINT2:',IINTL,IINTR,IINT2
        END IF

*---------------------------------------------------------------------*
* read the single excitation part of the response vectors and 
* lagrangian multipliers and calculate the response Lambda matrices:
*---------------------------------------------------------------------*
        DTIME = SECOND()

        KZETA1  = KEND2
        KT1AMPA = KZETA1  + NT1AM(ISYCTR)
        KLAMDPA = KT1AMPA + NT1AM(ISYAMP)
        KLAMDHA = KLAMDPA + NGLMDT(ISYAMP)
        KEND2   = KLAMDHA + NGLMDT(ISYAMP)
        LWRK2   = LWORK - KEND2

        IF (LWRK2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_FMAT. (8)')
        END IF
      

        IOPT = 1

        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                WORK(KZETA1),WORK(KDUM))

        CALL CC_RDRSP(LISTR,IDLSTR,ISYAMP,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMDPA), WORK(KLAMH0),
     &                   WORK(KLAMDHA),WORK(KT1AMPA),ISYAMP)

        TIMPRE = TIMPRE + SECOND() - DTIME

*---------------------------------------------------------------------*
*       do some memory allocation and read (ia|jb) integrals:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
           KXBARA = KEND2
           KYBARA = KXBARA + NMATIJ(ISYAMP)
           KA2IM  = KYBARA + NMATAB(ISYAMP)
           KEND2  = KA2IM  + NT1AM(ISYRES)
  
           KLIAJB = KEND2
           KXIAJB = KLIAJB + NT2SQ(ISYM0)
           KEND3  = KXIAJB + NT2AM(ISYM0)
        ELSE IF (LISTL(1:2).EQ.'L0') THEN
           KXIAJB = KEND2
           KEND3  = KXIAJB + NT2AM(ISYM0)
        ELSE
           KEND3  = KEND2
        END IF

        LWRK3  = LWORK - KEND3

        IF (LWRK3 .LE. 0) THEN
           CALL QUIT('Insufficient work space in CC_FMAT.')
        END IF

C       -----------------------
C       read (ia|jb) integrals:
C       -----------------------
        IF ((.NOT.CCS) .OR. LISTL(1:2).EQ.'L0') THEN
           CALL CCG_RDIAJB(WORK(KXIAJB),NT2AM(ISYM0))
        END IF

C       --------------------------------------------------------
C       for CCSD resort (ia|jb) integrals to full square matrix:
C       --------------------------------------------------------
        IF (.NOT.CCS) THEN
           CALL CC_T2SQ(WORK(KXIAJB),WORK(KLIAJB),ISYM0)
        END IF

*---------------------------------------------------------------------*
*       calculate the projection on <HF|:
*---------------------------------------------------------------------*
        IF (LISTL(1:2).EQ.'L0') THEN

C          -----------------------------------------------------------
C          calculate packed L(ia,jb) integrals and evaluate the 
C          projection contribution <HF|[[H,T],t_mu1]|CC>
C          -----------------------------------------------------------
           CALL CCSD_TCMEPK(WORK(KXIAJB),ONE,ISYM0)

           IOPT = 0
           CALL CCG_LXD(WORK(KTHETA1),ISYRES,WORK(KT1AMPA),ISYAMP,
     &                  WORK(KXIAJB),ISYM0,IOPT)
           CALL DSCAL(NT1AM(ISYRES),TWO,WORK(KTHETA1),1)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) MSGDBG,
     *            'norm of THETA1 after <HF| contrib.:',
     &          DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
c             WRITE (LUPRI,*) MSGDBG,
c    *            'result vector after <HF| contribution:'
c             CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
           END IF
        END IF

*---------------------------------------------------------------------*
*       for CCSD calculate second 21G contribution
*---------------------------------------------------------------------*
        IF (CCSD) THEN

C          ----------------------------------------------------------
C          memory allocation, reusing space for packed (ia|jb) array:
C          ----------------------------------------------------------
           KXIDJL  = KLIAJB  + NT2SQ(ISYM0)
           KXJLID  = KXIDJL  + NTRAOC(ISYAMP)
           KEND3   = KXJLID  + NTRAOC(ISYAMP)
           LWRK3   = LWORK - KEND3

           IF (LWRK3 .LE. 0) THEN
             CALL QUIT('Insufficient work space in CC_FMAT. (21G)')
           END IF

C          ----------------------------------------------------------
C          calculate g_iljd = (ia|jd) * T_al and compute 21G contrib.
C          ----------------------------------------------------------
           CALL CCG_TRANS4(WORK(KXIDJL),ISYAMP,WORK(KLIAJB),ISYM0,
     &                     WORK(KT1AMPA),ISYAMP)

           IOPT = 5
           CALL CCG_SORT1(WORK(KXIDJL),WORK(KXJLID),ISYAMP,IOPT)

           CALL CC_21GMO(WORK(KTHETA1),ISYRES,WORK(KXJLID),ISYAMP,
     &                   LUPQMO,FNPQMO,IADRPQMO(1,IINTL),ISYCTR,
     &                   WORK(KEND3),LWRK3)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'Theta1 after 21GMO contribution:'
              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
           END IF
        END IF

*---------------------------------------------------------------------*
*       for CC2 & CCSD precalculate the YBARA and A2 intermediates
*       for CC2 precalculate also the XBARA intermediate
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN

           KT2AMPA = KLIAJB  + NT2SQ(ISYM0)
           KEND3   = KT2AMPA + MAX(NT2AM(ISYAMP),NT2AM(ISYM0))
           LWRK3   = LWORK - KEND3

           IF (LWRK3 .LE. 0) THEN
             CALL QUIT(
     &          'Insufficient work space in CC_FMAT. (YBARA/A2IM)')
           END IF

C          ----------------------------------------------
C          calculate Liajb and read response amplitudes:
C          ----------------------------------------------
           CALL CCRHS_T2TR(WORK(KLIAJB),WORK(KEND3),LWRK3,ISYM0)
         
           IOPT = 2
           CALL CC_RDRSP(LISTR,IDLSTR,ISYAMP,IOPT,MODEL,
     &                   WORK(KDUM),WORK(KT2AMPA))
           CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,ISYAMP)

C          ----------------------------------------------
C          XBARA_kj = sum_dem T_djem L_mekd
C          ----------------------------------------------
           IF (CC2) THEN
              CALL CC_XBAR(WORK(KXBARA),WORK(KLIAJB),ISYM0,
     &                     WORK(KT2AMPA),ISYAMP,WORK(KEND3),LWRK3)
           END IF

C          ----------------------------------------------
C          YBARA_ba = sum_dlk T_dlbk L_ldka
C          ----------------------------------------------
           CALL CC_YBAR(WORK(KYBARA),WORK(KLIAJB),ISYM0,
     &                  WORK(KT2AMPA),ISYAMP,WORK(KEND3),LWRK3)

C          ----------------------------------------------
C          A2IM = sum_bj (2T_aibj - T_ajbi) Zeta_bj
C          ----------------------------------------------
           CALL CCSD_TCMEPK(WORK(KT2AMPA),ONE,ISYAMP)
           CALL CCG_LXD(WORK(KA2IM),ISYRES,WORK(KZETA1),ISYCTR,
     &                  WORK(KT2AMPA),ISYAMP,0)

        END IF

*---------------------------------------------------------------------*
* allocate & initialize doubles result vector:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
          KTHETA2 = KEND2
          KEND2   = KTHETA2 + NT2AM(ISYRES)
          CALL DZERO(WORK(KTHETA2),NT2AM(ISYRES))
        END IF

*---------------------------------------------------------------------*
* allocate work space some `small' intermediates:
*---------------------------------------------------------------------*
        IF (CCS.OR.CC2) THEN
          KFOCKA  = KEND2
          KGZETA  = KFOCKA + N2BST(ISYAMP)
          KEND2   = KGZETA + N2BST(ISYRES)
        END IF

        IF (.NOT.CCS) THEN
          KRIMA   = KEND2
          KXINT   = KRIMA + NEMAT1(ISYAMP)
          KYINT   = KXINT + NMATIJ(ISYRES)
          KFINT   = KYINT + NMATAB(ISYRES)
          KEND2   = KFINT + NT1AO(ISYRES) 
        END IF

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LE. NT1AM(ISYRES)) THEN
          CALL QUIT('Insufficient work space in CC_FMAT. (8)')
        END IF

*---------------------------------------------------------------------*
* for CCS/CC2 read the response Fock matrix and the Fock matrix like
* CCS/CC2 GZeta intermediate and calculate the CC_11A contribution
*---------------------------------------------------------------------*
        IF (CCS.OR.CC2) THEN
           CALL CC_RVEC(LUFK,FILFCK,LENFK,N2BST(ISYAMP),IINTR,
     &                  WORK(KFOCKA))
           CALL CC_FCKMO(WORK(KFOCKA),WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK(KEND2),LWRK2,ISYAMP,ISYM0,ISYM0)

           CALL CC_RVEC(LUGIM,GIMFIL,LENGIM,N2BST(ISYRES),IINT2,
     &                  WORK(KGZETA))
           CALL CC_11A(WORK(KTHETA1),WORK(KGZETA),ISYRES,WORK(KLAMH0),
     &                 WORK(KLAMP0),WORK(KEND2),LWRK2)
          
           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'CCS GZETA intermediate:'
              CALL CC_PRFCKAO(WORK(KGZETA),ISYRES)
           END IF
        END IF
      
*---------------------------------------------------------------------*
* for CC2 and CCSD read precalculated X and Y intermediates from file,
* and add the precalculated contributions from Rho1 to result vector:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
           CALL CC_RVEC(LUXIM,FILXIM,LENX,NMATIJ(ISYRES),IINT2,
     &                  WORK(KXINT))

           CALL CC_RVEC(LUYIM,FILYIM,LENY,NMATAB(ISYRES),IINT2,
     &                  WORK(KYINT))
           
           CALL CC_RVEC(LURHO,RHOFIL,LENRHO,NT1AM(ISYRES),IINT2,
     &                  WORK(KEND2))

           CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KEND2),1,WORK(KTHETA1),1)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'restored X and Y intermediates:'
              CALL CC_PREI(WORK(KYINT),WORK(KXINT),ISYRES,1)
              WRITE (LUPRI,*)
     *               LUXIM,FILXIM,LENX,NMATIJ(ISYRES),IINT2
              WRITE (LUPRI,*)
     *               LUYIM,FILYIM,LENY,NMATAB(ISYRES),IINT2

              WRITE (LUPRI,*) 'Rho1 read from file:'
              CALL CC_PRP(WORK(KEND2),WORK(KEND2),ISYRES,1,0)
              WRITE (LUPRI,*) 'result vector after Rho1 was added:'
              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)

              CALL FLSHFO(LUPRI) 
           END IF

        END IF

*---------------------------------------------------------------------*
* for CCSD read also the F and R intermediates;
* the F intermediate is transformed to MO and added to result vector
*---------------------------------------------------------------------*
        IF (.NOT.(CCS.OR.CC2)) THEN
           CALL CC_RVEC(LUFIM,FIMFIL,LENFIM,NT1AO(ISYRES),IINT2,
     &                  WORK(KFINT))
           CALL CC_T1AM(WORK(KTHETA1),ISYRES,WORK(KFINT),ISYRES,
     &                  WORK(KLAMH0),ISYM0,ONE)

           CALL CC_RVEC(LURIM,FILRIM,LENR,NEMAT1(ISYAMP),IINTR,
     &                  WORK(KRIMA))
        END IF

*---------------------------------------------------------------------*
* calculate the remaining contributions:
*---------------------------------------------------------------------*

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,'norm of THETA1 before CCFMAT2:',
     &       DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           IF (.NOT.CCS) 
     &      WRITE (LUPRI,*) MSGDBG,'norm of THETA2 before CCFMAT2:',
     &       DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
           WRITE (LUPRI,*) 'result vector before entering CCFMAT2:'
           CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
        END IF

        CALL CCFMAT2SC(WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &               LISTL,IDLSTL, WORK(KZETA1), ISYCTR,
     &               LISTR,IDLSTR, WORK(KT1AMPA),ISYAMP,
     &               WORK(KLAMDPA),WORK(KLAMDHA),
     &               WORK(KLAMP0), WORK(KLAMH0),
     &               WORK(KXINT),  WORK(KYINT),
     &               WORK(KXBARA), WORK(KYBARA), 
     &               WORK(KRIMA),  WORK(KA2IM),
     &               WORK(KFOCKA), WORK(KFOCK0AO), WORK(KFCKC0),
     &               LUBFI,FNBFI,IADRBFI(1,IINT2),
     &               LUBF, BFFIL,LENBF,IINTR,
     &               LUC,  CTFIL, LUCBAR, CBAFIL,
     &               LUD,  DTFIL, LUDBAR, DBAFIL,
     &               IOFFCD(IINTR), WORK(KEND2),LWRK2)

*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KTHETA0 = -999999

      IF (IOPTRES .EQ. 0  .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct access file, 
*       store start address in IFTRAN(3,ITRAN)

        IFTRAN(3,ITRAN) = IADRTH

        CALL PUTWA2(LUFMAT,FILFMA,WORK(KTHETA1),IADRTH,NT1AM(ISYRES))
        IADRTH = IADRTH + NT1AM(ISYRES)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUFMAT,FILFMA,WORK(KTHETA2),IADRTH,NT2AM(ISYRES))
          IADRTH = IADRTH + NT2AM(ISYRES)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 
     &        'F matrix transformation nb. ',ITRAN,' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IFTRAN(3,ITRAN),IADRTH-IFTRAN(3,ITRAN)
         XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &           DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('F matrix transformation written to file:')
         NDBLE = 1
         IF (CCS) NDBLE = 0
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,NDBLE)
        END IF

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN

*        write to a sequential file by a call to CC_WRRSP/CC_WARSP,
*        use FILFMA as LIST type and IFTRAN(3,ITRAN) as index

         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
         ELSE IF (IOPTRES.EQ.4) THEN
           CALL CC_WARSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 
     &              'Write ',LISTL,' * B * ',LISTR,' transformation',
     &              ' as ',FILFMA,' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. left  vector:',IFTRAN(1,ITRAN)
           WRITE (LUPRI,*) 'index of inp. right vector:',IFTRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of     result vector:',IFTRAN(3,ITRAN)
           XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           IF (.NOT.CCS) XNORM = XNORM +
     &       DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
         END IF
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KTHETA2),TWO,ISYRES)
         CALL CCDOTRSP(IFDOTS,FCONS,IOPTW,FILFMA,ITRAN,NFTRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                 WORK(KEND2),LWRK2)
      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_FMAT.')
      END IF
      
      TIMIO = TIMIO + SECOND() - DTIME

      TIMTRN = SECOND() - TIMTRN

      IF (IPRINT.GT.0) THEN

         IF (IOPTRES.EQ.5) THEN
            IVEC = 1
            DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
               IVEC = IVEC + 1
            END DO    
            WRITE (LUPRI,'(1X,2(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTL,
     &        '    | ',IDLSTR,'    | ',IVEC-1,'       | ',TIMTRN,'  |'
         ELSE
            WRITE (LUPRI,
     &        '(1X,2(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTL, '    | '
     &        ,IDLSTR,'    | ',IFTRAN(3,ITRAN),'       | ',TIMTRN,'  |'
         END IF 

      END IF

*---------------------------------------------------------------------*
* End of loop over F matrix transformations
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* close & remove scratch files:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUBF,   BFFIL,  'DELETE')
         CALL WCLOSE2(LUC,    CTFIL,  'DELETE')
         CALL WCLOSE2(LUD,    DTFIL,  'DELETE')
         CALL WCLOSE2(LUCBAR, CBAFIL, 'DELETE')
         CALL WCLOSE2(LUDBAR, DBAFIL, 'DELETE')
         CALL WCLOSE2(LURIM,  FILRIM, 'DELETE')
         CALL WCLOSE2(LUFIM,  FIMFIL, 'DELETE')
         CALL WCLOSE2(LUXIM,  FILXIM, 'DELETE')
         CALL WCLOSE2(LUYIM,  FILYIM, 'DELETE')
         CALL WCLOSE2(LUPQMO, FNPQMO, 'DELETE')
         CALL WCLOSE2(LUBFI,  FNBFI,  'DELETE')
      ELSE
         CALL WCLOSE2(LUGIM,  GIMFIL, 'DELETE')
         CALL WCLOSE2(LUFK,   FILFCK, 'DELETE')
      END IF

      IF (.NOT.CCS)  CALL WCLOSE2(LURHO,RHOFIL, 'DELETE')
      IF (CCSD)      CALL WCLOSE2(LUBFD, FNBFD, 'DELETE')


      TIMIO = TIMIO + SECOND() - DTIME

*---------------------------------------------------------------------*
* if IOPTRES=1 and enough work space available, read result
* vectors back into memory:
*---------------------------------------------------------------------*
      DTIME = SECOND()

* check size of work space:
      IF (IOPTRES .EQ. 1) THEN
        LENALL = IADRTH-1
        IF (LENALL .GT. LWORK) IOPTRES = 0
      END IF

* read the result vectors back into memory:
      IF (IOPTRES .EQ. 1) THEN

        CALL GETWA2(LUFMAT,FILFMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NFTRAN
            IF (ITRAN.LT.NFTRAN) THEN
              LEN     = IFTRAN(3,ITRAN+1)-IFTRAN(3,ITRAN)
            ELSE
              LEN     = IADRTH-IFTRAN(3,NFTRAN)
            END IF
            KTHETA1 = IFTRAN(3,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read F matrix transformation nb. ',NFTRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IFTRAN(3,NFTRAN),LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

      TIMIO = TIMIO + SECOND() - DTIME
*---------------------------------------------------------------------*
* close F matrix file, print timings & return
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUFMAT, FILFMA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUFMAT, FILFMA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_FMAT.')
      END IF

      TIMIO  = TIMIO + SECOND() - DTIME
      TIMALL = SECOND() - TIMALL
      CONVRT = 100.0D0/TIMALL

      IF (IPRINT.GT.0) THEN
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         WRITE (LUPRI,'(1X,A,I4,A,F10.2,A)')
     &     '| total time for',NFTRAN,' F transforms.:',TIMALL,' secs.|'
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         WRITE (LUPRI,
     &      '(1X,"|  % of time used in ",A,":",F10.2,"      |")')
     &      'start up org.', TIMPRE*CONVRT,
     &      'Fock interm. ', TIMFCK*CONVRT,
     &      'ERI          ', TIMINT*CONVRT,
     &      'CCRDAO       ', TIMRDAO*CONVRT,
     &      '(**|k del)   ', TIMTRBT*CONVRT,
     &      'BF term      ', TIMBF*CONVRT,
     &      'C term       ', TIMC*CONVRT,
     &      'D term       ', TIMD*CONVRT,
     &      'I/O          ', TIMIO*CONVRT

         WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'
      END IF
*---------------------------------------------------------------------*
* that's it; return:
*---------------------------------------------------------------------*

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCFMAT2 */
*=====================================================================*
       SUBROUTINE CCFMAT2SC(THETA1, THETA2,  ISYRES,
     &                    LISTL,  IDLSTL,  ZETA1,  ISYCTR,  
     &                    LISTR,  IDLSTR,  T1AMPA, ISYAMP,   
     &                    XLAMPA, XLAMHA,  XLAMP0, XLAMH0,  
     &                    XINT,   YINT,    XBARA,  YBARA, 
     &                    RIMA,   A2IM,    FOCKA,  FOCK0AO, FCKC0,
     &                    LUBFI,  FNBFI,   IADRBFI,
     &                    LUBF,   BFFIL,   LENBF,  IINTR,
     &                    LUC,    CTFIL,   LUCBAR, CBAFIL,
     &                    LUD,    DTFIL,   LUDBAR, DBAFIL,
     &                    IOFFCD, WORK,    LWORK                    )
*---------------------------------------------------------------------*
*
*    Purpose: calculate final contributions to F matrix transformation
*             from precalculated intermediates
*
*     Written by Christof Haettig, November 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <dummy_dec.h>
#  include <ccsdsym_dec.h>
#  include <maxorb_dec.h>
#  include <ccorb_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccfield.h>
#include <dummy.h>

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CCFMAT2> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0, LUBF0
      PARAMETER (ISYM0 = 1)
      
      CHARACTER*(*) LISTL, LISTR, FNBFI, BFFIL
      CHARACTER*(*) CTFIL, CBAFIL, DTFIL, DBAFIL
      INTEGER IDLSTL, IDLSTR, ISYRES, ISYCTR, ISYAMP, LWORK
      INTEGER LUBFI, IADRBFI(MAXORB), IOFFCD
      INTEGER LUC, LUCBAR, LUD, LUDBAR, LUBF, LENBF, IINTR


#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL THETA1(*), THETA2(*), ZETA1(*), T1AMPA(*)
      REAL XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      REAL XINT(*), YINT(*), XBARA(*), YBARA(*)
      REAL RIMA(*), A2IM(*), FOCKA(*), FOCK0AO(*), FCKC0(*)
      REAL XNORM, DDOT, ZERO, ONE, HALF, FACB, FF
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION THETA1(*), THETA2(*), ZETA1(*), T1AMPA(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      DOUBLE PRECISION XINT(*), YINT(*), XBARA(*), YBARA(*), RIMA(*)
      DOUBLE PRECISION A2IM(*), FOCKA(*), FOCK0AO(*), FCKC0(*)
      DOUBLE PRECISION XNORM, DDOT, ZERO, ONE, HALF, FACB, FF
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0)

      LOGICAL LGAMMA, LBFZETA, LO3BF, LRCON, LGCON, FCKCON
      CHARACTER*(10) MODEL
      INTEGER KGZETA, KGIMA, KEND0, LWRK0, ISYMB, ISYMJ, KFOCKA
      INTEGER IOPTG, ICON, KEND1, LWRK1, KBZETA, KFOCK
      INTEGER IOPT, KGAMMA, KBF, KZETA2, KOFF1, KOFF2, KLIAJB
      INTEGER KYPS, KCHI, ISYMA, ISYMD, ISYMK, NVIRA, NVIRD, KOFF3
      INTEGER IOPTE, IOPTB, KCDBAR, KEMAT1, KEMAT2
      INTEGER KONEHR, KONEHG, KONEH, KA2CON
      INTEGER NRHFK, NVIRC, ISYMC, ISYMI, IF
      INTEGER KFIELD, KFIELDR, KFIELDG

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (CCS) THEN
         KONEHR  = 1
         KONEHG  = KONEHR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEMAT1  = KONEHG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEMAT2  = KEMAT1 + NMATAB(ISYAMP)
         KEND0   = KEMAT2 + NMATIJ(ISYAMP)
      ELSE 
         KGZETA  = 1
         KGIMA   = KGZETA + NT1AO(ISYRES)
         KFOCKA  = KGIMA  + NT1AO(ISYAMP)
         KYPS    = KFOCKA + NT1AM(ISYAMP)
         KEMAT1  = KYPS   + NMATAB(ISYRES)
         KEMAT2  = KEMAT1 + MAX(NMATAB(ISYAMP),NEMAT1(ISYAMP))
         KONEH   = KEMAT2 + NMATIJ(ISYAMP)
         KONEHR  = KONEH  + N2BST(ISYM0)
         KONEHG  = KONEHR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KA2CON  = KONEHG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEND0   = KA2CON + NT1AM(ISYRES)
      END IF

      IF (CC2.AND.NONHF) THEN
         KFIELD  = KEND0
         KFIELDR = KFIELD  + N2BST(ISYM0)
         KFIELDG = KFIELDR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEND0   = KFIELDG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
      END IF

      LWRK0 = LWORK - KEND0
      IF (LWRK0.LT.0) THEN
         CALL QUIT('Insufficient memory in CCFMAT2. (0)')
      END IF

*---------------------------------------------------------------------*
* for CCS complete here the calculation of the FOCK^A matrix and
* put it into ONEHG and ONEHR arrays used for construction of E1 & E2:
*---------------------------------------------------------------------*
      IF (CCS.OR.CC2) THEN

         CALL DCOPY(N2BST(ISYM0),FOCK0AO,1,WORK(KONEHR),1)
         CALL DCOPY(N2BST(ISYM0),FOCK0AO,1,WORK(KONEHG),1)

         CALL CC_FCKMO(WORK(KONEHR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         
         CALL CC_FCKMO(WORK(KONEHG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)

         CALL DAXPY(N2BST(ISYAMP),ONE,WORK(KONEHG),1,WORK(KONEHR),1)
         CALL DAXPY(N2BST(ISYAMP),ONE,FOCKA,       1,WORK(KONEHR),1)

         CALL DCOPY(N2BST(ISYAMP),WORK(KONEHR),1,WORK(KONEHG),1)
         
      END IF

*---------------------------------------------------------------------*
* get one-electron hamiltonian and transform to MO:
* if frozen core used, add core-only zeroth-order Fock matrix
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN
         CALL CCRHS_ONEAO(WORK(KONEH),WORK(KEND0),LWRK0)
         DO IF = 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KONEH),WORK(KEND0),LWRK0,FF,1,LFIELD(IF))
         END DO
         IF (FROIMP.OR.FROEXP) THEN
           CALL DAXPY(N2BST(ISYM0),ONE,FCKC0,1,WORK(KONEH),1)
         END IF

         CALL DCOPY(N2BST(ISYM0),WORK(KONEH),1,WORK(KONEHR),1)
         CALL DCOPY(N2BST(ISYM0),WORK(KONEH),1,WORK(KONEHG),1)
    
         CALL CC_FCKMO(WORK(KONEH),XLAMP0,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYM0)
         CALL CC_FCKMO(WORK(KONEHR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         CALL CC_FCKMO(WORK(KONEHG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)
      END IF

*---------------------------------------------------------------------*
* get external fields added after the SCF and transform to MO:
*---------------------------------------------------------------------*
      IF (CC2.AND.NONHF) THEN
         CALL DZERO(WORK(KFIELD),N2BST(ISYM0))
         DO IF = 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KFIELD),WORK(KEND0),LWRK0,FF,1,LFIELD(IF))
         END DO

         CALL DCOPY(N2BST(ISYM0),WORK(KFIELD),1,WORK(KFIELDR),1)
         CALL DCOPY(N2BST(ISYM0),WORK(KFIELD),1,WORK(KFIELDG),1)
    
         CALL CC_FCKMO(WORK(KFIELD),XLAMP0,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYM0)
         CALL CC_FCKMO(WORK(KFIELDR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         CALL CC_FCKMO(WORK(KFIELDG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)
      END IF

*---------------------------------------------------------------------*
* calculate BZeta contribution and GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD) THEN
         KBZETA = KEND0
         KEND1  = KBZETA + NT2AO(ISYRES)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
 
         CALL CC_BFIFSORT(WORK(KBZETA),ISYRES,LUBFI,FNBFI,IADRBFI,
     &                    WORK(KEND1),LWRK1)

         CALL DZERO(WORK(KGZETA),NT1AO(ISYRES))

         ICON    = 1
         IOPTG   = 1
         LGAMMA  = .FALSE.
         LO3BF   = .FALSE.
         LBFZETA = .TRUE.
         CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBZETA),
     &                 THETA2,DUMMY,WORK(KGZETA),DUMMY,
     &                 XLAMH0,ISYM0,XLAMH0,ISYM0,
     &                 WORK(KEND1),LWRK1,ISYRES,
     &                 ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
      
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,'norm of THETA1 after B+E term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,'norm of THETA2 after B+E term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C           CALL AROUND('THETA2 after B+E term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the E1' contribution from the GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD) THEN

         CALL CC_T1AM(THETA1,ISYRES,WORK(KGZETA),ISYRES,
     &                XLAMH0,ISYM0,ONE)

         IF (LOCDBG) THEN
            CALL AROUND("THETA1 after E1' term:")
            CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate D2 contribution to Theta1 results vector: 
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN

         CALL CC_21EFM(THETA1,WORK(KONEH),ISYM0,XINT,YINT,ISYRES,
     &                 WORK(KEND0),LWRK0)

      ELSE IF (CC2.AND.NONHF) THEN

         CALL CC_21EFM(THETA1,WORK(KFIELD),ISYM0,XINT,YINT,ISYRES,
     &                 WORK(KEND0),LWRK0)

      END IF

c      IF (LOCDBG .AND. .NOT.CCS) THEN
         CALL AROUND('THETA after singles excit. D2 term:')
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
c      END IF

*---------------------------------------------------------------------*
* for CC2/NONHF and CCSD read double excitation part of the Lagrangian 
* multipliers and expand to a full square matrix:
*---------------------------------------------------------------------*
      IF ((CC2.AND.NONHF) .OR. CCSD) THEN
        KZETA2 = KEND0
        KEND1  = KZETA2 + NT2SQ(ISYCTR)
        LWRK1  = LWORK  - KEND1
      
        IF (LWRK1.LT.NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KZETA2),ISYCTR)

      ELSE
        KEND1 = KEND0
        LWRK1 = LWORK  - KEND1
      END IF

*---------------------------------------------------------------------*
* for CC2/NONHF calculate additional FF term to Theta2:
*---------------------------------------------------------------------*
      IF (CC2.AND.NONHF) THEN

         IOPT   = 1
         CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &               DUMMY,DUMMY,DUMMY,DUMMY,
     &               WORK(KFIELDR),WORK(KFIELDG),
     &               DUMMY,DUMMY,IDUMMY,DUMMY,DUMMY,IDUMMY,
     &               .TRUE.,.FALSE.,.FALSE.,.FALSE.,IOPT,ISYAMP)

         CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),WORK(KEND1),LWRK1,
     &                ISYAMP)

         CALL CCRHS_E(THETA2,WORK(KZETA2),WORK(KEMAT1),WORK(KEMAT2),
     &                WORK(KEND1),LWRK1,ISYCTR,ISYAMP)

      END IF
*---------------------------------------------------------------------*
* calculate the first E2 term contribution from the zero-order BF:
*---------------------------------------------------------------------*
      IF (CCSD) THEN
        KZETA2 = KEND0
        KBF    = KZETA2 + NT2SQ(ISYCTR)
        KEND1  = KBF    + 2 * NT2ORT(ISYM0)
        LWRK1  = LWORK  - KEND1
     
        IF (LWRK1.LT.0) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        LUBF0 = 0
        CALL GPOPEN(LUBF0,'CC_BFIM','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        READ(LUBF0) (WORK(KBF-1+I),I=1,2*NT2ORT(ISYM0))
        CALL GPCLOSE(LUBF0,'KEEP')


        ICON    = 3 
        IOPTG   = 0
        LGAMMA  = .FALSE.
        LO3BF   = .FALSE.
        LBFZETA = .FALSE.
        CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBF),
     &                DUMMY,DUMMY,DUMMY,DUMMY,
     &                XLAMP0,ISYM0,XLAMPA,ISYAMP,
     &                WORK(KEND1),LWRK1,ISYM0,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
     
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,
     &       'norm(THETA1) after first E2 (LT21C) term:',
     &       DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
           IF (.NOT.CCS) 
     &      WRITE (LUPRI,*) MSGDBG,
     &       'norm(THETA2) after first E2 (LT21C) term:',
     &       DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C          CALL AROUND('THETA after first E2 (LT21C) term:')
C          CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
        END IF

      END IF

*---------------------------------------------------------------------*
* calculate second E2 term contribution and the GAMMA and G 
* intermediates from the response BF intermediate: 
*---------------------------------------------------------------------*
      IF (CCSD) THEN
        KZETA2 = KEND0
        KBF    = KZETA2 + NT2SQ(ISYCTR)
        KGAMMA = KBF    + NT2AOIJ(ISYAMP)
        KEND1  = KGAMMA + NGAMMA(ISYAMP)
        LWRK1  = LWORK  - KEND1
      
        IF (LWRK1.LT.NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        CALL CC_RVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYAMP),IINTR,WORK(KBF))

        CALL DZERO(WORK(KGAMMA),NGAMMA(ISYAMP))
        CALL DZERO(WORK(KGIMA), NT1AO(ISYAMP))


        ICON    = 3
        IOPTG   = 2
        LGAMMA  = .TRUE.
        LO3BF   = .TRUE.
        LBFZETA = .FALSE.
        CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBF),
     &                DUMMY,WORK(KGAMMA),WORK(KGIMA),DUMMY,
     &                XLAMP0,ISYM0,XLAMP0,ISYM0,
     &                WORK(KEND1),LWRK1,ISYAMP,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
      
cs        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,
     &      'norm(THETA1) after second E2 (LT21C) term:',
     &      DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
          WRITE (LUPRI,*) MSGDBG,
     &      'norm(THETA2) after second E2 (LT21C) term:',
     &      DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) MSGDBG,'the GIMA intermediate:'
          WRITE (LUPRI,'(5F12.6)') (WORK(KGIMA-1+I),I=1,NT1AO(ISYAMP))
          CALL AROUND('THETA after second E2 (LT21C) term:')
          CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
          XNORM = DDOT(NGAMMA(ISYAMP),WORK(KGAMMA),1,WORK(KGAMMA),1)
          WRITE (LUPRI,*) 'Norm(GAMMA) = ',XNORM
cs        END IF

      END IF

*---------------------------------------------------------------------*
* for CCS/CCSD calculate here E1/E1* and E2 intermediates as nedded 
* for the B2 contribution from the G and YBAR intermediates:
*---------------------------------------------------------------------*
      LRCON  = .FALSE.
      LGCON  = .FALSE.
      FCKCON = .TRUE.
      IOPT   = 1

      IF (CCSD) LGCON  = .TRUE.

      CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &            DUMMY,DUMMY,WORK(KGIMA),DUMMY,
     &            WORK(KONEHR),WORK(KONEHG),
     &            XLAMH0,XLAMP0,ISYM0,DUMMY,DUMMY,IDUMMY,
     &            FCKCON,LRCON,LGCON,.FALSE.,IOPT,ISYAMP)

      IF (CC2) THEN
         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)
         CALL DAXPY(NMATIJ(ISYAMP),ONE,XBARA,1,WORK(KEMAT2),1)
      END IF

      IF (CCSD) THEN
         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)
      END IF

      IF (LOCDBG) THEN
         CALL AROUND('E-intermediates out from CC_EIM:')
         CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYAMP,1)
      END IF

*---------------------------------------------------------------------*
* calculate B2 contribution to Theta1 results vector: 
*---------------------------------------------------------------------*
      CALL CCLR_E1C1(THETA1,ZETA1,WORK(KEMAT1),WORK(KEMAT2),
     &               WORK(KEND1),LWRK1,ISYCTR,ISYAMP,'T')

c      IF (LOCDBG) THEN
         CALL AROUND('THETA1 B2 term:')
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
c      END IF

*---------------------------------------------------------------------*
* for CCSD calculate here the E intermed. from G, R and YBAR intermed.:
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN

         LRCON  = .TRUE.
         LGCON  = .TRUE.
         FCKCON = .TRUE.
         IOPT   = 1
         CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &               RIMA,DUMMY,WORK(KGIMA),DUMMY,
     &               WORK(KONEHR),WORK(KONEHG),
     &               XLAMH0,XLAMP0,ISYM0,DUMMY,DUMMY,IDUMMY,
     &               FCKCON,LRCON,LGCON,.FALSE.,IOPT,ISYAMP)

         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)

         IF (LOCDBG) THEN
            CALL AROUND('E-intermediates out from CC_EIM:')
            CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYAMP,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the A term and add to Theta2 result vector:
*---------------------------------------------------------------------*
      IF (CCSD) THEN
         CALL CC_GAMMA2(WORK(KGAMMA),WORK(KEMAT2),ISYAMP)

         IOPT = 2
         CALL CCRHS_A(THETA2,WORK(KZETA2),WORK(KGAMMA),
     &                WORK(KEND1),LWRK1,ISYAMP,ISYCTR,IOPT)
      
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA1) after doubles excit A term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA2) after doubles excit A term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C            CALL AROUND('THETA after doubles excit. A term:')
C            CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* read (ia|jb) integrals and calculate L(ia|jb) in place:
* (needed for H term and for A2 term)
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
         KLIAJB = KEND0
         KEND1  = KLIAJB + NT2AM(ISYM0)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
      
         CALL CCG_RDIAJB(WORK(KLIAJB),NT2AM(ISYM0))
         CALL CCSD_TCMEPK(WORK(KLIAJB),ONE,ISYM0)
      END IF

*---------------------------------------------------------------------*
* calculate A2 contribution:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
     
         CALL CCG_LXD(WORK(KA2CON),ISYRES,A2IM,ISYRES,
     &                WORK(KLIAJB),ISYM0,0)

         CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KA2CON),1,THETA1,1)

         IF (LOCDBG) THEN
            CALL AROUND('THETA1 A2 term:')
            CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
            CALL AROUND('A2IM intermediate:')
            CALL CC_PRP(A2IM,THETA2,ISYRES,1,0)
            CALL AROUND('A2CON contribution:')
            CALL CC_PRP(WORK(KA2CON),THETA2,ISYRES,1,0)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate H term:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
         KFOCK  = KEND1
         KCHI   = KFOCK  + N2BST(ISYAMP)
         KEND1  = KCHI   + NMATIJ(ISYRES)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
      
C        ------------------------
C        calculate F^A_jb matrix:
C        ------------------------
         IOPT = 0
         CALL CCG_LXD(WORK(KFOCKA),ISYAMP,T1AMPA,ISYAMP,
     &                WORK(KLIAJB),ISYM0,IOPT)

         IF (LOCDBG) THEN
            CALL AROUND('FOCKA_jb matrix:')
            CALL CC_PRP(WORK(KFOCKA),WORK,ISYRES,1,0)
         END IF

CTST
         CALL DZERO(WORK(KFOCK),N2BST(ISYAMP))
CTST

C        -----------------------------------------
C        resort Fock^A_jb into a full Fock matrix:
C        -----------------------------------------
         DO ISYMB = 1, NSYM
            ISYMJ = MULD2H(ISYAMP,ISYMB)

            DO J = 1, NRHF(ISYMJ)
            DO B = 1, NVIR(ISYMB)
               KOFF1 = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B-1) + J
               KOFF2 = IT1AM(ISYMB,ISYMJ)  + NVIR(ISYMB)*(J-1) + B
               WORK(KFOCK-1+KOFF1) = WORK(KFOCKA-1+KOFF2)
            END DO
            END DO

         END DO

C        ---------------------------------------
C        calculate the simple Fock contribution:
C        ---------------------------------------
         CALL CC_L1FCK(THETA2,ZETA1,WORK(KFOCK),ISYCTR,ISYAMP,
     &                 WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA1) after doubles excit H term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA2) after doubles excit H term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C            CALL AROUND('THETA after doubles exci. H term (Fock part):')
C            CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
           XNORM = DDOT(NT1AM(ISYCTR),ZETA1,1,ZETA1,1)
           WRITE (LUPRI,*) 'Norm of ZETA1:',XNORM
         END IF

C        ----------------------------------------------
C        calculate the transposed Ypsilon intermediate:
C        (for CC2 exclude the contribution from YINT)
C        ----------------------------------------------
         IF (CC2) THEN
            CALL DZERO(WORK(KYPS),NMATAB(ISYRES))
         ELSE
            CALL DCOPY(NMATAB(ISYRES),YINT,1,WORK(KYPS),1)
         END IF

         DO ISYMA = 1, NSYM

            ISYMD = MULD2H(ISYMA,ISYRES)
            ISYMK = MULD2H(ISYMA,ISYCTR)

            KOFF1 = IT1AM(ISYMD,ISYMK) + 1
            KOFF2 = IT1AM(ISYMA,ISYMK) + 1
            KOFF3 = KYPS + IMATAB(ISYMD,ISYMA) 

            NVIRD = MAX(NVIR(ISYMD),1)
            NVIRA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('N','T',NVIR(ISYMD),NVIR(ISYMA),NRHF(ISYMK),
     &                 ONE,T1AMPA(KOFF1),NVIRD,ZETA1(KOFF2),NVIRA,
     &                 ONE,WORK(KOFF3),NVIRD)
         END DO    
      
         CALL DSCAL(NMATAB(ISYRES),HALF,WORK(KYPS),1)

C        ------------------------------------------------------
C        calculate the Chi intermediate:
C        (for CCSD this contribution is included in the B term)
C        ------------------------------------------------------
         CALL DZERO(WORK(KCHI),NMATIJ(ISYRES))

         IF (CC2) THEN

           DO ISYMK = 1, NSYM
              ISYMC = MULD2H(ISYMK,ISYAMP)
              ISYMI = MULD2H(ISYMC,ISYCTR)

              KOFF1 = IT1AM(ISYMC,ISYMK) + 1
              KOFF2 = IT1AM(ISYMC,ISYMI) + 1
              KOFF3 = KCHI + IMATIJ(ISYMK,ISYMI)

              NRHFK = MAX(NRHF(ISYMK),1)
              NVIRC = MAX(NVIR(ISYMC),1)

              CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *                   ONE,T1AMPA(KOFF1),NVIRC,ZETA1(KOFF2),NVIRC,
     *                   ONE,WORK(KOFF3),NRHFK)

           END DO

           CALL DSCAL(NMATIJ(ISYRES),HALF,WORK(KCHI),1)

         END IF
         
C        --------------------------------
C        contract with (jb|id) integrals:
C        --------------------------------
         CALL CC_22EC(THETA2,WORK(KLIAJB),WORK(KCHI),WORK(KYPS),
     &                ISYRES,WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA1) after doubles excit H term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA2) after doubles excit H term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C           CALL AROUND('THETA after doubles excit. H term (2. part):')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
            CALL AROUND('the Ypsilon/Chi intermediates:')
            CALL CC_PREI(WORK(KYPS),WORK(KCHI),ISYRES,1)
            CALL AROUND('the Y and X intermediates:')
            CALL CC_PREI(YINT,XINT,ISYRES,1)
C           CALL AROUND('the integrals:')
C           CALL CC_PRSQ(WORK,WORK(KLIAJB),ISYM0,0,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the D & C term contributions:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS.OR.CC2)) THEN

        KCDBAR = KEND0
        KZETA2 = KCDBAR  + NT2SQ(ISYAMP)
        KEND1  = KZETA2  + NT2AM(ISYCTR)
        LWRK1  = LWORK   - KEND1

        IF (LWRK1 .LT. NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                DUMMY,WORK(KZETA2))

C       -----------------------------------------------
C       complete the calculation of the D intermediate:
C       (including the E1 contribution to the diagonal)
C       -----------------------------------------------
        IOPT  = 3
        IOPTB = 1
        IOPTE = 1
        FACB  = ONE
        CALL CCRHS_DIO3(WORK(KCDBAR),DUMMY,XLAMH0,
     &                  WORK(KEND1),LWRK1,1,ISYAMP,
     &                  LUD,DTFIL,LUC,CTFIL,IINTR,IOPT,
     &                  IOPTB,LUDBAR,DBAFIL,IOFFCD,FACB,
     &                  IOPTE,WORK(KEMAT1))

C       -----------------------------------------------
C       contract with the lagrangian multiplier vector:
C       -----------------------------------------------
        CALL CC_CD('D',-1,THETA2,ISYRES,WORK(KZETA2),ISYCTR,
     &             WORK(KCDBAR),ISYAMP,WORK(KEND1),LWRK1)

        IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA1) after doubles excit D term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA2) after doubles excit D term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(D int.) after doubles excit D term:',
     &        DDOT(NT2SQ(ISYAMP),WORK(KCDBAR),1,WORK(KCDBAR),1)
C           CALL AROUND('THETA2 after D term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
        END IF

C       -----------------------------------------------
C       complete the calculation of the C intermediate:
C       (including the E1 contribution to the diagonal)
C       -----------------------------------------------
        IOPT  = 3
        IOPTB = 1
        IOPTE = 1
        FACB  = ONE
        CALL CCRHS_CIO3(WORK(KCDBAR),DUMMY,XLAMH0,
     &                  WORK(KEND1),LWRK1,1,ISYAMP,
     &                  LUC,CTFIL,IINTR,IOPT,
     &                  IOPTB,LUCBAR,CBAFIL,IOFFCD,FACB,
     &                  IOPTE,WORK(KEMAT1))
       
C       -----------------------------------------------
C       contract with the lagrangian multiplier vector:
C       -----------------------------------------------
        CALL CC_CD('C',-1,THETA2,ISYRES,WORK(KZETA2),ISYCTR,
     &             WORK(KCDBAR),ISYAMP,WORK(KEND1),LWRK1)

        IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA1) after doubles excit C term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(THETA2) after doubles excit C term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
            WRITE (LUPRI,*) MSGDBG,
     &        'norm(C int.) after doubles excit C term:',
     &        DDOT(NT2SQ(ISYAMP),WORK(KCDBAR),1,WORK(KCDBAR),1)
C           CALL AROUND('THETA2 after C term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
        END IF

      END IF

*---------------------------------------------------------------------*
* that's it, return:
*---------------------------------------------------------------------*

      RETURN

      END
*=====================================================================*
