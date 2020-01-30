*=====================================================================*
c /* deck ccxopa_eomsetup */
*=====================================================================*
      SUBROUTINE CCXOPA_EOMSETUP(IFTRAN, IFDOTS, FCONS,  
     &                       NFTRAN, MXFTRAN, MXFVEC,
     &                       IEATRAN, IEADOTS, EACONS, 
     &                       NXE1TRAN,MXATRAN, MXAVEC,
     &                       IXE2TRAN,IX2DOTS, X2CONS, 
     &                       NXE2TRAN,MXXTRAN, MXXVEC,
     &                       IEOMXE2TRAN,IEOMX2DOTS, EOMX2CONS, 
     &                       NEOMXE2TRAN,MXEOMXTRAN, MXEOMXVEC,
     &                       IEOML0TRAN,IEOML0DOTS,EOML0CONS,
     &                       NEOML0TRAN,MXEOML0TRAN,MXEOML0VEC,
     &                       RESULT, MXOPA, LADD, WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC first-order transition moments
*         - list of B matrix transformations with eigenvectors
*         - list of A{X} matrix transformations with eigenvectors
*         - list of XKSI vector contractions with Nbar multipliers
*
*     Written by Christof Haettig, Oct 2003 
*
*     For EOM XOPA
*         - list of XKSI vector contractions with LE vectors
*           and of dot products tbar*RE)
*     Sonia Coriani, Nov 2015
*
*=====================================================================*
      IMPLICIT NONE  
#include "priunit.h"
#include "cclists.h"
#include "ccxopainf.h"
#include "ccroper.h"
#include "ccexci.h"
#include "ccsdinp.h"
#include "ccorb.h"

* local parameters:
      CHARACTER*(22) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCXOPA_EOMSETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)

      LOGICAL LADD
      INTEGER MXOPA,MXFTRAN,MXFVEC,MXATRAN,MXAVEC,MXXTRAN,MXXVEC
      INTEGER MXEOMXTRAN,MXEOMXVEC
      INTEGER MXEOML0TRAN,MXEOML0VEC

      INTEGER IFTRAN(MXDIM_FTRAN,MXFTRAN)
      INTEGER IFDOTS(MXFVEC,MXFTRAN)
      INTEGER IEATRAN(MXDIM_XEVEC,MXATRAN)
      INTEGER IEADOTS(MXAVEC,MXATRAN)
      INTEGER IXE2TRAN(MXDIM_XEVEC,MXXTRAN)
      INTEGER IX2DOTS(MXXVEC,MXXTRAN)
!
! EOM is one of these
!
      INTEGER IEOMXE2TRAN(MXDIM_XEVEC,MXEOMXTRAN)
      INTEGER IEOMX2DOTS(MXEOMXVEC,MXEOMXTRAN)
      INTEGER IEOML0TRAN(MXEOML0TRAN)
      INTEGER IEOML0DOTS(MXEOML0VEC,MXEOML0VEC)

      INTEGER NFTRAN, NXE1TRAN, NXE2TRAN, LWORK
      INTEGER NEOMXE2TRAN, NEOML0TRAN

#if defined (SYS_CRAY)
      REAL RESULT(MXOPA)
      REAL FCONS(MXFVEC,MXFTRAN)
      REAL EACONS(MXAVEC,MXATRAN)
      REAL X2CONS(MXXVEC,MXXTRAN)
      REAL EOMX2CONS(MXEOMXVEC,MXEOMXTRAN)     !(csi*LE)
      REAL EOML0CONS(MXEOML0VEC,MXEOML0VEC) !(Tbar*RE)
      !REAL EOMCONS(MXXVEC,MXEOMXTRAN)     !(csi*LE).(Tbar*RE)
      REAL WORK(LWORK)
      REAL ZERO, SIGN, EIGVI, EIGVF
      REAL WIAF, WXINIF, WIBF, WLiXi, WL0Rf, WLiXiRf
#else
      DOUBLE PRECISION RESULT(MXOPA)
      DOUBLE PRECISION FCONS(MXFVEC,MXFTRAN)
      DOUBLE PRECISION EACONS(MXAVEC,MXATRAN)
      DOUBLE PRECISION X2CONS(MXXVEC,MXXTRAN)
      DOUBLE PRECISION EOMX2CONS(MXEOMXVEC,MXEOMXTRAN)!(csi*LE)
      DOUBLE PRECISION EOML0CONS(MXEOML0VEC,MXEOML0TRAN) !(Tbar*RE)
      !DOUBLE PRECISION EOMCONS(MXEOMXVEC,MXEOMXTRAN)  !(csi*LE).(Tbar*RE)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, SIGN, EIGVI, EIGVF
      DOUBLE PRECISION WIAF, WXINIF, WIBF, WLiXi, WL0Rf, WLiXiRf
#endif
      PARAMETER (ZERO = 0.0D0)

      CHARACTER LABEL*(8)
      LOGICAL LORX, LPDBS
      INTEGER ITRAN, I, IRSD, IRSDX, ISTATEI, ISTATEF, ISYMI, ISYMF,
     &        ISTISY, ISTFSY, IOP, IOPER, ISYMO, ISYME, ITURN,
     &        IKAP, MXEAVEC, MXE2VEC, IN2VEC, IR1VEC, MFVEC, 
     &        ITMIF, IVEC, NBOPA, IDUM,
     &        IMULI, IMULF

      INTEGER MXEOME2VEC, MXEOM2
      INTEGER NL0TRAN

* external functions:
      INTEGER IR1TAMP
      INTEGER IN2AMP

*---------------------------------------------------------------------*
* initializations:
* initialize for EOM as well ....
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXATRAN
       IEATRAN(1,ITRAN)  = 0
       IEATRAN(2,ITRAN)  = 0
       IEATRAN(3,ITRAN)  = -1
       IEATRAN(4,ITRAN)  = -1
       IEATRAN(5,ITRAN)  = 0
       DO IVEC = 1, MXAVEC
        IEADOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      DO ITRAN = 1, MXXTRAN
       IXE2TRAN(1,ITRAN)  = 0
       IXE2TRAN(2,ITRAN)  = 0
       IXE2TRAN(3,ITRAN)  = -1
       IXE2TRAN(4,ITRAN)  = -1
       IXE2TRAN(5,ITRAN)  = 0
       DO IVEC = 1, MXXVEC
        IX2DOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      DO ITRAN = 1, MXEOMXTRAN
       IEOMXE2TRAN(1,ITRAN)  = 0
       IEOMXE2TRAN(2,ITRAN)  = 0
       IEOMXE2TRAN(3,ITRAN)  = -1
       IEOMXE2TRAN(4,ITRAN)  = -1
       IEOMXE2TRAN(5,ITRAN)  = 0
       DO IVEC = 1, MXEOMXVEC
        IEOMX2DOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      !megaredundant but I am too tired to make it smarter...      
      DO ITRAN = 1, MXEOML0TRAN
         IEOML0TRAN(ITRAN) = 0
         DO IVEC = 1, MXEOML0VEC
            IEOML0DOTS(IVEC,ITRAN) = 0
         END DO
      END DO

      DO ITRAN = 1, MXFTRAN
       DO I = 1, 3
        IFTRAN(I,ITRAN)  = 0
       END DO
       DO IVEC = 1, MXFVEC
        IFDOTS(IVEC,ITRAN)  = 0
       END DO
      END DO

      NFTRAN   = 0
      NXE1TRAN = 0
      NXE2TRAN = 0
      NEOMXE2TRAN = 0
      NEOML0TRAN = 0

      NBOPA   = 0
      MFVEC   = 0
      MXE2VEC = 0
      MXEAVEC = 0

      MXEOME2VEC = 0
      MXEOM2 = 0
      NL0TRAN = 0

! mi manca qualcosa qua sopra...

*---------------------------------------------------------------------*
* start loop over all requested transition moments:
*---------------------------------------------------------------------*
      DO IRSDX  = 1, 2*NXQR2ST
       ITURN = 1 + (IRSDX-1)/NXQR2ST
       IRSD  = IRSDX - (ITURN-1)*NXQR2ST

       IF (ITURN.EQ.1) THEN
         ISTATEI = IQR2ST(IRSD,1)
         ISTATEF = IQR2ST(IRSD,2)
       ELSE IF (ITURN.EQ.2) THEN
         ! switch state indices (and thereby also the sign of the freqs)
         ! to get the conjugated transition moments
         ISTATEI = IQR2ST(IRSD,2)
         ISTATEF = IQR2ST(IRSD,1)
       ELSE
         CALL QUIT('Error in CCXOPA_EOMSETUP')
       END IF

       ISYMI   = ISYEXC(ISTATEI)
       ISYMF   = ISYEXC(ISTATEF)
       ISYME   = MULD2H(ISYMI,ISYMF)
       ISTISY  = ISTATEI - ISYOFE(ISYMI)
       ISTFSY  = ISTATEF - ISYOFE(ISYMF)
       EIGVI   = EIGVAL(ISTATEI)
       EIGVF   = EIGVAL(ISTATEF)
CRF    We need to know multiplicities as well.
       IMULI   = IMULTE(ISTATEI)
       IMULF   = IMULTE(ISTATEF)

       IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'CCXOPA_EOMSETUP:'
         WRITE(LUPRI,*) 'ITURN,IRSD:',ITURN,IRSD
         WRITE(LUPRI,*) 'ISTATEI,ISTATEF:',ISTATEI,ISTATEF
         WRITE(LUPRI,*) 'ISYMI,ISYMF:',ISYMI,ISYMF
         WRITE(LUPRI,*) 'ISTISY,ISTFSY:',ISTISY,ISTFSY
         WRITE(LUPRI,*) 'EIGVI,EIGVF:',EIGVI,EIGVF
         WRITE(LUPRI,*) 'IMULI,IMULF:',IMULI,IMULF
       END IF

       IF (IMULI.NE.IMULF) THEN
          WRITE(LUPRI,*) 'Only singlet operators are currently'
     &                   //' implemented!'
          CYCLE
       END IF

       DO IOP = 1, NQR2OP
        IOPER = IQR2OP(IOP)
        LORX  = .FALSE.
        ISYMO = ISYOPR(IOPER)
        LABEL = LBLOPR(IOPER)
        LPDBS = LPDBSOP(IOPER)
        IKAP  = 0

        IF (LPDBS) CALL QUIT('perturbation-dependent basis sets not '//
     &              'implemented in CCXOPA_EOMSETUP.')

        IF (ISYMO.EQ.ISYME) THEN

          NBOPA = NBOPA + 1

          IF (NBOPA.GT.MXOPA) THEN
             CALL QUIT('NBOPA out of range in CCXOPA_EOMSETUP.')
          END IF

*---------------------------------------------------------------------*
*         in all cases we need LE x A{X} x RE
*---------------------------------------------------------------------*
          !write(lupri,*) "Call CC_SETXE('Eta')"
          !call flshfo(lupri)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXATRAN,MXAVEC,
     &                  ISTATEI,IOPER,IKAP,0,0,0,ISTATEF,ITRAN,IVEC)
          NXE1TRAN = MAX(NXE1TRAN,ITRAN)
          MXEAVEC  = MAX(MXEAVEC, IVEC)
          WIAF     = EACONS(IVEC,ITRAN)

*---------------------------------------------------------------------*
*         add N2 * Xksi{X} or LE * B * RE * R1, depending on QR22N1
*---------------------------------------------------------------------*
          WXINIF  = ZERO
          WIBF    = ZERO
          WLiXi   = ZERO
          WL0Rf   = ZERO
          WLiXiRf = ZERO

          IF (.NOT.CIS) THEN
            IF (QR22N1) THEN
              IN2VEC=IN2AMP(ISTATEI,-EIGVI,ISYMI,ISTATEF,+EIGVF,ISYMF)
              CALL CC_SETXE('Xi ',IXE2TRAN,IX2DOTS,MXXTRAN,MXXVEC,
     &                      0,IOPER,IKAP,0,0,0,IN2VEC,ITRAN,IVEC)
              NXE2TRAN = MAX(NXE2TRAN,ITRAN)
              MXE2VEC  = MAX(MXE2VEC, IVEC)
              WXINIF   = X2CONS(IVEC,ITRAN)
            ELSE IF (LEOMXOPA.AND.(IMULF.EQ.1).AND.(ISYMF.EQ.1)) then
              ! Only do this term if the F state is totally symmetric in
              ! both spin and space!
              !
              ! recover index of LE and do dot prod with Xi(X)
              !
              CALL CC_SETXE('Xi ',IEOMXE2TRAN,IEOMX2DOTS,
     &                            MXEOMXTRAN,MXEOMXVEC,
     &             0,IOPER,IKAP,0,0,0,ISTATEI,ITRAN,IVEC)
              NEOMXE2TRAN = MAX(NEOMXE2TRAN,ITRAN)
              MXEOME2VEC  = MAX(MXEOME2VEC, IVEC)
              WLiXi  = EOMX2CONS(IVEC,ITRAN)
              !write(lupri,*) "WLiXi:", WLiXi
              !
              !now we fill in the list required for the dot products
              !
              CALL CC_SETDOT(IEOML0TRAN,IEOML0DOTS,
     &                       MXEOML0TRAN,MXEOML0VEC,
     &                       0,ISTATEF,ITRAN,IVEC)
              NEOML0TRAN  = MAX(NEOML0TRAN,ITRAN)
              MXEOM2      = MAX(MXEOM2,IVEC)

              WL0Rf    = EOML0CONS(IVEC,ITRAN)
              !write(lupri,*) "WL0Rf:", WL0Rf
              WLiXiRf = WLiXi*WL0Rf
            ELSE IF(.NOT.LEOMXOPA) THEN
              IR1VEC = IR1TAMP(LABEL,LORX,EIGVI-EIGVF,IDUM)
              CALL CC_SETF12(IFTRAN,IFDOTS,MXFTRAN,MXFVEC,
     &                       ISTATEI,ISTATEF,IR1VEC,ITRAN,IVEC)
              NFTRAN = MAX(NFTRAN,ITRAN)
              MFVEC  = MAX(MFVEC, IVEC)
              WIBF   = FCONS(IVEC,ITRAN)
            END IF
            !end if
          END IF

*---------------------------------------------------------------------*
*          add contributions together:
*---------------------------------------------------------------------*
           IF (LADD) THEN

              ITMIF = (NQR2OP*(IRSD-1) + IOP-1)*2 + ITURN

              RESULT(ITMIF) = WIAF + WXINIF + WIBF - WLIXiRF

              IF (LOCDBG) THEN
                WRITE (LUPRI,*) '----- Summary after add ------'
                WRITE (LUPRI,*) 'IDX of result = ',ITMIF
                WRITE (LUPRI,*) 'ISTATEI, EIGVI:',ISTATEI,EIGVI
                WRITE (LUPRI,*) 'ISTATEF, EIGVF:',ISTATEF,EIGVF
                WRITE (LUPRI,*) 'OPERATOR:',LABEL
                WRITE (LUPRI,*) 'L^i A{X} x R^f :',WIAF
                WRITE (LUPRI,*) 'N^if x Xksi{X}:',WXINIF
                WRITE (LUPRI,*) 'L^i x B x R^f x R^X:',WIBF
                WRITE (LUPRI,*) 'L^i x Xksi{X}:',WLiXi
                WRITE (LUPRI,*) 'R^f x L0:',WL0Rf
                WRITE (LUPRI,*) '(L^i x Xksi{X})(R^f x L0):',WLiXiRF
                WRITE (LUPRI,*) 'Total result:',RESULT(ITMIF)
                WRITE (LUPRI,*) '------------------------------'
              END IF

           END IF

*---------------------------------------------------------------------*
*       end loop over transition moments
*---------------------------------------------------------------------*

        END IF
       END DO
      END DO

      IF      (MFVEC.GT.MXFVEC) THEN
         CALL QUIT('MFVEC has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (MXEAVEC.GT.MXAVEC) THEN
         CALL QUIT('MXEAVEC has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (MXE2VEC.GT.MXXVEC) THEN
         CALL QUIT('MXE2VEC has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (MXEOME2VEC.GT.MXEOMXVEC) THEN
         CALL QUIT('MXEOME2VEC has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (MXEOM2.GT.MXEOML0VEC) THEN
         CALL QUIT('MXEOM2 has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (NFTRAN.GT.MXFTRAN) THEN
         CALL QUIT('NFTRAN has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (NXE1TRAN.GT.MXATRAN) THEN
         CALL QUIT('NXE1TRAN has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (NXE2TRAN.GT.MXXTRAN) THEN
         CALL QUIT('NXE2TRAN has been out of bounds CCXOPA_EOMSETUP.')
      ELSE IF (NEOMXE2TRAN.GT.MXEOMXTRAN) THEN
         CALL QUIT('NEOMXE2TRAN has been out bounds CCXOPA_EOMSETUP.')
      ELSE IF (NEOML0TRAN.GT.MXEOML0TRAN) THEN
         CALL QUIT('NEOML0TRAN has been out bounds CCXOPA_EOMSETUP.')
      END IF

*---------------------------------------------------------------------*
* print the lists:
*---------------------------------------------------------------------*
* general statistics:
      IF ((.NOT.LADD) .OR. LOCDBG) THEN
       WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NBOPA,
     &      ' transition moments'
       WRITE(LUPRI,'((8X,A,I3,A))') 
     & ' - ',NFTRAN,  ' F matrix transformations with RE vectors',
     & ' - ',NXE1TRAN,' A{X} matrix transformations with LE vectors',
     & ' - ',NXE2TRAN,' extra XKSI vector calculations ',
     & ' - ',NEOMXE2TRAN,' extra XKSI vector calculations (EOM)',
     & ' - ',NEOML0TRAN,' extra L0 vector (EOM)'
       WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF

      IF (LOCDBG) THEN

         ! F matrix transformations:
         WRITE(LUPRI,*)'List of F matrix transformations:'
         DO ITRAN = 1, NFTRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MFVEC)
         END DO
         WRITE(LUPRI,*)

         ! LE x A{X} vector calculations:
         WRITE(LUPRI,*) 'List of A{O} matrix transformations:'
         DO ITRAN = 1, NXE1TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IEATRAN(I,ITRAN),I=1,5),(IEADOTS(I,ITRAN),I=1,MXEAVEC)
         END DO
         WRITE(LUPRI,*)

         ! extra Xi{O} vector calculations:
         WRITE(LUPRI,*) 'List of extra Xi{O} vector calculations:'
         DO ITRAN = 1, NXE2TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXE2TRAN(I,ITRAN),I=1,5),(IX2DOTS(I,ITRAN),I=1,MXE2VEC)
         END DO
         WRITE(LUPRI,*)

         ! extra Xi{O} vector calculations for EOM:
         WRITE(LUPRI,*) 'List of EOM Xi{O} vector calculations:'
         DO ITRAN = 1, NEOMXE2TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IEOMXE2TRAN(I,ITRAN),I=1,5),
     &             (IEOMX2DOTS(I,ITRAN),I=1,MXEOME2VEC)
         END DO
         WRITE(LUPRI,*)


         !extra L0 x RE vector dot products:
         WRITE (LUPRI,*) 'List of L0 x RE dot products:'
         DO ITRAN = 1, NEOML0TRAN
           WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &     IEOML0TRAN(ITRAN),(IEOML0DOTS(I,ITRAN),I=1,MXEOML0VEC)
         END DO
         WRITE (LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCXOPA_EOMSETUP                      *
*---------------------------------------------------------------------*
