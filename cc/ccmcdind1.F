*=====================================================================*
       SUBROUTINE CC_MCDIND1(WORK,LWORK) 
*---------------------------------------------------------------------*
*  Purpose: Determine which vectors are needed in magnetic circular
*           dichroism calculations 
*           Flags are set for: 2nd-order rhs vectors for T^AB, 
*                              1st-order T^X (w_X) response amplitudes 
*                              M^f(w_f) lagrangian vectors, 
*                              eigenvectors responses E^fX, Ebar^fX
*                              1st order rhs vectors for Tbar^A (eta part)
*                              projected Tbar^A (PL1)
*
*  Written by Sonia Coriani
*  Version: 04/04-2000
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccmcdinf.h>
#include <ccrspprp.h>
#include <ccexcinf.h>
#include <ccexci.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABSOP
      INTEGER IOPA,IOPB,IOPC
      LOGICAL LORXA,LORXB,LORXC, LPDBSA,LPDBSB,LPDBSC, SKIP_IT, LRELAX
      INTEGER ISYMA, ISYMB, ISYMC, ISYMAB, ISYMS_F, ISYMS, ISTATE
      INTEGER IEIGV_F, ISTAT_F, IEXCI_F, INUM, IOPER, IDX, IDXS
      INTEGER ISGNSOP,ISYSOP,NLORX,LWORK
      LOGICAL LPROJ
      
#if defined (SYS_CRAY)
      REAL EIGVA_F, ZERO, WORK(LWORK) 
#else
      DOUBLE PRECISION EIGVA_F, ZERO, WORK(LWORK)
#endif
      PARAMETER ( ZERO = 0.0d0 )

* external functions:
      INTEGER IROPER
      INTEGER IRHSR1
      INTEGER IRHSR2
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILRMAMP
      INTEGER IER1AMP
      INTEGER IEL1AMP
      INTEGER IETA1
      INTEGER IPL1ZETA

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./

*----------------------------------------------------------------------*
* Begin
*----------------------------------------------------------------------*

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_MCDIND> NMCDOPER = ',NMCDOPER
      END IF

*----------------------------------------------------------------------*
* test if operators are available and translate IAMCDOP,IBMCDOP,ICMCDOP
* arrays from the PRPLBL list to the new list maintained by IROPER.
* Note that NMCDOPER is the number of operator-triples (r,L,r)
*----------------------------------------------------------------------*
  
      LPROJ = .FALSE.

      IF (FIRSTCALL) THEN

        IOPER = 1
        DO WHILE (IOPER .LE. NMCDOPER)

          SKIP_IT = .FALSE.
          LABELA = PRPLBL(IAMCDOP(IOPER))
          LABELB = PRPLBL(IBMCDOP(IOPER))
          LABELC = PRPLBL(ICMCDOP(IOPER))
          LORXA  = LAMCDRX(IOPER)
          LORXB  = LBMCDRX(IOPER)
          LORXC  = LCMCDRX(IOPER)
          IOPA   = IROPER(LABELA,ISYMA)
          IOPB   = IROPER(LABELB,ISYMB)
          IOPC   = IROPER(LABELC,ISYMC)

          WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &     'CHECK TRIPLET:',LABELA, LABELB, LABELC
          CALL FLSHFO(LUPRI)


          IF ( (IOPA.LT.0) .OR. (IOPB.LT.0) .OR. (IOPC.LT.0) ) THEN

             WRITE(LUPRI,'(/2X,7A,/2X,2A)')
     &     ' WARNING: ONE OF THE OPERATORS WITH THE LABELS "',
     &       LABELA,'", "', LABELB,'", "', LABELC,'" IS NOT AVAILABLE.',
     &     ' MAGNE.CIRCUL.DICHR. CALCULATION IS CANCELED FOR THIS',
     &     ' OPERATOR TRIPLET.'
 
             SKIP_IT = .TRUE.
          END IF

          NLORX = 0
          IF (LORXA .OR. LPDBSOP(IOPA)) NLORX = NLORX + 1
          IF (LORXB .OR. LPDBSOP(IOPB)) NLORX = NLORX + 1
          IF (LORXC .OR. LPDBSOP(IOPC)) NLORX = NLORX + 1

          IF (NLORX.GT.1) THEN
            WRITE(LUPRI,'(/2X,8A,/2X,A,/2X,A)')
     &       ' WARNING: OPERATOR TRIPLET "',
     &         LABELA,'", "', LABELB,'", "', LABELC,'"',
     &       ' WITH MORE THAN ONE FIELD WHICH',
     &       ' INVOKES ORBITAL RELAXATION OR A PERTUR.-DEP. BASIS SET.',
     &       ' CALCULATION IS CANCELED FOR THIS OPERATOR TRIPLE.'
          END IF
          IF (.NOT. SKIP_IT) THEN
             ! if we have field-dependent basis sets, we need also
             ! to check, if the second-derivative integrals for this
             ! perturbation pair are available
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPB)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (LPDBSOP(IOPA) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (LPDBSOP(IOPB) .OR. LPDBSOP(IOPC)) THEN
                CALL CC_FIND_SO_OP(LABELB,LABELC,LABSOP,ISYSOP,
     &                             ISGNSOP,INUM,WORK,LWORK)
                IF (INUM.LT.0) SKIP_IT = .TRUE.
             END IF
             IF (SKIP_IT) THEN
               WRITE(LUPRI,'(/2X,7A,/2X,A,/2X,A)')
     &          ' WARNING: FOR THE OPERATOR TRIPLET "',
     &            LABELA,'", "', LABELB,'", "', LABELC,'"',
     &         ' A SEC. ORD. OPERATOR IS MISSING.',
     &         ' CALCULATION IS IGNORED.'
             END IF
          END IF



          IF (SKIP_IT) THEN
            DO IDX = IOPER, NMCDOPER-1
              IAMCDOP(IDX) = IAMCDOP(IDX+1)
              IBMCDOP(IDX) = IBMCDOP(IDX+1)
              ICMCDOP(IDX) = ICMCDOP(IDX+1)
              LAMCDRX(IDX) = LAMCDRX(IDX+1)
              LBMCDRX(IDX) = LBMCDRX(IDX+1)
              LCMCDRX(IDX) = LCMCDRX(IDX+1)
            END DO
            NMCDOPER = NMCDOPER - 1        !decrease # of triplets
          ELSE 
            WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &      'PUT TRIPLET:',LABELA, LABELB, LABELC,' ONTO THE LIST.'

            IAMCDOP(IOPER) = IROPER(LABELA,ISYMA)
            IBMCDOP(IOPER) = IROPER(LABELB,ISYMB)
            ICMCDOP(IOPER) = IROPER(LABELC,ISYMC)
            IOPER = IOPER + 1
          END IF

        END DO

        FIRSTCALL = .FALSE.

      END IF                             ! end if (FIRSTCALL)
*
*--------------------------------------------------------------------*
* set/check now symmetries and indices of the excited states
*--------------------------------------------------------------------*
*
      IF (SELMCDST) THEN

* check if all required states available, if not remove them from the list:

        IDXS = 1
        DO WHILE (IDXS .LE. NMCDST)
          IF ( IMCDSTNR(IDXS).GT.NCCEXCI(IMCDSTSY(IDXS),1)) THEN
            WRITE(LUPRI,'(/2X,A,I2,A,I2,A,/2X,A)')
     &       ' WARNING: THE STATE WITH SYMMETRY ',IMCDSTSY(IDXS),
     &       ' AND INDEX ',IMCDSTNR(IDXS) ,
     &       ' IS NOT AVAILABLE.',
     &       ' B TERM CALCULATION IS CANCELED FOR THIS STATE.'
            DO IDX = IDXS, NMCDST-1
              IMCDSTNR(IDX) = IMCDSTNR(IDX+1)     !move next index one step back
              IMCDSTNR(IDX) = IMCDSTNR(IDX+1)
            END DO
            NMCDST = NMCDST - 1
          ELSE
            IDXS = IDXS + 1
          END IF
        END DO
      
      ELSE

* Use default: MCD for all states specified in *CCEXCI

        DO ISYMS = 1, NSYM
          DO ISTATE = 1, NCCEXCI(ISYMS,1)
            IF (NMCDST.LT.MXMCDST) THEN
              NMCDST = NMCDST + 1
              IMCDSTSY(NMCDST) = ISYMS
              IMCDSTNR(NMCDST) = ISTATE
            END IF
          END DO
        END DO


      END IF

*
*--------------------------------------------------------------------*
* set list entries for all the required response vectors:
* NMCDOPER is # of operator triples (A,B,C)
*--------------------------------------------------------------------*
*
      DO 100 IOPER = 1, NMCDOPER
 
        LPROJ = .FALSE.

        LABELA = LBLOPR(IAMCDOP(IOPER))             !get labels back
        LABELB = LBLOPR(IBMCDOP(IOPER))
        LABELC = LBLOPR(ICMCDOP(IOPER))
        
        LPDBSA = LPDBSOP(IAMCDOP(IOPER))
        LPDBSB = LPDBSOP(IBMCDOP(IOPER))
        LPDBSC = LPDBSOP(ICMCDOP(IOPER))

        LORXA  = LAMCDRX(IOPER)
        LORXB  = LBMCDRX(IOPER)
        LORXC  = LCMCDRX(IOPER)

        ISYMA  = ISYOPR(IAMCDOP(IOPER))             !get symmetries back
        ISYMB  = ISYOPR(IBMCDOP(IOPER))
        ISYMC  = ISYOPR(ICMCDOP(IOPER))
        
        ISYMAB = MULD2H(ISYMA,ISYMB)

        LRELAX = LORXA.OR.LORXB.OR.LORXC.OR.LPDBSA.OR.LPDBSB.OR.LPDBSC

        WRITE(LUPRI,'(/2X,A,3(1X,A),A)')
     &   'require responses for op. triplet:',LABELA, LABELB, LABELC
          call flshfo(6)

        IF (ISYMAB.EQ.ISYMC) THEN

           DO 101 IDX = 1, NMCDST
              ISYMS_F = IMCDSTSY(IDX)        !symmetry of excited state
              ISTAT_F = IMCDSTNR(IDX)        !index of exc.state within symmetry
              !absolute index of the exc. state (pointer)
              IEXCI_F = ISYOFE(ISYMS_F) + ISTAT_F
              EIGVA_F = EIGVAL(IEXCI_F)              !excitation energy

              IF (ISYMS_F.EQ.ISYMC) THEN

                 IF (LOCDBG) THEN
                    WRITE (LUPRI,*) 'CC_MCDIND> put onto the list:',
     &               LABELA,'(',-EIGVA_F,'),  ', LABELB,'(',ZERO,'),  ',
     &               IEXCI_F,EIGVA_F
                 END IF

                 INUM = IR1TAMP(LABELA,LORXA,-EIGVA_F,ISYMA)
                 INUM = IR1TAMP(LABELB,LORXB,ZERO,ISYMB)
                 INUM = ILRMAMP(IEXCI_F,EIGVA_F,ISYMC)
                 INUM = IER1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                          LABELA,-EIGVA_F,ISYMA,.FALSE.)
                 INUM = IETA1(LABELB,LORXB,ZERO,ISYMB)
                 IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                 INUM = IEL1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                          LABELB, ZERO,ISYMB,LORXB,LPROJ)
                 IF (.NOT.LUSE2N1) THEN
                   INUM = IR1TAMP(LABELC,LORXC,-EIGVA_F,ISYMC)      
                 END IF
                 INUM = IRHSR1(LABELC,LORXC,EIGVA_F,ISYMC)
                 INUM = IETA1(LABELC,LORXC,EIGVA_F,ISYMC)

                 IF (LUSEPL1) THEN
                    IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                    INUM = IPL1ZETA(LABELA,LORXA,-EIGVA_F,ISYMA,
     &                              LPROJ,IEXCI_F, EIGVA_F,ISYMC)
                 ELSE

                    INUM = IRHSR2(LABELA,LORXA,-EIGVA_F,ISYMA,
     &                            LABELB,LORXB,ZERO,ISYMB) 
                    IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                    INUM = IER1AMP(IEXCI_F,EIGVA_F,ISYMC,
     &                             LABELB,ZERO,ISYMB,LPROJ)
                    INUM = IETA1(LABELA,LORXA,-EIGVA_F,ISYMA)
                 END IF

              END IF
 101       CONTINUE
        END IF
 100  CONTINUE                  
      call flshfo(LUPRI)

      RETURN
      END
*---------------------------------------------------------------------*
