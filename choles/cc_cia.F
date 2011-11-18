C  /* Deck cc_cia */
      SUBROUTINE CC_CIA(XINT,WORK,LWORK,SCLDIA,MEM,NPASS)
C
C     Coupled Cluster Cholesky Integral Assembler.
C     ============================================
C     Written by Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate (ai|bj)-type integrals from MO Cholesky vectors,
C
C              XINT(#ai,#bj) = XINT(#ai,#bj) + sum(J) L1(#ai,J) L2(#bj,J)
C
C              for batches of virtual orbitals. 
C              Note:
C              If the integrals are not initialized locally [INXINT = .FALSE.]
C              and if symmetrization is requested [SYMTRZ = .TRUE.], then the
C              XINT array is assumed symmetric on input.
C
C     Input from argument list:
C
C        WORK  : Scratch array, dimension LWORK.
C        SCLDIA: Scale (ai,bj) integral diagonal with SCLDIA.
C        NPASS : Counter of the number of passes through the Cholesky
C                files to be updated i this routine.
C
C     Input from ciarc.h:
C
C        IOF*1 : Offset to first virtual in each symmetry
C                of the a- and b-batches.
C        LVIR* : Local versions of NVIR array corresponding to
C                the a- and b-batches.
C        NX1AM*: Local versions of NT1AM array corresponding to
C                the a- and b-batches.
C        IX1AM*: Local versions of IT1AM array corresponding to
C                the a- and b-batches.
C        NX2SQ : Total dimension of the integral array XINT.
C                Somewhat redundant, but useful for initial decision
C                regarding algorithm.
C        IX2SQ : Local version of IT2SQ array corresponding to
C                the a- and b-batches.
C        NTOVEC: Total number of Cholesky vectors in each symmetry.
C        ISYCH*: Symmetry of L1 and L2.
C        ITYP* : Type specifiers for reading L1 and L2 from files as
C                specified in CHO_MOP.
C        LIDEN : Logical for L1=L2.
C                LIDEN = .TRUE. : L1 = L2 (redundant, but useful for
C                                 debugging).
C        SYMTRZ: Logical for symmetrization. Ignored if LIDEN = .TRUE.
C                SYMTRZ = .TRUE. : Symmetrize integrals in L1,L2 (equivalent
C                                  to ai,bj symmetrization for full matrix).
C                NOTE: the integral array is assumed symmetric on input for
C                      this option in CC_CIA3!
C        CIAMIO: Logical for minimum I/O. Relevant only for LIDEN = .FALSE.
C                and SYMTRZ = .TRUE. (I.e. in cc_cia4).
C                CIAMIO = .TRUE. : Minimum I/O at the expense of the size
C                                  of the Cholesky batch.
C        GETMNM: Logical for calculation of minimum memory requirement,
C                including memory needed for allocating XINT.
C                GETMNM = .TRUE. : calculate min. mem. requirement and
C                                  do nothing else (min. mem. requirement
C                                  is then returned in MEM). Only information
C                                  needed is ISYCH1, ISYCH2, LIDEN, SYMTRZ,
C                                  and CIAMIO.
C        INXINT: Logical for initialization of XINT array (to zero).
C        IPRCIA: Local print level.
C
C     Output (through argument list):
C
C        XINT  : The NX2SQ integrals.
C        MEM   : The memory (in words) actually used to assemble the integrals.
C        NPASS : Counter of passes through the Cholesky file(s).
C
C******* NOTICE *********************** NOTICE *************** NOTICE **********
C
C     If any changes are made to the allocations in CC_CIA routines,
C     please remember to make the corresponding changes in CC_CIAMNM.
C
C******* NOTICE *********************** NOTICE *************** NOTICE **********
C
C     TODO: The symmetrization does not take advantage of overlapping a- and
C           b-distributions, but calculates the integrals twice.....
C           Potentially, as the integral assembly is the most expensive step,
C           a relatively large amount of CPU time (as well as I/O) may be
C           saved by fixing this.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*6 SECNAM, STRAST
      PARAMETER (SECNAM = 'CC_CIA')
      PARAMETER (STRAST = '******')

      LOGICAL LOCDBG, FULLIN
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (TINY = 1.0D-14)
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.00D0)

C     If requested, calculate absolute min. mem. requirement.
C     That is, the minimum memory needed if there are one a
C     and one b for which the integrals are to be assembled.
C     Return without doing anything else...
C     ------------------------------------------------------

      IF (GETMNM) THEN
         CALL CC_CIAMNM(MEM)
         RETURN
      ENDIF

C     Initialize internal commons.
C     ----------------------------

      TIMCIA(1) = SECOND()
      DO ITMCIA = 2,NTMCIA
         TIMCIA(ITMCIA) = ZERO
      ENDDO

      DO ISYM = 1,NSYM
         NCHBAT(ISYM) = 0
      ENDDO

C     Debug: print input parameters.
C     ------------------------------

      IF (LOCDBG) THEN
         CALL HEADER('Debug: Input Parameters in '//SECNAM,-1)
         CALL CC_CIAIPR(XINT,LWORK,SCLDIA)
      ENDIF

C     Input tests.
C     ------------

      ISYMX = MULD2H(ISYCH1,ISYCH2)
      XX2SQ = ONE*NX2SQ

      IF (NX2SQ .EQ. 0) THEN
         WRITE(LUPRI,'(//,10X,A,A,/,10X,A,A,/,10X,A,I10,/,10X,A,A,//)')
     &   '************************************',STRAST,
     &   'WARNING: No integrals calculated in ',SECNAM,
     &   'Integral dimension from input:  ',NX2SQ,
     &   '************************************',STRAST
         RETURN
      ELSE IF ((XX2SQ.LT.ZERO) .OR. (XX2SQ.GT.XT2SQ(ISYMX))) THEN
         WRITE(LUPRI,'(//,10X,A,A,A)')
     &   'Error detected in ',SECNAM,':'
         WRITE(LUPRI,'(10X,A,5X,I10,/,10X,A,F15.1,/,10X,A,14X,I1)')
     &   'Integral dimension from input          : ',NX2SQ,
     &   'Full integral dimension from dccsdsym.h: ',XT2SQ(ISYMX),
     &   'Symmetry of integrals                  : ',ISYMX
         WRITE(LUPRI,'(10X,A,I1,/,10X,A,I1,/)')
     &   'Symmetry of first  Cholesky vectors: ',ISYCH1,
     &   'Symmetry of second Cholesky vectors: ',ISYCH2
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      IF (LIDEN .AND. (ISYMX.NE.1)) THEN
         WRITE(LUPRI,'(//,10X,A,A,A)')
     &   'Error detected in ',SECNAM,':'
         WRITE(LUPRI,'(10X,A,I1,A,/,10X,A)')
     &   'Integrals are non-tot. symmetric, ISYMX = ',ISYMX,',',
     &   'but Cholesky vectors are specified as identical!'
         WRITE(LUPRI,'(10X,A,I1,/,10X,A,I1,/)')
     &   'Symmetry of first  Cholesky vectors: ',ISYCH1,
     &   'Symmetry of second Cholesky vectors: ',ISYCH2
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Save NPASS.
C     -----------

      NPSIN = NPASS

C     Initialize integral array if requested.
C     ---------------------------------------

      IF (INXINT) THEN
         CALL DZERO(XINT,NX2SQ)
      ENDIF

C     Calculate initial norm.
C     -----------------------

      IF (IPRCIA .GE. IPRLVL) THEN
         XSTART = DSQRT(DDOT(NX2SQ,XINT,1,XINT,1))
      ENDIF

C     Is this the full integral matrix?
C     ---------------------------------

      DIFF   = XT2SQ(ISYMX) - XX2SQ
      FULLIN = (ABS(DIFF) .LE. TINY)

C     Call appropriate routine to do the job.
C     ---------------------------------------

      IF (LIDEN) THEN
         IF (FULLIN) THEN
            CALL CC_CIA1(XINT,WORK,LWORK,MEM,NPASS)
         ELSE
            CALL CC_CIA2(XINT,WORK,LWORK,MEM,NPASS)
         ENDIF
      ELSE
         IF (FULLIN) THEN
            IF (.NOT. INXINT) THEN
               CALL DSCAL(NX2SQ,HALF,XINT,1)
            ENDIF
            CALL CC_CIA3(XINT,WORK,LWORK,MEM,NPASS)
         ELSE
            CALL CC_CIA4(XINT,WORK,LWORK,MEM,NPASS)
         ENDIF
      ENDIF

C     Scale diagonal (if meaningful).
C     -------------------------------

      IF (SCLDIA .NE. ONE) THEN
         DTIME = SECOND()
         CALL CC_CIADSCL(XINT,ISYMX,SCLDIA,FULLIN)
         DTIME = SECOND() - DTIME
         TIMCIA(6) = TIMCIA(6) + DTIME
      ENDIF

C     Print summary.
C     --------------

      IF (IPRCIA .GE. IPRLVL) THEN

         XEND = DSQRT(DDOT(NX2SQ,XINT,1,XINT,1))

         TIMCIA(1) = SECOND() - TIMCIA(1)

         MNU = LWORK - MEM
         XLW = ONE*LWORK
         XMM = ONE*MEM
         XMN = ONE*MNU
         XMB = XLW*WTOMB
         XUS = D100*XMM/XLW
         XNU = D100*XMN/XLW

         CALL AROUND('Output from '//SECNAM)

         WRITE(LUPRI,'(5X,A,A,/,5X,A)')
     &   'Integrals were assembled in subroutine ',CIAROU,
     &   'Characteristics:'
         IF (INXINT) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      '- integral array was initialized by ',SECNAM
         ELSE
            WRITE(LUPRI,'(5X,A)')
     &      '- integral array was initialized by calling routine'
         ENDIF
         IF (LIDEN) THEN
            IF (FULLIN) THEN
               WRITE(LUPRI,'(5X,A)')
     &         '- full square integrals calculated with L1 = L2'
            ELSE
               WRITE(LUPRI,'(5X,A)')
     &         '- integral batch calculated with L1 = L2'
            ENDIF
            WRITE(LUPRI,'(5X,A)')
     &      '- integrals inherently symmetric in ai,bj'
            WRITE(LUPRI,'(5X,A,I1)')
     &      '- symmetry of L1 = L2: ',ISYCH1
            WRITE(LUPRI,'(5X,A,I4)')
     &      '- type specifier for L1 = L2: ',ITYP1
         ELSE
            IF (FULLIN) THEN
               WRITE(LUPRI,'(5X,A)')
     &         '- full square integrals calculated with L1 != L2'
               IF (SYMTRZ) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '- integrals symmetrized in ai,bj'
               ELSE
                  WRITE(LUPRI,'(5X,A)')
     &            '- integrals not symmetrized in ai,bj'
               ENDIF
            ELSE
               WRITE(LUPRI,'(5X,A)')
     &         '- integral batch calculated with L1 != L2'
               IF (SYMTRZ) THEN
                  WRITE(LUPRI,'(5X,A)')
     &            '- integrals symmetrized in L1,L2'
               ELSE
                  WRITE(LUPRI,'(5X,A)')
     &            '- integrals not symmetrized in L1,L2'
               ENDIF
            ENDIF
            WRITE(LUPRI,'(5X,A,I1)')
     &      '- symmetry of L1: ',ISYCH1
            WRITE(LUPRI,'(5X,A,I1)')
     &      '- symmetry of L2: ',ISYCH2
            WRITE(LUPRI,'(5X,A,I4)')
     &      '- type specifier for L1: ',ITYP1
            WRITE(LUPRI,'(5X,A,I4)')
     &      '- type specifier for L2: ',ITYP2
         ENDIF
         IF ((SCLDIA.NE.ONE) .AND. (ISYMX.EQ.1)) THEN
            WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &      '- ai,bj diagonal scaled with ',SCLDIA
         ENDIF
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Initial norm of integral array: ',XSTART,
     &   'Final   norm of integral array: ',XEND

         CALL HEADER('Memory Usage',-1)
         WRITE(LUPRI,'(5X,A,I10,A,F7.2,A,/,5X,A,I10,A,F7.2,A)')
     &   'Memory used     : ',MEM,' words (',XUS,' % )',
     &   'Memory not used : ',MNU,' words (',XNU,' % )'
         WRITE(LUPRI,'(5X,A,I10,A,F7.1,A)')
     &   'Memory available: ',LWORK,' words (',XMB,' Mb)'

         CALL HEADER('Cholesky Batch Info',-1)
         IPASS = NPASS - NPSIN
         IF (IPASS .GT. 1) THEN
            WRITE(LUPRI,'(5X,I4,A,A,A)')
     &      IPASS,' passes through the Cholesky files in this ',
     &      SECNAM,' call.'
         ENDIF
         WRITE(LUPRI,'(5X,A,I5,A,/)')
     &   'This is info for Cholesky file pass number',NPASS,':'
         WRITE(LUPRI,'(5X,A,/,5X,A)')
     &   'Sym. J   #Vectors   #Batches',
     &   '----------------------------'
         DO ISYCHO = 1,NSYM
            WRITE(LUPRI,'(8X,I1,3X,I10,1X,I10)')
     &      ISYCHO,NTOVEC(ISYCHO),NCHBAT(ISYCHO)
         ENDDO
         WRITE(LUPRI,'(5X,A)')
     &   '----------------------------'

         CALL HEADER('Timing',-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors     : ',TIMCIA(2),
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting Cholesky vectors    : ',TIMCIA(3),
     &   ' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for contracting Cholesky vectors: ',TIMCIA(4),
     &   ' seconds'
         IF (TIMCIA(5) .GT. ZERO) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for symmetrizing integrals      : ',TIMCIA(5),
     &      ' seconds'
         ENDIF
         IF (SCLDIA .NE. ONE) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for scaling integral diagonal   : ',TIMCIA(6),
     &      ' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A)')
     &  '--------------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                        : ',TIMCIA(1),
     &   ' seconds'

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
C  /* Deck cc_ciadscl */
      SUBROUTINE CC_CIADSCL(XINT,ISYMX,SCLDIA,FULLIN)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Scale diagonal of XINT(#ai,#bj).
C
#include <implicit.h>
      DIMENSION XINT(*)
      LOGICAL FULLIN
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>
#include <priunit.h>

      INTEGER A1, B1, AI, BJ

C     Check that XINT is tot. sym.
C     ----------------------------

      IF (ISYMX .NE. 1) RETURN

C     Scale.
C     ------

      IF (FULLIN) THEN

         DO ISYMAI = 1,NSYM
            DO AI = 1,NT1AM(ISYMAI)
               KAIAI = IT2SQ(ISYMAI,ISYMAI) + NT1AM(ISYMAI)*(AI - 1)
     &               + AI
               XINT(KAIAI) = SCLDIA*XINT(KAIAI)
            ENDDO
         ENDDO

      ELSE

         DO ISYMB = 1,NSYM

            ISYMA = ISYMB

            NA = LVIRA(ISYMA)
            NB = LVIRB(ISYMB)

            IF ((NA.GT.0) .AND. (NB.GT.0)) THEN

               LA1 = IOFA1(ISYMA)
               LA2 = LA1 + NA - 1

               LB1 = IOFB1(ISYMB)
               LB2 = LB1 + NB - 1

               IB1 = MAX(LB1,LA1)
               IB2 = MIN(LB2,LA2)
               NIB = IB2 - IB1 + 1

               IF (NIB .GT. 0) THEN

                  A1 = IB1 - LA1 + 1
                  B1 = IB1 - LB1 + 1

                  DO ISYMJ = 1,NSYM

                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = ISYMBJ
                     ISYMI  = ISYMJ

                     NJ = NRHF(ISYMJ)

                     DO J = 1,NJ

                        I = J

                        DO IB = 1,NIB

                           B = B1 + IB - 1
                           A = A1 + IB - 1

                           BJ = IX1AMB(ISYMB,ISYMJ)
     &                        + NB*(J - 1) + B
                           AI = IX1AMA(ISYMA,ISYMI)
     &                        + NA*(I - 1) + A

                           KDIA = IX2SQ(ISYMAI,ISYMBJ)
     &                          + NX1AMA(ISYMAI)*(BJ - 1) + AI

                           XINT(KDIA) = SCLDIA*XINT(KDIA)

                        ENDDO

                     ENDDO

                  ENDDO

               ENDIF

            ENDIF

         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_ciaipr */
      SUBROUTINE CC_CIAIPR(XINT,LWORK,SCLDIA)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Print input information for cc_cia.
C
#include <implicit.h>
      DIMENSION XINT(*)
#include <priunit.h>
#include <ccorb.h>
#include <ciarc.h>

      PARAMETER (ONE = 1.0D0)

      ISYMX = MULD2H(ISYCH1,ISYCH2)

      WRITE(LUPRI,'(A,/,A,/)')
     & 'Information about integral assembly:',
     & '------------------------------------'
      WRITE(LUPRI,'(A,I1,1X,I10)')
     & 'Symmetry and dimension of integral batch: ',ISYMX,NX2SQ
      IF (INXINT) THEN
         WRITE(LUPRI,'(A)')
     &   'Integrals will be initialized'
      ELSE
         WRITE(LUPRI,'(A)')
     &   'Integrals will not be initialized'
      ENDIF
      IF (.NOT. LIDEN) THEN
         IF (SYMTRZ) THEN
            WRITE(LUPRI,'(A)')
     &      'Integrals will be symmetrized in L1,L2.'
         ELSE
            WRITE(LUPRI,'(A)')
     &      'Integrals will not be symmetrized in L1,L2.'
         ENDIF
      ENDIF

      IF (SCLDIA .NE. ONE) THEN
         WRITE(LUPRI,'(A,1P,D15.6)')
     &    'Diagonal scale factor: ',SCLDIA
      ELSE
         WRITE(LUPRI,'(A,1P,D15.6,A)')
     &    'Diagonal scale factor: ',SCLDIA,' (ignored)'
      ENDIF

      WRITE(LUPRI,'(/,A,/,A,/)')
     & 'Information about Cholesky vectors:',
     & '-----------------------------------'
      IF (LIDEN) THEN
         WRITE(LUPRI,'(A)')
     &   'First and second Cholesky vectors are identical.'
         WRITE(LUPRI,'(A,I4)')
     &   'Type specifier for L1 Cholesky vectors: ',ITYP1
         WRITE(LUPRI,'(A,I4)')
     &   'Symmetry of Cholesky vectors          : ',ISYCH1
      ELSE
         WRITE(LUPRI,'(A)')
     &   'First and second Cholesky vectors are different.'
         WRITE(LUPRI,'(A,I4)')
     &   'Type specifier for L1 Cholesky vectors: ',ITYP1
         WRITE(LUPRI,'(A,I4)')
     &   'Type specifier for L2 Cholesky vectors: ',ITYP2
         WRITE(LUPRI,'(A,I4,/,A,I4)')
     &   'Symmetry of L1 Cholesky vectors: ',ISYCH1,
     &   'Symmetry of L2 Cholesky vectors: ',ISYCH2
      ENDIF
      WRITE(LUPRI,*)
      WRITE(LUPRI,1) 'NTOVEC',(NTOVEC(I),I=1,NSYM)

      WRITE(LUPRI,'(/,A,/,A,/)')
     & 'Index arrays for a- and b-batches:',
     & '----------------------------------'
      WRITE(LUPRI,1) 'IOFA1 ',(IOFA1(I),I=1,NSYM)
      WRITE(LUPRI,1) 'IOFB1 ',(IOFB1(I),I=1,NSYM)
      WRITE(LUPRI,*)
      WRITE(LUPRI,1) 'LVIRA ',(LVIRA(I),I=1,NSYM)
      WRITE(LUPRI,1) 'LVIRB ',(LVIRB(I),I=1,NSYM)
      WRITE(LUPRI,*)
      WRITE(LUPRI,1) 'NX1AMA',(NX1AMA(I),I=1,NSYM)
      WRITE(LUPRI,1) 'NX1AMB',(NX1AMB(I),I=1,NSYM)
      WRITE(LUPRI,*)
      DO I = 1,NSYM
         WRITE(LUPRI,1) 'IX1AMA',(IX1AMA(I,J),J=1,NSYM)
      ENDDO
      WRITE(LUPRI,*)
      DO I = 1,NSYM
         WRITE(LUPRI,1) 'IX1AMB',(IX1AMB(I,J),J=1,NSYM)
      ENDDO
      WRITE(LUPRI,*)
      DO I = 1,NSYM
         WRITE(LUPRI,1) 'IX2SQ ',(IX2SQ(I,J),J=1,NSYM)
      ENDDO
      WRITE(LUPRI,*)

      RETURN
    1 FORMAT(A6,': ',8I9)
      END
C  /* Deck cc_cia1 */
      SUBROUTINE CC_CIA1(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(ai,bj) = XINT(ai,bj) + sum(J) L(ai,J) * L(bj,J)
C
C              for full square XINT.
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_CIA1')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMC = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         IF (NTOVEC(ISYCHO) .GT. 0) THEN
            NEED = MAX(NEED,NT1AM(ISYMAI))
         ENDIF
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory stuff.
C     ------------------------

      MEM = 0

C     Loop over Cholesky symmetries.
C     ------------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = ISYMAI

C        Skip this symmetry if there's nothing to do.
C        --------------------------------------------

         IF (NTOVEC(ISYCHO) .LE. 0) GO TO 999
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 999

C        Open Cholesky file.
C        -------------------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

C        Set up batch.
C        -------------

         NVEC = MIN(LWORK/NT1AM(ISYMAI),NTOVEC(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai (= symmetry of bj): ',ISYMAI
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWORK,
     &      'Minimum required (pref. multiple): ',NT1AM(ISYMAI)
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NTOVEC(ISYCHO) - 1)/NVEC + 1

C        Store batch info for later print.
C        ---------------------------------

         NCHBAT(ISYCHO) = NBATCH

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NTOVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = 1
            KBLST = KCHAI + NT1AM(ISYMAI)*NUMV - 1
            MEM   = MAX(MEM,KBLST)

C           Read MO Cholesky vectors.
C           -------------------------

            DTIME = SECOND()
            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYMAI),NUMV,
     &                      JVEC1,LUCHO1)
            DTIME = SECOND() - DTIME
            TIMR  = TIMR     + DTIME

C           Calculate contribution to integrals.
C           ------------------------------------

            DTIME = SECOND()

            NTOAI = MAX(NT1AM(ISYMAI),1)
            NTOBJ = NTOAI

            KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + 1

            CALL DGEMM('N','T',NT1AM(ISYMAI),NT1AM(ISYMBJ),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHAI),NTOBJ,
     &                 ONE,XINT(KOFF1),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close Cholesky file.
C        --------------------

         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

C        Escape point for empty symmetry.
C        --------------------------------

  999    CONTINUE

      ENDDO

C     Update global timing.
C     ---------------------

      TIMCIA(2) = TIMCIA(2) + TIMR
      TIMCIA(4) = TIMCIA(4) + TIMC

C     Update Cholesky pass counter.
C     -----------------------------

      NPASS = NPASS + 1

C     Identify assembler routine.
C     ---------------------------

      CIAROU = SECNAM

      RETURN
      END
C  /* Deck cc_cia2 */
      SUBROUTINE CC_CIA2(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(#ai,#bj) = XINT(#ai,#bj) + sum(J) L(#ai,J) * L(#bj,J)
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_CIA2')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO

C     Initial test of memory.
C     -----------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = ISYMAI
         IF ((NX1AMA(ISYMAI).GT.0) .AND. (NX1AMB(ISYMBJ).GT.0)) THEN
            NEEDS = NT1AM(ISYMAI) + NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
            NEED  = MAX(NEED,NEEDS)
         ENDIF
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory stuff.
C     ------------------------

      MEM = 0

C     Loop over Cholesky symmetries.
C     ------------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = ISYMAI

         IF (NTOVEC(ISYMAI) .LE. 0) GO TO 999
         IF (NX1AMA(ISYMAI) .LE. 0) GO TO 999
         IF (NX1AMB(ISYMBJ) .LE. 0) GO TO 999

C        Open Cholesky file.
C        -------------------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

C        Allocation for 1 full Cholesky vector.
C        --------------------------------------

         KCHOL = 1
         KEND1 = KCHOL + NT1AM(ISYMAI)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            MINADD = NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
            NEED   = KEND1 - 1 + MINADD
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry         : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai (= sym. bj): ',ISYMAI
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required   : ',NEED,
     &      'Total memory available    : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
         NVEC = MIN(LWRK1/MINMEM,NTOVEC(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai (= symmetry of bj): ',ISYMAI
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWRK1,
     &      'Minimum required (pref. multiple): ',MINMEM
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NTOVEC(ISYCHO) - 1)/NVEC + 1

C        Store batch info for later print.
C        ---------------------------------

         NCHBAT(ISYCHO) = NBATCH

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            TIBAT = SECOND()

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NTOVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = KEND1
            KCHBJ = KCHAI + NX1AMA(ISYMAI)*NUMV
            KBLST = KCHBJ + NX1AMB(ISYMBJ)*NUMV - 1
            MEM   = MAX(MEM,KBLST)

C           Read MO Cholesky vectors and copy out subblocks.
C           ------------------------------------------------

            DO IVEC = 1,NUMV

               DTIME = SECOND()
               JVEC  = JVEC1 + IVEC - 1
               CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYMAI),1,
     &                         JVEC,LUCHO1)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

               DTIME = SECOND()
               CALL CC_CIACPS(WORK(KCHOL),WORK(KCHAI),
     &                        IVEC,IOFA1,LVIRA,NX1AMA,
     &                        IX1AMA,ISYMAI)
               CALL CC_CIACPS(WORK(KCHOL),WORK(KCHBJ),
     &                        IVEC,IOFB1,LVIRB,NX1AMB,
     &                        IX1AMB,ISYMBJ)
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

            ENDDO

C           Calculate integral contribution.
C           --------------------------------

            DTIME = SECOND()

            NTOAI = MAX(NX1AMA(ISYMAI),1)
            NTOBJ = MAX(NX1AMB(ISYMBJ),1)

            KOFFX = IX2SQ(ISYMAI,ISYMBJ) + 1

            CALL DGEMM('N','T',NX1AMA(ISYMAI),NX1AMB(ISYMBJ),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                 ONE,XINT(KOFFX),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close Cholesky file.
C        --------------------

         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

  999    CONTINUE

      ENDDO

C     Update global timings.
C     ----------------------

      TIMCIA(2) = TIMCIA(2) + TIMR
      TIMCIA(3) = TIMCIA(3) + TIMS
      TIMCIA(4) = TIMCIA(4) + TIMC

C     Update Cholesky pass counter.
C     -----------------------------

      NPASS = NPASS + 1

C     Identify assembler routine.
C     ---------------------------

      CIAROU = SECNAM

      RETURN
      END
C  */ Deck cc_cia3 */
      SUBROUTINE CC_CIA3(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(ai,bj) = XINT(ai,bj) + sum(J) L1(ai,J) * L2(bj,J)
C
C              for full square XINT.
C
C     NB!!! The final result is symmetrized in (ai,bj) indices
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_CIA3')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMC = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = MULD2H(ISYCHO,ISYCH2)
         IF ((NT1AM(ISYMAI).GT.0) .AND. (NT1AM(ISYMBJ).GT.0) .AND.
     &       (NTOVEC(ISYCHO) .GT. 0)) THEN
            NEEDS = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
            NEED  = MAX(NEED,NEEDS)
         ENDIF
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory stuff.
C     ------------------------

      MEM = 0

C     Loop over Cholesky symmetries.
C     ------------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = MULD2H(ISYCHO,ISYCH2)

C        Skip if nothing to do.
C        ----------------------

         IF (NTOVEC(ISYCHO) .LE. 0) GO TO 999
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 999
         IF (NT1AM(ISYMBJ)  .LE. 0) GO TO 999

C        Open Cholesky files.
C        --------------------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)
         CALL CHO_MOP(-1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
         NVEC   = MIN(LWORK/MINMEM,NTOVEC(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai                   : ',ISYMAI
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of bj                   : ',ISYMBJ
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWORK,
     &      'Minimum required (pref. multiple): ',MINMEM
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NTOVEC(ISYCHO) - 1)/NVEC + 1

C        Store batch info for later print.
C        ---------------------------------

         NCHBAT(ISYCHO) = NBATCH

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NTOVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = 1
            KCHBJ = KCHAI + NT1AM(ISYMAI)*NUMV
            KBLST = KCHBJ + NT1AM(ISYMBJ)*NUMV - 1
            MEM   = MAX(MEM,KBLST)

C           Read MO Cholesky vectors.
C           -------------------------

            DTIME = SECOND()
            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYMAI),NUMV,
     &                      JVEC1,LUCHO1)
            CALL CHO_MOREAD(WORK(KCHBJ),NT1AM(ISYMBJ),NUMV,
     &                      JVEC1,LUCHO2)
            DTIME = SECOND() - DTIME
            TIMR  = TIMR     + DTIME

C           Calculate contribution to integrals.
C           ------------------------------------

            DTIME = SECOND()

            NTOAI = MAX(NT1AM(ISYMAI),1)
            NTOBJ = MAX(NT1AM(ISYMBJ),1)

            KOFFX = IT2SQ(ISYMAI,ISYMBJ) + 1

            CALL DGEMM('N','T',NT1AM(ISYMAI),NT1AM(ISYMBJ),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                 ONE,XINT(KOFFX),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close Cholesky files.
C        ---------------------

         CALL CHO_MOP(1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)
         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

  999    CONTINUE

      ENDDO

C     Symmetrize the integrals if requested.
C     --------------------------------------

      IF (SYMTRZ) THEN
         TIMM = SECOND()
         ISYMX = MULD2H(ISYCH1,ISYCH2)
         CALL CC_CIASMT(XINT,ISYMX)
         TIMM = SECOND() - TIMM
      ENDIF

C     Update global timings.
C     ----------------------

      TIMCIA(2) = TIMCIA(2) + TIMR
      TIMCIA(4) = TIMCIA(4) + TIMC
      IF (SYMTRZ) TIMCIA(5) = TIMCIA(5) + TIMM

C     Update Cholesky pass counter.
C     -----------------------------

      NPASS = NPASS + 2

C     Identify assembler routine.
C     ---------------------------

      CIAROU = SECNAM

      RETURN
      END
C  */ Deck cc_cia4 */
      SUBROUTINE CC_CIA4(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(#ai,#bj) = XINT(#ai,#bj) + sum(J) L1(#ai,J) * L2(#bj,J)
C
C              with or without symmetrization in L1,L2 [specified through
C              SYMTRZ].
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_CIA4')

      PARAMETER (TINY = 1.0D-14)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Call appropriate assembler.
C     SYMTRZ: logical for symmetrization.
C     CIAMIO: logical for minimum I/O [only relevant for symmetrizing].
C     -----------------------------------------------------------------

      IF (SYMTRZ) THEN
         IF (CIAMIO) THEN
            CALL CC_CIA4_1(XINT,WORK,LWORK,MEM,NPASS)
         ELSE
            CALL CC_CIA4_2(XINT,WORK,LWORK,MEM,NPASS)
         ENDIF
      ELSE
         CALL CC_CIA4_3(XINT,WORK,LWORK,MEM,NPASS)
      ENDIF

      RETURN
      END
C  /* Deck cc_cia4_1 */
      SUBROUTINE CC_CIA4_1(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(#ai,#bj) = XINT(#ai,#bj)
C                            + P(L1,L2) sum(J) L1(#ai,J) * L2(#bj,J)
C
C              where P(L1,L2) symmetrizes with respect to L1,L2 according to
C
C              P(L1,L2) M(L1,L2) = M(L1,L2) + M(L2,L1)
C
C     The symmetrization is done with minimal I/O (i.e. 1 pass through
C     Cholesky files) at the expense of the batch size.
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CIA4_1')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         ISYAI1 = MULD2H(ISYCHO,ISYCH1)
         ISYBJ2 = MULD2H(ISYCHO,ISYCH2)
         ISYAI2 = ISYBJ2
         ISYBJ1 = ISYAI1
         MINMEM = NX1AMA(ISYAI1) + NX1AMA(ISYAI2)
     &          + NX1AMB(ISYBJ2) + NX1AMB(ISYBJ1)
         IF ((MINMEM.GT.0) .AND. (NTOVEC(ISYCHO).GT.0)) THEN
            NEEDS = NT1AM(ISYAI1) + NT1AM(ISYBJ2) + MINMEM
            NEED  = MAX(NEED,NEEDS)
         ENDIF
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory stuff.
C     ------------------------

      MEM = 0

C     Start loop over Cholesky symmetries (i.e. sym. of J-index).
C     -----------------------------------------------------------

      DO ISYCHO = 1,NSYM

         ISYAI1 = MULD2H(ISYCHO,ISYCH1)
         ISYBJ2 = MULD2H(ISYCHO,ISYCH2)
         ISYAI2 = ISYBJ2
         ISYBJ1 = ISYAI1

C        Skip if nothing to do.
C        ----------------------

         MINMEM = NX1AMA(ISYAI1) + NX1AMA(ISYAI2)
     &          + NX1AMB(ISYBJ2) + NX1AMB(ISYBJ1)

         IF (NTOVEC(ISYCHO) .LE. 0) GO TO 999
         IF (MINMEM         .LE. 0) GO TO 999

C        Open Cholesky files.
C        --------------------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)
         CALL CHO_MOP(-1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)

C        Allocation for 1 full L1 and L2.
C        --------------------------------

         KCHO1 = 1
         KCHO2 = KCHO1 + NT1AM(ISYAI1)
         KEND1 = KCHO2 + NT1AM(ISYBJ2)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            NEED   = KEND1 - 1 + MINMEM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry      : ',ISYCHO
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',NEED,
     &      'Total memory available : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         NVEC = MIN(LWRK1/MINMEM,NTOVEC(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWRK1,
     &      'Minimum required (pref. multiple): ',MINMEM
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NTOVEC(ISYCHO) - 1)/NVEC + 1

C        Store batch info for later print.
C        ---------------------------------

         NCHBAT(ISYCHO) = NBATCH

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NTOVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI1 = KEND1
            KCHBJ2 = KCHAI1 + NX1AMA(ISYAI1)*NUMV
            KCHAI2 = KCHBJ2 + NX1AMB(ISYBJ2)*NUMV
            KCHBJ1 = KCHAI2 + NX1AMA(ISYAI2)*NUMV
            KBLST  = KCHBJ1 + NX1AMB(ISYBJ1)*NUMV - 1
            MEM    = MAX(MEM,KBLST)

C           Read MO Cholesky vectors and copy out subblocks.
C           ------------------------------------------------

            DO IVEC = 1,NUMV

               DTIME = SECOND()
               JVEC = JVEC1 + IVEC - 1
               CALL CHO_MOREAD(WORK(KCHO1),NT1AM(ISYAI1),1,
     &                         JVEC,LUCHO1)
               CALL CHO_MOREAD(WORK(KCHO2),NT1AM(ISYBJ2),1,
     &                         JVEC,LUCHO2)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

               DTIME = SECOND()
               CALL CC_CIACPS(WORK(KCHO1),WORK(KCHAI1),
     &                        IVEC,IOFA1,LVIRA,NX1AMA,
     &                        IX1AMA,ISYAI1)
               CALL CC_CIACPS(WORK(KCHO1),WORK(KCHBJ1),
     &                        IVEC,IOFB1,LVIRB,NX1AMB,
     &                        IX1AMB,ISYBJ1)
               CALL CC_CIACPS(WORK(KCHO2),WORK(KCHBJ2),
     &                        IVEC,IOFB1,LVIRB,NX1AMB,
     &                        IX1AMB,ISYBJ2)
               CALL CC_CIACPS(WORK(KCHO2),WORK(KCHAI2),
     &                        IVEC,IOFA1,LVIRA,NX1AMA,
     &                        IX1AMA,ISYAI2)
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

            ENDDO

C           Calculate integral contributions.
C           ---------------------------------

            DTIME = SECOND()

            NTOAI1 = MAX(NX1AMA(ISYAI1),1)
            NTOBJ2 = MAX(NX1AMB(ISYBJ2),1)
            KOFFX  = IX2SQ(ISYAI1,ISYBJ2) + 1
            CALL DGEMM('N','T',NX1AMA(ISYAI1),NX1AMB(ISYBJ2),NUMV,
     &                 ONE,WORK(KCHAI1),NTOAI1,WORK(KCHBJ2),NTOBJ2,
     &                 ONE,XINT(KOFFX),NTOAI1)

            NTOAI2 = MAX(NX1AMA(ISYAI2),1)
            NTOBJ1 = MAX(NX1AMB(ISYBJ1),1)
            KOFFX  = IX2SQ(ISYAI2,ISYBJ1) + 1
            CALL DGEMM('N','T',NX1AMA(ISYAI2),NX1AMB(ISYBJ1),NUMV,
     &                 ONE,WORK(KCHAI2),NTOAI2,WORK(KCHBJ1),NTOBJ1,
     &                 ONE,XINT(KOFFX),NTOAI2)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close Cholesky files.
C        ---------------------

         CALL CHO_MOP(1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)
         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

  999    CONTINUE

      ENDDO

C     Update timing info.
C     -------------------

      TIMCIA(2) = TIMCIA(2) + TIMR
      TIMCIA(3) = TIMCIA(3) + TIMS
      TIMCIA(4) = TIMCIA(4) + TIMC

C     Update Cholesky pass counter.
C     -----------------------------

      NPASS = NPASS + 2

C     Identify assembler routine.
C     ---------------------------

      CIAROU = SECNAM

      RETURN
      END
C  /* Deck cc_cia4_2 */
      SUBROUTINE CC_CIA4_2(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(#ai,#bj) = XINT(#ai,#bj)
C                            + P(L1,L2) sum(J) L1(#ai,J) * L2(#bj,J)
C
C              where P(L1,bj) symmetrizes with respect to ai,bj according to
C
C              P(L1,L2) M(L1,L2) = M(L1,L2) + M(L2,L1)
C
C     The symmetrization is done with 2 passes through the each of 
C     the Cholesky files, thus maximizing batch size at the expense of I/O.
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CIA4_2')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Calculate L1 * L2 contribution.
C     -------------------------------

      CALL CC_CIA4_3(XINT,WORK,LWORK,M12,NPASS)

C     Calculate L2 * L1 contribution.
C     -------------------------------

      CALL CC_CIASWP
      CALL CC_CIA4_3(XINT,WORK,LWORK,M21,NPASS)
      CALL CC_CIASWP

C     Get MEM.
C     --------

      MEM = MAX(M12,M21)

      RETURN
      END
C  /* Deck cc_cia4_3 */
      SUBROUTINE CC_CIA4_3(XINT,WORK,LWORK,MEM,NPASS)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate
C
C              XINT(#ai,#bj) = XINT(#ai,#bj) + sum(J) L1(#ai,J) * L2(#bj,J)
C
C     Specific information about which integrals are calculated is contained
C     in the configuration file ciarc.h. See documentation in cc_cia for more
C     detailed info.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ciarc.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CIA4_3')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO

C     Initial test of memory.
C     -----------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = MULD2H(ISYCHO,ISYCH2)
         IF ((NX1AMA(ISYMAI).GT.0) .AND. (NX1AMB(ISYMBJ).GT.0) .AND.
     &       (NTOVEC(ISYCHO).GT.0)) THEN
            NEEDS = NT1AM(ISYMAI) + NX1AMA(ISYMAI)
     &            + NT1AM(ISYMBJ) + NX1AMB(ISYMBJ)
            NEED  = MAX(NEED,NEEDS)
         ENDIF
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize memory stuff.
C     ------------------------

      MEM = 0

C     Loop over Cholesky symmetries (i.e. sym. of vector index J).
C     ------------------------------------------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH1)
         ISYMBJ = MULD2H(ISYCHO,ISYCH2)

C        Skip if nothing to do.
C        ----------------------

         IF (NTOVEC(ISYCHO) .LE. 0) GO TO 999
         IF (NX1AMA(ISYMAI) .LE. 0) GO TO 999
         IF (NX1AMB(ISYMBJ) .LE. 0) GO TO 999

C        Open Cholesky files.
C        --------------------

         CALL CHO_MOP(-1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)
         CALL CHO_MOP(-1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)

C        Allocation for full L1 and L2.
C        ------------------------------

         KCHO1 = 1
         KCHO2 = KCHO1 + NT1AM(ISYMAI)
         KEND1 = KCHO2 + NT1AM(ISYMBJ)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Cholesky'
            MINADD = NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
            NEED   = KEND1 - 1 + MINADD
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry      : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai         : ',ISYMAI
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of bj         : ',ISYMBJ
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',NEED,
     &      'Total memory available : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
         NVEC   = MIN(LWRK1/MINMEM,NTOVEC(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of ai                   : ',ISYMAI
            WRITE(LUPRI,'(5X,A,9X,I1)')
     &      'Symmetry of bj                   : ',ISYMBJ
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWRK1,
     &      'Minimum required (pref. multiple): ',MINMEM
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NTOVEC(ISYCHO) - 1)/NVEC + 1

C        Store batch info for later print.
C        ---------------------------------

         NCHBAT(ISYCHO) = NBATCH

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NTOVEC(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = KEND1
            KCHBJ = KCHAI + NX1AMA(ISYMAI)*NUMV
            KBLST = KCHBJ + NX1AMB(ISYMBJ)*NUMV - 1
            MEM   = MAX(MEM,KBLST)

C           Read MO Cholesky vectors and copy out subblocks.
C           ------------------------------------------------

            DO IVEC = 1,NUMV

               DTIME = SECOND()
               JVEC  = JVEC1 + IVEC - 1
               CALL CHO_MOREAD(WORK(KCHO1),NT1AM(ISYMAI),1,
     &                         JVEC,LUCHO1)
               CALL CHO_MOREAD(WORK(KCHO2),NT1AM(ISYMBJ),1,
     &                         JVEC,LUCHO2)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

               DTIME = SECOND()
               CALL CC_CIACPS(WORK(KCHO1),WORK(KCHAI),
     &                        IVEC,IOFA1,LVIRA,NX1AMA,
     &                        IX1AMA,ISYMAI)
               CALL CC_CIACPS(WORK(KCHO2),WORK(KCHBJ),
     &                        IVEC,IOFB1,LVIRB,NX1AMB,
     &                        IX1AMB,ISYMBJ)
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

            ENDDO

C           Calculate integral contribution.
C           --------------------------------

            DTIME = SECOND()

            NTOAI = MAX(NX1AMA(ISYMAI),1)
            NTOBJ = MAX(NX1AMB(ISYMBJ),1)

            KOFFX = IX2SQ(ISYMAI,ISYMBJ) + 1

            CALL DGEMM('N','T',NX1AMA(ISYMAI),NX1AMB(ISYMBJ),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                 ONE,XINT(KOFFX),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

         ENDDO

C        Close Cholesky files.
C        ---------------------

         CALL CHO_MOP(1,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)
         CALL CHO_MOP(1,ITYP1,ISYCHO,LUCHO1,1,ISYCH1)

  999    CONTINUE

      ENDDO

C     Update global timings.
C     ----------------------

      TIMCIA(2) = TIMCIA(2) + TIMR
      TIMCIA(3) = TIMCIA(3) + TIMS
      TIMCIA(4) = TIMCIA(4) + TIMC

C     Update Cholesky pass counter.
C     -----------------------------

      NPASS = NPASS + 2

C     Identify assembler routine.
C     ---------------------------

      CIAROU = SECNAM

      RETURN
      END
C  /* Deck cc_ciasmt */
      SUBROUTINE CC_CIASMT(XINT,ISYMX)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Symmetrize XINT(ai,bj) stored as full square.
C
#include <implicit.h>
      DIMENSION XINT(*)
#include <ccsdsym.h>
#include <ccorb.h>

      INTEGER AI, BJ, AIBJ, BJAI

      PARAMETER (ONE = 1.00D0)

      IF (ISYMX .EQ. 1) THEN

         DO ISYMBJ = 1,NSYM

            ISYMAI = ISYMBJ

            DO BJ = 1,NT1AM(ISYMBJ)
               DO AI = 1,BJ

                  AIBJ = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(BJ - 1)
     &                 + AI
                  BJAI = IT2SQ(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(AI - 1)
     &                 + BJ

                  XINT(AIBJ) = XINT(AIBJ) + XINT(BJAI)
                  XINT(BJAI) = XINT(AIBJ)

               ENDDO
            ENDDO

         ENDDO

      ELSE

         DO ISYMBJ = 1,NSYM

            ISYMAI = MULD2H(ISYMBJ,ISYMX)

            IF (ISYMAI .LT. ISYMBJ) THEN

               DO BJ = 1,NT1AM(ISYMBJ)

                  KOFF1 = IT2SQ(ISYMBJ,ISYMAI) + BJ
                  KOFF2 = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(BJ - 1)
     &                  + 1

                  CALL DAXPY(NT1AM(ISYMAI),ONE,
     &                       XINT(KOFF1),NT1AM(ISYMBJ),XINT(KOFF2),1)

                  CALL DCOPY(NT1AM(ISYMAI),
     &                       XINT(KOFF2),1,XINT(KOFF1),NT1AM(ISYMBJ))

               ENDDO

            ENDIF

         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_ciacps */
      SUBROUTINE CC_CIACPS(CHOL,CHAI,IVEC,IOFA1,LVIRA,NX1AMA,IX1AMA,
     &                     ISYMAI)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Copy out subblocks of the Cholesky vectors for batches
C              of virtuals.
C
#include <implicit.h>
      DIMENSION CHOL(*), CHAI(*)
      INTEGER   IOFA1(8), LVIRA(8), NX1AMA(8), IX1AMA(8,8)
#include <ccsdsym.h>
#include <ccorb.h>

      IF (NX1AMA(ISYMAI) .EQ. NT1AM(ISYMAI)) THEN

         KOFF1 = NT1AM(ISYMAI)*(IVEC - 1) + 1

         CALL DCOPY(NT1AM(ISYMAI),CHOL,1,CHAI(KOFF1),1)

      ELSE

         DO ISYMI = 1,NSYM
            IF (NRHF(ISYMI) .GT. 0) THEN

               ISYMA = MULD2H(ISYMI,ISYMAI)

               IF (LVIRA(ISYMA) .EQ. NVIR(ISYMA)) THEN

                  NAI = NVIR(ISYMA)*NRHF(ISYMI)

                  KOFF1 = IT1AM(ISYMA,ISYMI) + 1
                  KOFF2 = NX1AMA(ISYMAI)*(IVEC - 1)
     &                  + IX1AMA(ISYMA,ISYMI) + 1

                  CALL DCOPY(NAI,CHOL(KOFF1),1,CHAI(KOFF2),1)

               ELSE IF (LVIRA(ISYMA) .GT. 0) THEN

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = IT1AM(ISYMA,ISYMI)
     &                     + NVIR(ISYMA)*(I - 1) + IOFA1(ISYMA)
                     KOFF2 = NX1AMA(ISYMAI)*(IVEC - 1)
     &                     + IX1AMA(ISYMA,ISYMI)
     &                     + LVIRA(ISYMA)*(I - 1) + 1

                     CALL DCOPY(LVIRA(ISYMA),CHOL(KOFF1),1,
     &                                       CHAI(KOFF2),1)

                  ENDDO

               ENDIF

            ENDIF
         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_ciaswp */
      SUBROUTINE CC_CIASWP
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Swap L1 and L2 information in ciarc.h
C
#include <implicit.h>
#include <ccorb.h>
#include <ciarc.h>

C     Swap symmetry info.
C     -------------------

      ISTMP  = ISYCH1
      ISYCH1 = ISYCH2
      ISYCH2 = ISTMP

C     Swap vector type info.
C     ----------------------

      ITMP  = ITYP1
      ITYP1 = ITYP2
      ITYP2 = ITMP

      RETURN
      END
C  /* Deck cc_ciamnm */
      SUBROUTINE CC_CIAMNM(MEM)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate min. mem. requirement for running
C              CC_CIA with one a and one b in each call.
C              The result includes the memory needed for the
C              integral array!
C              
C     NOTICE: THIS ROUTINE MUST BE CHANGED WHENEVER CHANGES ARE
C             MADE IN ALLOCATIONS IN CC_CIA ROUTINES.
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ciarc.h>

C     Initialize.
C     -----------

      NEED = 0

      IF (LIDEN) THEN

C        Calculate need for L1 = L2.
C        ===========================

         DO ISYMB = 1,NSYM
            DO ISYMA = 1,ISYMB
               ISYMIJ = MULD2H(ISYMA,ISYMB)
               NEEDC  = 0
               DO ISYCHO = 1,NSYM
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  ISYMBJ = ISYMAI
                  ISYMI  = MULD2H(ISYMAI,ISYMA)
                  ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                  NTST   = NT1AM(ISYMAI) + NRHF(ISYMI) + NRHF(ISYMJ)
                  NEEDC  = MAX(NEEDC,NTST)
               ENDDO
               NTST = NMATIJ(ISYMIJ) + NEEDC
               NEED = MAX(NEED,NTST)
            ENDDO
         ENDDO

      ELSE

C        Calculate need for L1 != L2.
C        ============================

         ISYINT = MULD2H(ISYCH1,ISYCH2)

         IF (SYMTRZ) THEN

C           Symmetrize option.
C           ------------------

            IF (CIAMIO) THEN

C              Min. I/O symmetrization requires more memory.
C              ---------------------------------------------

               DO ISYMB = 1,NSYM
                  DO ISYMA = 1,NSYM
                     ISYMAB = MULD2H(ISYMA,ISYMB)
                     ISYMIJ = MULD2H(ISYMAB,ISYINT)
                     NEEDC  = 0
                     DO ISYCHO = 1,NSYM
                        ISYMAI = MULD2H(ISYCHO,ISYCH1)
                        ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                        ISYMI  = MULD2H(ISYMAI,ISYMA)
                        ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                        NTST   = NT1AM(ISYMAI) + NRHF(ISYMI)
     &                         + NT1AM(ISYMBJ) + NRHF(ISYMJ)
                        ISYMAI = MULD2H(ISYCHO,ISYCH2)
                        ISYMBJ = MULD2H(ISYCHO,ISYCH1)
                        ISYMI  = MULD2H(ISYMAI,ISYMA)
                        ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                        NTST   = NTST + NRHF(ISYMI) + NRHF(ISYMJ)
                        NEEDC  = MAX(NEEDC,NTST)
                     ENDDO
                     NTST = NMATIJ(ISYMIJ) + NEEDC
                     NEED = MAX(NEED,NTST)
                  ENDDO
               ENDDO

            ELSE

C              Symmetrization with 2 passes through Cholesky files.
C              (In fact the same as no symmetrization in terms of
C              memory requirement).
C              ----------------------------------------------------

               DO ISYMB = 1,NSYM
                  DO ISYMA = 1,NSYM
                     ISYMAB = MULD2H(ISYMA,ISYMB)
                     ISYMIJ = MULD2H(ISYMAB,ISYINT)
                     NEEDC  = 0
                     DO ISYCHO = 1,NSYM
                        ISYMAI = MULD2H(ISYCHO,ISYCH1)
                        ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                        ISYMI  = MULD2H(ISYMAI,ISYMA)
                        ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                        NTST   = NT1AM(ISYMAI) + NRHF(ISYMI)
     &                         + NT1AM(ISYMBJ) + NRHF(ISYMJ)
                        NEEDC  = MAX(NEEDC,NTST)
                     ENDDO
                     NTST = NMATIJ(ISYMIJ) + NEEDC
                     NEED = MAX(NEED,NTST)
                  ENDDO
               ENDDO
            ENDIF

         ELSE

C           No symmetrization.
C           ------------------

            DO ISYMB = 1,NSYM
               DO ISYMA = 1,NSYM
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  ISYMIJ = MULD2H(ISYMAB,ISYINT)
                  NEEDC  = 0
                  DO ISYCHO = 1,NSYM
                     ISYMAI = MULD2H(ISYCHO,ISYCH1)
                     ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                     ISYMI  = MULD2H(ISYMAI,ISYMA)
                     ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                     NTST   = NT1AM(ISYMAI) + NRHF(ISYMI)
     &                      + NT1AM(ISYMBJ) + NRHF(ISYMJ)
                     NEEDC  = MAX(NEEDC,NTST)
                  ENDDO
                  NTST = NMATIJ(ISYMIJ) + NEEDC
                  NEED = MAX(NEED,NTST)
               ENDDO
            ENDDO

         ENDIF

      ENDIF

C     Return MEM.
C     -----------

      MEM = NEED

      RETURN
      END
