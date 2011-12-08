C  /* Deck iorder3 */
      SUBROUTINE IORDER3(IVAL,N,IJOB)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Sort IVAL (dimension: N) in ascending (IJOB>0)
C              or descending (IJOB<0) order.
C
C     Integer version of ORDER3 by T. Saue.
C     See the file gp/gpsaue.F for further details.
C
#include <implicit.h>
      INTEGER IVAL(*)
      PARAMETER (MONE = -1)

      IF ((N.LE.1) .OR. (IJOB.EQ.0)) RETURN

      IF (IJOB .LT. 0) CALL ISCAL(N,MONE,IVAL,1)

      DO I = 1,N-1

         IVMIN = IVAL(I)
         IMIN  = I

         DO J = I+1,N
            IF (IVAL(J) .LT. IVMIN) THEN
               IVMIN = IVAL(J)
               IMIN  = J
            ENDIF
         ENDDO

         IF (IMIN .NE. I) THEN
            IVAL(IMIN) = IVAL(I)
            IVAL(I)    = IVMIN
         ENDIF

      ENDDO

      IF (IJOB .LT. 0) CALL ISCAL(N,MONE,IVAL,1)

      RETURN
      END
C  /* Deck op_mp2tam */
      SUBROUTINE OP_MP2TAM(IOPT,LSYM,LUNIT,MSYM)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C
C        IOPT < 0: Open direct access file(s) for MP2 amplitudes.
C
C        IOPT = 0: Close and delete the files.
C
C        IOPT > 0: Close and keep the files.
C
#include <implicit.h>
      INTEGER LSYM(MSYM), LUNIT(MSYM)
#include <ccorb.h>
#include <ccsdsym.h>

      LOGICAL OLDDX

      CHARACTER*12 FILNAM(8)
      DATA FILNAM /'CHOMP2_TAM_1','CHOMP2_TAM_2','CHOMP2_TAM_3',
     &             'CHOMP2_TAM_4','CHOMP2_TAM_5','CHOMP2_TAM_6',
     &             'CHOMP2_TAM_7','CHOMP2_TAM_8'/

      IF (IOPT .LT. 0) THEN

         DO ISM = 1,MSYM

            ISYM  = LSYM(ISM)

            LRLEN = 2*NT1AM(ISYM)
            LRLEN = MAX(LRLEN,1)

            LUTAM = -1
            CALL GPOPEN(LUTAM,FILNAM(ISYM),'UNKNOWN','DIRECT',
     &                  'UNFORMATTED',LRLEN,OLDDX)
            LUNIT(ISM) = LUTAM

         ENDDO

      ELSE IF (IOPT .EQ. 0) THEN

         DO ISM = 1,MSYM
            ISYM  = LSYM(ISM)
            LUTAM = LUNIT(ISM)
            CALL GPCLOSE(LUTAM,'DELETE')
            LUNIT(ISM) = -1
         ENDDO

      ELSE

         DO ISM = 1,MSYM
            ISYM  = LSYM(ISM)
            LUTAM = LUNIT(ISM)
            CALL GPCLOSE(LUTAM,'KEEP')
            LUNIT(ISM) = -1
         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck wr_mp2tam */
      SUBROUTINE WR_MP2TAM(TAM,LTAM,BJ,LUNIT)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Write MP2 amplitudes.
C
C     Notes:
C       The first record on the direct access file is reserved
C       for the diagonal; to write to the diagonal record, specify
C       BJ=0.
C
C       The file is assumed open on entry, with logical unit number
C       LUNIT.
C
#include <implicit.h>
      DIMENSION TAM(LTAM)
      INTEGER BJ
#include <ccsdsym.h>

      IREC = BJ + 1
      WRITE(LUNIT,REC=IREC) (TAM(I),I=1,LTAM)

      RETURN
      END
C  /* Deck rd_mp2tam */
      SUBROUTINE RD_MP2TAM(TAM,LTAM,BJ,LUNIT)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Read MP2 amplitudes.
C
C     Notes:
C       The first record on the direct access file is reserved
C       for the diagonal; to read the diagonal record, specify
C       BJ=0.
C
C       The file is assumed open on entry, with logical unit number
C       LUNIT.
C
#include <implicit.h>
      DIMENSION TAM(LTAM)
      INTEGER BJ
#include <ccsdsym.h>

      IREC = BJ + 1
      READ(LUNIT,REC=IREC) (TAM(I),I=1,LTAM)

      RETURN
      END
C  /* Deck cho_mop */
      SUBROUTINE CHO_MOP(IOPT,ITYP,LSYM,LUNIT,MSYM,ISYMCH)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C
C        IOPT < 0: Open direct access file(s).
C        IOPT = 0: Close and delete the files.
C        IOPT > 0: Close and keep the files.
C
C        ITYP =  1: L(ia,J)   - from AO Cholesky vectors.
C        ITYP =  2: M(ia,J)   - from (ia|jb) decomposition.
C        ITYP =  3: L(ai,J)   - from AO Cholesky vectors.
C        ITYP =  4: L(ij,J)   - from AO Cholesky vectors.
C        ITYP =  5: M(ai,J)   - from (ai|bj) decomposition.
C        ITYP =  6: M(ai,J)   - from CC2 amplitude decomposition.
C        ITYP =  7: LD1(ai,J) - Dummy file.
C        ITYP =  8: LT(ai,J)  - perturbed MO Cholesky vectors
C                               for left-hand Jacobian transformations.
C        ITYP =  9: LB(ai,J)  - perturbed MO Cholesky vectors
C                               for right-hand Jacobian transformations.
C        ITYP = 10: LB(ij,J)  - perturbed MO Cholesky vectors
C                               for right-hand Jacobian transformations.
C                               (Used in F-matrix section).
C        ITYP = 11: LD2(ai,J) - Dummy file.
C
C        LSYM  : List of symmetries of the vector index (J).
C        LUNIT : List of unit numbers.
C        MSYM  : Dimension of LSYM and LUNIT lists.
C        ISYMCH: Overall symmetry of the Cholesky vectors to be written
C                (= 1 for ground state vectors). I.e. such that
C                   ISYMAI = MULD2H(ISYMJ,ISYMCH) [for L(ai,J) vectors]
C                          or
C                   ISYMIJ = MULD2H(ISYMJ,ISYMCH) [for L(ij,J) vectors]
C                where ISYMJ refer to the symmetry of the vector index J.
C
C     The symmetry index on the file names refer to the symmetry of
C     the J-index (vector index); this is, of course, irrelevant for
C     total symmetric Cholesky vectors.
C
#include <implicit.h>
      INTEGER LSYM(MSYM),LUNIT(MSYM)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <cyit.h>

      LOGICAL OLDDX

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CHO_MOP')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (NTYP = 11)

      CHARACTER*7 FILNAM(8,NTYP)
      DATA FILNAM /'CHOIA_1','CHOIA_2','CHOIA_3','CHOIA_4',
     &             'CHOIA_5','CHOIA_6','CHOIA_7','CHOIA_8',
     &             'CIAJB_1','CIAJB_2','CIAJB_3','CIAJB_4',
     &             'CIAJB_5','CIAJB_6','CIAJB_7','CIAJB_8',
     &             'CHOAI_1','CHOAI_2','CHOAI_3','CHOAI_4',
     &             'CHOAI_5','CHOAI_6','CHOAI_7','CHOAI_8',
     &             'CHOIJ_1','CHOIJ_2','CHOIJ_3','CHOIJ_4',
     &             'CHOIJ_5','CHOIJ_6','CHOIJ_7','CHOIJ_8',
     &             'CAIBJ_1','CAIBJ_2','CAIBJ_3','CAIBJ_4',
     &             'CAIBJ_5','CAIBJ_6','CAIBJ_7','CAIBJ_8',
     &             'CC2AM_1','CC2AM_2','CC2AM_3','CC2AM_4',
     &             'CC2AM_5','CC2AM_6','CC2AM_7','CC2AM_8',
     &             'CHLD1_1','CHLD1_2','CHLD1_3','CHLD1_4',
     &             'CHLD1_5','CHLD1_6','CHLD1_7','CHLD1_8',
     &             'CHLIA_1','CHLIA_2','CHLIA_3','CHLIA_4',
     &             'CHLIA_5','CHLIA_6','CHLIA_7','CHLIA_8',
     &             'CHRAI_1','CHRAI_2','CHRAI_3','CHRAI_4',
     &             'CHRAI_5','CHRAI_6','CHRAI_7','CHRAI_8',
     &             'CHRIJ_1','CHRIJ_2','CHRIJ_3','CHRIJ_4',
     &             'CHRIJ_5','CHRIJ_6','CHRIJ_7','CHRIJ_8',
     &             'CHLD2_1','CHLD2_2','CHLD2_3','CHLD2_4',
     &             'CHLD2_5','CHLD2_6','CHLD2_7','CHLD2_8'/

      TIMT = SECOND()

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   'Entering ',SECNAM,' - input variables:'
         WRITE(LUPRI,'(5X,A,I10)') 'IOPT  : ',IOPT
         WRITE(LUPRI,'(5X,A,I10)') 'ITYP  : ',ITYP
         WRITE(LUPRI,'(5X,A,I10)') 'MSYM  : ',MSYM
         WRITE(LUPRI,'(5X,A,I10)') 'ISYMCH: ',ISYMCH
         WRITE(LUPRI,'(5X,A)') 'LSYM: '
         WRITE(LUPRI,'(8I10)') (LSYM(I),I=1,MSYM)
         WRITE(LUPRI,'(5X,A)') 'LUNIT: '
         WRITE(LUPRI,'(8I10)') (LUNIT(I),I=1,MSYM)
         WRITE(LUPRI,'(A)') ' '
      ENDIF

C     Check ITYP.
C     -----------

      IF ((ITYP.LT.1) .OR. (ITYP.GT.NTYP)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,': specifier ITYP not recognized.'
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ITYP = ',ITYP
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      IF (IOPT .LT. 0) THEN

C        Open direct-access files.
C        -------------------------

         DO ISM = 1,MSYM

            ISYMJ  = LSYM(ISM)
            ISYMAI = MULD2H(ISYMJ,ISYMCH)

            IF ((ITYP.EQ.4) .OR. (ITYP.EQ.10)) THEN
               LRLEN = 2*NMATIJ(ISYMAI)
            ELSE
               LRLEN = 2*NT1AM(ISYMAI)
            ENDIF
            LRLEN = MAX(LRLEN,1)

            LUVEC = -1
            CALL GPOPEN(LUVEC,FILNAM(ISYMJ,ITYP),'UNKNOWN','DIRECT',
     &                  'UNFORMATTED',LRLEN,OLDDX)
            LUNIT(ISM) = LUVEC

         ENDDO

      ELSE IF (IOPT .EQ. 0) THEN

C        Close and delete files.
C        -----------------------

         DO ISM = 1,MSYM
            LUVEC = LUNIT(ISM)
            CALL GPCLOSE(LUVEC,'DELETE')
            LUNIT(ISM) = -1
         ENDDO

      ELSE

C        Close and keep files.
C        ---------------------

         DO ISM = 1,MSYM
            LUVEC = LUNIT(ISM)
            CALL GPCLOSE(LUVEC,'KEEP')
            LUNIT(ISM) = -1
         ENDDO

      ENDIF

      TIMT = SECOND() - TIMT
      TMCYIO(1) = TMCYIO(1) + TIMT

      RETURN
      END
C  /* Deck cho_mowrite */
      SUBROUTINE CHO_MOWRITE(VEC,NAI,NVEC,JVEC1,LUNIT)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Write MO Cholesky vectors.
C
C     Notes:
C       The file is assumed open on entry, with logical unit number
C       LUNIT.
C
#include <implicit.h>
      DIMENSION VEC(NAI,NVEC)
#include <cyit.h>

      TIMT = SECOND()
      DO JVEC = 1,NVEC
         IREC = JVEC1 + JVEC - 1
         CALL CHO_WRDA(LUNIT,IREC,VEC(1,JVEC),NAI)
      ENDDO
      TIMT = SECOND() - TIMT
      TMCYIO(1) = TMCYIO(1) + TIMT

      RETURN
      END
C  /* Deck cho_moread */
      SUBROUTINE CHO_MOREAD(VEC,NAI,NVEC,JVEC1,LUNIT)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Read MO Cholesky vectors.
C
C     Notes:
C       The file is assumed open on entry, with logical unit number
C       LUNIT.
C
#include <implicit.h>
      DIMENSION VEC(NAI,NVEC)
#include <cyit.h>

      TIMT = SECOND()
      DO JVEC = 1,NVEC
         IREC = JVEC1 + JVEC - 1
         CALL CHO_RDDA(LUNIT,IREC,VEC(1,JVEC),NAI)
      ENDDO
      TIMT = SECOND() - TIMT
      TMCYIO(1) = TMCYIO(1) + TIMT

      RETURN
      END
C  /* Deck cho_rdda */
      SUBROUTINE CHO_RDDA(LUNIT,IREC,DATA,LEN)
C
C     Purpose:
C        Read record IREC from direct access file.
C
#include <implicit.h>
      DIMENSION DATA(LEN)

      READ(LUNIT,REC=IREC) (DATA(I),I=1,LEN)

      RETURN
      END
C  /* Deck cho_wrda */
      SUBROUTINE CHO_WRDA(LUNIT,IREC,DATA,LEN)
C
C     Purpose:
C        Write record IREC from direct access file.
C
#include <implicit.h>
      DIMENSION DATA(LEN)

      WRITE(LUNIT,REC=IREC) (DATA(I),I=1,LEN)

      RETURN
      END
C  /* Deck cho_readn */
      SUBROUTINE CHO_READN(CHOVEC,IFIRST,NREAD,IND1,IND2,ISYCHO,ISTORE,
     &                     WORK,LWORK)
C
C     - renamed, TBP July 2005.
C
C     Purpose: Read NREAD Cholesky vectors (AO) starting with IFIRST
C              of symmetry ISYCHO. Calls CC_GETCHO for the reading
C              process.
C
C     ISTORE requests the storage format on output:
C
C        ISTORE = 1 --- not implemented; reserved for compressed vectors.
C
C        ISTORE = 2 --- Full square, incl. zeros.
C                       IND1 and IND2: dummies.
C
C        ISTORE = 3 --- Packed, incl. zeros.
C                       IND1 and IND2: dummies.
C
C     IF (REDUCE): IND1 must contain reduce index array on input.
C
#include <implicit.h>
#include <maxorb.h>
      DIMENSION CHOVEC(*), WORK(LWORK)
      INTEGER IND1(*), IND2(*)
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_READN')

C     Read the vectors.
C     -----------------

      CALL CC_GETCHO(CHOVEC,NREAD,ISYCHO,IFIRST,IND1,WORK,LWORK)
c     xnorm = ddot(nnbst(isycho),chovec,1,chovec,1)
c     write(lupri,'(a,i4,a,i2,a,d12.4)') 'Vector',jchol,
c    &     ' of symmetry',isycho,'. Norm after GETCHO: ',xnorm

C     Change storage format.
C     ======================

      IF (ISTORE .EQ. 2) THEN

C        Square.
C        -------

         IF (LWORK .LT. NREAD*NNBST(ISYCHO)) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for squaring vectors in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',NREAD*NNBST(ISYCHO),
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL DCOPY(NREAD*NNBST(ISYCHO),CHOVEC,1,WORK,1)
         DO IREAD = 1,NREAD
            KOFF1 = NNBST(ISYCHO)*(IREAD - 1) + 1
            KOFF2 = N2BST(ISYCHO)*(IREAD - 1) + 1
            CALL CCSD_SYMSQ(WORK(KOFF1),ISYCHO,CHOVEC(KOFF2))
         ENDDO

      ELSE IF (ISTORE .EQ. 3) THEN

C        Keep the format from CC_GETCHO.
C        -------------------------------

         RETURN

      ELSE

C        Unknown ISTORE.
C        ---------------

         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Illegal ISTORE request in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ISTORE = ',ISTORE
         CALL QUIT('Error in '//SECNAM)

      ENDIF

      RETURN
      END
C  /* Deck cc_getind1 */
      SUBROUTINE CC_GETIND1(WORK,LWORK,LIND1)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Read the reduce index vector from the Cholesky decomposition
C        into work space. On exit, LIND1 is the length of the index array.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <iratdef.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_GETIND1')

      IF (REDUCE) THEN

C        Determine length of index array.
C        --------------------------------

         LREDU = 0
         DO ISYM = 1,NSYM
            LREDU = MAX(LREDU,NNBST(ISYM))
         ENDDO
         LIND1 = NSYM*((LREDU - 1)/IRAT + 1)

         IF (LWORK .LT. LIND1) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Needed for index array: ',LIND1,
     &      'Available             : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read index array.
C        -----------------

         LREAD = NSYM*LREDU
         CALL CC_RDIND1(WORK,LREAD)

      ELSE

C        Cholesky vectors not reduced.
C        -----------------------------

         LIND1 = 0

      ENDIF

      RETURN
      END
C  /* Deck cc_rdind1 */
      SUBROUTINE CC_RDIND1(IND1,LEN)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose: Read reduce index array.
C
#include <implicit.h>
      INTEGER IND1(LEN)

      LOGICAL LOPEN

      INQUIRE(FILE='CHOLESKY.RST',OPENED=LOPEN)

      IF (.NOT. LOPEN) THEN

         LURST = -1
         CALL GPOPEN(LURST,'CHOLESKY.RST','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LURST)
         READ(LURST)
         READ(LURST) (IND1(I),I=1,LEN)
         CALL GPCLOSE(LURST,'KEEP')

      ELSE

         INQUIRE(FILE='CHOLESKY.RST',NUMBER=LURST)
         REWIND(LURST)
         READ(LURST)
         READ(LURST) (IND1(I),I=1,LEN)

      ENDIF

      RETURN
      END
C  /* Deck memrd */
      INTEGER FUNCTION MEMRD(NUMVEC,ISYCHO,ISTORE)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Function MEMRD returns the memory needed in
C        CHO_READN (incl. call to CC_GETCHO).
C
#include <implicit.h>
#include <iratdef.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>

      MEMRD = 1

      IF ((COMP) .AND. (.NOT. REDUCE)) THEN

C        AO Cholesky vectors are compressed (excl. zeros).
C        Memory needed to "pad" zeros and read index arrays.
C        ===================================================

         IF ((ISTORE.EQ.2) .OR. (ISTORE.EQ.3)) THEN
            NNVEC = NNBST(ISYCHO) + (NNBST(ISYCHO) - 1)/IRAT + 2
            MEMRD = NUMVEC*NNVEC + NNBST(ISYCHO)
         ENDIF

      ELSE IF (REDUCE) THEN

C        AO Cholesky vectors are stored with zeros.
C        Memory needed for expanding vectors and/or squaring.
C        ====================================================

         IF (ISTORE.EQ.3) THEN
            MEMRD = NUMVEC*NREDUC(ISYCHO)
         ELSE IF (ISTORE .EQ. 2) THEN
            MEMRD = NUMVEC*NNBST(ISYCHO)
         ENDIF

      ELSE

C        AO Cholesky vectors are stored with zeros.
C        Memory only needed for squaring in CHO_READN.
C        =============================================

         IF (ISTORE .EQ. 2) THEN
            MEMRD = NUMVEC*NNBST(ISYCHO)
         ENDIF

      ENDIF

      RETURN
      END
C  /* Deck onel_op */
      SUBROUTINE ONEL_OP(IOPT,ITYP,LUVEC)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C
C        IOPT < 0: Open direct access file.
C        IOPT = 0: Close and delete the file.
C        IOPT > 0: Close and keep the file.
C
C        ITYP = 1: h(ia)
C        ITYP = 2: h(alpha i)
C        ITYP = 3: F(ia)
C
C        LUVEC: Unit number.
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL OLDDX

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'ONEL_OP')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (NTYP = 3)

      CHARACTER*5 FILNAM(NTYP)
      DATA FILNAM /'ONEIA','ONALI','FCKIA'/

C     Check ITYP.
C     -----------

      IF ((ITYP.LT.1) .OR. (ITYP.GT.NTYP)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,': specifier ITYP not recognized.'
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ITYP = ',ITYP
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      IF (IOPT .LT. 0) THEN

C        Open direct-access files.
C        Record length is machine dependent.
C        -----------------------------------

         IF (ITYP .EQ. 2) THEN
            LRLEN = 2*NT1AO(1)
         ELSE
            LRLEN = 2*NT1AMX
         ENDIF
         LRLEN = MAX(LRLEN,1)

         LUVEC = -1
         CALL GPOPEN(LUVEC,FILNAM(ITYP),'UNKNOWN','DIRECT',
     &               'UNFORMATTED',LRLEN,OLDDX)

      ELSE IF (IOPT .EQ. 0) THEN

C        Close and delete file.
C        ----------------------

         CALL GPCLOSE(LUVEC,'DELETE')
         LUVEC = -1

      ELSE

C        Close and keep file.
C        --------------------

         CALL GPCLOSE(LUVEC,'KEEP')
         LUVEC = -1

      ENDIF

      RETURN
      END
C  /* Deck cc_cyitrnsl */
      SUBROUTINE CC_CYITRNSL(JSIDE,ISIDE)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Cheat the CC_CYI* I/O routines.
C
#include <implicit.h>

      JSIDE = ISIDE

      IF (ISIDE .EQ. -3) THEN
         JSIDE = -1
      ELSE IF (ISIDE .EQ. -2) THEN
         JSIDE = -1
      ELSE IF (ISIDE .EQ. 3) THEN
         JSIDE = 1
      ELSE IF (ISIDE .EQ. 33) THEN
         JSIDE = 3
      ENDIF

      RETURN
      END
C  /* Deck cc_cyiop */
      SUBROUTINE CC_CYIOP(IOPT,ISYCHO,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: 
C
C        IOPT < 0: Open file for Y(ai,J) intermediates.
C        IOPT = 0: Close and delete file for Y(ai,J) intermediates.
C        IOPT > 0: Close and keep file for Y(ai,J) intermediates.
C
C     ISYCHO is the symmetry of the vector index J, and ISIDE specifies
C     which transformation we are doing.
C
#include <implicit.h>
#include <chocc2.h>
#include <chocc2lhs.h>
#include <chocc2rhs.h>
#include <chocc2ftr.h>
#include <priunit.h>
#include <cyit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CYIOP')

      TIMT = SECOND()

      CALL CC_CYITRNSL(JSIDE,ISIDE)
      IF (IOPT .LT. 0) THEN

         IF (JSIDE .EQ. -1) THEN
            CALL WOPEN2(LCC2YL,FCC2YL(ISYCHO),64,0)
         ELSE IF (JSIDE .EQ. 0) THEN
            CALL WOPEN2(LCC2YM,FCC2YM(ISYCHO),64,0)
         ELSE IF (JSIDE .EQ. 1) THEN
            CALL WOPEN2(LCC2YR,FCC2YR(ISYCHO),64,0)
         ELSE IF (JSIDE .EQ. 2) THEN
            CALL WOPEN2(LCC2Y1,FCC2Y1(ISYCHO),64,0)
         ELSE IF (JSIDE .EQ. 3) THEN
            CALL WOPEN2(LCC2Y2,FCC2Y2(ISYCHO),64,0)
         ELSE IF (JSIDE .EQ. 4) THEN
            CALL WOPEN2(LCC2Y3,FCC2Y3(ISYCHO),64,0)
         ELSE
            WRITE(LUPRI,'(//,10X,A,A,A,I10)')
     &      '*** WARNING: Illegal ISIDE in ',SECNAM,': ',ISIDE
            WRITE(LUPRI,'(10X,A,I2/,10X,A,//)')
     &      '             No file has been opened for ISYCHO =',ISYCHO,
     &      '             Program continues nevetheless....'
         ENDIF

      ELSE IF (IOPT .EQ. 0) THEN

         IF (JSIDE .EQ. -1) THEN
            CALL WCLOSE2(LCC2YL,FCC2YL(ISYCHO),'DELETE')
         ELSE IF (JSIDE .EQ. 0) THEN
            CALL WCLOSE2(LCC2YM,FCC2YM(ISYCHO),'DELETE')
         ELSE IF (JSIDE .EQ. 1) THEN
            CALL WCLOSE2(LCC2YR,FCC2YR(ISYCHO),'DELETE')
         ELSE IF (JSIDE .EQ. 2) THEN
            CALL WCLOSE2(LCC2Y1,FCC2Y1(ISYCHO),'DELETE')
         ELSE IF (JSIDE .EQ. 3) THEN
            CALL WCLOSE2(LCC2Y2,FCC2Y2(ISYCHO),'DELETE')
         ELSE IF (JSIDE .EQ. 4) THEN
            CALL WCLOSE2(LCC2Y3,FCC2Y3(ISYCHO),'DELETE')
         ELSE
            WRITE(LUPRI,'(//,10X,A,A,A,I10)')
     &      '*** WARNING: Illegal ISIDE in ',SECNAM,': ',ISIDE
            WRITE(LUPRI,'(10X,A,I2,/,10X,A,//)')
     &      '             No file has been deleted for ISYCHO =',ISYCHO,
     &      '             Program continues nevetheless....'
         ENDIF

      ELSE

         IF (JSIDE .EQ. -1) THEN
            CALL WCLOSE2(LCC2YL,FCC2YL(ISYCHO),'KEEP')
         ELSE IF (JSIDE .EQ. 0) THEN
            CALL WCLOSE2(LCC2YM,FCC2YM(ISYCHO),'KEEP')
         ELSE IF (JSIDE .EQ. 1) THEN
            CALL WCLOSE2(LCC2YR,FCC2YR(ISYCHO),'KEEP')
         ELSE IF (JSIDE .EQ. 2) THEN
            CALL WCLOSE2(LCC2Y1,FCC2Y1(ISYCHO),'KEEP')
         ELSE IF (JSIDE .EQ. 3) THEN
            CALL WCLOSE2(LCC2Y2,FCC2Y2(ISYCHO),'KEEP')
         ELSE IF (JSIDE .EQ. 4) THEN
            CALL WCLOSE2(LCC2Y3,FCC2Y3(ISYCHO),'KEEP')
         ELSE
            WRITE(LUPRI,'(//,10X,A,A,A,I10)')
     &      '*** WARNING: Illegal ISIDE in ',SECNAM,': ',ISIDE
            WRITE(LUPRI,'(10X,A,I2,/,10X,A,//)')
     &      '             No file has been closed for ISYCHO =',ISYCHO,
     &      '             Program continues nevetheless....'
         ENDIF

      ENDIF

      TIMT = SECOND() - TIMT
      TMCYIO(2) = TMCYIO(2) + TIMT

      RETURN
      END
C  /* Deck cc_cyirdf */
      SUBROUTINE CC_CYIRDF(YIM,NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Read Y(ai,#J) intermediates from disk.
C     
C     ISYCHO is the symmetry of J-index.
C     ISYMY  is symmetry of the Y intermediates (overall).
C     ISIDE  determines which files are used (according to
C            type of Y intermediates, i.e. LH trf., gr.state,
C            or RH trf.)
C
#include <implicit.h>
      DIMENSION YIM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <chocc2lhs.h>
#include <chocc2rhs.h>
#include <chocc2ftr.h>
#include <priunit.h>
#include <cyit.h>

      INTEGER IADR, ICOL, LROW

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIRDF')

      TIMT = SECOND()

      CALL CC_CYITRNSL(JSIDE,ISIDE)

      ISYMAI = MULD2H(ISYCHO,ISYMY)
      ICOL   = JVEC1 - 1
      LROW   = NT1AM(ISYMAI)
      IADR   = LROW*ICOL + 1
      LEN    = NT1AM(ISYMAI)*NUMV

      IF (LEN .GT. 0) THEN
         IF (JSIDE .EQ. -1) THEN
            CALL GETWA2(LCC2YL,FCC2YL(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 0) THEN
            CALL GETWA2(LCC2YM,FCC2YM(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 1) THEN
            CALL GETWA2(LCC2YR,FCC2YR(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 2) THEN
            CALL GETWA2(LCC2Y1,FCC2Y1(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 3) THEN
            CALL GETWA2(LCC2Y2,FCC2Y2(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 4) THEN
            CALL GETWA2(LCC2Y3,FCC2Y3(ISYCHO),YIM,IADR,LEN)
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &      'Error in ',SECNAM,': Unknown option, ISIDE = ',ISIDE
            CALL QUIT('Error in '//SECNAM)
         ENDIF
      ENDIF
      TIMT = SECOND() - TIMT
      TMCYIO(2) = TMCYIO(2) + TIMT

      RETURN
      END
C  /* Deck cc_cyiwrf */
      SUBROUTINE CC_CYIWRF(YIM,NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Write Y(ai,#J) intermediates to disk.
C     
C     ISYCHO is the symmetry of J-index.
C     ISYMY  is symmetry of the Y intermediates (overall).
C     ISIDE  determines which files are used (according to
C            type of Y intermediates, i.e. LH trf., gr.state,
C            or RH trf.)
C
#include <implicit.h>
      DIMENSION YIM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <chocc2lhs.h>
#include <chocc2rhs.h>
#include <chocc2ftr.h>
#include <priunit.h>
#include <cyit.h>

      INTEGER IADR, ICOL, LROW

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIWRF')

      TIMT = SECOND()

      CALL CC_CYITRNSL(JSIDE,ISIDE)

      ISYMAI = MULD2H(ISYCHO,ISYMY)
      ICOL   = JVEC1 - 1
      LROW   = NT1AM(ISYMAI)
      IADR   = LROW*ICOL + 1
      LEN    = NT1AM(ISYMAI)*NUMV

      IF (LEN .GT. 0) THEN
         IF (JSIDE .EQ. -1) THEN
            CALL PUTWA2(LCC2YL,FCC2YL(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 0) THEN
            CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 1) THEN
            CALL PUTWA2(LCC2YR,FCC2YR(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 2) THEN
            CALL PUTWA2(LCC2Y1,FCC2Y1(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 3) THEN
            CALL PUTWA2(LCC2Y2,FCC2Y2(ISYCHO),YIM,IADR,LEN)
         ELSE IF (JSIDE .EQ. 4) THEN
            CALL PUTWA2(LCC2Y3,FCC2Y3(ISYCHO),YIM,IADR,LEN)
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &      'Error in ',SECNAM,': Unknown option, ISIDE = ',ISIDE
            CALL QUIT('Error in '//SECNAM)
         ENDIF
      ENDIF
      TIMT = SECOND() - TIMT
      TMCYIO(2) = TMCYIO(2) + TIMT

      RETURN
      END
C  /* Deck cc_cyird */
      SUBROUTINE CC_CYIRD(YIM,IOFA1,LVIRA,NX1AMA,IX1AMA,NUMV,JVEC1,
     &                    ISYCHO,ISYMY,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Read Y(#ai,#J) intermediates from disk.
C     
C     ISYCHO is the symmetry of J-index.
C     ISYMY  is symmetry of the Y intermediates (overall).
C     ISIDE  determines which files are used (according to
C            type of Y intermediates, i.e. LH trf., gr.state,
C            or RH trf.)
C
C     Note: Use CC_CYIRDF for full read without using "local" index
C           arrays.
C
#include <implicit.h>
      DIMENSION YIM(*)
      INTEGER   IOFA1(8), LVIRA(8), NX1AMA(8)
      INTEGER   IX1AMA(8,8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <chocc2lhs.h>
#include <chocc2rhs.h>
#include <chocc2ftr.h>
#include <priunit.h>
#include <cyit.h>

      INTEGER IADR, ICOL, LROW

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CYIRD')

      ISYMAI = MULD2H(ISYCHO,ISYMY)

      IF (NX1AMA(ISYMAI) .EQ. NT1AM(ISYMAI)) THEN

C        Full read.
C        ----------

         CALL CC_CYIRDF(YIM,NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)

      ELSE

C        Batch read.
C        -----------

         TIMT = SECOND()

         CALL CC_CYITRNSL(JSIDE,ISIDE)

         LROW = NT1AM(ISYMAI)

         DO IVEC = 1,NUMV

            JVEC = JVEC1 + IVEC - 1
            ICOL = JVEC  - 1

            DO ISYMI = 1,NSYM

               IF (NRHF(ISYMI) .GT. 0) THEN

                  ISYMA = MULD2H(ISYMI,ISYMAI)

                  IF (LVIRA(ISYMA) .EQ. NVIR(ISYMA)) THEN

                     KOFFY = NX1AMA(ISYMAI)*(IVEC - 1)
     &                     + IX1AMA(ISYMA,ISYMI) + 1
                     IADR  = LROW*ICOL + IT1AM(ISYMA,ISYMI) + 1
                     LEN   = NVIR(ISYMA)*NRHF(ISYMI)

                     IF (JSIDE .EQ. -1) THEN
                        CALL GETWA2(LCC2YL,FCC2YL(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 0) THEN
                        CALL GETWA2(LCC2YM,FCC2YM(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 1) THEN
                        CALL GETWA2(LCC2YR,FCC2YR(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 2) THEN
                        CALL GETWA2(LCC2Y1,FCC2Y1(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 3) THEN
                        CALL GETWA2(LCC2Y2,FCC2Y2(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 4) THEN
                        CALL GETWA2(LCC2Y3,FCC2Y3(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE
                        WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &                  'Error in ',SECNAM,': Unknown option, ISIDE = ',
     &                  ISIDE
                        CALL QUIT('Error in '//SECNAM)
                     ENDIF

                  ELSE IF (LVIRA(ISYMA) .GT. 0) THEN

                     LEN = LVIRA(ISYMA)

                     DO I = 1,NRHF(ISYMI)

                        KOFFY = NX1AMA(ISYMAI)*(IVEC - 1) 
     &                        + IX1AMA(ISYMA,ISYMI)
     &                        + LVIRA(ISYMA)*(I - 1) + 1
                        IADR  = LROW*ICOL + IT1AM(ISYMA,ISYMI)
     &                        + NVIR(ISYMA)*(I - 1) + IOFA1(ISYMA)


                        IF (JSIDE .EQ. -1) THEN
                           CALL GETWA2(LCC2YL,FCC2YL(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 0) THEN
                           CALL GETWA2(LCC2YM,FCC2YM(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 1) THEN
                           CALL GETWA2(LCC2YR,FCC2YR(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 2) THEN
                           CALL GETWA2(LCC2Y1,FCC2Y1(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 3) THEN
                           CALL GETWA2(LCC2Y2,FCC2Y2(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 4) THEN
                           CALL GETWA2(LCC2Y3,FCC2Y3(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE
                           WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &                     'Error in ',SECNAM,
     &                     ': Unknown option, ISIDE = ',ISIDE
                           CALL QUIT('Error in '//SECNAM)
                        ENDIF

                     ENDDO

                  ENDIF

               ENDIF

            ENDDO

         ENDDO

         TIMT = SECOND() - TIMT
         TMCYIO(2) = TMCYIO(2) + TIMT

      ENDIF

      RETURN
      END
C  /* Deck cc_cyiwr */
      SUBROUTINE CC_CYIWR(YIM,IOFA1,LVIRA,NX1AMA,IX1AMA,NUMV,JVEC1,
     &                    ISYCHO,ISYMY,ISIDE)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Write Y(#ai,#J) intermediates to disk.
C     
C     ISYCHO is the symmetry of J-index.
C     ISYMY  is symmetry of the Y intermediates (overall).
C     ISIDE  determines which files are used (according to
C            type of Y intermediates, i.e. LH trf., gr.state,
C            or RH trf.)
C
C     Note: Use CC_CYIWRF for full write without using "local" index
C           arrays.
C
#include <implicit.h>
      DIMENSION YIM(*)
      INTEGER   IOFA1(8), LVIRA(8), NX1AMA(8)
      INTEGER   IX1AMA(8,8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <chocc2lhs.h>
#include <chocc2rhs.h>
#include <chocc2ftr.h>
#include <priunit.h>
#include <cyit.h>

      INTEGER IADR, ICOL, LROW

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CYIWR')

      ISYMAI = MULD2H(ISYCHO,ISYMY)

      IF (NX1AMA(ISYMAI) .EQ. NT1AM(ISYMAI)) THEN

C        Full write.
C        -----------

         CALL CC_CYIWRF(YIM,NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)

      ELSE

C        Batch write.
C        ------------

         TIMT = SECOND()

         CALL CC_CYITRNSL(JSIDE,ISIDE)

         LROW = NT1AM(ISYMAI)

         DO IVEC = 1,NUMV

            JVEC = JVEC1 + IVEC - 1
            ICOL = JVEC  - 1

            DO ISYMI = 1,NSYM

               IF (NRHF(ISYMI) .GT. 0) THEN

                  ISYMA = MULD2H(ISYMI,ISYMAI)

                  IF (LVIRA(ISYMA) .EQ. NVIR(ISYMA)) THEN

                     KOFFY = NX1AMA(ISYMAI)*(IVEC - 1)
     &                     + IX1AMA(ISYMA,ISYMI) + 1
                     IADR  = LROW*ICOL + IT1AM(ISYMA,ISYMI) + 1
                     LEN   = NVIR(ISYMA)*NRHF(ISYMI)

                     IF (JSIDE .EQ. -1) THEN
                        CALL PUTWA2(LCC2YL,FCC2YL(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 0) THEN
                        CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 1) THEN
                        CALL PUTWA2(LCC2YR,FCC2YR(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 2) THEN
                        CALL PUTWA2(LCC2Y1,FCC2Y1(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 3) THEN
                        CALL PUTWA2(LCC2Y2,FCC2Y2(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE IF (JSIDE .EQ. 4) THEN
                        CALL PUTWA2(LCC2Y3,FCC2Y3(ISYCHO),YIM(KOFFY),
     &                               IADR,LEN)
                     ELSE
                        WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &                  'Error in ',SECNAM,': Unknown option, ISIDE = ',
     &                  ISIDE
                        CALL QUIT('Error in '//SECNAM)
                     ENDIF

                  ELSE IF (LVIRA(ISYMA) .GT. 0) THEN

                     LEN = LVIRA(ISYMA)

                     DO I = 1,NRHF(ISYMI)

                        KOFFY = NX1AMA(ISYMAI)*(IVEC - 1) 
     &                        + IX1AMA(ISYMA,ISYMI)
     &                        + LVIRA(ISYMA)*(I - 1) + 1
                        IADR  = LROW*ICOL + IT1AM(ISYMA,ISYMI)
     &                        + NVIR(ISYMA)*(I - 1) + IOFA1(ISYMA)

                        IF (JSIDE .EQ. -1) THEN
                           CALL PUTWA2(LCC2YL,FCC2YL(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 0) THEN
                           CALL PUTWA2(LCC2YM,FCC2YM(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 1) THEN
                           CALL PUTWA2(LCC2YR,FCC2YR(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 2) THEN
                           CALL PUTWA2(LCC2Y1,FCC2Y1(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 3) THEN
                           CALL PUTWA2(LCC2Y2,FCC2Y2(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE IF (JSIDE .EQ. 4) THEN
                           CALL PUTWA2(LCC2Y3,FCC2Y3(ISYCHO),
     &                                  YIM(KOFFY),IADR,LEN)
                        ELSE
                           WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &                     'Error in ',SECNAM,
     &                     ': Unknown option, ISIDE = ',ISIDE
                           CALL QUIT('Error in '//SECNAM)
                        ENDIF

                     ENDDO

                  ENDIF

               ENDIF

            ENDDO

         ENDDO

         TIMT = SECOND() - TIMT
         TMCYIO(2) = TMCYIO(2) + TIMT

      ENDIF

      RETURN
      END
C  /* Deck cc_cyini */
      SUBROUTINE CC_CYINI(WORK,LWORK,ISYMY,ISIDE,NUMCHO)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Initialize Y intermediates on disk.
C
C     NOTE: The number of Y intermediates must be passed in NUMCHO !!!
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      INTEGER   NUMCHO(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <cyit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_CYINI')

      PARAMETER (IOPEN = -1, IKEEP = 1)

C     Get a zero vector.
C     ------------------

      MAXT1 = NT1AM(1)
      DO ISYM = 2,NSYM
         MAXT1 = MAX(MAXT1,NT1AM(ISYM))
      ENDDO
      IF (MAXT1 .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',MAXT1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      CALL DZERO(WORK,MAXT1)

C     Write zero vectors to disk.
C     ---------------------------

      DO ISYCHO = 1,NSYM
         CALL CC_CYIOP(IOPEN,ISYCHO,ISIDE)
         DO JVEC = 1,NUMCHO(ISYCHO)
            CALL CC_CYIWRF(WORK,1,JVEC,ISYCHO,ISYMY,ISIDE)
         ENDDO
         CALL CC_CYIOP(IKEEP,ISYCHO,ISIDE)
      ENDDO

      RETURN
      END
C  /* Deck cho_imop */
      SUBROUTINE CHO_IMOP(IOPT,ITYP,LUNIT,ISYM)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose:
C
C        IOPT < 0: Open direct access file(s).
C        IOPT = 0: Close and delete the files.
C        IOPT > 0: Close and keep the files.
C
C        ITYP =  1: E(ij)    - ground state (global) intermediate.
C        ITYP =  2: E(ab)    - ground state (global) intermediate.
C        ITYP =  3: C(ai)    - C intermediates for left-hand Jacobian
C                              transformations.
C        ITYP =  4: ETA(ai)  - right-hand side for 0th order multiplier
C                              equations ('L0').
C        ITYP =  5: XI(ai)   - eff. right-hand side for 1st order multiplier
C                              equations ('R1').
C
C        LUNIT: Unit number.
C        ISYM : Symmetry of the intermediates.
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL OLDDX

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_IMOP')

      PARAMETER (NTYP = 5)

      CHARACTER*7 FILNAM(NTYP)
      DATA FILNAM /'CHO_EIJ','CHO_EAB','CHO_CIM','CHO_ETA',
     &             'CHO_XIE'/

C     Check ITYP.
C     -----------

      IF ((ITYP.LT.1) .OR. (ITYP.GT.NTYP)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,': specifier ITYP not recognized.'
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ITYP = ',ITYP
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      IF (IOPT .LT. 0) THEN

C        Open direct-access files.
C        -------------------------

         LRLEN = 2*NT1AM(ISYM)

         IF (ITYP .EQ. 1) THEN
            LRLEN = 2*NMATIJ(ISYM)
         ELSE IF (ITYP .EQ. 2) THEN
            LRLEN = 2*NMATAB(ISYM)
         ENDIF

         LRLEN = MAX(LRLEN,1)

         LUNIT = -1
         CALL GPOPEN(LUNIT,FILNAM(ITYP),'UNKNOWN','DIRECT',
     &               'UNFORMATTED',LRLEN,OLDDX)

      ELSE IF (IOPT .EQ. 0) THEN

C        Close and delete files.
C        -----------------------

         CALL GPCLOSE(LUNIT,'DELETE')
         LUNIT = -1

      ELSE

C        Close and keep files.
C        ---------------------

         CALL GPCLOSE(LUNIT,'KEEP')
         LUNIT = -1

      ENDIF

      RETURN
      END
C  /* Deck cho_wrrst */
      SUBROUTINE CHO_WRRST(ECC2,T1AM,OMEGA1)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Write Cholesky CC2 restart file.
C
#include <implicit.h>
      DIMENSION T1AM(*), OMEGA1(*)
#include <ccorb.h>
#include <ccsdsym.h>

      CHARACTER*10 FILRST
      PARAMETER (FILRST = 'CHOCC2.RST')

      LURST = -1
      CALL GPOPEN(LURST,FILRST,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURST)
      WRITE(LURST) NT1AMX
      WRITE(LURST) ECC2
      WRITE(LURST) (T1AM(I),   I = 1,NT1AMX)
      WRITE(LURST) (OMEGA1(I), I = 1,NT1AMX)
      CALL GPCLOSE(LURST,'KEEP')

      RETURN
      END
C  /* Deck cho_rdrst */
      SUBROUTINE CHO_RDRST(ECC2,T1AM,OMEGA1,RDE,RDT,RDO,IFAIL)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Read Cholesky CC2 restart file (MUST exist on entry).
C
C     RDE = .TRUE. : Read CC2 energy, ECC2
C     RDT = .TRUE. : Read CC2 singles amplitudes, T1AM
C     RDO = .TRUE. : Read CC2 singles vector function, OMEGA1
C
C     On entry:
C
C        IFAIL = -1: stop if inconsistent data lengths (see below).
C        IFAIL anything else: set IFAIL flag and return (see below).
C
C     On exit:
C
C        IFAIL = 0: everything seems OK. T1AM and OMEGA1 read from file.
C        IFAIL = 1: there is a problem with the amount of data on file:
C                   The length of T1AM and OMEGA1 is not equal to NT1AMX.
C                   Nothing is read in this case: T1AM and OMEGA1 undefined!
C
#include <implicit.h>
      DIMENSION T1AM(*), OMEGA1(*)
      LOGICAL RDE, RDT, RDO
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_RDRST')

      CHARACTER*10 FILRST
      PARAMETER (FILRST = 'CHOCC2.RST')

      LOGICAL ERRSTP

C     Set IFAIL.
C     ----------

      ERRSTP = IFAIL .EQ. -1
      IFAIL  = 0

C     Open file: note that the file MUST exist.
C     -----------------------------------------

      LURST = -1
      CALL GPOPEN(LURST,FILRST,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LURST)

C     Check that the right amount of data is on file.
C     -----------------------------------------------

      READ(LURST) LENTOT

      IF (LENTOT .NE. NT1AMX) THEN
         IF (ERRSTP) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,A)')
     &      SECNAM,': file ',FILRST,' contains incorrect data!'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'From file: LENTOT = ',LENTOT,
     &      'Expected : NT1AMX = ',NT1AMX
            CALL QUIT('Error in '//SECNAM)
         ELSE
            IFAIL = 1
            CALL GPCLOSE(LURST,'KEEP')
            RETURN
         ENDIF
      ENDIF

C     ECC2.
C     -----

      IF (RDE) THEN
         READ(LURST) ECC2
      ELSE
         READ(LURST)
      ENDIF

C     T1AM.
C     -----

      IF (RDT) THEN
         READ(LURST) (T1AM(I), I = 1,NT1AMX)
      ELSE
         READ(LURST)
      ENDIF

C     OMEGA1.
C     -------

      IF (RDO) THEN
         READ(LURST) (OMEGA1(I), I = 1,NT1AMX)
      ENDIF

C     Close file.
C     -----------

      CALL GPCLOSE(LURST,'KEEP')

      RETURN
      END
C  /* Deck cho_rdsir */
      SUBROUTINE CHO_RDSIR(POTNUC,ESCF,FOCKD,CMO,WORK,LWORK,RDE,RDF,RDC)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Read information from SIRIUS interface file SIRIFC:
C
C              POTNUC - Nuclear potential energy.
C              ESCF   - SCF energy.
C              FOCKD  - orbital energies, CC ordering, frozen orbitals deleted.
C              CMO    - MO coefficients,  CC ordering, frozen orbitals deleted.
C
C     RDE = .TRUE. : read POTNUC and ESCF
C     RDF = .TRUE. : read FOCKD
C     RDC = .TRUE. : read CMO
C
#include <implicit.h>
      DIMENSION CMO(*), FOCKD(*), WORK(LWORK)
      LOGICAL RDE, RDF, RDC
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*6 FSIRFC
      PARAMETER (FSIRFC = 'SIRIFC')

C     Open file.
C     ----------

      LUSIRG = -1
      CALL GPOPEN(LUSIRG,FSIRFC,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUSIRG)

C     ESCF (among other stuff).
C     -------------------------

      IF (RDE) THEN
         CALL MOLLAB('SIR IPH ',LUSIRG,LUPRI)
         READ(LUSIRG) POTNUC,EMY,EACTIV,ESCF,ISTATE,ISPIN,NACTEL,LSYM
      ENDIF

C     Position file.
C     ----------------

      CALL MOLLAB('TRCCINT ',LUSIRG,LUPRI)
      READ(LUSIRG)

C     Orbital energies.
C     -----------------

      IF (RDF) THEN
         READ(LUSIRG) (FOCKD(I), I = 1,NORBTS)
      ELSE
         READ(LUSIRG)
      ENDIF

C     MO coefficients.
C     ----------------

      IF (RDC) THEN
         READ(LUSIRG) (CMO(I),I=1,NLAMDS)
      ENDIF
     

C     Close file.
C     -----------

      CALL GPCLOSE(LUSIRG,'KEEP')

C     Reorder.
C     --------

      IF (RDF) THEN
         IF (FROIMP .OR. FROEXP) CALL CCSD_DELFRO(FOCKD,WORK,LWORK)
         CALL FOCK_REORDER(FOCKD,WORK,LWORK)
      ENDIF

      IF (RDC) THEN
         CALL CMO_REORDER(CMO,WORK,LWORK)
      ENDIF

      RETURN
      END
C  /* Deck rd_decmos */
      SUBROUTINE RD_DECMOS(THRDCM,SPADCM,MXQUAL,MXREAD,NTOVEC,NCOLUM,
     &                     XDIANL,NSYM,THRZER,FILSAV,SECALL)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Read information for MO decomposition.
C
#include <implicit.h>
#include <priunit.h>

      INTEGER NTOVEC(NSYM), NCOLUM(NSYM)
      DIMENSION XDIANL(3,NSYM)
      CHARACTER*(*) FILSAV, SECALL

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'RD_DECMOS')

      LUSAV = -1
      CALL GPOPEN(LUSAV,FILSAV,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUSAV)
      READ(LUSAV) LSYM
      IF (LSYM .NE. NSYM) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,A,/,5X,A,I10,/,5X,A,I10,/)')
     &   'Error detected in ',SECNAM,':',
     &   'Symmetry mismatch for skip option in ',SECALL,
     &   'Number of irreps read from disk: ',LSYM,
     &   'Actual number of irreps is     : ',NSYM
         CALL QUIT('Symmetry mismatch for skip option in '//SECNAM)
      ENDIF
      READ(LUSAV) THRDCM
      READ(LUSAV) SPADCM
      READ(LUSAV) MXQUAL
      READ(LUSAV) MXREAD
      READ(LUSAV) THRZER
      READ(LUSAV) (NTOVEC(ISYM),ISYM=1,NSYM)
      READ(LUSAV) (NCOLUM(ISYM),ISYM=1,NSYM)
      DO ISYM = 1,NSYM
         READ(LUSAV) (XDIANL(I,ISYM),I=1,3)
      ENDDO
      CALL GPCLOSE(LUSAV,'KEEP')

      RETURN
      END
C  /* Deck wr_decmos */
      SUBROUTINE WR_DECMOS(THRDCM,SPADCM,MXQUAL,MXREAD,NTOVEC,NCOLUM,
     &                     XDIANL,NSYM,THRZER,FILSAV)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Save information for MO decomposition.
C
#include <implicit.h>
      INTEGER NTOVEC(NSYM), NCOLUM(NSYM)
      DIMENSION XDIANL(3,NSYM)
      CHARACTER*(*) FILSAV

      LUSAV = -1
      CALL GPOPEN(LUSAV,FILSAV,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUSAV)
      WRITE(LUSAV) NSYM
      WRITE(LUSAV) THRDCM
      WRITE(LUSAV) SPADCM
      WRITE(LUSAV) MXQUAL
      WRITE(LUSAV) MXREAD
      WRITE(LUSAV) THRZER
      WRITE(LUSAV) (NTOVEC(ISYM),ISYM=1,NSYM)
      WRITE(LUSAV) (NCOLUM(ISYM),ISYM=1,NSYM)
      DO ISYM = 1,NSYM
         WRITE(LUSAV) (XDIANL(I,ISYM),I=1,3)
      ENDDO
      CALL FLSHFO(LUSAV)
      CALL GPCLOSE(LUSAV,'KEEP')

      RETURN
      END
