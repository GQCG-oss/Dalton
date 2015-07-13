C  /* Deck so_eres */
      SUBROUTINE SO_ERES(MODEL,  NOLDTR, NNEWTR,  DENSIJ,  LDENSIJ, 
     &                   DENSAB, LDENSAB, T2MP,    LT2MP,   FOCKD,
     &                   LFOCKD, DENSAI,  LDENSAI, NIT,     ISYMTR,
     &                   WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Frederik Beyer & Stephan P. A. Sauer: 27.08.2013: call to ERI corrected
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian matricx E[2]. 
C              The trial vector consists of four parts TR1E, TR1D, 
C              TR2E, and TR2D. E refers to excitations and D to 
C              de-excitations. 1 refer to the one-particle part and
C              2 to the two-particle part. The linear transformed 
C              trialvector is refered to as the resultvector and is
C              kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D, RES2E, and RES2D.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] is not constructed explicitly.
C

#include "implicit.h"

#ifdef VAR_MPI
#include "mpif.h"
#endif

#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "eribuf.h"

C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION FOCKD(LFOCKD)   
      DIMENSION DENSAI(LDENSAI) !intent(out) 
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*5
      integer :: inewtr, nnewtr, isymd1, ntosym
C   
      integer, save :: inew=0
C
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "soppinf.h"

#ifdef VAR_MPI
#include "iprtyp.h"
#include "infpar.h"
#include "parsoppa.h"
      integer, save :: numprocs
      forceupdate = .true.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_ERES')
C
C
C------------------------------------------------------------------
C     Determine the symmetri of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO   = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO  + LCMO
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_ERES.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ERES.1',' ',KEND1,LWORK)
C
C-------------------------------------------------------
C     Get the matrix which contains the MO coefficients.
C-------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(1)  = SOTIME(1) + DTIME
C
C---------------------------------
C     Work space allocation no. 2.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRES1E  = NT1AM(ISYMTR)
      LRES1D  = NT1AM(ISYMTR)
      LFOCK   = N2BST(ISYRES)
      LDENS   = N2BST(ISYMTR)
      LBTR1E  = NT1AO(ISYMTR)
      LBTR1D  = NT1AO(ISYMTR)
      LBTJ1E  = NMATAV(ISYMTR)
      LBTJ1D  = NMATAV(ISYMTR)

      LTR2E   = NT2AM(ISYMTR)
      LTR2D   = NT2AM(ISYMTR)
      LRES2E  = NT2AM(ISYMTR)
      LRES2D  = NT2AM(ISYMTR)
      LSIGAI1 = NT1AO(ISYRES)
      LSIGAI2 = NT1AO(ISYRES)
      LSIGDA1 = NMATAV(ISYRES)
      LSIGDA2 = NMATAV(ISYRES)
      LAIJ    = NRHFT*NRHFT
      LAAB    = NVIRT*NVIRT
C
      KTR1E   = KEND1
      KTR1D   = KTR1E   + LTR1E
      KTR2E   = KTR1D   + LTR1D
      KTR2D   = KTR2E   + LTR2E

      KRES1E  = KTR2D   + LTR2D
      KRES1D  = KRES1E  + LRES1E
      KRES2E  = KRES1D  + LRES1D
      KRES2D  = KRES2E  + LRES2E

      KFOCK   = KRES2D  + LRES2D
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E

      KSIGAI1 = KBTJ1D  + LBTJ1D
      KSIGAI2 = KSIGAI1 + LSIGAI1
      KSIGDA1 = KSIGAI2 + LSIGAI2
      KSIGDA2 = KSIGDA1 + LSIGDA1

      KAIJ    = KSIGDA2 + LSIGDA2
      KAAB    = KAIJ    + LAIJ
      KEND2   = KAAB    + LAAB
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('SO_ERES.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_ERES.2',' ',KEND2,LWORK)
C
C----------------------------
C     Initialize AIJ and AAB.
C----------------------------
C
      CALL DZERO(WORK(KAIJ),LAIJ)
      CALL DZERO(WORK(KAAB),LAAB)
C
C----------------------------------------------
C     Open files with trial and result vectors.
C----------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
C
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
      CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
      CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
C

      IF ( IPRSOP. GE. 7 ) THEN !Only printing related.
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
         DO 50 INEWTR = 1,NNEWTR
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
            INEW = NOLDTR + INEWTR
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,INEW)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
            CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,INEW)
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,'. new trial vector'
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
   50    CONTINUE
      END IF
C
C================================================
C     Loop over number of excitations considered.
C================================================
C
      DO 100 INEWTR = 1,NNEWTR
C
C-------------------------------------------------
C        Determine pointer to INEWTR trial vector.
C-------------------------------------------------
C
         INEW = NOLDTR + INEWTR
C
C--------------------------------------------------------------
C        Initialize RES1E, RES1D, RES2E, RES2D, SIGAI1, SIGAI2,
C                   SIGDA1, SIGDA2 and FOCK
C--------------------------------------------------------------
C
         CALL DZERO(WORK(KRES1E),LRES1E)
         CALL DZERO(WORK(KRES1D),LRES1D)
         CALL DZERO(WORK(KRES2E),LRES2E)
         CALL DZERO(WORK(KRES2D),LRES2D)
C
         CALL DZERO(WORK(KSIGAI1),LSIGAI1)
         CALL DZERO(WORK(KSIGAI2),LSIGAI2)
         CALL DZERO(WORK(KSIGDA1),LSIGDA1)
         CALL DZERO(WORK(KSIGDA2),LSIGDA2)
C
         CALL DZERO(WORK(KFOCK),LFOCK)
C
C--------------------------
C        Read trial vector.
C--------------------------
C
         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
         CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,INEW)
         CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,INEW)
C
C---------------------------------------------------
C        Calculate RPA-density matrices in AO basis.
C---------------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_AODENS(WORK(KDENS),LDENS,WORK(KCMO),LCMO,
     &                  WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,ISYMTR,
     &                  WORK(KEND2),LWORK2)
         DTIME     = SECOND()  - DTIME
         SOTIME(6) = SOTIME(6) + DTIME
C
C--------------------------------------------
C        Backtransformation of trial vectors.
C--------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_BCKTR(WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,WORK(KBTR1E),
     &                 LBTR1E,WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                 WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME
C
C=======================================================
C        Start the loop over distributions of integrals.
C=======================================================
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME  = SECOND()
            IF (HERDIR) THEN
               CALL HERDI1(WORK(KEND2),LWRK2,IPRINT)
            ELSE
               KCCFB1 = KEND2
               KINDXB = KCCFB1 + MXPRIM*MXCONT
               KEND2  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
               LWORK2 = LWORK  - KEND2

               CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                     KODPP1,KODPP2,KRDPP1,KRDPP2,KFREE,LFREE,
     &                     KEND2,WORK(KCCFB1),WORK(KINDXB),WORK(KEND2),
     &                     LWORK2,IPRINT)

               KEND2  = KFREE
               LWORK2 = LFREE
               DTIME     = SECOND()  - DTIME
               SOTIME(8) = SOTIME(8) + DTIME
            ENDIF
         ELSE
            NTOSYM = NSYM
         ENDIF
C
         KENDSV  = KEND2
         LWORKSV = LWORK2
         ICDEL1  = 0
C
         DO 210 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
               NTOT = MXCALL
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
C
#ifdef VAR_MPI
C***********************************************************
C     If more than 1 process is available, calculate 
C     in parallel with par_so_eres.
C
C     Note, there is a strong coherence between the order of 
C     these bcasts and the initialization of the slaves in 
C     the beginning of par_so_eres. 
C***********************************************************
      thisinewtr = inewtr
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      parallelsoppa: if (numprocs.gt.1) then

         ! make every slave join the parallel calculation
         ! Signal for PARA_SO_ERES=87 in iprtyp.h, maybe insert macro
         ! directly 
         call mpixbcast(87, 1, 'INTEGE', 0)
         ! Second broadcast -> printlevel for the slaves
         call mpixbcast( 0, 1, 'INTEGE', 0)

         ! Ready slaves for massive amount of updates.
         call mpi_bcast(forceupdate,1,mpi_logical,0,mpi_comm_world,ierr)

         ! These are input parameters to so_eres. To copy them to par_so_eres I rename then temporarily and store then in parsoppa.h
         ! Once inside par_so_eres the parameters are copied back from the common block to the master process.
         copymodel   = model
         copyldensai = ldensai
         copynit     = nit
         copyisymtr  = isymtr
         copyinewtr  = inewtr

         ! Copy the common block to all slaves
         call getbytespan(ntot, parsoppaLAST, bytesize)
         call mpi_bcast(ntot, bytesize, mpi_byte, 0,mpi_comm_world,ierr)

         ! Transfer mp2 amplitudes to all slaves
         call mpi_bcast(lt2mp,1, mpi_integer, 0,mpi_comm_world, ierr)
         call mpi_bcast(t2mp, lt2mp, mpi_real8, 0,mpi_comm_world,ierr)
         
         ! join the parallel calculation along w. slaves
         call par_so_eres(WORK, LWORK, 0)

         ! receive the newly assembled densai array from process 1 before continuing
         if ( getdensai )  then
            call mpi_recv(densai, ldensai, mpi_real8, 1, 0, 
     &      mpi_comm_world, mpistatus, ierr)
            getdensai = .false.
         endif

         goto 888  ! Skips over the do loop that would normally assemble the Fock/density matrix sequentially

      endif parallelsoppa
#endif !!endif VAR_MPI

            DO 220 ILLL = 1,NTOT
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
               IF (DIRECT) THEN
C
                  KEND2  = KENDSV
                  LWORK2 = LWORKSV
C
                  DTIME  = SECOND()
                  IF (HERDIR) THEN
                    CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
                  ELSE
C
                     CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                           WORK(KODCL1),WORK(KODCL2),
     &                           WORK(KODBC1),WORK(KODBC2),
     &                           WORK(KRDBC1),WORK(KRDBC2),
     &                           WORK(KODPP1),WORK(KODPP2),
     &                           WORK(KRDPP1),WORK(KRDPP2),
     &                           WORK(KCCFB1),WORK(KINDXB), 
     &                           WORK(KEND2),LWORK2,IPRINT)
C
                     DTIME     = SECOND()  - DTIME
                     SOTIME(9) = SOTIME(9) + DTIME
                  ENDIF
C
                  LRECNR = ( (NBUFX(0) -1) / IRAT ) + 1
                  KRECNR  = KEND2
                  KEND2   = KRECNR + LRECNR
                  LWORK2  = LWORK  - KEND2
C
                  CALL SO_MEMMAX ('SO_ERES.2B',LWORK2)
                  IF (LWORK2 .LT. 0) 
     &                CALL STOPIT('SO_ERES.2B',' ',KEND2,LWORK)
C
               ELSE
                  NUMDIS = 1
               ENDIF
C
C-------------------------------------------------------------------
C   Loop over number of distributions in disk.
C   In the case of ERI there are more than one distribution and IDEL2
C   loops over them and the actual index of the delta orbital IDEL is 
C   then obtained from the array INDEXA. In the case of a not direct 
C   calculation there is only one distribution on the disk, which 
C   implies that IDEL2 is always 1 and that IDEL is systematically 
C   incremented by one each time.
C--------------------------------------------------------------------
C
               DO 230 IDEL2 = 1,NUMDIS
C

                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL) !keeps track of current symmetry
                  ELSE
                     IDEL  = IBAS(ISYMD1) + ILLL
                     ISYMD = ISYMD1
                  ENDIF
C
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
C
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)
C
C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
C
                  LXINT  = NDISAO(ISYDIS)
C
                  KXINT   = KEND2
                  KEND3   = KXINT + LXINT
                  LWORK3  = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &                CALL STOPIT('SO_ERES.3',' ',KEND3,LWORK)
C
C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
C
                  DTIME      = SECOND()
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(10) = SOTIME(10) + DTIME
C
C---------------------------------------------
C                 Work space allocation no. 4.
C---------------------------------------------
C
                  ISAIJ = MULD2H(ISYMD,1)
C
                  IF (NVIR(ISYMD) .GT. 0) THEN
                     LT2M1 = NT2BCD(ISAIJ)
                  ELSE
                     LT2M1 = 0
                  END IF
C
                  KT2M1   = KEND3
                  KEND4   = KT2M1  + LT2M1
                  LWORK4  = LWORK  - KEND4
C
                  CALL SO_MEMMAX ('SO_ERES.4',LWORK4)
                  IF (LWORK4 .LT. 0) 
     &                CALL STOPIT('SO_ERES.4',' ',KEND4,LWORK)
C
C------------------------------------------------------------
C                 Construct the partially back-transformed T2
C                 MP-amplitudes.
C------------------------------------------------------------
C
                  DTIME      = SECOND()
                  CALL SO_T2M1(WORK(KT2M1),LT2M1,T2MP,LT2MP,
     &                         WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                         WORK(KEND4),LWORK4)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(12) = SOTIME(12) + DTIME
C
C---------------------------------------------
C                 Work space allocation no. 5.
C---------------------------------------------
C
                  LDSRHF = NDSRHF(ISYMD)
C
                  KDSRHF  = KEND4
                  KEND5   = KDSRHF + LDSRHF
                  LWORK5  = LWORK  - KEND5
C
                  CALL SO_MEMMAX ('SO_ERES.5',LWORK5)
                  IF (LWORK5 .LT. 0) 
     &                CALL STOPIT('SO_ERES.5',' ',KEND5,LWORK)
C
C----------------------------------------------------------------
C                 Transform one index in the integral batch to an
C                 occupied index.
C----------------------------------------------------------------
C
                  DTIME      = SECOND()
                  ISYMLP = 1
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                        ISYMLP,WORK(KEND5),LWORK5,ISYDIS)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(13) = SOTIME(13) + DTIME
C
C-------------------------------------------------------------------
C                 Calculate part of the second order density matrix.
C-------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  IF ( (MODEL.NE.'AOSOC') .AND. (NIT.EQ.1) .AND. 
     &               (  INEWTR.EQ.1) ) THEN
                     CALL SO_DENSAI1(DENSAI,LDENSAI,WORK(KDSRHF),LDSRHF,
     &                               WORK(KCMO),LCMO,WORK(KT2M1),LT2M1,
     &                               ISYMD,ISYDIS,WORK(KEND5),
     &                               LWORK5)

                  END IF
                  DTIME      = SECOND()   - DTIME
                  SOTIME(41) = SOTIME(41) + DTIME
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  DTIME      = SECOND()
                  CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                           WORK(KEND5),LWORK5,IDEL,ISYMD,.FALSE.,
     &                           'crashifTrue',ISYMTR)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(11) = SOTIME(11) + DTIME
C
C----------------------------------------------------------------------
C                 Calculate part of the result vectors RES1E and RES1D,
C                 specifically the first and the second term in eqs.
C                 (34,35). Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C


                  DTIME      = SECOND()
                  CALL SO_RES_A(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                          WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                          WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                          WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(14) = SOTIME(14) + DTIME
C
C-------------------------------------------------------------------
C                 Calculate the part of the result vectors RES1E and
C                 RES1D which originate from the C matrices. See 
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  CALL SO_RES_TCB(WORK(KRES1E),LRES1E,WORK(KRES1D),
     &                            LRES1D,WORK(KTR2E),LTR2E,
     &                            WORK(KTR2D),LTR2D,WORK(KDSRHF),LDSRHF,
     &                            WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                            ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(29) = SOTIME(29) + DTIME
C
C----------------------------------------------------------------------
C                 Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  CALL SO_RES_CB(WORK(KRES2E),LRES2E,WORK(KRES2D),
     &                          LRES2D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KBTR1E),LBTR1E,
     &                          WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                          WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                          IDEL,ISYMD,ISYDIS,ISYMTR,WORK(KEND5),
     &                          LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(15) = SOTIME(15) + DTIME
C
C--------------------------------------------------------------------
C                 Construct SIGMAI1(ALFA,I) and SIGMAI2(ALFA,I) which
C                 are used in SO_RES_B.
C--------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  CALL SO_SIGAI(WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),
     &                          LSIGAI2,WORK(KT2M1),LT2M1,WORK(KXINT),
     &                          LXINT,WORK(KBTR1E),LBTR1E,WORK(KBTR1D),
     &                          LBTR1D,WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),
     &                          LBTJ1D,WORK(KCMO),LCMO,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(16) = SOTIME(16) + DTIME
C
C--------------------------------------------------------------------
C                 Construct SIGDA1(DELTA,A) and SIGDA2(DELTA,A) which
C                 are used SO_RES_C.
C--------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  CALL SO_SIGDA(WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),
     &                          LSIGDA2,T2MP,LT2MP,WORK(KDSRHF),LDSRHF,
     &                          WORK(KBTR1E),LBTR1E,WORK(KBTR1D),LBTR1D,
     &                          WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),LBTJ1D,
     &                          WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(17) = SOTIME(17) + DTIME
C
  230          CONTINUE  ! End of IDEL2 loop
C
  220       CONTINUE !End of ILLL loop
C
  210    CONTINUE !end of ISYMD1 loop
C
  888    continue !Continuation point in case of parallel run.


C====================================================
C        End of loop over distributions of integrals.
C====================================================
C
C---------------------------------------------
C        Transform AO Fock matrix to MO basis.
C---------------------------------------------
C
         DTIME      = SECOND()
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     &                    WORK(KEND2),LWORK2,ISYRES,1,1)
         DTIME      = SECOND()   - DTIME
         SOTIME(24) = SOTIME(24) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA two-particle parts to the result 
C        vectors.
C------------------------------------------------------------------
C

         DTIME      = SECOND()
         CALL SO_TWOFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KFOCK),LFOCK,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(25) = SOTIME(25) + DTIME
C
C----------------------------------------------------------------------
C        Add contribution from sigma1(alfa,i) to RES1E and from 
C        sigma2(alfa,i) to RES1D. I.e. the third terms in eqs. (34) and 
C        (35).
C----------------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_B(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                 WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),LSIGAI2,
     &                 WORK(KCMO),LCMO,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(18) = SOTIME(18) + DTIME
C
C--------------------------------------------------------------------
C        Add contribution from sigda1(delta,a) to RES1E and from
C        sigda2(delta,a) to RES1D. I.e. the fourth terms in eqs. (34)
C        and (35).
C--------------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_C(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                 WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),LSIGDA2,
     &                 WORK(KCMO),LCMO,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(19) = SOTIME(19) + DTIME

C--------------------------------------------------------------
C        Calculate and add the symmetry correcting term to A in
C        eq. (44).
C--------------------------------------------------------------
C
         DTIME      = SECOND()

         CALL SO_RES_SYM(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KAIJ),LAIJ,WORK(KAAB),LAAB,WORK(KTR1E),
     &                   LTR1E,WORK(KTR1D),LTR1D,ISYRES)



         DTIME      = SECOND()   - DTIME
         SOTIME(20) = SOTIME(20) + DTIME
C
C---------------------------------------------------------
C        Calculate and add the Fock-term to A in eq. (40).
C---------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_FCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KTR1E),LTR1E,WORK(KTR1D),
     &                   LTR1D,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,
     &                   LDENSAB,ISYRES,ISYMTR)

         DTIME      = SECOND()   - DTIME
         SOTIME(21) = SOTIME(21) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA one-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()

         CALL SO_ONEFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,FOCKD,
     &                   LFOCKD,WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                   ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(26) = SOTIME(26) + DTIME
C
C-----------------------------------------------------------------
C        Construct D-contribution to 2p2h result vectors RES2E and
C        RES2D.
C-----------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_CD(WORK(KRES2E),LRES2E,WORK(KRES2D),LRES2D,
     &                  WORK(KTR2E),LTR2E,WORK(KTR2D),LTR2D,
     &                  FOCKD,LFOCKD,ISYRES,WORK(KEND2),LWORK2)
         
         DTIME      = SECOND()   - DTIME
         SOTIME(30) = SOTIME(30) + DTIME

C
C----------------------------------------
C        Write new resultvectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
         CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
         CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,INEW)
         CALL SO_WRITE(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,INEW)

C
  100 CONTINUE 
C
C==================================
C     End of loop over excitations.
C==================================
C
C----------------------------------------------------------------
C     Calculate the last part of the second order density matrix.
C----------------------------------------------------------------
C
      DTIME      = SECOND()
      IF ( (MODEL .NE. 'AOSOC')   .AND.
     &     (NIT .EQ. 1)         )  THEN
         CALL SO_DENSAI2(DENSAI,LDENSAI,FOCKD,LFOCKD)
      END IF
C
      DTIME      = SECOND()   - DTIME
      SOTIME(41) = SOTIME(41) + DTIME
C
      IF ( IPRSOP. GE. 7 ) THEN
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
         DO 400 INEWTR = 1,NNEWTR
            INEW = NOLDTR + INEWTR
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &                '. new E[2] linear transformed trial vector'
            CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
            CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS1D,INEW)
            CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS2E,INEW)
            CALL SO_READ(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,INEW)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LRES1E)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KRES2E+I-1),WORK(KRES2D+I-1),I=1,LRES2E)
  400    CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C   
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
C
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
C
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ERES')
C
      RETURN
      END SUBROUTINE

#ifdef VAR_MPI
C   /*deck par_so_eres */
      subroutine PAR_SO_ERES(WORK, LWORK, MYID)
C
C PAR_SO_ERES is a Parallel subroutine for SO_ERES
C It calculates the linear transformed trial vectors in parallel by distributing the 
C two-electron integrals over several processes and assembling the 
C partial matrices from every process before returning.
C
C WORK is the array where master and every slave assembles the fock matrix 
C incrementally. This is the 'output' of this subroutine.
C
C The load balancing is done by examining the number of AOs associated 
C with every distribution index in the two-electron integrals, before the integrals are 
C performed. The total load is calculated and the work is spread as evenly as 
C possible among the master and all slaves. The AO indices with the biggest 
C associated load of distributions is assigned first; then the smaller indices 
C are assigned. During the first iteration the actual walltimes are recorded and
C the calculation is rebalanced for the final time.
C
C A large number of common blocks need to be transferred from the master to 
C the slaves by using the subroutine getbytespan and mpi_bcast. The frequency
C with which these common blocks need to be transferred is regulated with a logical
C keyword, set in so_eres. 
C
      use dyn_iadrpk

#ifdef VAR_IFORT
      use IFPORT, ONLY: SLEEPQQ
#endif
! Parameter-only include files
#include "implicit.h"
#include "mpif.h" !use mpi fails with some profilers...
#include "priunit.h" 
#include "maxorb.h" 
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "iprtyp.h" 
#include "maxaqn.h" 
#include "chrnos.h" 
#include "ibtpar.h"

!#include "ibtfun.h"

! These include files depend on previous include files
#include "infpar.h" 
#include "eribuf.h" 

! The rest...
#include "inftap.h"
#include "ccorb.h"
#include "infind.h" 
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h" 
#include "soppinf.h"
#include "parsoppa.h"
#include "aobtch.h"
#include "odclss.h"
#include "ccom.h" 
#include "ericom.h" 
#include "eridst.h" 
#include "erithr.h"   
#include "erimem.h"  
#include "odbtch.h"
#include "nuclei.h"
#include "symmet.h"
#include "r12int.h"
#include "hertop.h"
#include "cbirea.h"
#include "erisel.h"
#include "symsq.h"
#include "gnrinf.h"
#include "ccpack.h"
#include "ccinftap.h"

C     Intent of Dummy arguments
      integer, intent(in) :: myid, lwork
      dimension work(lwork)  !intent(inout)

C     Variables for updating the slaves' /common/ blocks.
      integer :: bytesize
      logical :: resetvectors

C     MPI related and load balancing variables
      integer :: numprocs, request, allocsize
      integer, dimension(ntot) :: temparray
      integer, dimension(MPI_STATUS_SIZE) :: mpistatus

C     Pre-sorting load balancing variables
      logical, save :: forceloadbalance        =.true.
      logical, save :: forcedynamicloadbalance =.false.
      logical, save :: slavesupdatetasks       =.false.
      integer, save :: maxnumjobs
      logical       :: assignnow, workreceived
      integer       :: getnumjobs, getindices
      integer, allocatable, dimension(:), save :: AssignedIndices

C     Timing variables used for load monitoring
      integer :: color, workindex
      integer :: counti, countf, count_rate
      real    :: timespent
      real, dimension(:), allocatable, save :: timings
      
C     Miscellaneous
      DIMENSION INDEXA(MXCORB)
C     The real*8 datatype is due to the defined implicit type in implicit.h
C     If this datatype changes, every instance of real8 in this subroutine should
C     be changed to the correct type.
      real*8, allocatable, dimension(:) :: t2mp
      real*8, allocatable, dimension(:) :: densai
      real*8, allocatable, dimension(:) :: temp
      real*8, allocatable, dimension(:) :: tempdens

      logical       :: aijreduce        !This array is only calculated under special circumstances. This variable regulates when it needs to be transferred back to so_eres
      integer, save :: inewtr
      integer       :: allocstatus, deallocstatus, isymtr, ldensai, nit

      CHARACTER MODEL*5
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)


C**************************************************
C     Initialization variables.
C**************************************************

      aijreduce = .false.
      ! This bcast determines how many variables and common blocks
      ! need to be updated before the parallel calculation can start.
      if (myid.ne.0) then
         call mpi_bcast(forceupdate,1,mpi_logical,0,mpi_comm_world,ierr)

         ! This common block contains copies of the dummy parameters 
         ! passed to SO_ERES to I need to copy those values back before 
         ! continuing
         call getbytespan(ntot, parsoppaLAST, bytesize)
         call mpi_bcast(ntot,bytesize,mpi_byte, 0, mpi_comm_world, ierr)
      endif


      if (inewtr.ne.copyinewtr) then
         resetvectors = .true.
      else
         resetvectors = .false.
      endif

      model   = copymodel
      ldensai = copyldensai
      nit     = copynit
      isymtr  = copyisymtr
      inewtr  = copyinewtr

C     The array with MP2 amplitudes is passed to SO_ERES as a dummy parameter so in order to move the array to the parallel routine the master has to broadcast the length and contents to all slaves before the master joins the parallel calculation. Then when the master joins, one of the slaves transfer the contents of the array to the master.
      if (myid.ne.0) then

         ! Get length of array with MP2 amplitudes
         call mpi_bcast(lt2mp, 1, mpi_integer, 0, mpi_comm_world, ierr)
         allocate(t2mp(lt2mp), stat=allocstatus)
         if (allocstatus.ne.0)  then
            call quit('Allocation error in par_so_eres (t2mp)')
         endif

         ! Get the contents of the array
         ! NOTE: The master cannot enter PAR_SO_ERES before this bcast
         ! completes.
         call mpi_bcast(t2mp, lt2mp, mpi_real8, 0, mpi_comm_world, ierr)
      endif


      ! let the first slave transfer the info back to the master
      if (myid.eq.1) then
         call mpi_send(lt2mp, 1, mpi_integer,0,0,mpi_comm_world,ierr) 
         call mpi_send(t2mp, lt2mp, mpi_real8,0,0,mpi_comm_world,ierr) 
      elseif (myid.eq.0) then
         call mpi_recv(lt2mp, 1, mpi_integer, 1, 0, mpi_comm_world, 
     &                 mpistatus, ierr)
         allocate(t2mp(lt2mp), stat=allocstatus) 
         if (allocstatus.ne.0)  then
            call quit('Allocation error in par_so_eres (t2mp)')
         endif
         call mpi_recv(t2mp, lt2mp, mpi_real8, 1, 0, mpi_comm_world, 
     &   mpistatus, ierr)
      else
         continue
      endif 


C**********************************************************
C     Update all necessary common blocks
C     This large sequence of bcasts will only execute
C     once per call to RP_ERES
C
C     The <name> of the common block can be read from the <name>LAST
C     variable in the getbytespan call.
C**********************************************************
      commonblocktransfers: if (forceupdate) then

      aijreduce = .true.

      call getbytespan(lbuf, eribufLAST, bytesize)
      call mpi_bcast(lbuf, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(luaorc, eritapLAST, bytesize)
      call mpi_bcast(luaorc, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(nsym, ccorbLAST, bytesize)
      call mpi_bcast(nsym, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(irow, infindLAST, bytesize)
      call mpi_bcast(irow, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(centsh, blocksLAST, bytesize)
      call mpi_bcast(centsh, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(skip, ccsdgninpLAST, bytesize)
      call mpi_bcast(skip, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(ccs, ccmodelsLAST, bytesize)
      call mpi_bcast(ccs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(etmp, etmpLAST, bytesize)
      call mpi_bcast(etmp, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nckijmax, ccsdmaxLAST, bytesize)
      call mpi_bcast(nckijmax, bytesize, mpi_byte,0,mpi_comm_world,ierr)

      call getbytespan(nt1amx, ccsdsymLAST, bytesize)
      call mpi_bcast(nt1amx, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(it2del, ccsdioLAST, bytesize)
      call mpi_bcast(it2del, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mxcall, distclLAST, bytesize)
      call mpi_bcast(mxcall, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(runeri, cbieriLAST, bytesize)
      call mpi_bcast(runeri, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(sotime, soppinfLAST, bytesize)
      call mpi_bcast(sotime, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nexci2, soppexcLAST, bytesize)
      call mpi_bcast(nexci2, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(soorwc, rwinfLAST, bytesize)
      call mpi_bcast(soorwc, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(expbt, aobtchLAST, bytesize)
      call mpi_bcast(expbt, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nitcl, odclssLAST, bytesize)
      call mpi_bcast(nitcl, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(thrs, ccomLAST, bytesize)
      call mpi_bcast(thrs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gtotyp, ccomcLAST, bytesize)
      call mpi_bcast(gtotyp, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(scrmab, ericomLAST, bytesize)
      call mpi_bcast(scrmab, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrsh, erithrLAST, bytesize)
      call mpi_bcast(thrsh, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(memadd, erimemLAST, bytesize)
      call mpi_bcast(memadd, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nodbch, odbtchLAST, bytesize)
      call mpi_bcast(nodbch, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(charge, nucleiLAST, bytesize)
      call mpi_bcast(charge, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(ndistr, eridstLAST, bytesize)
      call mpi_bcast(ndistr, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(namn, nuclecLAST, bytesize)
      call mpi_bcast(namn, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(maxrep, symmtiLAST, bytesize)
      call mpi_bcast(maxrep, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(fmult, symmtrLAST, bytesize)
      call mpi_bcast(fmult, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gamac, comr12LAST, bytesize)
      call mpi_bcast(gamac, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(mbas1, cmmmulLAST, bytesize)
      call mpi_bcast(mbas1, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(jtop, hertopLAST, bytesize)
      call mpi_bcast(jtop, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(zcmval, cbireaLAST, bytesize)
      call mpi_bcast(zcmval, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mulnam, cbirea_cLAST, bytesize)
      call mpi_bcast(mulnam, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nmulbs, cmmbasLAST, bytesize)
      call mpi_bcast(nmulbs, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(i2bst, symsqLAST, bytesize)
      call mpi_bcast(i2bst, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrpckint, ccpackLAST, bytesize)
      call mpi_bcast(thrpckint,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(luiajb, cc_tapLAST, bytesize)
      call mpi_bcast(luiajb,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(gradml, gnrinfLAST, bytesize) 
      call mpi_bcast(gradml, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call stupid_isao_bcast_routine()

      ! Common blocks with explicitly set sizes, no need to calculate bytespan
      ! Though we should probably still insert the actual parameter governing
      ! their length, eg. lbasdir+12 = 612
      call mpi_bcast(basdir,612, mpi_character, 0, mpi_comm_world, ierr)
      call mpi_bcast(fnvajkl, 10, mpi_character, 0, mpi_comm_world,ierr)
      call mpi_bcast(vclthr, 24, mpi_byte, 0, mpi_comm_world, ierr)

      ! The slaves need to create the iadrpk array by 
      ! calling the module functions rather than with a bcast
      ! iadrpk_dim is initialized in get_iadrpk, no need to send it 
!      call mpi_bcast(iadrpk_dim, 1, mpi_integer, 0,mpi_comm_world,ierr) 
      if (.not. allocated(iadrpk) ) then
         call get_iadrpk(lupri,nsym,muld2h,nbas,
     &                   nbast,i2bst,iaodis,iaodpk)
      endif


      CALL DZERO(WORK(KAIJ),LAIJ)
      CALL DZERO(WORK(KAAB),LAAB)

      forceupdate = .false. ! This avoids redundant updating of common blocks. SO_ERES sets this to true when necessary and this routine sets it to false when necessary.
      endif commonblocktransfers

C**************************************************
C     This data needs to be updated every time PAR_SO_ERES is called 
C**************************************************

C
C This broadcasts the trial-vectors, ao-to-mo matrix, as well as 
C the result vectors and temporary stuff. The latter need probably only
C be explicitly zeroed here.  
      call mpi_bcast(work(KCMO), kendsv, mpi_real8, 0, mpi_comm_world,
     &                ierr)

      call mpi_bcast( work(kindxb), (8*MXSHEL*MXCON+ 1)/IRAT, 
     & mpi_real8, 0, mpi_comm_world, ierr)

      call mpi_bcast(work(kdens),ldens,mpi_real8,0,mpi_comm_world, ierr)


C      if (resetvectors) then
      if(myid.ne.0) then
         call dzero( work(kfock),   lfock)
         call dzero( work(kres1e),  lres1e)
         CALL dzero( work(kres1d),  lres1d)
         call dzero( work(kres2e),  lres2e)
         call dzero( work(kres2d),  lres2d)
         call dzero( work(ksigai1), lsigai1)
         call dzero( work(ksigai2), lsigai2)
         call dzero( work(ksigda1), lsigda1)
         call dzero( work(ksigda2), lsigda2)
      endif

      if  ( (model.ne.'AOSOC') .and. 
     &      (nit.eq.1)         .and.
     &      (inewtr.eq.1)    ) then
         
         allocate(densai(ldensai))
         call dzero(densai, ldensai)
      endif


C**************************************************    
C   Start load balancing
C**************************************************

C   In case there are too many processes compared
C   to the number of ILLL indices, the surplus processes are
C   redirected. They will return from the polling barrier when
C   the master sends a non-blocking logical equal to .true.
      if (myid.ge.ntot) color = 1
      if (myid.lt.ntot) color = 0 !these are the ones performing timings, other processes stall in the polling barrier.
      call mpi_comm_split(mpi_comm_world, color, myid, soppa_comm_active
     &                    ,ierr)
      if (myid.ge.ntot) then
         call pollingbarrier(1)
         goto  800 !skip to assembly of matrices. 
      endif


C     Start the dynamic load balancing
      if(forcedynamicloadbalance) then
          forcedynamicloadbalance = .false. ! This is the last time work is redistributed in an AOSOPPA computation.
          slavesupdatetasks       = .true.  ! Tels the slaves they have to receive work from the master.
          if (myid.eq.0) then
             call dynloadbal_parsoppa(timings, allocsize, temparray)
             AssignedIndices = temparray(1:allocsize)
          endif
      endif

C     Start the presorting load balancing.
      if (forceloadbalance) then
         forceloadbalance        = .false. ! Don't use the pre-sorting load balancing more than once, next time use dynamic instead.
         forcedynamicloadbalance = .true.  ! Monitor the tasks and re-balance once the actual walltimes are known.
         slavesupdatetasks       = .true.

         !performs balancing and sends tasks to slaves.
         if (myid.eq.0) then 
            call presortloadbal_parsoppa(work,lwork,allocsize,temparray)
            maxnumjobs = allocsize
            if (.not. allocated(AssignedIndices) ) then
               allocate(AssignedIndices(maxnumjobs), stat=allocstatus)
            endif
            AssignedIndices = temparray(1:allocsize)
         endif 

         if (.not. allocated(timings) ) then
            allocate(timings(ntot) )
            timings = 0.0
         endif

      endif !forceloadbalance.


C************************************************** 
C   Slaves receive the indices they have to work with
C**************************************************
         if (slavesupdatetasks) then
         slavesupdatetasks = .false.
         IF (myid.ne.0) then

C           Get the number of tasks and allocate memory for the indices
C           Maybe these should be made blocking (since we basically
C           just pol untill the finish, with nothing i between.
C            call mpi_scatter(maxnumjobs, 1, mpi_integer, 0, myid,
C     &                  mpi_comm_world, getnumjobs, ierr)
C            call mpi_irecv(maxnumjobs, 1, mpi_integer, 0, myid,
C     &                  mpi_comm_world, getnumjobs, ierr)
C001         continue
C            call mpi_test(getnumjobs, assignnow, mpistatus, ierr)
C            if (assignnow)  then
               maxnumjobs = ntot - min(nodtot, ntot) + 1  
               if (.not. allocated(AssignedIndices) ) then
                  allocate(AssignedIndices(maxnumjobs),stat=allocstatus)
                  if(.not.(allocstatus.eq.0) ) then
                     call quit('Allocation error; AssignedIndices ')
                  endif
               endif
C            else 
C               goto 001
C            endif



C           Get the pre-sorted ILLL indices.
C           Recieve the scatter, send arguments should be ignored
            call mpi_scatter( AssignedIndices, maxnumjobs, mpi_integer,
     &                        AssignedIndices, maxnumjobs, mpi_integer,
     &                        0, soppa_comm_active, ierr )

C            call mpi_irecv(AssignedIndices, maxnumjobs, mpi_integer, 0, 
C     &                  myid, mpi_comm_world, getindices, ierr)
C002         continue 
C            call mpi_test(getindices, workreceived, mpistatus, ierr)
C            if (workreceived)  then
C               continue 
C            else 
C               goto 002
C            endif

         ENDIF 
         endif

C *******************************************
C           Embarrasingly parallel do-loop 
C *******************************************
      
      doilll: DO 300 workindex=1, maxnumjobs
         ILLL = AssignedIndices(workindex)
         if (ILLL.eq.0) exit


         if (forcedynamicloadbalance) then
            call system_clock(counti, count_rate)
         endif


         if (.not.direct) then
            call quit('AOSOPPA in parallel must use .DIRECT. ',
     &                'integral transformations')
         endif

         KEND2  = KENDSV
         LWORK2 = LWORKSV

         IF (HERDIR) THEN 
            CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
         ELSE

            CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                  WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  WORK(KCCFB1),WORK(KINDXB), 
     &                  WORK(KEND2),LWORK2,IPRINT)
         ENDIF
         LRECNR = ((NBUFX(0) - 1)/IRAT) + 1 
         KRECNR  = KEND2
         KEND2   = KRECNR + LRECNR
         LWORK2  = LWORK  - KEND2

         CALL SO_MEMMAX ('SO_ERES.2B',LWORK2)
         IF (LWORK2 .LT. 0) CALL STOPIT('SO_ERES.2B',' ',KEND2,LWORK)


C------------------------------------------------------------------------
C              Loop over number of distributions in disk.
C------------------------------------------------------------------------
               idel2loop: DO 310 IDEL2 = 1,NUMDIS

                  
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL) ! Keeps track of symmetry
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)

C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
                  LXINT  = NDISAO(ISYDIS)
                  KXINT   = KEND2
                  KEND3   = KXINT  + LXINT
                  LWORK3  = LWORK  - KEND3

                  CALL SO_MEMMAX ('SO_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &               CALL STOPIT('SO_ERES.3',' ',KEND3,LWORK)

C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
C
C---------------------------------------------
C                 Work space allocation no. 4.
C---------------------------------------------
C
                  ISAIJ = MULD2H(ISYMD,1)
C
                  IF (NVIR(ISYMD) .GT. 0) THEN
                     LT2M1 = NT2BCD(ISAIJ)
                  ELSE
                     LT2M1 = 0
                  END IF
C
                  KT2M1   = KEND3
                  KEND4   = KT2M1  + LT2M1
                  LWORK4  = LWORK  - KEND4
C
                  CALL SO_MEMMAX ('SO_ERES.4',LWORK4)
                  IF (LWORK4 .LT. 0) 
     &                CALL STOPIT('SO_ERES.4',' ',KEND4,LWORK)
C
C------------------------------------------------------------
C                 Construct the partially back-transformed T2
C                 MP-amplitudes.
C------------------------------------------------------------
C
C
                  CALL SO_T2M1(WORK(KT2M1),LT2M1,T2MP,LT2MP,
     &                         WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                         WORK(KEND4),LWORK4)
C
C---------------------------------------------
C                 Work space allocation no. 5.
C---------------------------------------------
C
                  LDSRHF = NDSRHF(ISYMD)
C
                  KDSRHF  = KEND4
                  KEND5   = KDSRHF + LDSRHF
                  LWORK5  = LWORK  - KEND5
C
                  CALL SO_MEMMAX ('SO_ERES.5',LWORK5)
                  IF (LWORK5 .LT. 0) 
     &                CALL STOPIT('SO_ERES.5',' ',KEND5,LWORK)
C
C----------------------------------------------------------------
C                 Transform one index in the integral batch to an
C                 occupied index.
C----------------------------------------------------------------
C
                  ISYMLP = 1
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                        ISYMLP,WORK(KEND5),LWORK5,ISYDIS)
C
C-------------------------------------------------------------------
C                 Calculate part of the second order density matrix.
C-------------------------------------------------------------------
C
                  IF ( (MODEL.NE.'AOSOC') .AND. (NIT.EQ.1) .AND. 
     &               (  INEWTR.EQ.1) ) THEN
                     CALL SO_DENSAI1(DENSAI,LDENSAI,WORK(KDSRHF),LDSRHF,
     &                               WORK(KCMO),LCMO,WORK(KT2M1),LT2M1,
     &                               ISYMD,ISYDIS,WORK(KEND5),
     &                               LWORK5)
                  END IF
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                           WORK(KEND5),LWORK5,IDEL,ISYMD,.FALSE.,
     &                           'crashifTrue',ISYMTR)
C----------------------------------------------------------------------
C                 Calculate part of the result vectors RES1E and RES1D,
C                 specifically the first and the second term in eqs.
C                 (34,35). Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C
                  CALL SO_RES_A(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                          WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                          WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                          WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
C
C-------------------------------------------------------------------
C                 Calculate the part of the result vectors RES1E and
C                 RES1D which originate from the C matrices. See 
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
C
                  CALL SO_RES_TCB(WORK(KRES1E),LRES1E,WORK(KRES1D),
     &                            LRES1D,WORK(KTR2E),LTR2E,
     &                            WORK(KTR2D),LTR2D,WORK(KDSRHF),LDSRHF,
     &                            WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                            ISYMTR,WORK(KEND5),LWORK5)
C
C----------------------------------------------------------------------
C                 Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
                  CALL SO_RES_CB(WORK(KRES2E),LRES2E,WORK(KRES2D),
     &                          LRES2D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KBTR1E),LBTR1E,
     &                          WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                          WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                          IDEL,ISYMD,ISYDIS,ISYMTR,WORK(KEND5),
     &                          LWORK5)
C
C--------------------------------------------------------------------
C                 Construct SIGMAI1(ALFA,I) and SIGMAI2(ALFA,I) which
C                 are used in SO_RES_B.
C--------------------------------------------------------------------
C
                  CALL SO_SIGAI(WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),
     &                          LSIGAI2,WORK(KT2M1),LT2M1,WORK(KXINT),
     &                          LXINT,WORK(KBTR1E),LBTR1E,WORK(KBTR1D),
     &                          LBTR1D,WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),
     &                          LBTJ1D,WORK(KCMO),LCMO,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
C
C--------------------------------------------------------------------
C                 Construct SIGDA1(DELTA,A) and SIGDA2(DELTA,A) which
C                 are used SO_RES_C.
C--------------------------------------------------------------------
C
                  CALL SO_SIGDA(WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),
     &                          LSIGDA2,T2MP,LT2MP,WORK(KDSRHF),LDSRHF,
     &                          WORK(KBTR1E),LBTR1E,WORK(KBTR1D),LBTR1D,
     &                          WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),LBTJ1D,
     &                          WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)

310      enddo idel2loop


         if (forcedynamicloadbalance) then
            call system_clock(countf)
            timespent =   ( countf-counti )/REAL(count_rate)
            timings(ILLL) = timespent
         endif

300   enddo doilll

      if (forcedynamicloadbalance ) then
              !assemble the array with timings 
              if (myid.eq.0) then
                 call mpi_reduce(mpi_in_place, timings, ntot, MPI_REAL,
     &                           MPI_SUM, 0, soppa_comm_active, ierr)
              else
                 call mpi_reduce(timings, temparray, ntot, MPI_REAL,
     &                           MPI_SUM, 0, soppa_comm_active, ierr)
              endif
      endif

C**************************************************  
C  END OF PARALLEL PORTION 
C**************************************************  


C     Relase all surplus processes from the spinning barrier.
C     The slaves are redirected to label 800
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      if (numprocs.gt.ntot) then
         if (myid.eq.0) then
            do i=ntot, (numprocs-1)
               call mpi_isend(.true., 1, MPI_LOGICAL, i,
     &         i, mpi_comm_world, request, ierr)
            enddo
         endif
      endif
800   continue 


C ************************************************** 
C Now assemble the Fock matrix, the result vectors and the sigma vectors
C in the master's work array.
C Assemble the densai matrix as well and make it ready to be
C transferred back to the master.
C ************************************************** 

      lres = lres1e +  lres1d +  lres2e +  lres2d 
      lsig = lsigai1 + lsigai2 + lsigda1 + lsigda2
      allocate( temp( max(lfock,lres,lsig, laab, laij) ) )

C     Assemble the Fock Matrix in the master's work array
      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, work(kfock), lfock, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce( work(kfock), temp, lfock, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif


C     Assemble the Result vector in the master's work array
      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, work(kres1e), lres, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce( work(kres1e), temp, lres, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif
        
 
C     Assemble the sigda vector in the master's work array
      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, work(ksigai1), lsig, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce( work(ksigai1), temp, lsig, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif


C     Assemble KAIJ and KAAB (from routine so_res_a) in the master's work array
      if (aijreduce) then
      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, work(kaab), laab, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call mpi_reduce(MPI_IN_PLACE, work(kaij), laij, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce( work(kaab), temp, laab, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call mpi_reduce( work(kaij), temp, laij, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif
            aijreduce = .false. 
      endif



C**************************************************C
C     Make use of auxiliary slave to
C     send info back to master.
C**************************************************C

C     Assemble the density matrix on process 1.
C     Process 1 will send this matrix to the master, once the master
C     has successfully exited PAR_SO_ERES and returned to SO_ERES.
         IF ( (MODEL.NE.'AOSOC') .AND. 
     &        (NIT.EQ.1)         .AND. 
     &        (INEWTR.EQ.1)  )    THEN

            allocate(tempdens(ldensai))
            tempdens = 0.0000

            call mpi_reduce(densai, tempdens, ldensai, MPI_REAL8,
     &      MPI_SUM, 1, MPI_COMM_WORLD, ierr)

            ! This array is given to SO_ERES as an input parameter. It needs
            ! to be transferred back from PAR_SO_ERES before process 1 returns.
            ! The master posts a recv call immediately after returning from PAR_SO ERES
            getdensai = .true. !signals the master that densai will be transferred.
            if (myid.eq.1) then
               call mpi_send(tempdens, ldensai, mpi_real8, 0, 0, 
     &         mpi_comm_world, ierr)
            endif

            deallocate(tempdens, stat=deallocstatus)
            deallocate(densai, stat=deallocstatus)

         endif


C************************************************** 
C      Clean up allocated arrays 
C************************************************** 
      deallocate(temp, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Deallocation error in PAR_SO_ERES (temp)')
      endif

      deallocate(t2mp, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
        call quit('t2mp was not deallocated in PAR_SO_ERES')
      endif


      return
      end subroutine
C
C
C
      subroutine dynloadbal_parsoppa(timings, allocsize, temparray)
C    Dynamic Load Balancing for the parallel SOPPA calculations (and parallel RPA). 
C    The routine assumed that the MPI processes that does the load balancing is the master. If any other routine enters, there will be issues with the updated ILLL indices in the AssignedIndices array.  Right now, any slave that enters is immediately evicted from the routine.
C
C The routine takes an array of timings as input. The timings are creatd on the fly by every parallel process and the index corresponds to the ILLL index in the parallel calculation. The time is given in integers from calling system_clock. The timings are not given in any human-readable form. This choice was made so that one can reuse the subroutine getallocsize, which relies on an integer input array. 
C The work is resorted once the actual time associated with an ILLL index is known and the indices are rebalanced once more and sent to all slaves.
C The improvement over the presorted balancing scheme is expected to be very small, but gets better the larger the basis set.
C 
C F.Beyer Oct. 2014.

#include"implicit.h"
#include"mpif.h"
#include"parsoppa.h" !fetches ntot

      ! Dummy parameters
      real, dimension(ntot), intent(inout) :: timings
      integer, intent(out) :: allocsize
      integer, dimension(ntot), intent(out) :: temparray

      ! Bookkeeping
      integer, dimension(:,:), allocatable :: sortedindices
      integer, dimension(:), allocatable   :: numassignjobs
      real, dimension(:), allocatable   :: sumofwork
      integer :: maxrows, maxcols, numprocs, numrecipients
      integer :: getnumjobs, targetID, myid
      integer :: colindex, rowindex, assignILLL, col

      ! Evict slaves from the routine
      call mpi_comm_rank(mpi_comm_world, myid, ierr)
      if (myid.ne.0) return

      ! In case the amount of work is smaller than
      ! the number of MPI processes...
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      maxcols = min(numprocs, ntot)

      ! Explicitly set the maximum number of jobs a single MPI process can be allocated.
      maxrows = ntot-maxcols+1

      allocate(sortedindices(maxrows, maxcols) )
      allocate(numassignjobs( maxcols) )
      allocate(sumofwork( maxcols) )
      sortedindices(:,:) = 0
      numassignjobs(:) = 0
      sumofwork(:) = 0


      DO i=1, ntot
              ! Find the largest chunk of available work
              assignILLL = maxloc(timings,dim=1, mask=timings.gt.0)

              ! Find the laziest slave
              colindex = minloc(sumofwork, dim=1)

              ! Write the ILLL index to the correct row in the sortedmatrix
              numassignjobs(colindex) = numassignjobs(colindex) + 1
              rowindex = numassignjobs(colindex)
              sortedindices(rowindex, colindex) = assignILLL
              
              ! Update the slave's expected work/walltime and
              ! remove the work from the timings array
              sumofwork(colindex) = sumofwork(colindex) 
     &                             +timings(assignILLL)
              timings(assignILLL) = 0.0

       ENDDO




      ! Tell the slaves how much space to allocate.
      ! Let the slave calculate this themself 
      ! (if the slaves doesn't even agree with master on this, we're in
      ! deep shit)
      allocsize = maxrows
C      numrecipients = min(numprocs, ntot)
C      do targetID=1, (numrecipients-1)
C         call mpi_isend(allocsize, 1, mpi_integer, targetID, 
C     &                  targetID, mpi_comm_world, getnumjobs, ierr) 
C      enddo
         
       
      ! Send the info to every slave that does computation (some slaves might be stalled in the polling barrier in case there are too many MPI processes compared to the number of tasks).
C      ! Use scatter instead?
C      numcols = min(numprocs, ntot)
C      DO col=2, numcols
C         TargetID = col-1 
C         ! Note this isend makes use of column major storage in FORTAN.
C         call mpi_isend(sortedindices(1,col),allocsize,mpi_integer, 
C     &   TargetID, TargetID, mpi_comm_world, getindices,  ierr) 
C      ENDDO

      call mpi_scatter( sortedindices, allocsize, mpi_integer, 
     &                  temparray, allocsize, mpi_integer,
     &                  0, soppa_comm_active, ierr) 

      ! Update the master's own array of ILLL indices
      temparray(1:allocsize) = sortedindices(:,1)

      ! Are we in anyway guaranteed that the above sends finish
      ! before deallocation, or do we need to store all getindices
      ! call mpi_waitall( numcols-1, getindices, mpi_statuses_ignore)
      ! 
      deallocate(sortedindices)
      deallocate(numassignjobs)
      deallocate(sumofwork)

      return
      end subroutine






      subroutine presortloadbal_parsoppa(work, lwork, 
     &           allocsize, temparray)
C     This subroutine load balances a parallel SOPPA/RPA 
C     calculation.
C
C     The routine makes a best guess of loadbalancing by giving 
C     every MPI process an equal number of distributions to handle.
C     Testing shows this is a better first guess that giving every
C     MPI process and equal number of ILLL indices to work with 
C     (since there is a different number of distributions 
C     associated with every ILLL index).
C
C     F.Beyer Oct. 2014.


#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
#include "parsoppa.h" ! to get the value for NTOT and KINDXB

C     Dummy parameters
      integer, intent(in) :: lwork
      integer, intent(out) :: allocsize
      dimension work(lwork) !intent in
      integer, dimension(ntot), intent(out) :: temparray

C     Pre-sorting load balancing variables
      integer :: allocstatus, maxnumjobs, TargetID, 
     &           getnumjobs, getindices, myid, numprocs, ierr
      integer, allocatable, dimension(:,:) :: presortarray, finalsorted
      integer :: col
    

      call mpi_comm_size(mpi_comm_world, numprocs, ierr) 
      call mpi_comm_rank(mpi_comm_world, myid, ierr) 
      call izero(temparray, ntot)

C     FIND THE AMOUNT OF WORK ASSOCIATED WITH EVERY AO INDEX
      allocate( presortarray(2, ntot), stat=allocstatus )
      if(.not.(allocstatus.eq.0) ) then
         call quit('Allocation error in presortarray')
      endif
      call presortaodist(ntot, work(kindxb), presortarray)


      maxnumjobs = ntot-min(numprocs,ntot)+1 
      !call getallocsize(ntot, presortarray, maxnumjobs) 

C     The master tells the slaves how much space to allocate ...
C      numrecipients = min(numprocs, ntot)
C      do targetID=1, (numrecipients-1)
C         call mpi_isend(maxnumjobs, 1, mpi_integer, targetID, 
C     &                  targetID, mpi_comm_world, getnumjobs, ierr) 
C      enddo
         

C     CREATE THE FINAL MATRIX OF PRE-SORTED AO INDICES
      allocate(finalsorted(maxnumjobs, numprocs)
     &               ,stat=allocstatus )
      if(.not.(allocstatus.eq.0) ) then
         call quit('The matrix finalsorted was not properly ', 
     &             'allocated')
      endif
      call partitionAOindices(ntot, maxnumjobs, numprocs, 
     &                        presortarray, finalsorted)


C     TRANSFER AO INDICES TO SLAVES
C      numcols = min(numprocs, ntot)
C      DO col=2, numcols
C         TargetID = col-1 
C         call mpi_isend(finalsorted(1,col),maxnumjobs,mpi_integer, 
C     &   TargetID, TargetID, mpi_comm_world, getindices,  ierr) 
C      ENDDO
C     Use scatter instead
      call mpi_scatter( finalsorted, maxnumjobs, mpi_integer,
     &                  temparray, maxnumjobs, mpi_integer,
     &                  0, soppa_comm_active, ierr )

C     Keep own work in a vector and don't leak memory...
      allocsize = maxnumjobs
!     temparray(1:maxnumjobs) = finalsorted(:,1)  
      deallocate(finalsorted)
      deallocate(presortarray)



      return
      end subroutine 






      
      subroutine pollingbarrier(pollinginterval)
C     This subroutine is implemented in case there is ever a need to 
C     remove certain processes from a calculation.
C
C     The processes will repeatedly poll until they receive a non-blocking 
C     send from the master with a logical value equal to .true.
C
C     F.Beyer Oct. 2014
      
#include "implicit.h"
#include "mpif.h"

      integer, dimension(MPI_STATUS_SIZE) :: mpistatus
      integer, intent(in)                 :: pollinginterval !in milliseconds
      integer                             :: myid, request
      logical                             :: flag, exitbarrier
      volatile                            :: exitbarrier

      call mpi_comm_rank(mpi_comm_world, myid, ierr)
      print *, 'I am in the polling barrier ', myid
           
      ! Initiate the polling variables and ask for the non-blocking update
      exitbarrier = .false. 
      call mpi_irecv(exitbarrier, 1, MPI_LOGICAL, 0, myid,
     &               mpi_comm_world, request, ierr) 

      ! Polling barrier, the process will cycle repeatedly until released.
130   continue 
#ifdef VAR_IFORT
      call sleepqq(pollinginterval)  
#endif
      call mpi_test(request, flag, mpistatus, ierr)
      if (.not. flag) goto 130
      ! Warning, seems that we'll go into an infinite loop,
      ! if exitbarrier is ever sent as .false.
      if (.not.exitbarrier) then
         goto 130 ! Cycle to the top of the barrier.
      else
         return
      endif

      return
      end subroutine

      subroutine stupid_isao_bcast_routine()
#include "implicit.h"
#include "mpif.h"
#include "maxorb.h"
#include "ccisao.h"
#include "infpar.h"
C
C  This is a simple (but stupid) work-around to the problem that
C  the array ISAO exist in both infind.h and ccisao.h and that it seem
C  pretty random which version is needed by which routine
C  (the content appear to be the same).
C  A better solution would maybe be to remove isao from common block
C  in infind.h and explicitly include ccisao.h there...
C  Or even just replace ISAO with the function
C
C  pure function isao ( aoindex ) result(isym)
C  #include ! something that sets nbas(8), nsym
C  integer, intent(in) :: aoindex
C  integer :: i, isym, popsum
C  sum = 0
C  do isym = 1, nsym
C     popsum = nbas(isym) + popsum
C     if ( aoindex .le. popsum ) return 
C  enddo
C  ! Some error statement, here would be in order
C  end function isao
C
      integer  :: bytesize, ierr

C    This size SHOULD be MXCORB x "integer size"
C    But since this is a cc block, someone may change that to MXCORB_CC       
      call getbytespan(isao, ccisaolast, bytesize )
      call mpi_bcast(isao, bytesize, mpi_byte, master, mpi_comm_world,
     &               ierr)

      end subroutine stupid_isao_bcast_routine


#endif !VAR_MPI at the beginning of PAR_SO_ERES

#ifdef NOTDEFINED
      subroutine stupid_checksum( nameofdata , checkdata, ldata )
C   For debugging purposes
C  Calculate a few "checksums" of the data in order to ensure 
C  that it is (relatively) intact
C
C  This is not a cryptografic checksum, just calculat the sum
C  in a few different ways (also since we want to keep the order
C  of the significant bits) 
         implicit none 
#include "priunit.h"
         integer, parameter  :: dp = kind(1.0D0)
         character(len=*),intent(in):: nameofdata
         integer, intent(in) :: ldata
         real(dp), intent(in):: checkdata(ldata)
         real(dp)            :: signed_sum, abs_sum, squared_sum, temp
         integer             :: i
         
         signed_sum = 0.0_dp
         abs_sum = 0.0_dp
         squared_sum = 0.0_dp

         do i = 1, ldata 
            temp = checkdata(i)
            signed_sum = temp + signed_sum
            abs_sum = abs(temp) + abs_sum
            squared_sum = temp*temp + squared_sum
         enddo

         write (lupri,90) nameofdata, squared_sum, abs_sum, signed_sum
90       format(A6,3F23.15)

      end subroutine stupid_checksum
#endif
