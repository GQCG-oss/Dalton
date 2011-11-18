!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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
#include <single.h>
C
* ==================================================================== *
      SUBROUTINE WFNAO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           NOINT,KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI,
     &           WORK,LWORK)
* -------------------------------------------------------------------- *
C
C     generate a *.wfn interface file for atoms-in-molecules programs
C     (PROAIM, PROMEGA, MURPHY, etc.) which contains the basis set
C     information in form of a list of exponents, functions types and
C     transformation matrix from the primitive cartesian GTO to the
C     symmetry adapted AO basis
C
* -------------------------------------------------------------------- *

#include <implicit.h>
#include <molde.h>

#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <codata.h>
      PARAMETER (D0 = 0.0D0)
C
#include <ccom.h>
#include <cbirea.h>
#include <nuclei.h>
#include <primit.h>
#include <maxorb.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <sphtrm.h>
#include <frame.h>


      CHARACTER*10 CHRSEG
      LOGICAL SEG, NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          SEG(KBLOCK),NOINT(KBLOCK)
      DIMENSION WORK(LWORK)

C     ----------------------------------------------------------------
C     count:  -  number of cartesian primitive GTO's
C             -  number of (spherical or cartesian) AO's
C             -  number of centers
C     ----------------------------------------------------------------
      NAO     = 0
      NPRIM   = 0
      NCENTER = 0

      ICENT  = 0
      IBLOCK = 0
      DO I = 1, NONTYP   ! atom types
      DO N = 1, NONT(I)  ! symmetry distinct atoms of type I
         ICENT   = ICENT + 1
         NCENTER = NCENTER + NUCDEG(ICENT)
         KBCH    = IBLOCK
         DO J = 1, IQM(I)    ! angular momentum + 1
         DO K = 1, JCO(J,I)  ! contraction blocks for this ang. mom. 
            KBCH  = KBCH + 1
            NPRIM = NPRIM + NUC(KBCH) * NUCDEG(ICENT) * KCK(J)
            NAO   = NAO   + NRC(KBCH) * NUCDEG(ICENT) * KHK(J)
            print *, I,ICENT,KBCH,NRC(KBCH),NUCDEG(ICENT),KHK(J),NAO
         END DO
         END DO
      END DO
      IBLOCK = IBLOCK + NBLCK(I)
      END DO

      IF (NAO .NE. NBASIS) THEN
         PRINT *, "Inconsistency in the number of AO's in WFNAO."
         STOP     "Inconsistency in the number of AO's in WFNAO."
      END IF

      IF (NCENTER .NE. NUCDEP) THEN
         PRINT *, "Inconsistency in the number of centers in WFNAO."
         STOP     "Inconsistency in the number of centers in WFNAO."
      END IF

C     ----------------------------------------------------------------
C     memory administration:
C     ----------------------------------------------------------------
      KCNTASS = 1
      KPTYPE  = KCNTASS + NPRIM
      KEXPON  = KPTYPE  + NPRIM
      KOCCNUM = KEXPON  + NPRIM
      KORBENE = KOCCNUM + NBASIS
      KPOSAO  = KORBENE + NBASIS
      KSCR    = KPOSAO  + NBASIS * NPRIM
      KEND    = KSCR    + NBASIS * NPRIM 

      IF (KEND .GT. LWORK) THEN
         PRINT *,'Insufficient memory in WFNAO.'
         STOP    'Insufficient memory in WFNAO.'
      END IF

C     ----------------------------------------------------
C     set occupation numbers and orbital energies to zero:
C     ----------------------------------------------------
      CALL DZERO(WORK(KOCCNUM),NBASIS)
      CALL DZERO(WORK(KORBENE),NBASIS)

C     ----------------------------------------------------------------
C     compute primitive to SAO transformation matrix:
C     1. step: generate primitive to AO transformation by patching
C     all contraction coefficients into one large matrix and 
C     incorporating thereby the cartesian -> spherical transformation
C     (if not the DOCART flag is set)
C     ----------------------------------------------------------------

      ICENT   = 0
      IBLOCK  = 0
      IPSTRT  = 0
      IOFFPR  = 0
      IOFFAO  = 0
      IOFFCNT = 0

      DO I = 1, NONTYP   ! atom types
        DO N = 1, NONT(I)  ! symmetry distinct atoms of type I
          ICENT = ICENT + 1
          NDEG  = NUCDEG(ICENT)
          KBCH  = IBLOCK
          DO J = 1, IQM(I)    ! angular momentum + 1
          DO K = 1, JCO(J,I)  ! contraction blocks for this ang. mom. 
            KBCH = KBCH + 1
            NNUC = NUC(KBCH)
            NNRC = NRC(KBCH)
            IF (NNUC .NE. 0) THEN

              IOFFWRK = KSCR + IOFFPR * NBASIS 
              LENGTH  = NBASIS * NDEG * KCK(J) * NNUC
              IF (IOFFWRK+LENGTH .GT. LWORK) THEN
                PRINT *, 'Insufficient memory in WFNAO.'
                STOP     'Insufficient memory in WFNAO.'
              END IF
              CALL DZERO(WORK(IOFFWRK+1),NBASIS*NDEG*KCK(J)*NNUC)

              DO LCENT = 1, NDEG     ! `degenerate' centers

                ICNT = IOFFCNT + LCENT 

              DO ISPHC = 1, KHK(J)   ! spherical components
              DO ICARC = 1, KCK(J)   ! cartesian components

                IF (DOCART) THEN
                  COEF = 0.0D0
                  IF (ISPHC.EQ.ICARC) COEF = 1.0D0
                ELSE
                  COEF = CSP(ISPADR(J)+(ICARC-1)*KHK(J)+ISPHC-1)
                END IF

                DO INNUC = 1, NNUC     ! uncontracted function
c
c                 primitives orderd as: prim.bf. / sph.comp. / deg.cent.
c
                  IPR = IOFFPR + ((INNUC-1)*KCK(J)+(ICARC-1))*NDEG+LCENT

                  WORK(KCNTASS-1 + IPR) = DBLE(ICNT)
                  WORK(KEXPON-1  + IPR) = PRIEXP(IPSTRT+INNUC)
                  WORK(KPTYPE-1  + IPR) = DBLE((J-1)*J*(J+1)/6 + ICARC)

                  DO INNRC = 1, NNRC     ! contracted function
c
c                   AO's ordered as: contr.bf. / car.comp. / deg.cent.
c
                    IAO = IOFFAO+((INNRC-1)*KHK(J)+(ISPHC-1))*NDEG+LCENT

                    IWAO = KSCR-1 + (IPR-1) * NBASIS + IAO
                    WORK(IWAO) = COEF * CPRIMU(INNUC,INNRC,KBCH)

                    WRITE (6,'(X,A3,X,2I5,F15.7,X,F15.10,X,2I5,F5.1)')
     &                   GTOTYP(NHKOFF(J)+ISPHC),IAO,IPR,
     &                   PRIEXP(IPSTRT+INNUC),CPRIMU(INNUC,INNRC,KBCH),
     &                   ISPHC,ICARC,COEF

                  END DO
                END DO

              END DO
              END DO
              END DO

              IOFFAO = IOFFAO + NNRC * NDEG * KHK(J)
              IOFFPR = IOFFPR + NNUC * NDEG * KCK(J)
              IPSTRT = IPSTRT + NNUC

            END IF ! (NNUC .NE. 0) 
          END DO
          END DO
          IOFFCNT = IOFFCNT + NDEG
        END DO
        IBLOCK = IBLOCK + NBLCK(I)
      END DO

      ! total number of primitives
      NPRIM = IOFFPR

      CALL AROUND('PRIMITIVE TO AO TRANSFORMATION:')
      CALL OUTPUT(WORK(KSCR),1,NBASIS,1,NPRIM,NBASIS,NPRIM,1,6)


      IF ( (KPOSAO + NBASIS*NPRIM) .GT. LWORK) THEN
        PRINT *, 'Insufficient memory in WFNAO.'
        STOP     'Insufficient memory in WFNAO.'
      END IF

C     ----------------------------------------------------------------
C     compute primitive to SAO transformation matrix:
C     2. step: multiply the primitive->AO transformation with the 
C     sparse matrix describing the symmetrization; at the same time 
C     we transpose the transformation matrix such the primitive index
C     becomes the leading index
C     ----------------------------------------------------------------
      ISAO = 0
      DO IREP = 1, MAXREP+1
        IF (NAOS(IREP).GT.0) THEN
           DO L = 1, NAOS(IREP)
             ISAO = ISAO + 1    
             ICENT = IPCEN(ISAO)
             J = NUCDEG(ICENT)
             
             WRITE (*,'(I5,3X,A6,3X,A4,5X,I3,7(2X,A,1X,I3))')
     &          ISAO,NAMN(ICENT),GTOTYP(IPTYP(ISAO)),ITRAN(ISAO,1),
     &          (CHRSGN(NINT(CTRAN(ISAO,K))),ITRAN(ISAO,K),K=2,J)

             IAO    = ITRAN(ISAO,1)
             IPOSAO = KPOSAO + NPRIM*(ISAO-1)
             CALL DCOPY(NPRIM,WORK(KSCR-1+IAO),NBASIS,WORK(IPOSAO),1)

             DO K = 2, J
               IAO = ITRAN(ISAO,K)
               FAC = CTRAN(ISAO,K)
 
               IPOSAO = KPOSAO + NPRIM*(ISAO-1)
               CALL DAXPY(NPRIM,FAC,WORK(KSCR-1+IAO),NBASIS,
     &                              WORK(IPOSAO),1)
             END DO

           END DO

        END IF
      END DO

      CALL AROUND('PRIMITIVE TO SAO TRANSFORMATION:')
      CALL OUTPUT(WORK(KPOSAO),1,NPRIM,1,NBASIS,NPRIM,NBASIS,1,6)

      CALL WRWFN('WFN.AO','AO',NBASIS,0,NPRIM,NCENTER,
     &           WORK(KCNTASS),WORK(KPTYPE),WORK(KEXPON),
     &           WORK(KOCCNUM),WORK(KORBENE),WORK(KPOSAO),
     &           'NUCLEAR REPULSION=',POTNUC,
     &           ' XXXXXXXXXXXXXXXXX',0.0D0)

      RETURN
      END

* ==================================================================== *
      SUBROUTINE WFNMO(FNWFNMO,MOTYPE,CMO,OCC,
     &                 STRENERGY,ENERGY,STRVIRIAL,VIRIAL,
     &                 WORK,LWORK)
* -------------------------------------------------------------------- *
C
C     generate a *.wfn interface file for atoms-in-molecules programs
C     (PROAIM, PROMEGA, MURPHY, etc.) which contains the basis set
C     information in form of a list of exponents, functions types and
C     transformation matrix from the primitive cartesian GTO to the
C     MO/NO basis
C
C     basis set information is read from an existing WFN.AO file
C
* -------------------------------------------------------------------- *
#include <implicit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <molinp.h>
#include <nuclei.h>
#include <symmet.h>
#include <molde.h>
#include <inforb.h>
     
      INTEGER LUWFN
      PARAMETER( LUWFN = 65 )
      CHARACTER*6 FNWFNAO
      PARAMETER( FNWFNAO = 'WFN.AO' )

      CHARACTER FNWFNMO*(6), MOTYPE*(2), STRENERGY*(18), STRVIRIAL*(18)
      CHARACTER LINE*(80)
      DIMENSION CMO(*), OCC(*), WORK(LWORK)

      INTEGER IWRK(20)

C     ---------------------------------------------------
C     open "WFN.AO" file and read the header information:
C     ---------------------------------------------------
      OPEN(UNIT=LUWFN,FILE=FNWFNAO,IOSTAT=IOS,ACCESS='SEQUENTIAL',
     &     STATUS='OLD',FORM='FORMATTED',ERR=990)

      READ(LUWFN,'(A80)') LINE
      READ(LUWFN,'(8X,I3,9X,3(I3,17X))') IDUM, NAO, NPRIM, NATOM

      IF (NAO .NE. NBAST) THEN
         PRINT *, "Inconsistency in the number of AO's in WFNMO."
         STOP     "Inconsistency in the number of AO's in WFNMO."
      END IF

C     ----------------------
C     memory administration:
C     ----------------------
      KCNTASS = 1
      KEXPON  = KCNTASS + NPRIM
      KPTYPE  = KEXPON  + NPRIM
      KPOSAO  = KPTYPE  + NPRIM
      KEND    = KPOSAO  + NPRIM * NAO

      IF (KEND .GT. LWORK) THEN
         PRINT *,'Insufficient memory in WFNMO.'
         STOP    'Insufficient memory in WFNMO.'
      END IF


      DO IATOM = 1, NATOM
        READ(LUWFN,*) 
      END DO

      DO ILINE = 1, (NPRIM+19)/20
        MPRIM = MIN(ILINE*20,NPRIM) - (ILINE-1)*20
        READ(LUWFN,'(20X,20I3)') (IWRK(I),I=1,MPRIM)
        DO I = 1, 20
          WORK(KCNTASS-1 + (ILINE-1)*20 + I) = DBLE(IWRK(I))
        END DO
        WRITE(*,'("CENTER ASSIGNMENTS  ",20I3)') (IWRK(I),I=1,MPRIM)
      END DO

      DO ILINE = 1, (NPRIM+19)/20
        MPRIM = MIN(ILINE*20,NPRIM) - (ILINE-1)*20
        READ(LUWFN,'(20X,20I3)') (IWRK(I),I=1,MPRIM)
        DO I = 1, 20
          WORK(KPTYPE-1 + (ILINE-1)*20 + I) = DBLE(IWRK(I))
        END DO
        WRITE(*,'("TYPE ASSIGNMENTS    ",20I3)') (IWRK(I),I=1,MPRIM)
      END DO

      READ(LUWFN,'(10X,5D14.7)') (WORK(KEXPON-1+I),I=1,NPRIM)

      DO IDXMO = 1, NBAST
       READ(LUWFN,*) 
       READ(LUWFN,'(5D16.8)',IOSTAT=IOS,ERR=990)
     &   (WORK(KPOSAO-1+(IDXMO-1)*NPRIM+I),I=1,NPRIM) 
      END DO

      CLOSE(LUWFN,IOSTAT=IOS,ERR=990)  

*----------------------------------------------------------------------*
*     patch symmetry blocked CMO coefficients into one large matrix,
*     and multiply with the PO -> SAO transformation matrix:
*----------------------------------------------------------------------*
      KAOMO = KEND
      KPOMO = KAOMO + NBAST * NBAST
      KEND  = KPOMO + NBAST * NPRIM

      IF (KEND .GT. LWORK) THEN
         PRINT *,'Insufficient memory in WFNMO (2).'
         STOP    'Insufficient memory in WFNMO (2).'
      END IF

      CALL DZERO(WORK(KAOMO),NBAST*NBAST)
   
      KCMO    = 1
      IOFFAO  = 0
      IOFAOMO = 0
      DO ISYM = 1, NSYM
         DO IDXMO = 1, NAOS(ISYM)
            IADR = KAOMO + IOFAOMO + IOFFAO
            CALL DCOPY(NAOS(ISYM),CMO(KCMO),1,WORK(IADR),1)
            IOFAOMO = IOFAOMO + NBAST
            KCMO    = KCMO + NAOS(ISYM)
         END DO
         IOFFAO = IOFFAO + NAOS(ISYM)
      END DO

      PRINT *,NPRIM,NBAST,WORK(KPOSAO),WORK(KAOMO),WORK(KPOMO)
      CALL DGEMM('N','N',NPRIM,NBAST,NBAST,
     &           1.0D0,WORK(KPOSAO),NPRIM,WORK(KAOMO),NBAST,
     &           0.0D0,WORK(KPOMO),NPRIM)

      CALL AROUND('PO TO SAO TRANSFORMATION:')
      CALL OUTPUT(WORK(KPOSAO),1,NPRIM,1,NBAST,NPRIM,NBAST,1,6)

      CALL AROUND('SAO TO MO TRANSFORMATION:')
      CALL OUTPUT(WORK(KAOMO),1,NBAST,1,NBAST,NBAST,NBAST,1,6)

      CALL AROUND('PO TO MO TRANSFORMATION:')
      CALL OUTPUT(WORK(KPOMO),1,NPRIM,1,NBAST,NPRIM,NBAST,1,6)

*----------------------------------------------------------------------*
*     write "WFN.MO" file:
*----------------------------------------------------------------------*
      CALL WRWFN(FNWFNMO,MOTYPE,NBAST,0,NPRIM,NATOM,
     &           WORK(KCNTASS),WORK(KPTYPE),WORK(KEXPON),
     &           OCC,OREN,WORK(KPOMO),
     &           STRENERGY,ENERGY,STRVIRIAL,VIRIAL)

      RETURN
*----------------------------------------------------------------------*
*     Handle I/O Error
*----------------------------------------------------------------------*
990   CONTINUE
      WRITE(*,'(/5x,a,a,/5x,a,i5)')
     .     'An error occured while reading from file ', FNWFNAO,
     .     'IOSTAT = ', IOS
      CLOSE(UNIT=LUWFN)
      STOP 'I/O ERROR IN WFNMO.'
      END 

* ==================================================================== *
      SUBROUTINE WRWFN(FNWFN,MOTYPE,NMO,NVO,NPRIM,NATOM,
     &                 CNTASS,PTYPE,EXPONENT,OCCNUM,ORBENE,POMO,
     &                 STRENERGY,ENERGY,STRVIRIAL,VIRIAL)
* -------------------------------------------------------------------- *
C
C     Purpose: write *.wfn interface file 
C
* -------------------------------------------------------------------- *
#include <implicit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <molinp.h>
#include <nuclei.h>
#include <symmet.h>

      INTEGER LUWFN
      PARAMETER( LUWFN = 65 )

      INTEGER NVO, NMO, NPRIM, NATOM
      CHARACTER MOTYPE*(2), FNWFN*(6), STRENERGY*(18), STRVIRIAL*(18)
      DIMENSION CNTASS(*), EXPONENT(*), PTYPE(*), OCCNUM(*), ORBENE(*)
      DIMENSION POMO(NPRIM,*)

      INTEGER IOS, ICENT, MCENT, LCENT, I, N
#include <ibtfun.h>

C     ----------------------------------------------------
C     open wfn file and write header:
C     ----------------------------------------------------
      OPEN(UNIT=LUWFN,FILE=FNWFN,IOSTAT=IOS,ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED',ERR=990)

      WRITE(LUWFN,'(A80)') MLINE(NMLAU-1)

      WRITE(LUWFN,'(A8,I3,A9,3(I3,A17))') 'DALTON  ',
     &      NVO,   'VIRTUALS ',
     &      NMO,   ' MOL ORBITALS    ',
     &      NPRIM, ' PRIMITIVES      ',
     &      NATOM, ' NUCLEI          '      

C     ----------------------------------------------------
C     write labels, coordinates and charge of the centers:
C     ----------------------------------------------------
      MCENT  = 0 
      DO ICENT = 1, NUCIND  ! symmetry distinct atoms of type I
         
         IF ( MULT(ISTBNU(ICENT)) .EQ. 1 ) THEN
           MCENT = MCENT + 1
           WRITE(LUWFN,'(1X,A4,3X,A11,I3,A2,3F12.8,A10,F5.1)')
     &         NAMN(ICENT),'    (CENTRE',MCENT,') ',
     &         CORD(1,ICENT), CORD(2,ICENT), CORD(3,ICENT),
     &         '  CHARGE =', CHARGE(ICENT)      

         ELSE
           LCENT = 0
           DO ISYMOP = 0, MAXOPR
             IF (IBTAND(ISYMOP,ISTBNU(ICENT)) .EQ. 0) THEN
               LCENT = LCENT + 1
               MCENT = MCENT + 1

               X = PT(IBTAND(ISYMAX(1,1),ISYMOP)) * CORD(1,ICENT)
               Y = PT(IBTAND(ISYMAX(2,1),ISYMOP)) * CORD(2,ICENT)
               Z = PT(IBTAND(ISYMAX(3,1),ISYMOP)) * CORD(3,ICENT)

               WRITE(LUWFN,'(1X,A4,A1,I1,1X,A11,I3,A2,3F12.8,A10,F5.1)')
     &             NAMN(ICENT),'#',LCENT,'    (CENTRE',MCENT,') ',
     &             X,Y,Z, '  CHARGE =', CHARGE(ICENT)      

             END IF
           END DO
         END IF

      END DO

      IF (MCENT.NE.NATOM) THEN
         PRINT *, "Inconsistency in the number of centers in WRWFN."
         STOP     "Inconsistency in the number of centers in WRWFN."
      END IF

C     ---------------------------------------------------------------
C     write center and type assignments and exponents for primitives:
C     ---------------------------------------------------------------
      WRITE(LUWFN,'("CENTER ASSIGNMENTS  ",20I3)',IOSTAT=IOS,ERR=990)
     &      (NINT(CNTASS(I)),I=1,NPRIM)

      WRITE(LUWFN,'("TYPE ASSIGNMENTS    ",20I3)',IOSTAT=IOS,ERR=990)
     &      (NINT(PTYPE(I)),I=1,NPRIM)

      WRITE(LUWFN,'("EXPONENTS ",5D14.7)',IOSTAT=IOS,ERR=990)
     &      (EXPONENT(I),I=1,NPRIM)                       

C     ---------------------------------------------------------------
C     write orbitals to wfn file:
C     ---------------------------------------------------------------
      DO IORB = 1, NMO + NVO
       WRITE(LUWFN,'(A2,I3,A7,A2,A11,A10,F12.7,A15,F12.7)')
     &      MOTYPE, IORB, '       ', MOTYPE, ' 0.0       ',
     &      ' OCC NO = ', OCCNUM(IORB), '  ORB. ENERGY =', ORBENE(IORB)
       WRITE(LUWFN,'(5D16.8)',IOSTAT=IOS,ERR=990)
     &      (POMO(I,IORB),I=1,NPRIM)                       
      END DO

C     ---------------------------------------------------------------
C     write end data record, energy and virial coefficient:
C     ---------------------------------------------------------------
      WRITE(LUWFN,'(A8)',IOSTAT=IOS,ERR=990) 'END DATA'
      WRITE(LUWFN,'(A18,F19.12,A18,1X,F12.8)',IOSTAT=IOS,ERR=990)
     &       STRENERGY,  ENERGY, STRVIRIAL, VIRIAL                  

C     ---------------------------------------------------------------
C     close file and return:
C     ---------------------------------------------------------------
      CLOSE(LUWFN,IOSTAT=IOS,ERR=990)  

      RETURN
*----------------------------------------------------------------------*
*     Handle I/O Error
*----------------------------------------------------------------------*
990   CONTINUE
      WRITE(*,'(/5x,a,a,/5x,a,i5)')
     .     'An error occured while writing to file ', FNWFN,
     .     'IOSTAT = ', IOS
      CLOSE(UNIT=LUWFN)
      STOP 'I/O ERROR IN WRWFN.'
                                   
      END
* ==================================================================== *
