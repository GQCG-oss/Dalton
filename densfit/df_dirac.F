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
      SUBROUTINE DENS_FIT (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IPRFCK,
     &                     WORK,LWORK)
C*****************************************************************************
C    
C     DENS_FIT: Entree routine for density fitting in Dirac & Dalton.
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbieri.h"
#include "odclss.h"
#include "inforb.h"
#include "inftap.h"
      DIMENSION DMAT(*), FMAT(*), WORK(LWORK)
      DIMENSION KOD(12)
      DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)

      LOGICAL FIRST_CALL
      SAVE FIRST_CALL
      DATA FIRST_CALL /.TRUE./
#include <memint.h> 
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C     Transfer information about the number of density matrices
C     to the eri common block
C
      NDMAT = NDMT
C
C     Allocate memory for fit coefficients, work arrays,
C     contraction coefficients and index array
C     TODO : check whether allocation of FTMP is really necessary
C
      CALL MEMGET ('REAL',KCFIT,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KWMAT,NBASISAUX*NBASISAUX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KWVEC,NBASISAUX*NDMT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KFTMP,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KCCFBT,MXPRIM*MXCONT,WORK,KFREE,LFREE)
      CALL MEMGET ('INTE',KINDXB,MXSHEL*MXCONT*8,WORK,KFREE,LFREE)
C
C     Set the CCFBT and INDBX arrays and get the ODC pointers
C
      CALL ERI_POINT (WORK(KFREE),LFREE,WORK(KCCFBT),WORK(KINDXB),
     &                KOD,IPRFCK)
C
C     Update KFREE and LFREE since we want to keep the data generated in
C     ERIPOINT. We will call that part of WORK WKOD from now on.
C
      KWKOD = KFREE
      KFREE = KFREE + KOD(11) - 1
      LFREE = KOD(12)
C
C     Generate the 2-index repulsion integrals over the fit set
C
      IF (FIRST_CALL) THEN
C        The 2-index integrals are written to disk
         WRTINT = .TRUE.
         FCKINT = .FALSE.
C
C        We have symmetry between density 1 and 2
C
         ODTR12 = .TRUE.
C
C        Set the active classes for the ERI calculation.
C        For each electron we take as first index the fit set functions and 
C        as second the null function.
C
C        Index 1
         ICLASSES(1) = 3
         NCENTS(1)   = 0
         ICENTS(1,1) = 0
C        Index 2
         ICLASSES(2) = 0
         NCENTS(2)   = 1
         ICENTS(1,2) = 0
C        Index 3
         ICLASSES(3) = 3
         NCENTS(3)   = 0
         ICENTS(1,3) = 0
C        Index 4
         ICLASSES(4) = 0
         NCENTS(4)   = 1
         ICENTS(1,4) = 0
C
         CALL PICK_ACT_BATCH (NCENTS,ICENTS,ICLASSES,IPRFCK)
C
C        Open the file to which these integrals are to be written
C
         LUINTA = -1
         CALL GPOPEN(LUINTA,'FIT_INTS','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
C        Generate and write the integrals
C
         CALL ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,0,IPRFCK,
     &                    WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                    WORK(KWKOD),WORK(KFREE),LFREE)
C
         FIRST_CALL = .FALSE.
      ENDIF
C
C     From now on we will produce (modified) Fock matrices in ERI
C
      WRTINT = .FALSE.
      FCKINT = .TRUE.
C
C     We will have no symmetry between density 1 and 2
C
      ODTR12 = .FALSE.
C
C     Get the fit coefficients by looping over atom pairs
C
      CALL DRIVE_FIT (FMAT,DMAT,WORK(KWMAT),WORK(KWVEC),NDMT,
     &                ISYMDM,IFCTYP,
     &                WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                WORK(KWKOD),WORK(KFREE),LFREE,WORK(KCFIT),
     &                IPRFCK)
C
C     Build the fock matrix with the fitted density
C     Now the first two indices are regular and the last two relate
C     to the fit set and the null function respectively.
C
C     Index 1
      ICLASSES(1) = 1
      NCENTS(1)   = 0
      ICENTS(1,1) = 0
C     Index 2
      ICLASSES(2) = 1
      NCENTS(2)   = 0
      ICENTS(1,2) = 0
C     Index 3
      ICLASSES(3) = 3
      NCENTS(3)   = 0
      ICENTS(1,3) = 0
C     Index 4
      ICLASSES(4) = 0
      NCENTS(4)   = 1
      ICENTS(1,4) = 0
C
      CALL PICK_ACT_BATCH (NCENTS,ICENTS,ICLASSES,IPRFCK)
C
C     Compute the integrals and build the Fock matrix
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Fit coefficients ',-1)
         CALL OUTPUT(WORK(KCFIT),1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
C
      CALL ERI_MKINTS (FMAT,WORK(KCFIT),NDMT,ISYMDM,IFCTYP,1,IPRFCK,
     &                 WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
     &                 WORK(KWKOD),WORK(KFREE),LFREE)
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Fitted Fock matrix',-1)
         CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C---->The code below restores the original Fock matrix. Should be removed
C     once the implementation is tested and full debugged !
C
*     odtr12 = .true.
*     call pickao(0)
*     call dzero(fmat,ndmt*n2basx)
*     CALL ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,0,IPRFCK,
*    &                 WORK(KCCFBT),WORK(KINDXB),KOD,WORK(KFTMP),
*    &                 WORK(KWKOD),WORK(KFREE),LFREE)
*     IF (IPRFCK.GT.2) THEN
*        CALL HEADER('Exact Fock matrix',-1)
*        CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
*     ENDIF
C----< End of test code.
C
C     Release the memory allocated for the Fock matrix build. This
C     will also check for errors due to memory bound failures
C
      CALL MEMREL('DENS_FIT',WORK,KCFIT,KCFIT,KFREE,LFREE)
      CALL TIMER('DENS_FIT',TIMSTR,TIMEND)
C
      RETURN
      END
C
C  /* Deck eri_point */
      SUBROUTINE ERI_POINT (WORK,LWORK,CCFBT,INDXBT,KOD,IPRFCK)
C*****************************************************************************
C    
C     ERI_POINT: Initialization of ERI pointers for Fock matrix calculation.
C        
C     Taken from ERIFCK by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <veclen.h>
#include <odbtch.h>
#include <symmet.h>
#include <infpar.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
      DIMENSION  WORK(*), CCFBT(*), INDXBT(*)
      DIMENSION KOD(12)
#include <memint.h> 
C
C     Initialization of ERI (if we are the master)
C
      IF (SLAVE) THEN
         IPRINT = IPRFCK
      ELSE
C        
C        Initialization in ER2INI
C
         CALL ER2INI
C
         IPRINT = MAX(IPRERI,IPRFCK)
      END IF
C
C     Give permutational symmetries
C
      PMS12 = .TRUE.
      PMSAB = .FALSE.
      PMSCD = .FALSE.
C
      THRSH  = MAX(THRS,1.00D-15)
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
C     AO batches
C     ==========
C
      CALL SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C     
      CALL ODBCHS(KOD(1),KOD(2),KOD(3),KOD(4),KOD(5),
     &            KOD(6),KOD(7),KOD(8),KOD(9),KOD(10),
     &            KFREE,LFREE,CCFBT,WORK,
     &            LWORK,IPRINT)
C
C     Save the current values of KFREE and LFREE, these are called 
C     KLAST and LWRK, resp. in HR2FCK
C
      KOD(11) = KFREE
      KOD(12) = LFREE
C     
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWORK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C     
      ICALL = 0
      CALL GETDST(ICALL,ICALL,IPRINT)
C
      RETURN
      END
C
C  /* Deck eri_mkints */
      SUBROUTINE ERI_MKINTS (FMAT,DMAT,NDMT,ISYMDM,IFCTYP,IFIT_DMAT,
     &                       IPRFCK,CCFBT,INDXBT,KOD,FCKTMP,WKOD,
     &                       WORK,LWORK)
C*****************************************************************************
C    
C     ERI_MKINTS: Make integrals by calling ERI. 
C     The active classes should be set outside of this routine.
C     Control of writing versus Fock matrix construction is also done
C     on the outside.
C        
C     Second part of ERIFCK, split by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include "implicit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include "nuclei.h"
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <veclen.h>
#include <odbtch.h>
#include <symmet.h>
#include <infpar.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*)
      DIMENSION KOD(12)
C
C     Transfer information about the type of Fock matrix to ERI common
C     Set print level for ERI
C
      IFITDM = IFIT_DMAT
      IPRINT = MAX(IPRERI,IPRFCK)
C     
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(INDXBT,IPRINT)
#if defined (VAR_VECTOR)
      IF (FCKINT) THEN
         ICHUNK = MAX(IVECLN/NDMT,1)
         CALL DZERO(FCKTMP,ICHUNK*NDMT*(NBASE + NODD)*NBASE)
      ENDIF
#endif
C
      KODCL1 = KOD(1)
      KODCL2 = KOD(2)
      KODBC1 = KOD(3)
      KODBC2 = KOD(4)
      KRDBC1 = KOD(5)
      KRDBC2 = KOD(6)
      KODPP1 = KOD(7)
      KODPP2 = KOD(8)
      KRDPP1 = KOD(9)
      KRDPP2 = KOD(10)
C
C     Calculate integrals
C     ===================
C
      IF (SLAVE) THEN
#if defined (VAR_VECTOR)
         CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &               WKOD(KODBC1),WKOD(KODBC2),
     &               WKOD(KRDBC1),WKOD(KRDBC2),
     &               WKOD(KODPP1),WKOD(KODPP2),
     &               WKOD(KRDPP1),WKOD(KRDPP2),
     &               FCKTMP,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK,LWORK,IPRINT)
#else
         CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &               WKOD(KODBC1),WKOD(KODBC2),
     &               WKOD(KRDBC1),WKOD(KRDBC2),
     &               WKOD(KODPP1),WKOD(KODPP2),
     &               WKOD(KRDPP1),WKOD(KRDPP2),
     &               FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &               INDXBT,WORK,LWORK,IPRINT)
#endif
      ELSE
         IF (.NOT.INTSKP) THEN
#if defined (VAR_VECTOR)
            CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &                  WKOD(KODBC1),WKOD(KODBC2),
     &                  WKOD(KRDBC1),WKOD(KRDBC2),
     &                  WKOD(KODPP1),WKOD(KODPP2),
     &                  WKOD(KRDPP1),WKOD(KRDPP2),
     &                  FCKTMP,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK,LWORK,IPRINT)
#else
            CALL ODCDRV(WKOD(KODCL1),WKOD(KODCL2),
     &                  WKOD(KODBC1),WKOD(KODBC2),
     &                  WKOD(KRDBC1),WKOD(KRDBC2),
     &                  WKOD(KODPP1),WKOD(KODPP2),
     &                  WKOD(KRDPP1),WKOD(KRDPP2),
     &                  FMAT,DMAT,NDMT,IFCTYP,DUMMY,IDUMMY,CCFBT,
     &                  INDXBT,WORK,LWORK,IPRINT)
#endif
C
C           Error message in case of insufficient memory
C
            IF (.NOT.MEMOK) THEN
               WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &            ' Not enough memory for this run of ERIFCK.',
     &            ' Available memory in ERIFCK:',LWORK,
     &            ' Required memory for ERIFCK:',LWORK + MEMADD,
     &            ' Increase memory (LWORK) by:',MEMADD
               WRITE (LUPRI,'(/,1X,A,2I5)')
     &            ' Memory requirements largest for OD classes :',
     &              MODAB,MODCD
               CALL QUIT('Insufficient memory in ERIFCK.')
            END IF
         END IF
C
C        Copy and/or print of Fock matrix if we are in that branch
C        This code will only work with regular 2-index Fock matrices,
C        not with the modified 1-index Fock matrix that is made if
C        IFITDM .EQ. 2. The IF statement below should catch this.
C
         IF (FCKINT.AND.ABS(IFITDM).LE.1) THEN
#if defined (VAR_VECTOR)
            IOFF = 0
            DO I = 1, ICHUNK
               DO J = 1, NDMT
                  DO L = 1, NBASE
                     DO K = 1, NBASE
C                       FMAT(K,L,J) = FMAT(K,L,J) + FCKTMP(IOFF + K)
                        FMAT(K+(L-1)*NBASE+(J-1)*NBASE*NBASE) =
     &                       FMAT(K+(L-1)*NBASE+(J-1)*NBASE*NBASE) +
     &                       FCKTMP(IOFF + K)
                     END DO
                     IOFF = IOFF + NBASE+NODD
                  END DO
               END DO
            END DO
#endif
C
C           Print densities and Fock matrix
C           ===============================
C
            IF (IPRINT.GT.4) THEN
               CALL HEADER('Density and Fock matrices in ERIFCK',-1)
               KSTR = 1
               DO I = 1, NDMT
                  WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
                  CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                        NBASE,1,LUPRI)
                  WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
                  CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                        NBASE,1,LUPRI)
                  KSTR = KSTR + NBASE*NBASE
               END DO
            END IF
C
C        End of IF Block for FMAT construction
C
         END IF
C
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END
C
C  /* Deck drive_fit */
       SUBROUTINE DRIVE_FIT (FMAT,DMAT,WMAT,WVEC,NDMT,ISYMDM,
     &                       IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                       WKOD,WORK,LWORK,CFIT,IPRINT)
C*****************************************************************************
C    
C     DRIVE_FIT: Driver for density fitting.
C        
C     Written by Luuk Visscher, november 2003
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
#include <denfit.h> 
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*), CFIT(*)
      DIMENSION KOD(12)
#include "nuclei.h"
#include "symmet.h"
C
#include <memint.h> 
C
C     Read the 2-index integrals from file
C
      CALL READ_FIT2INT(WMAT,WORK(KFREE),LFREE)
C
C     Allocate memory for the indexing of fit functions
C
      CALL MEMGET ('INTE',KINDAB,NBASISAUX,WORK,KFREE,LFREE)
C
C     Initialize the fit coefficients
C
      CALL DZERO (CFIT,NBASISAUX*NDMT)
C
C     **************************************************
C     ************* Carry out the fitting ************** 
C     **************************************************
C
      !If requested, loop over symmetry unique center pairs (diatom fit)
      IF (DIATOM) THEN
         NUM_OF_CENTS = NUCIND
      !else do global fitting
      ELSE
         NUM_OF_CENTS = 1
      END IF
C
      !Loop over symmetry unique centers
      DO ICENTB = 1, NUM_OF_CENTS
         !Loop over symmetry unique centers
         DO ICENTA = 1, ICENTB
            CALL FIT_DENSITY (ICENTA,ICENTB,FMAT,DMAT,WMAT,WVEC,NDMT,
     &                        ISYMDM,IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                        WKOD,WORK(KFREE),LFREE,CFIT,
     &                        WORK(KINDAB),IPRINT)
         END DO 
      END DO 
C
      RETURN
      END
C  /* Deck fit_density */
      SUBROUTINE FIT_DENSITY (ICENTA,ICENTB,FMAT,DMAT,WMAT,WVEC,NDMT,
     &                        ISYMDM,IFCTYP,CCFBT,INDXBT,KOD,FCKTMP,
     &                        WKOD,WORK,LWORK,CFIT,INDAB,IPRINT)
C*****************************************************************************
C    
C     FIT_DENSITY: Fit an atom pair density.
C        
C     Written by Luuk Visscher, november 2003
C     
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "priunit.h"
      DIMENSION FMAT(*),DMAT(*), WMAT(*), WVEC(*), WORK(LWORK), WKOD(*),
     &          IFCTYP(NDMT), ISYMDM(NDMT), CCFBT(*), INDXBT(*),
     &          FCKTMP(*), CFIT(*), INDAB(*)
      DIMENSION KOD(12)
#include "nuclei.h"
#include "symmet.h"
#include "inforb.h"
#include "denfit.h"
C
      DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)
C
#include <memint.h> 
C
C     We take as first index the fit set functions on either center A or B,
C     and as second index the null function. The remaining two are ordinary
C     basis functions on either center A or B.
C
      IF (DIATOM) THEN 
         !Index 1
         ICLASSES(1) = 3
         NCENTS(1)   = 2
         ICENTS(1,1) = ICENTA
         ICENTS(2,1) = ICENTB
         !Index 2
         ICLASSES(2) = 0
         NCENTS(2)   = 1
         ICENTS(1,2) = 0
         !Index 3
         ICLASSES(3) = 1
         NCENTS(3)   = 2
         ICENTS(1,3) = ICENTA
         ICENTS(2,3) = ICENTB
         !Index 4
         ICLASSES(4) = 1
         NCENTS(4)   = 2
         ICENTS(1,4) = ICENTA
         ICENTS(2,4) = ICENTB
      ELSE
         !Index 1
         ICLASSES(1) = 3
         NCENTS(1)   = 0
         ICENTS(1,1) = 0
         !Index 2
         ICLASSES(2) = 0
         NCENTS(2)   = 1
         ICENTS(1,2) = 0
         !Index 3
         ICLASSES(3) = 1
         NCENTS(3)   = 0
         ICENTS(1,3) = 0
         !Index 4
         ICLASSES(4) = 1
         NCENTS(4)   = 0
         ICENTS(1,4) = 0
      END IF
C
      CALL PICK_ACT_BATCH (NCENTS,ICENTS,ICLASSES,IPRFCK)
C
C     Initialize the weight vector and subsequently fill it by
C     contracting the active integrals with density matrix
C
      CALL DZERO (WVEC,NBASISAUX*NDMT)
      IOPT = 2
      !To prevent many fittings for a single atom in diatom scheme
      IF (DIATOM.AND.(ICENTA.EQ.ICENTB)) IOPT = -2
      CALL ERI_MKINTS (WVEC,DMAT,NDMT,ISYMDM,IFCTYP,IOPT,IPRFCK,
     &                 CCFBT,INDXBT,KOD,FCKTMP,
     &                 WKOD,WORK,LWORK)
C
C     Do some printing in case we're really desparate
C
      IF (IPRINT.GT.10) THEN
         IF (DIATOM) THEN
          WRITE (LUPRI,'(//8X,A,2I6//)') 'Fit atom pair:',ICENTA,ICENTB
         END IF
         CALL HEADER('Full weight matrix',-1)
         CALL OUTPUT(WMAT,1,NBASISAUX,1,NBASISAUX,NBASISAUX,NBASISAUX,
     &            1,LUPRI)
         CALL HEADER('Full weight vector',-1)
         CALL OUTPUT(WVEC,1,NBASISAUX,1,1,NBASISAUX,1,1,LUPRI)
      ENDIF
C
C     Determine how many fit functions we have (should sit on either
C     A or B) and make an index array to do gather/scatter. The indices
C     need to be shifted let them start at 1.
C
      CALL NFUN_ACT_BATCH (1,N_FIT,INDAB)
      CALL DF_INDEX_SHIFT(INDAB,N_FIT,INDAB,N_FIT,N_FIT,ICLASSES,1,1)
C
C     Allocate memory for temporary storage of A,B subblocks
C
      CALL MEMGET ('REAL',KJAB,N_FIT*N_FIT,WORK,KFREE,LFREE)
      CALL MEMGET ('REAL',KCFTAB,N_FIT,WORK,KFREE,LFREE)
C
C     Extract the weight vector for the pairs
C
      INDNULL = 1
      CALL EXTR_ACT_BATCH (WVEC,NBASISAUX,WORK(KCFTAB),N_FIT,
     &                     INDAB,N_FIT,INDNULL,1,1)
C
C     Extract the weight matrix from the full weight matrix.
C
      CALL EXTR_ACT_BATCH (WMAT,NBASISAUX,WORK(KJAB),N_FIT,
     &                     INDAB,N_FIT,INDAB,N_FIT,1)
C
C     Printing these matrices may sometimes be helpful
C
      IF (IPRINT.GT.4) THEN
         IF (DIATOM) THEN
          WRITE (LUPRI,'(//8X,A,2I6/)')'Fitting atom pair',ICENTA,ICENTB
         END IF
         CALL HEADER('Weight matrix',-1)
         CALL OUTPUT(WORK(KJAB),1,N_FIT,1,N_FIT,N_FIT,N_FIT,1,LUPRI)
         CALL HEADER('Weight vector',-1)
         CALL OUTPUT(WORK(KCFTAB),1,N_FIT,1,1,N_FIT,1,1,LUPRI)
      ENDIF
C
C     Solve the equation for the fit coefficients
C
      CALL DPOSV ( 'L', N_FIT, 1, WORK(KJAB), max(1,N_FIT),
     &             WORK(KCFTAB), max(1,N_FIT), INFO )
C
C     We again have some information that may be useful
C
      IF (IPRINT.GT.4) THEN
         WRITE (LUPRI,'(//8X,A,2I6//)') 'INFO code from DPOSV',INFO
         CALL HEADER('Fit coefficients',-1)
         CALL OUTPUT(WORK(KCFTAB),1,N_FIT,1,1,N_FIT,1,1,LUPRI)
      ENDIF
C
C     Add the fit coefficients to the full array
C
      CALL EXTR_ACT_BATCH (CFIT,NBASISAUX,WORK(KCFTAB),N_FIT,
     &                     INDAB,N_FIT,INDNULL,1,3)
C
      CALL MEMREL('FIT_DENSITY',WORK,KJAB,KJAB,KFREE,LFREE)
C
      RETURN
      END
C  /* Deck pick_act_batch */
      SUBROUTINE PICK_ACT_BATCH (NCENTS,ICENTS,ICLASSES,IPRINT)
C*****************************************************************************
C    
C     PICK_AOAB: Define active batches in density fitting call to ERI.
C        
C     Written by Luuk Visscher, october 2004
C     
C     Arrays needed to pick the active batches:
C     - NCENTS :   the number of active centers for a given index (maximally 2)
C     - ICENTS :   the corresponding center indices 
C     - ICLASSES : the active class (can only be one at a time).
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <mxcent.h>
#include <shells.h>
#include <erisel.h>
#include <aobtch.h>
      DIMENSION NCENTS(4),ICENTS(2,4),ICLASSES(4)
C
      DO I = 1, NAOBCH
         ICENT = NCNTBT(I)
         ICLASS = KCLSBT(I)
         DO J = 1, 4
C           The shell must be of the right class
            IF (ICLASS.EQ.ICLASSES(J)) THEN
C              The center must be active (ncents=0 means all centers active)
               IF (NCENTS(J).EQ.0) THEN
                  ACTVBT(I,J) = .TRUE.
               ELSEIF (NCENTS(J).EQ.1.OR.NCENTS(J).EQ.2) THEN
                  ACTVBT(I,J) = .FALSE.
                  DO K = 1, NCENTS(J)
                     IF (ICENTS(K,J).EQ.ICENT) ACTVBT(I,J) = .TRUE.
                  ENDDO
               ELSE
                  CALL QUIT ("Wrong value for NCENTS in PICK_ACT_BATCH")
               ENDIF
            ELSE
               ACTVBT(I,J) = .FALSE.
            ENDIF
         END DO
      END DO
C     
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Output from PICK_ACT_BATCH',-1)
         DO J = 1, 4
             IF (NCENTS(J).EQ.0) THEN
                WRITE (LUPRI,'(A,I1,A)')
     &         '  Index ',J,', all centers active'
             ELSE
                WRITE (LUPRI,'(A,I1,A,(10I5))')
     &         '  Index ',J,', active centers :',
     &           (ICENTS(K,J),K=1,NCENTS(J))
             ENDIF
         ENDDO
         WRITE (LUPRI,'(A,4I5,/)')
     &      '  Active classes for each index:',
     &      (ICLASSES(I),I=1,4)
         WRITE (LUPRI,'(A/)') '  ACTVBT(1:NAOBCH,1:4) in PICKAO:'
         DO I = 1, NAOBCH
            WRITE (LUPRI,'(I15,5X,4L5)') I, (ACTVBT(I,J),J=1,4)
         END DO
      END IF
C
      RETURN
      END
C  /* Deck nfun_act_batch */
      SUBROUTINE NFUN_ACT_BATCH (IND,N_FIT,INDEX)
C*****************************************************************************
C    
C     Counts number of functions in active batches for index IND
C     Sets up index array that can be used in gather/scatter
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <mxcent.h>
#include <shells.h>
#include <erisel.h>
#include <aobtch.h>
C
      DIMENSION INDEX(*)
C
      N_FIT = 0
      J = 0
      DO ISH = 1, NAOBCH
        DO I = 1, NORBBT(ISH)
           J = J + 1
           IF (ACTVBT(ISH,IND)) THEN
              N_FIT = N_FIT + 1
              INDEX(N_FIT) = J
           ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck extr_act_batch */
      SUBROUTINE EXTR_ACT_BATCH (FULL_MAT,N_FULL,RED_MAT,N_RED,
     &                           INDA,NA,INDB,NB,MODE)
C*****************************************************************************
C    
C     EXTR_ACTB: Extract/insert the active parts out of /in to a full matrix.
C     MODE 1 : Extraction
C          2 : Insertion
C          3 : Addition of reduced matrix to full matrix
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION FULL_MAT(N_FULL,*),RED_MAT(N_RED,*)
      DIMENSION INDA(*),INDB(*)
C
      IF (MODE.EQ.1) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               RED_MAT(I,J) = FULL_MAT(IF,JF)
            ENDDO
         ENDDO
      ELSEIF (MODE.EQ.2) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               FULL_MAT(IF,JF) = RED_MAT(I,J)
            ENDDO
         ENDDO
      ELSEIF (MODE.EQ.3) THEN
         DO J = 1, NB
            JF = INDB(J)
            DO I = 1, NA
               IF = INDA(I)
               FULL_MAT(IF,JF) = FULL_MAT(IF,JF) + RED_MAT(I,J)
            ENDDO
         ENDDO
      ELSE
         CALL QUIT ('Wrong MODE in EXTR_ACT_BATCH')
      ENDIF
C     
      RETURN
      END
C  /* Deck read_fit2int */
      SUBROUTINE READ_FIT2INT(WMAT,WORK,LWORK)
C
C     Wrapper routine for read of the full set of 2-index fit integrals
C
#include <implicit.h>
#include <iratdef.h>
#include <inforb.h>
#include <cbieri.h>
      DIMENSION WMAT(*),WORK(LWORK)
      LBUF = LBFINP
      KIBUF = 1
      KTOP  = KIBUF + 4*LBUF
      KLEFT = LWORK - KTOP + 1
      CALL READ_FIT2INT_2(WMAT,WORK(KIBUF),LBUF,
     &                    WORK(KTOP),KLEFT)
      RETURN
      END
C  /* Deck read_fit2int_2 */
      SUBROUTINE READ_FIT2INT_2(WMAT,IBUF,LBUF,WRK,LWORK)
C
C     Read the full set of 2-index fit integrals to memory.
C     The indices are shifted so that the first fit function is placed at position 1
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include "mxcent.h"
#include <nuclei.h>
#include <inforb.h>
      INTEGER A, B, C, D
      DIMENSION WMAT(NBASISAUX,NBASISAUX),
     &          IBUF(4,LBUF), MBAS(8), WRK(LWORK)
      DIMENSION ICLASSES(4)
#include <ibtfun.h>
C
      KFREE = 1
      LFREE = LWORK
C
C     We will read integrals of the type (nf|nf)
C     where |f> is a fit function and |n> the null function
      ICLASSES(1) = 0
      ICLASSES(2) = 3
      ICLASSES(3) = 0
      ICLASSES(4) = 3
C
C     Two-electron part
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'FIT_INTS','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM, (MBAS(I),I=1,8), KBUF, NIBUF, NBITS
      IF (KBUF .NE. LBUF) CALL QUIT('READ_FIT2INT_2: KBUF .ne. LBUF')
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('INTE',KINT,LENINT,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
C
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
      CALL DZERO(WMAT,NBASISAUX*NBASISAUX)
      JBUF = 608
  200 CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUINTA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LUINTA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,IBUF,NINT)
         IF (NINT .EQ. 0) GO TO 200
         IF (NINT .LT. 0) GO TO 400
         IF (NINT .GT. 0) THEN
C           Indices are relative to the fulllist, shift them to the aux. list
            CALL DF_INDEX_SHIFT(IBUF,LBUF,IBUF,LBUF,NINT,ICLASSES,4,2)
            DO 300 I = 1, NINT
               GINT  = WRK(KINT + I - 1)
               A = IBUF(1,I)
               B = IBUF(2,I)
               C = IBUF(3,I)
               D = IBUF(4,I)
               WMAT(B,D) = GINT
               WMAT(D,B) = GINT
  300       CONTINUE
         END IF
      GO TO 200
  400 CONTINUE
      CALL GPCLOSE(LUINTA,'KEEP')
      RETURN
      END
C  /* Deck df_index_shift */
      SUBROUTINE DF_INDEX_SHIFT(IBIN,LBIN,IBUF,LBUF,NINTS,ICLASSES,
     &                          NIND,MODE)
C*****************************************************************************
C    
C     Shifts the index for the auxilliary fit function such that the first fit
C     functions get index 1 and the last nbasisaux
C
C     MODE 1 : To be used in the Fock matrix construction where the labels are
C              organized as in IBIN
C     MODE 2 : To be used in the read of integrals were the labels are 
C              organized as in IBUF
C
C     ICLASSES : class type for ech of the four indices
C     IBIN/IBUF: integral labels that are to be shifted
C     LBIN/LBUF: length of the buffer
C     NIND:      number of indices that are to be shifted
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <aovec.h>
#include <mxcent.h>
#include <shells.h>
#include <nuclei.h>
#include <erisel.h>
#include <aobtch.h>
C
      DIMENSION IBIN(LBIN,NIND),IBUF(NIND,LBUF),IND_SHIFT(NIND)
      DIMENSION ICLASSES(NIND)
C
C     Determine the shifts for each index
C
      DO J = 1, NIND
         IF (ICLASSES(J).EQ.1) THEN
C           Nonrelativistic or large component functions
            IND_SHIFT(J) = 0
         ELSEIF (ICLASSES(J).EQ.2) THEN
C           Huckel or small component functions
            IND_SHIFT(J) = NLARGE
         ELSEIF (ICLASSES(J).EQ.3) THEN
C           Auxilliary functions for density fitting
            IND_SHIFT(J) = NLARGE + NSMALL
         ELSEIF (ICLASSES(J).EQ.0) THEN
C           The null function comes after the normal aux. functions
            IND_SHIFT(J) = NLARGE + NSMALL + NBASISAUX
         ENDIF
      END DO
C
      IF (MODE.EQ.1) THEN
         DO J = 1, NIND
            IF (IND_SHIFT(J).NE.0) THEN
               DO I = 1, NINTS
                  IBIN(I,J) = IBIN(I,J) - IND_SHIFT(J)
               END DO
            ENDIF
         END DO
      ELSE
         DO I = 1, NINTS
            DO J = 1, NIND
               IBUF(J,I) = IBUF(J,I) - IND_SHIFT(J)
            END DO
         END DO
      ENDIF
C
      RETURN
      END
C  /* Deck index_canon */
      SUBROUTINE INDEX_CANON(IBIN,LBIN,NINTS)
C*****************************************************************************
C    
C     Make pair of indices canonical (I>J)
C
C     IBIN: integral labels that are to be made canonical
C     LBIN: length of the buffer
C        
C     Written by Luuk Visscher, october 2004
C     
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION IBIN(LBIN,2)
C
      DO INT = 1, NINTS
         I = IBIN(INT,1)
         J = IBIN(INT,2)
         IBIN(INT,1) = MAX(I,J)
         IBIN(INT,2) = MIN(I,J)
      ENDDO
C
      RETURN
      END
C  /* Deck densfit_inp */
      SUBROUTINE DENSFIT_INP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <denfit.h>
      PARAMETER ( NTABLE = 5 )
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DATA TABLE /'.DIATOM','.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
C
      WORD1 = WORD
C
 100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)

            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5), I
                  END IF
 200          CONTINUE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FCK3INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FCK3INP.')

 1    CONTINUE
        !.DIATOM - fit density for pairs of atoms
        DIATOM     = .TRUE.
      GO TO 100
 2    CONTINUE
      GO TO 100
 3    CONTINUE
      GO TO 100
 4    CONTINUE
      GO TO 100
 5    CONTINUE
      GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized for *DENFIT.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 CALL QUIT('Illegal keyword for *DENFIT.')
            END IF
 300  CONTINUE

*     CALL HEADER('Settings for FCK3 calculation:',0)

*     IF (IPRINT .NE. IPRFDEF) THEN
*         WRITE (LUPRI,'(A,I5)') ' Print level in FCK3:',IPRFCK3
*     END IF

*     IF ( FIXTHR ) THEN
*       WRITE (LUPRI,'(A,I5)') 'Threshold is fixed.'
*     ELSE
*       WRITE (LUPRI,'(A,I5)') 'Threshold set dynamically.'
*     END IF

*     IF ( NOSKIPCL ) THEN
*       WRITE (LUPRI,'(A,I5)') 
*    &  'No skipping of the nonclassical integrals (only set to 0).'
*     END IF

      RETURN
      END
C  /* Deck densfit_ini */
      SUBROUTINE DENSFIT_INI
#include <implicit.h>
#include <denfit.h>

      !.DIATOM - fit density for pairs of atoms (default is the global fit)
      DIATOM  = .FALSE.

      RETURN
      END
