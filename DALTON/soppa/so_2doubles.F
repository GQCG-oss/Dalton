C  /* Deck so_2doubles */
      SUBROUTINE SO_2DOUBLES(OMEGA2,T2AM,T2AMSQ,WORK,LWORK,FOCKD,LFOCKD)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Written by Lilli Irene Ør Kristensen, Summer 2017.
C
C     PURPOSE: Calculate the 2nd order double excited correction to the
C               wavefunction.
C
C     LIST OF VARIABLES:
C       FOCK:   Orbital energies
C       LFOCK:  Lenght of the FOCK array
C       OMEGA2: contains the 2nd order doubles amplitudes.
C       T2AM:   Contains the 1st order doubles amplitudes.
C       T2AMSQ: Contains the squared 1st order doubles amplitudes.
C
C
#include "implicit.h"
#include "priunit.h"
C      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C      PARAMETER (FOUR = 4.0D0)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccisao.h"
#include "blocks.h"
#include "ccfield.h"
#include "ccsections.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "eribuf.h"
#include "ccnoddy.h"
#include "cbirea.h"
      DIMENSION OMEGA2(6*6*NT2SQ(1)),T2AM(NT2AMX),T2AMSQ(NT2SQ(1))
C      DIMENSION WORK(*),FOCKD(*)
C
C
C--------------------------------------------------------------------------
C       Calculate the contibutions in the amplitude.
C-------------------------------------------------------------------------

     CALL SO_2NDDECOEF(OMEGA2,T2AM,T2AMSQ,WORK,LWORK)

C-------------------------------------------------------------------------
C       Devide with the orbital energies
C--------------------------------------------------------------------

      CALL SO_MP2T2(FOCKD,LFOCKD,OMEGA2,36*NT2SQ(1))

      RETURN
      END
