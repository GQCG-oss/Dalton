C 
C  codata.h revised 1999/04/29 and combined with units.h 
C
C
C  From
C  "CODATA Recommended Values of the Fundamental Physical Constants: 1998"
C                     Peter J. Mohr and Barry N. Taylor
C  Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999
C
C -- Fundamental constants, atomic units, Avogadro's constant
#if defined (SYS_CRAY)
      REAL CVEL, ALPHAC, ALPHA2
      REAL XTANG,XFAMU,ECHARGE,HBAR,XFMOL
      REAL XTJ,XTKAYS,XTHZ,XTEV,XKJMOL,XKCMOL,XTKJML,
     &     XTKCML,XTNM,XAJOUL,XTANGM10,XPRTMAS 
      REAL XFSEC,XTKMML,TESLA,AUTK,DEBYE,PMASS,EMASS,CCM
#else
      DOUBLE PRECISION CVEL, ALPHAC, ALPHA2
      DOUBLE PRECISION XTANG,XFAMU,ECHARGE,HBAR,XFMOL
      DOUBLE PRECISION XTJ,XTKAYS,XTHZ,XTEV,XKJMOL,XKCMOL,XTKJML,
     &                 XTKCML,XTNM,XAJOUL,XTANGM10,XPRTMAS
      DOUBLE PRECISION XFSEC,XTKMML,TESLA,AUTK,DEBYE,PMASS,EMASS,CCM
#endif
#include <pi.h>
C
      PARAMETER ( XTANG = 0.529 177 2083 D0, XFAMU  = 1822.88848 D0,
     &            ECHARGE = 1.602 176 462D-19, HBAR = 1.054 571 596D-34)
      PARAMETER ( XFMOL = 6.022 141 99 D23 )
      PARAMETER ( XTANGM10 = XTANG*1.0D-10)
C
      PARAMETER ( PMASS = 1.007276470D0, EMASS = 9.10938188D-31)
C
#include <alphac.h>
C -- conversion from hartree (au) to:
      PARAMETER ( XTJ  = HBAR**2/(XTANGM10*XTANGM10*EMASS),
     &            XTHZ =  HBAR/(2.0D0*PI*XTANGM10*XTANGM10*EMASS),
     &            XTKAYS = 1.0D-2*XTHZ/CCM,
     &            XTEV = XTJ/ECHARGE,
     &            XKJMOL = XTJ*XFMOL*1.D-3, XKCMOL = XKJMOL/4.184 D0,
     &            XTKJML = XKJMOL,          XTKCML = XKCMOL,
     &            XTNM = 1.D7/XTKAYS,       XAJOUL = 1.0 D18*XTJ)
C
C -- other
      PARAMETER ( XFSEC = HBAR/XTJ)
      PARAMETER ( XTKMML = 974.864 D0)
      PARAMETER ( TESLA=(XTANG*XTANG*ECHARGE/HBAR)*1.D-20 )
      PARAMETER ( AUTK = 3.157 7465 D5 )
      PARAMETER ( DEBYE = ECHARGE*XTANG*CCM*1.D11 )
      PARAMETER ( XPRTMAS = 1836.152 6675 D0 )   
C
