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
#include <pi.h>
C
      PARAMETER ( XTANG = 0.529 177 2083 D0, XFAMU  = 1822.88848 D0,
     &            ECHARGE = 1.602 176 462D-19, HBAR = 1.054 571 60D-34)
      PARAMETER ( XFMOL = 6.022 141 99 D23 )
      PARAMETER ( XTANGM10 = XTANG*1.0D-10)
C
      PARAMETER ( PMASS = 1.007276470D0, EMASS = 9.10938188D-31)
C
#include <alphac.h>
C -- conversion from hartree (au) to:
      PARAMETER ( XTJ  = HBAR**2/(XTANGM10*XTANGM10*EMASS),
     &            XTKAYS = 1.0D-2*HBAR/(CCM*2.0*PI*XTANGM10**2*EMASS),
     *            XTHZ =  HBAR/(2.0*PI*XTANGM10*XTANGM10*EMASS),
     &            XTEV = 27.211 3834 D0,
     &            XKJMOL = XTJ*XFMOL*1.D-3, XKCMOL = XKJMOL/4.184 D0,
     &            XTKJML = XKJMOL,          XTKCML = XKCMOL,
     *            XTNM = 1.D7/XTKAYS,       XAJOUL = 1.0 D18*XTJ)
C
C -- other
      PARAMETER ( AUTIME = HBAR/XTJ, XFSEC = AUTIME)
      PARAMETER ( XTKMML = 974.864 D0)
      PARAMETER ( TESLA=(XTANG*XTANG*ECHARGE/HBAR)*1.D-20 )
      PARAMETER ( AUTK = 3.157 7465 D5 )
      PARAMETER ( DEBYE = ECHARGE*XTANG*CCM*1D11 )
C
