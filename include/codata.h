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
#include <alphac.h>
C
      PARAMETER ( XTANG = 0.529 177 2083 D0, XFAMU  = 1822.88848 D0,
     &            ECHARGE = 1.602 176 462D-19, HBAR = 1.054 571 60D-34)
      PARAMETER ( XFMOL = 6.022 141 99 D23 )
C -- conversion from hartree (au) to:
      PARAMETER ( XTJ  = 4.359 743 81 D-18,  XTKAYS = 219474.6313710D0,
     *            XTHZ = 6.579 683 920 735 D15, XTEV = 27.211 3834 D0,
     &            XKJMOL = XTJ*XFMOL*1.D-3, XKCMOL = XKJMOL/4.184 D0,
     &            XTKJML = XKJMOL,          XTKCML = XKCMOL,
     *            XTNM = 1.D7/XTKAYS,       XAJOUL = 1.0 D18*XTJ)
C
C -- other
      PARAMETER ( AUTIME = HBAR/XTJ, XFSEC = AUTIME)
      PARAMETER ( XTKMML = 974.864 D0)
      PARAMETER ( TESLA=(XTANG*XTANG*ECHARGE/HBAR)*1.D-20 )
      PARAMETER ( AUTK = 3.157 7465 D5 )
C
