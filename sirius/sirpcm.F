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
C FILE: sirius/sirpcm.F
C
#include <single.h>
C
C/* Deck pcminp */
      SUBROUTINE PCMINP(WORD,INPERR,ALLOPT)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
C
      PARAMETER (D0=0.0D+00)
C
C     ----- set up the NAMELIST $PCM -----
C
      PARAMETER (NTABLE=16)
      CHARACTER*8 RZSOL
      CHARACTER*6 SOLWD
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7
      LOGICAL ALLOPT
#include <nuclei.h>
#include <pcmnuclei.h>
#include <infinp.h>
#include <infpri.h>
#include <inftap.h>
      DATA SOLWD/'PCM   '/
      DATA TABLE/'.SOLVNT','.EPS   ','.EPSINF','.RSOLV ',
     *           '.RET   ','.NESFP ','.ICESPH','.ICOMPC',
     *           '.PRINT ','.NON-EQ','.NEQRSP','.NPCMMT',
     *           'xANISO ','.LOCFLD','.OUTFLD','.NEWQR '/
C
C
C     --- READ IN PARAMETERS TO PERFORM A PCM SOLVATION CALCULATION ---
C
C     The default solvent is assumed to be WATER at 25 C, and the
C     cavity is built as usual, with a sphere of scaled atomic radius
C     around each solute atom.
C
C       this group defines the solvent (no default)
C
      PCM    = .TRUE.
      NEWQR  = .FALSE.
      RZSOL  = 'INPUT   '
      EPSI   = D0
      EPSINI = D0
      RSOLVI = D0
      NONEQ  = .FALSE.
      NEQRSP = .FALSE.
C
C Correction level of the PCM matrices (see file sirief)
C
C NPCMMT = 0 No correction of the DI, SI and C matrices
C NPCMMT = 1 Correction of DI and SI  (default)
C NPCMMT = 2 Correction of DI, SI and C
C
      NPCMMT = 1
C
C       this group defines the cavity
C
      ICESPH = 0
      NESFP  = 0
      OMEGA  = 40.0D+00
      FRO    = 0.7D+00
      RET    = 0.2D+00
      NRWCAV = 0
      AREATS = 0.4D+00
Clf oldcen true means that we use the old tessera center definition
C     (average of the vertices).
      OLDCEN = .FALSE.
      ICOMPCM = 0
      IPRPCM = 0
C
C     Do the input processing the Dalton way
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108,109,110,
     &                   111,112,113,114,115,116), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
               GO TO 100
            END IF
            WRITE (6,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
            GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 999
         ELSE
            WRITE (6,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
         CALL QUIT(' ILLEGAL KEYWORD IN PCMINP ')
C
C *** Option 1 >RZSOL<  Solvent type
  101 READ(LUCMD,'(A8)') RZSOL
      GO TO 100
C
C *** Option 2 >EPSI<  Static Dielectric constant
  102 READ(LUCMD,*) EPSI
      EPS = EPSI
      GO TO 100
C
C *** Option 3 >EPSINI<  Optical Dielectric constant
  103 READ(LUCMD,*) EPSINI
      EPSINF = EPSINI
      GO TO 100
C
C *** Option 4 >RSOLV <  Solvent radius
  104 READ(LUCMD,*) RSOLVI
      IF (RSOLVI.LE.D0) THEN 
         CALL QUIT('PCM INPUT ERROR: RSOLV MUST BE POSITIVE')
      END IF
      GO TO 100
C
C *** Option 5 >RET<  Minimum Radius of added spheres
  105 READ(LUCMD,*) RET
      GO TO 100
C
C *** Option 6 >NESFP< number of spheres
  106 READ(LUCMD,*) NESFP
      GO TO 100
C
C *** Option 7 >ICESPH< option to define sphere's center
  107 READ(LUCMD,*) ICESPH
      GO TO 100
C
C *** Option 8 >ICOMPCM< Charge renormalization
  108 READ(LUCMD,*) ICOMPCM
      GO TO 100
C
C *** Option 9 >PRINT <  Print level
  109 READ(LUCMD,*) IPRPCM
      GO TO 100
C
C *** Option 10 >NON-EQ <  Non-equilibrium solvation (PCM)
 110  NONEQ = .TRUE.
      GO TO 100
C
C *** Option 11 >NEQRSP <  Non-equilibrium response  (PCM)
 111  NEQRSP = .TRUE.
      GO TO 100
C
C *** Option 12 >NPCMMT <  Numerical correction level of PCM matrices
 112  READ(LUCMD,*) NPCMMT
      GO TO 100
C
C *** Option 13 >ANISO <   Anisotropic solvent with dielectric tensor
 113  READ(LUCMD,*) EPSXX, EPSYY, EPSZZ
      call quit('anisotropic solvent not yet implemented')
Clf      LANISO = .TRUE.
      GO TO 100
C
C *** Option 14 >LOCFLD < Local field effect
 114  LOCFLD = .TRUE.
      GO TO 100
C
C *** Option 15 >OUTFLD < Local field effect for sum frequency generation
C Ref: Mizrahi and Sipe -- Phys. Rev. B Vol 34 n. 6 pp. 3700-3709 
C
 115  LOCFLD = .TRUE.
      OUTFLD = .TRUE.
      GO TO 100
C
C *** Option 16 >NEWQR <  Go for the new quadratic response routine
 116  NEWQR = .TRUE.
      GO TO 100
C
C
C     We need to create to PCM files
C
 999  CONTINUE
C
C        No longer an input value, DR is the distance away from the
C        representative point of the tessera along the surface normal
C
      DR = 1.0D-04
C
C        look up numerical parameters for the solvent, default is water
C        any values which were input should override stored values
C
      IF(RZSOL.NE.'INPUT   ') THEN
         CALL DATSOL(RZSOL,EPS,EPSINF,RSOLV,VMOL,TCE,STEN,DSTEN,CMF)
      END IF
      IF(EPSI  .NE.D0) EPS    = EPSI
      IF(EPSINI.NE.D0) EPSINF = EPSINI
      IF(RSOLVI.NE.D0) RSOLV  = RSOLVI
C
C        print what we've got for input
C
      WRITE(LUPRI,9020) ICOMPCM, RZSOL,EPS,EPSINF,RSOLV,
     *     ICESPH,NESFP,OMEGA,RET,FRO,IPRPCM,NONEQ,NEQRSP
C
C        Check that input for booboos
C
      NERR=0
      IF(RZSOL.EQ.'INPUT   ' .AND. (EPS.EQ.D0.OR.RSOLV.EQ.D0)) THEN
         WRITE(LUPRI,9010) 'PICK SOLVNT OR RSOLV/EPS/EPSINF...'
         NERR = NERR+1
      END IF
      IF(NESFP.GT.MXSP) THEN
         WRITE(LUPRI,9010) 'EXCESSIVE NUMBER OF RADII (>MXSP)'
         NERR=NERR+1
      END IF
      IF(ICESPH.GT.3) THEN
         WRITE(LUPRI,9010) 'WRONG VALUE FOR -ICESPH- (>3)'
         NERR=NERR+1
      END IF
      IF(ICOMPCM.GT.3) THEN
         WRITE(LUPRI,9010) 'WRONG VALUE FOR -ICOMPCM- (>3)'
         NERR=NERR+1
      END IF
C
      IF(NERR.GT.0) THEN
         WRITE(LUPRI,*) 'PLEASE FIX THE ABOVE ERROR(S) IN *PCM'
         CALL QUIT('Errors detected in PCM input')
      END IF
C
C        read information defining the cavity spheres
C
      CALL MAKCAV(WORD,INPERR,ALLOPT)
C
      RETURN
C
 9010 FORMAT(' *** INCONSISTENCY FOUND IN *PCM INPUT GROUP ***'/5X,A)
 9020 FORMAT(/5X,35('-')/
     *   5X,'Input for PCM solvation calculation '/5X,35('-')/
     *   5X,'ICOMPCM =',I8,  5X,
     *   5X,'SOLVNT=',A8,  5X,'EPS   =',F8.4,5X,'EPSINF=',F8.4/
     *   5X,'RSOLV =',F8.4//
     *   5X,'ICESPH =',I8,  5X,'NESFP =',I8/
     *   5X,'OMEGA =',F8.4,5X,'RET   =',F8.4,5X,'FRO   =',F8.4//
     *   5X,'IPRPCM=',I8//5X,'NON-EQ = ',L1,5X,'NEQRSP =',L1)
      END
C/* Deck pcmsolvnt */
      SUBROUTINE PCMSOLVNT(WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
ckr
#include <aovec.h>
#include <primit.h>
ckr
#include <pcmdef.h>
C
      LOGICAL SOME,FNDLAB
C
#include <pcm.h>
#include <pcmlog.h>
#include <inforb.h>
#include <infpri.h>
#include <infvar.h>
C
      DIMENSION WORK(*)
C
C     ----- set up a PCM-BEM calculation -----
C
C     The common blocks used by the PCM package all begin with PCM...
C       /PCMABC/ reaction field grid
C       /PCMCAV/ cavity sphere info       size=          2*mxsp
C       /PCMCHG/ surface charge info      size= 2*mxts
C       /PCMDBS/ dispersion basis set
C       /PCMDAT/ input control params
C       /PCMDIS/ dispersion info
C       /PCMPLY/ polyhedron definitions   size= 1*mxts + 5*mxsp
C       /PCMPAR/ input control params
C       /PCMPRT/ final print results
C       /PCMREP/ repulsion info
C       /PCMTES/ tessera data             size= 8*mxts
C       /PCMTHF/ time dependent HF data
C       /PCMUGG/ surface energy data
C     Those listed without a size are very small.
C     Total fixed memory demand is roughly 7*MXSP + 11*MXTS.
C
C     The disk files used by the PCM package are
C       PCMDATA - NFT26 - All runs have the D inverse matrix, while
C                         gradient runs also have the DERPUNT, DERTES,
C                         DERCENTR, and DERRAD arrays present.
C
C     The file 26 data cannot be stored on DAF, unfortunately, as the
C     number of tesserae on the cavity surface changes during the run.
C
C
C     The memory usage by the PCM package is controlled by 2 parameters,
C       MXSP = maximum number of spheres in the cavity.  This should be
C              similar to the number of atoms being treated, although a
C              methyl group, for example, might be in a single sphere.
C       MXTS = maximum number of cavity tessera.  This should be roughly
C              25 times MXSP, more or less.
C     The 'mung' script can be used to resize MXSP or MXTS, but the
C     other three require hand editing of the source routines.
C
      L2 = NNBASX
      L3 = NNBASX*NNBASX
      SOME = IPRPCM.LE.2
      IF(SOME) WRITE(LUPRI,900)
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C Memory allocation
C
      NTS3 = NTS*NTS
      LWRK = LWORK
      KFREE = 1
      CALL MEMGET('REAL',LDMATM,NTS3,    WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LSE1,  NTS3,    WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LDE1,  NTS3,    WORK,KFREE,LWRK)
Clf      IF(LANISO) THEN
      IF(.FALSE.) THEN
         CALL MEMGET('REAL',LSE2,  NTS3,  WORK,KFREE,LWRK)
         CALL MEMGET('REAL',LDE2,  NTS3,  WORK,KFREE,LWRK)
      ELSE
         LSE2=LSE1
         LDE2=LDE1
      END IF
      CALL MEMGET('REAL',LWRK1, NTS,     WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LIPVT, NTS,     WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LMEP,  NTS,     WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LCAM,  NTS,     WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LBMCHG,NNBASX,  WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LVEC0 ,NTS,     WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LVEC, NNBASX,WORK,KFREE,LWRK)
      CALL MEMGET('REAL',LVEC2,NNBASX,WORK,KFREE,LWRK)
C
      CALL DZERO(WORK,KFREE)
C

ClfC
ClfC     We always assume we will do gradient calculations
ClfC
Clf      NDER = 0
Clf      IF(NDER.EQ.1) THEN
Clf         READ (LUPCMD)
Clf         READ (LUPCMD)
Clf      END IF
C     CCVEBF writes the electric field integrals to disk file -NFT27-.
C
      NTNIRR = NTS * NTSIRR
      LUBKUP = LUPCMD
      IF (LUPCMD .LT. 0) CALL GPOPEN(LUPCMD,'PCMDATA','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUPCMD)
      IF (FNDLAB('Q-PCMMAT',LUPCMD)) THEN
         CALL READT(LUPCMD,NTNIRR,WORK(LDMATM))
      ELSE
         WRITE (LUPRI,'(/A)') ' Matrix label Q-PCMMAT not '//
     &        'found on file PCMDATA'
         CALL QUIT('Integral label not found in PCMSOLVNT')
      END IF
      CALL ICVEV(WORK(LDMATM),WORK(LSE1),WORK(LSE1),
     $     WORK(LBMCHG),WORK(LVEC0),WORK(LMEP),WORK(LCAM),NTS,NNBASX,
     $     SOME,WORK(KFREE),LWRK)
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C     4) Calculation of the interaction free energy between
C        solute's nuclei-nuclear induced charges.
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      CALL PCMVNN
      IF(SOME) THEN
         WRITE(LUPRI,'(/A)')
     &   ' ..... DONE WITH INDUCED NUCLEAR CHARGES .....'
         CALL TIMER('PCMVNN   ',TIMSTR,TIMEND)
      ELSE
C
         WRITE(LUPRI,'(/A)')
     &   ' ..... DONE SETTING UP PCM CALCULATION .....'
      END IF
      RETURN
C
  900 FORMAT(/10X,'-------------------------------------'/
     *        10X,'---- POLARISABLE CONTINUUM MODEL ----'/
     *        10X,'----      UNIVERSITY OF PISA     ----'/
     *        10X,'-------------------------------------')
      END
C/* Deck makcav */
      SUBROUTINE MAKCAV(WORD,INPERR,ALLOPT)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <nuclei.h>
#include <pcmnuclei.h>
#include <infpri.h>
#include <inftap.h>
ckr
#include <aovec.h>
#include <maxorb.h>
#include <primit.h>
ckr
C
C LASTAR: LAST ATOMIC RADIUS DEFINED
      PARAMETER (LASTAR=86)
      PARAMETER (D0=0.0D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
C
C     ----- set up the NAMELIST $PCMCAV -----
C
      PARAMETER (NTABLE=8)
      CHARACTER WORD*7, WORD1*7, PROMPT*1, TABLE(NTABLE)*7,
     &          CAVWD*7
      LOGICAL ALLOPT
      DATA CAVWD/'*PCMCAV'/
      DATA TABLE/'.CENTER','.RIN   ','.ALPHA ','.INA   ',
     & '.POLYG ','.RWCAV ','.AREATS','.OLDCEN'/
C
C        sphere centers XE,YE,ZE and radii RIN must be input in Angstrom
C
      CALL DZERO(XE,MXSP)
      CALL DZERO(YE,MXSP)
      CALL DZERO(ZE,MXSP)
      CALL DZERO(RIN,MXSP)
      CALL DZERO(ALPHA,MXSP)
      CALL IZERO(INA,MXSP)
      CALL IZERO(IAN,MXCENT)
      POLYG=60.0D0
C
      
Clf      IF(ICESPH.LE.0) THEN
Clf         DO J=1,NESFP
Clf            XE(J)=PCMCORD(1,J)*TOANGS
Clf            YE(J)=PCMCORD(2,J)*TOANGS
Clf            ZE(J)=PCMCORD(3,J)*TOANGS
Clf            NUCZ = INT(PCMCHG(J))
Clf            RIN(J) = RVDW(NUCZ)
Clf       	   IAN(J)=1
Clf         ENDDO
Clf      END IF
C
C
C     Do the input processing the Dalton way
C
      WORD1 = WORD
      IF (ALLOPT) CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
  100    READ (LUCMD, '(A7)') WORD
         CALL UPCASE(WORD)
         PROMPT = WORD(1:1)
         IF (PROMPT .EQ. '.') THEN
            DO 1102 II = 1, NTABLE
               IF (TABLE(II) .EQ. WORD) THEN
                  GO TO (101,102,103,104,105,106,107,108), II
               END IF
 1102       CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
               IF (.NOT.ALLOPT)
     &            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
               GO TO 100
            END IF
            WRITE (6,'(/4A/)') ' Keyword ',WORD,
     *         ' not recognized for ',WORD1
         ELSE IF (PROMPT .EQ. '#' .OR. PROMPT .EQ. '!') THEN
               GO TO 100
         ELSE IF (PROMPT .EQ. '*') THEN
            GOTO 999
         ELSE
            WRITE (6,'(/3A/2A/)')
     *         ' Keyword ',WORD,' does not begin with',
     *         ' one of the four characters ".*!#" for ',WORD1
         END IF
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',6)
         CALL QUIT(' ILLEGAL KEYWORD IN PCMCAV ')
C
C *** Option 1 >CENTER <  Sphere Centers
  101 IF(ICESPH.EQ.1)THEN
        DO II=1,NESFP
          READ(LUCMD,*) XE(II),YE(II),ZE(II)
        ENDDO
      ENDIF
      GO TO 100
C
C *** Option 2 >RIN <  Sphere radii
  102 DO II=1,NESFP
        READ(LUCMD,*) RIN(II)
      ENDDO
      GO TO 100
C
C *** Option 3 >ALPHA <  Scaling factor
  103 DO II=1,NESFP
        READ(LUCMD,*) ALPHA(II)
      ENDDO
      GO TO 100
C
C *** Option 4 >INA   <  Nuclei associated to spheres
  104 IF(ICESPH.EQ.2)THEN
      DO II=1,NESFP
        READ(LUCMD,*) INA(II)
      ENDDO
      END IF
      GO TO 100
C
C *** Option 5 >POLYG <  Choice of polyhedra for tessellation
  105 READ(LUCMD,*) POLYG
      GO TO 100
C
C     *** Option 6 >RWCAV <  Read cavity from or write cavity to file
 106  READ(LUCMD,*) NRWCAV
      IF (NRWCAV .LT. 0 .OR. NRWCAV .GE. 3) THEN
         WRITE(LUPRI,*) '.RWCAV=',NRWCAV,' VALUE NOT ALLOWED'
         WRITE(LUPRI,*) 'ALLOWED VALUES ARE 1 (READ CAVITY) OR ',
     $        '2 (WRITE CAVITY)'
         CALL QUIT('SIRPCM: INPUT ERROR')
      ENDIF
      GO TO 100
C
C     *** Option 7 >AREATS<  Define the area of a tessera in the initial tesselation of a sphere
 107  READ(LUCMD,*) AREATS
      GO TO 100
C
C     *** Option 8 >OLDCEN<  Ask for old style center of tessera definition
 108  OLDCEN = .TRUE.
      write(lupri,*) 'read oldcen', oldcen
      GO TO 100
C
C
 999  CONTINUE
Clf      IF(ICESPH.EQ.2)THEN
Clf      DO J=1,NESFP
Clf       	XE(J)=PCMCORD(1,INA(J))*TOANGS
Clf       	YE(J)=PCMCORD(2,INA(J))*TOANGS
Clf       	ZE(J)=PCMCORD(3,INA(J))*TOANGS
Clf       	IAN(INA(J))=1
Clf      ENDDO
Clf      ENDIF
      
C Positive numbers indicate the number of faces of the polyhedra to be used in the tessellation.
C    ES: 60, 120, 240, ...
C Negative numbers indicate the average area of the tessera according to which 
C the number of faces of the polyhedra is chosen :
C  Es: -0.4 (angstrom^2)
C (this second option is active only for calculations without symmetry, i.e. C1).
      IF(POLYG.GT.0.0D0)THEN
        IPOLYG=INT(POLYG)
      ELSE
        IPOLYG=INT(POLYG*1000)
      ENDIF
      WRITE(LUPRI,*)'    POLYG',IPOLYG
C

      WRITE(LUPRI,9010)
      NERR = 0
      IF(ICESPH .LT. 3) THEN
         DO J=1,NESFP
            WRITE(LUPRI,9020) J,XE(J),YE(J),ZE(J),RIN(J)
            IF(RIN(J).LE.0) NERR=NERR+1
         ENDDO
      END IF
      IF(NERR.GT.0) THEN
         WRITE(LUPRI,9030)
         CALL QUIT('Inconsistent input in PCMCAV')
      END IF
C
C        on exit, XE,YE,ZE should be converted to Bohr, but RIN stays A.
C
      DO J=1,NESFP
         XE(J)=XE(J)/TOANGS
         YE(J)=YE(J)/TOANGS
         ZE(J)=ZE(J)/TOANGS
      ENDDO
C
C     If parameter ALPHA(I) is greater than 0 the I-th radius is
C     multiplied by ALPHA(I).
C     This allows to use radii=R(van der Waals)*ALPHA in the calculation
C     of electrostatic (and eventually SCF disp-rep) contribution, and
C     radii = RvdW for the cavitation.
C     It is possible to give a different value of ALPHA to each atom:
C     if the first one is less than D0, it is fixed equal to 1,
C     otherwise it keeps its value. Each following value less than D0
C     is changed to 1.
C
      IF(ALPHA(1).LE.0.0D+00) ALPHA(1) = 1.2D+00
      DO I = 2, MXSP
        IF(ALPHA(I).LE.0.0D+00) ALPHA(I) = ALPHA(1)
      ENDDO
      RETURN
C
C9000 FORMAT(' **** ERROR IN $',A8,' INPUT')
 9010 FORMAT(/5X,'INPUT FOR CAVITY DEFINITION '/5X,27(1H-)/
     *   5X,'ATOM         COORDINATES           RADIUS ')
 9020 FORMAT(2X,I4,4(2X,F8.4))
 9030 FORMAT(' **** ERROR **** PCM SPHERE(S) MUST HAVE',
     *          ' A POSITIVE RADIUS')
      END
C/* Deck pcmvnn */
      SUBROUTINE PCMVNN
C
#include <implicit.h>
#include <pcmdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
C
#include <nuclei.h>
#include <pcmnuclei.h>
#include <symmet.h>
#include <pcm.h>
#include <pcmlog.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D0)
C
C     Routine for calculation of free energy of interaction between
C     solute nuclei and the polarization charges induced by the
C     nuclei themselves.
C
      PCMNN = D0
      DO ITS = 1, NTSIRR
        XI=XTSCOR(ITS)
        YI=YTSCOR(ITS)
        ZI=ZTSCOR(ITS)
        DO JATOM = 1, NUCDEP
           R2 = (PCMCORD(1,JATOM)-XI)**2 + (PCMCORD(2,JATOM)-YI)**2 + 
     *          (PCMCORD(3,JATOM)-ZI)**2
           R = DSQRT(R2)
           PCMNN = PCMNN + QSN(ITS) * PCMCHG(JATOM) / R
        ENDDO
      ENDDO
      PCMNN = (MAXREP + 1) * (PCMNN * 0.5D0)
      RETURN
      END
C/* Deck pcmjmat */
      SUBROUTINE PCMJMAT(VEC2,NUM2,WORK,LWORK)
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <inftap.h>
C
      CHARACTER*8 RTNLBL(2)
      LOGICAL FNDLAB
      DIMENSION VEC(NUM2),VEC2(NUM2)
C
#include <pcm.h>
#include <pcmlog.h>
C
C     This routine computes J matrix  (interaction between induced
C     nuclear charges and electronic potential due to the solute):
C         J(m,n) = sum[ V(m,n;i)qsn(i) ]
C     and if IREP=1 h_(rep) matrix (repulsive interaction between
C     solute and solvent):
C          h_(rep)(m,n) = gamma*(S(m,n)+sum_i[E_i(m,n)*a(i)]/4PI)
C     (S=overlap matrix, E_i=matrix of normal components of the field)
C     These matrices are used to correct 1-e integrals obtained without
C     solvent (they are stored in file 11, or for MCSCF in file 87).
C
      CALL DZERO(VEC2,NUM2)
C
      LUBKP = LUPROP
      CALL J1INT(QSN,.FALSE.,VEC2,1,.FALSE.,'NPETES ',1,WORK,LWORK)
      IF (LUPROP .LT. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
C        
      REWIND (LUPROP)
C
C     Write modified 1-e integrals
C
      RTNLBL(2) = 'SYMMETRI'
      CALL WRTPRO(VEC2,NUM2,'NE-PCMIN',RTNLBL,IPRPCM)
      IF (LUBKP .LT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
      RETURN
      END
C     /* Deck datsol */
      SUBROUTINE DATSOL(ZSOL,EPS,EPSINF,RSOLV,VMOL,TCE,STEN,DSTEN,CMF)
C     
#include <implicit.h>
#include <priunit.h>
#include <infpri.h>
C
C To add a new solvent:
C 1) Update the NSOL parameter
C 2) Insert the appropriate name and formula in the data structures
C 3) Add a new label in the goto statement 
C 4) Add a new section with the solvent parameters
C     
      PARAMETER (NSOL=18)
      CHARACTER*8 ZSOL, NAMSOL(NSOL), LABSOL(NSOL)
      DATA LABSOL    /'H2O     ','CH3OH   ','C2H5OH  ','CHCL3   ',
     $     'CH2CL2  ','C2H4CL2 ','CCL4    ','C6H6    ','C6H5CH3 ',
     $     'C6H5CL  ','CH3NO2  ','C7H16   ','C6H12   ','C6H5NH2 ',
     $     'CH3COCH3','THF     ','DMSO    ','CH3CN   '/
      DATA NAMSOL    /'WATER   ','METHANOL','ETHANOL ','CLFORM  ',
     $     'METHYLCL','12DICLET','TETRACLC','BENZENE ','TOLUENE ',
     $     'CLBENZ  ','NITROMET','N-EPTANE','CYCLOHEX','ANILINE ',
     $     'ACETONE ','TETHYDFU','DIMETSOX','ACETONIT'/
C     
C     Database of optical and physical data for various solvents
C     
      DO I = 1, NSOL
         IF ((LABSOL(I) .EQ. ZSOL) .OR. (NAMSOL(I) .EQ. ZSOL)) THEN
            GO TO (101,102,103,104,105,106,107,108,
     &           109,110,111,112,113,114,115,116,117,118), I
         
         END IF
      END DO
      GOTO 300
C     
C     Water solvent
C     
 101  CONTINUE
      EPS = 78.39D+00
      EPSINF = 1.776D+00
      RSOLV = 1.385D+00
      VMOL = 18.07D+00
      TCE = 2.57D-04
      STEN = 71.81D+00
      DSTEN = 0.650D+00
      CMF = 1.277D+00
      GOTO 200
C     
C     METHANOL
C     
 102  CONTINUE
      EPS = 32.63D+00
      EPSINF = 1.758D+00
      RSOLV = 1.855D+00
      VMOL = 40.7D+00
      TCE = 1.182D-03
      STEN = 22.12D+00
      DSTEN = 1.154D+00
      CMF = 1.776D+00
      GOTO 200
C     
C     ETHANOL       
C     
 103  CONTINUE
      EPS = 24.55D+00
      EPSINF = 1.847D+00
      RSOLV = 2.180D+00
      VMOL = 58.7D+00
      TCE = 1.103D-03
      STEN = 21.89D+00
      DSTEN = 1.146D+00
      CMF = 1.543D+00
      GOTO 200
C     
C     CHLOROFORM
C     
 104  CONTINUE
      EPS = 4.90D+00
      EPSINF = 2.085D+00
      RSOLV = 2.48D+00
      VMOL = 80.7D+00
      TCE = 1.255D-03
      STEN = 26.53D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     METHYLENECHLORIDE
C     
 105  CONTINUE
      EPS = 8.93D+00
      EPSINF = 2.020D+00
      RSOLV = 2.27D+00
      VMOL = 64.5D+00
      TCE = 1.367D-03
      STEN = 27.33D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     1,2-DICHLOROETHANE
C     
 106  CONTINUE
      EPS = 10.36D+00
      EPSINF = 2.085D+00
      RSOLV = 2.505D+00
      VMOL = 79.4D+00
      TCE = 1.156D-03
      STEN = 31.54D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     CARBON TETRACHLORIDE
C     
 107  CONTINUE
      EPS = 2.228D+00
      EPSINF = 2.129D+00
      RSOLV = 2.685D+00
      VMOL = 96.5D+00
      TCE = 1.270D-03
      STEN = 26.15D+00
      DSTEN = 1.436D+00
      CMF = 0.629D+00
      GOTO 200
C     
C     BENZENE
C     
 108  CONTINUE
      EPS = 2.247D+00
      EPSINF = 2.244D+00
      RSOLV = 2.63D+00
      VMOL = 88.91D+00
      TCE = 1.380D-03
      STEN = 28.18D+00
      DSTEN = 1.469D+00
      CMF = 0.629D+00
      GOTO 200
C     
C     TOLUENE
C     
 109  CONTINUE
      EPS = 2.379D+00
      EPSINF = 2.232D+00
      RSOLV = 2.82D+00
      VMOL = 106.3D+00
      TCE = 1.08D-03
      STEN = 27.92D+00
      DSTEN = 1.391D+00
      CMF = 0.679D+00
      GOTO 200
C     
C     CHLOROBENZENE
C     
 110  CONTINUE
      EPS = 5.621D+00
      EPSINF = 2.320D+00
      RSOLV = 2.805D+00
      VMOL = 101.79D+00
      TCE = 0.981D-03
      STEN = 32.69D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     NITROMETHANE
C     
 111  CONTINUE
      EPS = 38.20D+00
      EPSINF = 1.904D+00
      RSOLV = 2.155D+00
      VMOL = 53.68D+00
      TCE = 1.192D-03
      STEN = 36.47D+00
      DSTEN = 1.373D+00
      CMF = 0.808D+00
      GOTO 200
C     
C     N-EPTHANE
C     
 112  CONTINUE
      EPS = 1.92D+00
      EPSINF = 1.918D+00
      RSOLV = 3.125D+00
      VMOL = 146.56D+00
      TCE = 1.25D-03
      STEN = 19.80D+00
      DSTEN = 1.505D+00
      CMF = 0.687D+00
      GOTO 200
C     
C     CYCLOHEXANE
C     
 113  CONTINUE
      EPS = 2.023D+00
      EPSINF = 2.028D+00
      RSOLV = 2.815D+00
      VMOL = 108.10D+00
      TCE = 1.20D-03
      STEN = 24.38D+00
      DSTEN = 1.467D+00
      CMF = 0.621D+00
      GOTO 200
C     
C     ANILINE
C     
 114  CONTINUE
      EPS = 6.89D+00
      EPSINF = 2.506D+00
      RSOLV = 2.80D+00
      VMOL = 91.15D+00
      TCE = 0.85D-03
      STEN = 42.79D+00
      DSTEN = 0.731D+00
      CMF = 0.972D+00
      GOTO 200
C     
C     ACETONE
C     
 115  CONTINUE
      EPS = 20.7D+00
      EPSINF = 1.841D+00
      RSOLV = 2.38D+00
      VMOL = 73.52D+00
      TCE = 1.42D-03
      STEN = 22.67D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     TETRHYDROFURANE
C     
 116  CONTINUE
      EPS = 7.58D+00
      EPSINF = 1.971D+00
      RSOLV = 2.9D+00
      VMOL = 81.11D+00
      TCE = 1.142D-03
      STEN = 26.40D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     DIMETHYLSULFOXIDE
C     
 117  CONTINUE
      EPS = 46.7D+00
      EPSINF = 2.179D+00
      RSOLV = 2.455D+00
      VMOL = 70.94D+00
      TCE = 9.82D-02
      STEN = 42.86D+00
      DSTEN = 0.0D+00
      CMF = 0.0D+00
      GOTO 200
C     
C     ACETONITRILE
C     
 118  CONTINUE
      EPS = 36.64d0
      EPSINF = 1.806d0
      RSOLV = 2.155
      VMOL = 53.68d0
      TCE = 1.192d-3
      STEN = 36.47d0
      DSTEN = 1.373d0
      CMF = 0.808d0
      GOTO 200
C
C We arrive here if we don't match the solvent
C
 300  CONTINUE
      WRITE(LUPRI,9010) ZSOL
      DO I =1, NSOL
         WRITE(LUPRI,9015) LABSOL(I),NAMSOL(I)
      ENDDO
      CALL QUIT('Unrecognized solvent in DATSOL')

 200  CONTINUE
      WRITE(LUPRI,9020) ZSOL,EPS,EPSINF,RSOLV,VMOL,TCE,STEN,DSTEN,CMF
C     
 9010 FORMAT(/'*** ERROR: NO DATA TABULATED FOR SOLVENT = ',A8
     *      //'ALLOWED SOLVENTS ARE:')
 9015 FORMAT('      ',A10,' ( OR ',A10,' )')
C     
 9020 FORMAT(/' ** LOOKING UP INTERNALLY STORED DATA FOR SOLVENT = ',
     *     A8,' **'
     *    /' Optical and physical constants:'
     *    /' EPS=',F7.3,'; EPSINF=',F7.3,';',
     *     ' RSOLV=',F7.3,' A; VMOL=',F8.3,' ML/MOL;'
     *    /' TCE=',1P,D11.5,0P,' 1/K; STEN=',F7.3,
     *     ' DYN/CM;  DSTEN=',F8.4,'; CMF=',F8.4/)
      RETURN
      END
C/* Deck v2q */
      SUBROUTINE V2Q(CINVMT,TJ2,QCREAT,QET,NONEQQ)
C
C     Transform electronic potential to electronic charges
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      PARAMETER (D0 = 0.0D0)
Clf      DIMENSION CINVMT(NTS*NTS), TJ2(NTS)
      DIMENSION CINVMT(NTSIRR*NTS), TJ2(*), QCREAT(*)
      LOGICAL QCHECK, FNDLAB, NONEQQ
#include <pcm.h>
#include <pcmlog.h>
C
      LPCMBK = LUPCMD
      NTNIRR = NTS * NTSIRR
      IF (LUPCMD .LT. 0) CALL GPOPEN(LUPCMD,'PCMDATA','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUPCMD)
      IF (NONEQQ) THEN
         IF (FNDLAB('Q-DINMAT',LUPCMD)) THEN
            CALL READT(LUPCMD,NTNIRR,CINVMT)
         ELSE
            WRITE (LUPRI,'(/A)') ' Matrix label Q-DINMAT not '//
     &           'found on file PCMDATA'
            CALL QUIT('Integral label not found in V2Q')
         END IF
      ELSE
         IF (FNDLAB('Q-PCMMAT',LUPCMD)) THEN
            CALL READT(LUPCMD,NTNIRR,CINVMT)
         ELSE
            WRITE (LUPRI,'(/A)') ' Matrix label Q-PCMMAT not '//
     &           'found on file PCMDATA'
            CALL QUIT('Integral label not found in V2Q')
         END IF
      END IF
      SQTNOP = DBLE(MAXREP + 1)
      DO ISYM = 0, MAXREP
         ISTART = ISYM * NTSIRR ** 2 + 1
         JSTART = ISYM * NTSIRR + 1
         CALL DGEMV('N',NTSIRR,NTSIRR,1.0D0,CINVMT(ISTART),NTSIRR,
     $        TJ2(JSTART),1,D0,QCREAT(JSTART),1)
      ENDDO
      QET = D0
      DO ITS = 1, NTS
         QCREAT(ITS) = QCREAT(ITS)*AS(ITS)/SQTNOP
         QET = QET + QCREAT(ITS)
      END DO
C
C SYMMETRY RENORMALIZATION
C
      QET = SQTNOP * QET
C
      RETURN
      END
C
C/* Deck iefmat */
C
      SUBROUTINE IEFMAT(WORK,LWRK)
C Author: Luca Frediani
C
C Date: 09/29/2002
C
C Purpose: Driver of all cavity matrices needed for a PCM calculation
C
C C-INVMAT: Matrix to compute q = MV for equilibrium calculations
C C-INVDIN: Matrix to compute q = MV for non-equilibrium response calculations
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      PARAMETER (D0 = 0.0D0, D1=1.0D0)
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION WORK(*)
Clf      LUPCMBK = LUPCMD
      LUPCMD = -10009
      IF (LUPCMD .LT. 0) THEN
         CALL GPOPEN(LUPCMD,'PCMDATA','UNKNOWN','SEQUENTIAL',
     $        'UNFORMATTED',IDUMMY,.FALSE.)
      END IF
      REWIND (LUPCMD)
      KFRSAV = 1
      LFRSAV = LWRK
      KFREE  = KFRSAV
      LFREE  = LWRK

      LWORK  = LWRK
      CALL PEDRAM(WORK,LWORK)
Clf   IF(SOME) THEN
         WRITE(LUPRI,'(/A)')
     &   ' ..... DONE GENERATION CAVITY .....'
         CALL TIMER('PEDRAM',TIMSTR,TIMEND)
Clf   END IF
      NTS2 = NTS*NTS
      CALL FLSHFO(LUPRI)
      CALL MEMGET('REAL',LDMATM,NTS2,    WORK,KFREE,LWORK)
      CALL MEMGET('REAL',LSE1,  NTS2,    WORK,KFREE,LWORK)
      CALL MEMGET('REAL',LDE1,  NTS2,    WORK,KFREE,LWORK)
Clf      IF(LANISO) THEN
      IF(.FALSE.) THEN
         CALL MEMGET('REAL',LSE2,  NTS2,  WORK,KFREE,LWORK)
         CALL MEMGET('REAL',LDE2,  NTS2,  WORK,KFREE,LWORK)
      ELSE
         LSE2=LSE1
         LDE2=LDE1
      END IF
      CALL MEMGET('REAL',LWRK1, NTSIRR,  WORK,KFREE,LWORK)
      CALL MEMGET('REAL',LIPVT, NTSIRR,  WORK,KFREE,LWORK)
C
C     Calculation of D matrix and its inverse.
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL IEFCMMV(WORK(LDMATM),WORK(LDMATM),WORK(LSE1),WORK(LSE1),
     $     WORK(LDE1),WORK(LDE1),WORK(LWRK1),WORK(LIPVT),WORK(KFREE),
     $     LWORK)
C     
CLF      SUBROUTINE WRTIEF(AINT,LENGTH,LABEL,LUNIT)

      NTNIRR = NTS * NTSIRR
      CALL WRTIEF(WORK(LDMATM),NTNIRR,'Q-PCMMAT',LUPCMD)
Clf      IF(SOME) THEN
      IF(.TRUE.) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
     &        ' ..... DONE GENERATING -Q-  MATRIX .....'
         CALL TIMER('Q-MAT ',TIMSTR,TIMEND)
      END IF
      IF (NEQRSP) THEN
         EPSOLD=EPS
         EPS=EPSINF
         CALL IEFCMMV(WORK(LDMATM),WORK(LDMATM),WORK(LSE1),WORK(LSE1),
     $        WORK(LDE1),WORK(LDE1),WORK(LWRK1),WORK(LIPVT),WORK(KFREE),
     $        LWORK)
         EPS=EPSOLD
         CALL WRTIEF(WORK(LDMATM),NTNIRR,'Q-DINMAT',LUPCMD)
      ENDIF
C      IF (OUTFLD) THEN
C         CALL WRTIEF(WORK(LDMATM),NTNIRR,'D-MATRIX',LUPCMD)
C      END IF
      IF (LOCFLD.AND.(.NOT.OUTFLD)) THEN
         IF (NEQRSP) THEN
            FACT=(EPSINF+D1)/(EPSINF-D1)
         ELSE
            FACT=(EPS+D1)/(EPS-D1)
         END IF
         CALL DM1LOC(WORK(LDE1),WORK(LSE1),WORK(LWRK1),WORK(LIPVT),
     $        AS,FACT,NTSIRR,MAXREP,IPRPCM)
         CALL WRTIEF(WORK(LDE1),NTNIRR,'D-LOCINV',LUPCMD)
      END IF
      CALL GPCLOSE(LUPCMD,'KEEP')
      CALL MEMREL('IEFMAT',WORK,KFRSAV,KFRSAV,KFREE,LWORK)
      RETURN
      END
C
C/*DECK DM1LOC*/
C
      SUBROUTINE DM1LOC(DPAK,SPAK,WORK,IPVT,AS,FACT,NTSIRR,MAXREP,
     $     IPRPCM)
C
C 1/10/2002 L. Frediani
C Calculation and inversion of the matrix FOR local field calculations
C with PCM
C
#include <implicit.h>
#include <pi.h>
      PARAMETER (D0=0.0D0,DP5=0.50D0,FPI=4.0D0*PI)
#include <priunit.h>
#include <inftap.h>
      DIMENSION DPAK(NTSIRR,*),SPAK(NTSIRR,*),WORK(NTSIRR),IPVT(*)
      DIMENSION AS(*)
      NTS = (MAXREP + 1)*NTSIRR
      DO JSYM = 0, MAXREP
         JSTART = JSYM * NTSIRR
         DO I=1,NTSIRR
            DO J = 1,NTSIRR
               JTS = J + JSTART
C     
C     symmmetry adapted transpose matrix indeces
C     
               ITSTR = J
               JTSTR = I + JSTART
               
               DELTAIJ = D0
               IF (I.EQ.J) then 
                  DELTAIJ = AS(I) * FACT * DP5
               END IF
               SPAK(I,JTS) = FPI*(DELTAIJ - DPAK(ITSTR,JTSTR))/AS(I)
            END DO
         END DO
         DO I=1,NTSIRR
            DO J=1,NTSIRR
               JTS = J + JSTART
               DPAK(I,JTS) = -SPAK(I,JTS)
            END DO
         END DO
      END DO
      IF(IPRPCM.GE.10) THEN
         WRITE(LUPRI,*) '========== LF MATRIX PACKED   =========='
         CALL OUTPUT(DPAK,1,NTSIRR,1,NTS,NTSIRR,NTS,1,LUPRI)
      END IF
      DO ISYM = 0, MAXREP
         ISTART = NTSIRR * ISYM + 1
         INFO = 0
         CALL DGEFA(DPAK(1,ISTART),NTSIRR,NTSIRR,IPVT,INFO)
         IF(INFO.NE.0) THEN
            WRITE(LUPRI,*) ' THE DI MATRIX IS SINGULAR'
            CALL QUIT('Singular DI matrix in J4LOC')
         END IF
         CALL DGEDI(DPAK(1,ISTART),NTSIRR,NTSIRR,IPVT,DET,WORK,01)
      ENDDO
      IF(IPRPCM.GE.6) THEN
         WRITE(LUPRI,*) '========== LF MATRIX INV PACK =========='
         CALL OUTPUT(DPAK,1,NTSIRR,1,NTS,NTSIRR,NTS,1,LUPRI)
      END IF
      RETURN
      END
C
C/*Deck J1INT*/
C
      SUBROUTINE J1INT(EXPVAL,EXP1VL,DENMAT,NOSIM,TOFILE,INTLAB,
     &                 KSYMP,WORK,LWORK)
C
C
C     DENMAT : if (EXP1VL) then input density matrix; else output matrix
C     
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <codata.h>
#include <infpri.h>
#include <nuclei.h>
#include <inforb.h>
#include <pcmnuclei.h>
#include <pcm.h>
#include <pcmlog.h>
#include <orgcom.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include "infpar.h"
#include "inftap.h"
#include "krcom.h"
C
      CHARACTER*7 INTLAB
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT), EXPVAL(NTS,NOSIM)
      DIMENSION DENMAT(*), WORK(LWORK)
      LOGICAL TOFILE, TRIMAT, VCHECK, EXP1VL
C
Clf Do not insert other definitions below the following include
#include <ibtfun.h>
C
      NBASTOLD  = NBAST
      NNBASXOLD = NNBASX
      NBAST  = ISUM(MAXREP+1,NAOS,1)
      NNBASX = NBAST * (NBAST + 1)/2
      N2BASX = NBAST*NBAST
      IF (INTLAB .EQ. 'PCMBSOL') THEN
         MATDIM = N2BASX
         IPRTYP = 13
         NCOMP  = 3
         TRIMAT = .FALSE.
      ELSE
         MATDIM = NNBASX
         IPRTYP = 11
         NCOMP  = 1
         TRIMAT = .TRUE.
      END IF
C
C     We use as a quick way of transfering tessera coordinates to HERMIT
C     the dipole origin. Need to be restored.
C
      XDIPORG = DIPORG(1)
      YDIPORG = DIPORG(2)
      ZDIPORG = DIPORG(3)
C
C  2) Calculation of apparent charges generated by the solute's nuclei.
C
      IF (EXP1VL) THEN
         KDEN = 1
         KLAST = KDEN + NNBASX
         LWRK  = LWORK - KLAST
         ISYMD = KSYMP - 1
         IF (KSYMP .EQ. 1) THEN
            CALL PKSYM1(DENMAT,WORK(KLAST),NBAS,NSYM,1)
            CALL DSYM1(WORK(KDEN),DUMMY,WORK(KLAST),DUMMY,.FALSE.,
     &                 NBAST,IPRPCM)
         ELSE
            CALL DCOPY(NNBASX,DENMAT,1,WORK(KDEN),1)
         END IF
         IF (NOSIM .GT. 1) CALL QUIT('NOSIM .GT. 1 and EXP1VL not '//
     &        'permitted in J1INT')
      ELSE
         KDEN = 1
         IF (NODTOT .GE. 1) THEN
            KLAST = KDEN + MATDIM*NOSIM
            CALL DCOPY(MATDIM*NOSIM,DENMAT,1,WORK(KDEN),1)
         ELSE
            KLAST = KDEN
         END IF
         LWRK = LWORK - KLAST
      END IF
#if defined (VAR_MPI)
      IF (NODTOT .GE. 1) THEN
         CALL J1INTP(NBAST,NOSIM,KSYMP,WORK(KDEN),EXP1VL,EXPVAL,TOFILE,
     &               IPRTYP,MATDIM,WORK(KLAST),LWRK)
         IF (.NOT. EXP1VL) CALL DAXPY(MATDIM*NOSIM,1.0D0,WORK(KDEN),1,
     &                                DENMAT,1)
      ELSE
#endif
      DO  ITS = 1, NTSIRR
         DIPORG(1) = XTSCOR(ITS)
         DIPORG(2) = YTSCOR(ITS)
         DIPORG(3) = ZTSCOR(ITS)
         NTESP = 1
         KPATOM = 0
C
C        Calculates nuclear potential energy integrals (in AO basis) for
C        the given tessera
C
         L=1
         KTMP = KLAST
         IF (.NOT. TOFILE .AND. .NOT. EXP1VL) THEN
            KMAT = KTMP + 8
            IF (IPRTYP .EQ. 11) THEN
               KLAST = KMAT + (MAXREP + 1)*MATDIM
            ELSE
               KLAST = KMAT + (MAXREP + 1)*MATDIM*NCOMP
            END IF
            NCOMP = NSYM
         ELSE
            KMAT  = KTMP + 8
            KLAST = KMAT
            NCOMP = 0
         END IF
         CALL GET1IN(WORK(KMAT),INTLAB,NCOMP,WORK(KLAST),LWRK,LABINT,
     &               INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,WORK(KTMP),
     &               EXP1VL,WORK(KDEN),IPRPCM)
         IF (IPRTYP .EQ. 13) THEN
            JMAT = KMAT
            DO IOSIM = 1, NOSIM
               CALL DAXPY(MATDIM,EXPVAL(ITS,IOSIM),WORK(JMAT),1,
     &                    DENMAT(MATDIM*(IOSIM - 1) + 1),1)
               JMAT = JMAT + MATDIM
            END DO
         ELSE IF (EXP1VL) THEN
            DO I = 1, NCOMP
               EXPVAL(ITS+(I-1)*NTSIRR,1) = -WORK(KTMP+I-1)
            END DO
         ELSE IF (.NOT. TOFILE) THEN
            DO IOSIM = 1, NOSIM
               IADR = KMAT + (KSYMP - 1)*MATDIM
               CALL DAXPY(MATDIM,-EXPVAL(ITS,IOSIM),WORK(IADR),1,
     &                    DENMAT(MATDIM*(IOSIM - 1) + 1),1)
            END DO
         END IF
      ENDDO 

#if defined (VAR_MPI)
      END IF
#endif
C     
      DIPORG(1) = XDIPORG
      DIPORG(2) = YDIPORG
      DIPORG(3) = ZDIPORG
      NBAST  = NBASTOLD
      NNBASX = NNBASXOLD
      RETURN
      END
C
C/*Deck PCMSPHGEN*/
C
      SUBROUTINE PCMSPHGEN
C
C Generation of the initial set of spheres if not defined in input
C Luca Frediani 15 11 2004
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <nuclei.h>
#include <pcmnuclei.h>
#include <pcmdef.h>
#include <pcm.h>
#include <codata.h>
#include <symmet.h>
      PARAMETER (LASTAR=99)
      PARAMETER (D0=0.0D0)
      DIMENSION RVDW(LASTAR)
#include <ibtfun.h>
C
C       van der Waals radii taken from "The Elements", 2nd edition,
C       John Emsley, Clarendon Press, Oxford, 1991.  Unknown values
C       are simply set to D0, rather than trying to guess them.
C
C       A.Bondi, J.Phys.Chem. 68: 441-451(1964) gives alternate
C       values, and a few transition metals.
C
      DATA (RVDW(II),II=1,18)/1.20D+00,1.22D+00,
     *      0.00D+00,0.00D+00,2.08D+00,1.85D+00,
     *      1.54D+00,1.40D+00,1.35D+00,1.60D+00,
     *      2.31D+00,0.00D+00,2.05D+00,2.00D+00,
     *      1.90D+00,1.85D+00,1.81D+00,1.91D+00/
      DATA (RVDW(II),II=19,36)/2.31D+00,13*0.0D+00,
     *      2.00D+00,2.00D+00,1.95D+00,1.98D+00/
      DATA (RVDW(II),II=37,54)/2.44D+00,13*0.0D+00,
     *      2.20D+00,2.20D+00,2.15D+00,0.00D+00/
      DATA (RVDW(II),II=55,99)/2.62D+00,27*0.0D+00,2.40D+00,16*0.0D+00/
C
C     override the above table with U. Pisa's experience
C     as to what works best for singly bonded C,N,O
C
      RVDW(6) = 1.70D+00
      RVDW(7) = 1.60D+00
      RVDW(8) = 1.50D+00
Clf guessed Li and Be radii for some testing
Clf      RVDW(3) = 1.45D+00
Clf      RVDW(4) = 1.45D+00
C
C

      I=0
      DO N=1,NUCIND
         I=I+1
         PCMCORD(1,I)=CORD(1,N)
         PCMCORD(2,I)=CORD(2,N)
         PCMCORD(3,I)=CORD(3,N)
         PCMCHG(I)=CHARGE(N)
      END DO
      DO K = 1, MAXREP
         DO N = 1, NUCIND
            IF (IBTAND(K,ISTBNU(N)) .EQ. 0) THEN
               DS = (CORD(1,N)-CORD(1,N)*PT(IBTAND(ISYMAX(1,1),K)))**2
     &            + (CORD(2,N)-CORD(2,N)*PT(IBTAND(ISYMAX(2,1),K)))**2
     &            + (CORD(3,N)-CORD(3,N)*PT(IBTAND(ISYMAX(3,1),K)))**2
               DIST = SQRT(DS)
               IF (DIST.LT.0.1D0) GOTO 5000
               I=I+1
               IF (I.GT.MXCENT) GOTO 5010
               PCMCORD(1,I)=CORD(1,N)*PT(IBTAND(ISYMAX(1,1),K))
               PCMCORD(2,I)=CORD(2,N)*PT(IBTAND(ISYMAX(2,1),K))
               PCMCORD(3,I)=CORD(3,N)*PT(IBTAND(ISYMAX(3,1),K))
               PCMCHG(I)=CHARGE(N)
               RSPH(I) = RSPH(N)
               !print *, 'nucdep',i,rsph(i)
            END IF
         END DO
      END DO

      IF(I.NE.NUCDEP) CALL QUIT('PCMSPHGEN: WRONG NUMBER OF NUCLEI')

Clf      IPCM = 0
ClfC
ClfC GENERATION OF ALL COORDINATES 
ClfC
Clf      DO 100 ICENT = 1, NUCIND
Clf         MULCNT = ISTBNU(ICENT)
Clf         IF (MULT(MULCNT) .EQ. 1) THEN
Clf            IPCM = IPCM + 1
Clf            NATOMS = NATOMS + 1
Clf            PCMCORD(1,IPCM) = CORD(1,ICENT)
Clf            PCMCORD(2,IPCM) = CORD(2,ICENT)
Clf            PCMCORD(3,IPCM) = CORD(3,ICENT)
Clf            PCMCHG(IPCM) = CHARGE(ICENT)
Clf         ELSE
Clf            JATOM = 0
Clf            DO 200 ISYMOP = 0, MAXOPR
Clf               IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
Clf                  IPCM = IPCM + 1
Clf                  JATOM  = JATOM + 1
Clf                  NATOMS = NATOMS + 1
Clf                  PCMCORD(1,IPCM) = PT(IBTAND(ISYMAX(1,1),ISYMOP))
Clf     *                        *CORD(1,ICENT)
Clf                  PCMCORD(2,IPCM) = PT(IBTAND(ISYMAX(2,1),ISYMOP))
Clf     *                        *CORD(2,ICENT)
Clf                  PCMCORD(3,IPCM) = PT(IBTAND(ISYMAX(3,1),ISYMOP))
Clf     *                        *CORD(3,ICENT)
Clf                  PCMCHG(IPCM) = CHARGE(ICENT)
Clf               END IF
Clf  200       CONTINUE
Clf         END IF
Clf  100 CONTINUE
ClfC
ClfC GENERATION OF ALL SPHERES FROM COORDINATES
ClfC 
C ICESPH = 0 SPHERES ON ALL ATOMS WITH TABULATED RADII                  
C ICESPH = 1 SPHERES GIVEN FROM INPUT (X,Y,Z,R)                    
C ICESPH = 2 SPHERES ON SELECTED NUCLEI WITH RADII GIVEN FROM INPUT
C ICESPH = 3 Spheres from mol file (Sphere=1,2...)

      IF(ICESPH.EQ.0) THEN
         NESFP = I
         DO J=1,NESFP
            XE(J)=PCMCORD(1,J)
            YE(J)=PCMCORD(2,J)
            ZE(J)=PCMCORD(3,J)
            NUCZ = NINT(PCMCHG(J))
            IF(NUCZ.GT.LASTAR) CALL QUIT('PCMSPHGEN: NUCZ IS TOO BIG')
            RIN(J) = RVDW(NUCZ)
            IF(RIN(J).EQ.D0) THEN 
               CALL QUIT('PCMSPHGEN: SPHERE NOT DEFINED')
            end if
         ENDDO
      ELSE IF(ICESPH.EQ.1) THEN
         WRITE(LUPRI,'(A)') 'PCMSPHGEN: SPHERE CENTERS FROM INPUT'
      ELSE IF(ICESPH.EQ.2) THEN
         DO J=1,NESFP
            XE(J)=PCMCORD(1,INA(J))
            YE(J)=PCMCORD(2,INA(J))
            ZE(J)=PCMCORD(3,INA(J))
         ENDDO
      ELSE IF(ICESPH.EQ.3) THEN
         NESFP=0
         DO J=1,NUCDEP
            IF(RSPH(J).GT.0.1) THEN
               NESFP = NESFP + 1
               INA(NESFP) = J
               RIN(NESFP) = RSPH(J)
               ALPHA(NESFP) = 1.0D0
               XE(NESFP)  = PCMCORD(1,J)
               YE(NESFP)  = PCMCORD(2,J)
               ZE(NESFP)  = PCMCORD(3,J)
            END IF
         END DO
      ELSE 
         CALL QUIT('PCMSPHGEN: WRONG ICESPH VALUE')
      ENDIF
      NESF = NESFP
CLF   IF(IPRCAV.GE.5) THEN
         WRITE(LUPRI,*) '********SPHERES IN PCMSPHGEN************'
         WRITE(LUPRI,*) 
     $        'INDEX        X        Y         Z        R'
         DO N=1,NESF
            WRITE(LUPRI,1000) N,XE(N),YE(N),ZE(N),RIN(N)
         END DO
CLF   END IF
      RETURN

 5000 CONTINUE
      WRITE (LUPRI,'(2(A,I5),A,1P,D15.5)')
     &     ' The nucleus',N,' is too close to its',K,
     &     'th transformation DIST=',DIST
      CALL QUIT('PCMSPHGEN: Nuclei too close.')
 5010 CONTINUE
      CALL QUIT('PCMSPHGEN: MXCENT exceeded.')
 1000 FORMAT(I4,1P,4D20.10)
      END
C  --- end of sirius/sirpcm.F ---
