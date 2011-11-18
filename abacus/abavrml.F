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
#ifdef UNDEF
===========================================================================
/* Comdeck log */
960319-vebjornb: New module for visualization with VRML models.
===========================================================================
#endif
C  /* Deck mkvrml */
      SUBROUTINE MKVRML(LAST, ATMARR, IEDIM, EVEC, EVC1, EVC2)
C     
C     Punch out geometry to VRML-file
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#include <gnrinf.h>
#include <cbiwlk.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL LAST
      DIMENSION ATMARR(MXCENT,8), EVEC(IEDIM,IEDIM)
      DIMENSION EVC1(MXCOOR), EVC2(MXCOOR)
      CHARACTER*12 FILENM
#include <ibtfun.h>
C
C     Initialize the ATMARR array. The first index runs over all
C     atoms, the second marks the following properties:
C
C               1 - Element number
C               2 - X coordinate of atom
C               3 - Y coordinate of atom
C               4 - Z coordinate of atom
C               5 - Covalent radius
C
      CALL ATMINI(ATMARR,IATOM,.FALSE.)
C
C     We proceed to open the output-file
C
      LUVRML = -1
      FILENM = 'first.wrl'
      IF (LAST) THEN
         LUVRML = -1
         FILENM = 'last.wrl'
      END IF
      CALL VRINI(LUVRML,FILENM)
C
C     We create all the atoms
C
      INDX = 1
      CALL DRWATM(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
C
C     Then we draw bonds between the atoms if this is requested.
C
      IF (VRBOND) CALL DRWBND(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
C
C     If the coordinate axes are requested they are drawn
C
      IF (VRCORD) CALL DRWAXS(LUVRML,INDX,IATOM,ATMARR,EVEC)
      CALL VREND(LUVRML)
C
C     Finally all eigenvectors are visualized
C
      IF (VREIGV) THEN
         LUVRML = -1
         IEIG = 1
         IEVEC = 1
         IF (IPRINT .GT. 0) THEN
            WRITE(LUPRI,*)
            CALL HEADER('VRML Visualization of Eigenvectors',-1)
            WRITE(LUPRI,'(A)')
     &           ' Eig.Vec.     Filename      Symmetry      Eig.Value '
            WRITE(LUPRI,'(A)')
     &           '----------------------------------------------------'
         END IF
C     
C     We loop over all symmetries...
C
         DO 30 IREP = 0, MAXREP
            IF (DOREPW(IREP)) THEN
               II = 0
               NCR = NCRREP(IREP,1)
               DO 35 I = 0, IREP - 1
                  II = II + NCRREP(I,1)
 35            CONTINUE
C
C     ... and all vectors in each symmetry
C
               DO 40 IVEC = 1, NCR
C
C     Only eigenvectors with a non-zero eigenvalue is visualized,
C     that is only the eigenvectors associated with.
C
                  IF (EVAL(IEVEC+IVEC-1) .GT. 1.0D-3) THEN
                     DO 50 I = 1, NCR
                        EVC1(I) = EVEC(II+I,II+IVEC)
 50                  CONTINUE
                     INDX = 1
                     FILENM = 'eigv_XXX.wrl'
                     WRITE(FILENM(6:8),'(I3)') IEIG
                     IF (IEIG .LT. 100) WRITE(FILENM(6:6),'(A1)') '0'
                     IF (IEIG .LT. 10) WRITE(FILENM(7:7),'(A1)') '0'
                     CALL VRINI(LUVRML,FILENM)
                     CALL DRWATM(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
                     IF (VRBOND)
     &                    CALL DRWBND(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
                     CALL DRWEIG(LUVRML,INDX,ATMARR,EVC1,EVC2,IREP,NCR)
                     CALL VREND(LUVRML)
                     IF (IPRINT .GT. 0) THEN
                        WRITE(LUPRI,'(I5,A,A12,A,I1,A,F16.6)')
     &                       IEIG,'       ',FILENM,'       ',
     &                       IREP,'   ',EVAL(IEVEC+IVEC-1)
                     END IF
                     IEIG = IEIG + 1
                  END IF
 40            CONTINUE
            END IF
            IEVEC = IEVEC + NCRREP(IREP,1)
 30      CONTINUE
      END IF
      RETURN
      END

C  /* Deck mkvrvb */
      SUBROUTINE MKVRVB(NCORD,IATOM,GVEC,EVEC,ATCHRG,
     &     MODENR,FREQ,WORK,LWORK)
C     
C     Make VRML representation of normal modes.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
      DIMENSION GVEC(NCORD), EVEC(NCORD), ATCHRG(NCORD), WORK(LWORK)
      KATMAR = 1
      KLAST = KATMAR + 8*MXCENT
      IF (KLAST .GT. LWORK) CALL STOPIT('MKVRVB',' ',KLAST,LWORK) 
      CALL MKVRVB_1(NCORD,IATOM,GVEC,EVEC,ATCHRG,MODENR,FREQ,
     &     WORK(KATMAR))
      RETURN
      END

C  /* Deck MKVRVB_1 */
      SUBROUTINE MKVRVB_1(NCORD,IATOM,GVEC,EVEC,ATCHRG,
     &     MODENR,FREQ,ATMARR)
C     
C     See MKVRVB.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#include <gnrinf.h>
#include <cbiwlk.h>
#include <priunit.h>
#include <symmet.h>
#include <codata.h>
      DIMENSION GVEC(NCORD), EVEC(NCORD), ATCHRG(NCORD)
      DIMENSION ATMARR(MXCENT,8)
      CHARACTER*12 FILENM
#include <ibtfun.h>
C
C     This factor is a humble attempt to give the vectors a nice length
C
!     FAC = 0.5D0*SQRT(1.0D0*NUCDEP)
      FAC = 0.25D0*NUCDEP
C
C     Initialize the ATMARR array. The first index runs over all
C     atoms, the second marks the following properties:
C
C               1 - Element number
C               2 - X coordinate of atom
C               3 - Y coordinate of atom
C               4 - Z coordinate of atom
C               5 - Covalent radius
C               6 - X-component of normal mode vector
C               7 - Y-component of normal mode vector
C               8 - Z-component of normal mode vector
C
      DO 10 I = 1, IATOM
         ATMARR(I,1) = ATCHRG(I)
         DO 12 J = 1, 3
            ATMARR(I,J+1) = XTANG*GVEC((I-1)*3+J)
            ATMARR(I,J+5) = FAC*EVEC((I-1)*3+J)
 12      CONTINUE
         ATMARR(I,5) = RADIUS(NINT(ATMARR(I,1)))
 10   CONTINUE
C
C     We proceed to open the output-file
C
      LUVRML = -1
      FILENM = 'norm_XXX.wrl'
      WRITE(FILENM(6:8),'(I3)') MODENR
      IF (MODENR .LT. 100) WRITE(FILENM(6:6),'(A1)') '0'
      IF (MODENR .LT. 10) WRITE(FILENM(7:7),'(A1)') '0'
      CALL VRINI(LUVRML,FILENM)
C
C     We create all the atoms
C
      INDX = 1
      CALL DRWATM(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
C
C     Then we draw bonds between the atoms if this is requested.
C
      IF (VRBOND) CALL DRWBND(LUVRML,INDX,IATOM,.TRUE.,ATMARR)
C
C     Finally all the vectors are drawn
C
      IF ((MODENR .EQ. 1) .AND. (IPRINT .GT. 0)) THEN
            CALL HEADER('VRML Visualization of Normal Modes',-1)
            WRITE(LUPRI,'(A)')
     &           ' Nrm.Mode     Filename      Frequency (cm-1)'
            WRITE(LUPRI,'(A)')
     &           '--------------------------------------------'
      END IF
      IF (IPRINT .GT. 0) THEN
         WRITE(LUPRI,'(I5,A,A12,A,F16.6)')
     &        MODENR,'       ',FILENM,'       ',FREQ
      END IF
C
C     We draw all the white vectors
C      
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,1)
C
C     ... then we turn all the vectors before we draw the black ones.
C
      DO 20 I = 1, IATOM
         ATMARR(I,6) = -ATMARR(I,6)
         ATMARR(I,7) = -ATMARR(I,7)
         ATMARR(I,8) = -ATMARR(I,8)
 20   CONTINUE
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,0)
      CALL VREND(LUVRML)
      RETURN
      END

C  /* Deck mkvrsy */
      SUBROUTINE MKVRSY(ATM,DRTAXS,MAXAXS,DMRPLN,MAXMIR,WORK,LWORK)
C     
C     Make VRML representation of symmetry elements.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
      DIMENSION ATM(6,0:MXCENT)
      DIMENSION DRTAXS(5,0:MAXAXS), DMRPLN(5,0:MAXMIR)
      DIMENSION WORK(LWORK)
      KATMAR = 1
      KLAST  = KATMAR + 8*MXCENT
      IF (KLAST .GT. LWORK) CALL STOPIT('MKVRSY',' ',KLAST,LWORK)
      CALL MKVRSY_1(ATM,DRTAXS,MAXAXS,DMRPLN,MAXMIR,WORK(KATMAR))
      RETURN
      END

C  /* Deck MKVRSY_1 */
      SUBROUTINE MKVRSY_1(ATM,DRTAXS,MAXAXS,DMRPLN,MAXMIR,ATMARR)
C     
C     See MKVRSY.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#include <gnrinf.h>
#include <cbiwlk.h>
#include <priunit.h>
#include <symmet.h>
#include <codata.h>
      DIMENSION ATM(6,0:MXCENT)
      DIMENSION DRTAXS(5,0:MAXAXS), DMRPLN(5,0:MAXMIR)
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION VEC(3)
      CHARACTER*12 FILENM
      LOGICAL TURN
#include <ibtfun.h>
C
C     We copy the contents of ATM to ATMARR
C
      CALL DZERO(ATMARR,8*MXCENT)
      IATOM = NINT(ATM(1,0))
      NAXS  = NINT(DRTAXS(1,0))
      NPLN  = NINT(DMRPLN(1,0))
      DO 10 I = 1, IATOM
         DO 15 J = 1, 3
            ATMARR(I,J+1) = ATM(J,I)*XTANG
 15      CONTINUE
         ATMARR(I,1) = ATM(4,I)
         ATMARR(I,5) = RADIUS(NINT(ATMARR(I,1)))
 10   CONTINUE
C
C     We find the largest coordinate
C
      CRDMX = 0.0D0
      DO 17 IAT = 1, IATOM
         IF (ABS(ATMARR(IAT,2)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,2))
         IF (ABS(ATMARR(IAT,3)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,3))
         IF (ABS(ATMARR(IAT,4)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,4))
 17   CONTINUE
      CRDMX = CRDMX + MAX(0.20D0,MIN(1.0D0, 0.25D0*CRDMX))
C
C     We proceed to open the output-file
C
      LUVRML = -1
      FILENM = 'firstsym.wrl'
      CALL VRINI(LUVRML,FILENM)
C
C     We create all the atoms
C
      INDX = 1
      CALL DRWATM(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
C
C     Then we draw bonds between the atoms if this is requested.
C
      IF (VRBOND) CALL DRWBND(LUVRML,INDX,IATOM,.FALSE.,ATMARR)
C
C     We draw all the rotational axes. The colour is determined by
C     the order:
C                     2 - Red
C                     3 - Green
C                     4 - Blue
C                     5 - Orange
C                     6 - Yellow
C                     7 - Violet
C                    >7 - Black
C
      IF (NPLN .GT. 0) THEN
         CALL DZERO(ATMARR,8*MXCENT)
         DO 20 II = 1, NPLN
C
C     The normalvector is scaled to reflect the size it should have.
C            
            DO 25 I = 1, 3
               ATMARR(II,I+5) = CRDMX*DMRPLN(I,II)
 25         CONTINUE
 20      CONTINUE
         CALL DRWPLN(LUVRML,INDX,NPLN,ATMARR,.TRUE.)
      END IF
      IF (NAXS .GT. 0) THEN
         CRDMX = CRDMX + 0.25D0
         DO 30 IORD = NINT(DRTAXS(4,1)),2,-1
            CALL DZERO(ATMARR,8*MXCENT)
            NVEC = 0
            DO 32 II = 1, NAXS
               IF (DRTAXS(4,II) .EQ. IORD) THEN
                  NVEC = NVEC + 1
                  DO 34 I = 1, 3
                     VEC(I) = DRTAXS(I,II)
 34               CONTINUE
                  TURN = .FALSE.
C
C     All vectors are turned appropriately
C
                  IF (VEC(1) .LT. 0.0D0) THEN
                     TURN = .TRUE.
                  ELSE IF (ABS(VEC(1)) .LT. 1.0D-10) THEN
                     IF (VEC(2) .LT. 0.0D0) THEN
                        TURN = .TRUE.
                     ELSE IF (ABS(VEC(2)) .LT. 1.0D-10) THEN
                        IF (VEC(3) .LT. 0.0D0) TURN = .TRUE.
                     END IF
                  END IF
                  IF (TURN) THEN
                     VEC(1) = -VEC(1)
                     VEC(2) = -VEC(2)
                     VEC(3) = -VEC(3)
                  END IF
                  DO 35 I = 1, 3
                     ATMARR(NVEC,I+1) = -CRDMX*VEC(I)
                     ATMARR(NVEC,I+5) = 2.0D0*CRDMX*VEC(I)
 35               CONTINUE
               END IF
 32         CONTINUE
            CALL DRWVEC(LUVRML,INDX,NVEC,ATMARR,IORD)
 30      CONTINUE
      END IF
      CALL VREND(LUVRML)
      RETURN
      END

C  /* Deck vrini */
      SUBROUTINE VRINI(LUVRML,FILENM)
C     
C     Open unit LUVRML for output
C
#include <implicit.h>
#include <dummy.h>
      CHARACTER*12 FILENM
      CHARACTER*8 CHDATE, CHTIME
      CALL GETDAT(CHDATE,CHTIME)
      CALL GPOPEN(LUVRML,FILENM,'UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE(LUVRML,'(A)') '#VRML V1.0 ascii'
      WRITE(LUVRML,'(A)') 'Separator {'
      WRITE(LUVRML,'(A)') 'Info {'
      WRITE(LUVRML,'(A)') 'string "Created by DALTON'
      WRITE(LUVRML,'(A,A8,A8)') 'date: ',CHDATE,CHTIME
      WRITE(LUVRML,'(A)') '" }'
      WRITE(LUVRML,'(A)') 'Separator {'
      RETURN
      END

C  /* Deck vrend */
      SUBROUTINE VREND(LUVRML)
C     
C     End of output, close unit LUVRML
C
#include <implicit.h>
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') '}'
      CALL GPCLOSE(LUVRML,'KEEP')
      RETURN
      END

C  /* Deck drwatm */
      SUBROUTINE DRWATM(LUVRML,INDX,IATOM,VIBRA,ATMARR)
C     
C     Draws all the atoms as spheres with appropriate colour
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL VIBRA
      DIMENSION ATMARR(MXCENT,8), RGBCOL(3)      
      DIMENSION COLOUR(60)
#include <ibtfun.h>
C
C     Array with colours for the first 20 elements
C
      DATA (COLOUR(I), I = 1, 60)/
     & 0.750D0, 0.750D0, 0.750D0, 0.950D0, 0.950D0, 0.600D0,
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.400D0,
     & 0.400D0, 0.400D0, 0.300D0, 0.300D0, 0.300D0, 0.150D0, 0.150D0,
     & 1.000D0, 1.000D0, 0.100D0, 0.100D0, 0.100D0, 0.950D0, 0.100D0,
     & 0.700D0, 0.100D0, 0.100D0,
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0,
     & 0.650D0, 0.650D0, 0.400D0, 0.400D0, 0.400D0, 0.760D0, 0.600D0,
     & 0.000D0, 0.960D0, 0.800D0, 0.200D0, 0.800D0, 0.960D0, 0.600D0,
     & 0.950D0, 0.600D0, 0.950D0,
     & 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0, 0.650D0/
C
C     We create all the atoms as spheres with appropriate colour
C
      LSTCHG = 0
      DO 10 INUC = 1, IATOM
         ICHARG = NINT(ATMARR(INUC,1))
         IF (ICHARG .NE. LSTCHG) THEN
            IF (INDX .LT. 10) THEN
               WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
            ELSE
               WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
            END IF
            IF (ICHARG .LE. 20) THEN
               RGBCOL(1) = COLOUR(3*(ICHARG-1)+1)
               RGBCOL(2) = COLOUR(3*(ICHARG-1)+2)
               RGBCOL(3) = COLOUR(3*(ICHARG-1)+3)
            ELSE
               RGBCOL(1) = 0.500D0
               RGBCOL(2) = 0.500D0
               RGBCOL(3) = 0.500D0
            END IF
            WRITE(LUVRML,'(A,3F7.3)') ' ambientColor',
     &           RGBCOL(1), RGBCOL(2), RGBCOL(3)
            WRITE(LUVRML,'(A,3F7.3)') ' diffuseColor',
     &           RGBCOL(1), RGBCOL(2), RGBCOL(3)
            WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
            WRITE(LUVRML,'(A)') ' shininess 0.750 }'
            INDX= INDX + 1
            LSTCHG = ICHARG
         END IF
         WRITE(LUVRML,'(A)') 'Separator {'
         WRITE(LUVRML,'(A,3F7.3,A)')
     &        ' Translation { translation ', ATMARR(INUC,2),
     &        ATMARR(INUC,3), ATMARR(INUC,4), ' }'
         RAD = 0.5D0*ATMARR(INUC,5)
         IF (VIBRA) RAD = 0.15D0
         WRITE(LUVRML,'(A,F5.3,A)') ' Sphere { radius ',RAD,' } }'
 10   CONTINUE
      RETURN
      END

C  /* Deck drwbnd */
      SUBROUTINE DRWBND(LUVRML,INDX,IATOM,VIBRA,ATMARR)
C
C     Draws bonds between the atoms
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <pi.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
      LOGICAL VIBRA
      DIMENSION ATMARR(MXCENT,8), VEC1(3), VEC2(3), VEC3(3)
#include <ibtfun.h>
C
      RADIUS = 0.06D0
      IF (VIBRA) RADIUS = 0.015D0
C
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') 'Separator {'
      IF (INDX .LT. 10) THEN
         WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
      ELSE
         WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
      END IF
      WRITE(LUVRML,'(A)') ' ambientColor 0.600 0.600 0.600'
      WRITE(LUVRML,'(A)') ' diffuseColor 0.600 0.600 0.600'
      WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
      WRITE(LUVRML,'(A)') ' shininess 0.750 }'
      INDX = INDX + 1
      DO 10 I = 1, IATOM - 1
         DO 20 J = I + 1, IATOM
            RADI = ATMARR(I,5)
            RADJ = ATMARR(J,5)            
C     We find the bond vector ...
            VEC1(1) = ATMARR(I,2)-ATMARR(J,2)
            VEC1(2) = ATMARR(I,3)-ATMARR(J,3)
            VEC1(3) = ATMARR(I,4)-ATMARR(J,4)
            DIST = SQRT(DDOT(3,VEC1,1,VEC1,1))
            IF (DIST .LE. 1.2D0*(RADI+RADJ)) THEN
C     ... the center of the bond ...
               VEC2(1) = 0.5D0*(ATMARR(I,2)+ATMARR(J,2))
               VEC2(2) = 0.5D0*(ATMARR(I,3)+ATMARR(J,3))
               VEC2(3) = 0.5D0*(ATMARR(I,4)+ATMARR(J,4))
C     ... the angle to rotate ...
               VEC3(1) = 0.0D0
               VEC3(2) = 1.0D0
               VEC3(3) = 0.0D0
               ANG = VECANG(VEC3,VEC1)
C     ... and finally a vector to rotate the bond around ... 
               VEC3(1) = VEC1(3)
               VEC3(2) = 0.0D0
               VEC3(3) = -VEC1(1)
C
C     This procedure causes problems if VEC2 is parallel to VEC3.
C     The solution is to use another axis. to rotate around
C
               IF (DDOT(3,VEC3,1,VEC3,1) .LT. 1.0D-16) THEN
                  VEC3(1) = 1.0D0
                  VEC3(2) = 0.0D0
                  VEC3(3) = 0.0D0
                  ANG = 0.0D0
                  IF (VEC2(2) .LT. 0.0D0) ANG = PI
               END IF
C     For vibrational visualization, all radii are equal
               H = 0.25D0*(RADJ-RADI)
               IF (VIBRA) THEN
                  RADI = 0.300D0
                  RADJ = 0.300D0
                  H = 0.0D0
               END IF
C     We have to move the center of the bond in accordance with the radii
               VEC2(1) = VEC2(1) + H*(VEC1(1)/DIST)
               VEC2(2) = VEC2(2) + H*(VEC1(2)/DIST)
               VEC2(3) = VEC2(3) + H*(VEC1(3)/DIST)
C     We also have to shorten the bond, so that it only touches the spheres
               H = 1.05D0*(DIST - 0.5D0*(RADI + RADJ))
               IF (VIBRA) H = 0.97D0*H
C
               WRITE(LUVRML,'(A)') 'Separator {'
               WRITE(LUVRML,'(A,3F7.3,A)')
     &              ' Translation { translation ', 
     &              VEC2(1), VEC2(2), VEC2(3), ' }'
               WRITE(LUVRML,'(A,4F7.3,A)')
     &              ' Rotation { rotation ', 
     &              VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
               WRITE(LUVRML,'(A,F5.3,A,F5.3,A)')
     &              ' Cylinder { parts SIDES radius ',RADIUS,' height ',
     &              H,' } }'
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END

C  /* Deck drwaxs */
      SUBROUTINE DRWAXS(LUVRML,INDX,IATOM,ATMARR,TMPMAT)
C
C     Draws coordinate axes
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
      DIMENSION ATMARR(MXCENT,8), TMPMAT(MXCENT,8)
#include <ibtfun.h>
      CALL DZERO(TMPMAT,8*MXCENT)
C
C     We have to find the highest absolute value of any
C     Cartesian coordinate, then we make our coordinate vectors
C     slightly longer than this.
C
      CRDMX = 0.0D0
      DO 10 IAT = 1, IATOM
         IF (ABS(ATMARR(IAT,2)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,2))
         IF (ABS(ATMARR(IAT,3)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,3))
         IF (ABS(ATMARR(IAT,4)) .GT. CRDMX) CRDMX = ABS(ATMARR(IAT,4))
 10   CONTINUE
      CRDMX = CRDMX + MAX(0.35D0,MIN(1.0D0, 0.25D0*CRDMX))
C
C     x-axis
C
      TMPMAT(1,2) = -CRDMX
      TMPMAT(1,3) = 0.0D0
      TMPMAT(1,4) = 0.0D0
      TMPMAT(1,6) = 2.0D0*CRDMX
      TMPMAT(1,7) = 0.0D0
      TMPMAT(1,8) = 0.0D0
      CALL DRWVEC(LUVRML,INDX,1,TMPMAT,2)
C
C     y-axis
C
      TMPMAT(1,2) = 0.0D0
      TMPMAT(1,3) = -CRDMX
      TMPMAT(1,4) = 0.0D0
      TMPMAT(1,6) = 0.0D0
      TMPMAT(1,7) = 2.0D0*CRDMX
      TMPMAT(1,8) = 0.0D0
      CALL DRWVEC(LUVRML,INDX,1,TMPMAT,3)
C
C     z-axis
C
      TMPMAT(1,2) = 0.0D0
      TMPMAT(1,3) = 0.0D0
      TMPMAT(1,4) = -CRDMX
      TMPMAT(1,6) = 0.0D0
      TMPMAT(1,7) = 0.0D0
      TMPMAT(1,8) = 2.0D0*CRDMX
      CALL DRWVEC(LUVRML,INDX,1,TMPMAT,4)
      RETURN
      END

C  /* Deck drweig */
      SUBROUTINE DRWEIG(LUVRML,INDX,ATMARR,EVEC1,EVEC2,IREP,NCR)
C
C     Draws vectors to illustrate vibrational modes
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <priunit.h>
#include <symmet.h>
#include <trkoor.h>
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION EVEC1(MXCOOR), EVEC2(MXCOOR)
      DIMENSION VEC1(3), VEC2(3), VEC3(3)
#include <ibtfun.h>
C
      CALL DZERO(EVEC2,MXCOOR)
      CALL DAXPY(NCR,1.0D0,EVEC1,1,EVEC2,1)
      CALL DZERO(EVEC1,MXCOOR)
      DO 10 I = 1, NCR
         EVEC1(I) = EVEC2(I)*1.0D0
 10   CONTINUE
C
C     The eigenvector in symmetry basis is transformed to
C     cartesian coordinates.
C
      CALL DZERO(EVEC2,MXCOOR)
      DO 20 IAT = 1, NUCIND
         DO 30 ICO = 1, 3
            ICCOOR = 3*(IAT - 1) + ICO
            ISCOOR = IPTCNT(ICCOOR,IREP,1)
            IF (ISCOOR .GT. 0) THEN
               EVEC2(ICCOOR)=EVEC1(ICCOOR)/SQRT(FMULT(ISTBNU(IAT)))
             END IF
 30       CONTINUE
 20   CONTINUE
C
C     We add necessary information to the ATMARR array:
C
C               6 - X component of eigenvector
C               7 - Y component of eigenvector
C               8 - Z component of eigenvector
C
      IATOM = 1
      DO 40 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 45 ISYMOP = 0, MAXOPR
C
C     This factor is a humble attempt to give the vectors a nice length
C
            FAC = 0.5D0*SQRT(1.0D0*NUCDEP)
C
C     Vectors on the symmetry dependent centres, should have a direction
C     according to the symmetry.
C
C     *****************************************************************
C     NOTE!!!!! This test is probably _NOT_ correct!
C     *****************************************************************
C
            IF (IBTAND(ISYMOP,IREP) .GT. 0) FAC = -FAC
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               ATMARR(IATOM,6) =
     &      FAC*PT(IBTAND(ISYMAX(1,1),ISYMOP))*EVEC2((ICENT-1)*3+1)
               ATMARR(IATOM,7) =
     &      FAC*PT(IBTAND(ISYMAX(2,1),ISYMOP))*EVEC2((ICENT-1)*3+2)
               ATMARR(IATOM,8) =
     &      FAC*PT(IBTAND(ISYMAX(3,1),ISYMOP))*EVEC2((ICENT-1)*3+3)
               IATOM = IATOM + 1
            END IF
 45      CONTINUE
 40   CONTINUE
      IATOM = IATOM - 1
C
C     We draw the white vectors...
C
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,1)
C
C     ... then we turn all the vectors before we draw the black ones.
C
      DO 50 I = 1, IATOM
         ATMARR(I,6) = -ATMARR(I,6)
         ATMARR(I,7) = -ATMARR(I,7)
         ATMARR(I,8) = -ATMARR(I,8)
 50   CONTINUE
      CALL DRWVEC(LUVRML,INDX,IATOM,ATMARR,0)
      RETURN
      END

C  /* Deck drwvec */
      SUBROUTINE DRWVEC(LUVRML,INDX,IATOM,ATMARR,ICOLR)
C
C     Draws vectors
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <priunit.h>
#include <pi.h>
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION VEC1(3), VEC2(3), VEC3(3)
C
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') 'Separator {'
      IF (INDX .LT. 10) THEN
         WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
      ELSE
         WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
      END IF
C
C     The colour of the vector is chosen based on the variable ICOLR:
C          0 - Black
C          1 - White
C          2 - Red
C          3 - Green
C          4 - Blue
C          5 - Orange
C          6 - Yellow
C          7 - Violet
C
      IF (ICOLR .EQ. 1) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 0.900 0.900 0.900'
         WRITE(LUVRML,'(A)') ' diffuseColor 0.900 0.900 0.900'
      ELSE IF (ICOLR .EQ. 2) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 1.000 0.100 0.100'
         WRITE(LUVRML,'(A)') ' diffuseColor 1.000 0.100 0.100'
      ELSE IF (ICOLR .EQ. 3) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 0.100 1.000 0.100'
         WRITE(LUVRML,'(A)') ' diffuseColor 0.100 1.000 0.100'
      ELSE IF (ICOLR .EQ. 4) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 0.100 0.100 1.000'
         WRITE(LUVRML,'(A)') ' diffuseColor 0.100 0.100 1.000'
      ELSE IF (ICOLR .EQ. 5) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 1.000 0.600 0.100'
         WRITE(LUVRML,'(A)') ' diffuseColor 1.000 0.600 0.100'
      ELSE IF (ICOLR .EQ. 6) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 1.000 1.000 0.100'
         WRITE(LUVRML,'(A)') ' diffuseColor 1.000 1.000 0.100'
      ELSE IF (ICOLR .EQ. 7) THEN
         WRITE(LUVRML,'(A)') ' ambientColor 1.000 0.100 1.000'
         WRITE(LUVRML,'(A)') ' diffuseColor 1.000 0.100 1.000'
C
C     All others (including ICOLR = 0) are black (very dark grey)
C
      ELSE
         WRITE(LUVRML,'(A)') ' ambientColor 0.200 0.200 0.200'
         WRITE(LUVRML,'(A)') ' diffuseColor 0.200 0.200 0.200'
      END IF
      WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
      WRITE(LUVRML,'(A)') ' shininess 0.750 }'
      INDX = INDX + 1
      DO 10 I = 1, IATOM
C
C     Vec1 contains the position of the atom.
C
         VEC1(1) = ATMARR(I,2)
         VEC1(2) = ATMARR(I,3)
         VEC1(3) = ATMARR(I,4)
C
C     Vec2 contains the vector.
C
         VEC2(1) = ATMARR(I,6)
         VEC2(2) = ATMARR(I,7)
         VEC2(3) = ATMARR(I,8)
         VECNRM = SQRT(DDOT(3,VEC2,1,VEC2,1))
C
C     The vectors are only drawn if they have a certain length
C
         IF (VECNRM .GT. 0.1D0) THEN
            VEC3(1) = 0.0D0
            VEC3(2) = 1.0D0
            VEC3(3) = 0.0D0
C
C     We calculate a vector to rotate around, and the amount to rotate
C
            ANG = VECANG(VEC3,VEC2)
            VEC3(1) = VEC2(3)
            VEC3(2) = 0.0D0
            VEC3(3) = -VEC2(1)
C
C     This procedure causes problems if VEC2 is parallel to VEC3.
C     The solution is to use another axis. to rotate around
C
            IF (DDOT(3,VEC3,1,VEC3,1) .LT. 1.0D-6) THEN
               VEC3(1) = 1.0D0
               VEC3(2) = 0.0D0
               VEC3(3) = 0.0D0
               ANG = 0.0D0
               IF (VEC2(2) .LT. 0.0D0) ANG = PI
            END IF
C
C     The center of the vector is placed correctly
C
            VEC1(1) = VEC1(1)+0.5D0*VEC2(1)
            VEC1(2) = VEC1(2)+0.5D0*VEC2(2)
            VEC1(3) = VEC1(3)+0.5D0*VEC2(3)
            WRITE(LUVRML,'(A)') 'Separator {'
            WRITE(LUVRML,'(A,3F7.3,A)')
     &           ' Translation { translation ', 
     &           VEC1(1), VEC1(2), VEC1(3), ' }'
            WRITE(LUVRML,'(A,4F7.3,A)')
     &           ' Rotation { rotation ', 
     &           VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
            WRITE(LUVRML,'(A,F5.3,A,F5.3,A)')
     &           ' Cylinder { parts SIDES radius 0.03 height ',
     &           VECNRM,' } }'
            VEC1(1) = VEC1(1)+0.5D0*VEC2(1)
            VEC1(2) = VEC1(2)+0.5D0*VEC2(2)
            VEC1(3) = VEC1(3)+0.5D0*VEC2(3)
            WRITE(LUVRML,'(A)') 'Separator {'
            WRITE(LUVRML,'(A,3F7.3,A)')
     &           ' Translation { translation ', 
     &           VEC1(1), VEC1(2), VEC1(3), ' }'
            WRITE(LUVRML,'(A,4F7.3,A)')
     &           ' Rotation { rotation ', 
     &           VEC3(1), VEC3(2), VEC3(3), ANG, ' }'
            WRITE(LUVRML,'(A)')
     &        ' Cone { parts ALL bottomRadius 0.060 height 0.100 } }'
         END IF
 10   CONTINUE
      RETURN
      END

C  /* Deck drwpln */
      SUBROUTINE DRWPLN(LUVRML,INDX,IPLN,ATMARR,TRANSP)
C
C     Draws vectors
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <priunit.h>
      DIMENSION ATMARR(MXCENT,8)
      DIMENSION VEC1(3), VEC2(3), VEC3(3)
      LOGICAL TRANSP
C
      WRITE(LUVRML,'(A)') ' }'
      WRITE(LUVRML,'(A)') 'Separator {'
      IF (INDX .LT. 10) THEN
         WRITE(LUVRML,'(A,I1)') 'Material{ #',INDX
      ELSE
         WRITE(LUVRML,'(A,I2)') 'Material{ #',INDX
      END IF
      WRITE(LUVRML,'(A)') ' ambientColor 0.500 0.500 0.500'
      WRITE(LUVRML,'(A)') ' diffuseColor 0.500 0.500 0.500'
      IF (TRANSP) WRITE(LUVRML,'(A)') ' transparency 0.70'
      WRITE(LUVRML,'(A)') ' specularColor 0.800 0.800 0.800'
      WRITE(LUVRML,'(A)') ' shininess 0.750 }'
      INDX = INDX + 1
      DO 10 I = 1, IPLN
C
C     Vec1 contains the center of the plane
C
         VEC1(1) = ATMARR(I,2)
         VEC1(2) = ATMARR(I,3)
         VEC1(3) = ATMARR(I,4)
C
C     Vec2 contains the normalvector.
C     The norm of this vector defines the size of the plane
C     (each side of the square is twize the norm).
C
         VEC2(1) = ATMARR(I,6)
         VEC2(2) = ATMARR(I,7)
         VEC2(3) = ATMARR(I,8)
         VECNRM = SQRT(DDOT(3,VEC2,1,VEC2,1))
C
C     The planes are only drawn if they have a certain size
C
         IF (VECNRM .GT. 0.1D0) THEN
            VEC3(1) = 0.0D0
            VEC3(2) = 0.0D0
            VEC3(3) = 1.0D0
C
C     We calculate a vector to rotate around, and the amount to rotate
C
            ANG1 = VECANG(VEC3,VEC2)
            IF (ABS(VEC2(1)) .GT. 1.0D-8) THEN
               ANG2 = ATAN(VEC2(2)/VEC2(1))
            ELSE
               ANG2 = 0.0D0
            END IF
            VEC3(1) = VEC2(2)
            VEC3(2) = -VEC2(1)
            VEC3(3) = 0.0D0
C
C     The center of the vector is placed correctly
C
            WRITE(LUVRML,'(A)') 'Separator {'
            WRITE(LUVRML,'(A,3F7.3,A)')
     &           ' Translation { translation ', 
     &           VEC1(1), VEC1(2), VEC1(3), ' }'
            WRITE(LUVRML,'(A,4F7.3,A)')
     &           ' Rotation { rotation ', 
     &           VEC3(1), VEC3(2), VEC3(3), ANG1, ' }'
            WRITE(LUVRML,'(A,F7.3,A)')
     &           ' Rotation { rotation   0.000  0.000  1.000', 
     &           ANG2, ' }'
C
C     As VRML has no plane primitive, we use a cube with zero depth.
C
            WRITE(LUVRML,'(A,F5.3,A,F5.3,A)')
     &           ' Cube { width ',2.0D0*VECNRM,' height ',
     &           2.0D0*VECNRM,' depth 0.000 } }'
         END IF
 10   CONTINUE
      RETURN
      END

C  /* Deck atmini */
      SUBROUTINE ATMINI(ATMARR,IATOM,BOHR)

C     Expand all symmetry-dependent atoms.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <molinp.h>
#include <optinf.h>
#include <gnrinf.h>
#include <cbiwlk.h>
#include <priunit.h>
#include <codata.h>
#include <symmet.h>
      DIMENSION ATMARR(MXCENT,8)
      LOGICAL BOHR
#include <ibtfun.h>
C
C     We initialize the ATMARR array. The first index runs over all
C     atoms, the second marks the following properties:
C
C               1 - Element number
C               2 - X coordinate of atom
C               3 - Y coordinate of atom
C               4 - Z coordinate of atom
C               5 - Covalent radius
C
      FAC = XTANG
      IF (BOHR) FAC = 1.0D0
C
      IATOM = 0
      DO 10 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 20 ISYMOP = 0, MAXOPR
            IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               ATMARR(IATOM,1) = IZATOM(ICENT)
               ATMARR(IATOM,2) =
     &              FAC*PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
               ATMARR(IATOM,3) =
     &              FAC*PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
               ATMARR(IATOM,4) =
     &              FAC*PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               RAD = RADIUS( IZATOM(ICENT) )
               IF (RAD .LT. 0.0D0) RAD = 1.0D0
               ATMARR(IATOM,5) = RAD
            END IF
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
