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
C     Original: abamolden  alex a auer 5.2.99
C
C this file contains all output routines thar write parts of the molden.inp 
C file which will be read by MOLDEN.
C infos about MOLDEN can be obtained at :
C        http://www.caos.kun.nl/~schaft/molden/molden.html
C
C
       SUBROUTINE MOHEAD
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <inftap.h>
#include <molde.h>
C
      DONEIT=.FALSE.
      DONEIU=.FALSE.
      DONEIV=.FALSE.
      DONEIW=.FALSE.
      DONEIX=.FALSE.
      IXYZPRINTS = 0 
      CALL GPOPEN(LUMOLDEN,'molden.inp','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND LUMOLDEN
      WRITE(LUMOLDEN,'(A)') '[Molden Format]'
      CALL FLSHFO(LUMOLDEN)
      RETURN
      END
      SUBROUTINE MOGTO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &                 KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI)

#include <implicit.h>
#include <maxorb.h>
#include <molde.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <aosotr.h>
#include <maxaqn.h>
#include <codata.h>
      PARAMETER (D0 = 0.0D0)
C
#include <ccom.h>
#include <cbirea.h>
#include <nuclei.h>
#include <primit.h>


      CHARACTER CHRSEG*10, GTOTMP*1, SPDCAR*1
      LOGICAL SEG, NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          SEG(KBLOCK)
      IF (.NOT. DONEIT) THEN

         WRITE(LUMOLDEN,'(A)') '[GTO]'

         ICENT  = 0
c
cPRT   icento is the offset counter for all centres;
cPRT   icent  counts only symmetry-distinct centres
c
         ICENTO = 1
         IBLOCK = 0
         INON   = 0
         IPROLD = 0
         IPRIM = 0
         IBS = 0
         DO 100 I = 1, NONTYP
            DO 110 N = 1, NONT(I)
               ICENT = ICENT + 1
               NDEG  = NUCDEG(ICENT)
               ILL   = 0
               DO L = 1, NDEG
                  WRITE (LUMOLDEN,'(/I5,A2)') ICENTO+L-1,' 0'
                  KBCH  = IBLOCK
                  IF (L.EQ.1) THEN
                     IPROLD  = IPRIM
                  ELSE
                     IPRIM = IPROLD
                  END IF
                  DO 200 J = 1, IQM(I) ! l quantum number s, p, d, ...
                     GTOTMP = SPDCAR(J-1)
                     DO 200 K = 1, JCO(J,I) ! n of blocks for a given AO
                        KBCH = KBCH + 1
                        NNUC  = NUC(KBCH)
                        NNRC  = NRC(KBCH)
                        IF (NNUC .EQ. 0) GO TO 200
                        ITYP = NHKOFF(J)
                        IPSTRT = IPRIM
                        IPRIM = IPRIM + NNUC
                        ITYP = ITYP + 1  
                        DO 420 INNRC = 1, NNRC
                           
Clf                  DO 400 L = 1, NDEG
Clf                     IF (.NOT. INON .EQ. I .OR. .NOT. ILL .EQ. L )
Clf                     INON = I      
Clf                     ILL  = L
                           WRITE (LUMOLDEN,'(1X,A1,1X,I3,A5)')
     &                          GTOTMP, NNUC, ' 1.00'
                           DO 410 INNUC = 1, NNUC
                              WRITE (LUMOLDEN,'(1X,F15.7,1X,F15.10)')
     &                             PRIEXP(IPSTRT+INNUC),
     &                             CPRIMU(INNUC,INNRC,KBCH)
                              
 410                       CONTINUE
Clf 400              CONTINUE
 420                    CONTINUE
 200                 CONTINUE
                  END DO
                  ICENTO = ICENTO + NDEG
 110           CONTINUE
            IBLOCK = IBLOCK + NBLCK(I)
 100     CONTINUE
         
      ENDIF

      DONEIT=.TRUE.
      CALL FLSHFO(LUMOLDEN)
      RETURN
      END

      SUBROUTINE MOATOMS(WORD)
C          
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <codata.h>

      CHARACTER*6 NAME
      CHARACTER*4 WORD

C     include character ASYMB(103)*2 with atomic symbols
#include <asymb.h>

#include <molde.h>
#include <nuclei.h>
#include <symmet.h>
#include <pgroup.h>
#include <cbirea.h>
#include <inftap.h>
#include <chrxyz.h>
#include <chrsgn.h>
#include <chrnos.h>
#include <ibtfun.h>

       IF (.NOT. DONEIU .OR. WORD .EQ. 'XYZ' .OR. WORD .EQ. 'FREQ') THEN

          NCOOR = 3*NUCDEP
          IF (WORD .EQ. 'XYZ ') WRITE(LUMOLDEN,'(I5/)') NUCDEP 
          IF (WORD .EQ. 'ATOM') WRITE(LUMOLDEN,'(A)') '[Atoms] AU' 
          IF (WORD .EQ. 'FREQ') WRITE(LUMOLDEN,'(/A)')'[FR-COORD]' 
       
          IATOM = 0
          DO 100 ICENT = 1, NUCIND
             MULCNT = ISTBNU(ICENT)
             NAME   = '      '
             J = 0
             DO I = 1,4
                IF (NAMN(ICENT)(I:I) .NE. ' ') THEN
                   J = J + 1
                   NAME(J:J) = NAMN(ICENT)(I:I)
                END IF
             END DO
             IF (MULT(MULCNT) .EQ. 1) THEN
                ICHARGE = IZATOM(ICENT)
                IF (WORD .EQ. 'ATOM') THEN
C FIXME ?          hjaaj Oct 2003: should we remove point charges ??
                   IATOM = IATOM + 1
                   WRITE (LUMOLDEN,'(A,1X,I5,1X,I5,3(1X,F20.10))')
     &                    NAME,IATOM,ICHARGE,(CORD(K,ICENT),K=1,3)
                ELSE IF (WORD .EQ. 'FREQ') THEN
                   WRITE (LUMOLDEN,'(A,3(1X,F20.10))')
     &                    ASYMB(ICHARGE),(CORD(K,ICENT),K=1,3)
                ELSE IF (WORD .EQ. 'XYZ ') THEN
                   WRITE (LUMOLDEN,'(A,3(1X,F20.10))')
     &                    NAME,(XTANG*CORD(K,ICENT),K=1,3)
                END IF
             ELSE
                JATOM = 0
                J = J + 1
                NAME(J:J) = '_'
                J = J + 1
                DO 200 ISYMOP = 0, MAXOPR
                   IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
                      JATOM = JATOM + 1
                      NAME(J:J) = CHRNOS(JATOM)
                      CRX = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                      CRY = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                      CRZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                      ICHARGE = IZATOM(ICENT)
                      IF (WORD .EQ. 'ATOM') THEN
                         IATOM = IATOM + 1
                         WRITE (LUMOLDEN,'(A,1X,I5,1X,I5,3(1X,F20.10))')
     &                        NAME,IATOM,ICHARGE,CRX,CRY,CRZ
                      ELSE IF (WORD .EQ. 'FREQ') THEN
                         WRITE (LUMOLDEN,'(A,3(1X,F20.10))')
     &                        ASYMB(ICHARGE),CRX,CRY,CRZ
                      ELSE IF (WORD .EQ. 'XYZ ') THEN
                         WRITE (LUMOLDEN,'(A,3(1X,F20.10))')
     &                        NAME,XTANG*CRX,XTANG*CRY,XTANG*CRZ
                      END IF
                   END IF
 200            CONTINUE
             END IF        
 100      CONTINUE

       END IF
       DONEIU = .TRUE.
       CALL FLSHFO(LUMOLDEN)
       RETURN
       END

      SUBROUTINE MOMOS(ITASK,CMO,OCCUP,ORBTRA,UCMO,ORBVEC)
C
C CMO    =  contains MO coefficients, symmetry packed
C UCMO   =  contains MO coefficients, symmetry unpacked
C itask  =  1 print everything to file
C           2 save orbital energies in CMO

#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <inftap.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aosotr.h>
#include <molde.h>
#include <cbieri.h>
#include <inforb.h>
#include <nuclei.h>
#include <symmet.h>
#include <pgroup.h>
    
      LOGICAL WRTELEM,WRTZERO
      DIMENSION CMO(*), OCCUP(*), ORBTRA(*), UCMO(NBAST,NORBT)
      DIMENSION ORBVEC(*)

#include <chrsgn.h>

      MOPRIN = 0
      IF (ITASK .EQ. 1) THEN
         IF (.NOT. DONEIV) THEN 
            WRITE(LUMOLDEN,'(/A)') '[5D]' 
            WRITE(LUMOLDEN,'(/A)') '[MO]' 
C           call traorb to generate SO to AO transformation
            CALL TRAORB(ORBTRA,NBAST,-1,MOPRIN)
            CALL DZERO(UCMO,  NBAST * NORBT)
            CALL DZERO(ORBVEC,NBAST)
            CALL UPKCMO(CMO,UCMO)
            JOCC = 0
            jORB = 0
            DO ISYM = 1,NSYM
               NOCCI = NOCC(ISYM)
               DO I=1,NORB(ISYM)
                  jORB = jORB + 1
                  II = IORB(ISYM) + I         
                  WRITE(LUMOLDEN,'(2A)')     'Sym= ',REP(ISYM-1)
                  WRITE(LUMOLDEN,'(A,F9.4)') 'Ene= ',OREN(II)
                  WRITE(LUMOLDEN,'(A)')      'Spin= Alpha'
                  IF (I .LE. NOCCI) THEN
                     JOCC = JOCC + 1
                     WRITE(LUMOLDEN,'(A,F6.4)') 'Occup= ',OCCUP(JOCC)
                  ELSE
                     WRITE(LUMOLDEN,'(A     )') 'Occup= 0.0000'
                  END IF
                  ICENT = IPCEN(I)
                  CALL DGEMV('T',NBAST,NBAST,1.0D0,ORBTRA,NBAST,
     $                 UCMO(1,jORB),1,0.0D0,ORBVEC,1)
                  DO IB=1,NBAST
                     WRITE(LUMOLDEN,'(I5,1X,2F15.6)') 
     $                    IB,ORBVEC(JAOAO(IB,2))
                  END DO
               END DO
            END DO
            CALL MOATOMS('ATOM') 
         END IF
         DONEIV = .TRUE.
      END IF

      IF (ITASK .EQ. 2) THEN
         DO I = 1, NORBT
            OREN(I)=CMO(I)
         ENDDO
      END IF
      CALL FLSHFO(LUMOLDEN)
      RETURN 
      END

      SUBROUTINE MOFREQ(EVEC,NUMMOD,NCORD,FREQAU)
#include <implicit.h> 
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <molde.h>
#include <inftap.h>

      DIMENSION EVEC(NCORD,NCORD), FREQAU(NCORD)
      SXFAMU = SQRT(XFAMU)
      WRITE(LUMOLDEN,'(A)') '[FREQ]'

      DO IMODE = 1, NUMMOD
         WRITE(LUMOLDEN,'(F10.2)') FREQAU(IMODE)*XTKAYS
      END DO

      CALL MOATOMS('FREQ')       

      WRITE(LUMOLDEN,'(/A)') '[FR-NORM-COORD]'
      DO IMODE = 1, NUMMOD
         WRITE(LUMOLDEN,'(A,I10)') 'Vibration ',IMODE 
         WRITE(LUMOLDEN,'(3(1X,F20.8))') 
     &    (SXFAMU*EVEC(I,IMODE),I=1,NCORD)
      END DO
      CALL FLSHFO(LUMOLDEN)
      END

      SUBROUTINE MOSCFCON(ITER,EMCSCF,WRITENOW)
#include <implicit.h> 
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <molde.h>
#include <inftap.h>
        LOGICAL WRITENOW
        IF(WRITENOW)THEN
           IF (.NOT. DONEIW) THEN 
         
              WRITE(LUMOLDEN,'(A)') '[SCFCONV]'
              WRITE(LUMOLDEN,'(A,I3)') 'scf-first  1  THROUGH ',ITER
              WRITE(LUMOLDEN,'(F20.10)') (OROC(I),I=1,ITER)
              CALL FLSHFO(LUMOLDEN)

           END IF 
           DONEIW = .TRUE.

        ELSE
           OROC(ITER)=EMCSCF
        ENDIF

        RETURN
      END

      SUBROUTINE MOGECON(WRITENOW,EMCSCF)
#include <implicit.h> 
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <molde.h>
#include <inftap.h>
      LOGICAL WRITENOW
      IF(WRITENOW)THEN
         
         WRITE(LUMOLDEN,'(A)') '[GEOCONV]'
         WRITE(LUMOLDEN,'(A)') 'energy'

         DO 100 I=1,IXYZPRINTS
            WRITE(LUMOLDEN,'(F20.10)') EMCEN(I)
 100     CONTINUE
            
      ELSE
         IF (.NOT. DONEIX) THEN  
            WRITE(LUMOLDEN,'(A)') '[GEOMETRIES] XYZ'
            DONEIX = .TRUE.
         END IF

         CALL MOATOMS('XYZ ')
         IXYZPRINTS = IXYZPRINTS + 1
         EMCEN(IXYZPRINTS) = EMCSCF
      ENDIF
      CALL FLSHFO(LUMOLDEN)
      RETURN
      END
