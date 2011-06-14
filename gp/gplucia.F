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

***********************************************************************
*
*========================================================*
*  Common LUCITA and LUCIAREL routines                   *
*  LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig*
*                                                        *
*   collected by Timo Fleig, Jul 12 - Sep, 2002          *
*========================================================*
*
      subroutine abend(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) 'abend called with: ',CHARINP
      call quit (CHARINP)
      return
      end
      subroutine abend1(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) 'abend1 called with: ',CHARINP
      call quit (CHARINP)
      return
      end
      subroutine abend2(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) 'abend2 called with: ',CHARINP
      call quit (CHARINP)
      return
      end
      subroutine abend3(CHARINP)
      implicit none
      character*(*) CHARINP
      write(6,*) 'abend3 called with: ',CHARINP
      call quit (CHARINP)
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE ADADS1_GAS(NK,I1,XI1S,LI1,IORB,NIORB,JORB,NJORB,
     &                KSTR,NKEL,NKSTR,IREO,IZ,
     &                NOCOB,KMAX,KMIN,IEND,SCLFAC)
*
* Obtain I1(KSTR) = +/- A+ IORB A+ JORB !KSTR>
* Only orbital pairs IOB .gt. JOB are included
*
* KSTR is restricted to strings with relative numbers in the
* range KMAX to KMIN
* =====
* Input
* =====
* IORB : First I orbital to be added
* NIORB : Number of orbitals to be added : IORB to IORB-1+NIORB
*        are used. They must all be in the same TS group
* JORB : First J orbital to be added
* LORB : Number of orbitals to be added : JORB to JORB-1+NJORB
*        are used. They must all be in the same TS group
* KMAX : Largest allowed relative number for K strings
* KMIN : Smallest allowed relative number for K strings
*
* ======
* Output
* ======
*
* NK      : Number of K strings
* I1(KSTR,JORB) : ne. 0 =>  a+IORB a+JORB !KSTR> = +/-!ISTR>
* XI1S(KSTR,JORB) : above +/-
*          : eq. 0    a + JORB !KSTR> = 0
* Offset is KMIN
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER KSTR(NKEL,NKSTR)
      INTEGER IREO(*), IZ(NOCOB,*)
*.Output
      INTEGER I1(LI1,*)
      DIMENSION XI1S(LI1,*)
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' ADADS1_GAS speaking '
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' IORB,NIORB ', IORB,NIORB
       WRITE(6,*) ' JORB,NJORB ', JORB,NJORB
      END IF
*
      IORBMIN = IORB
      IORBMAX = IORB + NIORB - 1
*
      JORBMIN = JORB
      JORBMAX = JORB + NJORB - 1
*
      NIJ = NIORB*NJORB
*
      KEND = MIN(NKSTR,KMAX)
      IF(KEND.LT.NKSTR) THEN
        IEND = 0
      ELSE
        IEND = 1
      END IF
      NK = KEND-KMIN+1
*
      DO KKSTR = KMIN,KEND
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Occupation of string ', KKSTR
         CALL IWRTMA(KSTR(1,KKSTR),1,NKEL,1,NKEL)
       END IF
*. Loop over electrons after which JORB can be added
       DO JEL = 0, NKEL
*
         IF(JEL.EQ.0 ) THEN
           JORB1 = JORBMIN - 1
         ELSE
           JORB1 = MAX(JORBMIN-1,KSTR(JEL,KKSTR))
         END IF
         IF(JEL.EQ.NKEL) THEN
           JORB2 = JORBMAX + 1
         ELSE
           JORB2 = MIN(JORBMAX+1,KSTR(JEL+1,KKSTR))
         END IF
         IF(NTEST.GE.1000)
     &    WRITE(6,*) ' JEL JORB1 JORB2 ',JEL,JORB1,JORB2
*
         IF(JEL.GT.0.AND.JORB1.GE.JORBMIN.AND.
     &                   JORB1.LE.JORBMAX) THEN
*. vanishing for any IORB
           IJOFF = (JORB1-JORBMIN)*NIORB
           DO IIORB = 1, NIORB
             IJ = IJOFF + IIORB
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JORB1 IIORB' , JORB1,IIORB
               write(6,*) ' ijoff ', ijoff
               Call Abend('adads1_gas')
             end if
             I1(KKSTR-KMIN+1,IJ) = 0
             XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
           END DO
         END IF
*
         IF(JORB1.LT.JORBMAX.AND.JORB2.GT.JORBMIN) THEN
*. Orbitals JORB1+1 - JORB2-1 can be added after electron JEL
           SIGNJ = (-1) ** JEL * SCLFAC
*. reverse lexical number of the first JEL ELECTRONS
           ILEX0 = 1
           DO JJEL = 1, JEL
             ILEX0 = ILEX0 + IZ(KSTR(JJEL,KKSTR),JJEL)
           END DO
           DO JJORB = JORB1+1, JORB2-1
* And electron JEL + 1
             ILEX1 = ILEX0 + IZ(JJORB,JEL+1)
*. Add electron IORB
             DO IEL = JEL, NKEL
               IF(IEL.EQ.JEL) THEN
                 IORB1 = MAX(JJORB,IORBMIN-1)
               ELSE
                 IORB1 = MAX(IORBMIN-1,KSTR(IEL,KKSTR))
               END IF
               IF(IEL.EQ.NKEL) THEN
                 IORB2 = IORBMAX+1
               ELSE
                 IORB2 = MIN(IORBMAX+1,KSTR(IEL+1,KKSTR))
               END IF
               IF(NTEST.GE.1000)
     &          WRITE(6,*) ' IEL IORB1 IORB2 ',IEL,IORB1,IORB2
               IF(IEL.GT.JEL.AND.IORB1.GE.IORBMIN.AND.
     &                           IORB1.LE.IORBMAX) THEN
                 IJ = (JJORB-JORBMIN)*NIORB+IORB1-IORBMIN+1
             if(ij.gt.nij) then
               write(6,*) ' ij .gt. nij '
               write(6,*) ' JJORB IORB1' , JJORB,IORB1
               write(6,*) ' ijoff ', ijoff
               Call Abend('adads1_gas')
             end if
                 I1(KKSTR-KMIN+1,IJ) = 0
                 XI1S(KKSTR-KMIN+1,IJ) = 0.0D0
               END IF
               IF(IORB1.LT.IORBMAX.AND.IORB2.GT.IORBMIN) THEN
*. Orbitals IORB1+1 - IORB2 -1 can be added after ELECTRON IEL in KSTR
*. Reverse lexical number of the first IEL+1 electrons
                 ILEX2 = ILEX1
                 DO IIEL = JEL+1,IEL
                   ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+1)
                 END DO
*. add terms for the last electrons
                 DO IIEL = IEL+1,NKEL
                   ILEX2 = ILEX2 + IZ(KSTR(IIEL,KKSTR),IIEL+2)
                 END DO
                 IJOFF = (JJORB-JORBMIN)*NIORB
                 SIGNIJ =  SIGNJ*(-1.0D0) ** (IEL+1)
                 DO IIORB = IORB1+1, IORB2-1
                   IJ = IJOFF + IIORB - IORBMIN + 1
                   IF(IJ.LE.0.OR.IJ.GT.NIJ) THEN
                     WRITE(6,*) ' PROBLEMO ADADS1 : IJ : ', IJ
                     WRITE(6,*) ' IJOFF IORBMIN ', IJOFF,IORBMIN
                     WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
                     Call Abend('adads1_gas')
                   END IF
                   ILEX = ILEX2 + IZ(IIORB,IEL+2)
                   IACT = IREO(ILEX)
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' IIORB JJORB ', IIORB,JJORB
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' IJ ILEX,IACT',IJ, ILEX,IACT
                   IF(NTEST.GE.1000)
     &             WRITE(6,*) ' ILEX0 ILEX1 ILEX2 ILEX ',
     &                          ILEX0,ILEX1,ILEX2,ILEX
                   I1(KKSTR-KMIN+1,IJ) = IACT
                   XI1S(KKSTR-KMIN+1,IJ) = SIGNIJ
                   IF(IJ.LT.0) THEN
                     Call Abend2( ' NEGATIVE IJ in ADADS1 ' )
                   END IF
                 END DO
               END IF
             END DO
           END DO
         END IF
       END DO
      END DO
*
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADADST1_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          IJ = 0
          DO  JJORB = JORB,JORB+NJORB-1
            JJORBR = JJORB-JORB+1
            DO  IIORB = IORB, IORB + NIORB - 1
              IJ = IJ + 1
C?            WRITE(6,*) ' IJ = ', IJ
C?            IF(IIORB.GT.JJORB) THEN
                IIORBR = IIORB - IORB + 1
                WRITE(6,*)
     &          ' Info for orbitals (iorb,jorb) ', IIORB,JJORB
                WRITE(6,*) ' Excited strings and sign '
                CALL IWRTMA(I1(1,IJ),1,NK,1,NK)
                CALL WRTMT_LU(XI1S(1,IJ),1,NK,1,NK)
C?            END IF
            END DO
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE ADAST_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,
     *              ISTMAP,SGNMAP,SCLFAC,IAC,LROW_IN,IEC)
*
* Annihilation/Creation mappings from K-strings of given sym in each gasspace
*
* Input
* NSTB : Number of strings before active gasspace
* NSTA : Number of strings after active gasspace
* IOFFK : Offset for K group of strings in active gasspace, i.e. start of
*         this symmetry of active K group strings
* IOFFI : Offset for I group of strings in active gasspace, i.e. start of
*         this symmetry of active I group strings
* IOFFISP: Offset for this symmetrydistribution of active I supergroup strings
* IOFFKSP: Offset for this symmetrydistribution of active K supergroup strings
* ICREORB : Orbital part of creation map for active K groupstrings
* ICRESTR : String  part of creation map for active K groupstrings
* IORBTSF   : First active orbital ( first orbital in in active GASspace
*           with required sym)
* IORBTF   : First orbital in active gas space, (can have any sym)
* NORBTS  : Number of orbitals of given symmetry and type
* NSTAK : Number of K groupstrings with given correct symmetry
* NSTAKT: Total Number of K groupstrings in active group (all symmetries)
* NSTAKTS: Total Number of K supergroup strings with correct symmetry
* ISTAKTS: Offset for K supergroup strings with hiven symmetrydistribution
* NSTAI : Number of I groupstrings in active gasspace
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION ICREORB(LROW_IN,*), ICRESTR(LROW_IN,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
*
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
*
C?    WRITE(6,*) ' ADAST_GASSM: NSTAK, NORBTS  ',NSTAK,NORBTS
C?    WRITE(6,*) ' IAC, IEC , SCLFAC = ', IAC, IEC, SCLFAC
C?    WRITE(6,*) ' LROW_IN = ', LROW_IN
      SIGN0 = (-1)**NELB*SCLFAC
C?    WRITE(6,*) ' SIGN0 = ', SIGN0
      DO KSTR = IOFFK, NSTAK+IOFFK-1
C?      WRITE(6,*) ' KSTR ', KSTR
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
C?        WRITE(6,*) ' IORB = ', IORB
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
C?        write(6,*) ' IORBRTS IORBRT', IORBRTS,IORBRT
*. Change of active group
          I_AM_ACTIVE = 0
          IF(IAC.EQ.2) THEN
C?          WRITE(6,*) ' ICREORB = ', ICREORB(IORBRT,KSTR)
C?          WRITE(6,*) ' ICRESTR = ', ICRESTR(IORBRT,KSTR)
            IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Creation is nonvanishing
              I_AM_ACTIVE = 1
              IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                SIGN = SIGN0
                ISTR = ICRESTR(IORBRT,KSTR)
              ELSE
                SIGN = -SIGN0
                ISTR = -ICRESTR(IORBRT,KSTR)
              END IF
            END IF
          ELSE IF(IAC.EQ.1) THEN
             IF(IEC.EQ.1) THEN
*. Expanded map
               IF(ICREORB(IORBRT,KSTR) .LT. 0 ) THEN
*. Annihilation is non-vanishing
                 I_AM_ACTIVE = 1
                 IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
                   SIGN = SIGN0
                   ISTR = ICRESTR(IORBRT,KSTR)
                 ELSE
                   SIGN = -SIGN0
                   ISTR = -ICRESTR(IORBRT,KSTR)
                 END IF
               END IF
             ELSE
*. Compressed map
               DO IROW = 1, LROW_IN
C?              write(6,*) ' IROW, ICREORB(IROW,KSTR)' ,
C?   &                       IROW, ICREORB(IROW,KSTR)
COLD             IF(ICREORB(IROW,KSTR) .EQ. -IORBRT ) THEN
                 IF(ICREORB(IROW,KSTR) .EQ. -IORB   ) THEN
*. Annihilation is non-vanishing
                   I_AM_ACTIVE = 1
                   IF(ICRESTR(IROW,KSTR) .GT. 0 ) THEN
                     SIGN = SIGN0
                     ISTR = ICRESTR(IROW,KSTR)
                   ELSE
                     SIGN = -SIGN0
                     ISTR = -ICRESTR(IROW,KSTR)
                   END IF
                 END IF
               END DO
             END IF
*            ^ End of expanded/compact switch
           END IF
*          ^ End of Creation/annihilation switch

          IF(I_AM_ACTIVE .EQ. 1  ) THEN
*. Excitation is open, corresponding active I string
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
C?          WRITE(6,*) ' ISTR, relative = ', ISTR
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          WRITE(6,*) ' IADRK0 IOFFK IOFFKSP ',
C?   &                   IADRK0,IOFFK,IOFFKSP
C?          WRITE(6,*) ' ISTR, IADRI0, IOFFISP ',
C?   &                   ISTR, IADRI0, IOFFISP
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
*
C?          WRITE(6,*) ' ISTR NSTA NSTB ',ISTR,NSTA,NSTB
C?          WRITE(6,*) ' NSTAI,NSTAK',NSTAI,NSTAK
            DO IB = 1, NSTB
              DO IA = 1, NSTA
C               IBKA = IADRI0 + (IB-1)*NSTAI*NSTA+IA
C               KBKA = IADRK0 + (IB-1)*NSTAK*NSTA+IA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
          END IF
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADAST_GASSM '
        WRITE(6,*) ' ======================== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Mapped strings and sign '
            CALL IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
            CALL WRTMT_LU(SGNMAP(1,IORBR),1,NK,1,NK)
          END DO
        END IF
      END IF

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE ADDDIA(A,FACTOR,NDIM,IPACK)
*
* add factor to diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      DO 100 I = 1,NDIM
        IF(IPACK .EQ. 0 ) THEN
          II = (I-1)*NDIM + I
        ELSE
          II = I*(I+1)/2
        END IF
        A(II) = A(II) + FACTOR
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE ADD_SKAIIB(SB,NI,NIA,SKAIIB,NKA,NIB,
     &                          I,ISCA,SSCA)
*
* Update Transposed sigma block with contributions for given orbital 
* index j from the matrix S(Ka,i,Ib)
*
* S(Ib,Isca(Ka)) =  S(Ib,Isca(Ka)) + Ssca(Ka)*S(Ka,I,Ib)
*
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION SKAIIB(*),SSCA(*),ISCA(*)
*. Input and Output
       DIMENSION SB(NIB,NIA)
*
C     LBLK = 100
      LBLK = 40
      NBLK = NIB/LBLK
      IF(LBLK*NBLK.LT.NIB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NIB)
        ICONST = NKA*NI
        IADR0 =  (I-1)*NKA+(ICOFF-1-1)*NKA*NI
        IF(ICEND.GT.ICOFF) THEN
*. Use form with Inner loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (I-1)*NKA+(ICOFF-1-1)*NKA*NI
              IADR = IADR0 + KA
              DO IB = ICOFF,ICEND
*. Address of S(Ka,i,Ib)
                IADR = IADR + ICONST
                SB(Ib,IROW) = SB(Ib,IROW)+S*SKAIIB(IADR)
              END DO
            END IF
          END DO
        ELSE
*. Form with no loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA + ICONST
C             DO IB = ICOFF,ICEND
*. Address of S(Ka,i,Ib)
C               IADR = IADR + ICONST
                SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
C             END DO
            END IF
          END DO
        END IF
*       ^ End of test of ICOFF=ICEND
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Save an integer to a real array
*
      subroutine i_save_r(WORK,IVALUE,IELMNT)
*
      implicit none
      integer WORK(*),IVALUE,IELMNT
*
      WORK(IELMNT) = IVALUE
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE ADSTN_GASSM(NSTB,NSTA,IOFFK,IOFFI,IOFFISP,
     &              IOFFKSP,ICREORB,ICRESTR,
     &              IORBTSF,IORBTF,NORBTS,NSTAK,NSTAKT,NSTAI,
     &              NSTAKTS,ISTAKTS,NELB,NACGSOB,
     *              ISTMAP,SGNMAP,SCLFAC)
*
* Creation mappings from K-strings of given sym in each gasspace
*
* Input
* NSTB : Number of strings before active gasspace
* NSTA : Number of strings after accive gasspace
* IOFFK : Offset for K group of strings in active gasspace, i.e. start of
*         this symmetry of active K group strings
* IOFFI : Offset for I group of strings in active gasspace, i.e. start of
*         this symmetry of active I group strings
* IOFFISP: Offset for this symmetrydistribution of active I supergroup strings
* IOFFKSP: Offset for this symmetrydistribution of active K supergroup strings
* ICREORB : Orbital part of creation map for active K groupstrings
* ICRESTR : String  part of creation map for active K groupstrings
* IORBTSF   : First active orbital ( first orbital in in active GASspace
*           with required sym)
* IORBTF   : First orbital in active gas space, (can have any sym)
* NORBTS  : Number of orbitals of given symmetry and type
* NSTAK : Number of K groupstrings with given correct symmetry
* NSTAKT: Total Number of K groupstrings in active group (all symmetries)
* NSTAKTS: Total Number of K supergroup strings with correct symmetry
* ISTAKTS: Offset for K supergroup strings with hiven symmetrydistribution
* NSTAI : Number of I groupstrings in active gasspace
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION ICREORB(NACGSOB,*), ICRESTR(NACGSOB,*)
*. Output
      DIMENSION ISTMAP(NSTAKTS,*),SGNMAP(NSTAKTS,*)
*
C?    WRITE(6,*) ' Creation maps, ORB and STR '
C?    CALL IWRTMA(ICREORB,NACGSOB, NSTAK+IOFFK-1,NACGSOB,
C?   *            NSTAK+IOFFK-1)
C?    CALL IWRTMA(ICRESTR,NACGSOB, NSTAK+IOFFK-1,NACGSOB,
C?   *             NSTAK+IOFFK-1)
*
C?    WRITE(6,*) ' ADSTN_GASSM : NSTA, NSTB, NSTAK',NSTA,NSTB,NSTAK
C?    WRITE(6,*) ' IOFFISP,IOFFKSP', IOFFISP, IOFFKSP
C?    WRITE(6,*) ' IORBTSF IORBTF ', IORBTSF,IORBTF
      IMULTK = NSTAK*NSTB
      IMULTI = NSTAI*NSTB
C?    WRITE(6,*) ' NSTAKT ', NSTAKT
*
      SIGN0 = (-1)**NELB*SCLFAC
C?    WRITE(6,*) ' NELB sign0 = ', NELB, SIGN0
      DO KSTR = IOFFK, NSTAK+IOFFK-1
        DO IORB = IORBTSF, IORBTSF-1+NORBTS
*. Relative to Type-symmetry start
          IORBRTS = IORB-IORBTSF+1
*. Relative to type start
          IORBRT = IORB-IORBTF+1
C?         write(6,*) 'IORB IORBRT KSTR ', IORB,IORBRT, KSTR
C?         WRITE(6,*) 'ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)',
C?   &                 ICRESTR(IORBRT,KSTR),ICREORB(IORBRT,KSTR)
          IF(ICREORB(IORBRT,KSTR) .GT. 0 ) THEN
*. Excitation is open, corresponding active I string
            IF(ICRESTR(IORBRT,KSTR) .GT. 0 ) THEN
              SIGN = SIGN0
              ISTR = ICRESTR(IORBRT,KSTR)
            ELSE
              SIGN = -SIGN0
              ISTR = -ICRESTR(IORBRT,KSTR)
            END IF
* Relative to start of given symmetry for this group
            ISTR = ISTR - IOFFI+ 1
*. This Creation is active for all choices of strings in supergroup
*. before and after the active type. Store the corrsponding mappings
            IADRK0 = (KSTR-IOFFK)*NSTA +IOFFKSP-1
            IADRI0 = (ISTR-1)*NSTA     +IOFFISP-1
C?          write(6,*) ' ISTR IADRK0 IADRI0 = ', ISTR, IADRK0,IADRI0
*
            NSTAINSTA = NSTAI*NSTA
            NSTAKNSTA = NSTAK*NSTA
            DO IB = 1, NSTB
              DO IA = 1, NSTA
                ISTMAP(IADRK0+IA,IORBRTS) = IADRI0 + IA
                SGNMAP(IADRK0+IA,IORBRTS) = SIGN
              END DO
              IADRI0 = IADRI0 +  NSTAINSTA
              IADRK0 = IADRK0 +  NSTAKNSTA
            END DO
          END IF
*
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Output from ADSTN_GASSM '
        WRITE(6,*) ' ======================== '
        NK = NSTB*NSTAK*NSTA
        WRITE(6,*) ' Number of K strings accessed ', NK
        IF(NK.NE.0) THEN
          DO IORB = IORBTSF,IORBTSF + NORBTS  - 1
            IORBR = IORB-IORBTSF+1
            WRITE(6,*) ' Update Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(ISTMAP(1,IORBR),1,NK,1,NK)
            CALL WRTMT_LU(SGNMAP(1,IORBR),1,NK,1,NK)
          END DO
        END IF
      END IF

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE APRBLM2(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)

#include "priunit.h"

      WRITE(lupri,*) ' Blocked matrix '
      WRITE(lupri,*) '================'

      IBASE = 1
      DO 100 IBLK = 1, NBLK
        WRITE(lupri,'(/A,I3)') ' Block ... ',IBLK
        IF(ISYM.EQ.0) THEN
        IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
        CALL WRTMT_LU(A(IBASE),LROW(IBLK),LCOL(IBLK),
     &                  LROW(IBLK),LCOL(IBLK) )
       ELSE
        IF(IBLK .NE. 1 )
     &  IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
        CALL PRSYM(A(IBASE),LROW(IBLK))
       END IF

  100 CONTINUE
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
        SUBROUTINE BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
C
C       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      DOUBLE PRECISION E,F
C*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(NSIZE,1),EL(NSIZE,1)
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO 30 I=1,N
        DO 10 J=1,N
 10     EL(I,J)=0.0D0
 30     EL(I,I)=1.0D0
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO 50 J=1,N1
        DO 45 I=M,N
        IF(A(I,J).EQ.0.0D0)GO TO 45
        D=DSQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
        C=A(J,J)/D
        S=A(I,J)/D
 38     DO 39 K=J,N
        D=C*A(J,K)+S*A(I,K)
        A(I,K)=C*A(I,K)-S*A(J,K)
        A(J,K)=D
 39     CONTINUE
        DO 40 K=1,N
        D=C*EL(J,K)+S*EL(I,K)
        EL(I,K)=C*EL(I,K)-S*EL(J,K)
        EL(J,K)=D
 40     CONTINUE
 45     CONTINUE
 50     M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
 51     DETERP=A(1,1)
        DO 52 I=2,N
 52     DETERP=DETERP*A(I,I)
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
 520    F=A(1,1)
        E=A(1,1)
        DO 58 I=2,N
        IF(DABS(F).LT.DABS(A(I,I)))F=A(I,I)
        IF(DABS(E).GT.DABS(A(I,I)))E=A(I,I)
 58     CONTINUE
        EPSILP=EPSIL
        IF(EPSILP.LE.0)EPSILP=1.0E-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO 100 J1=1,N
C       CALL SLITE(2)
        I=J
        ISL2=1
        DO 90 I1=1,J
C       CALL SLITET(2,K000FX)
        IF(ISL2.EQ.0)K000FX=2
        IF(ISL2.EQ.1)K000FX=1
        IF(ISL2.EQ.1)ISL2=0
        GO TO (70,75),K000FX
 70     A(I,J)=1.0D0/A(I,I)
        GO TO 90
 75     KS=I+1
        D=0.0D0
        DO 80 K=KS,J
 80     D=D+A(I,K)*A(K,J)
        A(I,J)=-D/A(I,I)
 90     I=I-1
 100    J=J-1
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX

C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
 105    M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
 107    D=D+A(I,K)*EL(K,J)
        EL(I,J)=D
 118    CONTINUE
 120    M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
 123    DO 124 I=1,N
        DO 124 J=1,N
 124    A(I,J)=EL(I,J)
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CHECK_CPLX( I_VAR_IN, I_VAR_OUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      IONE = 1
      ITWO = 2
C
CSK      WRITE(6,*) 'I_VAR_IN',I_VAR_IN
CSK      WRITE(6,*) 'I_VAR_OUT',I_VAR_OUT
      IF( I_VAR_IN .eq. IONE ) I_VAR_OUT = IONE
      IF( I_VAR_IN .eq. ITWO ) I_VAR_OUT = ITWO
CSK      WRITE(6,*) 'I_VAR_OUT',I_VAR_OUT
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE CLASS_PROD3(VEC1,VEC2,IBLOCK_OFF,NBLOCK,IBLOCK,
     &                       IBLTOCLS,NCLS,CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions of product to each occupation class
*
* Current version (PROD3 !!) uses IBLTOCLS to give relation
* between blocks and classes
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK(8,*), IBLTOCLS(*)
*. Input/output
      DIMENSION CLSVEC(*)

*
      IOFF = 1
      DO JBLOCK = IBLOCK_OFF,IBLOCK_OFF-1+NBLOCK
      IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        JATP = IBLOCK(1,JBLOCK)
        JBTP = IBLOCK(2,JBLOCK)
        NELMNT = IBLOCK(8,JBLOCK)
        JCLS = IBLTOCLS(JBLOCK)
        XTERM = DDOT(NELMNT,VEC1(IOFF),1,VEC2(IOFF),1)
C?      WRITE(6,*)
C?   &  ' CLASS_PROD : CLASS and XTERM = ', JCLS,XTERM
          CLSVEC(JCLS) = CLSVEC(JCLS) + XTERM
COLD    END IF
        IOFF = IOFF + NELMNT
      END IF
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMT_LU(CLSVEC,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE CLS_TO_BLK(NBLOCK,IBLK_TO_CLS,ICLS_A,IBLK_A)
*
* an array ICLS_A is given for each class as well as
* a block to class array IBLK_TO_CLS.
* Obtain ICLS_A in block form (Well, I am bit tired and my
* pedagogical explanations can be pretty lousy even when I
* am awake !A)
*
* Jeppe Olsen, Jan. 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*.input
      DIMENSION ICLS_A(*),IBLK_TO_CLS(*)
*.output
      INTEGER IBLK_A(*)
*
      DO IBLK = 1, NBLOCK
        IBLK_A(IBLK) = ICLS_A(IBLK_TO_CLS(IBLK))
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CLS_TO BLK : IBLKS_a array '
        CALL IWRTMA(IBLK_A,1,NBLOCK,1,NBLOCK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
*
* Two lists of excitations/annihilations/creations are given.
* Compress to common nonvanishing entries
*
* Jeppe Olsen, November 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION I1(NKIN,N1),XI1(NKIN,N1)
      DIMENSION I2(NKIN,N2),XI2(NKIN,N2)
*
      NKOUT = 0
      DO K = 1, NKIN
        I1ACT  = 0
        DO I = 1, N1
          IF(I1(K,I).NE.0) I1ACT = 1
        END DO
        I2ACT = 0
        DO I = 1, N2
          IF(I2(K,I).NE.0) I2ACT = 1
        END DO
        IF(I1ACT.EQ.1.AND.I2ACT.EQ.1) THEN
          NKOUT = NKOUT + 1
          IF(NKOUT.NE.K) THEN
            DO I = 1, N1
               I1(NKOUT,I) = I1(K,I)
              XI1(NKOUT,I) =XI1(K,I)
            END DO
            DO I = 1, N2
               I2(NKOUT,I) = I2(K,I)
              XI2(NKOUT,I) =XI2(K,I)
            END DO
          END IF
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE COPDIA(A,VEC,NDIM,IPACK)
*
* Copy diagonal of matrix A into vector VEC
*
*   IPACK = 0 : Full matrix
*   IPACK = 1 : Lower triangular matrix
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION A(*),VEC(*)
*
      IF(IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          VEC(I) = A((I-1)*NDIM+I)
  100   CONTINUE
      ELSE
        DO 200 I = 1, NDIM
          VEC(I) = A(I*(I+1)/2)
  200   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
*
* VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
*
* VDSV = SUM(I) VECIN(I) ** 2 /( DIAG(I) + SHIFT )
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VECOUT(1),VECIN(1),DIAG(1)
*
      THRES=1.0D-10
      VDSV = 0.0D0
      DO 100 I=1,NDIM
*
        DIVIDE=DIAG(I)+SHIFT
        IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
*
        VDSV = VDSV + VECIN(I) ** 2 /DIVIDE
        VECOUT(I)=VECIN(I)/DIVIDE
*
  100 CONTINUE
*
      write(6,*) 'DIAVC3 : VECIN, DIAG,VECOUT '
      DO I = 1, NDIM
        WRITE(6,'(3E15.8)') VECIN(I),DIAG(I),VECOUT(I)
      END DO
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,NOBTP,IL,IR)
*
* Obtain types of I,J,K,L so
* <L!a+I a+K a L a J!R> is nonvanishing
* only combinations with type(I) .ge. type(K) and type(J).ge.type(L)
* are included
*
      INTEGER IL(NOBTP),IR(NOBTP)
      INTEGER ITP(*),JTP(*),KTP(*),LTP(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DXTYP_GAS in action '
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' Occupation of left string '
        CALL IWRTMA(IL,1,NOBTP,1,NOBTP)
        WRITE(6,*) ' Occupation of right string '
        CALL IWRTMA(IR,1,NOBTP,1,NOBTP)
      END IF
*
*. Number of differing occupations
      NANNI = 0
      NCREA = 0
      NDIFT = 0
*
      ICREA1 = 0
      ICREA2 = 0
      IANNI1 = 0
      IANNI2 = 0
      DO IOBTP = 1, NOBTP
        NDIFT = NDIFT + ABS(IL(IOBTP)-IR(IOBTP))
        NDIF = IL(IOBTP)-IR(IOBTP)
        IF(NDIF.EQ.2) THEN
*. two electrons of type IOBTP must be created
          ICREA1 = IOBTP
          ICREA2 = IOBTP
          NCREA = NCREA + 2
        ELSE IF (NDIF .EQ. -2 ) THEN
*. Two electrons of type IOBTP must be annihilated
          IANNI1 = IOBTP
          IANNI2 = IOBTP
          NANNI = NANNI + 2
        ELSE IF (NDIF.EQ.1) THEN
*. one electron of type IOBTP must be created
          IF(NCREA.EQ.0) THEN
            ICREA1 = IOBTP
          ELSE
            ICREA2 = IOBTP
          END IF
          NCREA = NCREA + 1
        ELSE IF (NDIF.EQ.-1) THEN
* One electron of type IOBTP must be annihilated
          IF(NANNI.EQ.0) THEN
            IANNI1 = IOBTP
          ELSE
            IANNI2 = IOBTP
          END IF
          NANNI = NANNI + 1
        END IF
      END DO
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)  ' NCREA, NANNI ', NCREA, NANNI
        WRITE(6,*)  ' ICREA2, IANNI2 ', ICREA2,IANNI2
C       WRITE(6,*)  ' ICREA11,IANNI11 ', ICREA11,IANNI11
C       WRITE(6,*)  ' ICREA21,IANNI21 ', ICREA21,IANNI21
      END IF
*
      NDXTP = 0
      IF(NDIFT.GT.4) THEN
        NDXTP = 0
      ELSE
        IF(NCREA.EQ.0.AND.NANNI.EQ.0) THEN
*. strings identical, include diagonal excitions  itp = jtp, ktp=ltp 
          DO IJTP = 1, NOBTP
            IF(IR(IJTP).GE.1) THEN
              DO KLTP = 1, IJTP 
                IF((IJTP.NE.KLTP.AND.IR(KLTP).GE.1).OR.
     &             (IJTP.EQ.KLTP.AND.IR(KLTP).GE.2)) THEN
                   NDXTP = NDXTP + 1
                   ITP(NDXTP) = IJTP
                   JTP(NDXTP) = IJTP
                   KTP(NDXTP) = KLTP
                   LTP(NDXTP) = KLTP
                END IF
              END DO
            END IF
          END DO
*. Strings differ by single excitation
        ELSE IF( NCREA.EQ.1.AND.NANNI.EQ.1) THEN
*. diagonal excitation plus creation in ICREA1 
*                   and annihilation in IANNI1
          DO IDIA = 1, NOBTP
            IF((IDIA.NE.IANNI1.AND.IR(IDIA).GE.1).OR.
     &         (IDIA.EQ.IANNI1.AND.IR(IDIA).GE.2)) THEN
               NDXTP = NDXTP + 1
               ITP(NDXTP) = MAX(ICREA1,IDIA)
               KTP(NDXTP) = MIN(ICREA1,IDIA)
               JTP(NDXTP) = MAX(IANNI1,IDIA)
               LTP(NDXTP) = MIN(IANNI1,IDIA)
            END IF
          END DO
        ELSE IF(NCREA.EQ.2.AND.NANNI.EQ.2) THEN
*. Strings differ by double excitation
          NDXTP = 1
          ITP(1) = ICREA2
          KTP(1) = ICREA1
          JTP(1) = IANNI2
          LTP(1) = IANNI1
        END IF
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,I4)')
     &  ' Number of connecting double excitations ', NDXTP
        IF(NDXTP.NE.0) THEN
          WRITE(6,*) '  ITYP KTYP LTYP JTYP '
          WRITE(6,*) '  ===================='
          DO  IDX = 1,NDXTP
            WRITE(6,'(1X,4I5)')ITP(IDX),KTP(IDX),LTP(IDX),JTP(IDX)
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE EIGEN_LUCI(A,R,N,MV,MFKR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*),R(*)
      DATA TESTIT/1.D-20/
      DATA TESTX/1.D-26/
      DATA TESTY/1.D-18/
C
C        PURPOSE
C           COMPUTE EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C           MATRIX
C
C        USAGE
C           CALL EIGEN_LUCI(A,R,N,MV,MFKR)
C
C        DESCRIPTION OF PARAMETERS
C           A - ORIGINAL MATRIX (SYMMETRIC), DESTROYED IN COMPUTATION.
C               RESULTANT EIGENVALUES ARE DEVELOPED IN DIAGONAL OF
C               MATRIX A IN ASSCENDING ORDER.
C           R - RESULTANT MATRIX OF EIGENVECTORS (STORED COLUMNWISE,
C               IN SAME SEQUENCE AS EIGENVALUES)
C           N - ORDER OF MATRICES A AND R
C           MV- INPUT CODE
C   0   COMPUTE EIGENVALUES AND EIGENVECTORS
C   1   COMPUTE EIGENVALUES ONLY (R NEED NOT BE
C       DIMENSIONED BUT MUST STILL APPEAR IN CALLING
C       SEQUENCE)
C           MFKR=0 NO SORT
C               =1 SORT
C
C        REMARKS
C           ORIGINAL MATRIX A MUST BE REAL SYMMETRIC (STORAGE MODE=1)
C           MATRIX A CANNOT BE IN THE SAME LOCATION AS MATRIX R
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           DIAGONALIZATION METHOD ORIGINATED BY JACOBI AND ADAPTED
C           BY VON NEUMANN FOR LARGE COMPUTERS AS FOUND IN ?MATHEMATICAL
C           METHODS FOR DIGITAL COMPUTERS?, EDITED BY A. RALSTON AND
C           H.S. WILF, JOHN WILEY AND SONS, NEW YORK, 1962, CHAPTER 7
C
C     ..................................................................
C
C
C        ...............................................................
C
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
C        STATEMENT WHICH FOLLOWS.
C
C     DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
C    1 COSX2,SINCS,RANGE
C
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
C        ROUTINE.
C
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTS
C        40, 68, 75, AND 78 MUST BE CHANGED TO DSQRT.  ABS IN STATEMENT
C        62 MUST BE CHANGED TO DABS. THE CONSTANT IN STATEMENT 5 SHOULD
C        BE CHANGED TO 1.0D-12.
C
C        ...............................................................
C
C        GENERATE IDENTITY MATRIX
C
    5 RANGE=1.0D-12
      IF(MV-1) 10,25,10
   10 IQ=-N
      DO 20 J=1,N
      IQ=IQ+N
      DO 20 I=1,N
      IJ=IQ+I
      R(IJ)=0.0D+00
      IF(I-J) 20,15,20
   15 R(IJ)=1.0D+00
   20 CONTINUE
C
C        COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANORMX)
C
   25 ANORM=0.0D+00
      DO 35 I=1,N
      DO 35 J=I,N
      IF(I-J) 30,35,30
   30 IA=I+(J*J-J)/2
      ANORM=ANORM+A(IA)*A(IA)
   35 CONTINUE
      IF(ANORM .LE. 0.0D0) GO TO 165
      ANORM=1.414D+00*DSQRT(ANORM)
      ANRMX=ANORM*RANGE/DFLOAT(N)
C
C        INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR
C
      IND=0
      THR=ANORM
   45 THR=THR/DFLOAT(N)
      IF(THR.LT.TESTY)THR=0.D0
   50 L=1
   55 M=L+1
C
C        COMPUTE SIN AND COS
C
   60 MQ=(M*M-M)/2
      LQ=(L*L-L)/2
      LM=L+MQ
      IF(DABS(A(LM)).LT.TESTY)A(LM)=0.D0
      IF(DABS(A(LM)).EQ.0.D0.AND.THR.EQ.0.D0)GO TO 130
   62 IF( DABS(A(LM))-THR) 130,65,65
   65 IND=1
      LL=L+LQ
      MM=M+MQ
      X=0.5D+00*(A(LL)-A(MM))
      AJUK=(A(LM)*A(LM)+X*X)
      AJUK=DSQRT(AJUK)
      IF(DABS(AJUK).LT.TESTIT)WRITE(6,3000)TESTIT,AJUK,A(LM)
 3000 FORMAT(/,'***DENOMINATOR LT ',D12.6,'. VALUE=',D14.8,
     &'. NUMERATOR=',D14.8)
      Y=0.D0
      IF(DABS(AJUK).LT.TESTIT)GO TO 67
      Y=-A(LM)/AJUK
   67 CONTINUE
   68 CONTINUE
C  68 Y=-A(LM)/ DSQRT(A(LM)*A(LM)+X*X)
      IF(X) 70,75,75
   70 Y=-Y
   75 AJUK=(1.D0-Y*Y)
      IF(AJUK.LT.0.D0)WRITE(6,3001) AJUK
 3001 FORMAT(/,'***DSQRT OF ',D14.8)
      IF(AJUK.LT.0.D0)AJUK=0.D0
      AJUK=DSQRT(AJUK)
      AJUK=2.D0*(1.D0+AJUK)
      AJUK=DSQRT(AJUK)
      SINX=Y/AJUK
   76 CONTINUE
C     SINX=Y/ DSQRT(2.0D+00*(1.0D+00+( DSQRT(1.0D+00-Y*Y))))
      SINX2=SINX*SINX
C  78 COSX= DSQRT(1.0D+00-SINX2)
   78 CONTINUE
      AJUK=1.D0-SINX2
      IF(AJUK.LT.TESTX)AJUK=0.D0
      COSX=DSQRT(AJUK)
      COSX2=COSX*COSX
      SINCS =SINX*COSX
C
C        ROTATE L AND M COLUMNS
C
      ILQ=N*(L-1)
      IMQ=N*(M-1)
      DO 125 I=1,N
      IQ=(I*I-I)/2
      IF(I-L) 80,115,80
   80 IF(I-M) 85,115,90
   85 IM=I+MQ
      GO TO 95
   90 IM=M+IQ
   95 IF(I-L) 100,105,105
  100 IL=I+LQ
      GO TO 110
  105 IL=L+IQ
  110 X=A(IL)*COSX-A(IM)*SINX
      A(IM)=A(IL)*SINX+A(IM)*COSX
      A(IL)=X
  115 IF(MV-1) 120,125,120
  120 ILR=ILQ+I
      IMR=IMQ+I
      X=R(ILR)*COSX-R(IMR)*SINX
      R(IMR)=R(ILR)*SINX+R(IMR)*COSX
      R(ILR)=X
  125 CONTINUE
      X=2.0D+00*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C        TESTS FOR COMPLETION
C
C        TEST FOR M = LAST COLUMN
C
  130 IF(M-N) 135,140,135
  135 M=M+1
      GO TO 60
C
C        TEST FOR L = SECOND FROM LAST COLUMN
C
  140 IF(L-(N-1)) 145,150,145
  145 L=L+1
      GO TO 55
  150 IF(IND-1) 160,155,160
  155 IND=0
      GO TO 50
C
C        COMPARE THRESHOLD WITH FINAL NORM
C
  160 IF(THR-ANRMX) 165,165,45
C
C        SORT EIGENVALUES AND EIGENVECTORS
C
  165 IQ=-N
      IF(MFKR.EQ.0)GO TO 186
  166 CONTINUE
      DO 185 I=1,N
      IQ=IQ+N
      LL=I+(I*I-I)/2
      JQ=N*(I-2)
      DO 185 J=I,N
      JQ=JQ+N
      MM=J+(J*J-J)/2
      IF(A(MM)-A(LL)) 170,185,185
  170 X=A(LL)
      A(LL)=A(MM)
      A(MM)=X
      IF(MV-1) 175,185,175
  175 DO 180 K=1,N
      ILR=IQ+K
      IMR=JQ+K
      X=R(ILR)
      R(ILR)=R(IMR)
  180 R(IMR)=X
  185 CONTINUE
186   CONTINUE
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*  Dummy routine for normal compilations
*
      subroutine errmsg
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
*  This (stupid) little subroutine merely checks whether an
*  integer is odd or even and passes this info to the calling
*  routine.
*   Timo Fleig
*
      subroutine evenodd(IEVOD,INTTEST)
      implicit real*8 (A-H,O-Z)
* add 1
      INTTESTP = INTTEST + 1
* put to real value
      A=INTTEST
      B=INTTESTP
* divide both values by 2
      FIRST = A/2
      SECOND = B/2
* cut off places right of point
      ICUT1 = FIRST
      ICUT2 = SECOND
* compare results
      if (ICUT1.eq.ICUT2) then
* INTTEST is even
         IEVOD = 2
      else
* INTTEST is odd
         IEVOD = 1
      end if
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      subroutine findlow(IARRAY,LENGTH,MVALUE,MPOS)
*
*  Find lowest integer on IARRAY
*  Positive integers assumed, smaller than 1 million
*
      implicit real*8(A-H,O-Z)
*
      dimension IARRAY(LENGTH)
*
      MVALUE = 1000000
      do I=LENGTH,1,-1
         if (IARRAY(I).lt.MVALUE) THEN
           MVALUE = IARRAY(I)
           MPOS = I
         end if
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      subroutine findhigh(IARRAY,LENGTH,MVALUE)
*
*  Find highest integer on IARRAY
*  Positive integers assumed
*
      implicit real*8(A-H,O-Z)
*
      dimension IARRAY(LENGTH)
*
      MVALUE = 0
      do I=1,LENGTH,1
         if (IARRAY(I).gt.MVALUE) MVALUE = IARRAY(I)
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_CKAJJB(CB,NJ,NJA,CKAJJB,NKA,NJB,
     &                          J,ISCA,SSCA)
*
* Obtain for given orbital index j the gathered matrix
*
* C(Ka,j,Jb) = SSCA(Ka)C(Jb,ISCA(Ka))
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CB(NJB,NJA), SSCA(*),ISCA(*)
*. Output
      DIMENSION CKAJJB(*)
*
C?    WRITE(6,*) ' From GET_CKAJJB'
C     LBLK = 100
      LBLK = 40
      NBLK = NJB/LBLK
      IF(LBLK*NBLK.LT.NJB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NJB)
        ICONST = NKA*NJ 
        IADR0 =  (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
        IF(ICEND.GT.ICOFF) THEN
*. Inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
                CKAJJB(IADR) = S*CB(JB,IROW)
              END DO
            ELSE  
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
C               IADR = KA + (J-1)*NKA+(JB-1)*NKA*NJ
                IADR = IADR + ICONST
                CKAJJB(IADR) = 0.0D0          
              END DO
            END IF
          END DO
        ELSE
*. No inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
              IADR = IADR0 + KA
C             DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
                CKAJJB(IADR) = S*CB(ICOFF,IROW)
C             END DO
            ELSE  
              IADR = IADR0 + KA
C             DO JB = ICOFF,ICEND
C               IADR = KA + (J-1)*NKA+(JB-1)*NKA*NJ
                IADR = IADR + ICONST
                CKAJJB(IADR) = 0.0D0          
C             END DO
            END IF
          END DO
        END IF
*       ^ End of test ICEND,ICOFF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      FUNCTION GIJKLL(IREOTS,IPNTR,ISL,XINT,ISMFTO,IBSO,NACOB,NSMOB,
     &         NOCOBS,I,J,K,L)
*
* Obtain (IJ!KL), Lucas order
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IREOTS(*),IPNTR(NSMOB,NSMOB,NSMOB)
      DIMENSION ISL(NSMOB,NSMOB,NSMOB)
      DIMENSION IBSO(*),NOCOBS(*) ,ISMFTO(*)
      DIMENSION XINT(*)
*
C?    write(6,*) ' Hi from GIJKLL'
      II = IREOTS(I)
C?    write(6,*) ' II ',II
      ISM = ISMFTO(I)
C?    write(6,*) ' ISM ',ISM
      NI = NOCOBS(ISM)
C?    write(6,*) ' NI ',NI
      IREL = II - IBSO(ISM) + 1
C?    write(6,*) ' IREL ',IREL
      JJ = IREOTS(J)
      JSM = ISMFTO(J)
      JREL = JJ - IBSO(JSM) + 1
      NJ = NOCOBS(JSM)
      IJ = (IREL-1)*NJ + JREL
      JI = (JREL-1)*NI + IREL
      NJI = NI * NJ
      IJSM = (ISM-1)*NSMOB + JSM
*
      KK = IREOTS(K)
      KSM = ISMFTO(K)
      KREL = KK - IBSO(KSM) + 1
      NK = NOCOBS(KSM)
      LL = IREOTS(L)
      LSM = ISMFTO(L)
      LREL = LL - IBSO(LSM) + 1
      NL = NOCOBS(LSM)
      LK = (LREL-1)*NK + KREL
      KL = (KREL-1)*NL + LREL
      NLK = NK * NL
      KLSM = (KSM-1)*NSMOB + LSM
C?    WRITE(6,*) ' IJSM KLSM ', IJSM,KLSM
C?    WRITE(6,*) ' ISM JSM KSM LSM ',ISM,JSM,KSM,LSM

      IF(  (IJSM.GE.KLSM.AND.LSM.NE.ISL(ISM,JSM,KSM))
     &.OR. (IJSM.LT.KLSM.AND.JSM.NE.ISL(KSM,LSM,ISM)) )   THEN
        GIJKLL = 0.0D0
        write(6,*) ' Symmetry zero returned '
      ELSE
*
        IF(IJSM.GT.KLSM) THEN
C         IJKLO = (IJ-1)*NKL + KL + IPNTR(ISM,JSM,KSM)-1
          IJKLO = (LK-1)*NJI + JI + IPNTR(ISM,JSM,KSM)-1
        ELSE IF(IJSM.LT.KLSM) THEN
C         IJKLO = (KL-1)*NIJ + IJ + IPNTR(KSM,LSM,ISM)-1
          IJKLO = (JI-1)*NLK + LK + IPNTR(KSM,LSM,ISM)-1
        ELSE IF( IJSM.EQ.KLSM) THEN
C         IF(IJ.GE.KL) THEN
          IF(JI.GE.LK) THEN
C           IJKLO = IJ*(IJ-1)/2+KL + IPNTR(ISM,JSM,KSM)-1
            IJKLO = JI*(JI-1)/2+LK + IPNTR(ISM,JSM,KSM)-1
          ELSE
C           IJKLO = KL*(KL-1)/2+IJ + IPNTR(ISM,JSM,KSM)-1
            IJKLO = LK*(LK-1)/2+JI + IPNTR(ISM,JSM,KSM)-1
          END IF
        END IF
        GIJKLL = XINT(IJKLO)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,5I4,3X,E18.12)')
     &  ' GIJKLL I J K L ,IJKLO,(IJ!KL) ', I,J,K,L,IJKLO,GIJKLL
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE GRAPW_LUCI(W,Y,MINEL,MAXEL,NORB,NEL,NTEST)
*
* A graph of strings has been defined from
*
*      MINEL(I) is the smallest allowed number of electrons in
*      orbitals 1 through I
*
*      MAXEL(I) is the largest allowed number of electrons in
*      orbitals 1 through I
*
* Set up vertex weights W
* Set up arc weights    Y
*
* Reverse lexical ordering is used with
* weights of unoccupied orbitals set to 0
*
* Jeppe Olsen
*
       IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
       INTEGER W(NORB+1,NEL+1)
       INTEGER Y(NORB,NEL)
       INTEGER MAXEL(NORB),MINEL(NORB)
*
C      NTEST = 0
       call izero(W,(NEL+1)*(NORB+1))
       call izero(Y,NEL*NORB)
*
*================
*  Vertex weights
*================
*
*. (Weight for vertex(IEL,IORB) is stored in W(IORB+1,IEL+1) )
      W(1,1) = 1
      DO 300 IEL = 0, NEL
        DO 200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            IF( IEL .GT. 0 ) THEN
              W(IORB+1,IEL+1) = W(IORB-1+1,IEL+1)
     &                        + W(IORB-1+1,IEL-1+1)
            ELSE
              W(IORB+1,1) = W(IORB-1+1,1)
            END IF
          END IF
  200   CONTINUE
  300 CONTINUE
*
*=============
* Arc weights
*=============
*
*. Weight for arc connecting vertices (IORB-1,IEL-1) and(IORB,IEL)
*. is stored in Y(IORB,IEL)
*. Y(IORB,IEL) = W(IORB-1,IEL)
      DO 1300 IEL = 1, NEL
        DO 1200 IORB = 1, NORB
          IF(MINEL(IORB).LE.IEL .AND. IEL .LE. MAXEL(IORB) ) THEN
            Y(IORB,IEL) = W(IORB-1+1,IEL+1)
          END IF
 1200   CONTINUE
 1300 CONTINUE
*
      IF( NTEST .GE.10 ) THEN
C       WRITE(6,'(A)') ' Matrix of vertex weights '
C       WRITE(6,'(A)') ' ========================'
C       CALL IWRTMA(W,NORB+1,NEL+1,NORB+1,NEL+1)
        WRITE(lupri,'(A)') '  Matrix for arc weights  '
        WRITE(lupri,'(A)') '  ======================'
        CALL IWRTMA(Y,NORB,NEL,NORB,NEL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE GSTTBL(C,CTT,IATP,IASM,IBTP,IBSM,IOCOC,
     &                  NOCTPA,NOCTPB,NSASO,NSBSO,PSSIGN,ICOOSC,IDC,
     &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
*. If ISCALE = 1, the routine scales and returns the block
*  in determinant normalization, and SCLFAC = 1.0D0
*
* If ISCALE = 0, the routine does not perform any overall
* scaling, and a scale factor is returned in SCLFAC
*
* IF ISCALE = 0, zero blocks are not set explicitly to zero,
* instead  zero is returned in SCLFAC
*
* ISCALE, SCLFAC added May 97
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION C(*),CTT(*),NSASO(NSMST, *),NSBSO(NSMST, *)
      DIMENSION IOCOC(NOCTPA,NOCTPB),ICOOSC(NOCTPA,NOCTPB,*)
      DIMENSION SCR(*)
*
      NTEST = 000
*
C?    write(6,*) ' GSTTBL  ,IATP,IASM,IBTP,IBSM,IOCOC'
C?    write(6,*)            IATP,IASM,IBTP,IBSM,IOCOC(1,1)
* =================
* Read in from disc
* =================
      IF(LUC.NE.0) THEN
        CALL IFRMDS(LBL,1,-1,LUC)
C?      write(6,*) ' LBL = ', LBL
        IF(ISCALE.EQ.1) THEN
          CALL FRMDSC_LUCI(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
        ELSE
          NO_ZEROING = 1
          CALL FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
        END IF
*
        IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
          SCLFAC = 0.0D0
        ELSE
          NAST = NSASO(IASM,IATP)
          NBST = NSBSO(IBSM,IBTP)
          IF(LBL.NE.0) THEN
            isgvst = 0
C           ... not used, just to silence ftnchek warning /Aug 2004
            CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &           IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
          ELSE
            SCLFAC = 0.0D0
          END IF
        END IF
*
C?      WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
C?   &  ISCALE,SCLFAC

*. ISGVST and PLSIGN missing to make it work for IDC = 3,4
      ELSE
        write(6,*) 'C array not initialized !!'
        Call Abend2( 'Fatal error in GSTTBL.' )
* =================
* Pack out from C
* =================
      IF(ISCALE.EQ.0) THEN
         WRITE(6,*) ' GSTTBL : LUC = 0 and ISCALE = 0'
         WRITE(6,*) ' I will scale as normal '
         SCLFAC = 1.0D0
      END IF
* Permutation sign
      IF(IDC.EQ.2) THEN
        PSIGN = PSSIGN
      ELSE IF(IDC .EQ. 3 ) THEN
        PSIGN = PLSIGN
      END IF
      PLSSGN = PLSIGN * PSSIGN
* check for different packing possibilities and unpack
      IF(IASM.GT.IBSM.OR.IDC.EQ.1
     &   .OR.(IDC.EQ.3.AND.IASM.GE.IBSM))THEN
**************
** IASM > IBSM
**************
        IF ( IDC.LT.4 ) THEN
*. Simple copy
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
          CALL COPVEC(C(IBASE),CTT,NELMNT)
*
          if (NTEST.ge.100) then
            write(6,*) ' simple copy IBASE NELMNT ',IBASE,NELMNT
            CALL WRTMT_LU(CTT,NSASO(IASM,IATP),NSBSO(IBSM,IBTP),
     &      NSASO(IASM,IATP),NSBSO(IBSM,IBTP))
          end if
*
        ELSE IF( IDC.EQ.4 ) THEN
*. MLMS packed
          IF(IATP.GT.IBTP) THEN
            IBASE = ICOOSC(IATP,IBTP,IASM)
            NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
            CALL COPVEC(C(IBASE),CTT,NELMNT)
          ELSE IF(IATP.EQ.IBTP) THEN
            IBASE = ICOOSC(IATP,IATP,IASM)
            NAST = NSASO(IASM,IATP)
            CALL TRIPK3(CTT,C(IBASE),2,NAST,NAST,PLSIGN*PSSIGN)
          ELSE IF( IATP.LT.IBTP) THEN
            IBASE = ICOOSC(IBTP,IATP,IASM)
            NROW  = NSASO(IASM,IBTP)
            NCOL  = NSBSO(IBSM,IATP)
            CALL TRPMT3(C(IBASE),NROW,NCOL,CTT)
            NELMNT = NROW*NCOL
            CALL SCALVE(CTT,PLSIGN*PSSIGN,NELMNT)
          END IF
        END IF
      ELSE IF( IASM.EQ.IBSM) THEN
**************
** IASM = IBSM
**************
        IF(IATP.GT.IBTP.OR.IDC.EQ.3) THEN
*.. simple copying
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
          CALL COPVEC(C(IBASE),CTT,NELMNT)
        ELSE IF( IATP.EQ.IBTP) THEN
*.. expand triangular packed matrix
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NAST = NSASO(IASM,IATP)
          CALL TRIPK3(CTT,C(IBASE),2,NAST,NAST,PSSIGN)
        ELSE IF( IATP .LT. IBTP) THEN
*.. transpose ibtp iasm iatp ibsm block
          IBASE = ICOOSC(IBTP,IATP,IASM)
          NRI = NSASO(IASM,IBTP)
          NCI = NSBSO(IASM,IATP)
          CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
          IF(PSSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
        END IF
      ELSE IF( IASM .LT. IBSM ) THEN
**************
** IASM < IBSM
**************
*.. transpose ibtp ibsm iatp iasm block
        IF(IDC.LT.4) THEN
          IBASE = ICOOSC(IBTP,IATP,IBSM)
          NRI = NSASO(IBSM,IBTP)
          NCI = NSBSO(IASM,IATP)
          IF( IDC.EQ.2) THEN
            CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
          ELSE IF( IDC.EQ.3) THEN
            CALL COPVEC(C(IBASE),CTT,NRI*NCI)
          END IF
          IF(PSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
        ELSE IF ( IDC .EQ. 4 ) THEN
          IF(IBTP.GT.IATP) THEN
            IBASE = ICOOSC(IBTP,IATP,IBSM)
            NRI = NSASO(IBSM,IBTP)
            NCI = NSBSO(IASM,IATP)
            CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
            IF(PSSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
          ELSE IF (IBTP.EQ.IATP) THEN
            IBASE = ICOOSC(IBTP,IATP,IBSM)
            NRI   = NSASO(IBSM,IATP)
            NCI   = NSBSO(IASM,IATP)
            CALL TRIPK3(CTT,C(IBASE),2,NRI,NCI,PLSSGN)
            IF(PLSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
          ELSE IF( IBTP.LT.IATP) THEN
            IBASE = ICOOSC(IATP,IBTP,IBSM)
            NELMNT = NSASO(IBSM,IATP)*NSBSO(IASM,IBTP)
            CALL COPVEC(C(IBASE),CTT,NELMNT)
            IF(PLSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NELMNT)
          END IF
        END IF
      END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      FUNCTION GTH1ES(IREOTS,IPNT,H,IBSO,MXPNGAS,
     &           IBTSOB,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,IJSM)
*
* one electron integral between orbitals (iorb,itp,ism,jorb,jsm,jtp)
* correct combination of row and column symmetry is assumed
* IJSM = 1 => Lower triangular packed
*      else=> No triangular packing
*
* Last Revision January 98 (IJSM added )
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER IREOTS(*),IPNT(*),IBTSOB(MXPNGAS,*),IBSO(*)
      INTEGER NACOBS(*)
      DIMENSION H(*)
*
      IABS = IORB+IBTSOB(ITP,ISM)-1
      IREO = IREOTS(IABS)
      JABS = JORB+IBTSOB(JTP,JSM)-1
      JREO = IREOTS(JABS)
C?    write(6,*) ' GTH1ES : IREO JREO ',IREO,JREO
*
C?    write(6,*) ' GTH1ES : IBSO ', IBSO(ISM)
      IF(IJSM.EQ.1) THEN
        IF(ISM.GT.JSM) THEN
          NI = NACOBS(ISM)
          IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
        ELSE IF(ISM.EQ.JSM) THEN
          IJMAX = MAX(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
          IJMIN = MIN(IREO-IBSO(ISM)+1,JREO-IBSO(JSM)+1)
          IJ = IPNT(ISM)-1+IJMAX*(IJMAX-1)/2+IJMIN
        ELSE IF (ISM.LT.JSM) THEN
          NJ = NACOBS(JSM)
          IJ = IPNT(JSM)-1+(IREO-IBSO(ISM))*NJ+JREO-IBSO(JSM)+1
        END IF
      ELSE
         NI = NACOBS(ISM)
         IJ = IPNT(ISM)-1+(JREO-IBSO(JSM))*NI+IREO-IBSO(ISM)+1
      END IF
*
      GTH1ES = H(IJ)
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' one electron integral '
        WRITE(6,*) ' IORB ITP ISM ',IORB,ITP,ISM
        WRITE(6,*) ' JORB JTP JSM ',JORB,JTP,JSM
        WRITE(6,*) ' IJ and H(IJ) ', IJ,H(IJ)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*                                                                     *
***********************************************************************
      SUBROUTINE IAIBCM_GAS(LCMBSPC,ICMBSPC,
     &                      MNMXOC,NOCTPA,NOCTPB,IOCA,IOCB,NELFTP,
     &                      MXPNGAS,NGAS,IOCOC,IPRNT)
*
* Allowed combinations of alpha and beta types, GAS version
*
*
* =====
*.Input
* =====
*
* LCMBSPC : Number of GAS spaces included in this expnasion
* ICMBSPC : Gas spaces included in this expansion
*
* MXMNOC(IGAS,1,IGASSPC) : 
*      Min accumulated occ for AS 1-IGAS for space IGASSPC
* MXMNOC(IGAS,2,IGASSPC) : 
*      Max accumulated occ for AS 1-IGAS for space IGASSPC
*
* NOCTPA : Number of alpha types
* NOCTPB : Number of beta types
*
* IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
* IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
*
* MXPNGAS : Largest allowed number of gas spaces
* NGAS    : Actual number of gas spaces

*
* ======
*.Output
* ======
*
* IOCOC(IATP,IBTP)  = 1 =>      allowed combination
* IOCOC(IATP,IBTP)  = 0 => not allowed combination
*
*.Input
      INTEGER ICMBSPC(LCMBSPC)
      INTEGER MNMXOC(MXPNGAS,2,*)
C     INTEGER MNOCC(NGAS),MXOCC(NGAS)
      INTEGER IOCA(MXPNGAS,NOCTPA),IOCB(MXPNGAS,NOCTPB)
      INTEGER NELFTP(*)
*.Output
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
#include "parluci.h"
*
      NTEST = 00000
      NTEST = MAX(NTEST,IPRNT)

*
      IF(NTEST.GE.100) THEN
        IF(luci_MYPROC.EQ.luci_MASTER) THEN
          WRITE(luwrt,*) ' IAIBCM_GAS entered '
          WRITE(luwrt,*) ' ==================='
          WRITE(luwrt,*)
          WRITE(luwrt,*) ' Number of GAS spaces included ', LCMBSPC
          WRITE(luwrt,*) ' GAS spaces included ',
     &                     (ICMBSPC(II),II=1,LCMBSPC)
          WRITE(luwrt,*)
          IF(NTEST.GE.400) THEN
            WRITE(luwrt,*) ' IOCA and IOCB '
            CALL IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NOCTPA)
            CALL IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NOCTPB)
          END IF
        END IF
      END IF
*
      CALL IZERO(IOCOC,NOCTPA*NOCTPB)
      DO 100 IATP = 1, NOCTPA
         DO 90 IBTP = 1, NOCTPB
*. is this combination allowed in any of the GAS spaces included
CMI        INCLUDE = 0
           INCLUD = 0
           DO JJCMBSPC = 1, LCMBSPC
             JCMBSPC = ICMBSPC(JJCMBSPC)
             IEL = 0
             IAMOKAY = 1
             DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
               IF(IEL.LT.MNMXOC(IGAS,1,JCMBSPC).OR.
     &            IEL.GT.MNMXOC(IGAS,2,JCMBSPC))
     &         IAMOKAY = 0
             END DO
CMI          IF(IAMOKAY.EQ.1) INCLUDE = 1
             IF(IAMOKAY.EQ.1) INCLUD = 1
           END DO
CMI        IF(INCLUDE .EQ.1) THEN
           IF(INCLUD .EQ.1) THEN
*. Congratulations , you are allowed
              IOCOC(IATP,IBTP) = 1
          END IF
   90   CONTINUE
  100 CONTINUE
*
      IF( NTEST .GE. 300 ) THEN
        WRITE(luwrt,*) ' Matrix giving allowed combinations of types'
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ICOPVE2(IIN,IOFF,NDIM,IOUT)
*
* IOUT(I) = IIN(IOFF-1+I),I = 1, NDIM
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input 
      DIMENSION IIN(*)
*. Output
      DIMENSION IOUT(*)
*
      DO I = 1, NDIM
        IOUT(I) = IIN(IOFF-1+I)
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ICOPMT(MATI,NRI,NCI,MATO,NRO,NCO)
*
* Copy integer matrix MATI to MATO
*
*. Input
      INTEGER  MATI(NRI,NCI)
*. Output
      INTEGER MATO(NRO,NCO) 
*
      NREFF = MIN(NRI,NRO)
      NCEFF = MIN(NCI,NCO)
*
      DO IC = 1, NCEFF
        DO IR = 1, NREFF
          MATO(IR,IC) = MATI(IR,IC)
        END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ICPMT2(AIN,AOUT,NINR,NINC,NOUTR,NOUTC,IZERO)
*
* Copy INTEGER matrix AIN to AOUT . Dimensions can differ
* 
* If IZERO .ne. 0 , AOUT is zeroed  first
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input 
      INTEGER AIN(NINR,NINC)
*. Output
      INTEGER AOUT(NOUTR,NOUTC)
*
      IF(IZERO.NE.0) CALL ISETVC(AOUT,0,NOUTR*NOUTC)
      DO 100 J = 1, NINC
       CALL ICOPVE(AIN(1,J),AOUT(1,J),NINR)
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION IELSUM(IVEC,NELMNT)
*
* Sum elements of integer vector IVEC
*
      DIMENSION IVEC(*)
*
      ISUM = 0
      DO 100 IELMNT = 1, NELMNT
        ISUM = ISUM + IVEC(IELMNT)
  100 CONTINUE
*
      IELSUM = ISUM
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      INTEGER FUNCTION IGIVE_I_B(x_in)
      implicit none
      real(8), intent(in)  :: x_in
*
* return integer value of real*8 element
*
*
      IGIVE_I_B = INT( X_IN )
*
      END
***********************************************************************
      INTEGER FUNCTION IGET_ACTIVE_IRREP(IXSYMLIST,IXSYMACT,ISM)
      
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION IXSYMLIST(IXSYMACT)
C
C     return active irrep # corresponding to symmetry irrep ISM
C
      IGET_ACTIVE_IRREP = 0
      IRETURN           = 0
      DO I = 1, IXSYMACT
        IF( IXSYMLIST(I) .eq. ISM ) IRETURN = I
      END DO
      IGET_ACTIVE_IRREP = IRETURN
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION IFREQ(IVEC,IVAL,NDIM)
*
* Number of times IVAL occurs in IVEC
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION IVEC(*)
*
      NTIME = 0
      DO 100 I = 1, NDIM
        IF(IVEC(I).EQ.IVAL) NTIME = NTIME + 1
  100 CONTINUE
*
      IFREQ = NTIME
*
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INI_BLOCKL(IBLOCKI,NDIM,IBATS)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     initialize vector containing block length
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBLOCKI(NDIM), IBATS(8,*)
C
      CALL ICOPY(NDIM,IBATS(8,1),8,IBLOCKI,1)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INI_BLOCKL_PRP(IBLOCKI,IMAXDIM,IXSYM,NDIM,IBATS)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     initialize vector containing block length - property run
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IBLOCKI(IMAXDIM,*),IBATS(8,*)
C
      CALL ICOPY(NDIM,IBATS(8,1),8,IBLOCKI(1,IXSYM),1)
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      REAL*8 FUNCTION INPRODB(VEC1,VEC2,NBLK,LBLK,I0BLK)
*
* Inner products between blocked vectors with check of
* zero blocks. Zero blocks are flagged by a unit entry
* in I0BLK or a negative block length
*
*Jeppe Olsen, May 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION VEC1(*),VEC2(*), LBLK(NBLK),I0BLK(NBLK)
*
      REAL*8 INPROD
*
      X = 0.0D0
      IOFF = 1
      DO IBLK = 1, NBLK
        LLBLK = LBLK(IBLK)
C?      WRITE(6,*) ' INPRODB IBLK LLBLK IOFF ',IBLK,LLBLK,IOFF
        IF(I0BLK(IBLK).EQ.0.AND.LLBLK.GT.0) THEN
          X = X + DDOT(LLBLK,VEC1(IOFF),1,VEC2(IOFF),1)
C?        WRITE(6,*) ' Vec1 and Vec2 blocks '
C?        CALL WRTMT_LU(VEC1(IOFF),1,LLBLK,1,LLBLK)
C?        CALL WRTMT_LU(VEC2(IOFF),1,LLBLK,1,LLBLK)
C?        write(6,*) ' Updated x', X
        END IF
        IF(LLBLK.GT.0) IOFF = IOFF + LLBLK
      END DO
*
      INPRODB = X
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM)
C FIND INVERSE OF MATRIX A
C INPUT :
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      ITEST=0
      IF(NDIM.EQ.1)THEN
        IF(A(1,1) .NE. 0.0D0 ) THEN
           A(1,1) = 1.0D0/A(1,1)
        ELSE
           ITEST = 1
        END IF
      ELSE
        DETERM=0.0D0
        EPSIL=0.0D0
        CALL BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
      END IF
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (6,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMT_LU(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IVCSUM(IA,IB,IC,IFACB,IFACC,NDIM)
*
* Add two (scaled) integer vectors
*
*        IA(*) = IFACB*IB(*) + IFACC*IC(*)
*
      DIMENSION IA(*),IB(*),IC(*)
*
      DO 100 I = 1, NDIM
        IA(I) = IFACB * IB(I) + IFACC * IC(I)
  100 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      FUNCTION ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
*
* Adress of string IOCC
*
* version of Winter 1990 , Jeppe Olsen
*
      INTEGER Z
      DIMENSION IOCC(*),NEWORD(*),Z(NORB,*)
*
      IZ = 1
      DO 100 I = 1,NEL
        IZ = IZ + Z(IOCC(I),I)
  100 CONTINUE 
*
      IF(IREORD.EQ.0) THEN
        ISTRNM = IZ
      ELSE
        ISTRNM = NEWORD(IZ)
      END IF
*
      NTEST = 0
      IF ( NTEST .GT. 1 ) THEN
        WRITE(6,*) ' STRING'
        CALL IWRTMA(IOCC,1,NEL,1,NEL)
        WRITE(6,*) ' Z matrix ' 
        CALL IWRTMA(Z,NORB,NEL,NORB,NEL)
        WRITE(6,*) ' ADRESS OF STRING ',ISTRNM
        WRITE(6,*) ' REV LEX number : ', IZ
      END IF 
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IWRTMA10(IMAT,NROW,NCOL,MAXROW,MAXCOL)
* I10 format
      DIMENSION IMAT(MAXROW,MAXCOL)
C
      DO 100 I = 1, NROW
        WRITE(6,1110) (IMAT(I,J),J= 1,NCOL)
  100 CONTINUE
 1110 FORMAT(/,1X,8I10,/,(1X,8I10))
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE LFTPOS(CARD,LENGTH)
*
* Eliminate blanks and left position character string CARD
*
      CHARACTER*72 CARD
*
      IEFF = 0
C     WRITE(6,*) ' INPUT string to LFTPOS ' 
C     WRITE(6,'(1X,A)') CARD

      DO 100 IPOS = 1, LENGTH
       IF(CARD(IPOS:IPOS).NE.' ') THEN
         IEFF = IEFF + 1
         IF(IEFF.NE.IPOS) CARD(IEFF:IEFF) = CARD(IPOS:IPOS)
       END IF
  100 CONTINUE
*.Fill end with trailing blanks
      DO 200 IPOS = IEFF+1,LENGTH
        CARD(IPOS:IPOS) = ' ' 
  200 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Left adjusted character string '
        WRITE(6,'(1X,A)') CARD
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE MATCAS(CIN,COUT,NROWI,NROWO,IROWO1,NGCOL,ISCA,SCASGN)
*
* COUT(IR+IROWO1-1,ISCA(IC)) =
* COUT(IR+IROWO1-1,ISCA(IC)) + CIN(IR,IC)*SCASGN(IC)
* (if IGAT(IC).ne.0)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CIN(NROWI,*),COUT(NROWO,*)
      INTEGER ISCA(*)
      DIMENSION SCASGN(*)
*
      MAXCOL = 0
      DO 100 IC = 1, NGCOL
        IF(ISCA(IC).NE.0) THEN
          ICEXP = ISCA(IC)
          MAXCOL = MAX(MAXCOL,ICEXP)
          SIGN = SCASGN(IC)
          DO 50 IR = 1,NROWI
            COUT(IR+IROWO1-1,ICEXP) =
     &      COUT(IR+IROWO1-1,ICEXP) + SIGN*CIN(IR,IC)
   50     CONTINUE
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Output from MATCAS '
        CALL WRTMT_LU(COUT,NROWO,MAXCOL,NROWO,MAXCOL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE MATCG(CIN,COUT,NROWI,NROWO,IROWI1,NGCOL,IGAT,GATSGN)
*
* Gather columns of CIN with phase
*
* COUT(IR,IC) = GATSGN(IC)*CIN(IR+IROWI1-1,IGAT(IC)) if IGAT(IC) .ne.0
* COUT(IR,IC) = 0                           if IGAT(IC) .ne.0
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IGAT(*)
      DIMENSION GATSGN(*)
      DIMENSION CIN(NROWI,*),COUT(NROWO,*)
*
C?    write(6,*) ' MATCG NROWI,NROWO,IROWI1,NGCOL '
C?    write(6,*)         NROWI,NROWO,IROWI1,NGCOL
      DO 100 IG = 1, NGCOL
C?      write(6,*) ' igat,sign ',IGAT(IG),GATSGN(IG)
        IF(IGAT(IG).EQ.0) THEN
          DO 20 IR = 1, NROWO
            COUT(IR,IG)=0.0D0
   20     CONTINUE
        ELSE
         IGFRM = IGAT(IG)
         SIGN = GATSGN(IG)
         DO 30 IR = 1, NROWO
           COUT(IR,IG) = SIGN*CIN(IROWI1-1+IR,IGFRM)
   30    CONTINUE
        END IF
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Column gathered matrix '
        CALL WRTMT_LU(COUT,NROWO,NGCOL,NROWO,NGCOL)
      END IF
*
      RETURN
      END
***********************************************************************

      SUBROUTINE MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION   MATRIX(MATDIM,MATDIM),VECIN(2),VECOUT(2)
C
C     VECOUT=MATRIX*VECIN FOR ITRNSP=0
C     VECOUT=MATRIX(TRANSPOSED)*VECIN FOR ITRNSP .NE. 0
C
      DO 10 I=1,NDIM
   10 VECOUT(I)=0.0D0
      IF(ITRNSP.EQ.0) THEN
C
       DO 100 J=1,NDIM
        VECINJ=VECIN(J)
        DO 90 I=1,NDIM
         VECOUT(I)=VECOUT(I)+MATRIX(I,J)*VECINJ
   90   CONTINUE
  100  CONTINUE
      END IF
C
      IF(ITRNSP.NE.0) THEN
       DO 200 I=1,NDIM
        X=0.0D0
        DO 190 J=1,NDIM
         X=X+MATRIX(J,I)*VECIN(J)
  190   CONTINUE
        VECOUT(I)=X
  200  CONTINUE
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE MATVCC(A,VIN,VOUT,NROW,NCOL,ITRNS)
*
* ITRNS = 0 : VOUT(I) = A(I,J)*VIN(J)
* ITRNS = 1 : VOUT(I) = A(J,I)*VIN(J)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NROW,NCOL)
      DIMENSION VIN(*),VOUT(*)

      if(itrns.eq.0)then
        call dgemv('N',nrow,ncol,1.0d0,a,nrow,vin,1,1.0d0,vout,1)
      else
        call dgemv('T',nrow,ncol,1.0d0,a,nrow,vin,1,1.0d0,vout,1)
      end if
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      FUNCTION IOFF_SYM_DIST(ISYM,NGASL,IOFF,MAXVAL,MINVAL)
*
* A ts block of string is given and the individual
* symmetrydisrtributions has been obtained ( for example
* by TS_SYM_PNT)
*
* Obtain offset for symmetrycombination defined by ISYM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISYM(*),IOFF(*),MAXVAL(*),MINVAL(*)
* Address in IOFF is
*     1
*     + (ISM1-MINVAL(1))
*     + (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     + (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
C     write(6,*) ' Isym and minval '
C     call iwrtma(isym,1,ngasl,1,ngasl)
C     call iwrtma(minval,1,ngasl,1,ngasl)
      I = 1
      IMULT = 1
      DO IGAS = 1, NGASL-1
        I = I + (ISYM(IGAS)-MINVAL(IGAS)) * IMULT
        IMULT = IMULT*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
C       write(6,*) ' igas i imult ',igas,i,imult
      END DO
      IOFF_SYM_DIST=IOFF(I)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IOFF_SYM_DIST'
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Address and offset ',I,IOFF_SYM_DIST
        WRITE(6,*) ' Symmetry distribution : ', (ISYM(J),J=1,NGASL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE MXMNOC_SPGP(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,
     &                       NTESTG)
*
* Construct accumulated MAX and MIN arrays for a GAS supergroup
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Output
      DIMENSION MINEL(*),MAXEL(*)
*. Input
      INTEGER NORBFTP(*),NELFTP(*)
*
      NTESTL = 00
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==========='
        WRITE(6,*) ' MXMNOC_SPGP'
        WRITE(6,*) ' ==========='
        WRITE(6,*)
C?      WRITE(6,*) ' NORBFTP : '
C?      CALL IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
      END IF
*
      DO IORBTP = 1, NORBTP
*. Max and min at start of this type and at end of this type
        IF(IORBTP.EQ.1) THEN
          IORB_START = 1
          IORB_END = NORBFTP(1)
          NEL_START = 0
          NEL_END   = NELFTP(1)
        ELSE
          IORB_START =  IORB_START + NORBFTP(IORBTP-1)
          IORB_END   =  IORB_START + NORBFTP(IORBTP)-1
          NEL_START = NEL_END
          NEL_END   = NEL_START + NELFTP(IORBTP)
        END IF
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' IORBTP,IORB_START-IORB_END,NEL_START,NEL_END '
          WRITE(6,*)   IORBTP,IORB_START-IORB_END,NEL_START,NEL_END
        END IF
*
        DO IORB = IORB_START, IORB_END
          MAXEL(IORB) = MIN(IORB,NEL_END)
          MINEL(IORB) = NEL_START
          IF(NEL_END-MINEL(IORB).GT. IORB_END-IORB)
     &    MINEL(IORB) = NEL_END - ( IORB_END - IORB )
        END DO
      END DO
*
      IF( NTEST .GE. 100 ) THEN
        NORB = IELSUM(NORBFTP,NORBTP)
        WRITE(6,*) ' MINEL : '
        CALL IWRTMA(MINEL,1,NORB,1,NORB)
        WRITE(6,*) ' MAXEL : '
        CALL IWRTMA(MAXEL,1,NORB,1,NORB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      FUNCTION NDXFSM(NSMOB,NSMSX,MXPOBS,NO1PS,NO2PS,NO3PS,NO4PS,
     &         IDXSM,ADSXA,SXDXSX,IS12,IS34,IS1234,IPRNT)
*
* Number of double excitations with total symmetry IDXSM
*
* IS12 (0,1,-1)   : Permutational symmetry between index 1 and 2
* IS34 (0,1,-1)   : Permutational symmetry between index 3 and 3
* IS1234 (0,1,-1) : permutational symmetry between index 12 and 34
*
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*. Specific input
      INTEGER NO1PS(*),NO2PS(*),NO3PS(*),NO4PS(*)
*
      MDX = 0
      DO 200 I12SM = 1, NSMSX
        DO 190 I1SM = 1, NSMOB
          I2SM = ADSXA(I1SM,I12SM)
          IF(IS12.NE.0.AND.I1SM.LT.I2SM) GOTO 190
          IF(IS12.EQ.0) THEN
           I12NUM = (I1SM-1)*NSMSX+I2SM
          ELSE
           I12NUM =  I1SM*(I1SM+1)/2+I2SM
          END IF
          IF(IS12.EQ.0.OR.I1SM.NE.I2SM) THEN
            N12 = NO1PS(I1SM)*NO2PS(I2SM)
          ELSE IF(IS12.EQ.1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)+1)/2
          ELSE IF(IS12.EQ.-1.AND.I1SM.EQ.I2SM) THEN
            N12 = NO1PS(I1SM)*(NO1PS(I1SM)-1)/2
          END IF
          I34SM = SXDXSX(I12SM,IDXSM)
          DO 90 I3SM = 1, NSMOB
            I4SM = ADSXA(I3SM,I34SM)
            IF(IS34.NE.0.AND.I3SM.LT.I4SM) GOTO 90
            IF(IS34.EQ.0) THEN
             I34NUM = (I3SM-1)*NSMSX+I4SM
            ELSE
             I34NUM =  I3SM*(I3SM+1)/2+I4SM
            END IF
            IF(IS1234.NE.0.AND.I12NUM.LT.I34NUM) GOTO 90
            IF(IS34.EQ.0.OR.I3SM.NE.I4SM) THEN
            N34 = NO3PS(I3SM)*NO4PS(I4SM)
            ELSE IF(IS34.EQ.1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)+1)/2
            ELSE IF(IS34.EQ.-1.AND.I3SM.EQ.I4SM) THEN
              N34 = NO3PS(I3SM)*(NO3PS(I3SM)-1)/2
            END IF
            IF(IS1234.EQ.0.OR.I12NUM.NE.I34NUM) THEN
              MDX = MDX + N12 * N34
            ELSE IF( IS1234.EQ.1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12+1)/2
              ELSE IF( IS1234.EQ.-1.AND.I12NUM.EQ.I34NUM) THEN
              MDX =  MDX + N12*(N12-1)/2
            END IF
C?          WRITE(6,*) ' I1SM I2SM I3SM I4SM MDX '
C?          WRITE(6,*)   I1SM,I2SM,I3SM,I4SM,MDX
   90       CONTINUE
C 100     CONTINUE
  190   CONTINUE
  200 CONTINUE
*
      NDXFSM = MDX
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Number of double excitations obtained ', MDX
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE NEWTYPS(INSPGP,IACOP,ITPOP,NOP,
     &           NGAS,ISPGPAN,ISPGPCR,OUTSPGP)
*
* Strinf of operators X supergroup => new supergroup
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ISPGPAN(NGAS,*),ISPGPCR(NGAS,*)
*. Specific input
      DIMENSION IACOP(NOP),ITPOP(NOP)
*. Output
      INTEGER OUTSPGP
*
      OUTSPGP = INSPGP
      DO IOP = 1, NOP
        IF(IACOP(IOP).EQ.1) THEN
          OUTSPGP = ISPGPAN(ITPOP(IOP),OUTSPGP)
        ELSE
          OUTSPGP = ISPGPCR(ITPOP(IOP),OUTSPGP)
        END IF
*
        IF(OUTSPGP.EQ.0) THEN
C         WRITE(6,*) ' NEWTYPS, cul de sac'
C         WRITE(6,*) ' undefined supergroup type '
C         WRITE(6,*) ' String of operator : IAC,ITPOP'
C         CALL IWRTMA(IACOP,1,NOP,1,NOP)
C         CALL IWRTMA(ITPOP,1,NOP,1,NOP)
          GOTO 1001
        END IF
      END DO
*
 1001 CONTINUE
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' NEWTYPS : IAC and ITPOP'
        CALL IWRTMA(IACOP,1,NOP,1,NOP)
        CALL IWRTMA(ITPOP,1,NOP,1,NOP)
        WRITE(6,*) ' Input and Output group ',INSPGP,OUTSPGP
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      FUNCTION NSXFSM(NSMOB,MXPOBS,NO1PS,NO2PS,ISXSM,ADSXA,
     &                ISYM,IPRNT)
*
* Number of single excitations of symmetry ISXSM
*
* ISYM = 0 : All symmetry allowed excitations
* ISYM = 1 : Only excitations a+iaj with I.ge.J
* ISYM =-1 : Only excitations a+iaj with I.gt.J
      INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER NO1PS(*),NO2PS(*)
*
      MSXFSM = 0
C?    WRITE(6,*) ' NSMOB ',NSMOB
      DO 100 IO1SM = 1,NSMOB
        IO2SM = ADSXA(IO1SM,ISXSM)
C?      WRITE(6,*) ' IO1SM,IO2SM',IO1SM,IO2SM
        IF(ISYM.EQ.0.OR.IO1SM.GT.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*NO2PS(IO2SM)
        ELSE IF( ISYM.EQ. 1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)+1)/2
        ELSE IF( ISYM.EQ.-1 .AND. IO1SM.EQ.IO2SM) THEN
          MSXFSM = MSXFSM + NO1PS(IO1SM)*(NO1PS(IO1SM)-1)/2
        END IF
  100 CONTINUE
*
      NSXFSM = MSXFSM
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)

      IF(NTEST.NE.0) THEN
        WRITE(6,*)
     &  ' Number of single excitations of symmetry ',ISXSM,',',NSXFSM
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Output
      SUBROUTINE NXTIJ(I,J,NI,NJ,IJSM,NONEW)
*
* An ordered pair (I,J) is given ,I.LE.NI,J.LE.NJ
*
* Find next pair, if IJSM .ne. 0 ,I .ge. J
*
      NONEW = 0
  100 CONTINUE
      IF(I.LT.NI) THEN
        I = I + 1
      ELSE
        IF(J.LT.NJ) THEN
          I = 1
          J = J+1
        ELSE
          NONEW = 1
          GOTO 101
        END IF
      END IF
      IF(IJSM.NE.0.AND.I.LT.J) GOTO 100
  101 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' next (i,j) pair ', I,J
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* A set of numbers INUM(I),I=1,NELMNT is
* given. Find next compound number.
* Digit I must be in the range MINVAL,MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen Oct 1994
*
*. Input
      DIMENSION MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTNUM '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IF(NELMNT.EQ.0) THEN
       NONEW = 1
       GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = 1
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE NXTNUM3(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* A set of numbers INUM(I),I=1,NELMNT is given.
* Find next compound number.
* Digit I must be in the range MINVAL(I),MAXVAL(I).
*
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen, Oct 1994
*
*. Input
      DIMENSION MINVAL(*),MAXVAL(*)
*. Input and output
      DIMENSION INUM(*)
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input number in NXTNUM '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      IF(NELMNT.EQ.0) THEN
        NONEW = 1
        GOTO 1001
      END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF(INUM(IPLACE).LT.MAXVAL(IPLACE)) THEN
          INUM(IPLACE) = INUM(IPLACE) + 1
          NONEW = 0
          GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
          DO JPLACE = 1, IPLACE
            INUM(JPLACE) = MINVAL(JPLACE)
          END DO
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
          NONEW = 1
          GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
* Find next alpha/beta occupation type and symmetry.
* Reducing redundant code
*   Timo Fleig, Sep. 16, 2002
*
      subroutine nxt_tts(IORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
*
      implicit real*8(A-H,O-Z)
*
      IF(IORD.EQ.1) THEN
*. Old order : IB, IA, ISM  (leftmost inner loop )
        IF(IB.LT.NOCTPB) THEN
          IB = IB + 1
        ELSE
          IB = 1
          IF(IA.LT.NOCTPA) THEN
            IA = IA+ 1
          ELSE
            IA = 1
            IF(ISM.LT.NSMST) THEN
              ISM = ISM + 1
            ELSE
              IFINI = 1
            END IF
          END IF
        END IF
      ELSE IF(IORD.EQ.2) THEN
*. Else : New order : ISM,IB,IA (leftmost inner loop )
        IF(ISM.LT.NSMST) THEN
          ISM = ISM + 1
        ELSE
          ISM = 1
          IF(IB.LT.NOCTPB) THEN
            IB = IB + 1
          ELSE
            IB = 1
            IF(IA.LT.NOCTPA) THEN
              IA = IA + 1
            ELSE
              IFINI = 1
            END IF
          END IF
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
*     Find back to the first TT-1 block to start next batch from.
*     Stefan Knecht, July 02, 2007.
*
      subroutine bck_tts(IORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,NBLOCK,
     &                   LENGTH,LENGTHP,NBATCH,LUWRT)
*
      IMPLICIT REAL*8           (A-H,O-Z)
      dimension LBATCH(*), LEBATCH(*), IBATCH(8,*)
*
      IF( IORD .eq. 2 ) THEN
*
 100    CONTINUE
*
        IF( ISM .gt. 1 ) THEN
          ISM = ISM - 1 
          LBATCH(NBATCH)  = LBATCH(NBATCH) - 1
          LEBATCH(NBATCH) = LEBATCH(NBATCH) - IBATCH(8,IBLOCK)
          IBLOCK = IBLOCK - 1
          NBLOCK = NBLOCK - 1
          GOTO 100 
        ELSE
*
CSK        WRITE(LUWRT,*) ' reset was successful'
CSK        WRITE(LUWRT,*) ' IBLOCK, NBLOCK are now',IBLOCK,NBLOCK
*
        END IF
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE NXTBLK(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,IBLTP,IDC,
     &                  NONEW,IOCOC,ITTSS_ORD)
*
* Obtain allowed block following IATP IBTP IASM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IBLTP(*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
*
*.Initialize
      ISM = IASM
      IA = IATP
      IB = IBTP
      NONEW = 0
* NEW START
*. Loop over blocks in batch
 1000 CONTINUE
      call nxt_tts(ITTSS_ORD,IA,IB,ISM,NONEW,NOCTPA,NOCTPB,NSMST)
* NEW END
*.Next block
      IF(NONEW.EQ.1) GOTO 1001
*. Should this block be included
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.0) GOTO 1000
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
 1001 CONTINUE
*
      IATP = IA
      IBTP = IB
      IASM = ISM
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)')
     &  ' NXTBLK : ISM IA IB NONEW ', IASM,IA,IB,NONEW
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE NXTORD_LUCI(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
*
* An ordered set of numbers INUM(I),I=1,NELMNT is
* given in strictly ascending order. Values of INUM(*) is
* restricted to the interval MINVAL,MAXVAL .
*
* Find next higher number.
*
* NONEW = 1 on return indicates that no additional numbers
* could be obtained.
*
* Jeppe Olsen May 1989
*
      DIMENSION INUM(*)
*
       NTEST = 0
       IF( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' Initial number in NXTORD '
         CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
       END IF
*
      IPLACE = 0
 1000 CONTINUE
        IPLACE = IPLACE + 1
        IF( IPLACE .LT. NELMNT .AND.
     &      INUM(IPLACE)+1 .LT. INUM(IPLACE+1)
     &  .OR.IPLACE.EQ. NELMNT .AND.
     &      INUM(IPLACE)+1.LE.MAXVAL) THEN
              INUM(IPLACE) = INUM(IPLACE) + 1
              NONEW = 0
              GOTO 1001
        ELSE IF ( IPLACE.LT.NELMNT) THEN
              IF(IPLACE .EQ. 1 ) THEN
                INUM(IPLACE) = MINVAL
              ELSE
                INUM(IPLACE) = INUM(IPLACE-1) + 1
              END IF
        ELSE IF ( IPLACE. EQ. NELMNT ) THEN
              NONEW = 1
              GOTO 1001
        END IF
      GOTO 1000
 1001 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' New number '
        CALL IWRTMA(INUM,1,NELMNT,1,NELMNT)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
       SUBROUTINE ORBINH1(IORBINH1,NTOOBS,NTOOB,NSMOB)
*
* Obtain array of 2 orbital indeces,
* for symmetry packed matrices, stored as lower halfs
*
* resulting indeces are with respect to start of given symmetry block
* while input orbital indeces are absolute and in symmetry order
*
* Jeppe Olsen, March 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION NTOOBS(NSMOB)
*. output
      DIMENSION IORBINH1(NTOOB,NTOOB)
*
C?    WRITE(6,*) ' ORBINH1 speaking '
C?    WRITE(6,*) ' NSMOB NTOOB ',NSMOB,NTOOB
C?    WRITE(6,*) ' NTOOBS '
C?    CALL IWRTMA(NTOOBS,1,NSMOB,1,NSMOB)

*. Loop over symmetries of orbitals

      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)
        END IF
        DO JSM = 1, NSMOB
          IF(JSM.EQ.1) THEN
            JOFF = 1
          ELSE
            JOFF = JOFF + NTOOBS(JSM-1)
          END IF
C?        WRITE(6,*) ' ISM JSM IOFF JOFF', ISM,JSM,IOFF,JOFF
          DO IORB = 1, NTOOBS(ISM)
            IABS = IOFF -1 + IORB
            DO JORB = 1, NTOOBS(JSM)
              JABS = JOFF -1 + JORB
C?            write(6,*) ' IORB JORB IABS JABS ',IORB,JORB,IABS,JABS
              IF(ISM.GT.JSM) THEN
                INDEX = (IORB-1)*NTOOBS(JSM) + JORB
              ELSE IF(ISM.EQ.JSM) THEN
                IF(IORB.GE.JORB) THEN
                  INDEX = IORB*(IORB-1)/2 + JORB
                ELSE
                  INDEX = JORB*(JORB-1)/2 + IORB
                END IF
              ELSE IF(ISM.LT.JSM) THEN
                INDEX = (JORB-1)*NTOOBS(ISM) + IORB
              END IF
              IORBINH1(IABS,JABS) = INDEX
            END DO
          END DO
*. End of loops over orbital indeces
        END DO
      END DO
*. End of loop over orbital symmetries
*
      NTEST = 000
      IF(NTEST .GE. 100 ) THEN
        WRITE(6,*) ' IORBINH1 matrix delivered from ORBINH1'
        CALL IWRTMA(IORBINH1,NTOOB,NTOOB,NTOOB,NTOOB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRMBLK(IDC,ISGV,IASM,IBSM,IATP,IBTP,PS,PL,
     &                  JATP,JBTP,JASM,JBSM,ISGN,ITRP,NPERM)
*
* A block of CI coefficients defined by by IATP,IASM,IBTP,IBSM is given
*
* Obtain the number of other blocks that can be obtained by spin
* and relection symmetry.
*
* Jeppe Olsen, July 1993
*
* =====
* Output
* =====
* JATP(I),JASM(I),JBTP(I),JBSM(I) Indeces for Block I
* NPERM : Number of blocks  that can be obtained
* ITRP(I) = 1 => block should     be transposed
*         = 0 => block should not be transposed
* ISGN   : Sign to multiply previous block with to getnew sign
*
*
* There are four types of permutations
*
*    operation   *      JASM  *      JBSM  * JATP * JBTP * Iperm * Sign *
*   *********************************************************************
*   * Identity   *      IASM  *      IBSM  * IATP * IBTP *   0   * 1    *
*   * Ml         * ISGV(IASM) * ISGV(IBSM) * IATP * IBTP *   0   * PL   *
*   * Ms         *      IBSM  *      IASM  * IBTP * IATP *   1   * PS   *
*   * Ms+Ml      * ISGV(IBSM) * ISGV(IASM) * IBTP * IATP *   1   * PS PL*
*   *********************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*.Input
      DIMENSION ISGV(*)
*.Output
      DIMENSION JATP(4),JBTP(4),JASM(4),JBSM(4),ISGN(4),ITRP(4)
*
      SAVE LSIGN, LTRP
      DATA LSIGN, LTRP/0,0/
*
      NPERM = 0
      DO 100 IPERM = 1, 4
        ISET = 0
        IF(IPERM.EQ.1) THEN
*
* Identity operation
*
          KASM = IASM
          KBSM = IBSM
          KATP = IATP
          KBTP = IBTP
          KSIGN = 1
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.2.AND.(IDC.EQ.3.OR.IDC.EQ.4)) THEN
*
* Ml reflection
*
          KASM = ISGV(IASM)
          KBSM = ISGV(IBSM)
          KATP = IATP
          KBTP = IBTP
          IF(PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PL .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.3.AND.(IDC.EQ.2.OR.IDC.EQ.4)) THEN
*
* Ms reflection
*
          KASM = IBSM
          KBSM = IASM
          KATP = IBTP
          KBTP = IATP
          IF(PS.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        ELSE IF(IPERM.EQ.4 .AND. IDC.EQ.4) THEN
*
* Ms Ml  reflection
*
          KASM = ISGV(IBSM)
          KBSM = ISGV(IASM)
          KATP = IBTP
          KBTP = IATP
          IF(PS*PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        END IF
*
        IF(ISET.EQ.1) THEN
*. A new permutation was found, check and see if it was obtained previously
          INEW = 1
          DO 50 LPERM = 1, NPERM
            IF(JATP(LPERM).EQ.KATP  .AND. JASM(LPERM).EQ.KASM .AND.
     &         JBTP(LPERM).EQ.KBTP  .AND. JBSM(LPERM).EQ.KBSM) INEW = 0
   50     CONTINUE
          IF(INEW.EQ.1) THEN
*. The permutation was new, add it to the list
            NPERM = NPERM + 1
            JASM(NPERM) = KASM
            JBSM(NPERM) = KBSM
            JATP(NPERM) = KATP
            JBTP(NPERM) = KBTP
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KSIGN.EQ.LSIGN))THEN
              ISGN(NPERM) = 1
            ELSE
              ISGN(NPERM) = -1
            END IF
            LSIGN = KSIGN
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KTRP.EQ.LTRP))THEN
              ITRP(NPERM) = 0
            ELSE
              ITRP(NPERM) = 1
            END IF
            LTRP = KTRP
          END IF
        END IF
  100 CONTINUE
*
*. Should the block be trnasposed or scaled to return to initial form
      ITRP(NPERM+1) = LTRP
      ISGN(NPERM+1) = LSIGN
      IFNSGN = LSIGN
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)') ' Blocks obtained from IASM IBSM IATP IBTP ',
     &  IASM,IBSM,IATP,IBTP
        WRITE(6,*)
        WRITE(6,'(A)') ' JASM JBSM JATP JBTP Isgn Itrp  '
        WRITE(6,*)
        DO 10 IPERM = 1, NPERM
          WRITE(6,'(2x,6I4)') JASM(IPERM),JBSM(IPERM),JATP(IPERM),
     &                        JBTP(IPERM),ISGN(IPERM),ITRP(IPERM)
   10   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRTITL(LINES)
*     
* Print title cards
*
      CHARACTER*72 LINES
      DIMENSION LINES(3)
      CHARACTER*80 STARS
*
      STARS(1:1) = ' '
      DO 80 I = 2, 80
        STARS(I:I) = '*'
   80 CONTINUE
      WRITE(6,'(A)') STARS
      WRITE(6,'(3A)') ' *  ',LINES(1),'  *'
      WRITE(6,'(3A)') ' *  ',LINES(2),'  *'
      WRITE(6,'(3A)') ' *  ',LINES(3),'  *'
      WRITE(6,'(A)') STARS
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRTSTR(ISTR,NEL,NSTR)
*
* Print NSTR strings each containing NEL electrons
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ISTR(NEL,NSTR)
*
      DO JSTR = 1, NSTR
        WRITE(6,'(/,A,I6,A,4X,10(2X,I4),/,(20X,10(2X,I4)))' )
     &   ' String ',JSTR,' : ',(ISTR(IEL,JSTR),IEL=1,NEL)
      END DO
*
      RETURN 
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PRTRCCTOS(IMAT,IDIM1,IDIM2,LUPRI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IMAT(IDIM1,IDIM2)
*
      DO I = 1, IDIM2
         WRITE(LUPRI,*) ' priniting IMAT=RCCTOS for pair #I',I
         CALL IWRTMAMN(IMAT(1,I),1,IDIM1,1,IDIM1,LUPRI)
      END DO
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE RFTTS(BLOCKSI,BLOCKSO,IBLOCK,NBLOCK,
     &                 ICOPY,NSMST,NOCTPA,NOCTPB,
     &                 NSASO,NSBSO,IDC,PS,IWAY,IPRNT)
*
* Reformat between determinant and combination form of
* matrices. No scaling is performed .
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* Combination storage mode is defined BY IDC
*
*. Jeppe Olsen, August 1995
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKSI(*),BLOCKSO(*)
      INTEGER IBLOCK(8,NBLOCK)
*

      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
*
      LENGTH = 0
      IF(IWAY.EQ.1) THEN
        ISCI = 1
        ISCO = 2
      ELSE
        ISCI = 2
        ISCO = 1
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTTTS(BLOCKSI,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,ISCI)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF(IWAY.EQ.1) THEN
          IOFFI = IBLOCK(5,JBLOCK)
          IOFFO = IBLOCK(6,JBLOCK)
        ELSE
          IOFFO = IBLOCK(5,JBLOCK)
          IOFFI = IBLOCK(6,JBLOCK)
        END IF
*. Is this block diagonal in packed form
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
*. Number of elements in output block
        IF(IPACK .EQ. 1 .AND. ISCO.EQ.2 ) THEN
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
C?     write(6,*)
C?   & ' RFTTS : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
C?     WRITE(6,*)
C?   & ' RFTTS : NIA NIB IOFFI,IOFFO',NIA,NIB,IOFFI,IOFFO
*
        IF(IPACK.EQ.0) THEN
*. Just copy
          CALL COPVEC(BLOCKSI(IOFFI),BLOCKSO(IOFFO),NELMNT)
        ELSE
          IF(IWAY.EQ.1) THEN
*. unpacked => packed
C TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
            CALL TRIPK3(BLOCKSI(IOFFI),BLOCKSO(IOFFO),1,NIA,NIA,PS)
          ELSE
*. Packed => unpacked
            CALL TRIPK3(BLOCKSO(IOFFO),BLOCKSI(IOFFI),2,NIA,NIA,PS)
          END IF
        END IF
        LENGTH = LENGTH + NELMNT
        END IF
      END DO
*
      IF(ICOPY.NE.0) THEN
        CALL COPVEC(BLOCKSO,BLOCKSI,LENGTH)
      END IF
*
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*) ' Information from RFTTS  '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Output vector '
        CALL WRTTTS(BLOCKSO,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,ISCO)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SCDTTS(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
*
* Scale batch of
* blocks between determinant and combination form
*
*
* IWAY = 1 : dets to combs
* IWAY = 2 : combs to dets
*
* The blocks are assumed to be in packed form !!
*
*. Jeppe Olsen, August 1995
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
*
C?    LOGICAL DIAGBL
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRNT)
      IF( NTEST .GT. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from SCDTTS '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input vector '
        CALL WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,2)
      END IF
*
      SQ2 = SQRT(2.0D0)
      SQ2I = 1.0D0/SQ2
*
      IBASE = 1
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IOFFP= IBLOCK(6, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*. Is this block diagonal in packed form
        IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA   = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
        IF(IPACK .EQ. 1 ) THEN
          NELMNT =  NIA*(NIA+1)/2
        ELSE
          NELMNT =  NIA*NIB
        END IF
*Ms combinations
        IF(IDC.EQ.2) THEN
          IF(IWAY.EQ.1) THEN
            FACTOR = SQ2
          ELSE
            FACTOR = SQ2I
          END IF
          CALL SCALVE(BLOCKS(IOFFP),FACTOR,NELMNT)
          IF(IPACK.EQ.1 ) THEN
            FACTOR = 1.0D0/FACTOR
            CALL SCLDIA(BLOCKS(IOFFP),FACTOR,NIA,1)
          END IF
        END IF
*
        END IF
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Output vector '
        CALL WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,2)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SCLDIA(A,FACTOR,NDIM,IPACK)
*
* scale diagonal of square matrix A
*
* IPACK = 0 : full matrix
* IPACK .NE. 0 : Lower triangular packed matrix
*                assumed packed columnwise !!!!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*)
*
      IF( IPACK .EQ. 0 ) THEN
        DO 100 I = 1,NDIM
          II = (I-1)*NDIM + I
          A(II) = A(II) * FACTOR
  100   CONTINUE
      ELSE
        II = 1
        DO 200 I = 1, NDIM
          A(II) = A(II) * FACTOR
          II = II + NDIM - I + 1
  200   CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SDCMRF(CSD,CCM,IWAY,IATP,IBTP,IASM,IBSM,NA,NB,
     &                  IDC,PS,PL,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
*
* Change a block of coefficients between combination format and
* Slater determinant format
*
*     IWAY = 1 : SD => Combinations
*     IWAY = 2 : Combinations => SD
*
* Input
* =====
* CSD : Block in determinant form
* CCM : Block in combination  form
* IWAY : as above
* IATP,IBTP : type of alpha- and beta- string
* NA,NB : Number of alpha- and beta- strings
* IDC  : Combination type
* PS   : Spin combination sign
* PL   : Ml   combination sign
* ISGVST : Ml reflection of strings
*
*
* If ISCALE .EQ. 0, no overall scaling is performed,
*                   the overall scale factor is returned
*                   as SCLFAC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CSD(*),CCM(*),ISGVST(*)
*
      SQRT2  = SQRT(2.0D0)
      SQRT2I = 1.0D0/SQRT2
*
*. Is combination array packed ?
*
      SCLFAC = 1.0D0
      IPACK = 0
      FACTOR = 1.0D0
*
      IF(IDC.EQ.2.OR.IDC.EQ.4) THEN
         SIGN_T = PS
         FACTOR = SQRT2
         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) IPACK = 1
      ELSE IF( IDC.EQ.4.AND.IASM.EQ.ISGVST(IBSM)) THEN
        IF(IATP.EQ.IBTP) IPACK = 1
        SIGN_T = PS*PL
        FACTOR = 2.0D0
      END IF
*
      LDET = NA * NB
      IF( IPACK.EQ.0) THEN
        LCOMB = LDET
      ELSE
        LCOMB = NA*(NA+1)/2
      END IF
      IF(IDC.EQ.4.AND.IPACK.EQ.0) FACTOR = SQRT2
      IF(IWAY.EQ.2) FACTOR = 1.0D0/FACTOR
*
*. SD => combination transformation
*
      IF(IWAY .EQ. 1 ) THEN
        IF(IPACK.EQ.1) THEN
*. Pack to triangular form
          CALL TRIPK3(CSD,CCM,1,NA,NA,SIGN_T)
C              TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN_T)
        ELSE
          CALL DCOPY(NA*NB,CSD,1,CCM,1)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL DSCAL(LCOMB,FACTOR,CCM,1)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1 ) THEN
            CALL SCLDIA(CCM,SQRT2I,NA,1)
          END IF
        END IF
      END IF
*
*. Combination => SD transformation
*
      IF(IWAY.EQ.2) THEN
        IF(IPACK.EQ.1) THEN
*. Unpack from triangular form
          CALL TRIPK3(CSD,CCM,2,NA,NA,SIGN_T)
        ELSE
          CALL DCOPY(NA*NB,CCM,1,CSD,1)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL DSCAL(LDET,FACTOR,CSD,1)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1) THEN
             CALL SCLDIA(CSD,SQRT2,NA,0)
          END IF
        END IF
      END IF
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
C     IF(NTEST.NE.0.AND.IWAY.EQ.1) THEN
        WRITE(6,*) ' Information from SDCMRF '

        WRITE(6,'(A,6I4)') ' IWAY IATP IBTP IASM IBSM IDC ',
     &                   IWAY,IATP,IBTP,IASM,IBSM,IDC
        WRITE(6,'(A,I4,3X,2E13.8)') ' IPACK FACTOR SIGN_T',
     &  IPACK,FACTOR,SIGN_T
        IF(NTEST.GE. 100 ) THEN
          WRITE(6,*) ' Slater determinant block '
          CALL WRTMT_LU(CSD,NA,NB,NA,NB)
          WRITE(6,*)
          WRITE(6,*) ' Combination block '
          IF(IPACK.EQ.1) THEN
            CALL PRSM2(CCM,NA)
          ELSE
            CALL WRTMT_LU(CCM,NA,NB,NA,NB)
          END IF
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXST2,JEXSTV,IPRT)
C
C PURPOSE: FIND THE KEXSTV LOWEST VALUES IN WRK(KZVAR)
C          CHECK IF THE LAST ELEMENTS ARE DEGENERATE
C          JEXSTV IS THE NUMBER OF SORTED ELEMENTS WHERE
C          NO DEGENERACIES OCCUR AMONG THE HIGHEST ELEMENTS
* INPUT
*======
* WRK : ARRAY TO BE SORTED
* KEXST2 : NUMBER OF ELEMENTS TO BE OBTAINED  + 1
* KZVAR : LENGTH OG WRK
* OUTPUT
*=======
* STVAL  : STVAL(I) IS VALUE OF SORTED ELEMENT I
* ISTART : SCATTER POINTER ISTART(I) IS ADRESS IN FULL LIST OF
*          SORTED ELEMENT I
* KZVAR : LENGTH OG WRK
* JEXSTV : NUMBER OF ELEMENTS RETURNED , CAN BE LESS THAN
*          KEXSTV-1 IF DEGENERNCIES OCCURS  AMONG THE LAST ELEMENTS
*
C
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION WRK(*),STVAL(*),ISTART(*)
C
      PARAMETER ( TOLEQL=1.0D-6 )
      PARAMETER ( D0=0.0D0 , D1=1.0D0 , DM1 = -1.0D0 )
      LOGICAL FULL
C
      IF(KEXST2.GE.KZVAR) THEN
         FULL = .TRUE.
         KEXSTV = KZVAR
      ELSE
         FULL = .FALSE.
         KEXSTV = KEXST2+1
      END IF
*
      DO 100 K=1,KEXSTV
         ISTART(K) = K
         STVAL(K)  = WRK(K)
 100  CONTINUE
      KK=KEXSTV
C
      DO 210 I=1,KEXSTV
         DO 220 J=I+1,KEXSTV
            IF ((STVAL(J)-STVAL(I)) .LT. D0) THEN
               X =STVAL(I)
               II=ISTART(I)
               STVAL(I) =STVAL(J)
               ISTART(I)=ISTART(J)
               STVAL(J) =X
               ISTART(J)=II
            ENDIF
 220     CONTINUE
 210  CONTINUE
      GO TO 115
C     REPEAT UNTIL ...
 105  CONTINUE
         DO 110 I = KEXSTV,2,-1
            J = I - 1
            IF ((STVAL(J)-STVAL(I)) .GT. D0) THEN
               X =STVAL(I)
               II=ISTART(I)
               STVAL(I) =STVAL(J)
               ISTART(I)=ISTART(J)
               STVAL(J) =X
               ISTART(J)=II
            ELSE
               GO TO 115
            ENDIF
 110     CONTINUE
 115     CONTINUE
         STMAX=STVAL(KEXSTV)
C
 125     CONTINUE
            KK=KK+1
            IF (KK.LE.KZVAR) THEN
               IF ((STMAX-WRK(KK)).GT.D0) THEN
                  ISTART(KEXSTV) = KK
                  STVAL(KEXSTV)  = WRK(KK)
                  GO TO 105
C     ^--------------------
               ENDIF
               GO TO 125
C        ^--------------
            END IF
C
C     Check for degeneracy among diagonal elements
C
      JEXSTV = KEXSTV
      IF(.NOT.FULL) THEN

  160   JEXSTV = JEXSTV - 1
        IF (
     *  (ABS(STVAL(JEXSTV+1)-STVAL(JEXSTV))).LE.TOLEQL) GO TO 160
      END IF
      IF((IPRT.GT.8).AND.(KEXST2.NE.JEXSTV)) WRITE(6,1600)KEXST2,JEXSTV

1600  FORMAT(/' NUMBER OF START VECTORS IS DIMINISHED TO',I5,' FROM',I5)
      IF (IPRT.GE.100) THEN
         WRITE(6,*) '(I,(ISTART(I),STVAL(I)),I=1,JEXSTV)'
      DO 170 I = 1,JEXSTV
            WRITE(6,*) I,ISTART(I),STVAL(I)
  170    CONTINUE
         WRITE(6,*) 'THE FIRST',JEXSTV,' ELEMENTS ARE SELECTED.'
      END IF
C
C     END OF SORVAL
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
*
* Some GAS routines
*
* Nomenclature
*
*. Group of strings : set of strings with a given number of orbitals
*                    in a given GASspace
*
*. Supergroup of strings  : product of NGAS groups, i.e. a string with a
*                    given numb er of electrons in each GAS space
*
*. Type of string : Type is defined by the total number of electrons
*                   in the string. A type will therefore in general
*                   consists of several supergroups
*
      SUBROUTINE SPGP_AC(INSPGRP,NINSPGRP,IOUTSPGRP,NOUTSPGRP,
     &           NGAS,MXPNGAS,IAC,ISPGRP_AC,IBASEIN,IBASEOUT)
*
* Annihilation/Creation mapping of strings
*
* Jeppe Olsen, April 1, 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input : Number of electrons in each gasspace
      INTEGER INSPGRP(MXPNGAS,*),IOUTSPGRP(MXPNGAS,*)
*. Output
      INTEGER ISPGRP_AC(NGAS,*)
*. Check first that supergroups + IAC information is consistent
      NELIN = 0
      NELOUT = 0
      DO IGAS = 1, NGAS
        NELIN  = NELIN + INSPGRP(IGAS,IBASEIN)
        NELOUT = NELOUT + IOUTSPGRP(IGAS,IBASEOUT)
      END DO
      IF(.NOT.((IAC.EQ.1.AND.NELIN.EQ.NELOUT+1).OR.
     &         (IAC.EQ.2.AND.NELIN.EQ.NELOUT-1))) THEN
        WRITE(6,*) ' Inconsistent data provided to SPGRP_AC'
        WRITE(6,*) ' NELIN NELOUT IAC=',NELIN,NELOUT,IAC
        Call Abend2( ' Inconsistent data provided to SPGRP_AC' )
      END IF
*
      DO ISPGRP = IBASEIN, IBASEIN+NINSPGRP-1
        DO IGAS = 1, NGAS
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          END IF
*. Find corresponding supergroup
          ITO = 0
          DO JSPGRP = IBASEOUT, IBASEOUT+NOUTSPGRP-1
            IAMOKAY = 1
            DO JGAS = 1, NGAS
              IF( INSPGRP(JGAS,ISPGRP).NE.IOUTSPGRP(JGAS,JSPGRP))
     &        IAMOKAY=0
            END DO
            IF(IAMOKAY.EQ.1) ITO = JSPGRP
          END DO
          ISPGRP_AC(IGAS,ISPGRP) = ITO
*. And clean up
          IF(IAC.EQ.1) THEN
            INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) + 1
          ELSE IF (IAC.EQ.2) THEN
             INSPGRP(IGAS,ISPGRP) = INSPGRP(IGAS,ISPGRP) - 1
          END IF
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Input supergroups '
        CALL IWRTMA(INSPGRP(1,IBASEIN),NGAS,NINSPGRP,
     &  MXPNGAS,NINSPGRP)
        WRITE(6,*) ' Output supergroups '
        CALL IWRTMA(IOUTSPGRP(1,IBASEOUT),NGAS,NOUTSPGRP,
     &  MXPNGAS,NOUTSPGRP)
      END IF
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Output from SPGP_AC '
       WRITE(6,*) ' ===================='
       WRITE(6,*)
       WRITE(6,*) ' IAC = ', IAC
       WRITE(6,*) ' Mapping : '
       CALL IWRTMA(ISPGRP_AC(1,IBASEIN),NGAS,NINSPGRP,NGAS,NINSPGRP)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,
     &                  ISPGRPCON,IPRNT)
*
* FInd connection matrix for string types
*
* ISPGRPCON(ISPGP,JSPGRP) = 0 => spgrps are identical
*                         = 1 => spgrps are connected by single 
*                                excitation
*      .                  = 2 => spgrps are connected by double 
*                                excitation
*       .              . ge.3 => spgrps are connected by triple or
*        .                       higher excitation
*
* Jeppe Olsen, September 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IELFSPGRP(MXPNGAS,*)
*. output
      DIMENSION ISPGRPCON(NSPGRP,NSPGRP)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
      DO ISPGRP = 1, NSPGRP
        ISPGRPA = IOFSPGRP-1+ISPGRP
        DO JSPGRP = 1, ISPGRP
          JSPGRPA = IOFSPGRP-1+JSPGRP
          IDIF = 0
          DO IGAS = 1, NGAS
            IDIF = IDIF
     &    + ABS(IELFSPGRP(IGAS,ISPGRPA)-IELFSPGRP(IGAS,JSPGRPA))
          END DO
          NEXC = IDIF/2
          ISPGRPCON(ISPGRP,JSPGRP) = NEXC
          ISPGRPCON(JSPGRP,ISPGRP) = NEXC
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) '==================== '
        WRITE(6,*) 'output from SPGRPCON '
        WRITE(6,*) '==================== '
        WRITE(6,*)
        NEXC1 = 0
        NEXC2 = 0
        DO ISPGRP=1, NSPGRP
          DO JSPGRP = 1, NSPGRP
            IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.1) THEN
              NEXC1 = NEXC1 + 1
            ELSE IF(ISPGRPCON(ISPGRP,JSPGRP).EQ.2) THEN
              NEXC2 = NEXC2 + 1
            END IF
          END DO
        END DO
*
        WRITE(6,*)
     &  ' single excitation interactions',NEXC1,
     &   '( ',(NEXC1*100.0D0)/NSPGRP**2,' % ) '
        WRITE(6,*)
     &  ' double excitation interactions',NEXC2,
     &   '( ',(NEXC2*100.0D0)/NSPGRP**2,' % ) '
*
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Supergroup connection matrix '
         CALL IWRTMA(ISPGRPCON,NSPGRP,NSPGRP,NSPGRP,NSPGRP)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
*
* Two supergroups are given. Find single excitations that connects
* these supergroups
*
* Jeppe Olsen, July 1995
*
* Dec 97 : IPHGAS added :
*          Occupations of particle spaces (IPHGAS=2) are allowed to
*          have occupations less than zero in intermidiate steps
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ILTP(NGAS),IRTP(NGAS),IPHGAS(*)
*. Output
      DIMENSION ITP(*),JTP(*)
*. Differences between occupations :
      NCREA = 0
      NANNI = 0
      DO IAS = 1, NGAS
        IF(ILTP(IAS).GT.IRTP(IAS)) THEN
         NCREA = NCREA + ILTP(IAS) - IRTP(IAS)
         ICREA = IAS
        ELSE IF(IRTP(IAS).GT.ILTP(IAS)) THEN
         NANNI = NANNI + IRTP(IAS) - ILTP(IAS)
         IANNI = IAS
        END IF
      END DO
*
      IF(NCREA.GT.1) THEN
*. Sorry : No single excitation connects
        NSXTYP = 0
      ELSE IF(NCREA .EQ. 1 ) THEN
*. Supergroups differ by one single excitation.
        NSXTYP = 1
        ITP(1) = ICREA
        JTP(1) = IANNI
      ELSE IF (NCREA.EQ.0 ) THEN
*. Supergroups are identical, connects with all
*  diagonal excitations.
        NSXTYP = 0
        DO IAS = 1, NGAS
          IF(IRTP(IAS).NE.0.OR.IPHGAS(IAS).EQ.2) THEN
            NSXTYP = NSXTYP + 1
            ITP(NSXTYP) = IAS
            JTP(NSXTYP) = IAS
          END IF
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from SXTYP_GAS : '
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Input left  supergroup '
        CALL IWRTMA(ILTP,1,NGAS,1,NGAS)
        WRITE(6,*) ' Input right supergroup '
        CALL IWRTMA(IRTP,1,NGAS,1,NGAS)
        WRITE(6,*)
     &  ' Number of connecting single excitations ', NSXTYP
        IF(NSXTYP.NE.0) THEN
          WRITE(6,*) ' Connecting single excitations '
          DO ISX = 1, NSXTYP
            WRITE(6,*) ITP(ISX),JTP(ISX)
          END DO
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SWAPVE(VEC1,VEC2,NDIM)
C
C     SWAP ELEMENTS OF VECTORS VEC1 AND VEC2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      CALL DSWAP(NDIM,VEC1,1,VEC2,1)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SYMCM1(ITASK,IOBJ,I1,I2,I12)
*
* Symmetries I1,I2,I12 are related as
* I1*I2 = 12
* IF(ITASK = 1 ) I2 and I12 are known, find I1
* IF(ITASK = 2 ) I1 and I12 are known, find I1
* IF(ITASK = 3 ) I1 and I2 are known , find I12
*
* D2h version , written for compatibility with general symmetry
*
      INTEGER SYMPRO(8,8)
      DATA  SYMPRO/1,2,3,4,5,6,7,8,
     &             2,1,4,3,6,5,8,7,
     &             3,4,1,2,7,8,5,6,
     &             4,3,2,1,8,7,6,5,
     &             5,6,7,8,1,2,3,4,
     &             6,5,8,7,2,1,4,3,
     &             7,8,5,6,3,4,1,2,
     &             8,7,6,5,4,3,2,1 /
*
      IF(ITASK.EQ.1) THEN
        I1 = SYMPRO(I2,I12)
      ELSE IF(ITASK.EQ.2) THEN 
        I2 = SYMPRO(I1,I12) 
      ELSE IF (ITASK.EQ.3) THEN
        I12 = SYMPRO(I1,I2)
      END IF
* 
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SYMCM2(IDUM1,IDUM2,IDUM3,IDUM4,IDUM5)
      WRITE(6,*) ' Entering dummy SYMCM2 is fatal for me'
      Call Abend2( ' Dummy SYMCM2' )
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE TRIPAK_LUCI(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
C   Added Sign multiplied to elements when packing/unpacking
C      Timo Fleig, Dec. 17, 2002
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 I = 1,NDIM
          DO 50  J = 1, I
           APAK(IJ+J) = SIGN * AUTPAK(J,I)
   50     CONTINUE
          IJ = IJ + I
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 I = 1,NDIM
          DO 150  J = 1, I
           AUTPAK(I,J) = SIGN * APAK(IJ+J)
           AUTPAK(J,I) = SIGN * APAK(IJ+J)
  150     CONTINUE
          IJ = IJ + I
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPAK '
        CALL WRTMT_LU(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
*
* Some considerations on cache minimization used for IMET = 2 Loop
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
*
*. Packing : No problem with cache misses
*
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO J = 1,NDIM
          DO I = J , NDIM
            APAK(IJ+I) = AUTPAK(I,J)
          END DO
          IJ = IJ +NDIM-J
        END DO
      END IF
*
* Unpacking : cache misses can occur so two routes
*
      IF( IWAY .EQ. 2 ) THEN
*. Use blocked algorithm
      IMET = 2
C     IMET =2
      IF(IMET.EQ.1) THEN
*. No blocking
        IJ = 0
        DO J = 1,NDIM
          DO I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
          END DO
          IJ = IJ + NDIM-J
        END DO
      ELSE IF (IMET .EQ. 2 ) THEN
*. Blocking
C       LBLK = 40
        LBLK = 40
        NBLK = MATDIM/LBLK
        IF(LBLK*NBLK.LT.MATDIM) NBLK = NBLK + 1
        DO JBLK = 1, NBLK
          IF(JBLK.EQ.1) THEN
            JOFF = 1
          ELSE
            JOFF = JOFF + LBLK
          END IF
          JEND = MIN(JOFF+LBLK-1,MATDIM)
          DO IBLK = JBLK, NBLK
            IF(IBLK.EQ.JBLK) THEN
              IOFF = JOFF
            ELSE
              IOFF = IOFF + LBLK
            END IF
            IEND = MIN(IOFF+LBLK-1,MATDIM)
              DO J = JOFF,JEND
                IF(IBLK.EQ.JBLK) THEN
                  IOFF2 = J
                ELSE
                  IOFF2 = IOFF
                END IF
                IJOFF = (J-1)*MATDIM-J*(J-1)/2
                DO I = IOFF2,IEND
                  AUTPAK(J,I) = SIGN*APAK(IJOFF+I)
                  AUTPAK(I,J) = APAK(IJOFF+I)
                END DO
              END DO
*. End of loop over I and J
            END DO
          END DO
*. End of loop over blocks of I and J
        END IF
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK3 '
        CALL WRTMT_LU(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE TRPAD3(MAT,FACTOR,NDIM)
*
*  MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
*
*. With some considerations of effective cache use for large
*  matrices
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL*8            MAT(NDIM,NDIM)
              FAC2 = 1.0D0 - FACTOR**2
*
C     IWAY = 1
      IWAY = 2
      IF(IWAY.EQ.1) THEN
*
*. No blocking
*
*. Lower half
        DO  J = 1, NDIM
          DO  I = J, NDIM
            MAT(I,J) =MAT(I,J) + FACTOR * MAT(J,I)
          END DO
        END DO
*. Upper half
        IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
          FAC2 = 1.0D0 - FACTOR**2
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
            END DO
          END DO
        ELSE IF(FACTOR .EQ. 1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) = MAT(I,J )
            END DO
          END DO
        ELSE IF(FACTOR .EQ. -1.0D0) THEN
          DO I = 1, NDIM
            DO J = 1, I - 1
              MAT(J,I) =-MAT(I,J )
            END DO
          END DO
        END IF
      ELSE IF(IWAY .EQ. 2 ) THEN
*. Simple blocking of matrix
        LBLK = 40
        NBLK = NDIM/LBLK
        IF(NBLK*LBLK.LT.NDIM) NBLK = NBLK + 1
        IOFF = 1-LBLK
C?      write(6,*) 'NBLK ',nblk
        DO IBLK = 1, NBLK
          IF(IBLK.EQ.-1) write(6,*) 'IBLK = ',IBLK
          IOFF = IOFF + LBLK
          IEND = MIN(IOFF+LBLK-1,NDIM)
          JOFF = 1 - LBLK
          DO JBLK = 1, IBLK
            JOFF = JOFF + LBLK
            JEND = MIN(JOFF+LBLK-1,NDIM)
*. Lower half
            DO  I = IOFF,IEND
              IF(IBLK.EQ.JBLK) JEND = I
              DO J = JOFF,JEND
                MAT(I,J) = MAT(I,J) + FACTOR*MAT(J,I)
              END DO
            END DO
*. Upper half
            IF( ABS(FACTOR) .NE. 1.0D0 ) THEN
              FAC2 = 1.0D0 - FACTOR**2
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = FACTOR*MAT(I,J ) + FAC2 * MAT(J,I)
                 END DO
               END DO
            ELSE IF(FACTOR .EQ. 1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I -1
                DO J = JOFF, JEND
                  MAT(J,I) = MAT(I,J )
                END DO
              END DO
            ELSE IF(FACTOR .EQ. -1.0D0) THEN
              DO I = IOFF, IEND
                IF(IBLK.EQ.JBLK) JEND = I
                DO J = JOFF, JEND
                  MAT(J,I) = - MAT(I,J )
                END DO
              END DO
            END IF
*. ENd of loop over blocks
          END DO
        END DO
*. End of IWAY branching
      END IF
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE TRPMT3(XIN,NROW,NCOL,XOUT)
*
* XOUT(I,J) = XIN(J,I)
*
*. With a few considerations for large scale cases with cache 
*  minimization
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XIN(NROW,NCOL),XOUT(NCOL,NROW)
*
*!!!!!!!!!!!!!!!!!!!!!!!!
      IWAY = 2
*!!!!!!!!!!!!!!!!!!!!!!!!
*
      IF(IWAY.EQ.1) THEN
*. Straightforward, no blocking
        DO IROW =1, NROW
          DO ICOL = 1, NCOL
            XOUT(ICOL,IROW) = XIN(IROW,ICOL)
          END DO
        END DO
      ELSE IF(IWAY.EQ.2) THEN
*. Simple blocking of matrix
        LRBLK = 40
        LCBLK = 40
        NRBLK = NROW/LRBLK
        NCBLK = NCOL/LCBLK
        IF(LRBLK*NRBLK.NE.NROW) NRBLK = NRBLK + 1
        IF(LCBLK*NCBLK.NE.NCOL) NCBLK = NCBLK + 1
*
        DO IRBLK = 1,NRBLK
          IF(IRBLK.EQ.1) THEN
            IROFF = 1
          ELSE
            IROFF = IROFF + LRBLK
          END IF
          IREND = MIN(NROW,IROFF+LRBLK-1)
          DO ICBLK = 1, NCBLK
            IF(ICBLK.EQ.1) THEN
              ICOFF = 1
            ELSE
              ICOFF = ICOFF + LCBLK
            END IF
            ICEND = MIN(NCOL,ICOFF+LCBLK-1)
*
            DO IROW = IROFF,IREND
              DO ICOL = ICOFF,ICEND
                XOUT(ICOL,IROW) = XIN(IROW,ICOL)
              END DO
            END DO
*
          END DO
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE UPPCAS(LINE,LENGTH)
*
* Convert letters in character string LINE to upper case
*
* very stupid and not vectorized !
*
      CHARACTER*(*) LINE
      PARAMETER (NCHAR = 41)
      CHARACTER*1 LOWER(NCHAR)
      CHARACTER*1 UPPER(NCHAR)
*
      DATA LOWER/'a','b','c','d','e',
     &           'f','g','h','i','j',
     &           'k','l','m','n','o',
     &           'p','q','r','s','t',
     &           'u','v','w','x','y',
     &           'z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
      DATA UPPER/'A','B','C','D','E',
     &           'F','G','H','I','J',
     &           'K','L','M','N','O',
     &           'P','Q','R','S','T',
     &           'U','V','W','X','Y',
     &           'Z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
* 
      DO 100 ICHA = 1, LENGTH
        DO 50 I = 1,NCHAR
          IF(LINE(ICHA:ICHA).EQ.LOWER(I)) 
     &    LINE(ICHA:ICHA) = UPPER(I)
   50   CONTINUE
  100 CONTINUE 
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,
     &                       ISCR,NTEST)
*
* construct vertex weights for given supergroup
*
* Reverse lexical ordering is used
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NELFTP(NORBTP),NORBFTP(NORBTP)
*. Ouput
      INTEGER Z(*)
*. Scratch length : 2 * NORB + (NEL+1)*(NORB+1)
      INTEGER ISCR(*)
*
      NORB = IELSUM(NORBFTP,NORBTP)
      NEL  = IELSUM(NELFTP,NORBTP)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Subroutine WEIGHT_SPGP in action '
        WRITE(6,*) ' ================================='
        WRITE(6,*) 'NELFTP '
        CALL IWRTMA(NELFTP,1,NORBTP,1,NORBTP)
      END IF
*
      KLFREE = 1
      KLMAX = KLFREE
      KLFREE = KLFREE + NORB
*
      KLMIN = KLFREE
      KLFREE = KLFREE + NORB
*
      KW = KLFREE
      KLFREE = KW + (NEL+1)*(NORB+1)
*.Max and min arrays for strings
      CALL MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP,
     &                 NTEST)
*. Arc weights
      CALL GRAPW_LUCI(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL,NTEST)
*
      RETURN
      END
***********************************************************************
      SUBROUTINE WRTMT_LU(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
#include "priunit.h"
C
      DO 100 I=1,NROW
        WRITE(lupri,1010) I,(A(I,J),J=1,NCOL)
  100 CONTINUE
 1010 FORMAT(/I4,':',5(1X,E15.8),/,(5X,5(1X,E15.8)))
      RETURN
      END
***********************************************************************
      SUBROUTINE WRTMATMN(A,NROW,NCOL,NMROW,NMCOL,LUWRT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
C
!     XXX = 0.0D0
!     XXX = DDOT(NMCOL,A,1,A,1) <- hjaaj: this cannot be correct /dec 2010
!     IF( XXX .eq. 0.0D0 ) GOTO 200
C
      DO 100 I=1,NROW
        WRITE(LUWRT,1010) I,(A(I,J),J=1,NCOL)
  100 CONTINUE
 1010 FORMAT(/I4,':',5(1X,E15.8),/,(5X,5(1X,E15.8)))
! 200 CONTINUE
      END
***********************************************************************
      SUBROUTINE WRTMATMN8(A,NROW8,NCOL8,NMROW8,NMCOL8,LUWRT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*8 NROW8, NCOL8
      DIMENSION A(NMROW8,NMCOL8)
C
      DO 100 I=1,NROW8
         WRITE(LUWRT,1010) I,(A(I,J),J=1,NCOL8)
  100 CONTINUE
 1010 FORMAT(/I4,':',5E15.8,/,(5X,5E15.8))
      END
***********************************************************************
      SUBROUTINE IWRTMA8(MAT,NROW,NCOL,NMROW,NMCOL)
C
      DIMENSION MAT(NMROW,NMCOL)
#include "priunit.h"
C
      DO 100 I=1,NROW
         WRITE(lupri,1010) I,(MAT(I,J),J=1,NCOL)
  100 CONTINUE
 1010 FORMAT(/I4,':',5I15,/,(5X,5I15))
      END
***********************************************************************
      SUBROUTINE IWRTMAMN(MAT,NROW,NCOL,NMROW,NMCOL,LUWRT)
C
      DIMENSION MAT(NMROW,NMCOL)
C
      DO 100 I=1,NROW
         WRITE(LUWRT,1010) I,(MAT(I,J),J=1,NCOL)
  100 CONTINUE
 1010 FORMAT(/I4,':',5I15,/,(5X,5I15))
      END
***********************************************************************
      SUBROUTINE IWRTMAMN8(MAT,NROW8,NCOL8,NMROW8,NMCOL8,LUWRT)
C
      INTEGER*8 NROW8, NCOL8
      DIMENSION MAT(NMROW8,NMCOL8)
C
      DO 100 I=1,NROW8
         WRITE(LUWRT,1010) I,(MAT(I,J),J=1,NCOL8)
  100 CONTINUE
 1010 FORMAT(/I4,':',5I15,/,(5X,5I15))
      END
***********************************************************************
      SUBROUTINE IWRTMAS(MAT,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAT(NMROW,NMCOL)
      INTEGER ITEST,INACTIVE
#include "priunit.h"

      INACTIVE = 0
C
      DO 100 I=1,NROW
         ITEST = 0
        DO II = 1, NCOL
         ITEST = ITEST + MAT(I,II)
        ENDDO
        IF(ITEST.EQ.0) THEN
         INACTIVE = INACTIVE + 1
        ELSE
         WRITE(lupri,1010) I,(MAT(I,J),J=1,NCOL)
        END IF
  100 CONTINUE
 1010 FORMAT(I4,': ',60I2,/,(6X,60I2))
      WRITE(lupri,'(/1X,A,I4 )')'TOTAL    NUMBER:',NCOL
      WRITE(lupri,'( 1X,A,I4 )')'ACTIVE   BLOCKS:',NCOL-INACTIVE
      WRITE(lupri,'( 1X,A,I4/)')'INACTIVE BLOCKS:',INACTIVE 
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BLOCK_LENGTH(IBLOCKL,N_BLK,ILEN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IBLOCKL(*)
C
      ILEN = IBLOCKL(N_BLK)
C      WRITE(6,*) 'found block length',ILEN,N_BLK
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BLOCK_PROC(IBLOCKD,N_BLK,IPROC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IBLOCKD(*)
C
      IPROC = IBLOCKD(N_BLK)
C      WRITE(6*) 'found proc <--> block ',IPROC,N_BLK
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTTTS(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC)
*
* Print a batch of TTS blocks as given by IBLOCK
*
*
* ISC = 1 : In slater determinant form
* ISC = 2 : In Combination        form
*
*. Jeppe Olsen, August 1995
*
#include "implicit.h"
#include "priunit.h"
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
*
*
      WRITE(lupri,*) '  Batch of blocks '
      WRITE(lupri,*) ' ================= '
      WRITE(lupri,*)
      WRITE(lupri,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF (ISC.EQ.1 ) THEN
          IOFF = IBLOCK(5,JBLOCK)
        ELSE
          IOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(IPACK.EQ.1) THEN
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(lupri,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(lupri,'(A)')
     &      '  ============================'
            CALL PRSM2(BLOCKS(IOFF) ,NIA)
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(lupri,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(lupri,'(A)')
     &      '  ============================'
            CALL WRTMT_LU(BLOCKS(IOFF) ,NIA,NIB,NIA,NIB)
          END IF
        END IF
*
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE ZOOS(ISMOST,IBLTP,MAXSYM,IOCOC,NSSOA,NSSOB,
     &                NOCTPA,NOCTPB,IDC,IOOS,NOOS,NCOMB,IXPND)
*
      implicit real*8 (A-H,O-Z)
*
* Generate offsets for CI vector for RAS CI expansion of given symmetry
* Combination type is defined by IDC
* Total number of combinations NCOMB is also obtained
*
* Symmetry is defined through ISMOST
*
* ICBLTP gives typo of symmetry block :
* = 0 : symmetry block is not included
* = 1 : symmetry block is included , all OO types
* = 2 : symmetry block is included , lower OO types
*
* If IXPND .ne. 0 , the diagonal blocks are always choosen expanded
*
* ========
*  Output
* ========
*
* IOOS(IOCA,IOCB,ISYM) : Start of block with alpha strings of
*                        symmetry ISYM and type IOCA, and
*                        betastrings of type IOCB
* NOOS(IOCA,IOCB,ISYM) : Number of combinations
* The ordering used for the CI vector is
*
*    SYMMETRY  FOR ALPHA STRINGS..(GIVES SYMMETRY OF BETA STRING )
*         OCCUPATION TYPE  FOR ALPHA STRING
*            OCCUPATION TYPE FOR    BETA STRING
*                BETA STRING ( COLUMN INDEX)
*                ALPHA STRINGS ( ROW INDEX )
*    END OF LOOPS
*
*
*. Input
      DIMENSION IOCOC(NOCTPA,NOCTPB),ISMOST(*)
      DIMENSION NSSOA(MAXSYM,NOCTPA),NSSOB(MAXSYM,NOCTPB)
      DIMENSION IBLTP(*)
*. output
      DIMENSION IOOS(NOCTPA,NOCTPB,MAXSYM)
      DIMENSION NOOS(NOCTPA,NOCTPB,MAXSYM)
*
      CALL ISETVC(IOOS,0,NOCTPA*NOCTPB*MAXSYM)
      CALL ISETVC(NOOS,0,NOCTPA*NOCTPB*MAXSYM)
C?    CALL ISETVC(ICBLTP,0,MAXSYM)
      NCOMB = 0
C?    WRITE(6,*) ' ZOOS : IDC  = ', IDC
      DO 100 IASYM = 1, MAXSYM
*
        IBSYM = ISMOST(IASYM)
        IF(IBSYM .EQ. 0 ) GOTO 100
*. Allowed combination symmetry block ?
        IF(IDC.NE.1.AND.IBLTP(IASYM).EQ.0) GOTO 100
*. Allowed occupation combinations
        DO  95 IAOCC = 1,NOCTPA
          IF(IBLTP(IASYM).EQ.1) THEN
            MXBOCC = NOCTPB
            IREST1 = 0
          ELSE
            MXBOCC = IAOCC
            IREST1 = 1
          END IF
          DO 90 IBOCC = 1, MXBOCC
*.Is this block allowed
            IF(IOCOC(IAOCC,IBOCC).EQ.1) THEN
              IOOS(IAOCC,IBOCC,IASYM) = NCOMB+1
              IF(IXPND.EQ.0 .AND. IREST1.EQ.1 .AND. IAOCC.EQ.IBOCC)THEN
                NCOMB = NCOMB
     &      +   (NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
                NOOS(IAOCC,IBOCC,IASYM) =
     &          (NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
              ELSE
                NCOMB = NCOMB
     &      +   NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
                NOOS(IAOCC,IBOCC,IASYM) =
     &          NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
              END IF
            END IF
C?      write(6,*) ' NOOS(IA,IB,ISM) ',NOOS(IAOCC,IBOCC,IASYM)
   90     CONTINUE
   95   CONTINUE
  100 CONTINUE
*
      NTEST = 0 
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) 
         WRITE(6,*) ' ==============='
         WRITE(6,*) ' ZOOS reporting '
         WRITE(6,*) ' ==============='
         WRITE(6,*) 
         WRITE(6,*) ' NSSOA, NSSOB ( input ) '
         CALL IWRTMA(NSSOA,MAXSYM,NOCTPA,MAXSYM,NOCTPA)
         CALL IWRTMA(NSSOB,MAXSYM,NOCTPB,MAXSYM,NOCTPB)
         WRITE(6,*) 
         WRITE(6,*) ' Number of combinations obtained ',NCOMB
         WRITE(6,*) ' Offsets for combination OOS blocks '
         DO 111 IASYM = 1,MAXSYM
           WRITE(6,'(A,I2)') '  Symmetry ',IASYM
           CALL IWRTMA(IOOS(1,1,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  111    CONTINUE
         WRITE(6,*) ' Number of  combinations per OOS blocks '
         DO 112 IASYM = 1,MAXSYM
           WRITE(6,'(A,I2)') '  Symmetry ',IASYM
           CALL IWRTMA(NOOS(1,1,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  112    CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE ZSPGPIB(NSTSO,ISTSO,NSPGP,NSMST)
*
*  Offset for supergroups of strings with given sym.
*. Each supergroup is given start address 1
*
* Jeppe Olsen, Still summer of 95
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER NSTSO(NSMST,NSPGP)
*. Output
      INTEGER ISTSO(NSMST,NSPGP)
*
      DO ISPGP = 1, NSPGP
        ISTSO(1,ISPGP) = 1
        DO ISMST = 2, NSMST
          ISTSO(ISMST,ISPGP) = ISTSO(ISMST-1,ISPGP) + NSTSO(ISMST,ISPGP)
        END DO
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ZSPGPIB '
        WRITE(6,*) ' =================== '
        WRITE(6,*)
        CALL IWRTMA(ISTSO,NSMST,NSPGP,NSMST,NSPGP)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                     NSMST,MXLNG,IOCOC,ISMOST,
     &                     NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
     &                     ITTSS_ORD)
*
* Partition a CI vector into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete CI vector is constructed
*
*
* Compared to PART_CIV, the NOOS arrays have been eliminated.
* They are becoming the size defining arrays - atleast at
* the laptop
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
C     INTEGER NOOS(NOCTPA,NOCTPB,NSMST)
C     INTEGER NOOSP(NOCTPA,NOCTPB,NSMST)
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '     PART_CIV2     '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
*
*. block zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV2'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV2. ' )
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches ', NBATCH
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, by Timo Fleig and Jeppe Olsen                             *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PART_CIV2_SPC(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD)
*
* Partition a CI vector into batches of blocks.
* The length of a batch must be atmost MXLNG
*
* IF ICOMP. eq. 1 the complete CI vector is constructed
*
*
* Compared to PART_CIV, the NOOS arrays have been eliminated.
* They are becoming the size defining arrays - atleast at
* the laptop
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, August 1995
*
*
* IMPORTANT change compared to PART_CIV2:
*
* MXLNG is no longer assumed to be equals to the total number
* of dets. Therefore we have to partition the CI vector being aware
* of putting symmetry blocks of same TT into different batches.
*
*
* REASON: SIGMA-vector / DENSITY matrix calculation is done wrt loops
* over NSMST!
*
* S. Knecht - July 02 2007
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
C     INTEGER NOOS(NOCTPA,NOCTPB,NSMST)
C     INTEGER NOOSP(NOCTPA,NOCTPB,NSMST)
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      I_SAVE_IBLOCK = 0
      N_TRY = 0
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) '   PART_CIV2_SPC   '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
*
*. block zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV2_SPC'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV2_SPC. ' )
      ELSE
*
*     this batch is finished, goto next batch
*
*       is this the TT-1 block?
        IF( ISM .ne. 1 )THEN
*
*         go back to ISM == 1 since we loop over NSMST in GNSIDE_REL
*
          N_TRY = N_TRY + 1
*
          IF( N_TRY .gt. 1 .and. I_SAVE_IBLOCK .eq. IBLOCK ) THEN
            WRITE(6,*) ' Not enough scratch space to include a '//
     &                     ' single Block'
            WRITE(6,*) ' Since I cannot procede I will stop '
            WRITE(6,*) ' Insufficient buffer detected in '//
     &                     ' PART_CIV2_SPC'
            WRITE(6,*) '  LENGTH, LBLOCK ',LENGTH
            WRITE(6,*) ' Alter GAS space of raise Buffer from',MXLNG
            call abend2(' In PART_CIV2_SPC because of N_TRY .ge. 1. ')
          END IF
*
          I_SAVE_IBLOCK = IBLOCK
*         set batch back to first available TT-1 block
          call bck_tts(ITTSS_ORD,ISM,LBATCH,LEBATCH,IBATCH,IBLOCK,
     &                 NBLOCK,LENGTH,LENGTHP,NBATCH,6)
*
          GOTO 2000
*
        ELSE
          N_TRY = 0
          I_SAVE_IBLOCK = 0
          GOTO 2000
        END IF
*       ^ ISM check
      END IF
*     ^ length of batch < MXLNG?
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches ', NBATCH
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
      SUBROUTINE DIAVC2_TRUNC(VECOUT,VECIN,DIAG,SHIFT,NDIM,
     &                        THR_RTRUNC,THR_ETRUNC)
C
C VECOUT(I)=VECIN(I)/(DIAG(I)+SHIFT)
C
C Set VECOUT = 0 if VECIN = 0
C   Timo Fleig, 18-12-2002
C Implemented THR_RTRUNC and THR_ETRUNC /HJAaJ 28-Mar-2007
C
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (THRES=1.0D-4)
      DIMENSION VECOUT(*),VECIN(*),DIAG(*)
C
      n_zero_1 = 0
      n_zero_2 = 0
      DO I=1,NDIM
         IF(ABS(VECIN(I)) .LE. THR_RTRUNC) THEN
           VECOUT(I) = 0.0D0
           n_zero_1 = n_zero_1 + 1
         ELSE
           DIVIDE=DIAG(I)+SHIFT
           IF(ABS(DIVIDE).LE.THRES) DIVIDE=THRES
           VECOUT(I)=VECIN(I)/DIVIDE
           IF (THR_ETRUNC .GT. 0.0D0) THEN
              IF (ABS(VECOUT(I)*VECIN(I)) .LE. THR_ETRUNC) THEN
                 VECOUT(I) = 0.0D0
                 n_zero_2 = n_zero_2 + 1
              END IF
           END IF
         END IF
      END DO
      x_nonzero = ndim-n_zero_1-n_zero_2
      x_nonzero = x_nonzero / (0.01D0*ndim)
!     write (lupri,'(a,i3,a,4i10,f7.2,a)')
!    & 'Node ',mytid,' DIAVC2_TRUNC, NDIM,zero1/2,nonzero =',
!    & NDIM,n_zero_1, n_zero_2, ndim-n_zero_1-n_zero_2,x_nonzero,'%'
      END
