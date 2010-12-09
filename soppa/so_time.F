#include <single.h>
C
C  /* Deck so_time */
      SUBROUTINE SO_TIME(TIMTOT,TIMWTO,PARRAY,LPARRA)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, February 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Sort and print information about timings in SOPPA
C              calculation.
C
#include <implicit.h>
#include <priunit.h>
C
#include <soppinf.h>
C
C
      CHARACTER*14 SUBNM(LSOTIM)
      DATA SUBNM  /' LAMMAT       ',' SO_INIT      ',' SO_MOENERGY  ',
     &             ' SO_DENS      ',' SO_T2MP      ',' SO_AODENS    ',
     &             ' SO_BCKTR     ',' ERIDI1       ',' ERIDI2       ',
     &             ' CCRDAO       ',' CC_AOFOCK    ',' SO_T2M1      ',
     &             ' CCTRBT       ',' SO_RES_A     ',' SO_RES_CB    ',
     &             ' SO_SIGAI     ',' SO_SIGDA     ',' SO_RES_B     ',
     &             ' SO_RES_C     ',' SO_RES_SYM   ',' SO_RES_FCK   ',
     &             ' SO_RES_OVL   ',' SO_RES_OVLR  ',' CC_FCKMO     ',
     &             ' SO_TWOFOCK   ',' SO_ONEFOCK   ',' SO_CHECK     ',
     &             ' SO_INCRED    ',' SO_RES_TCB   ',' SO_RES_CD    ',
     &             ' SO_DIAG      ',' SO_TRIAL1    ',' SO_REDEV     ',
     &             ' SO_TRIAL2    ',' SO_ERES      ',' SO_ORDEIG    ',
     &             ' SO_OPTVEC    ',' SO_ROPT      ',' SO_ORTH_TRN  ',
     &             ' SO_SRES      ',' SO_DENSAI    ',' ???????????  ',
     &             ' ?????????    '/
      CHARACTER*14 SUBNMW(LSOWTI)
      DATA SUBNMW /' SO_ERES      ',' SO_SRES      ',' SO_REDEV     ',
     &             ' SO_TRIAL2    '/
      CHARACTER*14 RWOCI(LORWCI)
      DATA RWOCI  /' SO_READ      ',' SO_WRITE     ',' SO_OPEN      ',
     &             ' SO_CLOSE     '/
C
      INTEGER      PARRAY(LPARRA), SOHOUR, SOMIN, SOSEC
      CHARACTER*56 LINE1,LINE2
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TIME')
C
C---------------------------------------------------------------
C     Subtract from time in SO_RES all time spend in subroutines
C     which are timed separately.
C---------------------------------------------------------------
C
      SOTIME(35) = SOTIME(35) 
     &           - SOTIME(1)  - SOTIME(6)  - SOTIME(7)  - SOTIME(8)
     &           - SOTIME(9)  - SOTIME(10) - SOTIME(11) - SOTIME(12)
     &           - SOTIME(13) - SOTIME(14) - SOTIME(15) - SOTIME(16)
     &           - SOTIME(17) - SOTIME(18) - SOTIME(19) - SOTIME(20)
     &           - SOTIME(21) - SOTIME(24)
     &           - SOTIME(25) - SOTIME(26) - SOTIME(29) - SOTIME(30)
     &           - SOTIME(41)
C
C-----------------------------------------------------------------
C     Subtract from time in SO_REDEV all time spend in subroutines
C     which are timed separately.
C-----------------------------------------------------------------
C
      SOTIME(33) = SOTIME(33)
     &           - SOTIME(28) - SOTIME(36) - SOTIME(37)
C
C------------------------------------------------------------------
C     Subtract from time in SO_TRIAL2 time spend in SO_ORTH_TRO and
C     SO_ORTH_TRN.
C------------------------------------------------------------------
C
      SOTIME(34) = SOTIME(34) - SOTIME(38) - SOTIME(39)
C
C-----------------------------------------------------------
C     Subtract from time in SO_SIGAI time spend in subparts.
C-----------------------------------------------------------
C
      SOTIME(16) = SOTIME(16) 
     &           - SOTIME(40) - SOTIME(41) - SOTIME(42) - SOTIME(43)
C
C------------------------------
C     Initialize pointer-array.
C------------------------------
C
      DO 100 I = 1,LSOTIM
         PARRAY(I) = I
  100 CONTINUE
C
C--------------------------------------------------------------------
C     Arrange pointer-array to have first pointer to point at largest
C     time, second pointer to second largest and so on.
C--------------------------------------------------------------------
C
      DO 200 J = 1,LSOTIM
C
         DO 300 I = 1,LSOTIM-J
C
            IF ( SOTIME(PARRAY(I+1)) .GT. SOTIME(PARRAY(I)) ) THEN
C
               PI          = PARRAY(I+1)
               PARRAY(I+1) = PARRAY(I)
               PARRAY(I)   = PI
C
            ENDIF
C
  300    CONTINUE
C
  200 CONTINUE
C
      DTIME  = SECOND()
      TIMTOT = DTIME - TIMTOT
      CALL GETTIM (DUMMY,WTIME)
      TIMWTO = WTIME - TIMWTO
C
      LINE1 = '|=====================================================|'
      LINE2 = '|-----------------------------------------------------|'
C
      WRITE(LUPRI,'(1X)')
      WRITE(LUPRI,'(1X,A)') LINE1
      WRITE(LUPRI,'(1X,A)')
     &   '|           Timings for AO-SOPPA subroutines          |'
      WRITE(LUPRI,'(1X,A)') LINE1
      WRITE(LUPRI,'(1X,A,A)') '|   Routine            cpu time   ',
     &                     'of total        sum |'
      WRITE(LUPRI,'(1X,A)') LINE2
      TPRCEN =0.D0
C
      DO 400 I = 1,LSOTIM
C
C--------------------------------------------------------------
C        Calculate hours, minutes, and seconds used in various
C        subroutines. Also calculate the time in % of total and
C        the running sum in %.
C--------------------------------------------------------------
C
         SOHOUR = INT ( SOTIME(PARRAY(I)) / 3600.D0 )
         SOMIN  = INT ( MOD ( SOTIME(PARRAY(I)), 3600.D0 ) / 60.D0 )
         SOSEC  = NINT ( MOD ( MOD ( SOTIME(PARRAY(I)), 3600.D0 ),
     &                   60.D0 ) )
         PRCENT = ( SOTIME(PARRAY(I)) / TIMTOT ) * 100
         TPRCEN = TPRCEN + PRCENT
C
C-------------------------------------
C        Print the calculated timings.
C-------------------------------------
C
         WRITE(LUPRI,
     &         '(1X,A3,A14,3X,I5,A1,I2,A1,I2,5X,F4.1,A2,4X,F5.1,A4)')
     &         '|  ',SUBNM(PARRAY(I)),SOHOUR,':',SOMIN,':',SOSEC,PRCENT,
     &         ' %',TPRCEN,' % |'
C
  400 CONTINUE
C
C-----------------------------------------------------------------
C     Calculate and print information about total timing in SOPPA.
C-----------------------------------------------------------------
C
      SOHOUR = INT ( TIMTOT / 3600.D0 )
      SOMIN  = INT ( MOD ( TIMTOT, 3600.D0 ) / 60.D0 )
      SOSEC  = NINT ( MOD ( MOD ( TIMTOT, 3600.D0 ), 60.D0 ) )
C
      WRITE(LUPRI,'(1X,A)') LINE2
      WRITE(LUPRI,'(1X,A19,1X,I5,A1,I2,A1,I2,23X,A1)')
     &       '|   Total CPU time ',SOHOUR,':',SOMIN,':',SOSEC,'|'
      WRITE(LUPRI,'(1X,A)') LINE1
C
C-----------------------------------------------------
C     Print separate timing for input output in SOPPA.
C-----------------------------------------------------
C
      WRITE(LUPRI,'(/,1X,A)') LINE2
C
      DO 500 I = 1,LORWCI
C
         SOHOUR = INT ( SOORWC(I) / 3600.D0 )
         SOMIN  = INT ( MOD ( SOORWC(I), 3600.D0 ) / 60.D0 )
         SOSEC  = NINT ( MOD ( MOD ( SOORWC(I), 3600.D0 ), 60.D0 ) )
         PRCENT = ( SOORWC(I) / TIMTOT ) * 100
C
         WRITE(LUPRI,'(1X,A3,A14,3X,I5,A1,I2,A1,I2,5X,F4.1,A2,12X,A1)')
     &         '|  ',RWOCI(I),SOHOUR,':',SOMIN,':',SOSEC,PRCENT,
     &         ' %','|'
C
  500 CONTINUE
C
      WRITE(LUPRI,'(1X,A)') LINE2
C
C-----------------------
C     Print wall timing.
C-----------------------
C
      WRITE(LUPRI,'(1X)')
      WRITE(LUPRI,'(1X,A)') LINE1
      WRITE(LUPRI,'(1X,A)')
     &   '|        Wall timings for AO-SOPPA subroutines        |'
      WRITE(LUPRI,'(1X,A)') LINE1
      WRITE(LUPRI,'(1X,A,A)') '|   Routine           wall time   ',
     &                     'of total       sum  |'
      WRITE(LUPRI,'(1X,A)') LINE2
      TPRCEN =0.D0
C
      DO 600 I = 1,LSOWTI
C
         SOHOUR = INT ( SOWTIM(I) / 3600.D0 )
         SOMIN  = INT ( MOD ( SOWTIM(I), 3600.D0 ) / 60.D0 )
         SOSEC  = NINT ( MOD ( MOD ( SOWTIM(I), 3600.D0 ), 60.D0 ) )
         PRCENT = ( SOWTIM(I) / TIMWTO ) * 100
         TPRCEN = TPRCEN + PRCENT
C
         WRITE(LUPRI,
     &         '(1X,A3,A14,3X,I5,A1,I2,A1,I2,5X,F4.1,A2,4X,F5.1,A4)')
     &         '|  ',SUBNMW(I),SOHOUR,':',SOMIN,':',SOSEC,PRCENT,
     &         ' %',TPRCEN,' % |'
C
  600 CONTINUE
C
C----------------------------------------------------------------------
C     Calculate and print information about total wall timing in SOPPA.
C----------------------------------------------------------------------
C
      SOHOUR = INT ( TIMWTO / 3600.D0 )
      SOMIN  = INT ( MOD ( TIMWTO, 3600.D0 ) / 60.D0 )
      SOSEC  = NINT ( MOD ( MOD ( TIMWTO, 3600.D0 ), 60.D0 ) )
C
      WRITE(LUPRI,'(1X,A)') LINE2
      WRITE(LUPRI,'(1X,A19,1X,I5,A1,I2,A1,I2,23X,A1)')
     &      '|   Total wall time',SOHOUR,':',SOMIN,':',SOSEC,'|'
      WRITE(LUPRI,'(1X,A)') LINE1
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TIME')
C
C
      RETURN
      END
