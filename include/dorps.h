C
C     Parameters NDORL must be updated after changes (for parallelization)
C
C     NOTE: There must be no variables between DOREPS and DOCOOR.
C
      PARAMETER (NDORL = 8 + 3*MXCENT)
      LOGICAL DOREPS, DOCOOR, DCORD, DCORGD, DOPERT
      COMMON /DORPS/ DOREPS(0:7), DOCOOR(3,MXCENT), NDCORD(2),
     &               DCORD(MXCENT,3,2), DCORGD(0:MXCENT,3,2),
     &               DOPERT(0:3*MXCENT,2)
