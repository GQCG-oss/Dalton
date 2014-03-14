C
C
C
C     SLESLA = Sleeping Slaves
C     control module for the subroutines
C     sleepslaves, sleepslaves2 and wakeslaves
C     to make sure that repeated calls to sleepslaves do not broadcast
C     uncessessary messages from the master to the slaves.



      logical:: awake=.true.
      common /DOSLEEP/ awake
