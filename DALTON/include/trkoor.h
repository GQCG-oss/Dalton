! FILE : trkoor.h
C
C Those marked "%" are not used in current Dalton /Oct.2013 hjaaj
C
C     NCOOR  - total number of Cartesian coordinates
C %   NCOORS - number of Cartesian coordinates of each symmetry
C     DEPEND - true for dependent symmetry coordinates
C %   NCDEP  - total number of dependent/trarot coordinates
C %   NCDEPS - number of dependent/trarot coordinates of each symmetry
C %   NCIND
C %   NCINDS
C
C     IPTTRO - identifies symmetry-ordered trarot coordinate as
C              Tx, Ty, Tz, Rx, Ry, or Rz
C     NTRREP - number of trarot coordinates in each symmetry
C     NPRREP - number of trarot coordinates to project out in each symmetry
C
      INTEGER NCOOR, NCOORS, NCDEP, NCDEPS, NCIND, NCINDS
      INTEGER IPTTRO, NTRREP, NPRREP
      LOGICAL DEPEND
      COMMON /TRKOOR/ NCOOR, NCOORS(8),
     &                NCDEP, NCDEPS(8), NCIND, NCINDS(8),
     &                IPTTRO(6,8), NTRREP(0:7), NPRREP(0:7),
     &                DEPEND(MXCOOR)
! -- end of trkoor.h --
