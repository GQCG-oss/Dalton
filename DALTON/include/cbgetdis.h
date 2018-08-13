!     -----------------------------------------------------------------
!     Purpose of /CBGETD/: transfer information to GETIN2 about
!     how the "H2AC" integrals are available.  The distribution
!     type DISTYP tells about packing and index symmetry.
!     IADINT .lt. 0 means integrals are in-core.
!     IADINT .ge. 0 is offset for reading an active-active distribution
!     from disk.  IADH2, IADH2X, IADH2D are used for off-sets for
!     H2AC, H2XAC, and H2DAC, resp. (H2DAC from ABACUS).
      INTEGER         DISTYP, IADINT,IADH2,IADH2X,IADH2D
      COMMON /CBGETD/ DISTYP, IADINT,IADH2,IADH2X,IADH2D
!     -----------------------------------------------------------------
