!
! VERSION : $Revision: 1.1.2.1 $
! DATE    : $Date: 2009-11-13 10:35:26 $
! FILE    : parluci.h
!
!                     PARALLEL INFORMATION
!
!     common blocks which are general to LUCITA and LUCIAREL
!     ======================================================
!
C
C     definitions, communication handles and group information
C
      INTEGER    MSLVOUT_REL
      PARAMETER (MSLVOUT_REL = 7)
C
      INTEGER         LUCI_MASTER,LUCI_MYPROC,LUCI_NMPROC,
     &                MY_GROUPN,MYNEW_COMM,MYNEW_ID,NEWCOMM_PROC,
     &                JCOMCOM,ICOMM,ICOMM_ID,ICOMM_SIZE,
     &                IBLOCKCOMM,N_MASTER,
     &                MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,
     &                N_MASTER_SM, MYNEW_ID_SM_C, MYNEW_COMM_SM_C,
     &                NEWCOMM_PROC_SM_C, N_MASTER_SM_C, IXCOMM, 
     &                IXCOMM_SZ,IXCOMM_RK, NUM_BLOCKS, NUM_BLOCKS2,
     &                NPTEST_VAR, IRC_SAVE, LUCIWT, JVEC_SF,
     &                LBLOCK_LUCI, IAM_NOT_INV, IRUNPA, IODENSEPAR, 
     &                NCBLOCKS_USED, IIOMOD, NFLGRPS, NROOT_INFO, 
     &                IDISTROUTE, IPRINTSTAT, LUWRT, IIOMOD_REL, 
     &                NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE, L_COMBI, 
     &                I_NZERO_LEN, ISHARE_DENSM, IPARALLELIO, 
     &                LF2_ZERO, ISMEMFAC, IIJKL_ROD_LL, LEVEL_SM, 
     &                IT_SHL, IC_SHL, ICOUNTABLK_C, ISI_CALC_BL,
     &                IBI_MULT_BL, LMEMFREE_PTR
C
      COMMON/LUCIPAR/ LUCI_MASTER,LUCI_MYPROC,LUCI_NMPROC,
     &                MY_GROUPN,MYNEW_COMM,MYNEW_ID,NEWCOMM_PROC,
     &                JCOMCOM,ICOMM,ICOMM_ID,ICOMM_SIZE,
     &                IBLOCKCOMM,N_MASTER,
     &                MYNEW_ID_SM, MYNEW_COMM_SM, NEWCOMM_PROC_SM,
     &                N_MASTER_SM, MYNEW_ID_SM_C, MYNEW_COMM_SM_C,
     &                NEWCOMM_PROC_SM_C, N_MASTER_SM_C, IXCOMM, 
     &                IXCOMM_SZ,IXCOMM_RK, NUM_BLOCKS, NUM_BLOCKS2,
     &                NPTEST_VAR, IRC_SAVE, LUCIWT, JVEC_SF,
     &                LBLOCK_LUCI, IAM_NOT_INV, IRUNPA, IODENSEPAR, 
     &                NCBLOCKS_USED, IIOMOD, NFLGRPS, NROOT_INFO, 
     &                IDISTROUTE, IPRINTSTAT, LUWRT, IIOMOD_REL, 
     &                NFLGRPS_REL, NROOT_SAVE, NVEC_SAVE, L_COMBI, 
     &                I_NZERO_LEN, ISHARE_DENSM, IPARALLELIO, 
     &                LF2_ZERO, ISMEMFAC, IIJKL_ROD_LL, LEVEL_SM, 
     &                IT_SHL, IC_SHL, ICOUNTABLK_C, ISI_CALC_BL,
     &                IBI_MULT_BL, LMEMFREE_PTR
C
C     MPI file lists and MPI file handles
C
      INTEGER IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4, IALL_LU5,
     &        IALL_LU6, IALL_LU7, IALL_LUC, IIJKL_ROD
      COMMON/LUCIAPFILE/ IDIA,ILU1,ILU2,ILU3,ILU4,ILU5,ILU6,ILU7,
     &                   ILUC,
     &                   IALL_LU1, IALL_LU2, IALL_LU3, IALL_LU4,
     &                   IALL_LU5, IALL_LU6, IALL_LU7, IALL_LUC, 
     &                   IIJKL_ROD
#if defined (VAR_MPI)
      INTEGER MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL
      INTEGER*8 MY_VEC1_IOFF, MY_VEC2_IOFF, MY_DIA_OFF, MY_LU1_OFF,
     &          MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF, MY_LU5_OFF,
     &          MY_LU6_OFF, MY_LU7_OFF, MY_LUC_OFF
      COMMON/LUCIFILEOFF/ MY_VEC1_IOFF, MY_VEC2_IOFF,
     &                    MY_DIA_OFF, MY_LU1_OFF,
     &                    MY_LU2_OFF, MY_LU3_OFF, MY_LU4_OFF,
     &                    MY_LU5_OFF, MY_LU6_OFF, MY_LU7_OFF, 
     &                    MY_LUC_OFF,
     &                    MY_ACT_BLK1, MY_ACT_BLK2, MY_ACT_BLK_ALL
#endif
!     integer*8 block
      INTEGER*8           LEN_ALL_INT, IS_LENGTH_TT, NINT_TOTAL,
     &                    LEN_T_BUFF_NZ, MY_T_LEN, MY_T_LEN_D
      COMMON/LUCIPAR_I8/  LEN_ALL_INT, IS_LENGTH_TT, NINT_TOTAL,
     &                    LEN_T_BUFF_NZ, MY_T_LEN, MY_T_LEN_D
!     double precision block
      DOUBLE PRECISION  RNORM_FAC, TRUNC_FAC, mem_cp_time, xreadtimebi,
     &                  xcomputesi
      COMMON/LUCI_REAL/ RNORM_FAC, TRUNC_FAC, mem_cp_time, xreadtimebi,
     &                  xcomputesi
!     logical block
      LOGICAL CSCREEN, SHARED_M, REORD_IJKL, CHECK_TC, USE_EX_IDIA
!
      COMMON/LUCI_LOG/
     &        CSCREEN, SHARED_M, REORD_IJKL, CHECK_TC, USE_EX_IDIA
!     character block
      CHARACTER*9 NIIJKL_ROD
      COMMON/LUCI_CHAR/ NIIJKL_ROD
