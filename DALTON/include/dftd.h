!AMT for now a dirty comon block to pass parameters
!AMT this one needs MXCOOR defined...
      LOGICAL DO_DFTD2
      LOGICAL DO_DFTD3
      LOGICAL DO_BJDAMP
      LOGICAL DO_3BODY

      LOGICAL L_INP_D2PAR
      LOGICAL L_INP_D3PAR

!Initialized in DFTINI in dft_aux.F
      COMMON /DFTDLOG/ DO_DFTD2, DO_DFTD3, DO_BJDAMP, DO_3BODY,         &
     &                  L_INP_D2PAR, L_INP_D3PAR

      DOUBLE PRECISION CN
      COMMON /DFTD3/ CN(MXCOOR)

      DOUBLE PRECISION D2_s6_inp, D2_alp_inp, D2_rs6_inp
      DOUBLE PRECISION D3_s6_inp, D3_alp_inp, D3_rs6_inp
      DOUBLE PRECISION D3_rs18_inp, D3_s18_inp

      COMMON /DFTDINP/ D2_s6_inp, D2_alp_inp, D2_rs6_inp,               &
     &                  D3_s6_inp, D3_alp_inp, D3_rs6_inp,              &
     &                  D3_rs18_inp, D3_s18_inp


