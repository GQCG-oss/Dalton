
      integer :: MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV, ICISTR
      integer :: NOCSF,IDIAG,NOINT,IDMPIN,MXINKA,ICJKAIB, INIREF,IRESTRF
      integer :: MOCAA, MOCAB,IPERT,NPERT, IAPRREF,IAPRZER,IIDUM
      integer :: NSEQCI, ISEQCI, IEXTKOP, IE0AVEX,IC1DSC, IH0SPC,NPTSPC
      integer :: IPTSPC, IOCPTSPC, IH0INSPC, IRFROOT
      integer :: INIDEG, LCSBLK,IFINMO,NPSSPC, NPSSH, ICLSSEL,IDENSI
      integer :: IPTEKT, NPTEKT,LPTEKT, IH0ROOT,IRST2,IPAD, ispnden
      integer :: IXYZSYM, NPROP, ITRAPRP,NEXCSTATE,IEXCSYM,IRESPONS 
      integer :: NRESP, MAXORD_OP, MXITLE,N_AVE_OP,IROOTHOMING,IUSE_PH
      integer :: IADVICE,ITRACI, IUSE_PA,IPTFOCK,ITRANS ,NSXE,NDXEE 
      integer :: N_CC_AMP, I_RESTRT_CC,LAST_CC_SPC,LAST_CC_RUN
      integer :: IRELAX, LBLOCK, L0BLOCK, lsingle_resolution_block
      real(8) :: THRES_E, THRES_E_aux, E_THRE,C_THRE, E_CONV,C_CONV
      real(8) :: E0READ, RESP_W, XLAMBDA
      COMMON/CRUN/MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,       &
     &            ICISTR,NOCSF,IDIAG,NOINT,IDMPIN,MXINKA,ICJKAIB,       &
     &            INIREF,IRESTRF,MOCAA,MOCAB,IPERT,NPERT,               &
     &            IAPRREF,IAPRZER,IIDUM,                                &
     &            NSEQCI(10),ISEQCI(10,10),IEXTKOP,                     &
     &            IE0AVEX,IC1DSC,IH0SPC,NPTSPC,IPTSPC(MXPPTSPC),        &
     &            IOCPTSPC(2,MXPNGAS,MXPPTSPC),IH0INSPC(MXPPTSPC),      &
     &            IRFROOT,                            INIDEG,           &
     &            LCSBLK,IFINMO,NPSSPC,                                 &
     &            NPSSH(MXPIRR,MXPNGAS),ICLSSEL,IDENSI,IPTEKT,          &
     &            NPTEKT,LPTEKT(20,2),IH0ROOT,IRST2,        IPAD,       &
     &            IXYZSYM(3),NPROP,ispnden,                             &
     &            ITRAPRP,NEXCSTATE,IEXCSYM,IRESPONS,NRESP,             &
     &            MAXORD_OP(2,20), lsingle_resolution_block,            &
     &            MXITLE,N_AVE_OP,IROOTHOMING,IUSE_PH,IADVICE,ITRACI ,  &
     &            IUSE_PA,IPTFOCK,ITRANS ,NSXE,NDXEE,N_CC_AMP,          &
     &            I_RESTRT_CC,LAST_CC_SPC,LAST_CC_RUN,                  &
     &            IRELAX,LBLOCK,L0BLOCK,THRES_E, XLAMBDA, E_THRE,C_THRE,&
     &            E_CONV,C_CONV, THRES_E_aux, E0READ, RESP_W(20)
      CHARACTER*6 PROPER,ENVIRO
      CHARACTER*2 PADDY
      CHARACTER*8 RESP_OP, CSEQCI, AVE_OP, ITRACI_CR,ITRACI_CN
      COMMON /CRUNC/ PROPER(20) ,ENVIRO, PADDY,                         &
     &            RESP_OP(2,20),CSEQCI(10,10), AVE_OP(20),              &
     &            ITRACI_CR,ITRACI_CN
