*fordeck SysCtl.INC $Revision: 1.1.2.1 $
*----------------------------------------------------------------------*
*     General definitions                                              *
*----------------------------------------------------------------------*
      Integer ON,OFF
CStart Molpro
cStart cray 1-i8
c;      Data ON/8hON      /,OFF/8hOFF     /
cElse
c;      Data ON/4hON      /,OFF/4hOFF     /
cEnd
CElse

      Data ON/8hON      /,OFF/8hOFF     /

CEnd
*
      Parameter ( ipStat   = 1          )
      Parameter ( ipTrace  = ipStat+1   )
      Parameter ( ipQuery  = ipTrace+1  )
      Parameter ( ipSysOut = ipQuery+1  )
*----------------------------------------------------------------------*
*     Definitions specific to the subroutines SetMem and GetMem        *
*     and IniMem                                                       *
*----------------------------------------------------------------------*
      Parameter ( mxFld    = 1024       )
*
      Parameter ( ipClear  = ipSysOut+1 )
      Parameter ( ipCheck  = ipClear+1  )
      Parameter ( ipFldCnt = ipCheck+1  )
      Parameter ( ipFldDat = ipFldCnt+1 )
      Parameter ( lMemCtl  = ipFldDat+6*mxFld)
*
      Integer MemCtl(lMemCtl)
      Common / MOLCAS_GetMem / MemCtl
*----------------------------------------------------------------------*
*     Definitions specific to the subroutines ErrMsg, SysMsg,          *
*     IniMsg and SetMsg                                                *
*----------------------------------------------------------------------*
*
      Integer LuErrDB,LuSysDB
      Character*8 FnErrDB,FnSysDB
      Data LuErrDB / 90 / , FnErrDB / 'ERRDB   ' /
      Data LuSysDB / 91 / , FnSysDB / 'SYSDB   ' /
*
      Integer Itrap,Wtrap,Etrap,Strap,Xtrap
CStart Molpro
cStart cray 1-i8
c;      Data Itrap /8hI       / , Wtrap /8hW       / , Etrap /8hE       /
c;      Data Strap /8hS       / , Xtrap /8hX       /
cElse
c;      Data Itrap /4hI       / , Wtrap /4hW       / , Etrap /4hE       /
c;      Data Strap /4hS       / , Xtrap /4hX       /
cEnd
CElse

      Data Itrap /8hI       / , Wtrap /8hW       / , Etrap /8hE       /
      Data Strap /8hS       / , Xtrap /8hX       /

CEnd
*
      Parameter ( ipAtrI   = ipSysOut+1 )
      Parameter ( ipAtrW   = ipAtrI+1   )
      Parameter ( ipAtrE   = ipAtrW+1   )
      Parameter ( ipAtrS   = ipAtrE+1   )
      Parameter ( lErrCtl  = ipAtrS+1   )
*
      Integer ErrCtl(lErrCtl)
      Common / MOLCAS_ErrHnd / ErrCtl
*----------------------------------------------------------------------*
*     Definitions specific to the subroutines qEnter, qExit, qStat,    *
*     IniQue and SetQue                                                *
*----------------------------------------------------------------------*
      Parameter ( mxQTbl   = 256        )
      Parameter ( mxQStk   = 128        )
*
      Parameter ( iplTbl   = ipSysOut+1 )
      Parameter ( ipTbl    = iplTbl  +1 )
      Parameter ( iplStk   = ipTbl   + 4*mxQTbl )
      Parameter ( ipStk    = iplStk  +1 )
      Parameter ( lQueCtl  = ipStk   + 3*mxQStk )
*
      Integer QueCtl(lQueCtl)
      Common / MOLCAS_Query  / QueCtl
*----------------------------------------------------------------------*
*     Definitions system buffer which is referenced by the subroutines *
*     Remark, StkI, StkR, StkA and DmpMsg                              *
*----------------------------------------------------------------------*
      Parameter ( iZone1 = 1 )
      Parameter ( iZone2 = 2048 )
*
      Parameter ( lSysBuf = 128*128  )
      Parameter ( lHeuer  = 8        )
*
      Integer Heuer(lHeuer)
      Integer StkIni
      Character*1 SysBuf(lSysBuf)
      Common / MOLCAS_SysBuf / Heuer,StkIni,SysBuf
