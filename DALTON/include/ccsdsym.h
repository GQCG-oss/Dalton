      INTEGER NCKIJMAX, NCKAMAX

      INTEGER CCSDMAXLAST, CCSDSYMLAST !parallelization synchronization target variables

      INTEGER NT1AMX, NT2AMX, NT1AM, NT2AM, NT2AMA, NT2AMT, NNBST,      &
     &        NDISAO, NDSRHF, IDSAOG, ILMRHF, NDISAOSQ, NDSRHFSQ,       &
     &        ILMVIR, IT1AM, IT2AM, IAODIS, IDSAOGSQ,                   &
     &        NLAMDT, NLMRHF, IRHF, IVIR , NEMAT1,                      &
     &        IEMAT1, NGAMMA, NGAMSQ,NT1AO, NT2AO,                      &
     &        IT1AO, IT2AO, N2BST, NT2BCD,                              &
     &        IT2BCD, NT2BGD, IT2BGD, IT2SQ,                            &
     &        NMATIJ, IMATIJ, IGAMMA,IDSRHF,IDSRHFSQ,IGAMSQ,            &
     &        IT1AOT, IT1AMT, IT2BCT,                                   &
     &        IT2BGT, IFCVIR, IFCRHF,                                   &
     &        IMATAB, NMATAB, NT2AOS, IT2AOS,                           &
     &        NT2SQ, NMIJP, IMIJP, NT2ORT,                              &
     &        IT2ORT,NT2AOIJ,IT2AOIJ,NT2ORT3,IT2ORT3,                   &
     &        ICKID,NCKI,ICKI,ICKITR,                                   &
     &        NTOTOC,NTRAOC,ICKASR,NCKIJ,                               &
     &        ICKAD,NCKA,ICKA,ICKALP,                                   &
     &        ICKATR,NCKATR,ISAIK,ISAIKJ,                               &
     &        NMAJIK,ISJIK,ISJIKA,ISAIKL,                               &       
     &        IMATAV,NMATAV,ICKDAO,NT1AOX,                              &
     &        ICKBD,ISYMOP,IPRCC,NGLMDT,NGLMRH,                         &
     &        IGLMRH,IGLMVI,NT2MMO,NT2MAO,                              &
     &        IT2SP,ISEC,ICOR,                                          &
     &        ILRHSI,ILVISI,NLRHSI,NLRHFR,                              &
     &        I3ODEL,I3ORHF,I3OVIR,                                     &
     &        N3ODEL,N3ORHF,N3OVIR,IT2AIJ,                              &
     &        NT2AIJ,NMAIJK,IMAIJK,NCKASR,                              &
     &        IMAIJA,NMAIJA,ID2IJG,ID2AIG,                              &
     &        ID2ABG,ND2IJG,ND2AIG,ND2ABG,                              &
     &        IFCKDO,IFCKDV, NDSGRH, IT2VO,                             &
     &        I3VOOO,N3VOOO,                                            &
     &        NU2AMX, NU2AM, IU2AM,                                     &
     &        ICDKAO,ICDKVI,N3VDEL, N3VVIR,                             &
     &        I3VDEL, I3VVIR,IMAABC, NMAABC,                            &
     &        IMAABCI,NMAABCI,IMAABI,NMAABI,                            &
     &        IMAAB_CI,NMAAB_CI,IMAJBAI,                                &
     &        IMAAOBCI,NMAAOBCI,IMAJBAIT,                               &
     &        IMAIAB,NMAIAB,IMAIAJ,NMAIAJ,                              &
     &        IMAAOBI,NMAAOBI,IG1AM,NG1AM,                              &
     &        NH1AMX,NH1AM,IH1AM,NH2AMX,NH2AM,IH2AM,                    &
     &        NVAJKL,IVAJKL,NVABKL,IVABKL,NT1VM,NT1VMX,                 &
     &        I3AORHF,N3AORHF,I3AO,N3AO,IRHF3O,NRHF3O,                  &
     &        IRHFA,IRHFB,                                              &
     &        NMATKL,NMATKI,NTR12AM,NR12R12P,NTR12SQ,NR12R12SQ,         &
     &        IMATKL,IMATKI,ITR12AM,IR12R12P,ITR12SQ,ITR12SQT,          &
     &        IR12R12SQ,NT2R12,IT2R12,NTG2SQ,ITG2SQ 

      LOGICAL OMEGSQ,T2TCOR,OMEGOR,CC3LR,RSPIM,LSEC,LCOR,NEWGAM,INTTR,  &
     &        TRIPIM



      COMMON /CCSDMAX/ NCKIJMAX, NCKAMAX


      COMMON /CCSDMAX/ CCSDMAXLAST
      !  Very important !!!
      !  Always keep CCSDMAXLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.



      COMMON /CCSDSYM/ NT1AMX, NT2AMX, NT1AM(8), NT2AM(8), NT2AMA(8),   &
     &                 NT2AMT(8), NNBST(8),                             &
     &                 NDISAO(8),NDISAOSQ(8),IDSAOG(8,8),IDSAOGSQ(8,8), &
     &                 NDSRHF(8),NDSRHFSQ(8),IDSRHF(8,8),IDSRHFSQ(8,8), &
     &                 ILMRHF(8),ILMVIR(8), IT1AM(8,8), IT2AM(8,8),     &
     &                 NLAMDT, NLMRHF, IRHF(8), IVIR(8) , NEMAT1(8),    &
     &                 IEMAT1(8,8), NGAMMA(8), NT1AO(8), NT2AO(8),      &
     &                 NGAMSQ(8),                                       &
     &                 IT1AO(8,8), IT2AO(8,8), N2BST(8), NT2BCD(8),     &
     &                 IT2BCD(8,8), NT2BGD(8), IT2BGD(8,8), IT2SQ(8,8), &
     &                 NMATIJ(8), IMATIJ(8,8), IGAMMA(8,8),IGAMSQ(8,8), &
     &                 IT1AOT(8,8),IT1AMT(8,8),IT2BCT(8,8),IT2BGT(8,8), &
     &                 IFCVIR(8,8), IFCRHF(8,8),IAODIS(8,8),            &
     &                 IMATAB(8,8), NMATAB(8), NT2AOS(8), IT2AOS(8,8),  &
     &                 NT2SQ(8), NMIJP(8), IMIJP(8,8), NT2ORT(8),       &
     &                 IT2ORT(8,8),NT2AOIJ(8),IT2AOIJ(8,8),             &
     &                 NT2ORT3(8),IT2ORT3(8,8),                         &
     &                 ICKID(8,8),NCKI(8),ICKI(8,8),ICKITR(8,8),        &
     &                 NTOTOC(8),NTRAOC(8),ICKASR(8,8),NCKIJ(8),        &
     &                 ICKAD(8,8),NCKA(8),ICKA(8,8),ICKALP(8,8),        &
     &                 ICKATR(8,8),NCKATR(8),ISAIK(8,8),ISAIKJ(8,8),    &
     &                 NMAJIK(8),ISJIK(8,8),ISJIKA(8,8),ISAIKL(8,8),    &
     &                 IMATAV(8,8),NMATAV(8),ICKDAO(8,8),NT1AOX,        &
     &                 ICKBD(8,8),ISYMOP,IPRCC,NGLMDT(8),NGLMRH(8),     &
     &                 IGLMRH(8,8),IGLMVI(8,8),NT2MMO(8,8),NT2MAO(8,8), &
     &                 IT2SP(8,8),ISEC(8),ICOR(8),                      &
     &                 ILRHSI(8),ILVISI(8),NLRHSI,NLRHFR(8),            &
     &                 I3ODEL(8,8),I3ORHF(8,8),I3OVIR(8,8),             &
     &                 N3ODEL(8),N3ORHF(8),N3OVIR(8),IT2AIJ(8,8),       &
     &                 NT2AIJ(8),NMAIJK(8),IMAIJK(8,8),NCKASR(8),       &
     &                 IMAIJA(8,8),NMAIJA(8),ID2IJG(8,8),ID2AIG(8,8),   &
     &                 ID2ABG(8,8),ND2IJG(8),ND2AIG(8),ND2ABG(8),       &
     &                 ICDKAO(8,8),ICDKVI(8,8),                         &
     &                 IFCKDO(8),IFCKDV(8),NDSGRH(8,8), IT2VO(8,8),     &
     &                 IMAABC(8,8), NMAABC(8), I3VDEL(8,8),             &
     &                 I3VVIR(8,8), N3VDEL(8), N3VVIR(8),               &
     &                 OMEGSQ,T2TCOR,OMEGOR,CC3LR,RSPIM,LSEC,LCOR,      &
     &                 NEWGAM,INTTR,TRIPIM,                             &
     &                 I3VOOO(8,8),N3VOOO(8),                           &
     &                 NU2AMX,NU2AM(8),IU2AM(8,8),                      &
     &                 IMAABCI(8,8),NMAABCI(8),IMAABI(8,8),NMAABI(8),   &
     &                 IMAAB_CI(8,8),NMAAB_CI(8),IMAJBAI(8,8),          &
     &                 IMAAOBCI(8,8),NMAAOBCI(8),IMAJBAIT(8,8),         &
     &                 IMAIAB(8,8),NMAIAB(8),IMAIAJ(8,8),NMAIAJ(8),     &
     &                 IMAAOBI(8,8),NMAAOBI(8),NH1AMX,NH1AM(8),         &
     &                 IH1AM(8,8),NH2AMX,NH2AM(8),IH2AM(8,8),           &
     &                 IG1AM(8,8),NG1AM(8),                             &
     &                 NVAJKL(8),IVAJKL(8,8),NVABKL(8),IVABKL(8,8),     &
     &                 NT1VM(8),NT1VMX,I3AORHF(8,8),N3AORHF(8),         &
     &                 I3AO(8,8),N3AO(8),IRHF3O(8,8),NRHF3O(8),IRHFA(8),&
     &                 IRHFB(8),NMATKL(8),NMATKI(8),NTR12AM(8),         &
     &                 NR12R12P(8),NTR12SQ(8),NR12R12SQ(8),             &
     &                 IMATKL(8,8),IMATKI(8,8),ITR12AM(8,8),            &
     &                 IR12R12P(8,8),ITR12SQ(8,8),ITR12SQT(8,8),        &
     &                 IR12R12SQ(8,8),NT2R12(8),IT2R12(8,8),            &
     &                 NTG2SQ(8),ITG2SQ(8,8)                           
!

      COMMON /CCSDSYM/ CCSDSYMLAST
      !  Very important !!!
      !  Always keep CCSDSYMLAST as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.


      INTEGER A,B,C,D,E,F,G,P,Q,R,S,I,J,K,L,M,N
