************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     This code or parts thereof may not be copied or redistributed,   *
*     without the written permission of the author. The use is re -    *
*     stricted to research purposes only, and the matrial may not be   *
*     included in any commercial product.                              *
*                                                                      *
************************************************************************
*                                                                      *
*     This file defines the common blocks /FIO1/ - /FIO4/              *
*     including all entries required for fast I/O system.              *
*                                                                      *
*     Following the list of entries and their usage:                   *
*     Name    : Name of the file                                       *
*     Open    : open/close flag                                        *
*     FSCB    : file descriptors (C language)                          *
*     Addr    : pointer to the current position                        *
*     MxAddr  : largest address written                                *
*     Count   : counters for the I/O statistics                        *
*               (# of I/O startups, amount of data transferred)        *
*     Trace   : enable/disable debugging output                        *
*     Query   : enable/disable the traceback facilities                *
*     FirstCall : enable/disable multifile partitiong                  *
*     RefBuf  : anchor of the gather/scatter process                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L Serrano-Andres,University of Lund, Sweden, 1996                *
*                                                                      *
************************************************************************

      Integer Open,FSCB,Addr,MxAddr,RefBuf,MaxFileSize
      Common /FIO1/ Open(99),
     &              FSCB(99),
     &              Addr(99),
     &              MxAddr(99),
     &              RefBuf(8),
     &              MaxFileSize

      Logical     Trace,Query,FirstCall
      Common /FIO2/ Trace,
     &              Query,
     &              FirstCall

      Real*8      Count
      Common /FIO3/ Count(4,99)

      Character*8 Name
      Common /FIO4/ Name(99)


