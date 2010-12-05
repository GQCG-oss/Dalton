/*
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras, T. Saue, 
C...   S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren. 
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
*/
/* extra routines needed by intel's ifc compiler which does not
 * provide them itself.
 * Pawel Salek, pawsa@theochem.kth.se, 2002.01.31
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

void getenv_(const char* str, char* val, int strl, int vall)
{
  char name[32];
  char *envv;
  int i;

  if(strl>sizeof(name)-1) {
    printf("oops\n");
    return;
  }
  strncpy(name, str, strl);
  name[strl] = '\0';

  /*printf("getenv(%s) called with strl=%d, vall=%d '%s'\n",
    name, strl, vall, getenv(name));*/
  if( (envv=getenv(name)) != NULL)
    strncpy(val, envv, vall);
  else val[0] = '\0';
  /*printf("getenv(%s) returned '%s' (%d, %d)\n", name, val, strl, vall);*/
  for(i=strlen(val); i< vall; i++)
    val[i] = ' ';
}

void system_(char* str, int len)
{
    char buff[256];
    if(len>sizeof(buff)-1) 
        printf("Bloody crayio. I ignore your commands.\n");
    else {
        int sz = sizeof(buff)-1>len ? len : sizeof(buff)-1;
        strncpy(buff, str, sizeof(buff)-1);        
        buff[sz] = '\0';
        system(buff);
    }
}

void flush_(int* lu)
{
  /* no flushing */
}

float etime_(float* et)
{
  static int ticksclk = -1;
  struct tms tm;

  times(&tm);
  if(ticksclk<0) ticksclk = sysconf(_SC_CLK_TCK);  
  et[0] = tm.tms_utime/(double)ticksclk;
  et[1] = tm.tms_stime/(double)ticksclk;
  return et[0] + et[1];
}

#ifdef VAR_G77
char* fdate_(void)
{
  static char buf[24];
  time_t curr_time;
  time(&curr_time);
  strncpy(buf, ctime(&curr_time), 24);
  return buf;
}
#else
void fdate_(char* dt, int dtlen)
{
  time_t curr_time;
  time(&curr_time);
  strncpy(dt, ctime(&curr_time), 24);
}
#endif

int time_(void)
{ 
  return (int)time(NULL);
}

void hostnm_(char* buf, int buflen)
{
  int i;
  gethostname(buf, buflen);
  /* printf("Hostname: '%s' (%d)\n", buf, buflen); */
  for(i=strlen(buf); i<buflen; i++)
    buf[i] = ' ';
}

double* malloc_(int cnt)
{
  return malloc(cnt);
}
