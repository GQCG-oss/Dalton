/*
       Copyright (c) 2001 by the authors of Dalton (see below).
       All Rights Reserved.
    
       The source code in this file is part of
       "Dalton, a molecular electronic structure program, Release 1.2
       (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
       J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
       O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
       T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
       H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
       H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
       T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
       B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
       and O. Vahtras"
    
       This source code is provided under a written licence and may be
       used, copied, transmitted, or stored only in accord with that
       written licence.
    
       In particular, no part of the source code or compiled modules may
       be distributed outside the research group of the licence holder.
       This means also that persons (e.g. post-docs) leaving the research
       group of the licence holder may not take any part of Dalton,
       including modified files, with him/her, unless that person has
       obtained his/her own licence.
    
       For questions concerning this copyright write to:
          dalton-admin@kjemi.uio.no
    
       For information on how to get a licence see:
          http://www.kjemi.uio.no/software/dalton/dalton.html
*/

/* 
   this version of CC was generated from version used 1998/99 on
   jensen by a merge with the version used on Linux PCs
   Christof Haettig, April 1999

   fixed for 64 bit mode on IBM AIX with VAR_INT64 and SYS_AIX set
   Christof Haettig, Mar 19 2003
*/

/* Simulate the cray 64 bit word addressable I/O routines and
   allow for large buffering in core.


   CALL WOPEN(UNIT, NAME, LENNAM, BLOCKS, STATS, IERR)
                    ----  ------
   
   CALL WCLOSE(UNIT, IERR)

   CALL GETWA(UNIT, RESULT, ADDR, COUNT, IERR)

   CALL PUTWA(UNIT, SOURCE, ADDR, COUNT, IERR) 

   Currently the I/O is syncronous and unbuffered */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>

#define max_file 99

/* for machines, which don´t have lseek64, we overwrite here
   lseek64 by lseek and off64_t by off_t
*/
#if defined (SYS_LINUX) || defined (SYS_DEC) || defined (SYS_HPUX)
#define lseek64 lseek
#define off64_t off_t
#endif

/* define here the integer type, which has the same lengths
   as the integers used under fortran */
#if defined (VAR_INT64) && defined (SYS_AIX)
#define INTEGER long
#else
#define INTEGER int
#endif

#if defined (SYS_FUJITSU)
#ifdef __uxp__
#define L_XTND SEEK_END
#define lseek64 lseek
#define _LLTYPES
#endif
#endif

/* for some machines we have to add an underscore to the routine names:  */
#if defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_AIX) || defined (SYS_SUN)
#define WOPEN wopen_
#define WCLOSE wclose_
#define GETWA getwa_
#define PUTWA putwa_
#else
#define WOPEN wopen
#define WCLOSE wclose
#define GETWA getwa
#define PUTWA putwa
#endif

static INTEGER first_call = 1;   /* need to do stuff on the first call */

static struct w_file {
  int fds;                     /* file descriptor */
  long long length;            /* file length in bytes */
  long long position;          /* current file position in bytes a la lseek */
  char *path;                  /* file name */
  INTEGER stats;                   /* boolean flag to collect statistics */
  double words_write;          /* total no. of words written */
  double words_read;           /* total no. of words read */
  double time_write;           /* total wall time writing */
  double time_read;            /* total wall time reading */
  INTEGER n_read;                  /* no. of read requests */
  INTEGER n_write;                 /* no. of write reqeusts */
  INTEGER seek_read;               /* no. of seeks on read */
  INTEGER seek_write;              /* no. of seeks on write */
} file_array[max_file];

void walltm_(ai)
    double *ai;
/* return ai with the wall clock time in seconds as a double.
   it might be accurate to about 0.01s at best */
{
  struct timeval tp;
  struct timezone tzp;

  (void) gettimeofday(&tp,&tzp);
  *ai = (double) tp.tv_sec + ((double) tp.tv_usec) * 1.0e-6;
}

static INTEGER CheckUnit(unit)
     INTEGER unit;
{
  if ( (unit < 0) || (unit >= max_file) )
    return -1;
  
  if ( file_array[unit].fds == -1 )
    return -1;

  return 0;
}

static INTEGER CheckAddr(addr)
     INTEGER addr;
{
  if (addr <= 0)
    return -4;
  else
    return 0;
}

static INTEGER CheckCount(count)
     INTEGER count;
{
  if (count < 0)
    return -4;
  else
    return 0;
}

void InitFileStats(file)
     struct w_file *file;
{
  file->stats = 1;
  file->words_write = 0.0e0;
  file->words_read = 0.0e0;
  file->time_read = 0.0e0;
  file->time_write = 0.0e0;
  file->n_read = 0;
  file->n_write = 0;
  file->seek_read = 0;
  file->seek_write = 0;
}

void PrintFileStats(unit, file)
     struct w_file *file;
     INTEGER unit;
{
  double ave_read=0.0e0, ave_write=0.0e0;
  double rate_read=0.0e0, rate_write=0.0e0;

  (void) fflush(stdout);

  if (file->n_read) {
    ave_read = file->words_read / (double) file->n_read;
    if (file->time_read > 0.0e0)
      rate_read = file->words_read / (1000000.0e0 * file->time_read);
  }

  if (file->n_write) {
    ave_write = file->words_write / (double) file->n_write;
    if (file->time_write > 0.0e0)
      rate_write = file->words_write / (1000000.0e0 * file->time_write);
  }

  (void) fflush(stdout);
  (void) fprintf(stderr,"CRAYIO: Statistics for unit %d, file '%s', length=%ld bytes.\n",
		 unit, file->path, file->length);
  (void) fprintf(stderr,"CRAYIO: oper :  #req.  :  #seek  :   #words  :");
  (void) fprintf(stderr," #w/#req : time(s) :  MW/s \n");
  (void) fprintf(stderr,"CRAYIO: read : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
		 file->n_read, file->seek_read, (INTEGER) file->words_read, 
		 (INTEGER) ave_read, file->time_read, rate_read);
  (void) fprintf(stderr,"CRAYIO:write : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
		 file->n_write, file->seek_write, (INTEGER) file->words_write, 
		 (INTEGER) ave_write, file->time_write, rate_write);
}

void InitFileData(file)
     struct w_file *file;
{
  file->fds = -1;
  file->length = (long long) -1;
  file->path = (char *) NULL;
  file->position = (long long) -1;
}

void FirstCall()
     /* Initialization on first call to anything */
{
  INTEGER i;

  for (i=0; i<max_file; i++) {

    InitFileData(&file_array[i]);

    InitFileStats(&file_array[i]);
  }

  first_call = 0;
}

void WCLOSE(unit, ierr)
     INTEGER *unit, *ierr;
{
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = CheckUnit(*unit))
    return;

  file = file_array + *unit;

  *ierr = close(file->fds);

  if (file->stats)
    PrintFileStats(*unit, file);
  
  InitFileData(file);

  InitFileStats(file);
}

/* ARGSUSED */
void WOPEN(unit, name, lennam, blocks, stats, ierr)
     INTEGER *unit, *lennam, *blocks, *stats, *ierr;
     char *name;
{
  struct w_file *file;

  *ierr = (INTEGER) 0;

  if (first_call)
    FirstCall();

  if ( (*unit < 0) || (*unit >= max_file) ) {
    *ierr = -1;
    return;
  }
    
  file = file_array + *unit;

  file->stats = *stats;

  if (*lennam > 0) {
    file->path = malloc((size_t) (*lennam + 1));
    (void) strncpy(file->path,name,(INTEGER) *lennam);
    /* file->path[*lennam] = NULL; */
    file->path[*lennam] = 0;
    }
  else {
    file->path = malloc((size_t) 8);
    (void) sprintf(file->path,"fort.%.2d",*unit);
  }
  if (( file->fds = open(file->path, (O_RDWR|O_CREAT), 0660))
      == -1) {
    *ierr = -6;
    return;
  }

  file->length = lseek64(file->fds, (off64_t) 0, SEEK_END);
  file->position = lseek64(file->fds, (off64_t) 0, SEEK_SET);

}

void GETWA(unit, result, addr, count, ierr)
     INTEGER *unit, *addr, *count, *ierr;
     double *result;
{
  long nbytes, con2;
  long long where, con1;
  double start, end;
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = CheckUnit(*unit))
    return;

  if (*ierr = CheckAddr(*addr))
    return;

  if (*ierr = CheckCount(*count))
    return;

  file = file_array + *unit;

  con1 = *addr;
  con2 = *count;

  nbytes = con2 * (long) 8;
  where = (con1 - (long long) 1) * (long long) 8;

  if ( (where+nbytes) > file->length ) {
    *ierr = -5;
    (void) fflush(stdout);
    (void) fprintf(stderr,"GETWA: where %ld \n",where);
    (void) fprintf(stderr,"GETWA: nbytes %ld \n",nbytes);
    (void) fprintf(stderr,"GETWA: file->length %ld \n",file->length);
    (void) PrintFileStats(*unit, file);
    (void) fflush(stdout);
    return;
  }

  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    file->seek_read++;
    if ( (file->position = lseek64(file->fds, (off64_t) where, SEEK_SET)) == (long long) -1) {
      *ierr = -4;
      return;
    }
  }

  if ((long) read(file->fds, (char *) result, (long) nbytes) != nbytes) {
    *ierr = -6;
    return;
  }
  
  file->position += nbytes;

  if (file->stats) {
    walltm_(&end);
    file->n_read++;
    file->words_read += (double) *count;
    file->time_read +=  end - start;
  }

  *ierr = 0;
}
  
void PUTWA(unit, source, addr, count, ierr)
     INTEGER *unit, *addr, *count, *ierr;
     double *source;
{
  size_t nbytes,con2;
  long long where, con1;
  double start, end;
/*  long long jerr; */
  struct w_file *file;

  if (first_call)
    FirstCall();

  if ( *ierr = CheckUnit(*unit))
    return;

  if (*ierr = CheckAddr(*addr))
    return;

  if (*ierr = CheckCount(*count))
    return;

  file = file_array + *unit;

  con1 = *addr;
  con2 = *count;

  nbytes = con2 * (long) 8;
  where = (con1 - (long long) 1) * (long long) 8;
  
  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    file->seek_write++;
    if ( (file->position = lseek64(file->fds, (off64_t) where, SEEK_SET)) == (long long) -1) {
      *ierr = -4;
      return;
    }
  }

  if ( (*ierr=write(file->fds, (char *) source, nbytes)) != nbytes) {
    printf("\n write returned %d \n",*ierr);
    *ierr = -6;
    return;
  }

  where += nbytes;
  file->position += nbytes;
  if (file->length < where)
    file->length = where;
  

  if (file->stats) {
    walltm_(&end);
    file->n_write++;
    file->words_write += (double) *count;
    file->time_write +=  end - start;
  }

  *ierr = 0;
}

