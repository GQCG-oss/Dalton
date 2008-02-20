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
/* 
   This version of CC was generated from version used 1998/99 on
   jensen by a merge with the version used on Linux PCs
   Christof Haettig, April 1999

   fixed for 64 bit mode on IBM AIX with VAR_INT64 and SYS_AIX set
   Christof Haettig, Mar 19 2003
   
   General interfacing with fortran code compiled with non-default
   integer size done by Pawel Salek, Feb 2008.
*/

/* Simulate the cray 64 bit word addressable I/O routines and
   allow for large buffering in core.


   CALL WOPEN(UNIT, NAME, LENNAM, BLOCKS, STATS, IERR)
                    ----  ------
   
   CALL WCLOSE(UNIT, IERR)

   CALL GETWA(UNIT, RESULT, ADDR, COUNT, IERR)

   CALL PUTWA(UNIT, SOURCE, ADDR, COUNT, IERR) 

   Currently the I/O is syncronous and unbuffered */
#define _BSD_SOURCE 1
#define _LARGEFILE64_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>

#define max_file 99

/* define here the integer type, which has the same lengths
   as the integers used under fortran */
#if defined (VAR_INT64) || defined (SYS_AIX)
#include <stdint.h>
typedef int64_t INTEGER;
#else
typedef int INTEGER;
#endif

/* Mark Fortran-callable API with FSYM */
#if defined(NO_UNDERSCORE) || (defined(FUNDERSCORE) &&FUNDERSCORE == 0)
#define FSYM(a) a
#define FSYM2(a) a
#else
#define FSYM(a) a ## _
#if (defined(FUNDERSCORE) && FUNDERSCORE == 2)
#define FSYM2(a) a ## __
#else
#define FSYM2(a) a ## _
#endif
#endif

/* Disable old cruft. */
#if defined(OLD_CRUFT)
/* for machines, which don´t have lseek64, we overwrite here
   lseek64 by lseek and off64_t by off_t
*/
#if defined (SYS_DEC) || defined (SYS_HPUX)
#define lseek64 lseek
#define off64_t off_t
#endif

#if defined (SYS_FUJITSU)
#ifdef __uxp__
#define L_XTND SEEK_END
#define lseek64 lseek
#define _LLTYPES
#endif

#endif /* OLD_CRUFT */

#endif 

static int first_call = 1;   /* need to do stuff on the first call */

static struct w_file {
  int fds;                     /* file descriptor */
  off64_t length;              /* file length in bytes */
  off64_t position;            /* current file position in bytes a la lseek */
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

/** returns ai with the wall clock time in seconds as a double.  The
    actual accuracy is OS-dependent. */
void
FSYM(walltm)(double *ai)
{
  struct timeval tp;
  struct timezone tzp;

  (void) gettimeofday(&tp,&tzp);
  *ai = (double) tp.tv_sec + ((double) tp.tv_usec) * 1.0e-6;
}


static INTEGER
isUnitValidAndOpen(INTEGER unit)
{
  if ( (unit < 0) || (unit >= max_file) )
    return -1;
  
  if ( file_array[unit].fds == -1 )
    return -1;

  return 0;
}

static INTEGER
isAddressValid(INTEGER addr)
{
  return (addr <= 0) ?  -4 : 0;
}

static INTEGER
isCountValid(INTEGER count)
{
  return (count < 0) ? -4 : 0;
}

static void
InitFileStats(struct w_file* file)
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

static void
PrintFileStats(INTEGER unit, const struct w_file *file)
{
  double ave_read=0.0e0, ave_write=0.0e0;
  double rate_read=0.0e0, rate_write=0.0e0;

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

  fflush(stdout);
  fprintf(stderr,"CRAYIO: Statistics for unit %d, file '%s', length=%lu bytes.\n",
          unit, file->path, (unsigned long)file->length);
  fprintf(stderr,
          "CRAYIO: oper :  #req.  :  #seek  :   #words  :"
          " #w/#req : time(s) :  MW/s \n"
          "CRAYIO: read : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
          file->n_read, file->seek_read, (int) file->words_read, 
          (int) ave_read, file->time_read, rate_read);
  fprintf(stderr,"CRAYIO:write : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
          file->n_write, file->seek_write, (int) file->words_write, 
          (int) ave_write, file->time_write, rate_write);
}

static void
InitFileData(struct w_file *file)
{
  file->fds = -1;
  file->length = (off64_t) -1;
  file->path = NULL;
  file->position = (off64_t) -1;
}

static void
FirstCall()
     /* Initialization on first call to anything */
{
  INTEGER i;

  for (i=0; i<max_file; i++) {

    InitFileData(&file_array[i]);
    InitFileStats(&file_array[i]);
  }

  first_call = 0;
}

void
FSYM(wclose)(const INTEGER *unit, INTEGER *ierr)
{
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = isUnitValidAndOpen(*unit))
    return;

  file = file_array + *unit;

  *ierr = close(file->fds);

  if (file->stats)
    PrintFileStats(*unit, file);
  
  InitFileData(file);

  InitFileStats(file);
}

/* ARGSUSED */
void FSYM(wopen)(const INTEGER *unit, const char *name, const INTEGER *lennam,
		 const INTEGER* blocks, const INTEGER *stats, INTEGER *ierr)
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

void
FSYM(getwa)(const INTEGER *unit, double *result, const INTEGER *addr, 
	    const INTEGER *count, INTEGER *ierr)
{
  size_t nbytes, con2;
  off64_t where, con1;
  double start, end;
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = isUnitValidAndOpen(*unit))
    return;

  if (*ierr = isAddressValid(*addr))
    return;

  if (*ierr = isCountValid(*count))
    return;

  file = file_array + *unit;

  con1 = *addr;
  con2 = *count;

  nbytes = con2 * 8;
  where = (con1 - (off64_t) 1) * (off64_t) 8;

  if ( (where+nbytes) > file->length ) {
    *ierr = -5;
    fflush(stdout);
    fprintf(stderr,
            "GETWA: where %lu \n"
            "GETWA: nbytes %lu \n"
            "GETWA: file->length %lu \n",
            (unsigned long)where, (unsigned long)nbytes,
	    (unsigned long)file->length);
    PrintFileStats(*unit, file);
    return;
  }

  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    file->seek_read++;
    if ( (file->position = lseek64(file->fds, where, SEEK_SET)) == (off64_t) -1) {
      *ierr = -4;
      return;
    }
  }

  if ( read(file->fds, result, nbytes) != nbytes) {
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
  
void
FSYM(putwa)(const INTEGER *unit, const double *source, const INTEGER *addr,
	    const INTEGER *count, INTEGER *ierr)
{
  size_t nbytes,con2;
  off64_t where, con1;
  double start, end;
  struct w_file *file;

  if (first_call)
    FirstCall();

  if ( *ierr = isUnitValidAndOpen(*unit))
    return;

  if (*ierr = isAddressValid(*addr))
    return;

  if (*ierr = isCountValid(*count))
    return;

  file = file_array + *unit;

  con1 = *addr;
  con2 = *count;

  nbytes = con2 * 8;
  where = (con1 - (off64_t) 1) * (off64_t) 8;
  
  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    file->seek_write++;
    if ( (file->position = lseek64(file->fds, where, SEEK_SET)) == (off64_t) -1) {
      *ierr = -4;
      return;
    }
  }

  if ( (*ierr=write(file->fds, source, nbytes)) != nbytes) {
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
