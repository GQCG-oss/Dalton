/*

!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!

!
*/
/** scatter-io.c:

    efficient(?) implementation of a scatter-gather io to be used from
    Fortran routines.

    Pawel Salek, 2003.12.03

    General comments: 
    - the file position counters point to "words" (i.e. integers).
    Public Routines:
    
*/

/* XOPEN_SOURCE needed for fseeko */
#define _XOPEN_SOURCE 500
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "../dft/general.h"

#if defined(SYS_DEC)
/* DEC C V5.9-005 on Digital UNIX V4.0 (Rev. 1229) has no fseeko,
 * fall back to old interface. */
#define fseeko fseek
#endif

/* define max number of open files at the same time so that we can
 * simplify the code and allocate the table of open files statically
 * instead of using realloc(). */
#define DA_MAXOPEN 10


/* block size (in bytes): space on the disk is allocated in chunks of
   this size to be able to address large amounts of data with 32-bit
   addresses/integers.  Generally, the smaller block size the better
   on the other hand having too small block size may lead to lot of
   wasted space if blocks are not filled... Original code has blocks
   of 4096 bytes.  Excessive block size can increase the file size for
   eg HSOINT.  Too small block size on the other hand may wrap "record
   numbers" on 32-bit architectures - and some people learnt it the
   hard way... We can have no blocks (ie. block size=1) on 64-bit
   architectures with current disk sizes.
*/
#define DA_BLOCKSZ 1024

static struct file_data {
    FILE *file;
    char *name;
} open_files[DA_MAXOPEN] = { {NULL,NULL} };
static int da_get_first_free_entry();

#define CHECK_ARGS(expr) \
do {                                                               \
  if(!(expr))                                                      \
      dalton_quit("file %s: line %d: argument test '%s' failed\n", \
                  __FILE__, __LINE__, #expr);                      \
} while(0)

/** fastio_() is used just for setting various options, tracing etc.
 * and is not used any more. We could remove just calls to it. */
void
FSYM(fastio)(void)
{
    printf("fastio called.\n");
}

/** daopen(lu,name) opens given file for reading/writing. If the file
 * does not exist, it is created. */
void
FSYM(daopen)(integer *lu, const char *name, int len)
{
    int idx;
    *lu = idx = da_get_first_free_entry();
    (*lu)++;
    if(*lu<1) /* abort!? */
        dalton_quit("daopen: too many open files!\n");

    open_files[idx].name = dal_malloc(len+1);
    memcpy(open_files[idx].name, name, len);
    open_files[idx].name[len] = '\0';
    open_files[idx].file = fopen(open_files[idx].name, "r+");
    if(!open_files[idx].file)
        open_files[idx].file = fopen(open_files[idx].name, "w+");
    if(!open_files[idx].file) {
        *lu = -1;
        dalton_quit("daopen: could not open file '%s'!\n", 
                    open_files[idx].name);
    }
}

/** darmov_() closes file identified by lu and removes it. */
void
FSYM(darmov)(integer *lu)
{
    int idx = *lu-1;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    fclose(open_files[idx].file);
    unlink(open_files[idx].name);
    free  (open_files[idx].name);
    open_files[idx].file = NULL;
    open_files[idx].name = NULL;
    *lu = -1;
}

/** daclos_() closes file identified by lu */
void
FSYM(daclos)(integer *lu)
{
    int idx = *lu-1;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    fclose(open_files[idx].file);
    free  (open_files[idx].name); 
    open_files[idx].file = NULL;
    open_files[idx].name = NULL;
    *lu = -1;
}

/** dawrite() writes wrdlen words (integers) starting at buffer to
 * file lu at position pos (in integers). pos is updated on exit. */
void
FSYM(dawrite)(integer *lu, void *buffer, integer *wrdlen, integer *pos)
{
    int idx = *lu-1;
    off_t off;
    integer occupied_blocks = 1 + ((*wrdlen)*sizeof(integer)-1)/DA_BLOCKSZ;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    off = *pos; fseeko(open_files[idx].file, off*DA_BLOCKSZ, SEEK_SET);
    if(fwrite(buffer, sizeof(integer), 
              *wrdlen, open_files[idx].file) != *wrdlen)
        dalton_quit("dawrite: writing error encountered.\n");
    *pos += occupied_blocks;
}

/** dawrite() reads wrdlen words (integers) to the buffer from the
 * file lu at position pos (in integers). pos is updated on exit. */
void
FSYM(daread)(integer *lu, void *buffer, integer *wrdlen, integer *pos)
{
    int idx = *lu-1;
    off_t off;
    integer occupied_blocks = 1 + ((*wrdlen)*sizeof(integer)-1)/DA_BLOCKSZ;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    off = *pos; fseeko(open_files[idx].file, off*DA_BLOCKSZ, SEEK_SET);
    if(fread(buffer, sizeof(integer),
             *wrdlen, open_files[idx].file) != *wrdlen) {
        fort_print("problematic position: %ld offset %ld\n", (long)*pos, off);
        dalton_quit("daread: error at pos %ld words (%ld words to read) ferror=%d.",
                    (long)off, (long)*wrdlen, ferror(open_files[idx].file));
    }
    *pos += occupied_blocks;
}

/* darelist_() reads from a file lu, at position pos, a list
 * consisting of NLIST elements. It assumes that the list elements
 * follow in pairs (arr, len) and that array consists of WORDs defined
 * as INTEGERs.
 */

void
FSYM(darelist)(integer *lu, integer *pos, integer *nlist,...)
{
    va_list alist;
    int i;
    int idx = *lu-1;
    off_t off;
    integer occupied_blocks, wrdlen = 0;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    va_start(alist, nlist);
    off = *pos; fseeko(open_files[idx].file, off*DA_BLOCKSZ, SEEK_SET);
    for(i=0; i<*nlist; i++) {
        integer *arr = va_arg(alist, integer*);
        integer *len = va_arg(alist, integer*);
        if( fread(arr, sizeof(integer), *len, open_files[idx].file) != *len)
            dalton_quit("darelist: failed on read.\n");
        wrdlen += *len;
    }
    va_end(alist);
    occupied_blocks = 1 + (wrdlen*sizeof(integer)-1)/DA_BLOCKSZ;
    *pos += occupied_blocks;
}

/* dawrlist_() writes to a file lu, at position pos, a list
 * consisting of NLIST elements. It assumes that the list elements
 * follow in pairs (arr, len) and that array consists of WORDs defined
 * as INTEGERs.
 */
void
FSYM(dawrlist)(integer *lu, integer *pos, integer *nlist,...)
{
    va_list alist;
    int i;
    int idx = *lu-1;
    off_t off;
    integer occupied_blocks, wrdlen = 0;
    CHECK_ARGS(idx>=0 && idx <DA_MAXOPEN);
    CHECK_ARGS(open_files[idx].file != NULL);

    va_start(alist, nlist);
    off = *pos; fseeko(open_files[idx].file, off*DA_BLOCKSZ, SEEK_SET);
    for(i=0; i<*nlist; i++) {
        integer *arr = va_arg(alist, integer*);
        integer *len = va_arg(alist, integer*);
        if( fwrite(arr, sizeof(integer), *len, open_files[idx].file) != *len)
            dalton_quit("darelist: failed on read.\n");
        wrdlen += *len;
    }
    va_end(alist);
    occupied_blocks = 1 + (wrdlen*sizeof(integer)-1)/DA_BLOCKSZ;
    *pos += occupied_blocks;
}

/* daskip updates pos taking into account block size */
void
FSYM(daskip)(integer *lu, integer *wrdlen, integer *pos)
{
    integer occupied_blocks = 1 + ((*wrdlen)*sizeof(integer)-1)/DA_BLOCKSZ;
    *pos += occupied_blocks;
}

/* =================================================================== *
 *                  LOCAL AUXILLIARY ROUTINES                          *
 * =================================================================== */
static int
da_get_first_free_entry()
{
    int i;
    for(i=0; i<DA_MAXOPEN && open_files[i].file; i++)
        ;
    return i<DA_MAXOPEN ? i : -1;
}

