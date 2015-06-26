/*!
\file
\brief Various functions for transfering data to/from slave/master memory
\date Started 6/3/2013
\author George
*/


#include "common.h"


/* the working directory; initialized during init */
static char xfer_wdir[BDMPI_WDIR_LEN];

/* disk I/O will be done in chunks of this number of bytes */
#define BDMPI_DISK_CHUNK      131072


/*************************************************************************/
/*! Sets the static working directory. */
/*************************************************************************/
void xfer_setwdir(char *wdir)
{
  strcpy(xfer_wdir, wdir);
}


/*************************************************************************/
/*! Returns a globally unique file number. */
/*************************************************************************/
ssize_t xfer_getfnum()
{
  static ssize_t nextfnum=1;

  nextfnum++;

  return ((ssize_t)getpid())*10000000 + nextfnum;
}


/*************************************************************************/
/*! Unlinks the file associated with fnum. */
/*************************************************************************/
void xfer_unlink(ssize_t fnum)
{
  char *fname;

  if (asprintf(&fname, "%s/%zd", xfer_wdir, fnum) == -1)
    errexit("xfer_unlink: Failed to create filename.\n");

  unlink(fname);

  free(fname);

  return;
}


/*************************************************************************/
/*! Copies data into remote buffer space via scb. */
/*************************************************************************/
void xfer_out_scb(bdscb_t *scb, void *vbuf, size_t count, BDMPI_Datatype datatype)
{
  size_t i, dtsize, chunk, len;
  char *buf = (char *)vbuf;

  dtsize = bdmp_sizeof(datatype);
  chunk  = scb->size/dtsize;

  /* get into a loop copying and storing the data */
  for (i=0; i<count; i+=chunk) {
    len = (i+chunk < count ? chunk : count - i);
    bdscb_wait_empty(scb);
    //bdprintf("out_scb: %zd bytes from %p [%zu %zu]\n", len*dtsize, buf+i*dtsize, count, i);
    memcpy(scb->buf, buf+i*dtsize, len*dtsize);
    bdscb_post_full(scb);
  }

  return;
}


/*************************************************************************/
/*! Copies data into local buffer space via scb. */
/*************************************************************************/
void xfer_in_scb(bdscb_t *scb, void *vbuf, size_t count, BDMPI_Datatype datatype)
{
  size_t i, dtsize, chunk, len;
  char *buf = (char *)vbuf;

  dtsize = bdmp_sizeof(datatype);
  chunk  = scb->size/dtsize;

  /* get into a loop copying and storing the data */
  for (i=0; i<count; i+=chunk) {
    len = (i+chunk < count ? chunk : count - i);
    bdscb_wait_full(scb);
    //bdprintf("in_scb: %zd bytes from %p [%zu %zu]\n", len*dtsize, buf+i*dtsize, count, i);
    memcpy(buf+i*dtsize, scb->buf, len*dtsize);
    bdscb_post_empty(scb);
  }

  return;
}


/*************************************************************************/
/*! Copies data into the specified file. */
/*************************************************************************/
void xfer_out_disk(ssize_t fnum, char *buf, size_t count, BDMPI_Datatype datatype)
{
  char *fname;
  int fd;
  size_t len, size = bdmp_msize(count, datatype);
  ssize_t ret;

  if (asprintf(&fname, "%s/%zd", xfer_wdir, fnum) == -1)
    errexit("xfer_out_disk: Failed to create filename.\n");

  if ((fd = open(fname, O_CREAT|O_WRONLY|O_TRUNC, S_IRUSR|S_IWUSR)) == -1)
    errexit("xfer_out_disk: Failed to open file %s: %s\n", fname, strerror(errno));

  //bdprintf("out_disk: %zu bytes from %p [%zu %zu] %s\n", size, buf, size,
  //  fnum, fname);
  do {
    len = (size > BDMPI_DISK_CHUNK ? BDMPI_DISK_CHUNK : size);
    if ((ret=gk_write(fd, buf, len)) != len)
      errexit("[%5d] xfer_out_disk: Write size does not match: %s %zd %zu\n",
        (int)getpid(), strerror(errno), ret, len);
    buf += len;
    size -= len;
  } while (size > 0);

  close(fd);

  free(fname);

  return;
}


/*************************************************************************/
/*! Copies data from the specified file. */
/*************************************************************************/
void xfer_in_disk(ssize_t fnum, char *buf, size_t count, BDMPI_Datatype datatype,
         int rmfile)
{
  char *fname;
  int fd;
  size_t len, rsize, size = bdmp_msize(count, datatype);

  if (asprintf(&fname, "%s/%zd", xfer_wdir, fnum) == -1)
    errexit("xfer_in_disk: Failed to create filename.\n");

  if ((rsize = gk_getfsize(fname)) == -1)
    errexit("xfer_in_disk: Error with file: %s [%s]\n", fname, strerror(errno));

  if (rsize > size)
    errexit("xfer_in_disk: Size of file %s is larger than supplied buffer space: [%zu %zu]\n",
        fname, rsize, size);

  if ((fd = open(fname, O_RDONLY)) == -1)
    errexit("xfer_in_disk: Failed to open file %s: %s\n", fname, strerror(errno));

  //bdprintf("in_disk: %zu bytes from %p [%zu %zu]\n", rsize, buf, size, fnum);

  size = rsize;
  do {
    ssize_t retval;
    len = (size > BDMPI_DISK_CHUNK ? BDMPI_DISK_CHUNK : size);
    if ((retval=gk_read(fd, buf, len)) != len)
      errexit("[%5d] xfer_in_disk: Read size does not match: %s (%zd, %zu) " \
        "(%d) %p\n", (int)getpid(), strerror(errno), retval, len, errno, buf);
    buf  += len;
    size -= len;
  } while (size > 0);

  close(fd);

  if (rmfile)
    unlink(fname);

  free(fname);

  return;
}
