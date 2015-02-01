/*
 * bdmq.c
 *
 * The BDMP's functions for providing message queues.
 *
 * Started 4/5/2013
 * George
 *
 */


#include "common.h"


/*************************************************************************/
/*! Creates the neccessary components of a message queue.

    \param tag is a string to be used as a prefix of the MQ
    \param num is used to uniquely identify the MQ been created
           (eg., based on master or slave pid)
    \return the allocated MQ object.
*/
/*************************************************************************/
bdmq_t *bdmq_create(char *tag, int num)
{
  char name[256];
  bdmq_t *mq;
  struct mq_attr attr;

  mq = (bdmq_t *)gk_malloc(sizeof(bdmq_t), "bdmq_create: mq");

  /* set some deafult attributes */
  attr.mq_flags   = 0;
  attr.mq_maxmsg  = 512;
  attr.mq_msgsize = sizeof(bdmsg_t);
  attr.mq_curmsgs = 0;

  /* create the message queue */
  sprintf(name, "/%s%d", tag, num);
  mq->name = gk_strdup(name);
  mq_unlink(name);

  mq->mqdes = mq_open(name, O_CREAT|O_RDWR, S_IRUSR|S_IWUSR, &attr);
  //mq->mqdes = mq_open(name, O_CREAT|O_RDWR, S_IRUSR|S_IWUSR, NULL);
  if (mq->mqdes == -1)
    errexit("!! Failed on mq_open(mq->mqdes) for %s: %s\n", name, strerror(errno));

  mq_getattr(mq->mqdes, &attr);
  mq->msgsize = attr.mq_msgsize;
  mq->buf = gk_cmalloc(mq->msgsize, "mq->buf");

  //printf("mq_flags: %ld, mq_maxmsg: %ld, mq_msgsize: %ld, mq_curmsgs: %ld\n",
  //    attr.mq_flags, attr.mq_maxmsg, attr.mq_msgsize, attr.mq_curmsgs);

  return mq;
}


/*************************************************************************/
/*! Opens a previously created message queue.

    \param tag is a string to be used as a prefix of the MQ
    \param num is used to uniquely identify the MQ been created
           (eg., based on master or slave pid)
    \return the allocated MQ object.
*/
/*************************************************************************/
bdmq_t *bdmq_open(char *tag, int num)
{
  char name[256];
  bdmq_t *mq;
  struct mq_attr attr;

  mq = (bdmq_t *)gk_malloc(sizeof(bdmq_t), "bdmq_create: mq");

  /* open the message queue */
  sprintf(name, "/%s%d", tag, num);
  mq->name = gk_strdup(name);

  mq->mqdes = mq_open(name, O_RDWR);
  if (mq->mqdes == -1)
    errexit("Failed on mq_open(mq->mqdes): %s %s\n", name, strerror(errno));

  mq_getattr(mq->mqdes, &attr);
  mq->msgsize = attr.mq_msgsize;
  mq->buf = gk_cmalloc(mq->msgsize, "mq->buf");

  return mq;
}


/*************************************************************************/
/*! Closes a message queue. This should be called from the slaves. */
/*************************************************************************/
void bdmq_close(bdmq_t *mq)
{
  if (bdmq_length(mq) > 0)
    printf("Closing a non-empty message queue %s.\n", mq->name);

  if (mq_close(mq->mqdes) == -1)
    errexit("Failed on mq_close(mq->mqdes): %s\n", strerror(errno));

  gk_free((void **)&mq->name, &mq->buf, &mq, LTERM);
}


/*************************************************************************/
/*! Closes and destroys a message queue. This should be called from the
    master. */
/*************************************************************************/
void bdmq_destroy(bdmq_t *mq)
{
  if (bdmq_length(mq) > 0)
    printf("Destroying a non-empty message queue.\n");

  if (mq_close(mq->mqdes) == -1)
    errexit("Failed on mq_close(mq->mqdes): %s\n", strerror(errno));
  if (mq_unlink(mq->name) == -1)
    errexit("Failed on mq_unlink(mq->name): %s\n", strerror(errno));

  gk_free((void **)&mq->name, &mq->buf, &mq, LTERM);
}


/*************************************************************************/
/*! Sends a message to a message queue. */
/*************************************************************************/
int bdmq_send(bdmq_t *mq, void *buf, size_t size)
{
  if (size > mq->msgsize)
    errexit("bdmq_send: Message size %zu exceeds max limit of %zu.\n", size, mq->msgsize);

  return mq_send(mq->mqdes, buf, size, 0);
}


/*************************************************************************/
/*! Receives a message from a message queue.
    \returns the number of bytes received or -1. */
/*************************************************************************/
ssize_t bdmq_recv(bdmq_t *mq, void *buf, size_t size)
{
  unsigned int priority;
  ssize_t rsize;

  if ((rsize = mq_receive(mq->mqdes, mq->buf, mq->msgsize, &priority)) != -1) {
    if (size != rsize)
      errexit("bdmq_recv: Received message size %zu is not the same as requested size %zu.\n", rsize, size);

    if (size < rsize)
      return -1;
    memcpy(buf, mq->buf, rsize);
    return rsize;
  }
  return -1;
}


/*************************************************************************/
/*! Receives a message from a message queue with a timeout.
    \param dt is the number of nanoseconds from now!
    \returns the number of bytes received or -1. */
/*************************************************************************/
ssize_t bdmq_timedrecv(bdmq_t *mq, void *buf, size_t size, long dt)
{
  unsigned int priority;
  struct timespec abs_timeout;
  ssize_t rsize;

  //printf("[%d] timedrecv-in\n", (int)time(0));
  if (clock_gettime(CLOCK_REALTIME, &abs_timeout) == -1)
    return -1;

  if (abs_timeout.tv_nsec + dt < 1000000000)
    abs_timeout.tv_nsec += dt;
  else {
    abs_timeout.tv_sec++;
    abs_timeout.tv_nsec = (abs_timeout.tv_nsec + dt)%1000000000;
  }

  //abs_timeout.tv_nsec = gk_min(1000000000-1, abs_timeout.tv_nsec+dt);

  if ((rsize = mq_timedreceive(mq->mqdes, mq->buf, mq->msgsize, &priority, &abs_timeout)) != -1) {
    if (size < rsize)
      rsize = -1;
    else
      memcpy(buf, mq->buf, rsize);
  }
  //printf("[%d] timedrecv-out\n", (int)time(0));

  return rsize;
}


/*************************************************************************/
/*! Registers a function to be executed when a message is pending in the
    queue. */
/*************************************************************************/
void bdmq_register_function(bdmq_t *mq, void (*function)(union sigval))
{
  struct sigevent action;

  memset(&action, 0, sizeof(struct sigevent));
  action.sigev_notify            = SIGEV_THREAD;
  action.sigev_notify_function   = function;
  action.sigev_notify_attributes = NULL;
  action.sigev_value.sival_ptr   = mq;

  if (mq_notify(mq->mqdes, &action) == -1)
    errexit("Failed on mq_notify(mq->mqdes, &action): %s\n", strerror(errno));

}


/*************************************************************************/
/*! Returns the number of current messages in a queue.
*/
/*************************************************************************/
int bdmq_length(bdmq_t *mq)
{
  struct mq_attr attr;

  mq_getattr(mq->mqdes, &attr);
  return attr.mq_curmsgs;

}
