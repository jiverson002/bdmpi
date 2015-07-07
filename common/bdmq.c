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

  mq = (bdmq_t *)bd_malloc(sizeof(bdmq_t), "bdmq_create: mq");

  /* set some deafult attributes */
  attr.mq_flags   = 0;
  attr.mq_maxmsg  = 512;
  attr.mq_msgsize = sizeof(bdmsg_t);
  attr.mq_curmsgs = 0;

  /* create the message queue */
  sprintf(name, "/%s%d", tag, num);
  mq->name = bd_strdup(name);
  mq_unlink(name);

  mq->mqdes = mq_open(name, O_CREAT|O_RDWR, S_IRUSR|S_IWUSR, &attr);
  //mq->mqdes = mq_open(name, O_CREAT|O_RDWR, S_IRUSR|S_IWUSR, NULL);
  if (mq->mqdes == -1)
    errexit("!! Failed on mq_open(mq->mqdes) for %s: %s\n", name, strerror(errno));

  mq_getattr(mq->mqdes, &attr);
  mq->msgsize = attr.mq_msgsize;
  mq->buf = bd_malloc(mq->msgsize, "mq->buf");

  //bdprintf("create:%s: %zu\n", name, attr.mq_msgsize);

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

  mq = (bdmq_t *)bd_malloc(sizeof(bdmq_t), "bdmq_create: mq");

  /* open the message queue */
  sprintf(name, "/%s%d", tag, num);
  mq->name = bd_strdup(name);

  mq->mqdes = mq_open(name, O_RDWR);
  if (mq->mqdes == -1)
    errexit("Failed on mq_open(mq->mqdes): %s %s\n", name, strerror(errno));

  mq_getattr(mq->mqdes, &attr);
  mq->msgsize = attr.mq_msgsize;
  mq->buf = bd_malloc(mq->msgsize, "mq->buf");

  //bdprintf("open:%s: %zu\n", name, attr.mq_msgsize);

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

  bd_free((void **)&mq->name, &mq->buf, &mq, LTERM);
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

  bd_free((void **)&mq->name, &mq->buf, &mq, LTERM);
}


/*************************************************************************/
/*! Sends a message to a message queue. */
/*************************************************************************/
ssize_t bdmq_send(bdmq_t *mq, void *buf, size_t size)
{
  ssize_t ret;
  struct mq_attr attr;
  mq_getattr(mq->mqdes, &attr);

  if (size > mq->msgsize) {
    errexit("[%5d] bdmq_send: Message size %zu exceeds max limit of %zu (%s).\n",
      (int)getpid(), size, mq->msgsize, mq->name);
  }
  for (;;) {
    ret = mq_send(mq->mqdes, buf, size, 0);
    if (-1 == ret && EINTR != errno)
      break;
    else if (-1 != ret)
      break;
  }
  return ret;
}


/*************************************************************************/
/*! Receives a message from a message queue.
    \returns the number of bytes received or -1. */
/*************************************************************************/
ssize_t bdmq_recv(bdmq_t *mq, void *buf, size_t size)
{
  unsigned int priority;
  ssize_t ret;

  for (;;) {
    ret = mq_receive(mq->mqdes, mq->buf, mq->msgsize, &priority);
    if (-1 == ret && EINTR != errno)
      break;
    else if (-1 != ret)
      break;
  }

  if (-1 != ret) {
    if (size != ret)
      errexit("bdmq_recv: Received message size %zu is not the same as "
        "requested size %zu.\n", ret, size);

    if (size < ret)
      ret = -1;
    else
      memcpy(buf, mq->buf, ret);
  }

  return ret;
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
  ssize_t ret;

  /*printf("[%d] timedrecv-in\n", (int)time(0));*/
  if (clock_gettime(CLOCK_REALTIME, &abs_timeout) == -1)
    return -1;

  if (abs_timeout.tv_nsec + dt < 1000000000)
    abs_timeout.tv_nsec += dt;
  else {
    abs_timeout.tv_sec++;
    abs_timeout.tv_nsec = (abs_timeout.tv_nsec + dt)%1000000000;
  }
  //abs_timeout.tv_nsec = gk_min(1000000000-1, abs_timeout.tv_nsec+dt);

  for (;;) {
    ret = mq_timedreceive(mq->mqdes, mq->buf, mq->msgsize, &priority,
      &abs_timeout);
    if (-1 == ret && EINTR != errno)
      break;
    else if (-1 != ret)
      break;
  }

  if (-1 != ret) {
    if (size != ret)
      errexit("bdmq_timedrecv: Received message size %zu is not the same as "
        "requested size %zu.\n", ret, size);

    if (size < ret)
      ret = -1;
    else
      memcpy(buf, mq->buf, ret);
  }
  /*printf("[%d] timedrecv-out\n", (int)time(0));*/

  return ret;
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
