#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <stdio.h>

//#ifndef NDEBUG
//#define printd(args...)\
//  fprintf(stderr, "%s:%d: ",__FILE__,__LINE__);\
//  fprintf(stderr, args)
//#endif

#ifndef NDEBUG
#define printd(args...)\
  fprintf(stderr, args)
#endif

#ifdef NDEBUG
#define printd(args...)
#endif


#endif // _DEBUG_H_
