#ifndef _CHECK_PTR_H_
#define _CHECK_PTR_H_

#include <stdio.h>

#define check_ptr(ptr)\
  if(ptr == 0) {\
    fprintf(stderr, "Could not allocate memory successfully %s, %d",\
            __FILE__, __LINE__);\
    exit(1);\
  }

#endif // _CHECK_PTR_H_
