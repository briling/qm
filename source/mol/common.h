#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define   BA    0.5291772
#define   AB    1.88972616356109068947

#define MAX(x,y) ( ((x) > (y)) ? (x) : (y) )
#define MIN(x,y) ( ((x) < (y)) ? (x) : (y) )

#define GOTOHELL { \
      fprintf(stderr, "%s:%d %s() -- ", \
      __FILE__, __LINE__, __FUNCTION__); \
      abort(); }

#define printalive printf("alive @ %s:%d\n", __FILE__, __LINE__);

#define fprintf_if(F, ...) if(F) fprintf(F, __VA_ARGS__);

typedef char styp[8];
void     g_print(int n, double * G, const char s[], FILE * f);

