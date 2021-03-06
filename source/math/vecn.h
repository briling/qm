#include <stdlib.h>
#include <math.h>

void vecset (size_t n, double * u, double s);
void veccp  (size_t n, double * u, double * v);
void vecsum (size_t n, double * w, double * u, double * v);
void vecscal(size_t n, double * u, double s);
void vecadd (size_t n, double * u, double * v);
void vecadds(size_t n, double * u, double * v, double s);
double vecabsmax(size_t n, double * u);

