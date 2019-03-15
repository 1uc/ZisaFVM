#if ZISA_HAS_OPENMP != 0
#include <omp.h>

#define ZISA_OMP_FOR_SCHEDULE_DEFAULT schedule(dynamic, 64)

#endif
