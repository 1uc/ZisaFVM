// SPDX-License-Identifier: MIT
// Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval

#if ZISA_HAS_OPENMP != 0
#include <omp.h>

#define ZISA_OMP_FOR_SCHEDULE_DEFAULT schedule(static, 8)

#endif
