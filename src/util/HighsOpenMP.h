/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsOpenMP.h
 * @brief Wrapper for Open MP methods that may not be available
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSOPENMP_H_
#define UTIL_HIGHSOPENMP_H_

#ifdef OPENMP
#include "omp.h"
#endif

int ompGetMaxThreads() {
  int omp_max_threads = 1;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
  assert(omp_max_threads > 0);
#ifdef HiGHSDEV
  if (omp_max_threads <= 0)
    printf("WARNING: omp_get_max_threads() returns %d\n", omp_max_threads);
  printf("Running with %d OMP thread(s)\n", omp_max_threads);
#endif
#endif
  return omp_max_threads;
}

int ompGetThreadNum() {
  int omp_thread_num = 0;
#ifdef OPENMP
  omp_thread_num = omp_get_thread_num();
  assert(omp_thread_num >= 0);
#ifdef HiGHSDEV
  if (omp_thread_num < 0)
    printf("WARNING: omp_get_thread_num() returns %d\n", omp_thread_num);
#endif
#endif
  return omp_thread_num;
}

#endif /* UTIL_HIGHSOPENMP_H_ */
