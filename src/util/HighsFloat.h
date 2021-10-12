/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsFloat.h
 * @brief The definition for the float type to use
 */

#ifndef UTIL_HIGHS_FLOAT_H_
#define UTIL_HIGHS_FLOAT_H_

#include <util/HighsCD0uble.h>

// Settings for double precision
// typedef double HighsFloat;
// const HighsFloat kHighsFloatTiny = 1e-14; // = kHighsTiny;

// Settings for psudo-quad precision
typedef HighsCD0uble HighsFloat;
const HighsFloat kHighsFloatTiny = 1e-30;

#endif
