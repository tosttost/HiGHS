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
/**@file simplex/HEkkInterface.cpp
 * @brief
 */

#include "lp_data/HighsLpUtils.h"
#include "simplex/HEkk.h"

void HEkk::appendColsToVectors(const HighsInt num_new_col,
                               const vector<HighsFloat>& colCost,
                               const vector<HighsFloat>& colLower,
                               const vector<HighsFloat>& colUpper) {
  appendColsToLpVectors(lp_, num_new_col, colCost, colLower, colUpper);
}

void HEkk::appendRowsToVectors(const HighsInt num_new_row,
                               const vector<HighsFloat>& rowLower,
                               const vector<HighsFloat>& rowUpper) {
  appendRowsToLpVectors(lp_, num_new_row, rowLower, rowUpper);
}
