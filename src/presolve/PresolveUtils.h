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
#ifndef PRESOLVE_PRESOLVE_UTILS_H_
#define PRESOLVE_PRESOLVE_UTILS_H_

#include <vector>

#include "util/HighsInt.h"

namespace presolve {

void printRowwise(
    const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsFloat>& colCost, const std::vector<HighsFloat>& colLower,
    const std::vector<HighsFloat>& colUpper, const std::vector<HighsFloat>& rowLower,
    const std::vector<HighsFloat>& rowUpper, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<HighsFloat>& ARvalue);

void printRow(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<HighsFloat>& rowLower, const std::vector<HighsFloat>& rowUpper,
    const std::vector<HighsFloat>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<HighsFloat>& ARvalue);

void printCol(
    const HighsInt col, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<HighsFloat>& colLower, const std::vector<HighsFloat>& colUpper,
    const std::vector<HighsFloat>& values, const std::vector<HighsInt>& Astart,
    const std::vector<HighsInt>& Aend, const std::vector<HighsInt>& Aindex,
    const std::vector<HighsFloat>& Avalue);

void printCol(
    const HighsInt col, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<HighsFloat>& colLower, const std::vector<HighsFloat>& colUpper,
    const std::vector<HighsFloat>& values, const std::vector<HighsInt>& Astart,
    const std::vector<HighsInt>& Aindex, const std::vector<HighsFloat>& Avalue);

void printRowOneLine(
    const HighsInt row, const HighsInt numRow, const HighsInt numCol,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    const std::vector<HighsFloat>& rowLower, const std::vector<HighsFloat>& rowUpper,
    const std::vector<HighsFloat>& values, const std::vector<HighsInt>& ARstart,
    const std::vector<HighsInt>& ARindex, const std::vector<HighsFloat>& ARvalue);

void printAR(const HighsInt numRow, const HighsInt numCol,
             const std::vector<HighsFloat>& colCost,
             const std::vector<HighsFloat>& rowLower,
             const std::vector<HighsFloat>& rowUpper,
             const std::vector<HighsInt>& ARstart,
             const std::vector<HighsInt>& ARindex,
             std::vector<HighsFloat>& ARvalue);

void printA(const HighsInt numRow, const HighsInt numCol,
            const std::vector<HighsFloat>& colCost,
            const std::vector<HighsFloat>& rowLower,
            const std::vector<HighsFloat>& rowUpper,
            const std::vector<HighsFloat>& colLower,
            const std::vector<HighsFloat>& colUpper,
            const std::vector<HighsInt>& Astart,
            const std::vector<HighsInt>& Aindex, std::vector<HighsFloat>& Avalue);

}  // namespace presolve

#endif
