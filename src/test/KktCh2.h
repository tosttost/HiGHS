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
/**@file test/KktChStep.h
 * @brief
 */
#ifndef TEST_KKTCH2_H_
#define TEST_KKTCH2_H_

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stack>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "test/DevKkt.h"
#include "util/HighsInt.h"

namespace presolve {

namespace dev_kkt_check {

class KktCheck;

class KktChStep {
 public:
  KktChStep() {}
  virtual ~KktChStep() {}

  std::vector<HighsFloat> RcolCost;
  std::vector<HighsFloat> RcolLower;
  std::vector<HighsFloat> RcolUpper;
  std::vector<HighsFloat> RrowLower;
  std::vector<HighsFloat> RrowUpper;

  int print = 1;

  std::stack<std::vector<std::pair<HighsInt, HighsFloat> > > rLowers;
  std::stack<std::vector<std::pair<HighsInt, HighsFloat> > > rUppers;
  std::stack<std::vector<std::pair<HighsInt, HighsFloat> > > cLowers;
  std::stack<std::vector<std::pair<HighsInt, HighsFloat> > > cUppers;
  std::stack<std::vector<std::pair<HighsInt, HighsFloat> > > costs;

  // full matrix
  void setBoundsCostRHS(const std::vector<HighsFloat>& colUpper_,
                        const std::vector<HighsFloat>& colLower_,
                        const std::vector<HighsFloat>& cost,
                        const std::vector<HighsFloat>& rowLower_,
                        const std::vector<HighsFloat>& rowUpper_);
  void addChange(int type, HighsInt row, HighsInt col, HighsFloat valC,
                 HighsFloat dualC, HighsFloat dualR);
  void addCost(HighsInt col, HighsFloat value);

  dev_kkt_check::State initState(
      const HighsInt numCol_, const HighsInt numRow_,
      const std::vector<HighsInt>& Astart_, const std::vector<HighsInt>& Aend_,
      const std::vector<HighsInt>& Aindex_, const std::vector<HighsFloat>& Avalue_,
      const std::vector<HighsInt>& ARstart_,
      const std::vector<HighsInt>& ARindex_,
      const std::vector<HighsFloat>& ARvalue_,
      const std::vector<HighsInt>& flagCol_,
      const std::vector<HighsInt>& flagRow_,
      const std::vector<HighsFloat>& colValue_, const std::vector<HighsFloat>& colDual_,
      const std::vector<HighsFloat>& rowValue_, const std::vector<HighsFloat>& rowDual_,
      const std::vector<HighsBasisStatus>& col_status_,
      const std::vector<HighsBasisStatus>& row_status_);
};

}  // namespace dev_kkt_check

}  // namespace presolve
#endif /* TEST_KKTCHSTEP_H_ */
