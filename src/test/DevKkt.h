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
/**@file test/DevKkt.h
 * @brief
 */
#ifndef TEST_DEV_KKT_H_
#define TEST_DEV_KKT_H_

#include <map>
#include <vector>

#include "lp_data/HConst.h"

namespace presolve {
namespace dev_kkt_check {

struct State {
  State(
      const HighsInt numCol_, const HighsInt numRow_,
      const std::vector<HighsInt>& Astart_, const std::vector<HighsInt>& Aend_,
      const std::vector<HighsInt>& Aindex_, const std::vector<HighsFloat>& Avalue_,
      const std::vector<HighsInt>& ARstart_,
      const std::vector<HighsInt>& ARindex_,
      const std::vector<HighsFloat>& ARvalue_, const std::vector<HighsFloat>& colCost_,
      const std::vector<HighsFloat>& colLower_,
      const std::vector<HighsFloat>& colUpper_,
      const std::vector<HighsFloat>& rowLower_,
      const std::vector<HighsFloat>& rowUpper_,
      const std::vector<HighsInt>& flagCol_,
      const std::vector<HighsInt>& flagRow_,
      const std::vector<HighsFloat>& colValue_, const std::vector<HighsFloat>& colDual_,
      const std::vector<HighsFloat>& rowValue_, const std::vector<HighsFloat>& rowDual_,
      const std::vector<HighsBasisStatus>& col_status_,
      const std::vector<HighsBasisStatus>& row_status_)
      : numCol(numCol_),
        numRow(numRow_),
        Astart(Astart_),
        Aend(Aend_),
        Aindex(Aindex_),
        Avalue(Avalue_),
        ARstart(ARstart_),
        ARindex(ARindex_),
        ARvalue(ARvalue_),
        colCost(colCost_),
        colLower(colLower_),
        colUpper(colUpper_),
        rowLower(rowLower_),
        rowUpper(rowUpper_),
        flagCol(flagCol_),
        flagRow(flagRow_),
        colValue(colValue_),
        colDual(colDual_),
        rowValue(rowValue_),
        rowDual(rowDual_),
        col_status(col_status_),
        row_status(row_status_) {}

  const HighsInt numCol;
  const HighsInt numRow;

  const std::vector<HighsInt>& Astart;
  const std::vector<HighsInt>& Aend;
  const std::vector<HighsInt>& Aindex;
  const std::vector<HighsFloat>& Avalue;

  const std::vector<HighsInt>& ARstart;
  const std::vector<HighsInt>& ARindex;
  const std::vector<HighsFloat>& ARvalue;

  const std::vector<HighsFloat>& colCost;
  const std::vector<HighsFloat>& colLower;
  const std::vector<HighsFloat>& colUpper;
  const std::vector<HighsFloat>& rowLower;
  const std::vector<HighsFloat>& rowUpper;

  const std::vector<HighsInt>& flagCol;
  const std::vector<HighsInt>& flagRow;

  // solution
  const std::vector<HighsFloat>& colValue;
  const std::vector<HighsFloat>& colDual;
  const std::vector<HighsFloat>& rowValue;
  const std::vector<HighsFloat>& rowDual;

  // basis
  const std::vector<HighsBasisStatus>& col_status;
  const std::vector<HighsBasisStatus>& row_status;
};

enum class KktCondition {
  kColBounds,
  kPrimalFeasibility,
  kDualFeasibility,
  kComplementarySlackness,
  kStationarityOfLagrangian,
  kBasicFeasibleSolution,
  kUnset,
};

struct KktConditionDetails {
  KktConditionDetails() {}
  KktConditionDetails(KktCondition type_) : type(type_) {}

  KktCondition type = KktCondition::kUnset;
  HighsFloat max_violation = 0.0;
  HighsFloat sum_violation_2 = 0.0;
  HighsInt checked = 0;
  HighsInt violated = 0;
};

struct KktInfo {
  std::map<KktCondition, KktConditionDetails> rules;
  bool pass_col_bounds = false;
  bool pass_primal_feas_matrix = false;
  bool pass_dual_feas = false;
  bool pass_st_of_L = false;
  bool pass_comp_slackness = false;
  bool pass_bfs = false;
};

KktInfo initInfo();

bool checkKkt(const State& state, KktInfo& info);

void checkPrimalBounds(const State& state, KktConditionDetails& details);
void checkPrimalFeasMatrix(const State& state, KktConditionDetails& details);
void checkDualFeasibility(const State& state, KktConditionDetails& details);
void checkComplementarySlackness(const State& state,
                                 KktConditionDetails& details);
void checkStationarityOfLagrangian(const State& state,
                                   KktConditionDetails& details);
void checkBasicFeasibleSolution(const State& state,
                                KktConditionDetails& details);

}  // namespace dev_kkt_check
}  // namespace presolve

#endif /* TEST_KKTCHSTEP_H_ */
