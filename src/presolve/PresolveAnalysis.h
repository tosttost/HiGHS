/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_ANALYSIS_H_
#define PRESOLVE_PRESOLVE_ANALYSIS_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsTimer.h"

namespace presolve {

using std::max;
using std::min;

constexpr double inf = std::numeric_limits<double>::infinity();

enum PresolveRule {
  // Presolve rules.
  EMPTY_ROW,
  FIXED_COL,
  SING_ROW,
  DOUBLETON_EQUATION,
  REMOVE_FORCING_CONSTRAINTS,
  FORCING_ROW,
  REDUNDANT_ROW,
  DOMINATED_ROW_BOUNDS,
  REMOVE_COLUMN_SINGLETONS,
  FREE_SING_COL,
  SING_COL_DOUBLETON_INEQ,
  IMPLIED_FREE_SING_COL,
  REMOVE_DOMINATED_COLUMNS,
  DOMINATED_COLS,
  WEAKLY_DOMINATED_COLS,
  DOMINATED_COL_BOUNDS,
  EMPTY_COL,
  // HTICK_PRE_DUPLICATE_ROWS,
  // HTICK_PRE_DUPLICATE_COLUMNS,

  // For timing.
  MATRIX_COPY,
  RESIZE_MATRIX,

  RUN_PRESOLVERS,
  REMOVE_ROW_SINGLETONS,
  REMOVE_DOUBLETON_EQUATIONS,
  REMOVE_EMPTY_ROW,

  TOTAL_PRESOLVE_TIME,
  // Number of presolve rules.
  PRESOLVE_RULES_COUNT,

  // Items required by postsolve
  DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE,
  DOUBLETON_EQUATION_X_ZERO_INITIALLY,
  DOUBLETON_EQUATION_NEW_X_NONZERO,
  DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE,
  DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE,
  SING_COL_DOUBLETON_INEQ_SECOND_SING_COL,
  FORCING_ROW_VARIABLE
};

enum presolveNumerics {
  PRESOLVE_INCONSISTENT_BOUNDS,
  PRESOLVE_FIXED_COLUMN,
  PRESOLVE_DOUBLETON_EQUATION_BOUND,
  PRESOLVE_DOUBLETON_INEQUALITY_BOUND,
  PRESOLVE_SMALL_MATRIX_VALUE,
  PRESOLVE_EMPTY_ROW_BOUND,
  PRESOLVE_DOMINATED_COLUMN,
  PRESOLVE_WEAKLY_DOMINATED_COLUMN,
  PRESOLVE_NUMERICS_COUNT
};

enum postsolveNumerics {
  POSTSOLVE_INCONSISTENT_BOUNDS,
  POSTSOLVE_DOUBLETON_INEQUALITY_BASIC,
  POSTSOLVE_DOUBLETON_INEQUALITY_LO_UP,
  POSTSOLVE_DOUBLETON_INEQUALITY_INFEAS,
  POSTSOLVE_BOUND_ON_LBYZJ_FIXED,
  POSTSOLVE_BOUND_ON_LBYZJ_AT_BOUND,
  POSTSOLVE_DUALS_SING_ROW_COL_BASIC,
  POSTSOLVE_DUALS_SING_ROW_ROW_BASIC,
  POSTSOLVE_DUALS_SING_ROW_ROW_BELOW_LB,
  POSTSOLVE_DUALS_SING_ROW_ROW_ABOVE_UB,
  POSTSOLVE_DUALS_SING_ROW_ROW_DUAL,
  POSTSOLVE_DUALS_DOUBLETON_EQUALITY_X_VALUE,
  POSTSOLVE_DUALS_DOUBLETON_EQUALITY_X_BOUND,
  POSTSOLVE_DUALS_DOUBLETON_EQUALITY_Y0,
  POSTSOLVE_DUALS_DOUBLETON_EQUALITY_Y1,
  POSTSOLVE_NUMERICS_COUNT
};

struct PresolveRuleInfo {
  PresolveRuleInfo(PresolveRule id, std::string name, std::string name_ch3)
      : rule_id(id),
        rule_name(std::move(name)),
        rule_name_ch3(std::move(name_ch3)) {}
  PresolveRule rule_id;

  std::string rule_name;
  std::string rule_name_ch3;

  int count_applied = 0;
  int rows_removed = 0;
  int cols_removed = 0;

  int clock_id = 0;
  double total_time = 0;
};

struct numericsRecord {
  std::string name;
  double tolerance;
  int num_test;
  int num_negative;
  int num_zero_true;
  int num_tol10_true;
  int num_tol_true;
  int num_10tol_true;
  int num_clear_true;
  double max_below_tolerance;
  double min_above_tolerance;
};

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules);

class PresolveTimer {
 public:
  PresolveTimer(HighsTimer& timer) : timer_(timer) {
    initializePresolveRuleInfo(rules_);
    for (PresolveRuleInfo& rule : rules_) {
      int clock_id =
          timer_.clock_def(rule.rule_name.c_str(), rule.rule_name_ch3.c_str());
      rule.clock_id = clock_id;
    }
  }

  std::vector<numericsRecord> presolve_numerics;
  std::vector<numericsRecord> postsolve_numerics;

  void recordStart(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    timer_.start(rules_[rule].clock_id);
  }

  void recordFinish(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    timer_.stop(rules_[rule].clock_id);

    if (rule == TOTAL_PRESOLVE_TIME)
      total_time_ = timer_.read(rules_[rule].clock_id);
  }

  void addChange(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    rules_[rule].count_applied++;
  }

  void increaseCount(bool row_count, PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    if (row_count)
      rules_[rule].rows_removed++;
    else
      rules_[rule].cols_removed++;
  }

  void reportClocks() {
    std::vector<int> clocks;
    for (int id = 0; id < PRESOLVE_RULES_COUNT - 1; id++) {
      assert(rules_[id].rule_id == id);
      if (id == RUN_PRESOLVERS) continue;
      if (id == REMOVE_ROW_SINGLETONS) continue;
      if (id == REMOVE_DOUBLETON_EQUATIONS) continue;
      if (id == REMOVE_EMPTY_ROW) continue;
      clocks.push_back(rules_[id].clock_id);
    }
    int ideal_time_rule;
    double ideal_time;
    ideal_time_rule = TOTAL_PRESOLVE_TIME;
    ideal_time = getRuleTime(ideal_time_rule);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    clocks.push_back(rules_[RUN_PRESOLVERS].clock_id);
    clocks.push_back(rules_[RESIZE_MATRIX].clock_id);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = RUN_PRESOLVERS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[REMOVE_ROW_SINGLETONS].clock_id);
    clocks.push_back(rules_[REMOVE_FORCING_CONSTRAINTS].clock_id);
    clocks.push_back(rules_[REMOVE_COLUMN_SINGLETONS].clock_id);
    clocks.push_back(rules_[REMOVE_DOUBLETON_EQUATIONS].clock_id);
    clocks.push_back(rules_[REMOVE_DOMINATED_COLUMNS].clock_id);
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_FORCING_CONSTRAINTS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[REMOVE_EMPTY_ROW].clock_id);
    clocks.push_back(rules_[FORCING_ROW].clock_id);
    clocks.push_back(rules_[REDUNDANT_ROW].clock_id);
    clocks.push_back(rules_[DOMINATED_ROW_BOUNDS].clock_id);
    timer_.report_tl("grep--RmFrcCs", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_COLUMN_SINGLETONS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[FREE_SING_COL].clock_id);
    clocks.push_back(rules_[SING_COL_DOUBLETON_INEQ].clock_id);
    clocks.push_back(rules_[IMPLIED_FREE_SING_COL].clock_id);
    timer_.report_tl("grep-RmColSng", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_DOMINATED_COLUMNS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[DOMINATED_COLS].clock_id);
    clocks.push_back(rules_[WEAKLY_DOMINATED_COLS].clock_id);
    timer_.report_tl("grep-RmDomCol", clocks, ideal_time, 0);
    std::cout << std::endl;
  }

  void initialiseNumericsRecord(std::string name, const double tolerance,
                                numericsRecord& numerics_record) {
    // Make sure that the tolerance has been set to a positive value
    assert(tolerance > 0);
    numerics_record.name = name;
    numerics_record.tolerance = tolerance;
    numerics_record.num_test = 0;
    numerics_record.num_negative = 0;
    numerics_record.num_zero_true = 0;
    numerics_record.num_tol10_true = 0;
    numerics_record.num_tol_true = 0;
    numerics_record.num_10tol_true = 0;
    numerics_record.num_clear_true = 0;
    numerics_record.max_below_tolerance = 0;
    numerics_record.min_above_tolerance = HIGHS_CONST_INF;
  }

  void initialisePresolveNumericsRecord(int record, std::string name,
                                        const double tolerance) {
    initialiseNumericsRecord(name, tolerance, presolve_numerics[record]);
  }

  void initialisePostsolveNumericsRecord(int record, std::string name,
                                         const double tolerance) {
    initialiseNumericsRecord(name, tolerance, postsolve_numerics[record]);
  }

  void updateNumericsRecord(numericsRecord& numerics_record,
                            const double value) {
    double tolerance = numerics_record.tolerance;
    numerics_record.num_test++;
    if (value < 0) {
      numerics_record.num_negative++;
      return;
    }
    if (value == 0) {
      numerics_record.num_zero_true++;
    } else if (10 * value <= tolerance) {
      numerics_record.num_tol10_true++;
    } else if (value <= tolerance) {
      numerics_record.num_tol_true++;
    } else if (value <= 10 * tolerance) {
      numerics_record.num_10tol_true++;
    } else {
      numerics_record.num_clear_true++;
    }
    if (value <= tolerance)
      numerics_record.max_below_tolerance =
          max(value, numerics_record.max_below_tolerance);
    if (value >= tolerance)
      numerics_record.min_above_tolerance =
          min(value, numerics_record.min_above_tolerance);
  }

  void updatePresolveNumericsRecord(int record, const double value) {
    updateNumericsRecord(presolve_numerics[record], value);
  }

  void updatePostsolveNumericsRecord(int record, const double value) {
    updateNumericsRecord(postsolve_numerics[record], value);
  }

  void reportNumericsRecord(const numericsRecord& numerics_record) {
    if (!numerics_record.num_test) return;
    assert(numerics_record.num_test ==
           numerics_record.num_negative + numerics_record.num_zero_true +
               numerics_record.num_tol10_true + numerics_record.num_tol_true +
               numerics_record.num_10tol_true + numerics_record.num_clear_true);
    printf("%-30s %6.1g %9d %9d %9d %9.2g %9d %9.2g %9d %9d\n",
           numerics_record.name.c_str(), numerics_record.tolerance,
           numerics_record.num_negative, numerics_record.num_zero_true,
           numerics_record.num_tol10_true, numerics_record.max_below_tolerance,
           numerics_record.num_tol_true, numerics_record.min_above_tolerance,
           numerics_record.num_10tol_true, numerics_record.num_clear_true);
  }

  void reportNumericsCsvRecord(const numericsRecord& numerics_record,
                               const bool header = false) {
    if (header) {
      printf(",Tol/10,MaxBw,Tol,MinAb,10Tol");
    } else {
      printf(",");
      if (numerics_record.num_tol10_true)
        printf("%d", numerics_record.num_tol10_true);
      printf(",");
      if (numerics_record.max_below_tolerance > 0)
        printf("%g", numerics_record.max_below_tolerance);
      printf(",");
      if (numerics_record.num_tol_true)
        printf("%d", numerics_record.num_tol_true);
      printf(",");
      if (numerics_record.min_above_tolerance < HIGHS_CONST_INF)
        printf("%g", numerics_record.min_above_tolerance);
      printf(",");
      if (numerics_record.num_10tol_true)
        printf("%d", numerics_record.num_10tol_true);
    }
  }

  void reportNumericsRecords(
      const std::string type,
      const std::vector<numericsRecord>& numerics_record) {
    int num_record = numerics_record.size();
    printf("%s numerics analysis for %s\n", type.c_str(), model_name.c_str());
    printf(
        "Rule                              Tol  Negative      Zero    Tol/10   "
        "  "
        "MaxBw       Tol     MinAb     10Tol     Clear\n");

    for (int record = 0; record < num_record; record++)
      reportNumericsRecord(numerics_record[record]);
    printf("grep_%sNumerics:,0,%s", type.c_str(), model_name.c_str());
    for (int record = 0; record < num_record; record++)
      printf(",%s,,", numerics_record[record].name.c_str());
    printf("\n");

    printf("grep_%sNumerics:,1,%s", type.c_str(), model_name.c_str());
    for (int record = 0; record < num_record; record++)
      reportNumericsCsvRecord(numerics_record[record], true);
    printf("\n");

    printf("grep_%sNumerics:,2,%s", type.c_str(), model_name.c_str());
    for (int record = 0; record < num_record; record++)
      reportNumericsCsvRecord(numerics_record[record]);
    printf("\n");
  }

  void reportPresolveNumericsRecords() {
    reportNumericsRecords("Presolve", presolve_numerics);
  }
  void reportPostsolveNumericsRecords() {
    reportNumericsRecords("Postsolve", postsolve_numerics);
  }

  void updateInfo();
  double getTotalTime() { return total_time_; }

  HighsTimer& timer_;

  double getRuleTime(const int rule_id) {
    return timer_.read(rules_[rule_id].clock_id);
  }

  inline double getTime() { return timer_.readRunHighsClock(); }

  inline bool reachLimit() {
    if (time_limit == inf || time_limit <= 0) return false;
    if (getTime() < time_limit) return false;
    return true;
  }

  double start_time = 0.0;
  double time_limit = 0.0;
  std::string model_name;

 private:
  std::vector<PresolveRuleInfo> rules_;

  double total_time_ = 0.0;
};

}  // namespace presolve

#endif
