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
#ifndef MIP_HIGHS_MIP_SOLVER_H_
#define MIP_HIGHS_MIP_SOLVER_H_

#include "Highs.h"
#include "lp_data/HighsOptions.h"

struct HighsMipSolverData;
class HighsCutPool;
struct HighsPseudocostInitialization;
class HighsCliqueTable;
class HighsImplications;

class HighsMipSolver {
 public:
  const HighsOptions* options_mip_;
  const HighsLp* model_;
  const HighsLp* orig_model_;
  HighsModelStatus modelstatus_;
  std::vector<HighsFloat> solution_;
  HighsFloat solution_objective_;
  HighsFloat bound_violation_;
  HighsFloat integrality_violation_;
  HighsFloat row_violation_;
  HighsFloat dual_bound_;
  HighsFloat primal_bound_;
  int64_t node_count_;

  bool submip;
  const HighsBasis* rootbasis;
  const HighsPseudocostInitialization* pscostinit;
  const HighsCliqueTable* clqtableinit;
  const HighsImplications* implicinit;

  std::unique_ptr<HighsMipSolverData> mipdata_;

  void run();

  HighsInt numCol() const { return model_->num_col_; }

  HighsInt numRow() const { return model_->num_row_; }

  HighsInt numNonzero() const { return model_->a_matrix_.numNz(); }

  const HighsFloat* colCost() const { return model_->col_cost_.data(); }

  HighsFloat colCost(HighsInt col) const { return model_->col_cost_[col]; }

  const HighsFloat* rowLower() const { return model_->row_lower_.data(); }

  HighsFloat rowLower(HighsInt col) const { return model_->row_lower_[col]; }

  const HighsFloat* rowUpper() const { return model_->row_upper_.data(); }

  HighsFloat rowUpper(HighsInt col) const { return model_->row_upper_[col]; }

  bool isSolutionFeasible(const std::vector<HighsFloat>& solution) const;

  const HighsVarType* variableType() const {
    return model_->integrality_.data();
  }

  HighsVarType variableType(HighsInt col) const {
    return model_->integrality_[col];
  }

  HighsMipSolver(const HighsOptions& options, const HighsLp& lp,
                 const HighsSolution& solution, bool submip = false);

  ~HighsMipSolver();

  void setModel(const HighsLp& model) {
    model_ = &model;
    solution_objective_ = kHighsInf;
  }

  mutable HighsTimer timer_;
  void cleanupSolve();
};

#endif
