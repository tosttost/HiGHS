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
#ifndef HIGHS_LP_PROPAGATOR_H_
#define HIGHS_LP_PROPAGATOR_H_

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsCD0uble.h"

/// propagates domains as part of LP presolve
/// final propagated bounds are relaxed by a wide enough margin
/// so that they cannot be used in any basic feasible solution
class HighsLpPropagator {
  std::vector<HighsCD0uble> activitymin_;
  std::vector<HighsCD0uble> activitymax_;
  std::vector<HighsInt> activitymininf_;
  std::vector<HighsInt> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<HighsInt> propagateinds_;

  std::vector<HighsFloat>& Avalue_;
  std::vector<HighsInt>& Aindex_;
  std::vector<HighsInt>& Astart_;
  std::vector<HighsInt>& Aend_;

  std::vector<HighsFloat>& ARvalue_;
  std::vector<HighsInt>& ARindex_;
  std::vector<HighsInt>& ARstart_;

  const std::vector<HighsInt>& flagRow;
  const std::vector<HighsInt>& flagCol;
  std::vector<HighsFloat>& rowLower_;
  std::vector<HighsFloat>& rowUpper_;
  const std::vector<HighsVarType>& integrality_;

  HighsInt infeasible_ = 0;
  HighsInt numBoundChgs_ = 0;

  void computeMinActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const HighsFloat* ARvalue, HighsInt& ninfmin,
                          HighsCD0uble& activitymin);

  void computeMaxActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const HighsFloat* ARvalue, HighsInt& ninfmax,
                          HighsCD0uble& activitymax);

  HighsInt propagateRowUpper(const HighsInt* Rindex, const HighsFloat* Rvalue,
                             HighsInt Rlen, HighsFloat Rupper,
                             const HighsCD0uble& minactivity, HighsInt ninfmin,
                             HighsDomainChange* boundchgs);

  HighsInt propagateRowLower(const HighsInt* Rindex, const HighsFloat* Rvalue,
                             HighsInt Rlen, HighsFloat Rlower,
                             const HighsCD0uble& maxactivity, HighsInt ninfmax,
                             HighsDomainChange* boundchgs);

  void updateActivityLbChange(HighsInt col, HighsFloat oldbound, HighsFloat newbound);

  void updateActivityUbChange(HighsInt col, HighsFloat oldbound, HighsFloat newbound);

  HighsFloat doChangeBound(const HighsDomainChange& boundchg);

 public:
  std::vector<HighsFloat> colLower_;
  std::vector<HighsFloat> colUpper_;

  HighsLpPropagator(
      const std::vector<HighsFloat>& colLower, const std::vector<HighsFloat>& colUpper,
      const std::vector<HighsVarType>& integrality_,
      std::vector<HighsFloat>& Avalue_, std::vector<HighsInt>& Aindex_,
      std::vector<HighsInt>& Astart_, std::vector<HighsInt>& Aend_,
      std::vector<HighsFloat>& ARvalue_, std::vector<HighsInt>& ARindex_,
      std::vector<HighsInt>& ARstart_, const std::vector<HighsInt>& flagRow,
      const std::vector<HighsInt>& flagCol, std::vector<HighsFloat>& rowLower_,
      std::vector<HighsFloat>& rowUpper_);

  void markPropagate(HighsInt row);

  void computeRowActivities();

  bool infeasible() const { return infeasible_ != 0; }

  void changeBound(HighsDomainChange boundchg);

  HighsInt propagate();

  HighsInt tightenCoefficients();

  HighsInt getNumChangedBounds() const { return numBoundChgs_; }
};

#endif
