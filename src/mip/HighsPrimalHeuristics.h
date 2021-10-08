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
#ifndef HIGHS_PRIMAL_HEURISTICS_H_
#define HIGHS_PRIMAL_HEURISTICS_H_

#include <vector>

#include "lp_data/HStruct.h"
#include "lp_data/HighsLp.h"
#include "util/HighsRandom.h"

class HighsMipSolver;

class HighsPrimalHeuristics {
 private:
  HighsMipSolver& mipsolver;
  size_t lp_iterations;

  HighsFloat successObservations;
  HighsInt numSuccessObservations;
  HighsFloat infeasObservations;
  HighsInt numInfeasObservations;

  HighsRandom randgen;

  std::vector<HighsInt> intcols;

 public:
  HighsPrimalHeuristics(HighsMipSolver& mipsolver);

  void setupIntCols();

  bool solveSubMip(const HighsLp& lp, const HighsBasis& basis,
                   HighsFloat fixingRate, std::vector<HighsFloat> colLower,
                   std::vector<HighsFloat> colUpper, HighsInt maxleaves,
                   HighsInt maxnodes, HighsInt stallnodes);

  HighsFloat determineTargetFixingRate();

  void rootReducedCost();

  void RENS(const std::vector<HighsFloat>& relaxationsol);

  void RINS(const std::vector<HighsFloat>& relaxationsol);

  void feasibilityPump();

  void centralRounding();

  void flushStatistics();

  bool tryRoundedPoint(const std::vector<HighsFloat>& point, char source);

  bool linesearchRounding(const std::vector<HighsFloat>& point1,
                          const std::vector<HighsFloat>& point2, char source);

  void randomizedRounding(const std::vector<HighsFloat>& relaxationsol);
};

#endif
