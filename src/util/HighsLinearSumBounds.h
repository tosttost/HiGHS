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
/**@file util/HighsLinearSumBounds.h
 * @brief Data structure to compute and update bounds on a linear sum of
 * variables with finite or infinite bounds
 */

#ifndef HIGHS_LINEAR_SUM_BOUNDS_H_
#define HIGHS_LINEAR_SUM_BOUNDS_H_

#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsCD0uble.h"

class HighsLinearSumBounds {
  std::vector<HighsCD0uble> sumLowerOrig;
  std::vector<HighsCD0uble> sumUpperOrig;
  std::vector<HighsInt> numInfSumLowerOrig;
  std::vector<HighsInt> numInfSumUpperOrig;
  std::vector<HighsCD0uble> sumLower;
  std::vector<HighsCD0uble> sumUpper;
  std::vector<HighsInt> numInfSumLower;
  std::vector<HighsInt> numInfSumUpper;
  const HighsFloat* varLower;
  const HighsFloat* varUpper;
  const HighsFloat* implVarLower;
  const HighsFloat* implVarUpper;
  const HighsInt* implVarLowerSource;
  const HighsInt* implVarUpperSource;

 public:
  void setNumSums(HighsInt numSums) {
    numInfSumLower.resize(numSums);
    numInfSumUpper.resize(numSums);
    sumLower.resize(numSums);
    sumUpper.resize(numSums);
    numInfSumLowerOrig.resize(numSums);
    numInfSumUpperOrig.resize(numSums);
    sumLowerOrig.resize(numSums);
    sumUpperOrig.resize(numSums);
  }

  void setBoundArrays(const HighsFloat* varLower, const HighsFloat* varUpper,
                      const HighsFloat* implVarLower, const HighsFloat* implVarUpper,
                      const HighsInt* implVarLowerSource,
                      const HighsInt* implVarUpperSource) {
    this->varLower = varLower;
    this->varUpper = varUpper;
    this->implVarLower = implVarLower;
    this->implVarUpper = implVarUpper;
    this->implVarLowerSource = implVarLowerSource;
    this->implVarUpperSource = implVarUpperSource;
  }

  void sumScaled(HighsInt sum, HighsFloat scale) {
    sumLowerOrig[sum] *= scale;
    sumUpperOrig[sum] *= scale;
    sumLower[sum] *= scale;
    sumUpper[sum] *= scale;

    if (scale < 0) {
      std::swap(sumLower[sum], sumUpper[sum]);
      std::swap(sumLowerOrig[sum], sumUpperOrig[sum]);
      std::swap(numInfSumLower[sum], numInfSumUpper[sum]);
      std::swap(numInfSumLowerOrig[sum], numInfSumUpperOrig[sum]);
    }
  }

  void add(HighsInt sum, HighsInt var, HighsFloat coefficient);

  void remove(HighsInt sum, HighsInt var, HighsFloat coefficient);

  void updatedVarUpper(HighsInt sum, HighsInt var, HighsFloat coefficient,
                       HighsFloat oldVarUpper);

  void updatedVarLower(HighsInt sum, HighsInt var, HighsFloat coefficient,
                       HighsFloat oldVarLower);

  void updatedImplVarUpper(HighsInt sum, HighsInt var, HighsFloat coefficient,
                           HighsFloat oldImplVarUpper,
                           HighsInt oldImplVarUpperSource);

  void updatedImplVarLower(HighsInt sum, HighsInt var, HighsFloat coefficient,
                           HighsFloat oldImplVarLower,
                           HighsInt oldImplVarLowerSource);

  HighsFloat getResidualSumLower(HighsInt sum, HighsInt var,
                             HighsFloat coefficient) const;

  HighsFloat getResidualSumUpper(HighsInt sum, HighsInt var,
                             HighsFloat coefficient) const;

  HighsFloat getResidualSumLowerOrig(HighsInt sum, HighsInt var,
                                 HighsFloat coefficient) const;

  HighsFloat getResidualSumUpperOrig(HighsInt sum, HighsInt var,
                                 HighsFloat coefficient) const;

  HighsFloat getSumLowerOrig(HighsInt sum) const {
    return numInfSumLowerOrig[sum] == 0 ? HighsFloat(sumLowerOrig[sum])
                                        : -kHighsInf;
  }

  HighsFloat getSumUpperOrig(HighsInt sum) const {
    return numInfSumUpperOrig[sum] == 0 ? HighsFloat(sumUpperOrig[sum]) : kHighsInf;
  }

  HighsFloat getSumLower(HighsInt sum) const {
    return numInfSumLower[sum] == 0 ? HighsFloat(sumLower[sum]) : -kHighsInf;
  }

  HighsFloat getSumUpper(HighsInt sum) const {
    return numInfSumUpper[sum] == 0 ? HighsFloat(sumUpper[sum]) : kHighsInf;
  }

  HighsFloat getSumLower(HighsInt sum, HighsFloat offset) const {
    return numInfSumLower[sum] == 0 ? HighsFloat(sumLower[sum] + offset)
                                    : -kHighsInf;
  }

  HighsFloat getSumUpper(HighsInt sum, HighsFloat offset) const {
    return numInfSumUpper[sum] == 0 ? HighsFloat(sumUpper[sum] + offset)
                                    : kHighsInf;
  }

  HighsFloat getSumLower(HighsInt sum, HighsCD0uble offset) const {
    return numInfSumLower[sum] == 0 ? HighsFloat(sumLower[sum] + offset)
                                    : -kHighsInf;
  }

  HighsFloat getSumUpper(HighsInt sum, HighsCD0uble offset) const {
    return numInfSumUpper[sum] == 0 ? HighsFloat(sumUpper[sum] + offset)
                                    : kHighsInf;
  }

  HighsInt getNumInfSumLower(HighsInt sum) const { return numInfSumLower[sum]; }

  HighsInt getNumInfSumUpper(HighsInt sum) const { return numInfSumUpper[sum]; }

  HighsInt getNumInfSumLowerOrig(HighsInt sum) const {
    return numInfSumLowerOrig[sum];
  }

  HighsInt getNumInfSumUpperOrig(HighsInt sum) const {
    return numInfSumUpperOrig[sum];
  }

  void shrink(const std::vector<HighsInt>& newIndices, HighsInt newSize);
};

#endif
