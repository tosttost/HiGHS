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
/**@file model/HighsModel.h
 * @brief
 */
#ifndef MODEL_HIGHS_MODEL_H_
#define MODEL_HIGHS_MODEL_H_

#include <vector>

#include "lp_data/HighsLp.h"
#include "model/HighsHessian.h"

class HighsModel;

class HighsModel {
 public:
  HighsLp lp_;
  HighsHessian hessian_;
  bool isQp() const { return this->hessian_.dim_; }
  bool isMip() const { return this->lp_.isMip(); }
  void clear();
  HighsFloat objectiveValue(const std::vector<HighsFloat>& solution) const;
  void objectiveGradient(const std::vector<HighsFloat>& solution,
                         std::vector<HighsFloat>& gradient) const;
};

#endif
