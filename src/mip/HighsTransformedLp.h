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
/**@file mip/HighsTransformedLp.h
 * @brief LP transformations useful for cutting plane separation. This includes
 * bound substitution with simple and variable bounds, handling of slack
 * variables, flipping the complementation of integers.
 */

#ifndef MIP_HIGHS_TRANSFORMED_LP_H_
#define MIP_HIGHS_TRANSFORMED_LP_H_

#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsImplications.h"
#include "mip/HighsSparseVectorSum.h"
#include "util/HighsCD0uble.h"
#include "util/HighsInt.h"

class HighsLpRelaxation;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTransformedLp {
 private:
  const HighsLpRelaxation& lprelaxation;

  std::vector<const std::pair<const HighsInt, HighsImplications::VarBound>*>
      bestVub;
  std::vector<const std::pair<const HighsInt, HighsImplications::VarBound>*>
      bestVlb;
  std::vector<HighsFloat> simpleLbDist;
  std::vector<HighsFloat> simpleUbDist;
  std::vector<HighsFloat> lbDist;
  std::vector<HighsFloat> ubDist;
  std::vector<HighsFloat> boundDist;
  enum class BoundType : uint8_t {
    kSimpleUb,
    kSimpleLb,
    kVariableUb,
    kVariableLb,
  };
  std::vector<BoundType> boundTypes;
  HighsSparseVectorSum vectorsum;

 public:
  HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                     HighsImplications& implications);

  HighsFloat boundDistance(HighsInt col) const { return boundDist[col]; }

  bool transform(std::vector<HighsFloat>& vals, std::vector<HighsFloat>& upper,
                 std::vector<HighsFloat>& solval, std::vector<HighsInt>& inds,
                 HighsFloat& rhs, bool& integralPositive, bool preferVbds = false);

  bool untransform(std::vector<HighsFloat>& vals, std::vector<HighsInt>& inds,
                   HighsFloat& rhs, bool integral = false);
};

#endif
