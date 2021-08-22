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

#include "mip/HighsTransformedLp.h"

#include "mip/HighsMipSolverData.h"
#include "pdqsort/pdqsort.h"
#include "util/HighsCDouble.h"

HighsTransformedLp::HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                                       HighsCutPool& cutpool,
                                       HighsImplications& implications)
    : lprelaxation(lprelaxation), cutpool(cutpool) {
  assert(lprelaxation.scaledOptimal(lprelaxation.getStatus()));
  const HighsMipSolver& mipsolver = implications.mipsolver;
  const HighsSolution& lpSolution = lprelaxation.getLpSolver().getSolution();

  HighsInt numTransformedCol = lprelaxation.numCols() + lprelaxation.numRows();

  boundDist.resize(numTransformedCol);
  simpleLbDist.resize(numTransformedCol);
  simpleUbDist.resize(numTransformedCol);
  lbDist.resize(numTransformedCol);
  ubDist.resize(numTransformedCol);
  bestVlb.resize(numTransformedCol);
  bestVub.resize(numTransformedCol);
  boundTypes.resize(numTransformedCol);

  for (HighsInt col : mipsolver.mipdata_->continuous_cols) {
    mipsolver.mipdata_->implications.cleanupVarbounds(col);
    if (mipsolver.mipdata_->domain.infeasible()) return;

    if (mipsolver.mipdata_->domain.isFixed(col)) continue;

    double bestub = mipsolver.mipdata_->domain.col_upper_[col];
    simpleUbDist[col] = bestub - lpSolution.col_value[col];
    if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
      simpleUbDist[col] = 0.0;

    double minbestub = bestub;
    size_t bestvubnodes = 0;

    double bestlb = mipsolver.mipdata_->domain.col_lower_[col];
    simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
    if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
      simpleLbDist[col] = 0.0;
    double maxbestlb = bestlb;
    size_t bestvlbnodes = 0;

    for (const auto& vub : implications.getVUBs(col)) {
      if (vub.second.coef == kHighsInf) continue;
      if (mipsolver.mipdata_->domain.isFixed(vub.first)) continue;
      assert(mipsolver.mipdata_->domain.isBinary(vub.first));
      double vubval = lpSolution.col_value[vub.first] * vub.second.coef +
                      vub.second.constant;

      assert(vub.first >= 0 && vub.first < mipsolver.numCol());
      if (vubval <= bestub + mipsolver.mipdata_->feastol) {
        size_t vubnodes =
            vub.second.coef > 0
                ? mipsolver.mipdata_->nodequeue.numNodesDown(vub.first)
                : mipsolver.mipdata_->nodequeue.numNodesUp(vub.first);
        double minvubval = vub.second.minValue();
        if (bestVub[col] == nullptr || vubnodes > bestvubnodes ||
            (vubnodes == bestvubnodes &&
             minvubval < minbestub - mipsolver.mipdata_->feastol)) {
          bestub = vubval;
          minbestub = minvubval;
          bestVub[col] = &vub;
          bestvubnodes = vubnodes;
        }
      }
    }

    for (const auto& vlb : implications.getVLBs(col)) {
      if (vlb.second.coef == -kHighsInf) continue;
      if (mipsolver.mipdata_->domain.isFixed(vlb.first)) continue;
      assert(mipsolver.mipdata_->domain.isBinary(vlb.first));
      assert(vlb.first >= 0 && vlb.first < mipsolver.numCol());
      double vlbval = lpSolution.col_value[vlb.first] * vlb.second.coef +
                      vlb.second.constant;

      if (vlbval >= bestlb - mipsolver.mipdata_->feastol) {
        size_t vlbnodes =
            vlb.second.coef > 0
                ? mipsolver.mipdata_->nodequeue.numNodesUp(vlb.first)
                : mipsolver.mipdata_->nodequeue.numNodesDown(vlb.first);
        double maxvlbval = vlb.second.maxValue();
        if (bestVlb[col] == nullptr || vlbnodes > bestvlbnodes ||
            (vlbnodes == bestvlbnodes &&
             maxvlbval > maxbestlb + mipsolver.mipdata_->feastol)) {
          bestlb = vlbval;
          maxbestlb = maxvlbval;
          bestVlb[col] = &vlb;
          bestvlbnodes = vlbnodes;
        }
      }
    }

    lbDist[col] = lpSolution.col_value[col] - bestlb;
    if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
    ubDist[col] = bestub - lpSolution.col_value[col];
    if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;

    boundDist[col] = std::min(lbDist[col], ubDist[col]);
  }

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    double bestub = mipsolver.mipdata_->domain.col_upper_[col];
    double bestlb = mipsolver.mipdata_->domain.col_lower_[col];

    if (bestub - bestlb > 1.5 && bestub != kHighsInf && bestlb != -kHighsInf &&
        lprelaxation.getLpSolver().getBasis().col_status[col] !=
            HighsBasisStatus::kBasic &&
        lpSolution.col_value[col] > bestlb + mipsolver.mipdata_->feastol &&
        lpSolution.col_value[col] < bestub - mipsolver.mipdata_->feastol) {
      if (lprelaxation.getLpSolver().getBasis().col_status[col] ==
          HighsBasisStatus::kLower) {
        extendedCols.emplace_back(lprelaxation.colLower(col), col);
        extendedCols.back().dynamicColIndex = cutpool.acquireExtendedCol(
            extendedCols.back().splitCol, extendedCols.back().splitPos);
        boundDist[col] = 0.0;
        lbDist[col] = 0.0;
        ubDist[col] = bestub - lprelaxation.colLower(col);
        simpleLbDist[col] = lbDist[col];
        simpleUbDist[col] = ubDist[col];
        boundTypes[col] = BoundType::kLocalLb;
      } else {
        extendedCols.emplace_back(std::ceil(lprelaxation.colUpper(col) + 0.5),
                                  col);
        extendedCols.back().dynamicColIndex = cutpool.acquireExtendedCol(
            extendedCols.back().splitCol, extendedCols.back().splitPos);
        boundDist[col] = 0.0;
        lbDist[col] = lprelaxation.colUpper(col) - bestlb;
        ubDist[col] = 0.0;
        simpleLbDist[col] = lbDist[col];
        simpleUbDist[col] = ubDist[col];
        boundTypes[col] = BoundType::kLocalUb;
      }
      continue;
    }

    // todo: use binary variable bounds on integers?
    if (true || bestub - bestlb < 100.5) {
      if (bestlb == bestub) continue;
      lbDist[col] = lpSolution.col_value[col] - bestlb;
      if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
      simpleLbDist[col] = lbDist[col];
      ubDist[col] = bestub - lpSolution.col_value[col];
      if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;
      simpleUbDist[col] = ubDist[col];
      boundDist[col] = std::min(lbDist[col], ubDist[col]);
    } else {
      mipsolver.mipdata_->implications.cleanupVarbounds(col);
      if (mipsolver.mipdata_->domain.infeasible()) return;
      simpleUbDist[col] = bestub - lpSolution.col_value[col];
      if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
        simpleUbDist[col] = 0.0;

      double minbestub = bestub;
      size_t bestvubnodes = 0;

      simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
      if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
        simpleLbDist[col] = 0.0;
      double maxbestlb = bestlb;
      size_t bestvlbnodes = 0;

      for (const auto& vub : implications.getVUBs(col)) {
        if (vub.second.coef == kHighsInf) continue;
        if (mipsolver.mipdata_->domain.isFixed(vub.first)) continue;
        assert(mipsolver.mipdata_->domain.isBinary(vub.first));
        double vubval = lpSolution.col_value[vub.first] * vub.second.coef +
                        vub.second.constant;

        assert(vub.first >= 0 && vub.first < mipsolver.numCol());
        if (vubval <= lpSolution.col_value[col] + mipsolver.mipdata_->feastol) {
          size_t vubnodes =
              vub.second.coef > 0
                  ? mipsolver.mipdata_->nodequeue.numNodesDown(vub.first)
                  : mipsolver.mipdata_->nodequeue.numNodesUp(vub.first);
          double minvubval = vub.second.minValue();
          if (bestVub[col] == nullptr || vubnodes > bestvubnodes ||
              (vubnodes == bestvubnodes &&
               minvubval < minbestub - mipsolver.mipdata_->feastol)) {
            bestub = vubval;
            minbestub = minvubval;
            bestVub[col] = &vub;
            bestvubnodes = vubnodes;
          }
        }
      }

      for (const auto& vlb : implications.getVLBs(col)) {
        if (vlb.second.coef == -kHighsInf) continue;
        if (mipsolver.mipdata_->domain.isFixed(vlb.first)) continue;
        assert(mipsolver.mipdata_->domain.isBinary(vlb.first));
        assert(vlb.first >= 0 && vlb.first < mipsolver.numCol());
        double vlbval = lpSolution.col_value[vlb.first] * vlb.second.coef +
                        vlb.second.constant;

        if (vlbval >= lpSolution.col_value[col] - mipsolver.mipdata_->feastol) {
          size_t vlbnodes =
              vlb.second.coef > 0
                  ? mipsolver.mipdata_->nodequeue.numNodesUp(vlb.first)
                  : mipsolver.mipdata_->nodequeue.numNodesDown(vlb.first);
          double maxvlbval = vlb.second.maxValue();
          if (bestVlb[col] == nullptr || vlbnodes > bestvlbnodes ||
              (vlbnodes == bestvlbnodes &&
               maxvlbval > maxbestlb + mipsolver.mipdata_->feastol)) {
            bestlb = vlbval;
            maxbestlb = maxvlbval;
            bestVlb[col] = &vlb;
            bestvlbnodes = vlbnodes;
          }
        }
      }

      lbDist[col] = lpSolution.col_value[col] - bestlb;
      if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
      ubDist[col] = bestub - lpSolution.col_value[col];
      if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;

      boundDist[col] = std::min(lbDist[col], ubDist[col]);
    }
  }

  // setup information of slackVariables
  HighsInt numLpRow = lprelaxation.numRows();
  HighsInt indexOffset = mipsolver.numCol();
  for (HighsInt row = 0; row != numLpRow; ++row) {
    HighsInt slackIndex = indexOffset + row;
    double bestub = lprelaxation.slackUpper(row);
    double bestlb = lprelaxation.slackLower(row);

    if (bestlb == bestub) continue;

    lbDist[slackIndex] = lpSolution.row_value[row] - bestlb;
    if (lbDist[slackIndex] <= mipsolver.mipdata_->feastol)
      lbDist[slackIndex] = 0.0;
    simpleLbDist[slackIndex] = lbDist[slackIndex];

    ubDist[slackIndex] = bestub - lpSolution.row_value[row];
    if (ubDist[slackIndex] <= mipsolver.mipdata_->feastol)
      ubDist[slackIndex] = 0.0;
    simpleUbDist[slackIndex] = ubDist[slackIndex];

    boundDist[slackIndex] = std::min(lbDist[slackIndex], ubDist[slackIndex]);
  }

  vectorsum.setDimension(numTransformedCol, -cutpool.minExtendedColIndex());
}

HighsTransformedLp::~HighsTransformedLp() {
  for (const ExtendedCol& extCol : extendedCols)
    cutpool.releaseExtendedCol(extCol.dynamicColIndex);
}

bool HighsTransformedLp::transform(std::vector<double>& vals,
                                   std::vector<double>& upper,
                                   std::vector<double>& solval,
                                   std::vector<HighsInt>& inds, double& rhs,
                                   bool& integersPositive, bool preferVbds) {
  HighsCDouble tmpRhs = rhs;

  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = lprelaxation.numCols();

  HighsInt numNz = inds.size();
  bool removeZeros = false;

  for (HighsInt i = 0; i != numNz; ++i) {
    HighsInt col = inds[i];

    double lb;
    double ub;

    if (col < 0) {
      auto split = mip.mipdata_->cutpool.getExtendedColSplit(col);
      if (mip.mipdata_->domain.col_lower_[split.first] >= split.second) {
        // extended column is globally fixed to 1.0
        tmpRhs -= vals[i];
        vals[i] = 0.0;
        removeZeros = true;
        continue;
      } else if (mip.mipdata_->domain.col_upper_[split.first] <
                 split.second - 0.5) {
        // extended column is globally fixed to 0.0
        vals[i] = 0.0;
        removeZeros = true;
        continue;
      } else {
        // extended column is globally unfixed
        integersPositive = false;
      }
      continue;
    } else if (col < slackOffset) {
      lb = mip.mipdata_->domain.col_lower_[col];
      ub = mip.mipdata_->domain.col_upper_[col];
    } else {
      HighsInt row = col - slackOffset;
      lb = lprelaxation.slackLower(row);
      ub = lprelaxation.slackUpper(row);
    }

    if (ub - lb < mip.options_mip_->small_matrix_value) {
      tmpRhs -= std::min(lb, ub) * vals[i];
      vals[i] = 0.0;
      removeZeros = true;
      continue;
    }

    if (lb == -kHighsInf && ub == kHighsInf) return false;

    if (lprelaxation.isColIntegral(col)) {
      if (lb == -kHighsInf || ub == kHighsInf) integersPositive = false;
      if (boundTypes[col] != BoundType::kLocalLb &&
          boundTypes[col] != BoundType::kLocalUb) {
        bool useVbd = false;
        if (ub - lb > 1.5) {
          if (vals[i] < 0 && ubDist[col] == 0.0 &&
              simpleUbDist[col] > mip.mipdata_->feastol) {
            boundTypes[col] = BoundType::kVariableUb;
            useVbd = true;
          } else if (vals[i] > 0.0 && lbDist[col] == 0.0 &&
                     simpleLbDist[col] > mip.mipdata_->feastol) {
            boundTypes[col] = BoundType::kVariableLb;
            useVbd = true;
          }
        }

        if (!useVbd) continue;
      }
    } else {
      if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        if (!bestVlb[col])
          boundTypes[col] = BoundType::kSimpleLb;
        else if (preferVbds || vals[i] > 0 ||
                 simpleLbDist[col] > lbDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        if (!bestVub[col])
          boundTypes[col] = BoundType::kSimpleUb;
        else if (preferVbds || vals[i] < 0 ||
                 simpleUbDist[col] > ubDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      } else if (vals[i] > 0) {
        if (bestVlb[col])
          boundTypes[col] = BoundType::kVariableLb;
        else if (preferVbds && bestVub[col])
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else {
        if (bestVub[col])
          boundTypes[col] = BoundType::kVariableUb;
        else if (preferVbds && bestVlb[col])
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      }
    }

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb:
        if (vals[i] > 0) {
          tmpRhs -= lb * vals[i];
          vals[i] = 0.0;
          removeZeros = true;
        }
        break;
      case BoundType::kSimpleUb:
        if (vals[i] < 0) {
          tmpRhs -= ub * vals[i];
          vals[i] = 0.0;
          removeZeros = true;
        }
        break;
      case BoundType::kVariableLb:
        tmpRhs -= bestVlb[col]->second.constant * vals[i];
        vectorsum.add(bestVlb[col]->first, vals[i] * bestVlb[col]->second.coef);
        if (vals[i] > 0) vals[i] = 0;
        break;
      case BoundType::kVariableUb:
        tmpRhs -= bestVub[col]->second.constant * vals[i];
        vectorsum.add(bestVub[col]->first, vals[i] * bestVub[col]->second.coef);
        vals[i] = -vals[i];
        if (vals[i] > 0) vals[i] = 0;
        break;
      case BoundType::kLocalLb: {
        integersPositive = false;
        double localLb = lprelaxation.colLower(col);
        HighsInt extCol = cutpool.getExtendedColIndex(col, localLb);

        double constant = lb;
        double coef = (localLb - lb);
        assert(coef > 0);
        tmpRhs -= constant * vals[i];
        vectorsum.add(extCol, vals[i] * coef);
        break;
      }
      case BoundType::kLocalUb: {
        integersPositive = false;
        double localUb = lprelaxation.colUpper(col);
        HighsInt extCol =
            cutpool.getExtendedColIndex(col, std::ceil(localUb + 0.5));

        double constant = localUb;
        double coef = (ub - localUb);
        assert(coef > 0);
        tmpRhs -= constant * vals[i];
        vectorsum.add(extCol, vals[i] * coef);
        vals[i] = -vals[i];
      }
    }
  }

  if (!vectorsum.getNonzeros().empty()) {
    for (HighsInt i = 0; i != numNz; ++i) {
      if (vals[i] != 0.0) vectorsum.add(inds[i], vals[i]);
    }

    double maxError = 0.0;
    auto IsZero = [&](HighsInt col, double val) {
      double absval = std::abs(val);
      if (absval <= mip.options_mip_->small_matrix_value) return true;

      return false;
    };

    vectorsum.cleanup(IsZero);
    if (maxError > mip.mipdata_->feastol) return false;

    inds = vectorsum.getNonzeros();
    numNz = inds.size();

    vals.resize(numNz);
    for (HighsInt j = 0; j != numNz; ++j) vals[j] = vectorsum.getValue(inds[j]);

    vectorsum.clear();
  } else if (removeZeros) {
    for (HighsInt i = numNz - 1; i >= 0; --i) {
      if (vals[i] == 0) {
        --numNz;
        vals[i] = vals[numNz];
        inds[i] = inds[numNz];
        std::swap(vals[i], vals[numNz]);
        std::swap(inds[i], inds[numNz]);
      }
    }

    vals.resize(numNz);
    inds.resize(numNz);
  }

  if (integersPositive) {
    // complement integers to make coefficients positive
    for (HighsInt j = 0; j != numNz; ++j) {
      HighsInt col = inds[j];
      if (!lprelaxation.isColIntegral(inds[j])) continue;

      if (vals[j] > 0)
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    }
  } else {
    // complement integers with closest bound
    for (HighsInt j = 0; j != numNz; ++j) {
      HighsInt col = inds[j];

      if (col < 0) continue;

      if (!lprelaxation.isColIntegral(inds[j])) continue;

      // do not touch integers that use a local bound as they are
      // already complemented with their closest bound and we need to keep
      // their bound type as local for correct retransformation
      if (boundTypes[col] == BoundType::kLocalLb ||
          boundTypes[col] == BoundType::kLocalUb)
        continue;

      if (lbDist[col] < ubDist[col])
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    }
  }

  upper.resize(numNz);
  solval.resize(numNz);

  for (HighsInt j = 0; j != numNz; ++j) {
    HighsInt col = inds[j];

    double lb;
    double ub;

    if (col < 0) {
      auto split = mip.mipdata_->cutpool.getExtendedColSplit(col);
      upper[j] = 1.0;

      if (lprelaxation.colLower(split.first) >= split.second) {
        // extended column is locally fixed to 1.0
        solval[j] = 1.0;
      } else if (lprelaxation.colUpper(split.first) < split.second - 0.5) {
        // extended column is locally fixed to 0.0
        solval[j] = 0.0;
      } else {
        // extended column is locally unfixed, set solution to most feasible
        // value
        solval[j] = vals[j] < 0;
      }

      if (solval[j] == 1.0) {
        extColIsComplemented.insert(col);
        tmpRhs -= vals[j];
        vals[j] = -vals[j];
        solval[j] = 0.0;
      } else {
        extColIsComplemented.erase(col);
      }

      continue;
    }

    if (col < slackOffset) {
      lb = mip.mipdata_->domain.col_lower_[col];
      ub = mip.mipdata_->domain.col_upper_[col];
    } else {
      HighsInt row = col - slackOffset;
      lb = lprelaxation.slackLower(row);
      ub = lprelaxation.slackUpper(row);
    }

    upper[j] = ub - lb;

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb: {
        assert(lb != -kHighsInf);
        tmpRhs -= lb * vals[j];
        solval[j] = lbDist[col];
        break;
      }
      case BoundType::kSimpleUb: {
        assert(ub != kHighsInf);
        tmpRhs -= ub * vals[j];
        vals[j] = -vals[j];
        solval[j] = ubDist[col];
        break;
      }
      case BoundType::kVariableLb:
        solval[j] = lbDist[col];
        break;
      case BoundType::kVariableUb:
        solval[j] = ubDist[col];
        break;
      case BoundType::kLocalLb:
        upper[j] = ub - lprelaxation.colLower(col);
        solval[j] = 0.0;
        break;
      case BoundType::kLocalUb:
        upper[j] = lprelaxation.colUpper(col) - lb;
        solval[j] = 0.0;
    }
  }

  rhs = double(tmpRhs);

  if (numNz == 0 && rhs >= -mip.mipdata_->feastol) return false;

  return true;
}

bool HighsTransformedLp::untransform(std::vector<double>& vals,
                                     std::vector<HighsInt>& inds, double& rhs,
                                     bool integral) {
  HighsCDouble tmpRhs = rhs;
  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = mip.numCol();

  HighsInt numNz = inds.size();

  for (HighsInt i = 0; i != numNz; ++i) {
    if (vals[i] == 0.0) continue;
    HighsInt col = inds[i];
    if (col < 0) {
      if (extColIsComplemented.find(col)) {
        tmpRhs -= vals[i];
        vectorsum.add(col, -vals[i]);
      } else
        vectorsum.add(col, vals[i]);

      continue;
    }

    switch (boundTypes[col]) {
      case BoundType::kVariableLb: {
        tmpRhs += bestVlb[col]->second.constant * vals[i];
        vectorsum.add(bestVlb[col]->first,
                      -vals[i] * bestVlb[col]->second.coef);
        vectorsum.add(col, vals[i]);
        break;
      }
      case BoundType::kVariableUb: {
        tmpRhs -= bestVub[col]->second.constant * vals[i];
        vectorsum.add(bestVub[col]->first, vals[i] * bestVub[col]->second.coef);
        vectorsum.add(col, -vals[i]);
        break;
      }
      case BoundType::kSimpleLb: {
        if (col < slackOffset) {
          tmpRhs += vals[i] * mip.mipdata_->domain.col_lower_[col];
          vectorsum.add(col, vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs += vals[i] * lprelaxation.slackLower(row);

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
        break;
      }
      case BoundType::kSimpleUb: {
        if (col < slackOffset) {
          tmpRhs -= vals[i] * mip.mipdata_->domain.col_upper_[col];
          vectorsum.add(col, -vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs -= vals[i] * lprelaxation.slackUpper(row);
          vals[i] = -vals[i];

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
        break;
      }
      case BoundType::kLocalLb: {
        // printf("got nonzero coefficient of integer sitting at local
        // bound\n");
        double localLb = lprelaxation.colLower(col);
        HighsInt extCol = cutpool.getExtendedColIndex(col, localLb);

        double constant = mip.mipdata_->domain.col_lower_[col];
        double coef = (localLb - mip.mipdata_->domain.col_lower_[col]);

        tmpRhs += constant * vals[i];
        vectorsum.add(extCol, -vals[i] * coef);
        vectorsum.add(col, vals[i]);
        break;
      }
      case BoundType::kLocalUb: {
        // printf("got nonzero coefficient of integer sitting at local
        // bound\n");
        double localUb = lprelaxation.colUpper(col);
        HighsInt extCol =
            cutpool.getExtendedColIndex(col, std::ceil(localUb + 0.5));

        double constant = localUb;
        double coef = mip.mipdata_->domain.col_upper_[col] - localUb;

        tmpRhs -= constant * vals[i];
        vectorsum.add(extCol, vals[i] * coef);
        vectorsum.add(col, -vals[i]);
        break;
      }
    }
  }

  // if (!localBoundStrengthenings.empty()) {
  //   printf("have cut with %ld local bound indicators\n",
  //          localBoundStrengthenings.size());
  // }

  if (integral) {
    // if the cut is integral, we just round all coefficient values and the
    // right hand side to the nearest integral value, as small deviation
    // only come from numerical errors during resubstitution of slack variables

    auto IsZero = [&](HighsInt col, double val) {
      assert(col < mip.numCol());
      return std::round(val) == 0.0;
    };

    vectorsum.cleanup(IsZero);
    rhs = std::round(double(tmpRhs));
  } else {
    bool abort = false;
    auto IsZero = [&](HighsInt col, double val) {
      double absval = std::abs(val);
      if (absval <= mip.options_mip_->small_matrix_value) return true;

      if (absval <= mip.mipdata_->feastol) {
        if (val > 0) {
          double lb = col >= 0 ? mip.mipdata_->domain.col_lower_[col] : 0;
          if (mip.mipdata_->domain.col_lower_[col] == -kHighsInf)
            abort = true;
          else
            tmpRhs -= val * mip.mipdata_->domain.col_lower_[col];
        } else {
          double ub = col >= 0 ? mip.mipdata_->domain.col_upper_[col] : 1;
          if (mip.mipdata_->domain.col_upper_[col] == kHighsInf)
            abort = true;
          else
            tmpRhs -= val * mip.mipdata_->domain.col_upper_[col];
        }
        return true;
      }
      return false;
    };

    vectorsum.cleanup(IsZero);
    if (abort) {
      vectorsum.clear();
      return false;
    }
    rhs = double(tmpRhs);
  }

  inds = vectorsum.getNonzeros();
  numNz = inds.size();
  vals.resize(numNz);

  if (integral)
    for (HighsInt i = 0; i != numNz; ++i)
      vals[i] = std::round(vectorsum.getValue(inds[i]));
  else
    for (HighsInt i = 0; i != numNz; ++i) vals[i] = vectorsum.getValue(inds[i]);

  vectorsum.clear();

  return true;
}
