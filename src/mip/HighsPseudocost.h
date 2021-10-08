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

#ifndef HIGHS_PSEUDOCOST_H_
#define HIGHS_PSEUDOCOST_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include "util/HighsInt.h"

class HighsMipSolver;
namespace presolve {
class HighsPostsolveStack;
}

class HighsPseudocost;

struct HighsPseudocostInitialization {
  std::vector<HighsFloat> pseudocostup;
  std::vector<HighsFloat> pseudocostdown;
  std::vector<HighsInt> nsamplesup;
  std::vector<HighsInt> nsamplesdown;
  std::vector<HighsFloat> inferencesup;
  std::vector<HighsFloat> inferencesdown;
  std::vector<HighsInt> ninferencesup;
  std::vector<HighsInt> ninferencesdown;
  std::vector<HighsFloat> conflictscoreup;
  std::vector<HighsFloat> conflictscoredown;
  HighsFloat cost_total;
  HighsFloat inferences_total;
  HighsFloat conflict_avg_score;
  int64_t nsamplestotal;
  int64_t ninferencestotal;

  HighsPseudocostInitialization(const HighsPseudocost& pscost,
                                HighsInt maxCount);
  HighsPseudocostInitialization(
      const HighsPseudocost& pscost, HighsInt maxCount,
      const presolve::HighsPostsolveStack& postsolveStack);
};
class HighsPseudocost {
  friend struct HighsPseudocostInitialization;
  std::vector<HighsFloat> pseudocostup;
  std::vector<HighsFloat> pseudocostdown;
  std::vector<HighsInt> nsamplesup;
  std::vector<HighsInt> nsamplesdown;
  std::vector<HighsFloat> inferencesup;
  std::vector<HighsFloat> inferencesdown;
  std::vector<HighsInt> ninferencesup;
  std::vector<HighsInt> ninferencesdown;
  std::vector<HighsInt> ncutoffsup;
  std::vector<HighsInt> ncutoffsdown;
  std::vector<HighsFloat> conflictscoreup;
  std::vector<HighsFloat> conflictscoredown;

  HighsFloat conflict_weight;
  HighsFloat conflict_avg_score;
  HighsFloat cost_total;
  HighsFloat inferences_total;
  int64_t nsamplestotal;
  int64_t ninferencestotal;
  int64_t ncutoffstotal;
  HighsInt minreliable;
  HighsFloat degeneracyFactor;

 public:
  HighsPseudocost() = default;
  HighsPseudocost(const HighsMipSolver& mipsolver);

  void subtractBase(const HighsPseudocost& base) {
    HighsInt ncols = pseudocostup.size();

    for (HighsInt i = 0; i != ncols; ++i) {
      pseudocostup[i] -= base.pseudocostup[i];
      pseudocostdown[i] -= base.pseudocostdown[i];
      nsamplesup[i] -= base.nsamplesup[i];
      nsamplesdown[i] -= base.nsamplesdown[i];
    }
  }

  void increaseConflictWeight() {
    conflict_weight *= 1.02;

    if (conflict_weight > 1000.0) {
      HighsFloat scale = 1.0 / conflict_weight;
      conflict_weight = 1.0;
      conflict_avg_score *= scale;

      HighsInt numCol = conflictscoreup.size();
      for (HighsInt i = 0; i < numCol; ++i) {
        conflictscoreup[i] *= scale;
        conflictscoredown[i] *= scale;
      }
    }
  }

  void setDegeneracyFactor(HighsFloat degeneracyFactor) {
    assert(degeneracyFactor >= 1.0);
    this->degeneracyFactor = degeneracyFactor;
  }

  void increaseConflictScoreUp(HighsInt col) {
    conflictscoreup[col] += conflict_weight;
    conflict_avg_score += conflict_weight;
  }

  void increaseConflictScoreDown(HighsInt col) {
    conflictscoredown[col] += conflict_weight;
    conflict_avg_score += conflict_weight;
  }

  void setMinReliable(HighsInt minreliable) { this->minreliable = minreliable; }

  HighsInt getMinReliable() const { return minreliable; }

  HighsInt getNumObservations(HighsInt col) const {
    return nsamplesup[col] + nsamplesdown[col];
  }

  HighsInt getNumObservationsUp(HighsInt col) const { return nsamplesup[col]; }

  HighsInt getNumObservationsDown(HighsInt col) const {
    return nsamplesdown[col];
  }

  void addCutoffObservation(HighsInt col, bool upbranch) {
    ++ncutoffstotal;
    if (upbranch)
      ncutoffsup[col] += 1;
    else
      ncutoffsdown[col] += 1;
  }

  void addObservation(HighsInt col, HighsFloat delta, HighsFloat objdelta) {
    assert(delta != 0.0);
    assert(objdelta >= 0.0);
    if (delta > 0.0) {
      HighsFloat unit_gain = objdelta / delta;
      HighsFloat d = unit_gain - pseudocostup[col];
      nsamplesup[col] += 1;
      pseudocostup[col] += d / nsamplesup[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / nsamplestotal;
    } else {
      HighsFloat unit_gain = -objdelta / delta;
      HighsFloat d = unit_gain - pseudocostdown[col];
      nsamplesdown[col] += 1;
      pseudocostdown[col] += d / nsamplesdown[col];

      d = unit_gain - cost_total;
      ++nsamplestotal;
      cost_total += d / nsamplestotal;
    }
  }

  void addInferenceObservation(HighsInt col, HighsInt ninferences,
                               bool upbranch) {
    HighsFloat d = ninferences - inferences_total;
    ++ninferencestotal;
    inferences_total += d / ninferencestotal;
    if (upbranch) {
      d = ninferences - inferencesup[col];
      ninferencesup[col] += 1;
      inferencesup[col] += d / ninferencesup[col];
    } else {
      d = ninferences - inferencesdown[col];
      ninferencesdown[col] += 1;
      inferencesdown[col] += d / ninferencesdown[col];
    }
  }

  bool isReliable(HighsInt col) const {
    return std::min(nsamplesup[col], nsamplesdown[col]) >= minreliable;
  }

  bool isReliableUp(HighsInt col) const {
    return nsamplesup[col] >= minreliable;
  }

  bool isReliableDown(HighsInt col) const {
    return nsamplesdown[col] >= minreliable;
  }

  HighsFloat getAvgPseudocost() const { return cost_total; }

  HighsFloat getPseudocostUp(HighsInt col, HighsFloat frac, HighsFloat offset) const {
    HighsFloat up = std::ceil(frac) - frac;
    HighsFloat cost;

    if (nsamplesup[col] == 0 || nsamplesup[col] < minreliable) {
      HighsFloat weightPs = nsamplesup[col] == 0
                            ? 0
                            : 0.9 + 0.1 * nsamplesup[col] / (HighsFloat)minreliable;
      cost = weightPs * pseudocostup[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostup[col];
    return up * (offset + cost);
  }

  HighsFloat getPseudocostDown(HighsInt col, HighsFloat frac, HighsFloat offset) const {
    HighsFloat down = frac - std::floor(frac);
    HighsFloat cost;

    if (nsamplesdown[col] == 0 || nsamplesdown[col] < minreliable) {
      HighsFloat weightPs = nsamplesdown[col] == 0 ? 0
                                               : 0.9 + 0.1 * nsamplesdown[col] /
                                                           (HighsFloat)minreliable;
      cost = weightPs * pseudocostdown[col];
      cost += (1.0 - weightPs) * getAvgPseudocost();
    } else
      cost = pseudocostdown[col];

    return down * (offset + cost);
  }

  HighsFloat getPseudocostUp(HighsInt col, HighsFloat frac) const {
    HighsFloat up = std::ceil(frac) - frac;
    if (nsamplesup[col] == 0) return up * cost_total;
    return up * pseudocostup[col];
  }

  HighsFloat getPseudocostDown(HighsInt col, HighsFloat frac) const {
    HighsFloat down = frac - std::floor(frac);
    if (nsamplesdown[col] == 0) return down * cost_total;
    return down * pseudocostdown[col];
  }

  HighsFloat getConflictScoreUp(HighsInt col) const {
    return conflictscoreup[col] / conflict_weight;
  }

  HighsFloat getConflictScoreDown(HighsInt col) const {
    return conflictscoredown[col] / conflict_weight;
  }

  HighsFloat getScore(HighsInt col, HighsFloat upcost, HighsFloat downcost) const {
    HighsFloat costScore = std::max(upcost, 1e-6) * std::max(downcost, 1e-6) /
                       std::max(1e-6, cost_total * cost_total);
    HighsFloat inferenceScore = std::max(inferencesup[col], 1e-6) *
                            std::max(inferencesdown[col], 1e-6) /
                            std::max(1e-6, inferences_total * inferences_total);

    HighsFloat cutOffScoreUp =
        ncutoffsup[col] /
        std::max(1.0, HighsFloat(ncutoffsup[col] + nsamplesup[col]));
    HighsFloat cutOffScoreDown =
        ncutoffsdown[col] /
        std::max(1.0, HighsFloat(ncutoffsdown[col] + nsamplesdown[col]));
    HighsFloat avgCutoffs =
        ncutoffstotal / std::max(1.0, HighsFloat(ncutoffstotal + nsamplestotal));

    HighsFloat cutoffScore = std::max(cutOffScoreUp, 1e-6) *
                         std::max(cutOffScoreDown, 1e-6) /
                         std::max(1e-6, avgCutoffs * avgCutoffs);

    HighsFloat conflictScoreUp = conflictscoreup[col] / conflict_weight;
    HighsFloat conflictScoreDown = conflictscoredown[col] / conflict_weight;
    HighsFloat conflictScoreAvg =
        conflict_avg_score / (conflict_weight * conflictscoreup.size());
    HighsFloat conflictScore = std::max(conflictScoreUp, 1e-6) *
                           std::max(conflictScoreDown, 1e-6) /
                           std::max(1e-6, conflictScoreAvg * conflictScoreAvg);

    auto mapScore = [](HighsFloat score) { return 1.0 - 1.0 / (1.0 + score); };
    return mapScore(costScore) / degeneracyFactor +
           degeneracyFactor *
               (1e-2 * mapScore(conflictScore) +
                1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  HighsFloat getScore(HighsInt col, HighsFloat frac) const {
    HighsFloat upcost = getPseudocostUp(col, frac);
    HighsFloat downcost = getPseudocostDown(col, frac);

    return getScore(col, upcost, downcost);
  }

  HighsFloat getScoreUp(HighsInt col, HighsFloat frac) const {
    HighsFloat costScore = getPseudocostUp(col, frac) / std::max(1e-6, cost_total);
    HighsFloat inferenceScore =
        inferencesup[col] / std::max(1e-6, inferences_total);

    HighsFloat cutOffScoreUp =
        ncutoffsup[col] /
        std::max(1.0, HighsFloat(ncutoffsup[col] + nsamplesup[col]));
    HighsFloat avgCutoffs =
        ncutoffstotal / std::max(1.0, HighsFloat(ncutoffstotal + nsamplestotal));

    HighsFloat cutoffScore = cutOffScoreUp / std::max(1e-6, avgCutoffs);

    HighsFloat conflictScoreUp = conflictscoreup[col] / conflict_weight;
    HighsFloat conflictScoreAvg =
        conflict_avg_score / (conflict_weight * conflictscoreup.size());
    HighsFloat conflictScore = conflictScoreUp / std::max(1e-6, conflictScoreAvg);

    auto mapScore = [](HighsFloat score) { return 1.0 - 1.0 / (1.0 + score); };

    return mapScore(costScore) +
           (1e-2 * mapScore(conflictScore) +
            1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  HighsFloat getScoreDown(HighsInt col, HighsFloat frac) const {
    HighsFloat costScore =
        getPseudocostDown(col, frac) / std::max(1e-6, cost_total);
    HighsFloat inferenceScore =
        inferencesdown[col] / std::max(1e-6, inferences_total);

    HighsFloat cutOffScoreDown =
        ncutoffsdown[col] /
        std::max(1.0, HighsFloat(ncutoffsdown[col] + nsamplesdown[col]));
    HighsFloat avgCutoffs =
        ncutoffstotal / std::max(1.0, HighsFloat(ncutoffstotal + nsamplestotal));

    HighsFloat cutoffScore = cutOffScoreDown / std::max(1e-6, avgCutoffs);

    HighsFloat conflictScoreDown = conflictscoredown[col] / conflict_weight;
    HighsFloat conflictScoreAvg =
        conflict_avg_score / (conflict_weight * conflictscoredown.size());
    HighsFloat conflictScore = conflictScoreDown / std::max(1e-6, conflictScoreAvg);

    auto mapScore = [](HighsFloat score) { return 1.0 - 1.0 / (1.0 + score); };

    return mapScore(costScore) +
           (1e-2 * mapScore(conflictScore) +
            1e-4 * (mapScore(cutoffScore) + mapScore(inferenceScore)));
  }

  HighsFloat getAvgInferencesUp(HighsInt col) const { return inferencesup[col]; }

  HighsFloat getAvgInferencesDown(HighsInt col) const {
    return inferencesdown[col];
  }
};

#endif
