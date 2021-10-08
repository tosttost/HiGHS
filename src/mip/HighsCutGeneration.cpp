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
#include "mip/HighsCutGeneration.h"

#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "pdqsort/pdqsort.h"
#include "util/HighsIntegers.h"

HighsCutGeneration::HighsCutGeneration(const HighsLpRelaxation& lpRelaxation,
                                       HighsCutPool& cutpool)
    : lpRelaxation(lpRelaxation),
      cutpool(cutpool),
      randgen(lpRelaxation.getMipSolver().options_mip_->random_seed +
              lpRelaxation.getNumLpIterations() + cutpool.getNumCuts()),
      feastol(lpRelaxation.getMipSolver().mipdata_->feastol),
      epsilon(lpRelaxation.getMipSolver().mipdata_->epsilon) {}

bool HighsCutGeneration::determineCover(bool lpSol) {
  if (rhs <= 10 * feastol) return false;

  cover.clear();
  cover.reserve(rowlen);

  for (HighsInt j = 0; j != rowlen; ++j) {
    if (!isintegral[j]) continue;

    if (lpSol && solval[j] <= feastol) continue;

    cover.push_back(j);
  }

  HighsInt maxCoverSize = cover.size();
  HighsInt coversize = 0;
  HighsInt r = randgen.integer();
  coverweight = 0.0;
  if (lpSol) {
    // take all variables that sit at their upper bound always into the cover
    coversize = std::partition(cover.begin(), cover.end(),
                               [&](HighsInt j) {
                                 return solval[j] >= upper[j] - feastol;
                               }) -
                cover.begin();

    for (HighsInt i = 0; i != coversize; ++i) {
      HighsInt j = cover[i];

      assert(solval[j] >= upper[j] - feastol);

      coverweight += vals[j] * upper[j];
    }

    const auto& nodequeue = lpRelaxation.getMipSolver().mipdata_->nodequeue;
    // sort the remaining variables by the contribution to the rows activity in
    // the current solution
    pdqsort(cover.begin() + coversize, cover.begin() + maxCoverSize,
            [&](HighsInt i, HighsInt j) {
              if (upper[i] < 1.5 && upper[j] > 1.5) return true;
              if (upper[i] > 1.5 && upper[j] < 1.5) return false;

              HighsFloat contributionA = solval[i] * vals[i];
              HighsFloat contributionB = solval[j] * vals[j];

              if (contributionA > contributionB + feastol) return true;
              if (contributionA < contributionB - feastol) return false;
              // for equal contributions take the larger coefficients first
              // because this makes some of the lifting functions more likely
              // to generate a facet
              // if the value is equal too, choose a random tiebreaker based
              // on hashing the column index and the current number of pool
              // cuts
              if (std::abs(vals[i] - vals[j]) <= feastol)
                return HighsHashHelpers::hash(std::make_pair(inds[i], r)) >
                       HighsHashHelpers::hash(std::make_pair(inds[j], r));
              return vals[i] > vals[j];
            });
  } else {
    // the current solution
    const auto& nodequeue = lpRelaxation.getMipSolver().mipdata_->nodequeue;

    pdqsort(cover.begin() + coversize, cover.begin() + maxCoverSize,
            [&](HighsInt i, HighsInt j) {
              if (solval[i] > feastol && solval[j] <= feastol) return true;
              if (solval[i] <= feastol && solval[j] > feastol) return false;

              int64_t numNodesA;
              int64_t numNodesB;

              numNodesA = complementation[i] ? nodequeue.numNodesDown(inds[i])
                                             : nodequeue.numNodesUp(inds[i]);

              numNodesB = complementation[j] ? nodequeue.numNodesDown(inds[j])
                                             : nodequeue.numNodesUp(inds[j]);

              if (numNodesA > numNodesB) return true;
              if (numNodesA < numNodesB) return false;

              return HighsHashHelpers::hash(std::make_pair(inds[i], r)) >
                     HighsHashHelpers::hash(std::make_pair(inds[j], r));
            });
  }

  const HighsFloat minlambda =
      std::max(10 * feastol, feastol * std::abs(HighsFloat(rhs)));

  for (; coversize != maxCoverSize; ++coversize) {
    HighsFloat lambda = HighsFloat(coverweight - rhs);
    if (lambda > minlambda) break;

    HighsInt j = cover[coversize];
    coverweight += vals[j] * upper[j];
  }
  if (coversize == 0) return false;

  coverweight.renormalize();
  lambda = coverweight - rhs;

  if (lambda <= minlambda) return false;

  cover.resize(coversize);
  assert(lambda > feastol);

  return true;
}

void HighsCutGeneration::separateLiftedKnapsackCover() {
  const HighsFloat feastol = lpRelaxation.getMipSolver().mipdata_->feastol;

  const HighsInt coversize = cover.size();

  std::vector<HighsFloat> S;
  S.resize(coversize);
  std::vector<int8_t> coverflag;
  coverflag.resize(rowlen);
  pdqsort_branchless(cover.begin(), cover.end(),
                     [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; });

  HighsCD0uble abartmp = vals[cover[0]];
  HighsCD0uble sigma = lambda;
  for (HighsInt i = 1; i != coversize; ++i) {
    HighsCD0uble delta = abartmp - vals[cover[i]];
    HighsCD0uble kdelta = HighsFloat(i) * delta;
    if (HighsFloat(kdelta) < HighsFloat(sigma)) {
      abartmp = vals[cover[i]];
      sigma -= kdelta;
    } else {
      abartmp -= sigma * (1.0 / i);
      sigma = 0.0;
      break;
    }
  }

  if (HighsFloat(sigma) > 0) abartmp = HighsCD0uble(rhs) / HighsFloat(coversize);

  HighsFloat abar = HighsFloat(abartmp);

  HighsCD0uble sum = 0.0;
  HighsInt cplussize = 0;
  for (HighsInt i = 0; i != coversize; ++i) {
    sum += std::min(abar, vals[cover[i]]);
    S[i] = HighsFloat(sum);

    if (vals[cover[i]] > abar + feastol) {
      ++cplussize;
      coverflag[cover[i]] = 1;
    } else
      coverflag[cover[i]] = -1;
  }
  assert(std::abs(HighsFloat(sum - rhs) / HighsFloat(rhs)) <= 1e-14);
  bool halfintegral = false;

  /* define the lifting function */
  auto g = [&](HighsFloat z) {
    HighsFloat hfrac = z / abar;
    HighsFloat coef = 0.0;

    HighsInt h = std::floor(hfrac + 0.5);
    if (h != 0 && std::abs(hfrac - h) * std::max(1.0, abar) <= epsilon &&
        h <= cplussize - 1) {
      halfintegral = true;
      coef = 0.5;
    }

    h = std::max(h - 1, HighsInt{0});
    for (; h < coversize; ++h) {
      if (z <= S[h] + feastol) break;
    }

    return coef + h;
  };

  rhs = coversize - 1;

  for (HighsInt i = 0; i != rowlen; ++i) {
    if (vals[i] == 0.0) continue;
    if (coverflag[i] == -1) {
      vals[i] = 1;
    } else {
      vals[i] = g(vals[i]);
    }
  }

  if (halfintegral) {
    rhs *= 2;
    for (HighsInt i = 0; i != rowlen; ++i) vals[i] *= 2;
  }

  // resulting cut is always integral
  integralSupport = true;
  integralCoefficients = true;
}

bool HighsCutGeneration::separateLiftedMixedBinaryCover() {
  HighsInt coversize = cover.size();
  std::vector<HighsFloat> S;
  S.resize(coversize);
  std::vector<uint8_t> coverflag;
  coverflag.resize(rowlen);

  if (coversize == 0) return false;

  for (HighsInt i = 0; i != coversize; ++i) coverflag[cover[i]] = 1;

  pdqsort_branchless(cover.begin(), cover.end(),
                     [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; });
  HighsCD0uble sum = 0;

  HighsInt p = coversize;
  for (HighsInt i = 0; i != coversize; ++i) {
    if (vals[cover[i]] - lambda <= epsilon) {
      p = i;
      break;
    }
    sum += vals[cover[i]];
    S[i] = HighsFloat(sum);
  }
  if (p == 0) return false;
  /* define the lifting function */
  auto phi = [&](HighsFloat a) {
    for (HighsInt i = 0; i < p; ++i) {
      if (a <= S[i] - lambda) return HighsFloat(i * lambda);

      if (a <= S[i]) return HighsFloat((i + 1) * lambda + (HighsCD0uble(a) - S[i]));
    }

    return HighsFloat(p * lambda + (HighsCD0uble(a) - S[p - 1]));
  };

  rhs = -lambda;

  integralCoefficients = false;
  integralSupport = true;
  for (HighsInt i = 0; i != rowlen; ++i) {
    if (!isintegral[i]) {
      if (vals[i] < 0)
        integralSupport = false;
      else
        vals[i] = 0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = std::min(vals[i], HighsFloat(lambda));
      rhs += vals[i];
    } else {
      vals[i] = phi(vals[i]);
    }
  }

  return true;
}

bool HighsCutGeneration::separateLiftedMixedIntegerCover() {
  HighsInt coversize = cover.size();

  HighsInt l = -1;

  std::vector<uint8_t> coverflag;
  coverflag.resize(rowlen);
  for (HighsInt i : cover) coverflag[i] = 1;

  auto comp = [&](HighsInt a, HighsInt b) { return vals[a] > vals[b]; };
  pdqsort_branchless(cover.begin(), cover.end(), comp);

  std::vector<HighsCD0uble> a;
  std::vector<HighsCD0uble> u;
  std::vector<HighsCD0uble> m;

  a.resize(coversize);
  u.resize(coversize + 1);
  m.resize(coversize + 1);

  HighsCD0uble usum = 0.0;
  HighsCD0uble msum = 0.0;
  // set up the partial sums of the upper bounds, and the contributions
  for (HighsInt c = 0; c != coversize; ++c) {
    HighsInt i = cover[c];

    u[c] = usum;
    m[c] = msum;
    a[c] = vals[i];
    HighsFloat ub = upper[i];
    usum += ub;
    msum += ub * a[c];
  }

  u[coversize] = usum;
  m[coversize] = msum;

  // determine which variable in the cover we want to create the MIR inequality
  // from which we lift we try to select a variable to have the highest chance
  // of satisfying the facet conditions for the superadditive lifting function
  // gamma to be satisfied.
  HighsInt lpos = -1;
  HighsInt bestlCplusend = -1;
  HighsFloat bestlVal = 0.0;
  bool bestlAtUpper = true;

  for (HighsInt i = 0; i != coversize; ++i) {
    HighsInt j = cover[i];
    HighsFloat ub = upper[j];

    bool atUpper = solval[j] >= ub - feastol;
    if (atUpper && !bestlAtUpper) continue;

    HighsFloat mju = ub * vals[j];
    HighsCD0uble mu = mju - lambda;

    if (mu <= 10 * feastol) continue;
    if (std::abs(vals[j]) < 1000 * feastol) continue;

    HighsFloat mudival = HighsFloat(mu / vals[j]);
    if (std::abs(std::round(mudival) - mudival) <= feastol) continue;
    HighsFloat eta = ceil(mudival);

    HighsCD0uble ulminusetaplusone = HighsCD0uble(ub) - eta + 1.0;
    HighsCD0uble cplusthreshold = ulminusetaplusone * vals[j];

    HighsInt cplusend =
        std::upper_bound(cover.begin(), cover.end(), HighsFloat(cplusthreshold),
                         [&](HighsFloat cplusthreshold, HighsInt i) {
                           return cplusthreshold > vals[i];
                         }) -
        cover.begin();

    HighsCD0uble mcplus = m[cplusend];
    if (i < cplusend) mcplus -= mju;

    HighsFloat jlVal = HighsFloat(mcplus + eta * vals[j]);

    if (jlVal > bestlVal || (!atUpper && bestlAtUpper)) {
      lpos = i;
      bestlCplusend = cplusend;
      bestlVal = jlVal;
      bestlAtUpper = atUpper;
    }
  }

  if (lpos == -1) return false;

  l = cover[lpos];
  HighsCD0uble al = vals[l];
  HighsFloat upperl = upper[l];
  HighsCD0uble mlu = upperl * al;
  HighsCD0uble mu = mlu - lambda;

  a.resize(bestlCplusend);
  cover.resize(bestlCplusend);
  u.resize(bestlCplusend + 1);
  m.resize(bestlCplusend + 1);

  if (lpos < bestlCplusend) {
    a.erase(a.begin() + lpos);
    cover.erase(cover.begin() + lpos);
    u.erase(u.begin() + lpos + 1);
    m.erase(m.begin() + lpos + 1);
    for (HighsInt i = lpos + 1; i < bestlCplusend; ++i) {
      u[i] -= upperl;
      m[i] -= mlu;
    }
  }

  HighsInt cplussize = a.size();

  assert(mu > 10 * feastol);

  HighsFloat mudival = HighsFloat(mu / al);
  HighsFloat eta = ceil(mudival);
  HighsCD0uble r = mu - floor(mudival) * HighsCD0uble(al);
  // we multiply with r and it is important that it does not flip the sign
  // so we safe guard against tiny numerical errors here
  if (r < 0) r = 0;

  HighsCD0uble ulminusetaplusone = HighsCD0uble(upperl) - eta + 1.0;
  HighsCD0uble cplusthreshold = ulminusetaplusone * al;

  HighsInt kmin = floor(eta - upperl - 0.5);

  auto phi_l = [&](HighsFloat a) {
    assert(a < 0);

    int64_t k = std::min(int64_t(a / HighsFloat(al)), int64_t(-1));

    for (; k >= kmin; --k) {
      if (a >= k * al + r) {
        assert(a < (k + 1) * al);
        return HighsFloat(a - (k + 1) * r);
      }

      if (a >= k * al) {
        assert(a < k * al + r);
        return HighsFloat(k * (al - r));
      }
    }

    assert(a <= -lambda + epsilon);
    return HighsFloat(kmin * (al - r));
  };

  int64_t kmax = floor(upperl - eta + 0.5);

  auto gamma_l = [&](HighsFloat z) {
    assert(z > 0);
    for (HighsInt i = 0; i < cplussize; ++i) {
      HighsInt upperi = upper[cover[i]];

      for (HighsInt h = 0; h <= upperi; ++h) {
        HighsCD0uble mih = m[i] + h * a[i];
        HighsCD0uble uih = u[i] + h;
        HighsCD0uble mihplusdeltai = mih + a[i] - cplusthreshold;
        if (z <= mihplusdeltai) {
          assert(mih <= z);
          return HighsFloat(uih * ulminusetaplusone * (al - r));
        }

        int64_t k = ((int64_t)(HighsFloat)((z - mihplusdeltai) / al)) - 1;
        for (; k <= kmax; ++k) {
          if (z <= mihplusdeltai + k * al + r) {
            assert(mihplusdeltai + k * al <= z);
            return HighsFloat((uih * ulminusetaplusone + k) * (al - r));
          }

          if (z <= mihplusdeltai + (k + 1) * al) {
            assert(mihplusdeltai + k * al + r <= z);
            return HighsFloat((uih * ulminusetaplusone) * (al - r) + z - mih -
                          a[i] + cplusthreshold - (k + 1) * r);
          }
        }
      }
    }

    int64_t p = ((int64_t)(HighsFloat)((z - m[cplussize]) / al)) - 1;
    for (;; ++p) {
      if (z <= m[cplussize] + p * al + r) {
        assert(m[cplussize] + p * al <= z);
        return HighsFloat((u[cplussize] * ulminusetaplusone + p) * (al - r));
      }

      if (z <= m[cplussize] + (p + 1) * al) {
        assert(m[cplussize] + p * al + r <= z);
        return HighsFloat((u[cplussize] * ulminusetaplusone) * (al - r) + z -
                      m[cplussize] - (p + 1) * r);
      }
    }
  };

  rhs = (HighsCD0uble(upperl) - eta) * r - lambda;
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt i = 0; i != rowlen; ++i) {
    if (vals[i] == 0.0) continue;
    if (!isintegral[i]) {
      if (vals[i] < 0.0)
        integralSupport = false;
      else
        vals[i] = 0.0;
      continue;
    }

    if (coverflag[i]) {
      vals[i] = -phi_l(-vals[i]);
      rhs += vals[i] * upper[i];
    } else {
      vals[i] = gamma_l(vals[i]);
    }
  }

  return true;
}

static HighsFloat fast_floor(HighsFloat x) { return (int64_t)x - (x < (int64_t)x); }

bool HighsCutGeneration::cmirCutGenerationHeuristic(HighsFloat minEfficacy,
                                                    bool onlyInitialCMIRScale) {
  using std::abs;
  using std::floor;
  using std::max;
  using std::sqrt;

  HighsFloat continuouscontribution = 0.0;
  HighsFloat continuoussqrnorm = 0.0;

  deltas.clear();
  deltas.reserve(rowlen + 3);
  integerinds.clear();
  integerinds.reserve(rowlen);

  HighsFloat maxabsdelta = 0.0;
  constexpr HighsFloat maxCMirScale = 1e6;
  constexpr HighsFloat f0min = 0.005;
  constexpr HighsFloat f0max = 0.995;

  complementation.resize(rowlen);

  for (HighsInt i = 0; i != rowlen; ++i) {
    if (isintegral[i]) {
      integerinds.push_back(i);

      if (upper[i] < 2 * solval[i]) {
        complementation[i] = 1 - complementation[i];
        rhs -= upper[i] * vals[i];
        vals[i] = -vals[i];
        solval[i] = upper[i] - solval[i];
      }

      if (onlyInitialCMIRScale) continue;

      if (solval[i] > feastol) {
        HighsFloat delta = abs(vals[i]);
        if (delta <= 1e-4 || delta == maxabsdelta) continue;
        maxabsdelta = max(maxabsdelta, delta);
        deltas.push_back(delta);
      }
    } else {
      continuouscontribution += vals[i] * solval[i];
      continuoussqrnorm += vals[i] * vals[i];
    }
  }

  if (continuoussqrnorm == 0 && deltas.size() > 1) {
    HighsFloat intScale = HighsIntegers::integralScale(deltas, feastol, kHighsTiny);

    if (intScale != 0.0 && intScale <= 1e4) {
      HighsFloat scalrhs = HighsFloat(rhs) * intScale;
      HighsFloat downrhs = fast_floor(scalrhs);

      HighsFloat f0 = scalrhs - downrhs;
      if (f0 >= f0min && f0 <= f0max) deltas.push_back(1.0 / intScale);
    }
  }

  deltas.push_back(std::min(1.0, initialScale));
  if (!onlyInitialCMIRScale)
    deltas.push_back(maxabsdelta + std::min(1.0, initialScale));

  if (deltas.empty()) return false;

  pdqsort(deltas.begin(), deltas.end());
  HighsFloat curdelta = deltas[0];
  for (size_t i = 1; i < deltas.size(); ++i) {
    if (deltas[i] - curdelta <= 10 * feastol)
      deltas[i] = 0.0;
    else
      curdelta = deltas[i];
  }

  deltas.erase(std::remove(deltas.begin(), deltas.end(), 0.0), deltas.end());
  HighsFloat bestdelta = -1;
  HighsFloat bestefficacy = minEfficacy;

  for (HighsFloat delta : deltas) {
    HighsFloat scale = 1.0 / delta;
    HighsFloat scalrhs = HighsFloat(rhs) * scale;
    HighsFloat downrhs = fast_floor(scalrhs);

    HighsFloat f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) continue;
    HighsFloat oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) continue;

    HighsFloat sqrnorm = scale * scale * continuoussqrnorm;
    HighsFloat viol = scale * continuouscontribution * oneoveroneminusf0 - downrhs;

    for (HighsInt j : integerinds) {
      HighsFloat scalaj = vals[j] * scale;
      HighsFloat downaj = fast_floor(scalaj + kHighsTiny);
      HighsFloat fj = scalaj - downaj;
      HighsFloat aj = downaj + std::max(0.0, fj - f0);

      viol += aj * solval[j];
      sqrnorm += aj * aj;
    }

    HighsFloat efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  /* try if multiplying best delta by 2 4 or 8 gives a better efficacy */
  for (HighsInt k = 1; !onlyInitialCMIRScale && k <= 3; ++k) {
    HighsFloat delta = bestdelta * (1 << k);
    HighsFloat scale = 1.0 / delta;
    HighsFloat scalrhs = HighsFloat(rhs) * scale;
    HighsFloat downrhs = fast_floor(scalrhs);
    HighsFloat f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) continue;

    HighsFloat oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) continue;

    HighsFloat sqrnorm = scale * scale * continuoussqrnorm;
    HighsFloat viol = scale * continuouscontribution * oneoveroneminusf0 - downrhs;

    for (HighsInt j : integerinds) {
      HighsFloat scalaj = vals[j] * scale;
      HighsFloat downaj = fast_floor(scalaj + kHighsTiny);
      HighsFloat fj = scalaj - downaj;
      HighsFloat aj = downaj + std::max(0.0, fj - f0);

      viol += aj * solval[j];
      sqrnorm += aj * aj;
    }

    HighsFloat efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestdelta = delta;
      bestefficacy = efficacy;
    }
  }

  if (bestdelta == -1) return false;

  // try to flip complementation of integers to increase efficacy
  for (HighsInt k : integerinds) {
    if (upper[k] == kHighsInf) continue;
    if (solval[k] <= feastol) continue;

    complementation[k] = 1 - complementation[k];
    solval[k] = upper[k] - solval[k];
    rhs -= upper[k] * vals[k];
    vals[k] = -vals[k];

    HighsFloat delta = bestdelta;
    HighsFloat scale = 1.0 / delta;
    HighsFloat scalrhs = HighsFloat(rhs) * scale;
    HighsFloat downrhs = fast_floor(scalrhs);

    HighsFloat f0 = scalrhs - downrhs;
    if (f0 < f0min || f0 > f0max) {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];

      continue;
    }

    HighsFloat oneoveroneminusf0 = 1.0 / (1.0 - f0);
    if (oneoveroneminusf0 > maxCMirScale) {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
      continue;
    }

    HighsFloat sqrnorm = scale * scale * continuoussqrnorm;
    HighsFloat viol = scale * continuouscontribution * oneoveroneminusf0 - downrhs;

    for (HighsInt j : integerinds) {
      HighsFloat scalaj = vals[j] * scale;
      HighsFloat downaj = fast_floor(scalaj + kHighsTiny);
      HighsFloat fj = scalaj - downaj;
      HighsFloat aj = downaj + std::max(0.0, fj - f0);

      viol += aj * solval[j];
      sqrnorm += aj * aj;
    }

    HighsFloat efficacy = viol / sqrt(sqrnorm);
    if (efficacy > bestefficacy) {
      bestefficacy = efficacy;
    } else {
      complementation[k] = 1 - complementation[k];
      solval[k] = upper[k] - solval[k];
      rhs -= upper[k] * vals[k];
      vals[k] = -vals[k];
    }
  }

  HighsCD0uble scale = 1.0 / HighsCD0uble(bestdelta);
  HighsCD0uble scalrhs = rhs * scale;
  HighsFloat downrhs = floor(HighsFloat(scalrhs));

  HighsCD0uble f0 = scalrhs - downrhs;
  HighsCD0uble oneoveroneminusf0 = 1.0 / (1.0 - f0);

  rhs = downrhs * bestdelta;
  integralSupport = true;
  integralCoefficients = false;
  for (HighsInt j = 0; j != rowlen; ++j) {
    if (vals[j] == 0.0) continue;
    if (!isintegral[j]) {
      if (vals[j] > 0.0)
        vals[j] = 0.0;
      else {
        vals[j] = HighsFloat(vals[j] * oneoveroneminusf0);
        integralSupport = false;
      }
    } else {
      HighsCD0uble scalaj = scale * vals[j];
      HighsFloat downaj = floor(HighsFloat(scalaj + kHighsTiny));
      HighsCD0uble fj = scalaj - downaj;
      HighsCD0uble aj = downaj;
      if (fj > f0) aj += fj - f0;

      vals[j] = HighsFloat(aj * bestdelta);
    }
  }

  return true;
}

bool HighsCutGeneration::postprocessCut() {
  // right hand sides slightly below zero are likely due to numerical errors and
  // can cause numerical troubles with scaling, so set them to zero
  if (rhs < 0 && rhs > -epsilon) rhs = 0;

  if (integralSupport && integralCoefficients) {
    // if the cut is known to be integral no postprocessing is needed and we
    // simply remove zero coefficients
    for (HighsInt i = rowlen - 1; i >= 0; --i) {
      if (vals[i] == 0.0) {
        --rowlen;
        inds[i] = inds[rowlen];
        vals[i] = vals[rowlen];
      }
    }
    return true;
  }

  HighsDomain& globaldomain = lpRelaxation.getMipSolver().mipdata_->domain;
  // determine maximal absolute coefficient
  HighsFloat maxAbsValue = 0.0;
  for (HighsInt i = 0; i != rowlen; ++i)
    maxAbsValue = std::max(std::abs(vals[i]), maxAbsValue);

  // determine minimal allowed coefficient
  HighsFloat minCoefficientValue = 100 * feastol * std::max(maxAbsValue, 1e-3);

  // remove small coefficients and check whether the remaining support is
  // integral
  integralSupport = true;
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0) continue;
    if (std::abs(vals[i]) <= minCoefficientValue) {
      if (vals[i] < 0) {
        HighsFloat ub = globaldomain.col_upper_[inds[i]];
        if (ub == kHighsInf)
          return false;
        else
          rhs -= ub * vals[i];
      } else {
        HighsFloat lb = globaldomain.col_lower_[inds[i]];
        if (lb == -kHighsInf)
          return false;
        else
          rhs -= lb * vals[i];
      }

      vals[i] = 0.0;
      continue;
    }

    if (integralSupport && !lpRelaxation.isColIntegral(inds[i]))
      integralSupport = false;
  }

  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    }
  }

  if (integralSupport) {
    // integral support -> determine scale to make all coefficients integral
    HighsFloat intscale =
        HighsIntegers::integralScale(vals, rowlen, feastol, epsilon);

    bool scaleSmallestValToOne = true;

    if (intscale != 0.0 &&
        intscale * std::max(1.0, maxAbsValue) <= (HighsFloat)(uint64_t{1} << 53)) {
      // A scale to make all value integral was found. The scale is only
      // rejected if it is in a range where not all integral values are
      // representable in HighsFloat precision anymore. Otherwise we want to always
      // use the scale to adjust the coefficients and right hand side for
      // numerical safety reasons. If the resulting integral values are too
      // large, however, we scale the cut down by shifting the exponent.
      rhs.renormalize();
      rhs *= intscale;
      maxAbsValue = std::round(maxAbsValue * intscale);
      for (HighsInt i = 0; i != rowlen; ++i) {
        HighsCD0uble scaleval = intscale * HighsCD0uble(vals[i]);
        HighsCD0uble intval = round(scaleval);
        HighsFloat delta = HighsFloat(scaleval - intval);

        vals[i] = (HighsFloat)intval;

        // if the coefficient would be strengthened by rounding, we add the
        // upperbound constraint to make it exactly integral instead and
        // therefore weaken the right hand side
        if (delta < 0.0) {
          HighsFloat ub = globaldomain.col_upper_[inds[i]];
          if (ub == kHighsInf) return false;

          rhs -= delta * ub;
        } else {
          HighsFloat lb = globaldomain.col_lower_[inds[i]];
          if (lb == -kHighsInf) return false;

          rhs -= delta * lb;
        }
      }

      // finally we can round down the right hand side. Therefore in most cases
      // small errors for which the upper bound constraints where used and the
      // right hand side was weakened, do not weaken the final cut.
      rhs = floor(rhs + feastol);

      if (intscale * maxAbsValue * feastol < 0.5) {
        // integral scale leads to small enough values, accept scale
        scaleSmallestValToOne = false;
        integralCoefficients = true;
      }
    }

    if (scaleSmallestValToOne) {
      // integral scale lead to very large coefficient values. We now shift the
      // exactly integral values down such that the smallest coefficient is
      // around 1
      HighsFloat minAbsValue = kHighsInf;
      for (HighsInt i = 0; i != rowlen; ++i)
        minAbsValue = std::min(std::abs(vals[i]), minAbsValue);

      int expshift;
      std::frexp(minAbsValue - epsilon, &expshift);
      expshift = -expshift;

      rhs = std::ldexp((HighsFloat)rhs, expshift);

      for (HighsInt i = 0; i != rowlen; ++i)
        vals[i] = std::ldexp(vals[i], expshift);
    }
  } else {
    // the support is not integral, scale cut to have the largest coefficient
    // around 1.0
    int expshift;
    std::frexp(maxAbsValue, &expshift);
    expshift = -expshift;
    rhs = std::ldexp((HighsFloat)rhs, expshift);

    for (HighsInt i = 0; i != rowlen; ++i)
      vals[i] = std::ldexp(vals[i], expshift);
  }

  return true;
}

bool HighsCutGeneration::preprocessBaseInequality(bool& hasUnboundedInts,
                                                  bool& hasGeneralInts,
                                                  bool& hasContinuous) {
  // preprocess the inequality before cut generation
  // 1. Determine the maximal activity to check for trivial redundancy and
  // tighten coefficients
  // 2. Check for presence of continuous variables and unbounded integers as not
  // all methods for cut generation are applicable in that case
  // 3. Remove coefficients that are below the feasibility tolerance to avoid
  // numerical troubles, use bound constraints to cancel them and
  // reject base inequalities where that is not possible due to unbounded
  // variables
  hasUnboundedInts = false;
  hasContinuous = false;
  hasGeneralInts = false;
  HighsInt numZeros = 0;

  HighsFloat maxact = -feastol;
  HighsFloat maxAbsVal = 0;
  for (HighsInt i = 0; i < rowlen; ++i)
    maxAbsVal = std::max(std::abs(vals[i]), maxAbsVal);

  int expshift = 0;
  std::frexp(maxAbsVal, &expshift);
  expshift = -expshift;
  initialScale = std::ldexp(1.0, expshift);
  rhs *= initialScale;
  for (HighsInt i = 0; i < rowlen; ++i) vals[i] = std::ldexp(vals[i], expshift);

  isintegral.resize(rowlen);
  for (HighsInt i = 0; i != rowlen; ++i) {
    // we do not want to have integral variables with small coefficients as this
    // may lead to numerical instabilities during cut generation
    // Therefore we relax integral variables with small coefficients to
    // continuous ones because they still might have a non-negligible
    // contribution e.g. when they come from integral rows with coefficients on
    // the larger side. When we relax them to continuous variables they will be
    // complemented so that their solution value is closest to zero and then
    // will be relaxed if their value is positive or their maximal contribution
    // is below feasibility tolerance.
    isintegral[i] =
        lpRelaxation.isColIntegral(inds[i]) && std::abs(vals[i]) > 10 * feastol;

    if (!isintegral[i]) {
      if (upper[i] - solval[i] < solval[i]) {
        if (complementation.empty()) complementation.resize(rowlen);

        complementation[i] = 1 - complementation[i];
        rhs -= upper[i] * vals[i];
        vals[i] = -vals[i];
      }

      // relax positive continuous variables and those with small contributions
      if (vals[i] > 0 || std::abs(vals[i]) * upper[i] <= 10 * feastol) {
        // printf("remove: vals[i] = %g  upper[i] = %g\n", vals[i], upper[i]);
        if (vals[i] < 0) {
          if (upper[i] == kHighsInf) return false;
          rhs -= vals[i] * upper[i];
        }

        ++numZeros;
        vals[i] = 0.0;
        continue;
      }

      hasContinuous = true;
      // if (lpRelaxation.isColIntegral(inds[i]))
      //   printf("vals[i] = %g  upper[i] = %g\n", vals[i], upper[i]);

      if (vals[i] > 0) {
        if (upper[i] == kHighsInf)
          maxact = kHighsInf;
        else
          maxact += vals[i] * upper[i];
      }
    } else {
      if (upper[i] == kHighsInf) {
        hasUnboundedInts = true;
        hasGeneralInts = true;
      } else if (upper[i] != 1.0) {
        hasGeneralInts = true;
      }

      if (vals[i] > 0) maxact += vals[i] * upper[i];
    }
  }

  HighsInt maxLen = 100 + 0.15 * (lpRelaxation.numCols());

  if (rowlen - numZeros > maxLen) {
    HighsInt numCancel = rowlen - numZeros - maxLen;
    std::vector<HighsInt> cancelNzs;

    for (HighsInt i = 0; i != rowlen; ++i) {
      HighsFloat cancelSlack = vals[i] > 0 ? solval[i] : upper[i] - solval[i];
      if (cancelSlack <= feastol) cancelNzs.push_back(i);
    }

    if ((HighsInt)cancelNzs.size() < numCancel) return false;
    if ((HighsInt)cancelNzs.size() > numCancel)
      std::partial_sort(cancelNzs.begin(), cancelNzs.begin() + numCancel,
                        cancelNzs.end(), [&](HighsInt a, HighsInt b) {
                          return std::abs(vals[a]) < std::abs(vals[b]);
                        });

    for (HighsInt i = 0; i < numCancel; ++i) {
      HighsInt j = cancelNzs[i];

      if (vals[j] < 0) {
        rhs -= vals[j] * upper[j];
      } else
        maxact -= vals[j] * upper[j];

      vals[j] = 0.0;
    }

    numZeros += numCancel;
  }

  if (numZeros != 0) {
    // remove zeros in place
    if (complementation.empty()) {
      for (HighsInt i = rowlen - 1; i >= 0; --i) {
        if (vals[i] == 0.0) {
          --rowlen;
          inds[i] = inds[rowlen];
          vals[i] = vals[rowlen];
          upper[i] = upper[rowlen];
          solval[i] = solval[rowlen];
          isintegral[i] = isintegral[rowlen];
          if (--numZeros == 0) break;
        }
      }
    } else {
      for (HighsInt i = rowlen - 1; i >= 0; --i) {
        if (vals[i] == 0.0) {
          --rowlen;
          inds[i] = inds[rowlen];
          vals[i] = vals[rowlen];
          upper[i] = upper[rowlen];
          solval[i] = solval[rowlen];
          isintegral[i] = isintegral[rowlen];
          complementation[i] = complementation[rowlen];
          if (--numZeros == 0) break;
        }
      }
    }
  }

  return maxact > rhs;
}

#if 0
static void checkNumerics(const HighsFloat* vals, HighsInt len, HighsFloat rhs) {
  HighsFloat maxAbsCoef = 0.0;
  HighsFloat minAbsCoef = kHighsInf;
  HighsCD0uble sqrnorm = 0;
  for (HighsInt i = 0; i < len; ++i) {
    sqrnorm += vals[i] * vals[i];
    maxAbsCoef = std::max(std::abs(vals[i]), maxAbsCoef);
    minAbsCoef = std::min(std::abs(vals[i]), minAbsCoef);
  }

  HighsFloat norm = HighsFloat(sqrt(sqrnorm));

  // printf("length: %" HIGHSINT_FORMAT
  //       ", minCoef: %g, maxCoef, %g, norm %g, rhs: %g, dynamism=%g\n",
  //       len, minAbsCoef, maxAbsCoef, norm, rhs, maxAbsCoef / minAbsCoef);
}
#endif

bool HighsCutGeneration::generateCut(HighsTransformedLp& transLp,
                                     std::vector<HighsInt>& inds_,
                                     std::vector<HighsFloat>& vals_, HighsFloat& rhs_,
                                     bool onlyInitialCMIRScale) {
#if 0
  if (vals_.size() > 1) {
    std::vector<HighsInt> indsCheck_ = inds_;
    std::vector<HighsFloat> valsCheck_ = vals_;
    HighsFloat tmprhs_ = rhs_;
    bool intsPositive = true;
    if (!transLp.transform(valsCheck_, upper, solval, indsCheck_, tmprhs_,
                           intsPositive))
      return false;

    rowlen = indsCheck_.size();
    this->inds = indsCheck_.data();
    this->vals = valsCheck_.data();
    this->rhs = tmprhs_;
    complementation.clear();
    bool hasUnboundedInts = false;
    bool hasGeneralInts = false;
    bool hasContinuous = false;
    // printf("before preprocessing of base inequality:\n");
    checkNumerics(vals, rowlen, HighsFloat(rhs));
    if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                  hasContinuous))
      return false;
    // printf("after preprocessing of base inequality:\n");
    checkNumerics(vals, rowlen, HighsFloat(rhs));

    tmprhs_ = (HighsFloat)rhs;
    valsCheck_.resize(rowlen);
    indsCheck_.resize(rowlen);
    if (!transLp.untransform(valsCheck_, indsCheck_, tmprhs_)) return false;

    // printf("after untransform of base inequality:\n");
    checkNumerics(vals, rowlen, HighsFloat(rhs));

    // finally check whether the cut is violated
    rowlen = indsCheck_.size();
    inds = indsCheck_.data();
    vals = valsCheck_.data();
    lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
        inds, vals, rowlen, tmprhs_);
  }
#endif

  bool intsPositive = true;
  if (!transLp.transform(vals_, upper, solval, inds_, rhs_, intsPositive))
    return false;

  rowlen = inds_.size();
  this->inds = inds_.data();
  this->vals = vals_.data();
  this->rhs = rhs_;
  complementation.clear();
  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;
  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  // it can happen that there is an unbounded integer variable during the
  // transform call so that the integers are not tranformed to positive values.
  // Now the call to preprocessBaseInequality may have removed the unbounded
  // integer, e.g. due to a small coefficient value, so that we can still use
  // the lifted inequalities instead of cmir. We need to make sure, however,
  // that the cut values are transformed to positive coefficients first, which
  // we do below.
  if (!hasUnboundedInts && !intsPositive) {
    complementation.resize(rowlen);

    for (HighsInt i = 0; i != rowlen; ++i) {
      if (vals[i] > 0 || !isintegral[i]) continue;

      complementation[i] = 1 - complementation[i];
      rhs -= upper[i] * vals[i];
      vals[i] = -vals[i];
      solval[i] = upper[i] - solval[i];
    }
  }

  const HighsFloat minEfficacy = 10 * feastol;

  if (hasUnboundedInts) {
    if (!cmirCutGenerationHeuristic(minEfficacy, onlyInitialCMIRScale))
      return false;
  } else {
    // 1. Determine a cover, cover does not need to be minimal as neither of
    // the
    //    lifting functions have minimality of the cover as necessary facet
    //    condition
    std::vector<HighsFloat> tmpVals(vals, vals + rowlen);
    std::vector<HighsInt> tmpInds(inds, inds + rowlen);
    HighsCD0uble tmpRhs = rhs;
    bool success = false;
    do {
      if (!determineCover()) break;

      // 2. use superadditive lifting function depending on structure of base
      //    inequality:
      //    We have 3 lifting functions available for pure binary knapsack sets,
      //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
      if (!hasContinuous && !hasGeneralInts) {
        separateLiftedKnapsackCover();
        success = true;
      } else if (hasGeneralInts) {
        success = separateLiftedMixedIntegerCover();
      } else {
        assert(hasContinuous);
        assert(!hasGeneralInts);
        success = separateLiftedMixedBinaryCover();
      }
    } while (false);

    HighsFloat minMirEfficacy = minEfficacy;
    if (success) {
      HighsFloat violation = -HighsFloat(rhs);
      HighsFloat sqrnorm = 0.0;

      for (HighsInt i = 0; i < rowlen; ++i) {
        sqrnorm += vals[i] * vals[i];
        violation += vals[i] * solval[i];
      }

      HighsFloat efficacy = violation / std::sqrt(sqrnorm);
      if (efficacy <= minEfficacy) {
        success = false;
        rhs = tmpRhs;
      } else {
        minMirEfficacy += efficacy;
        if (!complementation.empty()) {
          // remove the complementation if it exists, so that the values stored
          // values are uncomplemented
          for (HighsInt i = 0; i != rowlen; ++i) {
            if (complementation[i]) {
              rhs -= upper[i] * vals[i];
              vals[i] = -vals[i];
              solval[i] = upper[i] - solval[i];
            }
          }
        }
        std::swap(tmpRhs, rhs);
      }
    }

    inds = tmpInds.data();
    vals = tmpVals.data();

    bool cmirSuccess =
        cmirCutGenerationHeuristic(minMirEfficacy, onlyInitialCMIRScale);

    if (cmirSuccess) {
      // take the cmir cut as it is better
      inds_.swap(tmpInds);
      vals_.swap(tmpVals);
      inds = inds_.data();
      vals = vals_.data();
    } else if (success) {
      // take the previous lifted cut as cmir could not improve
      // as we already removed the complementation we simply clear
      // the vector if altered by the cmir routine and restore the old
      // right hand side and values
      rhs = tmpRhs;
      complementation.clear();
      inds = inds_.data();
      vals = vals_.data();
    } else
      // neither cmir nor lifted cut successful
      return false;
  }

  if (!complementation.empty()) {
    // remove the complementation if exists
    for (HighsInt i = 0; i != rowlen; ++i) {
      if (complementation[i]) {
        rhs -= upper[i] * vals[i];
        vals[i] = -vals[i];
      }
    }
  }

  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    }
  }

  // transform the cut back into the original space, i.e. remove the bound
  // substitution and replace implicit slack variables
  rhs_ = (HighsFloat)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);
  if (!transLp.untransform(vals_, inds_, rhs_)) return false;

  rowlen = inds_.size();
  inds = inds_.data();
  vals = vals_.data();
  rhs = rhs_;

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(inds, vals,
                                                               rowlen, rhs_);
  // apply cut postprocessing including scaling and removal of small
  // coeffiicents
  if (!postprocessCut()) return false;
  rhs_ = (HighsFloat)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), rowlen, rhs_);

  // finally determine the violation of the cut in the original space
  HighsCD0uble violation = -rhs_;
  const auto& sol = lpRelaxation.getSolution().col_value;
  for (HighsInt i = 0; i != rowlen; ++i) violation += sol[inds[i]] * vals_[i];

  if (violation <= 10 * feastol) return false;

  lpRelaxation.getMipSolver().mipdata_->domain.tightenCoefficients(
      inds, vals, rowlen, rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                     vals_.data(), inds_.size(), rhs_,
                                     integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::generateConflict(HighsDomain& localdomain,
                                          std::vector<HighsInt>& proofinds,
                                          std::vector<HighsFloat>& proofvals,
                                          HighsFloat& proofrhs) {
  this->inds = proofinds.data();
  this->vals = proofvals.data();
  this->rhs = proofrhs;
  rowlen = proofinds.size();

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds, vals, rowlen, proofrhs);

  complementation.assign(rowlen, 0);

  upper.resize(rowlen);
  solval.resize(rowlen);

  HighsDomain& globaldomain = lpRelaxation.getMipSolver().mipdata_->domain;
  for (HighsInt i = 0; i != rowlen; ++i) {
    HighsInt col = inds[i];

    upper[i] = globaldomain.col_upper_[col] - globaldomain.col_lower_[col];

    solval[i] = vals[i] < 0 ? std::min(globaldomain.col_upper_[col],
                                       localdomain.col_upper_[col])
                            : std::max(globaldomain.col_lower_[col],
                                       localdomain.col_lower_[col]);
    if (vals[i] < 0 && globaldomain.col_upper_[col] != kHighsInf) {
      rhs -= globaldomain.col_upper_[col] * vals[i];
      vals[i] = -vals[i];
      complementation[i] = 1;

      solval[i] = globaldomain.col_upper_[col] - solval[i];
    } else {
      rhs -= globaldomain.col_lower_[col] * vals[i];
      complementation[i] = 0;
      solval[i] = solval[i] - globaldomain.col_lower_[col];
    }
  }

  bool hasUnboundedInts = false;
  bool hasGeneralInts = false;
  bool hasContinuous = false;

  if (!preprocessBaseInequality(hasUnboundedInts, hasGeneralInts,
                                hasContinuous))
    return false;

  if (hasUnboundedInts) {
    if (!cmirCutGenerationHeuristic(-kHighsInf)) return false;
  } else {
    // 1. Determine a cover, cover does not need to be minimal as neither of
    // the
    //    lifting functions have minimality of the cover as necessary facet
    //    condition
    std::vector<HighsFloat> tmpVals(vals, vals + rowlen);
    std::vector<HighsInt> tmpInds(inds, inds + rowlen);
    std::vector<uint8_t> tmpComplementation(complementation);
    HighsCD0uble tmpRhs = rhs;
    bool success = false;
    do {
      if (!determineCover(false)) break;

      // 2. use superadditive lifting function depending on structure of base
      //    inequality:
      //    We have 3 lifting functions available for pure binary knapsack sets,
      //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
      if (!hasContinuous && !hasGeneralInts) {
        separateLiftedKnapsackCover();
        success = true;
      } else if (hasGeneralInts) {
        success = separateLiftedMixedIntegerCover();
      } else {
        assert(hasContinuous);
        assert(!hasGeneralInts);
        success = separateLiftedMixedBinaryCover();
      }
    } while (false);

    HighsFloat minEfficacy = -kHighsInf;
    if (success) {
      HighsFloat violation = -HighsFloat(rhs);
      HighsFloat sqrnorm = 0.0;

      for (HighsInt i = 0; i < rowlen; ++i) {
        violation += vals[i] * solval[i];
        sqrnorm += vals[i] * vals[i];
      }

      minEfficacy = violation / std::sqrt(sqrnorm);
      minEfficacy += feastol;
      std::swap(tmpRhs, rhs);
    }

    inds = tmpInds.data();
    vals = tmpVals.data();

    bool cmirSuccess = cmirCutGenerationHeuristic(minEfficacy);

    if (cmirSuccess) {
      // take the cmir cut as it is better
      proofinds.swap(tmpInds);
      proofvals.swap(tmpVals);
      inds = proofinds.data();
      vals = proofvals.data();
    } else if (success) {
      // take the previous lifted cut as cmir could not improve
      // we restore the old complementation vector, right hand side, and values
      rhs = tmpRhs;
      complementation.swap(tmpComplementation);
      inds = proofinds.data();
      vals = proofvals.data();
    } else
      // neither cmir nor lifted cut successful
      return false;
  }

  // remove the complementation
  if (!complementation.empty()) {
    for (HighsInt i = 0; i != rowlen; ++i) {
      if (complementation[i]) {
        rhs -= globaldomain.col_upper_[inds[i]] * vals[i];
        vals[i] = -vals[i];
      } else
        rhs += globaldomain.col_lower_[inds[i]] * vals[i];
    }
  }

  // apply cut postprocessing including scaling and removal of small
  // coefficients
  if (!postprocessCut()) return false;

  proofvals.resize(rowlen);
  proofinds.resize(rowlen);
  proofrhs = (HighsFloat)rhs;

  bool cutintegral = integralSupport && integralCoefficients;

  lpRelaxation.getMipSolver().mipdata_->domain.tightenCoefficients(
      proofinds.data(), proofvals.data(), rowlen, proofrhs);

  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(),
                                     proofinds.data(), proofvals.data(), rowlen,
                                     proofrhs, cutintegral, true, true, true);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}

bool HighsCutGeneration::finalizeAndAddCut(std::vector<HighsInt>& inds_,
                                           std::vector<HighsFloat>& vals_,
                                           HighsFloat& rhs_) {
  complementation.clear();
  rowlen = inds_.size();
  inds = inds_.data();
  vals = vals_.data();
  rhs = rhs_;

  integralSupport = true;
  integralCoefficients = false;
  // remove zeros in place
  for (HighsInt i = rowlen - 1; i >= 0; --i) {
    if (vals[i] == 0.0) {
      --rowlen;
      inds[i] = inds[rowlen];
      vals[i] = vals[rowlen];
    } else {
      integralSupport &= lpRelaxation.isColIntegral(inds[i]);
    }
  }

  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(inds, vals,
                                                               rowlen, rhs_);
  // apply cut postprocessing including scaling and removal of small
  // coeffiicents
  if (!postprocessCut()) return false;
  rhs_ = (HighsFloat)rhs;
  vals_.resize(rowlen);
  inds_.resize(rowlen);

  lpRelaxation.getMipSolver().mipdata_->debugSolution.checkCut(
      inds_.data(), vals_.data(), rowlen, rhs_);

  // finally determine the violation of the cut in the original space
  HighsCD0uble violation = -rhs_;
  const auto& sol = lpRelaxation.getSolution().col_value;
  for (HighsInt i = 0; i != rowlen; ++i) violation += sol[inds[i]] * vals_[i];

  if (violation <= 10 * feastol) return false;

  lpRelaxation.getMipSolver().mipdata_->domain.tightenCoefficients(
      inds, vals, rowlen, rhs_);

  // if the cut is violated by a small factor above the feasibility
  // tolerance, add it to the cutpool
  HighsInt cutindex = cutpool.addCut(lpRelaxation.getMipSolver(), inds_.data(),
                                     vals_.data(), inds_.size(), rhs_,
                                     integralSupport && integralCoefficients);

  // only return true if cut was accepted by the cutpool, i.e. not a duplicate
  // of a cut already in the pool
  return cutindex != -1;
}
