#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Prototypes of functions defined in fast.cpp (same shared library).
double coalAlphaSumSortedC(const NumericVector& leaves, const NumericVector& nodes);
double coalDeltaC(const NumericVector& leaves, const NumericVector& nodes, double alpha,
                  double oldval, double newval);
double localTermsC(const NumericMatrix& tab, IntegerVector nodes1, double mu, double sigma,
                   double minbralen, int model, bool useRec);
double localTermsVecC(const NumericMatrix& tab, const std::vector<int>& nodes1,
                      double mu, double sigma, double minbralen,
                      int model, bool useRec);
double localTermsArrC(const NumericMatrix& tab, const int* nodes1, int n1,
                      double mu, double sigma, double minbralen,
                      int model, bool useRec);
double fullLikC(const NumericMatrix& tab, double mu, double sigma,
                double minbralen, int model, bool useRec);
int findCurrootC(const NumericMatrix& edge, NumericVector rootchildren, int rootnode);

// Prototypes of functions defined in probs.cpp (same shared library).
double coalpriorC(const NumericVector& leaves, const NumericVector& nodes, double alpha);
void changeinorderedvec(NumericVector vec, double old, double n);

// ---------------------------------------------------------------------------
// Whole-MCMC-loop implementation of the main.R iteration, in three modes:
//   mode 0 = "slow"  : bit-identical replica of the original R loop
//   mode 1 = "exact" : bit-identical to the original R loop, using the
//                      prefix-sum prior cache and children-structure lookups
//   mode 2 = "fast"  : mathematically identical but incrementally accumulated
//                      (localized likelihood terms, incremental prior delta)
// RNG consumption matches the R loop draw for draw (R::rnorm/R::runif/...
// call R's own RNG), so identical seeds give identical streams.
// mode>0 requires modelcode>0. showProgress runs must use the R loop.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List mcmcLoopC(NumericMatrix tab, const NumericMatrix& edge, int n, int nbIts, int thin,
               NumericVector orderedleafdates, NumericVector orderednodedates,
               IntegerVector misDates, const NumericMatrix& rangedate,
               double mu, double sigma, double alpha, double l, double p,
               double sdMu, double sdSigma, double sdDates,
               bool updateMu, bool updateSigma, bool updateAlpha,
               bool updateRootAll, bool updateRootBranch, bool mixedModel,
               bool tuning, bool useCoalPrior, int model, bool useRec,
               double minbralen, int mode, NumericMatrix record) {
  int nrowtab = tab.nrow();
  int nint = nrowtab - n;
  double acceptance_target = 0.234;
  double delta = 1.0 / acceptance_target / (1.0 - acceptance_target);
  bool exact = (mode == 1);
  bool fast = (mode == 2);

  // children structure (1-based node -> list of 1-based children)
  std::vector<std::vector<int> > kids(nrowtab + 1);
  for (int x = 1; x <= nrowtab; ++x) {
    double f = tab(x - 1, 3);
    if (!NumericVector::is_na(f)) kids[(int)f].push_back(x);
  }
  auto rebuildKids = [&]() {
    for (int x = 1; x <= nrowtab; ++x) kids[x].clear();
    for (int x = 1; x <= nrowtab; ++x) {
      double f = tab(x - 1, 3);
      if (!NumericVector::is_na(f)) kids[(int)f].push_back(x);
    }
  };
  std::vector<int> rootch = kids[n + 1];
  NumericVector rootchN = NumericVector::create((double)rootch[0], (double)rootch[1]);
  int curroot = findCurrootC(edge, rootchN, n + 1);

  // Prefix-sum cache for the bit-identical coalescent prior (exact mode):
  // walk state after t merge steps, matching coalpriorC step for step.
  const double* LV = orderedleafdates.begin();
  const double* ND = orderednodedates.begin();
  int ne = orderedleafdates.size();
  int msteps = 2 * ne - 2;             // merge steps after the initial leaf
  std::vector<double> Ssum, PrevS;
  std::vector<int> I1s, I2s;
  std::vector<double> inc;             // P4: per-step increments of the current state
  std::vector<double> winInc;          // P4: proposal-state window increments (scratch)
  int X = 0;
  auto countGtL = [&](double val) {
    int lo = 0, hi = ne;
    while (lo < hi) { int mid = (lo + hi) / 2; if (LV[mid] > val) lo = mid + 1; else hi = mid; }
    return lo;
  };
  auto countGtN = [&](double val) {
    int lo = 0, hi = ne - 1;
    while (lo < hi) { int mid = (lo + hi) / 2; if (ND[mid] > val) lo = mid + 1; else hi = mid; }
    return lo;
  };
  // P4: counts of events >= val (decreasing order), for the identical-tail
  // boundary pos_lo.
  auto countGeL = [&](double val) {
    int lo = 0, hi = ne;
    while (lo < hi) { int mid = (lo + hi) / 2; if (LV[mid] >= val) lo = mid + 1; else hi = mid; }
    return lo;
  };
  auto countGeN = [&](double val) {
    int lo = 0, hi = ne - 1;
    while (lo < hi) { int mid = (lo + hi) / 2; if (ND[mid] >= val) lo = mid + 1; else hi = mid; }
    return lo;
  };
  // P1: lookup table of the coalescent weight k*(k-1)/(2*alpha). Every entry
  // is computed with the exact same expression coalpriorC evaluates inline,
  // so a table read is bit-identical to recomputing it, and the hot walk
  // drops one multiplication + one division per step. Rebuilt whenever alpha
  // changes (Gibbs update); k stays within [0, ne] for any valid walk over a
  // binary tree's dates. Only built/used in exact mode (COAL_STEP users).
  std::vector<double> coalw(ne + 1);
  auto buildCoalw = [&]() {
    for (int k = 0; k <= ne; ++k) coalw[k] = k * (k - 1.0) / (2.0 * alpha);
  };
  // P2: cache entries are written only every COAL_CK steps ("checkpoints").
  // X always sits on a checkpoint; a read at an arbitrary pos_hi restarts
  // from the checkpoint at or below it (at most COAL_CK-1 recomputed steps)
  // instead of storing all four arrays every step. This cuts the cache write
  // traffic ~COAL_CK-fold while the walked arithmetic - and therefore every
  // bit of p2 - is unchanged.
  // P4 changed the invalidation story: the proposal walk writes NO
  // checkpoints (the identical tail is served from the increment cache inc[]
  // instead), so every checkpoint <= c lies in the prefix identical to both
  // states. On accept X = c; on reject X is left untouched.
  #define COAL_CK 32
  #define COAL_CKPT(sumvar, t)                                                 \
    if ((t) % COAL_CK == 0) { Ssum[t] = sumvar; I1s[t] = i1; I2s[t] = i2; PrevS[t] = prev; }
  // One merge step of the coalpriorC walk (inlined, same arithmetic/order).
  // Consumes the newer of LV[i1], ND[i2] into the running sum *pp.
  #define COAL_STEP(pp)                                                        \
    do {                                                                       \
      int k = i1 - i2;                                                         \
      if (i1 < ne && LV[i1] > ND[i2]) {                                        \
        pp = pp - (coalw[k] * (prev - LV[i1]));                                \
        prev = LV[i1++];                                                       \
      } else {                                                                 \
        pp = pp - (coalw[k] * (prev - ND[i2]));                                \
        prev = ND[i2++];                                                       \
      }                                                                        \
    } while (0)
  // P4: same step, but additionally stores the per-step increment
  // inc[t] = coalw[k]*(prev - event) into the increment cache.
  #define COAL_STEP_INC(pp, t)                                                 \
    do {                                                                       \
      int k = i1 - i2;                                                         \
      if (i1 < ne && LV[i1] > ND[i2]) {                                        \
        double w = coalw[k] * (prev - LV[i1]);                                 \
        inc[(t)] = w;                                                          \
        pp = pp - w;                                                           \
        prev = LV[i1++];                                                       \
      } else {                                                                 \
        double w = coalw[k] * (prev - ND[i2]);                                 \
        inc[(t)] = w;                                                          \
        pp = pp - w;                                                           \
        prev = ND[i2++];                                                       \
      }                                                                        \
    } while (0)
  // P4: same step, stashing the increment into winInc (proposal-state window).
  #define COAL_STEP_WIN(pp, t)                                                 \
    do {                                                                       \
      int k = i1 - i2;                                                         \
      if (i1 < ne && LV[i1] > ND[i2]) {                                        \
        double w = coalw[k] * (prev - LV[i1]);                                 \
        winInc[(t)] = w;                                                       \
        pp = pp - w;                                                           \
        prev = LV[i1++];                                                       \
      } else {                                                                 \
        double w = coalw[k] * (prev - ND[i2]);                                 \
        winInc[(t)] = w;                                                       \
        pp = pp - w;                                                           \
        prev = ND[i2++];                                                       \
      }                                                                        \
    } while (0)
  if (exact) {
    Ssum.resize(msteps + 1); PrevS.resize(msteps + 1);
    I1s.resize(msteps + 1); I2s.resize(msteps + 1);
    inc.resize(msteps + 1); winInc.resize(msteps + 1);
    buildCoalw();
    double pp = -log(alpha) * (ne - 1);
    int i1 = 1, i2 = 0; double prev = LV[0];
    Ssum[0] = pp; I1s[0] = 1; I2s[0] = 0; PrevS[0] = prev;
    for (int t = 1; t <= msteps; ++t) {
      COAL_STEP_INC(pp, t);
      COAL_CKPT(pp, t);
    }
    X = (msteps / COAL_CK) * COAL_CK;
  }

  // Slow-mode prior: full likelihood on explicitly built 3/4-row mintabs,
  // exactly as the original R loop computes likelihood(mintab,...).
  auto slowNodeLocal = [&](int j1, double rowDate) {
    int ntab = (j1 > n + 1) ? 4 : 3;
    NumericMatrix mt(ntab, tab.ncol());
    const std::vector<int>& ch = kids[j1];
    for (int c = 0; c < (int)ch.size(); ++c)
      for (int cc = 0; cc < tab.ncol(); ++cc) mt(c, cc) = tab(ch[c] - 1, cc);
    int rr = (int)ch.size();
    if (j1 > n + 1) {
      int fa = (int)tab(j1 - 1, 3);
      for (int cc = 0; cc < tab.ncol(); ++cc) mt(rr, cc) = tab(fa - 1, cc);
      rr++;
    }
    for (int cc = 0; cc < tab.ncol(); ++cc) mt(rr, cc) = tab(j1 - 1, cc);
    mt(rr, 2) = rowDate;
    for (int q = 0; q < ntab; ++q) mt(q, 3) = (q < 2) ? 4.0 : R_NaReal;
    if (ntab == 4) mt(3, 3) = 3.0;
    mt(2, 1) = minbralen;
    return fullLikC(mt, mu, sigma, minbralen, model, useRec);
  };
  auto slowLeafLocal = [&](int j1, double rowDate) {
    NumericMatrix mt(2, tab.ncol());
    for (int cc = 0; cc < tab.ncol(); ++cc) mt(0, cc) = tab(j1 - 1, cc);
    mt(0, 2) = rowDate;
    int fa = (int)tab(j1 - 1, 3);
    for (int cc = 0; cc < tab.ncol(); ++cc) mt(1, cc) = tab(fa - 1, cc);
    mt(0, 3) = 2.0; mt(1, 3) = R_NaReal;
    mt(1, 1) = minbralen;
    return fullLikC(mt, mu, sigma, minbralen, model, useRec);
  };

  std::vector<int> terms;              // reused index set for local terms
  NumericVector rn(nint);              // preallocated per-iteration proposal noise
  for (int i = 1; i <= nbIts; ++i) {
    //Record
    if (i % thin == 0) {
      int row = i / thin - 1;
      for (int q = 0; q < nrowtab; ++q) record(row, q) = tab(q, 2);
      for (int q = 0; q < nrowtab; ++q) record(row, nrowtab + q) = tab(q, 3);
      for (int q = 0; q < nrowtab; ++q) record(row, 2 * nrowtab + q) = tab(q, 1);
      record(row, 3 * nrowtab) = l;
      record(row, 3 * nrowtab + 1) = mu;
      record(row, 3 * nrowtab + 2) = sigma;
      record(row, 3 * nrowtab + 3) = alpha;
      record(row, 3 * nrowtab + 4) = p;
      record(row, 3 * nrowtab + 5) = (double)curroot;
    }

    if (updateMu) {
      //MH move using Gamma(1e-3,1e3) prior
      double mu2 = fabs(R::rnorm(mu, sdMu));
      double omu = mu; mu = mu2;
      double l2 = fullLikC(tab, mu, sigma, minbralen, model, useRec);
      mu = omu;
      double mh = l2 - l + R::dgamma(mu2, 1e-3, 1e3, 1) - R::dgamma(mu, 1e-3, 1e3, 1);
      if (log(R::runif(0.0, 1.0)) < mh) { l = l2; mu = mu2; }
      if (tuning) sdMu = sdMu * exp(delta / (i + 1.0) * (std::min(1.0, exp(mh)) - acceptance_target));
    }

    if (updateSigma && sigma > 0) {
      //MH move using Gamma(1e-3,1e3) prior
      double sigma2 = fabs(R::rnorm(sigma, sdSigma));
      double osig = sigma; sigma = sigma2;
      double l2 = fullLikC(tab, mu, sigma, minbralen, model, useRec);
      sigma = osig;
      double mh = l2 - l + R::dgamma(sigma2, 1e-3, 1e3, 1) - R::dgamma(sigma, 1e-3, 1e3, 1);
      if (log(R::runif(0.0, 1.0)) < mh) { l = l2; sigma = sigma2; }
      if (tuning) sdSigma = sdSigma * exp(delta / (i + 1.0) * (std::min(1.0, exp(mh)) - acceptance_target));
    }

    if (mixedModel) {
      //Reversible-jump move
      double sigma2, qratio, pratio;
      if (sigma == 0) {
        sigma2 = R::rexp(1.0);
        qratio = sigma2;
        pratio = R::dgamma(sigma2, 1e-3, 1e3, 1);
      } else {
        sigma2 = 0;
        qratio = -sigma;
        pratio = -R::dgamma(sigma, 1e-3, 1e3, 1);
      }
      double osig = sigma; sigma = sigma2;
      double l2 = fullLikC(tab, mu, sigma, minbralen, model, useRec);
      sigma = osig;
      if (log(R::runif(0.0, 1.0)) < l2 - l + pratio + qratio) { l = l2; sigma = sigma2; }
    }

    if (updateAlpha) {
      //Gibbs move using inverse-gamma prior
      double su = coalAlphaSumSortedC(orderedleafdates, orderednodedates);
      alpha = 1.0 / R::rgamma(n + 0.001 - 1.0, 2000.0 / (su * 1000.0 + 2.0));
      p = useCoalPrior ? coalpriorC(orderedleafdates, orderednodedates, alpha) : 0.0;
      if (exact) {
        //rebuild the weight table and prefix cache for the new alpha
        buildCoalw();
        double pp = -log(alpha) * (ne - 1);
        int i1 = 1, i2 = 0; double prev = LV[0];
        Ssum[0] = pp; I1s[0] = 1; I2s[0] = 0; PrevS[0] = prev;
        for (int t = 1; t <= msteps; ++t) {
          COAL_STEP_INC(pp, t);
          COAL_CKPT(pp, t);
        }
        X = (msteps / COAL_CK) * COAL_CK;
      }
    }

    //MH to update internal dates
    for (int q = 0; q < nint; ++q) rn[q] = R::rnorm(0.0, sdDates);
    for (int j1 = n + 1; j1 <= nrowtab; ++j1) {
      double r = rn[j1 - n - 1];
      int row = j1 - 1;
      double old = tab(row, 2);
      double nw = old + r;
      bool fatherExists = !NumericVector::is_na(tab(row, 3));
      double fatherDate = fatherExists ? tab((int)tab(row, 3) - 1, 2) : 0.0;
      double minChildDate = R_PosInf;
      const std::vector<int>& ch = kids[j1];
      for (size_t k = 0; k < ch.size(); ++k)
        if (tab(ch[k] - 1, 2) < minChildDate) minChildDate = tab(ch[k] - 1, 2);
      bool boundary1 = (r < 0) && fatherExists && (nw < fatherDate);   //older than father
      bool boundary2 = (r > 0) && (nw > minChildDate);                 //younger than sons
      double mh;
      if (boundary1 || boundary2) {
        mh = R_NegInf;
      } else {
        double l2;
        if (mode == 0) {
          double oldlocal = slowNodeLocal(j1, old);
          double newlocal = slowNodeLocal(j1, nw);
          l2 = l - oldlocal + newlocal;
        } else {
          if (ch.size() < 16) {
            int termarr[16];
            int nterm = (int)ch.size();
            for (int k = 0; k < nterm; ++k) termarr[k] = ch[k];
            if (j1 != n + 1) termarr[nterm++] = j1;
            double oldlocal = localTermsArrC(tab, termarr, nterm, mu, sigma, minbralen, model, useRec);
            tab(row, 2) = nw;
            double newlocal = localTermsArrC(tab, termarr, nterm, mu, sigma, minbralen, model, useRec);
            l2 = l - oldlocal + newlocal;
          } else {
            terms = ch;
            if (j1 != n + 1) terms.push_back(j1);
            double oldlocal = localTermsVecC(tab, terms, mu, sigma, minbralen, model, useRec);
            tab(row, 2) = nw;
            double newlocal = localTermsVecC(tab, terms, mu, sigma, minbralen, model, useRec);
            l2 = l - oldlocal + newlocal;
          }
          //date stays set to nw; reverted below on rejection
        }
        changeinorderedvec(orderednodedates, old, nw);
        double p2;
        int pos_hi = 0;
        int c = 0;                       //checkpoint at or below pos_hi (exact path)
        int winEnd = 0;                  //P4: last walked (non-cached) step
        if (!useCoalPrior) {
          p2 = 0.0;
        } else if (mode == 0 || fast) {
          p2 = (mode == 0) ? coalpriorC(orderedleafdates, orderednodedates, alpha)
                           : p + coalDeltaC(orderedleafdates, orderednodedates, alpha, old, nw);
        } else {
          // Events strictly newer than hi are identical in the current and
          // proposed states, so the cached partial sum at that point is exact
          // for both; the first leaf (newest event) is consumed before step 1.
          double hi = nw > old ? nw : old;
          pos_hi = countGtL(hi) + countGtN(hi);
          if (pos_hi > 0) pos_hi--;
          // P2: advance the valid boundary X to the checkpoint at or below
          // pos_hi (that prefix is identical in both states), writing
          // checkpoints only.
          c = (pos_hi / COAL_CK) * COAL_CK;
          if (c > X) {
            double pp = Ssum[X];
            int i1 = I1s[X], i2 = I2s[X]; double prev = PrevS[X];
            for (int t = X + 1; t <= c; ++t) {
              COAL_STEP(pp);
              COAL_CKPT(pp, t);
            }
            X = c;
          }
          // Restart from that checkpoint and walk up to pos_hi (no writes).
          p2 = Ssum[c];
          int i1 = I1s[c], i2 = I2s[c]; double prev = PrevS[c];
          for (int t = c + 1; t <= pos_hi; ++t) COAL_STEP(p2);
          // P4: identical-tail boundary. Let lo = min(old,nw). The state-
          // dependent events are exactly those with date in [lo,hi] (the ones
          // involving old/nw); events with date < lo are identical in both
          // states and occupy positions countGe(lo)+1 .. end. Since
          // inc[t] = coalw[k_t]*(e_t - e_{t+1}) and k_t is a leaf-minus-node
          // count (old,nw are the same node, so counts match), inc[t] is
          // state-independent iff both e_t and e_{t+1} have date < lo, i.e.
          // iff t >= countGe(lo)+1. So the window must be walked through
          // t = countGe(lo) and only steps >= countGe(lo)+1 may use the cache.
          double lo = nw < old ? nw : old;
          winEnd = countGeL(lo) + countGeN(lo);   // = countGe(lo)
          if (winEnd > msteps) winEnd = msteps;
          // Window walk over the proposed arrays (stashing increments).
          for (int t = pos_hi + 1; t <= winEnd; ++t) COAL_STEP_WIN(p2, t);
          // Identical tail: subtract cached current-state increments in the
          // original left-to-right order (bit-identical accumulation).
          for (int t = winEnd + 1; t <= msteps; ++t) p2 -= inc[t];
        }
        mh = l2 - l + p2 - p;
        if (log(R::runif(0.0, 1.0)) < mh) {
          l = l2; p = p2;                 //accept: keep nw
          if (exact) {
            // P4: refresh the window increments for the new state; tail
            // increments are already correct. No checkpoint was written during
            // the proposal and every checkpoint <= c lies in the prefix
            // identical to both states, so X = c.
            for (int t = pos_hi + 1; t <= winEnd; ++t) inc[t] = winInc[t];
            X = c;
          }
        } else {
          changeinorderedvec(orderednodedates, nw, old);
          tab(row, 2) = old;              //reject: revert date
          // P4: proposal wrote no checkpoints and left inc untouched.
        }
      }
      if (tuning) {
        double ii = (double)(i - 1) * (double)nint + (double)j1;
        sdDates = sdDates * exp(delta / (ii + 1.0) * (std::min(1.0, exp(mh)) - acceptance_target));
      }
    }

    //MH to update missing leaf dates
    for (int q = 0; q < misDates.size(); ++q) {
      int j1 = misDates[q];               //1-based leaf row
      int row = j1 - 1;
      double old = tab(row, 2);
      double nw = R::rnorm(old, (rangedate(j1 - 1, 1) - rangedate(j1 - 1, 0)) * 0.05);
      if (nw - old < 0 && (!NumericVector::is_na(tab(row, 3)) && nw < tab((int)tab(row, 3) - 1, 2))) continue; //older than father
      if (nw > rangedate(j1 - 1, 1) || nw < rangedate(j1 - 1, 0)) continue; //stay within prior range
      double l2;
      if (mode == 0) {
        l2 = l - slowLeafLocal(j1, old) + slowLeafLocal(j1, nw);
      } else {
        terms.assign(1, j1);
        double oldlocal = localTermsVecC(tab, terms, mu, sigma, minbralen, model, useRec);
        tab(row, 2) = nw;
        double newlocal = localTermsVecC(tab, terms, mu, sigma, minbralen, model, useRec);
        l2 = l - oldlocal + newlocal;
      }
      changeinorderedvec(orderedleafdates, old, nw);
      double p2;
      int pos_hi = 0;
      int c = 0;                     //checkpoint at or below pos_hi (exact path)
      int winEnd = 0;                //P4: last walked (non-cached) step
      if (!useCoalPrior) p2 = 0.0;
      else if (mode == 0) p2 = coalpriorC(orderedleafdates, orderednodedates, alpha);
      else if (fast) p2 = p + coalDeltaC(orderedleafdates, orderednodedates, alpha, old, nw);
      else { //exact: bit-identical prefix-cache tail walk (same as node dates)
        double hi = nw > old ? nw : old;
        pos_hi = countGtL(hi) + countGtN(hi);
        if (pos_hi > 0) pos_hi--;
        // P2: advance the valid boundary X to the checkpoint at or below
        // pos_hi (identical in both states), writing checkpoints only.
        c = (pos_hi / COAL_CK) * COAL_CK;
        if (c > X) {
          double pp = Ssum[X];
          int i1 = I1s[X], i2 = I2s[X]; double prev = PrevS[X];
          for (int t = X + 1; t <= c; ++t) {
            COAL_STEP(pp);
            COAL_CKPT(pp, t);
          }
          X = c;
        }
        // Restart from that checkpoint and walk up to pos_hi (no writes).
        p2 = Ssum[c];
        int i1 = I1s[c], i2 = I2s[c]; double prev = PrevS[c];
        for (int t = c + 1; t <= pos_hi; ++t) COAL_STEP(p2);
        // P4: identical-tail boundary (same reasoning as the node branch).
        double lo = nw < old ? nw : old;
        winEnd = countGeL(lo) + countGeN(lo);   // = countGe(lo)
        if (winEnd > msteps) winEnd = msteps;
        for (int t = pos_hi + 1; t <= winEnd; ++t) COAL_STEP_WIN(p2, t);
        for (int t = winEnd + 1; t <= msteps; ++t) p2 -= inc[t];
      }
      if (log(R::runif(0.0, 1.0)) < l2 - l + p2 - p) {
        l = l2; p = p2;
        if (exact) {
          for (int t = pos_hi + 1; t <= winEnd; ++t) inc[t] = winInc[t];
          X = c;
        }
      } else {
        changeinorderedvec(orderedleafdates, nw, old);
        tab(row, 2) = old;
        // P4: proposal wrote no checkpoints and left inc untouched.
      }
    }

    if (updateRootAll || updateRootBranch) {
      //Move root on current branch
      std::vector<int> sidesv = kids[n + 1];
      if (mode == 0) {
        //original R: root=which(is.na(tab[,4])); sides=which(tab[,4]==root)
        int root = -1;
        for (int q = 0; q < nrowtab; ++q) if (NumericVector::is_na(tab(q, 3))) { root = q + 1; break; }
        sidesv.clear();
        for (int q = 0; q < nrowtab; ++q) if (!NumericVector::is_na(tab(q, 3)) && (int)tab(q, 3) == root) sidesv.push_back(q + 1);
      }
      double o0 = tab(sidesv[0] - 1, 1), o1 = tab(sidesv[1] - 1, 1);
      double oldlocal = fast ? localTermsVecC(tab, sidesv, mu, sigma, minbralen, model, useRec) : 0.0;
      double r = R::runif(0.0, 1.0);
      tab(sidesv[0] - 1, 1) = (o0 + o1) * r;
      tab(sidesv[1] - 1, 1) = (o0 + o1) * (1 - r);
      double l2 = (mode == 0 || exact) ? fullLikC(tab, mu, sigma, minbralen, model, useRec)
        : l - oldlocal + localTermsVecC(tab, sidesv, mu, sigma, minbralen, model, useRec);
      if (log(R::runif(0.0, 1.0)) < l2 - l) l = l2;
      else { tab(sidesv[0] - 1, 1) = o0; tab(sidesv[1] - 1, 1) = o1; }
    }

    if (updateRootAll) {
      //Move root branch
      std::vector<int> sidesv = kids[n + 1];
      if (mode == 0) {
        int root = -1;
        for (int q = 0; q < nrowtab; ++q) if (NumericVector::is_na(tab(q, 3))) { root = q + 1; break; }
        sidesv.clear();
        for (int q = 0; q < nrowtab; ++q) if (!NumericVector::is_na(tab(q, 3)) && (int)tab(q, 3) == root) sidesv.push_back(q + 1);
      }
      int left, right;
      if (tab(sidesv[0] - 1, 2) < tab(sidesv[1] - 1, 2)) { left = sidesv[0]; right = sidesv[1]; }
      else { left = sidesv[1]; right = sidesv[0]; }
      if (left > n) {
        std::vector<int> abv = (mode == 0) ? std::vector<int>() : kids[left];
        if (mode == 0) {
          for (int q = 0; q < nrowtab; ++q)
            if (!NumericVector::is_na(tab(q, 3)) && (int)tab(q, 3) == left) abv.push_back(q + 1);
        }
        int a = (R::runif(0.0, 1.0) < 0.5) ? abv[0] : abv[1];
        double olda2 = tab(a - 1, 1), oldleft2 = tab(left - 1, 1), oldright2 = tab(right - 1, 1);
        int olda4 = (int)tab(a - 1, 3), oldright4 = (int)tab(right - 1, 3);
        double oldleft5 = useRec ? tab(left - 1, 4) : 0.0;
        if (mode == 0) {
          NumericMatrix oldtab = clone(tab);
          tab(a - 1, 3) = n + 1;
          tab(right - 1, 3) = left;
          double rr = R::runif(0.0, 1.0);
          tab(a - 1, 1) = oldtab(a - 1, 1) * rr;
          tab(left - 1, 1) = oldtab(a - 1, 1) * (1 - rr);
          tab(right - 1, 1) = oldtab(right - 1, 1) + oldtab(left - 1, 1);
          if (useRec) tab(left - 1, 4) = tab(a - 1, 4);
          double l2 = fullLikC(tab, mu, sigma, minbralen, model, useRec);
          double hast = log((oldtab(a - 1, 1) == 0 ? 1.0 : oldtab(a - 1, 1)) /
                            (tab(right - 1, 1) == 0 ? 1.0 : tab(right - 1, 1)));
          if (log(R::runif(0.0, 1.0)) < l2 - l + hast) {
            l = l2;
            rebuildKids();
            std::vector<int> rc = kids[n + 1];
            curroot = findCurrootC(edge, NumericVector::create((double)rc[0], (double)rc[1]), n + 1);
          } else tab = oldtab;
        } else {
          std::vector<int> aff; aff.push_back(a); aff.push_back(left); aff.push_back(right);
          std::sort(aff.begin(), aff.end());
          double oldlocal = fast ? localTermsVecC(tab, aff, mu, sigma, minbralen, model, useRec) : 0.0;
          tab(a - 1, 3) = n + 1;
          tab(right - 1, 3) = left;
          double rr = R::runif(0.0, 1.0);
          tab(a - 1, 1) = olda2 * rr;
          tab(left - 1, 1) = olda2 * (1 - rr);
          tab(right - 1, 1) = oldright2 + oldleft2;
          if (useRec) tab(left - 1, 4) = tab(a - 1, 4);
          double l2 = exact ? fullLikC(tab, mu, sigma, minbralen, model, useRec)
            : l - oldlocal + localTermsVecC(tab, aff, mu, sigma, minbralen, model, useRec);
          double hast = log((olda2 == 0 ? 1.0 : olda2) / (tab(right - 1, 1) == 0 ? 1.0 : tab(right - 1, 1)));
          if (log(R::runif(0.0, 1.0)) < l2 - l + hast) {
            l = l2;
            rebuildKids();
            std::vector<int> rc = kids[n + 1];
            curroot = findCurrootC(edge, NumericVector::create((double)rc[0], (double)rc[1]), n + 1);
          } else {
            tab(a - 1, 3) = olda4; tab(right - 1, 3) = oldright4;
            tab(a - 1, 1) = olda2; tab(left - 1, 1) = oldleft2; tab(right - 1, 1) = oldright2;
            if (useRec) tab(left - 1, 4) = oldleft5;
          }
        }
      }
    }
  } //End of MCMC loop

  return List::create(_["tab"] = tab, _["l"] = l, _["p"] = p,
                      _["mu"] = mu, _["sigma"] = sigma, _["alpha"] = alpha,
                      _["sdMu"] = sdMu, _["sdSigma"] = sdSigma, _["sdDates"] = sdDates,
                      _["curroot"] = curroot, _["record"] = record);
}
