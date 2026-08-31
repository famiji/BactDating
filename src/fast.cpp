#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Prototypes of functions defined in probs.cpp (same shared library).
double coalpriorC(const NumericVector& leaves, const NumericVector& nodes, double alpha);
void changeinorderedvec(NumericVector vec, double old, double n);

// ---------------------------------------------------------------------------
// Per-node likelihood term, replicating the body of the per-row accumulation in
// src/probs.cpp exactly (same expressions, same evaluation order) so that the
// local likelihood computed here is bit-identical to the full likelihood terms.
//
//   i        : 0-based row index in tab (node = i+1)
//   mu,sigma : clock parameters
//   minbralen: minimum branch length threshold (used by some models)
//   model    : 1=poisson 2=negbin 3=strictgamma 4=relaxedgamma 5=arc 6=carc
//   useRec   : whether tab has a 5th (unrecombined) column
// ---------------------------------------------------------------------------
static inline double termC(const NumericMatrix& tab, int i, double mu, double sigma,
                           double minbralen, int model, bool useRec) {
  double unrec = 1.0;
  double l;
  switch (model) {
    case 1: { // poisson (likelihoodPoissonC)
      if (useRec) unrec = tab(i, 4);
      return R::dpois(round(unrec * tab(i, 1)),
                      unrec * mu * (tab(i, 2) - tab((int)tab(i, 3) - 1, 2)), 1);
    }
    case 2: { // negbin (likelihoodNegbinC)
      double k = mu * mu / sigma / sigma;
      double theta = sigma * sigma / mu;
      if (useRec) unrec = tab(i, 4);
      l = unrec * (tab(i, 2) - tab((int)tab(i, 3) - 1, 2));
      return R::dnbinom(round(unrec * tab(i, 1)), k, 1 / (1.0 + theta * l), 1);
    }
    case 3: { // strictgamma (likelihoodGammaC)
      if (useRec) unrec = tab(i, 4);
      if (unrec * tab(i, 1) <= minbralen)
        return R::pgamma(minbralen, unrec * mu * (tab(i, 2) - tab((int)tab(i, 3) - 1, 2)), 1, 1, 1);
      else
        return R::dgamma(unrec * tab(i, 1), unrec * mu * (tab(i, 2) - tab((int)tab(i, 3) - 1, 2)), 1, 1);
    }
    case 4: { // relaxedgamma (likelihoodRelaxedgammaC)
      double ratevar = sigma * sigma;
      l = tab(i, 2) - tab((int)tab(i, 3) - 1, 2);
      if (useRec) { unrec = tab(i, 4); l = l * unrec; }
      if (unrec * tab(i, 1) <= minbralen)
        return R::pgamma(minbralen, l * mu * mu / (mu + l * ratevar), 1.0 + l * ratevar / mu, 1, 1);
      else
        return R::dgamma(unrec * tab(i, 1), l * mu * mu / (mu + l * ratevar), 1.0 + l * ratevar / mu, 1);
    }
    case 5: { // arc (likelihoodArcC)
      l = tab(i, 2) - tab((int)tab(i, 3) - 1, 2);
      if (useRec) { unrec = tab(i, 4); l = l * unrec; }
      return R::dnbinom(round(unrec * tab(i, 1)), mu * l / sigma, 1.0 - sigma / (1.0 + sigma), 1);
    }
    case 6: { // carc (likelihoodCarcC)
      if (useRec) unrec = tab(i, 4);
      l = unrec * (tab(i, 2) - tab((int)tab(i, 3) - 1, 2));
      if (unrec * tab(i, 1) <= minbralen)
        return R::pgamma(minbralen, mu * l / (1.0 + sigma), 1.0 + sigma, 1, 1);
      else
        return R::dgamma(unrec * tab(i, 1), mu * l / (1.0 + sigma), 1.0 + sigma, 1);
    }
  }
  return 0.0; // unreachable
}

// ---------------------------------------------------------------------------
// Local log-likelihood of the set of terms that change when node j's date is
// modified: the term of j itself (unless j is the global root n+1, which carries
// no branch above it and is skipped by the full likelihood) plus the terms of all
// children of j (their branch length depends on j's date). This exactly equals
// likelihood(mintab, ...) from main.R for the corresponding node.
//   j1 : 1-based node index
// ---------------------------------------------------------------------------
static double localLik(const NumericMatrix& tab, int j1, int n,
                       const std::vector<std::vector<int> >& children,
                       double mu, double sigma, double minbralen,
                       int model, bool useRec) {
  double s = 0.0;
  const std::vector<int>& ch = children[j1];
  // Children first (matches mintab row order: children rows precede j's row).
  for (size_t k = 0; k < ch.size(); ++k)
    s += termC(tab, ch[k] - 1, mu, sigma, minbralen, model, useRec);
  // j's own term, unless j is the global root (1-based n+1) which is skipped.
  if (j1 != n + 1)
    s += termC(tab, j1 - 1, mu, sigma, minbralen, model, useRec);
  return s;
}

// ---------------------------------------------------------------------------
// Sum of the likelihood terms of an explicit set of nodes (1-based indices),
// in the given order. Used to localize the root-move likelihood updates.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double localTermsC(const NumericMatrix& tab, IntegerVector nodes1, double mu, double sigma,
                   double minbralen, int model, bool useRec) {
  double s = 0.0;
  for (int i = 0; i < nodes1.size(); ++i)
    s += termC(tab, nodes1[i] - 1, mu, sigma, minbralen, model, useRec);
  return s;
}

// Same as localTermsC but without the SEXP-wrapping cost: used from the C++
// MCMC loop where the index set is already a plain integer vector.
static inline double localTermsVec(const NumericMatrix& tab, const std::vector<int>& nodes1,
                                   double mu, double sigma, double minbralen,
                                   int model, bool useRec) {
  double s = 0.0;
  for (size_t i = 0; i < nodes1.size(); ++i)
    s += termC(tab, nodes1[i] - 1, mu, sigma, minbralen, model, useRec);
  return s;
}

// Non-static alias for cross-file use (mcmcLoopC.cpp).
double localTermsVecC(const NumericMatrix& tab, const std::vector<int>& nodes1,
                      double mu, double sigma, double minbralen,
                      int model, bool useRec) {
  return localTermsVec(tab, nodes1, mu, sigma, minbralen, model, useRec);
}

// Same, taking a plain C array + length: lets the MCMC loop pass a fixed-size
// stack buffer instead of copying a std::vector per node update.
double localTermsArrC(const NumericMatrix& tab, const int* nodes1, int n1,
                      double mu, double sigma, double minbralen,
                      int model, bool useRec) {
  double s = 0.0;
  for (int i = 0; i < n1; ++i)
    s += termC(tab, nodes1[i] - 1, mu, sigma, minbralen, model, useRec);
  return s;
}

// ---------------------------------------------------------------------------
// Full-table log-likelihood via termC, replicating the row order of the
// likelihoodXC functions in probs.cpp (all rows in order, root row n+1
// skipped). Bit-identical to them.
//
// The per-row terms are mutually independent (termC only reads tab and calls
// pure Rmath density functions - no RNG, no global state), so on OpenMP-capable
// builds (Linux toolchains; enabled via src/Makevars) the terms are evaluated
// in parallel into a buffer and then summed strictly in row order, which
// reproduces the original sequential accumulation bit for bit. On builds
// without OpenMP (e.g. macOS) the pragmas are ignored and this stays serial.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double fullLikC(const NumericMatrix& tab, double mu, double sigma,
                double minbralen, int model, bool useRec) {
  int n = (tab.nrow() + 1) / 2;
  int nr = tab.nrow();
#ifdef _OPENMP
  // Parallel path: only when OpenMP is active, more than one thread is
  // available, and the table is large enough to amortize buffer allocation
  // and team startup. Each worker writes a disjoint, contiguous index range
  // (schedule(static)); the serial row-order sum then reads every element
  // once, reproducing the original accumulation bit for bit (verified
  // identical to the serial path on the same inputs).
  int maxthreads = omp_get_max_threads();
  if (nr > 512 && maxthreads > 1) {
    // P3: persistent shared buffer (function-static). The MCMC loop is serial
    // outside these OpenMP regions, so only one fullLikC is ever active and
    // no locking is needed; all nr entries are written before any is read,
    // so stale contents are harmless. This avoids a fresh ~nr*8-byte
    // allocation (and zero-fill) on each of the thousands of calls per run.
    static std::vector<double> buf;
    if ((int)buf.size() < nr) buf.resize(nr);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nr; ++i)
      buf[i] = (i == n) ? 0.0 : termC(tab, i, mu, sigma, minbralen, model, useRec);
    double s = 0.0;
    for (int i = 0; i < nr; ++i) {
      if (i == n) continue;
      s += buf[i];
    }
    return s;
  }
#endif
  // Serial path: direct accumulation, no buffer - bit-identical to the
  // original likelihoodXC implementations in probs.cpp.
  double s = 0.0;
  for (int i = 0; i < nr; ++i) {
    if (i == n) continue;
    s += termC(tab, i, mu, sigma, minbralen, model, useRec);
  }
  return s;
}

// Number of elements >= val in a vector sorted in decreasing order.
static int countGe(const NumericVector& v, double val) {
  int lo = 0, hi = v.size();
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (v[mid] >= val) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// Number of elements > val in a vector sorted in decreasing order.
static int countGt(const NumericVector& v, double val) {
  int lo = 0, hi = v.size();
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    if (v[mid] > val) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// ---------------------------------------------------------------------------
// Fast equivalent of main.R's curroot search:
//   for (j in 1:nrow(tree$edge)) if (setequal(rootchildren,tree$edge[j,])) ...
//   if (is.na(curroot)) curroot=min(which(tree$edge[,1]==(n+1)))
// Returns the same 1-based edge index, without the per-edge setequal/sort
// overhead. rootchildren holds the two children of the global root n+1.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
int findCurrootC(const NumericMatrix& edge, NumericVector rootchildren, int rootnode) {
  double r1 = rootchildren[0], r2 = rootchildren[1];
  for (int j = 0; j < edge.nrow(); ++j) {
    double e1 = edge(j, 0), e2 = edge(j, 1);
    if ((e1 == r1 && e2 == r2) || (e1 == r2 && e2 == r1)) return j + 1;
  }
  for (int j = 0; j < edge.nrow(); ++j)
    if (edge(j, 0) == (double)rootnode) return j + 1;
  return NA_INTEGER;
}

// ---------------------------------------------------------------------------
// Sort-free equivalent of the sum-of-sorted-differences statistic used by the
// alpha Gibbs move in main.R:
//   su = sum(k*(k-1)*difs) over the merged decreasing sequence of leaf and
//   node dates, where k is the running leaf-minus-node count.
// Merging two already-sorted vectors reproduces the R sort-based computation
// with identical arithmetic (tie order is immaterial: the gap between tied
// values is zero and the count after the tied block is order-independent).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double coalAlphaSumC(const NumericVector& leaves, const NumericVector& nodes) {
  std::vector<double> lv(leaves.begin(), leaves.end());
  std::vector<double> nd(nodes.begin(), nodes.end());
  std::sort(lv.begin(), lv.end(), std::greater<double>());
  std::sort(nd.begin(), nd.end(), std::greater<double>());
  int nl = (int)lv.size(), nn = (int)nd.size();
  int i1 = 0, i2 = 0, k = 0;
  double su = 0.0, prev = 0.0;
  bool first = true;
  while (i1 < nl || i2 < nn) {
    double e; bool isLeaf;
    if (i2 >= nn || (i1 < nl && lv[i1] >= nd[i2])) { e = lv[i1]; isLeaf = true; }
    else { e = nd[i2]; isLeaf = false; }
    if (!first) su += (double)k * (double)(k - 1) * (prev - e);
    k += isLeaf ? 1 : -1;
    prev = e;
    first = false;
    if (isLeaf) i1++; else i2++;
  }
  return su;
}

// Same statistic without the defensive copies and sorts: requires both
// inputs to already be sorted in decreasing order (as maintained by
// changeinorderedvec / the sort() calls in main.R). Sorting an already
// sorted vector is an identity operation, so this is bit-identical to
// coalAlphaSumC on such inputs. Used by the C++ MCMC loop's Gibbs move.
double coalAlphaSumSortedC(const NumericVector& leaves, const NumericVector& nodes) {
  const double* lv = leaves.begin();
  const double* nd = nodes.begin();
  int nl = (int)leaves.size(), nn = (int)nodes.size();
  int i1 = 0, i2 = 0, k = 0;
  long double su = 0.0;
  double prev = 0.0;
  bool first = true;
  while (i1 < nl || i2 < nn) {
    double e; bool isLeaf;
    if (i2 >= nn || (i1 < nl && lv[i1] >= nd[i2])) { e = lv[i1]; isLeaf = true; }
    else { e = nd[i2]; isLeaf = false; }
    if (!first) su += (double)k * (double)(k - 1) * (prev - e);
    k += isLeaf ? 1 : -1;
    prev = e;
    first = false;
    if (isLeaf) i1++; else i2++;
  }
  return (double)su;
}

// ---------------------------------------------------------------------------
// Incremental coalescent prior update.
//
// coalpriorC is O(n) and was called once per internal node per iteration.
// Moving a single node date from oldval to newval only changes the lineage
// count between the two dates, so the prior difference can be computed over
// that interval alone:
//
//   With P = -log(alpha)(n-1) - S, S = sum k(k-1)/(2alpha) * dt, letting
//   I = integral of k_new(t) over (lo, hi) in the UPDATED state
//   (nodes already contains newval instead of oldval):
//     newval < oldval : dP = ((hi - lo) - I) / alpha
//     newval > oldval : dP =  I / alpha
//
// Cost is O(log n + events crossed), typically O(1) for tuned MH proposals.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double coalDeltaC(const NumericVector& leaves, const NumericVector& nodes, double alpha,
                  double oldval, double newval) {
  if (newval == oldval) return 0.0;
  double lo = std::min(oldval, newval), hi = std::max(oldval, newval);
  int nl = leaves.size(), nn = nodes.size();
  int i1 = countGe(leaves, hi), i2 = countGe(nodes, hi);
  int k = i1 - i2;                 // lineages just below hi in the new state
  double I = 0.0, prev = hi;
  while (true) {
    double nextL = (i1 < nl) ? leaves[i1] : R_NegInf;
    double nextN = (i2 < nn) ? nodes[i2] : R_NegInf;
    double e = std::max(nextL, nextN);
    if (!(e > lo)) break;          // events exactly at lo lie outside (lo,hi)
    I += (double)k * (prev - e);
    if (nextL >= nextN) { k++; i1++; } else { k--; i2++; }
    prev = e;
  }
  I += (double)k * (prev - lo);
  if (newval < oldval) return ((hi - lo) - I) / alpha;
  return I / alpha;
}

// ---------------------------------------------------------------------------
// Fast replacement for the internal-node-date MH loop of main.R (lines ~199-224).
// Mutates tab (dates column) and orderednodedates in place; returns updated
// l, p and sdDates. RNG consumption matches the original R loop exactly:
// the per-node proposals rn are drawn in R beforehand, and one uniform is drawn
// (via R::runif) per node that passes the boundary checks, in the same order.
//
// When incrementalPrior is false the coalescent prior is still evaluated
// bit-identically to coalpriorC, but via a prefix-sum cache: coalpriorC is a
// single sequential accumulation from the newest event to the oldest, and a
// proposal only relocates one node event, so every event strictly newer than
// max(old,new) yields the same partial sum as a full recomputation. The sweep
// caches those partial sums and only walks the tail (events at or older than
// max(old,new)) for each proposal - the same operations, in the same order,
// as coalpriorC on the proposed arrays, hence bit-identical results at a
// fraction of the cost.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List nodeDatesUpdateC(NumericMatrix tab, int n,
                      NumericVector orderedleafdates, NumericVector orderednodedates,
                      NumericVector rn, double mu, double sigma, double l, double p,
                      double alpha, double minbralen, double sdDates, int outerit,
                      bool tuning, bool useCoalPrior, int model, bool useRec,
                      bool incrementalPrior) {
  int nrowtab = tab.nrow();
  int nint = nrowtab - n;               // number of internal nodes (incl root)
  double acceptance_target = 0.234;
  double delta = 1.0 / acceptance_target / (1.0 - acceptance_target);

  // Build children structure from the father column (fixed during this sweep).
  std::vector<std::vector<int> > children(nrowtab + 1);
  for (int x = 1; x <= nrowtab; ++x) {
    double f = tab(x - 1, 3);
    if (!NumericVector::is_na(f))
      children[(int)f].push_back(x);
  }

  // Prefix-sum cache for the exact (bit-identical) prior evaluation.
  // Walk state after t merge steps: leaves/nodes consumed (I1s/I2s), last
  // event value (PrevS) and accumulated p (Ssum), matching coalpriorC's
  // arithmetic step for step. Valid for the current state up to step X.
  const double* LV = orderedleafdates.begin();
  const double* ND = orderednodedates.begin();
  int ne = orderedleafdates.size();
  int msteps = 2 * ne - 2;             // merge steps after the initial leaf
  std::vector<double> Ssum, PrevS;
  std::vector<int> I1s, I2s;
  int X = 0;
  bool exactPrior = useCoalPrior && !incrementalPrior;
  if (exactPrior) {
    Ssum.resize(msteps + 1); PrevS.resize(msteps + 1);
    I1s.resize(msteps + 1); I2s.resize(msteps + 1);
    double pp = -log(alpha) * (ne - 1);
    int i1 = 1, i2 = 0; double prev = LV[0];
    Ssum[0] = pp; I1s[0] = 1; I2s[0] = 0; PrevS[0] = prev;
    for (int t = 1; t <= msteps; ++t) {
      int k = i1 - i2;
      if (i1 < ne && LV[i1] > ND[i2]) {
        pp = pp - (k * (k - 1.0) / (2.0 * alpha) * (prev - LV[i1]));
        prev = LV[i1++];
      } else {
        pp = pp - (k * (k - 1.0) / (2.0 * alpha) * (prev - ND[i2]));
        prev = ND[i2++];
      }
      Ssum[t] = pp; I1s[t] = i1; I2s[t] = i2; PrevS[t] = prev;
    }
    X = msteps;
  }

  for (int j1 = n + 1; j1 <= nrowtab; ++j1) {
    double r = rn[j1 - n - 1];          // rn is 0-based in C++, rn[j-n] in R
    int row = j1 - 1;
    double old = tab(row, 2);
    double nw = old + r;

    bool fatherExists = !NumericVector::is_na(tab(row, 3));
    double fatherDate = fatherExists ? tab((int)tab(row, 3) - 1, 2) : 0.0;
    double minChildDate = R_PosInf;
    const std::vector<int>& ch = children[j1];
    for (size_t k = 0; k < ch.size(); ++k)
      if (tab(ch[k] - 1, 2) < minChildDate) minChildDate = tab(ch[k] - 1, 2);

    double mh;
    bool boundary1 = (r < 0) && fatherExists && (nw < fatherDate);   // older than father
    bool boundary2 = (r > 0) && (nw > minChildDate);                 // younger than sons
    if (boundary1) {
      mh = R_NegInf;
    } else if (boundary2) {
      mh = R_NegInf;
    } else {
      double oldlocal = localLik(tab, j1, n, children, mu, sigma, minbralen, model, useRec);
      tab(row, 2) = nw;
      double newlocal = localLik(tab, j1, n, children, mu, sigma, minbralen, model, useRec);
      double l2 = l - oldlocal + newlocal;
      changeinorderedvec(orderednodedates, old, nw);
      double p2;
      int pos_hi = 0;                  // merge steps consuming events > max(old,nw)
      if (!useCoalPrior) {
        p2 = 0.0;
      } else if (incrementalPrior) {
        p2 = p + coalDeltaC(orderedleafdates, orderednodedates, alpha, old, nw);
      } else {
        // Events strictly newer than hi are identical in the current and
        // proposed states, so the cached partial sum at that point is exact
        // for both; the first leaf (newest event) is consumed before step 1.
        double hi = nw > old ? nw : old;
        pos_hi = countGt(orderedleafdates, hi) + countGt(orderednodedates, hi);
        if (pos_hi > 0) pos_hi--;
        if (pos_hi > X) {
          // Extend the cache over events > hi (identical in both states).
          double pp = Ssum[X];
          int i1 = I1s[X], i2 = I2s[X]; double prev = PrevS[X];
          for (int t = X + 1; t <= pos_hi; ++t) {
            int k = i1 - i2;
            if (i1 < ne && LV[i1] > ND[i2]) {
              pp = pp - (k * (k - 1.0) / (2.0 * alpha) * (prev - LV[i1]));
              prev = LV[i1++];
            } else {
              pp = pp - (k * (k - 1.0) / (2.0 * alpha) * (prev - ND[i2]));
              prev = ND[i2++];
            }
            Ssum[t] = pp; I1s[t] = i1; I2s[t] = i2; PrevS[t] = prev;
          }
          X = pos_hi;
        }
        // Walk the tail over the proposed arrays: the same operations, in
        // the same order, as coalpriorC from that point on.
        p2 = Ssum[pos_hi];
        int i1 = I1s[pos_hi], i2 = I2s[pos_hi]; double prev = PrevS[pos_hi];
        for (int t = pos_hi + 1; t <= msteps; ++t) {
          int k = i1 - i2;
          if (i1 < ne && LV[i1] > ND[i2]) {
            p2 = p2 - (k * (k - 1.0) / (2.0 * alpha) * (prev - LV[i1]));
            prev = LV[i1++];
          } else {
            p2 = p2 - (k * (k - 1.0) / (2.0 * alpha) * (prev - ND[i2]));
            prev = ND[i2++];
          }
          Ssum[t] = p2; I1s[t] = i1; I2s[t] = i2; PrevS[t] = prev;
        }
      }
      mh = l2 - l + p2 - p;
      if (log(R::runif(0.0, 1.0)) < mh) {
        l = l2; p = p2;                 // accept: keep tab(row,2)=nw
        if (exactPrior) X = msteps;     // tail sums now describe the new state
      } else {
        changeinorderedvec(orderednodedates, nw, old);
        tab(row, 2) = old;              // reject: revert date
        if (exactPrior) X = pos_hi;     // cache valid only up to the prefix
      }
    }

    if (tuning) {
      double ii = (double)(outerit - 1) * (double)nint + (double)j1;
      sdDates = sdDates * exp(delta / (ii + 1.0) * (std::min(1.0, exp(mh)) - acceptance_target));
    }
  }

  return List::create(_["l"] = l, _["p"] = p, _["sdDates"] = sdDates);
}
