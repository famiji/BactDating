#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// C++ reimplementation of leafDates() from R/roototip.R.
//
// The original R version walks each leaf up to the root, using which() to
// find the row of the edge whose child (=edge[,2]) is the current node.
// That which() is an O(n) scan per step, making the whole function O(n^2
// log n) - it dominated >90% of bactdate() runtime at n=1000 (Rprof).
//
// This version pre-builds a direct-mapped lookup table (child node number
// -> edge row index) so each step is O(1), bringing the total to O(n).
//
// The floating-point accumulation order is IDENTICAL to the R version:
//   dates[i] = rootdate
//   dates[i] = dates[i] + edge.length[r]   (leaf -> parent)
//   dates[i] = dates[i] + edge.length[r]   (parent -> grandparent)
//   ...
// so the result is bit-identical to leafDates() on the same tree, which
// keeps the initial rate (mu0) unchanged and thus preserves the MCMC
// trajectory bit-for-bit (roottotip runs after set.seed, but leafDates
// itself consumes no RNG - only the permutation test does).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector leafDatesC(const IntegerMatrix& edge, const NumericVector& elen,
                         int ntip, double rootdate) {
  int ne = edge.nrow();
  int ntotal = 2 * ntip - 1;            // max node number
  std::vector<int> child2row(ntotal + 1, -1);
  for (int r = 0; r < ne; ++r)
    child2row[edge(r, 1)] = r;          // 1-based child -> 0-based row
  NumericVector dates(ntip);
  for (int i = 0; i < ntip; ++i) dates[i] = rootdate;
  for (int i = 0; i < ntip; ++i) {
    int w = i + 1;                       // 1-based leaf
    while (true) {
      int r = child2row[w];
      if (r < 0) break;                  // reached root
      dates[i] += elen[r];
      w = edge(r, 0);                    // parent node
    }
  }
  return dates;
}

// [[Rcpp::export]]
double coalpriorC(const NumericVector& leaves, const NumericVector& nodes, double alpha) {
  int n=leaves.length();
  //NumericVector nodes = clone(intnodes);
  //std::sort(nodes.begin(), nodes.end(), std::greater<double>());
  double p = -log(alpha) * (n - 1);
  int i1=1,i2=0,k=1;
  double prev=leaves[0];
  for (int i=1;i<(n+n-1);i++) {
    if (i1<n&&leaves[i1]>nodes[i2]) {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-leaves[i1]));
      prev=leaves[i1++];
      k++;
    } else {
      p = p - (k * (k - 1.0) / (2.0 * alpha) * (prev-nodes[i2]));
      prev=nodes[i2++];
      k--;
    }
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodGammaC(NumericMatrix tab, double mu) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  double minbralen=tab(n,1);
  //if (NumericVector::is_na(minbralen)) minbralen=0.1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    if (unrec*tab(i,1)<=minbralen)
      p+=R::pgamma(minbralen,unrec*mu*(tab(i,2)-tab(tab(i,3)-1,2)),1,1,1);
    else
      p+=R::dgamma(unrec*tab(i,1),unrec*mu*(tab(i,2)-tab(tab(i,3)-1,2)),1,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodRelaxedgammaC(NumericMatrix tab, double mu, double sigma) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double l=0;
  double unrec=1;
  double minbralen=tab(n,1);
  //if (NumericVector::is_na(minbralen)) minbralen=0.1;
  double ratevar=sigma*sigma;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    l=tab(i,2)-tab(tab(i,3)-1,2);
    if (tab.ncol()==5) {unrec=tab(i,4);l=l*unrec;}
    if (unrec*tab(i,1)<=minbralen)
      p+=R::pgamma(minbralen,l*mu*mu/(mu+l*ratevar),1.0+l*ratevar/mu,1,1);
    else
      p+=R::dgamma(unrec*tab(i,1),l*mu*mu/(mu+l*ratevar),1.0+l*ratevar/mu,1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodPoissonC(NumericMatrix tab, double mu) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    p+=R::dpois(round(unrec*tab(i,1)),unrec*mu*(tab(i,2)-tab(tab(i,3)-1,2)),1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodNegbinC(NumericMatrix tab, double mu, double sigma) {
  double k=mu*mu/sigma/sigma;
  double theta=sigma*sigma/mu;
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  double l;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    l=unrec*(tab(i,2)-tab(tab(i,3)-1,2));
    p+=R::dnbinom(round(unrec*tab(i,1)),k,1/(1.0+theta*l),1);
  }
  return(p);
}

// [[Rcpp::export]]
void changeinorderedvec(NumericVector vec,double old,double n) {
  //NumericVector res = clone(vec);
  //Binary search for the first occurrence of old (the vector is sorted in
  //decreasing order): same index as the original linear scan, but O(log n).
  int lo=0, hi=vec.length();
  while (lo<hi) {
    int mid=(lo+hi)/2;
    if (vec[mid]>old) lo=mid+1; else hi=mid;
  }
  int i=lo;
  vec(i)=n;
  while (1) {
    if (i>0               &&vec(i-1)<vec(i)) {vec(i)=vec(i-1);vec(i-1)=n;i--;continue;}
    if (i<(vec.length()-1)&&vec(i+1)>vec(i)) {vec(i)=vec(i+1);vec(i+1)=n;i++;continue;}
    break;
  }
  //  return(res);
}

// [[Rcpp::export]]
double likelihoodArcC(NumericMatrix tab, double mu, double sigma) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double l=0;
  double unrec=1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    l=tab(i,2)-tab(tab(i,3)-1,2);
    if (tab.ncol()==5) {unrec=tab(i,4);l=l*unrec;}
    p+=R::dnbinom(round(unrec*tab(i,1)),mu*l/sigma,1.0-sigma/(1.0+sigma),1);
  }
  return(p);
}

// [[Rcpp::export]]
double likelihoodCarcC(NumericMatrix tab, double mu, double sigma) {
  int n = (tab.nrow()+1)/2;
  double p=0;
  double unrec=1;
  double l;
  double minbralen=tab(n,1);
  //if (NumericVector::is_na(minbralen)) minbralen=0.1;
  for (int i=0;i<tab.nrow();i++) {
    if (i==n) continue;
    if (tab.ncol()==5) unrec=tab(i,4);
    l=unrec*(tab(i,2)-tab(tab(i,3)-1,2));
    if (unrec*tab(i,1)<=minbralen)
      p+=R::pgamma(minbralen,mu*l/(1.0+sigma),1.0+sigma,1,1);
    else
      p+=R::dgamma(unrec*tab(i,1),mu*l/(1.0+sigma),1.0+sigma,1);
  }
  return(p);
}
