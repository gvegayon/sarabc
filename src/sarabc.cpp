// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// This function generates p matrix of the form W^1, ..., W^p
// to be used later for computing approximate (In - rho*W)
// [[Rcpp::export]]
std::vector< arma::sp_mat > Powers(const arma::sp_mat & X, int p=5) {
  std::vector< arma::sp_mat > Xs(p);
  Xs[0] = X;
  int i=0;
  while(++i < p) 
    Xs.at(i)=Xs[i-1] * X;
  
  return Xs;
}

// Simulates y in the sar model by approximation
// y = (I_n - rho*W)^(-1)(X * beta)
// where (In - rho*W)^(-1) = I_n + rho W + rho^2W^2 + ...
void pred_sar_fast(
    arma::sp_colvec & y,
    arma::sp_mat & TMP,
    const std::vector< arma::sp_mat > & Ws,
    const arma::sp_mat & X,
    double rho,
    arma::sp_colvec & beta
  ) {
  
  int n = y.n_rows;
  TMP = arma::speye(n,n) + rho*Ws[0];
  
  for (int i=1;i<Ws.size();i++)
    TMP += pow(rho, i+1)*Ws[i];
  
  y = TMP * (X*beta);
    
  return;
}

// Function that can be called from R. This is a wrapper of sim_sar_fast
// which requires passing a vector of Gs (power graphs) a container for the
// output y and a sparse temp matrix.
// [[Rcpp::export]]
arma::colvec pred_sar(
    arma::sp_mat W, arma::mat X, double rho, arma::colvec beta, int p=5) {
  
  // Initializing values
  int n=X.n_rows;
  std::vector< arma::sp_mat > Ws = Powers(W,p);
  arma::sp_mat TMP(n,n);
  arma::sp_colvec ysp(n);
  arma::sp_mat Xsp(X);
  arma::sp_colvec betasp(beta);
  
  pred_sar_fast(ysp, TMP, Ws, Xsp, rho, betasp);
  
  arma::colvec y(ysp);
  return y;
  
}

// Comparative Stats -----------------------------------------------------------

// Computes fitness aprox
void fitness(double & val, arma::sp_colvec y,
             arma::sp_colvec & Wy, arma::sp_mat & X,
             double rho, arma::sp_colvec & beta) {

  val = arma::norm(y - rho*Wy - X*beta);
  
  return;
  
}

void pred_sar_naive(
    arma::sp_colvec & val,
    const arma::sp_colvec & Wy, 
    const arma::sp_mat & X,
    double rho, arma::sp_colvec & beta) {
  val = rho*Wy + X*beta;
  return;
}

// Computes moran's I
void moran(
    double & val,
    const arma::sp_colvec & x,
    const arma::sp_mat & W,
    double & Wsum) {
  
  arma::colvec xcent(x);
  xcent-=arma::mean(x);
  double numer = accu((xcent * xcent.t()) % W);
  double denom = accu(arma::pow(xcent, 2.0));
  
  val = (x.size()/Wsum)*(numer/(denom + 1e-15));
  return;
}

// Function to compute distances -----------------------------------------------

void dist(
    double & val,
    std::vector<double> & y,
    std::vector<double> & yhat,
    const std::vector<double> & sweights) {
  
  int n = y.size();
  
  val= 0;
  for(int i=0;i<n;i++)
    val += pow(y[i] - yhat[i],2.0)*sweights[i];
  
  val = sqrt(val);
  
  return;
}

void distArma(
  double & val,
  const arma::sp_colvec & y,
  arma::sp_colvec & yhat) {
  int n = y.size();
  val=0;
  for (int i=0;i<n;i++)
    val += pow(y[i] - yhat[i],2.0)/n;
  
  val = sqrt(val);
  
  return;
}


// SAR ABC ---------------------------------------------------------------------
// [[Rcpp::export]]
List sar_abc_cpp(const arma::colvec & y, const arma::mat & X, const arma::sp_mat & W,
             const std::vector<double> & rho,
             const arma::mat & beta, 
             const std::vector<double> & sweights,
             int N=1e3, int p = 5,
             bool no_inv = false, bool no_moran=true) {
  
  // Phase 1: Generate variables -----------------------------------------------
  int n = X.n_rows;
  int k = X.n_cols;

  // Common variables
  arma::sp_colvec ysp(y);
  arma::sp_mat Xsp(X);
  arma::sp_mat betasp(beta.t());

  // Computing aprox inverse
  const std::vector< arma::sp_mat > Ws = Powers(W,p);
  arma::sp_mat TMP(n,n);
  arma::sp_colvec yhat_sp(y);
  
  // Computing MSE
  const arma::sp_colvec Wy = W*ysp;
  
  // Computing summary stats:
  //  - moran's I
  //  - mean
  //  - sd
  std::vector< double > stats0(3), stats1(3);
  double Wsum = accu(W);
  
  // Computing stats0
  NumericMatrix stats(N,3);
  stats0[0] = 0;
  stats0[1] = arma::mean(y);
  if (!no_moran) stats0[2] = 0;
  else stats[2] = 0;
  
  // Outputs
  NumericVector D(N);
  
  // Phase 2: Randomize and fit the data
  arma::sp_colvec bsp(k);
  for (int i=0;i<N;i++) {
    bsp = betasp.col(i);
    
    // Prediction
    if (!no_inv) pred_sar_fast(yhat_sp, TMP,  Ws, Xsp, rho[i], bsp);
    else        pred_sar_naive(yhat_sp, Wy, Xsp, rho[i], bsp);
    
    // Prediction error, mean and morans I
    distArma(stats1[0], ysp, yhat_sp);
    stats1[1] =     arma::mean(yhat_sp);
    
    if (!no_moran) moran(stats1[2], ysp - yhat_sp, W, Wsum);
    else stats1[2] = 0; 
    
    // Computing distances
    dist(D[i], stats0, stats1, sweights);
    stats(i,0) = stats1[0];
    stats(i,1) = stats1[1];
    stats(i,2) = stats1[2];
  }
  
  // Phase 3: Returning
  return List::create(
    _["rho"] = rho,
    _["beta"] = beta,
    _["distance"] = D,
    _["stats0"] = stats0,
    _["stats1"] = stats
  );
  
}

