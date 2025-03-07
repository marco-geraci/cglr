//* C++ code for the cglr package by Marco Geraci, Sapienza University of Rome*//

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace Rcpp;
using namespace Eigen;

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
List gaussQuad(int n, std::string dist = "uniform", double l = 0, double u = 1, double mu = 0, double sigma = 1, double alpha = 1, double beta = 1) {
	if (n < 0) {
        throw std::invalid_argument("need non-negative number of nodes");
    }
    if (n == 0) {
        return List::create(Named("nodes") = NumericVector(0), 
                            Named("weights") = NumericVector(0));
    }
    std::vector<std::string> valid_dists = {"uniform", "beta1", "beta2", "normal", "beta", "gamma"};
    if (std::find(valid_dists.begin(), valid_dists.end(), dist) == valid_dists.end()) {
        throw std::invalid_argument("invalid distribution");
    }
    if (n == 1) {
        double x;
        if (dist == "uniform") {
            x = (l + u) / 2;
        } else if (dist == "beta1" || dist == "beta2" || dist == "beta") {
            x = alpha / (alpha + beta);
        } else if (dist == "normal") {
            x = mu;
        } else if (dist == "gamma") {
            x = alpha * beta;
        }
        return List::create(Named("nodes") = NumericVector::create(x), Named("weights") = NumericVector::create(1));
    }
    if (dist == "beta" && alpha == 0.5 && beta == 0.5) {
        dist = "beta1";
    }
    if (dist == "beta" && alpha == 1.5 && beta == 1.5) {
        dist = "beta2";
    }
    VectorXd a(n), b(n - 1);
    if (dist == "uniform") {
        a.setZero();
        for (int i = 0; i < n - 1; ++i) {
            b(i) = (i + 1) / std::sqrt(4 * std::pow(i + 1, 2) - 1);
        }
    } else if (dist == "beta1") {
        a.setZero();
        b.setConstant(0.5);
        b(0) = std::sqrt(0.5);
    } else if (dist == "beta2") {
        a.setZero();
        b.setConstant(0.5);
    } else if (dist == "normal") {
        a.setZero();
        for (int i = 0; i < n - 1; ++i) {
            b(i) = std::sqrt((i + 1) / 2.0);
        }
    } else if (dist == "beta") {
        double ab = alpha + beta;
        a(0) = (alpha - beta) / ab;
        for (int i = 1; i < n; ++i) {
            double abi = ab - 2 + 2 * (i + 1);
            a(i) = ((alpha - 1) * (alpha - 1) - (beta - 1) * (beta - 1)) / ((abi - 2) * abi);
        }
        b(0) = std::sqrt(4 * alpha * beta / (ab * ab * (ab + 1)));
        for (int i = 1; i < n - 1; ++i) {
            double abi = ab - 2 + 2 * (i + 1);
            b(i) = std::sqrt(4 * (i + 1) * (i + 1 + alpha - 1) * (i + 1 + beta - 1) * (i + 1 + ab - 2) / ((abi * abi - 1) * abi * abi));
        }
    } else if (dist == "gamma") {
        for (int i = 0; i < n; ++i) {
            a(i) = 2 * (i + 1) + alpha - 2;
        }
        for (int i = 0; i < n - 1; ++i) {
            b(i) = std::sqrt((i + 1) * (i + 1 + alpha - 1));
        }
    }
    b.conservativeResize(n);
    b(n - 1) = 0;

    // Eigenvalue decomposition
    MatrixXd T = MatrixXd::Zero(n, n);
    T.diagonal() = a;
    T.diagonal(1) = b.head(n - 1);
    T.diagonal(-1) = b.head(n - 1);
    SelfAdjointEigenSolver<MatrixXd> solver(T);
    VectorXd x = solver.eigenvalues();
    VectorXd w = solver.eigenvectors().row(0).array().square();

    if (dist == "uniform") {
        x = l + (u - l) * (x.array() + 1) / 2;
    } else if (dist == "beta1" || dist == "beta2" || dist == "beta") {
        x = (x.array() + 1) / 2;
    } else if (dist == "normal") {
        x = mu + std::sqrt(2) * sigma * x.array();
    } else if (dist == "gamma") {
        x = beta * x.array();
    }

    return List::create(Named("nodes") = x, Named("weights") = w);	
}

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
double dmvnrm_arma(arma::rowvec const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
    int xdim = x.n_cols;
    double out;
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;

    arma::rowvec z = x - mean;
    inplace_tri_mat_mult(z, rooti);
	out = other_terms - 0.5 * arma::dot(z, z);     
      
    if (logd)
      return out;
    return exp(out);
}

// [[Rcpp::export]]
double dmvpn_arma(arma::rowvec const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
    double out;
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
 
 	// z%*%Sigmainv%*%z
	arma::rowvec z = x;
    inplace_tri_mat_mult(z, rooti);
	double f1 = arma::dot(z, z);
	// z%*%Sigmainv%*%m
	arma::rowvec m = mean;
    inplace_tri_mat_mult(m, rooti);
	double f2 = arma::dot(m, z);
	// m%*%Sigmainv%*%m
	double f3 = arma::dot(m, m);
	// f2/sqrt(f1)
	double f4 = f2/sqrt(f1);
	
	out = -0.5*f3 - log2pi - log(f1) - 0.5*log(arma::det(sigma)) + log(1 + f4*arma::normcdf(f4)/arma::normpdf(f4));
	
    if (logd)
      return out;
    return exp(out);
}

//[[Rcpp::export]]
arma::vec C_loglik_proj_cglr(arma::mat w, arma::vec theta, arma::mat xb, arma::mat sigma, arma::mat muw, arma::vec alpha, int Hv){

	int n = w.n_rows;
	int K = w.n_cols;
	double LOGC = -std::numeric_limits<double>::infinity();
	bool single_gq = (alpha.size() == 1);

	arma::vec ans(n, arma::fill::zeros);
	arma::mat tmp(n, Hv, arma::fill::zeros);
	arma::mat sigmasw = sigma.submat(K, 0, K+1, K-1);
	arma::mat sigmaws = sigma.submat(0, K, K-1, K+1);
	arma::mat sigmaw = sigma.submat(0, 0, K-1, K-1);
	arma::mat sigmas = sigma.submat(K, K, K+1, K+1);
	arma::mat SS = sigmasw * arma::inv(sigmaw);
	arma::mat mw = xb.submat(0,0,n-1,K-1);
	arma::mat ms = xb.submat(0,K,n-1,K+1);
	//arma::mat ew = w - mw;
	//ms = ms + ew*SS.t();
	
	sigmas = sigmas - SS * sigmaws;
	arma::vec ct = cos(theta);
	arma::vec st = sin(theta);
	arma::mat z = arma::join_rows(ct, st);

	List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[0], 1);
	NumericVector rv = GQ["nodes"];
	NumericVector wv = GQ["weights"];
	//Rcpp::Rcout << alpha << std::endl;
	//Rcpp::Rcout << rv << std::endl;
	//Rcpp::Rcout << wv << std::endl;
	
	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;
		if (!single_gq) {
			List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[i], 1);
			NumericVector rv = GQ["nodes"];
			NumericVector wv = GQ["weights"];
		}
		
		arma::rowvec zi = z.row(i);
		arma::rowvec mis = ms.row(i);

		arma::rowvec wi = w.row(i);
		arma::rowvec miw = mw.row(i);
	
		//Rcpp::Rcout << mis << std::endl;
		for(int gq = 0; gq < Hv; ++gq){
			miw = rv(gq)*muw.row(i) + miw;
			arma::rowvec ew = wi - miw;
			mis = mis + ew*SS.t();
			// Integrated log-likelihood
			tmp(i,gq) = dmvnrm_arma(wi, miw, sigmaw*rv(gq), true) + dmvpn_arma(zi, mis, sigmas*rv(gq), true) + log(fabs(wv(gq)));
			if(tmp(i,gq) > LOGC) LOGC = tmp(i,gq); // constant for OVERFLOW protection
		} // GQ
	}  // i
	
	for(int i = 0; i < n; ++i){
		double val = 0;
		for(int gq = 0; gq < Hv; ++gq){
			val += copysign(1, wv(gq)) * exp(tmp(i,gq) - LOGC);
		}
		//Rcpp::Rcout << LOGC + log(val) << std::endl;
		ans(i) = LOGC + log(val);
	}

	return(ans);
}

//[[Rcpp::export]]
arma::vec C_loglik_cond_cglr(arma::mat y, arma::mat xb, arma::mat sigma, arma::mat mu, arma::vec alpha, int Hv){

	int n = y.n_rows;
	double LOGC = -std::numeric_limits<double>::infinity();
	bool single_gq = (alpha.size() == 1);
	
	arma::vec ans(n, arma::fill::zeros);
	arma::mat tmp(n, Hv, arma::fill::zeros);

	List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[0], 1);
	NumericVector rv = GQ["nodes"];
	NumericVector wv = GQ["weights"];

	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;
		if (!single_gq) {
			List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[i], 1);
			NumericVector rv = GQ["nodes"];
			NumericVector wv = GQ["weights"];
		}

		arma::rowvec yi = y.row(i);

		for(int gq = 0; gq < Hv; ++gq){

			arma::rowvec mi = rv(gq)*mu.row(i) + xb.row(i);

			// Integrated log-likelihood
			tmp(i,gq) = dmvnrm_arma(yi, mi, rv(gq)*sigma, true) + log(fabs(wv(gq)));
			if(tmp(i,gq) > LOGC) LOGC = tmp(i,gq); // constant for OVERFLOW protection
		} // GQ
	} // i

	for(int i = 0; i < n; ++i){
		double val = 0;
		for(int gq = 0; gq < Hv; ++gq){
			val += copysign(1, wv(gq)) * exp(tmp(i,gq) - LOGC);
		}
		ans(i) = LOGC + log(val);
	}
	//Rcpp::Rcout << ans << std::endl;

	return(ans);
}

//[[Rcpp::export]]
double C_loglik_proj_cpnr(arma::mat w, arma::vec theta, arma::mat m, arma::mat sigma){

	double ans = 0;
	int n = w.n_rows;
	int K = w.n_cols;

	arma::mat sigmasw = sigma.submat(K, 0, K+1, K-1);
	arma::mat sigmaws = sigma.submat(0, K, K-1, K+1);
	arma::mat sigmaw = sigma.submat(0, 0, K-1, K-1);
	arma::mat sigmas = sigma.submat(K, K, K+1, K+1);
	arma::mat SS = sigmasw * arma::inv(sigmaw);
	arma::mat mw = m.submat(0,0,n-1,K-1);
	arma::mat ms = m.submat(0,K,n-1,K+1);
	
	sigmas = sigmas - SS * sigmaws;
	arma::vec ct = cos(theta);
	arma::vec st = sin(theta);
	arma::mat z = arma::join_rows(ct, st);

	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;
		
		arma::rowvec zi = z.row(i);
		arma::rowvec mis = ms.row(i);

		arma::rowvec wi = w.row(i);
		arma::rowvec miw = mw.row(i);
	
		ans += dmvnrm_arma(wi, miw, sigmaw, true) + dmvpn_arma(zi, mis, sigmas, true);
	}  // i
	
	return(ans);
}

//[[Rcpp::export]]
double C_loglik_cond_cpnr(arma::mat y, arma::mat m, arma::mat sigma){

	double ans = 0;
	int n = y.n_rows;

	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;

		arma::rowvec yi = y.row(i);
		arma::rowvec mi = m.row(i);
		ans += dmvnrm_arma(yi, mi, sigma, true);
	} // i

	//Rcpp::Rcout << ans << std::endl;

	return(ans);
}


//* Densities *//


//[[Rcpp::export]]
double C_dpgl(double theta, arma::rowvec m, arma::mat sigma, double alpha, int Hv, bool const logd = false){

	double LOGC = -std::numeric_limits<double>::infinity();
	double ans = 0;
	arma::vec tmp(Hv, arma::fill::zeros);
	
	arma::rowvec z = {std::cos(theta), std::sin(theta)};

	List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha, 1);
	NumericVector rv = GQ["nodes"];
	NumericVector wv = GQ["weights"];
	
	//Rcpp::Rcout << i << std::endl;
	for(int gq = 0; gq < Hv; ++gq){
		// Integrated log-likelihood of PN to get projected GL
		tmp(gq) = dmvpn_arma(z, m, sigma*rv(gq), true) + log(fabs(wv(gq)));
		if(tmp(gq) > LOGC) LOGC = tmp(gq); // constant for OVERFLOW protection
	} // GQ
	
	for(int gq = 0; gq < Hv; ++gq){
		ans += copysign(1, wv(gq)) * exp(tmp(gq) - LOGC);
	}
	//Rcpp::Rcout << LOGC + log(val) << std::endl;
	ans = LOGC + log(ans);

    if (logd)
      return ans;
    return exp(ans);
}

//[[Rcpp::export]]
arma::vec C_dpgl_vv(arma::vec theta, arma::mat m, arma::mat sigma, arma::vec alpha, int Hv, bool const logd = false){

	int n = theta.n_elem;
	double LOGC = -std::numeric_limits<double>::infinity();
	bool single_gq = (alpha.size() == 1);

	arma::vec ans(n, arma::fill::zeros);
	arma::mat tmp(n, Hv, arma::fill::zeros);
	arma::vec ct = cos(theta);
	arma::vec st = sin(theta);
	arma::mat z = arma::join_rows(ct, st);

	List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[0], 1);
	NumericVector rv = GQ["nodes"];
	NumericVector wv = GQ["weights"];
	
	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;
		if (!single_gq) {
			List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha[i], 1);
			NumericVector rv = GQ["nodes"];
			NumericVector wv = GQ["weights"];
		}

		arma::rowvec zi = z.row(i);
		arma::rowvec mi = m.row(i);
		
		for(int gq = 0; gq < Hv; ++gq){
			// Integrated log-likelihood of PN to get projected GL
			tmp(i,gq) = dmvpn_arma(zi, mi, sigma*rv(gq), true) + log(fabs(wv(gq)));
			if(tmp(i,gq) > LOGC) LOGC = tmp(i,gq); // constant for OVERFLOW protection
			//Rcpp::Rcout << tmp(i,gq) << std::endl;
		} // GQ
	}  // i
	
	for(int i = 0; i < n; ++i){
		double val = 0;
		for(int gq = 0; gq < Hv; ++gq){
			val += copysign(1, wv(gq)) * exp(tmp(i,gq) - LOGC);
		}
		//Rcpp::Rcout << LOGC + log(val) << std::endl;
		ans(i) = LOGC + log(val);
	}

    if (logd)
      return ans;
    return exp(ans);
}

//[[Rcpp::export]]
arma::vec C_dpn_vv(arma::vec theta, arma::mat m, arma::mat sigma, bool const logd = false){

	int n = theta.n_elem;
	
	arma::vec ans(n, arma::fill::zeros);
	arma::vec ct = cos(theta);
	arma::vec st = sin(theta);
	arma::mat z = arma::join_rows(ct, st);

	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;
		arma::rowvec zi = z.row(i);
		arma::rowvec mi = m.row(i);
		ans(i) = dmvpn_arma(zi, mi, sigma, true);

	}  // i
	
    if (logd)
      return ans;
    return exp(ans);
}


//[[Rcpp::export]]
arma::vec C_dmgl_vv(arma::mat y, arma::rowvec theta, arma::mat sigma, arma::rowvec mu, double alpha, int n, int Hv, bool const logd = false){

	double LOGC = -std::numeric_limits<double>::infinity();

	arma::vec ans(n, arma::fill::zeros);
	arma::mat tmp(n, Hv, arma::fill::zeros);

	List GQ = gaussQuad(Hv, "gamma", 0, 1, 0, 1, alpha, 1);
	NumericVector rv = GQ["nodes"];
	NumericVector wv = GQ["weights"];
	
	for(int i = 0; i < n; ++i){
		//Rcpp::Rcout << i << std::endl;

		arma::rowvec yi = y.row(i);
		
		for(int gq = 0; gq < Hv; ++gq){
			// Integrated log-likelihood of PN to get projected GL
			tmp(i,gq) = dmvnrm_arma(yi, theta + mu*rv(gq), sigma*rv(gq), true) + log(fabs(wv(gq)));
			if(tmp(i,gq) > LOGC) LOGC = tmp(i,gq); // constant for OVERFLOW protection
			//Rcpp::Rcout << tmp(i,gq) << std::endl;
		} // GQ
	}  // i
	
	for(int i = 0; i < n; ++i){
		double val = 0;
		for(int gq = 0; gq < Hv; ++gq){
			val += copysign(1, wv(gq)) * exp(tmp(i,gq) - LOGC);
		}
		//Rcpp::Rcout << LOGC + log(val) << std::endl;
		ans(i) = LOGC + log(val);
	}

    if (logd)
      return ans;
    return exp(ans);
}
