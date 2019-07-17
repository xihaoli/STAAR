// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
List matsp(arma::mat G) {

	int i;
	int n = G.n_rows;
	int p = G.n_cols;
	
	arma::vec maf;
	maf.zeros(p);
	
	// maf and flip
	for(i = 0; i < p; i++)
	{
		maf(i) = arma::mean(G(arma::span(0,n - 1),i))/2;
	}

	return List::create(Named("Geno") = arma::sp_mat(G), Named("maf") = maf);

}

