// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matrix_impute(arma::mat G) {

	int i,j;
	int n = G.n_rows;
	int p = G.n_cols;

	arma::vec AF;
	AF.zeros(p);

	double num = 0;

	// Calculate AF
	for(i = 0; i < p; i++)
	{
		num = 0;

		for(j = 0; j < n; j++)
		{
			if(G(j,i) > -1)
			{
				AF(i) = AF(i) + G(j,i);
				num = num + 1;
			}
		}

		AF(i) = AF(i)/2/num;
	}

	// Genotype Imputation
	for(i = 0; i < p; i++)
	{
		if(AF(i) <= 0.5)
		{
			for(j = 0; j < n; j++)
			{
				if(!(G(j,i) > -1))
				{
					G(j,i) = 0;
				}
			}
		}else
		{
			for(j = 0; j < n; j++)
			{
				if(!(G(j,i) > -1))
				{
					G(j,i) = 2;
				}
			}
		}
	}

	return G;
}

