// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// declare K
double K(double x, arma::vec egvalues);
// declare K1 (first derivative)
double K1(double x, arma::vec egvalues, double q);
// declare K2 (second derivative)
double K2(double x, arma::vec egvalues);
// declare bisection
double Bisection(arma::vec egvalues, double q, double xmin, double xmax);
// declare saddlepoint
double Saddle(double q, arma::vec egvalues);
// declare CCT_pval
double CCT_pval(arma::vec x, arma::vec weights);

// [[Rcpp::export]]
arma::vec STAAR_O_SMMAT(arma::sp_mat G, arma::mat P, arma::vec residuals, arma::mat weights_B, arma::mat weights_S, arma::mat weights_A, arma::vec mac, int mac_thres=10)
{
	int i,k;
	int j;
	double sum0 = 0.0;
	double sumw = 0.0;
	double sumx = 0.0;
	
	bool lower = false;
	bool logp = false;
	

	// int vn = G.n_rows;
	int un = G.n_cols;
	int wn = weights_B.n_cols;
	
	arma::vec res;
	res.zeros(3*wn);
	
	arma::mat Cov;
	Cov.zeros(un,un);
	
	arma::mat Covw;
	Covw.zeros(un,un);
	
	arma::rowvec x = trans(residuals)*G;

	arma::vec eigenvals;
	eigenvals.zeros(un);
	

	arma::mat Wleft;
	Wleft.zeros(un,un);
	
	arma::mat Wright;
	Wright.zeros(un,un);
	
	double c1 = 0.0;
	double c2 = 0.0;
	double c4 = 0.0;
	double l = 0.0;
	
	int ii;
	
	Cov = trans(G)*P*G;
	

	//ACAT
	arma::uvec id_veryrare = arma::find(mac <= mac_thres); 
	arma::uvec id_common = arma::find(mac > mac_thres); 
	
	
	int n0 = id_veryrare.size();
	int n1 = id_common.size();
	
	arma::vec pseq;
	pseq.zeros(un);
	
	arma::vec wseq;
	wseq.zeros(un);
	
	// double SSR = 0.0;
	// double SST = 0.0;
	
	for(k = 0; k < n1; k++)
	{
		pseq(k) = pow(x(id_common(k)),2)/Cov(id_common(k),id_common(k));
		pseq(k) = R::pchisq(pseq(k),1,lower,logp);
	}
	
	
	
	
	int n = x.size();
	for(i = 0; i < wn; i++)
	{
		// SKAT
		Wright.each_row() = trans(weights_S.col(i));
		Wleft.each_col() = weights_S.col(i);
		
		Covw = Wleft%Cov%Wright;
		
		sum0 = 0.0;
		for (k = 0; k < n; k++)
		{
			sum0 = sum0 + pow(x(k), 2) *pow(weights_S(k,i), 2);
		}
	
		eigenvals = arma::eig_sym(Covw);
		for(j = 0; j < un; j++)
		{
			if(eigenvals(j) < 1e-8)
			{
				eigenvals(j) = 0.0;
			}
		}

	
		res(i) = Saddle(sum0,eigenvals);
		
		if(res(i)== 2)
		{
			c1 = 0.0;
			c2 = 0.0;
			c4 = 0.0;
			
			for(ii = 0; ii < un; ii++)
			{
				c1 = c1 + Covw(ii, ii);
			}
		
			Covw = Covw*Covw;
		
			for(ii = 0; ii < un; ii++)
			{
				c2 = c2 + Covw(ii, ii);
			}
		
			Covw = Covw*Covw;

			for(ii = 0; ii < un; ii++)
			{
				c4 = c4 + Covw(ii, ii);
			}
		
			sum0 = (sum0 - c1)/sqrt(2*c2);
			l = pow(c2,2)/c4;
			res(i) = R::pchisq(sum0*sqrt(2*l)+l,l,lower,logp);
		}
		
		// Burden
		Wright.each_row() = trans(weights_B.col(i));
		Wleft.each_col() = weights_B.col(i);
		
		Covw = Wleft%Cov%Wright;
		
		sumw = arma::accu(Covw);
		
		sum0 = 0.0;
		for (k = 0; k < n; k++)
		{
			sum0 = sum0 + x(k)*weights_B(k,i);
		}
		
		sumx = pow(sum0, 2) / sumw;
	
		res(wn + i) = R::pchisq(sumx,1,lower,logp);
		
		// ACAT
		for(k = 0; k < n1; k++)
		{
			wseq(k) = weights_A(id_common(k),i);
		}
		
		if(n0 == 0)
		{
			res(2*wn + i) = CCT_pval(pseq(arma::span(0,n1-1)),wseq(arma::span(0,n1-1)));
		}else
		{
			sum0 = 0.0;
			sumw = 0.0;
			for (k = 0; k < n0; k++)
			{
				sum0 = sum0 + x(id_veryrare(k))*weights_B(id_veryrare(k),i);
				sumw = sumw + weights_A(id_veryrare(k),i);
			}

			sumx = pow(sum0, 2) / arma::accu(Covw.submat(id_veryrare,id_veryrare));
			pseq(n1) = R::pchisq(sumx,1,lower,logp);
			wseq(n1) = sumw/n0;
			res(2*wn + i) = CCT_pval(pseq(arma::span(0,n1)),wseq(arma::span(0,n1)));
		}	
		
		
	}
	
	
	return res;
}


double K(double x, arma::vec egvalues)
{
    double res = 0.0;
    const int en = egvalues.size();
    
    for(int i = 0; i < en; i++)
    {
        res = res + log(1-2*egvalues(i)*x);
    }
    
    res = res*(-0.5);
    
    return res;
}

double K1(double x, arma::vec egvalues, double q)
{
    double res = 0.0;
    const int en = egvalues.size();
    
    for(int i = 0; i < en; i++)
    {
        res = res + egvalues(i)/(1-2*egvalues(i)*x);
    }
    
    res = res - q;
    
    return res;
}

double K2(double x, arma::vec egvalues)
{
    double res = 0.0;
    const int en = egvalues.size();
    
    for(int i = 0; i < en; i++)
    {
        res = res + pow(egvalues(i),2)/pow(1-2*egvalues(i)*x,2.0);
    }
    
    res = res*2;
    
    return res;
}

double Bisection(arma::vec egvalues, double q, double xmin, double xmax)
{
   
    // the range of x to search
    double xupper = xmax;
    double xlower = xmin;
    
    double x0 = 0.0;
    double K1x0 = 1.0;
	
   
    while (fabs(xupper-xlower) > 1e-08)
    {
		x0 = (xupper + xlower)/2.0;
		K1x0 = K1(x0,egvalues,q);
		
		if(K1x0 == 0){
			break;
        }else if (K1x0 > 0){
            xupper = x0;
        }else{
            xlower = x0;
        }
    }
    
    return x0;
}


double Saddle(double q, arma::vec egvalues)
{
	bool lower = false;
	bool logp = false;
	double xmin = 0.0;
	double xmax = 0.0;
	double v = 0.0;
	double res = 0.0;
	
	int i;
	const int en = egvalues.size();
	double lambdamax = max(egvalues);
	q = q/lambdamax;
	
	for(i = 0; i < en; i++)
	{
		egvalues[i] = egvalues[i]/lambdamax;
	}
	lambdamax = 1.0;
	
    if (q > arma::sum(egvalues)) 
	{
        xmin = -0.01;
	}else
	{
        xmin = -en/(2 * q);
    }
	
    xmax = 1/(2*lambdamax) * 0.99999;
	
    double xhat = Bisection(egvalues,q,xmin,xmax);
    double w = sqrt(2*(xhat*q-K(xhat,egvalues)));
    if (xhat < 0){
        w = -w;
    }
	
	v = xhat*sqrt(K2(xhat,egvalues));
	if(fabs(xhat)<1e-04)
	{
		res = 2;
	}else
	{
		res =  R::pnorm(w+log(v/w)/w,0.0,1.0,lower,logp);
	}
    return res;


}



double CCT_pval(arma::vec x, arma::vec weights)
{
	double cct_stat = 0.0;
	double pval = 0.0;
	
	bool lower = false;
	bool logp = false;
	
	weights = weights/sum(weights);

	int n = x.size();
	int k ;


	for(k = 0;k < n;k++)
	{
		if(x(k) < 1e-16)
		{
			cct_stat = cct_stat + weights(k)/x(k)/PI;
		}else
		{
			cct_stat = cct_stat + weights(k)*tan((0.5-x(k))*PI);
		}
	}
		
	if (cct_stat > 1e+15)
	{
        pval = (1/cct_stat)/PI;
    }else
	{
		pval = R::pcauchy(cct_stat,0.0,1.0,lower,logp);
    }
	
	return pval;
		
}







