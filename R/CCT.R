#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT.pval} function takes in a vector of p-values, and return the
#' aggregated p-value using Cauchy method.
#' @param Pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @return \code{pval}: the aggregated p-value combining p-values from the vector
#' \code{Pvals}.
#' @examples p.values <- c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples CCT.pval(Pvals=p.values)
#' @references Liu, Y. and Xie, J. (2019). Cauchy combination test: a powerful test with analytic
#' p-value calculation under arbitrary dependency structures.
#' \emph{Journal of American Statistical Association (In press)}.
#' (\href{https://amstat.tandfonline.com/doi/abs/10.1080/01621459.2018.1554485?journalCode=uasa20}{pub})
#' @export

CCT.pval<-function(Pvals){
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-mean(tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((1/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(tan((0.5-Pvals[!is.small])*pi))
    cct.stat<-cct.stat/length(Pvals)
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}

