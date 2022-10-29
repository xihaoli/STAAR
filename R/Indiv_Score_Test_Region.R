#' Score test for individual variants in a given variant-set
#'
#' The \code{Indiv_Score_Test_Region} function takes in genotype and the object from fitting the null
#' model to analyze the associations between a quantitative/dichotomous phenotype and
#' all individual variants in a given variant-set by using score test.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a data frame with p rows corresponding to the p genetic variants in the given variant-set
#' and three columns: \code{Score} (the score test statistic), \code{SE} (the standard error associated
#' with the score test statistic), and \code{pvalue} (the score test p-value).
#' If a variant in the given variant-set has minor allele frequency = 0 or
#' greater than \code{rare_maf_cutoff}, the corresponding row will be \code{NA}. If a variant in
#' the given variant-set has standard error equal to 0, the p-value will be set as 1.
#' @export

Indiv_Score_Test_Region <- function(genotype,obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2){

  if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  results <- data.frame(Score = rep(NA, dim(genotype)[2]),
                        SE = rep(NA, dim(genotype)[2]),
                        pvalue = rep(NA, dim(genotype)[2]))

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype,MAF)
  gc()

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    rm(Geno_rare)
    gc()

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P
        P_scalar <- sqrt(dim(P)[1])
        P <- P*P_scalar

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)

        results[RV_label,] <- do.call(cbind,Indiv_Score_Test_SMMAT(G,P,residuals.phenotype))
      }else{
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals

        results[RV_label,] <- do.call(cbind,Indiv_Score_Test_SMMAT_sparse(G,Sigma_i,Sigma_iX,cov,residuals.phenotype))
      }
    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        fam <- 1
      }else if(obj_nullmodel$family[1] == "gaussian"){
        fam <- 0
      }

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      results[RV_label,] <- do.call(cbind,Indiv_Score_Test(G,X,working,sigma,fam,residuals.phenotype))
    }

    return(results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

