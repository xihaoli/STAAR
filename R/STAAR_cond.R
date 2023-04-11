#' STAAR procedure for conditional analysis using omnibus test
#'
#' The \code{STAAR_cond} function takes in genotype, the genotype of variants to be
#' adjusted for in conditional analysis, the object from fitting the null
#' model, and functional annotation data to analyze the conditional association between a
#' quantitative/dichotomous phenotype and a variant-set by using STAAR procedure,
#' adjusting for a given list of variants. For each variant-set, the conditional
#' STAAR-O p-value is a p-value from an omnibus test that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1)
#' together with conditional p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param genotype_adj an n*p_adj genotype matrix (dosage matrix) of the target
#' sequence, where n is the sample size and p_adj is the number of genetic variants
#' to be adjusted for in conditional analysis (or a vector of a single variant with length n
#' if p_adj is 1).
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, STAAR will perform the
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), ACAT-V(1,1)
#' and ACAT-O tests (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{genotype_adj}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{genotype_adj}
#' and taking the residuals (default = \code{optimal}).
#' @return A list with the following members:
#' @return \code{num_variant}: the number of variants with minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set using STAAR.
#' @return \code{cMAC}: the cumulative minor allele count of variants with
#' minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{RV_label}: the boolean vector indicating whether each variant in the given
#' variant-set has minor allele frequency > 0 and less than \code{rare_maf_cutoff}.
#' @return \code{results_STAAR_O_cond}: the conditional STAAR-O p-value that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) together
#' with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_ACAT_O_cond}: the conditional ACAT-O p-value that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) using Cauchy method.
#' @return \code{results_STAAR_S_1_25_cond}: a vector of conditional STAAR-S(1,25) p-values,
#' including conditional SKAT(1,25) p-value weighted by MAF, the conditional SKAT(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-S(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_S_1_1_cond}: a vector of conditional STAAR-S(1,1) p-values,
#' including conditional SKAT(1,1) p-value weighted by MAF, the conditional SKAT(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-S(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_B_1_25_cond}: a vector of conditional STAAR-B(1,25) p-values,
#' including conditional Burden(1,25) p-value weighted by MAF, the conditional Burden(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-B(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_B_1_1_cond}: a vector of conditional STAAR-B(1,1) p-values,
#' including conditional Burden(1,1) p-value weighted by MAF, the conditional Burden(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-B(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_A_1_25_cond}: a vector of conditional STAAR-A(1,25) p-values,
#' including conditional ACAT-V(1,25) p-value weighted by MAF, the conditional ACAT-V(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-A(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_A_1_1_cond}: a vector of conditional STAAR-A(1,1) p-values,
#' including conditional ACAT-V(1,1) p-value weighted by MAF, the conditional ACAT-V(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-A(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @references Li, Z., Li, X., et al. (2020). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @export

STAAR_cond <- function(genotype,genotype_adj,obj_nullmodel,annotation_phred=NULL,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       method_cond=c("optimal","naive")){

  method_cond <- match.arg(method_cond) # evaluate choices
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }

  if(inherits(genotype_adj, "numeric")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)

      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
      w_A <- cbind(w_A_1,w_A_2)
      w_A <- as.matrix(w_A)
    }

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P
        P_scalar <- sqrt(dim(P)[1])
        P <- P*P_scalar

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)
        PX_adj <- P%*%X_adj
        P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
          PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
          X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
        rm(P)
        gc()

        pvalues <- STAAR_O_SMMAT(G,P_cond,residuals.phenotype,
                                 weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                 mac=as.integer(round(MAF*2*dim(G)[1])))
      }else{
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)

        pvalues <- STAAR_O_SMMAT_sparse_cond(G,Sigma_i,Sigma_iX,cov,X_adj,residuals.phenotype,
                                             weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                             mac=as.integer(round(MAF*2*dim(G)[1])))
      }
    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        P <- diag(working) - X%*%solve(t(X)%*%diag(working)%*%X)%*%t(X)
      }else if(obj_nullmodel$family[1] == "gaussian"){
        P <- diag(length(working)) - X%*%solve(t(X)%*%X)%*%t(X)
      }

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      if(method_cond == "optimal"){
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+model.matrix(obj_nullmodel)-1)
      }else{
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
      }
      residuals.phenotype <- residuals.phenotype.fit$residuals
      X_adj <- model.matrix(residuals.phenotype.fit)
      PX_adj <- P%*%X_adj
      P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
        PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
        X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
      rm(P)
      gc()

      pvalues <- STAAR_O_SMMAT(G,P_cond,residuals.phenotype,
                               weights_B=w_B,weights_S=w_S,weights_A=w_A,
                               mac=as.integer(round(MAF*2*dim(G)[1])))
    }

    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2]+1
    results_STAAR_O <- CCT(pvalues)
    results_ACAT_O <- CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
    pvalues_STAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
    pvalues_STAAR_S_1_1 <- CCT(pvalues[(num_annotation+1):(2*num_annotation)])
    pvalues_STAAR_B_1_25 <- CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
    pvalues_STAAR_B_1_1 <- CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
    pvalues_STAAR_A_1_25 <- CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
    pvalues_STAAR_A_1_1 <- CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])

    results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_S_1_25)
    results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))

    results_STAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_S_1_1)
    results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))

    results_STAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

    results_STAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

    results_STAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_STAAR_A_1_25)
    results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))

    results_STAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_STAAR_A_1_1)
    results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))

    if(dim(annotation_phred)[2] == 0){
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)","STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)","STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)","STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)","STAAR-A(1,1)")
    }else{
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                          paste0("SKAT(1,25)-",colnames(annotation_phred)),
                                          "STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)",
                                         paste0("SKAT(1,1)-",colnames(annotation_phred)),
                                         "STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-",colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                          paste0("ACAT-V(1,25)-",colnames(annotation_phred)),
                                          "STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                         paste0("ACAT-V(1,1)-",colnames(annotation_phred)),
                                         "STAAR-A(1,1)")
    }

    return(list(num_variant = num_variant,
                cMAC = cMAC,
                RV_label = RV_label,
                results_STAAR_O_cond = results_STAAR_O,
                results_ACAT_O_cond = results_ACAT_O,
                results_STAAR_S_1_25_cond = results_STAAR_S_1_25,
                results_STAAR_S_1_1_cond = results_STAAR_S_1_1,
                results_STAAR_B_1_25_cond = results_STAAR_B_1_25,
                results_STAAR_B_1_1_cond = results_STAAR_B_1_1,
                results_STAAR_A_1_25_cond = results_STAAR_A_1_25,
                results_STAAR_A_1_1_cond = results_STAAR_A_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

