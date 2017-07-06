
cat(" ################################# \n # \n # \n # \n # \n # Metabolomics Phenotype Association: \n # \n # Yann Ilboudo, Guillaume Lettre 2016 \n # \n #\n # \n ################################# \n")


minp_permutations <- function(phenotype,input_filename, testStat = "linear", NFML, NLML, Nperm,output_filename, use.cov =F, cov_names, inv.norm=F)
{
  start.time <- proc.time()[3]
  set.seed(4)
  
  Metabo_Pheno_File<- read.table(input_filename, h=T, sep="\t")

  Metabo_Pheno_File <- Metabo_Pheno_File[complete.cases(Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)]),]

  invnorm <- function(x) 
  {norm <- (x - mean(x ))/ sd(x) 
   znorm<- qnorm((rank(norm,na.last="keep")-0.5)/sum(!is.na(norm)))
   return(znorm)
  }

  #Define Test Statistic and Phenotype of Interest. Pick which test statistic to run linear or logistic regression

  if (inv.norm==T) { 
    
    PofI = invnorm(Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)])
  
  } else {
    
    PofI= Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)]
  }
  
  #Name of first and last metabolite in file
  
  fmlist <- which(colnames(Metabo_Pheno_File)==NFML)
  lmlist <- which(colnames(Metabo_Pheno_File)==NLML)

  
  
  if (use.cov==T && testStat =="linear") {

    
    linearRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,4]
    linearRegcov_stderr <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,2]
    linearRegcov_effect <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,1]
    linearRegcov_samplesize <- function(y,...) nrow(model.frame(glm(y ~ . , data.frame(y,...) ,family="gaussian")))
    
    covariates <- strsplit(cov_names, split=",")[[1]]
    cov=Metabo_Pheno_File[,covariates]

    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_p, y=PofI, data=cov)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_stderr, y=PofI, data=cov)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_effect, y=PofI, data=cov)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_samplesize, y=PofI, data=cov)

    } else if (use.cov == T && testStat=="logistic") {

    logisticRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="binomial")))[2,4]
    logisticRegcov_effect <- function(y,...) exp(coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,1])
    logisticRegcov_confint_upper   <- function(y,...) exp(coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,1] + (1.96*coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,2])) 
    logisticRegcov_confint_lower <- function(y,...) exp(coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,1] - (1.96*coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,2])) 
    logisticRegcov_samplesize <- function(y,...) nrow(model.frame(glm(y ~ . , data.frame(y,...),family="binomial")))
    
    covariates <- strsplit(cov_names, split=",")[[1]]
    cov=Metabo_Pheno_File[,covariates]

    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_p, y=PofI, data=cov)
    confint_upper <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_confint_upper, y=PofI, data=cov)
    confint_lower <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_confint_lower, y=PofI, data=cov)
    odds_ratio <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_effect, y=PofI, data=cov)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_samplesize, y=PofI, data=cov)
      
    } else if (use.cov == F && testStat=="linear") {

    linearReg_p <- function(y,x) coef(summary(glm(y ~ x ,family="gaussian")))[2,4]
    linearReg_stderr <- function(y,x) coef(summary(glm(y ~ x  ,family="gaussian")))[2,2]
    linearReg_effect <- function(y,x) coef( summary(glm(y ~ x ,family="gaussian" )) ) [2,1]
    linearReg_samplesize <- function(y,x) nrow(model.frame(glm(y ~ x ,family="gaussian")))

    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=PofI)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_stderr, y=PofI )
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_effect, y=PofI)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_samplesize, y=PofI)

    } else if (use.cov == F && testStat=="logistic") {

    logisticReg_p <- function(y,x) coef(summary(glm(y ~ x ,family="binomial")))[2,4]
    logisticReg_confint_upper <- function(y,x) exp(coef(summary(glm(y ~ x ,family="binomial")))[2,1] + (1.96*coef(summary(glm(y ~ x , family="binomial")))[2,2]))
    logisticReg_confint_lower <- function(y,x) exp(coef(summary(glm(y ~ x , family="binomial")))[2,1] - (1.96* coef(summary(glm(y ~ x , family="binomial")))[2,2]))
    logisticReg_effect <- function(y,x) exp(coef(summary(glm(y ~ x , family="binomial")))[2,1])
    logisticReg_samplesize <- function(y,x) nrow(model.frame(glm(y ~ x ,family="binomial")))
    
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=PofI)
    confint_upper <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_confint_upper, y=PofI)
    confint_lower <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_confint_lower, y=PofI)
    odds_ratio <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_effect, y=PofI)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_samplesize, y=PofI)

    } else {

      print("NA")
    }
  




  #Define Number of metabolites and number of permutations
  Nmetabo=length(fmlist:lmlist)
  permutations_pval <- matrix(NA,Nmetabo,Nperm)
  
  #Calculate test statistic for each permutations one for logistic one for linear
  cat(" \n \n ################################# \n # \n # Computing Permutation Matrix...\n # \n ################################# \n")
  if (testStat=="linear" && use.cov==F) {
    linearReg_p <- function(y,x) coef(summary(glm(y ~ x ,family="gaussian")))[2,4]

    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=sample(PofI, replace=F, size=length(PofI)))}
 
  } else if (testStat=="linear" && use.cov==T) {

    linearRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,4]    
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_p, y=sample(PofI, replace=F, size=length(PofI)), data=cov)}

  } else if (testStat=="logistic" && use.cov==F) {
    logisticReg_p <- function(y,x) coef(summary(glm(y ~ x ,family="binomial")))[2,4]

    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=sample(PofI, replace=F, size=length(PofI)))}
  
  } else if (testStat=="logistic" && use.cov==T) {

    logisticRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="binomial")))[2,4]
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_p, y=sample(PofI, replace=F, size=length(PofI)), data=cov)}
  
  } else {

    print("Zero. Nothing to do.")
  }

  cat("\n")
  cat(" ################################# \n # \n # Finished!\n # \n ################################# \n\n")
  #Store smallest pvalue for each rounds
  smallest_pval_per_round <- matrix(NA, 1,Nperm)
  for(i in 1:Nperm) {smallest_pval_per_round[,i] <- min(permutations_pval[,i])}
  
  pval_perm <- matrix(NA,1,length(fmlist:lmlist))
  for (i in 1:Nmetabo) { pval_perm[,i] <- mean(smallest_pval_per_round <= observed_pval[i])}
  
  #Empirical Pvalue and Permutation Pvalue
  emp_pval <- sapply(seq_along(fmlist:lmlist), function(x) mean(permutations_pval[x,] <= observed_pval[x])) 
  
  if (testStat=="logistic") {
  
  all_values<- cbind(sample_size, odds_ratio, confint_upper,confint_lower, observed_pval, emp_pval,t(pval_perm))

    } else if (testStat=="linear") {

  all_values<- cbind(sample_size, effect_size, std_err, observed_pval, emp_pval,t(pval_perm))
  
  } else {

    print("NA")
  }

  all_values <- apply(all_values,2, signif, digits=3)
  #format(all_metrics, scientific = F)
  
  #Write to file
  Outputfile <- file(output_filename,'w')
  writeLines(paste("\t\t\t",phenotype), Outputfile)
  
  if (testStat=="logistic") {

    writeLines(paste("Metabolites","Sample Size","Odds Ratio","CI-U","CI-L","Observed_pvalue","Empirical_pvalue","Permutation_pvalue",sep="\t"), Outputfile)

    
  } else if (testStat=="linear") {

    writeLines(paste("Metabolites","Sample Size","Effect Size","Standard Error","Observed_pvalue","Empirical_pvalue","Permutation_pvalue",sep="\t"), Outputfile)
  
  }  else {
    print ("NA")
  }
  write.table(all_values, file = Outputfile, quote = FALSE, row.names = TRUE, col.names = FALSE,sep="\t")
  close(Outputfile)
  cat(" ################################# \n # \n # Results Ready!\n # \n ################################# \n")


    running.time <- proc.time()[3] - start.time 
    if (running.time < 60){
      out.run.time <- round(running.time, digits = 2)
      out.time.units <- "seconds"
    }
    if (running.time > 60 && running.time < 3600){
      out.run.time <- round((running.time/60), digits = 2)
      out.time.units <- "minutes"
    }
    if (running.time > 3600){
      out.run.time <- round((running.time/3600), digits = 2)
      out.time.units <- "hours"
    }
    
      print(paste("RUNNING TIME: ",   out.run.time, out.time.units, sep = " "))
    
    q(status = 0)

}
