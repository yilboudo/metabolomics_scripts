start.time <- proc.time()[3]
options(echo = TRUE)
options(warn=-1)
library(batch)

cat(" ################################# \n # \n # \n # \n # \n # Metabolomics Phenotype Association: \n # \n # Yann Ilboudo, Guillaume Lettre 2016 \n # \n #\n # \n ################################# \n")


minp_permutations <- function(phenotype,input_filename, testStat = "linear", NFML, NLML, Nperm,output_filename, use.cov = F, cov_names, inv.norm=F)
{
  set.seed(4)
  
  Metabo_Pheno_File<- read.table(input_filename, h=T, sep="\t")
  #Define Test Statistic and Phenotype of Interest. Pick which test statistic to run linear or logistic regression
  

  PofI= Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)]

  
  #Name of first and last metabolite in file
  
  fmlist <- which(colnames(Metabo_Pheno_File)==NFML)
  lmlist <- which(colnames(Metabo_Pheno_File)==NLML)

  covariates <- strsplit(cov_names, split=",")[[1]]

  cov=Metabo_Pheno_File[,covariates]
  
  if (use.cov==T && testStat =="linear") {

    
    linearRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,4]
    linearRegcov_stderr <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,2]
    linearRegcov_effect <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="gaussian")))[2,1]
    linearRegcov_samplesize <- function(y,...) nrow(model.frame(glm(y ~ . , data.frame(y,...) ,family="gaussian")))
      
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_p, y=PofI, data=cov)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_stderr, y=PofI, data=cov)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_effect, y=PofI, data=cov)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_samplesize, y=PofI, data=cov)

    } else if (use.cov == T && testStat=="logistic") {

    logisticRegcov_p <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...) ,family="binomial")))[2,4]
    logisticRegcov_stderr <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,2]
    logisticRegcov_effect <- function(y,...) coef(summary(glm(y ~ . , data.frame(y,...),family="binomial")))[2,1]
    logisticRegcov_samplesize <- function(y,...) nrow(model.frame(glm(y ~ . , data.frame(y,...),family="binomial")))
      
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_p, y=PofI, data=cov)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_stderr, y=PofI, data=cov)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_effect, y=PofI, data=cov)
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
    logisticReg_stderr <- function(y,x) coef(summary(glm(y ~ x ,family="binomial")))[2,2]
    logisticReg_effect <- function(y,x) coef(summary(glm(y ~ x ,family="binomial")))[2,1]
    logisticReg_samplesize <- function(y,x) nrow(model.frame(glm(y ~ x ,family="binomial")))
    
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=PofI)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_stderr, y=PofI)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_effect, y=PofI)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_samplesize, y=PofI)

    } else {

      print("Zero. Nothing to do.")
    }
  




  #Define Number of metabolites and number of simulations
  Nmetabo=length(fmlist:lmlist)
  permutations_pval <- matrix(NA,Nmetabo,Nperm)
  
  #Calculate test statistic for each permutations one for logistic one for linear
  cat(" \n \n ################################# \n # \n # Computing Permutation Matrix...\n # \n ################################# \n")
  if (testStat=="linear" && use.cov==F) {

    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=sample(PofI, replace=F, size=length(PofI)))}
 
  } else if (testStat=="linear" && use.cov==T) {
    
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearRegcov_p, y=sample(PofI, replace=F, size=length(PofI)), data=cov)}

  } else if (testStat=="logistic" && use.cov==F) {

    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=sample(PofI, replace=F, size=length(PofI)))}
  
  } else if (testStat=="logistic" && use.cov==T) {

    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticRegcov_p, y=sample(PofI, replace=F, size=length(PofI)), data=cov)}
  
  } else {

    print("Zero. Nothing to do.")
  }

  cat("\n")
  cat(" Finished!\n\n")
  #Store smallest pvalue for each rounds
  smallest_pval_per_round <- matrix(NA, 1,Nperm)
  for(i in 1:Nperm) {smallest_pval_per_round[,i] <- min(permutations_pval[,i])}
  
  pval_perm <- matrix(NA,1,length(fmlist:lmlist))
  for (i in 1:Nmetabo) { pval_perm[,i] <- mean(smallest_pval_per_round <= observed_pval[i])}
  
  #Empirical Pvalue and Permutation Pvalue
  emp_pval <- sapply(seq_along(fmlist:lmlist), function(x) mean(permutations_pval[x,] <= observed_pval[x]))
  
  all_values<- cbind(sample_size,effect_size,std_err,observed_pval,emp_pval,t(pval_perm))
  #format(all_metrics, scientific = F)
  
  #Write to file
  Outputfile <- file(output_filename,'w')
  writeLines(paste("\t\t\t",phenotype), Outputfile)
  writeLines(paste("Metabolites","Sample Size","Effect Size","Standard Error","Observed_pvalue","Empirical_pvalue","Permutation_pvalue",sep="\t"), Outputfile)
  write.table(all_values, file = Outputfile, quote = FALSE, row.names = TRUE, col.names = FALSE,sep="\t")
  close(Outputfile)
  cat(" ################################# \n # \n # Results Ready!\n # \n ################################# \n")


  #if(print.time){
  #cat(" ################################# \n # \n #   Print time \n # \n ################################# \n")
    #}
    running.time <- proc.time()[3] - start.time 
    if (running.time < 60){
      out.run.time <- round(running.time, digits = 2)
      out.time.units <- "seconds"
    }
    if (running.time > 60 & running.time < 3600){
      out.run.time <- round((running.time/60), digits = 2)
      out.time.units <- "minutes"
    }
    if (running.time > 3600){
      out.run.time <- round((running.time/3600), digits = 2)
      out.time.units <- "hours"
    }
    #if(print.time){
      print(paste("RUNNING TIME: ",   out.run.time, out.time.units, sep = " "))
    #}
    quit()
  #}
}
