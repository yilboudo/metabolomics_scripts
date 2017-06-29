minp_permutations <- function(phenotype,input_filename, testStat = "linear", NFML, NLML, Nperm,output_filename, sex, age, ...)
{
  set.seed(4)
  
  Metabo_Pheno_File<- read.table(input_filename, h=T, sep="\t")
  #Define Test Statistic and Phenotype of Interest. Pick which test statistic to run linear or logistic regression
  
  
  linearReg_p <- function(y,x,w,z) coef(summary(lm(y ~ x + w + z)))[2,4]
  linearReg_stderr <- function(y,x,w,z) coef(summary(lm(y ~ x + w + z)))[2,2]
  linearReg_effect <- function(y,x,w,z) coef( summary(lm(y ~ x + w + z)) ) [2,1]
  linearReg_samplesize <- function(y,x, w, z) nrow(model.frame(lm(y ~ x + w + z)))

  
  logisticReg_p <- function(y,x, w, z) coef(summary(glm(y ~ x + w + z,family=binomial(link='logit'))))[2,4]
  logisticReg_stderr <- function(y,x, w, z) coef(summary(glm(y ~ x + w + z,family=binomial(link='logit'))))[2,2]
  logisticReg_effect <- function(y,x, w, z) coef(summary(glm(y ~ x + w + z,family=binomial(link='logit'))))[2,1]
  logisticReg_samplesize <- function(y,x, w, z) nrow(model.frame(glm(y ~ x + w + z,family=binomial(link='logit'))))



  PofI= Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)]
  if missing(sex){
     Sex = NA} else {
  Sex = Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==sex)]
  }
 
  if missing(age){
     Age = NA} else {
  
  Age = Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==age)]
  }

  #Name of first Metabolite in list and last
  
  fmlist <- which(colnames(Metabo_Pheno_File)==NFML)
  lmlist <- which(colnames(Metabo_Pheno_File)==NLML)
  
  
  if(testStat =="linear"){

    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=PofI, w=Sex, z=Age)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_stderr, y=PofI, w=Sex, z=Age)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_effect, y=PofI, w=Sex, z=Age)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_samplesize, y=PofI, w=Sex, z=Age)

  } else if(testStat=="logistic") {
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=PofI, w=Sex, z=Age)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_stderr, y=PofI, w=Sex, z=Age)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_effect, y=PofI, w=Sex, z=Age)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_samplesize, y=PofI, w=Sex, z=Age)

  } else {
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=PofI, w=Sex, z=Age)
    std_err <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_stderr, y=PofI, w=Sex, z=Age)
    effect_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_effect, y=PofI, w=Sex, z=Age)
    sample_size <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_samplesize, y=PofI, w=Sex, z=Age)
  }
  
  #Define Number of metabolites and number of simulations
  Nmetabo=length(fmlist:lmlist)
  permutations_pval <- matrix(NA,Nmetabo,Nperm)
  
  #Calculate test statistic for each permutations one for logistic one for linear
  cat(" Calculating Permutation Matrix...\n")
  if(testStat=="linear"){
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=sample(PofI, replace=F, size=length(PofI)), w=Sex, z=Age)}
  }
  else if(testStat=="logistic"){
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg_p, y=sample(PofI, replace=F, size=length(PofI)),w=Sex, z=Age)}
  } else {
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg_p, y=sample(PofI, replace=F, size=length(PofI)),w=Sex, z=Age)}
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
  cat(" Results ready!\n")
}
