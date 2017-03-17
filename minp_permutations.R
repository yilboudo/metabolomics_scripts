minp_permutations <- function(phenotype,input_filename, testStat = "linear", NFML, NLML, Nperm,output_filename)
{
  set.seed(4)
  
  Metabo_Pheno_File<- read.table(input_filename, h=T, sep="\t")
  #Define Test Statistic and Phenotype of Interest. Pick which test statistic to run linear or logistic regression
  
  linearReg <- function(y, x) coef(summary(lm(y ~ x)))[2,4]
  logisticReg <- function(y,x) coef(summary(glm(y ~ x ,family=binomial(link='logit'))))[2,4]
  PofI=Metabo_Pheno_File[,which(colnames(Metabo_Pheno_File)==phenotype)]
  
  #Name of first Metabolite in list and last
  
  fmlist <- which(colnames(Metabo_Pheno_File)==NFML)
  lmlist <- which(colnames(Metabo_Pheno_File)==NLML)
  
  
  if(testStat =="linear"){
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg, y=PofI)
  } else if(testStat=="logistic") {
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg, y=PofI)
  } else {
    observed_pval <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg, y=PofI)
  }
  
  #Define Number of metabolites and number of simulations
  Nmetabo=length(fmlist:lmlist)
  permutations_pval <- matrix(NA,Nmetabo,Nperm)
  
  #Calculate test statistic for each permutations one for logistic one for linear
  cat(" Calculating Permutation Matrix...\n")
  if(testStat=="linear"){
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg, x=sample(PofI, replace=F, size=length(PofI)))}
  }
  else if(testStat=="logistic"){
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], logisticReg, y=sample(PofI, replace=F, size=length(PofI)))}
  } else {
    for (i in 1:Nperm) {permutations_pval[,i] <- sapply(Metabo_Pheno_File[,fmlist:lmlist], linearReg, x=sample(PofI, replace=F, size=length(PofI)))}
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
  
  all_pval<- cbind(observed_pval,emp_pval,t(pval_perm))
  format(all_pval, scientific = F)
  
  #Write to file
  Outputfile <- file(output_filename,'w')
  writeLines(paste("\t\t",phenotype), Outputfile)
  writeLines(paste("Metabolites","Observed_pvalue","Empirical_pvalue","Permutation_pvalue",sep="\t"), Outputfile)
  write.table(all_pval, file = Outputfile, quote = FALSE, row.names = TRUE, col.names = FALSE,sep="\t")
  close(Outputfile)
  cat(" Results ready!\n")
}