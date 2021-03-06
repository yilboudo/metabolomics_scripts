####################################################################
####################################################################
# This script calculates association statistics  
# between metabolites and a given phenotypes
# using linear or logistic regression providing 
# ajusted p-values in permutation procedure. 
#
#
# The input data should contain both the metabolites data and 
# phenotype(s) of interest. For the script to work, ensure that
# all metabolites data follow each other in the table. In such a
# way that you can give the script the name of the first metabolite
# (called NFML in the script) and that of the last metabolite 
# (called NLML in the script).
#
# Example run:
# source("minp_permutations.R")
# 
# For Linear Model
#
#
# minp_permutations(
# input_filename = "Data_with_Phenotypes_and_Metabolites.txt",
# output_filename = "HbF_Metabolites_Association.txt"
# phenotype = "HbF"
# testStat = "linear", 
# NFML = "phenylalanine.d8", 
# NLML = "C26.carnitine", 
# Nperm = 100,
# output_filename = "output.txt")
#
#
# For Logistic Model
#
#
# minp_permutations(
# phenotype = "Survival",
# input_filename = "Data_with_Phenotypes_and_Metabolites.txt",
# output_filename = "Survival_Metabolites_Association.txt",
# testStat = "logistic",
# NFML = "phenylalanine.d8", 
# NLML = "C26.carnitine", 
# Nperm = 100 
# output_filename = "output.txt")
#
####################################################################
####################################################################