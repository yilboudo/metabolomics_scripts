Brief explanation of the variables used in the function minp_permutations.R


input_filename: Name of the your input file

output_filename: Name of your output file

phenotype: Phenotype of interest, the column header corresponding 
           to the phenotype you want to run your association tests

testStat: It stands for test statistic. Depending on whether your phenotype is continuous or categorical,
          you can choose the appropriate test statistic, logistic or linear. The default is linear.

NFML: Name of First Metabolite in List, this is the column header corresponding
      to the first metabolite in your table.

NLML: Name of Last Metabolite in List, this is the column header corresponding
      to the last metabolite in your table.

Nperm: Number of permutations applied to the test statistic. 

Note: All the metabolites columns in your table should follow each other. They don't have to be in any particular order. 
