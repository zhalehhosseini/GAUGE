function [SC_lambda,HC_lambda]=SC_HC_calculator(model,PCC,F2C2_Solver,SC_PCC,HC_PCC)
%calculates slightly correlated and highly correlated reaction pairs in
%metabolic model and corresponding flux ratios
%INPUT:
%   model: metabolic model in COBRA format with known gene_reaction
%       associations in model.rules and model.rxnGeneMat field
%   PCC: a matrix containing pearson correlation coefficients of gene pairs
%       number of columns=number of rows=number of genes in the model
%       PCC(i,j)=pearson correlation coefficient of gene pair (i,j)
%   F2C2_solver: the solver which is used in F2C2 function, (glpk,Lindo,clp,SoPlex,linprog)
%   SC_PCC: cutoff for choosing slightly correlated reaction pairs
%   HC_PCC: cutoff for choosing highly correlated reaction pairs
%OUTPUT:
%   SC_lambda: slightly correlated reaction pairs in the model with their
%   flux ratio
%   HC: highly correlated reaction pairs in the model with their
%   flux ratio      

%calculating gene coupling relations
GCM=GCA(model);
%slightly correlated and highly correlated reaction pairs
[SC,HC]=SCHC_pairs(model,GCM,PCC,F2C2_Solver,SC_PCC,HC_PCC);
%calculating flux ratios for SC and HC
[SC_lambda]=lambda_calculation(SC,model);
[HC_lambda]=lambda_calculation(HC,model);
end