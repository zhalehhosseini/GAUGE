function [alternatives,f_max,f_min]=GAUGE(model,PCC,F2C2_Solver,SC_PCC,HC_PCC,U)
% INPUTS:
%   model: metabolic model in COBRA format with known gene_reaction
%       associations in model.rules and model.rxnGeneMat field
%   PCC: a matrix containing pearson correlation coefficients of gene pairs
%       number of columns=number of rows=number of genes in the model
%       PCC(i,j)=pearson correlation coefficient of gene pair (i,j)
%   F2C2_solver: the solver which is used in F2C2 function, (glpk,Lindo,clp,SoPlex,linprog)
%   SC_PCC: cutoff for choosing slightly correlated reaction pairs
%   HC_PCC: cutoff for choosing highly correlated reaction pairs
%   U: structure representing universal dataset of reactions with 3 fields
%       U.S: stoichiometric matrix of the reactions in the database
%       U.lb: lower bounds of reactions
%       U.ub: upper bounds of reactions
% OUTPUTS:
%   alternatives: all of the alternative solutions. each column shows a
%       possible solution
%   f_min: minimum number of reactions which should be added to the model
%   f_max: maximum number of inconsistancies which can be resolved

    [SC_lambda,HC_lambda]=SC_HC_calculator(model,PCC,F2C2_Solver,SC_PCC,HC_PCC);
    f_max = first_step(U,model,SC_lambda,HC_lambda);
    [alternatives,f_min]=alternative_solutions(U,model,S,H,f_k);
    

end