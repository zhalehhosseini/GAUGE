function [A_lambda]=lambda_calculation(A,model)
% calculating flux ratio
% INPUTS:
%   A: fully coupled reaction pairs
%   model: metabolic model in COBRA format
% OUTPUTS:
%   A_lambda: a n*3 matrix with n=number of slightly correlated fully coupled reaction pairs. 
%       in every row, first and second entries are the number of reactions in the pair in model. 
%       the third entry is the flux ratio of the reactions in the pair 
A_lambda=A;
for i=1:size(A_lambda,1)
    new=model;
    new.lb(A_lambda(i,1))=1;
    new.ub(A_lambda(i,1))=1;
    new.c(:)=0;
    new.c(A_lambda(i,2))=1;
    [sol]=optimizecbmodel(new);
    [sol2]=optimizecbmodel(new,'min');
    if sol.stat==1 & sol.f==sol2.f
        A_lambda(i,3)=1/sol.f;
    else
        new.lb(A_lambda(i,1))=-1;
        new.ub(A_lambda(i,1))=-1;
        [sol]=optimizecbmodel(new);
        [sol2]=optimizecbmodel(new,'min');
        A_lambda(i,3)=-1/sol.f;
    end
end
for i=1:size(A_lambda,1)
    if abs(A_lambda(i,3))<1
        A_lambda(i,3)=1/A_lambda(i,3);
        temp=A_lambda(i,2);
        A_lambda(i,2)=A_lambda(i,1);
        A_lambda(i,1)=temp;
    end
end

end