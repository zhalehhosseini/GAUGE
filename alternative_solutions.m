function [x_total,f_k1]=alternative_solutions(U,model,S,H,max_obj)

% calculating alternative solutions for adding reactions to the model
%INPUTS:
%   U: structure representing universal dataset of reactions with 3 fields
%       U.S: stoichiometric matrix of the reactions in the database
%       U.lb: lower bounds of reactions
%       U.ub: upper bounds of reactions
%   model: original metabolic model with 3 fields
%       model.S: stoichiometric matrix of the reactions in the database
%       model.lb: lower bounds of reactions
%       model.ub: upper bounds of reactions
%   S: a n*3 matrix with n=number of slightly correlated fully coupled reaction pairs. 
%       in every row, first and second entries are the number of reactions in the pair in model. 
%       the third entry is the flux ratio of the reactions in the pair 
%   H: a m*3 matrix with m=number of highly correlated fully coupled reaction pairs. 
%       in every row, first and second entries are the number of reactions in the pair in model. 
%       the third entry is the flux ratio of the reactions in the pair 
%   max_obj: objective value of the first step. maximum number of resolved
%       inconsistencies
%OUTPUTS:
%   x_total: all of the alternative solutions. each column shows a possible solution 
%   f_k1: minimum number of reactions which should be added to the model

x_total=[];
   
[f_k1, x,c, Ain, x_L, x_U, b_L, b_U,IntVars] = second_step(U,model,S,H,max_obj);
x_total(:,1)=x;
       
        f_k=f_k1;
q=size(Ain,1);
i=1;

while abs(f_k-f_k1)<1e-6
    i=i+1;
    
       
        y=x(size(model.S,2)+size(U.S,2)+1:size(model.S,2)+2*size(U.S,2));
        t= abs(y)<1e-6;y(t)=0;
        t= abs(y)>0.999999;y(t)=1;
        for j=1:size(U.S,2)
            if ~y(j)
                Ain(q+i-1,size(model.S,2)+size(U.S,2)+j)=-1;
            else
                Ain(q+i-1,size(model.S,2)+size(U.S,2)+j)=1;
            end
        end
        
    temp=find(y);
    b_U(end+1)=length(temp)-1;
    b_L(end+1)=-inf;
    [x, slack, v, rc, f_k, ninf, sinf, Inform] = cplex(c, Ain, x_L, x_U, b_L, b_U,[],[],[],[],IntVars);
    x_total(:,i)=x;
    
end
    
end