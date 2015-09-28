function f_k = first_step(U,model,S,H)
% first step of GAUGE. finding maximum number of inconssistencies which can
% be resolved
% INPUTS:
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
% note that model.S and U.S should have the same number of rows.
% common metabolites in model and U should be at the same rows. for
% metabolites from U(model) which are not present in model(U) rows
% with all zero elements should be added to model.S(U.S).
% OUTPUTS:
%   f_k=maximum number of inconsistancies which can be resolved


n_ori=size(model.S,2);
n_new=size(U.S,2);
Aeq=[model.S,U.S];




Ain(1:n_new,n_ori+1:n_ori+n_new)=sparse(1:n_new,1:n_new,-1);
Ain(1:n_new,n_ori+n_new+1:n_ori+2*n_new)=sparse(1:n_new,1:n_new,U.lb);
b_U(1:n_new)=0;
b_L(1:n_new)=U.lb-U.ub;

Ain(n_new+1:2*n_new,n_ori+1:n_ori+n_new)=sparse(1:n_new,1:n_new,1);
Ain(n_new+1:2*n_new,n_ori+n_new+1:n_ori+2*n_new)=sparse(1:n_new,1:n_new,-U.ub);
b_U(n_new+1:2*n_new)=0;
b_L(n_new+1:2*n_new)=U.lb-U.ub;

for i=1:size(S,1)
    Ain(2*n_new+i,S(i,1))=-1;
    Ain(2*n_new+i,S(i,2))=S(i,3);
    Ain(2*n_new+i,n_ori+2*n_new+i)=-1001;
    Ain(2*n_new+i,n_ori+2*n_new+size(S,1)+i)=-1001;
    b_U(2*n_new+i)=-1e-6;
    if S(i,3)>0
        b_L(2*n_new+i)=-model.ub(S(i,1))+S(i,3)*model.lb(S(i,2))-1001*2;
    else
        b_L(2*n_new+i)=-model.ub(S(i,1))+S(i,3)*model.ub(S(i,2))-1001*2;
    end
end

for i=1:size(S,1)
    Ain(2*n_new+size(S,1)+i,n_ori+2*n_new+3*size(S,1)+i)=1;
    Ain(2*n_new+size(S,1)+i,n_ori+2*n_new+size(S,1)+i)=1001;
    b_U(2*n_new+size(S,1)+i)=1001;
    b_L(2*n_new+size(S,1)+i)=0;
end

for i=1:size(S,1)
    Ain(2*n_new+2*size(S,1)+i,S(i,1))=1;
    Ain(2*n_new+2*size(S,1)+i,S(i,2))=-S(i,3);
    Ain(2*n_new+2*size(S,1)+i,n_ori+2*n_new+i)=1001;
    Ain(2*n_new+2*size(S,1)+i,n_ori+2*n_new+2*size(S,1)+i)=-1001;
    b_U(2*n_new+2*size(S,1)+i)=1001-10e-6;
    if S(i,3)>0
        b_L(2*n_new+2*size(S,1)+i)=model.lb(S(i,1))-S(i,3)*model.ub(S(i,2))-1001;
    else
        b_L(2*n_new+2*size(S,1)+i)=model.lb(S(i,1))-S(i,3)*model.lb(S(i,2))-1001;
    end
end

for i=1:size(S,1)
    Ain(2*n_new+3*size(S,1)+i,n_ori+2*n_new+3*size(S,1)+i)=1;
    Ain(2*n_new+3*size(S,1)+i,n_ori+2*n_new+2*size(S,1)+i)=1001;
    b_U(2*n_new+3*size(S,1)+i)=1001;
    b_L(2*n_new+3*size(S,1)+i)=0;
end

for i=1:size(H,1)
    Ain(2*n_new+4*size(S,1)+i,H(i,1))=1;
    Ain(2*n_new+4*size(S,1)+i,H(i,2))=-H(i,3);
    b_U(2*n_new+4*size(S,1)+i)=0;
    b_L(2*n_new+4*size(S,1)+i)=0;
end




Aeq=[Aeq,sparse(size(model.S,1),size(Ain,2)-n_ori-n_new)];
beq=zeros(size(model.S,1),1);
            
Ain=[Ain;Aeq];
b_U=[b_U,beq'];
b_L=[b_L,beq'];




c=zeros(1,n_ori+2*n_new+4*size(S,1));

c(n_ori+2*n_new+3*size(S,1)+1:n_ori+2*n_new+4*size(S,1))=1;
x_L=[model.lb;U.lb;zeros(size(Ain,2)-n_ori-n_new,1)];

x_U=[model.ub;U.ub;ones(size(Ain,2)-n_ori-n_new,1)];

IntVars=false(length(x_L),1);
IntVars(n_new+n_ori+1:length(x_L))=1;

c=c';

[x, slack, v, rc, f_k, ninf, sinf, Inform] = cplex(-c, Ain, x_L, x_U, b_L, b_U,[],[],[],[],IntVars);
f_k=-f_k;

end

