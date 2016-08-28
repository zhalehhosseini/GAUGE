function [SC,HC]=SCHC_pairs(model,GCM,PCC,F2C2_solver,PCC_SC,PCC_HC)
% calculating slightly correlated and highly correlated reaction pairs
% INPUTS:
%   model: metabolic network model in COBRA format including model.rxnGeneMat
%       field
%   GCM: gene coupling matrix; number of columns=number of rows=number of
%       genes in the model
%   PCC: a matrix containing pearson correlation coefficients of gene pairs
%       number of columns=number of rows=number of genes in the model
%       PCC(i,j)=pearson correlation coefficient of gene pair (i,j)
%   F2C2_solver: the solver which is used in F2C2 function, (glpk,Lindo,clp,SoPlex,linprog)
%   PCC_SC: cutoff for choosing slightly correlated reaction pairs
%   PCC_HC: cutoff for choosing highly correlated reaction pairs
% OUTPUTS:
%   SC: slightly correlated reaction pairs in the model
%   HC: highly correlated reaction pairs in the model

F2C2_model=CobraToF2C2(model);
[FCM2,blocked]=F2C2(F2C2_solver,F2C2_model);
x=[];y=[];GCM=triu(GCM);
GCM(1:size(GCM,1)+1:size(GCM,1)*size(GCM,1))=0;
FCM = zeros(length(model.rxns), length(model.rxns));
FCM(~blocked, ~blocked) = FCM2;
FCM=triu(FCM);
FCM(1:size(FCM,1)+1:size(FCM,1)*size(FCM,1))=0;
HC=[];
SC=[];
for i=1:size(GCM,1)
    for j=1:size(GCM,1)
        if GCM(i,j)==1 
            if abs(PCC(i,j))<PCC_SC
           
            a=find(model.rxnGeneMat(:,i));
            b=find(model.rxnGeneMat(:,j));
            for k=1:length(a)
                for l=1:length(b)
                    
                    if FCM(a(k),b(l))==1 || FCM(b(l),a(k))==1 
                        if a(k)<b(l)
                            x=(b(l)-1)*size(FCM,1)+a(k);
                        SC=[SC;[a(k),b(l)]];
                        
                        else
                            x=(a(k)-1)*size(FCM,1)+b(l);
                        SC=[SC;[b(l),a(k)]];
                        end
                    end
                    
                end
            end
            end
            if abs(PCC(i,j))>PCC_HC
            
            a=find(model.rxnGeneMat(:,i));
            b=find(model.rxnGeneMat(:,j));
            for k=1:length(a)
                for l=1:length(b)
                    if FCM(a(k),b(l))==1 || FCM(b(l),a(k))==1 
                        
                        if a(k)<b(l)
                             y=(b(l)-1)*size(FCM,1)+a(k);
                        HC=[HC;[a(k),b(l)]];
                        else
                            y=(a(k)-1)*size(FCM,1)+b(l);
                        HC=[HC;[b(l),a(k)]];
                        end
                    end
                    
                end
            end
            end
        end
    end
end

if ~isempty(x)
    for i=1:length(x)
        if i<=length(x)
            a=find(x==x(i));
            if length(a)>1
                for j=length(a):-1:2
                x(a(j))=[];
                SC(a(j),:)=[];
                end
            end
        end
    end
end

if ~isempty(y)
    for i=1:length(y)
        if i<=length(y)
            a=find(y==y(i));
            if length(a)>1
                for j=length(a):-1:2
                    y(a(j))=[];
                    SC(a(j),:)=[];
                end
            end
        end
    end
end

end
                    
                    
                