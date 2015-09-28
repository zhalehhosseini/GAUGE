function GCM=GCA(model)
% calculating gene couplings for every pair of model genes. 
%       INPUT:
%        model: metabolic model in COBRA format with known gene_reaction
%        associations in model.rules field
%       OUTPUT:
%        GCM: the resulting gene coupling table
%       Interpretation for element (i, j):
%           0 - uncoupled
%           1 - fully coupled
%           2 - gene i is directionally coupled to j
%           4 - gene j is directionally coupled to i

GCM=-ones(length(model.genes),length(model.genes));
for i=1:length(model.rxns)
    if model.rev(i)
        model.lb(i)=-1000;
        model.ub(i)=1000;
    else
        model.lb(i)=0;
        model.ub(i)=1000;
    end
end
for i=1:length(model.genes)
    
    expressed_rxns=find(model.rxnGeneMat(:,i));
    
    for j=1:length(model.genes)
        if j==i
            continue
        end
        expressed_rxns_2=expressed_rxns;
        x=ones(length(model.genes),1);
        x(j)=0;
        deleted=-ones(length(model.rxns));
        for k=1:length(model.rxns)
            a=strcmp(model.rules(k),'');
            if a
                deleted(k)=-1;
            else
            deleted(k)=eval(model.rules{k});
            end
        end
    
        deleted_rxns=find(deleted==0);
        
        model2=model;
        model2.lb(deleted_rxns)=0;
        model2.ub(deleted_rxns)=0;
        all_del_rxns=[];
        for k=1:length(expressed_rxns_2)
            model2.c(:)=0;
            model2.c(expressed_rxns_2(k))=1;
            [sol]=optimizeCbModel(model2,'max');
            max=sol.f;
            [sol]=optimizeCbModel(model2,'min');
            min=sol.f;
            if abs(min)<10e-7 && abs(max)<10e-7
                all_del_rxns=[all_del_rxns;k];
            end
        end
        
        expressed_rxns_2(all_del_rxns)=[];
        if isempty(expressed_rxns_2)
            GCM(i,j)=1;
        end
    end
end
GCM(1:length(model.genes)+1:length(model.genes)*length(model.genes))=1;
for i=1:length(model.genes)
    for j=1:length(model.genes)
        if GCM(i,j)==1 && GCM(j,i)==-1
            GCM(i,j)=2;GCM(j,i)=3;
        end
        if GCM(i,j)==-1 && GCM(j,i)==1
            GCM(i,j)=3;GCM(j,i)=2;
        end
        if GCM(i,j)==-1 && GCM(j,i)==-1
            GCM(i,j)=0;GCM(j,i)=0;
        end
    end
end

end