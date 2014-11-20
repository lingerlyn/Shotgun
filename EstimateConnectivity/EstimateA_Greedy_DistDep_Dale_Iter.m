function [EW,idents]=EstimateA_Greedy_DistDep_Dale_Iter(Cxx,Cxy,rates,spar,lambda,M,PC,lambda2,idenTol)

    maxIter=5;

    N=length(Cxx);
    idents=zeros(N,1);
    oldIdents=inf(N,1);
    iter=0;

    while (sum(oldIdents~=idents)/N)>idenTol && iter<=maxIter
        
        iter=iter+1;
        oldIdents=idents;
        
        EW=EstimateA_Greedy_DistDep_Dale(Cxx,Cxy,rates,spar,lambda,M,PC,lambda2,idents);
        
        idents=getIdentities(EW,idenTol);
        
    end

    
end

