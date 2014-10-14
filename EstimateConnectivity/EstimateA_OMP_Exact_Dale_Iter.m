function [EW,idents]=EstimateA_OMP_Exact_Dale_Iter(Cxx,Cxy,rates,spar,lambda,M,idenTol)

    maxIter=20;

    N=length(Cxx);
    idents=zeros(N,1);
    oldIdents=inf(N,1);
    iter=0;

    while (sum(oldIdents~=idents)/N)>idenTol && iter<=maxIter
        
        iter=iter+1;
        oldIdents=idents;
        
        EW=EstimateA_OMP_Exact_Dale(Cxx,Cxy,rates,spar,lambda,M,idents);
        
        idents=getIdentities(EW,idenTol);
        
    end

    
end

