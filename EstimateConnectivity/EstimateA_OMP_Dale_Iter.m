function [EW,idents]=EstimateA_OMP_Dale_Iter(Cxx,Cxy,rates,spar,lambda,M,idenTol)

    maxIter=20;

    N=length(Cxx);
    idents=zeros(N,1);
    oldIdents=inf(N,1);
    iter=0;

    while (sum(oldIdents~=idents)/N)>idenTol && iter<=maxIter
        
        iter=iter+1;
        oldIdents=idents;
        
        EW=EstimateA_OMP_Dale(Cxx,Cxy,spar,lambda,M,rates,idents);
        
        idents=getIdentities(EW,idenTol);
        
    end

    
end

