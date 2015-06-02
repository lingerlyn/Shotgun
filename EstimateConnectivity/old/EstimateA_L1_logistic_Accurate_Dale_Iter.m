function [EW,idents]=EstimateA_L1_logistic_Accurate_Dale_Iter(Cxx,Cxy,rates,spar,N_stim,pen_diag,warm,idenTol)

    maxIter=20;

    N=length(Cxx);
    idents=zeros(N,1);
    oldIdents=nan(N,1);
    iter=0;

    while (sum(oldIdents~=idents)/N)>idenTol && iter<=maxIter
        
        iter=iter+1;
        oldIdents=idents;
        
        EW=EstimateA_L1_logistic_Accurate_Dale(Cxx,Cxy,rates,spar,N_stim,pen_diag,warm,idents);
        
        idents=getIdentities(EW,idenTol);
        
    end

    
end

