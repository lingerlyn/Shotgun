function [R,correlation, zero_matching,sign_matching,TPR_p,FPR_p,TPR_n,FPR_n] = GetWeightsErrors( W,EW )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%     MSE=sqrt(mean((EA(:)-true_A(:)).^2)/mean(true_A(:).^2));    
    R=sqrt(1-mean((EW(:)-W(:)).^2)/mean((W(:)-mean(W(:))).^2));   
    if imag(R)>0
        R=0;
    end
%     SIGN_ERROR=mean(abs(sign(EA(:))-sign(true_A(:))));    
    correlation=corr(W(:),EW(:));
    zero_matching=1-0.5*sum( (EW(:)~=0).*(W(:)==0)+ (EW(:)==0).*(W(:)~=0))/sum(W(:)~=0) ;
    sign_matching=1-0.5*sum( (EW(:)~=0).*(W(:)~=0).*abs(sign(EW(:))-sign(W(:))))/sum((EW(:)~=0).*(W(:)~=0)) ;
    
    N=size(W,1);
    if ismatrix(W)==2
        W_save=W(eye(N)<0.5);
        EW_save=EW(eye(N)<0.5);
    else
        W_save=W;
        EW_save=EW;
    end 
    for sgn=[-1, 1]        
        W=W_save;
        W(sign(W)~=sgn)=0;
        W=abs(W);
        EW=EW_save;
        EW(sign(EW)~=sgn)=0;
        EW=abs(EW);
        positive=abs(W(:))>0;
        prediction=EW(:)>0;
        FP=mean((prediction).*(~positive));
        TP=mean((prediction).*(positive));
        TN=mean((~prediction).*(~positive));
        FN=mean((~prediction).*(positive));
        if sgn==1
            TPR_p=TP./(TP+FN);
            FPR_p=FP./(FP+TN);
        else
            TPR_n=TP./(TP+FN);
            FPR_n=FP./(FP+TN);
        end
    end
end

