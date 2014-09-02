function [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( true_A,EA )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%     MSE=sqrt(mean((EA(:)-true_A(:)).^2)/mean(true_A(:).^2));    
    R=sqrt(1-mean((EA(:)-true_A(:)).^2)/mean((true_A(:)-mean(true_A(:))).^2));   
    if imag(R)>0
        R=0;
    end
%     SIGN_ERROR=mean(abs(sign(EA(:))-sign(true_A(:))));    
    correlation=corr(true_A(:),EA(:));
    zero_matching=1-0.5*sum( (EA(:)~=0).*(true_A(:)==0)+ (EA(:)==0).*(true_A(:)~=0))/sum(true_A(:)~=0) ;
    sign_matching=1-0.5*sum( (EA(:)~=0).*(true_A(:)~=0).*abs(sign(EA(:))-sign(true_A(:))))/sum((EA(:)~=0).*(true_A(:)~=0)) ;
end

