function [MSE,correlation, sparsity_error] = GetWeightsErrors( true_A,EA )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    MSE=sqrt(mean((EA(:)-true_A(:)).^2)/mean(true_A(:).^2));    
%     SIGN_ERROR=mean(abs(sign(EA(:))-sign(true_A(:))));
    sparsity_error=mean( (EA(:)~=0).*(true_A(:)==0)+ (EA(:)==0).*(true_A(:)~=0)) ;
    correlation=corr(true_A(:),EA(:));
end

