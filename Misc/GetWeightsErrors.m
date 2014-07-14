function [MSE,correlation,SIGN_ERROR] = GetWeightsErrors( true_A,EA )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    MSE=mean((EA(:)-true_A(:)).^2);    
    SIGN_ERROR=mean(abs(sign(EA(:))-sign(true_A(:))));
    correlation=corr(true_A(:),EA(:));
end

