function bias = GetBias( W,Cxx,rates)
%GETBIAS Summary of this function goes here
% input:
% W - weight matrix
% Cxx - covariance
% rates- mean firing rates
% ouput:
% bias
%   Detailed explanation goes here
% According to Poisson approximation
mat=0.5*W*Cxx*W';
bias=sqrt(1+(pi*8)*mat)*log(rates./(1-rates))-W*rates;

end

