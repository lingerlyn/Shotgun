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
bias=log(rates)-W*rates-mat(eye(size(mat))>0.5);

end

