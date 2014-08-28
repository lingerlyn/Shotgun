function bias = GetBias2( EW,Cxx,rates)
% input:
% W - weight matrix
% Cxx - covariance
% rates- mean firing rates
% ouput:
% bias

%   Detailed explanation goes here
% According to linearized approximation of logistic

N=length(rates);
bias=zeros(N,1);

% options = optimset('Display','off');
options = optimset('GradObj','on','Display','off');

for kk=1:N
  C=EW(kk,:)*Cxx*EW(kk,:)';
  [bias(kk),~,exitflag]=fminunc(@(x) bias_func(x,rates(kk),C),-1,options)
end

end

function [objective,grad] = bias_func( x,C,rate)
%GETBIAS Summary of this function goes here
% %   Detailed explanation goes here


objective=-rate*x+log(1+exp(x))+0.5*C./((1+exp(-x)).*(1+exp(x)));
grad=-rate+1./(1+exp(-x))-0.5*C*(exp(x)-exp(-x))./((1+exp(-x)).*(1+exp(x))).^2;

% objective=-rate*x+log(1+exp(x));
% grad=-rate+1./(1+exp(-x));

end

