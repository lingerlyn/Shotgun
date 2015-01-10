function [ log_lik ] = logistic_lik( Eta,X,W_n,Theta_n,thres,n)
%LOGISTIC_LIK Returns likelihood of Eta(n,:) for a given parameter
[N T] = size(Eta);
prob = 1./(1+exp(-(Theta_n'*X(:,2:T)+W_n*Eta(:,1:T-1)+thres)));
log_lik = sum(log(prob).*Eta(n,2:T)+log(1-prob).*(1-Eta(n,2:T)));


end

