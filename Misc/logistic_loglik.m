function [ log_lik ] = logistic_loglik( Eta,W,B)
%LOGISTIC_LIK Returns likelihood of Eta(n,:) for a given parameter
prob = 1./(1+exp(-(W*Eta(:,1:end-1)+B(:,1:end-1))));
log_lik = sum(log(prob).*Eta(:,2:end)+log(1-prob).*(1-Eta(:,2:end)),2);


end

