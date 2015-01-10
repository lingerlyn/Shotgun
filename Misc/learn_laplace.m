function [old_w_nk,w_nk_var] = learn_laplace(Eta,X,w_n,Theta_n,thres_n,n,k)
%LEARN_LAPLACE learns the mean and variance of the laplace approximation of
% likelihood function of W(n,k) for the slab part

% OUTPUT
% old_w_nk --mean of the approximated lieklihood
% w_nk_var -- corresponding variance
[N T] = size(Eta);
old_w_nk = 0;
change = Inf;
i=1;
%Eta(k,randsample(1:T-1,ceil(0.01*T)))=rand(1,ceil(0.01*T))<0.5;
while change>0.01 && i<3    
    w_n(k) = old_w_nk;
    prob = 1./(1+exp(-(Theta_n'*X(:,2:T)+w_n*Eta(:,1:T-1)+thres_n)));
    dl_dw = sum((Eta(n,2:T)-prob).*Eta(k,1:T-1));
    d2l_dw2 = -sum((1-prob).*prob.*((Eta(k,1:T-1).^2)));
    new_w_nk = old_w_nk - dl_dw/d2l_dw2;
    change = abs(new_w_nk-old_w_nk);
    old_w_nk = new_w_nk;
    i=i+1;
end
w_n(k) = old_w_nk;
prob = 1./(1+exp(-(Theta_n'*X(:,2:T)+w_n*Eta(:,1:T-1)+thres_n)));
w_nk_var = 1./sum((1-prob).*prob.*(Eta(k,1:T-1).^2));