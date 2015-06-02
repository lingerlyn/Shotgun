function [EW,iterFlag]=EstimateA_L1_logistic_distdep(CXX,CXY,sparsity,V,f,lambdafunc)
% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1
% linear regression.
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% OUTPUTS:
% EW: ML estimate of weights

%params
iterations=30;
% Tol_sparse=0.1; %tolerance for sparsity level
% Tol_sparse=.05;
Tol_sparse=.01;

%initialize FISTA
L=2*max(eig(CXX)); %lipshitz constant

%initialize binary search
alpha_high=1e2*L; %maximum bound for alpha
% alpha_low=1e-4;  %minimum bound for alpha
% alpha_high=1;
alpha_low=1e-20;
loop_cond=1;  %flag for while llop
iter=0;
maxIter=1000;
iterFlag=0;

sparsity=mean(sparsity);

while  loop_cond %binary search for lambda that give correct sparsity level
    iter=iter+1;
    if iter>maxIter
        EW=x; iterFlag=1; break; 
    end
        
    if sparsity==1
        alpha=0;
    else
        alpha=(alpha_high+alpha_low)/2;
    end
    
    lambda=alpha*lambdafunc(f);


%%% FISTA
t=1;
x=0*CXY';
y=x;

for kk=1:iterations
    x_prev=x;
    u=y-(2/L)*(V*y*CXX-CXY');
    x=ThresholdOperator(u,lambda/L);
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);  
    t=t_next;
end
%%% 


    sparsity_measure=mean(~~x(:));
    cond=sparsity_measure<sparsity;
    loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
    if sparsity==1
        loop_cond=0;
    end

    if cond
        alpha_high=alpha;
    else
        alpha_low=alpha;
    end
end
        
EW=x;
end

function y = ThresholdOperator( x , lambda )
%THRESHOLDOPERATOR Summary of this function goes here
% This function implements the threshold operator for the l1/l2 group norm
% as specified on "convex optimization with sparsity-inducing norms", page 12, last equation

%   Detailed explanation goes here
% x - input
% y - output
% lambda - regularization constant

    temp=(1-lambda./abs(x));
    y=temp.*x;
    y(temp<0)=0;
    



end

