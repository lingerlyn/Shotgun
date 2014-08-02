function EW=EstimateA_L1_proj(CXX,CXY,S,sparsity,l2,Pinit,xinit)
% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% S - matrix of slab means
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% l2 - projected L2 penalty term
% Pinit - initial projection matrix
% xinit - initial starting point
% OUTPUTS:
% EW: ML estimate of weights

%params
iterations=30;
Tol_sparse=0.1; %tolerance for sparsity level
maxIter=500;

%initialize FISTA
if isempty(xinit)
    x=0*CXY';
    y=x;
else
    x=xinit;
    y=x;
end
L=2*max(eig(CXX))+l2; %lipshitz constant. compensate for the added term

%initialize binary search
lambda_high=1e2*L; %maximum bound for lambda
lambda_low=1e-4;  %minimum bound for lambda
loop_cond=1;  %flag for while llop

iter=0;
while  loop_cond %binary search for lambda that give correct sparsity level
    iter=iter+1;
    if iter>maxIter; EW=x; break; end
    if sparsity==1
        lambda=0;
    else
        lambda=(lambda_high+lambda_low)/2;
    end

%%% FISTA
t=1;

for kk=1:iterations
    if kk==1
        P=Pinit;
    else
        P=(~~x)';
    end
        
    x_prev=x;
%     u=y-(2/L)*(y*CXX-CXY');
    u=y-(2/L)*(y*(CXX+diag(l2*P))-(CXY+l2*S.*P)');
    x=ThresholdOperator(u,lambda/L);
    if any(isnan(x))|| any(isinf(x)); keyboard; end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);    
end
%%% 


    sparsity_measure=mean(~~x(:));
    cond=sparsity_measure<sparsity;
    loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
    if sparsity==1
        loop_cond=0;
    end

    if cond
        lambda_high=lambda;
    else
        lambda_low=lambda;
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

