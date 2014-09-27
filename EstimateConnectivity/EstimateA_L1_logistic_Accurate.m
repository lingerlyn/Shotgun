function EW=EstimateA_L1_logistic_Accurate(CXX,CXY,rates,sparsity,N_stim,pen_diag,warm)
% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% initial - initial conditions
% OUTPUTS:
% EW: ML estimate of weights

%params
% iterations=300;

N=length(rates);
Tol_sparse=0.01; %tolerance for sparsity level
Tol_FISTA=1/N^2; %toleratnce for fista
max_iterations=100;

%initialize FISTA
% if initial
%     x=initial;
% else
    x=0*CXY';
% end
y=x;


V=-diag(pi*(rates.*log(rates)+(1-rates).*log(1-rates))/8); %for the gradient
V(isnan(V))=0; %take care of 0*log(0) cases...
L=2*max(V(:))*max(eig(CXX)); %lipshitz constant

%initialize binary search
lambda_high=max(1e8,1e2*L); %maximum bound for lambda
lambda_low=min(1e-8,1e-2*L);  %minimum bound for lambda
loop_cond=1;  %flag for while llop

while  loop_cond %binary search for lambda that give correct sparsity level
    if sparsity==1
        lambda=0;
    else
        lambda=(lambda_high+lambda_low)/2;
    end

%%% FISTA
if ~warm
    x=0*CXY';
    y=x;
end
    
if pen_diag
    mask=ones(N);
else  
    mask=ones(N)-eye(N); % used to remove L1 penalty from diagonal
end
t_next=1;
FISTA_cond=1;
iteration=0;

while FISTA_cond
    t=t_next;
    x_prev=x;
    
    A=1./sqrt(1+pi*(y*CXX*y')/8);
    B=diag(A(eye(N)>0.5));
    u=y-(2/L)*(V*B*y*CXX-CXY');
      
    x=ThresholdOperator(u,mask.*lambda/L);
    if any(~isfinite(u(:)))
        error('non finite x!')
    end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);   
    iteration=iteration+1;    
    FISTA_cond=(mean(abs(x(:)-x_prev(:)))>Tol_FISTA);
    if (iteration>max_iterations);
        break
        warning('FISTA max iteration reached');
    end
end
%%% 

    temp=~~x(1:(end-N_stim),1:(end-N_stim));
    sparsity_measure=mean(temp(:));
    cond=sparsity_measure<sparsity;

    
    loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
    if sparsity==1
        loop_cond=0;
    end

    if cond
        lambda_high=lambda
    else
        lambda_low=lambda
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
    temp(isnan(temp))=1;
    y=temp.*x;
    y(temp<0)=0;
    



end

