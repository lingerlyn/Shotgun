function [EW, Eb]=EstimateA_L1_logistic_fullyobserved(CXX,CXY,rates,spikes,sparsity,N_stim,pen_diag,warm)
% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% OUTPUTS:
% EW: MAP estimate of weights
% Eb: MAP estimate of weights

%params
max_iterations=5000;
Tol_sparse=0.02; %tolerance for sparsity level
N=length(rates);
T=size(spikes,2);
Tol_FISTA=1e-6; %toleratnce for fista

%initialize FISTA
x=0*CXY';
y=x;
L=10*2*max(rates)*max(eig(CXX)); %lipshitz constant - and educated guess. make sure this is large enough to ensure convergence

%initialize binary search
lambda_high=1e-2; %max(1e8,1e8*L); %maximum bound for lambda
lambda_low=1e-4;%min(1e-8,1e-8*L);  %minimum bound for lambda
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
Eb=0;
z=0;
MXY=(spikes(:,1:end-1)*spikes(:,2:end)')';
FISTA_cond=1;
iteration=0;

while FISTA_cond
    
    t=t_next;
    x_prev=x;
    Eb_prev=Eb;
    f=1./(1+exp(-bsxfun(@plus,y*spikes(:,1:end-1),z)));
    u=y-(2/L)*( f*(spikes(:,1:end-1))'-MXY)/(T-1);
    Eb=z-(2/L)*(mean(f,2)-rates);
      
    x=ThresholdOperator(u,mask.*lambda/L);
    if any(~isfinite(u(:)))
        error('non finite x!')
    end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);    
    z=Eb+((t-1)/t_next)*(Eb-Eb_prev);   
    FISTA_cond=(mean(abs(x(:)-x_prev(:)))>Tol_FISTA);
    iteration=iteration+1;
    
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

    sparsity_measure
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

