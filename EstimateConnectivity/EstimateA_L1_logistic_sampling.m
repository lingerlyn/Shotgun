function [EW, Eb]=EstimateA_L1_logistic_sampling(CXX,CXY,rates,sparsity,N_stim,pen_diag,warm,is_spikes)
% Following derivation in GLM_GradUsingSampling 13.1.2015.lyx

% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% is_spikes - if equal 1 we use a cavity-style averages
% OUTPUTS:
% EW: MAP estimate of weights
% Eb: MAP estimate of weights

%params
Tol_sparse=0.1; %tolerance for sparsity level
N=length(rates);
Tol_FISTA=1/N^2; %toleratnce for fista
samples_num=1e3;
max_iterations=1e3;

%initialize FISTA
x=0*CXY';
y=x;
Eb=0;
z=0;
L=20*max(eig(CXX)); %lipshitz constant - and educated guess. make sure this is large enough to ensure convergence


%initialize binary search
lambda_high=max(1e4,1e1*L); %maximum bound for lambda
lambda_low=min(1e-8,1e-8*L);  %minimum bound for lambda
% lambda_high=1e-1;
% lambda_low=1e-5;
loop_cond=1;  %flag for while llop
cond=1;

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
    Eb=0;
    z=0;
end
    
if pen_diag
    mask=ones(N);
else  
    mask=ones(N)-eye(N); % used to remove L1 penalty from diagonal
end
t_next=1;

MYX=CXY'+rates*rates';
Ef_cavity=zeros(N,N);

FISTA_cond=1;
iteration=0;

while FISTA_cond
    t=t_next;
    x_prev=x;
    Eb_prev=Eb;
    samples=randn(N,samples_num);
    spikes=bsxfun(@plus,rates,chol(CXX)*samples);
    Ef=mean(1./(1+exp(-bsxfun(@plus,y*spikes,z))),2);
    if is_spikes==1
        for nn=1:N        
            rates_cavity=rates+CXX(:,nn)*(1-rates(nn))/CXX(nn,nn);        
            CXX_cavity=CXX-CXX(:,nn)*CXX(nn,:)/CXX(nn,nn);
            ind=1:N;
            ind(nn)=[];
            spikes(ind,:)=bsxfun(@plus,rates_cavity(ind),chol(CXX_cavity(ind,ind))*samples(ind,:));
            spikes(nn,:)=1;
            Ef_cavity(:,nn)=mean(1./(1+exp(-bsxfun(@plus,y*spikes,z))),2);
            u=y-(2/L)*( Ef_cavity*diag(rates)-MYX);
        end
    else
        Efspikes=(1./(1+exp(-bsxfun(@plus,y*spikes,z))))*spikes'/size(spikes,2);
        u=y-(2/L)*(Efspikes-MYX);
    end    
    Eb=z-(2/L)*(Ef-rates);
      
    x=ThresholdOperator(u,mask.*lambda/L);
    if any(~isfinite(u(:)))
        error('non finite x!')
    end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);    
    z=Eb+((t-1)/t_next)*(Eb-Eb_prev); 
    
    iteration=iteration+1;    
    FISTA_cond=((mean(abs(x(:)-x_prev(:)))>Tol_FISTA))||(iteration<30);
    if nnz(x(:))==0
        break        
    end
    if (iteration>max_iterations);
        break
        warning('FISTA max iteration reached');
    end
end
%%% 

    temp=~~x(1:(end-N_stim),1:(end-N_stim));
    sparsity_measure=mean(temp(:));
    cond_old=cond;
    cond=sparsity_measure<sparsity;
    % decrease noise and step size near convergence point
%     if cond_old~=cond
%         samples_num=round(samples_num*1.1);
%         L=L*2;
%         lambda=lambda*2;
%         lambda_high=2*lambda_high;
%         lambda_low=2*lambda_low;
%     end
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

