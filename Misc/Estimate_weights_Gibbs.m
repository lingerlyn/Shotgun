function W_next=Estimate_weights_Gibbs(h0,mu0,sigma0,b,spikes,W_prev)
% input:    
% h_prev,mu_prev,sigma_prev - previous estimate of Slab probability, mean and variance in spike&slab prior of connectivty matrix W
% h_prev,mu_prev,sigma_prev - Prior on Slab probability, mean and variance in spike&slab prior of connectivty matrix W
% b - bias
% spikes - spikes
% W_prev - previous sample of W
use_MH=1; %should we use Metropolis-Hastings algorithm to correct proposal?

% output:
% W - sample from weight matrix

%% initialize
N=size(W_prev,1);
T=size(spikes,2);
spikes_stat=mean(spikes,2);
func=@(x) 1./(1+exp(-x));

%% VB loop - weights
%update weights
B=repmat(b,[1,T]);
spikes_delayed=[spikes_stat spikes(:,1:(end-1))];
W_next=W_prev;
W=W_prev;
rand_indices=randperm(N);

for nn=rand_indices
    W_nn=W; W_nn(:,nn)=0;           

    U=W_nn*spikes_delayed+B;
    f=func(U);
    omega=(spikes-f)*spikes_delayed(nn,:)';
    epsilon=(f-f.^2)*spikes_delayed(nn,:)';
   
    ss=(sigma0(:,nn).^(-2)+epsilon).^(-1); %sigma square
    mu=(ss./(sigma0(:,nn).^2)).*(mu0(:,nn)+omega.*sigma0(:,nn).^2);
    h=h0(:,nn)+0.5*log(ss./(sigma0(:,nn).^2))+0.5*(mu.^2)./ss-0.5*(mu0(:,nn)./sigma0(:,nn)).^2;
        
   slab_occur=(rand(N,1)<func(h));
   slab_value=mu+sqrt(ss).*randn(N,1);
   W_next(:,nn)=slab_occur.*slab_value; %sample from proposal density   
  
   %% Metropolis-Hastings part
   if use_MH
    W_temp=W;
    W_temp(:,nn)=W_next(:,nn);
    log_lik_next=logistic_loglik(spikes,W_temp,B);
    W_temp=W;
    log_lik_prev=logistic_loglik(spikes,W_temp,B);
    log_prior_next=zeros(N,1);
    ind_0=W_next(:,nn)==0;
    log_prior_next(ind_0)=func(-h0(ind_0,nn));
    log_prior_next(~ind_0)=func(h0(~ind_0,nn)).*((2*pi*sigma0(~ind_0,nn)).^(-0.5)).*exp(-0.5*(W_next(~ind_0,nn)-mu0(~ind_0,nn)).^2./sigma0(~ind_0,nn).^2);
    log_prior_prev=zeros(N,1);
    ind_0=W_prev(:,nn)==0;
    log_prior_prev(ind_0)=func(-h0(ind_0,nn));
    log_prior_prev(~ind_0)=func(h0(~ind_0,nn)).*((2*pi*sigma0(~ind_0,nn)).^(-0.5)).*exp(-0.5*(W_prev(~ind_0,nn)-mu0(~ind_0,nn)).^2./sigma0(~ind_0,nn).^2);
    log_proposal_forward=zeros(N,1);
    ind_0=W_next(:,nn)==0;
    log_proposal_forward(ind_0)=func(-h(ind_0));
    log_proposal_forward(~ind_0)=func(h(~ind_0)).*((2*pi*ss(~ind_0)).^(-0.5)).*exp(-0.5*(W_next(~ind_0,nn)-mu(~ind_0)).^2./ss(~ind_0).^2);
    log_proposal_reverse=zeros(N,1);
    ind_0=W_prev(:,nn)==0;
    log_proposal_reverse(ind_0)=func(-h(ind_0));
    log_proposal_reverse(~ind_0)=func(h(~ind_0)).*((2*pi*ss(~ind_0)).^(-0.5)).*exp(-0.5*(W_prev(~ind_0,nn)-mu(~ind_0)).^2./ss(~ind_0).^2);
    
    p_accept=min(1,exp(log_lik_next-log_lik_prev).*log_prior_next.*log_proposal_reverse.*(log_prior_prev.^(-1)).*(log_proposal_forward.^(-1))) ;
    accept=rand(N,1)<p_accept;
    W(accept,nn)=W_next(accept,nn);
   else
       W(:,nn)=W_next(:,nn);
   end
    
end



end

