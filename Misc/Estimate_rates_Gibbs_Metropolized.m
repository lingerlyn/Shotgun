function [eta_out, beta]=Estimate_rates_Gibbs_Metropolized(W,b,eta_obs,eta_sample)
% input:    
% W - connectivty matrix
% b - bias
% eta_obs - observed spikes
% eta_sample - previous sample

% output:
% eta -sampled spikes
% beta -  predicted average firing rates
%parameters
Gibbs_Reps_start=30; %sample to start from
Gibbs_Reps=Gibbs_Reps_start+1;
disp(['Sampling unobserved spikes...'])

%initialize
N=size(eta_sample,1);
T=size(eta_sample,2); %we assume T is even
if mod(T,2)
    error('T must be even!!');
end
beta_stat=mean(eta_sample,2);
beta=eta_sample*0;

%Gibbs loop
f=@(x) (1+exp(-x)).^(-1);
g=@(x) log(f(x));
eta_odd=eta_sample(:,1:2:end);
eta_even=eta_sample(:,2:2:end);
eta_obs_odd=eta_obs(:,1:2:end);
eta_obs_even=eta_obs(:,2:2:end);
eta_out=zeros(N,T);
W_ext=repmat(W,[1 1 T/2]);
B=repmat(b,1,T/2);

%Maybe metropolized version will be better, but what is the proposal density?
for rep=1:Gibbs_Reps
    if ~mod(rep,floor(Gibbs_Reps/10))
        disp([num2str(100*rep/Gibbs_Reps,2) '%'])
    end
    eta_stat=rand(N,1)<beta_stat;
    alpha_odd1=W*[eta_stat eta_even(:,1:(end-1))]+B-W.'*(1-eta_even); 
    beta_odd=alpha_odd1*0;
    
    %odd times
    for nn=1:N
        W_nn=W; W_nn(:,nn)=0;
        temp=W_nn*eta_odd+B;
        alpha_odd2=sum(g(temp+squeeze(W_ext(:,nn,:)))-g(temp),1);
        
        beta_odd(nn,:)=f(alpha_odd1(nn,:)+alpha_odd2);
        eta_odd(nn,:)=rand(1,T/2)<beta_odd(nn,:);
        eta_odd(nn,eta_obs_odd(nn,:)==1)=1;
        eta_odd(nn,eta_obs_odd(nn,:)==0)=0;
    end
    
    %even times
    eta_stat=rand(N,1)<beta_stat;
    alpha_even1=W*eta_odd+B-W.'*(1-[eta_odd(:,2:end) eta_stat]);
    beta_even=alpha_even1*0;    
    
    for nn=1:N
        W_nn=W; W_nn(:,nn)=0;
        temp=W_nn*eta_even+B;
        alpha_even2=sum(g(temp+squeeze(W_ext(:,nn,:)))-g(temp),1);
        
        beta_even(nn,:)=f(alpha_even1(nn,:)+alpha_even2);
        eta_even(nn,:)=rand(1,T/2)<beta_even(nn,:);        
        eta_even(nn,eta_obs_even(nn,:)==1)=1;
        eta_even(nn,eta_obs_even(nn,:)==0)=0;
    end
    
    if rep>Gibbs_Reps_start
        beta(:,1:2:end)=beta_odd+beta(:,1:2:end);
        beta(:,2:2:end)=beta_even+beta(:,2:2:end);
        beta(eta_obs==0)=0;
        beta(eta_obs==1)=1;
    end
   
end

beta=beta/(Gibbs_Reps-Gibbs_Reps_start);
beta(eta_obs==1)=1;
beta(eta_obs==0)=0;
eta_out(:,1:2:end)=eta_odd;
eta_out(:,2:2:end)=eta_even;

end
