function [W_est, accept_array]= EstimateA_Gibbs( b,spikes,obs_neurons,p_0, mu_0, std_0,rates)
% This function estimates the weight matrix A using gibbs sampling

%% INPUT:
% b : N x 1 vector , Base stimulus for the neurons
% spikes - spiking matrix, N x T
% a = prior zero weight probability
% mu_0 = prior slab mean
% var_0 = prior slab variance
% rates: mean firing rates, measure empirically
% obs_neurons : N x T binary matrix indicating observed neurons

% iter: Number of sampling iterations to learn W 

%% OUTPUT
% W - N x N sampled connectivity matrix
% W_rb - N x N posterior mean of connectivity 
% accept_array - array of accept ratios


%% Algorithm

% Internal parameters
use_old=1; %use old code, with lapalace approximation
L = 1; 
burnin = 10;
iterations=10;

[N,T]=size(spikes);
spikes=full(spikes);

%initialization
W=zeros(N); 
sampled_W=zeros(N,N,iterations);
temp=repmat(rates,1,T)>rand(N,T); 
spikes(obs_neurons==0)=temp(obs_neurons==0);
spikes_obs=spikes;
spikes_obs(~obs_neurons)=0.5;

if ~use_old
    h_0=log(p_0./(1-p_0)); %in new code we use h instead of p for numerical accuracy
else
    Theta=zeros(1,N); %D x N matrix, the external stimulus filter parameter
    X=zeros(1,T);  % D x T matrix of external stimulus
end

accept_array=[];
inner_reps=1;

for ii=1:(iterations+burnin)
    if use_old
        % Sample spikes matrix given W
%         if any(obs_neurons(:)==0)
%             cell_spikes = GibbsSampleSpikes(obs_neurons, Theta, b, X, L, 0, spikes, W);
%             spikes=cell_spikes{1};
%         end
         if any(obs_neurons(:)==0)
            [spikes, ~]=Estimate_rates_Gibbs_Metropolized(W,b,spikes_obs,spikes);         
         end
        % Sample W given the sampled spike matrix
        [ W, ~ ,accept_ratio] = GibbsSampleW(W, Theta, b, spikes, X, p_0, mu_0, std_0);
    else
         if any(obs_neurons(:)==0)
            [spikes, ~]=Estimate_rates_Gibbs_Metropolized(W,b,spikes_obs,spikes);         
         end
         for kk=1:inner_reps
            [W,accept_ratio]=Estimate_weights_Gibbs(h_0, mu_0, std_0,b,spikes,W);
         end
    end

sampled_W(:,:,ii)=W;
disp([num2str(100*ii/(iterations+burnin)) '%']);
accept_array(end+1)=accept_ratio; %#ok

end

sampled_W(:,:,1:burnin)=[];
W_est=W*0;
for ii=1:N
    for jj=1:N
        temp=sampled_W(ii,jj,:);
        W_est(ii,jj)=mean(temp(~~temp));
        if mean(~temp)>0.5
            W_est(ii,jj)=0;
        end
    end
end
end

