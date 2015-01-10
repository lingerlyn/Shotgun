function W_est= EstimateA_Gibbs( b,spikes,obs_neurons,p_0, mu_0, std_0)
% This function estimates the weight matrix A using gibbs sampling

%% INPUT:
% b : N x 1 vector , Base stimulus for the neurons
% spikes - spiking matrix, N x T
% a = prior zero weight probability
% mu_0 = prior slab mean
% var_0 = prior slab variance
% obs_neurons : N x T binary matrix indicating observed neurons

% iter: Number of sampling iterations to learn W 

%% OUTPUT
% W - N x N sampled connectivity matrix
% W_rb - N x N posterior mean of connectivity matrix


%% Algorithm

% Internal variables
use_old=1; %use old code, with lapalace approximation
L = 1; 
burnin = 30;
iterations=100;


[N,T]=size(spikes);

W=zeros(N); %initial estimate for W
sampled_W=zeros(N,N,iterations);

if ~use_old
    spikes_obs=spikes;
    spikes_obs(~obs_neurons)=0.5;
    spikes=repmat(mean(spikes,2),1,T); %initial guess
    spikes(spikes_obs==1)=1;
    spikes(spikes_obs==0)=0;    
    h_0=log(p_0./(1-p_0)); %in new code we use h instead of p for numerical accuracy
else
    p_0=mean(p_0(:));
    mu_0=mean(mu_0(:));
    std_0=mean(std_0(:));    
    Theta=zeros(1,N); %D x N matrix, the external stimulus filter parameter
    X=zeros(1,T);  % D x T matrix of external stimulus
end

for ii=1:(iterations+burnin)
    if use_old
        % Sample spikes matrix given W
        cell_spikes = GibbsSampleSpikes(obs_neurons, Theta, b, X, L, 0, spikes, W);
        spikes=cell_spikes{1};
        % Sample W given the sampled spike matrix
        [ W, ~ ] = GibbsSampleW(W, Theta, b, spikes, X, p_0, mu_0, std_0);
    else
         W=Estimate_weights_Gibbs(h_0, mu_0, std_0,b,spikes,W);
         [spikes ~]=Estimate_rates_Gibbs_Metropolized(W,b,spikes_obs,spikes);         
    end


sampled_W(:,:,ii)=W;
disp([num2str(100*ii/iterations) '%']);
end

sampled_W(:,:,1:burnin)=[];
W_est=median(sampled_W,3);

end

