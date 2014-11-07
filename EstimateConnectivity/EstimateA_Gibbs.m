function [ W, W_rb ] = EstimateA_Gibbs( Theta, thres, Eta, X, a, mu_0, var_0, obs_neurons, W, iter)
% This function estimates the weight matrix A using gibbs sampling

%% INPUT:
% Theta : D x N matrix, the external stimulus filter parameter
% thres : N x 1 vector , Base stimulus for the neurons
% Eta - spiking matrix, N x T
% X : D x T matrix of external stimulus
% a = prior zero weight probability
% mu_0 = prior slab mean
% var_0 = prior slab variance
% obs_neurons : N x T binary matrix indicating observed neurons
% W : initial estimate of W matrix
% iter: Number of sampling iterations to learn W 

%% OUTPUT
% W - N x N sampled connectivity matrix
% W_rb - N x N posterior mean of connectivity matrix


%% Algorithm

% Sample unobserved spikes given W
L = 1; 
burnin = 0;

ii = 0;
while (ii<iter)
[ Eta ] = sample_spikes(obs_neurons, Theta, thres, X, L, burnin, Eta, W);

% Sample W given the sampled spike matrix
[ W, W_rb ] = sample_W(W, Theta, thres, Eta, X, a, mu_0, var_0);
ii=ii+1;
end

end

