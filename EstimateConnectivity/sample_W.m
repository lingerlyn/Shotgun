function [ W, W_rb ] = sample_W( W, Theta, thres, Eta, X, a, mu_0, var_0)
%% SAMPLE_W 
%  Gibbs samples the weight matrix given the spiking data and prior parameters

%% INPUT
% W - previously sampled weight matrix
% Theta, thres (or b) -- other params
% Eta - spiking matrix, N x T
% X - stimulus matrix, P x T
% obs_neurons - N x T observation matrix
% a = prior zero weight probability
% mu_0 = prior slab mean
% var_0 = prior slab variance

%% OUTPUT
% W - N x N sampled connectivity matrix
% W_rb - N x N posterior mean of connctivity matrix (for rao blackwell, we
% are not using it right now)

%% Algorithm

[N T] = size(Eta);
W_rb = zeros(N);

% Sampling each row of W in parallel
parfor i = 1:N
    W_i = W(i,:);
    W_i_rb = W_rb(i,:);
    rand_neuron_order = randperm(N);
    % Sampling each element of row in random order
    for n = rand_neuron_order
        w_prev = W_i(n);
        [w_mean,w_var] = learn_laplace(Eta,X,W_i,Theta(:,i),thres(i),i,n);
        W_temp = W_i;
        W_temp(n) = w_mean;
        c = logistic_lik(Eta,X,W_temp,Theta(:,i),thres(i),i);
        w_var_post = 1./(1./var_0 + 1./w_var);
        w_mean_post = w_var_post*(w_mean/w_var+mu_0/var_0); 
        temp = (w_mean_post.^2 - w_mean.^2)/(0.5*w_var) + (w_mean_post.^2 - mu_0.^2)/(0.5*var_0);
        log_prob_1_prop = log(a)+c+temp+log(sqrt(w_var_post/var_0));
        W_temp = W_i;
        W_temp(n) = 0;
        temp = logistic_lik( Eta,X,W_temp,Theta(:,i),thres(i),i);
        log_prob_0_prop = log(1-a)+(temp);
        prob_1 = exp(log_prob_1_prop-logsumexp([log_prob_1_prop,log_prob_0_prop]));
        % sampling z by metropolized gibbs
        z_eq_0 = rand>prob_1;
        if z_eq_0==1
            prop_w = 0;
        else
            prop_w = normrnd(w_mean_post,sqrt(w_var_post));
        end
        W_temp = W_i;
        W_temp(n) = prop_w;
        log_term1_num = logistic_lik(Eta,X,W_temp,Theta(:,i),thres(i),i);
        W_temp = W_i;
        W_temp(n) = w_prev;
        log_term1_den = logistic_lik(Eta,X,W_temp,Theta(:,i),thres(i),i);
        
        log_term2 = 0.5*((prop_w-w_mean_post).^2/w_var_post - (w_prev-w_mean_post).^2/w_var_post);
        accept_prob = min(1,exp(log_term1_num - log_term1_den+log_term2));
        accept = rand<accept_prob;
        if accept
            W_i(n) = prop_w;
        end
        
        if prop_w~=0 && accept
            W_i_rb(n) = w_mean_post;
        elseif w_prev~=0 && ~accept
            W_i_rb(n) = w_prev;
        end
        
    end
    W(i,:) = W_i;
    W_rb(i,:) = W_i_rb;
end
end

