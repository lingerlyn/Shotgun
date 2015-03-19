function [ W, W_rb,accept_ratio ] = GibbsSampleW( W, Theta, thres, Eta, X, p_0, mu_0, std_0)
%% SAMPLE_W 
%  Gibbs samples the weight matrix given the spiking data and prior parameters

%% INPUT
% W - previously sampled weight matrix
% Theta, thres (or b) -- other params
% Eta - spiking matrix, N x T
% X - stimulus matrix, P x T
% obs_neurons - N x T observation matrix
% p_0 = prior non-zero weight probability
% mu_0 = prior slab mean
% std_0 = prior slab variance

%% OUTPUT
% W - N x N sampled connectivity matrix
% W_rb - N x N posterior mean of connctivity matrix (for rao blackwell, we
% are not using it right now)
% accept_ratio - fraction of samples accepts in MH scheme

%% Algorithm

[N, T] = size(Eta);
W_rb = zeros(N);
accept_ratio=0;

% Sampling each row of W in parallel
for ii = 1:N
    W_i = W(ii,:);
    W_i_rb = W_rb(ii,:);
    rand_neuron_order = randperm(N);
    % Sampling each element of row in random order
    for nn = rand_neuron_order
        w_prev = W_i(nn);
        [w_mean,w_var] = learn_laplace(Eta,X,W_i,Theta(:,ii),thres(ii),ii,nn);
        W_temp = W_i;
        W_temp(nn) = w_mean;
        c = logistic_lik(Eta,X,W_temp,Theta(:,ii),thres(ii),ii);
        w_var_post = 1./(1./std_0(ii,nn) + 1./w_var);
        w_mean_post = w_var_post*(w_mean/w_var+mu_0(ii,nn)/std_0(ii,nn)); 
        temp = (w_mean_post.^2 - w_mean.^2)/(0.5*w_var) + (w_mean_post.^2 - mu_0(ii,nn).^2)/(0.5*std_0(ii,nn));
        log_prob_1_prop = log(p_0(ii,nn))+c+temp+log(sqrt(w_var_post/std_0(ii,nn)));
        W_temp = W_i;
        W_temp(nn) = 0;
        temp = logistic_lik( Eta,X,W_temp,Theta(:,ii),thres(ii),ii);
        log_prob_0_prop = log(1-p_0(ii,nn))+(temp);
        prob_1 = exp(log_prob_1_prop-logsumexp([log_prob_1_prop,log_prob_0_prop]));
        % sampling z by metropolized gibbs
        z_eq_0 = rand>prob_1;
        if z_eq_0==1
            prop_w = 0;
        else
            prop_w = normrnd(w_mean_post,sqrt(w_var_post));
        end
        W_temp = W_i;
        W_temp(nn) = prop_w;
        log_term1_num = logistic_lik(Eta,X,W_temp,Theta(:,ii),thres(ii),ii);
        W_temp = W_i;
        W_temp(nn) = w_prev;
        log_term1_den = logistic_lik(Eta,X,W_temp,Theta(:,ii),thres(ii),ii);
        
        log_term2 = 0.5*((prop_w-w_mean_post).^2/w_var_post - (w_prev-w_mean_post).^2/w_var_post);
        accept_prob = min(1,exp(log_term1_num - log_term1_den+log_term2));
        accept = rand<accept_prob;
        if accept
            W_i(nn) = prop_w;
        end
        
        if prop_w~=0 && accept
            W_i_rb(nn) = w_mean_post;
        elseif w_prev~=0 && ~accept
            W_i_rb(nn) = w_prev;
        end
        
        accept_ratio=accept_ratio+accept/N^2;
    end
    W(ii,:) = W_i;
    W_rb(ii,:) = W_i_rb;    
end
end

