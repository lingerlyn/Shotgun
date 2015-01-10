function [ Etas ] = GibbsSampleSpikes(obs_neurons, Theta, thres, X, L, burnin, initial_Eta, W)
%% GIBBS_LOGISTIC 
%  Gibbs samples the unobserved spikes given current estimates of parametes

%% INPUT:
% obs_neurons : N x T binary matrix indicating observed neurons
% Theta : D x N matrix, the external stimulus filter parameter
% thres : N x 1 vector , Base stimulus for the neurons
% W : N x N weight matrix
% X : D x T matrix of external stimulus
% L : Number of samples of Eta you want to draw
% burnin: Burn-in, discards initial burnin samples
% initial_Eta : N x T matrix, Initial sample of Eta

%% OUTPUT
% Etas : The cell of length L, each containing a sample of Eta of size N x T
%        after neglecting the burn-in samples

%% Algorithm

[N T] = size(obs_neurons);
Etas = {};
last_Eta = initial_Eta;
non_zero_weights = W~=0;

for samples = 1:L
    
    Eta = logical(last_Eta);
    
    %% Updating eta in parallel for even time points
    Eta_even_sample = zeros(N,T/2);
    Eta_shift_back = Eta(:,1:2:T);
    Eta_shift_for = Eta(:,3:2:T);
    Eta_shift_for(:,T/2) = 0;
    Eta_even = Eta(:,2:2:T);
    X_even = X(:,2:2:T);
    X_odd = X(:,3:2:T);
    X_odd(:,T/2) = 0;
    obs_neurons_even = obs_neurons(:,2:2:T);
    parfor t = 1:(T/2)
        obs_neurons_t = obs_neurons_even(:,t);
        Eta_t = Eta_even(:,t); % we want to update this
        for n = 1:N
            % sample only for unobserved neurons
            if obs_neurons_t(n)==0
                conn_neurons = non_zero_weights(:,n); % neurons taking input from neuron n
                
                if Eta_t(n)==1
                    % log of prior probability for the flipped sign
                    term_1_1 = -log((1+exp((Theta(:,n)'*X_even(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    % log of prior probability for the same sign
                    term_0_1 = -log((1+exp(-(Theta(:,n)'*X_even(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                else
                    % log of prior probability for the flipped sign
                    term_1_1 = -log((1+exp(-(Theta(:,n)'*X_even(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    % log of prior probability for the same sign
                    term_0_1 = -log((1+exp((Theta(:,n)'*X_even(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                end
                
                % for neurons at time T, we don't have a likelihood term
                % calcating probability of accepting flipped sign
                if t~=(T/2)
                    % calculate the log likelihood of observing neurons at
                    % next time
                    if ~isempty(conn_neurons)
                        eta_temp = Eta_t;
                        eta_temp(n) = 1-Eta_t(n);
                        phi = 1./(1+exp(-(Theta(:,conn_neurons)'*X_odd(:,t)+W(conn_neurons,:)*eta_temp+thres(conn_neurons))));
                        temp=(Eta_shift_for(conn_neurons,t));
                        term_1_2 = sum(log(phi(temp)));
                        term_1_3 = sum(log((1-phi(~temp))));
                    else
                        term_1_2 = 0;
                        term_1_3 = 0;
                    end
                    % log probability of the flipped sign (upto
                    % proportionality)
                    log_prob_1_prop = term_1_1+term_1_2+term_1_3;
                else
                    log_prob_1_prop = term_1_1;
                end
                
                % calcating probability of not flipped sign
                if t~=(T/2)
                    % calculate the log likelihood of observing neurons at
                    % next time
                    if ~isempty(conn_neurons)
                        % calculating likelihood by taking into account
                        % only the neurons receiving input from neuron n
                        eta_temp = Eta_t; % no flip of sign
                        % probablities of positive spike for neurons at
                        % next time conditional on eta_temp
                        phi = 1./(1+exp(-(Theta(:,conn_neurons)'*X_odd(:,t)+W(conn_neurons,:)*eta_temp+thres(conn_neurons))));
                        %phi(find(phi==1)) = 0.99999;
                        %phi(find(phi==0)) = 0.00001;
                        % likelihood terms without flipped sign
                        temp=(Eta_shift_for(conn_neurons,t));
                        term_0_2 = sum(log(phi(temp)));
                        term_0_3 = sum(log((1-phi(~temp))));
                        %term_0_2 = sum(log(phi.^Eta_shift_for(conn_neurons,t)));
                        %term_0_3 = sum(log((1-phi).^(1-Eta_shift_for(conn_neurons,t))));
                    else
                        term_0_2 = 0;
                        term_0_3 =0;
                    end
                    % log probability of not flipped sign (upto
                    % proportionality)
                    log_prob_0_prop = term_0_1+term_0_2+term_0_3;
                else
                    log_prob_0_prop = term_0_1;
                end
                % calculate acceptance probability of flipped sign
                %acc_prob = min(1,exp(log_prob_1_prop-log_prob_0_prop));
                %accept = binornd(1,acc_prob);
                if log_prob_1_prop > log_prob_0_prop
                    Eta_t(n)= (1-Eta_t(n)); %accept=1; 
                else
                    accept = (rand<exp(log_prob_1_prop-log_prob_0_prop));
                    if accept
                        Eta_t(n)= (1-Eta_t(n));
                    end 
                end;
                %Eta_t(n)= (1-Eta_t(n))*accept+Eta_t(n)*(1-accept);

            end
        end
        Eta_even_sample(:,t) = Eta_t;
    end
    Eta(:,2:2:T) = Eta_even_sample;
    %% update eta in parallel for odd times
    
    
    % since we can't refer to t+c while parforing over t, we store the
    % matrix to used for the odd times
    Eta_odd_sample = zeros(N,T/2);
    obs_neurons_odd = obs_neurons(:,1:2:T);
    Eta_odd = Eta(:,1:2:T);
    Eta_shift_for = Eta(:,2:2:T);
    Eta_shift_back = zeros(N,T/2);
    Eta_shift_back(:,2:end) = Eta(:,2:2:(T-2));
    X_odd = X(:,1:2:T);
    X_even = X(:,2:2:T);
    
    parfor t = 1:(T/2)
        obs_neurons_t = obs_neurons_odd(:,t);
        % we want to update this
        Eta_t = Eta_odd(:,t);
        for n = 1:N
            % sample only for unobserved neurons
            if obs_neurons_t(n)==0
                conn_neurons = non_zero_weights(:,n); % neurons taking input from neuron n
                
                % calcating prior probability of flipped sign
                if ~isempty(conn_neurons)
                    eta_temp = Eta_t;
                    eta_temp(n) = 1-Eta_t(n); % eta after flipping the sign
                    % probability of positive spike at next time
                    % conditional on eta_temp
                    phi = 1./(1+exp(-(Theta(:,conn_neurons)'*X_even(:,t)+W(conn_neurons,:)*eta_temp+thres(conn_neurons))));
                    %phi(find(phi==1)) = 0.99999;
                    %phi(find(phi==0)) = 0.00001;
                    % likelihood terms for observed spke at next time
  %                  term_1_2 = sum(log(phi.^Eta_shift_for(conn_neurons,t)));
   %                 term_1_3 = sum(log((1-phi).^(1-Eta_shift_for(conn_neurons,t))));
                       temp=(Eta_shift_for(conn_neurons,t));
                        term_1_2 = sum(log(phi(temp)));
                        term_1_3 = sum(log((1-phi(~temp))));
                 else
                    term_1_2 = 0;
                    term_1_3 = 0;
                end
                % prior probability has different expressions for t = 1 and
                % rest
                if t~=1
                    % log of prior probability for the flipped sign
                    if Eta_t(n)==1
                        term_1_1 = -log((1+exp((Theta(:,n)'*X_odd(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    else
                        term_1_1 = -log((1+exp(-(Theta(:,n)'*X_odd(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    end
                    % log probability of accepting the flipped sign (upto
                    % proportionality)
                    log_prob_1_prop = term_1_1+term_1_2+term_1_3;
                else
                    % log of prior probability for the flipped sign
                    if Eta_t(n)==1
                        term_1_1 = -log((1+exp((Theta(:,n)'*X_odd(:,t)+thres(n)))));
                    else
                        term_1_1 = -log((1+exp(-(Theta(:,n)'*X_odd(:,t)+thres(n)))));
                    end
                    % log probability of accepting the flipped sign (upto
                    % proportionality)
                    log_prob_1_prop = term_1_1+term_1_2+term_1_3;
                end
                
                
                
                % calcating probability of not flipped sign
                if ~isempty(conn_neurons)
                    eta_temp = Eta_t; % no sign flip
                    % probability of positive spike at next time
                    % conditional on eta_temp
                    phi = 1./(1+exp(-(Theta(:,conn_neurons)'*X_even(:,t)+W(conn_neurons,:)*eta_temp+thres(conn_neurons))));
                    %phi(find(phi==1)) = 0.99999;
                    %phi(find(phi==0)) = 0.00001;
                    % likelihood terms for observed spke at next time
%                    term_0_2 = sum(log(phi.^Eta_shift_for(conn_neurons,t)));
 %                   term_0_3 = sum(log((1-phi).^(1-Eta_shift_for(conn_neurons,t))));
                                            temp=(Eta_shift_for(conn_neurons,t));
                        term_0_2 = sum(log(phi(temp)));
                        term_0_3 = sum(log((1-phi(~temp))));

                else
                    term_0_2 = 0;
                    term_0_3 = 0;
                end
                % prior probability has different expressions for t = 1 and
                % rest
                if t~=1
                    % log of prior probability for the unflipped sign
                    if Eta_t(n)==1
                        term_0_1 = -log((1+exp(-(Theta(:,n)'*X_odd(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    else
                        term_0_1 = -log((1+exp((Theta(:,n)'*X_odd(:,t)+W(n,:)*Eta_shift_back(:,t)+thres(n)))));
                    end
                    % log probability of not flipped sign (upto
                    % proportionality)
                    log_prob_0_prop = term_0_1+term_0_2+term_0_3;
                else
                    if Eta_t(n)==1
                        term_0_1 = -log((1+exp(-(Theta(:,n)'*X_odd(:,t)+thres(n)))));
                    else
                        term_0_1 = -log((1+exp((Theta(:,n)'*X_odd(:,t)+thres(n)))));
                    end
                    % log probability of not flipped sign (upto
                    % proportionality)
                    log_prob_0_prop = term_0_1+term_0_2+term_0_3;
                end
                
                
                
                % calculate acceptance probability of flipped sign
                %acc_prob = min(1,exp(log_prob_1_prop-log_prob_0_prop));
                %accept = binornd(1,acc_prob);
                if log_prob_1_prop > log_prob_0_prop
                    Eta_t(n)= (1-Eta_t(n)); %accept=1; 
                else
                    accept = (rand<exp(log_prob_1_prop-log_prob_0_prop));
                    if accept
                        Eta_t(n)= (1-Eta_t(n));
                    end 
                end;
                %Eta_t(n)= (1-Eta_t(n))*accept+Eta_t(n)*(1-accept);
            end
        end
        Eta_odd_sample(:,t) = Eta_t;
    end
    Eta(:,1:2:T) = Eta_odd_sample;
    % only save the samples after the burnin
    if samples>burnin
        Etas{samples-burnin} = Eta;
    end
    last_Eta = Eta;
end

end