function [ W ] = construct_weights_common_input(N, W_spar,lr_conn)
%CONSTRUCT_WEIGHTS returna the weight matrix such that
% first ceil(N/3) neurons are not connected directly but are cross
% connected, and every other neuron is connected directly to ceil(0.15*N) neurons
% out of which half of the connections are with nearby neurons

%INPUT
% N = total number of neurons


% OUTPUT
% N x N weight matrix

%rnd(seed)
total_conn = ceil(W_spar*N);
%W = diag(unifrnd(-2,0,N,1)); % take diagonal wieghts to be negative
% W = -1*eye(N);
W=diag(unifrnd(-1,1,N,1));
% W=diag(randn(N,1));
obs_subset = (ceil(N/3)+1):N;
unobs_subset = 1:ceil(N/3);
nearby_conn = ceil(total_conn*(1-lr_conn));
for ii = obs_subset
    for j = 1:nearby_conn
        if ii+j<=N
            %            W(i,i+j)=unifrnd(-2,2,1,1);
            W(ii,ii+j)=normrnd(-1/total_conn,1);
        end
    end
    left_conn = total_conn - nearby_conn;
    other_conn = randsample(find(~W(ii,:)),left_conn);
    %    W(i,other_conn)=unifrnd(-2,2,length(other_conn),1);
    W(ii,other_conn)=2*ones(length(other_conn),1);
end

for ii = unobs_subset
    conn_unobs_neu = randsample(obs_subset,2);
    W(ii,conn_unobs_neu(1)) = normrnd(-1,0.1);
    W(conn_unobs_neu(1),conn_unobs_neu(2)) = normrnd(-1,0.1);
    W(ii+1,conn_unobs_neu(2)) = normrnd(-1,0.1);
end

for ii = 1:N
    if ii>unobs_subset(end)
        other_conn = randsample(find(~W(ii,:)),1);
    else
        other_conn = randsample(find(~W(ii,obs_subset))+unobs_subset(end),1);
    end
    %    W(i,other_conn)=unifrnd(-2,2,length(other_conn),1);
    W(ii,other_conn)=-sum(W(ii,:));
    W(ii,:) = W(ii,:)/std(W(ii,find(W(ii,:))));
  
    
end

