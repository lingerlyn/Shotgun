function [ W ] = construct_weights_combi(N, spar,lr_conn)
%CONSTRUCT_WEIGHTS returna the weight matrix such that
% first ceil(N/3) neurons are not connected directly but are cross
% connected, and every other neuron is connected directly to ceil(0.15*N) neurons
% out of which half of the connections are with nearby neurons

%INPUT
% N = total number of neurons
% spar = sparsity level

% OUTPUT
% N x N weight matrix

%rnd(seed)
total_conn = ceil(spar*N);
%W = diag(unifrnd(-2,0,N,1)); % take diagonal wieghts to be negative
W = -2*eye(N);
not_obs_subset = (ceil(N/3)+1):N;
nearby_conn = ceil(total_conn*(1-lr_conn));
for i = not_obs_subset
    for j = 1:nearby_conn
        if i+j<=N
            %            W(i,i+j)=unifrnd(-2,2,1,1);
            W(i,i+j)=2*sign(unifrnd(-1,1));
        end
    end
    left_conn = total_conn - nearby_conn;
    other_conn = randsample(find(~W(i,:)),left_conn);
    %    W(i,other_conn)=unifrnd(-2,2,length(other_conn),1);
    W(i,other_conn)=2*ones(length(other_conn),1);
end

for i = 1:(ceil(N/3))
    conn_unobs_neu = randsample((ceil(N/3)+1):N,2);
    sgn = rand<0.5;
    if sgn ==0
        sgn = -1;
    end
    W(i,conn_unobs_neu(1)) = sgn*2;
    W(conn_unobs_neu(1),conn_unobs_neu(2)) = 2;
    W(i+1,conn_unobs_neu(2)) = sgn*2;
end
