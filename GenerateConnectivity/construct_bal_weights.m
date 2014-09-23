function [ W ] = construct_bal_weights(N, spar,lr_conn)
%CONSTRUCT_WEIGHTS returna the weight matrix such that
% first N_unobs neurons are not connected directly but are cross
% connected, and every other neuron is connected directly to N_unobs/2 neurons
% out of which half of the connections are with nearby neurons

%INPUT
% N = total number of neurons
% spar = sparsity level
% lr_conn = fraction of long range connections

% OUTPUT
% N x N weight matrix

alpha=0.5; % control decay of NN weights

total_conn = ceil(spar*N);
W = -1*eye(N);

nearby_conn = ceil(total_conn*(1-lr_conn))-1;
for ii = 1:N
    for jj = 1:floor(nearby_conn/2)            
            W(ii,mod(ii+jj-1,N)+1)=unifrnd(-1,1,1,1)/jj^alpha;
    end
    for jj = 1:ceil(nearby_conn/2)
            W(ii,mod(ii-jj-1,N)+1)=unifrnd(-1,1,1,1)/jj^alpha;
    end
    left_conn = total_conn - nearby_conn-1;
    other_conn = randsample(find(~W(ii,:)),left_conn);
    %    W(i,other_conn)=unifrnd(-2,2,length(other_conn),1);
    W(ii,other_conn)=unifrnd(-1,1,1,1);
end

%% make network balanced and with zero mean
for ii = 1:N
    ind=~~W(ii,:);
    mu=sum(W(ii,ind))/(total_conn-1);
    ind(ii)=~1;
    W(ii,ind)=W(ii,ind)-mu;
    W(ii,:) = W(ii,:)/std(W(ii,~~(W(ii,:))));
end

%% Add "common hidden connections"
N_subset=ceil(N/3);
shift=floor(N_subset/2);
subset=1:N_subset;
amp=2;
num=2;
for ii = subset 
    temp=(N_subset+1):N;    
    sgn = amp*(2*(rand<0.5)-1);
    for jj=1:num
        conn_unobs_neu = randsample(temp,2);        
        W(ii,conn_unobs_neu(1)) = normrnd(sgn,0.1);
        W(conn_unobs_neu(2),conn_unobs_neu(1)) = normrnd(sgn,0.1);
        W(mod(ii+shift,N_subset)+1,conn_unobs_neu(2)) = normrnd(sgn,0.1);
    end
end

%% Add "common input"
N_subset=0;
shift=N_subset/2;
subset=1:N_subset;
amp=2;
num=2;
for ii = subset 
    temp=(N_subset+1):N;    
    sgn = amp;%*(2*(rand<0.5)-1);
    for jj=1:num
        conn_unobs_neu = randsample(temp,1);        
        W(ii,conn_unobs_neu) = normrnd(sgn,0.1);
        W(mod(ii+shift,N_subset)+1,conn_unobs_neu) = normrnd(sgn,0.1);
    end
end

