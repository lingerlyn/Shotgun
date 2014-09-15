function [ W ] = construct_weights_realistic(N, spar,lr_conn)
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

IE_ratio=5; % ratio of excitatory conenctions over inhbitory connections
if mod(N/IE_ratio,1)>0
    error('N/IE_ratio must be integer')
end
N_excite=(1-1/IE_ratio)*N;
N_inhib=(1/IE_ratio)*N;
total_conn = ceil(spar*N)-1;
%W = diag(unifrnd(-2,0,N,1)); % take diagonal weights to be negative
W = -1*eye(N);

nearby_conn = ceil(total_conn*(1-lr_conn))-1;
excite_subset = 1:N_excite;
inhib_subset = (N_excite+1):N;

%% Interpyramid spike transmission stabilizes the sparseness of recurrent network activity,  by Ikegaya
m=10^(-0.31);
v=10^(-0.3);
mu=log(m^2+(m^2+v));
si=sqrt(log(1+v/m^2));          
            
for jj = excite_subset    
    for ii = 1:nearby_conn   
            %            W(i,i+j)=unifrnd(-2,2,1,1);
            if ~mod(ii,IE_ratio)                
                W(N_excite+mod(jj+ii/IE_ratio-1,N_inhib)+1,jj)=unifrnd(0,1);%exp(mu+randn*si);
            else
                index_e=mod(ii+1+nearby_conn/2,nearby_conn+1)-nearby_conn/2-1;
                W(mod(jj+index_e-1,N_excite)+1,jj)=unifrnd(0,1)/sqrt(abs(index_e)); %exp(mu+randn*si)/sqrt(abs(index_e));
            end
    end
    left_conn = total_conn - nearby_conn;
    other_conn = randsample(find(~W(:,jj)),left_conn);
    W(other_conn,jj)=exp(mu+randn(size(other_conn))*si);
end

for jj = inhib_subset 
    for ii = 1:total_conn                
            if ~mod(ii,IE_ratio)   
                temp=inhib_subset;
                temp(temp==jj)=[];
                index_i=randsample(inhib_subset ,1);
                W(index_i,jj)=-unifrnd(0,1); %exp(mu+randn*si);
            else
                index_i=randsample(excite_subset,1);
                W(index_i,jj)=-IE_ratio*unifrnd(0,1); %exp(mu+randn*si);
            end
    end
end

for jj = 1:N
    W(jj,:) = W(jj,:)/std(W(jj,~~(W(jj,:))));
end


