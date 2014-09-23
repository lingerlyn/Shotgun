function s=network_simulation_logistic_with_delays(A,b,T,T0,seed)
% This function simulates a network with parameters
% A - network connectivity (NxN)
% b - bias (NxT)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% and outputs 
% s - network spikes (NxT)

N=size(A,1);
s=zeros(N,T);
f = @(x) 1./(1+exp(-x));

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

% T0=1e2; %time to wait so network activity becomes stationary
s0=rand(N,1)<0.5;

for tt=1:T0    
    s0=rand(N,1)<f(A*s0+b(:,1));
end

s(:,1)=s0;
s0=rand(N,1)<f(A*s0+b(:,1));
s(:,2)=s0;

delay=rand(1,N)<0.5;
for tt=2:(T-1)
%     delay=rand(1,N)<0.5;
    temp=s(:,tt-1:tt)';
    ind=(1:2:(2*N))+delay;
    temp=temp(ind)';
    s(:,tt+1)=rand(N,1)<f(A*temp+b(:,tt));
end

end