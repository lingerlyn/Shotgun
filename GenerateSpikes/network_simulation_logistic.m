function s=network_simulation_logistic(A,b,T,seed)
% This function simulates a network with parameters
% A - network connectivity (NxN)
% b - bias (Nx1)
% T - simulation duration (scalar)
% seed - random seed
% and outputs 
% s - network spikes (NxT)

N=size(A,1);
s=zeros(N,T);
f = @(x) 1./(1+exp(-x));

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

T0=1e2; %time to wait so network activity becomes stationary
s0=rand(N,1)<0.5;
for tt=1:T0
    s0=rand(N,1)<f(A*s0+b);
end

s(:,1)=s0;

for tt=1:(T-1)
    s(:,tt+1)=rand(N,1)<f(A*s(:,tt)+b);
end

end