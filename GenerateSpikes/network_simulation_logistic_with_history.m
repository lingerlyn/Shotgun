function s=network_simulation_logistic_with_history(A,b,T,T0,seed,timescale)
% This function simulates a network with parameters
% A - network connectivity (NxN)
% b - bias (NxT)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% timescale - of neuronal integration
% and outputs 
% s - network spikes (NxT)

gamma=1/timescale; %rate

N=size(A,1);
s=zeros(N,T);
f = @(x) 1./(1+exp(-x));

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

% T0=1e2; %time to wait so network activity becomes stationary
s0=rand(N,1)<0.5;
v=0; % "voltage" - integation variable

for tt=1:T0
    v=v*(1-gamma)+gamma*f(A*s0+b(:,1));
    s0=rand(N,1)<v;
end

s(:,1)=s0;

for tt=1:(T-1)
    v=(1-gamma)*v+gamma*f(A*s(:,tt)+b(:,tt));
    s(:,tt+1)=rand(N,1)<v;
end

end