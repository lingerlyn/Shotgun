function s=network_simulation_linear(A,b,T,T0,seed)
% This function simulates a network with parameters
% A - network connectivity (NxN)
% b - bias (NxT)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% and outputs 
% s - network activity (NxT)

sigma_noise=1;

N=size(A,1);
s=zeros(N,T);

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

% T0=1e2; %time to wait so network activity becomes stationary
s0=rand(N,1)<0.5;
for tt=1:T0
    s0=A*s0+b(:,1)+sigma_noise*randn(N,1);
end

s(:,1)=s0;

for tt=1:(T-1)
    s(:,tt+1)=A*s(:,tt)+b(:,tt)+sigma_noise*randn(N,1);
end

end