function s=network_simulation_logistic_with_history(A,b,T,T0,seed,timescale,s0,verbos)
% This function simulates a network with parameters
% A - network connectivity (NxN)
% b - bias (NxT)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% timescale - of neuronal integration
% s0 - (Nx1) initial condition. can be empty
% outputs 
% s - network spikes (NxT)

gamma=1/timescale; %rate

N=size(A,1);
s=zeros(N,T);
f = @(x) 1./(1+exp(-x));

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

% T0=1e2; %time to wait so network activity becomes stationary
v=0; % "voltage" - integation variable

if isempty(s0)
    s0=rand(N,1)<0.5;   

    for tt=1:T0
        v=v*(1-gamma)+gamma*f(A*s0+b(:,1));
        s0=rand(N,1)<v;
    end
end

s(:,1)=s0;
if verbos
    disp(['Generating spikes...'])
end
for tt=1:(T-1)
    v=(1-gamma)*v+gamma*f(A*s(:,tt)+b(:,tt));
    s(:,tt+1)=rand(N,1)<v;
    
    if verbos
        if ~mod(tt,floor(T/10))
            disp([num2str(100*tt/T,2) '%'])
        end
    end
end

end