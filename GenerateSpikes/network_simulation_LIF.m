function s=network_simulation_LIF(A,b,T,T0,seed,timescale,s0,verbos)
% This function simulates a sign thresholded network with parameters
% A - network connectivity (NxN)
% b - bias (NxT)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% timescale - of neuronal integration
% s0 - (Nx1) initial condition. can be empty
% and outputs 
% s - network activity (NxT)

sigma_noise=1;
gamma=1/timescale; %rate;

N=size(A,1);
s=zeros(N,T);

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

% T0=1e2; %time to wait so network activity becomes stationary
v=0;
if isempty(s0)
    s0=rand(N,1)<0.5;   
    
for tt=1:T0
    Input=A*s0+sigma_noise*randn(N,1);
    v=v*(1-gamma)+gamma*Input;
    s0=(sign(v+b(:,1)))/2;
    v(s0>0.5)=0; %reset
end

end

s(:,1)=s0;

for tt=1:(T-1)
        Input=A*s(:,tt)+sigma_noise*randn(N,1);
        v=v*(1-gamma)+gamma*Input;
        s(:,tt+1)=(sign(v+b(:,tt))+1)/2;
        v(s(:,tt+1)>0.5)=0;     %reset    
        
        if verbos
            if ~mod(tt,floor(T/10))
                disp([num2str(100*tt/T,2) '%'])
            end
        end
end

end