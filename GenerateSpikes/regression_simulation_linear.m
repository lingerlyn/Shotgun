function s=regression_simulation_linear(A,b,T,seed)
% This function simulates a linear regression in each 2 timebins
% A - network connectivity (NxN)
% b - bias (Nx1)
% T - simulation duration (scalar)
% seed - random seed
% and outputs 
% s - network activity (NxT)

sigma_noise=1;

N=size(A,1);
s=zeros(N,T);

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

for tt=1:T
    if mod(tt,2)==1
        s(:,tt)=b+sigma_noise*randn(N,1);
    elseif mod(tt,2)==0
        s(:,tt)=A*s(:,tt-1)+b+sigma_noise*randn(N,1);
    end
end

end