function [ W ] = construct_weights_probabilistic(N, spar,inhib_frac)
%CONSTRUCT_WEIGHTS returns the weight matrix such that
% connectivity is randomly drawn from some diftribution

%INPUT
% N = total number of neurons
% spar = sparsity level
% inhib_frac = fraction of  inhbitory neurons in the population
% OUTPUT
% N x N weight matrix

% addpath('../Misc') % so we can use the GetProb function


W = -1*eye(N); % diagonal negative

sgn_array=ones(N,1);
sgn_array(1:(1/inhib_frac):end)=-1;

for ii=1:N
%     sgn=sign(rand-inhib_frac);
    sgn=sgn_array(ii);
    amp=1;
%     if sgn==1
%         amp=inhib_frac/(1-inhib_frac);
%     else
%         amp=1;
%     end
    
     p=GetProb(N,spar,ii);
    conn=rand(N,1)<p;
%     conn=0.5<p;
%     conn(mod((ii-1):(ii+1)-1,N)+1)=~1; % since we already have the diagonal
    conn(ii)=~1; % since we already have the diagonal
    W(conn,ii)=amp*sgn*unifrnd(0,1,sum(conn),1);
end

end



