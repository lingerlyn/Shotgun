function [ W ] = construct_weights_probabilistic(N, spar,inhib_frac)
%CONSTRUCT_WEIGHTS returns the weight matrix such that
% connectivity is randomly drawn from some diftribution

%INPUT
% N = total number of neurons
% spar = sparsity level
% inhib_frac = fraction of  inhbitory neurons in the population
% OUTPUT
% N x N weight matrix

W = -1*eye(N); % diagonal negative
D=1; %number of dimentions
% centers=rand(N,D); %locations of each neuron in a 3D box
centers=linspace(0,1,N)';
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
    
    V=(pi^(D/2))/gamma(1+D/2);
    dist=sqrt(sum((mod(bsxfun(@plus,centers(ii,:),-centers)+0.5,1)-0.5).^2,2))/(spar*V)^(1/D); %distance metric - assumes all neurons are on a 3D box with cyclic boundary conditions
    p=GetProb(dist);
    conn=rand(N,1)<p;
%     conn=0.5<p;
%     conn(mod((ii-1):(ii+1)-1,N)+1)=~1; % since we already have the diagonal
    conn(ii)=~1; % since we already have the diagonal
    W(conn,ii)=amp*sgn*unifrnd(0,1,sum(conn),1);
end

end

function p=GetProb(dist)
    p=exp(-3*dist);
end

