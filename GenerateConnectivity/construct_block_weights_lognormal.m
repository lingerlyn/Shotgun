function [W,labels]=construct_block_weights_lognormal(N,K,str_mean,str_var,pconn,seed)
% Constructs a block matrix of different spike-and-slab distributions.
% INPUTS:
% N - number of neurons
% K - block fractions
% str_mean - slab means for each block
% str_var - slab variance for each block
% pconn - inclusion probability for each block pair
% seed - seed
% OUPUTS:
% W - weight matrix
% labels - block identity indices for each neuron.

    numK=length(K);

    if mod(N,numK)
        error('N must be divisible into the right number of blocks')
    end

    stream = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(stream);

    % Construct labels of network first
    labels=zeros(N,1);
    labels( 1:K(1)*N )=1;
    for k=2:numK
        labels( sum(K(1:k-1))*N+1: sum(K(1:k-1))*N+K(k)*N) = k;
    end

    %create weight matrix
    W=zeros(N);
    for n=1:N
        for nn=1:N
            m=str_mean(labels(n),labels(nn));
            v=str_var(labels(n),labels(nn));
            
            mu=log(m^2+(m^2+v));
            si=sqrt(log(1+v/m^2));
            W(n,nn)=sign(m)*exp(mu+randn*si);
            W(n,nn)=W(n,nn)*(rand<pconn(labels(n),labels(nn)));
        end
    end
    
    
    

end
        

        
        
        
        
        
        

