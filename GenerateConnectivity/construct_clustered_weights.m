function W=construct_clustered_weights(N,K,spar)
% Constructs of random weights divides into K-sized clusters, where each cluster is not connectec to the other clusters
% INPUTS:
% N - number of neurons
% K - size of blocks
% spar - sparsity level

if mod(N,K)
    error('N must be divisible into the right number of blocks')
end

NumK=N/K;

    %create weight matrix
    W=zeros(N);
    for nn=1        
            W((nn-1)*NumK+(1:K),(nn-1)*NumK+(1:K))=randn(K).*(rand(K)<spar);
    end
    
    

end
        

        
        
        
        
        
        

