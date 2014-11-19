function W = GetWeights( N,network_type,spar, seed,scale,N_stim,sbmparams)
%GetWeights Summary of this function goes here
% N - number of neurons
% network_type - what kind of connectivity to use (note each connectivity has additional parameters)
% spar - sparsity level
% seed - seed used to generate the random connections
% scale - weights are multiplied by this constant

if ~isempty(seed)
    stream = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(stream);
end

switch network_type
    case 'combi'
        lr_conn = 0.2;
        A=construct_weights_combi(N, spar,lr_conn);
    case 'prob'
        inhib_frac = 0.5;
        A=construct_weights_probabilistic(N, spar,inhib_frac);
    case 'realistic'
        lr_conn = 0.1;
        A=construct_weights_realistic(N, spar,lr_conn);
    case 'balanced'   
        lr_conn = 0.1;
        A = construct_bal_weights(N,spar,lr_conn);      
    case 'common_input'   
        lr_conn = 0.1;        
        A=construct_weights_common_input(N,spar,lr_conn);
    case 'circular'
        NN_range=5;
        A = construct_weights_circ_NN(N,NN_range);
    case 'rand'
        A=randn(N).*(rand(N)<spar);
    case 'cluster'
        K=10;
        A=construct_clustered_weights(N,K,spar);
    case 'block'
        nTypes=length(sbmparams.blockFracs);
        str_var=sbmparams.str_var*ones(nTypes);
        seed=randi(1000);
        if sbmparams.DistDep
            A=construct_block_weights(N,sbmparams.blockFracs,sbmparams.block_means,str_var,ones(sbmparams.nblocks),seed);
            %connectivity
            A(1:N,1:N)=A(1:N,1:N).*(sbmparams.fd>rand(N) | ~~eye(N)*sbmparams.Realistic );
        else
            A=construct_block_weights(N,sbmparams.blockFracs,sbmparams.block_means,str_var,sbmparams.pconn,seed);
        end
        
        %fix the diagonal
        if sbmparams.Realistic
            A(~~eye(N))=scale*sbmparams.self_inhibition;
        end
    otherwise
        error('unknown connectivity type!')
end

    A=A*scale; 
%     G=scale*(rand(N,N_stim)<spar);
    G=ones(N,N_stim);
%     G=zeros(N,N_stim);
%     G(1,1)=1;
%     G(2,2)=1;   
    W=[A, G; zeros(N_stim,N_stim+N)]; 

end

