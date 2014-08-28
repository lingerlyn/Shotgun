function W = GetWeights( N,network_type,spar, seed,scale,N_stim,params)
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
        lr_conn = 0.01;
        A=construct_weights_combi(N, spar,lr_conn);
    case 'balanced'   
        lr_conn = 0.01;
        A = construct_bal_weights(N,spar,lr_conn);
    case 'circular'
        NN_range=5;
        A = construct_weights_circ_NN(N,NN_range);
    case 'rand'
        A=randn(N).*(rand(N)<spar)/(sqrt(N*spar));
    case 'block'
        nTypes=length(params.sbm.blockFracs);
        str_mean=params.sbm.abs_mean*ones(nTypes)-2*params.sbm.abs_mean*eye(nTypes); %hard-coded pattern
        str_var=params.sbm.str_var*ones(nTypes);
        seed=randi(1000);
        A=construct_block_weights(N,params.sbm.blockFracs,str_mean,str_var,params.sbm.pconn,seed);
    otherwise
        error('unknown connectivity type!')
end

    A=A*scale;
    G=rand(N,N_stim);
    W=[A, G; zeros(N_stim,N_stim+N)]; 

end

