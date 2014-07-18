function A = GetWeights( N,network_type,spar, seed,scale )
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
end

A=A*scale;

end

