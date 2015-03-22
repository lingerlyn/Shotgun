function [W,centers] = GetWeights( N,network_type,spar,inhib_frac,weight_dist,seed,scale,N_stim,stim_type,sbmparams)
%GetWeights Summary of this function goes here
% INPUT
% N - number of neurons
% network_type - what kind of connectivity to use (note each connectivity has additional parameters)
% spar - sparsity level
% seed - seed used to generate the random connections
% scale - weights are multiplied by this constant
% N_stim - number of stimuli
% stim_type - stimulus type
% inhib_frac - fraction of inhibitory neurons
% weight_dist - type of weight distribution 
% OUTPUT
% W = N x N weight matrix
% centers = position of all neurons

if ~isempty(seed)
    stream = RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(stream);
end
centers=[];

switch network_type
    case 'combi'
        lr_conn = 0.2;
        A=construct_weights_combi(N, spar,lr_conn);
    case 'prob'        
        A=construct_weights_probabilistic(N, spar,inhib_frac,weight_dist);
    case 'realistic'        
        [A,centers]=construct_weights_realistic(N,inhib_frac,spar);
    case 'balanced'   
        lr_conn = 0.2;
        A = construct_bal_weights(N,spar,lr_conn);     
    case 'balanced2'   
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
    
    if strcmp(stim_type,'Markov')
        % use clusterd inputs
        G=(rand(N,N_stim)-0.5);
        for nn=1:N_stim
            kk=randi(round(1/spar));
            ind_low=(1:N)<max(1,round((kk-1)*N*spar));
            ind_high=(1:N)>min(N,round(kk*N*spar));        
            G(ind_low,nn)=0;
            G(ind_high,nn)=0;
        end
        
        % use soft clustered structure
        for kk=1:round(1/spar)
                ind_low=(1:N)<max(1,round((kk-1)*N*spar)+1);
                ind_high=(1:N)>min(N,round(kk*N*spar)-1);
                ind_range=((~ind_low).*(~ind_high))>0.5;
                A(ind_range,ind_range)=5*A(ind_range,ind_range);
        end
    else
%         G=rand(N,N_stim)<spar;
        G=ones(N,N_stim);
    end
%     G=zeros(N,N_stim);
%     G(1,1)=1;
%     G(2,2)=1;   
    W=[A, G; zeros(N_stim,N_stim+N)]; 

end

