clear all
% close all
clc

%% Set params - later write a function for several default values
addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')

%Network parameters
N=100; %number of neurons
N_stim=0; %number of stimulation sources
spar =0.2; %sparsity level; 
bias=-1.2*ones(N,1)+0.1*randn(N,1); %bias  - if we want to specify a target rate and est the bias from that instead
target_rates=[]; %set as empty if you want to add a specific bias.
seed_weights=1; % random seed
weight_scale=1;%1/sqrt(N*spar*2); % scale of weights  
conn_type='block';
connectivity=v2struct(N,spar,bias,seed_weights, weight_scale, conn_type,N_stim);

% Spike Generation parameters
T=1e5; %timesteps
T0=1e2; %burn-in time 
sample_ratio=0.2; %fraction of observed neurons per time step
neuron_type='logistic'; %'logistic' or 'linear' or 'sign' or 'linear_reg'
sample_type_set={'continuous','fixed_subset','spatially_random','prob'};
sample_type=sample_type_set{3};
stim_type_set={'pulses','delayed_pulses','sine'};
stim_type=stim_type_set{1};
seed_spikes=1;
seed_sample=1;
spike_gen=v2struct(T,T0,sample_ratio,sample_type,seed_spikes,N_stim,stim_type, neuron_type);

% SBM parameters
if strcmp(conn_type,'block')
    DistDep=1; %distance dependent connectivity
    Realistic=1; %Adhere to Dale's law and add a negative diagonal
    str_var=.005; %variance of block weights
    blockFracs=[1/2;1/2];
    
    if DistDep
        sigm=@(z) 1./(1+exp(-z));
        distfun=@(a,b,x) sigm(a*x+b);
        neuron_positions=rand(N,1);
        distdep_a=-15;
        distdep_b=1.4;
        pconn=[];
        
    else
        pconn=spar*ones(length(blockFracs));
        distfun=[];
        neuron_positions=[];
        distdep_a=[];
        distdep_b=[];
    end
    
    if Realistic
        self_inhibition=-1*weight_scale; %weight of the diagonal
        idenTol=.25; %Tolerance for classifying neurons as exc or inh
    else
        self_inhibition=[];
        idenTol=[];
    end
    
end

% Sufficeint Statistics Estimation flags
glasso=0; %use glasso?
restricted_penalty=0; % use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)?
pos_def=1; % use positive semidefinite projection?
est_spar=spar;%spar; % estimated sparsity level. If empty, we "cheat". We just use the prior (not it is not accurate in some types of matrices, due to long range connections), and increase it a bit to reduce shrinkage
stat_flags=v2struct(glasso,pos_def,restricted_penalty,est_spar); %add more...

% Connectivity Estimation Flags
pen_diag=0; %penalize diagonal entries in fista
warm=0; %use warm starts in fista
conn_est_flags=v2struct(pen_diag,warm);

% Set up the SBM mean matrix
if strcmp(conn_type,'block')
    nblocks=length(blockFracs);
    block_means=ones(nblocks).*[ones(nblocks,1) -2*ones(nblocks,1)];

    %increase the rank of MeanMatrix by scaling the values a little
    block_means(1,1)=block_means(1,1)*1.1;
    block_means(1,2)=block_means(1,2)*.9;

    MeanMatrix=GetBlockMeans(N,blockFracs,block_means)*weight_scale;
    if Realistic
        MeanMatrix(~~eye(N))=self_inhibition; 
    end
    
    if DistDep
        neuron_distances=abs(repmat(neuron_positions,[1,N])-repmat(neuron_positions',[N,1]));
        fd=distfun(distdep_a,distdep_b,neuron_distances); %distance-dependent connection probs
    else
        neuron_distances=[];
        fd=[];
    end
    
    sbm=v2struct(Realistic,DistDep,blockFracs,nblocks,str_var,block_means,MeanMatrix,self_inhibition,neuron_positions,...
        neuron_distances,fd,distfun,distdep_a,distdep_b,pconn,idenTol);

else
    sbm=[];
    MeanMatrix=eye(N+N_stim);
end

% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,conn_est_flags,sbm);

% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
tic
W=GetWeights(N,conn_type,spar,seed_weights,weight_scale,N_stim,sbm);
RunningTime.GetWeights=toc;

if ~isempty(target_rates)
    bias=SetBiases(W,target_rates*ones(N,1));
end

%%
%Generate Spikes
addpath('GenerateSpikes');
tic
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type);
RunningTime.GetSpikes=toc;
tic
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;
RunningTime.SampleSpikes=toc;

% spikes=sparse(GetSpikes(W,bias,T,T0,seed_spikes,neuron_type));
% observations=sparse(SampleSpikes(N,T,sample_ratio,sample_type,seed_sample+1));
% sampled_spikes=sparse(observations.*spikes);
%% Handle case of fix obsereved subset
if strcmp(sample_type,'fixed_subset')||strcmp(sample_type,'random_fixed_subset')
    ind=any(observations,2);        
    W=W(ind,ind);    
    spikes=spikes(ind,1:end);
    sampled_spikes=sampled_spikes(ind,1:end);
    observations=observations(ind,1:end);
    ind(end+1-N_stim:end)=[];
    N=sum(ind);% note change in N!!!
    bias=bias(ind);
    est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
end
est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
%% Estimate sufficeint statistics
addpath('EstimateStatistics');
