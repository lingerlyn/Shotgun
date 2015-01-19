function params=SetParams()

%Network parameters
N=50; %number of neurons
N_stim=0; %number of stimulation sources
spar =0.2; %sparsity level; 
bias=-1.2*ones(N,1)+0.1*randn(N,1); %bias  - if we want to specify a target rate and est the bias from that instead
target_rates=[]; %set as empty if you want to add a specific bias.
seed_weights=1; % random seed
weight_scale=1;%1/sqrt(N*spar*2); % scale of weights  
conn_type='prob';
connectivity=v2struct(N,spar,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates);
% Spike Generation parameters
T=1e5; %timesteps
T0=1e2; %burn-in time 
sample_ratio=1; %fraction of observed neurons per time step
neuron_type_set={'logistic','logistic_with_history','linear','linear_reg', 'sign'};
neuron_type=neuron_type_set{1}; 
sample_type_set={'continuous','fixed_subset','spatially_random','prob'};
sample_type=sample_type_set{3};
stim_type_set={'pulses','delayed_pulses','sine','Markov'};
stim_type=stim_type_set{1};
timescale=1; %timescale of filter in neuronal type 'logistic_with_history' - does not affect anything in other neuron models
seed_spikes=1;
seed_sample=1;
spike_gen=v2struct(T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale);

% Sufficeint Statistics Estimation flags
glasso=0; %use glasso?
restricted_penalty=0; % use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)?
pos_def=0; % use positive semidefinite projection?
est_spar=[];% estimated sparsity level. If empty, we "cheat" - we just use the true sparsity level
stat_flags=v2struct(glasso,pos_def,restricted_penalty,est_spar); %add more...

% Connectivity Estimation Flags
pen_diag=0; %penalize diagonal entries in fista
warm=0; %use warm starts in fista
est_type_set={'ELL','Gibbs','FullyObservedGLM'};
est_type=est_type_set{1};
conn_est_flags=v2struct(pen_diag,warm,est_type);

% SBM parameters
if strcmp(conn_type,'block')
    sbm=SetSbmParams(N,weight_scale); %does not work - please correct this
else
    sbm=[];
    MeanMatrix=eye(N+N_stim);
    DistDep=0;
end

% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,conn_est_flags,sbm);
end