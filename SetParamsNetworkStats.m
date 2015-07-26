function params=SetParamsNetworkStats()

%% Network parameters
N=1e3; %number of neurons
N_stim=0; %number of stimulation sources
N_unobs=round(0.2*N); %number of neurons completely unobserved 
N=N_unobs+N;
spar =0.2; %sparsity level - set as empty for default value in realistic conn_type; 
bias=-0.5*ones(N,1)+0.1*randn(N,1);  %bias  - if we want to specify a target rate and est the bias from that instead
target_rates=[0.05]; %set as empty if you want to add a specific bias.
seed_weights=1; % random seed
weight_scale=1;%1/sqrt(N*spar*2); % scale of weights  
conn_types={'realistic', 'realistic+1','rand','common_input', 'balanced', 'balanced2','NIPS_network'};
conn_type=conn_types{1};
inhib_frac=0.2;
weight_dist_types={ 'lognormal','uniform'};
weight_dist=weight_dist_types{1}; %
connectivity=v2struct(N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates,N_unobs);

%% Spike Generation parameters
T=1e4; %timesteps
T0=1e2; %burn-in time 
sample_ratio=1; %fraction of observed neurons per time step
neuron_type_set={'logistic','logistic_with_history','linear','linear_reg', 'sign','Poisson','LIF'};
neuron_type=neuron_type_set{2}; 
sample_type_set={'continuous','spatially_random','prob','double_continuous'};
sample_type=sample_type_set{2};
stim_type_set={'pulses','delayed_pulses','sine','Markov','white'};
stim_type=stim_type_set{5};
timescale=1; %timescale of filter in neuronal type 'logistic_with_history' - does not affect anything in other neuron models
seed_spikes=1;
seed_sample=1e6;
obs_duration=100; %duration we observe each neurons
CalciumObs=0; %use a calcium observation model

spike_gen=v2struct(T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration,CalciumObs);

%% Calcium model parameters- only relevant if  CalciumObs=1
sn=0.2;  % observation noise level
amp=1; % spike amplitude
baseline=0.3;  % baseline
g=0.97; % calcium rate (AR model poles - can be a vector
est_type_set={'FOOPSI','GreedyAccurate','GreedyFast'};
est_type=est_type_set{2};

% default values I used in paper:
% sn=0.2;  % observation noise level
% amp=1; % spike amplitude
% baseline=0.3;  % baseline (does not matter)
% g=0.97; % calcium rate (AR model poles - can be a vector)
% est_type = GreedyAccurate
calcium_model=v2struct(sn,amp,baseline,g,est_type);

%% Connectivity Estimation Flags
pen_diag=0; %penalize diagonal entries in fista
pen_dist=0; %penalize distance in fista
warm=1; %use warm starts in fista
est_type_set={'ELL','Cavity','Gibbs','FullyObservedGLM','EM'};
est_type=est_type_set{2};
est_spar=1;% estimated sparsity level. If empty, we just use the true sparsity level
conn_est_flags=v2struct(pen_diag,pen_dist,warm,est_type,est_spar);

%% Sufficeint Statistics Estimation flags - obsolete
glasso=0; %use glasso?
restricted_penalty=0; % use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)?
pos_def=0; % use positive semidefinite projection?
bin_num=1e3; %number of bins in marginal estimation of fitlered spikes (only relevant if timescale>1)
stat_flags=v2struct(glasso,pos_def,restricted_penalty,est_spar,bin_num); %add more...

%% SBM parameters - obsolete
if strcmp(conn_type,'block')
    sbm=SetSbmParams(N,weight_scale); %does not work - please correct this
else
    sbm=[];
    MeanMatrix=eye(N+N_stim);
    DistDep=0;
end

%% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,conn_est_flags,calcium_model,sbm);
% end