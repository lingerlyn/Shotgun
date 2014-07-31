clear all
close all
clc

% T_array=1e4*(1:100);

%% Set params - later write a function for several default values
addpath('Misc')

%Network parameters
N=100; %number of neurons
spar =0.2; %sparsity level; 
bias=0; %bias 
seed_weights=1; % random seed
weight_scale=0.1; % scale of weights
conn_type='balanced';
connectivity=v2struct(N,spar,bias,seed_weights);

% Spike Generation parameters
T=1e3; %timesteps
T0=0; %burn-in time 
sample_ratio=1; %fraction of observed neurons per time step
neuron_type='linear'  ; %'logistic'or 'linear' or 'sign' or 'linear_reg'
sample_type='fully_random';
seed_spikes=1;
seed_sample=1;
spike_gen=v2struct(T,T0,sample_ratio,sample_type,seed_spikes);

% Sufficeint Statistics Estimation flags
glasso=1; %use glasso?
restricted_penalty=0; % use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)?
pos_def=0; % use positive semidefinite projection?
est_spar=spar; % estimated sparsity level. If empty, we "cheat". We just use the prior (not it is not accurate in some types of matrices, due to long range connections), and increase it a bit to reduce shrinkage
stat_flags=v2struct(glasso,pos_def,restricted_penalty,est_spar); %add more...

% Connectivity Estimation prior parameters
est_priors.eta=zeros(N); %means
est_priors.ss2=.05*ones(N); %slab variance
est_priors.noise_var=.09*ones(N,1); %noise variance
est_priors.a=spar*ones(N);%.2*ones(N); %sparsity

% sbm=1; %use sbm?
% est=v2struct(sbm); %add more...
% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,est_priors);

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
W=GetWeights(N,conn_type,spar, seed_weights,weight_scale );
bias=bias*ones(N,1); %set bias

%% Generate Spikes
addpath('GenerateSpikes');
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type);
sampled_spikes=SampleSpikes(spikes,sample_ratio,sample_type,seed_sample);

%% Estimate sufficeint statistics
addpath('EstimateStatistics')
% temp=W;
% temp(eye(N)>0.5)=0;
% est_spar=mean(~~temp(:));
tic
[Cxx, Cxy,EW,eye_mat,rates,spike_count] = GetStat(sampled_spikes,glasso,restricted_penalty,pos_def,est_spar,W);
% [Cxx, Cxy,EW,eye_mat,rates,spike_count] = GetStat_reg(sampled_spikes,glasso,restricted_penalty,pos_def,est_spar,W);
% EW=eye_mat\EW*(mean(eye_mat(eye(N)>0.5)));
t_elapsed=toc
params.t_elapsed=t_elapsed;
Ebias=0;

%% Estimate Connectivity
addpath('EstimateConnectivity');
% [EW2,alpha, rates_A, s_sq]=EstimateA(Cxx,Cxy,rates,spike_count,est_priors);
% EW2=EstimateA_L1(Cxx,Cxy,lambda)
% EW2=Cxy'/Cxx;
EW2=EstimateA_L1(Cxx,Cxy,est_spar);
%% Save Results
file_name=GetName(params);  %need to make this a meaningful name
save(file_name,'W','bias','EW','EW2','Cxx','Cxy','Ebias','params');

%% Plot
Plotter

