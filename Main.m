clear all
close all
clc

%% Set params - later write a function for several default values
addpath('Misc')

%Network parameters
N=100; %number of neurons
spar = 0.1; %sparsity level; 
bias=0; %bias 
seed_weights=1; % random seed
conn_type='balanced';
connectivity=v2struct(N,spar,bias,seed_weights);

% Spike Generation parameters
T=1e4; %timesteps
sample_ratio=0.8; %fraction of observed neurons per time step
sample_type='fully_random';
seed_spikes=1;
seed_sample=1;
spike_gen=v2struct(T,sample_ratio,sample_type,seed_spikes);

% Sufficeint Statistics Estimation parameters
glasso=1; %use glasso?
stat=v2struct(glasso); %add more...

% Connectivity Estimation parameters
sbm=1; %use sbm?
est=v2struct(sbm); %add more...

% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat,est);

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
A=GetWeights(N,conn_type,spar, seed_weights );
bias=bias*ones(N,1); %set bias

%% Generate Spikes
addpath('GenerateSpikes');
spikes=network_simulation_logistic(A,bias,T,seed_spikes);
sampled_spikes=SampleSpikes(spikes,sample_ratio,sample_type,seed_sample);

%% Estimate sufficeint statistics
addpath('EstimateStatistics')
[Cxx, Cxy] = GetStat(sampled_spikes,stat);

%% Estimate Connectivity
addpath('EstimateConnectivity');
[EA, Eb]=EstimateA(Cxx, Cxy,est); %estimate A and b

%% Save Results
name=GetName(params);  %need to make this a meaningful name
Save(['Results\' name],'A','b','EA','Eb','params');


