clear all
close all
clc

% T_array=1e4*(1:100);

%% Set params - later write a function for several default values
addpath('Misc')

%Network parameters
N=99; %number of neurons
spar =0.2; %sparsity level; 
bias=0; %bias 
seed_weights=1; % random seed
weight_scale=0.05; % scale of weights
conn_type='block';
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


% SBM parameters
if strcmp(conn_type,'block')
    blockFracs=[1/3;1/3;1/3];
    abs_mean=1;
    str_var=.5;
    noise_var=1;
    pconn=spar*ones(length(blockFracs));
    sbm=v2struct(blockFracs,abs_mean,str_var,noise_var,pconn);
else
    sbm=[];
end

if ~isempty(sbm)

% Connectivity Estimation prior parameters
    naive=0; %use correct mean prior or zero mean prior
    if naive
        est_priors.eta=zeros(N); 
    else
        str_mean=(abs_mean*ones(length(blockFracs))-2*abs_mean*eye(length(blockFracs)))*weight_scale; %this structure is hard-coded into the sbm for now
        est_priors.eta=GetBlockMeans(N,blockFracs,str_mean); 
    end
    est_priors.ss2=sbm.str_var*ones(N)*weight_scale^2;
    est_priors.noise_var=sbm.noise_var*ones(N,1);
    est_priors.a=spar*ones(N);
end

% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,est_priors,sbm);

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
W=GetWeights(N,conn_type,spar, seed_weights,weight_scale,params);
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

