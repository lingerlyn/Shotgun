clear all
close all
clc

%% Set params - later write a function for several default values
addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')

params=SetParams;

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
[N,spar,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates]=v2struct(params.connectivity);

tic
W=GetWeights(N,conn_type,spar,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
RunningTime.GetWeights=toc;

%% Generate Spikes
addpath('GenerateSpikes');
addpath('CalciumMeasurements');
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale]=v2struct(params.spike_gen);

tic
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale);
RunningTime.GetSpikes=toc;

tic
noise_std=0.1;
amp=0.3;
b=0.2;
gamma=1/timescale;
Y = Spikes2Calcium(spikes,gamma,b,amp,noise_std);
RunningTime.GetSpikes=toc;

tt=1:T;
plotyy(tt,spikes,tt,Y)
