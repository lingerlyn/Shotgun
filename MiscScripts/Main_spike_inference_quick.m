clear all
close all
clc

addpath(genpath('CalciumMeasurements')); %adds folder with subfolders
data_types={'Tolias','Tim','Manolis','Sim'};
data=data_types(2);
order_arpfit=1;
order_sysid=1;
plot_stuff=1;

if strcmp(data,'Sim')

%% Set params - later write a function for several default values
addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')

params=SetParams;

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
[N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates]=v2struct(params.connectivity);

tic
W=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
RunningTime.GetWeights=toc;

%% Generate Spikes
addpath('GenerateSpikes');
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration]=v2struct(params.spike_gen);

tic
s0=[]; verbose=1; iter=0;
 true_spikes=GetSpikes(W,bias,T,T0,seed_spikes+iter,neuron_type,N_stim,stim_type,timescale,s0,verbose); 
true_spikes=true_spikes(1:N,:);
RunningTime.GetSpikes=toc;

tic
sn=0.2;
amp=1;
b=3; 
g=0.95;
noise_mode=0;
Y = Spikes2Calcium(true_spikes,g,b,amp,sn,noise_mode);
tt=1:length(Y);
g=[g ; zeros(order_sysid-1,1)];
dt=1;
end
%% Use real data
data_sets_num=1; %for Tim's data
K=6; %=size(spike_est,3)
CR_cell=cell(data_sets_num,K); 
KL_cell=CR_cell;

for pp=data_sets_num

if strcmp(data,'Tolias')    
    target_folder='C:\Users\Daniel\Copy\Columbia\Research\Data\ToliasWithEphys';
    load(fullfile(target_folder,'data_before_filtering.mat'),'data');
    % tt=1:length(data{ind}.calcium);
    tt=(1:2e3)+2e4;        
    Y=double(data{pp}.calcium(tt)); %ind in 1 to 3
    true_spikes=double(data{pp}.spikes(tt));
    dt=1/data{pp}.fps;
elseif strcmp(data,'Tim')   
     [ Y,true_spikes,dt ] = GetTimData(pp); %ind in 1 to..?
     tt=1:size(Y,2);
 elseif strcmp(data,'Manolis')   
     [ Y,~,dt ] = GetManolisData('New'); %ind in 1 to..?
     tt=1:size(Y,2);
end

if ~strcmp(data,'Sim')   
    N=size(Y,1);
    g=zeros(max(order_arpfit,order_sysid),1);
    sn=0;
    b=0;
    T=length(tt);
    % Pre-process data
    for nn=1:N
        C=max(Y(nn,:));
        Y(nn,:)=2*Y(nn,:)/C; %re-normalize Y
    end
end

%% Get Spikes
P= GetParams(Y,order_arpfit,'psd','arpfit');
[spikes,relative_std_cell] = Calcium2Spikes_Greedy(Y,P);
% [spikes,b] = Calcium2Spikes(Y,P);
end
%%
CalciumInferenceQuality
