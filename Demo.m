% Run this script for a simple demo of the algorithm

clear all
close all
clc

% dimension
N=100; 
% # time instants
M=1e5; 


%% Generate network connectivity
addpath('GenerateConnectivity')
% connectivity parameters (for details see "SetParams.m")
conn_type='realistic'; spar=0.1;inhib_frac=0.2;weight_dist='lognormal';seed_weights=2;weight_scale=1;N_stim=0;
% connectivity matrix
[W,~]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,[],[]);
%bias
b=-3+0.1*randn(N,1); 


%% Generate spikes
addpath('GenerateSpikes');
% spiking parameters (for details see "SetParams.m")
T0=1e2;seed_spikes=1;neuron_type='logistic_with_history';N_stim=0; stim_type='pulses';timescale=1;s0=[],verbose=1;
% Get Spikes
S=GetSpikes(W,b,M,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbose);

%observation matrix
O=1+0*S;  % currently, assume full observations
 %sampled spikes
SS=O.*S;

%% Calcuate sufficient statistics
mY=sum(SS,2); mYn=sum(O,2);
rates=mY./(mYn+eps); %estimate the mean firing rates

XX=SS*SS'; XXn=O*O';
Cxx=XX./(XXn+eps)-rates*rates'; %estimate the covariance (not including stim terms for now)

XY=SS(:,1:(end-1))*(SS(:,2:end))'; XYn=O(:,1:(end-1))*(O(:,2:end))';
Cxy=XY./(XYn+eps)-rates*rates'; %estimate cross-covariance


%% Estimate connectivity
est_spar=spar; %sparsity target (set here to correct sparisty level -in general this should estiamted from data)
N_stim=0;  %number of external stimuli
pen_diag=0;pen_dist=0; warm=1;W_obs=[]; centers=[];  % (for details see "SetParams.m")

%Main algorithm of the dpaper
[EW,Ebias,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,pen_dist,warm,W_obs,centers);     

%% Plot
colormap('jet')
figure(1)
a=2;b=2;
subplot(a,b,1)
imagesc(W); colorbar; title('W')
subplot(a,b,2)
imagesc(S); title('S')
subplot(a,b,3)
imagesc(EW); colorbar; title('Inferred W')
subplot(a,b,4)
plot(W(:),EW(:),'.'); grid on
xlabel('W'); ylabel('Inferred W')