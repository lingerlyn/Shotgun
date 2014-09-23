clear all
close all
clc
tic

%% Set params - later write a function for several default values
addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')

%Network parameters
N=50; %number of neurons
N_stim=0; %number of stimulation sources
spar =0.2; %sparsity level; 
bias=-0.5*ones(N,1)+0.1*randn(N,1); %bias  - if we want to specify a target rate and est the bias from that instead
target_rates=[]; %set as empty if you want to add a specific bias.
seed_weights=1; % random seed
weight_scale=1; % scale of weights  1/sqrt(N*spar*2)
conn_type='balanced';
connectivity=v2struct(N,spar,bias,seed_weights, weight_scale, conn_type);

% Spike Generation parameters
T=2e6; %timesteps
T0=1e2; %burn-in time 
sample_ratio=0.3; %fraction of observed neurons per time step
neuron_type='logistic'; %'logistic' or 'linear' or 'sign' or 'linear_reg'
sample_type_set={'continuous','fixed_subset','spatially_random'};
sample_type=sample_type_set{3};
stim_type='sine';
seed_spikes=1;
seed_sample=1;
spike_gen=v2struct(T,T0,sample_ratio,sample_type,seed_spikes,N_stim,stim_type, neuron_type);

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

% SBM parameters
if strcmp(conn_type,'block')
    Realistic=1; %Adhere to Dale's law and add a negative diagonal
    DistDep=1;
    blockFracs=[1/3;1/3;1/3];
    nblocks=length(blockFracs);
    abs_mean=10^(-0.31);
    str_var=0.1;%10^(-0.6);
    noise_var=0.1;
    pconn=spar*ones(nblocks);
    
    %Dale's law
    if Realistic
        block_means=abs_mean*ones(nblocks).*[ones(nblocks,1) -1*ones(nblocks,1) ones(nblocks,1)];
%         c=-2*abs_mean*weight_scale; %self-inhibition.
        c=-1*weight_scale;%self-inhibition.
        MeanMatrix=GetBlockMeans(N,blockFracs,block_means)*weight_scale;
        MeanMatrix(~~eye(N))=c;
        sbm=v2struct(Realistic,DistDep,blockFracs,nblocks,abs_mean,str_var,noise_var,pconn,block_means,MeanMatrix,c);
    else
        block_means=abs_mean*(ones(nblocks)-2*eye(nblocks)); %default blockmodel
        MeanMatrix=GetBlockMeans(N,blockFracs,block_means)*weight_scale;
        sbm=v2struct(Realistic,DistDep,blockFracs,nblocks,abs_mean,str_var,noise_var,pconn,block_means,MeanMatrix);
    end
    
else
    sbm=[];
end

% Combine all parameters 
params=v2struct(connectivity,spike_gen,stat_flags,conn_est_flags,sbm);
%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
W=GetWeights(N,conn_type,spar,seed_weights,weight_scale,N_stim,sbm);

if ~isempty(sbm) && sbm.Realistic %add self-inhibition diag
    temp=W(1:N,1:N);
    temp(~~eye(N))=c+randn(N,1)*sqrt(str_var)*weight_scale*.1; %10x less variance
    W(1:N,1:N)=temp;
end

if ~isempty(target_rates)
    bias=SetBiases(W,target_rates*ones(N,1));
end

%% Distance-depdendent Connectivity
sigm=@(z) 1./(1+exp(-z));
distfun=@(a,b,x) sigm(a*x+b);

if ~isempty(sbm) && sbm.DistDep
    a=-5;
    b=GetDistDepBias(a,spar);

    %Generate some uniform random distances
    D=triu(rand(N),1); D=D+D'; DD=distfun(a,b,D);
%     figure; plot(D(:),DD(:)); %look at dist func
    W(1:N,1:N)=W(1:N,1:N).*(DD>rand(N));
    
end

%% Generate Spikes
addpath('GenerateSpikes');
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type);
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;

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

%% Estimate sufficeint statistics
addpath('EstimateStatistics')
[Cxx, Cxy,EW3,rates,obs_count] = GetStat(sampled_spikes,observations,glasso,restricted_penalty,pos_def,est_spar,W);
% Ebias=GetBias( EW,Cxx,rates);
%% Estimate Connectivity
addpath('EstimateConnectivity');
% [EW2,alpha, rates_A, s_sq]=EstimateA(Cxx,Cxy,rates,obs_count,est_priors);
if strcmp(neuron_type,'logistic')
    V=diag(rates);
else
    V=eye(N);
end
   
EW3=Cxy'/(V*Cxx);
EW=EstimateA_L1(Cxx,Cxy,est_spar);
% EW=EstimateA_L1_logistic(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);
Ebias=GetBias( EW,Cxx,rates);

if strcmp(neuron_type,'logistic')
    [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
else
    amp=1;
    Ebias2=Ebias;
end

EW2=diag(amp)*EW;
% EW2=median(amp)*EW;  %somtimes this works better...
% EW=EstimateA_L1_logistic_known_b(Cxx,Cxy,bias,est_spar);

%% Remove stimulus parts
if N_stim>0
    W=W(1:N,1:N);
    EW=EW(1:N,1:N);
    EW2=EW2(1:N,1:N);
    Ebias=Ebias(1:N);
    Ebias2=Ebias2(1:N);
    spikes=spikes(1:N,:);
end
%% Save Results
t_elapsed=toc
params.t_elapsed=t_elapsed;
file_name=GetName(params);  %need to make this a meaningful name
save(file_name,'W','bias','EW','EW2','Cxx','Cxy','Ebias','Ebias2','params');

%% Plot
Plotter

