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

if ~isempty(target_rates)
    bias=SetBiases(W,target_rates*ones(N,1));
end

if isempty(params.stat_flags.est_spar)
    params.stat_flags.est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
end

%sorted W
% if DistDep
%     [~,idx]=sort(neuron_positions);
%     sortedW=W(idx,idx);
% else
%     sortedW=[];
% end

%% Generate Spikes
addpath('GenerateSpikes');
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale]=v2struct(params.spike_gen);

tic
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale);
RunningTime.GetSpikes=toc;
tic
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;
RunningTime.SampleSpikes=toc;

% spikes=sparse(GetSpikes(W,bias,T,T0,seed_spikes,neuron_type));
% observations=sparse(SampleSpikes(N,T,sample_ratio,sample_type,seed_sample+1));
% sampled_spikes=sparse(observations.*spikes);
%% Handle case of fixed observed subset
if strcmp(sample_type,'fixed_subset')||strcmp(sample_type,'random_fixed_subset')
    ind=any(observations,2);        
    W=W(ind,ind);    
    spikes=spikes(ind,1:end);
    sampled_spikes=sampled_spikes(ind,1:end);
    observations=observations(ind,1:end);
    ind(end+1-N_stim:end)=[];
    N=sum(ind);% note change in N!!!
    bias=bias(ind);
    params.stat_flags.est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
end

%% Estimate sufficeint statistics
addpath('EstimateStatistics')
[glasso,pos_def,restricted_penalty,est_spar]=v2struct(params.stat_flags);

tic
if timescale==1
 [Cxx, Cxy,~,rates,obs_count] = GetStat(sampled_spikes,observations,glasso,restricted_penalty,pos_def,est_spar,W);
else 
    filter_list=cell(2,N);
    gamma=1/timescale;
    for nn=1:N
        filter_list{1,nn}=[1 -(1-gamma)];
        filter_list{2,nn}=gamma;
    end
    [Cxx, Cxy,~,rates,obs_count] = GetStat2(sampled_spikes,observations,filter_list,glasso,restricted_penalty,pos_def,est_spar,W);
end
RunningTime.GetStat=toc;
% Ebias=GetBias( EW,Cxx,rates);
%% Estimate Connectivity
addpath('EstimateConnectivity');
[pen_diag,warm,est_type]=v2struct(params.conn_est_flags);

tic

EW3=Cxy'/Cxx;

if strncmpi(neuron_type,'logistic',8)
    [amp, Ebias2]=logistic_ELL(rates,EW3,Cxx,Cxy);    
else
    amp=eye(N);
end
EW3=diag(amp)*EW3;

switch est_type
    case 'Gibbs'
        p_0=est_spar*ones(N);
        if ~pen_diag
            p_0(eye(N)>0.5)=1;
        end
%         mu_0=-eye(N);
        mu_0=zeros(N);
        std_0=std(W(~~W(~eye(N))))*ones(N);
        EW = EstimateA_Gibbs( bias,spikes,observations,p_0, mu_0, std_0);
        Ebias=GetBias( EW,Cxx,rates);
        if strncmpi(neuron_type,'logistic',8)
            [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        else
            amp=1;
            Ebias2=Ebias;
        end
        EW2=diag(amp)*EW;
    case 'ELL'
        EW=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);              
        Ebias=zeros(N,1);
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        EW=diag(amp)*EW;
        if timescale==1
            is_spikes=1;
        else
            is_spikes=0;
        end
        EW2=EstimateA_L1_logistic_sampling(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm,is_spikes);       
    case 'FullyObservedGLM'
%         temp=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);      
%         [amp, Ebias2]=logistic_ELL(rates,temp,Cxx,Cxy);
%         EW=diag(amp)*temp;
        EW=EstimateA_L1_logistic_sampling(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);             
        [EW2, Ebias]=EstimateA_L1_logistic_fullyobserved(Cxx,Cxy,rates,spikes,est_spar,N_stim,pen_diag,warm);
end



RunningTime.EstimateWeights=toc;

CheckMore=0;
if CheckMore
   
    %OMP
    omp_lambda=0; %#ok
    EW_OMP=EstimateA_OMP(Cxx,Cxy,spar,omp_lambda,MeanMatrix,rates);
    if strncmpi(neuron_type,'logistic',8)
        [amp, ~]=logistic_ELL(rates,EW_OMP,Cxx,Cxy);
    else
        amp=1;
    end
    EW_OMP_ELL=diag(amp)*EW_OMP;
    if Realistic

        %Dale's Law L1
        [EW_DL1,idents]=EstimateA_L1_logistic_Accurate_Dale_Iter(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm,idenTol);
        if strncmpi(neuron_type,'logistic',8)
            [amp, Ebias2]=logistic_ELL(rates,EW_DL1,Cxx,Cxy);
        else
            amp=1;
            Ebias2=Ebias;
        end
        EW_DL1_ELL=diag(amp)*EW_DL1;

        lambda=0; %Dale's law OMP (linear)
        [EW_DOMP,identsDOMP]=EstimateA_OMP_Dale_Iter(Cxx,Cxy,rates,est_spar,lambda,MeanMatrix,idenTol);
        if strncmpi(neuron_type,'logistic',8)
            [amp, Ebias2]=logistic_ELL(rates,EW_DOMP,Cxx,Cxy);
        else
            amp=1;
            Ebias2=Ebias;
        end
        EW_DOMP_ELL=diag(amp)*EW_DOMP;

        %Dale's law OMP Exact
        lambda=0;
        [EW_DOMP_Exact,identsDOMP_Exact]=EstimateA_OMP_Exact_Dale_Iter(Cxx,Cxy,rates,est_spar,lambda,MeanMatrix,idenTol);
        if strncmpi(neuron_type,'logistic',8)
            [amp, Ebias2]=logistic_ELL(rates,EW_DOMP_Exact,Cxx,Cxy);
        else
            amp=1;
            Ebias2=Ebias;
        end
        EW_DOMP_Exact_ELL=diag(amp)*EW_DOMP_Exact;

    end

    %Exact OMP
    lambda=0;
    EW_OMP_Exact=EstimateA_OMP_Exact(Cxx,Cxy,rates,est_spar,lambda,MeanMatrix);
    if strncmpi(neuron_type,'logistic',8)
        [amp, Ebias2]=logistic_ELL(rates,EW_OMP_Exact,Cxx,Cxy);
    else
        amp=1;
        Ebias2=Ebias;
    end
    EW_OMP_Exact_ELL=diag(amp)*EW_OMP_Exact;

    disp(['L1: ' num2str(corr(EW(:),W(:)))]);
    disp(['L1+ELL: ' num2str(corr(EW2(:),W(:)))]);
    disp(['OMP: ' num2str(corr(EW_OMP(:),W(:)))]);
    disp(['OMP+ELL: ' num2str(corr(EW_OMP_ELL(:),W(:)))]);
    disp(['OMP Exact: ' num2str(corr(EW_OMP_Exact(:),W(:)))]);
    disp(['OMP Exact+ELL: ' num2str(corr(EW_OMP_Exact_ELL(:),W(:)))]);
    if Realistic
        disp(['Dales Law L1: ' num2str(corr(EW_DL1(:),W(:)))]);
        disp(['Dales Law L1+ELL: ' num2str(corr(EW_DL1_ELL(:),W(:)))]);
        disp(['Dales Law OMP: ' num2str(corr(EW_DOMP(:),W(:)))]);
        disp(['Dales Law OMP+ELL: ' num2str(corr(EW_DOMP_ELL(:),W(:)))]);
        disp(['Dales Law OMP Exact: ' num2str(corr(EW_DOMP_Exact(:),W(:)))]);
        disp(['Dales Law OMP Exact+ELL: ' num2str(corr(EW_DOMP_Exact_ELL(:),W(:)))]);
    end

end


%% Remove stimulus parts
if N_stim>0
    W=W(1:N,1:N);
    EW=EW(1:N,1:N);
    EW2=EW2(1:N,1:N);
    EW3=EW3(1:N,1:N);
    Ebias=Ebias(1:N);
    Ebias2=Ebias2(1:N);
    spikes=spikes(1:N,:);
end
%% Save Results

params.RunningTime=RunningTime;
file_name=GetName(params);  %need to make this a meaningful name
% save(file_name,'W','bias','EW','EW2','V','Cxx','Cxy','rates','Ebias','Ebias2','params');
% end
%% Plot
Plotter

% CommonInputPlotC