clear all
% close all
clc

addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')
addpath('GenerateSpikes');
addpath('CalciumMeasurements');

sample_ratio_array=[1];
for kk=1:length(sample_ratio_array)
%% Set params - later write a function for several default values
params=SetParams;
params.spike_gen.sample_ratio=sample_ratio_array(kk);

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
[N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates,N_unobs]=v2struct(params.connectivity);

tic
[W,centers]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
centers=centers(N_unobs+1:end,:); % remove unobserved part
RunningTime.GetWeights=toc;

    if ~isempty(target_rates)
        bias=SetBiases(W,target_rates,params.spike_gen);
    end


if isempty(params.stat_flags.est_spar)
    params.stat_flags.est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
    params.conn_est_flags.est_spar=nnz(W(1:N,1:N))/N^2; %correct sparsity estimation. Cheating????
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
addpath('CalciumMeasurements');
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration,CalciumObs]=v2struct(params.spike_gen);

memory_threshold=1e10;  %maximial size of array to allow
splits=ceil((T*N^2)/memory_threshold)
mY=0; mYn=0;
XX=0; XXn=0;
XY=0; XYn=0;
t_start=0;

for iter=1:splits
    if iter<splits
        T_split=floor(T/splits);
    else
        T_split=T-(splits-1)*floor(T/splits);
    end
    
    
    verbose=1;
    
    if iter==1
        s0=[];
    else
        s0=full(spikes(:,end));
    end
    tic
    spikes=GetSpikes(W,bias,T_split,T0,seed_spikes+iter,neuron_type,N_stim,stim_type,timescale,s0,verbose);
    RunningTime.GetSpikes=toc;
    tic
    
    if CalciumObs==1
         true_spikes=spikes;
        [ Y,spikes] = Spikes2Calcium2Spikes( true_spikes );
    else
        true_spikes=0*spikes;
        Y=spikes*0;
    end
    
    observations=SampleSpikes(N,T_split,sample_ratio,sample_type,obs_duration,N_stim,seed_sample+iter,t_start);
    t_start=t_start+T_split;
    sampled_spikes=sparse(observations.*spikes);
    RunningTime.SampleSpikes=toc;

    %% Estimate sufficeint statistics
%     addpath('EstimateStatistics')
%     [glasso,pos_def,restricted_penalty,est_spar,bin_num]=v2struct(params.stat_flags);
% 
%     tic
%     if timescale==1
%      [Cxx_split, Cxy_split,~,rates_split,obs_count_split] = GetStat(sampled_spikes,observations,glasso,restricted_penalty,pos_def,est_spar,W);
%          filter_list=cell(2,N);
%         for nn=1:N
%             filter_list{1,nn}=1;
%             filter_list{2,nn}=1;
%         end
%         p_x_split=0;
%     else 
%         filter_list=cell(2,N);
%         gamma=1/timescale;
%         for nn=1:N
%             filter_list{1,nn}=[1 -(1-gamma)];
%             filter_list{2,nn}=gamma;
%         end
%         [Cxx_split, Cxy_split,~,rates_split,obs_count_split,p_x_split] = GetStat2(sampled_spikes,observations,filter_list,glasso,restricted_penalty,pos_def,est_spar,W,bin_num);
%     end
%     RunningTime.GetStat=toc;
%     % Ebias=GetBias( EW,Cxx,rates);
        
        sampled_spikes_obs=sampled_spikes(N_unobs+1:end,:);
        observations_obs=observations(N_unobs+1:end,:);
        mY=mY+full(sum(sampled_spikes_obs,2));
        mYn=mYn+full(sum(observations_obs,2));
        XX=XX+sampled_spikes_obs*sampled_spikes_obs';
        XXn=XXn+observations_obs*observations_obs';
        XY=XY+sampled_spikes_obs(:,1:(end-1))*(sampled_spikes_obs(:,2:end))';
        XYn=XYn+observations_obs(:,1:(end-1))*(observations_obs(:,2:end))';
        
        imagesc(XYn);
        pause(1e-6);
    
end

rates=mY./(mYn+eps); %estimate the mean firing rates
Cxx=full(XX./(XXn+eps))-rates*rates'; %estimate the covariance (not including stim terms for now)
Cxy=full(XY./(XYn+eps))-rates*rates'; %estimate cross-covariance

addpath('EstimateStatistics')
 [Cxx,Cxy ] = PosProj( Cxx,Cxy );

%% Estimate Connectivity
addpath('EstimateConnectivity');
[pen_diag,pen_dist,warm,est_type,est_spar]=v2struct(params.conn_est_flags);

tic

EW3=Cxy'/Cxx;

if strncmpi(neuron_type,'logistic',8)
    [amp, Ebias2]=logistic_ELL(rates,EW3,Cxx,Cxy);    
else
    amp=ones(size(EW3,1),1);
end
EW3=diag(amp)*EW3;
EW2=EW3;
EW=EW3;

N_obs=N-N_unobs;
switch est_type
    case 'Gibbs'
        p_0=est_spar*ones(N_obs);
        if ~pen_diag
            p_0(eye(N_obs)>0.5)=1;
        end
%         mu_0=-eye(N);
        mu_0=zeros(N_obs);
        std_0=std(W(~~W(~eye(N_obs))))*ones(N_obs);
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
        Ebias=zeros(N_obs,1);
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        EW2=diag(amp)*EW;
    case 'Cavity'
%         if timescale==1
            is_spikes=1;
%         else
%             is_spikes=0;
%         end

        Cxx=full(Cxx);
        Cxy=full(Cxy);
        rates=full(rates);
        
        W_now=W(N_unobs+1:end,N_unobs+1:end);
        [EW,Ebias2,quality]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,pen_dist,warm,W_now,centers);     
         if strncmpi(neuron_type,'LIF',8)            
            EW=bsxfun(@times,EW,(std(W(~~W(:))/std(EW(~~EW(:))))));            
         end
%         EW2=EstimateA_MLE_cavity(Cxx,Cxy,rates);        %MLE
%         mask=~~EW;
%         [EW2,Ebias2,RMSE]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,1,N_stim,pen_diag,warm,is_spikes,W_now);                       
%         EW2=EW2.*mask;
%         [amp, Ebias2]=logistic_ELL(rates,EW2,Cxx,Cxy);
%         EW2=diag(amp)*EW2;
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


%% Save Results
% Remove stimulus parts
if N_stim>0
    W=W(1:N,1:N);
    EW=EW(1:N,1:N);
    EW2=EW2(1:N,1:N);
    EW3=EW3(1:N,1:N);
    Ebias=Ebias(1:N);
    Ebias2=Ebias2(1:N);
    spikes=spikes(1:N,:);
    true_spikes=true_spikes(1:N,:);
end
% Remove unobvserved parts
W_full=W;
bias_full=bias;
W=W(N_unobs+1:end,N_unobs+1:end);
bias=bias(N_unobs+1:end);
spikes=spikes(N_unobs+1:end,:);
true_spikes=true_spikes(N_unobs+1:end,:);
Y=Y(N_unobs+1:end,:);

params.connectivity.N=N-N_unobs;
params.RunningTime=RunningTime;

file_name=GetName(params);  %need to make this a meaningful name
save(file_name,'W_full','bias_full','W','bias','centers','EW','EW2','quality','Cxx','Cxy','rates','params'); %,'Ebias','Ebias2'?


end
%% Plot
Plotter

% CommonInputPlotC