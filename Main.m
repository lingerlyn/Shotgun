% Run this script to explore more options

clear all
% close all
clc

addpath('Misc')
addpath('EstimateConnectivity')
addpath('GenerateConnectivity')
addpath('GenerateSpikes');
addpath(genpath('CalciumMeasurements'));

parameter_scan_array=1;

for kk=1:length( parameter_scan_array)
%% Change params inside SetParams function
params=SetParams;

%% Uncomment for parameter scans 
% params.spike_gen.sample_ratio= parameter_scan_array(kk);
% params.connectivity.N_stim=parameter_scan_array(kk);
% params.spike_gen.N_stim=parameter_scan_array(kk);

%% Generate Connectivity - a ground truth N x N glm connectivity matrix, and bias
addpath('GenerateConnectivity')
[N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates,N_unobs]=v2struct(params.connectivity);
N_obs=N-N_unobs;

tic
[W,centers]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
centers=centers(N_unobs+1:end,:); % remove unobserved part

% W=fliplr(flipud(W));

RunningTime.GetWeights=toc;
W=sparse(W);
    if ~isempty(target_rates)
        bias=SetBiases(W,target_rates,params.spike_gen);
    end

if isempty(params.stat_flags.est_spar)
    W_nostim=W(1:N,1:N);
    params.stat_flags.est_spar=nnz(W_nostim(eye(N)<0.5))/(N^2-N); %correct sparsity estimation to actual sparsity
    params.conn_est_flags.est_spar=nnz(W_nostim(eye(N)<0.5))/(N^2-N); 
%     params.stat_flags.est_spar=nnz(W_nostim)/(N^2-N); %correct sparsity estimation to actual sparsity
%     params.conn_est_flags.est_spar=nnz(W_nostim)/(N^2-N); 
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
addpath(genpath('CalciumMeasurements'));

[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration,CalciumObs]=v2struct(params.spike_gen);

memory_threshold=1e8;  %maximial size of array to allow
splits=ceil((T*N^2)/memory_threshold)
mY=0; mYn=0;
XX=0; XXn=0;
XY=0; XYn=0;
t_start=0;
RunningTime.SampleSpikes=0;
RunningTime.GetSpikes=0;
RunningTime.GetStat=0;

for iter=1:splits
    disp(['iter=' num2str(iter) '/' num2str(splits)])
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
    true_spikes=GetSpikes(W,bias,T_split,T0,seed_spikes+iter,neuron_type,N_stim,stim_type,timescale,s0,verbose);
    RunningTime.GetSpikes=RunningTime.GetSpikes+toc;
    
    if CalciumObs==1        
        [ Y,spikes,relative_std_cell] = Spikes2Calcium2Spikes( true_spikes, params.calcium_model);

    else
        spikes=true_spikes;
        Y=spikes*0;
    end
    
    tic
    observations=SampleSpikes(N,T_split,sample_ratio,sample_type,obs_duration,N_stim,seed_sample+iter,t_start);
    if strcmp(params.connectivity.conn_type,'realistic+1') %"+1" neuron is always observed
        observations(N,:)=1;
    end
    t_start=t_start+T_split;
    sampled_spikes=sparse(observations.*spikes);
    RunningTime.SampleSpikes=RunningTime.SampleSpikes+toc;

    %% Estimate sufficeint statistics - previous version (without chopping up data series)
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
        
%% Estimate sufficeint statistics 
        tic
        sampled_spikes_obs=sampled_spikes(N_unobs+1:end,:);
        observations_obs=observations(N_unobs+1:end,:);
        mY=mY+full(sum(sampled_spikes_obs,2));
        mYn=mYn+full(sum(observations_obs,2));
        XX=XX+sampled_spikes_obs*sampled_spikes_obs';
        XXn=XXn+observations_obs*observations_obs';
        XY=XY+sampled_spikes_obs(:,1:(end-1))*(sampled_spikes_obs(:,2:end))';
        XYn=XYn+observations_obs(:,1:(end-1))*(observations_obs(:,2:end))';
        RunningTime.GetStat=RunningTime.GetStat+toc;
        
        imagesc(XYn);
        pause(1e-6);
    
end

rates=mY./(mYn+eps); %estimate the mean firing rates
Cxx=full(XX./(XXn+eps))-rates*rates'; %estimate the covariance (not including stim terms for now)
Cxy=full(XY./(XYn+eps))-rates*rates'; %estimate cross-covariance

% addpath('EstimateStatistics')
%  [Cxx,Cxy ] = PosProj( Cxx,Cxy );

%% Estimate Connectivity
addpath('EstimateConnectivity');
[pen_diag,pen_dist,warm,est_type,est_spar]=v2struct(params.conn_est_flags);
W_obs=W(N_unobs+1:end,N_unobs+1:end);     %observed W   
bias_obs=bias(N_unobs+1:end);  
quality=[];
error_rates=[];

tic

EW3=Cxy'/Cxx;
% 
% if strncmpi(neuron_type,'logistic',8)
%     [amp, Ebias2]=logistic_ELL(rates,EW3,Cxx,Cxy);    
% else
    Ebias2=bias*0;
    amp=ones(size(EW3,1),1);
% end
Ebias=Ebias2;
EW3=diag(amp)*EW3;
EW2=EW3;
EW=EW3;


switch est_type
%     case 'Linear'
        %         EW=EstimateA_L1(Cxy,Cxy,est_spar);
    case 'Gibbs'
        diag_mat=eye(N_obs)>0.5;
        p_0=est_spar*ones(N_obs);
        mu_0=zeros(N_obs);
        std_0=std(W_obs(find(~~W_obs(~diag_mat))))*ones(N_obs);   %#ok
        
        if ~pen_diag
            p_0(diag_mat)=1;
            mu_0(diag_mat)=mean(W_obs(diag_mat));
            std_0(diag_mat)=std(W_obs(diag_mat));
        end
        
        [EW,  quality] = EstimateA_Gibbs( bias_obs,sampled_spikes_obs,observations_obs,p_0, mu_0, std_0,rates); %quality here is accept array
        Ebias=0*bias_obs;%GetBias( EW,Cxx,rates);
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

   [EW,Ebias,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,pen_dist,warm,W_obs,centers);     
    EW2=EW3;
              
%         EW3=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);              
%         [amp, Ebias2]=logistic_ELL(rates,EW3,Cxx,Cxy);
%         EW2=diag(amp)*EW3;
%         EW2=EstimateA_MLE_cavity(Cxx,Cxy,rates);        %MLE
%         mask=~~EW;
      
%         EW2=EW2.*mask;
%         [amp, Ebias2]=logistic_ELL(rates,EW2,Cxx,Cxy);
%         EW2=diag(amp)*EW2;
    case 'FullyObservedGLM'
%         temp=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);      
%         [amp, Ebias2]=logistic_ELL(rates,temp,Cxx,Cxy);
%         EW=diag(amp)*temp;
        EW=EstimateA_L1_logistic_sampling(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);             
        [EW2, Ebias]=EstimateA_L1_logistic_fullyobserved(Cxx,Cxy,rates,spikes,est_spar,N_stim,pen_diag,warm);
     case 'EM'     
         tic
         [EW,Ebias,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,pen_dist,warm,W_obs,centers); 
         RunningTime.Cavity=toc;
         sampled_spikes_obs=full(sampled_spikes_obs);
            
        %Sample missing spikes       
        if sample_ratio<1
            for pp=1:2                
                spikes0=sampled_spikes_obs;
                if pp==1
                    EW2=EW;Ebias2=Ebias;
                    temp=repmat(rates,1,T)>rand(N_obs,T); 
                    spikes0(observations_obs==0)=temp(observations_obs==0);
                end
                sampled_spikes_obs(~observations_obs)=0.5;
                tic
                [sampled_spikes_obs, ~]=Estimate_rates_Gibbs_Metropolized(EW2,Ebias2,sampled_spikes_obs,spikes0);
                RunningTime.GibbsSpikes(pp)=toc;
                
%             [EW_out, Ebias2]=EstimateA_L1_logistic_fullyobserved(Cxx,Cxy,rates,sampled_spikes_obs,est_spar,N_stim,pen_diag,warm);
                XX=sampled_spikes_obs*sampled_spikes_obs';                
                XY=sampled_spikes_obs(:,1:(end-1))*(sampled_spikes_obs(:,2:end))';            
                mY=full(sum(sampled_spikes_obs,2));        
                rates=mY/T; %estimate the mean firing rates
                
                Cxx=full(XX/T)-rates*rates'; %estimate the covariance (not including stim terms for now)
                Cxy=full(XY/T)-rates*rates'; %estimate cross-covariance
%                 [Cxx,Cxy ] = PosProj( Cxx,Cxy );
                 
                [EW_out,Ebias2,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,pen_dist,warm,W_obs,centers); 
 
                if pp==1
                    EW2=EW_out;
                elseif pp==2                   
                    EW3=EW_out;
                end
            end
        else 
            tic
                [EW2, Ebias2]=EstimateA_L1_logistic_fullyobserved(Cxx,Cxy,rates,sampled_spikes_obs,est_spar,N_stim,pen_diag,warm);
            RunningTime.logistic_fullyobserved=toc;
        end

    otherwise 
        error('unknown inference method')
end

% correct weight scaling if model is LIF
 if strncmpi(neuron_type,'LIF',8)            
    EW=bsxfun(@times,EW,(std(W_obs(~~W_obs(:))/std(EW_obs(~~EW_obs(:))))));             
 end

%% Save Results
% Remove stimulus parts
if N_stim>0
    W=W(1:N,1:N);
    EW=EW(1:N_obs,1:N_obs);
    EW2=EW2(1:N_obs,1:N_obs);
    EW3=EW3(1:N_obs,1:N_obs);
    Ebias=Ebias(1:N_obs);
    Ebias2=Ebias2(1:N_obs);
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

% Spike_rec_correlation=zeros(size(spikes,1),1);
% for ii=1:size(spikes,1)
%     Spike_rec_correlation(ii)=corr(spikes(ii,:)',true_spikes(ii,:)');
% end
dt=0.01;
bin_min=1;bin_max=30;
DF = bin_min:bin_max;
timebins=DF*dt;
Spike_rec_correlation = zeros(N_obs,length(DF));
for nn = 1:N_obs
    [Spike_rec_correlation(nn,:),~] = ca_metrics(squeeze(spikes(nn,:)'),true_spikes(nn,:)',DF);
end

sample_time=1:min(T_split,1e4);
sample_traces.Y=Y(:,sample_time);
sample_traces.spikes=spikes(:,sample_time);
sample_traces.true_spikes=true_spikes(:,sample_time);
params.connectivity.N=N_obs;
params.RunningTime=RunningTime;

file_name=GetName(params);  %need to make this a meaningful name
save(file_name,'W_full','bias_full','W','bias','centers','EW','EW2','EW3','quality','error_rates'...
    ,'Spike_rec_correlation','timebins','sample_traces','Cxx','Cxy','rates','params'); %,'Ebias','Ebias2'?


end
%% Plot
Plotter