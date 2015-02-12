clear all
close all
clc

load_data=1;
sample_type_set={'continuous','fixed_subset','spatially_random','prob'};
sample_type=sample_type_set{3};

if load_data
    load('zebrafish_data_small','data');
else 
    N=50;
    subsets={'all','front','back'};
    subset=subsets(2);

    addpath(genpath('CalciumMeasurements'))
    [data_original, cent_merged] = GetZebraFishData(subset,N);

    % Pre-processing
    data=data_original(1:end-1,:);
    % data=bsxfun(@times,data,1./mean(data,2));
    % data=diff(data,[],2);
    order=1;
    P= GetParams(data,order,'psd','arpfit');
    [data,b] = Calcium2Spikes(data,P);
    data=bsxfun(@times,data,1./max(data,2));

    components=0; %number of NN components to remove
    for kk=1:components
        [v_s,v_t] = GreedyNNPCA(data,[]);
        data=data-v_s*v_t;
    end

    data=[data ; data_original(end,:)];
end

[N,T]=size(data); % Includes the stimulus now. Also the original N decreases if it is larger then subset

%% observations
sample_ratio_array=[0.1,0.2,0.4,0.8,1];
est_spar_array=[0.05,0.1,0.2,0.4,0.8,1];
AUROC_cell=cell(length(sample_ratio_array),length(est_spar_array));
EW_cell=cell(length(sample_ratio_array),length(est_spar_array));
Cxx_cell=cell(length(sample_ratio_array),1);
Cxy_cell=cell(length(sample_ratio_array),1);
actual_sparsity=zeros(length(sample_ratio_array),length(est_spar_array));

for pp=1:length(sample_ratio_array)
    sample_ratio=sample_ratio_array(pp);
    N_stim=0;
    seed_obs=1;
    addpath('GenerateSpikes')
    observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_obs);
    full_observations=ones(size(data))>0.5;
    sampled_data=observations.*data;
    %% estimate stats
    addpath('EstimateStatistics')
    glasso=0;
    restricted_penalty=1;
    pos_def=1;
    est_spar=1;
    W=0;

    [Cxx, Cxy,~,rates,obs_count] = GetStat(sampled_data,observations,glasso,restricted_penalty,pos_def,est_spar,W);
%     [Cxx_full, Cxy_full,~,rates_full,obs_count_full] = GetStat(data,full_observations,glasso,restricted_penalty,pos_def,est_spar,W);
    Cxx_cell(pp)={Cxx};
    Cxy_cell(pp)={Cxy};
    %% estimate connectivity
    addpath('EstimateConnectivity')
    pen_diag=0;
    warm=1;
    
    for kk=1:length(est_spar_array)
        est_spar=est_spar_array(kk);

        EW=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);              
        % [amp, Ebias]=logistic_ELL(rates,EW,Cxx,Cxy);
        % EW=diag(amp)*EW;
        
        EW_cell(pp,kk)={EW};    

        %% 
    end
end

addpath('Misc')
for pp=1:length(sample_ratio_array)
            Cxy=Cxy_cell{pp};
            Cxx=Cxx_cell{pp};            
    for kk=1:length(est_spar_array)
        if (est_spar_array(kk)==1)&(sample_ratio_array(pp)==1)
            EW=Cxy'/Cxx;
        else
            EW=EW_cell{pp,kk};
        end
        actual_sparsity(pp,kk)=nnz(EW(:))/numel(EW(:));
        
        Ebias=0;
        predicted_data=1./(1+exp(-bsxfun(@plus,EW*data,Ebias)));
        AUROC=zeros(N-1,1);
        for ii=1:N-1
            AUROC(ii) = GetAUROC(predicted_data(ii,1:end-1),data(ii,2:end));
        end
         AUROC_cell(pp,kk)={mean(AUROC)};
    end
end

save('zebrafish_results','EW_cell','AUROC_cell','Cxx_cell','Cxy_cell','actual_sparsity');
%%
% figure(1)
% subplot(211)
% % EW_full=Cxy_full'/Cxx_full;
% EW_full_visualize=EW_full(1:end-1,1:end-1);
% ma=max(EW_full_visualize(:));mi=min(EW_full_visualize(:));
% imagesc(EW_full_visualize,[mi ma])
% title('EW full')
% subplot(212)
% % EW=Cxy'/Cxx;
% EW_visualize=EW(1:end-1,1:end-1);
% imagesc(EW_visualize,[mi ma])
% title('EW')

figure(2000)
AUROC_mat=cell2mat(AUROC_cell);
imagesc(AUROC_mat,[0.5 1]);
set(gca,'XTick',1:length(est_spar_array),'XTickLabel',est_spar_array,...
    'YTick',1:length(sample_ratio_array),'YTickLabel',sample_ratio_array);
ylabel('sample ratio')
xlabel('target sparsity')
colorbar
target_folder='Figures';
Export2Folder([' zebrafish_AUROC.pdf'],target_folder)          

figure(2001)
plot(est_spar_array,actual_sparsity)
%% 
figure(1)
a=5; b=1;
dt=1/2.1; %2.1Hz
tt=(1:T)*dt;
subplot(a,b,[1 3])
imagesc(tt,1:N,data(1:end-1,:))
colormap('gray')
% xlabel('t [sec]')
ylabel('cell #')
% colorbar('location','northoutside')
subplot(a,b,4)
plot(tt,mean(data(1:end-1,:),1))
xlim([tt(1) tt(end)])
% xlabel('t [sec]')
title('mean population rate')
subplot(a,b,5)
plot(tt,data(end,:))
xlim([tt(1) tt(end)])
xlabel('t [sec]')
title('mean unobserved population rate')
target_folder='Figures';
Export2Folder([' zebrafish_data.pdf'],target_folder)          


    % figure(3)
    % subplot(221)
    % imagesc(Cxx_full)
    % title('full Cxx')
    % subplot(222)
    % imagesc(Cxy_full)
    % title('full Cxy')
    % subplot(223)
    % imagesc(Cxx)
    % title('Cxx')
    % subplot(224)
    % imagesc(Cxy)
    % title('Cxy')