clear all
close all
clc

addpath(genpath('CalciumMeasurements')); %adds folder with subfolders
data_types={'Tolias','Tim','Manolis','Sim'};
data=data_types(4);
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
 spikes=GetSpikes(W,bias,T,T0,seed_spikes+iter,neuron_type,N_stim,stim_type,timescale,s0,verbose); 
spikes=spikes(1:N,:);
RunningTime.GetSpikes=toc;

tic
sn=0.2;
amp=1;
b=3; 
g=0.95;
noise_mode=0;
Y = Spikes2Calcium(spikes,g,b,amp,sn,noise_mode);
tt=1:length(Y);
g=[g ; zeros(order_sysid-1,1)];
dt=1;
end
%% Use real data
data_sets_num=1; %for Tim's data
K=6; %=size(spike_est,3)
CR_cell=cell(data_sets_num,K); 
KL_cell=CR_cell;

for pp=1:data_sets_num

if strcmp(data,'Tolias')    
    target_folder='C:\Users\Daniel\Copy\Columbia\Research\Data\ToliasWithEphys';
    load(fullfile(target_folder,'data_before_filtering.mat'),'data');
    % tt=1:length(data{ind}.calcium);
    tt=(1:2e3)+2e4;        
    Y=double(data{pp}.calcium(tt)); %ind in 1 to 3
    spikes=double(data{pp}.spikes(tt));
    dt=1/data.fps;
elseif strcmp(data,'Tim')   
     [ Y,spikes,dt ] = GetTimData(pp); %ind in 1 to..?
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
P_arpfit=GetParams(Y,order_arpfit,'keep','arpfit');
P_psd= GetParams(Y,order_arpfit,'psd','arpfit');
P_sysid = GetParams(Y,order_sysid,'keep','sysid');
P_sysid2 = GetParams(Y,order_sysid+1,'keep','sysid');

[spikes_arpfit,b_arpfit] = Calcium2Spikes(Y,P_arpfit);
[spikes_psd,b_psd] = Calcium2Spikes(Y,P_psd);
[spikes_sysid,b_sysid] = Calcium2Spikes(Y,P_sysid);
[spikes_sysid2,b_sysid2] = Calcium2Spikes(Y,P_sysid2);

Reps=3;
spikes_sysid_iterated=spikes_psd;
b_sysid_iterated=b_psd;
for ii=1:Reps
    P_sysid_iterated = GetParams(Y,order_sysid,'keep','sysid',[spikes_sysid_iterated, b_sysid_iterated]);
    [spikes_sysid_iterated,b_sysid_iterated] = Calcium2Spikes(Y,P_sysid_iterated);
end

spikes_sysid_iterated2=spikes_psd;
b_sysid_iterated2=b_psd;
for ii=1:Reps
    P_sysid_iterated2 = GetParams(Y,order_sysid+1,'keep','sysid',[spikes_sysid_iterated2, b_sysid_iterated2]);
    [spikes_sysid_iterated2,b_sysid_iterated2] = Calcium2Spikes(Y,P_sysid_iterated2);
end

% P_sysid2 = GetParams(Y,order_sysid,'keep','sysid',[spikes, [b;b]]);
spike_est=zeros(N,T,4);
spike_est(:,:,1)=spikes_psd;
spike_est(:,:,2)=spikes_arpfit;
spike_est(:,:,3)=spikes_sysid;
spike_est(:,:,4)=spikes_sysid_iterated;
spike_est(:,:,5)=spikes_sysid2;
spike_est(:,:,6)=spikes_sysid_iterated2;

RunningTime.GetSpikes=toc;

%% Plot estimated params
if plot_stuff
%     ind=1;
%     P_arpfit{ind}.g=[P_arpfit{ind}.g ; zeros(order_sysid-order_arpfit,1)];

    sn_vect=[sn,P_arpfit{ind}.sn,P_psd{ind}.sn,P_sysid{ind}.sn,P_sysid_iterated{ind}.sn];
    g_vect=[g,P_arpfit{ind}.g,P_sysid{ind}.g,P_sysid_iterated{ind}.g]';
    b_vect=[b,P_arpfit{ind}.Cb,P_psd{ind}.Cb,P_sysid{ind}.Cb,P_sysid_iterated{ind}.Cb]';
    roots_vec=g_vect*0;
    for ii=1:size(g_vect,1)
        if ~any(g_vect(ii,:))
            roots_vec(ii,:)=g_vect(ii,:)*0;
        else
            roots_vec(ii,:)=roots([1 -g_vect(ii,:)]);
        end
    end
    figure(1)
    subplot(4,1,1)
    bar(sn_vect);
    set(gca,'XTickLabel',{'true', 'arpfit','psd', 'sysid','sysid_iterated'})
    ylabel('Sn')
    subplot(4,1,2)
    bar(g_vect);
    ylabel('g')
    set(gca,'XTickLabel',{'true', 'arpfit', 'sysid','sysid_iterated'})
    subplot(4,1,3)
    bar(abs(roots_vec));
    ylabel('g roots')
    set(gca,'XTickLabel',{'true', 'arpfit', 'sysid','sysid_iterated'})
    subplot(4,1,4)
    bar(b_vect);
    ylabel('b')
    set(gca,'XTickLabel',{'true', 'arpfit','psd', 'sysid','sysid2'})

    roots_vec
    %% plot
    xmin=min(tt);
    xmax=max(tt);
    % xmin=1e3;
    % xmax=2e3;
    a=7;
    b=1;

    for ind=1:N
        figure(1+ind)
        subplot(a,b,1)
        Gx =G_mat(spikes_psd(ind,:),1,T,P_psd{ind}.g,0,P_psd{ind}.z)+b_psd(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_psd(ind,:)/max(spikes_psd(ind,:)+eps));    
        ylabel('psd')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,2)
        Gx =G_mat(spikes_arpfit(ind,:),1,T,P_arpfit{ind}.g,0,P_arpfit{ind}.z)+b_arpfit(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_arpfit(ind,:)/max(spikes_arpfit(ind,:)+eps));    
        ylabel('arpfit')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,3)
        Gx =G_mat(spikes_sysid(ind,:),1,T,P_sysid{ind}.g,0,P_sysid{ind}.z)+b_sysid(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_sysid(ind,:)/max(spikes_sysid(ind,:)+eps));    
        ylabel('sysid')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,4)
        Gx =G_mat(spikes_sysid_iterated(ind,:),1,T,P_sysid_iterated{ind}.g,0,P_sysid_iterated{ind}.z)+b_sysid_iterated(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_sysid_iterated(ind,:)/max(spikes_sysid_iterated(ind,:)+eps));    
        ylabel('sysid_iterated')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,5)
        Gx =G_mat(spikes_sysid2(ind,:),1,T,P_sysid2{ind}.g,0,P_sysid2{ind}.z)+b_sysid2(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_sysid2(ind,:)/max(spikes_sysid2(ind,:)+eps));    
        ylabel('sysid2')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,6)
        Gx =G_mat(spikes_sysid2(ind,:),1,T,P_sysid2{ind}.g,0,P_sysid2{ind}.z)+b_sysid2(ind);
        plot(tt,Y(ind,:),tt,Gx,tt,spikes_sysid2(ind,:)/max(spikes_sysid2(ind,:)+eps));    
        ylabel('sysid2_iterated')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);
        subplot(a,b,7)    
        plot(tt,Y(ind,:),tt,spikes(ind,:)/max(spikes(ind,:)));    
        ylabel('true spikes')
        xlabel('t')
        xlim([xmin,xmax]);
        ylim([min(Y(ind,:)) max(Y(ind,:))]);

    target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Code\Shotgun\Figures';
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

        if strcmp( data, 'Tim')
            Export2Folder([' Machado - neuron #' num2str(ind) '.pdf'],target_folder)    
        elseif strcmp( data, 'Tolias')
            Export2Folder([' Tolias - neuron #' num2str(ind) '.pdf'],target_folder)  
        elseif strcmp( data, 'Sim')
            Export2Folder([' Simulation - neuron #' num2str(ind) '.pdf'],target_folder)              
        end

    end
end

%% 

    bin_min=1;bin_max=50;
    DF = bin_min:bin_max;
for nn=1:size(spike_est,3)        
    CR = zeros(size(N,1),length(DF));
    KL = CR;
    for ii = 1:N
        [CR(ii,:),KL(ii,:)] = ca_metrics(squeeze(spike_est(ii,:,nn))',spikes(ii,:)',DF);
    end
%     CR(isnan(CR))=0;
%     KL(isnan(CR))=1;
    CR_cell(pp,nn)={CR};
    KL_cell(pp,nn)={KL};
    

end

end
%%
CR_mat=cell2mat(CR_cell(:,1));
KL_mat=cell2mat(KL_cell(:,1));
  
for ii=2:size(CR_cell,2)
        CR_mat(:,:,end+1)=cell2mat(CR_cell(:,ii));
        KL_mat(:,:,end+1)=cell2mat(KL_cell(:,ii));
end

%%
figure(1000)
a=2;
b=1;
method_names={'psd (p=1)','arpfit (p=1)','sysid (z=1,p=1)','iterated sysid (z=1,p=1)','sysid (z=2,p=2)',...
        'iterated sysid (z=2,p=2)'};
    
timebins=DF*dt;
for nn=1:size(spike_est,3)      
    subplot(a,b,1)
%     colorOrder = get(gca, 'ColorOrder');
%     current_color=colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :);
%     shadedErrorBar(DF*dt,mean(CR,1),std(CR,1)/sqrt(N),{'color', current_color},1)
    errorbar(timebins,mean(CR_mat(:,:,nn),1),std(CR_mat(:,:,nn),1)/sqrt(size(CR_mat,1)))
    ylabel('Correlation')
     xlabel('Time bin size [sec]')   
    xlim([bin_min bin_max]*dt)
    hold all   
    subplot(a,b,2)
    errorbar(timebins,mean(KL_mat(:,:,nn),1),std(KL_mat(:,:,nn),1)/sqrt(size(CR_mat,1)))
    ylabel('-KL divergence')
    xlabel('Time bin size [sec]')    
    xlim([bin_min bin_max]*dt)
    hold all

    
        legend(method_names,'location','southeast')
    set(gcf,'units','normalized','outerposition',[0 0 0.4 0.8])
        target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Code\Shotgun\Figures';
%         if strcmp( data, 'Tim')
%             Export2Folder([' Machado - Quality Measures.pdf'],target_folder)    
%         elseif strcmp( data, 'Tolias')
%             Export2Folder([' Tolias - Quality Measures.pdf'],target_folder)  
%         elseif strcmp( data, 'Sim')
%             Export2Folder([' Simulation - Quality Measures.pdf'],target_folder)              
%         end
end

a=2;
b=3;
%%
for nn=1:size(spike_est,3)      
        subplot(a,b,nn)
    imagesc(timebins,1:size(spike_est,3),CR_mat(:,:,nn),[0 1])        
   xlabel('Time bin size [sec]')    
   ylabel('Trials')    
    title(method_names{nn})
%         if strcmp( data, 'Tim')
%             Export2Folder([' Machado - Correlations.pdf'],target_folder)    
%         elseif strcmp( data, 'Tolias')
%             Export2Folder([' Tolias - Correlations.pdf'],target_folder)  
%         elseif strcmp( data, 'Sim')
%             Export2Folder([' Simulation - Correlations.pdf'],target_folder)              
%         end
end
for nn=1:size(spike_est,3)      
    subplot(a,b,nn)
    imagesc(timebins,1:size(spike_est,3)  ,KL_mat(:,:,nn),[-5 0]) 
    title(method_names{nn})
%             if strcmp( data, 'Tim')
%             Export2Folder([' Machado - -KL.pdf'],target_folder)    
%         elseif strcmp( data, 'Tolias')
%             Export2Folder([' Tolias - -KL.pdf'],target_folder)  
%         elseif strcmp( data, 'Sim')
%             Export2Folder([' Simulation - -KL.pdf'],target_folder)              
%             end
   xlabel('Time bin size [sec]')    
   ylabel('Trials')    
end

% 
% 
% save('results','CR_cell', 'KL_cell');