clear all
% close all
clc

% Generate stim figure
addpath('Results')
addpath('Misc')
addpath('GenerateSpikes');
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',14)
set(0,'defaultlinelinewidth',2.5)
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.08], [0.05 0.05], [0.1 0.01]);

%%
T=1e5;
N=50;
p_obs=0.2;

%% Figure 1 - Toy model 

L=9; %number of stimuli + 1
K=5; %width of subplots
x_ticks={'R','C','Z','S'};
fontsize=12;
fontsize2=1.3*fontsize; 
if N<100
    dot_size=100; %for scatter plots
else
    dot_size=0.01;
end
    
% load(['Run_N=50_obs=' num2str(observations_ratios(1)) '_T=200000.mat']);
% neuron_type='logistic'  ;
% spikes=GetSpikes(W,bias,params.spike_gen.T,params.spike_gen.T0,params.spike_gen.seed_spikes,neuron_type,params.spike_gen.N_stim,params.spike_gen.stim_type);

figure(1)

%%
corr_array=zeros(L,1);
mean_rate=zeros(L,1);

for ii=1:L
    
    if ii==1
        load(['Run_N=' num2str(N) '_obs=' num2str(p_obs) '_T=' num2str(T) '.mat'],'W','EW2','rates','params');
    else
        load(['Run_N=' num2str(N) '_obs=' num2str(p_obs) '_T=' num2str(T) '_N_stim=' num2str(ii-1) '.mat'],'W','EW2','rates','params');
    end
%         Run_N=50_obs=0.2_T=100000_N_stim=1
    
%     mi=min([ W(:); EW2(:)] );ma=max([W(:); EW2(:)]);
%     mi2=min([ EW2(:)] );ma2=max([ EW2(:)]);
    mi=-1;ma=1;
    mi2=mi; ma2=ma;
    subplot(L+1,K,[1 2])
    imagesc(W,[mi ma]); h=colorbar;
%     colormap('gray')
    ylabel('True W','fontsize',fontsize2)
    set(h, 'ylim', [mi ma])
    subplot(L+1,K,[3 4])
    scatter(W(:),W(:),dot_size,'b.')
    xlabel('W','fontsize',fontsize2)
    ylabel('W','fontsize',fontsize2)
    axis([mi ma mi2 ma2])
    subplot(L+1,K,5)
    y_ticks=[1,1,1,1];   
    bar(y_ticks);    
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
    
    subplot(L+1,K,K*ii+[3 4])
    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind,'r-');
    hold all
    scatter(W(:),EW2(:),dot_size,'b.')
    axis([mi ma mi2 ma2])
    hold off
%     legend('x=y','EW2','EW22')
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    
   subplot(L+1,K,K*ii+5)
     [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW2 );
     corr_array(ii)=correlation; 
     mean_rates(ii)=mean(rates);      
%     [MSE2,correlation2,SE2] = GetWeightsErrors( W,EW22 );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
%     title({[' EW2 corr =' num2str(correlation) ', EW22 corr =' num2str(correlation2)]; ...
%          [' EW2 MSE =' num2str(MSE) ', EW22 MSE =' num2str(MSE2)]; ...
%          [' EW2 SE =' num2str(SE) ', EW22 SE =' num2str(SE2) ]});
%     hold off

    
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(EW2,[mi ma]); h=colorbar;
    set(h, 'ylim', [mi ma])
    if ii==1    
        ylabel('No stimulus','fontsize',fontsize2);
    else
        ylabel('With stimulus','fontsize',fontsize2);
    end

    
end
%%
figure(2)
fontsize=17;
stim_strength=(1:L)-1;
[AX,H1,H2] =plotyy(stim_strength,corr_array,stim_strength,mean_rates);
xlabel('Stimulus amplitude','fontsize',fontsize)
 set(get(AX(1),'Ylabel'),'String','C (correlation)','fontsize',fontsize) 
 set(get(AX(2),'Ylabel'),'String','m (mean spike probability)','fontsize',fontsize) 

%%
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript'
Export2Folder(['Stim_N=' num2str(N) '.eps'],target_folder) 