clear all
% close all
clc

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('GenerateSpikes');
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',10)
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.07], [0.06 0.1], [0.08 0.01]);
T=2e6;
T_view=3e2;
N=50;
dt=1e-2; %100 Hz imaging frame rate

if N==50
    observations_ratios= [1,0.2,0.1,0.04];
    dot_size=3.8;
    ind_array=1:N;
elseif N==1e3;
    observations_ratios= [1,0.2,0.1];
    dot_size=3.8;
    M=50;
    ind_array=randi(N,M,1);
end

fontsize=10;
fontsize2=12; 
L_bins=30; %number of bins in the weight histogram
%% Figure 1 - Toy model 
% In this figure:
% N=50, 
% T=2e5 (~1 hour experiment, assuming 200ms bins)
% firing rate~0.2   (~1Hz firing rate)
% balanced network (no Dale's law yet)

L=length(observations_ratios);
K=7; %width of subplots
x_ticks={'R','C','Z','S'};

%% Plot Network statistics
figure(1)
load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(2)) '_T=' num2str(T) '_Cavity.mat'],'W','W_full','bias_full','EW','bias','rates','params','centers');
W_temp=W(ind_array,ind_array);
mi=min(W_temp(:));
ma=max(W_temp(:));
% mi=-1;ma=1;
mi2=min(W(:))*1.2;
ma2=max(W(:))*1.2;

subplot(L+1,K,[1 2])
imagesc(W(ind_array,ind_array),[mi ma]);
% h=colorbar;
% set(h, 'ylim', [mi ma])
ylabel('true W','fontsize',fontsize2)
freezeColors
% h=colorbar;
% cbfreeze(h)

subplot(L+1,K,[3 4])
[hist_W,bins]=hist(W(~~W(:)),L_bins);
plot(bins,hist_W)
xlabel('W');
ylabel('count');

addpath('GenerateSpikes');
subplot(L+1,K,[ 5 6 7])
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration]=v2struct(params.spike_gen);
s0=[];verbose=1;
spikes=GetSpikes(W_full,bias_full,T_view,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbose);
 
tt=(1:T_view)*dt; nn=1:length(bias_full);
imagesc(tt,nn,spikes)
% line()
colormap('gray')
freezeColors
xlabel('time [sec]');
ylabel('neuron');
line([-1 2*tt(end)],N*[1 1],'color','red')


%%
for ii=1:L
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity.mat'],'W','EW','rates');
        %     high_firing_rate=rates>0.05;
%     W=W(high_firing_rate,high_firing_rate);
%     EW=EW(high_firing_rate,high_firing_rate);
       
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(EW(ind_array,ind_array),[mi ma]);
%     h=colorbar;
%     set(h, 'ylim', [mi ma])
    ylabel(['$p_{\mathrm{obs}}=$' num2str(observations_ratios(ii))],'fontsize',fontsize2)
    colormap('jet')

    subplot(L+1,K,K*ii+[5 6])
    A_ind=linspace(mi2,ma2,100);
    plot(A_ind,A_ind,'r-');
    hold all
    EW_temp=EW(ind_array,:);
    W_temp=W(ind_array,:);
    plot(W_temp(:),EW_temp(:),'b.','MarkerSize',dot_size)
%     cloudPlot(W(:),EW(:),[mi ma mi ma],1,[30 30])
    axis([mi2 ma2 mi2 ma2])
    hold off
%     legend('x=y','EW','EW2')
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    subplot(L+1,K,K*ii+[3 4])
    mi3=-2;%min([W(:) ;EW(:)]);
    ma3=2;%max([W(:) ;EW(:)]);
    bins=linspace(mi3,ma3,L_bins);
%     hist_W=hist(W(~~W(:)),bins);
%     hist_EW=hist(EW(~~EW(:)),bins);
    hist_W=hist(W(:),bins);
    hist_EW=hist(EW(:),bins);
    semilogy(bins,hist_W,'b-',bins,hist_EW,'r.');
    xlim([mi3 ma3]);
    xlabel('W','fontsize',fontsize2)
    ylabel('count','fontsize',fontsize2)
%     if ii==1
%         h=legend('W','$\hat{W}$')
%         set(h,'Interpreter','latex','fontsize',fontsize,'location','northwest','orientation','horizontal')
%     end
    
   subplot(L+1,K,K*ii+7)
   [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    ylabel('quality','fontsize',fontsize)
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
    

end

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript'
Export2Folder(['Sparsity_N=' num2str(N) '.pdf'],target_folder) 