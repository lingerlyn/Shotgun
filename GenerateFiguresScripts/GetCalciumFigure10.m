clear all
close all
clc


% Generate all figures for the paper
% addpath('Results\Calcium')
addpath('Results')
addpath('Misc')
addpath('Plotting')
% addpath('GenerateSpikes');
% addpath('GenerateConnectivity')
% set(0,'DefaultTextInterpreter', 'latex');
% set(0,'DefaultAxesFontSize',13)
SetDefaultGraphicSettings( )

fontsize=12;
x_ticks={'R','C','Z','S'};
h=figure(1003);
% set(h,'units','normalized','outerposition',[0.3859    0.0620    0.6042    0.7593])

units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[2   2   18    15])
set(gcf, 'Position',[2   2   18   15])


dt=1e-2; %100 Hz imaging frame rate

N=50;
T=2e6;
observations_ratios=[1 0.2 0.1];
title_pos=[-0.15 1.1];
title_pos2=[-0.025 1.1];
title_pos3=[-0.07 1.03];
title_pos4=[-0.07 1.08];

a=4; b=12;
spike_mag=2; %height of spikes in subfigure A

color_set={[0    0    1],[0.9290    0.6940    0.1250],[0.6350    0.0780    0.1840]};
%% spike reconstructions
subplot = @(p) subtightplot (a,b, p, [0.12 0.08], [0.08 0.04], [0.07 0.01]);

for kk=1:2
    if kk==1
        load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(1)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'W','EW', 'Spike_rec_correlation','timebins','sample_traces');
    else
         load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(1)) '_T=' num2str(T) '_Cavity_CalciumObsHighNoise.mat'],'W','EW', 'Spike_rec_correlation','timebins','sample_traces');
    end


    ind=1:5e2;
%      ind=1:1e4;
    tt=ind*dt;
    ii=14;
    subplot([1 2 3 4 5 6 7 8 9])
    plot(tt,sample_traces.Y(ii,ind)/max(sample_traces.Y(ii,ind))+0.2+0.6*(1-(kk-1)),'color',color_set{kk+1},'linewidth',1);
    xlim([tt(1) tt(end)])
    hold all
    title('(A)','color', 'k',  'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold')


    hold all
    if kk==1
        marker='v';
        markersize=4;
        shift=1.05;
    else        
        marker='^';
        markersize=4;
        shift=0.95;
    end
%         temp=sample_traces.true_spikes(ii,ind).*(1-abs(sample_traces.spikes(ii,ind)-sample_traces.true_spikes(ii,ind)));
        plot(tt,shift*spike_mag*sample_traces.spikes(ii,ind),marker,'markersize',markersize,'color',color_set{kk+1})
        hold all     
    %s,'MarkerEdgeColor','r','MarkerFaceColor','r')
    plot(tt,spike_mag*sample_traces.true_spikes(ii,ind),'+','markersize',6,'color',color_set{1})
    ylim(0.05+1.1*[0, spike_mag]);
    
    xlbl=xlabel('time [sec]','fontsize',fontsize,'units','normalized');
%     xpos=get(xlbl,'position');
%     set(xlbl,'units','normalized','position',xpos+[0.05 0.1 0])
    set(gca,'ytick',[]);
    % title('A','color', 'k', 'fontweight', 'bold', 'Units','normalized', 'interpreter','none','position',[-0.03 0.95],'fontsize',fontsize)

    % spike reconstruction correlations
    T_view=size(sample_traces.Y,2);
    N_view=size(sample_traces.Y,1);
    tt=(1:T_view)*dt;

    subplot([10 11 12])
    mCR=mean(Spike_rec_correlation,1);
    sCR=std(Spike_rec_correlation,[],1);
    hold all
    shadedErrorBar(timebins*1e3,mCR,sCR,{'color',color_set{kk+1}},0)
    title('(B)','color', 'k',  'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')
    xlabel('bin size [msec]' ,'fontsize',fontsize)
    ylabel('correlation','fontsize',fontsize);
    xlim([timebins(1) timebins(end/3)]*1e3);
    set(gca,'xtick',[10, 30, 50, 70, 90]);
    % title('B','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)
end
%% Matrix reconstruction
subplot = @(p) subtightplot (a,b, p, [0.09 0.03], [0.05 0.06], [0.07 0.01]);

for ii=1:length(observations_ratios)
    
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'EW','W','params');
%     ind1=EW<0; ind2=EW>0;
%     g1=std(EW(ind1))/std(W(ind1));
%     g1=0.5;
%     EW(ind1)=EW(ind1)/g1;
    
    EW_c=EW;
    W_c=W;

    [R_c,correlation_c, zero_matching_c,sign_matching_c] = GetWeightsErrors( W_c,EW_c ); 
    
     load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_CalciumObsHighNoise.mat'],'EW','W','params');
%         ind1=EW<0; ind2=EW>0;
%     g1=0.25;
%     EW(ind1)=EW(ind1)/g1;
%     ind1=EW<0; ind2=EW>0;
%     g1=(1-mean(Spike_rec_correlation));
%     g2=1-g1;
%     EW(ind1)=EW(ind1)./g1;
%     EW(ind2)=EW(ind2)./g2;
    
% %      g1=norm(EW(ind1))/norm(W(ind1));
%     g1=std(EW(ind1))/std(W(ind1));
%     EW(ind1)=EW(ind1)/g1;
%     EW(ind1)=EW(ind1)+mean(W(ind1)-EW(ind1));
%     g2=std(diag(EW))/std(diag(W));
%     EW(eye(N)>0.5)=EW(eye(N)>0.5)/g2;
%     EW(eye(N)>0.5)=EW(eye(N)>0.5)+mean(diag(W-EW));
%     
    EW_c2=EW;
    W_c2=W;

    [R_c2,correlation_c2, zero_matching_c2,sign_matching_c2] = GetWeightsErrors( W_c2,EW_c2 ); 


load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_noCalcium.mat'],'W','EW');

subplot([1 2 3 4 1+b 2+b 3+b 4+b]+(ii-1)*(b/3)+b)
mi=min([W(:); EW(:)])*1.2;
ma=max([W(:); EW(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(W(:),EW(:),'.','color',color_set{1});
hold all
plot(W_c(:),EW_c(:),'.','color',color_set{2});
hold all
plot(W_c2(:),EW_c2(:),'.','color',color_set{3});
grid on
axis([mi ma mi ma])
xlabel('true weights','fontsize',fontsize)
if ii==1
    ylabel('inferred weights','fontsize',fontsize);
end
text(0.02, 0.8,['p_{\rm{obs}}= ' num2str(observations_ratios(ii))],'Units', 'normalized','fontsize',fontsize+2)
letter=['(' char(ii+1+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos3,'fontweight','bold')

subplot([1 2 3 4 ]+(ii-1)*(b/3)+3*b)
[R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW ); 
h=bar(1:4,[ [R,correlation, zero_matching,sign_matching]; [R_c,correlation_c, zero_matching_c,sign_matching_c];  [R_c2,correlation_c2, zero_matching_c2,sign_matching_c2] ]',2);    
set(h(1),'FaceColor',color_set{1});
set(h(2),'FaceColor',color_set{2});
set(h(3),'FaceColor',color_set{3});
ylim([0 1])
xlim([0.5 4.5])
set(gca, 'YTick', [0 0.25 0.5 0.75 1]);
if ii==1
    ylabel('quality','fontsize',fontsize)
else
    set(gca, 'YTickLabel', []);
end
set(gca, 'XTickLabel', x_ticks,'ygrid','on');
letter=['(' char(ii+4+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos4,'fontweight','bold')
% title('F','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)
end
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\FinalProduction';
Export2Folder(['Fig10.tif'],target_folder) 