clear all
close all
clc


% Generate all figures for the paper
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
h=figure(1001);
% set(h,'units','normalized','outerposition',[0.3859    0.0620    0.6042    0.7593])

units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[2   2   18    12])
set(gcf, 'Position',[2   2   18   12])


dt=1e-2; %100 Hz imaging frame rate

N=50;
T=2e6;
observations_ratios=[1 0.2 0.1];
title_pos=[-0.15 1.1];
title_pos2=[-0.025 1.1];
title_pos3=[-0.07 1.03];
title_pos4=[-0.07 1.08];


load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(1)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'W','EW', 'Spike_rec_correlation','sample_traces');

a=4; b=12;
%% spike reconstructions
subplot = @(p) subtightplot (a,b, p, [0.13 0.08], [0.08 0.06], [0.07 0.01]);

% T_view=size(sample_traces.Y,2);
T_view=4e3;
tt=(1:T_view)*dt;
ii=6;
subplot([1 2 3 4 5 6 7 8 9])
plot(tt,sample_traces.Y(ii,1:T_view)/max(sample_traces.Y(ii,1:T_view))+0.2,'k')
xlim([0 max(tt)])
hold all
title('(A)','color', 'k',  'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold')

plot(tt,sample_traces.spikes(ii,1:T_view),'ro')%s,'MarkerEdgeColor','r','MarkerFaceColor','r')
hold all
plot(tt,sample_traces.true_spikes(ii,1:T_view),'.b')

ylim([0.1 1.1]);
xlbl=xlabel('time [sec]','fontsize',fontsize,'units','normalized');
xpos=get(xlbl,'position');
set(xlbl,'units','normalized','position',xpos+[0.05 0.1 0])
set(gca,'ytick',[]);
% title('A','color', 'k', 'fontweight', 'bold', 'Units','normalized', 'interpreter','none','position',[-0.03 0.95],'fontsize',fontsize)

% spike reconstruction correlations
T_view=size(sample_traces.Y,2);
N_view=size(sample_traces.Y,1);
tt=(1:T_view)*dt;
ii=6;
subplot([10 11 12])
plot(Spike_rec_correlation,'k')

title('(B)','color', 'k',  'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')
xlabel('neuron #','fontsize',fontsize)
ylabel('correlation','fontsize',fontsize);
xlim([1 N_view]);
% title('B','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)

%% Matrix reconstruction
subplot = @(p) subtightplot (a,b, p, [0.09 0.03], [0.05 0.06], [0.07 0.01]);

for ii=1:length(observations_ratios)
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'EW','W');
    EW_c=EW;
    W_c=W;

    [R_c,correlation_c, zero_matching_c,sign_matching_c] = GetWeightsErrors( W,EW_c ); 


load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_noCalcium.mat'],'W','EW');

subplot([1 2 3 4 1+b 2+b 3+b 4+b]+(ii-1)*(b/3)+b)
mi=min([W(:); EW(:)])*1.2;
ma=max([W(:); EW(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(W(:),EW(:),'.b');
hold all
plot(W_c(:),EW_c(:),'.r');
axis([mi ma mi ma])
xlabel('true weights','fontsize',fontsize)
if ii==1
    ylabel('inferred weights','fontsize',fontsize);
end
text(0.05, 0.9,['p_{\rm{obs}}= ' num2str(observations_ratios(ii))],'Units', 'normalized','fontsize',fontsize+2)
letter=['(' char(ii+1+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos3,'fontweight','bold')

subplot([1 2 3 4 ]+(ii-1)*(b/3)+3*b)
[R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW ); 
h=bar(1:4,[ [R,correlation, zero_matching,sign_matching];  [R_c,correlation_c, zero_matching_c,sign_matching_c]]',2);    
set(h(1),'FaceColor','b');
set(h(2),'FaceColor','r');
ylim([0 1])
xlim([0.5 4.5])
if ii==1
    ylabel('quality','fontsize',fontsize)
else
    set(gca, 'YTick', []);
end
set(gca, 'XTickLabel', x_ticks);
letter=['(' char(ii+4+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos4,'fontweight','bold')
% title('F','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)
end
% target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2';
% Export2Folder(['Fig10.eps'],target_folder) 