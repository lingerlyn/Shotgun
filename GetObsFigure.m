clear all
% close all
clc

% Generate all figures for the paper
addpath('GenerateSpikes');
addpath('Misc')

set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',11)
set(0,'defaultAxesFontName', 'Times')
set(0,'defaultTextFontName', 'Times')

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.1], [0.2 0.2], [0.1 0.1]);
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.08], [0.07 0.04], [0.05 0.01]);

title_height=1;
title_horz=-0.2;
title_font=11;
figure(1)
% set(h,'units','normalized','outerposition',[0 0 0.6 1])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',1.8*[0.1,1, 9, 12])
set(gcf, 'Position',1.8*[0.1,1, 9, 12])


N=10;
T=1e5*N;
t_show=1:3e3;
dt=0.01; %100 Hz frame rate
tt=t_show*dt;
sample_ratio=0.2;
N_stim=0;
sample_type_set={'continuous','fixed_subset','spatially_random','prob','double_continuous','fully_random'};
obs_duration=100;
t_start=1;

seed_sample=1;
L=100; %number of histogram bins
M=5; %show M cases
a=M; b=3; %subplot grid

for kk=1:M
switch kk
    case 1
        sample_type=sample_type_set{2};
        name='Fixed';
    case 2
        sample_type=sample_type_set{1};
        name='Serial';
    case 3
        sample_type=sample_type_set{6};
        name='Fullly Random';        
    case 4
        sample_type=sample_type_set{3};
        name='Random Blocks';
    case 5
        sample_type=sample_type_set{5};
        name='Doubled Serial';
end 

observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample,t_start);
% inputs:
subplot(a,b,(kk-1)*b+[1 2])
imagesc(tt,[],observations(:,t_show))
% title('(A)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
xlabel('t [sec]','fontweight','bold');
ylabel(name,'fontweight','bold');
% ylabel('cell number');
% title(name, 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
colormap('gray')
freezeColors

subplot(a,b,(kk-1)*b+3)
sample_type=sample_type_set{3};
% observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
XYn=observations(:,1:(end-1))*(observations(:,2:end))'/T;
ma=max(XYn(:));
imagesc(XYn,[0 ma]);
% title('(B)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
% title('results in all spike pairs being observed', 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
% xlabel('cell number');
% ylabel('cell number');
caxis([0 0.05]);
% ['random scanning' -   ; 'all spike pairs observed (good)']
hc=colorbar;
colorbar_labels = get(hc,'YTickLabel');
colorbar_labels=[colorbar_labels, [repmat(' ',1,size(colorbar_labels,1)-1), '<']'];
set(hc,'YTickLabel',colorbar_labels);
colormap('jet')
freezeColors

% subplot(a,b,(kk-1)*b+5)
% 
% bins=linspace(0,1,L);
% [hist_O,bins]=hist(XYn(:),bins);
% hist_O=hist_O/sum(hist_O);
% hist_O(hist_O==0)=nan;
% bar(bins,hist_O,'k')
% ylim([0 1]);
% switch kk
%     case 1
%         xlim([-0.05 1.05]);
%     otherwise
%         xlim([-0.01 0.25]);
% end
% 
% ylabel('frequency','fontweight','bold')
% xlabel('pair obs. frequency','fontweight','bold')
end

% %
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision';
Export2Folder(['Observations.eps'],target_folder) 