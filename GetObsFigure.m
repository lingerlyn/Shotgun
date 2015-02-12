clear all
% close all
clc

% Generate all figures for the paper
addpath('GenerateSpikes');
addpath('Misc')

set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',18)
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.1], [0.05 0.05]);

colormap('gray')
title_height=1;
title_horz=-0.2;
title_font=18;
h=figure(1)
set(h,'units','normalized','outerposition',[0 0 0.6 1])

N=10;
T=1e5*N;
t_show=1:3e3;
sample_ratio=0.2;
N_stim=0;
sample_type_set={'continuous','fixed_subset','spatially_random','prob','double_continuous'};
obs_duration=100;

seed_sample=1;
a=2; b=3;
L=30; %number of histogram bins

subplot(a,b,1)
sample_type=sample_type_set{3};
observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
imagesc(observations(:,t_show))
% title('(A)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
xlabel('t');
ylabel('cell number');
title('Random scanning', 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
colorbar

subplot(a,b,2)
sample_type=sample_type_set{3};
% observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
XYn=observations(:,1:(end-1))*(observations(:,2:end))'/T;
ma=max(XYn(:));
imagesc(XYn,[0 ma]);
% title('(B)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
% title('results in all spike pairs being observed', 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
xlabel('cell number');
ylabel('cell number');
% ['random scanning' -   ; 'all spike pairs observed (good)']
colorbar

subplot(a,b,3)
[hist_O,bins]=hist(XYn(:),L);
stem(bins,hist_O,'k')
hist_max=100;
ylim([0 hist_max]);

ylabel('count')
xlabel('pair obs. frequency')

subplot(a,b,4)
sample_type=sample_type_set{5};
observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
imagesc(observations(:,t_show))
% title('(C)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
title('Continuous scanning', 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
xlabel('t');
ylabel('cell number');
colorbar


subplot(a,b,5)
sample_type=sample_type_set{5};
% observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
XYn=observations(:,1:(end-1))*(observations(:,2:end))'/T;
ma=max(XYn(:));
imagesc(XYn,[0 ma])
% title('(D)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
% title('results in some spike pairs being unobserved', 'Units', 'normalized', 'fontweight','bold','fontsize',title_font) 
xlabel('cell number');
ylabel('cell number');
colorbar

subplot(a,b,6)
[hist_O,bins]=hist(XYn(:),L);
stem(bins,hist_O,'k')
ylabel('count')
xlabel('pair obs. frequency')
ylim([0 hist_max]);

%%
% target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript';
% Export2Folder(['Observations.eps'],target_folder) 