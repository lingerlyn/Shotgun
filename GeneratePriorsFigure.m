close all
clear all
clc

addpath('Misc')
addpath('Plotting')

set(0,'defaultlinelinewidth',2.5)
set(0,'DefaultAxesFontSize',24)
load(fullfile('SBM_results','Full_model','results.mat'));

N=params.connectivity.N;

V=params.sbm.MeanMatrix;
% M(eye(N)>0.5)=-1; %diagonal

E=sqrt(0.05)*randn(N);
% spar=0.2;
% p=GetProb(N,spar,1:N);
% S=rand(N)<p;
S=sortedW~=0;

mi=min(V(:)+E(:));
ma=max(V(:)+E(:));

% a=3; b=3;
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\PriorFigures\';


% subplot(a,b,5)
figure(1)
imagesc(S)
colormap(gray)
cbr=colorbar
set(gca,'xtick',[],'ytick',[])
set(cbr,'YTick',[0 1])
Export2Folder(['S.eps'],target_folder) 

% subplot(a,b,1)
figure(2)
colormap('jet')
imagesc(V,[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['V.eps'],target_folder) 
% subplot(a,b,4)
figure(3)
imagesc(V+E,[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['VE.eps'],target_folder) 
% subplot(a,b,6)
figure(4)
imagesc(sortedW,[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['W.eps'],target_folder) 
% subplot(a,b,8)
figure(5)
% p=GetProb(N,spar,1);
d=(0:0.01:N/2)/N;
p=1./(1+exp(15*d-1.3));
plot(d,p);
xlim([0 d(end)])
xlabel('distance')
ylabel('conn. probability')
Export2Folder(['P(d).eps'],target_folder) 