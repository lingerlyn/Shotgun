close all
clear all
clc

addpath('Misc')
addpath('Plotting')

set(0,'defaultlinelinewidth',2.5)
set(0,'DefaultAxesFontSize',24)

N=33;

B=ones(N/3);
M=[ 0.3*B -0.4*B B ; B -B 0.8*B ; 0.7*B -0.6*B 0.6*B ]; %mean matrix
M(eye(N)>0.5)=-1; %diagonal

E=sqrt(0.02)*randn(N);
spar=0.2;
p=GetProb(N,spar,1:N);
S=rand(N)<p;

mi=min(M(:)+E(:));
ma=max(M(:)+E(:));


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
imagesc(M,[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['M.eps'],target_folder) 
% subplot(a,b,4)
figure(3)
imagesc(M+E,[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['ME.eps'],target_folder) 
% subplot(a,b,6)
figure(4)
imagesc(S.*(M+E),[mi ma])
colorbar
set(gca,'xtick',[],'ytick',[])
Export2Folder(['SME.eps'],target_folder) 
% subplot(a,b,8)
figure(5)
p=GetProb(N,spar,1);
d=1:floor(N/2);
plot(d,p(d));
xlim([0 d(end)])
xlabel('distance')
ylabel('conn. probability')
Export2Folder(['P(d).eps'],target_folder) 