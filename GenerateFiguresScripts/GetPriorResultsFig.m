clear all
% close all
clc

fontsize=12;
fontsize2=2*fontsize;

addpath('Misc')
addpath('Plotting')
load(fullfile('SBM_results','Full_model','results.mat'));
addpath('GenerateSpikes');
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',fontsize)

% subplot(3,3,1)
% imagesc(params.sbm.MeanMatrix)
% title('V','fontsize',fontsize2)
% colorbar
% subplot(3,3,2)
% imagesc(W)
% title('W','fontsize',fontsize2)
% colorbar
% subplot(3,3,3)
% % plot()
% title('$f(d)$','fontsize',fontsize2)
% 
% subplot(3,2,[3 4 5 6])
[lassoR,lassoC,lassoZ,lassoS] = GetWeightsErrors(W,lassoEW_ell);

[lasso_dale_R,lasso_dale_C,lasso_dale_Z,lasso_dale_S] = GetWeightsErrors(W,lasso_dale_ell);

[greedy_R,greedy_C,greedy_Z,greedy_S] = GetWeightsErrors(W,greedy_EW_ell);

[greedy_dale_R,greedy_dale_C,greedy_dale_Z,greedy_dale_S] = GetWeightsErrors(W,greedy_EW_dale_ell);

[greedy_infer_m_R,greedy_infer_m_C,greedy_infer_m_Z,greedy_infer_m_S] = GetWeightsErrors(W,greedy_infer_m);

[greedy_infer_d_R,greedy_infer_d_C,greedy_infer_d_Z,greedy_infer_d_S] = GetWeightsErrors(W,greedy_infer_d);

[greedy_infer_both_R,greedy_infer_both_C,greedy_infer_both_Z,greedy_infer_both_S] = GetWeightsErrors(W,greedy_infer_both);

[greedy_known_R,greedy_known_C,greedy_known_Z,greedy_known_S] = GetWeightsErrors(W,greedy_known_ell);

%create bar plot
% x_ticks={'R','C','Z','S'};
x_ticks={'L1','L1,D','L0','L0,D','L0,D,sbm','L0,D,dd','L0,D,sbm,dd','L0,D,Ksbm,Kdd'};
% x_ticks={'L1',{'L1' ,'D'},'L0',{'L0,Dale'},{'L0','sbm'},{'L0','dd'},{'L0','sbm','dd'},{'L0','known sbm','known dd'}};

bins=1*(1:length(x_ticks));


% bar( [R,correlation, zero_matching,sign_matching] );  
bar(bins,[lassoC,lasso_dale_C,greedy_C,greedy_dale_C,greedy_infer_m_C,greedy_infer_d_C,greedy_infer_both_C,greedy_known_C]);
ylim([0 1]); 
% my_xticklabels(gca,bins,x_ticks);
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
ylabel('C (correlation)','fontsize',fontsize2)

% target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\';
% Export2Folder(['PriorResults.eps'],target_folder) 