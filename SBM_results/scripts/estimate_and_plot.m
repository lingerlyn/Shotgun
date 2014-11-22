


%% First determine distance-dependent penalty

ts=T; %looping over T optional - not done here.
obs=[.05;.1;.2;1];
%dist-dep penalties
% sls=[0;1e-4;2e-4;5e-4;7e-4;1e-3;.005;.01;.05;.1];
sls=[0;1e-8;1e-5;1e-4;2e-4;5e-4;7e-4;1e-3;.005;.01;.05;.1];

nts=length(ts);
nobs=length(obs);
nsls=length(sls);

EW_d=cell(nts,nobs,nsls,1);
EW_d_ell=cell(nts,nobs,nsls,1);
EW_d_corr=zeros(nts,nobs,nsls,1);
EW_d_ell_corr=zeros(nts,nobs,nsls,1);


for oo=1
disp(['                           obs ' num2str(oo) ' of ' num2str(nobs)])

    for tt=1:nts
        disp(['                tt ' num2str(tt) ' of ' num2str(nts)])

        sample_ratio=obs(oo);
        observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
        sampled_spikes=observations.*spikes;

        [Cxx, Cxy,EW,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);

        for l=1:nsls
            disp(l)
            EW_d{tt,oo,l}=EstimateA_Greedy_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,0,zeros(N),sbm.fd,sls(l),idenTol);
%             EW_d{tt,oo,l}=EstimateA_OMP2_Exact_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,0,zeros(N),sbm.fd,sls(l),idenTol);
            [amp, Ebias2]=logistic_ELL(rates,EW_d{tt,oo,l},Cxx,Cxy);
            EW_d_ell{tt,oo,l}=diag(amp)*EW_d{tt,oo,l};
            EW_d_corr(tt,oo,l)=corr(EW_d{tt,oo,l}(:),W(:));
            EW_d_ell_corr(tt,oo,l)=corr(EW_d_ell{tt,oo,l}(:),W(:));
        end
    
    end
end
%
% figure; plot(squeeze(EW_d_ell_corr(1,1,:)),'.')
[~,sll]=max(squeeze(EW_d_ell_corr(1,1,:)));

% Now find correct mean penalty parameter
%
ts=T;
nts=length(ts);

% ls=[0;.0001;.001;.005;.01;.05;.1];
ls=[0;.0001;.001;.005;.01;.05;.1;.15;.2;.25;.3;.5;1];
nls=length(ls);
obs=[.05;.1;.2;1];
nobs=length(obs);

greedy_dale_ests=zeros(N,N,nts,nobs,nls);
greedy_dale_mses=zeros(nts,nobs,nls);
greedy_dale_corrs=zeros(nts,nobs,nls);

greedy_dale_ell_ests=zeros(N,N,nts,nobs,nls);
greedy_dale_ell_mses=zeros(nts,nobs,nls);
greedy_dale_ell_corrs=zeros(nts,nobs,nls);

GreedyIdents=zeros(N,nts,nobs,nls);

%
for oo=1
disp(['                           obs ' num2str(oo) ' of ' num2str(nobs)])

    for tt=1:nts
        disp(['                tt ' num2str(tt) ' of ' num2str(nts)])

        sample_ratio=obs(oo);
        observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
        sampled_spikes=observations.*spikes;

        [Cxx, Cxy,EW,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);


         for l=1:nls
             tic;
             disp(['l ' num2str(l) ' of ' num2str(nls)])

            %Dale

            [EW,GreedyIndents(:,tt,oo,l)]=EstimateA_Greedy_Dale_Iter(Cxx,Cxy,rates,est_spar,ls(l),MeanMatrix,idenTol);
%             [EW,GreedyIndents(:,tt,oo,l)]=EstimateA_OMP2_Exact_Dale_Iter(Cxx,Cxy,rates,est_spar,ls(l),MeanMatrix,idenTol);


            greedy_dale_ests(:,:,tt,oo,l)=EW;
            greedy_dale_mses(tt,oo,l)=norm(EW-W,'fro')/norm(W,'fro');
            greedy_dale_corrs(tt,oo,l)=corr(EW(:),W(:));

            [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
            EW_ELL=diag(amp)*EW;

            greedy_dale_ell_ests(:,:,tt,oo,l)=EW_ELL;
            greedy_dale_ell_mses(tt,oo,l)=norm(EW_ELL-W,'fro')/norm(W,'fro');
            greedy_dale_ell_corrs(tt,oo,l)=corr(EW_ELL(:),W(:));
            
            toc
         end

    end
end
% figure; hold on;
% plot(squeeze(greedy_dale_ell_corrs(1,:,:))')

[~,ll]=max(squeeze(greedy_dale_ell_corrs(1,1,:)));


% now loop to infer unknown mean and distdep

%which observation level / amount of time to use
oo=1;
tt=1;

% parameters picked from 'cross-validation'
% ll=4; % mean penalty
% sll=9; % distance dependence

% ll=4; %from the 'high' var stuff
% sll=8;
%%
opt_lambda_M=ls(ll)/2;
opt_lambda_D=sls(sll)/2;
%

maxRank=2;
nLambda=50;
lambda=0;
niter=10;

allEWs=cell(niter,1);
allfhats=cell(niter,1);
loop_mses=zeros(niter,1);
loop_corrs=zeros(niter,1);
fhat_mses=zeros(niter,1);

%using only distance-dependence
allEWs_distonly=cell(niter,1);
loop_distonly_mses=zeros(niter,1);
loop_distonly_corrs=zeros(niter,1);
allfhats_distonly=cell(niter,1);
fhat_mses_distonly=zeros(niter,1);

%using only low-rank means
allEWs_meanonly=cell(niter,1);
loop_meanonly_mses=zeros(niter,1);
loop_meanonly_corrs=zeros(niter,1);
allMs_m=cell(niter,1);
M_meanonly_mses=zeros(niter,1);
M_meanonly_corrs=zeros(niter,1);


allMs=cell(niter,1);
M_mses=zeros(niter,1);
M_corrs=zeros(niter,1);


sample_ratio=obs(oo);
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;

[Cxx, Cxy,EW,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);

%
for i=1:niter
    tic;
    disp(i)
    if i==1
        EW=EstimateA_Greedy_Dale_Iter(Cxx,Cxy,rates,est_spar,0,zeros(N),idenTol);
%         EW=EstimateA_OMP2_Exact_Dale_Iter(Cxx,Cxy,rates,est_spar,0,zeros(N),idenTol);
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        EW=diag(amp)*EW;
%         
        EW_d=EW; %all the same on the first iteration.
        EW_m=EW;
        
        
    else
        EW=EstimateA_Greedy_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,M,fhat,opt_lambda_D,idenTol);
%         EW=EstimateA_OMP2_Exact_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,M,fhat,opt_lambda_D,idenTol);
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        EW=diag(amp)*EW;
        
        EW_d=EstimateA_Greedy_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,0,M,fhat_d,opt_lambda_D,idenTol);
%         EW_d=EstimateA_OMP2_Exact_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,0,M,fhat_d,opt_lambda_D,idenTol);
        [amp, Ebias2]=logistic_ELL(rates,EW_d,Cxx,Cxy);
        EW_d=diag(amp)*EW_d;
        
        EW_m=EstimateA_Greedy_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,M_m,idenTol);
%         EW_m=EstimateA_OMP2_Exact_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,M_m,idenTol);
        [amp, Ebias2]=logistic_ELL(rates,EW_m,Cxx,Cxy);
        EW_m=diag(amp)*EW_m;
    end
    
    
    allEWs{i}=EW;
    loop_mses(i)=norm(EW-W,'fro')/norm(W,'fro');
    loop_corrs(i)=corr(W(:),EW(:));
    
    allEWs_distonly{i}=EW_d;
    loop_distonly_mses(i)=norm(EW_d-W,'fro')/norm(EW_d,'fro');
    loop_distonly_corrs(i)=corr(EW_d(:),W(:));
    
    allEWs_meanonly{i}=EW_m;
    loop_meanonly_mses(i)=norm(EW_m-W,'fro')/norm(EW_m,'fro');
    loop_meanonly_corrs(i)=corr(EW_m(:),W(:));
    
    
    %dist-dep
    x=sbm.neuron_distances(:); y=~~EW(:);
    options = optimset('GradObj','on','Display','off','LargeScale','off');
    [new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
    fhat=reshape(sigm(new_ab(1)*sbm.neuron_distances(:)+new_ab(2)),[N,N]);
    allfhats{i}=fhat;
    
    fhat_mses(i)=sum( (fhat(:)-sbm.fd(:)).^2)/sum(sbm.fd(:).^2);
    
    %soft-impute
    MM=EW; MM(~~eye(N))=0;
    [M,ranks]=softImputePlus(MM,maxRank,nLambda,lambda);
    idx=~~diag(EW); %select non-zeros - should be all of them...
    M(~~eye(N))=diag(mean(EW( false(N)|diag(idx))));
    allMs{i}=M;
    
    M_mses(i)=norm(M-MeanMatrix,'fro')/norm(MeanMatrix,'fro');
    M_corrs(i)=corr(M(:),MeanMatrix(:));
    
    
    %dist-dep for distonly
    x=sbm.neuron_distances(:); y=~~EW_d(:);
    options = optimset('GradObj','on','Display','off','LargeScale','off');
    [new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
    fhat_d=reshape(sigm(new_ab(1)*sbm.neuron_distances(:)+new_ab(2)),[N,N]);
    allfhats_distonly{i}=fhat_d;
    fhat_mses_distonly(i)=sum( (fhat_d(:)-sbm.fd(:)).^2)/sum(sbm.fd(:).^2);
    
    %SI for means only
    MM_m=EW_m; MM_m(~~eye(N))=0;
    [M_m,ranks]=softImputePlus(MM_m,maxRank,nLambda,lambda);
    idx=~~diag(EW_m); %select non-zeros - should be all of them...
    M_m(~~eye(N))=diag(mean(EW_m( false(N)|diag(idx))));
    allMs_m{i}=M_m;
    
    M_meanonly_mses(i)=norm(M_m-MeanMatrix,'fro')/norm(MeanMatrix,'fro');
    M_meanonly_corrs(i)=corr(M_m(:),MeanMatrix(:));

    toc
end
%
%estimate with both known
greedy_known=EstimateA_Greedy_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,MeanMatrix,sbm.fd,opt_lambda_D,idenTol);
% greedy_known=EstimateA_OMP2_Exact_DistDep_Dale_Iter(Cxx,Cxy,rates,est_spar,opt_lambda_M,MeanMatrix,sbm.fd,opt_lambda_D,idenTol);
[amp, ~]=logistic_ELL(rates,greedy_known,Cxx,Cxy);
greedy_known_ell=diag(amp)*greedy_known;

greedy_known_corr=corr(greedy_known(:),W(:));
greedy_known_ell_corr=corr(greedy_known_ell(:),W(:));


figure; hold on; 
plot(0:niter-1,loop_corrs); title('corr'); xlabel('iteration')
plot(0:niter-1,loop_meanonly_corrs,'r');
plot(0:niter-1,loop_distonly_corrs,'g');
plot(0,greedy_known_ell_corr,'*r');

xlabel('iteration'); ylabel('correlation')
legend 'full model' 'inferring mean' 'inferring dist dep' 'known mean and dist dep'

%% Estimate other stuff

%lasso estimation
pen_diag=0; warm=0; N_stim=0;
lassoEW=EstimateA_L1_logistic_Accurate(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);
[amp, Ebias2]=logistic_ELL(rates,lassoEW,Cxx,Cxy);
lassoEW_ell=diag(amp)*lassoEW;

%lasso with dale
lasso_dale=EstimateA_L1_logistic_Accurate_Dale_Iter(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm,idenTol);
[amp, Ebias2]=logistic_ELL(rates,lasso_dale,Cxx,Cxy);
lasso_dale_ell=diag(amp)*lasso_dale;

%greedy without means
greedy_EW=EstimateA_Greedy(Cxx,Cxy,rates,est_spar,0,zeros(N));
[amp, Ebias2]=logistic_ELL(rates,greedy_EW,Cxx,Cxy);
greedy_EW_ell=diag(amp)*greedy_EW;

%greedy without means with dale
greedy_EW_dale=EstimateA_Greedy_Dale_Iter(Cxx,Cxy,rates,est_spar,0,M,idenTol);
[amp, Ebias2]=logistic_ELL(rates,greedy_EW_dale,Cxx,Cxy);
greedy_EW_dale_ell=diag(amp)*greedy_EW_dale;


%greedy with inferring means dale
greedy_infer_m=allEWs_meanonly{end};

%greedy with inferring distdep dale
greedy_infer_d=allEWs_distonly{end};

%greedy with inferring means dale distdep
greedy_infer_both=allEWs{end};

% make bar figure
%lassoEW greedy_EW greedy_EW_dale lasso_dale greedy_EW_infer_m greedy_EW_infer_d opm2_known
% [R,C,Z,S] = GetWeightsErrors(true_A,EA);

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
x_ticks={'lasso','lasso+dale','OMP','OMP+dale','OMP+infer means','OMP+infer dd','OMP+infer both','OMP+known means & dd'};
fontsize=12;
fontsize2=1.5*fontsize;


% bar( [R,correlation, zero_matching,sign_matching] );  
figure; subplot(4,1,1);
bar([lassoR,lasso_dale_R,greedy_R,greedy_dale_R,greedy_infer_m_R,greedy_infer_d_R,greedy_infer_both_R,greedy_known_R]);
ylim([0 1]); title('R')
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

subplot(4,1,2);
bar([lassoC,lasso_dale_C,greedy_C,greedy_dale_C,greedy_infer_m_C,greedy_infer_d_C,greedy_infer_both_C,greedy_known_C]);
ylim([0 1]); title('C')
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

subplot(4,1,3);
bar([lassoZ,lasso_dale_Z,greedy_Z,greedy_dale_Z,greedy_infer_m_Z,greedy_infer_d_Z,greedy_infer_both_Z,greedy_known_Z]);
ylim([0 1]); title('Z')
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

subplot(4,1,4);
bar([lassoS,lasso_dale_S,greedy_S,greedy_dale_S,greedy_infer_m_S,greedy_infer_d_S,greedy_infer_both_S,greedy_known_S]);
ylim([0 1]); title('S')
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

%% color + scatter plots

set(0,'DefaultTextInterpreter', 'latex');

cc=[-2,2]; %have to pick after looking after them once

%names: 'lasso+dale','OMP','OMP+dale','OMP+infer means','OMP+infer dd','OMP+infer both','OMP+known means & dd'
figure; 
subplot(4,2,1); imagesc(lassoEW_ell); colorbar; title('lasso'); caxis(cc);
subplot(4,2,2); plot(W(:),lassoEW_ell(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('lasso');

subplot(4,2,3); imagesc(lasso_dale_ell); colorbar; title('lasso+Dale'); caxis(cc);
subplot(4,2,4); plot(W(:),lasso_dale_ell(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('lasso+Dale');

subplot(4,2,5); imagesc(greedy_EW_ell); colorbar; title('OMP'); caxis(cc);
subplot(4,2,6); plot(W(:),greedy_EW_ell(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP');

subplot(4,2,7); imagesc(greedy_EW_dale_ell); colorbar; title('OMP+Dale'); caxis(cc);
subplot(4,2,8); plot(W(:),greedy_EW_dale_ell(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP+Dale');


figure;
subplot(4,2,1); imagesc(greedy_infer_m); colorbar; title('OMP+Dale+Infer M'); caxis(cc);
subplot(4,2,2); plot(W(:),greedy_infer_m(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP+Dale+Infer M');

subplot(4,2,3); imagesc(greedy_infer_d); colorbar; title('OMP+Dale+Infer M'); caxis(cc);
subplot(4,2,4); plot(W(:),greedy_infer_d(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP+Dale+Infer M');

subplot(4,2,5); imagesc(greedy_infer_both); colorbar; title('OMP+Dale+Infer M and D'); caxis(cc);
subplot(4,2,6); plot(W(:),greedy_infer_both(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP+Dale+Infer M and D');

subplot(4,2,7); imagesc(greedy_known_ell); colorbar; title('OMP+Dale+Known M and D'); caxis(cc);
subplot(4,2,8); plot(W(:),greedy_known_ell(:),'.'); xlabel('W'); ylabel('$\hat{W}$'); title('OMP+Dale+Known M and D');

%% Save estimates

save('results.mat','params','W','lassoEW_ell','lasso_dale_ell','greedy_EW_ell','greedy_EW_dale_ell','greedy_infer_m','greedy_infer_d',...
    'greedy_infer_both','greedy_known_ell','allEWs','allEWs_distonly','allEWs_meanonly','allMs','allMs_m','opt_lambda_M','opt_lambda_D');