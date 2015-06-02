%%%%%%%%%%%%% KNOWN MU %%%%%%%%%%%%%%%%%%%%%%%%%%

% ts=[.01 .05 .1 .5 1]*T;
ts=[.1,.5,1]*T;
nts=length(ts);

ls=[0;.0001;.001;.01;.05;.1];
nls=length(ls);
tol=.01;

obs=[.05;.1;.2;1];
nobs=length(obs);


lasso_ests=zeros(N,N,nts,nobs);
lasso_mse=zeros(nts,nobs,1);
lasso_corr=zeros(nts,nobs,1);

lasso_ell_ests=zeros(N,N,nts,nobs);
lasso_ell_mse=zeros(nts,nobs,1);
lasso_ell_corr=zeros(nts,nobs,1);

omp_ests=zeros(N,N,nts,nobs,nls);
omp_mses=zeros(nts,nobs,nls);
omp_corrs=zeros(nts,nobs,nls);

omp_ell_ests=zeros(N,N,nts,nobs,nls);
omp_ell_mses=zeros(nts,nobs,nls);
omp_ell_corrs=zeros(nts,nobs,nls);

glasso_ests=zeros(N,N,nts,nobs);
glasso_mses=zeros(nts,nobs);
glasso_corrs=zeros(nts,nobs);

%%
for oo=1:nobs
disp(['                           obs ' num2str(oo) ' of ' num2str(nobs)])

for tt=1:nts
    disp(['                tt ' num2str(tt) ' of ' num2str(nts)])
    
    sample_ratio=obs(oo);
    observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
    sampled_spikes=observations.*spikes;
    
    [Cxx, Cxy,EW,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);
    
    V=-diag(pi*(rates.*log(rates)+(1-rates).*log(1-rates))/8); %for the gradient
    V(isnan(V))=0; %take care of 0*log(0) cases...
    
    glasso_ests(:,:,nts,nobs)=EW;
    glasso_mses(tt,oo)=norm(EW-W,'fro')/norm(W,'fro');
    glasso_corrs(tt,oo)=corr(EW(:),W(:));
    
    %%%%%% omp %%%%%%%%%%%
     for l=1:nls
        disp(['l ' num2str(l) ' of ' num2str(nls)])
        
        EW=EstimateA_OMP(Cxx,Cxy,spar,tol,ls(l),MeanMatrix,rates);
        omp_ests(:,:,tt,oo,l)=EW;
        omp_mses(tt,oo,l)=norm(EW-W,'fro')/norm(W,'fro');
        omp_corrs(tt,oo,l)=corr(EW(:),W(:));
        
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
        EW2=diag(amp)*EW;
        omp_ell_ests(:,:,tt,oo,l)=EW2;
        omp_ell_mses(tt,oo,l)=norm(EW2-W,'fro')/norm(W,'fro');
        omp_ell_corrs(tt,oo,l)=corr(EW2(:),W(:));
     end
     
     
    
    %%%%%% lasso %%%%%%%%%%%
    EW=EstimateA_L1_logistic(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);
    
    lasso_ests(:,:,tt,oo)=EW;
    lasso_mse(tt,oo)=norm(EW-W,'fro')/norm(W,'fro');
    lasso_corr(tt,oo)=corr(EW(:),W(:));
    
    if strcmp(neuron_type,'logistic')
    [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
    else
        amp=1;
        Ebias2=Ebias;
    end

    EW2=diag(amp)*EW;
    lasso_ell_ests(:,:,tt,oo)=EW2;
    lasso_ell_mse(tt,oo)=norm(EW2-W,'fro')/norm(W,'fro');
    lasso_ell_corr(tt,oo)=corr(EW2(:),W(:));

end

end


%% Looping unknown mu %%%%%%%%%%%%%%%%%%%%%%%%%

tt=2;
ll=3;
oo=4;

maxRank=3;
nLambda=50;
lambda=0;

niter=10;
allEWs=cell(niter,1);
allMs=cell(niter,1);
loop_mses=zeros(niter,1);
loop_corrs=zeros(niter,1);
M_mses=zeros(niter,1);
M_corrs=zeros(niter,1);

sample_ratio=obs(oo);
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;
    
[Cxx, Cxy,EW,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),gg,restricted_penalty,pos_def,est_spar,W);    

for iter=1:niter
    disp(iter)
    if iter==1
        EW=squeeze(omp_ests(:,:,tt,oo,1));
    else
        EW=EstimateA_OMP(Cxx,Cxy,spar,tol,ls(ll),M,rates);
    end
    allEWs{iter}=EW;
    loop_mses(iter)=norm(EW-W,'fro')/norm(W,'fro');
    loop_corrs(iter)=corr(W(:),EW(:));

    [M,ranks]=softImputePlus(EW,maxRank,nLambda,lambda);
    allMs{iter}=M;
    
    M_mses(iter)=norm(M-MeanMatrix,'fro')/norm(MeanMatrix,'fro');
    M_corrs(iter)=corr(M(:),MeanMatrix(:));
    
end

%% distance-dependent connectivity






