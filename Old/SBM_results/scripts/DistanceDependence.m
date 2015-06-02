%% Using L1 methods.

%set up the weird penalty thingie
betadist=@(a,b,x) x.^(1-a).*(1-x).^(b-1);
aaa=1;
bbb=2;
f=.01:.01:1;
% figure; plot(f,betadist(aaa,bbb,f));

lambdafunc=@(x)betadist(aaa,bbb,x);



ts=T; %for now
obs=1;
% tt=3; oo=4;
tt=1; oo=1;
%%
niter=10;

sample_ratio=obs(oo);
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;
[Cxx, Cxy,EW3,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);

V=diag(rates);

EW=EstimateA_L1_logistic(Cxx,Cxy,rates,est_spar,N_stim,pen_diag,warm);
if strcmp(neuron_type,'logistic')
    [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
else amp=1; 
end
EW=diag(amp)*EW;


allEWs=cell(niter+1,1);
allfhats=cell(niter+1,1);
distdep_mse=zeros(niter+1,1);
distdep_corrs=zeros(niter+1,1);
distdep_sel=zeros(niter+1,1);

allEWs{1}=EW;
x=D(:); y=~~EW(:);

options = optimset('GradObj','on','Display','off','LargeScale','off');
[new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
fhat=reshape(sigm(new_ab(1)*D(:)+new_ab(2)),[N,N]);
allfhats{1}=fhat;

distdep_mse(1)=norm(EW-W,'fro')/norm(W,'fro');
distdep_corrs(1)=corr(EW(:),W(:));
distdep_sel(1)=(sum(~~EW(:) & ~~W(:))+sum(~EW(:) & ~W(:)))/numel(W);


for i=1:niter
    disp(i)
    
    [EW,iterflag]=EstimateA_L1_mult(Cxx,Cxy,mean(fhat(:)),V,fhat,lambdafunc);
    if strcmp(neuron_type,'logistic')
        [amp, Ebias2]=logistic_ELL(rates,EW,Cxx,Cxy);
    else amp=1; 
    end
    EW=diag(amp)*EW;
       
    allEWs{i+1}=EW;
    
    x=D(:); y=~~EW(:);
    [new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
    fhat=reshape(sigm(new_ab(1)*D(:)+new_ab(2)),[N,N]);
    allfhats{i+1}=fhat;
    
    distdep_mse(i+1)=norm(EW-W,'fro')/norm(W,'fro');
    distdep_corrs(i+1)=corr(EW(:),W(:));
    distdep_sel(i+1)=(sum(~~EW(:) & ~~W(:))+sum(~EW(:) & ~W(:)))/numel(W);

    
end

figure; plot(distdep_mse); title('rmse'); xlabel('iteration')
figure; plot(distdep_corrs); title('corr'); xlabel('iteration')
figure; plot(distdep_sel); title('percent correct selection'); xlabel('iteration');

%plot real f vs estimate
figure; hold on;
plot(D(:),DD(:),'.');
plot(D(:),allfhats{end}(:),'.r');
legend 'true' 'estimated'
xlabel('distance'); ylabel('P(connection)')


%% Using OMP



ts=T; %for now
obs=1;
% tt=3; oo=4;
tt=1; oo=1;

niter=10;

sample_ratio=obs(oo);
observations=SampleSpikes(N,T,sample_ratio,sample_type,N_stim,seed_sample+1);
sampled_spikes=observations.*spikes;
[Cxx, Cxy,EW3,rates,obs_count] = GetStat(sampled_spikes(:,1:ts(tt)),observations(:,1:ts(tt)),glasso,restricted_penalty,pos_def,est_spar,W);
V=diag(rates);

omp_lambda=0;
tol=0.01;
EW=EstimateA_OMP(Cxx,Cxy,spar,tol,omp_lambda,zeros(N),rates);

allEWs=cell(niter+1,1);
allfhats=cell(niter+1,1);
distdep_mse=zeros(niter+1,1);
distdep_corrs=zeros(niter+1,1);
distdep_sel=zeros(niter+1,1);

allEWs{1}=EW;
x=D(:); y=~~EW(:);

options = optimset('GradObj','on','Display','off','LargeScale','off');
[new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
fhat=reshape(sigm(new_ab(1)*D(:)+new_ab(2)),[N,N]);
allfhats{1}=fhat;

distdep_mse(1)=norm(EW-W,'fro')/norm(W,'fro');
distdep_corrs(1)=corr(EW(:),W(:));
distdep_sel(1)=(sum(~~EW(:) & ~~W(:))+sum(~EW(:) & ~W(:)))/numel(W);

lambda2=1e-6; %NO idea...
% lambda2=5e-6;
% lambda2=10;
for i=1:niter
    disp(i)
    
    EW=EstimateA_OMP_DistDep(Cxx,Cxy,spar,tol,omp_lambda,zeros(N),rates,fhat,lambda2);

    
    x=D(:); y=~~EW(:);
    [new_ab,junk,exitflag]=fminunc(@(ab)logistic_opt(ab,x,y),[0;0],options);
    fhat=reshape(sigm(new_ab(1)*D(:)+new_ab(2)),[N,N]);
    allfhats{i+1}=fhat;
    
    distdep_mse(i+1)=norm(EW-W,'fro')/norm(W,'fro');
    distdep_corrs(i+1)=corr(EW(:),W(:));
    distdep_sel(i+1)=(sum(~~EW(:) & ~~W(:))+sum(~EW(:) & ~W(:)))/numel(W);

    
end

figure; plot(distdep_mse); title('rmse'); xlabel('iteration')
figure; plot(distdep_corrs); title('corr'); xlabel('iteration')
figure; plot(distdep_sel); title('percent correct selection'); xlabel('iteration');

%plot real f vs estimate
figure; hold on;
plot(D(:),DD(:),'.');
plot(D(:),allfhats{end}(:),'.r');
legend 'true' 'estimated'
xlabel('distance'); ylabel('P(connection)')


