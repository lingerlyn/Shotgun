clear;
load demo_data.mat

niter=10;
maxRank=3;
lambda=.1; %these params seem to do fine with softimpute.
nLambda=200;
EAs=cell(niter,1);
Fs=cell(niter,1);

N=length(CXX);
priors.eta=zeros(N);
priors.ss2=.05*ones(N);
priors.noise_var=.09*ones(N,1);
priors.a=.2*ones(N);

params.n_eff=mYn;
params.m=mY./(mYn+eps);
params.options = optimset('GradObj','on','Display','off');

for iter=1:niter
    
    if iter==1
        %use specified priors
        [EA,alpha, rates_A, s_sq]=EstimateA(CXX,CXY,priors,params);
        
    else
        %use priors from matrix completion
        priors.eta=F;
        [EA,alpha, rates_A, s_sq]=EstimateA(CXX,CXY,priors,params);
        
    end
    
    [F,ranks]=softImputePlus(EA,maxRank,nLambda,lambda);
    
    %save some more information
    EAs{iter}=EA;
    Fs{iter}=F;

end

    figure; imagesc(EA); colorbar; caxis([-1,1]);