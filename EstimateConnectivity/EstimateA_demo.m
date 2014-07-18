load EstimateA_demo_data.mat

N=length(CXX);
priors.eta=zeros(N);
priors.ss2=.05*ones(N);
priors.noise_var=.09*ones(N,1);
priors.a=.2*ones(N);

params.n_eff=mYn;
params.m=mY./(mYn+eps);
params.options = optimset('GradObj','on','Display','off');

[EA,alpha, rates_A, s_sq]=EstimateA(CXX,CXY,priors,params);

figure; imagesc(EA); colorbar; caxis([-1,1]);