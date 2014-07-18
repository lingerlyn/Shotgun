function [EA,alpha, rates_A, s_sq]=EstimateA(CXX,CXY,mY,mYn,priors)
% Calculate VB posterior distribution of spike-and-slab distributed weights from Carbonetto and Stephens. 
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% priors.eta - prior slab means (for each weight, so an NXN matrix)
% priors.a - prior inclusion probability (probability of being non-zero)
% priors.ss2 - slab variance (NXN)
% priors.noise_var - noise variance of y_k=W(k,:)*x+noise (NX1)
% params:
% n_eff - number of effective observations
% m - mean firing rate for each neuron
% options - optimization options for fminunc
% OUTPUTS:
% EA: ML estimate of weights, with a gain and bias correction found by
% logistic ELL maximization
% alpha: NXN matrix of posterior inclusion probabilities
% rates_A: posterior slab means
% s_sq: posterior slab variances

    addpath('Misc')
    addpath('EstimateConnectivity\Carbonetto')

    params.n_eff=mYn;
    params.m=mY./(mYn+eps);
    params.options = optimset('GradObj','on','Display','off');

    N=length(CXX);
    
    alpha=zeros(N);
    rates_A=zeros(N);
    s_sq=zeros(N);
    
    eta=priors.eta;
    eta0=eta;
    a=priors.a;
    a0=a;
    ss2=priors.ss2;
    noise_var=priors.noise_var;
    
    EA=zeros(N);
    m=params.m;
    z=randn(1e4,1);

	for k=1:N	
            
		[alpha(k,:), rates_A(k,:), s_sq(k,:)]=...
			varbvs_general_ss(CXX*params.n_eff(k),CXY(:,k)*params.n_eff(k),noise_var(k),ss2(k,:)',logit(a(k,:)'),eta(k,:)',a0(k,:)',eta0(k,:)',1e-6);

        %Use ELL to fit bias and gain terms:
        ea=(alpha(k,:)>0.5).*rates_A(k,:);
        En=m(k);
        Enx=ea*(CXY(:,k)+m*m(k));
        x=sqrt(ea*CXX*ea')*z+ea*m;
        [new_ab,~,~]=fminunc(@(ab) twod_logistic_ELL(ab,Enx,En,x),[0;0],params.options);
        EA(k,:)=ea*new_ab(1);
        
        
	end




end

