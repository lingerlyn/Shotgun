function [EW, Eb,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity(CXX,CXY,rates,sparsity,N_stim,pen_diag,pen_dist,warm,true_W,centers)
% Following derivation in GLM_GradUsingSampling 13.1.2015.lyx

% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% pen_dist - 1 if we should we adjust penalty according to distance; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% true_W - for calculating Error
% centers - NxD with locations of neurons in D-dimensional space
% OUTPUTS:
% EW: MAP estimate of weights
% Eb: MAP estimate of weights
% quality: quality estimates during convergence
% error_rates: 4xL array erros rates, for positive and negative weights estimated during during exponential search (on regularization path)
% lambda_path: regularization constants used during binary search


%internal params
Tol_sparse=0.02; %tolerance for sparsity level
N=length(rates);
Tol_FISTA=1e-6; %toleratnce for fista
max_iterations=1000;
last_max_iterations=max_iterations;
lambda0=0.002; %initial value  of regularization constant
rho=2; % parameter for exponential search. should be larger then 1
use_sampling=0; %use sampling to calculate Gaussian integral. Typically slower, so don't use
show_progress=0; %show figures with progress
show_w=0; %#ok show figures with W (estimate vs true) 
reweighted_L1_Reps=1;
reweighted_L1_eps=1e-3;

CXX=full(CXX);
CXY=full(CXY);
rates=full(rates);

V=-diag(pi*(rates.*log(rates)+(1-rates).*log(1-rates))/8); %for the gradient
V(isnan(V))=0; %take care of 0*log(0) cases...
L=200*max(V(:))*max(eig(CXX)); %lipshitz constant - an educated guess. make sure this is large enough to ensure convergence
% L=5*max(eig(CXX)); %lipshitz constant

%initialize FISTA
% x=zeros(N);%*sparse(N,N);
x=sparse(CXY'/CXX);
y=x;
Eb=0;
z=0;
quality=[]; %MSE from true W        
lambda_path=[];
error_rates=[];

if pen_diag
    mask0=ones(N);
else  
    mask0=ones(N)-eye(N); % used to remove L1 penalty from diagonal
end

if pen_dist
    p=GetDistProb(centers);
    mask0=mask0./(p/max(p(:)));
end

for kk=1:reweighted_L1_Reps
    if kk>1
        mask=mask0./(abs(x)+reweighted_L1_eps);
    else
        mask=mask0;
    end
    
    %initialize exponential binary search
    lambda_low=-1;
    lambda_high=-1;

    if sparsity==1
        relative_sparsity=0;
        lambda=0;
    else
        relative_sparsity=inf;  
        lambda=lambda0;
    end
    spar_array=[];
    loop_cond=1;

    while  loop_cond %binary search for lambda that give correct sparsity level

    loop_cond=( relative_sparsity>  Tol_sparse);

    %%% FISTA
    if ~warm
%         x=sparse(N,N);
        x=sparse(CXY'/CXX);
        y=x;
        Eb=0;
        z=0;
    end

    t_next=1;

    MYX=CXY'+rates*rates';

    FISTA_cond=1;
    iteration=0;
    MAE=[]; %relative mean absolute error of the current estimate, in comparsion to the previous estimate

    while FISTA_cond
      
        t=t_next;
        x_prev=x;
        Eb_prev=Eb;       
        mean_U=y*rates+z;
        var_U=diag(y*CXX*y');
%         var_U(var_U<0)=0; %correct for numerical errors;
        CXX_diag=CXX(eye(N)>0.5);
        mean_U_cavity=bsxfun(@plus,mean_U,y*(CXX*diag((1-rates)./CXX_diag)));
        var_U_cavity=bsxfun(@plus,var_U,-((y*CXX).^2)*diag(1./CXX_diag));
        var_U_cavity(var_U_cavity<0)=0; %correct for numerical errors;

        if use_sampling %slower, so I don't use it
            L_z=1e4;
            samples=randn(1,1,L_z);
            u=bsxfun(@plus,mean_U,bsxfun(@times,sqrt(var_U),samples));        
            u_cavity=bsxfun(@plus,mean_U_cavity,bsxfun(@times,sqrt(var_U_cavity),samples));    
            Ef=squeeze(mean(1./(1+exp(-u)),3));
            Ef_cavity=squeeze(mean(1./(1+exp(-u_cavity)),3));
        else
            Ef=1./(1+exp(-(mean_U./sqrt(1+pi*var_U/8)))); % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
    %         ind=Ef<approx_thresh;
    %         Ef(ind)=exp(mean_U(ind)+0.5*var_U(ind));

            Ef_cavity=1./(1+exp(-(mean_U_cavity./sqrt(1+pi*var_U_cavity/8)))); % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
    %         ind=Ef_cavity<approx_thresh;
    %         Ef_cavity(ind)=exp(mean_U_cavity(ind)+0.5*var_U_cavity(ind));
        end
    %         for nn=1:N        
    %             rates_cavity=rates+CXX(:,nn)*(1-rates(nn))/CXX(nn,nn);        
    %             CXX_cavity=CXX-CXX(:,nn)*CXX(nn,:)/CXX(nn,nn); %can remove cavity here since it does not seem to  hurt perfomance (?)
    %             A_cavity=1./sqrt(1+pi*y*CXX_cavity*y'/8);
    %             Ef_cavity(:,nn)=1./(1+exp(-((y*rates_cavity+z).*A_cavity(eye(N)>0.5))));    % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
    %         end
            u=y-(2/L)*( Ef_cavity*diag(rates)-MYX);
            Eb=z-(2/L)*(Ef-rates);

        x=ThresholdOperator(u,mask.*lambda/L);
        if any(~isfinite(u(:)))
            error('non finite x!')
        end
        t_next=(1+sqrt(1+4*t^2))/2;
        y=x+((t-1)/t_next)*(x-x_prev);    
        z=Eb+((t-1)/t_next)*(Eb-Eb_prev); 

%         Adaptive restart from Donoghue2012    
        do_restart=(sum(sum((u-y).*(x-x_prev))))+sum((Eb-z).*(Eb-Eb_prev))>0;
        if do_restart
            t_next=1;
            y=x;        
            z=Eb;
        end

        iteration=iteration+1;    
        MAE(end+1)=mean(abs(x(:)-x_prev(:))); %#ok    
        [R,correlation, zero_matching,sign_matching,~,~,~,~] = GetWeightsErrors( true_W,x );
        quality(end+1,:)=[R,correlation, zero_matching,sign_matching]; %#ok  
        lambda_path(end+1)=lambda; %#ok  
            
        if show_progress
            if ~show_w
                disp(['MAE=' num2str(MAE(end))]);
            end
            figure(1000)

            a=3; b=2;
             if show_w
                 mi=min([true_W(:); x(:)]);
                ma=max([true_W(:); x(:)]);
                subplot(a,b,1)
                imagesc(true_W,[mi ma])
                subplot(a,b,2)
                imagesc(x,[mi ma])
            end
            subplot(a,b,3)
            plot(MAE)
            ylabel('MAE');
            xlabel('iteration');
            subplot(a,b,4)
            plot(spar_array)
            hold all
            plot(sparsity+spar_array*0)
            hold off
            ylabel('Sparsity');
            xlabel('iteration');
            if show_w
                subplot(a,b,5)   
                foo=linspace(mi,ma,100);
                plot(foo,foo);
                hold all
                plot(true_W(:),x(:),'.');
                hold off
            end
            pause(1e-6)
            subplot(a,b,6)   
            bar([lambda_low,lambda, lambda_high])
        end

       FISTA_cond=(MAE(end)>Tol_FISTA)||(iteration<30);
        
        if loop_cond
            if (iteration>max_iterations);
                break    
            end
        else %last round        
            if (iteration>last_max_iterations); 
                warning('FISTA max iteration reached');
                break    
            end     
                   
        end

        if nnz(x(eye(N)<0.5))==0 %if everything is zero except diagonal, quite
            break        
        end

    end
    %%% 
        
        [~,~,~,~,TPR_p,FPR_p,TPR_n,FPR_n]= GetWeightsErrors( true_W,x );
        error_rates(end+1,:)=[TPR_p,FPR_p,TPR_n,FPR_n]; %#ok  
        
    

        temp=~~x(1:(end-N_stim),1:(end-N_stim));
        sparsity_measure=mean(temp(eye(size(temp,1))<0.5))    
        spar_array(end+1)=sparsity_measure; %#ok
        relative_sparsity=abs(sparsity_measure-sparsity)/sparsity;

        if sparsity_measure<sparsity
            lambda_high=lambda
        else
            lambda_low=lambda
        end
        
        %exponential binary search
        if lambda_high==-1
            lambda=rho*lambda;
        elseif lambda_low==-1
            lambda=lambda/rho;
        else
            lambda=(lambda_high+lambda_low)/2;
        end
  
    end
end
        
EW=x;
end

function y = ThresholdOperator( x , lambda )
%THRESHOLDOPERATOR Summary of this function goes here
% This function implements the threshold operator for the l1/l2 group norm
% as specified on "convex optimization with sparsity-inducing norms", page 12, last equation

%   Detailed explanation goes here
% x - input
% y - output
% lambda - regularization constant

    temp=(1-lambda./abs(x));
    temp(isnan(temp))=1;
    y=temp.*x;
    y(temp<0)=0;    
end

