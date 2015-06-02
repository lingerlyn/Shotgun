function [EW, Eb,quality,error_rates,lambda_path]=EstimateA_L1_cavity_forLogReg(XX,YX,mX,mY,sparsity,warm,true_W,use_sampling )
% Following derivation in GLM_GradUsingSampling 13.1.2015.lyx

% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% XX - inputs second moment matrix
% YX - output-input cross moment
% mX - mean inputs
% mY - mean outputs
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% true_W - for calculating Error
% use_sampling - use sampling to calculate Gaussian integral. Typically slower, so don't use

% OUTPUTS:
% EW: MAP estimate of weights
% Eb: MAP estimate of weights
% quality: quality estimates during convergence
% error_rates: 4xL array erros rates, for positive and negative weights estimated during during exponential search (on regularization path)
% lambda_path: regularization constants used during binary search

%internal params
Tol_sparse=0.02; %tolerance for sparsity level
N=length(mX);
Tol_FISTA=1e-6; %toleratnce for fista
max_iterations=2000;
last_max_iterations=max_iterations;
lambda0=0.002; %initial value  of regularization constant
rho=2; % parameter for exponential search. should be larger then 1
show_progress=0; %show figures with progress
show_w=0; %#ok show figures with W (estimate vs true) 
reweighted_L1_Reps=1;
reweighted_L1_eps=1e-3;

CXX=XX-mX*mX';
CXY=YX'-mX*mY';

V=-diag(pi*(mY.*log(mY)+(1-mY).*log(1-mY))/8); %for the gradient
V(isnan(V))=0; %take care of 0*log(0) cases...
L=500*max(V(:))*max(eig(CXX)); %lipshitz constant - an educated guess. make sure this is large enough to ensure convergence
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


for kk=1:reweighted_L1_Reps
    if kk>1
        mask=1./(abs(x)+reweighted_L1_eps);
    else
        mask=1;
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

    FISTA_cond=1;
    iteration=0;
    MAE=[]; %relative mean absolute error of the current estimate, in comparsion to the previous estimate

    while FISTA_cond
      
        t=t_next;
        x_prev=x;
        Eb_prev=Eb;       
        mean_U=y*mX+z;
        var_U=diag(y*CXX*y');
        CXX_diag=CXX(eye(N)>0.5);
        mean_U_cavity=bsxfun(@plus,mean_U,y*(CXX*diag((1-mX)./CXX_diag)));
        var_U_cavity=bsxfun(@plus,var_U,-((y*CXX).^2)*diag(1./CXX_diag));
        var_U_cavity(var_U_cavity<0)=0; %correct for numerical errors;

        if use_sampling %slower, so I don't use it
            L_z=1e3;
            samples=randn(1,1,L_z);
            u=bsxfun(@plus,mean_U,bsxfun(@times,sqrt(var_U),samples));        
            u_cavity=bsxfun(@plus,mean_U_cavity,bsxfun(@times,sqrt(var_U_cavity),samples));    
            Ef=squeeze(mean(1./(1+exp(-u)),3));
            Ef_cavity=squeeze(mean(1./(1+exp(-u_cavity)),3));
        else
            Ef=1./(1+exp(-(mean_U./sqrt(1+pi*var_U/8)))); % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
            Ef_cavity=1./(1+exp(-(mean_U_cavity./sqrt(1+pi*var_U_cavity/8)))); % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
        end
            u=y-(2/L)*(Ef_cavity*diag(mX)-YX);
            Eb=z-(2/L)*(Ef-mY);

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
            elseif iteration>last_max_iterations/2
                EW=EW+x/(last_max_iterations/2);
            else
                EW=0;
            end
                   
        end

        if nnz(x)==0 %if everything is zero except diagonal, quite
            break        
        end

    end
    %%% 
        EW=x;
        [~,~,~,~,TPR_p,FPR_p,TPR_n,FPR_n]= GetWeightsErrors( true_W,x );
        error_rates(end+1,:)=[TPR_p,FPR_p,TPR_n,FPR_n]; %#ok  
        
    
        sparsity_measure=mean(~~x);
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

