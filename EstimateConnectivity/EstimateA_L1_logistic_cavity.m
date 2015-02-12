function [EW, Eb,MSE]=EstimateA_L1_logistic_cavity(CXX,CXY,rates,sparsity,N_stim,pen_diag,warm,is_spikes,true_W)
% Following derivation in GLM_GradUsingSampling 13.1.2015.lyx

% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% is_spikes - if equal 1 we use a cavity-style averages
% true_W - for calculating Error
% use_sampling - 1 if we should we use sampling to calculate Gaussian Integral
% OUTPUTS:
% EW: MAP estimate of weights
% Eb: MAP estimate of weights


%internal params
Tol_sparse=0.05; %tolerance for sparsity level
N=length(rates);
Tol_FISTA=1e-8; %toleratnce for fista
max_iterations=300;
last_max_iterations=5*max_iterations;
use_sampling=0; %use sampling to calculate Gaussian integral
show_progress=0; %show figures with progress
show_w=0; %#ok show figures with W (estimate vs true) 
reweighted_L1_Reps=1;
reweighted_L1_eps=1e-3;

L=5*max(eig(CXX)); %lipshitz constant - and educated guess. make sure this is large enough to ensure convergence
lambda_high0=1e2;%max(1e4,1e1*L); %maximum bound for lambda
lambda_low0=min(1e-8,1e-8*L);  %minimum bound for lambda
% approx_thresh=0.05; %approximation threshold

%initialize FISTA
% x=zeros(N);%*sparse(N,N);
x=sparse(CXY'/CXX);
y=x;
Eb=0;
z=0;


if pen_diag
    mask0=ones(N);
else  
    mask0=ones(N)-eye(N); % used to remove L1 penalty from diagonal
end

for kk=1:reweighted_L1_Reps
    if kk>1
        mask=mask0./(abs(x)+reweighted_L1_eps);
    else
        mask=mask0;
    end
    
    %initialize binary search
    lambda_high=lambda_high0;%max(1e4,1e1*L); %maximum bound for lambda
    lambda_low=lambda_low0;  %minimum bound for lambda
    % lambda_high=1e-1;
    % lambda_low=1e-5;
    loop_cond=1;  %flag for while llop
    cond=1;
    spar_array=[];

    while  loop_cond %binary search for lambda that give correct sparsity level
        if sparsity==1
            lambda=0;
        else
            lambda=(lambda_high+lambda_low)/2;
        end

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
    MSE=[];

     if is_spikes~=1
            error('non-spike case is not written yet')
    end    

    while FISTA_cond
        t=t_next;
        x_prev=x;
        Eb_prev=Eb;       
        mean_U=y*rates+z;
        var_U=diag(y*CXX*y');    
        CXX_diag=CXX(eye(N)>0.5);
        mean_U_cavity=bsxfun(@plus,mean_U,y*(CXX*diag((1-rates)./CXX_diag)));
        var_U_cavity=bsxfun(@plus,var_U,-((y*CXX).^2)*diag(1./CXX_diag));
        var_U_cavity(var_U_cavity<0)=0; %correct for numerical errors;

        if use_sampling
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
        MSE(end+1)=mean(abs(x(:)-x_prev(:))); %#ok    

        if show_progress
            disp(['MSE=' num2str(MSE(end))]);
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
            plot(MSE)
            ylabel('MSE');
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

    %     if length(MSE)>20
    %         if MSE(end-1)<MSE(end)
    %             disp('Liphshitz constant doubled')
    %             L=2*L;
    %             lambda_high=lambda_high*2;
    %             lambda_low=lambda_low*2;
    %             lambda=lambda*2;
    %         end
    %     end

        FISTA_cond=(MSE(end)>Tol_FISTA)||(iteration<30);

        temp=~~x(1:(end-N_stim),1:(end-N_stim));
        sparsity_measure=mean(temp(:));    
        cond=sparsity_measure<sparsity;
        relative_sparsity=abs(sparsity_measure-sparsity)/sparsity;
        loop_cond=( relative_sparsity>  Tol_sparse);
        if sparsity==1
            loop_cond=0;
        end

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


        sparsity_measure
        spar_array(end+1)=sparsity_measure; %#ok
        if cond
            lambda_high=lambda
        else
            lambda_low=lambda
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

