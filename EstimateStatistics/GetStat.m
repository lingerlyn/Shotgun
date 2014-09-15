function [ CXX, CXY,W,rates,obs_count] = GetStat( sampled_spikes,observations,glasso,restricted_penalty,pos_def,sparsity,true_W)
% inputs:
% sampled_spikes -NxT observed spikes (with NaNs in unobsereved samples)
% glasso - flag that indicate whether or not use glasso
% restricted_penalty - flag that indicates whether or not to use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)
% sparsity-  the average nnz in the inv(CXX) matrix.  if empty, test using true matrix
% true_W - for  testing purposes

% outputs:
% Cxx  - NxN matrix of estimated spike covariance at a single timestep
% Cxy - NxN matrix of estimated spike cross-covariance  between two adjacent timesteps
% W - infered weights matrix
% rates - mean firing rates
% obs_count - observation count for each neuron


% internal parameters
Tol=1e-6; %numerical tolerance for glasso
msg=1; %verbosity of glasso algorithm
maxIter=100;  %glasso max iterations
Tol_sparse=1e-2; % tolerance for sparsity of W
thresh=5e6; %max size of spikes array for which we use effiecient computation of sufficient statistics

[N, T] = size(sampled_spikes);

%initialize the sufficient statistics arrays
observations=double(observations);
if N*T<thresh
    XX=sampled_spikes*sampled_spikes';
    XXn=observations*observations';
    XY=sampled_spikes(:,1:(end-1))*(sampled_spikes(:,2:end))';
    XYn=observations(:,1:(end-1))*(observations(:,2:end))';
else

    XX=zeros(N);
    XXn=zeros(N);
    XY=zeros(N);
    XYn=zeros(N);
    g = find(observations(:,1));
    sg = double(sampled_spikes(g,1));
    disp(['Gathering sufficient statistics...'])
    for tt = 2:T
        f = g;
        sf = sg;
        g = find(observations(:,tt));
        sg = double(sampled_spikes(g,tt));

        XX(f,f)=XX(f,f)+sf*sf';
        XXn(f,f)=XXn(f,f)+1;

        XY(f,g)=XY(f,g)+sf*sg';
        XYn(f,g)=XYn(f,g)+1;
        if ~mod(tt,T/10)
            disp([num2str(100*tt/T,2) '%'])
        end
    end
end

mY=sum(sampled_spikes,2);
mYn=sum(observations,2);

%% Calculate sufficient stats

rates=mY./(mYn+eps); %estimate the mean firing rates
CXX=XX./(XXn+eps)-rates*rates'; %estimate the covariance (not including stim terms for now)
% CXX((XXn<10))=0;%set elements to zero that haven't been observed sufficiently
CXY=XY./(XYn+eps)-rates*rates'; %estimate cross-covariance
% CXY((XYn<10))=0;
COV = [CXX CXY; CXY' CXX];

if(any(eig(COV)<0))
    disp('COV is not positive semidefinite;')
    if pos_def %positive semidefinite projection, so we won't have problems with glasso
        disp('correcting...');
        [v,d]=eig(COV);
        X0=v*spdiags(max(diag(d),0),0,2*N,2*N)*v';
        COV=X0;
    end
end

inv_COV=COV\eye(2*N);

%% Use glasso to sparsify stats
if glasso==1    
    disp('starting glasso...')



    
    addpath(fullfile('EstimateStatistics','QUIC')) %mex files for QUIC glasso implementation
    
    lambda_high=1; %maximum bound for lambda
    lambda_low=1e-4;  %minimum bound for lambda
    if restricted_penalty %more accurate, but much slower
        regularization_mat=ones(2*N);
        regularization_mat(1:N,1:N)=0; %remove penalty with upper left block
        regularization_mat(eye(2*N)>0.5)=0;  %remove penalty from diagonal
        regularization_mat([zeros(N), eye(N); zeros(N), zeros(N)]>0.5)=0;  %remove penalty from diagonal of upper right block
        regularization_mat([zeros(N), zeros(N); eye(N), zeros(N)]>0.5)=0;  %remove penalty from diagonal of lower left block         
    else
         regularization_mat=1;
    end
    
    
    if isempty(sparsity) % in case we "cheat" choouse optimal lambda
        func= @(x) OptFun( x,COV,regularization_mat, Tol, msg, maxIter ,true_W);
        lambda_opt=fminbnd(func,lambda_low,lambda_high);
        lambda_mat=lambda_opt*regularization_mat;
        [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter);
        
    else %set lambda according to estimated sparsity level of W, using binary search
        loop_cond=1;  %flag for while llop
        flag_first=1; % flag for first iteration of loop    
         
        while  loop_cond
            if sparsity==1
                lambda=0;
            else
                lambda=(lambda_high+lambda_low)/2;
            end
                lambda_mat=lambda*regularization_mat;
                if flag_first
                    [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter);
                else %use warm start
                    flag_first=0;
                    [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter,inv_COV_res, COV_res);
                end
            
                % W block
                temp=inv_COV_res((N+1):end,1:N);
                temp(eye(N)>0.5)=0;

                sparsity_measure=mean(~~temp(:)); 
                cond=sparsity_measure<sparsity;
                loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
                 if sparsity==1
                    loop_cond=0;
                 end
                 
            if cond
                lambda_high=lambda
            else
                lambda_low=lambda
            end
        end
    end
    
        CXY = COV_res(1:N,N+1:end);%%%%%
        CXX = COV_res(1:N,1:N);%%%%%
        W=-inv_COV_res((N+1):end,1:N);            
    else
        CXY = COV(1:N,N+1:end);%%%%%
        CXX = COV(1:N,1:N);%%%%%
        W=-inv_COV((N+1):end,1:N); % W estimate without glasso        
end

obs_count=mYn;
end

