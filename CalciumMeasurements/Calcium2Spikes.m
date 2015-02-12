function [spikes, b ] = Calcium2Spikes(Y,P)
% input:
% Y - NxT calcium trace
% P - parameter struct of size N:
% P.g - poles polynomial coeffieints 
% P.Cb - bias
% P.sn - noise stdev
% P.z - zeros

% output:
% spikes - NxT de-convolved spikes
% b-  Nx1 estimated bias

[N,T] = size(Y);
spikes=zeros(N,T);
b=zeros(N,1);

for nn=1:N
    %% initialize parameters    
    
        Gx = @(x,mode) G_mat(x,mode,T,P{nn}.g,1,P{nn}.z);
%     Gx = @(x,mode) G_mat(x,mode,T,P.g,1);

    %% Re-weighted L1
%     addpath('utilities/spgl1-1.8'); 
    resparse_eps = 1e-6;
    N_resparse=3; % number of times to resparse the solution using re-weighted l1-norm

    for ii = 1:N_resparse
        if ii>1        
            wt = [1./(spikes(nn,:)'+resparse_eps);1];
            options = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0,'weights',wt);
        else
            options = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0);
        end    
        tic
        noise_residual=P{nn}.sn*sqrt(T);
        [sp_temp,~,~,~] = spg_bpdn( Gx, Y(nn,:)', noise_residual,options);
        toc
        spikes(nn,:) = sp_temp(1:T);
        b(nn)=sp_temp(end);
    end

%     spikes(nn,:)=spikes(nn,:)/mean(spikes(nn,~~spikes(nn,:))); %make maximum equl 1 so we  can interpet this as spikes
    end

end