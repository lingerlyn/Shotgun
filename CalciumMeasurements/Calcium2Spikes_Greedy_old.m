function [spikes, b ] = Calcium2Spikes_Greedy(Y,P)
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
 % initialize filter   
        Gx = @(x,mode) G_mat(x,mode,T,P{nn}.g,0,P{nn}.z);
        residual=Gx(Y(nn,:),2); %= G'y - multipling data with matched 
        m=mean(Gx(ones(T,1),1));
        num_spikes=round(T*mean(Y(nn,:),2)/m); %later modify this with exponential search
        
        % get delta response for G_trans
        delta=zeros(1,T);
        delta(round(T/2))=1;
        G_trans_delta=Gx(Gx(delta,1),2); 
        cutoff=round(5/(1-P{nn}.g)); %cut off range for delta response support
        support=round(T/2)+(-cutoff:cutoff);
        G_trans_delta=G_trans_delta(support);
        
        for ss=1:num_spikes
            [~,ind]=max(residual);    
            spikes(nn,ind)=1;
            support_shifted=ind+support-round(T/2);
            ind_cut=or(support_shifted<1,support_shifted>T);
            support_shifted(ind_cut)=[];
            G_trans_delta_temp=G_trans_delta;
            G_trans_delta_temp(ind_cut)=[];
            residual(support_shifted)=residual(support_shifted)-G_trans_delta_temp;
            plot(residual)
        end
        if ~mod(nn,floor(N/10))
            disp([num2str(100*nn/N,2) '%'])
        end
end
end