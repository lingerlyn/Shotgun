function [spikes,relative_std_cell] = Calcium2Spikes_GreedyFast(Y,P)
% input:
% Y - NxT calcium trace
% P - parameter struct of size N:
% P.g - poles polynomial coeffieints 
% P.Cb - bias
% P.sn - noise stdev
% P.z - zeros

% output:
% spikes - NxT de-convolved spikes

% internal parameters - stop conditions
std_tol=0; % if reliatve std of residual is below this tolerance, stop
max_errors=1; %how many errors (and increase in std_res) are we allowed to have before algorithms stops 
show_progress=0;
amp=1; %amplitude of a single spike. estimate visually from calcium trace. Is there a way to get this automatically?

[N,T] = size(Y);
tt=1:T;
spikes=zeros(N,T);
baseline_est=0; %do not estimate basline
relative_std_cell={};


for nn=1:N
 % initialize filter   
        Gx = @(x,mode) amp*G_mat(x,mode,T,P{nn}.g,baseline_est,P{nn}.z);
        residual=Gx(Y(nn,:),2); %= G'y - multipling data with matched filter
        target_std=std(Gx(P{nn}.sn*randn(1,T),2)); 
%         m=mean(Gx(ones(T,1),1));
%         num_spikes=round(T*mean(Y(nn,:),2)/m); %later modify this with exponential search
        
        % get delta response for G_trans
        delta=zeros(1,T);
        delta(round(T/2))=1;
        G_trans_delta=Gx(Gx(delta,1),2); 
        cutoff=round(5/(1-P{nn}.g)); %currently assume AR(1) - need to generalize
        support=round(T/2)+(-cutoff:cutoff);
        G_trans_delta=G_trans_delta(support);
        std_res=std(residual);
        count=max_errors; %how many errors (and increase in std_res) are we allowed to have before algorithms stops 
        relative_std=[];
        loop_cond=1;

        while loop_cond
            if show_progress
                figure(1)
                a=2;b=1;
                subplot(a,b,1)
                plot(tt,Y(nn,:)./max(Y(nn,:)),tt,spikes(nn,:),'o');
                subplot(a,b,2)
                plot(tt,residual);
                pause
            end
            [~,ind]=max(residual);   
            support_shifted=ind+support-round(T/2);
            ind_cut=or(support_shifted<1,support_shifted>T);
            support_shifted(ind_cut)=[];
            G_trans_delta_temp=G_trans_delta;
            G_trans_delta_temp(ind_cut)=[];
%             var_residual_support=var(residual(support_shifted));
            residual(support_shifted)=residual(support_shifted)-G_trans_delta_temp;
%             var_residual_diff=var(residual(support_shifted))-var_residual_support;
            std_res_prev=std_res;
%             std_res=sqrt(std_res.^2+var_residual_diff*(length(support_shifted)/T));
            std_res=std(residual);  
            relative_std(end+1)=abs(1-std_res/target_std);
            if std_res_prev>std_res
                spikes(nn,ind)=spikes(nn,ind)+1;
                count=max_errors;
            else
                count=count-1;
            end
            
            if or(count<0,relative_std(end)<std_tol)
                loop_cond=0;
            end
        end
        relative_std_cell{end+1}=relative_std;
        
        
end
end