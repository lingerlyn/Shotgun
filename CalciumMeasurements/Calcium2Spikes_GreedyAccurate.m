function [spikes,relative_var_cell] = Calcium2Spikes_GreedyAccurate(Y,P)
% this alternative version of the greedy algorithm does not work well.
% here we have a differet stopping condition: if the MSE of the residual
% does not decrease we go to the next height point 

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
var_tol=0; % if reliatve var of residual is below this tolerance, stop
max_num_peaks=inf; %how many peaks should we check before algorithms stops 
show_progress=0;
amp=1; %amplitude of a single spike. estimate visually from calcium trace. Is there a way to get this automatically?
baseline_quantile=0.4; % should somewhat smaller then 0.5; 
baseline_time=10; % re-calculate baseline value after baseline_time peaks removed - larger then 1 only to improve speed (this is the most time consuming step in the algorithm if done at each iteration)

[N,T] = size(Y);
tt=1:T;
spikes=zeros(N,T);
baseline_est=0; %do not estimate basline
relative_var_cell={};


for nn=1:N
 % initialize filter   
        Gx = @(x,mode) amp*G_mat(x,mode,T,P{nn}.g,baseline_est,P{nn}.z);
        residual=Gx(Y(nn,:),2); %= G'y - multipling data with matched filter
        target_var=var(Gx(P{nn}.sn*randn(1,T),2)); 
%         m=mean(Gx(ones(T,1),1));
%         num_spikes=round(T*mean(Y(nn,:),2)/m); %later modify this with exponential search
        
        % get delta response for G_trans
        delta=zeros(1,T);
        delta(round(T/2))=1;
        G_trans_delta=Gx(Gx(delta,1),2); 
        cutoff=round(5/(1-P{nn}.g)); %currently assume AR(1) - need to generalize
        support=round(T/2)+(-cutoff:cutoff);
        G_trans_delta=G_trans_delta(support);
        relative_var=[];
        loop_cond=1;
        var_res_total=var(residual);  
        [Maxima,MaxIdx] = findpeaks(residual,'SORTSTR','descend');
        count=0;

        while loop_cond
            if show_progress
                figure(1)
                a=2;b=1;
                subplot(a,b,1)
                plot(tt,Y(nn,:)./max(Y(nn,:)),tt,spikes(nn,:),'o');
                subplot(a,b,2)
                plot(tt,residual);
                pause(1e-6)
            end
            
            if count==0
                mu=quantile(residual,baseline_quantile);
                count=baseline_time;
            else
                count=count-1;
            end
            num_peaks  =min(max_num_peaks ,length(MaxIdx));

            for jj = 1:num_peaks  
                ind=MaxIdx(jj);
                support_shifted=ind+support-round(T/2);
                ind_cut=or(support_shifted<1,support_shifted>T);
                support_shifted(ind_cut)=[];
                G_trans_delta_temp=G_trans_delta;
                G_trans_delta_temp(ind_cut)=[];
                residual_check=residual(support_shifted);
                var_res_prev=mean((residual_check-mu).^2);
                residual_check=residual_check-G_trans_delta_temp;                
                var_res=mean((residual_check-mu).^2);  
                if var_res<var_res_prev
                    residual(support_shifted)=residual_check;
                    spikes(nn,ind)=spikes(nn,ind)+1;
%                     [Maxima_local,MaxIdx_local] = findpeaks(residual_check);

                    temp=diff(residual_check);
                    MaxIdx_local=find((temp(2:end)<=0).*(temp(1:end-1)>=0))+1;
                    Maxima_local=residual_check(MaxIdx_local);

                    MaxIdx_local=MaxIdx_local+support_shifted(1)-1;
                    ind_remove=ismember(MaxIdx,support_shifted);
                    MaxIdx(ind_remove)=[]; Maxima(ind_remove)=[]; 
                    [Maxima,sort_ind]=sort([Maxima,Maxima_local],'descend');
                    temp=[MaxIdx,MaxIdx_local];
                    MaxIdx=temp(sort_ind);
                    break
                else
                    if jj==num_peaks
                        loop_cond=0;
                    end
                end
            end
            
            var_res_total=var_res_total+var_res-var_res_prev;
            relative_var(end+1)=abs(1-var_res_total/target_var); %#ok 
            
%             if or(relative_var(end)<var_tol,any(residual<0))
%                 loop_cond=0;
%             end
        end
        relative_var_cell{end+1}=relative_var;%#ok 
        
        
end
end