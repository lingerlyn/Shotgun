function [ MPE, MSE, MLE ] = GetSpikeError(estimated_rate,eta,ind_check )
%GETMEANERROR Summary of this function goes here
 % estimated_rates - self-explantory
 % eta - spikes
 % ind - indicator for where to checkerror
 
 %out:
 % MPE - mean pred error
 % MSE - mean square error
 % MLE - mean loglikelihood error
%   Detailed explanation goes here

        pred_error=abs(round(estimated_rate)-eta);
        square_error=(estimated_rate-eta).^2;
        loglikelihood_error=-(eta.*log(estimated_rate)+(1-eta).*log(1-estimated_rate)); %minus loglikelihood      
                       
        MPE=mean(pred_error(ind_check));
        MSE=mean(square_error(ind_check));
        MLE=mean(loglikelihood_error(ind_check));


end

