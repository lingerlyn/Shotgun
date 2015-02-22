function [ quality,d_bins ] = GetQualityDistance( W,EW,centers )
%GETQUALITYDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    L=100;
    [~,dist]=GetDistProb(centers);
    [~,d_bins]=hist(dist,L);
    d_bins_size=mean(diff(d_bins));
    quality=zeros(L,4);
    for ii=1:L
        ind=(dist-d_bins(ii))<d_bins_size/2;
        [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W(ind),EW(ind) )
        quality(ii,:)=[R,correlation, zero_matching,sign_matching];
    end
end

