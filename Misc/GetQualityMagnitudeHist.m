function [ quality,m_bins ] = GetQualityMagnitudeHist( W,EW)
%GETQUALITYDISTANCE Summary of this function goes here
%   Detailed explanation goes here
    L=30;
    [~,m_bins]=hist(W(:),L);
    m_bins_size=mean(diff(m_bins));
    quality=nan(L,2);
    for ii=1:L
        ind=abs(W-m_bins(ii))<=m_bins_size/2;
        if sum(ind(:))>30 %enough samples for histrogram
             [~, ~, zero_matching, sign_matching] = GetWeightsErrors( W(ind),EW(ind) );
%              correlation=mean(W(ind).*EW(ind))/sqrt(mean((W(ind).^2).*(EW(ind).^2)));
%              R=1-sqrt(mean((W(ind)-EW(ind)).^2)./std(W(:))^2);
%               correlation=mean((W(ind)-EW(ind)).^2);
             quality(ii,:)=[ zero_matching, sign_matching];        
        end
    end
end

