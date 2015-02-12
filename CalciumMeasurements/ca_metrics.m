function [CR,I] = ca_metrics(spikes_foopsi,spiketimes,DF)
% inputs:
% spikes_foopsi - estimated spikes
% spiketimes - true spike times
% DF - vector of time bins sizes
% outputs:
% CR - correlations
% I - mutual information

al = spiketimes'*spikes_foopsi/norm(spikes_foopsi)^2;

ln = length(DF);
I = zeros(ln,1);
CR = zeros(ln,1);
for i = 1:length(DF)
    df = DF(i);
    spikes_true = filter(ones(df,1)/df,1,full(spiketimes(:)));
    spikes_pred = filter(ones(df,1)/df,1,full(spikes_foopsi(:)))*al;
    I(i) = 1/length(spiketimes)*(spikes_true'*log2((spikes_pred+eps)./(spikes_true+eps))) + mean(spikes_true) - mean(spikes_pred);
    CR(i) = corr(spikes_true(:),spikes_pred(:));
end