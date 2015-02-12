function AUROC = GetAUROC( predicted_label,label )
% Calculate the Area Under an ROC Curve
% label - Nx1 long binary vector (0/1)
% predicted_label - Nx1 long vector of probabilities
%   Detailed explanation goes here

L=100; %number of points to sample the ROC curve
TP=zeros(1,L);
FP=zeros(1,L);
TN=zeros(1,L);
FN=zeros(1,L);
thresholds=linspace(min(predicted_label),max(predicted_label),L);
for kk=1:L
    FP(kk)=mean(((predicted_label-thresholds(kk))>0).*(label<0.5));
    TP(kk)=mean(((predicted_label-thresholds(kk))>0).*(label>0.5));
    TN(kk)=mean(((predicted_label-thresholds(kk))<0).*(label<0.5));
    FN(kk)=mean(((predicted_label-thresholds(kk))<0).*(label>0.5));
end

TPR=TP./(TP+FN);
FPR=FP./(FP+TN);
delta_FPR=-diff(FPR,1);
AUROC=sum(TPR(2:end).*delta_FPR);

end

