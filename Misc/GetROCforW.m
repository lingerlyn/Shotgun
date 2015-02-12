function [AUROC, ROC] = GetROCforW( W,EW,sgn)
% Calculate the Area Under an ROC Curve
% inputs: 
% label - Nx1 long binary vector (0/1)
% predicted_label - Nx1 long vector of probabilities
% outputs: 
% AUROC -  Total area under ROC curve
% ROC- an 2xL array for ROC curve: first row false positive rate,  second row true positive rate.
% sgn - for which sign to calculate the ROC (1 - excitatatory weights, 0 - inhibitory weights)

L=1e4; %number of points to sample the ROC curve

W(sign(W)~=sgn)=0;
W=abs(W);
EW(sign(EW)~=sgn)=0;
EW=abs(EW);

TP=zeros(1,L);
FP=zeros(1,L);
TN=zeros(1,L);
FN=zeros(1,L);
thresholds=linspace(-1e-15,max([W(:);EW(:)]),L);

for kk=1:L
    positive=abs(W(:))>0;
    prediction=EW(:)-thresholds(kk)>0;
    
    FP(kk)=mean((prediction).*(~positive));
    TP(kk)=mean((prediction).*(positive));
    TN(kk)=mean((~prediction).*(~positive));
    FN(kk)=mean((~prediction).*(positive));
end

TPR=TP./(TP+FN);
FPR=FP./(FP+TN);
ROC=zeros(2,L);
ROC(1,:)=FPR;
ROC(2,:)=TPR;
delta_FPR=-diff(FPR,1)+eps;
AUROC=sum(TPR(2:end).*delta_FPR);

end

