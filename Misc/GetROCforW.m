function [AUROC, ROC] = GetROCforW(error_rates)
% Get ROC and Calculate the Area Under an ROC Curve
% inputs: 
% quality - output from EstimateA_L1_logistic_cavity
% output: 
% AUROC - 2x1 Total area under ROC curve, for positive and negative weights
% ROC-  2x2xL array for ROC curve: dim1 - positive and negative weights,
% dim2 -  first row false positive rate,  second row true positive rate , dim3 - different values 

% semilogx(lambda_path,quality(:,5:6))

for ii=1:2
    kk=2*ii;
    x=error_rates(:,kk);
    y=error_rates(:,kk-1);
    [x,ind]=sort(x);
    y=y(ind);
    x=[0; x; 1]; %#ok
    y=[0; y; 1]; %#ok    
    ROC(ii,1,:)=x; %#ok
    ROC(ii,2,:)=y; %#ok
    AUROC(ii)=trapz(x,y); %#ok
end


end

