clear all
close all
clc

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Presentation&Figures\';

for kk=1:2
    if kk==1
        load('Run_N=50_obs=1_T=500000_FullyObservedGLM.mat');
        save_name='GLM_ELL_vs_Accurate=500000.pdf';
    elseif kk==2
        load('Run_N=50_obs=1_T=10000_FullyObservedGLM.mat');
        save_name='GLM_ELL_vs_Accurate=10000.pdf';
    end

    mi=min(W(:));ma=max(W(:));

    figure(kk)

    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind);
    hold all
    scatter(W(:),EW(:),'.')
    hold all
    scatter(W(:),EW2(:),'.')
    % hold all
    % scatter(W(:),glassoEW(:),'.')
    legend('x=y','Accurate','Approx.','location','northwest')
    xlabel('True weights')
    ylabel('Estimated weights')
    [R,C,Z,S] = GetWeightsErrors( W,EW );
    [R2,C2,Z2,S2] = GetWeightsErrors( W,EW2 );

    title({[' Accurate R =' num2str(R) ', Approx. R =' num2str(R2)]; ...
         [' Accurate C =' num2str(C) ', Approx. C =' num2str(C2)]; ...
         [' Accurate Z =' num2str(Z) ', Approx. Z =' num2str(Z2) ];...
         [' Accurate S =' num2str(S) ', Approx. S =' num2str(S2) ]});
    hold off    
Export2Folder(save_name,target_folder);
end
