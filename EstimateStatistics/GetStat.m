function [ CXX CXY ] = GetStat( sampled_spikes,glasso,sparsity,varargin)
% inputs:
% sampled_spikes - observed spikes
% outputs:
% Cxx  - NxN matrix of estimated spike covariance at a single timestep
% Cxy - NxN matrix of estimated spike cross-covariance  between two adjacent timesteps
% stat - struct with various parameters (add details later)
% sparsity: the average nnz in the inv(CXX) matrix. 
%          Calculate it using nnz(eye(N)-true_A*true_A')/(N^2); true_A
%          is true weight matrix


if length(varargin)>1
    err = MException('ResultChk:OutOfRange', ...
        'Resulting value is outside expected range');
    throw(err);
end




[N T] = size(sampled_spikes);

%initialize the sufficient statistics arrays
XX=zeros(N);
XXn=zeros(N);
XY=zeros(N);
XYn=zeros(N);
mY=zeros(N,1);
mYn=zeros(N,1);


%% Calculate sufficient stats

g = find(~isnan(sampled_spikes(:,1)));
sg = double(sampled_spikes(g,1));
for t = 2:T
    f = g;
    sf = sg;
    g = find(~isnan(sampled_spikes(:,t)));
    sg = double(sampled_spikes(g,t));
    
    XX(f,f)=XX(f,f)+sf*sf';
    XXn(f,f)=XXn(f,f)+1;
    
    XY(f,g)=XY(f,g)+sf*sg';
    XYn(f,g)=XYn(f,g)+1;
    
    mY(g)=mY(g)+sg;
    mYn(g)=mYn(g)+1;
end

m=mY./(mYn+eps); %estimate the mean firing rates
CXX=XX./(XXn+eps)-m*m'; %estimate the covariance (not including stim terms for now)
CXX(find(XXn<10))=0;%set elements to zero that haven't been observed sufficiently
CXY=XY./(XYn+eps)-m*m'; %estimate cross-covariance
CXY(find(XYn<10))=0;


if glasso==1
    COV = [CXX CXY; CXY' CXX];
    if(any(eig(COV)<0))
        disp('COV is not pos def; correcting...')
        [v,d]=eig(COV);
        X0=v*spdiags(max(diag(d),0),0,2*N,2*N)*v';
        COV=X0;
    end
    
    
    csvwrite('CXX.csv', COV);
    
    if length(varargin)==0
    csvwrite('sparsity.csv',sparsity);
    % Code to run the .R file. Will have to change it for different
    !unset DYLD_LIBRARY_PATH; /usr/bin/Rscript myscript.R   
    system('~/Dropbox/Min-Suraj/Research/final\ implementation/code_for_calcium_stuff/myscript.R') %Run Rcode
    end
    
    
    if length(varargin)==1
        rho = varargin{1};
        csvwrite('rho_vec.csv',rho);
        !unset DYLD_LIBRARY_PATH; /usr/bin/Rscript myscript2.R   
        system('~/Dropbox/Min-Suraj/Research/final\ implementation/code_for_calcium_stuff/myscript2.R') %Run Rcode
    end
        
    
    inv_COV = csvread('TEST_inv_cov_md_est_GL_mle.csv');%%%%%
    COV = inv_COV\eye(2*N);%%%%%%%
    CXY = COV(1:N,N+1:end);%%%%%
    CXX = COV(1:N,1:N);%%%%%
end

end

