function [EW, Eb]=LogisticRegression(X,Y,sparsity)
    m_x=mean(X,1);
    m_y=mean(Y,1);    
    if numel(X)<1e7
        CXX=(bsxfun(@plus,X,-m_x))'*(bsxfun(@plus,X,-m_x));    
    else
        CXX=zeros(size(X,2),size(X,2));
        for tt=1:size(X,1)
            CXX=CXX+(X(tt,:)-m_x)'*(X(tt,:)-m_x);        
        end
    end 
    CXY=(bsxfun(@plus,X,-m_x))'*(bsxfun(@plus,Y,-m_y));
    N_stim=0;
    pen_diag=1;
    warm=1;
    EW=EstimateA_L1_logistic_Accurate(CXX,CXY,m_y,sparsity,N_stim,pen_diag,warm);
     [amp, Eb]=logistic_ELL(m_y,m_x,EW,CXX,CXY);
     EW=EW*amp;

end

function EW=EstimateA_L1_logistic_Accurate(CXX,CXY,rates,sparsity,N_stim,pen_diag,warm)
% Algorithm Implements FISTA algorithm by Beck and Teboulle 2009 for L1 linear regression
% INPUTS: 
% CXX - covariance
% CXY - cross-covariance
% rates - mean firing rates
% sparsity - required sparsitiy level (percentage of non-zero enteries in EW)
% N_stim - number of stimuli
% pen_dial - 1 if we want to penalize diagonal entries; 0 otherwise
% warm - 1 if we want to do warm starts within FISTA; 0 otherwise
% initial - initial conditions
% OUTPUTS:
% EW: ML estimate of weights

%params
% iterations=300;
% tic

N=length(rates);
Tol_sparse=0.01; %tolerance for sparsity level
Tol_FISTA=1/N^2; %toleratnce for fista
max_iterations=1000;

%initialize FISTA
% if initial
%     x=initial;
% else
    x=0*CXY';
% end
y=x;

% this effectively removes weights in-going into a stimulus node
rates(rates<0)=0;
rates(rates>1)=0;
%%%

V=-diag(pi*(rates.*log(rates)+(1-rates).*log(1-rates))/8); %for the gradient
V(isnan(V))=0; %take care of 0*log(0) cases...
L=2*max(V(:))*max(eig(CXX)); %lipshitz constant

%initialize binary search
lambda_high=max(1e8,1e2*L); %maximum bound for lambda
lambda_low=min(1e-8,1e-2*L);  %minimum bound for lambda
loop_cond=1;  %flag for while llop

while  loop_cond %binary search for lambda that give correct sparsity level
    if sparsity==1
        lambda=0;
    else
        lambda=(lambda_high+lambda_low)/2;
    end

%%% FISTA
if ~warm
    x=0*CXY';
    y=x;
end
    
if pen_diag
    mask=ones(N);
else  
    mask=ones(N)-eye(N); % used to remove L1 penalty from diagonal
end
t_next=1;
FISTA_cond=1;
iteration=0;

while FISTA_cond
    t=t_next;
    x_prev=x;
    
    A=1./sqrt(1+pi*(y*CXX*y')/8);
    B=diag(A(eye(N)>0.5));
    u=y-(2/L)*(V*B*y*CXX-CXY');
      
    x=ThresholdOperator(u,mask.*lambda/L);
    if any(~isfinite(u(:)))
        error('non finite x!')
    end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);   
    iteration=iteration+1;    
    FISTA_cond=(mean(abs(x(:)-x_prev(:)))>Tol_FISTA);
    if (iteration>max_iterations);
        break
        warning('FISTA max iteration reached');
    end
end
%%% 

    temp=~~x(1:(end-N_stim),1:(end-N_stim));
    sparsity_measure=mean(temp(:));
    cond=sparsity_measure<sparsity;

    
    loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
    if sparsity==1
        loop_cond=0;
    end

    if cond
        lambda_high=lambda
    else
        lambda_low=lambda
    end
%     toc
end
        
EW=x;
end

function y = ThresholdOperator( x , lambda )
%THRESHOLDOPERATOR Summary of this function goes here
% This function implements the threshold operator for the l1/l2 group norm
% as specified on "convex optimization with sparsity-inducing norms", page 12, last equation

%   Detailed explanation goes here
% x - input
% y - output
% lambda - regularization constant

    temp=(1-lambda./abs(x));
    temp(isnan(temp))=1;
    y=temp.*x;
    y(temp<0)=0;
    



end


function [amp, bias]=logistic_ELL(m_y,m_x,EW,Cxx,Cxy)
% Given 
% m_y - mean output
% m_x -mean input
% CXX - covariance
% CXY - cross-covariance
% EW - estimated weight matrix
% This function outputs
% amp - amplitude correction to the rows of W
% bias - the bias

L_z=1e4; %length of sampling
N=size(Cxy,2);
amp=zeros(N,1);
bias=amp;
options = optimset('GradObj','on','Display','off','LargeScale','off');

disp('starting Logistic ELL step')
for kk=1:N %for each row in EW, find corrected amplitude and bias
    En=m_y;
    Enx=EW(kk,:)*(Cxy(:,kk)+m_x'*m_y);
    z=randn(L_z,1);
    x=sqrt(EW*Cxx*EW')*z+EW*m_x';
    [new_ab,~,exitflag]=fminunc(@(ab) twod_logistic_ELL_func(ab,Enx,En,x),[0;0],options); %OK    
%     disp(exitflag)    
    amp(kk)=new_ab(1); %gain
    bias(kk)=new_ab(2); %offset
%     if exitflag~=1
%        amp(kk)=1;
%        bias(kk)=NaN             
%     end
end

disp('Finished Logistic ELL step')
end

function [L,grad]=twod_logistic_ELL_func(ab,Enx,En,x)

a=ab(1); b=ab(2);
eaxb=exp(a*x+b);
L=a*Enx+b*En-mean(log(1+eaxb));
v=eaxb./(1+eaxb);
grad(1)=Enx-mean(x.*v);
grad(2)=En-mean(v);

L=-L;
grad=-grad;
end
