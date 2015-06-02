clear all
% close all
clc

N=100;
T=1e6;
rho_x=3;
Sigma=eye(N)+rho_x*randn(N);
mu=randn(N,1);
x=double(bsxfun(@plus,mu,Sigma*randn(N,T))>0);
b=0;
Sigma_w=0.3;
% w=Sigma_w*exprnd(1,[D 1]);
w=Sigma_w*randn(N,1);
z=w'*x-mean(w'*x)+b;
y=double((1./(1+exp(-z))>rand(1,T)));

hist(z,100)

%% Load external data
% x=1*(data>0)';
% y=(labels'+1)/2;

%% Estimate weights
addpath('EstimateConnectivity')
addpath('Misc')
mX=mean(x,2);
mean(mX)
mY=mean(y,2)
XX=(x*x')/T;
YX=(y*x')/T;
any(XX(:)==0)
any(YX(:)==0)

tic
Ew=EstimateA_MLE_cavity_forLogReg(XX,YX,mX,mY);
% sparsity=1;warm=1; use_sampling=0;
% [Ew, Eb,quality,error_rates,lambda_path]=EstimateA_L1_cavity_forLogReg(XX,YX,mX,mY,sparsity,warm,w,use_sampling);
time1=toc
tic
Ew2=(YX-mY*mX')/(XX-mX*mX');
% B= mnrfit(x',y'+1);
% Ew2=-B(2:end);
time2=toc

%% Plot
mi=min([w(:) Ew(:)]);
ma=max([w(:) Ew(:)]);
diagonal=mi:0.01:ma;
plot(diagonal,diagonal,'-',w,Ew,'.b',w,Ew2,'.r');
[R C]= GetWeightsErrors( w,Ew )
[R2 C2]= GetWeightsErrors( w,Ew2 )

title({[' EW 1-R^2 =' num2str(1-R^2) ', EW2 1-R^2  =' num2str(1-R2^2)]; ...
     [' EW C =' num2str(C) ', EW2 C =' num2str(C2)];...
     [' EW T =' num2str(time1) ', EW2 T =' num2str(time2)]});
