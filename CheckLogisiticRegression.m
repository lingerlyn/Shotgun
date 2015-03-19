clear all
% close all
clc

D=100;
N=1e5;
rho_x=1;
Sigma=eye(D)+rho_x*randn(D);
mu=randn(D,1);
x=double(bsxfun(@plus,mu,Sigma*randn(D,N))>0);
b=-4;
Sigma_w=1;
w=Sigma_w*randn(D,1);
z=w'*x+b;
y=double((1./(1+exp(-z))>rand(1,N)));

hist(z,100)

%% Estimate weights
addpath('EstimateConnectivity')
addpath('Misc')
mX=mean(x,2);
mean(mX)
mY=mean(y,2)
XX=(x*x')/N;
YX=(y*x')/N;
any(XX(:)==0)
any(YX(:)==0)

Ew=EstimateA_MLE_cavity_forLogReg(XX,YX,mX,mY);
Ew2=(YX-mY*mX')/(XX-mX*mX');
% B= mnrfit(x',y'+1);
% Ew2=-B(2:end);
mi=min([w(:) Ew(:)]);
ma=max([w(:) Ew(:)]);
diagonal=mi:0.01:ma;
plot(diagonal,diagonal,'-',w,Ew','.g',w,Ew2,'.r');
[R C]= GetWeightsErrors( w,Ew )
[R2 C2]= GetWeightsErrors( w,Ew2 )
