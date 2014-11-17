x=-5:0.001:5;
c=sqrt(pi/8);

% c=sqrt(2*pi)*exp(b^2/2)/(2+exp(b)+exp(-b));
y=normcdf(c*x,0,1);
c=0.36;
y=(x>0)-((x>0)-normcdf(c*x,0,1)).*(exp(-abs(c*x)));
z=1./(1+exp(-x));

plot(x,y,x,z,':r')

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript';
Export2Folder(['WangApprox.eps'],target_folder) 

%% 
% mu=-10:0.01:10;
% s=1;
mu=1;
s=0:0.01:100;
W=bsxfun(@times,s,randn(100000,1));
X=bsxfun(@plus,mu,W);
y=mean(log(1./(1+exp(-X))),1);
z=sqrt(1+(pi*s.^2)/8).*log(1./(1+exp(-mu./(sqrt(1+(pi*s.^2)/8)))));

% plot(mu,y,mu,z)
plot(s,y,s,z)