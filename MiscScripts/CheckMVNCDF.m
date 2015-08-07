L=1e4;
M=100;
Sigma=1;
y=zeros(M,1);
z=y;
mu=linspace(-10,10,M);
b=0.1;
xl=[-b 0];
xu=[b inf];
for ii=1:M
    x=mu(ii)+Sigma*randn(L,1);
    z(ii)=mean(normcdf(x*sqrt(pi/8)).*(x<b).*(x>-b));
    Sigma_mat=[Sigma.^2,Sigma.^2; Sigma.^2, Sigma^2+pi/8];
    y(ii)=mvncdf(xl,xu,mu(ii),Sigma_mat);   
end
plot(mu,z,mu,y)