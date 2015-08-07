L=1e3;
M=100;
b=2;
Sigma=1;
y=zeros(M,1);
z=y;
mu=linspace(-10,10,M);

Erf_p=erf((mu+b)./(sqrt(2)*Sigma));
Erf_m1=erf((mu-b)./(sqrt(2)*Sigma));
Erf_m2=erf((-mu+b)./(sqrt(2)*Sigma));
Exp_p=exp(-(mu+b).^2./(2*Sigma.^2));
Exp_m=exp(-(mu-b).^2./(2*Sigma.^2));

M1=Sigma.*(Exp_p-Exp_m)/(sqrt(2*pi))+0.5*mu.*(Erf_p-Erf_m);
M3=(2*Sigma.*Exp_p.*(mu.^2+2*Sigma.^2-mu.*b+b.^2-(mu.^2 ...
    +2*Sigma.^2+mu.*b+b.^2).*exp(2*mu.*b./(Sigma.^2))) ...
    +mu.*(mu.^2+3*Sigma.^2).*sqrt(2*pi).*(Erf_p+Erf_m2))/(2*sqrt(2*pi));
M5=(exp(-(mu.^2+b^2)./(Sigma.^2)).*(2*Sigma.*(mu.^4+8*Sigma.^4 ...
    -mu.^3*b+4*Sigma.^2*b^2+b^4+mu.^2.*(9*Sigma.^2+b^2)...
        -mu.*(7*Sigma.^2*b+b^3)).*exp((mu-b).^2./(2*Sigma.^2)) ...
        -2*Sigma.*(mu.^4+8*Sigma.^4+mu.^3*b+4*Sigma.^2*b^2 ...
    +b.^4+mu.^2.*(9*Sigma.^2+b^2)+mu.*(7*Sigma.^2*b+b^3 ...
    )).* exp((mu+b).^2./(2*Sigma.^2)) )+mu.*(mu.^4+10*mu.^2.*Sigma.^2 ...
    +15*Sigma.^4)*sqrt(2*pi).*(Erf_p+Erf_m2))/(2*sqrt(2*pi));
    
y=0.5+M1/4-M3/48+M5/480;

for ii=1:M
    x=mu(ii)+Sigma*randn(L,1);
    z(ii)=mean(1./(1+exp(-x)).*(x>=-b).*(x<=b));    
end
plot(mu,z,mu,y);