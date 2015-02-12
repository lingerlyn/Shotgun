
function y=sigmoid_int(mu,Sigma)
% This function approximate the integral of (1./(1+exp(-x)))*N(x|mu,Sigma^2)
% by using log(1+exp(-x))= exp(x) for x<-b,  exp(-x) for x>b,
% and some 2D gaussian integral otherwise, and b 

% note: we can make this function faster by further dividing into cases.
% for example, if both mu and Sigma are small, we can just use Phi(mu/sqrt(Sigma.^2+8/pi)); 


%% The 
% there is some numeric problem with the products of diverging and decaying
% exponents - did not solve this yet
b=1.5;  % some parameter that was chosen by hand to minmize boundary discontinuity

I1=exp(mu+0.5*Sigma.^2).*Phi(-(b+mu+Sigma.^2)./Sigma);
I1(~isfinite(I1))=0;
I2_a=-exp(-mu+0.5*Sigma.^2).*Phi(-(b-mu+Sigma.^2)./Sigma);
I2_a(~isfinite(I2_a))=0;
I2_b=Phi(-(b-mu)./Sigma);
I2=I2_a+I2_b;

%% Using 5th order Taylor expansion of 1/(1+exp(-x)) around 0
Erf_p=erf((mu+b)./(sqrt(2)*Sigma));
Erf_m=erf((mu-b)./(sqrt(2)*Sigma));
Erf_m2=erf((-mu+b)./(sqrt(2)*Sigma));
Exp_p=exp(-(mu+b).^2./(2*Sigma.^2));
Exp_m=exp(-(mu-b).^2./(2*Sigma.^2));

M0=Phi((b-mu)./Sigma)-Phi((-b-mu)./Sigma);
M0(~isfinite(M0))=0;
M1=Sigma.*(Exp_p-Exp_m)/(sqrt(2*pi))+0.5*mu.*(Erf_p-Erf_m);
M1(~isfinite(M1))=0;
M3=(2*Sigma.*Exp_p.*(mu.^2+2*Sigma.^2-mu.*b+b.^2-(mu.^2 ...
    +2*Sigma.^2+mu.*b+b.^2).*exp(2*mu.*b./(Sigma.^2))) ...
    +mu.*(mu.^2+3*Sigma.^2).*sqrt(2*pi).*(Erf_p+Erf_m2))/(2*sqrt(2*pi));
M3(~isfinite(M3))=0;
M5=(exp(-(mu.^2+b^2)./(Sigma.^2)).*(2*Sigma.*(mu.^4+8*Sigma.^4 ...
    -mu.^3*b+4*Sigma.^2*b^2+b^4+mu.^2.*(9*Sigma.^2+b^2)...
        -mu.*(7*Sigma.^2*b+b^3)).*exp((mu-b).^2./(2*Sigma.^2)) ...
        -2*Sigma.*(mu.^4+8*Sigma.^4+mu.^3*b+4*Sigma.^2*b^2 ...
    +b.^4+mu.^2.*(9*Sigma.^2+b^2)+mu.*(7*Sigma.^2*b+b^3 ...
    )).* exp((mu+b).^2./(2*Sigma.^2)) )+mu.*(mu.^4+10*mu.^2.*Sigma.^2 ...
    +15*Sigma.^4)*sqrt(2*pi).*(Erf_p+Erf_m2))/(2*sqrt(2*pi));
M5(~isfinite(M5))=0;
    
I3=0.5*M0+M1/4-M3/48+M5/480;

ind=Sigma<1e-2;
if any(ind)
    I3(ind)=Phi(mu(ind)./sqrt(Sigma(ind).^2+8/pi)).*(mu(ind)>-b).*(mu(ind)<b);    
end

  
%% using the 2D gaussian integral calculation - here better b=2.5
% N=length(mu);
% xl=[-b 0];
% xu=[b inf];
% for nn=1:N
%     if abs(mu(nn)./Sigma(nn))>3*b
%         I3(nn)=0;
%     elseif Sigma(nn)<1e-2 %in this case mvncdf will be very slow, or fail
%          I3(nn)=Phi(sqrt(pi/8)*mu(nn))*(mu(nn)>-b)*(mu(nn)<b);    
%     elseif Sigma(nn)>1e2 %if Sigma is too large, mvncdf will fail. This is also faster
%         I3(nn)=Phi((b-mu(nn))/Sigma(nn))-Phi(-mu(nn)/Sigma(nn));
%     else
%         Sigma_mat=[Sigma(nn).^2,Sigma(nn).^2; Sigma(nn).^2, Sigma(nn)^2+8/pi];
%         assert(all(eig(Sigma_mat)>0),'negative eigenvalues!');
%         I3(nn)=mvncdf(xl,xu,mu(nn),Sigma_mat);
%     end
% end

y=I1+I2+I3;

ind=Sigma>1e2;
if any(ind)
    y(ind)=Phi(mu(ind)./Sigma(ind));    
end
end

function y=Phi(x)
thresh=8;
y=x;
ind1=abs(x)<thresh;
ind2=x>=thresh;
ind3=x<=-thresh;

y(ind1)=normcdf(x(ind1));
y(ind2)=1-exp(-x(ind2).^2/2)./(x(ind2)*sqrt(2*pi));
y(ind3)=-exp(-x(ind3).^2/2)./(x(ind3)*sqrt(2*pi));

end

