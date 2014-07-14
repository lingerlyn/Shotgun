function [L,grad]=twod_logistic_ELL(ab,Enx,En,x);

a=ab(1); b=ab(2);
eaxb=exp(a*x+b);
L=a*Enx+b*En-mean(log(1+eaxb));
v=eaxb./(1+eaxb);
grad(1)=Enx-mean(x.*v);
grad(2)=En-mean(v);

L=-L;
grad=-grad;