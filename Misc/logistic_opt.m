function [L,grad]=logistic_opt(ab,x,y)

    a=ab(1); b=ab(2);

    eaxb=exp(a*x+b);
    v=eaxb./(1+eaxb);

    L=(a*x+b)'*y-sum(log(1+eaxb));

    grad(1)=x'*(y-v);
    grad(2)=sum(y-v);

    L=-L;
    grad=-grad;



end