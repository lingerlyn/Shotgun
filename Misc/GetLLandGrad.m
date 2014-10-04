function [L,grad]=GetLLandGrad(ww,params)

    N=length(params.Cxx);
    Cxx=params.Cxx;
    Cxy=params.Cxy;
    h=params.h; %entropy
    supp=params.supp; %support

%     ww=zeros(N,1);
%     ww(supp)=w;
    
    Cxx=Cxx(supp,supp); %do it faster this way.
    Cxy=Cxy(supp);



    L=ww'*Cxy-h*(1+pi/8*ww'*Cxx*ww)^(1/2);
    
    grad=Cxy-1/2*h*(pi/8*Cxx*ww)*(1+pi/8*ww'*Cxx*ww)^(-1/2);
%     grad=grad(supp); %select the weights being changed.
    
    L=-L;
    grad=-grad;

end

