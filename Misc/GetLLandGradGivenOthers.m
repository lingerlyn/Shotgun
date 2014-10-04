function [L,dd]=GetLLandGradGivenOthers(wj,params)
%Can speed this up, i.e., input products with fixed w entries.

    Cxx=params.Cxx;
    Cxy=params.Cxy;
    h=params.h;
    w=params.w;
    idx=params.idx; %index of w_j

    ww=w;
    ww(idx)=wj;

    L=ww'*Cxy-h*(1+pi/8*ww'*Cxx*ww)^(1/2);
    dd=Cxy(idx)-pi/8*h*Cxx(idx,:)*ww*(1+pi/8*ww'*Cxx*ww)^(-1/2);
    
    L=-L;
    dd=-dd;
    
end
