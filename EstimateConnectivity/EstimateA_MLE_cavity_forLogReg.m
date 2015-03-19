function EW=EstimateA_MLE_cavity_forLogReg(XX,YX,mX,mY)

D=size(XX,1);
beta=1:D;
alpha=1;

finv=@(x)log(x./(1-x));
% remove zeros from Myx
% temp = unique(YX(:));
% min2 = temp(2); %second smallest value in Myx
% YX(YX==0)=min2;

Cxx=XX-mX*mX';
Cxx_inv=inv(Cxx);

Fas=bsxfun(@plus,(1-mX(beta)').^2,(pi/8)*((finv(bsxfun(@times,YX(alpha,beta),1./mX(beta)'))).^2)*diag(diag(Cxx(beta,beta))));
Bas=(1./Fas).*(finv(mY(alpha))*(1-mX(beta)'))*diag(diag(Cxx(beta,beta)));
Cas=(1./Fas).*(bsxfun(@plus,(finv(mY(alpha))).^2,-(finv(bsxfun(@times,YX(alpha,beta),1./mX(beta)'))).^2))*diag(diag(Cxx(beta,beta)).^2);
Das=-Bas+sqrt(Bas.^2-Cas);
Aas=Das./sqrt(1-(pi/8)*Das*(Cxx_inv(beta,beta))*Das');
if any(imag(Aas(:))~=0)
    Das=-Bas-sqrt(Bas.^2-Cas);
    Aas=Das./sqrt(1-(pi/8)*Das*(Cxx_inv(beta,beta))*Das');    
end
EW=Aas*Cxx_inv(beta,beta);
end
% imagesc(EW)
% colorbar