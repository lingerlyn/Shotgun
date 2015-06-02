function EW=EstimateA_MLE_cavity_forLogReg(XX,YX,mX,mY)
% XX - inputs second moment matrix
% YX - output-input cross moment
% mX - mean inputs
% mY - mean outputs

D=size(XX,1);
beta=1:D;
alpha=1;

finv=@(x)log(x./(1-x));
% remove zeros from Myx
temp = unique(YX(:));
min2 = temp(2); %second smallest value in Myx
YX(YX==0)=min2;

Cxx=XX-mX*mX';
Cxx_inv=pinv(Cxx);
Cyx=YX-mY*mX';

Fas=bsxfun(@plus,(1-mX(beta)').^2,(pi/8)*((finv(bsxfun(@times,YX(alpha,beta),1./mX(beta)'))).^2)*diag(diag(Cxx(beta,beta))));
Bas=(1./Fas).*(finv(mY(alpha))*(1-mX(beta)'))*diag(diag(Cxx(beta,beta)));
Cas=(1./Fas).*(bsxfun(@plus,(finv(mY(alpha))).^2,-(finv(bsxfun(@times,YX(alpha,beta),1./mX(beta)'))).^2))*diag(diag(Cxx(beta,beta)).^2);
Das=-Bas-sqrt(abs(Bas.^2-Cas)); %added abs value here for numerical stability
ind=or(imag(Das)>0,sign(Cyx.*Das)<0);
Das(ind)=-Bas(ind)+sqrt(Bas(ind).^2-Cas(ind));
Aas=Das./sqrt(abs(1-(pi/8)*Das*(Cxx_inv(beta,beta))*Das')); %added abs value here for numerical stability
EW=(Aas*Cxx_inv(beta,beta));
end
% imagesc(EW)
% colorbar