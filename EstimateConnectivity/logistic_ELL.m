function [amp, bias]=logistic_ELL(rates,EW,Cxx,Cxy)
% Given 
% CXX - covariance
% CXY - cross-covariance
% EW - estimated weight matrix
% This function outputs
% amp - amplitude correction to the rows of W
% bias - the bias

L_z=1e4; %length of sampling
N=length(rates);
amp=zeros(N,1);
bias=amp;
options = optimset('GradObj','on','Display','off','LargeScale','off');

disp('starting Logistic ELL step')
for kk=1:N %for each row in EW, find corrected amplitude and bias
    En=rates(kk);
    Enx=EW(kk,:)*(Cxy(:,kk)+rates*rates(kk));
    z=randn(L_z,1);
    x=sqrt(EW(kk,:)*Cxx*EW(kk,:)')*z+EW(kk,:)*rates;
    [new_ab,~,exitflag]=fminunc(@(ab) twod_logistic_ELL_func(ab,Enx,En,x),[1;-1],options); %OK    
%     disp(exitflag)    
    amp(kk)=new_ab(1); %gain
    bias(kk)=new_ab(2); %offset
    if exitflag~=1
       amp(kk)=1;
       bias(kk)=NaN;
    else
        
    end
end

disp('Finished Logistic ELL step')
end

function [L,grad]=twod_logistic_ELL_func(ab,Enx,En,x)

a=ab(1); b=ab(2);
eaxb=exp(a*x+b);
L=a*Enx+b*En-mean(log(1+eaxb));
v=eaxb./(1+eaxb);
grad(1)=Enx-mean(x.*v);
grad(2)=En-mean(v);

L=-L;
grad=-grad;
end

%% visualize objective (run from inside twod_logistic_ELL_func)
% a_v=-5:0.1:5;
% b_v=-5:0.1:0;
% 
% [A B]=meshgrid(a_v,b_v);
% y=zeros(1,1,length(x));
% y(1,1,:)=x;
% eaxb=exp(bsxfun(@plus,bsxfun(@times,A,y),B));
% L_mat=-A*Enx-B*En+mean(log(1+eaxb),3);
% 
% mesh(A,B,L_mat)

