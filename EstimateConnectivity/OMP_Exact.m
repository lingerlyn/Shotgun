function EW=OMP_Exact(Cxx,Cxy,spar,rates)%,M,lambda)
%add penalties later.
N=length(Cxx);

h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy


EW=zeros(N);

for a=1:N
    w=zeros(N,1);
    S=[]; %support
    Sc=(1:N)'; %complement of support   
   
    
    while (nnz(w)/N)<spar
 
        wSc=zeros(length(Sc),1);
        LSc=zeros(length(Sc),1);
        
        %find optimal ws and evaluate the log-likelihoods
        for bb=1:length(Sc)
            b=Sc(bb);

            aa=(Cxy(b,a)/(pi/8*h(a)))^2;

            a2=aa*pi/8*Cxx(b,b)-Cxx(b,b)^2; %quadratic coeff


            if ~isempty(S)
                a1=aa*pi/4*w(S)'*Cxx(S,b)-2*w(S)'*Cxx(S,b)*Cxx(b,b); %linear coeff
                a0=aa*(1+pi/8*w(S)'*Cxx(S,S)*w(S))-(w(S)'*Cxx(S,b))^2; %constant
            else
                a1=0;
                a0=aa;
            end

    %       %we have two solutions
            wSc(bb,1)=(-a1 + sqrt(a1^2-4*a2*a0))/(2*a2);
            wSc(bb,2)=(-a1 - sqrt(a1^2-4*a2*a0))/(2*a2);


            %plus solution
            ww1=w;
            ww1(b)=wSc(bb,1);
            LSc(bb,1)=ww1'*Cxy(:,a)-h(a)*sqrt(1+pi/8*ww1'*Cxx*ww1);

            %minus solution
            ww2=w;
            ww2(b)=wSc(bb,2);
            LSc(bb,2)=ww2'*Cxy(:,a)-h(a)*sqrt(1+pi/8*ww2'*Cxx*ww2);

        end
    
    [~,idx]=max(LSc(:));
    
    if idx>size(LSc,1) %a negative soln
        idx=idx-size(LSc,1);
    end
    

    %update support
    S=sort([S;Sc(idx)]);
    Sc(idx)=[];
    params.supp=S;
        
    A=sqrt( (pi/8*h(a))^2 - (Cxy(S,a)'*(Cxx(S,S)\Cxy(S,a)))*pi/8);
    
    %find MLE restricted to support S
    www=1/A*Cxx(S,S)\Cxy(S,a);
    w(S)=www;
    
    end
    EW(a,:)=w;
end


end


