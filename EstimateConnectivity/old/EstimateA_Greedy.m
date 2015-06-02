function EW=EstimateA_Greedy(Cxx,Cxy,rates,spar,lambda,M)
%add penalties later.
N=length(Cxx);

h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy


EW=zeros(N);

for a=1:N
    
    w=zeros(N,1);
    S=[]; %support
    Sc=(1:N)'; %complement of support   
   
    
    while (nnz(w)/N)<spar
 
        wSc=zeros(length(Sc),4);
        LSc=zeros(length(Sc),4);
        
        %find optimal ws and evaluate the log-likelihoods
        for bb=1:length(Sc)
            b=Sc(bb);

%         c=(pi/8*h(a))^2;
%         cc=Cxy(b,a)^2+4*lambda*Cxy(b,a)*M(a,b)+4*lambda^2*M(a,b);

%         if ~isempty(S)
%             
%             a4=pi/2*lambda^2*Cxx(b,b);
%             a3=pi*lambda^2*w(S)'*Cxx(S,b)+(-4*lambda*Cxy(b,a)-8*lambda^2*M(a,b))*pi/8*Cxx(b,b);
%             
%             a2=pi/2*lambda^2*w(S)'*Cxx(S,S)*w(S)+4*lambda^2+(-4*lambda*Cxy(b,a)-8*lambda^2*M(a,b))*pi/4*w(S)'*Cxx(S,b)+...
%                 cc*pi/8*Cxx(b,b)-c*Cxx(b,b);
%             
%             a1=(-4*lambda*Cxy(b,a)-8*lambda^2*M(a,b))*(pi/8*w(S)'*Cxx(S,S)*w(S)+1)+cc*pi/4*w(S)'*Cxx(S,b)-2*c*w(S)'*Cxx(S,b)*Cxx(b,b);
%             
%             a0=cc*(pi/8*w(S)'*Cxx(S,S)*w(S)+1)-c*(w(S)'*Cxx(S,b))^2;
%             
%         else
%             a4=pi/2*lambda^2*Cxx(b,b);
%             
%             a3=(4*lambda*Cxy(b,a)-8*lambda^2*M(a,b))*pi/8*Cxx(b,b);
%             
%             a2=cc*pi/8*Cxx(b,b)-c*Cxx(b,b);
%             
%             a1=4*lambda*Cxy(b,a)-8*lambda^2*M(a,b);
%             
%             a0=cc;
%             
%         end

%         wSc(bb,:)=roots([a4 a3 a2 a1 a0]);
    if ~isempty(S)
        a_=(Cxy(b,a)+2*lambda*M(a,b))/(pi/8*h(a));
        b_=-2*lambda/(pi/8*h(a));
        c=1+pi/8*(w(S)'*Cxx(S,S)*w(S));
        d=pi/4*w(S)'*Cxx(S,b);
        e=pi/8*Cxx(b,b);
        f=w(S)'*Cxx(S,b);
        g=Cxx(b,b);
    else
        a_=(Cxy(b,a)+2*lambda*M(a,b))/(pi/8*h(a));
        b_=-2*lambda/(pi/8*h(a));
        c=1;
        d=0;
        e=pi/8*Cxx(b,b);
        f=0;
        g=Cxx(b,b);
    end
        
        d0=a_^2*c-f^2;
        d1=2*a_*b_*c+2*a_^2*d-2*f*g;
        d2=b_^2*c+2*a_*b_*d+a_^2*e-g^2;
        d3=b_^2*d+2*a_*b_*e;
        d4=b_^2*e;
        
        if isnan(d0+d1+d2+d3+d4)|| isinf(d0+d1+d2+d3+d4)
            keyboard;
        end
        
        temp=roots([d4 d3 d2 d1 d0]);       
        %could have 2 or 4 solutions...
        if length(temp)==2
            temp=[temp;temp];
        end
        wSc(bb,:)=temp;
            
            for ii=1:4
                
                if ~isreal(wSc(bb,ii))
                    LSc(bb,ii)=-inf;
                else
                    
                    ww=w;
                    ww(b)=wSc(bb,ii);
                    LSc(bb,ii)=ww'*Cxy(:,a)-h(a)*sqrt(1+pi/8*ww'*Cxx*ww)-lambda*(ww(b)-M(a,b))^2;
                
                end
            end
                 
        end
                
        %now form a penalized objective
    obj=LSc;
    maxval=max(obj(:));

    [row,col]=find(obj==maxval);
    idx=row;
    idx=idx(randi(length(idx))); %if the objs are equal pick a random one
    
    %update support
    S=sort([S;Sc(idx)]);
    Sc(idx)=[];
        
    A=sqrt( (pi/8*h(a))^2 - (Cxy(S,a)'*(Cxx(S,S)\Cxy(S,a)))*pi/8);
    
    %find MLE restricted to support S
    www=(1/A*Cxx(S,S)+2*lambda*eye(length(S)))\(Cxy(S,a)+2*lambda*M(a,S)');

    
    w(S)=www;
    
    end
    EW(a,:)=w;
end


end


