function EW=EstimateA_Greedy_Dale_Cavity(Cxx,Cxy,rates,spar,lambda,M,idents)
%add penalties later.
N=length(Cxx);

h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy

EW=zeros(N);

m=rates;

finv=@(x)log(x./(1-x));

for a=1:N
    disp(a)
    w=zeros(N,1);
    S=[]; %support
    Sc=(1:N)'; %complement of support   
    
    
    while (nnz(w)/N)<spar
 
        wSc=zeros(length(Sc),2);
        LSc=zeros(length(Sc),2);
        
        %find optimal ws and evaluate the log-likelihoods
        for bb=1:length(Sc)
            b=Sc(bb);
            
            const=(log( m(b)/(Cxy(b,a)+rates(a)*rates(b))-1))^2;
            
            if isnan(const)
                disp('negative');
            end
            
            constt=((1-m(b))/Cxx(b,b))^2+const*pi/(8*Cxx(b,b));
                        
            notb=(1:N)'~=b;
            wm=w(notb)'*m(notb);
            wc=w(notb)'*Cxx(notb,b);
            
            a_=m(b)^2+constt*Cxx(b,b)^2+2*(1-m(b))*m(b)-const*pi/8*Cxx(b,b);
            b_=2*m(b)*wm+2*constt*Cxx(b,b)*wc+2*(1-m(b))/Cxx(b,b)*(m(b)*wc+Cxx(b,b)*wm)-const*pi/4*wc;
            c_=wm^2+constt*wc^2+2*(1-m(b))/Cxx(b,b)*(wm*wc)-const*pi/8* w(notb)'*Cxx(notb,notb)*w(notb)-const;      
            
            temp=roots([a_,b_,c_]);

            if length(temp)==1 %if there's only one solution, duplicate
                wSc(bb,:)=[temp temp];
            else
                wSc(bb,:)=temp;
            end

            %look at resulting likelihoods of solutions found
            for ii=1:2
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
        obj=LSc-repmat(abs(idents(Sc)),[1,2]).*(10^50*(sign(wSc)~=repmat(idents(Sc),[1,2])));
        maxval=max(obj(:));

        [row,col]=find(obj==maxval);
        idx=row; %care about the row element

        if ismember(a,Sc) %if a diagonal weight remains, pick it instead
            idx=find(a==Sc); %since we know it's negative (and on).
        end

        idx=idx(randi(length(idx))); %if the objs are equal pick a random one

        %update support
        S=sort([S;Sc(idx)]);
        Sc(idx)=[];

        Fas=(1-m(S)).^2+pi/8*diag(Cxx(S,S)).*(finv( (Cxy(S,a)+m(a)*m(S))./m(S))).^2;
        Bas=diag(Cxx(S,S))./Fas.*(1-m(S))*finv(m(a));
        Cas=diag(Cxx(S,S)).^2./Fas.*((finv(m(a)))^2-(finv( (Cxy(S,a)+m(a)*m(S))./m(S))).^2);
        Das=-Bas+sqrt(Bas.^2-Cas);
        Aas=Das./sqrt(1-pi/8*Das'*(Cxx(S,S)\Das));

        if any(~isreal(Aas)); keyboard; end
            
        w(S)=Cxx(S,S)\Aas;
    end
    EW(a,:)=w;
end
    
end