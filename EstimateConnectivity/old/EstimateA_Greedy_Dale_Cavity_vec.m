function EW=EstimateA_Greedy_Dale_Cavity_vec(Cxx,Cxy,rates,spar,lambda,M,idents)


%add penalties later.
N=length(Cxx);

h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy
m=rates;
finv=@(x)log(x./(1-x));

EW=zeros(N);

%just add diagonal entries immediately
S=(1:N)';
Sc=ones(N,1)*(1:N); Sc=Sc';
Sc=Sc(~eye(N)); Sc=reshape(Sc,[N,N-1]);

for j=1:N %vectorize...?
    ew=EstimateA_MLE_cavity2(j,S(j,:),Cxx,Cxy,rates);
    ew(ew>0)=0; 
    EW(j,S(j,:))=ew;
end

while (nnz(EW)/N^2)<spar
    nS=size(S,2);
    disp(nS);
    
    wSc=zeros(N,size(Sc,2),2);
    LSc=zeros(N,size(Sc,2),2);
    
    for bb=1:size(Sc,2)
         
        b=Sc(:,bb);
        
        %get indices of W that correspond to the weights we have
        temp=[repmat( (1:N)',[nS,1]) S(:)];
        Sidx=N.*(temp(:,2)-1)+temp(:,1);
        
        %retrieve (b,a) index...
        temp=[b (1:N)'];
        cxyidx=N*(temp(:,2)-1)+temp(:,1);
        
%         const=(log( m(b)/(Cxy(b,a)+rates(a)*rates(b))-1))^2;       
        const=(log( m(b)./(Cxy(cxyidx)+rates.*rates(b))-1)).^2;
        
        temp=[b,b];
        bdiagidx=N.*(temp(:,2)-1)+temp(:,1);        

%         constt=((1-m(b))/Cxx(b,b))^2+const*pi/(8*Cxx(b,b));
        constt=((1-m(b))./Cxx(bdiagidx)).^2+const*pi./(8*Cxx(bdiagidx));
        
        %looping this part might be hard...
        %get W part that has current support
        WW=reshape(EW(Sidx),[N,nS]);
        
        
        wm=diag(WW*m(S)');
        
        temp=[repmat( b,[nS,1]) S(:)];
        Cxxidx=N.*(temp(:,2)-1)+temp(:,1);
        Cxxx=reshape(Cxx(Cxxidx),[N,nS]);
        
        
        wc=diag(WW*Cxxx');
        
        quadprod=zeros(N,1);
        for j=1:N
            quadprod(j)=EW(j,S(j,:))*Cxx(S(j,:),S(j,:))*EW(j,S(j,:))';
        end

        b_=2*m(b).*wm+2*constt.*Cxx(bdiagidx).*wc+2*(1-m(b))./Cxx(bdiagidx).*(m(b).*wc+Cxx(bdiagidx).*wm)-const*pi/4.*wc;
%         c_=wm.^2+constt.*wc.^2+2*(1-m(b))./Cxx(bdiagidx).*(wm.*wc)-const*pi/8.* w(notb)'*Cxx(notb,notb)*w(notb)-const;   
        c_=wm.^2+constt.*wc.^2+2*(1-m(b))./Cxx(bdiagidx).*(wm.*wc)-const*pi/8.* quadprod-const; 
        a_=m(b).^2+constt.*Cxx(bdiagidx).^2+2*(1-m(b)).*m(b)-const.*pi/8.*Cxx(bdiagidx);
        
%         solns= (-b_*ones(1,2) + [sqrt(b_^2-4*a_*c_) -sqrt(b_^2-4*a_*c_)])/(2*a_);
        solns= (-b_*ones(1,2)+ [sqrt(b_.^2-4*a_.*c_) -sqrt(b_.^2-4*a_.*c_)])./(2*a_*ones(1,2));
        wSc(:,bb,:)=solns; %make sure this lines up... X

        %look at resulting likelihoods of solutions found
        for ii=1:2
            
            ww=EW;
%             ww(bidx)=wSc(:,ii,bb);
            
            %need to index the bs
            temp=[(1:N)' b];
            widx=N.*(temp(:,2)-1)+temp(:,1);
            ww(widx)=solns(:,ii); %simpler...
%             wnonzero=ww(Sidx);
            
            %have to loop this part...
            quadprods=zeros(N,1);
            for j=1:N
                quadprods(j)=ww(j,[S(j,:) b(j)])*Cxx([S(j,:) b(j)],[S(j,:) b(j)])*ww(j,[S(j,:) b(j)])';
            end
            

            
%             LSc(:,ii,bb)=wnonzero'*Cxy(Sidx,a)-h(a)*sqrt(1+pi/8*wnonzero'*Cxx(Sidx,Sidx)*wnonzero)-lambda*(wSc(bb,ii)-M(a,b))^2;
            temp=[(1:N)' b];
            midx=N.*(temp(:,2)-1)+temp(:,1);
%             MM=reshape(M(midx),[N,1]);
            MM=M(midx);
            
            temp=[repmat( (1:N)',[nS+1,1]) [S(:);b]];
            Sbidx=N.*(temp(:,2)-1)+temp(:,1);
            WW=reshape(ww(Sbidx),[N,nS+1]);
            
            Cxyt=Cxy';
            Cxytt=reshape(Cxyt(Sbidx),[N,nS+1]);
            
%             disp(bb)

            LSc(:,bb,ii)=diag(WW*Cxytt')-h.*sqrt(1+pi/8*quadprods)-lambda*(squeeze(wSc(:,bb,ii))-MM).^2;
            
        end
    end
    
    
    %select only from correct + or - idents
%     identss=repmat(idents(Sc)',[N,1,2]);
    identss=repmat(idents(Sc),[1,1,2]);
    obj=LSc-identss.*(10^50*(sign(wSc)~=identss));
%     obj=LSc-repmat(abs(idents(Sc)),[1,2]).*(10^50*(sign(wSc)~=repmat(idents(Sc),[1,2])));
%     maxval=max(obj(:));
    [~,rows]=max(max(obj,[],3),[],2);
    
    temp=[(1:N)' rows];
    Scidx=N.*(temp(:,2)-1)+temp(:,1);
    
    S=sort([S Sc(Scidx)],2);
    for j=1:N %vectorize...?
        EW(j,S(j,:))=EstimateA_MLE_cavity2(j,S(j,:),Cxx,Cxy,rates);
    end

end


end