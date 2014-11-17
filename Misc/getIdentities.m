function idents=getIdentities(EW,tol)

    N=length(EW);
    idents=zeros(N,1);


     for ii=1:N
        nzs=EW(~~EW([1:ii-1 ii+1:N],ii),ii);
        fracEx=sum(nzs>0)/numel(nzs);
        fracIn=sum(nzs<0)/numel(nzs);

%         if fracIn==0 && fracEx>0
%             idents(ii)=1;
%         end
% 
%         if fracEx==0 && fracIn>0
%             idents(ii)=-1;
%         end

        if (fracEx-fracIn)>tol
            idents(ii)=1;
        end

        if (fracIn-fracEx)>tol
            idents(ii)=-1;
        end

    end         

         
end