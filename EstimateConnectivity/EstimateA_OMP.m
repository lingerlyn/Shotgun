function X=EstimateA_OMP(A_,B,spar,tol,lambda,M,rates)
% Does orthogonal matching pursuit until sparsity tolerance is reached

    
    N=size(A_,2);
    X=zeros(N);
    
    %loop over each row of W, i.e. x
    
    for i=1:N
        disp(i)
        x=zeros(N,1);
        S=[];
        Sc=(1:N)';

        b=B(:,i);
        r=b; %residuals
        
        A=A_*rates(i);
        while (mean(~~x)-spar)<tol

            AA=A(:,Sc);
            col_norms=sqrt(diag(AA'*AA));
            AAn=AA*diag(1./col_norms);

            prods=(AAn'*r).^2;
            zs=diag(col_norms.^(-2))*AA'*r;
            obj=prods-lambda*(zs-M(Sc)).^2;

            [~,idx]=max(obj);
            val=AA(:,idx)'*r/col_norms(idx)^2;

            S=sort([S;Sc(idx)]); %add the new dimension to the list.
            x(Sc(idx))=val;

            Sc(idx)=[];
            r=b-A*x;

        end

        X(i,:)=x;
    end



end
