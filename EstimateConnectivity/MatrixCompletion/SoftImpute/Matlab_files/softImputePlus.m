function [Mfinal,ranks]=softImputePlus(Minit,maxRank,nLambda,lambda)
% Uses soft-impute to find an svd that minimizes ||P(A)-B||-L*||A||_*, where P(.) is a projection that selects the 'non-missing' value 
% coordinates of B, the data matrix that contains missing values. This minimization is done for progressively smaller values of L until the 
% maxRank constraint is broken. The correct penalty value, lambda, is difficult to know a priori. This code starts with a low value 
% (giving a high rank) and increments the penalty until the lowest rank is achieved. For low (<10) ranks, setting lambda=.1 seems to work well.
% The 'plus' part does ordinary linear regression on the singular values to minimize the MSE between the low-rank approximation and the data matrix.
%
% Inputs
% Minit: observed data matrix
% maxRank: maximum rank of minimizing approximation
% nLambda: length of lambda sequence
% lambda: smallest value of lambda used (usually is set to 0.1, and the code adjusts lambda as needed).
%
% Outputs
% Mfinal: low-rank approximation
% ranks: sequence of ranks for each solution to the sequence of minimization problems

	if any(sum(Minit.^2,1)==0) || any(sum(Minit.^2,2)==0) %one row or column has all zeros
		error('All zero column or row');
	end

	if isempty(lambda)
    	lambda=.1; %pick a small lambda.
	end

    OPTS.NUMBER_GRID=nLambda;
    OPTS.MAX_RANK=maxRank;
    loopFailed=1;

    while loopFailed
        
        [Mat,Glr_mat_u,Glr_mat_d,Glr_mat_v,ranks,lambda,loopFailed]=...
            soft_impute_path_ps(sparse(Minit),lambda,OPTS);
        %soft_impute spits out the first lambda that gives a matrix exceeding maxRank, so just increment that lambda
        lambda=lambda+.05;   

    end


%Post-processing: Do ordinary linear regression on the singular values
    
    [u,s,v]=svd(double(Mat));
    
    estRank=sum(diag(Glr_mat_d)~=0);
 
    nonzeros=Minit~=0;
    %b is design matrix
    b=zeros(numel(Minit),estRank);
    for i=1:estRank
        temp=u(:,i)*v(:,i)'.*nonzeros;
        b(:,i)=temp(:);
    end
    
    %can reduce b
    B=b(nonzeros(:),:);
    
    %and y
    y=Minit(nonzeros(:));
    
    alpha=(B'*B)\(B'*y);
   
    Mfinal=u*diag([alpha;zeros(length(Minit)-estRank,1)])*v';
    

end
