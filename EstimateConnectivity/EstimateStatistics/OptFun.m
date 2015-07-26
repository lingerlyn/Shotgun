function quality = OptFun( lambda,COV,regularization_mat, Tol, msg, maxIter ,true_W)
% This is a wrapper function for optimzation purpuses
%   Detailed explanation goes here
 [inv_COV_res, ~, ~,~, ~, ~] = QUIC('default', COV, lambda*regularization_mat, Tol, msg, maxIter);
 N=size(COV,1)/2;
 temp=inv_COV_res((N+1):end,1:N);
 W=-temp;
 quality=-corr(W(:),true_W(:));
 if isnan(quality)
     quality=inf;
 end

end

