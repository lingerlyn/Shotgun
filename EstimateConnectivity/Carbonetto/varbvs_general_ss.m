
function [alpha, mu, s_sq] = varbvs_general_ss (XtX, Xty, sigma_n_sq, sigma_s_sq, logodds, eta, alpha0, mu0,tolerance)

  % XtX - just what you would think
  % Xty - ditto
  % sigma_n_sq - response noise variance (scalar)
  % sigma_s_sq - slab prior variances (K x 1)
  % logodds - logit(a) where a is prior inclusions probability (Kx1)
  % eta - slab prior mean (K x 1)
  % alpha0 - starting value for posterior inclusion probs (for coord
  % ascent)
  % mu0 - starting value for posterior slab mean (for coord ascent)
  % tolerance - convergence criterion, a good choice is 1e-4
 
  K = size(XtX,1);

  % CHECK INPUTS.
  

  % X must be single precision.
  if ~isa(XtX,'single')
    XtX = single(XtX);
  end

  
  if isscalar(sigma_s_sq)
      sigma_s_sq = sigma_s_sq*ones(K,1);
  end
  
  if isscalar(eta)
      eta = eta*ones(K,1);
  end
      

  % LOGODDS must be a double precision column vector of length P.
  if isscalar(logodds)
    logodds = repmat(logodds,K,1);
  end

  d  = diag(XtX);
  
  alpha  = alpha0;
  mu     = mu0;
  
  s_sq = sigma_n_sq./(diag(XtX) + 1./sigma_s_sq);
  
  iter = 0;


  while true

    % Go to the next iteration.
    iter = iter + 1;
    
    
    % Save the current variational parameters and lower bound.
    params0 = [ alpha; alpha .* mu ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a forward or backward pass of the coordinate ascent updates.
    if isodd(iter)
      I = 1:K;
    else
      I = K:-1:1;
    end
    [alpha mu] = varbvsupdatematlab_general_ss(XtX,double(Xty),double(sigma_n_sq),...
        double(sigma_s_sq),double(logodds),...
        double(eta),double(d),double(alpha),double(mu),double(I-1));

    
    
    % CHECK CONVERGENCE.
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha .* mu ];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    
    if max(err) < tolerance
        disp(['converged after ' num2str(iter) ' iterations'])
      break
    end
    
  end
