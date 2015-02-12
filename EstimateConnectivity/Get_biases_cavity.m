    Eb_prev=Eb;       
    mean_U=y*rates+z;
    var_U=diag(y*CXX*y');        
    Ef=1./(1+exp(-(mean_U./sqrt(1+pi*var_U/8)))); % use sigmoid_int here for possibly higher accuracy (but much slower computation. also asymptotic cases not working well)
    Eb=z-(2/L)*(Ef-rates);
      
    x=ThresholdOperator(u,mask.*lambda/L);
    if any(~isfinite(u(:)))
        error('non finite x!')
    end
    t_next=(1+sqrt(1+4*t^2))/2;
    y=x+((t-1)/t_next)*(x-x_prev);    
    z=Eb+((t-1)/t_next)*(Eb-Eb_prev); 