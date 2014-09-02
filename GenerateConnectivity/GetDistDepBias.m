function b=GetDistDepBias(a,spar)
% Finds the correct bias for a sigmoid distance dependent function given a
% in p(conn)=1/(1+exp(-(ax+b)))
% Inputs:
% a - gain of sigmoid
% spar - desired sparsity

    %requires solving a transcendental eqn so have to do it numerically
    syms bb;
    aa=num2str(a);
    as=num2str(a*spar);

    prblm=['exp(' aa '+bb)+1=exp(' as '+bb)+exp(' as ')'];
    b=double(solve(prblm));
    b=round(b*100)/100; %round

end
