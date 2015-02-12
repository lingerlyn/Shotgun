function b = G_mat(x,mode,T,g,bas_est,varargin)

% bas_est       flag for estimating baseline
 %G = spdiags(ones(T,1)*[-flipud(g(:))',1],-length(g):0,T,T);
 switch length(varargin)
     case 0
         z=1;
     case 1      
         if isempty(varargin{1})
            z=1;
         else
            z=varargin{1};
         end
     otherwise
         error('too many arguments!')
 end

if mode == 1
    %b = G\x(1:T) + x(end);
    if bas_est
        b = filter(z,[1;-g(:)],x(1:T)) + x(end);
    else
        b = filter(z,[1;-g(:)],x(1:T));
    end
elseif mode == 2
    if bas_est
        b = [flipud(filter(z,[1;-g(:)],flipud(x)));sum(x)];
    else
        b = flipud(filter(z,[1;-g(:)],flipud(x)));
    end
    %b = [G'\x;sum(x)] ;
end