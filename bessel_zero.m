function x = bessel_zero(m, k, kind, arg0)
% Determine the first k zeroes of the m:th-order Bessel function of three
% different kinds,
% kind = 1: first kind J_m
% kind = 2: second kind Y_m
% kind = 3: L_m, the superposition of J_m and Y_m
%
% arg0 is only needed for case 3, L_m.
    arguments
        m
        k
        kind = 1;
        arg0 = 0;
    end

    k3 = 3*k;
    x = zeros(k3,1);
    for j=1:k3
        
        % Initial guess of zeros 
        x0 = 1 + sqrt(2) + (j-1)*pi + m + m^0.4;
        
        % Do Halley's method
        x(j) = findzero(m, x0, kind, arg0);
        if x(j) == inf
            error('Bad guess.');
        end
    end

    x = sort(x);
    dx = [1;abs(diff(x))];
    x = x(dx>1e-8);
    x = x(1:k);
end

function x = findzero(n, x0, kind, arg0)
    n1 = n+1;     
    %n2 = n*n;  
    tol = 1e-12; % Tolerance
    MAXIT = 100; % Maximum number of times to iterate
    err = 1; % Initial error
    iter = 0;
    while abs(err) > tol && iter < MAXIT
        switch kind
            case 1
                a = besselj(n,x0);    
                b = besselj(n1,x0);   
            case 2
                a = bessely(n,x0);
                b = bessely(n1,x0);
            case 3
                a = bessel_L(n, x0, arg0);
                b = bessel_L(n1, x0, arg0);
        end
                
        x02 = x0*x0;
        err = 2*a*x0*(n*a - b*x0)/(2*b*b*x02 - a*b*x0*(4*n + 1) + (n*n1+x02)*a*a);
        
        x = x0 - err;
        x0 = x;
        iter = iter+1;
    end
    if iter > MAXIT - 1
        warning("Failed to converge to within tolerance. Try a different initial guess");
        x = inf;    
    end
end