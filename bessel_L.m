function L_m = bessel_L(m, n, rho, r_i, root_0n, options)
% Calculates L_m( xi_mn,k * r_1 / r_2 )
% where xi_mn,k is the n:th root of L of order m in waveguide k={1,2}
% where k={1,2} represent left side (1) or right side (2) of the junction.
%
% Typically,    r_i = r_i, i.e. radius of inner conductor in coax.
%               rho = {r_{o,k}, r_i}, i.e. either radius of outer wall of waveguide k or r_i.       
%
% Limitation of this function: 
% Since all occurrences of xi_{mn} in our problem are the roots of order
% m=0, this function will assume this mode when calculating the roots.
    arguments
        m % order
        n % n:th mode 
        rho
        r_i
        root_0n % n:th root of L_0
        options.error = 1e-14
    end
    %roots = bessel_L_root(0, n, rho, r_i); % m = 0
    J_0 = besselj(0, root_0n);
    Y_0 = bessely(0, root_0n);
    K = J_0 / Y_0;  % constant
    
    arg = root_0n * rho / r_i;
    J_m = besselj(m, arg);
    Y_m = bessely(m, arg);
    
    L_m = J_m - (K * Y_m);
    if abs(L_m) < options.error
        L_m = 0;
    end  
end

