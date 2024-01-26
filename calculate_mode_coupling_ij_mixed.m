function C_ij = calculate_mode_coupling_ij_mixed(n, q, r_o, r_i, d_t, d_z, root_0n, root_0q, case_area)
% Calculates the mode coupling element between modes i=TM_{0n} and j=TM_{0q}
% for a junction between two circular waveguides, one having an inner
% conductor (coaxial) while the other does not. The two waveguides are
% assumed to share a common outer radius.
%
% See chapter 2.2 in report for the theory.
%
% Assumes 
% case_area: cases (1) and (2) as described in chapter 2.2 of report (ver7)
% Have made use of the mathematical relations
%   * L_{-i}(x) = -L_{i}(x)
%   * J_{-i}(x) = -J_{i}(x)
%   where i=0,1,2,...
%
    arguments
        n       % mode n
        q       % mode q
        r_o     % single-valued; shared outer radii
        r_i     % single-valued; inner radius of coax conductor
        d_t     % matrix size N x 2, where N=[1,2,..., max amount of roots/modes]
        d_z     % matrix size N x 2, where N=[1,2,..., max amount of roots/modes]
        root_0n % single-valued; n:th root of L_0 (coax, case 1) or J_0 (circ. wg, case 2)
        root_0q % single-valued; q:th root of J_0 (circ. wg, case 1) or L_0 (coax, case 2) 
        case_area% single-valued
    end
    
    if root_0n == root_0q
        % Should never happen.
        error("Error in mode coupling calculation: modes n="+n+" and q="+q+". have equal roots, xi_0n="+root_0n+", xi_0q="+root_0q+".");
    end
    
    switch case_area
        case 1 % eq (2.14)
            % propagate from a coaxial to a circular waveguide       
            L_1_n1_ri = bessel_L(1, n, r_i ,r_i, root_0n);     
            L_1_n1_ro = bessel_L(1, n, r_o, r_i, root_0n);   
            J_1_q2 = besselj(1, root_0q);               
            prod3 = 1/sqrt( (-r_i^2 * L_1_n1_ri^2 + r_o^2 * L_1_n1_ro^2) * J_1_q2);
            
            L_2_n1_ro = bessel_L(2, n, r_o, r_i, root_0n);
            J_1_q2_ro = besselj(1, root_0q*r_o/r_i);
            J_2_q2_ro = besselj(2, root_0q);
            term1_1 = root_0n * L_2_n1_ro * J_1_q2_ro * r_o / r_i;
            term1_2 = root_0q * L_1_n1_ro * J_2_q2_ro; % (*r_o/r_o=1)
            term1 = term1_1 - term1_2;
            
            L_2_n1_ri = bessel_L(2, n, r_i, r_i, root_0n);
            J_1_q2_ri = besselj(1, root_0q*r_i/r_o);
            J_2_q2_ri = besselj(2, root_0q*r_i/r_o);
            term2_1 = root_0n * L_2_n1_ri * J_1_q2_ri; % (*r_i/r_i=1)
            term2_2 = root_0q * L_1_n1_ri * J_2_q2_ri * r_i / r_o;
            term2 = term2_1 - term2_2;
            
            denom = (root_0n * root_0n) / (r_i * r_i) - (root_0q * root_0q) / (r_o * r_o);
        case 2 % eq (2.16)
            L_1_q2_ri = bessel_L(1, q, r_i, r_i, root_0q);   
            L_1_q2_ro = bessel_L(1, q, r_o, r_i, root_0q);     
            J_1_n1 = besselj(1, root_0n);              
            prod3 = 1/sqrt( (-r_i^2 * L_1_q2_ri^2 + r_o^2 * L_1_q2_ro^2) * J_1_n1);
            
            L_2_q2_ro = bessel_L(2, q, r_o, r_i, root_0q);
            J_2_n1 = besselj(2, root_0n);
            term1_1 = root_0q * L_2_q2_ro * J_1_n1 * r_o / r_i;
            term1_2 = root_0n * L_1_q2_ro * J_2_n1; % (*r_o/r_o=1)
            term1 = term1_1 - term1_2;

            L_2_q2_ri = bessel_L(2, q, r_i, r_i, root_0q);
            J_1_n1_ri = besselj(1, root_0n*r_i/r_o);
            L_1_q2_ri = bessel_L(1, q, r_i, r_i, root_0q);
            J_2_n1_ri = besselj(2, root_0n*r_i/r_o);
            term2_1 = root_0q * L_2_q2_ri * J_1_n1_ri; % (*r_i/r_i=1)
            term2_2 = root_0n * L_1_q2_ri * J_2_n1_ri * r_i / r_o;
            term2 = term2_1 - term2_2;
            
            denom = (root_0q * root_0q) / (r_i * r_i) - (root_0n * root_0n) / (r_o * r_o);
        otherwise
            error("case_area can only take values 1 or 2.")
    end

    prod1 = 2*root_0n*root_0q / (r_i*r_o*r_o*d_t(n,1)*d_t(q,2)*conj(d_t(q,2)));  
    prod2 = sqrt(abs(d_z(n,1))*abs(d_z(q,2)));
    prod4 = (term1 - term2) / denom;
    
    C_ij = prod1 * prod2 * prod3 * prod4;   % eq (2.11) or (2.13) depending on case
end