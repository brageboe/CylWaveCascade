function C_ij = calculate_mode_coupling_ij_coaxials(n, q, r_o, r_i, d_t, d_z, root_0n, root_0q, case_area)
% Calculates the mode coupling element between modes i=TM_{0n} and j=TM_{0q}
% for a junction between two circular coaxial waveguides with different
% outer radii r_o(1) (left side) and r_o(2) (right side) and a common inner
% radius r_i.
%
% See chapter 2.1 in report for the theory.
%
% Assumes 
% size(r_o) = 1 x 2 with values r_o = [b R]
% d_t: transversal wavenumbers; left-side of junction in d_t(n,1), right-side in d_t(q,2)
% case_area: cases (1) and (2) as described in chapter 2.1 of report (ver6)
    arguments
        n       % mode n
        q       % mode q
        r_o     % array length 2
        r_i     % single-valued
        d_t     % matrix size N x 2, where N=[1,2,..., max amount of roots/modes]
        d_z     % matrix size N x 2, where N=[1,2,..., max amount of roots/modes]
        root_0n % single-valued; n:th root of L_0
        root_0q % single-valued; q:th root of L_0
        case_area% single-valued
    end
    
    if root_0n == root_0q
        % Should never happen.
        error("Error in mode coupling calculation: modes n="+n+" and q="+q+". have equal roots, xi_0n="+root_0n+", xi_0q="+root_0q+".");
    end
    
    switch case_area
        case 1 
            % left-side of junction has r_o_1=b, right-side r_o_1=R
            r_o_1 = r_o(1); % b
            r_o_2 = r_o(2); % R
            
            L_1_n1 = bessel_L(1,n,r_i,r_i,root_0n);
            L_1_n2 = bessel_L(1,n,r_o_1,r_i,root_0n);
            
            % calculate L_m exclusive to this case:
            L_2_n2 = bessel_L(2,n,r_o_1,r_i,root_0n);
            L_1_q3 = bessel_L(1,q,r_o_1,r_i,root_0q);
            L_2_q2 = bessel_L(2,q,r_o_1,r_i,root_0q);
            term1_1 = root_0n * L_2_n2 * L_1_q3 * r_o_1 / r_i;
            term1_2 = root_0q * L_1_n2 * L_2_q2 * r_o_1 / r_i;
            term1 = term1_1 - term1_2;
        case 2 
            % left-side of junction has r_o_1=R, right-side r_o_1=b
            r_o_1 = r_o(2); % R
            r_o_2 = r_o(1); % b
            
            L_1_n1 = bessel_L(1,n,r_i,r_i,root_0n);
            L_1_n2 = bessel_L(1,n,r_o_1,r_i,root_0n);
            
            % calculate L_m exclusive to this case here:
            L_2_n2 = bessel_L(2,n,r_o_2,r_i,root_0n);
            L_1_q3 = bessel_L(1,q,r_o_2,r_i,root_0q);
            L_2_q2 = bessel_L(2,q,r_o_2,r_i,root_0q);
            term1_1 = root_0n * L_2_n2 * L_1_q3 * r_o_2 / r_i;
            term1_2 = root_0q * L_1_n2 * L_2_q2 * r_o_2 / r_i;
            term1 = term1_1 - term1_2;
        otherwise
            error("case_area can only take values 1 or 2.")
    end
    
% 	prod1 = 2*root_0n*root_0q*d_z(q,2) / (pi*(r_i^2)*d_t(n,1)*d_t(q,2)*abs(d_z(q,2))); %can simply with sign(d_z(q,2))
%   prod1 = 2*root_0n*root_0q*sign(d_z(q,2)) / (pi*d_t(n,1)*d_t(q,2)); % 25/10 NEEDS UPDATE for ver7
    prod1 = 2*root_0n*root_0q / (r_i*r_i*d_t(n,1)*d_t(q,2)*conj(d_t(q,2)));   
    prod2 = sqrt(abs(d_z(n,1))*abs(d_z(q,2))); 
    prod3 = 1/sqrt( (r_o_1^2)*L_1_n2^2 - (r_i^2)*L_1_n1^2 );
    
    L_1_q1 = bessel_L(1,q,r_i,r_i,root_0q);
    L_1_q2 = bessel_L(1,q,r_o_2,r_i,root_0q);
    prod4 = 1/sqrt( (r_o_2^2)*L_1_q2^2 - (r_i^2)*L_1_q1^2 ); % sign switch?
    
    L_2_n1 = bessel_L(2,n,r_i,r_i,root_0n);
    L_2_q1 = bessel_L(2,q,r_i,r_i,root_0q);
    term2_1 = root_0n * L_2_n1 * L_1_q1;
    term2_2 = root_0q * L_1_n1 * L_2_q1;
    term2 = term2_1 - term2_2;
    denom = (root_0n^2 - root_0q^2) / r_i^2;
    prod5 = (term1 - term2) / denom;
    
    C_ij = prod1 * prod2 * prod3 * prod4 * prod5;   % eq (2.11) or (2.13) depending on case
end