function S = scattering_matrix_mixed(freq, radius_inner, radius_outer, length_1, length_2, number_of_modes, junction_case, options)
    % Calculates S-matrix between a junction of a coaxial waveguide and a
    % circular waveguide, including the propagator matrix from the junction
    % over the length of the second waveguide. Can optionally include the
    % propagator matrix up to the junction as well. 
    arguments
        freq
        radius_inner
        radius_outer
        length_1
        length_2
        number_of_modes
        junction_case % 1 or 2, equivalent to cases discussed in section [2.2]
        options.print_cutoff = false
        options.propagator_geometry1 = false
        options.propagator_geometry2 = true
    end
    %% Global parameters
    % Geometry
    r_i = radius_inner; % radius of inner conductor in coaxials
    r_o = radius_outer; % radius of outer waveguide walls R
    h_1 = length_1;
    h_2 = length_2;
    N = number_of_modes;

    % Physical constants
    mu_0 = 1.25663706212e-6;            % vacuum permeability
    eps_0 = 8.8541878128e-12;           % vacuum permittivity
    k2 = (2*pi*freq)^2 * mu_0 * eps_0;  % wavenumber squared

    %% Roots, wavenumbers and cutoff
    d_t = zeros(N,2); % transversal wavenumbers; left-side stored in d_t(:,1), right-side in d_t(:,2)
    switch junction_case
        case 1 % propagate from a coaxial to a circular waveguide.
            root_0n = bessel_L_root(0, N, r_o, r_i);    % geometry 1 = coax.
            root_0q = bessel_zero(0, N);                % geometry 2 = circ.
            d_t(:,1) = root_0n ./ r_i;
            d_t(:,2) = root_0q ./ r_o;
        case 2 % propagate from a circular waveguide to a coaxial.
            root_0n = bessel_zero(0, N);                % geometry 1 = circ.
            root_0q = bessel_L_root(0, N, r_o, r_i);    % geometry 2 = coax.
            d_t(:,1) = root_0n ./ r_o;
            d_t(:,2) = root_0q ./ r_i;
        otherwise
            error("Parameter junction_case must take value 1 or 2.")
    end
        
    % longitudinal wavenumbers; left-side stored in d_z(:,1), right-side in d_z(:,2)
    d_z = sqrt(k2 - d_t.^2);  

    % if 1st mode in either geometry is below cutoff, set S=0 and exit.
    if ~above_cutoff(1, d_z, 1) || ~above_cutoff(1, d_z, 2)
        S = zeros(2*N, 2*N);
        %disp("S=0 @ f = "+freq/1e9+" GHz") % only for testing/debugging.
        return
    end

    if options.print_cutoff
        for n=1:N % badly optimized for frequency sweeps; does not change wrt freq but prints for every freq iteration.
            f_c_1 = calculate_cutoff(n, d_t, 1); % cutoff freq for left-side geometry
            f_c_2 = calculate_cutoff(n, d_t, 2); % cutoff freq for right-side geometry
            disp("f_c_1("+n+") = "+f_c_1/1e9+"   ;   f_c_2("+n+") = "+f_c_2/1e9);
        end
    end
    
    %% 1st Propagator
    % Wave propagation from z=-length_1 to junction z=0
    zero_matrix = zeros(N,N);
    if options.propagator_geometry1
        P_1 = compile_propagator_matrix(N, 0, h_1, d_z, 1);
        P_1_macro = [zero_matrix P_1; P_1 zero_matrix];
        % Typically propagator up until junction is included in the
        % previous layer's scattering matrix.
        % Exception: direct comparisons with comsol may possibly require some finite distance before first junction. 
    end
    %% Junction
    % Junction between a circular wg and a circular coaxial wg.
    % See section [2.2] in theory document.

    % Diagonal matrix
    D_1 = compile_diagonal_matrix(d_z, 1);
    D_2 = compile_diagonal_matrix(d_z, 2); 

    % Coupling matrix
    C = zeros(N,N);
    for n=1:N
        for q=1:N
            if ~above_cutoff(n, d_z, 1) || ~above_cutoff(q, d_z, 2)
                % set to zero if either mode does not exist at this freq
                C(n,q) = 0;
%                 disp("C("+n+","+q+") = 0 @ f = "+freq/1e9+" GHz") % only for testing/debugging.
            else
                C(n,q) = calculate_mode_coupling_ij_mixed(n, q, r_o, r_i, d_t, d_z, root_0n(n), root_0q(q), junction_case);
%                 if n == q  % only for testing/debugging.
%                     disp("C("+n+","+q+") = "+C(n,q)+" @f="+freq/1e9+"GHz")
%                 end
            end
        end
    end

    % Scattering matrix
    S_junction = compile_scattering_matrix(D_1, D_2, C);

    %% 2nd Propagator
    % Wave propagation from junction (z=0) to end of geometry 2 (z=length_2)
    if options.propagator_geometry2
        P_2 = compile_propagator_matrix(N, h_1, h_1+h_2, d_z, 2);
        P_2_macro = [zero_matrix P_2; P_2 zero_matrix];
    end
    %% Total S-matrix
    S = S_junction;
    if options.propagator_geometry1
        S =  P_1_macro * S;
    elseif options.propagator_geometry2
        S = S * P_2_macro;      
    end

    if ~check_physical_realizability(S) 
        % If operating frequency < cutoff frequency then ignore,
        % otherwise print warning.
        % Start from highest mode 
        % (if freq>cutoff_highermode then also freq>cutoff_lowermode)
        for n=(N:-1:1)
            f_c_1 = calculate_cutoff(n, d_t, 1); % cutoff freq of mode n, for left-side geometry
            f_c_2 = calculate_cutoff(n, d_t, 2); % cutoff freq of mode n, for right-side geometry
            if freq > f_c_1 
                warning("Unphysical S-matrix @ frequency="+freq/1e9+"GHz. Cutoff in geometry 1 for mode "+n+" is "+f_c_1/1e9+"GHz.");
                break
            end
            if freq > f_c_2  
                warning("Unphysical S-matrix @ frequency="+freq/1e9+"GHz. Cutoff in geometry 2 for mode "+n+" is "+f_c_2/1e9+"GHz.");
                break
            end
        end
    end
end