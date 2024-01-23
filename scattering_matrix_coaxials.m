function S = scattering_matrix_coaxials(freq, radius_inner, radius_outer_1, radius_outer_2, length_1, length_2, number_of_modes, options)
%%
% Calculates the total scattering matrix of two coaxials connected by a
% junction at frequency freq.
%
% Waveguide i (i=1,2) has length longitudinal length_i and transversal radius radius_outer_i.
% Both coaxials must have equal conductor/inner radii.
    arguments
        freq
        radius_inner
        radius_outer_1
        radius_outer_2
        length_1
        length_2
        number_of_modes
        options.print_cutoff = false
        options.propagator_geometry1 = false
    end
    %% Global parameters
    % Geometry
    r_i = radius_inner; % radius of inner conductor in coaxials
    r_o = [radius_outer_1, radius_outer_2]; % radius of outer waveguide walls [b, R]
    h_1 = length_1;
    h_2 = length_2;
    N = number_of_modes;

    % Physical constants
    mu_0 = 1.25663706212e-6;            % vacuum permeability
    eps_0 = 8.8541878128e-12;           % vacuum permittivity
    k2 = (2*pi*freq)^2 * mu_0 * eps_0;  % wavenumber squared
    %eta = sqrt(mu_0/eps_0);             % wave impedance (not used atm)

    %% Roots, wavenumbers and cutoff
    root_0n = bessel_L_root(0, N, r_o(1), r_i); % geometry 1
    root_0q = bessel_L_root(0, N, r_o(2), r_i); % geometry 2

    d_t = zeros(N,2); % transversal wavenumbers; left-side stored in d_t(:,1), right-side in d_t(:,2)
    d_t(:,1) = root_0n ./ r_i;
    d_t(:,2) = root_0q ./ r_i;
    
    d_z = sqrt(k2 - d_t.^2);  % longitudinal wavenumbers; left-side stored in d_z(:,1), right-side in d_z(:,2)

    % if 1st mode in either geometry is below cutoff, set S=0 and exit.
    if ~above_cutoff(1, d_z, 1) || ~above_cutoff(1, d_z, 2)
        S = zeros(2*N, 2*N);
        %disp("S=0 @ f = "+freq/1e9+" GHz") % only for testing/debugging.
        return
    end
    
%     for n=1:N
%         if ~above_cutoff(n, d_z, 1) || ~above_cutoff(n, d_z, 2)
%             S = zeros(2*N, 2*N);
%             return
%         end
%     end
    if options.print_cutoff
        for n=1:N % badly optimized for frequency sweeps; does not change wrt freq but prints for every freq iteration.
            f_c_1 = calculate_cutoff(n, d_t, 1); % cutoff freq for left-side geometry
            f_c_2 = calculate_cutoff(n, d_t, 2); % cutoff freq for right-side geometry
            disp("f_c_1("+n+") = "+f_c_1/1e9+"   ;   f_c_2("+n+") = "+f_c_2/1e9);
        end
    end
    
    %% 1st Propagator
    % Wave propagation from (z=0) to junction (z=h_1)
    zero_matrix = zeros(N,N);
    if options.propagator_geometry1
        P_1 = compile_propagator_matrix(N, 0, h_1, d_z, 1);
        P_1_macro = [zero_matrix P_1; P_1 zero_matrix];
        % Typically propagator up until junction is included in the
        % previous layer's scattering matrix.
        % Exception: direct comparisons with comsol may possibly require some finite distance before first junction. 
    end
    %% Junction
    % Junction between two circular coaxial waveguides with different outer radii.
    % See section [2.1] in theory document.
    % **case 1** 
    % Left side r_o(1) = b
    % Right side r_o(2) = R
    % **case 2** 
    % Left side r_o(1) = R
    % Right side r_o(2) = b
    if r_o(1) < r_o(2)
        coax_case = 1;
    elseif r_o(1) > r_o(2)
        coax_case = 2;
    else
        error("Unrecognized value type.")
    end
    
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
                C(n,q) = calculate_mode_coupling_ij_coaxials(n, q, r_o, r_i, d_t, d_z, root_0n(n), root_0q(q), coax_case);
%                 if n == q  % only for testing/debugging.
%                     disp("C("+n+","+q+") = "+C(n,q)+" @f="+freq/1e9+"GHz")
%                 end
            end
        end
    end

    % Scattering matrix
    S_junction = compile_scattering_matrix(D_1, D_2, C); % Calculating this with P-matrix to get S_L for the next junction(correct?)

    %% 2nd Propagator
    % Wave propagation from junction (z=h_1) to (z=h_1+h_2)
    P_2 = compile_propagator_matrix(N, h_1, h_1+h_2, d_z, 2);
    P_2_macro = [zero_matrix P_2; P_2 zero_matrix];

    %% Total S-matrix
    if options.propagator_geometry1
        S =  P_1_macro * S_junction * P_2_macro;
    else
        S = S_junction * P_2_macro;
    end

    if ~check_physical_realizability(S, print_warning=false) 
        % If operating frequency < cutoff frequency then ignore,
        % otherwise print warning.
        % Start from highest mode; 
        % if freq>cutoff_highermode then also freq>cutoff_lowermode
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