function P = compile_propagator_matrix(N, z_1, z_2, d_z, i)
% Assumes d_z_1 is the longitudinal wavenumber in propagating medium    
    arguments
        N   % total amount of modes
        z_1 % longitudinal position 
        z_2 % longitudinal position 
        d_z % longitudinal wavenumber, size Nx1
        i   % ={1,2} specify left or right side of junction
    end
    if i ~= 1 && i~= 2
        warning("Parameter i should take values 1 or 2.")
    end
    P = zeros(N,N);
    for n=1:N
       P(n,n) = calculate_propagator_ij(n, z_1, z_2, d_z, i);
    end
end