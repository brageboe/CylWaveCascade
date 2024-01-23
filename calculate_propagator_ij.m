function P_ij = calculate_propagator_ij(n, z_1, z_2, d_z, i)
% Calculates the propagation coefficient of mode ij between two points
% between the cross-sectional waveguide junctions L and R in the scattering
% system.
% d_zi = d_z0n depend on the waveguide geometry between the points
% {z_1, z_2}
% Modes i=TM_{0n} and j=TM_{0q}
% Wave travelling from point z_1 to z_2, where z_2 is assumed to be the
% junction point between geometry 1 and geometry 2. I.e. the wave
% propagates in geometry 1 with wavenumber d_z0n.
    arguments
        n   % mode
        z_1 % longitudinal position 
        z_2 % longitudinal position 
        d_z % longitudinal wavenumber, size Nx2
        i   % ={1,2} specify left or right side of junction
    end
    P_ij = exp(-1i*d_z(n,i)*(z_1 - z_2));
end