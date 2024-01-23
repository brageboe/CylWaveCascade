function D = compile_diagonal_matrix(d_z, i)
% compile diagonal matrix D_i, where i={1,2} specifies which waveguide geometry
    arguments
        d_z % longitudinal wavenumbers; left-side stored in d_z(:,1), right-side in d_z(:,2)
        i   % 1 or 2, i.e. left-side or right-side of junction
    end
    D = diag(abs(d_z(:,i))./conj(d_z(:,i)));
end