function S = compile_scattering_matrix(D_1, D_2, C)
% For theory see end of section [2.2], eqs (2.19)-(2.22)
    arguments 
        D_1 % NxN matrix
        D_2 % NxN matrix
        C   % NxN matrix
    end
    S_11 = inv(conj(D_1) + conj(C)*inv(D_2)*transpose(C)) * (conj(D_1) - conj(C)*inv(D_2)*transpose(C));
    S_12 = 2 * inv(conj(D_1) + conj(C)*inv(D_2)*transpose(C)) * conj(C);
    S_21 = 2 * inv(transpose(C)*inv(conj(D_1))*conj(C) + D_2) * transpose(C);
    S_22 = inv(transpose(C)*inv(conj(D_1))*conj(C) + D_2) * (transpose(C)*inv(conj(D_1))*conj(C) - D_2);
    
    S = [S_11 S_12; S_21 S_22]; % 2Nx2N matrix
end