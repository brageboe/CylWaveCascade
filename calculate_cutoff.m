function freq_cutoff = calculate_cutoff(n, d_t, i)
    arguments
        n   % mode number
        d_t % transversal wavenumbers; left-side stored in d_t(:,1), right-side in d_t(:,2)
        i   % specify left side (1) or right side (2)
    end
    mu_0 = 1.25663706212e-6;            % vacuum permeability
    eps_0 = 8.8541878128e-12;           % vacuum permittivity
    freq_cutoff = d_t(n,i) / (2*pi*sqrt(mu_0*eps_0));
end