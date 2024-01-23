function freq_res = calculate_resonance_frequency_unloaded(d_t, n, d, p)
% Calculates the p:th resonance frequency of an UNLOADED cavity with length
% d and transversal wavenumber d_t of mode n
    mu_0 = 1.25663706212e-6;            % vacuum permeability
    eps_0 = 8.8541878128e-12;           % vacuum permittivity
    
    freq_res = sqrt( (d_t(n,1)/pi)^2 + (p/d)^2 ) / (2*eps_0*mu_0);
    % CARE: d_t(n,1) assumes left-side geometry (whatever the leftside is).
    % Should probably update this.
end