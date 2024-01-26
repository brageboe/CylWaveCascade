%% Used for testing together with COMSOL.
%
% A coaxial waveguide consisting of two geometries:
% Coax 1: outer radius b, coax length h_1
% Coax 2: outer radius R, coax length h_2
% Common inner conductor with radius a.
clear variables
%% User parameters
freq = linspace(1, 10, 200).*1e9; % operating frequency 
a = 1e-2;
b = 8e-2;
R = 12e-2;
h_1 = 0; % Coax 1 in Fig.1
h_2 = 6e-2; % parameter h in Fig.1
d = 18e-2;  % parameter d in Fig.1
N = 3; % number of modes (highest mode TM_{0N})

S1 = cell(length(freq),1);
S2 = S1;
S3 = S1;
S4 = S1;
S = S1;
R_mode11 = zeros(length(freq),1);
T_mode11 = zeros(length(freq),1);
R_mode12 = zeros(length(freq),1);
T_mode12 = zeros(length(freq),1);
R_mode1N = zeros(length(freq),1);
T_mode1N = zeros(length(freq),1);

for i=1:length(freq)
    S1{i} = scattering_matrix_coaxials(freq(i), a, b, R, 0, h_2, N);
    S2{i} = scattering_matrix_mixed(freq(i), a, R, 0, (d-2*h_2), N, 1); %h_1=0 cause previous S includes propagator up to current junction?
    S3{i} = scattering_matrix_mixed(freq(i), a, R, 0, h_2, N, 2); %h_1=0 cause previous S includes propagator up to current junction?
    S4{i} = scattering_matrix_coaxials(freq(i), a, R, b, 0, h_2, N, propagator_geometry2=false);
  	S{i} = S1{i} * S2{i} * S3{i} * S4{i};
   
    check_physical_realizability(S{i}, print_warning=false);
   
    R_mode11(i) = abs(S{i}(1,1))^2;
	T_mode11(i) = abs(S{i}(N+1,1))^2;
    if N > 1
        R_mode12(i) = abs(S{i}(2,1))^2;     % percentage power reflected of mode 2 when port 1 is excited by mode 1
        T_mode12(i) = abs(S{i}(N+2,1))^2;   % percentage power transmitted of mode 2 when port 1 is excited by mode 1 
        if N > 2
            R_mode1N(i) = abs(S{i}(N,1))^2;     % percentage power reflected of mode N when port 1 is excited by mode 1
            T_mode1N(i) = abs(S{i}(N+N,1))^2;   % percentage power transmitted of mode 1 when port 1 is excited by mode 1
        end
   end
end


figure(); hold on
lwidth = 1.5;
lwidth_k = 1;
fontsize = 12;

plot(freq/1e9, R_mode11, 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
plot(freq/1e9, T_mode11, 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
if N == 1
    plot(freq/1e9, R_mode11+T_mode11, 'k-', 'Linewidth', lwidth_k)
    legend("$|S11|^2$","$|S21|^2$","$|S11|^2+|S21|^2$ (TM$_{01}$)", 'Location', 'best', 'Interpreter', 'latex', 'FontSize', fontsize)
elseif N > 1
    plot(freq/1e9, R_mode12, '--', 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
    plot(freq/1e9, T_mode12, '--', 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
    if N == 2
        plot(freq/1e9, R_mode11+T_mode11+R_mode12+T_mode12, 'k-', 'Linewidth', lwidth_k)
        legend("$|S11|^2$ (TM$_{01}$)","$|S21|^2$ (TM$_{01}$)", ...
            "$|S11|^2$ (TM$_{02}$)","$|S21|^2$ (TM$_{02}$)","$|S11|^2+|S21|^2$ (all modes)", ...
            'Location', 'best', 'Interpreter', 'latex', 'FontSize', fontsize)
    end
    if N > 2
        plot(freq/1e9, R_mode1N, '-.', 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
        plot(freq/1e9, T_mode1N, '-.', 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
        if N == 3
            plot(freq/1e9, R_mode1N+T_mode1N+R_mode11+T_mode11+R_mode12+T_mode12, 'k-', 'Linewidth', lwidth_k)
            legend("$|S11|^2$ (TM$_{01}$)","$|S21|^2$ (TM$_{01}$)", ...
            "$|S11|^2$ (TM$_{02}$)","$|S21|^2$ (TM$_{02}$)", ...
            "$|S11|^2$ (TM$_{03}$)","$|S21|^2$ (TM$_{03}$)", ...
            "$|S11|^2+|S21|^2$ (all modes)", ...
            'Location', 'best', 'Interpreter', 'latex', 'FontSize', fontsize)
        end
    end
end

L = legend;
L.AutoUpdate = 'off'; 
%yline(1)
grid on
xlabel("Frequency (GHz)", 'interpreter', 'latex', 'FontSize', fontsize+2)
xlim([1.5, max(freq)/1e9])
set(gca, "FontSize", fontsize+2)
yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
hold off