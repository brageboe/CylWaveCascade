% A coaxial waveguide consisting of two geometries:
% Coax 1: outer radius b, coax length h_1
% Coax 2: outer radius R, coax length h_2
% Common inner conductor with radius a.
clear variables
%% User parameters
freq = linspace(1, 10, 100).*1e9; % operating frequency 
a = 1e-2;
b = 8e-2;
R = 12e-2;
h_1 = 6e-2;
h_2 = 6e-2;
N = 2; % number of modes (highest mode TM_{0N})

S = cell(length(freq),1);
R_port1_mode11 = zeros(length(freq),1);
T_port1_mode11 = zeros(length(freq),1);
R_port1_mode12 = zeros(length(freq),1);
T_port1_mode12 = zeros(length(freq),1);
R_port1_mode1N = zeros(length(freq),1);
T_port1_mode1N = zeros(length(freq),1);

%% Main calculation
for i=1:length(freq)
   S{i} = scattering_matrix_coaxials(freq(i), a, b, R, h_1, h_2, N, propagator_geometry1=false, print_cutoff=false);
   %S_b{i} = scattering_matrix_mixed(freq(i), a, R, h_2,
   
   
   R_port1_mode11(i) = abs(S{i}(1,1))^2;
   T_port1_mode11(i) = abs(S{i}(N+1,1))^2;
   if N > 1
        R_port1_mode12(i) = abs(S{i}(2,1))^2;     % percentage power reflected of mode 2 when port 1 is excited by mode 1
        T_port1_mode12(i) = abs(S{i}(N+2,1))^2;   % percentage power transmitted of mode 2 when port 1 is excited by mode 1 
        if N > 2
            R_port1_mode1N(i) = abs(S{i}(N,1))^2;     % percentage power reflected of mode N when port 1 is excited by mode 1
            T_port1_mode1N(i) = abs(S{i}(N+N,1))^2;   % percentage power transmitted of mode 1 when port 1 is excited by mode 1
        end
   end
end

%% Figure
figure(); hold on
lwidth = 1.25;

plot(freq/1e9, R_port1_mode11, 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
plot(freq/1e9, T_port1_mode11, 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
if N == 1
    plot(freq/1e9, R_port1_mode11+T_port1_mode11, 'k-', 'Linewidth', 0.5)
    legend("$|S11|^2$","$|S21|^2$","$|S11|^2+|S21|^2$ (TM$_{01}$)", 'Location', 'best', 'Interpreter', 'latex')
elseif N > 1
    plot(freq/1e9, R_port1_mode12, '--', 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
    plot(freq/1e9, T_port1_mode12, '--', 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
    if N == 2
        plot(freq/1e9, R_port1_mode11+T_port1_mode11+R_port1_mode12+T_port1_mode12, 'k-', 'Linewidth', 0.5)
        legend("$|S11|^2$ (TM$_{01}$)","$|S21|^2$ (TM$_{01}$)", ...
            "$|S11|^2$ (TM$_{02}$)","$|S21|^2$ (TM$_{02}$)","$|S11|^2+|S21|^2$ (all modes)", ...
            'Location', 'best', 'Interpreter', 'latex')
    end
    if N > 2
        plot(freq/1e9, R_port1_mode1N, '-.', 'Linewidth', lwidth, 'Color', [0 0.4470 0.7410])
        plot(freq/1e9, T_port1_mode1N, '-.', 'Linewidth', lwidth, 'Color', [0.8500 0.3250 0.0980])
        if N == 3
            plot(freq/1e9, R_port1_mode1N+T_port1_mode1N+R_port1_mode11+T_port1_mode11+R_port1_mode12+T_port1_mode12, 'k-', 'Linewidth', 0.5)
            legend("$|S11|^2$ (TM$_{01}$)","$|S21|^2$ (TM$_{01}$)", ...
            "$|S11|^2$ (TM$_{02}$)","$|S21|^2$ (TM$_{02}$)", ...
            "$|S11|^2$ (TM$_{03}$)","$|S21|^2$ (TM$_{03}$)", ...
            "$|S11|^2+|S21|^2$ (all modes)", ...
            'Location', 'best', 'Interpreter', 'latex')
        end
    end
end
L = legend;
L.AutoUpdate = 'off'; 
%yline(1)
grid on
xlabel("Frequency (GHz)", 'interpreter', 'latex')
hold off

