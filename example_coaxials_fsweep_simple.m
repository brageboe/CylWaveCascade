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
h_1 = 6e-2;
h_2 = 6e-2;
N = 1; % number of modes (highest mode TM_{0N})

S = cell(length(freq),1);

for i=1:length(freq)
   S{i} = scattering_matrix_coaxials(freq(i), a, b, R, h_1, h_2, N);
end