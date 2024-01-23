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
h_1 = 6e-2; % Coax 1 in Fig.1
h_2 = 6e-2; % parameter h in Fig.1
d = 18e-2;  % parameter d in Fig.1
N = 1; % number of modes (highest mode TM_{0N})

S1 = cell(length(freq),1);
S2 = S1;
S3 = S1;

for i=1:length(freq)
   S1{i} = scattering_matrix_coaxials(freq(i), a, b, R, h_1, h_2, N);
   S2{i} = scattering_matrix_mixed(freq(i), a, R, 0, (d-2*h_2), N); %h_1=0 cause previous S includes propagator up to current junction?
   % just read the fakin theory doc
end