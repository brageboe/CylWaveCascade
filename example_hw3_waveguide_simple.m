freq = linspace(1, 10, 200).*1e9; % operating frequency 
a = 1e-2;
b = 8e-2;
R = 12e-2;
h = 6e-2; % parameter h in Fig.1
d = 18e-2;  % parameter d in Fig.1
N = 3; % number of modes (highest mode TM_{0N})

S1 = cell(length(freq),1);
S2 = S1;
S3 = S1;
S4 = S1;
S = S1;

for i=1:length(freq)
   S1{i} = scattering_matrix_coaxials(freq(i), a, b, R, 0, h, N);
   S2{i} = scattering_matrix_mixed(freq(i), a, R, 0, d-2*h, N, 1); 
   S3{i} = scattering_matrix_mixed(freq(i), a, R, 0, h, N, 2); 
   S4{i} = scattering_matrix_coaxials(freq(i), a, R, b, 0, h, N, propagator_geometry2=false);
   S{i} = S1{i} * S2{i} * S3{i} * S4{i};
   check_physical_realizability(S{i}, print_warning=false);
end




