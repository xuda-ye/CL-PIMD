#!/bin/bash

var="A1" # choice of parameter

cat > "${var}.m" <<EOF
clc
file_name = mfilename('fullname');
% ---- method parameters ---- %
random_seed = str2double(file_name(2));
char = file_name(1);
switch char
    case 'A'; N = 5;
    case 'B'; N = 9;
    case 'C'; N = 17;
    case 'D'; N = 33;
    case 'E'; N = 65;
end
D = 128; % number of grid points
dt = 1/16; % time step
T = 1e6; % simulation time
T_warm = 1e3; % warm up time
ITER = round(T/dt); 
ITER_warm = round(T_warm/dt);
% ---- saved data ---- %
xi_data = zeros(5,ITER);
% warm up
fprintf("N: %d\n", N)
fprintf("D: %d\n", D)
pimd = CL_PIMD(N, D);
[xi, eta] = pimd.initialize();
for iter = 1:ITER_warm
    [xi, eta] = pimd.BCOCB(xi, eta, dt);
end
% ---- PIMD simulation ---- %
tic
for iter = 1:ITER
    [xi, eta] = pimd.BCOCB(xi, eta, dt);
    xi_data(:,iter) = xi(1:5);
    if mod(iter, ITER/100) == 0
        fprintf("%d\n", iter/(ITER/100))
    end
end
duration = toc;
fprintf("N: %d\n", N)
fprintf("D: %d\n", D)
fprintf("dt: %g\n", dt)
fprintf("time: %g\n", duration)
save(append(file_name,'_data.mat'));
EOF

cat > "${var}.sh" <<EOF
#!/bin/bash
#SBATCH -o ${var}.out
#SBATCH -J ${var}
#SBATCH --partition=C032M0128G   
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load matlab/R2022b
matlab -nodesktop -nosplash -nodisplay -r ${var}
EOF

sbatch "${var}.sh"