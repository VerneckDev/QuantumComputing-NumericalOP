clear, clc, close all;

I = [1,0; 0,1]; % Identity matrix
sig_y = [0,-1i; 1i,0]; % Pauli Y matrix
sig_z = [1,0; 0,-1]; % Pauli Z matrix

U1 = [1, 0, 0, 1;
0, 1, 1, 0;
0, 1, -1, 0;
1, 0, 0, -1];

U = (1/sqrt(2)) * U1;
Ut = U';

disp("Unitary matrix for the system");
disp(U);
disp("Transconjugated matrix");
disp(Ut);
disp("Ut * U = I");
disp(Ut*U);

st_up_up = kron([1;0],[1;0]);

disp("Initial state");
disp(st_up_up);

psi_theo = U * st_up_up;

disp("U Matrix applied to the initial state");
disp("Ket");
disp(psi_theo);
disp("Bra");
disp(psi_theo');

rho_theo = psi_theo*psi_theo';

disp("Theoretical density of states matrix");
disp(rho_theo);

% Experimental density of states matrix (obtained by the quantum computer) - GRAPE
rho_exp_r = [0.5799336410578559,0.004600340885815913,-0.012760127458909525,0.3646572264631376;
0.004600340885815906,0.002677036339008106,0.0016907274905727408,0.006744998574542948;
-0.012760127458909528,0.0016907274905727451,0.0019187023853373996,-0.005490992696460782;
0.3646572264631376,0.0067449985745429665,-0.005490992696460809,0.4154706202177984];

rho_exp_img = 1i*[-2.871175220871866E-17,0.025724850066156185,0.017693123741915742,-0.20630052229841542;
-0.025724850066156157,3.7947076036992655E-18,0.00150155798398627,-0.01763756050748353;
-0.0176931237419157,-0.0015015579839862668,1.5890338090490674E-18,-0.013364998253968353;
0.20630052229841545,0.017637560507483583,0.013364998253968405,2.0816681711721685E-17];

rho_exp = rho_exp_r + rho_exp_img;
disp("Rho experimental - GRAPE");
disp(rho_exp);

F = trace(sqrtm(sqrtm(rho_exp)*rho_theo*sqrtm(rho_exp)));

disp("Average fidelity of the experiment");
disp(F);

figure;
hold on;

h1 = bar3(abs(rho_exp_r),0.3);

for k = 1:length(h1)
    h1(k).FaceColor = [1 0 0];
    h1(k).XData = h1(k).XData - 0.15;
end

h2 = bar3(abs(rho_exp_img),0.3);

for k = 1:length(h2)
    h2(k).FaceColor = [0 0 1];
    h2(k).XData = h2(k).XData + 0.15;
end

xticks(1:4)
xticklabels({'00','01','10','11'})

yticks(1:4)
yticklabels({'00','01','10','11'})

view(135,30)
grid on
title('Experimental Density Matrix - GRAPE')