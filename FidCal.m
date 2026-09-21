clear, clc, close all;

I = [1,0; 0,1]; % Identity matrix
sig_y = [0,-1i; 1i,0]; % Pauli Y matrix
sig_z = [1,0; 0,-1]; % Pauli Z matrix

R1_y = cos(pi/8)*I - 1i*sin(pi/8)*sig_y; % Rotation matrix for the first qubit around the Y axis
R1_z = cos(pi/8)*I - 1i*sin(pi/8)*sig_z;

R2_y = cos(pi/8)*I + 1i*sin(pi/8)*sig_y; % Rotation matrix for the second qubit around the Y axis
R2_z = cos(pi/8)*I + 1i*sin(pi/8)*sig_z;

U1 = R1_z * R1_y;
U2 = R2_z * R2_y;

disp("Unitary matrix for the first qubit");
disp(U1);
disp("Unitary matrix for the second qubit");
disp(U2);

U = kron(U1,U2);
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

rho_theo = psi_theo * psi_theo'; % States density matrix
rho_theo_r = real(rho_theo);
rho_theo_img = imag(rho_theo);

disp("Theoretical density of states matrix");
disp(rho_theo);

% Experimental density of states matrix (obtained by the quantum computer)
rho_exp_r = [0.6297685024513513,-0.1844219529803129,0.06489216214469543,-0.10189119809348327;
-0.18442195298031305,0.12539411835333245,0.07883097796814893,0.0081303979540904;
0.06489216214469545,0.0788309779681488,0.20824844412297305,-0.05330967033964552;
-0.10189119809348324,0.00813039795409038,-0.053309670339645496,0.036588935072344245];

rho_exp_img = 1i*[2.370269108736562E-17,-0.20590168227980493,-0.33031031898967494,0.06778088859148035;
0.20590168227980515,-1.8431436932253575E-17,0.11448421745864856,-0.04595740693298753;
0.330310318989675,-0.11448421745864847,-1.5612511283791264E-17,-0.06406345418287152;
-0.06778088859148036,0.04595740693298746,0.06406345418287146,6.7220534694101275E-18];

rho_exp = rho_exp_r + rho_exp_img;
disp("Experimental density of states matrix");
disp(rho_exp);

F = trace(sqrtm(sqrtm(rho_exp)*rho_theo*sqrtm(rho_exp)));

disp("Average fidelity of the experiment");
disp(F);

figure(1);
hold on;

% Theoretical
h1 = bar3(abs(rho_theo_r),0.3);

for k = 1:length(h1)

    h1(k).FaceColor = [1 0 0];
    h1(k).XData = h1(k).XData - 0.15;

end

h2 = bar3(abs(rho_theo_img),0.3);

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
title('Theoretical Density Matrix')

% Experimental
figure(2);
hold on;

h3 = bar3(abs(rho_exp_r),0.3);

for k = 1:length(h3)

    h3(k).FaceColor = [1 0 0];
    h3(k).XData = h3(k).XData - 0.15;

end

h4 = bar3(abs(rho_exp_img),0.3);

for k = 1:length(h4)

    h4(k).FaceColor = [0 0 1];
    h4(k).XData = h4(k).XData + 0.15;

end

xticks(1:4)
xticklabels({'00','01','10','11'})

yticks(1:4)
yticklabels({'00','01','10','11'})

view(135,30)
grid on
title('Experimental Density Matrix')