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

rho_theo = psi_theo*psi_theo';

disp("Theoretical density of states matrix");
disp(rho_theo);

% Experimental density of states matrix (obtained by the quantum computer) - SMP
rho_exp_r = [0.6881552677287189,-0.20233054423729602,0.15986369097654937,-0.08473038548097143;
-0.20233054423729593,0.13206403778908254,0.0261080742726483,0.015039340693675816;
0.15986369097654932,0.02610807427264826,0.16556699449086698,-0.030137440085941884;
-0.08473038548097145,0.015039340693675832,-0.030137440085941887,0.014213699991331664];

rho_exp_img = 1i*[3.3203167759274493E-17,-0.22299367005419238,-0.21350197306502117,0.02890729139312426;
0.22299367005419246,-5.3017486234541167E-17,0.11253499455700894,-0.0366992256042837;
0.21350197306502117,-0.11253499455700892,1.0408340855860843E-17,-0.03211928309006394;
-0.028907291393124264,0.0366992256042837,0.032119283090063955,-6.505213034913027E-19];

rho_exp = rho_exp_r + rho_exp_img;
disp("Experimental density of states matrix");
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
title('Experimental Density Matrix - SMP')