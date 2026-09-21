clear, clc, close all;

fprintf('=== Starting GRAPE Optimization with particleswarm (multiple trials) ===\n');

J = 695.8; % Coupling constant (Hz)
w1_max = 6250; % Maximum RF amplitude (Hz)

I = eye(2);
sig_z = [1 0; 0 -1];
sig_y = [0 -1i; 1i 0];
sig_x = [0 1; 1 0];

H0_ij = 2*pi*J * kron(sig_z, sig_z) / 4; % Coupling Hamiltonian, used for free evolution between pulses

st_up_up = kron([1;0], [1;0]);
st_target = (kron([1;0], [1;0]) + kron([0;1], [0;1])) / sqrt(2); % Initial and target states (Bell)

T_total = 1000e-6;
N = 300;
dt = T_total / N; % Duration of each segment (time step)
fprintf('Total duration: %.2f μs, N segments: %d, dt: %.2f μs\n', T_total*1e6, N, dt*1e6);

nvars = 4*N; % Total number of variables: amplitudes and phases of the 2 channels, for each segment

lb = [zeros(1, N), -inf*ones(1, N), zeros(1, N), -inf*ones(1, N)]; % Lower bounds: [A_H, phi_H, A_P, phi_P]
ub = [ones(1, N),   inf*ones(1, N), ones(1, N),   inf*ones(1, N)]; % Upper bounds: [A_H, phi_H, A_P, phi_P]

options_ps = optimoptions('particleswarm', 'SwarmSize', min(100, 10*nvars), 'HybridFcn', @fmincon, 'Display', 'off', 'MaxIterations', 2000, 'FunctionTolerance', 1e-6, 'UseParallel', false);

num_trials = 5;
best_fval = inf;
best_pa = [];

fprintf('\n=== Starting %d attempts ===\n', num_trials);

for trial = 1:num_trials

    fprintf('\n--- Attempt %d ---\n', trial);
    
    tic;
    [pa_opt, fval] = particleswarm(@(p) cost_fun_grape(p, H0_ij, st_up_up, st_target, w1_max, N, dt), nvars, lb, ub, options_ps);
    tempo = toc;
    
    fprintf('Fidelity = %.6f, time = %.2f s\n', 1-fval, tempo);
    
    if fval < best_fval

        best_fval = fval;
        best_pa = pa_opt;
        fprintf('  >>> New best: %.6f <<<\n', 1-best_fval);

    end

end

fprintf('\n=== Best result ===\n');
fprintf('Final fidelity: %.8f\n', 1 - best_fval);

A_H_opt = best_pa(1:N); % Optimal amplitudes of channel H
phi_H_opt = best_pa(N+1:2*N); % Optimal phases of channel H
A_P_opt = best_pa(2*N+1:3*N); % Optimal amplitudes of channel P
phi_P_opt = best_pa(3*N+1:4*N); % Optimal phases of channel P

generate_json_grape(A_H_opt, phi_H_opt, A_P_opt, phi_P_opt, dt, 'GRAPE_Bell_pulse.json', 1-best_fval);

fprintf('\n=== Process complete ===\n');

function Q = cost_fun_grape(pa, H0_ij, st_up_up, st_target, w1_max, N, dt)

    dim = 4;
    I = eye(2);
    sig_x = [0 1; 1 0];
    sig_y = [0 -1i; 1i 0];
    
    A_H = pa(1:N);
    phi_H = pa(N+1:2*N);
    A_P = pa(2*N+1:3*N);
    phi_P = pa(3*N+1:4*N);
    
    A_H = max(0, min(1, A_H)); % Clip amplitudes to [0, 1]
    A_P = max(0, min(1, A_P)); % Clip amplitudes to [0, 1]
    
    U_total = eye(dim); % Identity as the starting propagator
    
    for k = 1:N

        H1_H = 2*pi * A_H(k) * w1_max * (cos(phi_H(k))*(sig_x/2) + sin(phi_H(k))*(sig_y/2)); % RF Hamiltonian on spin H (channel 1)
        H1_P = 2*pi * A_P(k) * w1_max * (cos(phi_P(k))*(sig_x/2) + sin(phi_P(k))*(sig_y/2)); % RF Hamiltonian on spin P (channel 2)
        H = H0_ij + kron(H1_H, I) + kron(I, H1_P); % Total Hamiltonian for segment k
        U_seg = expm(-1i * H * dt); % Propagator for segment k
        U_total = U_seg * U_total; % Accumulate propagators (time-ordered)

    end    

    st_final = U_total * st_up_up; % Final state obtained from the initial state
    F = abs(st_target' * st_final)^2; % Fidelity relative to the target state
    Q = 1 - F; % Cost to minimize

end

function generate_json_grape(A_H, phi_H, A_P, phi_P, dt, filename, fidelity)

    if nargin < 6, filename = 'GRAPE_pulse.json'; end
    if nargin < 7, fidelity = ''; end
    N = length(A_H);
    t_us = dt * 1e6; % Duration of each segment in microseconds

    pulse_data.description.TITLE = 'GRAPE';
    pulse_data.description.OWNER = 'User';
    pulse_data.description.DATE = datestr(now, 'dd-mm-yyyy');
    pulse_data.description.FIDELITY = num2str(fidelity);
    pulse_data.description.TOTALPULSEWIDTH = num2str(N * t_us);
    pulse_data.description.SLICES = num2str(N);
    
    pulse_data.parameters.offset{1}.channel1_pulsefree_offset = 0;
    pulse_data.parameters.offset{1}.channel1_framefree_offset = 0;
    pulse_data.parameters.offset{1}.channel2_pulsefree_offset = 0;
    pulse_data.parameters.offset{1}.channel2_framefree_offset = 0;
    
    for k = 1:N

        phase_H_deg = mod(phi_H(k) * 180/pi, 360); % Convert phase to degrees, in [0, 360)
        phase_P_deg = mod(phi_P(k) * 180/pi, 360); % Convert phase to degrees, in [0, 360)
        
        pulse_data.pulse.channel1_pulse{k}.detuning = 0;
        pulse_data.pulse.channel1_pulse{k}.phase = phase_H_deg;
        pulse_data.pulse.channel1_pulse{k}.amplitude = A_H(k) * 100; % Amplitude as a percentage
        pulse_data.pulse.channel1_pulse{k}.width = t_us;
        
        pulse_data.pulse.channel2_pulse{k}.detuning = 0;
        pulse_data.pulse.channel2_pulse{k}.phase = phase_P_deg;
        pulse_data.pulse.channel2_pulse{k}.amplitude = A_P(k) * 100; % Amplitude as a percentage
        pulse_data.pulse.channel2_pulse{k}.width = t_us;

    end   

    json_text = jsonencode(pulse_data, 'PrettyPrint', true);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_text);
    fclose(fid);
    fprintf('File %s generated.\n', filename);
    
end