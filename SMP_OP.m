clear, clc, close all;

fprintf('=== Starting SMP Optimization (without detuning) ===\n');

J = 695.8; % Coupling constant (Hz)
w1_max = 6250; % Maximum RF amplitude (Hz)

I = [1 0; 0 1];
sig_z = [1 0; 0 -1];
sig_y = [0 -1i; 1i 0];
sig_x = [0 1; 1 0];

H0_ij = 2*pi*J * kron(sig_z, sig_z) / 4; % Coupling Hamiltonian, used for free evolution between pulses

R1_y = cos(pi/8)*I - 1i*sin(pi/8)*sig_y;
R1_z = cos(pi/8)*I - 1i*sin(pi/8)*sig_z;
R2_y = cos(pi/8)*I + 1i*sin(pi/8)*sig_y;
R2_z = cos(pi/8)*I + 1i*sin(pi/8)*sig_z;
U1 = R1_z * R1_y;
U2 = R2_z * R2_y;
U = kron(U1,U2); % Target unitary gate

N = 6;
fprintf('Number of SMP segments: %d\n', N);

num_trials = 1;
best_fval = inf;
best_pa = [];

fprintf('\n=== Starting %d attempts ===\n', num_trials);

options = optimoptions('fminunc', 'Display', 'off', 'MaxIterations', 5000, 'MaxFunctionEvaluations', 500000, 'OptimalityTolerance', 1e-12, 'FunctionTolerance', 1e-12);

for trial = 1:num_trials

    fprintf('Attempt %d... ', trial);
    
    A_H_init = rand(1, N);
    phi_H_init = 2*pi*rand(1, N);
    A_P_init = rand(1, N);
    phi_P_init = 2*pi*rand(1, N);
    t_init = 20 + 10*rand(1, N);
    
    pa_init = [A_H_init, phi_H_init, A_P_init, phi_P_init, t_init]; % Random initial parameter vector
    
    tic;
    [pa_opt, fval] = fminunc(@(p) cost_fun(p, H0_ij, U, w1_max, N), pa_init, options); 
    tempo = toc;
    
    fprintf('fidelity = %.6f, time = %.2f s\n', 1-fval, tempo);
    
    if fval < best_fval
        best_fval = fval;
        best_pa = pa_opt;
        fprintf('  >>> New best fidelity: %.6f <<<\n', 1-best_fval);
    end
end

fprintf('\n=== Best result ===\n');
fprintf('Final fidelity: %.8f\n', 1 - best_fval);

[A_H_opt, phi_H_opt, A_P_opt, phi_P_opt, t_opt] = extr_p(best_pa, N);
t_us = t_opt * 1e6; % Convert segment durations from seconds to microseconds

generate_json_file(A_H_opt, phi_H_opt, A_P_opt, phi_P_opt, t_us, 'SMP_pulse_no_detuning.json', 1-best_fval);

fprintf('\n=== Process complete ===\n');

function Q = cost_fun(pa, H0_ij, U, w1_max, N)

    dim = 4;
    I = [1 0; 0 1];
    sig_x = [0 1; 1 0];
    sig_y = [0 -1i; 1i 0];
    
    [A_H, phi_H, A_P, phi_P, t] = extr_p(pa, N);
    
    A_H = max(0, min(1, A_H)); % Clip amplitudes to [0, 1]
    A_P = max(0, min(1, A_P)); % Clip amplitudes to [0, 1]
    
    U_smp = eye(dim); % Identity as the starting propagator

    for k = 1:N

        H1_H = 2*pi * A_H(k) * w1_max * (cos(phi_H(k))*(sig_x/2) + sin(phi_H(k))*(sig_y/2)); % RF Hamiltonian on spin H (channel 1)
        H1_P = 2*pi * A_P(k) * w1_max * (cos(phi_P(k))*(sig_x/2) + sin(phi_P(k))*(sig_y/2)); % RF Hamiltonian on spin P (channel 2)
        
        H_total = H0_ij + kron(H1_H, I) + kron(I, H1_P); % Total Hamiltonian for segment k
        
        U_seg = expm(-1i * H_total * t(k)); % Propagator for segment k
        U_smp = U_seg * U_smp; % Accumulate propagators (time-ordered)

    end    

    F = abs(trace(U' * U_smp)) / dim; % Gate fidelity
    Q = 1 - F; % Cost to minimize

end

function [A_H, phi_H, A_P, phi_P, t] = extr_p(pa, N)

    A_H = pa(1:N); % H-channel amplitudes
    phi_H = pa(N+1:2*N); % H-channel phases
    A_P = pa(2*N+1:3*N); % P-channel amplitudes
    phi_P = pa(3*N+1:4*N); % P-channel phases
    t = pa(4*N+1:5*N); % Segment durations
    
    A_H = max(0, min(1, A_H));
    A_P = max(0, min(1, A_P));
    
    t = abs(t) * 1e-6; % Convert to seconds (parameter is in microseconds)

end

function generate_json_file(A_H, phi_H, A_P, phi_P, t_us, filename, fidelity)

    if nargin < 6, filename = 'SMP_pulse.json'; end
    if nargin < 7, fidelity = ''; end
    N = length(t_us);
    
    pulse_data.description.TITLE = 'SMP';
    pulse_data.description.OWNER = 'User';
    pulse_data.description.DATE = datestr(now, 'dd-mm-yyyy');
    pulse_data.description.FIDELITY = num2str(fidelity);
    pulse_data.description.TOTALPULSEWIDTH = num2str(sum(t_us));
    pulse_data.description.SLICES = num2str(N);
    
    pulse_data.parameters.offset{1}.channel1_pulsefree_offset = 0;
    pulse_data.parameters.offset{1}.channel1_framefree_offset = 0;
    pulse_data.parameters.offset{1}.channel2_pulsefree_offset = 0;
    pulse_data.parameters.offset{1}.channel2_framefree_offset = 0;
    
    for k = 1:N

        phase_H_deg = mod(phi_H(k) * 180/pi, 360); % Convert phase to degrees in [0, 360)
        phase_P_deg = mod(phi_P(k) * 180/pi, 360); % Convert phase to degrees in [0, 360)
        
        pulse_data.pulse.channel1_pulse{k}.detuning = 0;
        pulse_data.pulse.channel1_pulse{k}.phase = phase_H_deg;
        pulse_data.pulse.channel1_pulse{k}.amplitude = A_H(k) * 100; % Amplitude as a percentage
        pulse_data.pulse.channel1_pulse{k}.width = t_us(k);
        
        pulse_data.pulse.channel2_pulse{k}.detuning = 0;
        pulse_data.pulse.channel2_pulse{k}.phase = phase_P_deg;
        pulse_data.pulse.channel2_pulse{k}.amplitude = A_P(k) * 100; % Amplitude as a percentage
        pulse_data.pulse.channel2_pulse{k}.width = t_us(k);
    end
    
    json_text = jsonencode(pulse_data, 'PrettyPrint', true);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s', json_text);
    fclose(fid);
    fprintf('File %s generated.\n', filename);
    
end