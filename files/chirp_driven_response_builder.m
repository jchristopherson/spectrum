% A short script to create a chirp driven second-order vibrating system.

% Time Parameters
time_duration = 5.0;
sample_rate_hz = 2048;

% Frequency Sweep Parameters
init_freq_hz = 1.0;
final_freq_hz = 500.0;

% Build the chirp signal
c = (final_freq_hz - init_freq_hz) / time_duration;
chirp = @(t) sin(2 * pi * (0.5 * c * t.^2 + init_freq_hz * t));
dydt = @(t) 2 * pi * (c * t + init_freq_hz) .* cos(2 * pi * (0.5 * c * t.^2 + init_freq_hz * t));

% Model Parameters
zeta = 0.01;
fn = 250.0;
wn = 2 * pi * fn;

% Define the equations of motion:
% x" + 2 * zeta * wn * x' + wn^2 * x = 2 * zeta * wn * y' + wn^2 * y
% Which can be rewritten:
% x" = 2 * zeta * wn * (y' - x') + wn^2 * (y - x)
eom = @(t, z) [z(2); 2 * zeta * wn * (dydt(t) - z(2)) + wn^2 * (chirp(t) - z(1))];
[t, z] = ode45(eom, [0.0:1/sample_rate_hz:time_duration], [0.0 0.0]);

% Write the results to a CSV file
csvwrite('chirp.csv', [t chirp(t) z]);
