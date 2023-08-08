% Danesh Abdollahi - 9723053 - Optical Communication Systems - Fall 2022
%%
clc; clear; close all;

%% Constants
p0 = 5e-3;  % watt
gamma = 1.2;  % 1/km
beta2 = -21.6e-24;  % s^2/km
l = 500;  % Km
t0_2 = abs(beta2) / (gamma * p0);  % s^2
t0 = sqrt(t0_2);  % s
Ts = 1e-12; 
t = -1e-9 : Ts : 1e-9 - Ts ;
fs = 1 / Ts ;

dz = 0.1;  % km
num_of_steps = round(l / dz);

C = 0 ; % chirp

%% Part 1
input_pulse = sqrt(p0) .* sech(t/t0);
df = fs / length( t ) ;
freq = -fs/2 : df : (fs/2) - df ;

output_symmetric = fftshift(fft(input_pulse)) ;
for i = 1 : num_of_steps 
    output_symmetric = D_operator(output_symmetric, beta2, dz/2, freq);
    output_symmetric = ifft(fftshift(output_symmetric));
    output_symmetric = N_operator(output_symmetric, gamma, dz) ;
    output_symmetric = fftshift(fft(output_symmetric));
    output_symmetric = D_operator(output_symmetric, beta2, dz/2, freq);
end
output_symmetric = ifft(fftshift(output_symmetric));

% Asymetric
output_asymmetric = fftshift(fft(input_pulse)) ;
for i = 1 : num_of_steps 
    output_asymmetric = D_operator(output_asymmetric, beta2, dz, freq);
    output_asymmetric = ifft(fftshift(output_asymmetric));
    output_asymmetric = N_operator(output_asymmetric, gamma, dz) ;
    output_asymmetric = fftshift(fft(output_asymmetric));
end
output_asymmetric = ifft(fftshift(output_asymmetric));

figure();
subplot(2,1,1);
plot(t * 1e9, input_pulse, LineWidth=1.3, Color="red", LineStyle="-");
hold on;
plot(t * 1e9, abs(output_symmetric), LineWidth=1.3, Color="blue", LineStyle="-.");
grid minor ;
title("Pulse Evolution (Soliton Pulse) (Symmetric)");
xlabel("Time (ns)");
ylabel("Magnitude");
legend("Input Pulse", "Output Pulse");

subplot(2,1,2);
plot(t * 1e9, input_pulse, LineWidth=1.3, Color="red", LineStyle="-");
hold on;
plot(t * 1e9, abs(output_asymmetric), LineWidth=1.3, Color="blue", LineStyle="-.");
grid minor;
title("Pulse Evolution (Soliton Pulse) (Asymmetric)");
xlabel("Time (ns)");
ylabel("Magnitude");
legend("Input Pulse", "Output Pulse");

%% Part 2
input_pulse = sqrt(p0) .* exp( -((1 + 1i * C) / (2*t0_2)) * t.^2 );
df = fs / length( t ) ;
freq = -fs/2 : df : (fs/2) - df ;

% Symmetric
output_symmetric = fftshift(fft(input_pulse)) ;
for i = 1 : num_of_steps 
    output_symmetric = D_operator(output_symmetric, beta2, dz/2, freq);
    output_symmetric = ifft(fftshift(output_symmetric));
    output_symmetric = N_operator(output_symmetric, gamma, dz) ;
    output_symmetric = fftshift(fft(output_symmetric));
    output_symmetric = D_operator(output_symmetric, beta2, dz/2, freq);
end
output_symmetric = ifft(fftshift(output_symmetric));

% Asymetric
output_asymmetric = fftshift(fft(input_pulse)) ;
for i = 1 : num_of_steps 
    output_asymmetric = D_operator(output_asymmetric, beta2, dz, freq);
    output_asymmetric = ifft(fftshift(output_asymmetric));
    output_asymmetric = N_operator(output_asymmetric, gamma, dz) ;
    output_asymmetric = fftshift(fft(output_asymmetric));
end
output_asymmetric = ifft(fftshift(output_asymmetric));

figure();
subplot(2,1,1);
plot(t * 1e9, input_pulse, LineWidth=1.3, Color="red", LineStyle="-");
hold on;
plot(t * 1e9, abs(output_symmetric), LineWidth=1.3, Color="blue", LineStyle="-.");
grid minor ;
title("Pulse Evolution (Gaussian Pulse) (Symmetric)");
xlabel("Time (ns)");
ylabel("Magnitude");
legend("Input Pulse", "Output Pulse");

subplot(2,1,2);
plot(t * 1e9, input_pulse, LineWidth=1.3, Color="red",LineStyle="-");
hold on;
plot(t * 1e9, abs(output_asymmetric), LineWidth=1.3, Color="blue", LineStyle="-.");
grid minor ;
title("Pulse Evolution (Gaussian Pulse) (Asymmetric)");
xlabel("Time (ns)");
ylabel("Magnitude");
legend("Input Pulse", "Output Pulse");


% Functions
function output = D_operator(input, beta2, dz, freq)
    D_freq = exp(1i * beta2 * dz * ((2*pi*freq).^2) / 2 );
    output = input .* D_freq ;
end

function output = N_operator(input, gamma, dz)
    N = exp( 1i * 8 * gamma * dz * (abs(input).^2) / 9 );
    output = input .* N ;
end
