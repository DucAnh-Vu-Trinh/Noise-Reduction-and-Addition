%%% Part 1
%% Task 1: Preliminary Analysis
% Load original and noisy audio files
[y_real, fs_real] = audioread("IceCream.mp4");
y_real = y_real(:,1);
[y_noise, fs_noise] = audioread("NNoisy_IceCream.wav");
y_noise = y_noise(1:length(y_real));
L_real = size(y_real, 1)-1;             % Length of signal (number of samples - 1)
L_noise = size(y_noise, 1)-1;

t_real = seconds(0:1/fs_real:L_real/fs_real);
t_noise = seconds(0:1/fs_noise:L_noise/fs_noise);

% Visualize waveforms4
figure("Name", "Waveforms")
subplot(2,1,1);
plot(t_real, y_real); title("Without Noise")
subplot(2, 1, 2);
plot(t_noise, y_noise); title("With Noise")

% disp('Click on the plot to select 2 points. Press ENTER when done.');
% [selected_x, selected_y] = ginput;
% start_time = selected_x(1);
% end_time = selected_x(2);
% Convert start and end time to sample indices
% start_index = round(start_time * fs_noise) + 1;  % Add 1 to ensure inclusive indexing
% end_index = round(end_time * fs_noise);
% Extract the desired portion of the audio data
% y_noise_extracted = y_noise(start_index:end_index, :);

% Visualize spectra
[pxx1,f1]=pwelch(y_real(:,1),1000,500,1000,fs_real); % power of the DFT
[pxx2,f2]=pwelch(y_noise,1000,500,1000,fs_noise);    % power of the DFT change real to noise
figure("Name", "Spectra")
plot(f1,pxx1); hold on
plot(f2,pxx2); hold off
legend('Without Noise','With Noise');
xlabel ('Frequency (Hz)'); ylabel('Power'); 

%% Task 2: Frequency Domain Analysis
% Perform FFT
Y_real = fft(y_real(:,1));
Y_noise = fft(y_noise);
figure("Name", "FFT")
subplot(2, 1, 1);
plot(fs_real/L_real*(-L_real/2:L_real/2),fftshift(abs(Y_real)), 'r')
title("Complex Magnitude of fft Spectrum for Clean Audio"); xlabel("f (Hz)")

subplot(2, 1, 2);
plot(fs_noise/L_noise*(-L_noise/2:L_noise/2),fftshift(abs(Y_noise)))
title("Complex Magnitude of fft Spectrum for Noisy Audio"); xlabel("f (Hz)")

%% Task 3: Noise Characterization
% Isolate the noise components.
noise_fft = (fftshift(abs(Y_noise)) - fftshift(abs(Y_real))) .* (fftshift(abs(Y_noise)) > fftshift(abs(Y_real)));

% Plot noise spectrum
L_noise_new = size(Y_noise, 1)-1;
figure("Name", "Noise")
plot(fs_noise/L_noise*(-L_noise_new/2:L_noise_new/2),abs(noise_fft));
title('Noise Spectrum'); xlabel('Frequency (Hz)'); ylabel('Magnitude');
% Analyze mean and variance of noise
mean_noise = mean(abs(noise_fft));
variance_noise = var(abs(noise_fft));

disp(['Mean of Noise: ', num2str(mean_noise)]);
disp(['Variance of Noise: ', num2str(variance_noise)]);

%% Task 4: Noise Reduction
% Design and apply filter
first_stop = [885.37 2560.84];
fpass = 3000;
y_filtered1 = bandstop(lowpass(y_noise, fpass, fs_noise), first_stop, fs_noise);
% sound([y_noise; y_filtered1], fs_noise)

%% Task 5: Reconstruction and Evaluation
% Compare visually
figure("Name", "Waveforms After Filtered")
subplot(3,1,1);
plot(t_real, y_real); title("Original")
subplot(3, 1, 2);
plot(t_noise, y_noise); title("With Noise")
subplot(3, 1, 3);
plot(t_noise, y_filtered1); title("Noise Filtered")

% Calculate SNR
r1 = snr(y_real, y_filtered1 - y_real);
r2 = snr(y_real, y_noise - y_real);

disp(['SNR_improved: ', num2str(r1), ' dB']);
disp(['SNR_noise: ', num2str(r2), ' dB']);

%%% Part 2
%% Task 1: Create Custom Noisy Signals
% Uniform Random Noise Addition
lower_bound = 1200;
upper_bound = 2500;
uniform_noise = 1/5 * (rand(size(y_real)) - 0.5); % Generate noise between -0.1 and 0.1
t = linspace(0, L_real/fs_real, length(y_real)); % Time vector
noise_frequency_period = 1; % Adjust the period of the noise frequency modulation
w_rand_range = 2 * pi * rand(1) * (upper_bound - lower_bound) * t;
f_noise_uniform = sin((w_rand_range+2*pi*lower_bound*t)/noise_frequency_period);
uniform_noise = uniform_noise .* f_noise_uniform'; % Apply frequency modulation

% Add uniform random noise to the original signal
y_noisy_signal_1 = y_real + uniform_noise;
figure("Name", "Uniform Noise")
subplot(2,1,1);
plot(t_real, y_real); title("Without Noise")
subplot(2, 1, 2);
plot(t_real, y_noisy_signal_1); title("With Noise")
audiowrite('noisy_signal_1.wav', y_noisy_signal_1, fs_real);

% Alternative Noise Addition Method
varn = 0.001;
avg = 0.01;
y_noisy_signal_2 = y_real + sqrt(varn)*randn(size(y_real))+avg; % Gaussian noise with mean avg and variance varn
figure("Name", "Gausian Noise")
subplot(2,1,1);
plot(t_real, y_real); title("Without Noise")
subplot(2, 1, 2);
plot(t_real, y_noisy_signal_2); title("With Noise")
audiowrite('noisy_signal_2.wav', y_noisy_signal_2, fs_real);

%% Task 2: Noise Reduction and Signal Reconstruction
% Normal Noise for a range of frequency
fn = fs_real/2;                                              % Nyquist Frequency (Hz)
Wp = lower_bound/fn;                                           % Passband Frequency (Normalised)
Ws = upper_bound/fn;                                           % Stopband Frequency (Normalised)
Rp =   2;                                               % Passband Ripple (dB)
Rs = 80;                                               % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
[z,p,k] = cheby2(n,Rs,[Wp Ws] ,'stop');                        % Filter Design
[soslp,glp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability                               % Filter Bode Plot
y_noisy_signal_1_filtered = filtfilt(soslp, glp, y_noisy_signal_1);
y_noisy_signal_1_filtered = smoothdata(y_noisy_signal_1_filtered);
fpass = 20000;
y_noisy_signal_1_filtered = lowpass(y_noisy_signal_1_filtered, fpass, fs_real);
% sound(y_noisy_signal_1_filtered, fs_real)
% Calculate SNR
r1 = snr(y_real, y_noisy_signal_1_filtered - y_real);
r2 = snr(y_real, y_noisy_signal_1 - y_real);

disp(['SNR_improved_Normal: ', num2str(r1), ' dB']);
disp(['SNR_noise_Normal: ', num2str(r2), ' dB']);
% Compare visually
figure("Name", "Normal Noise Waveforms After Filtered")
subplot(3,1,1);
plot(t_real, y_real); title("Original")
subplot(3, 1, 2);
plot(t_real, y_noisy_signal_1); title("With Noise")
subplot(3, 1, 3);
plot(t_real, y_noisy_signal_1_filtered); title("Noise Filtered")

% Gausian Noise                                     % Sampling Frequency (Hz)
y_noisy_signal_2_filtered = wdenoise(y_noisy_signal_2, "DenoisingMethod", 'Bayes', 'ThresholdRule', 'Soft', 'NoiseEstimate', 'LevelIndependent', 12, 'Wavelet', 'sym8');
% Compare visually
figure("Name", "Gausian Waveforms After Filtered")
subplot(3,1,1);
plot(t_real, y_real); title("Original")
subplot(3, 1, 2);
plot(t_real, y_noisy_signal_2); title("With Noise")
subplot(3, 1, 3);
plot(t_real, y_noisy_signal_2_filtered); title("Noise Filtered")
% sound(y_noisy_signal_2_filtered, fs_real)
% Calculate SNR
r1 = snr(y_real, y_noisy_signal_2_filtered - y_real);
r2 = snr(y_real, y_noisy_signal_2 - y_real);
disp(['SNR_improved_Gausian: ', num2str(r1), ' dB']);
disp(['SNR_noise_Gausian: ', num2str(r2), ' dB']);