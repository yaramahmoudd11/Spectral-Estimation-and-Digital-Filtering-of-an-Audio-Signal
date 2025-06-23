%% Main Script

% Read audio
[Y, fs] = audioread('music_test_fayrouz.mp3');
X = Y(3e5:4e5, 1);  % Extract approx. 3 seconds

% Time vector
n = 0:length(X)-1;
Ts = 1/fs;
t = n * Ts;

% Plot the audio segment
figure;
plot(t, X);
xlabel('Time (secs)');
ylabel('Amplitude');
title('Audio Signal (Approx. 3 seconds from middle)');

% Add sinusoidal interference
f2 = 15200;
A = 1.8;
in = A * cos(2 * pi * f2 * t).';
X_noisy = X + in;

% Play original and noisy audio
disp('Playing original audio...');
sound(X, fs);
pause(length(X)/fs + 1);

disp('Playing audio with interference...');
sound(X_noisy, fs);
pause(length(X_noisy)/fs + 1);



N = 1024;           % FFT size
N1=256;
N2=4096;
perc = 0.5;         % 50% overlap

%%  7 – Analyze Effects of  Parameters

signal = X_noisy;

%% A. Varying FFT Size (Frequency Resolution)
figure('Name', 'Effect of FFT Size');
subplot(3,1,1);
plot_welch_spectrum(signal, fs, 'hamming', perc, N1);
title('FFT Size = 256');

subplot(3,1,2);
plot_welch_spectrum(signal, fs, 'hamming', perc, N);
title('FFT Size = 1024');

subplot(3,1,3);
plot_welch_spectrum(signal, fs, 'hamming',perc, N2);
title('FFT Size = 4096');

%% B. Varying Window Type (Leakage Suppression)
figure('Name', 'Effect of Window Type');
subplot(3,2,1);
plot_welch_spectrum(X_noisy, fs, 'rectangular', perc, N);
title('Rectangular');

subplot(3,2,2);
plot_welch_spectrum(X_noisy, fs, 'hamming', perc, N);
title('Hamming');

subplot(3,2,3);
plot_welch_spectrum(X_noisy, fs, 'hann', perc, N);
title('Hann');

subplot(3,2,4);
plot_welch_spectrum(X_noisy, fs, 'blackman', perc, N);
title('Blackman');

subplot(3,2,5);
plot_welch_spectrum(X_noisy, fs, 'flattop', perc, N);
title('Flattop');


%% C. Varying Overlap Percentage (Variance Reduction)
figure('Name', 'Effect of Overlap Percentage');
subplot(3,1,1);
plot_welch_spectrum(signal, fs, 'hann', 0.0, N);
title('Overlap = 0%');

subplot(3,1,2);
plot_welch_spectrum(signal, fs, 'hann', perc, N);
title('Overlap = 50%');

subplot(3,1,3);
plot_welch_spectrum(signal, fs, 'hann', 0.75, N);
title('Overlap = 75%');

%% D. Varying Window Size (≤ FFT Size)
figure('Name', 'Effect of Window Size');

window_sizes = [256, 512, 1024];  % all ≤ N2 = 1024
for i = 1:length(window_sizes)
    subplot(length(window_sizes), 1, i);
    plot_welch_spectrum(signal, fs, 'hamming', perc, window_sizes(i));
    title(['Window Size = ', num2str(window_sizes(i))]);
end

%% E. Varying Sampling Frequency
fs_list = [fs, 22050, 8000,44000];
labels = {'Original fs', 'Downsampled to 22.05 kHz', 'Downsampled to 8 kHz','upsampled to 44 kHz'};

figure('Name', 'Effect of Sampling Frequency');
for i = 1:length(fs_list)
    % Resample signal
    fs_new = fs_list(i);
    X_rs = resample(signal, fs_new, fs);
    subplot(length(fs_list), 1, i);
    plot_welch_spectrum(X_rs, fs_new, 'hann', perc, N);
    title(labels{i});
end
%% ==========  filters uisng matlab  ==========

X_equimatlab = filter(fireqqqq,1, X_noisy);
X_firlsmatlab= filter(firls, 1, X_noisy);

X_iirelmatlab = filter(iir_elpp,1, X_noisy);
X_iirbuttmatlab= filter(iir_btt, 1, X_noisy);

disp('Playing filtered signal fir Equiripple');
sound(X_equimatlab, fs);
pause(length(X_equimatlab)/fs + 1);

disp('Playing filtered signalLS FIR');
sound(X_firlsmatlab, fs);
pause(length(X_firlsmatlab)/fs + 1);

disp('Playing Elliptic IIR filtered signal');
sound(X_iirelmatlab, fs);
pause(length(X_iirelmatlab)/fs + 1);

disp('Playing Butterworth IIR filtered signal');
sound(X_iirbuttmatlab, fs);
pause(length(X_iirbuttmatlab)/fs + 1);


%% ==========  IIR Filters ==========
% Elliptic Filter
f0 = 15200;  Q = 35; bw = f0 / Q;
Wn = [f0 - bw/2, f0 + bw/2] / (fs/2);
[b_elliptic, a_elliptic] = ellip(4, 1, 40, Wn, 'stop');

% Butterworth Filter
[b_butter, a_butter] = butter(6, Wn, 'stop');

%% ==========  FIR Filters ==========
% Frequency bands for FIR designs
f_nyq = fs/2;
bands = [0, f0-bw/2-500, f0-bw/2, f0+bw/2, f0+bw/2+500, f_nyq] / f_nyq;
desired = [1, 1, 0, 0, 1, 1]; 
weights = [1, 10, 1]; 

% 1. Equiripple FIR
N_equiripple = 100;
[b_equiripple, ~] = firpm(N_equiripple, bands, desired, weights);
filterDesigner
% 2.  LS FIR
N_hamming = 150;
[b_hamming, ~] = fircls1(N_hamming, (f0-bw/2)/f_nyq, (f0+bw/2)/f_nyq, 0.01, 0.99);

%% ========== Apply ALL Filters ==========
X_elliptic = filter(b_elliptic, a_elliptic, X_noisy);
X_butter = filter(b_butter, a_butter, X_noisy);
X_equiripple = filter(b_equiripple, 1, X_noisy);
X_hamming = filter(b_hamming, 1, X_noisy);

%% ========== Enhanced Visualization ==========
N_fft = 8192; % Higher resolution for clearer plots
f_axis = linspace(0, fs/2, N_fft/2);

% Compute all frequency responses
H_elliptic = 20*log10(abs(freqz(b_elliptic, a_elliptic, N_fft)));
H_butter = 20*log10(abs(freqz(b_butter, a_butter, N_fft)));
H_equiripple = 20*log10(abs(fft(b_equiripple, N_fft)));
H_hamming = 20*log10(abs(fft(b_hamming, N_fft)));

% Compute all signal spectra
X_orig_fft = 20*log10(abs(fft(X, N_fft)));
X_noisy_fft = 20*log10(abs(fft(X_noisy, N_fft)));
X_elliptic_fft = 20*log10(abs(fft(X_elliptic, N_fft)));
X_butter_fft = 20*log10(abs(fft(X_butter, N_fft)));
X_equiripple_fft = 20*log10(abs(fft(X_equiripple, N_fft)));
X_hamming_fft = 20*log10(abs(fft(X_hamming, N_fft)));

%% Plot 1: All Filter Responses (Zoomed)
figure('Position', [100 100 1200 800]);
subplot(2,1,1);
plot(f_axis, H_elliptic(1:N_fft/2), 'm', 'LineWidth', 2); hold on;
plot(f_axis, H_butter(1:N_fft/2), 'c', 'LineWidth', 2);
plot(f_axis, H_equiripple(1:N_fft/2), 'b', 'LineWidth', 2);
plot(f_axis, H_hamming(1:N_fft/2), 'r', 'LineWidth', 2);
xline(f0, '--k', 'Interference', 'LineWidth', 1.5);
xline([f0-bw/2, f0+bw/2], ':k', 'LineWidth', 1);
xlim([14000, 16400]); ylim([-100, 5]);
title('Filter Responses (15.2 kHz Region)', 'FontSize', 14);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Elliptic IIR', 'Butterworth IIR', 'Equiripple FIR', 'ls FIR', 'Location', 'southeast');
grid on;

%% Plot 2: Signal Spectra Comparison
subplot(2,1,2);
plot(f_axis, X_noisy_fft(1:N_fft/2), 'Color', [0.6 0.6 0.6], 'LineWidth', 1); hold on;
plot(f_axis, X_elliptic_fft(1:N_fft/2), 'm', 'LineWidth', 1.5);
plot(f_axis, X_butter_fft(1:N_fft/2), 'c', 'LineWidth', 1.5);
plot(f_axis, X_equiripple_fft(1:N_fft/2), 'b', 'LineWidth', 1.5);
plot(f_axis, X_hamming_fft(1:N_fft/2), 'r', 'LineWidth', 1.5);
xlim([0, 20000]); ylim([-20, 80]);
title('Filtered Signal Spectra', 'FontSize', 14);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Noisy', 'Elliptic', 'Butterworth', 'Equiripple', 'LS');
grid on;

%% ========== Audio Playback ==========
sound_types = {'Elliptic IIR', 'Butterworth IIR', 'Equiripple FIR', 'LS FIR'};
signals = {X_elliptic, X_butter, X_equiripple, X_hamming};

for i = 1:length(sound_types)
    fprintf('Playing %s filtered signal...\n', sound_types{i});
    sound(signals{i}, fs);
    pause(length(signals{i})/fs + 1); 
end

%% ========== Time Domain Analysis ==========
t = (0:1000)/fs; % First 1001 samples

figure('Position', [100 100 1000 800]);
for i = 1:4
    subplot(4,1,i);
    plot(t, X_noisy(1:length(t)), 'Color', [0.8 0.8 0.8]); hold on;
    plot(t, signals{i}(1:length(t)), 'LineWidth', 1.5);
    title(sprintf('%s Filtered Signal (Time Domain)', sound_types{i}));
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('Noisy', 'Filtered');
    grid on;
end

% Welch Power Spectrum of  Filter Output
figure;
pwelch(X_equiripple, blackman(1024), 512, 1024, fs);
title('Power Spectrum of Signal After  Filter');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

figure;
pwelch(X_noisy, hamming(1024), 512, 1024, fs);
title('Power Spectrum of Signal before Filter');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

figure;
pwelch(X, hamming(1024), 512, 1024, fs);
title('Power Spectrum of Signal before noise');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;



audiowrite('noise.wav', X_noisy, fs)
audiowrite('filtered_ausio.wav', X_equiripple, fs)
audiowrite('original.wav', X, fs)