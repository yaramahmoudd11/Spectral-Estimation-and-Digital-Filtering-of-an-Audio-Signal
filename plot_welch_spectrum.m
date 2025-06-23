function plot_welch_spectrum(signal, fs, window_type,perc,N)
    
    noverlap = N *perc;   %  overlap %
    nfft = N;           % FFT length

    % Select window
    switch lower(window_type)
        case 'rectangular'
            win = rectwin(N);
        case 'hamming'
            win = hamming(N);
        case 'hann'
            win = hann(N);
        case 'blackman'
            win = blackman(N);
        case 'flattop'
            win = flattopwin(N);
        case 'triangular'
            win = triang(N);
        otherwise
            error('Unknown window type. Choose: rectangular, hamming, hann, blackman, flattop, triangular.');
    end

    % Plot Welch power spectrum
    PSD=pwelch(signal, win, noverlap, nfft, fs);
    freq = 0:fs/nfft:fs/2;
    freq = freq/1e3;

    plot(freq, 10*log10(PSD), 'LineWidth', 1.5);
    title(['Power Spectrum using ', window_type, ' window']);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
end