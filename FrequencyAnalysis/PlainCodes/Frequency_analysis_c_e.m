close all; clear; clc

[y, SR] = audioread('single_tone.ogg');
y_L = y(:,1); %import left channel
n = length(y_L);

% Find the next power of 2
n_ext = 2^nextpow2(n);

% Extend the signal with zeros
y_ext = [y_L; zeros(n_ext-n, 1)];

%% FFT
Y_FFT = fft(y_ext);

%% (c)
S = 1/(n_ext^2)*(abs(Y_FFT)).^2;

f = (1:n_ext)*SR/n_ext; % frequencies in Hz

figure
% plot(f, S)
semilogy(f, S)
xlabel('Frequency (Hz)')
%xlim([0 SR/2])  % Limit the x-axis to [0, Nyquist frequency]
xlim([-0.3*10^(4) 4.7*10^4])
ylabel('Power Spectral Density')
title('Power Spectral Density vs Frequency')

grid on

%% (d)

% Find the fundamental frequency
[~, idx] = max(S);
fundamental_freq = f(idx);

%% (e)
% Consider only up to Nyquist frequency. Isn't that necessary, but allows
% to make a better approximation
nyquistIdx = SR/2;
S = S(1:nyquistIdx);
f = f(1:nyquistIdx);

% Identify top frequencies
[pks, locs] = findpeaks(S);
[~, sortIndex] = sort(pks, 'descend');
topFreqs = f(locs(sortIndex(1:6)));

% Generate signal
approxSignal = zeros(size(y_ext));
for freq = topFreqs
    approxSignal = approxSignal + sin(2*pi*freq*(0:length(y_ext)-1)/SR)';
end

% Normalize to [-1, 1]. Prevents ddistortion
approxSignal = approxSignal / max(abs(approxSignal));

audiowrite('approximation.ogg', approxSignal, SR);
