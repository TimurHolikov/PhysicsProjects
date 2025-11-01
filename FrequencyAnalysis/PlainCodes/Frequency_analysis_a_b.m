% Assignment 2. Task 1. Parts (a) - (b)
% Names: Timur Holikov, Lea Ebenschweiger

% Warning: The execution time of these two subtasks can exceed 15 minutes (at least this is the case on my computer)

[y, SR] = audioread('single_tone.ogg');
y_L = y(:,1); %import left channel
n = length(y_L);

% Find the next power of 2
n_ext = 2^nextpow2(n);

% Extend the signal with zeros
y_ext = [y_L; zeros(n_ext-n, 1)];

%% DFT
% tic; Measure the execution time
Y_DFT = zeros(n_ext, 1);
Y_DFT = DFT(n_ext, Y_DFT, y_ext);


%% FFT
Y_FFT = fft(y_ext);

% A negative exponent sign is used in Matlab to define FFT, so one need to
% take the complex conjugate for comparison
Y_FFT = conj(Y_FFT);

%% Comparison
A1 = abs(Y_FFT);
A2 = abs(Y_DFT);

figure
subplot(2,1,1)
plot(A1)
title('Amplitude spectrum FFT')
subplot(2,1,2)
plot(A2)
title('Amplitude spectrum DFT')
% t = toc; approx. 18 Minutes

%% (b):
m_values = linspace(128, n/5, 60); % segment lengths
time_values = zeros(size(m_values)); % computation times

for i = 1:length(m_values)
    % m = 2^nextpow2(m_values(i));
    m = round(m_values(i));
    y_segment = y_ext(1:m); % select a segment

    tic; % start timer
    Y_DFT_segment = zeros(m, 1);
    Y_DFT_segment = DFT(m, Y_DFT_segment, y_segment);

    time_values(i) = toc; % stop timer and record time
end

figure
semilogx(m_values, time_values, 'o-')
xscale log

xlabel('Segment length m')
ylabel('Computation time (s)')
title('Computation time of DFT as a function of m')

%% Functions
function Y_DFT = DFT(n, Y_DFT, y)
    for k = 0:n-1
        for j = 0:n-1
            Y_DFT(k+1) = Y_DFT(k+1) + y(j+1)*exp(1i*2*pi*k*j/n);
        end
    end
end