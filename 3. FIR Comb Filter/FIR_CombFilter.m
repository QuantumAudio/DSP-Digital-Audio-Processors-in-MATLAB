%% FIR Comb Filter

% Input signal
[s,fs] = audioread('original.wav');
% Initial parameters
D = 0.25*fs;
a = 0.45;
[N,k] = size(s);
% Internal delay buffer for s(n)
w = zeros(1, 3*D + 1);
% Delay buffer index variable
q = 1;
% Delay buffer taps
tap1 = D + 1;
tap2 = 2*D + 1;
tap3 = 3*D + 1;
% Loop through input signal
for n = 1:N
    % Read input into w.
    w(q) = s(n);
    % y(n) = s(n) + as(n-D) + a^2s(n-2D) + a^3s(n-3D)
    y(n) = w(q) + a*w(tap1) + a^2*w(tap2) + a^3*w(tap3);    
    q = q - 1;                  % Backshift index
    if q < 1                    % Circulate index
        q = 3*D + 1;
    end
    tap1 = tap1 - 1;              % Backshift tap1
    if tap1 < 1                  % Circulate tap1
        tap1 = 3*D + 1;
    end
        tap2 = tap2 - 1;         % Backshift tap2
    if tap2 < 1                  % Circulate tap2
        tap2 = 3*D + 1;
    end
        tap3 = tap3 - 1;         % Backshift tap3
    if tap3 < 1                  % Circulate tap1
        tap3 = 3*D + 1;
    end
end
% Normalize y(n)
ymax = max(y);
y = y/ymax;

% Output Results to Speakers
sound(y,fs);

%% Plot Results
close all
t = 1:160000;
subplot(2,1,1)
plot(t/fs,s); title('original')
subplot(2,1,2)
plot(t/fs,y); title('FIR Comb'), xlabel('t (sec)')

%% Output Results to File
audiowrite('FIRComb.wav',y,fs);