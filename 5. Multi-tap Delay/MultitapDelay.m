%% Multi-tap Delay

% Input signal
[s,fs] = audioread('original.wav');
% Initial parameters
D1 = 0.125*fs; D2 = 0.25*fs;
b0 = 1; b1 = 1; b2 = 1;
a1 = 0.2; a2 = 0.4;
[N,k] = size(s);

% Internal delay buffer for s(n)
w = zeros(1, D1 + D2 + 1);
% Delay buffer index variable
q = 1;
% Delay buffer taps
tap1 = D1 + 1;
tap2 = D1 + D2 + 1;
% Loop through input signal
for n = 1:N
    s1 = w(tap1);
    s2 = w(tap2);
    y(n) = b0*s(n) + b1*s1 + b2*s2;
    w(q) = s(n) + a1*s1 + a2*s2;
        q = q - 1;                 % Backshift index 1
    if q < 1                    % Circulate index 1
        q = D1 + D2 + 1;
    end
    tap1 = tap1 - 1;             % Backshift tap1
    if tap1 < 1                  % Circulate tap1
        tap1 = D1 + D2 + 1;
    end
        tap2 = tap2 - 1;         % Backshift tap2
    if tap2 < 1                  % Circulate tap2
        tap2 = D1 + D2 + 1;
    end
end
% Normalize y(n)
ymax = max(y);
y = y/ymax;

% Output results to speakers
sound(y,fs);

% Output results to file
audiowrite('MultitapDelay.wav',y,fs);

%% Plot Results
t = 1:160000;
close all
subplot(2,1,1)
plot(t/fs,s),title('original')
subplot(2,1,2)
plot(t/fs,y); title('multitap'), xlabel('t (sec)');
