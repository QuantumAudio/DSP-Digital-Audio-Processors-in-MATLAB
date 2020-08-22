%% Multi-Delay

% Input signal
[s,fs] = audioread('original.wav');
% Initial parameters
D1 = 0.25*fs; D2 = 0.125*fs;
b0 = 1; b1 = 1; b2 = 1;
a1 = 0.2; a2 = 0.4;
[N,k] = size(s);
% Internal delay buffer 1 for s(n)
w1 = zeros(1, D1 + 1);
% Internal delay buffer 2 for w1(n)
w2 = zeros(1, D2 + 1);
% Delay buffer 1 index variable
q1 = 1;
% Delay buffer 2 index variable
q2 = 1;
% Delay buffer taps
tap1 = D1 + 1;
tap2 = D2 + 1;
% Loop through input signal
for n = 1:N
    s1 = w1(tap1);
    s2 = w2(tap2);
    y(n) = b0*s(n) + b1*s1 + b2*s2;
    w2(q2) = s1 + a2*s2;
    w1(q1) = s(n) + a1*s1;
    q1 = q1 - 1;                 % Backshift index 1
    if q1 < 1                    % Circulate index 1
        q1 = D1 + 1;
    end
    tap1 = tap1 - 1;             % Backshift tap1
    if tap1 < 1                  % Circulate tap1
        tap1 = D1 + 1;
    end
        q2 = q2 - 1;             % Backshift index 2
    if q2 < 1                    % Circulate index 2
        q2 = D2 + 1;
    end
        tap2 = tap2 - 1;         % Backshift tap2
    if tap2 < 1                  % Circulate tap2
        tap2 = D2 + 1;
    end
end

% Normalize y(n)
ymax = max(y);
y = y/ymax;
% Output results to speakers
sound(y,fs);
% Output results to file
audiowrite('MultiDelay.wav',y,fs);