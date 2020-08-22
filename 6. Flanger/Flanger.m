%% Flanger

% Input signal
[x,fs] = audioread('noflange.wav');  
% Initial parameters
D = round(0.003*fs);
F = 2/fs;
% Internal delay buffer for x(n)
w = zeros(1, D + 1);
% Delay buffer index variable
q = 1;
a = 0.9;
[N,k] = size(x);

% Loop through input signal
for n = 1:N
    d = round((D/2)*(1 - cos(2*pi*F*n)));
    tap = q + d;
    if tap < 1
        tap = tap + (D + 1);
    end
    if tap > (D + 1)
        tap = tap - (D + 1);
    end
    y(n) = x(n) + a*w(tap);
    w(q) = x(n);
    q = q - 1;
    if q < 1
        q = D + 1;
    end
end

% Normalize y(n)
ymax = max(y);
y = y/ymax;

% Output results to speakers
sound(y,fs);

% Output results to file
audiowrite('FlangedFile.wav',y,fs);

%% Plot Results
t = 1:N;
close all
subplot(2,1,1)
plot(t/fs,x),title('no flange')
subplot(2,1,2)
plot(t/fs,y); title('flanged'), xlabel('t (sec)');
