%% Echo Remover

% Loading the file with the echo into a vector x.
[x,fs] = audioread('echo.wav'); % Sampling rate is 16kHz.
[s,fs] = audioread('original.wav');

% Compute the autocorrelation and find Rxx(0)
Rxx = xcorr(x);
[Rxx0, Rxx0_k] = max(Rxx);

%% Plot autocorrelation
close all
t = -159999:159999;
plot(t, Rxx)
grid on, xlabel('n (samples)'), ylabel('R_{xx}'), title('Autocorrelation of signal x(n)')

%% Finding Rxx(D) and D
range = Rxx((Rxx0_k + 8000):(Rxx0_k + 15000));
[RxxD,RxxD_k] = max(range);

D = 8000 + RxxD_k - 1;

%% Solving for a
syms A

quadratic = (A^2)*RxxD - A*Rxx0 + RxxD == 0;

a = double(solve(quadratic, A))

% Taking the first root
a = a(1);

%% Implement the echo-removing algorithm
D = 11000;
clear y;
[N,k] = size(x);
% Internal delay buffer for y(n)
w = zeros(1, D + 1);
% Delay buffer index variable
q = 1;
% Delay buffer tap
tap = D + 1;
% Loop through input signal
for n = 1:N
    % Read input sample directly into arithmetic unit to save space.
    w(q) = x(n) - a*w(tap);     % y(n) = x(n) - ay(n-D)
    y(n) = w(q);                % output y(n)
    q = q - 1;                  % Backshift index
    if q < 1                    % Circulate index
        q = D + 1;
    end
    tap = tap - 1;              % Backshift tap
    if tap < 1                  % Circulate tap
        tap = D + 1;
    end
end
%% Plot Results
close all
t = 1:160000;
subplot(3,1,1);
plot(t/fs,x); title('recorded, single echo');
subplot(3,1,2); 
plot(t/fs,y,'Color','r'); title('recovered');
subplot(3,1,3); 
plot(t/fs,s); title('original'); xlabel('t (sec)')

%% Output Results
audiowrite('EchoRemoved.wav',y,fs);