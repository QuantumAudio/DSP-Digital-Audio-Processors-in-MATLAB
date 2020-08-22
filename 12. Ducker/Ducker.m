%% Ducker

% Input signals
[xs, fs] = audioread('speech.wav');  % speech signal
[xm, fs] = audioread('music.wav');   % music signal
xms = xm + xs;                       % speech plus music
t = (1:length(xs))/fs;               % time in sec

% Audition input signal
% sound(xms,fs)

%% Plot Input Signals
subplot(3,1,1)
plot(t,xs),title('speech, x_s(t)'),grid on
subplot(3,1,2)
plot(t,xm),title('music, x_m(t)'),grid on
subplot(3,1,3)
plot(t,xms),xlabel('t(sec)'),
title('speech + music, x(t) = x_s(t) + x_m(t)'),grid on

%% Ducker Design
Aa = 0.997916; % lambda_a
Ar = 0.999826; % lambda_r
p = 1/5;     % rho
c0 = 0.006428;      % threshold
g1 = 0;     % for EMA gain-smoothing filter, g_n
g0 = 0;             % g_{n-1}
c = 0;          % Control signal
G = 0;          % Smoothed gain
n = 1;              % Index variable
N = length(xs);    % For limit on n

for n = 1:N
if abs(xs(n)) >= c                   % calculate control signal
    c = Aa*c + (1 - Aa)*abs(xs(n));
else
    c = Ar*c + (1 - Ar)*abs(xs(n));
end

cPlot(n) = c;    % For plotting

g0 = g1;    

if c == 0        % Calculate gain signal
    g1 = 1;
elseif c >= c0 
    g1 = (c/c0)^(p-1);
else
    g1 = 1;
end

gPlot(n) = g1;    % For plotting

if g1 >= g0     % EMA smoothing filter
    G = Aa*G + (1 - Aa)*g1;
else
    G = Ar*G + (1 - Ar)*g1;
end

GPlot(n) = G;    % For plotting

ym(n) = G*xm(n);

y(n) = xs(n) + ym(n);
end

% Output results to speakers
sound(y,fs);

% Output results to file
audiowrite('speech with ducked music.wav',y,fs);

%% Plot Results
close all
subplot(2,1,1)
plot(t,ym),title('ducked music, y_m(t)'),axis([0,7,-1,1]),grid on
subplot(2,1,2)
plot(t,y),title('speech + ducked music, y(t) = x_s(t) + G(t)x_m(t)'),
axis([0,7,-1,1]),grid on,xlabel('t (sec)')
%%
close all
subplot(3,1,1)
plot(t,cPlot,'r')
title('control signal, c(t)'),
grid on

subplot(3,1,2)
plot(t,gPlot),title('ducking gain, g(t)')
grid on

subplot(3,1,3)
plot(t,GPlot,'r')
xlabel('t (sec)'),title('smoothed gain, G(t)')
grid on
%%
close all
plot(t,xms,t,GPlot),grid on