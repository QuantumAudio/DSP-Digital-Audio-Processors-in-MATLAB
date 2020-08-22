%% Compressor

%% Defining Input Signal

fs = 8000; % sampling rate
Ts = 1/fs; % sampling period
ta = 0.002; % attack time
tr = 0.01; % release time
Aa = exp(-2.3*Ts/ta); % attack forgetting factor
Ar = exp(-2.3*Ts/tr); % release forgetting factor

% Defining three sinusoids to be played consecutively
A1 = 2; % Amplitude of first sinusoid
f1 = 300; % f1 = 0.3 kHz
t1 = [0:Ts:0.025 - Ts]; % t1 from 0 to 25 msec
A2 = 4; % Amplitude of second sinusoid
f2 = 600; % f2 = 0.6 kHz
t2 = [0.025:Ts:0.05 - Ts]; % t2 from 25 msec to 50 msec
A3 = 0.5; % Amplitude of third sinusoid
f3 = 1200; % f3 = 1.2 kHz
t3 = [0.05:Ts:0.075 - Ts]; % t3 from 50 msec to 75 msec
x = [A1*cos(2*pi*f1*t1),A2*cos(2*pi*f2*t2),A3*cos(2*pi*f3*t3)];

envInTop = [A1*ones(1,200),A2*ones(1,200),A3*ones(1,200)];
envInBot = [-A1*ones(1,200),-A2*ones(1,200),-A3*ones(1,200)];
t = [t1,t2,t3];

%% Plot Input Signal
plot(t*1000,x,t*1000,envInTop,'--g',t*1000,envInBot,'--g'); 
axis([0,75,-4.2,4.2]); legend('x(t)','envelope'), 
text(60,3.8,'f_1 = 0.3 kHz');text(60,3.3,'f_2 = 0.6 kHz'); 
text(60,2.8,'f_3 = 1.2 kHz'); xlabel('t (msec)'), title('input signal, x(t)')

%% Compressor Design
Aa = 0.8661; % lambda_a
Ar = 0.9717; % lambda_r
p = 1/3;     % rho
c0 = 1;      % threshold
L = ceil((1 + Aa)/(1 - Aa)); % Length of FIR gain-smoothing filter
g = zeros(1,L); % Internal buffer for FIR gain-smoothing filter
q = 1;              % Internal buffer index
c = 0;          % Control signal
G = 0;          % Smoothed gain
n = 1;              % Index variable
N = length(x);    % For limit on n

clear y
for n = 1:N
if abs(x(n)) >= c                   % calculate control signal
    c = Aa*c + (1 - Aa)*abs(x(n));
else
    c = Ar*c + (1 - Ar)*abs(x(n));
end

cPlot(n) = c;    % For plotting

if c == 0        % Calculate gain signal
    g(q) = 1;
elseif c >= c0 
    g(q) = (c/c0)^(p-1);
else
    g(q) = 1;
end

gPlot(n) = g(q);    % For plotting

q = q + 1;
if q > L
    q = 1;
end

G = (1/L)*sum(g); % FIR averaging filter

GPlot(n) = G;    % For plotting

y(n) = G*x(n); % Output
end

% Output results to speakers
sound(y,fs);

% Output results to file
audiowrite('Compressed.wav',y,fs);


%% Plot Results
subplot(2,2,1)
plot(t*1000,y,t*1000,envInTop,'--g',t*1000,envInBot,'--g',...
t*1000,c0*ones(1,N),'--r',t*1000, -c0*ones(1,N),'--r')
axis([0,75,-4.2,4.2]),xlabel('t (msec)'),
title('compressor, y(t) = G(t)x(t)'), text(60,3.8,'\rho = 1/3')


subplot(2,2,2)
plot(t*1000,cPlot,'r',t*1000, c0*ones(1,N),'--b')
axis([0,75,0,4]),xlabel('t (msec)'),title('control signal, c(t)'),
legend('c(t)','c_0')

subplot(2,2,3)
plot(t*1000,gPlot,'r')
axis([0,75,0,1.2]),xlabel('t (msec)'),title('gain, g(t) = F(c(t))')
grid on

subplot(2,2,4)
plot(t*1000,GPlot)
axis([0,75,0,1.2]),xlabel('t (msec)'),title('smoothed gain, G(t)')
grid on
