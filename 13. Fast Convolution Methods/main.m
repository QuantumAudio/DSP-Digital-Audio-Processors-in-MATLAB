%% Overlap-add & Overlap-save Fast Convolution Methods
x = 1:19; % Input signal

h = [1 2 3 4]; % Filter

M = length(h) - 1; % Order of filter

N = 2^ceil(log2(2*M + 1)); % N is the nearest power of 2 above 2*M

yconv = conv(h,x)

yovaddT = ovadd(h,x,N,'t')
yovaddF = ovadd(h,x,N,'f')
yovsaveT = ovsave(h,x,N,'t')
yovsaveF = ovsave(h,x,N,'f')

%% Matched Filter vs. Cross-Correlation
clear, clc
load xrec;
fs = 1e6; Ts = 1/fs;            % Sampling rate and period
t = 0:Ts*1e3:1;                 % Duration of pulse
s = @(t) sin(20*pi*t + 10*pi*t.^2); % Pulse signal
T = 1;                          % Max time value of pulse
h = s(T-t);                     % Matched filter impulse response
M = length(h) - 1;              % Order of matched filter
N = 2^ceil(log2(2*M + 1));      % N is the nearest power of 2 above 2*M

ymatch = ovsave(h,xrec,N,'t');           % Matched filter output
ymatch(length(xrec)+1:end) = [];         % Truncate convolution tail
ymatch = ymatch/max(ymatch);

t5 = 0:Ts*1e3:5;                         % Duration of received signal
s5 = [s(t),zeros(1,4000)];   % Extend pulse signal length to match received signal length 

xCorr = xcorr(xrec,s5);     % Compute the cross-correlation between the sent pulse and the received signal
xCorr(1:5000) = [];% Eliminate negative portion of xCorr to only look at positive values of delay
xCorr = xCorr/max(xCorr);

close all
plot(t5,s5,'k',t5,ymatch,t5,xCorr,'r--'),legend('transmitted','matched filter','correlator')
axis([0,5,-1.5,1.5]);title('processed received signal'),xlabel('t (msec)'), grid on

%% Just for Plots
close all
subplot(1,2,1)
plot(t,s(t)), axis([0,1,-1.5,1.5]), grid on, xlabel('t (msec)'),
title('chirped sinusoid')
subplot(1,2,2)
plot(t,h),axis([0,1,-1.5,1.5]), grid on, xlabel('t (msec)'),
title('matched filter')


close all
as = 0.5*s(t - 3);
as5 = [zeros(1,3000),as,zeros(1,1000)];
xrec = xrec/max(xrec);

subplot(2,1,1)
plot(t5,s5,'k',t5,as5,'r'),legend('transmitted','reflection')
axis([0,5,-1.5,1.5]);title('received noise-free signal'),xlabel('t (msec)'), grid on
subplot(2,1,2)
plot(t5,xrec,'r'),
axis([0,5,-1.5,1.5]);title('noisy received signal'),xlabel('t (msec)'), grid on

