%% 1. Phase Vocoder (phvoc function)
clear, clc
[x, fs] = audioread('04 Finger Eleven - One Thing.mp3', [44100*25, 44100*35 - 1]); % 10 sec audio sample
x = x(:,1); % as a column

%%
soundsc(x,fs)

%% 1.a Phase Vocoder, r > 1
r = 1.4;
Ra = 256;
N = 4096;
y = phvoc(x,r,Ra,N);
y = y/max(abs(y));

%%
audiowrite('phvoc_r_gt_1.wav', [x', zeros(1,fs), y'],fs);

%% 1.b Phase Vocoder, r < 1
r = 0.7;
Ra = 256;
N = 4096;
y = phvoc(x,r,Ra,N);
y = y/max(abs(y));

%%
audiowrite('phvoc_r_lt_1.wav', [x', zeros(1,fs),y'],fs);

%% 2. Pitch Modification (pitchmod function) 
[x,fs] = audioread('flute2.wav');
x = x(:,1);

%% 2.I Spectrum Plot, Original Signal
N = 4096;
Xmag = abs(fft(x,N));               % Take the FFT of x
Xmag = fftshift(Xmag'/max(Xmag));   % Normalize and FFT shift
Xmag(1:2048) = [];                  % Only look at positive frequencies
fmax = 1800;                        % Only look up to 1800 Hz.
kmax = round(fmax*N/fs);            % DTF index at around 1800 Hz
f = (0:1/kmax:(1-1/kmax))*fmax;     % frequencies from 0 to 1800 Hz
Xmag(kmax + 1:end) = [];            % frequencies from 0 to 1800 Hz
close all
plot(f,Xmag,398,1,'r*'); grid on; xlabel('f (Hz)'); ylabel('|X(f)|');
title('original spectrum'); axis([0,fmax,0,1.1])

%% 2.II Spectrogram, Original Signal
R = 256;        % Hop size
N = 4096;       % FFT length
X = stft(x,R,N);    % STFT of x
M = size(X,2) - 1;  
k = 0:N/2;      % positive-frequency half of Nyquist interval
f = k*fs/N;     % f from 0 to fs/2
t = (0:M)*R/fs; % M + 1 time frames spaced by R*Ts = R/fs

Xmag = abs(X(k+1,:)); % extract positive frequency half
Xmag = Xmag/max(max(Xmag));

S = 20*log10(Xmag);      % convert to dB

close all
surf(t,f,S,'edgecolor','none') % Create spectrogram with time on the 
% horizontal axis and frequency on the vertical axis
colormap(jet); colorbar; axis tight; view(0,90); % View from above
xlabel('{\itt} (sec)'); ylabel('{\itf} (Hz)'); 
ylim([0,1800]);             % Show frequency range from 0 1800 Hz
title('spectrogram, original signal')

%% 2.a Pitch Modification, r = 2
r = 2;
Ra = 256;
N = 2048;
y = pitchmod(x, r, Ra, N);
y = y/max(abs(y));
%%
soundsc([x', zeros(1,fs),y'],fs);
%%
audiowrite('pitchmod_r_equals_2.wav', [x', zeros(1,fs),y'],fs);

%% 2.a.I Spectrum Plot, r = 2
N = 4096;
Ymag = abs(fft(y,N));               % Take the FFT of x
Ymag = fftshift(Ymag'/max(Ymag));   % Normalize and FFT shift
Ymag(1:2048) = [];                  % Only look at positive frequencies
fmax = 1800;                        % Only look up to 1800 Hz.
kmax = round(fmax*N/fs);            % DTF index at around 1800 Hz
f = (0:1/kmax:(1-1/kmax))*fmax;     % frequencies from 0 to 1800 Hz
Ymag(kmax + 1:end) = [];            % frequencies from 0 to 1800 Hz
close all
plot(f,Ymag,784,1,'r*'); grid on; xlabel('f (Hz)'); ylabel('|X(f)|');
title('pitch-shifted spectrum, r = 2'); axis([0,fmax,0,1.1])

%% 2.a.II Spectrogram, r = 2
R = 256;        % Hop size
N = 4096;       % FFT length
Y = stft(y,R,N);    % STFT of y
M = size(Y,2) - 1;  
k = 0:N/2;      % positive-frequency half of Nyquist interval
f = k*fs/N;     % f from 0 to fs/2
t = (0:M)*R/fs; % M + 1 time frames spaced by R*Ts = R/fs

Ymag = abs(Y(k+1,:)); % extract positive frequency half
Ymag = Ymag/max(max(Ymag));

S = 20*log10(Ymag);      % convert to dB

close all
surf(t,f,S,'edgecolor','none') % Create spectrogram with time on the 
% horizontal axis and frequency on the vertical axis
colormap(jet); colorbar; axis tight; view(0,90); % View from above
xlabel('{\itt} (sec)'); ylabel('{\itf} (Hz)'); 
ylim([0,1800]);             % Show frequency range from 0 1800 Hz
title('spectrogram, pitch-shifted, r = 2')

%% 2.b Pitch Modification, r = 1/2
r = 1/2;
Ra = 256;
N = 2048;
y = pitchmod(x, r, Ra, N);
y = y/max(abs(y));

%%
soundsc([x', zeros(1,fs),y'],fs);
%%
audiowrite('pitchmod_r_equals_0.5.wav', [x', zeros(1,fs),y'],fs);

%% 2.b.I Spectrum Plot, r = 1/2
N = 4096;
Ymag = abs(fft(y,N));               % Take the FFT of y
Ymag = fftshift(Ymag'/max(Ymag));   % Normalize and FFT shift
Ymag(1:2048) = [];                  % Only look at positive frequencies
fmax = 1800;                        % Only look up to 1800 Hz.
kmax = round(fmax*N/fs);            % DTF index at around 1800 Hz
f = (0:1/kmax:(1-1/kmax))*fmax;     % frequencies from 0 to 1800 Hz
Ymag(kmax + 1:end) = [];            % frequencies from 0 to 1800 Hz
close all
plot(f,Ymag,185,1,'r*'); grid on; xlabel('f (Hz)'); ylabel('|X(f)|');
title('pitch-shifted spectrum, r = 1/2'); axis([0,fmax,0,1.1])

%% 2.b.II Spectrogram, r = 1/2
R = 256;        % Hop size
N = 4096;       % FFT length
Y = stft(y,R,N);    % STFT of y
M = size(Y,2) - 1;  
k = 0:N/2;      % positive-frequency half of Nyquist interval
f = k*fs/N;     % f from 0 to fs/2
t = (0:M)*R/fs; % M + 1 time frames spaced by R*Ts = R/fs

Ymag = abs(Y(k+1,:)); % extract positive frequency half
Ymag = Ymag/max(max(Ymag));

S = 20*log10(Ymag);      % convert to dB

close all
surf(t,f,S,'edgecolor','none') % Create spectrogram with time on the 
% horizontal axis and frequency on the vertical axis
colormap(jet); colorbar; axis tight; view(0,90); % View from above
xlabel('{\itt} (sec)'); ylabel('{\itf} (Hz)'); 
ylim([0,1800]);             % Show frequency range from 0 1800 Hz
title('spectrogram, pitch-shifted, r = 1/2')

%% 2.c A Closer Look at the Sample Level
fs = 44100; % Sampling Rate
load x4.mat;    % 4 msec sample of flute playing

t = (0:fs*0.004)*1000/ fs; % 4 msec in units of msec

close all
stem(t,x4); xlabel('t (msec)'), title('original at rate f_s')

% Speed up by r =2
r = 2;
[p,q] = rat(1/r);

y2 = resample(x4,p,q);          % resample at fs/r, r = 2
t2 = (0:fs/r*0.004)*r*1000/fs; % 4 msec in units of msec

close all
stem(t2,y2,'r'); xlabel('t (msec)'), 
title('resampled at rate f_s/r, r = 2')

% Resample at fs/r and playback at fs

y3 = [y2,zeros(1,length(y2)-1)]; % resmpled signal is contained in the
% first 2 milliseconds
close all
stem(t,y3,'r'); xlabel('t (msec)'), 
title('resampled at f_s/r, replayed at f_s')

% Extend in time by r with phase vocoder at 1/r
Ra = 4;
N = 32;
y4 = phvoc(y2,1/r,Ra,N);
y4 = y4/max(abs(y4));       % Normalize
y4(178:end) = [];           % Make duration equal to that of original
% by removing the zeros padded

close all
stem(t,y4); xlabel('t (msec)'), 
title('phase vocoder, time-stretched by factor r = 2')

%% 3. COLA Property
N = 128;    % Length of window
M = 10;     % Number of blocks
n = 0:N-1;  
w = 0.5 - 0.5*cos(2*pi*n/(N-1)); % Hanning window
t = 1:800;

%% 3.a R = N/2
R = N/2;

wbuff = repmat(w',1,M+1);   % Matrix of w(n) repeated M+1 times

wOla = zeros(1,800);        % Overlap-added window

for m = 0:M                     % Number of blocks (columns)     
    y = wbuff(:,m+1); 
    for n = 1:N                             % Overlap-add each block
        wOla(m*R + n) = wOla(m*R + n) + y(n);  
    end
end

close all
plot(t,wOla),axis([0,800,0,3]), title('window overlap-add, R = N/2'), 
xlabel('n');

r = 0:R-1;                  % DFT index
wr = 2*pi*r/R;              % DFT frequencies

% Calculate DFTs according to equation (7)
wDFT = zeros(1,R);
for r = 0:R-1
    for n = 0:N-1
        wDFT(r+1) = wDFT(r+1) + w(n+1)*exp(-1i*wr(r+1)*n);
    end
end

wDFT = abs(wDFT)/max(abs(wDFT));  % Normalize to unity

r = 0:R-1;                  % DFT index

close all
plot(r,wDFT, 'r*'); grid on; axis([-10,70,-0.1,1.1]);
xlabel('DFT index, k'), title('R-point DFT, W(k)'),

%% 3.b R = 3N/8
R = 3*N/8;

wbuff = repmat(w',1,M+1);   % Matrix of w(n) repeated M+1 times

wOla = zeros(1,800);        % Overlap-added window

for m = 0:M                     % Number of blocks (columns)     
    y = wbuff(:,m+1); 
    for n = 1:N                             % Overlap-add each block
        wOla(m*R + n) = wOla(m*R + n) + y(n);  
    end
end

close all
plot(t,wOla),axis([0,800,0,3]), title('window overlap-add, R = 3N/8'), 
xlabel('n');

r = 0:R-1;                  % DFT index
wr = 2*pi*r/R;              % DFT frequencies

% Calculate DFTs according to equation (7)
wDFT = zeros(1,R);
for r = 0:R-1
    for n = 0:N-1
        wDFT(r+1) = wDFT(r+1) + w(n+1)*exp(-1i*wr(r+1)*n);
    end
end

wDFT = abs(wDFT)/max(abs(wDFT)); % Normalize to unity

r = 0:R-1;                  % DFT index

close all
plot(r,wDFT, 'r*'); grid on; axis([-10,70,-0.1,1.1]);
xlabel('DFT index, k'), title('R-point DFT, W(k)'),

%% 3.c R = N/3
R = round(N/3);

wbuff = repmat(w',1,M+1);   % Matrix of w(n) repeated M+1 times

wOla = zeros(1,800);        % Overlap-added window

for m = 0:M                     % Number of blocks (columns)     
    y = wbuff(:,m+1); 
    for n = 1:N                             % Overlap-add each block
        wOla(m*R + n) = wOla(m*R + n) + y(n);  
    end
end

close all
plot(t,wOla),axis([0,800,0,3]), title('window overlap-add, R = N/3'), 
xlabel('n');

r = 0:R-1;                  % DFT index
wr = 2*pi*r/R;              % DFT frequencies

% Calculate DFTs according to equation (7)
wDFT = zeros(1,R);
for r = 0:R-1
    for n = 0:N-1
        wDFT(r+1) = wDFT(r+1) + w(n+1)*exp(-1i*wr(r+1)*n);
    end
end

wDFT = abs(wDFT)/max(abs(wDFT)); % Normalize to unity

r = 0:R-1;                  % DFT index

close all
plot(r,wDFT, 'r*'); grid on; axis([-10,70,-0.1,1.1]);
xlabel('DFT index, k'), title('R-point DFT, W(k)'),

%% 3.d R = N/4
R = N/4;

wbuff = repmat(w',1,M+1);   % Matrix of w(n) repeated M+1 times

wOla = zeros(1,800);        % Overlap-added window

for m = 0:M                     % Number of blocks (columns)     
    y = wbuff(:,m+1); 
    for n = 1:N                             % Overlap-add each block
        wOla(m*R + n) = wOla(m*R + n) + y(n);  
    end
end

close all
plot(t,wOla),axis([0,800,0,3]), title('window overlap-add, R = N/4'), 
xlabel('n');

r = 0:R-1;                  % DFT index
wr = 2*pi*r/R;              % DFT frequencies

% Calculate DFTs according to equation (7)
wDFT = zeros(1,R);
for r = 0:R-1
    for n = 0:N-1
        wDFT(r+1) = wDFT(r+1) + w(n+1)*exp(-1i*wr(r+1)*n);
    end
end

wDFT = abs(wDFT)/max(abs(wDFT)); % Normalize to unity

r = 0:R-1;                  % DFT index

close all
plot(r,wDFT, 'r*'); grid on; axis([-10,70,-0.1,1.1]);
xlabel('DFT index, k'), title('R-point DFT, W(k)'),