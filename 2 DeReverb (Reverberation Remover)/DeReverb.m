%% Appendix 1.b: DeReverb (Reverberation Removal)
% Import audio file
[x, fs] = audioread('reverb.wav');
% Import original to compare the results
[s, fs] = audioread('original.wav');
% Pulse signal
[p, fs] = audioread('pulse.wav');
% Pulse recorded in room acoustics
[r, fs] = audioread('pulse_rec.wav');
% Obtain the impulse response of the room acoustics 
h = filter(1,p,r);

%% Plot the signals
t = 1:24000;
subplot(3,1,1);
plot(t/fs,p); title('1 kHz pulse'); axis([0,1.5,-1.5,1.5]);
ylabel('p(t)');
subplot(3,1,2); 
plot(t/fs,r); title('recorded pulse'); axis([0,1.5,-1.5,1.5]);
ylabel('r(t)');
subplot(3,1,3); 
plot(t/fs,h,'Color','r'); title('impulse response'); xlabel('t (sec)');
axis([0,1.5,0,1.5]); ylabel('h(t)');

%% Remove the Room Reverberation
y = filter(1,h,x);

%% Plot the Results
t = 1:10*fs;
subplot(3,1,1);
plot(t/fs,x); title('recorded, reverb'); 
subplot(3,1,2); 
plot(t/fs,y,'Color','r'); title('recovered'); 
subplot(3,1,3); 
plot(t/fs,s); title('original'); xlabel('t (sec)');

%% Output the Recovered Signal
audiowrite('ReverbRemoved.wav',y,fs);
