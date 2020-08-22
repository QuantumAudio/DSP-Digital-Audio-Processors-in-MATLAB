function X = stft(x, Ra, N)
% Calculates the STFT of the input signal x with analysis hop size Ra and
% FFT length N

M = floor(length(x)/Ra);        % Determine how many blocks, (will be M+1)
Lext = M*Ra + N;                % Length of extended input signal
x(length(x)+1:Lext) = 0;        % Zero pad input signal
Xbuff = buffer(x, N, N-Ra,'nodelay'); % Separate signal into size N blocks with overlap N-Ra
n = 1:N;                            % For window of length N
w = 0.5 - 0.5*cos(2*pi*n/(N-1));    % Hanning window
w = w';                             % Transpose window
W = repmat(w,1,M+1);                % Window each block
X = fft(W.*Xbuff,N);                % Take the FFT of windowed blocks
end

