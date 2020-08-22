% Takes the inverse STFT of Y
function y = istft(Y, Rs, N)
[~,M] = size(Y);        % Determine how many columns, M+1, in Y
M = M-1;                % Determine the value of M
Ly = M*Rs + N;          % Determine length of the output signal
y = zeros(Ly,1);        % Initialize y(n) to zeros
n = 0:(N-1);                         % For window of length N
w = 0.5 - 0.5*cos(2*pi*n/(N-1)); % Hanning window
w = w';                          % Transpose window


Ybuff = real(ifft(Y,N));              % Take the IFFT of Y

for m = 0:M                     % Number of blocks (columns)          
    ym = Ybuff(:,m+1);          % Vectorize windowing
    ym = ym.*w;                 % Window each block of Ybuff
    for n = 1:N                             % Overlap-add each block
        y(m*Rs + n) = y(m*Rs + n) + ym(n);  % ym(n) is already windowed
    end
end

end