function y = ovsave(h,x,N,type)
% This function applies the overlap-save method of fast convolution
% h = impulse response, x = input signal, N = block size, type = t for
% time-domain, f for frequency-domain

n = 1;              % Counter variable for the number of output samples
M = length(h) - 1;      % Amount of overlap
L = N - M;              % To output the correct number of samples per block

xPad = [x(:); zeros(M,1)];  % Pad zeros for the correct number of blocks
xMat = buffer(xPad,N,M);    % Create matrix of overlapping blocks
[l,m] = size(xMat);         % To determine the number of blocks, m

if type == 't'          % Time-domain convolution

for j = 1:m
    y1 = conv(h,xMat(:,j));
    % Equivalent to wrapping modulo N and discarding the wrapped samples,
    % just less operations
    y1(1:M) = [];                %Discard first M samples
    y1(end-M+1:end) = [];        %Discard last M samples
        
    for k = 1:L
        if n <= (length(x)+M)       % Output correct number of samples
            y((j-1)*L + k) = y1(k);         % Output correct values
            n = n + 1;
        end
    end
end

elseif type == 'f'
    
    H = fft(h(:),N);    % FFT of impulse response
    
for j = 1:m
    X = fft(xMat(:,j),N);   % FFT of each block
    y1 = real(ifft(H.*X));  % output of convolution for each block
    % Equivalent to wrapping modulo N and discarding the wrapped samples,
    % just less operations
    y1(1:M) = [];           %Discard first M samples
    
    for k = 1:L
        if n <= (length(x)+M)       % Output correct number of samples
            y((j-1)*L + k) = y1(k);         % Output correct values
            n = n + 1;
        end
    end
end
end
end