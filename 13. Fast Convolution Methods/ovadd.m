function y = ovadd(h,x,N,type)
% This function applies the overlap-add method of fast convolution
% h = impulse response, x = input signal, N = block size, type = t for
% time-domain, f for frequency-domain

n = 1;                  % Counter variable for the number of output samples
M = length(h) - 1;      % Amount of overlap
L = N - M;              % Length of each block size

xPad = [x(:); zeros(M,1)];
xMat = buffer(xPad,L);      % Create matrix of nonoverlapping blocks
[l,m] = size(xMat);     % To determine the number of blocks, m

ytemp = zeros(M,1);     % Temporary buffer for convolution tails

if type == 't'          % Time-domain convolution

for j = 1:m
    y1 = conv(h,xMat(:,j)); % Convolve each block

    for k = 1:M
        y1(k) = y1(k) + ytemp(k);   % Add blocks with M overlaps
        ytemp(k) = y1(k+L);         % Input the overlaps for the next block
    end

    for k = 1:L
        if n <= (length(x)+M)       % Output correct number of samples
            y((j-1)*L + k) = y1(k);         % Output correct values
            n = n + 1;
        end
    end


end

elseif type == 'f'     % Frequency-domain convolution
    
    H = fft(h(:),N);    % FFT of impulse response
    
for j = 1:m
    X = fft(xMat(:,j),N);   % FFT of each block
    y1 = real(ifft(H.*X));  % output of convolution for each block

    for k = 1:M
        y1(k) = y1(k) + ytemp(k);   % Add blocks with M overlaps
        ytemp(k) = y1(k+L);         % Input the overlaps for the next block
    end

    for k = 1:L
        if n <= (length(x)+M)       % Output correct number of samples
            y((j-1)*L + k) = y1(k);         % Output correct values
            n = n + 1;
        end
    end

end
end
end