function y = phvoc(x, r, Ra, N)

Rs = round(Ra/r);       % Determine synthesis hop size

% This is necessary to determine the correct output length when expanding
% the signal in time
if r < 1
    M = floor(length(x)/Ra);
    Ly = M*Rs + N;
    x(end + 1:Ly) = 0;  % Zero pad input signal to be the length of the output signal
end

X = stft(x, Ra, N);     % Calculate the STFT
Y = phmap(X, r, Ra, N); % Map the phases to preserve them
y = istft(Y, Rs, N);    % Calculate the ISTFT

end