% Modifies pitch without changing duration
function y = pitchmod(x, r, Ra, N)

[p,q] = rat(1/r);           % Rationalize spead up factor for resampling

if p <= q                   % If r >= 1 resample at 1/r first, then phvoc
    x = resample(x, p, q);  % Resample at 1/r
    y = phvoc(x, 1/r, Ra, N);   % Map phases to 1/r
else                            % If r < 1 phvoc, then resample at 1/r
    x = phvoc(x, 1/r, Ra, N);   % Map phases to 1/r
    y = resample(x, p, q);      % Resample at 1/r
end

end
