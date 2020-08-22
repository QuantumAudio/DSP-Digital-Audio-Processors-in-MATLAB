% Maintains correct phase under time scale modifications
function Y = phmap(X,r,Ra,N)

mod2pi = @(x) mod(x + pi, 2*pi) - pi;

[~,M] = size(X);        % Determine how many columns, M+1, in X
M = M-1;                % Determine the value of M
Rs = round(Ra/r);       % Determine synthesis hop size

Ymag = abs(X);          % Preserve magnitude of input
Phi = angle(X);         % Transfer phase angles of input

k = 0:(N-1);
wk = 2*pi*k'/N;

Psi(:,1) = Phi(:,1);    % Psi(k,0) = Phi(k,0)

    for m = 1:M
        dw = (1/Ra)*mod2pi(Phi(:,m+1) - Phi(:,m) - Ra*wk); % Bring omega within the Nyquist interval
        Psi(:,m+1) = Psi(:,m) + Rs*(wk + dw);
    end
    Y = Ymag.*exp(1i.*Psi);
end
