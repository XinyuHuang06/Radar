function [S, f, alfa] = Analysis_CS_DFSM(signal,fs, df, M)
% Example:
% :param :
% :return :
% detailed description: Cyclostationary Estimator using the Direct Frequency Smoothing Method algorithm.
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 10/12/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------


x = signal(:); % Using real part of signal

N = (M*fs)/df;
N = pow2 (nextpow2(N)); % windowing record for FFT

X = fft(x,N);   % fft of the truncated (or zero padded) time series
X = fftshift(X);% shift components of fft
Xc = conj(X);                 % precompute the complex conjugate vector

S = zeros (N,N);              % size of the Spectral Correlation Density matrix
f = zeros (N,N);              % size of the frequency matrix;
alfa = zeros (N,N);           % size of the cycle frequency matrix
F = fs/(2*N);                 % precompute constants -  F = fs/(2*N);     
G = fs/N;                     % precompute constants -  G = fs/N;  
m = -M/2+1:M/2;               % set frequency smoothing window index

for k = 1:N                                % fix k
    % computes vectors of f and alfa,
    % store frequency and cycle frequency data for given k.
    k1 = 1:N;
    f(k,k1) = F*(k+k1-1) - fs/2;          % Computes f values and shift them to center in zero (f = (K+L)/2N) [1]
    alfa(k,k1) = G*(k-k1 + N-1) - fs;       % Computes alfa values and shift them to center in zero (alfa = (K-L)/N) [1]
    for k1 = 1:N % fix k1 = J

        % calculate X(K+m) & conj (X(J+m)) for arguments of X(1:N) only
        B = max(1-k, 1-k1);          % Largest min of 1 <= (K+m)| (J+m) <= N
        A = min (N-k, N-k1);         % Smallest max of 1 <= (K+m)| (J+m) <= N
        n = m((m<=A) & (m>=B));   %fix the index out of range problem by
                                                   % truncating the window
        if isempty(n)
            S(k,k1) = 0;
        else
            p = k+n;
            q = k1+n;
            Y = X(p).*Xc(q);
            S(k,k1) = sum(Y);
        end
    end
end

% S = abs(S./max(max(S)));% normalize output matrix
S = real(S./M);
mesh(alfa, f, S); grid;
xlabel('Cycle frequency (Hz)'); ylabel('Frequency (Hz)');
title (['Frequency Smoothing SCD ', ', df = ', int2str(df),', N = ', int2str(N)]);
colorbar;