function varargout = Analysis_CS_DFSM(varargin)
% Example: [S, f, alfa] = Analysis_CS_DFSM(fs, signal, df, M, 'bool_draw', 1)
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
    in_par = inputParser;
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'signal', 0);
    addOptional(in_par, 'df', 0);
    addOptional(in_par, 'M', 32);
    addParameter(in_par, 'bool_draw', 1);
    parse(in_par, varargin{:});
    fs = in_par.Results.fs;
    signal = in_par.Results.signal;
    M = in_par.Results.M;
    df = in_par.Results.df;
    bool_draw = in_par.Results.bool_draw;

    x = real(signal); % Using real part of signal
    x = x(:);
    N = length(signal);
    N = (M*fs)/df;
    N = pow2 (nextpow2(N)); % windowing record for FFT

    X = fft(x,N);               % fft of the truncated (or zero padded) time series
    X = fftshift(X);            % shift components of fft
    X = X*2/sqrt(N);
    Xc = conj(X);               % precompute the complex conjugate vector

    S = zeros (N,N);            % size of the Spectral Correlation Density matrix
    f = zeros (N,N);            % size of the frequency matrix;
    alfa = zeros (N,N);         % size of the cycle frequency matrix
    F = fs/(2*N);               % precompute constants -  F = fs/(2*N);     
    G = fs/N;                   % precompute constants -  G = fs/N;  
    m = -M/2+1:M/2;             % set frequency smoothing window index

    for k = 1:N                                % fix k
        % computes vectors of f and alfa,
        % store frequency and cycle frequency data for given k.
        k1 = 1:N;
        f(k,k1) = F*(k+k1-1) - fs/2;          % Computes f values and shift them to center in zero (f = (K+L)/2N) [1]
        alfa(k,k1) = G*(k-k1 + N-1) - fs;       % Computes alfa values and shift them to center in zero (alfa = (K-L)/N) [1]
        % f(k,k1) = k + k1;
        % alfa(k,k1) = k - k1;
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
    S = abs(real(S));
    S = S/max(S(:));
    % S = S/max(S(:));
    % S = abs(S./max(S(:)));
    f = f/fs;
    alfa = alfa/fs;
    if bool_draw
        % contour(alfa, f, S,'LevelStep',0.02); grid;
        mesh(f, alfa, S);
        xlabel('Frequency (Hz)');ylabel('Cycle frequency (Hz)');
        title (['Frequency Smoothing SCD ']);
        colorbar;
        SetDrawStyle;
    end
    if nargout == 1
        varargout{1} = S;
    elseif nargout == 3
        varargout{1} = S;
        varargout{2} = f;
        varargout{3} = alfa;
    end
end
