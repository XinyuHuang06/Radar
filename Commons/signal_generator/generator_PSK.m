function [signal_MPSK, t] = generator_PSK(varargin)
% Example: [MPSK, t] = PSK_generator(fs,fc,phase_seq,M,'noiseF',1,'SNR',10)
% :param fs: Required, the sample frequency(Hz)
% :param fc: Required, the carrier frequency(Hz)
% :param phase_seq: Required, the discrete phase seq, phase_seq(i) in {1,2,...,M}.
% :param M: the number of discrete phase
% :param noiseF: Optional Para, Default,0. 1->add awgn noise with SNR
% :param SNR: Optional Para, Default,0. the SNR of noise
% :return signal_MPSK: a Nx1 column vector
% :return t: the time sequence (Nx1 column vector)
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
% fs, fc, phase_seq, M, SNR
%------------------------------------------------------------------------------
    in_par = inputParser; % Intialization
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'fc', 0);
    addOptional(in_par, 'Rb', 0);
    addOptional(in_par, 'phase_seq', 0);
    addOptional(in_par, 'M', 0);
    addOptional(in_par, 'Rb', 0); % the rate of phase code transmission
    addParameter(in_par, 'noiseF', 0);
    addParameter(in_par, 'SNR', 0);
    parse(in_par,varargin{:});
    fs = in_par.Results.fs;
    fc = in_par.Results.fc;
    Rb = in_par.Results.Rb;
    phase_seq = in_par.Results.phase_seq;
    M = in_par.Results.M;
    noiseF = in_par.Results.noiseF;
    SNR = in_par.Results.SNR;
    % Generate the signal
    N = length(phase_seq); % the no. of phase code
    Nc = floor(fs/Rb);
    signal_MPSK = zeros(N*Nc,1);
    t = 0:1/fs:(N*Nc-1)/fs;
    for iN = 0:N-1
        signal_MPSK(iN*Nc+1:iN*Nc+Nc) = exp(1j*2*pi*fc*t(iN*Nc+1:iN*Nc+Nc)+...
        1j*2*pi*phase_seq(iN+1)/M);
    end
    if noiseF
        signal_MPSK = awgn(signal_MPSK, SNR, 'measured');
    end
end