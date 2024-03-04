function [signal_LFM, t] = generator_LFM(varargin)
% Example: [LFM, t] LFM_generator(fs,fc,B,T,'noiseF',1,'SNR',10,'type',1)
% :param fs: Required, the sample frequency(Hz)
% :param fc: Required, the carrier frequency(Hz)
% :param B: Required, the bandwidth(Hz)
% :param T: Required, the impulse time(s)
% :param noiseF: Optional Para, Default,0. 1->add awgn noise with SNR
% :param SNR: Optional Para, Default,0. the SNR of noise
% :param type: Optional Para, Default,0. 0->postive LFM, 1->negative LFM
% :return signal_LFM: a Nx1 column vector
% :return t: the time-line sequence (Nx1 column vector)
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser; % Intialization
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'fc', 0);
    addOptional(in_par, 'B', 0);
    addOptional(in_par, 'T', 0);
    addParameter(in_par, 'ratio', -1);
    addParameter(in_par, 'num_pulse', -1);
    addParameter(in_par, 'noiseF', 0);
    addParameter(in_par, 'SNR', 0);
    addParameter(in_par, 'type', 0);
    parse(in_par, varargin{:}); 
    fs = in_par.Results.fs;
    fc = in_par.Results.fc;
    B = in_par.Results.B;
    T = in_par.Results.T;
    ratio = in_par.Results.ratio;
    num_pulse = in_par.Results.num_pulse;
    noiseF = in_par.Results.noiseF;
    SNR = in_par.Results.SNR;
    type = in_par.Results.type;
    % Generate the signal
    N = ceil(T*fs); % The no. of data sample
    k = B/T; % The slope of LFM
    t = (0:1/fs:(N-1)/fs).';
    signal_carrier = exp(1j*2*pi*fc*t);
    if type
        signal_LFM = exp(1j*pi*(2*B*t-k*t.^2));
    else
        signal_LFM = exp(1j*pi*k*t.^2);
    end
    signal_LFM = signal_carrier.*signal_LFM;
    if ratio~=-1
        T_pulse = T/ratio;
        T_intra = T_pulse*num_pulse;
        N_pulse = ceil(T_intra*fs);
        N_inter = ceil(T_pulse*fs);
        t = (0:1/fs:(N_pulse-1)/fs).';
        signal_LFM_1 = zeros(N_pulse, 1);
        for i_np = 1:num_pulse
            signal_LFM_1(1+(i_np-1)*N_inter:N+(i_np-1)*N_inter) = signal_LFM;
        end
        signal_LFM = signal_LFM_1;
    end
    if noiseF
        signal_LFM = awgn(signal_LFM, SNR,'measured');
    end
end