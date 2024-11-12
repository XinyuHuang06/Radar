function [signal_OFDM, t] = generator_OFDM(varargin)
% 
% :param fs: sampling frequency
% :param fc: the frequency of carrier
% :param 
% :param T1: 
% :param M: 
% :return signal_OFDM:
% :return t:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 25/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    % Intialization
    in_par = inputParser;
    addOptional(in_par, 'fs', 0); % the sampling frequency
    addOptional(in_par, 'fc_1', 0); % 子载波起始频率
    addOptional(in_par, 'T1', 0); % OFDM符号持续时间
    addOptional(in_par, 'M', 0); % the number of sub-carriers
    addOptional(in_par, 'N', 0); % the number of OFDM character
    addParameter(in_par, 'm_type', 0); % the type of mudulation
    addParameter(in_par, 'code', -1); % the code sequence of modulation (for BPSK and QPSK)
    parse(in_par, varargin{:});
    fs = in_par.Results.fs;
    fc_1 = in_par.Results.fc_1;
    T1 = in_par.Results.T1;
    M = in_par.Results.M;
    N = in_par.Results.N;
    m_type = in_par.Results.m_type;
    code = in_par.Results.code;

    Num = T1*fs; % 每个符号采样点数
    t = (0:1/fs:(N*Num-1)/fs).'; % 时间坐标序列
    signal_OFDM = zeros(size(t)); % 信号序列
    inter_f = 1/T1;
    switch m_type
        case 'BPSK'
            if code == -1 % 若无数据传入，随机生成调制序列
                code = randi([0,1], M, N);
                code(code==0) = -1;
            end

            chouqu = zeros(M,N); % 随机抽取部分子载波归零，用于约束峰均比
            for iter_N = 0:N-1
                index = randperm(M,M/4);
                chouqu(index,iter_N+1) = 1;
            end

            for iter_N = 0:N-1
                sub_carrier = zeros(Num,M);
                for iter_M = 0:M-1
                    if chouqu(iter_M+1,iter_N+1)
                    sub_carrier(:,iter_M+1) = code(iter_M+1,iter_N+1)*exp(1j*2*pi*(iter_M*inter_f)*t(iter_N*Num+1:iter_N*Num+Num));
                    end
                end
                signal_OFDM(iter_N*Num+1:iter_N*Num+Num) = sum(sub_carrier,2);
            end
            signal_OFDM = signal_OFDM.*exp(1j*2*pi*fc_1*t);
        case 'QPSK'
            if code == -1
                code = randi([0,3], M, 1);
                code = code/pi;
            end
            for iter_M = 0:M-1
                temp_signal = exp(1j*2*pi*(iter_M*inter_f)*t + code(iter_M+1));
                signal_OFDM = signal_OFDM + temp_signal;
            end
        case 'LFM'
            chouqu = zeros(M,N);
            for iter_N = 0:N-1
                index = randperm(M,M/4);
                chouqu(index,iter_N+1) = 1;
            end
            for iter_N = 0:N-1
                sub_carrier = zeros(Num,M);
                for iter_M = 0:M-1
                    if chouqu(iter_M+1,iter_N+1)
                    sub_carrier(:,iter_M+1) = exp(1j*2*pi*(iter_M*inter_f)*t(iter_N*Num+1:iter_N*Num+Num)+1j*pi*(inter_f/T1).*t(iter_N*Num+1:iter_N*Num+Num).^2);
                    end
                end
                signal_OFDM(iter_N*Num+1:iter_N*Num+Num) = sum(sub_carrier,2);
            end
        case 'RM' % Random modulated 
        case 'None' % 不进行任何调制
            for iter_N = 0:N-1
                sub_carrier = zeros(Num,M);
                for iter_M = 0:M-1
                    sub_carrier(:,iter_M+1) = exp(1j*2*pi*(iter_M*inter_f)*t(iter_N*Num+1:iter_N*Num+Num));
                end
                signal_OFDM(iter_N*Num+1:iter_N*Num+Num) = sum(sub_carrier,2);
            end
            signal_OFDM = signal_OFDM.*exp(1j*2*pi*fc_1*t);            
    end
    signal_OFDM = signal_OFDM/max(abs(signal_OFDM));
end
