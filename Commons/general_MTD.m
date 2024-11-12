% Example: 
% :param M: a int number, The accumulation numbers
% :param signal: a vector, The signal
% :param fc: the carrier frequency
% :param PRI: the Pulse Repeation Inter
% :param TargetInfo: a M*3 matrix, the one-colum indict the location, the two-colum indict the velocity
%                      M is the numbers of target
% :return : The Range-Doppler Map
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 27/09/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function varargout = general_MTD(signal, M, fc, fs, PRI, TargetInfo, SNR)
    numTargets = size(TargetInfo, 1);
    c = 3*e8;
    N = round(PRI*fs);
    len = length(signal);
    echo_matrix = zero(len,M);
    for i = 1:numTargets
        range = TargetInfo(i,1);
        velocity = TargetInfo(i,2);
        RCS = TargetInfo(i,3);
        bias_time = 2*range/(c+velocity);
        bias_fre = 2*velocity/(c+velocity);
        temp_echo = zero(len,1);
        t_seq = 0:1/fs:(len-1)/fs;
        bias_freseq = exp(1j*2*pi*bias_fre*t_seq);
        temp_echo = signal.*exp(-j*2*pi*bias_time).*bias_freseq;
        echo_matrix(:,i) = circshift(temp_echo, round(bias_time*fs));
    end
    matchedFilter = conj(fliplr(signal));
    rangeDopplerMap = zeros(M, length(signal));
    % Range-Match Doppler-Match
    for i = 1:M
        pulseSegment = echo_matrix(:, i);
        rangeDopplerMap(:, i) = abs(ifft(fft(pulseSegment) .* fft(matchedFilter, length(pulseSegment))));
    end
    varargout{1} = rangeDopplerMap;
end

