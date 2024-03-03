% Example:
% :param fs:
% :param subN:
% :param N:
% :return signal_PRPC:
% detailed description: Pseudo-random phase code(PRPC) signal
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 28/02/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function [signal_PRPC, t] = generator_PRPC(fs, subN, N)
    t = (0:1/fs:N*subN/fs).';
    signal_PRPC = zeros(size(t));
    for i_N = 0:N-1
        phi = rand*2*pi;
        signal_PRPC(i_N*subN+1:i_N*subN+subN) = exp(1j*phi);
    end
end