% Example:
% :param s: a NxX matrix, the N is the lengths of signal , the X is the numbers of signal
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 28/02/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function [seq_eps] = Evaluate_eps(s, M, seq_SNR)
    N = size(s,1);
    X = size(s,2);
    num = 100;
    seq_eps = zeros(X, size(seq_SNR, 2));
    figure
    for i_X = 1:X
        temp_s0 = s(:,i_X);
        temp_curve = zeros(size(seq_SNR));
        temp_value = 0;
        i_index = 0;
        for SNR = seq_SNR
            for i_n = 1:num
                temp_s1 = awgn(temp_s0, SNR,'measured');
                Z = zeros(2*M,N-M+1);
                for i_c = 1:N-M+1
                    if rem(i_c,2) == 1
                        Z(1:M,i_c) = temp_s1(i_c:i_c+M-1);
                    else
                        Z(M+1:2*M,i_c) = temp_s1(i_c:i_c+M-1);
                    end
                end
            end
            temp_value = temp_value + norm(Z'*Z-M*eye(N-M+1));
            i_index = i_index +1;
            temp_curve(i_index) = temp_value/num;
        end
        semilogy(seq_SNR,temp_curve);
        hold on
    end
    xlabel('SNR(dB)');ylabel('\epsilon');
    hold off;
    SetDrawStyle;
end