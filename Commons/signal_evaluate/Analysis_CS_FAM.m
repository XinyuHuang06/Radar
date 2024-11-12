function [CS, f_plot_grid, a_plot_grid] = Analysis_CS_FAM(varargin)
% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 10/11/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser;
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'signal', 0);
    addOptional(in_par, 'N1', 32);
    addParameter(in_par, 'bool_draw', 1);
    parse(in_par, varargin{:});
    fs = in_par.Results.fs;
    signal = in_par.Results.signal;
    N1 = in_par.Results.N1;
    bool_draw = in_par.Results.bool_draw;

    
    signal = signal(:); % 转为列向量
    N = length(signal); % 数据采样点数
    L = N1/4; % 重叠区长度
    P = pow2(nextpow2(N/L)); % 第二次FFT的大小
 
    N = (P-1)*L+N1; % 总观测长度
    if length(signal)<N          
        signal(N) = 0; 
    elseif length(signal)>N
        signal = signal(1:N);
    end
    
    seg_signal = zeros(N1,P); % 信号分段
    for iter_seg = 0:P-1
        seg_signal(:,iter_seg+1) = signal(iter_seg*L+1:iter_seg*L+N1);
    end
    
    a = hamming(N1); % 信号加窗
    seg_signal_win = diag(a)*seg_signal;

    seg_signal_fre = fft(seg_signal_win); % 傅里叶变换
    seg_signal_fre = fftshift(seg_signal_fre);
    seg_signal_fre = [seg_signal_fre(:,P/2+1:P),seg_signal_fre(:,1:P/2)];

    % complex demodulation. seg_signal_fre: N1*P
    cd_coef = zeros(N1,P); % 下变频矩阵
    for k = -N1/2:N1/2-1
        for r = 0:P-1
            cd_coef(k+N1/2+1,r+1) = exp(-1j*2*pi*k*r/P);
        end
    end
    seg_signal_cd = seg_signal_fre.*cd_coef;
    seg_signal_cd = transpose(seg_signal_cd); % P*N1
    seg_signal_cd2 = zeros(P,N1^2); % 复解调
    for k = 1:N1
        for l = 1:N1
            seg_signal_cd2(:,(k-1)*N1+l) = seg_signal_cd(:,k).*conj(seg_signal_cd(:,l));
        end
    end
    seg_signal_cd3 = fft(seg_signal_cd2);
    seg_signal_cd3 = fftshift(seg_signal_cd3);
    seg_signal_cd3 = [seg_signal_cd3(:,N1^2/2+1:N1^2),seg_signal_cd3(:,1:N1^2/2)];
    seg_signal_cd3 = abs(seg_signal_cd3);
    seg_signal_cd3 = seg_signal_cd3/max(seg_signal_cd3(:));
    % drawing figure
    f_plot_grid = -fs/2:fs/N1:fs/2;
    a_plot_grid = -fs:fs/N:fs;
    CS = zeros(N1+1,2*N+1); 
    for k1 = 1:P/2+1
        for k2 = 1:N1^2
            % 计算l当前的频率位置[-N1/2,N1/2-1];
            if rem(k2,N1) == 0
                l = N1/2 - 1;
            else
                l = rem(k2,N1) - N1/2 - 1;
            end
            % 计算k当前的频率位置[-N1/2,N1/2-1];
            k = ceil(k2/N1)-N1/2-1;
            f = (k+l)/(2*N1); % 归一化谱频率位置((k+l)/2)*(fs/N1)
            
            q = k1-P/4-1; % q \in [-P/4,P/4+1]
            alpha = (k-l)/N1+(q-1)/N; % 归一化循环频率位置((k+l))*(fs/N1) + q*(fs/N)

            if (alpha<-1 || alpha>1) || (f<-0.5 || f>0.5) % 如果不满足频率支撑域要求，则跳过
                continue;
            else
                kk = 1+N1*(f + 0.5);
                ll = 1+N*(alpha + 1);
                CS(round(kk), round(ll)) = seg_signal_cd3(k1,k2);
            end
        end
    end
    if bool_draw
        % mesh(seg_signal_cd3);
        mesh(a_plot_grid, f_plot_grid , CS);
        grid off;
        xlabel('Cycle frequency (Hz)');
        ylabel('Frequency (Hz)');
    end
end