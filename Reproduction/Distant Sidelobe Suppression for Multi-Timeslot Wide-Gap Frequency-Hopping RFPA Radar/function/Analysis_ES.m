% Example: Analysis_ES(signal, fs, 'normalized', 0)
% :param signal: 信号样本
% :param fs: 信号采样频率
% :return :
% detailed description: 归一化能量谱计算函数
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function Analysis_ES(varargin)
    in_par = inputParser;
    addOptional(in_par, 'signal', 0);
    addOptional(in_par, 'fs', 0);
    addOptional(in_par, 'N', 128);
    addParameter(in_par, 'normalized', 0);
    parse(in_par, varargin{:});
    signal = in_par.Results.signal;
    N = in_par.Results.N;
    fs = in_par.Results.fs;
    flag_normalized = in_par.Results.normalized;
    fres = fftshift(fft(signal,N))/N; % 频域表示
    energyspectrum = (abs(fres)).^2; % 能量谱
    if flag_normalized
        energyspectrum = energyspectrum/max(energyspectrum(:)); % 归一化
    end
    energyspectrum = 10*log10(energyspectrum); % 对数表示
    fs = -fs/2+fs/N:fs/N:fs/2;
    plot(fs, energyspectrum);
    ylim([-140,20]);
end

