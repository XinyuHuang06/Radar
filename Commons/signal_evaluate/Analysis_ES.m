% Example: Analysis_ES(signal, fs, N, 'normalized', 0);
% :param signal: The sample sequence
% :param fs: The sample frequency
% :param N: The point numbers of FFT
% :return energyspectrum:
% :return fs_seq: 
% detailed description: The normalized energy spectrum.                                                                        
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function varargout = Analysis_ES(signal, fs, N, varargin)
    in_par = inputParser;
    addParameter(in_par, 'bool_draw', 0);
    addParameter(in_par, 'normalized', 0);
    parse(in_par, varargin{:});
    flag_normalized = in_par.Results.normalized;
    bool_draw = in_par.Results.bool_draw;
    fres = fftshift(fft(signal,N)); % 频域表示
    energyspectrum = (abs(fres)*1/N).^2; % 能量谱
    if flag_normalized
        energyspectrum = energyspectrum/max(energyspectrum(:)); % 归一化
    end
    energyspectrum = 10*log10(energyspectrum); % 对数表示 1e-10为0修正量
    fs_seq = (-fs/2+fs/N:fs/N:fs/2)/1e6;
    if bool_draw
    plot(fs_seq, energyspectrum);
    end
    ylim([-140,20]);
    if nargout == 1
        varargout{1} = energyspectrum;
    elseif nargout == 2
        varargout{1} = energyspectrum;
        varargout{2} = fs_seq;
    end
end