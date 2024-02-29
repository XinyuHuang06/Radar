function [out_freseq,out_FreSpec_Ampli,out_FreSpec_Phase] = FS_Analysis(f_s,signal,boo_draw)
%% generate the result of  frequecny analysis
% :param f_s: the sampling frequency
% :param signal: the signal to be analyzed, a 1xN vector
% :return: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 09/06/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    N = length(signal);
    s_fredo = fft(signal);
    % Amplitude and Phase spectrum
    out_FreSpec_Ampli = fftshift(abs(s_fredo));
    out_FreSpec_Ampli = out_FreSpec_Ampli/max(out_FreSpec_Ampli);
    out_FreSpec_Phase = fftshift(angle(s_fredo));
    % frequency sequence
    out_freseq = -f_s+(2*f_s)/N:(2*f_s)/N:f_s;
    if boo_draw
        plot(out_freseq,out_FreSpec_Ampli);
        xlabel('Frequency/Hz');
        ylabel('Normalized amplitude');
        title('Fre Spectrum')
        SetDrawStyle;
    end
end